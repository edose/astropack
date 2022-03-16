""" Module astropak.web:
    Web and internet utilities.
    Includes useful parsing utilities for certain astronomy data websites, notably:
        MPES, the Minor Planet Center's "minor planet ephemeris service".
    Removed facilities to download from minorplanet.info (e.g., OneAsteroidInfo pages),
        as for some reason known only to themselves, they now block all legitimate
        download attempts by all known packages. Bad show.
"""

__author__ = "Eric Dose, Albuquerque"

# Python core:
from datetime import datetime, timezone

# External packages:
# import pandas as pd
import astropy.units as u
from astropy.time import Time
# from astropy.coordinates import Latitude, Longitude
from astroquery.mpc import MPC

# From other modules, this package:
# from astropak.ini import Site


__________FUNCTIONS___________________________________________________________ = 0


def get_mp_ephem(mp_id, utc_start, step_hours, num_entries, site=None):
    """Return pandas DataFrame containing MPES information from one MPES web page.
       Step between entries is hourly.
    :param mp_id: MP identifier for which to query. [string, or int for numbered MP]
    Allowable formats:
        (3202)     numbered MP
        14442      numbered MP
        1997 XF11  unnumbered MP
        Badenia    MP name
    :param utc_start: Date and time to start. [py datetime, astropy Time object, or
    string convertible to astropy Time object]
    Will be truncated to most recent hour, zero minutes and seconds.
    If py datetime object in UTC (any timezone information will be overwritten, not converted).
    :param step_hours: hours between sequential entries. [int]
    :param num_entries: number of entries to be returned in dataframe. [int]
    :param site: Site object, or None for geocentric. [astropak Site object]
    :return: dataframe of MPES data for requested minor planet. [pandas DataFrame]
    Columns:

    """
    target = str(mp_id)
    if isinstance(utc_start, datetime):
        dt_as_utc = utc_start.replace(tzinfo=timezone.utc)
        start = Time(dt_as_utc)
    elif isinstance(utc_start, str):
        start = Time(utc_start, scale='utc')
    elif isinstance(utc_start, Time):
        start = utc_start
    else:
        raise TypeError('Parameter \'utc_start\' must be a py datetime, astropy Time, '
                        'or string convertible to astropy Time.')
    step = step_hours * u.hour
    number = num_entries
    if site is None:
        location = None
    else:
        location = site.mpc_code
    table = MPC.get_ephemeris(target=target, location=location, start=start, step=step, number=number,
                              unc_links=False)
    df = table.to_pandas()
    columns_to_keep = [c for c in df.columns if not c.startswith(('Uncertainty ', 'Unc. '))]
    df_mp = df.loc[:, columns_to_keep]
    return df_mp


def get_mp_info(mp_number=None, mp_name=None):
    """For one MP, retrieve selected MP info from MPC database.
    MP designation is not allowed, as MPC gives no data when designation passed in.
    :param mp_number: e.g., 333. [int]
    :param mp_name: e.g., 'Badenia'. Used if mp_number is not given. [string]
    :return: dictionary of selected MP data. [py dict]
    """
    result = None  # keep IDE happy.
    if mp_number is not None:
        if not isinstance(mp_number, int):
            raise TypeError('Parameter \'mp_number\' must be an integer.')
        result = MPC.query_object(target_type='asteroid', number=mp_number)
    elif mp_name is not None:
        if not isinstance(mp_name, str):
            raise TypeError('Parameter \'mp_name\' must be a string.')
        result = MPC.query_object(target_type='asteroid', name=mp_name)
    else:
        raise ValueError('No MP number, name, or designation given.')

    if result is None:
        return None
    if isinstance(result, list):
        if len(result) == 0:
            return None
    d = result[0]
    kept_keys = ['name', 'number', 'designation', 'period']
    new_keys = {'absolute_magnitude': 'H', 'phase_slope': 'G'}
    info_dict = {k: d[k] for k in kept_keys}
    for old_key in new_keys.keys():
        info_dict[new_keys[old_key]] = d[old_key]
    return info_dict
