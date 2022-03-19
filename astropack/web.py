""" Module astropack.web:
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


__all__ = ['get_mp_ephem',
           'get_mp_info']


__________FUNCTIONS___________________________________________________________ = 0


def get_mp_ephem(mp_id, utc_start, step_hours, num_entries, site=None):
    """Return pandas DataFrame containing MPES information from one MPES web page.

    Parameters
    ----------
    mp_id : int for numbered Minor Planet, or str
        Identifier for the minor planet.
        Allowable formats:
            (3202)     numbered MP
            14442      numbered MP
            1997 XF11  unnumbered MP
            Badenia    MP name

    utc_start : py datetime, |Time|, or str the |Time| can parse
        Date and time to start. Will be truncated to most recent hour (i.e.,
        minutes and seconds are set to zero).

    step_hours : int
        Number of hours between sequential entries in the result table.

    num_entries : int
        Number of entries wanted.

    site : `~.ini.Site` instance, or None, optional
        The earth location for which data are wanted, or
        None (default) for geocentric data.

    Returns
    -------
    df_mp : |DataFrame|
        Pandas DataFrame containing MPES data for minor planet `mp_id`.
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
    """Given either a Minor Planet number or name,
    return a python dictionary of relevant data retrieved from the Minor Planet
    Center database.

    Parameters
    ----------
    mp_number : int
        Number of target minor planet.
    mp_name : str
        Name of target minor planet, e.g., 'Badenia'. This may not be a designation,
        as the MPC database does not currently (March 2022) accept them.

    Returns
    -------
    info_dict : dict
        Dictionary of selected MP data, including orbital period, G, and H.
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
