""" Module astropack.web:
    Web and internet utilities.

    Includes useful parsing utilities for certain astronomy data websites, notably:
    MPES, the Minor Planet Center's "minor planet ephemeris service".

    We have removed all downloads from minorplanet.info (e.g., OneAsteroidInfo pages),
    as for some reason known only to themselves, they now block all legitimate
    download attempts by all known packages.
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


__all__ = ['get_mp_info', 'get_mp_ephem']


__________FUNCTIONS___________________________________________________________ = 0


def get_mp_ephem(mp_id, utc_start, step_hours, num_entries, site=None):
    """Return MPES ephemeris information from one MPES web page
    (i.e., for one asteroid/minor planet).

    Parameters
    ----------
    mp_id : str, or int
        Identifier for the minor planet. Allowable formats include:
          * ``(3202)`` for a numbered MP
          * ``14442`` for a numbered MP, no parentheses (may be :class:`python:int`)
          * ``1997 XF11`` for an unnumbered MP
          * ``Badenia`` MP name.

    utc_start : |py.datetime|, |Time|, or a str that |Time| can parse
        Date and time to start. Will be truncated to most recent hour (i.e.,
        minutes and seconds are set to zero).

    step_hours : int
        Number of hours between sequential ephemeris entries in the result table.

    num_entries : int
        Number of ephemeris entries wanted.

    site : |Site| instance, or None, optional
        The earth location for which data are wanted, or
        None (default) for geocentric data. Default is None.

    Returns
    -------
    df_mp : |DataFrame|
        Pandas DataFrame containing MPES ephemeris data for minor planet ``mp_id``.
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
    table = MPC.get_ephemeris(target=target, location=location, start=start,
                              step=step, number=number, unc_links=False)
    df = table.to_pandas()
    columns_to_keep = [c for c in df.columns
                       if not c.startswith(('Uncertainty ', 'Unc. '))]
    df_mp = df.loc[:, columns_to_keep]
    return df_mp


def get_mp_info(mp_number=None, mp_name=None):
    """Given either a Minor Planet number or name (not both),
    return a python dictionary of principal data as retrieved from the Minor Planet
    Center database.

    Parameters
    ----------
    mp_number : int, or str representing an int
        Number of target minor planet, e.g., 233 or '233'
    mp_name : str
        Name of target minor planet, e.g., 'Badenia'. This may *not* be a designation,
        as the MPC database does not currently (March 2022) accept them.

    Returns
    -------
    info_dict : dict, or None if no ``mp_number`` or ``mp_name`` given
        Dictionary of selected MP data, including orbital period, G, and H.

        Current (version 0.1beta) dict keys are:
            * name: MP name
            * number: MP number
            * designation: MPC designation
            * period: orbital period, in years
            * H: reduced magnitude, V band
            * G: phase slope, V band
    """
    query_result = None  # keep IDE happy.
    if mp_number is not None:
        if isinstance(mp_number, int):
            query_result = MPC.query_object(target_type='asteroid', number=mp_number)
        elif isinstance(mp_number, str):
            try:
                _ = int(mp_number)
            except ValueError:
                raise TypeError('Parameter \'mp_number\' must be None, an integer,'
                                ' or a string representing an integer.')
        query_result = MPC.query_object(target_type='asteroid', number=str(mp_number))
    elif mp_name is not None:
        if isinstance(mp_name, str):
            query_result = MPC.query_object(target_type='asteroid', name=mp_name)
        else:
            raise TypeError('Parameter \'mp_name\' must be a string.')
    else:
        raise ValueError('No MP number, name, or designation given.')
    if query_result is None:
        return None

    if isinstance(query_result, list):
        if len(query_result) == 0:
            return None
    d = query_result[0]
    kept_keys = ['name', 'number', 'designation', 'period']
    new_keys = {'absolute_magnitude': 'H', 'phase_slope': 'G'}
    info_dict = {k: d[k] for k in kept_keys}
    for old_key in new_keys.keys():
        info_dict[new_keys[old_key]] = d[old_key]
    return info_dict
