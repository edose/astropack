""" Module astropack.almanac:
    Astronight class + other almanac code. (py file)
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timezone, timedelta
from math import sin, cos, asin, acos, sqrt

# External packages:
from numpy import diff, flatnonzero, clip, sign
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from skyfield.api import load, wgs84, Star
from skyfield.almanac import risings_and_settings, fraction_illuminated, \
    meridian_transits
from skyfield.searchlib import find_discrete

# Astropack packages:
from .util import Timespan, ra_as_hours
from .reference import DEGREES_PER_RADIAN

THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')

HORIZON_USNO = -0.833  # USNO sun/moon effective horizon (radius & refraction), degrees.
HOURS_TO_ASSURE_HALF_NIGHT = 14

__all__ = ['Astronight',
           'NoDarkTimeError',
           'SunAlwaysUpError',
           'calc_timespan_no_sun',
           'find_target_up_down',
           'ha_dec_from_az_alt',
           'local_sidereal_time',
           'moon_illumination_pct',
           'moon_ra_dec',
           'moon_transit_time',
           'target_transit_time',
           'astropy_time_to_skyfield',
           'skyfield_time_to_astropy',
           'make_skyfield_observatory_from_site']


class SunAlwaysUpError(Exception):
    """Raised when new |Astronight| has no sun-down time for observing.

    User is responsible for catching this, if observing site's latitude allows for
    24 consecutive hours of sun-up daylight (only in earth's polar latitudes of about
    > +66 or <-66 degrees).
    """
    pass


class NoDarkTimeError(Exception):
    """Raised when new |Astronight| has no dark time for observing.

    User is responsible for catching this, if observing site's latitude allows for
    24 consecutive hours of daylight or twilight (typically only at earth's polar
    latitudes of about > +55 or <-55 degrees.
    """
    pass


__________CLASS_ASTRONIGHT_______________________________________________________ = 0


class Astronight:
    """**Astropack**'s main engine of astronomical planning.

    One need supply only a |Site| instance and an astronight date (yyyymmdd,
    representing the local date when nighttime starts). Astronight's constructor
    generates and supplies numerous properties and methods useful to planning
    observations at the location during that night.

    Parameters
    ----------
    site : |Site|
        Observing location.

    an_date : int of form yyyymmdd, or string of form 'yyyymmdd'
        The Astronight date.
        Formally, this is the local date (yyyymmdd) of the moment before midnight of
        the observing nighttime.

        This should always be the local date
        (format yyyymmdd) when the sun sets and observations can begin; it is
        equivalent to ACP's (Astronomer's Control Panel) definition of observing date.

    sun_altitude_dark : float, optional
        Altitude of sun's center, in degrees, when has just become sufficiently dark to
        begin observations, i.e., the end of twilight.
        Default is value of ``site``'s ``sun_altitude_dark`` attribute.

    Attributes
    ----------
    site : |Site|
        Observer's location info, from input parameter ``site``.
    site_name : str
        Name of the observer's location, from input parameter ``site``.
    an_date_string : str
        String representation of the astronight date, of form 'yyyymmdd'
    sun_altitude_dark : float
        Maximum sun altitude, in degrees, that is considered observably dark.
    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition object, corresponding to input parameter ``site``.

        See: |skyfield.GP|
    timespan_no_sun : |Timespan|
        Timespan when sun is down for this astronight, at this site.
    timespan_dark : |Timespan|
        Timespan when sky is observably dark for this astronight, at this site.
    utc_local_middark : |Time|
        Midpoint of ``timespan_dark``.
        More useful than clock midnight for almanac calculations, as it is practically
        guaranteed to be dark, unless that astronight has no dark time at all
        (which raises :class:`~.almanac.SunAlwaysUpError` exception).
    local_middark_lst_hour_string : str
        Local sidereal time as sexigesimal string of form 'hh:mm:ss'.
    moon_illumination : float
        Percent illumination of moon at local middark, in range [0-100].
    moon_ra : float
        Moon's Right Ascension at local middark, in degrees, in range [0, 360).
    moon_dec : float
        Moon's Declination at local middark, in degrees, in range [-90, +90].
    moon_skycoord : |SkyCoord|
        Moon's sky location at local middark.
    moon_transit : |Time|
        The time nearest local middark at which the moon crosses local meridian.
    moon_up_timespans : |Timespan|, or list of two |Timespan|
        Timespan when moon is above horizon, as limited by ``timespan_no_sun``.
        List of 2 |Timespan| if, for this astronight, the moon is above horizon at
        both sunset and sunrise and below horizon at some time between (rare).
    timespan_dark_no_moon : |Timespan|, or list of two |Timespan|
        Timespan when sun is below its user-chosen dark elevation, and moon is below
        horizon. The best observing timespan(s). May be a list of two |Timespan| if
        moon is up in middle of the night (very rare).

    Raises
    ------
    SunAlwaysUpError
        Raised when |Astronight| at this site has no sun-down time for observing.

    NoDarkTimeError
        Raised when |Astronight| at this site has no dark time for observing.

    Examples
    --------
    Note that the format for astronight date differs from ISO 8601 format.
    This is meant to emphasize the difference in meaning between the two date types.

      >>> from astropack.almanac import Astronight
      >>> from astropack.ini import Site
      >>> this_site = Site('My Dome')
      >>> this_date = 20220302  # or '20220302'
      >>> an = Astronight(this_site, this_date)

    The Astronight object then supplies relevant data about the local night.

      >>> mp = an.moon_illumination
      >>> dark = an.timespan_dark
      >>> very_dark = an.timespan_dark_no_moon

    The Astronight object can also supply target-observation data, typically taking a
    sky location and returning a |Time| or |Timespan| object.

      >>> from astropy.coordinates import SkyCoord
      >>> target = SkyCoord(10.4532, -25.32, unit="deg")
      >>> transit_time = an.transit(target)
      >>> observable = an.timespan_observable(target, min_alt=30, min_moon_dist=50)

    If a given site and date will have no nighttime or dark time (polar latitudes,
    in summer), Astronight will raise an exception:
    :class:`~.almanac.SunAlwaysUpError` or
    :class:`~.almanac.NoDarkTimeError`.
    Users planning for observations at polar latitudes are
    encouraged to handle these with try/catch blocks.
    """

    def __init__(self, site, an_date, sun_altitude_dark=None):
        # Handle inputs:
        self.site = site
        self.site_name = site.name
        self.an_date_string = str(an_date)
        if len(self.an_date_string) != 8:
            raise ValueError('Parameter \'.an_date_string\' appears '
                             'invalid (must be yyyymmdd).')
        if sun_altitude_dark is None:
            self.sun_altitude_dark = site.sun_altitude_dark
        else:
            self.sun_altitude_dark = sun_altitude_dark

        # Initiate skyfield package:
        self.master_eph = load('de440s.bsp')
        self.timescale = load.timescale()
        self.obs = make_skyfield_observatory_from_site(site)  # Skyfield observer.

        # Get utc_approx_midnight from an_date and longitude, to use as starting point:
        an_year = int(self.an_date_string[0:4])
        an_month = int(self.an_date_string[4:6])
        an_day = int(self.an_date_string[6:8])
        if site.utc_offset < 0:
            # West of Greenwich but east of Intl Date Line [Americas].
            longitude_hours = ((site.longitude % 360) - 360.0) / 15.0
        else:
            # East of Greenwich but west of Intl Date Line [Eurasia, Africa, ANZ].
            longitude_hours = (site.longitude % 360) / 15.0
        utc_approx_midnight = Time(datetime(an_year, an_month, an_day, 0, 0, 0),
                                   scale='utc') + \
            TimeDelta(-longitude_hours * 3600, format='sec') + \
            TimeDelta(24 * 3600, format='sec')

        # Refine utc_approx_midnight to utc_sun_lowest:
        # Nested function:
        def d_sun_alt_d_time(time):
            """Nested function to estimate first derivative of sun altitude vs time."""
            d_time = 4.0  # seconds.
            topos_at = (self.master_eph['earth'] + self.obs).at
            time_before = time - timedelta(seconds=d_time / 2)
            time_after = time + timedelta(seconds=d_time / 2)
            sun_alt_before = topos_at(time_before).observe(self.master_eph['sun']).\
                apparent().altaz()[0].degrees
            sun_alt_after = topos_at(time_after).observe(self.master_eph['sun']).\
                apparent().altaz()[0].degrees
            return sign((sun_alt_after - sun_alt_before) / d_time)

        d_sun_alt_d_time.rough_period = 0.25  # skyfield's VERY weird API requires this.
        start_time = astropy_time_to_skyfield(self.timescale, utc_approx_midnight) - \
                     timedelta(hours=3)
        end_time = astropy_time_to_skyfield(self.timescale, utc_approx_midnight) + \
                   timedelta(hours=3)
        event_times, _ = find_discrete(start_time, end_time, d_sun_alt_d_time)
        if len(event_times) >= 1:
            utc_sun_lowest = skyfield_time_to_astropy(event_times[0])
            sf_sun_lowest = event_times[0]
        else:
            utc_sun_lowest = utc_approx_midnight
            sf_sun_lowest = self.timescale.from_datetime(utc_approx_midnight.
                                                         to_datetime().
                                                         replace(tzinfo=timezone.utc))
        topos_at = (self.master_eph['earth'] + self.obs).at
        alt_sun_lowest = topos_at(sf_sun_lowest).observe(self.master_eph['sun']) \
            .apparent().altaz()[0].degrees
        if alt_sun_lowest > HORIZON_USNO:
            raise SunAlwaysUpError('Astronight date ' + self.an_date_string +
                                   ', latitude=' + '{0:.2f}'.format(site.latitude))
        if alt_sun_lowest > self.sun_altitude_dark:
            raise NoDarkTimeError('Astronight date ' + self.an_date_string +
                                  ', latitude=' + '{0:.2f}'.format(site.latitude))

        # No-sun (sun below horizon) timespan:
        self.timespan_no_sun = calc_timespan_no_sun(self.obs, self.master_eph,
                                                    self.timescale, utc_sun_lowest)
        # if self.timespan_no_sun is None:
        #     raise SunAlwaysUpError('Astronight date ' + self.an_date_string +
        #                            ', latitude=' + '{0:.2f}'.format(site.latitude))

        # Dark (sun far below horizon) timespan:
        self.timespan_dark = find_target_up_down(self.obs, self.master_eph,
                                                 self.master_eph['sun'],
                                                 self.timescale,
                                                 self.timespan_no_sun, 'down',
                                                 horizon=self.sun_altitude_dark)[0]
        # if self.timespan_dark is None:
        #     raise NoDarkTimeError('Astronight date ' + self.an_date_string +
        #                           ', latitude=' + '{0:.2f}'.format(site.latitude))
        self.utc_local_middark = self.timespan_dark.midpoint
        local_middark_lst_degrees = 15.0 * local_sidereal_time(site.longitude,
                                                               self.timescale,
                                                               self.utc_local_middark)
        self.local_middark_lst_hour_string = ra_as_hours(local_middark_lst_degrees)

        # Moon quantities and moon_down timespan:
        self.moon_illumination = moon_illumination_pct(self.master_eph, self.timescale,
                                                       self.utc_local_middark)
        self.moon_ra, self.moon_dec = moon_ra_dec(self.obs, self.master_eph,
                                                  self.timescale,
                                                  self.utc_local_middark)
        self.moon_skycoord = SkyCoord(self.moon_ra, self.moon_dec, unit="deg")
        self.moon_transit = moon_transit_time(self.obs, self.master_eph,
                                              self.timescale, self.utc_local_middark)
        self.moon_up_timespans = find_target_up_down(self.obs, self.master_eph,
                                                     self.master_eph['moon'],
                                                     self.timescale,
                                                     self.timespan_dark, 'up',
                                                     HORIZON_USNO)

        # Dark_no_moon timespan:
        dnm = self.timespan_dark.copy()
        for mdt in self.moon_up_timespans:
            dnm = dnm.subtract(mdt)
        self.timespan_dark_no_moon = dnm

    def transit(self, target_skycoord):
        """Return transit UTC datetime, closest to middark time, for given target.

        Parameters
        ----------
        target_skycoord : |SkyCoord|
            Sky position of the target to observe.

        Returns
        -------
        transit_time : |Time|
            Target's transit time (local meridian crossing) time at the astronight
            ``site``.
        """
        return target_transit_time(self.obs, self.master_eph, self.timescale,
                                   target_skycoord, self.utc_local_middark)

    # TODO: ensure that this works for a LIST of SkyCoords or array SkyCoords,
    #     as well as for a single SkyCoord.
    def timespan_observable(self, target_skycoord, min_alt=0.0, min_moon_dist=None):
        """ Return observable timespan for given target.

        Parameters
        ----------
        target_skycoord : |SkyCoord|
            Sky position of the target to observe.

        min_alt : float, or |Angle|
            Minimum altitude for target to be observable. In degrees, if a float.
            Default is zero (ideal horizon).

        min_moon_dist : float, or |Angle|
            Minimum distance from the moon for the target is observable (but only
            when the moon is up. In degrees, if a float. Set to zero to disable.
            Required.

        Returns
        -------
        timespan_observable : |Timespan|
            Timespan during which target is observable.

        Raises
        ------
        ValueError
            Raised when ``min_moon_dist`` is unspecified.
        """
        if min_moon_dist is None:
            raise ValueError('Astropy.timespan_observable requires value '
                             'for min_moon_dist.')
        # moon = self.master_eph['moon']
        # middark_sf = self.timescale.from_astropy(self.utc_local_middark)

        star_list = _skycoords_to_skyfield_star_list(target_skycoord)
        timespan_target_up = find_target_up_down(self.obs, self.master_eph, star_list,
                                                 self.timescale, self.timespan_dark,
                                                 'up', min_alt)
        moon_dist = self.moon_ra.separation_from(star_list)
        if moon_dist >= min_moon_dist:
            return timespan_target_up
        else:
            return [time_up.intersection(self.timespan_dark_no_moon)
                    for time_up in timespan_target_up]

    def __repr__(self):
        return 'astropack.almanac.Astronight(site=[' + self.site.name + ']'\
               ', an_date=' + self.an_date_string + \
               ', sun_altitude_dark=' + str(self.sun_altitude_dark) + ')'

    def __str__(self):
        return "Astronight '" + self.an_date_string + "' at site '" + \
               self.site_name + "'."


__________PUBLIC_FUNCTIONS______________________________________________________ = 0


def ha_dec_from_az_alt(latitude, az_alt):
    """Return hour angle and declination for a given azimuth and altitude,
    at a site's given latitude.

     Parameters
     ----------
     latitude : float
        Observing site latitude, in degrees.

     az_alt : tuple of 2 float, or list of tuple
        Tuple of (azimuth, altitude), in degrees.

    Returns
    -------
    ha_dec : tuple of 2 float, or list of tuple
       Tuple of (hourangle, declination), in degrees.
    """
    az_alt_is_list = isinstance(az_alt, list)
    if not az_alt_is_list:
        az_alt = [az_alt]

    # Formulae (1) and (2) from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm:
    # phi=latitude, delta=declination, H=hourangle, A=target azimuth, a=target altitude.
    cos_phi, sin_phi = cos(latitude / DEGREES_PER_RADIAN), \
                       sin(latitude / DEGREES_PER_RADIAN)
    results = []
    for (az, alt) in az_alt:
        cos_a, sin_a = cos(alt / DEGREES_PER_RADIAN), sin(alt / DEGREES_PER_RADIAN)
        cos_A = cos(az / DEGREES_PER_RADIAN)

        # (1) sin(δ) = sin(a) sin(φ) + cos(a) cos(φ) cos(A)
        sin_delta = (sin_a * sin_phi) + (cos_a * cos_phi * cos_A)
        cos_delta = sqrt(1.0 - sin_delta ** 2)  # cosine of declination >= 0 always.
        # Clip to correct any floating-point errors driving abs(sin(delta)) > 1:
        sin_delta = clip(sin_delta, -1.0, 1.0)
        delta = asin(sin_delta) * DEGREES_PER_RADIAN

        # (2) cos(H) = { sin(a) - sin(δ) sin(φ) } / { cos(δ) cos(φ) }
        cos_H = (sin_a - sin_delta * sin_phi) / (cos_delta * cos_phi)
        # Clip to correct any floating-point errors driving abs(cos(H)) > 1:
        cos_H = clip(cos_H, -1.0, 1.0)
        sign_H = 1.0 if 0.0 <= az <= 180.0 else -1.0  # because sign of H is ambiguous.
        H = sign_H * acos(cos_H) * DEGREES_PER_RADIAN
        results.append((H, delta))
    if not az_alt_is_list:
        results = results[0]
    return results


__________SUPPORT_FUNCTIONS______________________________________________________ = 0


def calc_timespan_no_sun(obs, master_eph, timescale, approx_midnight):
    """Return |Timespan| representing sunset to sunrise and containing
    ``approx_midnight``.

    Parameters
    ----------
    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition instance.

        See: |skyfield.GP|

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    timescale:  skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    approx_midnight : |Time|
        Approximate local midnight, UTC.

    Returns
    -------
    timespan_no_sun : |Timespan|, or None
        Timespan from sunset to sunrise.
        None if there is no time near ``approx_midnight`` when sun is down
        (polar summer).
    """
    midnight = timescale.from_astropy(approx_midnight)  # skyfield Time obj.
    time_before = midnight - timedelta(hours=HOURS_TO_ASSURE_HALF_NIGHT)
    time_after = midnight + timedelta(hours=HOURS_TO_ASSURE_HALF_NIGHT)
    sun_fn = risings_and_settings(master_eph, master_eph['sun'], obs, HORIZON_USNO)
    sun_fn.step_days = 1.0 / 12.0
    set_times, set_values = find_discrete(time_before, midnight, sun_fn)
    rise_times, rise_values = find_discrete(midnight, time_after, sun_fn)
    topos_at = (master_eph['earth'] + obs).at
    alt_midnight = topos_at(midnight).observe(master_eph['sun'])\
        .apparent().altaz()[0].degrees

    # Case: no rises or sets, sun is either up or down for the entire search period.
    if len(set_values) <= 0 and len(rise_values) <= 0:
        if alt_midnight <= HORIZON_USNO:
            return Timespan(skyfield_time_to_astropy(time_before),
                            skyfield_time_to_astropy(time_after))
        else:
            return None
    if len(set_values) >= 1:
        set_time = skyfield_time_to_astropy(set_times[-1])
    else:
        set_time = skyfield_time_to_astropy(time_before)
    if len(rise_values) >= 1:
        rise_time = skyfield_time_to_astropy(rise_times[0])
    else:
        rise_time = skyfield_time_to_astropy(time_after)
    return Timespan(set_time, rise_time)


def find_target_up_down(obs, master_eph, target_eph, timescale, timespan,
                        up_down, horizon=0.0):
    """Returns |Timespan| for which target is either up or down
    (above or below given horizon), within a given timespan.

    Parameters
    ----------
    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition object.

        See: |skyfield.GP|

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    target_eph : skyfield.vectorlib.VectorSum or skyfield.starlib.Star
        Target ephemeris; either a solar-system ephemeris (selected from
        ``master_eph``) or a skyfield Star object.

        See:
        https://github.com/skyfielders/python-skyfield/blob/master/skyfield/vectorlib.py
        or https://rhodesmill.org/skyfield/api-stars.html#skyfield.starlib.Star

    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    timespan : |Timespan|
        Timespan within which result will be constrained.
        If ``up_down`` is 'down', one convenient source for ``timespan``
        is ``.util.Astronight.timespan_dark``.

    up_down : {'up', 'down'}
        If 'up', results represent when target is above horizon.
        If 'down', results represent when target is below horizon. Required.

    horizon : float, optional
        Altitude above actual horizon that separates 'up' from 'down', in degrees.
        Default is zero (ideal horizon).

    Returns
    -------
    timespan_up_down : list of |Timespan|, or None.
        List of timespans for which target is either up or down as selected.
        Returns a list even if only one Timespan is returned.
        If ``up_down`` is satisfied throughout input timespan, the returned list will
        contain a copy of ``timespan``.
        If ``up_down`` is never satisfied during input timespan, return None.
    """
    up_down = up_down.lower().strip()
    if up_down not in ['up', 'down']:
        raise ValueError('Parameter up_down must equal \'up\' or \'down\'.')
    start = timescale.from_astropy(timespan.start)
    end = timescale.from_astropy(timespan.end)
    dark_fn = risings_and_settings(master_eph, target_eph, obs, horizon)
    event_times, event_values = find_discrete(start, end, dark_fn)
    # Case: no risings or settings found:
    if len(event_values) <= 0:
        topos_at = (master_eph['earth'] + obs).at
        mid_time = timescale.from_astropy(timespan.midpoint)
        alt = topos_at(mid_time).observe(target_eph).apparent().altaz()[0].degrees
        if (up_down == 'up' and alt > horizon) or \
            (up_down == 'down' and alt <= horizon):
            return [timespan.copy()]
        else:
            return None
    # Case: at least one rising or setting found:
    # First, convert times to astropy, so that start and end are passed through exactly.
    astropy_times = [timespan.start] + \
                    [skyfield_time_to_astropy(t) for t in event_times] + \
                    [timespan.end]
    values = [0.5] + list(event_values) + [0.5]
    diffs = diff(values)
    # Now choose events as up or down:
    if up_down == 'up':
        is_index_qualifying = diffs < 0
    else:
        is_index_qualifying = diffs > 0
    qualifying_indices = flatnonzero(diffs)[is_index_qualifying]
    timespans = [Timespan(astropy_times[i], astropy_times[i + 1])
                 for i in qualifying_indices]
    return timespans


def moon_ra_dec(obs, master_eph, timescale, time):
    """Returns moon's sky coordinates at given earth location and time.

    Parameters
    ----------
    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition object.

        See: |skyfield.GP|

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    time : |Time|
        Time UTC at which moon's sky location is wanted.

    Returns
    -------
    ra, dec : tuple of float
        (RA, Declination), in degrees, of moon at observer location and given time.
    """
    time_sf = timescale.from_astropy(time)
    ra, dec, _ = (master_eph['earth'] + obs).at(time_sf).\
        observe(master_eph['moon']).apparent().radec()
    return ra.hours * 15.0, dec.degrees


def moon_transit_time(obs, master_eph, timescale, time):
    """Return moon's transit time closest to given time.

    Parameters
    ----------
    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition object.

        See: |skyfield.GP|

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    time : |Time|
        Time UTC from which moon's nearest transit time is wanted.

    Returns
    -------
    transit_time : |Time|
        Moon's transit time UTC at observer's location and nearest to given time.
    """
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    time_start = time_sf - timedelta(hours=14)
    time_end = time_sf + timedelta(hours=14)
    moon_fn = meridian_transits(master_eph, master_eph['moon'], obs)
    times, values = find_discrete(time_start, time_end, moon_fn)
    return skyfield_time_to_astropy(_closest_transit(times, values, time_sf))


def target_transit_time(obs, master_eph, timescale, target_skycoord, time):
    """Return fixed target's transit time closest to a given time.

    Parameters
    ----------

    obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition object.

        See: |skyfield.GP|

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    target_skycoord : |SkyCoord|
        Target's fixed sky position.

    time : |Time|
        Time UTC nearest which target's nearest transit time is wanted.

    Returns
    -------
    transit_time : |Time|
        Moon's transit time UTC at observer's location and nearest to given time
    """
    star_list = _skycoords_to_skyfield_star_list(target_skycoord)
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    time_start = time_sf - timedelta(hours=14)
    time_end = time_sf + timedelta(hours=14)
    transit_list = []
    for star in star_list:
        star_fn = meridian_transits(master_eph, star, obs)
        times, values = find_discrete(time_start, time_end, star_fn)
        transit_list.append(skyfield_time_to_astropy(_closest_transit(times,
                                                                      values,
                                                                      time_sf)))
    if len(transit_list) == 1:
        return transit_list[0]  # return single time as scalar Time rather than list.
    return transit_list


def _closest_transit(times, values, time_sf):
    """For times [list of skyfield Time objects], values [list of ints], and
    a target time [Skyfield Time object], determine which transit is closest
    to target time, return the time as a Skyfield Time object."""
    transit_times = [t for (t, v) in zip(times, values) if v == 1]
    timedeltas = [abs(time - time_sf) for time in transit_times]
    transit_time = transit_times[timedeltas.index(min(timedeltas))]
    return transit_time


def local_sidereal_time(longitude, timescale, time):
    """Return local sidereal time, in hours, for a given time and earth longitude.

    Parameters
    ----------
    longitude : float
        Longitude of observer's earth location, in degrees, in range [-180, +180].

    timescale:  skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    time : |Time|
        Time UTC for which local sidereal time is wanted.

    Returns
    -------
    lst : float
        Local sidereal time at observer's longitude and time, in hours
    """
    from skyfield.earthlib import sidereal_time
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    gmst = sidereal_time(time_sf)  # in hours
    local_offset = longitude / 15.0  # in hours
    return gmst + local_offset


def moon_illumination_pct(master_eph, timescale, time):
    """Return moon illumination extent, as percentage, at a given time.

    Parameters
    ----------

    master_eph :  skyfield.jpllib.SpiceKernel
        Skyfield master ephemeris, available as Astronight.skyfield_eph

        See: |skyfield.ME|

    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    time : |Time|
        Time UTC for which moon's illumination is wanted.

    Returns
    -------
    moon_pct : float
        Moon's percent illumination at given time, in range [0, 100].
    """
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    fraction = 100.0 * fraction_illuminated(master_eph, 'moon', time_sf)
    return fraction


def make_skyfield_observatory_from_site(site):
    """ Return Skyfield observer (earth location) object made
    from astropack |Site| instance.

    Parameters
    ----------
    site : |Site|
        Earth location, as read from .ini file by |Site| class.

    Returns
    -------
    skyfield_obs : skyfield.toposlib.GeographicPosition
        Skyfield GeographicPosition instance.

        See: |skyfield.GP|
    """
    obs = wgs84.latlon(site.latitude, site.longitude, site.elevation)
    return obs


def astropy_time_to_skyfield(timescale, time):
    """Convert astropy time in UTC to skyfield time.

.. Warning:: Unfortunately, skyfield and astropy both use 'Time' as name
        for their main time class. Because you'll probably use astropy |Time| more
        often, if you also need skyfield Time, consider importing them both as,
        for example,

        >>> from astropy.time import Time
        >>> import skyfield.timelib.Time as sfTime

    Parameters
    ----------
    timescale : skyfield.timelib.Timescale
        Skyfield timescale object, available as Astronight.skyfield.ts

        See: |skyfield.TS|

    time : |Time|
        Astropy time instance, UTC

    Returns
    -------
    skyfield_time : skyfield.timelib.Time
        Skyfield time instance.

        See: |skyfield.Time|
    """
    return timescale.from_datetime(time.to_datetime().replace(tzinfo=timezone.utc))


def skyfield_time_to_astropy(skyfield_time):
    """Convert skyfield time to astropy time in UTC timescale and ISO format.
    Time arrays are OK.

.. Warning:: Unfortunately, skyfield and astropy both use 'Time' as name
        for their main time class. Because you'll probably use astropy |Time| more
        often, if you also need skyfield Time, consider importing it as, for example,

        >>> from astropy.time import Time
        >>> import skyfield.timelib.Time as sfTime

    Parameters
    ----------
    skyfield_time : skyfield.timelib.Time
        Skyfield time instance.

        See: |skyfield.Time|

    Returns
    -------
    astropy_time : |Time|
        Astropy time instance, UTC
    """
    astropy_time = skyfield_time.to_astropy().utc
    astropy_time.format = 'iso'
    return astropy_time


def _skycoords_to_skyfield_star_list(skycoord):
    """Convert astropy SkyCoord object or list of objects to a *list of* skyfield
    Star objects. Always returns a list, only one sky position per Star.
    """
    def stars_from_skycoord_object(skycoord_obj):
        """Nested function, always returns a list of Star objects."""
        if skycoord_obj.size == 1:
            return [Star(ra_hours=skycoord_obj.ra.degree / 15.0,
                         dec_degrees=skycoord_obj.dec.degree)]
        else:
            return [Star(ra_hours=sc_element.ra.degree / 15.0,
                         dec_degrees=sc_element.dec.degree)
                    for sc_element in skycoord_obj]

    if isinstance(skycoord, SkyCoord):
        return stars_from_skycoord_object(skycoord)
    elif isinstance(skycoord, list):
        star_list = []
        for sc_object in skycoord:
            star_list.extend(stars_from_skycoord_object(sc_object))
        return star_list
