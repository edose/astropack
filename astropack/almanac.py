""" Module astropack.almanac:
    Astronight class + other almanac code. (py file)
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timezone, timedelta
from math import sin, cos, asin, acos, sqrt
from abc import ABC, abstractmethod
from typing import TypeAlias, List, Tuple, Callable
from enum import Enum, auto

# External packages:
from numpy import diff, flatnonzero, clip, sign
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, Angle
import skyfield.vectorlib
import skyfield.starlib
import skyfield.timelib
import skyfield.toposlib
from skyfield.api import load, wgs84, Star
from skyfield.almanac import risings_and_settings, fraction_illuminated, \
    meridian_transits
from skyfield.searchlib import find_discrete

# Astropack packages:
from .util import Timespan, ra_as_hours
from .reference import DEGREES_PER_RADIAN
from .ini import Site

THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')

HORIZON_USNO = -0.833  # USNO sun/moon effective horizon (radius & refraction), degrees.
HOURS_TO_ASSURE_HALF_NIGHT = 14

__all__ = ['Astronight',
           'NoDarkTimeError',
           'SunAlwaysUpError',
           ]

AN_date_type: TypeAlias = 'AN_date'


class Invalid_ANDate_Error(Exception):
    """Raised when AN (Astronight) input cannot be converted to a valid and
    reasonable calendar date."""
    pass


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
    latitudes of about > +55 or <-55 degrees).
    """
    pass


class SkyCoordNotScalarError(Exception):
    """Raised when an array-based SkyCoord object is given
    but a scalar SkyCoord object is required."""
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

    sun_dark_alt : float, optional
        Altitude of sun's center, in degrees, when has just become sufficiently dark to
        begin observations, i.e., the end of twilight.
        Default is value of ``site``'s ``sun_altitude_dark`` attribute.

    Attributes
    ----------
    site : |Site|
        Observer's location info, from input parameter ``site``.
    site_name : str
        Name of the observer's location, from input parameter ``site``.
    sun_altitude_dark : float
        Maximum sun altitude, in degrees, that is considered observably dark.
    timespan_no_sun : |Timespan|
        Timespan when sun is down for this astronight, at this site.
    timespan_observable : |Timespan|
        Timespan when sky is observably dark for this astronight, at this site.
    local_middark_utc : |Time|
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
        horizon. The best observing timespan(s). This may be a list of two |Timespan|
        if moon is up in middle of the night (very rare).

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
      >>> dark = an.timespan_observable
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
    def __init__(self, site: Site, an_date: int | str,
                 sun_dark_alt: float | None = None) -> None:
        # Handle inputs:
        self.site = site
        self.site_name = site.name
        self.an_date = AN_date(an_date)
        if sun_dark_alt is None:
            self.sun_altitude_dark = site.sun_altitude_dark
        else:
            self.sun_altitude_dark = sun_dark_alt

        # Select the almanac engine to use (the ONLY time we're allowed to reference
        #   a subclass of abstract class AlmanacEngine):
        self.engine = SkyfieldEngine.from_site(site=site,
                                               sun_altitude_dark=self.sun_altitude_dark)
        engine = self.engine  # shorter alias.

        # Get utc_approx_midnight from an_date and longitude, to use as starting point:
        approx_midnight_utc = Time(datetime(self.an_date.year, self.an_date.month,
                                            self.an_date.day, 0, 0, 0), scale='utc') +\
            TimeDelta(-site.longitude_hours * 3600, format='sec') + \
            TimeDelta(24 * 3600, format='sec')

        # Get time when sun is farthest below horizon:
        sun_antitransit_utc = engine.sun_antitransit_utc(approx_midnight_utc)
        _, sun_antitransit_alt = engine.sun_azalt(time=sun_antitransit_utc)
        if sun_antitransit_alt > HORIZON_USNO:
            raise SunAlwaysUpError(f"Astronight date {self.an_date.an_str}, "
                                   f"latitude {site.latitude:.2f}")
        if sun_antitransit_alt > self.sun_altitude_dark:
            raise NoDarkTimeError(f"Astronight date {self.an_date.an_str}, "
                                  f"latitude {site.latitude:.2f}")

        # No-sun (sun below horizon) timespan:
        self.sunset_utc = engine.prev_sunset_utc(ref_time_utc=sun_antitransit_utc)
        self.sunrise_utc = engine.next_sunrise_utc(ref_time_utc=sun_antitransit_utc)
        self.timespan_no_sun = Timespan(self.sunset_utc, self.sunrise_utc)

        # Calculate dark (sun sufficiently below horizon) timespan,
        # but limited to 12 hours to each side of lowest sun angle:
        self.dark_start_utc = \
            engine.prev_dark_start_utc(ref_time_utc=sun_antitransit_utc)
        self.dark_end_utc = \
            engine.next_dark_end_utc(ref_time_utc=sun_antitransit_utc)
        self.observable_start_utc = max(self.dark_start_utc,
                                        sun_antitransit_utc - timedelta(hours=12))
        self.observable_end_utc = min(self.dark_end_utc,
                                      sun_antitransit_utc + timedelta(hours=12))
        self.timespan_observable = Timespan(self.observable_start_utc,
                                            self.observable_end_utc)

        # Other reference times for this Astronight:
        self.middark_utc = self.timespan_observable.midpoint
        self.middark_lst = engine.local_sidereal_time(time=self.middark_utc)

        # Moon quantities for this Astronight:
        self.moon_illumination = engine.moon_illumination(time=self.middark_utc)
        self.moon_skycoord = engine.moon_skycoord(time=self.middark_utc)
        self.moon_transit_utc = engine.moon_transit_time(ref_time_utc=self.middark_utc)
        moonsets_moonrises = engine.moonsets_and_moonrises(
            self.timespan_observable.start,
            self.timespan_observable.end)
        self.timespan_dark_no_moon = self._calc_dark_no_moon(moonsets_moonrises)

    def _calc_dark_no_moon(self, moonsets_moonrises: List[Tuple]) -> Timespan | None:
        """ Returns (longest if > 1) contiguous dark-and-no-moon timespan. """
        timespan_edge_times = \
            [self.timespan_observable.start] + \
            [t for (t, value) in moonsets_moonrises] + \
            [self.timespan_observable.end]
        timespans = [Timespan(timespan_edge_times[i], timespan_edge_times[i + 1])
                     for i in range(len(timespan_edge_times) - 1)]
        longest_moon_down_timespan = None
        longest_duration = 0  # in seconds
        for ts in timespans:
            if ts.duration > 0:
                if (self.engine.moon_azalt(ts.midpoint))[1] < 0:
                    if ts.seconds > longest_duration:
                        longest_moon_down_timespan = ts
                        longest_duration = ts.seconds
        return longest_moon_down_timespan

    def __repr__(self) -> str:
        return f'astropack.almanac.Astronight(site=[{self.site.name}],'\
               f' an_date={self.an_date.an_str},' \
               f' sun_altitude_dark={self.sun_altitude_dark})'

    def __str__(self) -> str:
        return f'Astronight "{self.an_date.an_str}" at site "{self.site_name}".'


__________SUPPORT_FUNCTIONS______________________________________________________ = 0


def calc_approx_midnight(an_date: AN_date_type, longitude_hours: float) -> Time:
    """Return an estimate of midnight for given AN_date, using longitude in hours
    which is measured from Greenwich meridian and which respects the International
    Date Line so that the correct midnight is chosen.

    Parameters
     ----------
    an_date : AN_date
        Astronight date for which approximate midnight is wanted.

    longitude_hours :
        longitude in hours as measured from Greenwich meridian and which
        respects the International Date Line so that the correct midnight is chosen.

    Returns
    -------
    approx_midnight : |Time|
       Approximate estimate of midnight for given AN date and given longitude in hours.
       May be used as starting point for determining that night's almanac times,
       especially sun-down and dark timespans, solar-midnight time, etc.
    """
    return Time(datetime(an_date.year, an_date.month,
                         an_date.day, 0, 0, 0), scale='utc') + \
           TimeDelta(-longitude_hours * 3600, format='sec') + \
           TimeDelta(24 * 3600, format='sec')


# def calc_timespan_no_sun(approx_midnight: Time) -> [Timespan, None]:
#     """Return |Timespan| representing sunset to sunrise and containing
#     ``approx_midnight``.
#
#     Parameters
#     ----------
#     approx_midnight : |Time|
#         Approximate local midnight, UTC.
#
#     Returns
#     -------
#     timespan_no_sun : |Timespan|, or None
#         Timespan from sunset to sunrise.
#         None if there is no time near ``approx_midnight`` when sun is down
#         (polar summer).
#     """
#
#     return skyfield_timespan_no_sun(approx_midnight)


# def skyfield_timespan_no_sun(approx_midnight: Time) -> [Timespan, None]:
#     """Return |Timespan| representing sunset to sunrise and containing
#     ``approx_midnight``.
#
#     Parameters
#     ----------
#     approx_midnight : |Time|
#         Approximate local midnight, UTC.
#
#     Returns
#     -------
#     timespan_no_sun : |Timespan|, or None
#         Timespan from sunset to sunrise.
#         None if there is no time near ``approx_midnight`` when sun is down
#         (polar summer).
#     """
#     midnight = sf_obs.sf_timescale.from_astropy(approx_midnight)  # skyfield Time obj.
#     time_before = midnight - timedelta(hours=HOURS_TO_ASSURE_HALF_NIGHT)
#     time_after = midnight + timedelta(hours=HOURS_TO_ASSURE_HALF_NIGHT)
#     sun_fn = risings_and_settings(sf_obs.sf_master_eph, sf_obs.sf_master_eph['sun'],
#                                   sf_obs.obs, HORIZON_USNO)
#     sun_fn.step_days = 1.0 / 12.0
#     set_times, set_values = find_discrete(time_before, midnight, sun_fn)
#     rise_times, rise_values = find_discrete(midnight, time_after, sun_fn)
#     topos_at = (sf_obs.sf_master_eph['earth'] + sf_obs.obs).at
#     alt_midnight = topos_at(midnight).observe(sf_obs.sf_master_eph['sun'])\
#         .apparent().altaz()[0].degrees
#
#     # Case: no rises or sets, sun is either up or down for the entire search period.
#     if len(set_values) <= 0 and len(rise_values) <= 0:
#         if alt_midnight <= HORIZON_USNO:
#             return Timespan(_skyfield_time_to_astropy(time_before),
#                             _skyfield_time_to_astropy(time_after))
#         else:
#             return None
#     if len(set_values) >= 1:
#         set_time = _skyfield_time_to_astropy(set_times[-1])
#     else:
#         set_time = _skyfield_time_to_astropy(time_before)
#     if len(rise_values) >= 1:
#         rise_time = _skyfield_time_to_astropy(rise_times[0])
#     else:
#         rise_time = _skyfield_time_to_astropy(time_after)
#     return Timespan(set_time, rise_time)


# def find_target_up_down(sf_obs: SkyfieldObserver,
#                         target_eph: [skyfield.vectorlib.VectorSum,
#                                      skyfield.starlib.Star],
#                         timespan: Timespan, up_down: str, horizon: float = 0.0)\
#                         -> [list[Timespan], None]:
#     """Returns |Timespan| for which target is either up or down
#     (above or below given horizon), within a given timespan.
#
#     Parameters
#     ----------
#     sf_obs : SkyfieldObserver
#         SkyfieldObserver instance.
#
#     target_eph : skyfield.vectorlib.VectorSum or skyfield.starlib.Star
#         Target ephemeris; either a solar-system ephemeris (selected from
#         ``master_eph``) or a skyfield Star object.
#
#         See:
#         https://github.com/skyfielders/python-skyfield/blob/master/skyfield/vectorlib.py
#         or https://rhodesmill.org/skyfield/api-stars.html#skyfield.starlib.Star
#
#     timespan : |Timespan|
#         Timespan within which result will be constrained.
#         If ``up_down`` is 'down', one convenient source for ``timespan``
#         is ``.util.Astronight.timespan_dark``.
#
#     up_down : {'up', 'down'}
#         If 'up', results represent when target is above horizon.
#         If 'down', results represent when target is below horizon. Required.
#
#     horizon : float, optional
#         Altitude above actual horizon that separates 'up' from 'down', in degrees.
#         Default is zero (ideal horizon).
#
#     Returns
#     -------
#     timespan_up_down : list of |Timespan|, or None.
#         List of timespans for which target is either up or down as selected.
#         Returns a list even if only one Timespan is returned.
#         If ``up_down`` is satisfied throughout input timespan, the returned list will
#         contain a copy of ``timespan``.
#         If ``up_down`` is never satisfied during input timespan, return None.
#     """
#     up_down = up_down.lower().strip()
#     if up_down not in ['up', 'down']:
#         raise ValueError('Parameter up_down must equal \'up\' or \'down\'.')
#     start = sf_obs.sf_timescale.from_astropy(timespan.start)
#     end = sf_obs.sf_timescale.from_astropy(timespan.end)
#     dark_fn = risings_and_settings(sf_obs.sf_master_eph, target_eph, sf_obs.obs,
#                                    horizon)
#     event_times, event_values = find_discrete(start, end, dark_fn)
#     # Case: no risings or settings found:
#     if len(event_values) <= 0:
#         topos_at = (sf_obs.sf_master_eph['earth'] + sf_obs.obs).at
#         mid_time = sf_obs.sf_timescale.from_astropy(timespan.midpoint)
#         alt = topos_at(mid_time).observe(target_eph).apparent().altaz()[0].degrees
#         if (up_down == 'up' and alt > horizon) or \
#            (up_down == 'down' and alt <= horizon):
#             return [timespan.copy()]
#         else:
#             return None
#     # Case: at least one rising or setting found:
#     # First, convert times to astropy, so that start and end are passed through exactly.
#     astropy_times = [timespan.start] + \
#                     [_skyfield_time_to_astropy(t) for t in event_times] + \
#                     [timespan.end]
#     values = [0.5] + list(event_values) + [0.5]
#     diffs = diff(values)
#     # Now choose events as up or down:
#     if up_down == 'up':
#         is_index_qualifying = diffs < 0
#     else:
#         is_index_qualifying = diffs > 0
#     qualifying_indices = flatnonzero(diffs)[is_index_qualifying]
#     timespans = [Timespan(astropy_times[i], astropy_times[i + 1])
#                  for i in qualifying_indices]
#     return timespans


# def moon_ra_dec(sf_obs: SkyfieldObserver, time: Time) -> tuple[float, float]:
#     """Returns moon's sky coordinates at given earth location and time.
#
#     Parameters
#     ----------
#     sf_obs : SkyfieldObserver
#         SkyfieldObserver instance.
#
#     time : |Time|
#         Time UTC at which moon's sky location is wanted.
#
#     Returns
#     -------
#     ra, dec : tuple of float
#         (RA, Declination), in degrees, of moon at observer location and given time.
#     """
#     time_sf = sf_obs.sf_timescale.from_astropy(time)
#     ra, dec, _ = (sf_obs.sf_master_eph['earth'] + sf_obs.obs).at(time_sf).\
#         observe(sf_obs.sf_master_eph['moon']).apparent().radec()
#     return ra.hours * 15.0, dec.degrees


# def moon_transit_time(sf_obs: SkyfieldObserver, time: Time) -> Time:
#     """Return moon's transit time closest to given time.
#
#     Parameters
#     ----------
#     sf_obs : SkyfieldObserver
#         SkyfieldObserver instance.
#
#     time : |Time|
#         Time UTC from which moon's nearest transit time is wanted.
#
#     Returns
#     -------
#     transit_time : |Time|
#         Moon's transit time UTC at observer's location and nearest to given time.
#     """
#     time_sf = sf_obs.sf_timescale.from_astropy(time)  # skyfield Time obj.
#     time_start = time_sf - timedelta(hours=14)
#     time_end = time_sf + timedelta(hours=14)
#     moon_fn = meridian_transits(sf_obs.sf_master_eph, sf_obs.sf_master_eph['moon'],
#                                 sf_obs.obs)
#     times, values = find_discrete(time_start, time_end, moon_fn)
#     return _skyfield_time_to_astropy(_closest_transit(times, values, time_sf))


# def target_transit_time(sf_obs: SkyfieldObserver,
#                         target_skycoord: SkyCoord,
#                         time: Time) -> Time:
#     """Return fixed target's transit time closest to a given time.
#
#     Parameters
#     ----------
#     sf_obs : SkyfieldObserver
#         SkyfieldObserver instance.
#
#     target_skycoord : |SkyCoord|
#         Target's fixed sky position.
#
#     time : |Time|
#         Time UTC nearest which target's nearest transit time is wanted.
#
#     Returns
#     -------
#     transit_time : |Time|
#         Moon's transit time UTC at observer's location and nearest to given time
#     """
#     star_list = _skycoords_to_skyfield_stars(target_skycoord)
#     time_sf = sf_obs.sf_timescale.from_astropy(time)  # skyfield Time obj.
#     time_start = time_sf - timedelta(hours=14)
#     time_end = time_sf + timedelta(hours=14)
#     transit_list = []
#     for star in star_list:
#         star_fn = meridian_transits(sf_obs.sf_master_eph, star, sf_obs.obs)
#         times, values = find_discrete(time_start, time_end, star_fn)
#         transit_list.append(_skyfield_time_to_astropy(_closest_transit(times,
#                                                                       values,
#                                                                       time_sf)))
#     if len(transit_list) == 1:
#         return transit_list[0]  # return single time as scalar Time rather than list.
#     return transit_list


# def _closest_transit(times: list[skyfield.timelib.Time],
#                      values: list[int],
#                      time_sf: skyfield.timelib.Time) -> skyfield.timelib.Time:
#     """For times [list of skyfield Time objects], values [list of ints], and
#     a target time [Skyfield Time object], determine which transit is closest
#     to target time, return the time as a Skyfield Time object."""
#     transit_times = [t for (t, v) in zip(times, values) if v == 1]
#     timedeltas = [abs(time - time_sf) for time in transit_times]
#     transit_time = transit_times[timedeltas.index(min(timedeltas))]
#     return transit_time


# def local_sidereal_time(longitude: float,
#                         timescale: skyfield.timelib.Timescale,
#                         time: Time) -> float:
#     """Return local sidereal time, in hours, for a given time and earth longitude.
#
#     Parameters
#     ----------
#     longitude : float
#         Longitude of observer's earth location, in degrees, in range [-180, +180].
#
#     timescale:  skyfield.timelib.Timescale
#         Skyfield timescale object, available as Astronight.skyfield.ts
#
#         See: |skyfield.TS|
#
#     time : |Time|
#         Time UTC for which local sidereal time is wanted.
#
#     Returns
#     -------
#     lst : float
#         Local sidereal time at observer's longitude and time, in hours
#     """
#     from skyfield.earthlib import sidereal_time
#     time_sf = timescale.from_astropy(time)  # skyfield Time obj.
#     gmst = sidereal_time(time_sf)  # in hours
#     local_offset = longitude / 15.0  # in hours
#     return gmst + local_offset


# def moon_illumination_pct(sf_obs: SkyfieldObserver, time: Time) -> float:
#     """Return moon illumination extent, as percentage, at a given time.
#
#     Parameters
#     ----------
#     sf_obs : SkyfieldObserver
#         SkyfieldObserver instance.
#
#     time : |Time|
#         Time UTC for which moon's illumination is wanted.
#
#     Returns
#     -------
#     moon_pct : float
#         Moon's percent illumination at given time, in range [0, 100].
#     """
#     time_sf = sf_obs.sf_timescale.from_astropy(time)  # skyfield Time obj.
#     fraction = 100.0 * fraction_illuminated(sf_obs.sf_master_eph, 'moon', time_sf)
#     return fraction


# def make_skyfield_observatory_from_site(site: Site) \
#                                         -> skyfield.toposlib.GeographicPosition:
#     """ Return skyfield observer (earth location) object made
#     from astropack |Site| instance.
#
#     Parameters
#     ----------
#     site : |Site|
#         Earth location, as read from .ini file by |Site| class.
#
#     Returns
#     -------
#     skyfield_obs : skyfield.toposlib.GeographicPosition
#         Skyfield GeographicPosition instance.
#
#         See: |skyfield.GP|
#     """
#     obs = wgs84.latlon(site.latitude, site.longitude, site.elevation)
#     return obs

# def _skycoords_to_skyfield_stars(skycoords: [SkyCoord, list[SkyCoord]]) -> list[Star]:
#     """Convert astropy SkyCoord object or list of SkyCoord objects
#     to a *list of* skyfield Star objects.
#     Always returns a list, only one sky position per Star.
#     """
#     # Nested function, always returns a list of Star objects.
#     def skycoord_object_to_skyfield_stars(skycoord_obj: SkyCoord):
#         """Nested function, always returns a list of Star objects."""
#         if skycoord_obj.size == 1:
#             return [Star(ra_hours=skycoord_obj.ra.degree / 15.0,
#                          dec_degrees=skycoord_obj.dec.degree)]
#         else:
#             return [Star(ra_hours=sc_element.ra.degree / 15.0,
#                          dec_degrees=sc_element.dec.degree)
#                     for sc_element in skycoord_obj]
#
#     # Build and return flattened list of skyfield star objects:
#     if isinstance(skycoords, SkyCoord):
#         return skycoord_object_to_skyfield_stars(skycoords)
#     elif isinstance(skycoords, list):
#         return [coord for sc_obj in skycoords for coord in sc_obj]


__________CLASS_AN_DATE_________________________________________________________ = 0


class AN_date:
    """Builds, holds, and delivers data based on an Astronight date format."""
    def __init__(self, an_input: [str, int]) -> None:
        self.an_str = str(an_input)
        try:
            self.an_int = int(self.an_str)
        except ValueError:
            raise Invalid_ANDate_Error(f'Input >{self.an_str}<.')
        today = datetime.utcnow()
        an_date_today = 10000 * today.year + 100 * today.month + today.day
        an_date_limits = 19700000, an_date_today + 50000  # 1970 to 5 years from now.
        if self.an_int < an_date_limits[0] or self.an_int > an_date_limits[1]:
            raise Invalid_ANDate_Error(f'Date {self.an_str} is unreasonable.')
        self.year = self.an_int // 10000
        self.month = (self.an_int - (10000 * self.year)) // 100
        self.day = self.an_int % 100
        if self.month == 0 or self.month > 12 or self.day == 0 or self.day > 31:
            raise Invalid_ANDate_Error(f'Date {self.an_str} does not represent '
                                       'a valid calendar date.')

__________ABSTRACT_CLASS_ENGINE_________________________________________________ = 0


class AlmanacEngine(ABC):
    """ Abstract class. Defines interface guaranteed to apply to implementations
    of an almanac engine (that is available to user code), e.g., skyfield package.

    All functions, input values, and return values are independent of engine used.
    """
    class Event(Enum):
        """ Type of (meridian or horizon) crossing event."""
        SETTING = auto()
        RISING = auto()
        TRANSIT = auto()
        ANTITRANSIT = auto()
        DARK_START = auto()
        DARK_END = auto()


    @abstractmethod
    def sun_azalt(self, time: Time) -> Tuple[float, float]:
        """
        Parameters
        ----------
        time : |Time|
            Time UTC for which sun's local azimuth and altitude is wanted.

        Returns
        -------
        az, alt : tuple of 2 floats
            Sun's local azimuth and altitude at given time.
        """
        pass

    @abstractmethod
    def moon_azalt(self, time: Time) -> Tuple[float, float]:
        """
        Parameters
        ----------
        time : |Time|
            Time UTC for which moon's local azimuth and altitude is wanted.

        Returns
        -------
        az, alt : tuple of 2 floats
            Moon's local azimuth and altitude at given time.
        """
        pass

    @abstractmethod
    def skycoord_azalt(self, target: SkyCoord, time: Time) -> Tuple[float, float]:
        """
        Parameters
        ----------
        target : |SkyCoord|
            SkyCoord for target (e.g., star) of constant RA and Declination.

        time : |Time|
            Time UTC for which moon's local azimuth and altitude is wanted.

        Returns
        -------
        az, alt : tuple of 2 floats
            Target's local azimuth and altitude at given time.
        """
        pass

    @abstractmethod
    def prev_sunset_utc(self, ref_time_utc: Time,
                        horizon_degrees=HORIZON_USNO) -> Time:
        """
        Parameters
        ----------
        ref_time_utc : |Time|
            Reference time for which latest previous sunset's time is to be found.

        horizon_degrees : float
            Horizon altitude, below which the sun is considered to have set, in degrees.
            Default = USNO's refraction-aware horizon angle (about -0.833 degrees).

        Returns
        -------
        prev_sunset_utc : |Time|
            Time (UTC) of latest sunset previous to reference time.
        """
        pass

    @abstractmethod
    def next_sunrise_utc(self, ref_time_utc: Time,
                         horizon_degrees=HORIZON_USNO) -> Time:
        """
        Parameters
        ----------
        ref_time_utc : |Time|
            Reference time for which earliest next sunset's time is to be found.

        horizon_degrees : float
            Horizon altitude, below which the sun is considered not yet to have risen,
            in degrees.
            Default = USNO's refraction-aware horizon angle (about -0.833 degrees).

        Returns
        -------
        next_sunrise_utc : |Time|
            Time (UTC) of earliest sunset after reference time.
        """
        pass

    @abstractmethod
    def prev_dark_start_utc(self, ref_time_utc: Time, sun_alt_dark: float=None) -> Time:
        """
        Parameters
        ----------
        ref_time_utc : |Time|
            Reference time for which earliest next sunset's time is to be found.

        sun_alt_dark : float
            Sun altitude below which the sky is considered to be dark.
            Typically, about -9 degrees.

        Returns
        -------
        prev_dark_start_utc : |Time|
            Time (UTC) of latest start of sky darkness (twilight)
            previous to reference time.
        """
        pass

    @abstractmethod
    def next_dark_end_utc(self, ref_time_utc: Time, sun_alt_dark: float=None) -> Time:
        """
        Parameters
        ----------
        ref_time_utc : |Time|
            Reference time for which earliest next sunset's time is to be found.

        sun_alt_dark : float
            Sun altitude below which the sky is considered to be dark.
            Typically, about -9 degrees.

        Returns
        -------
        next_dark_end_utc : |Time|
            Time (UTC) of earliest start of sky darkness (twilight)
            after reference time.
        """
        pass

    @abstractmethod
    def sun_antitransit_utc(self, ref_time_utc: Time) -> Time:
        """ Returns tuple of (sun antitransit time, sun antitransit altitude),
        for given site location, nearest the approximate midnight time given.

        Parameters
        ----------
        ref_time_utc : |Time|
            Approximate local midnight, UTC.

        Returns
        -------
        sun_antitransit_utc : |Time|
            Time of sun's antitransit that is nearest the given reference time.
        """
        pass

    @abstractmethod
    def moon_skycoord(self, time: Time) -> SkyCoord:
        """Returns moon's sky coordinates at site location and given time.

        Parameters
        ----------
        time : |Time|
            Time UTC at which moon's sky location is wanted.

        Returns
        -------
        moon_skycoord : |SkyCoord|
            Sky coordinates of moon at site location and given time.
        """
        pass

    @abstractmethod
    def moon_transit_time(self, ref_time: Time) -> Time:
        """Return moon's transit time nearest the given reference time.

        Parameters
        ----------
        ref_time : |Time|
            Time UTC from which moon's nearest transit time is wanted.

        Returns
        -------
        nearest_transit_time : |Time|
            Moon's transit time UTC at observer's location and nearest to given time.
        """
        pass

    @abstractmethod
    def moonsets_and_moonrises(self, start_utc: Time, end_utc: Time) \
            -> List[Tuple[Time, Event]]:
        """ XXX """

    @abstractmethod
    def moon_illumination(self, time: Time) -> float:
        """Return moon illumination extent, as fraction in range 0-1, at a given time.

        Parameters
        ----------
        time : |Time|
            Time UTC for which moon's illumination is wanted.

        Returns
        -------
        moon_illumination_fraction : float
            Moon's illumination fraction at given time, in range [0, 1].
        """
        pass

    @abstractmethod
    def local_sidereal_time(self, time: Time) -> float:
        """Return local sidereal time, in hours, at the site location, for a given time.

        Parameters
        ----------
        time : |Time|
            Time UTC for which local sidereal time is wanted.

        Returns
        -------
        lst : float
            Local sidereal time at observer's longitude and time, in hours.
        """
        pass


__________CLASS_SKYFIELDENGINE_________________________________________________ = 0


class SkyfieldEngine(AlmanacEngine):
    """ Holds skyfield data, makes skyfield-dependent facilities available
    to higher-level functions without *any* other reference to skyfield package. """

    # -----------------------------------------------------------------------
    # CONSTRUCTORS:

    def __init__(self,
                 site_latitude: float, site_longitude: float, site_elevation: float,
                 sun_altitude_dark: float):
        """
        Parameters
        ----------
        site_latitude : float
            Observer's latitude in degrees.

        site_longitude : float
            Observer's longitude in degrees, East positive, West negative.

        site_elevation : float
            Observer's elevation above sea level in meters.

        sun_altitude_dark : float
            Sun altitude relative to horizontal above which sky is considered light
                and below which sky is considered dark, in degrees. Typically about -9.
        """
        self.master_eph = load('de440s.bsp')
        self.timescale = load.timescale()
        self.sf_obs = wgs84.latlon(site_latitude, site_longitude, site_elevation)
        self.site_latitude = site_latitude
        self.site_longitude = site_longitude
        self.site_elevation = site_elevation
        self.topos_at = (self.master_eph['earth'] + self.sf_obs).at  # a function.
        self.sun_eph = self.master_eph['sun']
        self.moon_eph = self.master_eph['moon']
        # self.setting_value = 0
        # self.rising_value = 1
        self.sun_alt_dark = sun_altitude_dark

    @classmethod
    def from_site(cls, site: Site, sun_altitude_dark: float):
        """ Alternate constructor, using Site object. """
        return cls(site.latitude, site.longitude, site.elevation, sun_altitude_dark)

    # -----------------------------------------------------------------------
    # LOWEST-LEVEL, PRIVATE ENUMS and METHODS:

    class _HorizonCrossing(Enum):
        """ Distinguising rising from setting (horizon crossing, even if "horizon" is
        set to non-zero altitude). Values match skyfield convention.
        PRIVATE, for use internal to subclass SkyfieldEngine ONLY."""
        RISING = 1
        SETTING = 0

    class _MeridianCrossing(Enum):
        """ Distinguishing meridian transit from antitransit (antimeridian transit).
        Values match skyfield convention.
        PRIVATE, for use internal to subclass SkyfieldEngine ONLY."""
        TRANSIT = 1
        ANTITRANSIT = 0

    def _target_azalt(self, target: skyfield.vectorlib.VectorSum |
                                    skyfield.starlib.Star |
                                    skyfield.jpllib.ChebyshevPosition,
                      time: skyfield.timelib.Time) -> Tuple[float, float]:
        """ Given time (UTC), return altitude of one (stationary in RA,Dec) target,
        in degrees, at that time. """
        alt, az, _ = self.topos_at(time).observe(target).apparent().altaz()
        return az.degrees, alt.degrees

    def _horizon_crossings(self, target: skyfield.vectorlib.VectorSum |
                                         skyfield.starlib.Star |
                                         skyfield.jpllib.ChebyshevPosition,
                           start_utc: skyfield.timelib.Time,
                           end_time: skyfield.timelib.Time, horizon_degrees: float) \
            -> List[Tuple[skyfield.timelib.Time, _HorizonCrossing]]:
        """Engine/support method: return times and set/rise values when target crosses
        the given horizon angle."""
        fn = risings_and_settings(self.master_eph, target, self.sf_obs,
                                  horizon_degrees=horizon_degrees)
        times, values = find_discrete(start_utc, end_time, fn)
        crossing_types = [SkyfieldEngine._HorizonCrossing(v) for v in values]
        return list(zip(times, crossing_types))

    def _meridian_crossings(self, target: skyfield.vectorlib.VectorSum |
                                          skyfield.starlib.Star |
                                          skyfield.jpllib.ChebyshevPosition,
                            start_utc: skyfield.timelib.Time,
                            end_time: skyfield.timelib.Time) \
            -> List[Tuple[skyfield.timelib.Time, _MeridianCrossing]]:
        fn = meridian_transits(self.master_eph, target, self.sf_obs)
        times, values = find_discrete(start_utc, end_time, fn)
        crossing_types = [SkyfieldEngine._MeridianCrossing(v) for v in values]
        return list(zip(times, crossing_types))

    @staticmethod
    def _nearest_skyfield_time(ref_time: skyfield.timelib.Time,
                               times: List[skyfield.timelib.Time]) \
            -> skyfield.timelib.Time:
        """ From a list of times, return the time closest to ref_time. """
        diffs = [abs(t - ref_time) for t in times]
        return times[diffs.index(min(diffs))]

    # -----------------------------------------------------------------------
    # SKYFIELD CONVERSION METHODS (private):

    @staticmethod
    def _skycoord_to_skyfield_star(skycoord: SkyCoord) -> Star:
        """Convert one scalar-only astropy SkyCoord object to one skyfield Star object.
        Raises exception if SkyCoord object is array rather than scalar.

        Parameters
        ----------
        skycoord : |SkyCoord|
            Sky coordinates in astropy SkyCoord format. Scalar object only.

        Returns
        -------
        skyfield_star : |Star|
            Skyfield Star object, corresponding to input of astropy
            sky coordinate.

            See |skyfield.Star|
            """
        if not skycoord.isscalar:
            raise SkyCoordNotScalarError
        return Star(ra_hours=skycoord.ra.degree / 15.0,
                    dec_degrees=skycoord.dec.degree)

    @staticmethod
    def _astropy_time_to_skyfield(time: Time, sf_timescale=None) \
            -> skyfield.timelib.Time:
        """Convert one astropy time in UTC to skyfield time.

    .. Warning:: Unfortunately, skyfield and astropy both use 'Time' as name
            for their main time class.

        Parameters
        ----------
        time : |Time|
            Astropy time instance, UTC

        sf_timescale: skyfield.timelib.Timescale
            Skyfield timescale object.

        Returns
        -------
        skyfield_time : skyfield.timelib.Time
            Skyfield time instance.

            See: |skyfield.Time|
        """
        if sf_timescale is None:
            sf_timescale = load.timescale()
        return sf_timescale.from_datetime(time.to_datetime().
                                          replace(tzinfo=timezone.utc))

    @staticmethod
    def _skyfield_time_to_astropy(skyfield_time: skyfield.timelib.Time) -> Time:
        """Convert skyfield time to astropy time in UTC timescale and ISO format.
        Time arrays are OK.

    .. Warning:: Unfortunately, skyfield and astropy both use 'Time' as name
            for their main time class.

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

    # -----------------------------------------------------------------------
    # PUBLIC METHODS:

    def sun_azalt(self, time: Time | skyfield.timelib.Time) -> Tuple[float, float]:
        """ Given time (UTC), return sun azimuth, altitude in degrees at that time. """
        if isinstance(time, Time):
            time = self._astropy_time_to_skyfield(time)
        alt, az, _ = self.topos_at(time).observe(self.sun_eph).apparent().altaz()
        return az.degrees, alt.degrees

    def moon_azalt(self, time: Time | skyfield.timelib.Time) -> Tuple[float, float]:
        """ Given time (UTC), return moon azimuth, altitude in degrees at that time. """
        if isinstance(time, Time):
            time = self._astropy_time_to_skyfield(time)
        alt, az, _ = self.topos_at(time).observe(self.moon_eph).apparent().altaz()
        return az.degrees, alt.degrees

    def skycoord_azalt(self, target: SkyCoord, time: Time) \
            -> Tuple[float, float] | List[Tuple[float, float]]:
        """ Given time (UTC), return azimuth, altitude of (stationary in RA,Dec) target,
        in degrees, at that time.
        Handles only a scalar SkyCoord object, does not handle array-type SkyCoord objs.
        """
        if not target.isscalar:
            raise SkyCoordNotScalarError
        sf_target = self._skycoord_to_skyfield_star(target)
        print(f'\nsf_target: {sf_target}')
        sf_time = self._astropy_time_to_skyfield(time)
        return self._target_azalt(sf_target, sf_time)

    def prev_sunset_utc(self, ref_time_utc: Time, horizon_degrees=HORIZON_USNO) -> Time:
        """ Return sunset time nearest and previous to given reference time.
        Return None if no sunset found in prev. 25 hours (sun either always up or down).
        """
        sf_ref_time = self.timescale.from_astropy(ref_time_utc)
        sf_start_time = sf_ref_time - timedelta(hours=36)
        sf_end_time = sf_ref_time
        # Get all horizon crossings:
        crossings_list = self._horizon_crossings(self.sun_eph,
                                                 sf_start_time, sf_end_time,
                                                 horizon_degrees=horizon_degrees)
        # Select only crossings that are sunsets:
        sunset_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                        for (sf_time, value) in crossings_list
                        if value == SkyfieldEngine._HorizonCrossing.SETTING]
        return max(sunset_times)

    def next_sunrise_utc(self, ref_time_utc: Time,
                         horizon_degrees=HORIZON_USNO) -> Time:
        """ Return sunrise time nearest and after given reference time.
        Return None if no sunrise found in next 25 hours (sun either always up or down).
        """
        sf_ref_time = self.timescale.from_astropy(ref_time_utc)
        sf_start_time = sf_ref_time
        sf_end_time = sf_ref_time + timedelta(hours=36)
        # Get all horizon crossings:
        crossings_list = self._horizon_crossings(self.sun_eph,
                                                 sf_start_time, sf_end_time,
                                                 horizon_degrees=horizon_degrees)
        # Select only crossings that are sunrises:
        sunrise_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                         for (sf_time, value) in crossings_list
                         if value == SkyfieldEngine._HorizonCrossing.RISING]
        return min(sunrise_times)

    def prev_dark_start_utc(self, ref_time_utc: Time, sun_alt_dark: float=None) -> Time:
        """ Return going-dark time nearest and previous to given reference time,
        based on angle of sun below horizon that user has defined as twilight in
        sun_altitude_dark.
        """
        if sun_alt_dark is None:
            sun_alt_dark = self.sun_alt_dark
        return self.prev_sunset_utc(ref_time_utc=ref_time_utc,
                                    horizon_degrees=sun_alt_dark)

    def next_dark_end_utc(self, ref_time_utc: Time, sun_alt_dark: float=None) -> Time:
        """ Return going-light time nearest and after given reference time,
        based on angle of sun below horizon that user has defined as twilight in
        sun_altitude_dark.
        """
        if sun_alt_dark is None:
            sun_alt_dark = self.sun_alt_dark
        return self.next_sunrise_utc(ref_time_utc=ref_time_utc,
                                     horizon_degrees=sun_alt_dark)

    def sun_antitransit_utc(self, ref_time_utc: Time) -> Time:
        """ Returns sun antitransit time for given site location,
        nearest the approximate time given (ref_time is normally near midnight).
        *** See abstract method for parameters and return value.
        """
        sf_approx_midnight = self._astropy_time_to_skyfield(ref_time_utc)
        start_time = sf_approx_midnight - timedelta(hours=25)
        end_time = sf_approx_midnight + timedelta(hours=25)
        meridian_crossings = self._meridian_crossings(target=self.sun_eph,
                                                      start_utc=start_time,
                                                      end_time=end_time)
        antitransit_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                             for (sf_time, value) in meridian_crossings
                             if value == SkyfieldEngine._MeridianCrossing.ANTITRANSIT]
        diffs = [abs(t - ref_time_utc).sec for t in antitransit_times]
        nearest_antitransit_time = antitransit_times[diffs.index(min(diffs))]
        return nearest_antitransit_time

    def moon_skycoord(self, time: Time) -> SkyCoord:
        """Returns moon's sky coordinates at site location and given time.
        *** See abstract method for parameters and return value.
        """
        time_sf = self.timescale.from_astropy(time)
        ra, dec, _ = (self.master_eph['earth'] + self.sf_obs).at(time_sf). \
            observe(self.master_eph['moon']).apparent().radec()
        return SkyCoord(ra=ra.hours * 15.0, dec=dec.degrees, unit="deg")

    def moon_transit_time(self, ref_time_utc: Time) -> Time:
        """Return moon's transit time nearest the given time.
        *** See abstract method for parameters and return value.
        """
        time_sf = self.timescale.from_astropy(ref_time_utc)  # skyfield Time obj.
        time_start = time_sf - timedelta(hours=15)
        time_end = time_sf + timedelta(hours=15)
        meridian_crossings = self._meridian_crossings(self.moon_eph,
                                                      time_start, time_end)
        transit_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                         for (sf_time, value) in meridian_crossings
                         if value == SkyfieldEngine._MeridianCrossing.TRANSIT]
        diffs = [abs(t - ref_time_utc).sec for t in transit_times]
        nearest_transit_time = transit_times[diffs.index(min(diffs))]
        return nearest_transit_time

    def moonsets_and_moonrises(self, start_utc: Time, end_utc: Time) \
            -> List[Tuple[Time, AlmanacEngine.Event]]:
        """Return list of all moonsets and moonrises from start time to end time."""
        sf_start_time = self._astropy_time_to_skyfield(start_utc, self.timescale)
        sf_end_time = self._astropy_time_to_skyfield(end_utc, self.timescale)

        # Get all horizon crossings (with skyfield times and skyfield values):
        sf_crossings_list = self._horizon_crossings(self.moon_eph,
                                                    sf_start_time, sf_end_time,
                                                    horizon_degrees=0)
        # Convert to astropy times and AlamancEngine.Event values, return:
        astropy_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                         for (sf_time, sf_value) in sf_crossings_list]
        sf_setting = SkyfieldEngine._HorizonCrossing.SETTING
        event_values = [AlmanacEngine.Event.SETTING if sf_value == sf_setting
                        else AlmanacEngine.Event.RISING
                        for (sf_time, sf_value) in sf_crossings_list]
        return list(zip(astropy_times, event_values))

    def moon_illumination(self, time: Time) -> float:
        """Return moon illumination extent, as fraction in range 0-1, at a given time.
        *** See abstract method for parameters and return value.
        """
        time_sf = self.timescale.from_astropy(time)  # skyfield Time obj.
        # deprecated as of skyfield 1.42:
        #     return fraction_illuminated(self.master_eph, 'moon', time_sf)
        m = self.topos_at(time_sf).observe(self.moon_eph).apparent()
        return m.fraction_illuminated(self.sun_eph)

    def local_sidereal_time(self, time: Time) -> float:
        """Return local sidereal time in hours, at the site location, for a given time.
        *** See abstract method for parameters and return value.
        """
        from skyfield.earthlib import sidereal_time
        time_sf = self.timescale.from_astropy(time)  # skyfield Time obj.
        gmst = sidereal_time(time_sf)  # in hours
        local_offset = self.site_longitude / 15.0  # in hours
        return gmst + local_offset


