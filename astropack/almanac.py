""" Module astropack.almanac:
    Astronight class + other almanac code. (py file)
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timezone, timedelta, date
from math import sqrt
from abc import ABC, abstractmethod
from typing import TypeAlias, List, Tuple
from enum import Enum, auto

# External packages:
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, get_sun, EarthLocation, angular_separation
import skyfield.vectorlib
import skyfield.starlib
import skyfield.timelib
import skyfield.toposlib
from skyfield.api import load, wgs84, Star
from skyfield.almanac import risings_and_settings, meridian_transits
from skyfield.searchlib import find_discrete

# Astropack packages:
from .util import Timespan, nearest_time, hhmm, ra_as_hours
from .ini import Site

THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')

HORIZON_USNO = -0.833  # USNO sun/moon effective horizon (radius & refraction), degrees.
HOURS_TO_ASSURE_HALF_NIGHT = 14

DAYS_OF_WEEK = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday',
                'Saturday', 'Sunday']

__all__ = ['Astronight',
           'NoDarkTimeError',
           'SunAlwaysUpError',
           'AN_date',
           'calc_phase_angle_bisector'
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

    sun_dark_alt : float or None
        Altitude of sun's center, in degrees, when has just become sufficiently dark to
        begin observations, i.e., the end of twilight.
        If None, use default value of ``site``'s ``sun_altitude_dark`` attribute.

    Attributes
    ----------
    site : |Site|
        Observer's location info, from input parameter ``site``.
    site_name : str
        Name of the observer's location, from input parameter ``site``.
    an_date : AN_date
        Astronight date object.
    sun_altitude_dark : float
        Maximum sun altitude, in degrees, that is considered observably dark.
    engine : AlmanacEngine subclass
        AlmanacEngine subclass object, e.g., SkyfieldEngine object
    sun_antitransit_utc : |Time|
        Sun antitransit (lowest below the horizon) time for this night, UTC.
    sunset_utc : |Time|
        Sunset time for this night, UTC.
    sunrise_utc : |Time|
        Sunrise time for this night, UTC.
    timespan_no_sun : |Timespan|
        Timespan when sun is down for this astronight, at this site.
    dark_start_utc : |Time|
        Time when sun sets below given .sun_dark_alt, i.e., twilight, UTC.
    dark_end_utc : |Time|
        Time when sun rises above given .sun_dark_alt, i.e., twilight, UTC.
    timespan_dark : |Timespan|
        Timespan when sun is below .sun_dark_alt, i.e., sky is observably dark, UTC.
    moon_illumination : float
        Fraction illumination of moon at .sun_antitransit_utc, in range [0-1].
    moon_skycoord : |SkyCoord|
        SkyCoord object holding moon's sky location at .sun_antitransit_utc.
    moon_transit_utc : |Time|
        Moon's transit time nearest to .sun_antitransit_utc.
    timespan_dark_no_moon : |Timespan|
        Timespan when sun is below its user-chosen dark elevation, and moon is below
        horizon. The best observing timespan(s). If moon is up in the middle of
        the night and there are two such timespans, the longer is represented here.

    Raises
    ------
    Invalid_ANDate_Error
        Raised when AN_Date object cannot be constructed from given input.

    SunAlwaysUpError
        Raised when |Astronight| at this site has no sun-down time for observing.

    NoDarkTimeError
        Raised when |Astronight| at this site has no dark time for observing.

    SkyCoordNotScalarError
        Raised when a SkyCoord object is given but a scalar SkyCoord object is required.

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
                 sun_dark_alt: float | None = None):
        # Handle inputs:
        self.site = site
        self.site_name = site.name
        self.an_date = AN_date(an_date)
        if sun_dark_alt is None:
            self.sun_altitude_dark = site.sun_altitude_dark
        else:
            self.sun_altitude_dark = sun_dark_alt

        # Select the almanac engine to use (the *ONLY* code line allowed to reference
        #   a subclass of the abstract class AlmanacEngine):
        self.engine = SkyfieldEngine.from_site(site=site,
                                               sun_altitude_dark=self.sun_altitude_dark)

        # Get utc_approx_midnight from an_date and longitude, to use as starting point:
        approx_midnight_utc = Time(datetime(self.an_date.year, self.an_date.month,
                                            self.an_date.day, 0, 0, 0), scale='utc') +\
            TimeDelta(-site.longitude_hours * 3600, format='sec') + \
            TimeDelta(24 * 3600, format='sec')

        # Get time when sun is farthest below horizon:
        self.sun_antitransit_utc = self.engine.sun_antitransit_utc(approx_midnight_utc)
        _, sun_antitransit_alt = self.engine.sun_azalt(time=self.sun_antitransit_utc)
        if sun_antitransit_alt > HORIZON_USNO:
            raise SunAlwaysUpError(f"Astronight date {self.an_date.an_str}, "
                                   f"latitude {site.latitude:.2f}")
        if sun_antitransit_alt > self.sun_altitude_dark:
            raise NoDarkTimeError(f"Astronight date {self.an_date.an_str}, "
                                  f"latitude {site.latitude:.2f}")

        # No-sun (sun below horizon) timespan:
        self.sunset_utc = \
            self.engine.prev_sunset_utc(ref_time_utc=self.sun_antitransit_utc)
        self.sunrise_utc = \
            self.engine.next_sunrise_utc(ref_time_utc=self.sun_antitransit_utc)
        self.timespan_no_sun = Timespan(self.sunset_utc, self.sunrise_utc)

        # Calculate dark (sun sufficiently below horizon) timespan,
        # but limited to 12 hours to each side of lowest sun angle:
        self.dark_start_utc = \
            max(self.engine.prev_dark_start_utc(ref_time_utc=self.sun_antitransit_utc),
                self.sun_antitransit_utc - timedelta(hours=12))
        self.dark_end_utc = \
            min(self.engine.next_dark_end_utc(ref_time_utc=self.sun_antitransit_utc),
                self.sun_antitransit_utc + timedelta(hours=12))
        self.timespan_dark = Timespan(self.dark_start_utc, self.dark_end_utc)

        # Moon quantities for this Astronight:
        self.moon_illumination = \
            self.engine.moon_illumination(time=self.sun_antitransit_utc)
        self.moon_skycoord = \
            self.engine.moon_skycoord(time=self.sun_antitransit_utc)
        self.moon_transit_utc = \
            self.engine.moon_transit_time(ref_time_utc=self.sun_antitransit_utc)
        moonsets_moonrises = self.engine.moonsets_and_moonrises(
            self.timespan_dark.start,
            self.timespan_dark.end)
        self.timespan_dark_no_moon = self._calc_dark_no_moon(moonsets_moonrises)

    def moon_distance(self, target: SkyCoord) -> float:
        """ Returns distance of target from moon's position at sun antitransit,
        in degrees."""
        return angular_separation(self.moon_skycoord.ra, self.moon_skycoord.dec,
                                  target.ra, target.dec).to(u.deg).value

    def target_transit_utc(self, target: SkyCoord) -> Time:
        """ Returns astropy Time at which target crosses meridian. """
        return self.engine.target_transit_utc(target, self.timespan_dark.midpoint)

    def target_observable(self, target: SkyCoord, min_alt: float,
                          min_moon_dist: float) -> Timespan:
        """ Returns timespan within .timespan_dark during which
        target is observable (i.e., above min_alt in degrees,
        which is typically about 30). """
        moon_distance = self.moon_distance(target)
        if moon_distance < min_moon_dist:
            return Timespan(self.sun_antitransit_utc, self.sun_antitransit_utc)

        target_transit_utc = self.engine.target_transit_utc(target,
                                                            self.sun_antitransit_utc)
        target_rise_utc = self.engine.prev_target_rise_utc(target,
                                                           target_transit_utc, min_alt)
        if target_rise_utc is None:
            return None
        target_set_utc = self.engine.next_target_set_utc(target,
                                                         target_rise_utc, min_alt)
        timespan_target_up = Timespan(target_rise_utc, target_set_utc)
        timespan_target_observable = self.timespan_dark.intersection(timespan_target_up)
        return timespan_target_observable

    def time_from_hhmm(self, hhmm_string: str) -> Time:
        """ Returns astropy Time object corresponding to UTC time 'hhmm' or 'hh:mm'
            and occurring during this Astronight (or at least, nearest the night's
            time of sun antitransit).
        """
        hhmm = hhmm_string.replace(':', '')  # gives, e.g., '1234'
        if len(hhmm) != 4:
            raise ValueError(f'time_from_hhmm() cannot parse \'{hhmm_string}\'.')
        hhcmm = ':'.join([hhmm[0:2], hhmm[2:4]])  # gives, e.g., '12:34'
        central_candidate_time = \
            Time(' '.join([self.sun_antitransit_utc.iso.split()[0], hhcmm]))
        candidate_times = [central_candidate_time + TimeDelta(n * u.d)
                           for n in range(-2, 3)]  # i.e., 2-3 days before and after.
        return nearest_time(candidate_times, self.sun_antitransit_utc)

    def times_at_sun_alt(self, sun_alt: float) -> Tuple[Time | None, Time | None]:
        """ Returns list of two astropy Times during current Astronight
            at which sun is at the given sun altitude. If either time doesn't exist,
            returns None in that place in the list; normally zero place or two places
            will be None, only one being None would be very rare but possible.
        """
        prev_time = self.engine.prev_sunset_utc(self.sun_antitransit_utc, sun_alt)
        next_time = self.engine.next_sunrise_utc(self.sun_antitransit_utc, sun_alt)
        return prev_time, next_time

    @property
    def acp_header_string(self) -> List[str]:
        """ Return list of text lines (strings) suitable for helping form header
            text for ACP plan files and Astronight summary files.
            USAGE: lines = an.acp_header_string """
        lines = []

        # Sun line:
        dark_start_sidereal_hours = \
            self.engine.local_sidereal_time(self.dark_start_utc) % 24
        dark_start_sidereal_hhmm = ra_as_hours(15 * dark_start_sidereal_hours) \
            .rsplit(':', maxsplit=1)[0].replace(':', '')
        dark_end_sidereal_hours = self.engine.local_sidereal_time(self.dark_end_utc)
        if dark_end_sidereal_hours < 0:
            dark_end_sidereal_hours += 24.0
        dark_end_sidereal_hhmm = ra_as_hours(15 * dark_end_sidereal_hours) \
            .rsplit(':', maxsplit=1)[0].replace(':', '')
        sun_line = f'; sun --- down: {hhmm(self.sunset_utc)}-'\
                   f'{hhmm(self.sunrise_utc)} UTC,   '\
                   f'dark({round(self.sun_altitude_dark):+2d}\N{DEGREE SIGN}): '\
                   f'{hhmm(self.dark_start_utc)}-{hhmm(self.dark_end_utc)} UTC  = '\
                   f'{dark_start_sidereal_hhmm}-{dark_end_sidereal_hhmm} LST'
        lines.append(sun_line)

        # Moon line:
        if self.timespan_dark_no_moon is None or self.timespan_dark_no_moon.seconds < 1:
            dark_no_moon_string = 'Moon UP ALL NIGHT'
        elif self.timespan_dark_no_moon == self.timespan_dark:
            dark_no_moon_string = 'Moon DOWN ALL NIGHT'
        else:
            dark_no_moon_string = \
                f'no moon: {hhmm(self.timespan_dark_no_moon.start)}'\
                f'-{hhmm(self.timespan_dark_no_moon.end)} UTC'
        moon_line = f'; moon -- {round(100.0 * self.moon_illumination)}% '\
                    f'({self.moon_skycoord.ra.hour:.1f}h,'\
                    f'{int(round(self.moon_skycoord.dec.degree)):+d}' \
                    f'\N{DEGREE SIGN})   ' + \
                    dark_no_moon_string + '    '\
                    f'transit: {hhmm(self.moon_transit_utc)}'
        lines.append(moon_line)

        # LST line:
        lst_hours_at_sun_at = \
            self.engine.local_sidereal_time(self.sun_antitransit_utc)
        utc_as_datetime = self.sun_antitransit_utc.to_datetime()
        utc_hours_at_sun_at = utc_as_datetime.hour + \
            utc_as_datetime.minute / 60 + utc_as_datetime.second / 3600 + \
            utc_as_datetime.microsecond / (3600 * 1000000)
        lst_minus_utc_hours = (lst_hours_at_sun_at - utc_hours_at_sun_at) % 24.0
        utc_minus_lst_hours = (-lst_minus_utc_hours) % 24.0
        lst_minus_utc_minutes = int(round(60.0 * lst_minus_utc_hours))
        utc_minus_lst_minutes = int(round(60.0 * utc_minus_lst_hours))
        lst_minus_utc_hm = divmod(lst_minus_utc_minutes, 60.0)
        utc_minus_lst_hm = divmod(utc_minus_lst_minutes, 60.0)
        lst_minus_utc_hhmm = f'{int(lst_minus_utc_hm[0]):02d}{int(round(lst_minus_utc_hm[1])):02d}'
        utc_minus_lst_hhmm = f'{int(utc_minus_lst_hm[0]):02d}{int(round(utc_minus_lst_hm[1])):02d}'
        lst_line = f'; LST = UTC + {lst_minus_utc_hhmm}     '\
                   f'UTC = LST + {utc_minus_lst_hhmm}    '\
                   f'@sun antitransit = {hhmm(self.sun_antitransit_utc)} UTC'
        lines.append(lst_line)

        return lines

    def _calc_dark_no_moon(self, moonsets_moonrises: List[Tuple]) -> Timespan | None:
        """ Returns (longest if > 1) contiguous dark-and-no-moon timespan. """
        timespan_edge_times = \
            [self.timespan_dark.start] + \
            [t for (t, value) in moonsets_moonrises] + \
            [self.timespan_dark.end]
        timespans = [Timespan(timespan_edge_times[i], timespan_edge_times[i + 1])
                     for i in range(len(timespan_edge_times) - 1)]
        longest_moon_down_timespan = None
        longest_duration = 0  # in seconds
        for ts in timespans:
            if ts.duration > 0 * u.s:
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


def calc_phase_angle_bisector(times: List[Time], mp_skycoords: List[SkyCoord],
                              deltas: List[float], site: Site) -> List[SkyCoord]:
    """Calculate and return phase angle bisectors at time for a minor planet.
    :param times: list of times [list of astropy Time objects].
    :param mp_skycoords: list of MP sky locations corresponding to the times.
           [list of SkyCoord objects]
    :param deltas: list of distances site-to-MP, in AU. [list of floats]
    :param site: Site object for observing location. [Site object]
    :return: list of phase angle bisectors. [SkyCoord objects]
    """
    # Internally, all coordinates and vectors are in GCRS frame and
    # in meters (even if dimensionless).
    if len(times) != len(mp_skycoords):
        raise ValueError('parameters \'times\' and \'mp_skycoords\' '
                         'must agree in length.')

    # Get GCRS location of sun center, as ndarray [x,y,z] in meters:
    sun_xyz_gcrs = [get_sun(t).gcrs.cartesian.xyz.to(u.m) for t in times]

    # Get GCRS location of observing site, as ndarray [x,y,z] in meters:
    site_loc = EarthLocation.from_geodetic(lon=site.longitude, lat=site.latitude,
                                           height=site.elevation)
    site_xyz_gcrs = [site_loc.get_gcrs_posvel(t)[0].xyz.to(u.m) for t in times]

    # Get GCRS location of MP (NB: Delta is from MP to site (not to geocenter),
    # as ndarray [x,y,z] in meters:
    mp_xyz_from_geocenter = [SkyCoord(ra=sc.ra, dec=sc.dec, distance=delta * u.au,
                                      frame='gcrs').gcrs.cartesian.xyz.to(u.m)
                             for (sc, delta) in zip(mp_skycoords, deltas)]
    mp_xyz_gcrs = [(s[0] + m[0], s[1] + m[1], s[2] + m[2])
                   for (s, m) in zip(site_xyz_gcrs, mp_xyz_from_geocenter)]

    # Make vectors from MP to sun, MP to observing site, sun to site:
    vector_mp_to_sun = [(sun[0] - mp[0], sun[1] - mp[1], sun[2] - mp[2])
                        for (sun, mp) in zip(sun_xyz_gcrs, mp_xyz_gcrs)]
    vector_mp_to_site = [(site[0] - mp[0], site[1] - mp[1], site[2] - mp[2])
                         for (site, mp) in zip(site_xyz_gcrs, mp_xyz_gcrs)]
    vector_sun_to_site = [(site[0] - sun[0], site[1] - sun[1], site[2] - sun[2])
                          for (site, sun) in zip(site_xyz_gcrs, sun_xyz_gcrs)]

    # Getting vector lengths, renamed per angle bisector theorem:
    v_ab = vector_mp_to_sun
    v_ac = vector_mp_to_site
    v_bc = vector_sun_to_site
    len_ab = [sqrt(v[0].value ** 2 + v[1].value ** 2 + v[2].value ** 2) * u.m
              for v in v_ab]
    len_ac = [sqrt(v[0].value ** 2 + v[1].value ** 2 + v[2].value ** 2) * u.m
              for v in v_ac]
    len_bc = [sqrt(v[0].value ** 2 + v[1].value ** 2 + v[2].value ** 2) * u.m
              for v in v_bc]

    # Solve for point d along bc and on bac's angle bisector:
    fraction = [l_ab / (l_ac + l_ab)
                for (l_ab, l_ac, l_bc) in zip(len_ab, len_ac, len_bc)]
    d_xyz_gcrs = [(sun_loc[0] + f * vss[0],
                   sun_loc[1] + f * vss[1],
                   sun_loc[2] + f * vss[2])
                  for (f, sun_loc, vss) in
                  zip(fraction, sun_xyz_gcrs, vector_sun_to_site)]
    # Make vector ad, convert it to ecliptic coordinates:
    v_ad = [(d[0] - sun[0], d[1] - sun[1], d[2] - sun[2])
            for (sun, d) in zip(mp_xyz_gcrs, d_xyz_gcrs)]
    sc_list = [SkyCoord(x=-v[0], y=-v[1], z=-v[2], unit='m', frame='gcrs',
                        representation_type='cartesian').geocentricmeanecliptic
               for v in v_ad]
    return sc_list


__________CLASS_AN_DATE_________________________________________________________ = 0


class AN_date:
    """Builds, holds, and delivers data based on an Astronight date format."""
    def __init__(self, an_input: str | int) -> None:
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
        an_date = date(year=self.year, month=self.month, day=self.day)
        self.day_of_week = DAYS_OF_WEEK[an_date.weekday()]


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
    def prev_dark_start_utc(self, ref_time_utc: Time,
                            sun_alt_dark: float = None) -> Time:
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
    def next_dark_end_utc(self, ref_time_utc: Time,
                          sun_alt_dark: float = None) -> Time:
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
    def target_transit_utc(self, target: SkyCoord, ref_time_utc: Time) -> Time:
        """ Returns target's transit time nearest the given reference time.

        Parameters
        ----------
        target : |SkyCoord|
            Sky coordinates of target whose transit is wanted.

        ref_time_utc : |Time|
            Approximate local midnight, UTC.

        Returns
        -------
        nearest_transit_time : |Time|
            Time of target's transit nearest the given reference time.
        """

    @abstractmethod
    def prev_target_rise_utc(self, target: SkyCoord, ref_time_utc: Time,
                             horizon_degrees: float) -> Time:
        """ Return target rise time nearest and previous to given reference time.
        Ref. time is typically the target's transit time (see .target_transit_utc()).
        Rise is defined by horizon_degrees.

        Parameters
        ----------
        target : |SkyCoord|
            Sky coordinates of target whose rise time is wanted.

        ref_time_utc : |Time|
            Approximate local midnight, UTC.

        horizon_degrees : float
            Degrees above ideal horizon defining target rise (usually for determining
            observability), typically about 30.

        Returns
        -------
        prev_rise_time : |Time|
            Time of target's previous time of rise above horizon_degrees and
            nearest the given reference time.
        """

    @abstractmethod
    def next_target_set_utc(self, target: SkyCoord, ref_time_utc: Time,
                            horizon_degrees: float) -> Time:
        """ Return target set time nearest and after given reference time.
        Ref. time is typically the target's transit time (see .target_transit_utc()).
        Setting is defined by horizon_degrees.

        Parameters
        ----------
        target : |SkyCoord|
            Sky coordinates of target whose set time is wanted.

        ref_time_utc : |Time|
            Approximate local midnight, UTC.

        horizon_degrees : float
            Degrees above ideal horizon defining target set (usually for determining
            observability), typically about 30.

        Returns
        -------
        next_set_time : |Time|
            Time of target's next time of setting below horizon_degrees and
            nearest the given reference time.
        """

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
        sf_time = self._astropy_time_to_skyfield(time)
        return self._target_azalt(sf_target, sf_time)  # tuple (az, alt)

    def prev_sunset_utc(self, ref_time_utc: Time, horizon_degrees=HORIZON_USNO)\
            -> Time | None:
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
        if len(crossings_list) == 0:
            return None
        # Select only crossings that are sunsets:
        sunset_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                        for (sf_time, value) in crossings_list
                        if value == SkyfieldEngine._HorizonCrossing.SETTING]
        return max(sunset_times)

    def next_sunrise_utc(self, ref_time_utc: Time,
                         horizon_degrees=HORIZON_USNO) -> Time | None:
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
        if len(crossings_list) == 0:
            return None
        # Select only crossings that are sunrises:
        sunrise_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                         for (sf_time, value) in crossings_list
                         if value == SkyfieldEngine._HorizonCrossing.RISING]
        return min(sunrise_times)

    def prev_dark_start_utc(self, ref_time_utc: Time,
                            sun_alt_dark: float = None) -> Time:
        """ Return going-dark time nearest and previous to given reference time,
        based on angle of sun below horizon that user has defined as twilight in
        sun_altitude_dark.
        """
        if sun_alt_dark is None:
            sun_alt_dark = self.sun_alt_dark
        return self.prev_sunset_utc(ref_time_utc=ref_time_utc,
                                    horizon_degrees=sun_alt_dark)

    def next_dark_end_utc(self, ref_time_utc: Time,
                          sun_alt_dark: float = None) -> Time:
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
        # TODO: use util.nearest_time() here:
        diffs = [abs(t - ref_time_utc).sec for t in antitransit_times]
        nearest_antitransit_time = antitransit_times[diffs.index(min(diffs))]
        return nearest_antitransit_time

    def target_transit_utc(self, target: SkyCoord, ref_time_utc: Time) -> Time:
        """Return target's transit time nearest the given time.
        *** See abstract method for parameters and return value.
        """
        target_sf = self._skycoord_to_skyfield_star(target)
        time_sf = self.timescale.from_astropy(ref_time_utc)  # skyfield Time obj.
        time_start = time_sf - timedelta(hours=15)
        time_end = time_sf + timedelta(hours=15)
        meridian_crossings = self._meridian_crossings(target_sf,
                                                      time_start, time_end)
        transit_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                         for (sf_time, value) in meridian_crossings
                         if value == SkyfieldEngine._MeridianCrossing.TRANSIT]
        # TODO: use util.nearest_time() here:
        diffs = [abs(t - ref_time_utc).sec for t in transit_times]
        nearest_transit_time = transit_times[diffs.index(min(diffs))]
        return nearest_transit_time

    def prev_target_rise_utc(self, target: SkyCoord, ref_time_utc: Time,
                             horizon_degrees: float) -> Time | None:
        """ Return target rise time nearest and previous to given reference time.
        Ref. time is typically the target's transit time (see .target_transit_utc())."""
        target_sf = self._skycoord_to_skyfield_star(target)
        sf_ref_time = self.timescale.from_astropy(ref_time_utc)  # skyfield Time obj.
        sf_start_time = sf_ref_time - timedelta(hours=36)
        sf_end_time = sf_ref_time
        crossings_list = self._horizon_crossings(target_sf,
                                                 sf_start_time, sf_end_time,
                                                 horizon_degrees=horizon_degrees)
        if len(crossings_list) == 0:
            return None
        rise_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                      for (sf_time, value) in crossings_list
                      if value == SkyfieldEngine._HorizonCrossing.RISING]
        return max(rise_times)

    def next_target_set_utc(self, target: SkyCoord, ref_time_utc: Time,
                            horizon_degrees: float) -> Time:
        """ Return target set time nearest and after the given reference time.
        Ref. time is typically the target's transit time (see .target_transit_utc())."""
        target_sf = self._skycoord_to_skyfield_star(target)
        if ref_time_utc is None:
            iiii = 4
        sf_ref_time = self.timescale.from_astropy(ref_time_utc)  # skyfield Time obj.
        sf_start_time = sf_ref_time
        sf_end_time = sf_ref_time + timedelta(hours=36)
        crossings_list = self._horizon_crossings(target_sf,
                                                 sf_start_time, sf_end_time,
                                                 horizon_degrees=horizon_degrees)
        if len(crossings_list) == 0:
            return None
        set_times = [SkyfieldEngine._skyfield_time_to_astropy(sf_time)
                     for (sf_time, value) in crossings_list
                     if value == SkyfieldEngine._HorizonCrossing.SETTING]
        return min(set_times)

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


