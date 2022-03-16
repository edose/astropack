""" Module astropak.almanac:
    Astronight class + other almanac code.
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timezone, timedelta
from math import sin, cos, asin, acos, sqrt

# External packages:
from numpy import diff, flatnonzero, clip
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from skyfield.api import load, wgs84, Star
from skyfield.almanac import risings_and_settings, fraction_illuminated, meridian_transits
from skyfield.searchlib import find_discrete

# Author's packages:
from astropack.util import Timespan, ra_as_hours
from .reference import DEGREES_PER_RADIAN

THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')

HORIZON_USNO = -0.833  # USNO sun & moon effective horizon (radius & refraction), in degrees.


class SunAlwaysUpError(Exception):
    """Raised when sun-down time is required for observing but none exists that night.

    User is responsible for catching this, if there is a chance of the sun's being up for
    24 consecutive hours (only in earth's polar regions).
    """
    pass


class NoDarkTimeError(Exception):
    """Raised when dark time is required for observing but none exists that night.

    User is responsible for catching this, if there is a chance of daylight or twilight
    extending for 24 consecutive hours (typically only at earth's polar latitudes of about
    > +55 or <-55 degrees.
    """
    pass


__________CLASS_ASTRONIGHT_______________________________________________________ = 0


class Astronight:
    """ Container for observation-relevant data for one observing night at one site.
    New in astropy2022: Removed references to ephem package. Constructor uses new .ini file structure
    (see astropy.ini module).

    Examples
    --------
    One supplies a site object and an "astronight" date (yyyymmdd for local date that the nighttime starts),
    and Astronight generates and supplies numerous properties and methods, mostly useful for planning
    observations during that night. Format for the astronight date differs from ISO 8601 format to emphasize
    sharply the difference in meaning between the two date types.

      >>> from astropack.almanac import Astronight
      >>> from astropack.ini import Site
      >>> this_site = Site('My Dome')
      >>> this_date = 20220302  # or '20220302'
      >>> an = Astronight(this_site, this_date)

    The Astronight object then supplies relevant data about the local night.

      >>> mp = an.moon_illumination
      >>> dark = an.timespan_dark               # returns astropak.util.Timespan object
      >>> very_dark = an.timespan_dark_no_moon  # "

    The Astronight object can also supply target-observation data, typically as astropy.Time or
    astropak.util.Timespan objects.

      >>> from astropy.coordinates import SkyCoord
      >>> target = SkyCoord(10.4532, -25.32, unit="deg")
      >>> transit_time = an.transit(target)
      >>> available_timespan = an.timespan_observable(target, min_alt=30, min_moon_dist=50)

    If the site has no nighttime for the given date (summer at polar regions), Astronight will raise
    a SunAlwaysUpError exception, which the user is responsible for handling.

    If the site has no dark time (sun sufficiently below horizon) for the given date (summer at
    high or low latitudes), Astronight will raise a NoDarkTimeError exception,
    which the user is responsible for handling.

    Attributes
    ----------
    .site : the site object passed in. [astropak.ini.Site object]
    .site_name : the name of the site object passed in. [string]
    .an_date_string : yyyymmdd string representation of the Astronight date. [string, len=8]
    .sun_altitude_dark : the maximum sun altitude, in degrees, considered observably dark. [float]
    .obs : the skyfield observer object. [skyfield.Observer object]
    .timespan_no_sun : timespan when sun is down at this site. [astropak.util.Timespan object]
    .timespan_dark : timespan when sky is observably dark at this site. [astropak.util.Timespan object]
    .local_middark_utc : midpoint of dark timespan. [astropy.time.Time object]
    .local_middark_lst_hour_string : hex string representing local sidereal time in hours:minutes:seconds
        [string hh:mm:ss]
    .moon_illumination : percent illumination of moon at local middark, 0-100. [float]
    .moon_ra : moon right ascension at local middark, in degrees, 0-360. [float]
    .moon_dec : moon declination at local middark, in degrees, -90 - +90. [float]
    .moon_skycoord : SkyCoord object holding moon sky location at local middark.
        [astropy.coordinates.SkyCoord object]
    .moon_transit : the time (i.e., the nearest local middark) at which moon crosses meridian.
        [astropy.time.Time object]
    .moon_up_timespan : timespan when moon is above horizon,
        limited by .timespan_no_sun. List of two
        timespans if moon is up at sunset and at sunrise. [astropak.util.Timespan object, or list of 2]

    Parameters
    ----------
    site : `~astropak.ini.Site` class
        Observing location.
    an_date : int or string of exactly 8 digits or characters, respectively
        Local date (yyyymmdd) when the night starts.
    sun_altitude_dark : float, optional
        Altitude of sun's center, in degrees, when it is considered
        sufficiently dark to begin observations.
        If absent, uses site object's Sun Altitude Dark entry (as read from its .ini file).
    """

    def __init__(self, site, an_date, sun_altitude_dark=None):
        # Handle inputs:
        self.site = site
        self.site_name = site.name
        self.an_date_string = str(an_date)
        if len(self.an_date_string) != 8:
            raise ValueError('Parameter \'.an_date_string\' appears invalid (must be yyyymmdd).')
        if sun_altitude_dark is None:
            self.sun_altitude_dark = site.sun_altitude_dark
        else:
            self.sun_altitude_dark = sun_altitude_dark

        # Initiate skyfield package:
        self.master_eph = load('de440s.bsp')
        self.timescale = load.timescale()
        self.obs = make_skyfield_observatory_from_site(site)  # Skyfield observer.

        # Get UTC date and *approximate* local midnight time, as starting point:
        an_year = int(self.an_date_string[0:4])
        an_month = int(self.an_date_string[4:6])
        an_day = int(self.an_date_string[6:8])
        if site.utc_offset < 0:
            # West of Greenwich but east of Intl Date Line [Americas].
            longitude_hours = ((site.longitude % 360) - 360.0) / 15.0
        else:
            # East of Greenwich but west of Intl Date Line [Eurasia, Africa, ANZ].
            longitude_hours = (site.longitude % 360) / 15.0
        approx_midnight = \
            Time(datetime(an_year, an_month, an_day, 0, 0, 0).replace(tzinfo=timezone.utc)) + \
            TimeDelta(-longitude_hours * 3600, format='sec') + TimeDelta(24 * 3600, format='sec')

        # No-sun (sun below horizon) timespan:
        self.timespan_no_sun = find_no_sun(self.obs, self.master_eph, self.timescale, approx_midnight, 12)
        if self.timespan_no_sun is None:
            raise SunAlwaysUpError('Astronight date ' + self.an_date_string +
                                   ', latitude=' + '{0:.2f}'.format(site.latitude))

        # Dark (sun far below horizon) timespan:
        self.timespan_dark = find_target_up_down(self.obs, self.master_eph, self.master_eph['sun'],
                                                 self.timescale, self.timespan_no_sun, 'down',
                                                 horizon=self.sun_altitude_dark)
        if self.timespan_dark is None:
            raise NoDarkTimeError('Astronight date ' + self.an_date_string +
                                  ', latitude=' + '{0:.2f}'.format(site.latitude))
        self.timespan_dark = self.timespan_dark[0]

        self.local_middark_utc = self.timespan_dark.midpoint
        local_middark_lst_degrees = 15.0 * local_sidereal_time(site.longitude, self.timescale,
                                                               self.local_middark_utc)
        self.local_middark_lst_hour_string = ra_as_hours(local_middark_lst_degrees)

        # Moon quantities and moon_down timespan:
        self.moon_illumination = moon_illumination_pct(self.master_eph, self.timescale,
                                                       self.local_middark_utc)
        self.moon_ra, self.moon_dec = moon_ra_dec(self.obs, self.master_eph,
                                                  self.timescale, self.local_middark_utc)
        self.moon_skycoord = SkyCoord(self.moon_ra, self.moon_dec, unit="deg")
        self.moon_transit = moon_transit_time(self.obs, self.master_eph,
                                              self.timescale, self.local_middark_utc)
        self.moon_up_timespans = find_target_up_down(self.obs, self.master_eph,
                                                     self.master_eph['moon'], self.timescale,
                                                     self.timespan_dark, 'up', HORIZON_USNO)

        # Dark_no_moon timespan:
        dnm = self.timespan_dark.copy()
        for mdt in self.moon_up_timespans:
            dnm = dnm.subtract(mdt)
        self.timespan_dark_no_moon = dnm

    def transit(self, target_skycoord):
        """ Return transit UTC [py datetime], closest to middark time, for given target. """
        # TODO: rewrite to pass through a SkyCoord object rather than list of Star objects.
        return target_transit_time(self.obs, self.master_eph, self.timescale, target_skycoord,
                                   self.local_middark_utc)

    # TODO: ensure that this works for a LIST of SkyCoords or array SkyCoords,
    #     as well as for a single SkyCoord.
    def timespan_observable(self, target_skycoord, min_alt=None, min_moon_dist=None):
        """ Return observable timespan [Timespan object] for given target.
            :param target_skycoord: target's sky position, in astropy.coordinates.SkyCoord object
            :param min_alt: minimum altitude for target to be observable, in degrees (as int, float,
                or as astropy Quantity)
            :param min_moon_dist: minimum distance from the moon for target to be observable, in degrees
                (as int, float, or as astropy Quantity)
            :return: Timespan for which target is observable.
        """
        if min_alt is None or min_moon_dist is None:
            raise ValueError('Astropy.timespan_observable requires values for min_alt and min_moon_dist.')
        # moon = self.master_eph['moon']
        # middark_sf = self.timescale.from_astropy(self.local_middark_utc)

        star_list = _skycoords_to_skyfield_star_list(target_skycoord)
        timespan_target_up = find_target_up_down(self.obs, self.master_eph, star_list, self.timescale,
                                                 self.timespan_dark, 'up', min_alt)
        moon_dist = self.moon_ra.separation_from(star_list)
        if moon_dist >= min_moon_dist:
            return timespan_target_up
        else:
            return [time_up.intersection(self.timespan_dark_no_moon) for time_up in timespan_target_up]

    def __repr__(self):
        return 'astropak.almanac.Astronight(site=[' + self.site.name + ']'\
               ', an_date=' + self.an_date_string + \
               ', sun_altitude_dark=' + str(self.sun_altitude_dark) + ')  # astropak 2022.'

    def __str__(self):
        return "Astronight[2022] '" + self.an_date_string + "' at site '" + self.site_name + "'."


__________PUBLIC_FUNCTIONS______________________________________________________ = 0


def ha_dec_from_az_alt(latitude, az_alt):
    """Return (hourangle, declination) from known site latitude and target (azimuth, altitude).
       May pass in a single (az, alt) tuple, or a list of such tuples.
    :param latitude: site latitude, in degrees. [float]
    :param az_alt: 2-tuple (azimuth, altitude), in degrees [(float, float)] or list of such tuples.
    :return: target (hourangle, declination), in degrees. [2-tuple of (float, float)]
    """
    az_alt_is_list = isinstance(az_alt, list)
    if not az_alt_is_list:
        az_alt = [az_alt]

    # Formulae (1) and (2) from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm:
    # phi=latitude, delta=declination, H=hourangle, A=target azimuth, a=target altitude.
    cos_phi, sin_phi = cos(latitude / DEGREES_PER_RADIAN), sin(latitude / DEGREES_PER_RADIAN)
    results = []
    for (az, alt) in az_alt:
        cos_a, sin_a = cos(alt / DEGREES_PER_RADIAN), sin(alt / DEGREES_PER_RADIAN)
        cos_A = cos(az / DEGREES_PER_RADIAN)

        # (1) sin(δ) = sin(a) sin(φ) + cos(a) cos(φ) cos(A)
        sin_delta = (sin_a * sin_phi) + (cos_a * cos_phi * cos_A)
        cos_delta = sqrt(1.0 - sin_delta ** 2)  # cosine of declination is always non-negative.
        sin_delta = clip(sin_delta, -1.0, 1.0)  # prevent floating-point errors driving abs(sin(delta)) > 1
        delta = asin(sin_delta) * DEGREES_PER_RADIAN

        # (2) cos(H) = { sin(a) - sin(δ) sin(φ) } / { cos(δ) cos(φ) }
        cos_H = (sin_a - sin_delta * sin_phi) / (cos_delta * cos_phi)  # but sign of H is ambiguous.
        cos_H = clip(cos_H, -1.0, 1.0)  # prevent floating-point errors driving abs(cos(H)) > 1
        sign_H = 1.0 if 0.0 <= az <= 180.0 else -1.0
        H = sign_H * acos(cos_H) * DEGREES_PER_RADIAN
        results.append((H, delta))
    if not az_alt_is_list:
        results = results[0]
    return results


__________SUPPORT_FUNCTIONS______________________________________________________ = 0


def make_skyfield_observatory_from_site(site):
    """ Return Skyfield observer (earth location) object from astropak Site object. """
    obs = wgs84.latlon(site.latitude, site.longitude, site.elevation)
    return obs


def find_no_sun(obs, master_eph, timescale, approx_midnight, hours_each_side=12):
    """Return Timespan from sunset to sunrise within given Timespan.
       Special case required to handle 'nights' with no rise and/or no set.
    :param obs: observer earth location. [skyfield object]
    :param master_eph: skyfield master ephemeris. [skyfield object]
    :param timescale: skyfield timescale. [skyfield object]
    :param approx_midnight: approximate midnight. [astropy Time object]
    :param hours_each_side: hours on each side of midnight sufficient to capture
        both sunset and sunrise. [float]
    :return: Timespan from sunset to sunrise. [Timespan]
    """
    # print(approx_midnight)
    # print(type(approx_midnight))
    midnight = timescale.from_astropy(approx_midnight)  # skyfield Time obj.
    time_before = midnight - timedelta(hours=hours_each_side)
    time_after = midnight + timedelta(hours=hours_each_side)
    sun_fn = risings_and_settings(master_eph, master_eph['sun'], obs, HORIZON_USNO)
    sun_fn.step_days = 1.0 / 12.0
    set_times, set_values = find_discrete(time_before, midnight, sun_fn)
    rise_times, rise_values = find_discrete(midnight, time_after, sun_fn)
    topos_at = (master_eph['earth'] + obs).at
    alt_midnight = topos_at(midnight).observe(master_eph['sun']).apparent().altaz()[0].degrees
    # Case: no rises or sets, sun is either up or down for the entire search period.
    if len(set_values) <= 0 and len(rise_values) <= 0:
        if alt_midnight <= HORIZON_USNO:
            return Timespan(_skyfield_time_to_astropy(time_before), _skyfield_time_to_astropy(time_after))
        else:
            return None
    if len(set_values) >= 1:
        set_time = _skyfield_time_to_astropy(set_times[-1])
    else:
        set_time = _skyfield_time_to_astropy(time_before)
    if len(rise_values) >= 1:
        rise_time = _skyfield_time_to_astropy(rise_times[0])
    else:
        rise_time = _skyfield_time_to_astropy(time_after)
    return Timespan(set_time, rise_time)


def find_target_up_down(obs, master_eph, target_eph, timescale, timespan, up_down='up', horizon=0.0):
    """Return timespan for which target is either up or down (above or below given horizon),
       within a given timespan.
    :param obs: skyfield observatory object.
    :param master_eph: skyfield Ephmeris object.
    :param target_eph: either a solar-system ephemeris or a Star object.
    :param timescale: skyfield Timescale object.
    :param timespan:
    :param up_down:
    :param horizon:
    :return: list of timespans for which target is either up or down. If up_down is satisfied throughout
        input timespan, return list containing input timespan only. If up_down is never satisfied during
        input timespan, return None. [list of Timespan objects (always a list even if only one Timespan),
        or None]
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
        if (up_down == 'up' and alt > horizon) or (up_down == 'down' and alt <= horizon):
            return [timespan.copy()]
        else:
            return None
    # Case: at least one rising or setting found:
    # First, convert times to astropy, so that start and end are passed through exactly.
    astropy_times = [timespan.start] + \
                    [_skyfield_time_to_astropy(t) for t in event_times] + \
                    [timespan.end]
    values = [0.5] + list(event_values) + [0.5]
    diffs = diff(values)
    # Now choose events as up or down:
    if up_down == 'up':
        is_index_qualifying = diffs < 0
    else:
        is_index_qualifying = diffs > 0
    qualifying_indices = flatnonzero(diffs)[is_index_qualifying]
    timespans = [Timespan(astropy_times[i], astropy_times[i + 1]) for i in qualifying_indices]
    return timespans


def moon_ra_dec(obs, master_eph, timescale, time):
    """Return moon sky location, in degrees, at given earth location and time."""
    time_sf = timescale.from_astropy(time)
    ra, dec, _ = (master_eph['earth'] + obs).at(time_sf).observe(master_eph['moon']).apparent().radec()
    return ra.hours * 15.0, dec.degrees


def moon_transit_time(obs, master_eph, timescale, time):
    """Return moon's transit time closest to given time."""
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    time_start = time_sf - timedelta(hours=14)
    time_end = time_sf + timedelta(hours=14)
    moon_fn = meridian_transits(master_eph, master_eph['moon'], obs)
    times, values = find_discrete(time_start, time_end, moon_fn)
    return _skyfield_time_to_astropy(_closest_transit(times, values, time_sf))


def target_transit_time(obs, master_eph, timescale, target_skycoord, time):
    """Return fixed target's (SkyCoord) transit time [astropy Time]
    closest to target astropy Time."""
    star_list = _skycoords_to_skyfield_star_list(target_skycoord)
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    time_start = time_sf - timedelta(hours=14)
    time_end = time_sf + timedelta(hours=14)
    transit_list = []
    for star in star_list:
        star_fn = meridian_transits(master_eph, star, obs)
        times, values = find_discrete(time_start, time_end, star_fn)
        transit_list.append(_skyfield_time_to_astropy(_closest_transit(times, values, time_sf)))
    if len(transit_list) == 1:
        return transit_list[0]  # return single time as scalar Time rather than list.
    return transit_list


def _closest_transit(times, values, time_sf):
    """For times [list of skyfield Time objects], values [list of ints], and
       a target time [Skyfield Time object], determine which transit is closes to target time,
       return the time as a Skyfield Time object."""
    transit_times = [t for (t, v) in zip(times, values) if v == 1]
    timedeltas = [abs(time - time_sf) for time in transit_times]
    transit_time = transit_times[timedeltas.index(min(timedeltas))]
    return transit_time


def local_sidereal_time(longitude, timescale, time):
    """Return local sidereal time in hours for astropy Time at obs (earth location)."""
    from skyfield.earthlib import sidereal_time
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    gmst = sidereal_time(time_sf)  # in hours
    local_offset = longitude / 15.0  # in hours
    return gmst + local_offset


def moon_illumination_pct(master_eph, timescale, time):
    """Return moon illumination fraction at astropy time. """
    time_sf = timescale.from_astropy(time)  # skyfield Time obj.
    fraction = 100.0 * fraction_illuminated(master_eph, 'moon', time_sf)
    return fraction


def _skyfield_time_to_astropy(skyfield_time):
    """Convert skyfield time to astropy time in UTC timescale and ISO format. Time arrays are OK."""
    astropy_time = skyfield_time.to_astropy().utc
    astropy_time.format = 'iso'
    return astropy_time


def _skycoords_to_skyfield_star_list(skycoord):
    """Convert astropy SkyCoord object or list of objects
       to *list of* skyfield Star objects (always returns a list, only one sky position per Star).
       :param skycoord: astropy SkyCoord object, or list of objects, any of which can be array objects
           which will be flattened and converted one at a time.
       :return: *list* of skyfield Star objects.
       """
    def stars_from_skycoord_object(skycoord_obj):
        """Nested function, always returns a list of Star objects."""
        if skycoord_obj.size == 1:
            return [Star(ra_hours=skycoord_obj.ra.degree / 15.0, dec_degrees=skycoord_obj.dec.degree)]
        else:
            return [Star(ra_hours=sc_element.ra.degree / 15.0, dec_degrees=sc_element.dec.degree)
                    for sc_element in skycoord_obj]

    if isinstance(skycoord, SkyCoord):
        return stars_from_skycoord_object(skycoord)
    elif isinstance(skycoord, list):
        star_list = []
        for sc_object in skycoord:
            star_list.extend(stars_from_skycoord_object(sc_object))
        return star_list
