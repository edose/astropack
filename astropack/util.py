""" Module astropak.util
    Utility/support classes and functions. Not specific to observing targets or style.
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timedelta
from math import floor, pow, ceil

# External packages:
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy.stats import circmean


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


_____CLASSES________________________________________________ = 0


class Timespan:
    """ Holds one (start, end) span of time. Immutable.
        Input: 2 python datetimes (in UTC) or astropy.time.Time objects, defining start and end of timespan.
        Internal representation of times: in astropy.time.Time.
        methods:
        ts2 = ts.copy()  # (not too useful, as Timespan objects are immutable.)
        ts2 == ts  # only if both start and end are equal
        ts2 = ts.delay_seconds(120)  # returns new Timespan object, offset in both start and end.
        ts.intersect(other)  # returns True iff any overlap at all with other Timespan object.
        ts2 = ts.remove(other)  # returns new Timespan; longer of 2 possible spans if ambiguous.
        ts.contains_time(t)  # returns True iff ts.start <= t <= ts.end for some datetime object t.
        ts.contains_timespan(other)  # returns True iff ts wholly contains other Timespan.
        Timespan.longer(ts1, ts2)  # returns longer (in duration) of two Timespan objects.
        dt_list = ts.generate_events(jd_ref, period_days, max_events): generates list of up to
            max_events datetimes, all within the Timespan object ts, and the times beginning
            with JD_ref and spaced by period_days.
        str(ts)  # returns string describing Timespan's start, end, and duration in seconds.
    """
    def __init__(self, start_time, end_time):
        """ Constructor.
        :param start_time: start of timespan [py datetime (UTC) or Time]
        :param end_time: end of timespan [py datetime (UTC) or Time]
        """
        if not isinstance(start_time, (datetime, Time)) or not isinstance(end_time, (datetime, Time)):
            raise TypeError('Both inputs must be either datetime or Time objects.')
        if isinstance(start_time, datetime):
            if start_time.tzname() != 'UTC':
                raise ValueError('First input parameter (datetime object) must be in UTC.')
        if isinstance(end_time, datetime):
            if end_time.tzname() != 'UTC':
                raise ValueError('Second input parameter (datetime object) must be in UTC.')
        if isinstance(start_time, Time):
            if not start_time.isscalar:
                raise ValueError('First input parameter (Time object) must be scalar.')
        if isinstance(end_time, Time):
            if not end_time.isscalar:
                raise ValueError('Second input parameter(Time object) must be scalar.')
        self.start = Time(start_time)
        self.end = max(self.start, Time(end_time))
        self.duration = self.end - self.start
        self.seconds = self.duration.sec
        self.days = self.duration.jd
        self.midpoint = self.start + TimeDelta(self.duration / 2)

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def copy(self):
        """ Return new Timespan exactly equivalent to this Timespan."""
        return Timespan(self.start, self.end)

    def delay_by(self, delay):
        """ Return new Timespan with start and end delayed by param 'delay'.
        :param delay: delay to add to returned Timespan, as seconds [int or float], or
            as [astropy.time.TimeDelta], or as [datetime.timedelta].
            Delay may be negative to return earlier Timespan.
        :return: new Timespan delayed by parameter delay. [Timespan]
        """
        if isinstance(delay, (int, float)):
            delay = TimeDelta(delay, format='sec')
        elif isinstance(delay, (TimeDelta, timedelta)):
            delay = TimeDelta(delay)
        else:
            raise TypeError('Parameter \'delay\' must be int, float, py datetime, '
                            'or astropy.TimeDelta.')
        return Timespan(self.start + delay, self.end + delay)

    def expand_by(self, expansion):
        """ Return new Timespan object expanded or contracted at both start and end by param 'expansion'.
        :param expansion: Duration by which to extend both the start and end, as seconds [int or float] or
            as [astropy.time.TimeDelta] or as [datetime.timedelta] .
            May be negative to contract the timespan instead, but will not contract beyond zero duration.
        :return: new Timespan extended (or contracted) at each limit. [Timespan]
        """
        # Use negative seconds to contract Timespan. New Timespan will have non-negative duration.
        if isinstance(expansion, (int, float)):
            delay = TimeDelta(expansion, format='sec')
        elif isinstance(expansion, timedelta):
            delay = TimeDelta(expansion)
        elif isinstance(expansion, TimeDelta):
            delay = expansion
        else:
            raise TypeError('Parameter \'expansion\' must be int, float, py datetime, '
                            'or astropy.TimeDelta.')
        new_start = min(self.start - delay, self.midpoint)
        new_end = max(self.end + delay, self.midpoint)
        return Timespan(new_start, new_end)

    def intersection(self, other):
        """ Return Timespan containing only times common to this Timespan and the param Timespan.
            Return zero-duration Timespan if timespans do not intersect.
            :param other: Timespan object to intersect with this Timespan.
        """
        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)
        return Timespan(new_start, new_end)

    def union(self, other):
        """ Return union of this Timespan and other Timespan, or list of two Timespans if the two do not
            intersect.
        :param other: Timespan to add to this Timespan. [Timespan]
        :return: Timespan containing all times in either timespan, or list of two Timespans if the two do
            not intersect. [Timespan]
        """
        if (self.end < other.start) or (self.start > other.end):
            return [self, other]
        return Timespan(min(self.start, other.start), max(self.end, other.end))

    def add(self, other):
        """ Synonym for self.union(). """
        return self.union(other)

    def subtract(self, other):
        """ Return this timespan with other timespan removed, or list of two timespans
            if other timespan is wholly contained within this timespan.
            No effect if the two timespans do not intersect or if other has zero duration.
            Zero-length Timespan if
        :param other: timespan to remove from this timespan.
        :return: new timespan object of other timespan subtracted from this timespan. [Timespan]
        """
        if (self.intersection(other).seconds == 0) or (other.seconds == 0):
            return self.copy()
        if self.start < other.start and other.end < self.end:
            # Other timespan is wholly within this timespan, so return 2 timespans:
            return [Timespan(self.start, other.start), Timespan(other.end, self.end)]
        if other.start <= self.start and self.end <= other.end:
            # This timespan is wholly within other timespan, so return zero-length timespan at midpoint:
            return Timespan(self.midpoint, self.midpoint)
        if self.start < other.start:
            return Timespan(self.start, other.start)
        return Timespan(other.end, self.end)

    def contains(self, other):
        """ Returns True iff this timespan contains a Time or other timespan.
         :param other: the time or timespan to be tested as wholly contained within this timespan.
             [Time object or other Timespan object]
         :return: True iff other is contained within this timespan, else False. [bool]

         """
        if isinstance(other, (Time, datetime)):
            return self.start <= Time(other) <= self.end
        elif isinstance(other, Timespan):
            return (self.start <= other.start) & (self.end >= other.end)
        raise TypeError('Timespan.contains() requires Time, datetime, or Timespan object as parameter.')

    def split_at(self, time):
        """ Return list of two new timespans if time is within it, else return a copy of this timespan.
        :param time: time at which to split this timespan into two contiguous timespans.
        :return: Two contiguous timespans split at time, or a copy of this timespan if it does not contain
            the given time. [Timespan object, or list of 2 Timespan objects]
        """
        if not isinstance(time, (Time, datetime)):
            raise TypeError('Timespan.split_at() requires Time or datetime object as parameter.')
        time = Time(time)
        if self.start < time < self.end:
            return [Timespan(self.start, time), Timespan(time, self.end)]
        return self.copy()

    @staticmethod
    def longer(ts1, ts2):
        """ Returns Timespan with longer duration (larger .seconds).
            If equal duration, return earlier. If exactly coincident, return ts1.
        :param ts1: input Timespan object.
        :param ts2: input Timespan object.
        :return: the Timespan object with longer duration, or ts1 if durations are equal. [Timespan]
        """
        if ts1.seconds > ts2.seconds:
            return ts1
        elif ts2.seconds > ts1.seconds:
            return ts2
        elif ts1.start <= ts2.start:
            return ts1
        return ts2

    def periodic_events(self, ref_time, period, max_events=10):
        """ Returns a list of UTC times of period events within a given Timespan.
        :param ref_time: Time object of any occurence of the period event. [Time]
        :param period: cycle period of events. [TimeDelta]
        :param max_events: maximum number of events to return in list. [int]
        :return: Time object containing max_events periodic event times within this timespan.
           [Time object holding with list of times]
           Return Time object with no times if this timespan contains no such events.
        """
        if not isinstance(period, (timedelta, TimeDelta)):
            raise TypeError('Period must be a timedelta or TimeDelta object.')
        period = TimeDelta(period)
        if period <= 0.0:
            raise ValueError('Period must be positive.')
        if self.seconds == 0 or max_events <= 0:
            return []
        n_start = ceil((self.start - ref_time).sec / period.sec)
        first_time = ref_time + n_start * period
        if not self.contains(first_time):
            print('\n No events within Timespan.')
            return []
        n_events = min(max_events, 1 + floor(Timespan(first_time, self.end).seconds / period.sec))
        print('\nn_events', str(n_events))
        event_times = [first_time + i * period for i in range(n_events)]
        return event_times

    def __str__(self):
        return "Timespan '" + str(self.start) + "' to '" + str(self.end) + "' = " + \
               '{0:.3f}'.format(self.seconds) + " seconds."

    def __repr__(self):
        return 'Timespan(' + self.start.iso + ', ' + self.end.iso + ')'


__________TIME_and_DATE_FUNCTIONS_____________________________________________ = 0


def hhmm_from_datetime_utc(datetime_utc):
    """ For any datetime (UTC), return string 'hhmm' for the UTC time. TESTS OK 2022-03-11. """
    minutes_of_day = int(round(datetime_utc.hour*60  # NB: banker's rounding (nearest even value, if a tie)
                               + datetime_utc.minute
                               + datetime_utc.second/60
                               + datetime_utc.microsecond/(60*1000000))) % 1440
    hh, mm = divmod(minutes_of_day, 60)
    return '{0:0>4d}'.format(100 * hh + mm)


_____RA_and_DEC_FUNCTIONS_____________________________________ = 0


def ra_as_degrees(ra_string):
    """  Takes Right Ascension as string, returns degrees. TESTS OK 2022-03-11.
    :param ra_string: Right Ascension in either full hex ("12:34:56.7777" or "12 34 56.7777"),
               or degrees ("234.55") [string]
    :return: Right Ascension in degrees between 0 and 360 [float], or None if < 0 or > 360.
    Usage: ra_as_degrees('180.23')    # as degrees from 0 through 360.
           ra_as_degrees('11:16:30')  # as hex, from 0 hours through 24 hours.
    """
    ra_list = parse_hex(ra_string)
    if len(ra_list) == 1:
        ra_degrees = float(ra_list[0])  # input assumed to be in degrees.
    elif len(ra_list) == 2:
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1])/60.0)  # input assumed in hex.
    else:
        ra_seconds = float(ra_list[2].split()[0])  # in case trailing comments
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1]) / 60.0 +
                           ra_seconds/3600.0)  # input assumed in hex.
    if (ra_degrees < 0) | (ra_degrees > 360):
        ra_degrees = None
    return ra_degrees


def hex_as_degrees(hex_degrees_string):
    """ Takes angle in hex degrees string (general case) or degrees, returns degrees as float.
        TESTS OK 2020-10-24.
    :param hex_degrees_string: angle in either full hex ("-12:34:56.7777", or "-12 34 56.7777"),
           or degrees ("-24.55")
    :return degrees. [float]
    """
    # dec_list = hex_degrees_string.split(":")
    dec_list = parse_hex(hex_degrees_string)
    # dec_list = [dec.strip() for dec in dec_list]
    if dec_list[0].startswith("-"):
        sign = -1
    else:
        sign = 1
    if len(dec_list) == 1:
        dec_degrees = float(dec_list[0])  # input assumed to be in degrees.
    elif len(dec_list) == 2:
        dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1])/60.0)  # input is hex.
    else:
        deg_seconds = float(dec_list[2].split()[0])  # in case trailing comments
        dec_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1]) / 60.0 +
                              deg_seconds/3600.0)  # input is hex.
    return dec_degrees


def dec_as_degrees(dec_string):
    """ Takes Declination as string (hex or degrees), returns degrees as float. TESTS OK 2020-10-24.
    :param dec_string: declination in full hex ("-12:34:56.7777") or degrees ("-24.55"). [string]
    :return: degrees, limited to -90 to +90. [float, or None if outside Dec range]
    """
    dec_degrees = hex_as_degrees(dec_string)
    if (dec_degrees < -90) | (dec_degrees > +90):
        dec_degrees = None
    return dec_degrees


def ra_as_hours(ra_degrees, seconds_decimal_places=2):
    """ Takes Right Ascension degrees as float, returns RA string. TESTS OK 2020-10-24.
    :param ra_degrees: Right Ascension in degrees, limited to 0 through 360. [float]
    :param seconds_decimal_places: number of places at end of RA string (no period if zero). [int]
    :return: RA in hours/hex format. [string, or None if outside RA range]
    """
    if (ra_degrees < 0) | (ra_degrees > 360):
        return None
    seconds_decimal_places = int(max(0, seconds_decimal_places))  # ensure int and non-negative.
    total_ra_seconds = ra_degrees * (3600 / 15)
    # TODO: remove code duplication (next ca. 20 lines).
    int_hours = int(total_ra_seconds // 3600)
    remaining_seconds = total_ra_seconds - 3600 * int_hours
    int_minutes = int(remaining_seconds // 60)
    remaining_seconds -= 60 * int_minutes
    if seconds_decimal_places > 0:
        seconds, fract_seconds = divmod(remaining_seconds, 1)
        int_fract_seconds = int(round(fract_seconds * 10 ** seconds_decimal_places))
    else:
        seconds, fract_seconds, int_fract_seconds = round(remaining_seconds), 0, 0
    int_seconds = int(seconds)
    if seconds_decimal_places > 0:
        if int_fract_seconds >= 10 ** seconds_decimal_places:
            int_fract_seconds -= 10 ** seconds_decimal_places
            int_seconds += 1
    if int_seconds >= 60:
        int_seconds -= 60
        int_minutes += 1
    if int_minutes >= 60:
        int_minutes -= 60
        int_hours += 1
    if int_hours >= 24:
        int_hours -= 24
    if seconds_decimal_places > 0:
        format_string = '{0:02d}:{1:02d}:{2:02d}.{3:0' + str(int(seconds_decimal_places)) + 'd}'
    else:
        format_string = '{0:02d}:{1:02d}:{2:02d}'
    ra_string = format_string.format(int_hours, int_minutes, int_seconds, int_fract_seconds)
    return ra_string


def dec_as_hex(dec_degrees, arcseconds_decimal_places=0):
    """ Input: float of Declination in degrees. TESTS OK 2020-10-24.
        Returns: Declination in hex, to desired precision. [string]
    """
    if (dec_degrees < -90) | (dec_degrees > +90):
        return None
    dec_string = degrees_as_hex(dec_degrees, arcseconds_decimal_places)
    return dec_string


def degrees_as_hex(angle_degrees, arcseconds_decimal_places=2):
    """ Takes degrees, returns hex representation. TESTS OK 2020-10-24.
    :param angle_degrees: any angle as degrees. [float]
    :param arcseconds_decimal_places: dec. places at end of hex string (no period if zero). [int]
    :return: same angle in hex notation, with proper sign, unbounded. [string]
    """
    if angle_degrees < 0:
        sign = "-"
    else:
        sign = "+"
    abs_degrees = abs(angle_degrees)
    arcseconds_decimal_places = int(max(0, arcseconds_decimal_places))  # ensure int and non-negative.
    total_arcseconds = abs_degrees * 3600
    # TODO: remove code duplication (next ca. 20 lines).
    int_degrees = int(total_arcseconds // 3600)
    remaining_arcseconds = total_arcseconds - 3600 * int_degrees
    int_arcminutes = int(remaining_arcseconds // 60)
    remaining_arcseconds -= 60 * int_arcminutes
    if arcseconds_decimal_places > 0:
        arcseconds, fract_arcseconds = divmod(remaining_arcseconds, 1)
        int_fract_arcseconds = int(round(fract_arcseconds * 10 ** arcseconds_decimal_places))
    else:
        arcseconds, fract_arcseconds, int_fract_arcseconds = round(remaining_arcseconds), 0, 0
    int_arcseconds = int(arcseconds)
    if arcseconds_decimal_places > 0:
        if int_fract_arcseconds >= 10 ** arcseconds_decimal_places:
            int_fract_arcseconds -= 10 ** arcseconds_decimal_places
            int_arcseconds += 1
    if int_arcseconds >= 60:
        int_arcseconds -= 60
        int_arcminutes += 1
    if int_arcminutes >= 60:
        int_arcminutes -= 60
        int_degrees += 1
    if int_degrees >= 360:
        int_degrees -= 360
    if arcseconds_decimal_places > 0:
        format_string = '{0}{1:02d}:{2:02d}:{3:02d}.{4:0' + str(int(arcseconds_decimal_places)) + 'd}'
    else:
        format_string = '{0}{1:02d}:{2:02d}:{3:02d}'
    hex_string = format_string.format(sign, int(int_degrees), int(int_arcminutes), int_arcseconds,
                                      int_fract_arcseconds)
    return hex_string


def parse_hex(hex_string):
    """ Helper function for RA and Dec parsing, takes hex string, returns list of strings representing
        floats.
        Not normally called directly by user. TESTS OK 2020-10-24.
    :param hex_string: string in either full hex ("12:34:56.7777" or "12 34 56.7777"),
               or degrees ("234.55")
    :return: list of strings representing floats (hours:min:sec or deg:arcmin:arcsec).
    """
    colon_list = hex_string.split(':')
    space_list = hex_string.split()  # multiple spaces act as one delimiter
    if len(colon_list) >= len(space_list):
        return [x.strip() for x in colon_list]
    return space_list


def concatenate_skycoords(skycoord_list):
    """Take SkyCoord or list of SkyCOords and return one SkyCoord object with all coordinates
       concatenated, in original order. Handles scalar SkyCoords, but always returns a list SkyCoord.
       :param skycoord_list: list of SkyCoords, or a single scalar SkyCoord.
       :return master_skycoord: SkyCoord object with all coordinates, in order. [array-based SkyCoord obj]
       """
    # Combine all skycoords found into a single astropy SkyCoord object sc:
    if isinstance(skycoord_list, SkyCoord):
        if skycoord_list.shape == ():  # if a single scalar SkyCoord was passed in.
            return SkyCoord(ra=[skycoord_list.ra], dec=[skycoord_list.dec])
        else:
            return skycoord_list

    if isinstance(skycoord_list, list):
        if len(skycoord_list) <= 0:
            raise ValueError('Empty SkyCoord object is not allowed.')
        # Combine via RA and Dec since SkyCoord has no concatenate capability (odd):
        ra_list, dec_list = [], []
        for sc in skycoord_list:
            if sc.shape == ():  # if scalar is SkyCoord.
                ra_list.append(sc.ra)
                dec_list.append(sc.dec)
            else:
                ra_list.extend(sc.ra)
                dec_list.extend(sc.dec)
        return SkyCoord(ra=ra_list, dec=dec_list)

    raise TypeError('parameter \'skycoord_list\' must be a SkyCoord object or a list of them.')


def combine_ra_dec_bounds(skycoord_list, extension_percent=3):
    """Take SkyCoord or list of SkyCoords and return smallest RA, Dec bounding box that will cover
       them all. Extend by desired percentage if desired.
       Caution: ensure RA is consistent across 0 (or 360) degrees.
        :param skycoord_list: coordinates around which to find bounding box. Usually coordinates of
            corners of one or more images. [astropy SkyCoord object, or a list of them].
        :param extension_percent: to extend bounding box by x% beyond actual edges, enter x. [float]
        :return ra_min, ra_max, dec_min, dec_max, all in degrees. [4-tuple of floats]]
            RA values are within range 0 <= RA < 360.
    """
    sc = concatenate_skycoords(skycoord_list)

    ra_wrap_angle = (circmean(sc.ra) + 180.0 * u.deg) % (360.0 * u.deg)
    ra_list = list(sc.ra.wrap_at(ra_wrap_angle))
    ra_min = min(ra_list)
    ra_max = max(ra_list)
    ra_extension = (ra_max - ra_min) * (extension_percent / 100.0)
    ra_min = (ra_min - ra_extension).wrap_at(360.0 * u.deg)
    ra_max = (ra_max + ra_extension).wrap_at(360.0 * u.deg)

    from astropy.coordinates import Angle
    dec_list = list(sc.dec)
    dec_min = Angle(min(dec_list))  # Angle(), in case exceeds -90 or +90 deg.
    dec_max = Angle(max(dec_list))  # "
    dec_extension = (dec_max - dec_min) * (extension_percent / 100.0)
    dec_min = Angle(max(-90.0 * u.deg, dec_min - dec_extension))
    dec_max = Angle(min(+90.0 * u.deg, dec_max + dec_extension))
    return ra_min.degree, ra_max.degree, dec_min.degree, dec_max.degree


__________FILESYSTEM_FUNCTIONS___________________________________ = 0


def make_directory_if_not_exists(directory_path):
    """ As the name says. Will not touch existing directory with this path, no matter what.
    :param directory_path: path specification for desired new directory. [string]
    Usage: make_directory_if_not_exists('C:/Astro/ACP/AN20201111')
    :return: path_preexists, True iff path already existed and this fn did nothing, else False. [boolean]
    """
    path_preexists = (os.path.exists(directory_path) and os.path.isdir(directory_path))
    if not path_preexists:
        os.mkdir(directory_path)
    return path_preexists


def count_files_immediate(dir_fullpath):
    """Return number of files (but not subdirectories) in this immediate directory (not recursive)."""
    n_files = sum(1 for element in os.scandir(dir_fullpath) if element.is_file())  # fast.
    return n_files


__________GENERAL_UTILITY_FUNCTIONS__________________________________________ = 0


def pressure_from_elevation(elevation):
    """ Return standard atmospheric pressure in mbar (= hectoPascals) from elevation in meters.
        Based on NOAA formula.
    """
    pressure_mbar = 1013.25 * pow(1.0 - (elevation / 44307.694), 5.25530)
    return pressure_mbar
