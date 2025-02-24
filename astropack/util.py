""" Module astropack.util
    Utility/support classes and functions. Not specific to observing targets or style.
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os
from datetime import datetime, timedelta, timezone
from math import floor, pow, ceil
from typing import List

# External packages:
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from astropy.stats import circmean

# 'ok' in the following means: 'converted to astropy Time and tested'.
__all__ = ['Timespan',                      # ok
           'hhmm',                          # ok
           # 'jd_from_datetime_utc',        # DELETED (native support in astropy).
           # 'datetime_utc_from_jd',        # DELETED (native support in astropy).
           'ra_as_degrees',                 # ok
           'hex_as_degrees',                # ok
           'dec_as_degrees',                # ok
           'ra_as_hours',                   # ok
           'dec_as_hex',                    # ok
           'degrees_as_hex',                # ok
           'parse_hex',                     # ok
           'concatenate_skycoords',         # ok
           'combine_ra_dec_bounds',         # ok
           'make_directory_if_not_exists',  # ok
           'count_files_immediate',         # ok
           'pressure_from_elevation']       # ok


THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_____CLASSES________________________________________________ = 0


class Timespan:
    """ Represents a span of time and provides for operations on and between timespans.
    Inputs may be astropy Time or python datetime, but internal time representation
    and returned times are always astropy Time.

    |Timespan| instances are expected to be immutable. Attributes are computed at time
    of construction and will not be recomputed if ``start_time`` or ``end_time`` are
    changed by user. For that, use the various methods, or make a ``.copy()``.

        Parameters
        ----------
        start_time : |py.datetime| or |Time|, UTC timezone required.
            Start time of new |Timespan|.

        end_time :  |py.datetime| or |Time|, UTC timezone required.
            End time of new |Timespan|.

        Attributes
        ----------
        start : |Time|
            Start time for this timespan, from input parameter ``start_time``.

        end : |Time|
            End time for this timespan, from input parameter ``end_time``.

        duration : |TimeDelta|
            Time duration, from start time to end time.

        seconds : non-negative float
            Time duration, from start time to end time, in seconds.

        days : non-negative float
            Time duration, from start time to end time, in days.

        midpoint : |Time|
            Time that exactly divides this timespan into two equal durations.

        Raises
        ------
        TypeError
            Raised when ``start_time`` or ``end_time`` is not datetime or
            a scalar |Time|.

        ValueError
            Raised when either ``start_time`` or ``end_time`` has no
            timezone information, or has a timezone other than UTC.

        Examples
        --------

    """
    def __init__(self, start_time: datetime | Time, end_time: datetime | Time):
        if not isinstance(start_time, (datetime, Time)) or \
          not isinstance(end_time, (datetime, Time)):
            raise TypeError('Both inputs must be either datetime or Time objects.')
        if isinstance(start_time, datetime):
            if start_time.tzname() != 'UTC':
                raise ValueError('First input parameter (datetime object)'
                                 ' must be in UTC.')
        if isinstance(end_time, datetime):
            if end_time.tzname() != 'UTC':
                raise ValueError('Second input parameter (datetime object)'
                                 ' must be in UTC.')
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

    def __eq__(self, other: 'Timespan') -> bool:
        """Provides == operator."""
        return self.start == other.start and self.end == other.end

    def copy(self) -> 'Timespan':
        """Return new |Timespan| exactly equivalent to this timespan."""
        return Timespan(self.start, self.end)

    def delay_by(self, delay: float | timedelta | TimeDelta) -> 'Timespan':
        """ Return new |Timespan| with start and end delayed by param 'delay'.

        Parameters
        ----------
        delay : float, or |py.timedelta|, or |TimeDelta|
            A duration by which to delay the start and end of new Timespan, relative
            to the present instance. If a float, delay is in seconds.

        Returns
        -------
        new_timespan : |Timespan|
            A new Timespan made later by duration ``delay``.
        """
        if isinstance(delay, (int, float)):
            delay = TimeDelta(delay, format='sec')
        elif isinstance(delay, (TimeDelta, timedelta)):
            delay = TimeDelta(delay)
        else:
            raise TypeError('Parameter \'delay\' must be int, float, py datetime, '
                            'or astropy.TimeDelta.')
        return Timespan(self.start + delay, self.end + delay)

    def expand_by(self, expansion: float | timedelta | TimeDelta) -> 'Timespan':
        """ Return new |Timespan| expanded or contracted at both start and end.

        Parameters
        ----------
        expansion : float, or |py.timedelta|, or |TimeDelta|
            The extent to which to extend both the start and the end away from the
            midpoint, relative to this timespan.
            If a float, expansion is in seconds.
            If negative, new |Timespan| is contracted relative to this timespan,
            but it will not contract beyond zero duration, in which case both
            start and end of the new |Timespan| will equal this timespan's ``midpoint``.

        Returns
        -------
        new_timespan : |Timespan|
            A new |Timespan| expanded at each limit by duration ``expansion``.
        """
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

    def intersection(self, other: 'Timespan') -> 'Timespan':
        """Return new |Timespan| containing only times common to (overlapping with)
        this timespan and ``other``. Returns zero-duration |Timespan| if timespans
        do not intersect.

        Parameters
        ----------
        other : |Timespan|
            Timespan to intersect with this Timespan.

        Returns
        -------
        new_timespan : |Timespan|
            Intersection of timespans.
        """
        new_start = max(self.start, other.start)
        new_end = min(self.end, other.end)
        return Timespan(new_start, new_end)

    def union(self, other: 'Timespan') -> 'Timespan' | List['Timespan']:
        """Return a new |Timespan| representing all times in either timespan.
        If the present and other timespan do not intersect,
        return list of the two |Timespan|.

        Parameters
        ----------
        other : |Timespan|
            The timespan to combine with the present timespan.

        Returns
        -------
        new_timespan : |Timespan|, or list of 2 |Timespan|
            Union (combination) of the present timespan with ``other``.
        """
        if (self.end < other.start) or (self.start > other.end):
            return [self, other]
        return Timespan(min(self.start, other.start), max(self.end, other.end))

    def add(self, other: 'Timespan') -> 'Timespan':
        """Alias for :meth:`Timespan.union()`."""
        return self.union(other)

    def subtract(self, other: 'Timespan') -> 'Timespan' | List['Timespan']:
        """Return a new |Timespan| having all times in this timespan, but not
        in ``other``, in effect subtracting any times in ``other``.

        Parameters
        ----------
        other : |Timespan|
            Timespan to subtract from this timespan.

        Returns
        -------
        new_timespan : |Timespan| or list of 2 |Timespan|
            A new |Timespan| having all times in this timespan but not
            in ``other``.
            If the two timespans do not intersect,
            returns a copy of the present timespan.
            If this timespan is wholly contained within ``other``,
            returns zero-length |Timespan| instance.
            If ``other`` is wholly contained within the prsent timespan,
            returns a list of 2 |Timespan|.
        """
        if (self.intersection(other).seconds == 0) or (other.seconds == 0):
            return self.copy()
        if self.start < other.start and other.end < self.end:
            # Other timespan is wholly within this timespan, so return 2 timespans:
            return [Timespan(self.start, other.start), Timespan(other.end, self.end)]
        if other.start <= self.start and self.end <= other.end:
            # This timespan is wholly within other timespan,
            # so return zero-length timespan at midpoint:
            return Timespan(self.midpoint, self.midpoint)
        if self.start < other.start:
            return Timespan(self.start, other.start)
        return Timespan(other.end, self.end)

    def contains(self, other: 'Timespan') -> bool:
        """Returns True if this timespan wholly contains ``other``, else returns False.

        Parameters
        ----------
        other : |py.datetime|, or |Time|, or |Timespan|
            The time or timespan to be tested as wholly contained within the present
            timespan.

        Returns
        -------
        is_contained : bool
            True iff ``other`` is contained within this timespan, else False.
        """
        if isinstance(other, (Time, datetime)):
            return self.start <= Time(other) <= self.end
        elif isinstance(other, Timespan):
            return (self.start <= other.start) & (self.end >= other.end)
        raise TypeError('Timespan.contains() requires Time, datetime, '
                        'or Timespan object as parameter.')

    def split_at(self, time: datetime | Time) -> 'Timespan' | List['Timespan']:
        """Splits current timespan at ``time`` to give two contiguous |Timespan|,
        or the original timespan if ``time`` is not contained.

        Parameters
        ----------
        time : |py.datetime| or |Time|
            Time at which to split this timespan.

        Returns
        -------
        timespans : list of 2 |Timespan|, or |Timespan|
            Two contiguous timespans split at ``time``.
            If ``time`` is not contained within this timespan, return a copy
            this |Timespan| instance.
        """
        if not isinstance(time, (Time, datetime)):
            raise TypeError('Timespan.split_at() requires Time '
                            'or datetime object as parameter.')
        time = Time(time)
        if self.start < time < self.end:
            return [Timespan(self.start, time), Timespan(time, self.end)]
        return self.copy()

    @staticmethod
    def longer(ts1: 'Timespan', ts2: 'Timespan') -> 'Timespan':
        """Returns the |Timespan| with longer duration (larger ``.seconds``), or
        of equal duration, return ``ts1``.

        Parameters
        ----------
        ts1, ts2 : |Timespan|
            The two timespans, from which the longer is to be selected.

        Returns
        -------
        longer_timespan : |Timespan|
            Timespan having longer duration, or ``ts1`` if durations are equal.
        """
        if ts1.seconds > ts2.seconds:
            return ts1
        elif ts2.seconds > ts1.seconds:
            return ts2
        elif ts1.start <= ts2.start:
            return ts1
        return ts2

    def periodic_events(self, ref_time: datetime | Time, period: timedelta | TimeDelta,
                        max_events: int | None = 10) -> List[Time]:
        """Returns times of periodic events that occur within a given |Timespan|.

        Parameters
        ----------
        ref_time : |py.datetime| or |Time|
            Time of any occurence of the periodic event.
            Need not be contained within present timespan.

        period : |py.timedelta| or |TimeDelta|
            Length of time between periodic events.

        max_events : positive int, or None, optional
            Maximum number of event times to be returned.
            If None, number of returned event times is not limited.
            Default is 10.

        Returns
        -------
        events : list of |Time|
            List of periodic event times.
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
            return []
        n_events = min(max_events, 1 +
                       floor(Timespan(first_time, self.end).seconds / period.sec))
        event_times = [first_time + i * period for i in range(n_events)]
        return event_times

    def __str__(self) -> str:
        return f"Timespan '{str(self.start)}' to '{str(self.end)}' = " \
               f"{self.seconds:.3f} seconds."

    def __repr__(self):
        return f'Timespan({self.start.iso}, {self.end.iso})'


__________TIME_and_DATE_FUNCTIONS_____________________________________________ = 0


def hhmm(time: Time) -> str:
    """Return UTC time-of-day string of form 'hhmm' for a scalar astropy Time.

    If datetime falls on the half-minute, rounding is performed to the nearest even
    minute, per 'banker's rounding'.

    >>> from astropack.util import hhmm
    >>> from astropy.time import Time
    >>> dt = datetime(2016, 1, 31, 0, 0, 30, 0).replace(tzinfo=timezone.utc)
    >>> time = Time(dt)
    >>> print(hhmm(time))
    0000
    >>> dt = datetime(2016, 1, 31, 0, 1, 30, 0).replace(tzinfo=timezone.utc)
    >>> time = Time(dt)
    >>> print(hhmm(date))
    0002

    Parameters
    ----------
    time : |Time|, scalar
        Time to represent as 'hhmm'. A scalar astropy Time object.

    Returns
    -------
    hhmm : str
        String of form 'hhmm' (e.g., '1352') representing UTC time of day.
    """
    # Apply banker's rounding: for a tie, give the nearest even value:
    datetime_utc = time.to_datetime()
    minutes_of_day = int(round(datetime_utc.hour*60 +
                               datetime_utc.minute +
                               datetime_utc.second/60 +
                               datetime_utc.microsecond/(60*1000000))) % 1440
    hh, mm = divmod(minutes_of_day, 60)
    return '{0:0>4d}'.format(100 * hh + mm)


def nearest_time(time_list: List[Time], ref_time: Time) -> Time | None:
    """ Return the time from time_list that is nearest to ref_time. """
    if len(time_list) == 0:
        return None
    return min(time_list, key=lambda t: abs(t - ref_time))


_____RA_and_DEC_FUNCTIONS_____________________________________ = 0


def ra_as_degrees(ra_string: str) -> float | None:
    """Takes string representing Right Ascension as either sexigesimal hours or
    decimal degrees, returns value in degrees.

    Parameters
    ----------
    ra_string : str
         String representation of Right Ascension in hex (e.g.,
         "12:34:56.7777" or "12 34 56.7777") or in degrees (e.g., "234.55").

    Returns
    -------
    ra_degrees : float, or None
        Right Ascension in degrees in range [0, 360], or None if value is < 0 or > 360.
    """
    ra_list = parse_hex(ra_string)
    if len(ra_list) == 1:
        # input assumed to be in degrees:
        ra_degrees = float(ra_list[0])
    elif len(ra_list) == 2:
        # input assumed to be in sexigesimal:
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1])/60.0)
    else:
        ra_seconds = float(ra_list[2].split()[0])  # in case trailing comments
        ra_degrees = 15 * (float(ra_list[0]) + float(ra_list[1]) / 60.0 +
                           ra_seconds/3600.0)  # input assumed in hex.
    if (ra_degrees < 0) | (ra_degrees > 360):
        ra_degrees = None
    return ra_degrees


def hex_as_degrees(degrees_string: str) -> float:
    """Takes string representing an angle as sexigesimal or as decimal degrees,
    returns value in degrees.

    Parameters
    ----------
    degrees_string : string
        String, either in sexigesimal format Hex string representing an angle,
        (e.g., "-12:34:56.7777", or "-12 34 56.7777") or in decimal degrees
        (e.g., "-24.55").

    Returns
    -------
    angle_degrees : float
        Angle in degrees.
    """
    dec_list = parse_hex(degrees_string)
    # dec_list = [dec.strip() for dec in dec_list]
    if dec_list[0].startswith("-"):
        sign = -1
    else:
        sign = 1
    if len(dec_list) == 1:
        angle_degrees = float(dec_list[0])  # input assumed to be in degrees.
    elif len(dec_list) == 2:
        # input is hex.
        angle_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1])/60.0)
    else:
        deg_seconds = float(dec_list[2].split()[0])  # in case trailing comments
        angle_degrees = sign * (abs(float(dec_list[0])) + float(dec_list[1]) / 60.0 +
                                deg_seconds/3600.0)  # input is hex.
    return angle_degrees


def dec_as_degrees(dec_string: str) -> float | None:
    """Takes Declination as string (hex or degrees), returns degrees as float.

    Parameters
    ----------
    dec_string : str
        String representing declination in full sexigesimal format (e.g.,
        "-12:34:56.7777", or "-12 34 56.7777") or in decimal degrees (e.g., "-24.55").

    Returns
    -------
    dec_degrees : float, or None
        Declination in degrees, or None if ``dec_string`` is outside range [-90, 90].
    """
    dec_degrees = hex_as_degrees(dec_string)
    if (dec_degrees < -90) | (dec_degrees > +90):
        dec_degrees = None
    return dec_degrees


def ra_as_hours(ra_degrees: float, seconds_decimal_places: int = 2) -> str | None:
    """Takes Right Ascension in degrees, returns string representing RA in hours.

    Parameters
    ----------
    ra_degrees : float
        Right Ascension in degrees, limited to range [0, 360].

    seconds_decimal_places : non-negative int, optional
        Number of decimal places after the decimal point.
        If zero, decimal point is removed. Default is 2.

    Returns
    -------
    ra_string : str or None
        Right Ascension in hours, in sexigesimal format, as '23:59:57.60'.
        If ``ra_degrees`` is outside range [0, 360], return None.
    """
    if (ra_degrees < 0) | (ra_degrees > 360):
        return None
    seconds_decimal_places = int(max(0, seconds_decimal_places))
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
        format_string = '{0:02d}:{1:02d}:{2:02d}.{3:0' + \
                        str(int(seconds_decimal_places)) + 'd}'
    else:
        format_string = '{0:02d}:{1:02d}:{2:02d}'
    ra_string = format_string.format(int_hours, int_minutes, int_seconds,
                                     int_fract_seconds)
    return ra_string


def dec_as_hex(dec_degrees: float, arcseconds_decimal_places: int = 0) -> str | None:
    """Takes degrees Declination, returns sexigesimal string representation.

    Parameters
    ----------
    dec_degrees : float
        Declination in degrees, within range [-90, 90].

    arcseconds_decimal_places : non-negative int, optional
        In the arcseconds place, the number of decimal places after the decimal point.
        If zero, decimal point is removed. Default is zero.

    Returns
    -------
    dec_string : str, or None
        String representing declination in sexigesimal format, as '-13:59:57.60'.
        If ``dec_degrees`` is outside range [-90, 90], return None.
    """
    if (dec_degrees < -90) | (dec_degrees > +90):
        return None
    dec_string = degrees_as_hex(dec_degrees, arcseconds_decimal_places)
    return dec_string


def degrees_as_hex(angle_degrees: float, arcseconds_decimal_places: int = 2) -> str:
    """Takes angle in degrees, returns sexigesimal string representation.

    Parameters
    ----------
    angle_degrees : float
        Angle in degrees.

    arcseconds_decimal_places : non-negative int, optional
        Number of decimal places after the decimal point in the arcseconds field.
        If zero, the decimal point is removed. Default is 2.

    Returns
    -------
    hex_string : str
        Sexigesimal string representation of angle, e.g., '-23:44:51.22'.
    """
    sign = '-' if angle_degrees < 0 else '+'
    abs_degrees = abs(angle_degrees)
    arcseconds_decimal_places = int(max(0, arcseconds_decimal_places))
    total_arcseconds = abs_degrees * 3600
    # TODO: remove code duplication (next ca. 20 lines).
    int_degrees = int(total_arcseconds // 3600)
    remaining_arcseconds = total_arcseconds - 3600 * int_degrees
    int_arcminutes = int(remaining_arcseconds // 60)
    remaining_arcseconds -= 60 * int_arcminutes
    if arcseconds_decimal_places > 0:
        arcseconds, fract_arcseconds = divmod(remaining_arcseconds, 1)
        int_fract_arcseconds = \
            int(round(fract_arcseconds * 10 ** arcseconds_decimal_places))
    else:
        arcseconds, fract_arcseconds, int_fract_arcseconds = \
            round(remaining_arcseconds), 0, 0
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
        format_string = '{0}{1:02d}:{2:02d}:{3:02d}.{4:0' + \
                        str(int(arcseconds_decimal_places)) + 'd}'
    else:
        format_string = '{0}{1:02d}:{2:02d}:{3:02d}'
    hex_string = format_string.format(sign, int(int_degrees), int(int_arcminutes),
                                      int_arcseconds, int_fract_arcseconds)
    return hex_string


def parse_hex(hex_string: str) -> list[str, str, str]:
    """Takes sexigesimal string representing an angle,
    returns a list of 3 str representing int or float values.
    Helper function not normally called diretly by user code.

    Parameters
    ----------
    hex_string : str
        String representation of angle, sexagesimal (colon- or space-delimited) or
        decimal, e.g., ("12:34:56.7777" or "12 34 56.7777").

    Returns
    -------
    string_list : list of 3 str
        List of three strings parsed from ``hex_string``.
    """
    colon_list = hex_string.split(':')
    space_list = hex_string.split()  # multiple spaces act as one delimiter
    if len(colon_list) >= len(space_list):
        string_list = [x.strip() for x in colon_list]
    else:
        string_list = space_list
    if len(string_list) == 1:
        return string_list
    string_list = (string_list + ['0', '0'])[:3]  # pad with zeroes, clip to 3 elements.
    return string_list


def concatenate_skycoords(skycoord_input: SkyCoord | list[SkyCoord]) -> SkyCoord:
    """Compiles an astropy |SkyCoord| or list of |SkyCoord|, array-based or scalar
    in any combination, and returns one array-based |SkyCoord| object containing
    all input coordinates. Order of input skycoords is retained.

    Parameters
    ----------
    skycoord_input : |SkyCoord|, or list of |SkyCoord|
        The input astropy |SkyCoord| coordinates.
        If a list, each list entry must be a |SkyCoord| instance,
        but each such entry may be scalar- or array-based.

    Returns
    -------
    skycoord_result : |SkyCoord|, array-based
        One array-based |SkyCoord| instance with all input coordinates.
        Original order is retained.
    """
    # Combine all skycoords found into a single astropy SkyCoord object sc:
    if isinstance(skycoord_input, SkyCoord):
        if skycoord_input.shape == ():  # if a single scalar SkyCoord was passed in.
            return SkyCoord(ra=[skycoord_input.ra], dec=[skycoord_input.dec])
        else:
            return skycoord_input

    if isinstance(skycoord_input, list):
        if len(skycoord_input) <= 0:
            raise ValueError('Empty input list is not allowed.')
        # Combine via RA and Dec since SkyCoord has no concatenate capability (odd):
        ra_list, dec_list = [], []
        for sc in skycoord_input:
            if sc.shape == ():  # if SkyCoord instance is scalar.
                ra_list.append(sc.ra)
                dec_list.append(sc.dec)
            else:
                ra_list.extend(sc.ra)
                dec_list.extend(sc.dec)
        return SkyCoord(ra=ra_list, dec=dec_list)
    raise TypeError('Parameter \'skycoord_input\' must be a '
                    'SkyCoord object or a list of them.')


def combine_ra_dec_bounds(skycoord_input: SkyCoord | list[SkyCoord],
                          extension_percent: int = 3) \
        -> tuple[float, float, float, float] | None:
    """Computes and returns sky coordinates of the smallest RA, Dec bounding box
    that can cover all input sky coordinates.

    Parameters
    ----------
    skycoord_input : |SkyCoord| or list of |SkyCoord|
        Sky coordinates (at least one) around which to find bounding box.

        These are usually the coordinates of all four corners of one or more images.

    extension_percent : non-negative float, optional
        Percent of bounding box dimensions by which to extend the bounding box,
        at user's discretion. Default is 3.

    Returns
    -------
    bounding_box : tuple of 4 float, or None
        RA and Dec of bounding box edges, each in degrees, as
        (ra_min, ra_max, dec_min, dec_max).
        RA zero-crossing is handled gracefully.
        If no sky coordinates in ``skycoord_list``, return None.
    """
    sc = concatenate_skycoords(skycoord_input)
    if len(sc) == 0:
        return None

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


def make_directory_if_not_exists(directory_path: str) -> bool:
    """Makes new directory on the file system only if no directory already exists
    at ``directory_path``.

    Parameters
    ----------
    directory_path : str
        The full path to the new directory wanted.

    Returns
    -------
    path_preexists : bool
        True if a directory already existed at ``directory_path`` in which case no
        action is taken; else False and directory is created.
    """
    path_preexists = (os.path.exists(directory_path) and os.path.isdir(directory_path))
    if not path_preexists:
        os.mkdir(directory_path)
    return path_preexists


def count_files_immediate(dir_fullpath: str) -> int:
    """Rapidly counts number of files in the specified directory, but not in child
    directories (non-recursive count).

    Parameters
    ----------
    dir_fullpath : str
        Full path to directory.

    Returns
    -------
    n_files : int
        Count of files in directory ``dir_fullpath``.
    """
    n_files = sum(1 for element in os.scandir(dir_fullpath) if element.is_file())
    return n_files


__________GENERAL_UTILITY_FUNCTIONS__________________________________________ = 0


def pressure_from_elevation(elevation: float) -> float:
    """Return standard atmospheric pressure in mbar (= hectoPascals) from
    given elevation in meters. Derived from a formula from the U.S. NOAA at
    https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf

    Parameters
    ----------
    elevation : float
        Elevation in meters.

    Returns
    -------
    pressure_mbar : float
        Standard (nominal) atmospheric pressure in millibars.
    """
    pressure_mbar = 1013.25 * pow(1.0 - (elevation / 44307.694), 5.25530)
    return pressure_mbar
