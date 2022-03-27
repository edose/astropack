""" test_util.py """

__author__ = "Eric Dose, Albuquerque"

# Python core:
import os
from datetime import datetime, timezone, timedelta

# External packages:
import pytest
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord

# Test target:
from astropack import util


THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


__________TIME_and_DATE_FUNCTIONS____________________________________________ = 0


def test_hhmm_from_datetime_utc():
    dt = datetime(2016, 1, 1, 23, 34, 45, 454545).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '2335'
    dt = datetime(2016, 1, 1, 23, 34, 29, 999999).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '2334'
    dt = datetime(2016, 1, 1, 23, 59, 31, 454545).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    dt = datetime(2016, 1, 31, 0, 0, 0, 0).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    # Using banker's rounding to nearest minute:
    dt = datetime(2016, 1, 31, 0, 0, 30, 0).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0000'
    dt = datetime(2016, 1, 31, 0, 1, 30, 0).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0002'
    dt = datetime(2016, 1, 31, 0, 0, 30, 1).replace(tzinfo=timezone.utc)
    assert util.hhmm_from_datetime_utc(dt) == '0001'


__________RA_and_DEC_FUNCTIONS_____________________________________________________ = 0


def test_ra_as_degrees():
    # Case: pass in RA in degrees:
    assert util.ra_as_degrees("180") == 180.0
    assert util.ra_as_degrees("0") == 0.0
    assert util.ra_as_degrees("360") == 360.0
    assert util.ra_as_degrees("-0.1") is None
    assert util.ra_as_degrees("360.1") is None

    # Case: pass in hex RA hours:
    assert util.ra_as_degrees("12:00") == util.ra_as_degrees("12 00") == 180.0
    assert util.ra_as_degrees("12:00:00") == util.ra_as_degrees("12 00 00") == 180.0
    assert util.ra_as_degrees("0:00:00") == util.ra_as_degrees("0 00 00") == 0.0
    assert util.ra_as_degrees("11:16:30") == util.ra_as_degrees("11 16 30") == 169.125
    assert util.ra_as_degrees("11:16:30.5") == pytest.approx(169.127083, abs=0.000002)
    assert util.ra_as_degrees("24:00:01") is None
    assert util.ra_as_degrees("24 00 01") is None

    # Cases: unusual syntax:
    assert util.ra_as_degrees("11:16:30 # comment") == util.ra_as_degrees("11:16:30")
    assert util.ra_as_degrees("11:16:30.5 # comment") == \
           util.ra_as_degrees("11:16:30.5")
    assert util.ra_as_degrees("11:16:30:444:whatever") == util.ra_as_degrees("11:16:30")


def test_hex_as_degrees():
    # Case: pass in string(float) degrees (and retain as whole-number degrees):
    assert util.hex_as_degrees("12") == 12.0
    assert util.hex_as_degrees("-12") == -12.0
    assert util.hex_as_degrees("0") == 0.0
    assert util.hex_as_degrees("-0") == 0.0
    assert util.hex_as_degrees("90") == 90
    assert util.hex_as_degrees("+90") == 90
    assert util.hex_as_degrees("-90") == -90
    assert util.hex_as_degrees("90.125") == 90.125
    assert util.hex_as_degrees("-90.125") == -90.125

    # Pass in hex degrees:
    assert util.hex_as_degrees("88:45") == util.hex_as_degrees("88 45") == 88.75
    assert util.hex_as_degrees("-88:45") == util.hex_as_degrees("-88 45") == -88.75
    assert util.hex_as_degrees("12:34:30") == util.hex_as_degrees("12 34 30") == 12.575
    assert util.hex_as_degrees("-12:34:30") == \
           util.hex_as_degrees("-12 34 30") == -12.575
    assert util.hex_as_degrees("91:34:30") == util.hex_as_degrees("91 34 30") == 91.575
    assert util.hex_as_degrees("-91:34:30") == \
           util.hex_as_degrees("-91 34 30") == -91.575
    assert util.hex_as_degrees("91:45") == util.hex_as_degrees("91 45") == 91.75
    assert util.hex_as_degrees("-91:45") == util.hex_as_degrees("-91 45") == -91.75

    # Cases: unusual syntax:
    assert util.hex_as_degrees("12:34:30 # comment") == \
           util.hex_as_degrees("12:34:30")
    assert util.hex_as_degrees("12:34:30.5 # comment") == \
           util.hex_as_degrees("12:34:30.5")
    assert util.hex_as_degrees("12:34:30:444:whatever") == \
           util.hex_as_degrees("12:34:30")


def test_dec_as_degrees():
    # Case: pass in string(float):
    assert util.dec_as_degrees("+12") == 12.0
    assert util.dec_as_degrees("-12") == -12.0
    assert util.dec_as_degrees("0") == 0.0
    assert util.dec_as_degrees("-0") == 0.0
    assert util.dec_as_degrees("90") == 90
    assert util.dec_as_degrees("-90") == -90
    assert util.dec_as_degrees("90.1") is None
    assert util.dec_as_degrees("-90.1") is None

    # Case: pass in hex string:
    assert util.dec_as_degrees("88:45") == util.dec_as_degrees("88 45") == 88.75
    assert util.dec_as_degrees("-88:45") == util.dec_as_degrees("-88 45") == -88.75
    assert util.dec_as_degrees("12:34:30") == util.dec_as_degrees("12 34 30") == 12.575
    assert util.dec_as_degrees("-12:34:30") == \
           util.dec_as_degrees("-12 34 30") == -12.575
    assert util.dec_as_degrees("91:34:30") is None
    assert util.dec_as_degrees("91 34 30") is None
    assert util.dec_as_degrees("-91:34:30") is None
    assert util.dec_as_degrees("-91 34 30") is None
    assert util.dec_as_degrees("91:45") is None
    assert util.dec_as_degrees("91 45") is None
    assert util.dec_as_degrees("-91:45") is None
    assert util.dec_as_degrees("-91 45") is None

    # Cases: unusual syntax:
    assert util.dec_as_degrees("12:34:30 # comment") == \
           util.dec_as_degrees("12:34:30")
    assert util.dec_as_degrees("12:34:30.5 # comment") == \
           util.dec_as_degrees("12:34:30.5")
    assert util.dec_as_degrees("12:34:30:444:whatever") == \
           util.dec_as_degrees("12:34:30")
    assert util.dec_as_degrees("112:34:30:444:whatever") is None
    assert util.dec_as_degrees("-92:34:30:444:whatever") is None


def test_ra_as_hours():
    # Test with default number of decimal places = 2:
    assert util.ra_as_hours(0.0) == "00:00:00.00"
    assert util.ra_as_hours(20.0) == "01:20:00.00"
    assert util.ra_as_hours(169.125) == "11:16:30.00"
    assert util.ra_as_hours(180.0) == "12:00:00.00"
    assert util.ra_as_hours(360.0) == "00:00:00.00"
    assert util.ra_as_hours(359.99) == "23:59:57.60"
    assert util.ra_as_hours(359.999) == "23:59:59.76"
    assert util.ra_as_hours(359.9999) == "23:59:59.98"
    assert util.ra_as_hours(359.99999) == "00:00:00.00"
    assert util.ra_as_hours(359.999999) == "00:00:00.00"
    assert util.ra_as_hours(359.9999999) == "00:00:00.00"
    assert util.ra_as_hours(359.99999999) == "00:00:00.00"
    assert util.ra_as_hours(-0.01) is None
    assert util.ra_as_hours(360.01) is None
    assert util.ra_as_hours(-44) is None
    assert util.ra_as_hours(654) is None
    # Test decimal places other than 2:
    assert util.ra_as_hours(359.99, 0) == '23:59:58'
    assert util.ra_as_hours(359.999, 0) == '00:00:00'
    assert util.ra_as_hours(359.99, 4) == '23:59:57.6000'
    assert util.ra_as_hours(359.9999, 4) == '23:59:59.9760'


def test_dec_as_hex():
    # Test with default number of decimal places = 0:
    assert util.dec_as_hex(0.0) == "+00:00:00"
    assert util.dec_as_hex(+90.0) == "+90:00:00"
    assert util.dec_as_hex(-90.0) == "-90:00:00"
    assert util.dec_as_hex(0.001) == "+00:00:04"
    assert util.dec_as_hex(-69.125) == "-69:07:30"
    assert util.dec_as_hex(69.125) == "+69:07:30"
    assert util.dec_as_hex(90.001) is None
    assert util.dec_as_hex(-90.001) is None
    assert util.dec_as_hex(255) is None
    assert util.dec_as_hex(-255) is None
    # Test decimal places other than 0:
    assert util.dec_as_hex(0.0, 2) == "+00:00:00.00"
    assert util.dec_as_hex(+90.0, 2) == "+90:00:00.00"
    assert util.dec_as_hex(-90.0, 1) == "-90:00:00.0"
    assert util.dec_as_hex(0.001, 4) == "+00:00:03.6000"
    assert util.dec_as_hex(-69.125, 1) == "-69:07:30.0"
    assert util.dec_as_hex(69.125, 3) == "+69:07:30.000"
    assert util.dec_as_hex(89.999, 4) == '+89:59:56.4000'


def test_degrees_as_hex():
    assert util.degrees_as_hex(0.0) == "+00:00:00.00"
    assert util.degrees_as_hex(0.0) == util.degrees_as_hex(0.0, 2)  # default
    assert util.degrees_as_hex(+90.0) == "+90:00:00.00"
    assert util.degrees_as_hex(-90.0) == "-90:00:00.00"
    assert util.degrees_as_hex(0.001) == "+00:00:03.60"
    assert util.degrees_as_hex(-0.001) == "-00:00:03.60"
    assert util.degrees_as_hex(-69.125) == "-69:07:30.00"
    assert util.degrees_as_hex(69.125) == "+69:07:30.00"
    assert util.degrees_as_hex(90.001) == "+90:00:03.60"
    assert util.degrees_as_hex(-90.001, 4) == "-90:00:03.6000"
    assert util.degrees_as_hex(255) == "+255:00:00.00"
    assert util.degrees_as_hex(-255, 6) == "-255:00:00.000000"
    assert util.degrees_as_hex(359.9999, 4) == '+359:59:59.6400'
    assert util.degrees_as_hex(359.9999, 0) == "+00:00:00"


def test_parse_hex():
    assert util.parse_hex('00:00:00') == util.parse_hex('00 00 00') == \
           ['00', '00', '00']
    assert util.parse_hex('12:34:56.89') == \
           util.parse_hex('12 34 56.89') == ['12', '34', '56.89']
    assert util.parse_hex('-50:30') == \
           util.parse_hex('-50 30') == ['-50', '30', '0']
    assert util.parse_hex('34') == ['34']


def test_concatenate_skycoords():
    sc_scalar = SkyCoord(10.625, 41.25, frame='icrs', unit='deg')
    sc_1 = SkyCoord([11.625], [42.25], frame='icrs', unit='deg')
    sc_3 = SkyCoord([12.625, 11, 12], [43.25, 42, 43], frame='icrs', unit='deg')

    # Case: scalar SkyCoord objects:
    result = util.concatenate_skycoords(sc_scalar)
    assert result.shape == (1, )
    assert result.ra.degree == [sc_scalar.ra.degree]
    assert result.dec.degree == [sc_scalar.dec.degree]
    result = util.concatenate_skycoords([sc_scalar, sc_scalar])
    assert result.shape == (2, )
    assert all(result.ra.degree == [sc_scalar.ra.degree, sc_scalar.ra.degree])
    assert all(result.dec.degree == [sc_scalar.dec.degree, sc_scalar.dec.degree])

    # Case: list-based SkyCoord objects:
    result = util.concatenate_skycoords([sc_1, sc_3, sc_1])
    assert result.shape == (5, )
    assert list(result.ra.degree) == list(sc_1.ra.degree) + \
           list(sc_3.ra.degree) + list(sc_1.ra.degree)
    assert list(result.dec.degree) == list(sc_1.dec.degree) + \
           list(sc_3.dec.degree) + list(sc_1.dec.degree)

    # Case: mixed scalar and list-based SkyCoord objects:
    result = util.concatenate_skycoords([sc_1, sc_scalar, sc_3])
    assert result.shape == (5, )
    assert list(result.ra.degree) == list(sc_1.ra.degree) + \
           [sc_scalar.ra.degree] + list(sc_3.ra.degree)
    assert list(result.dec.degree) == list(sc_1.dec.degree) + \
           [sc_scalar.dec.degree] + list(sc_3.dec.degree)


def test_combine_ra_dec_bounds():
    sc_scalar = SkyCoord(10.625, 41.25, frame='icrs', unit='deg')
    sc_1 = SkyCoord([11.625], [42.25], frame='icrs', unit='deg')
    sc_3 = SkyCoord([12.625, 11, 12], [43.25, 42, 43], frame='icrs', unit='deg')

    # Cases: RA values not near 0 (or 360) degrees:
    bounds_a = util.combine_ra_dec_bounds(sc_scalar, extension_percent=0)
    assert bounds_a == (10.625, 10.625, 41.25, 41.25)
    bounds_b = util.combine_ra_dec_bounds(sc_3, extension_percent=0)
    assert bounds_b == (11.0, 12.625, 42, 43.25)
    bounds_c = util.combine_ra_dec_bounds([sc_scalar, sc_3, sc_1], extension_percent=3)
    assert bounds_c == pytest.approx((10.565, 12.685, 41.19, 43.31))

    # Cases: RA values that cross 0 (or 360) degrees:
    sc_cross = SkyCoord([359.75, 0.0, 359.875, 0.025],
                        [41.25, 42, 43, 42.5], frame='icrs', unit='deg')
    bounds = util.combine_ra_dec_bounds(sc_cross, extension_percent=0)
    assert bounds == (359.75, 0.025, 41.25, 43.0)
    bounds = util.combine_ra_dec_bounds(sc_cross, extension_percent=4)
    assert bounds == pytest.approx((359.739, 0.036, 41.18, 43.07))

    # Cases: Dec values near -90 or +90:
    sc_dec = SkyCoord([234, 235, 236], [-89, -89.5, -89.999], frame='icrs', unit='deg')
    bounds = util.combine_ra_dec_bounds(sc_dec, extension_percent=0)
    assert bounds == pytest.approx((234.0, 236.0, -89.999, -89.0))
    bounds = util.combine_ra_dec_bounds(sc_dec, extension_percent=2)
    assert bounds == pytest.approx((233.96, 236.04, -90.0, -88.98002))
    sc_dec = SkyCoord([234, 235, 236], [89, 89.5, 89.999], frame='icrs', unit='deg')
    bounds = util.combine_ra_dec_bounds(sc_dec, extension_percent=2)
    assert bounds == pytest.approx((233.96, 236.04, 88.98002, 90.0))

    # Case: empty input:
    sc_empty = SkyCoord([], [], frame='icrs', unit='deg')
    assert util.combine_ra_dec_bounds(sc_empty) is None


__________OTHER_FUNCTIONS__________________________________________________________ = 0


def test_make_directory_if_not_exists():
    fullpath = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, '$$dummy_directory$')
    # First, ensure directory does not exist at fullpath:
    try:
        os.rmdir(fullpath)
    except OSError:
        pass
    assert (os.path.exists(fullpath) and os.path.isdir(fullpath)) == False

    # Case: directory does not pre-exist:
    assert util.make_directory_if_not_exists(fullpath) == False
    assert (os.path.exists(fullpath) and os.path.isdir(fullpath)) == True

    # Case: directory does pre-exist:
    assert util.make_directory_if_not_exists(fullpath) == True
    assert (os.path.exists(fullpath) and os.path.isdir(fullpath)) == True

    # Cleanup: remove test directory:
    try:
        os.rmdir(fullpath)
    except OSError:
        pass
    assert (os.path.exists(fullpath) and os.path.isdir(fullpath)) == False


def test_count_files_immediate():
    # Case: many files:
    n_files = util.count_files_immediate('D:/Astro/Catalogs/ATLAS-refcat2/mag-0-16/')
    assert n_files == 64800

    # Case: no files:
    fullpath = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, '$$dummy_directory$')
    # First, ensure directory does not exist at fullpath:
    try:
        os.rmdir(fullpath)
    except OSError:
        pass
    util.make_directory_if_not_exists(fullpath)
    assert util.count_files_immediate(fullpath) == 0
    try:
        os.rmdir(fullpath)
    except OSError:
        pass


def test_pressure_from_elevation():
    assert util.pressure_from_elevation(0) == 1013.25
    assert util.pressure_from_elevation(200) == pytest.approx(989.443561)
    assert util.pressure_from_elevation(2000) == pytest.approx(794.87050397)


__________CLASSES_____________________________________________________________ = 0


def test_class_timespan():
    # Set up tests:
    dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
    dt2 = dt1 + timedelta(hours=1.5)
    t1, t2 = Time(dt1), Time(dt2)
    ts1 = util.Timespan(t1, t2)
    two_min = TimeDelta(120, format='sec')
    test_tolerance = TimeDelta(0.000001, format='sec')

    # Test constructor and attributes, normal case:
    assert isinstance(ts1.start, Time)
    assert isinstance(ts1.end, Time)
    assert isinstance(ts1.duration, TimeDelta)
    assert isinstance(ts1.seconds, float)
    assert isinstance(ts1.start, Time)
    assert ts1.start == t1
    assert ts1.end == t2
    assert ts1.duration == t2 - t1
    assert ts1.seconds == 1.5 * 3600
    assert ts1.days == 1 / 16
    assert ts1.midpoint == t1 + (t2-t1) / 2

    # Test exception on bad input:
    with pytest.raises(TypeError):
        _ = util.Timespan('this is a bad input', ts1.end)
    with pytest.raises(TypeError):
        _ = util.Timespan(ts1.start, 'this is a bad input')
    with pytest.raises(ValueError):
        _ = util.Timespan(datetime(2022, 1, 30, 13),
                          datetime(2022, 1, 30, 13, tzinfo=timezone.utc))
    with pytest.raises(ValueError):
        _ = util.Timespan(datetime(2022, 1, 30, 13, tzinfo=timezone.utc),
                          datetime(2022, 1, 30, 13))
    with pytest.raises(ValueError):
        _ = util.Timespan(Time([dt1, dt2]), Time([dt1, dt2]))

    # Test equality __eq__():
    assert util.Timespan(t1, t2) == ts1
    assert util.Timespan(t1, t1+timedelta(hours=1)) != ts1
    assert util.Timespan(t1+timedelta(hours=1), t1) != ts1
    assert util.Timespan(t1+timedelta(hours=1), t2+timedelta(hours=1)) != ts1
    assert util.Timespan(t1, t2) == util.Timespan(t1, t2)

    # Test equivalence and non-identity of Time and datetime inputs:
    assert util.Timespan(dt1, dt2) == util.Timespan(t1, t2)
    assert util.Timespan(dt1, dt2) is not util.Timespan(t1, t2)

    # Test copy():
    ts_copy = ts1.copy()
    assert ts_copy is not ts1
    assert ts_copy == ts1
    del ts_copy

    # Test .delay_by(), positive delay:
    target_positive = util.Timespan(ts1.start + two_min, ts1.end + two_min)
    ts1_delay = ts1.delay_by(120)  # input in positive seconds.
    assert ts1_delay == target_positive
    ts1_delay = ts1.delay_by(timedelta(seconds=120))  # input in positive timedelta.
    assert ts1_delay == target_positive
    ts1_delay = ts1.delay_by(TimeDelta(120, format='sec'))  # input in pos. TimeDelta.
    assert ts1_delay == target_positive
    # Test .delay_by(), negative delay:
    target_negative = util.Timespan(ts1.start - two_min, ts1.end - two_min)
    ts1_delay = ts1.delay_by(-120)  # input in negative seconds.
    assert ts1_delay == target_negative
    ts1_delay = ts1.delay_by(timedelta(seconds=-120))  # input in negative timedelta.
    assert abs(ts1_delay.start - target_negative.start) < test_tolerance
    assert abs(ts1_delay.end - target_negative.end) < test_tolerance
    ts1_delay = ts1.delay_by(TimeDelta(-120, format='sec'))  # input in neg. TimeDelta.
    assert abs(ts1_delay.start - target_negative.start) < test_tolerance
    assert abs(ts1_delay.end - target_negative.end) < test_tolerance
    # Test .delay_by(), zero delay:
    ts1_delay = ts1.delay_by(0)  # input in zero seconds.
    assert ts1_delay == ts1
    assert ts1_delay is not ts1
    del ts1_delay
    # Test type exception:
    with pytest.raises(TypeError):
        _ = ts1.delay_by('10 seconds')

    # Test .expand_by():
    target_positive = util.Timespan(ts1.start - two_min, ts1.end + two_min)
    ts1_expand = ts1.expand_by(two_min)  # input in positive seconds.
    assert abs(ts1_expand.start - target_positive.start) < test_tolerance
    assert abs(ts1_expand.end - target_positive.end) < test_tolerance
    # Test .expand_by() with small negative expansion (normal case):
    target_negative = util.Timespan(ts1.start + two_min, ts1.end - two_min)
    ts1_contract = ts1.expand_by(-two_min)
    assert abs(ts1_contract.start - target_negative.start) < test_tolerance
    assert abs(ts1_contract.end - target_negative.end) < test_tolerance
    # Test .expand_by() with large negative expansion (returns zero-length Timespan):
    target_negative = util.Timespan(ts1.midpoint, ts1.midpoint)
    ts1_contract = ts1.expand_by(-two_min * 1000)
    assert ts1_contract == target_negative
    assert ts1_contract is not ts1
    # Test .expand_by() with zero expansion (returns copy).
    ts1_contract = ts1.expand_by(0)
    assert ts1_contract == ts1
    assert ts1_contract is not ts1
    del ts1_expand, ts1_contract
    # Test type exception:
    with pytest.raises(TypeError):
        _ = ts1.expand_by('10 seconds')

    # Test .intersection():
    ts1_wider = ts1.expand_by(100)
    # Test .intersection() for one Timespan wholly contains the other:
    assert ts1_wider.intersection(ts1) == ts1
    assert ts1.intersection(ts1_wider) == ts1
    # Test .intersection() for partial overlap:
    ts1_earlier = ts1.delay_by(-120)
    # ts1_intersection = ts1_earlier.intersection(ts1)
    ts_target = util.Timespan(ts1.start, ts1_earlier.end)
    assert ts1.intersection(ts1_earlier) == ts_target
    assert ts1_earlier.intersection(ts1) == ts_target
    # Test .intersection() for no overlap at all:
    ts1_no_overlap = util.Timespan(ts1.end + timedelta(hours=1),
                                   ts1.end + timedelta(hours=2))
    assert ts1.intersection(ts1_no_overlap).seconds == 0
    assert ts1_no_overlap.intersection(ts1).seconds == 0
    del ts1_wider, ts1_earlier, ts1_no_overlap

    # Test .union() and .add() (synonyms):
    # Normal case:
    ts_other = util.Timespan(ts1.midpoint, ts1.end + timedelta(hours=1))
    assert ts1.union(ts_other) == util.Timespan(ts1.start, ts_other.end)
    assert ts_other.union(ts1) == ts1.union(ts_other)
    assert ts1.union(ts_other) == ts1.add(ts_other)
    # Test .union() for non-overlapping Timespans (return list of the two Timespans):
    ts_other = util.Timespan(ts1.end + timedelta(hours=1),
                             ts1.end + timedelta(hours=2))
    assert ts1.union(ts_other) == [ts1, ts_other]
    assert ts_other.union(ts1) == [ts_other, ts1]
    assert ts1.add(ts_other) == ts1.union(ts_other)
    # Test .union() where Timespans just touch:
    ts_other = util.Timespan(ts1.end, ts1.end + timedelta(hours=1))
    assert ts1.union(ts_other) == util.Timespan(ts1.start, ts_other.end)
    assert ts_other.union(ts1) == util.Timespan(ts1.start, ts_other.end)
    # Test .union() where one Timespan wholly contains the other:
    ts_other = ts1.expand_by(100)
    assert ts1.union(ts_other) == ts_other
    assert ts_other.union(ts1) == ts_other
    del ts_other

    # Test .subtract() for partial overlap:
    ts_later = ts1.delay_by(100)
    assert ts1.subtract(ts_later) == util.Timespan(ts1.start, ts_later.start)
    assert ts_later.subtract(ts1) == util.Timespan(ts1.end, ts_later.end)
    del ts_later

    # Test .subtract() for a Timespan of zero duration:
    ts_zero = util.Timespan(ts1.midpoint, ts1.midpoint)
    assert ts1.subtract(ts_zero) == ts1
    assert ts_zero.subtract(ts1).seconds == 0
    # Test .subtract() for no overlap (return zero-duration Timespan):
    ts_other = util.Timespan(ts1.end + timedelta(hours=1),
                             ts1.end + timedelta(hours=2))
    assert ts1.subtract(ts_other) == ts1
    assert ts_other.subtract(ts1) == ts_other
    # Test .subtract where one Timespan wholly contains the other:
    ts_expanded = ts1.expand_by(two_min)
    assert ts1.subtract(ts_expanded).seconds == 0
    assert ts_expanded.subtract(ts1) == [util.Timespan(ts_expanded.start, ts1.start),
                                         util.Timespan(ts1.end, ts_expanded.end)]
    del ts_zero, ts_other, ts_expanded

    # Test .contains():
    # Test exception on bad parameter type:
    with pytest.raises(TypeError):
        _ = ts1.contains(999)
    # Test .contains() for other as Time object:
    assert ts1.contains(ts1.start)
    assert ts1.contains(ts1.end)
    assert ts1.contains(ts1.midpoint)
    assert not ts1.contains(ts1.start - two_min)
    assert not ts1.contains(ts1.end + two_min)

    # Test.contains() for other as Timespan object:
    assert ts1.contains(ts1)  # edge case.
    ts_later = ts1.delay_by(two_min)
    assert not ts1.contains(ts_later)
    assert not ts_later.contains(ts1)
    ts_wider = ts1.expand_by(two_min)
    assert not ts1.contains(ts_wider)
    assert ts_wider.contains(ts1)
    ts_no_overlap = util.Timespan(ts1.end + two_min, ts1.end + 2 * two_min)
    assert not ts1.contains(ts_no_overlap)
    assert not ts_no_overlap.contains(ts1)
    del ts_later, ts_wider, ts_no_overlap

    # Test .split_at():
    dt1 = datetime(2016, 9, 10, 0, 0, 0, tzinfo=timezone.utc)
    dt2 = dt1 + timedelta(hours=1.5)
    dt_split = dt1 + timedelta(hours=1.0)
    splits = ts1.split_at(dt_split)
    assert splits == [util.Timespan(dt1, dt_split), util.Timespan(dt_split, dt2)]
    splits = ts1.split_at(dt1 - timedelta(hours=0.66))
    assert splits == util.Timespan(dt1, dt2)
    splits = ts1.split_at(dt1 + timedelta(hours=+22))
    assert splits == util.Timespan(dt1, dt2)

    # Test str():
    s = str(ts1)
    assert s.startswith("Timespan ")
    assert s.endswith(" = 5400.000 seconds.")
    del s

    # Test repr():
    r = repr(ts1)
    assert r == 'Timespan(2016-09-10 00:00:00.000, 2016-09-10 01:30:00.000)'

    # Test Timespan.longer():
    # Edge cases:
    assert util.Timespan.longer(ts1, ts1) == ts1
    ts_copy = ts1.copy()
    assert util.Timespan.longer(ts1, ts_copy) == ts1
    assert util.Timespan.longer(ts1, ts_copy) is ts1
    ts_zero = util.Timespan(ts1.start, ts1.start)
    assert util.Timespan.longer(ts_zero, ts_zero) == ts_zero
    assert util.Timespan.longer(ts1, ts_zero) == ts1
    assert util.Timespan.longer(ts_zero, ts1) == ts1
    # Equal length:
    ts_later = ts1.delay_by(two_min)
    assert util.Timespan.longer(ts1, ts_later) == ts1  # ts1 because it is earlier.
    assert util.Timespan.longer(ts_later, ts1) == ts1  # "
    # Normal cases:
    ts_wider = ts1.expand_by(two_min)
    assert util.Timespan.longer(ts1, ts_wider) == ts_wider
    assert util.Timespan.longer(ts_wider, ts1) == ts_wider
    del ts_copy, ts_zero, ts_later, ts_wider

    # Test Timespan.periodic_events():
    # Edge cases (but not subject to rounding error problems of astropy.time.Time):
    ts = util.Timespan(t1, t1)
    assert ts.periodic_events(t1 + two_min, TimeDelta(1000, format='sec')) == []
    assert ts1.periodic_events(ts1.start, 2 * ts1.duration) == [ts1.start]
    assert ts1.periodic_events(ts1.end, 2 * ts1.duration) == [ts1.end]
    assert ts1.periodic_events(ts1.start + two_min, 5 * ts1.duration) == \
           [ts1.start + two_min]
    ref = ts1.start - two_min
    period = 5 * ts1.duration
    assert ts1.periodic_events(ref, period) == []
    assert ts1.periodic_events(ts1.midpoint, 5 * ts1.duration) == [ts1.midpoint]

    # Test normal case:
    target_times = [ts1.start + TimeDelta(600, format='sec'),
                    ts1.start + TimeDelta(4600, format='sec'),
                    ts1.start + TimeDelta(8600, format='sec')]
    ts = util.Timespan(ts1.start, ts1.start + TimeDelta(10000, format='sec'))
    return_value = ts.periodic_events(ts.start - TimeDelta(3400, format='sec'),
                                      TimeDelta(4000, format='sec'))
    assert len(return_value) == len(target_times)
    assert all([rv.isclose(tt) for (rv, tt) in zip(return_value, target_times)])
