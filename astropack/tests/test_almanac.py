""" test_almanac.py
    All test functions and Astronight class appear to pass all tests as of 2022-02-09.
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os

# External packages:
import pytest
import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from skyfield.api import load, Star
import skyfield

# Author's packages:
from astropack import ini
from astropack.util import Timespan

# Test target:
from astropack import almanac

THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FOR_TEST_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'tests',
                                       '$data_for_test')

ONE_SECOND = TimeDelta(1, format='sec')
DSW_LON = -(105 + 31/60 + 44.32/3600)
DSW_LAT = 32 + 54/60 + 11.36/3600
DSW_ELEV = 2228  # meters


_____HELPERS_for_these_TESTS_________________________________________________ = 0


def make_new_site():
    """Returns a typical Site object."""
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    return site


def make_new_skyfield_engine():
    """Returns a SkyfieldEngine object."""
    site = make_new_site()
    sfe = almanac.SkyfieldEngine.from_site(site, -10)
    return sfe


def make_new_astronight(an_date=20220208):
    """Returns a typical Astronight object."""
    site = make_new_site()
    return almanac.Astronight(site, an_date, -9.0)


_____TEST_SUPPORT_CLASSES____________________________________________________ = 0


def test_class_an_date():
    a = almanac.AN_date('20221118')
    assert a.an_int == 20221118
    assert type(a.an_int) == int
    assert a.an_str == '20221118'
    assert type(a.an_str) == str
    assert (a.year, a.month, a.day) == (2022, 11, 18)

    a = almanac.AN_date(20221119)
    assert a.an_int == 20221119
    assert type(a.an_int) == int
    assert a.an_str == '20221119'
    assert type(a.an_str) == str
    assert (a.year, a.month, a.day) == (2022, 11, 19)

    assert almanac.AN_date(20221216).day_of_week.upper() == 'FRIDAY'
    assert almanac.AN_date(20221225).day_of_week.upper() == 'SUNDAY'

    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('not a date string')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('14920101')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('23910101')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('20220001')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('20221301')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('20220100')
    with pytest.raises(almanac.Invalid_ANDate_Error):
        _ = almanac.AN_date('20220188')


_____TEST_MODULE_LEVEL_FUNCTIONS_____________________________________________ = 0


def test_calc_approx_midnight():
    a = almanac.AN_date('20221118')
    approx_midnight = almanac.calc_approx_midnight(an_date=a, longitude_hours=-7)
    assert approx_midnight == Time('2022-11-19 07:00:00')
    approx_midnight = almanac.calc_approx_midnight(an_date=a, longitude_hours=0)
    assert approx_midnight == Time('2022-11-19 00:00:00')
    approx_midnight = almanac.calc_approx_midnight(an_date=a, longitude_hours=+2.25)
    assert approx_midnight == Time('2022-11-18 21:45:00')
    approx_midnight = almanac.calc_approx_midnight(an_date=a, longitude_hours=+14)
    assert abs(approx_midnight - Time('2022-11-18 10:00:00')) < ONE_SECOND
    approx_midnight = almanac.calc_approx_midnight(an_date=a, longitude_hours=-14)
    assert abs(approx_midnight - Time('2022-11-19 14:00:00')) < ONE_SECOND


def test_calc_phase_angle_bisector():
    site = make_new_site()
    # For MPs 3229 & 3166:
    times = [Time('2022-11-22 03:00:00'), Time('2022-11-24 05:00:00')]
    mp_skycoords = [SkyCoord(['01 25 39.7 +26 18 20.4'],
                             frame='icrs', unit=(u.hourangle, u.deg)),
                    SkyCoord(['04 31 32.2 +20 58 46.1'],
                             frame='icrs', unit=(u.hourangle, u.deg))]
    deltas = [1.291, 1.384]
    pabs = almanac.calc_phase_angle_bisector(times, mp_skycoords, deltas, site)
    pab_longitudes = [sc.lon.degree for sc in pabs]
    pab_latitudes = [sc.lat.degree for sc in pabs]
    expected_pab_longitudes = [36.4, 67.8]  # from minorplanet.info query
    expected_pab_latitudes = [12.8, -0.7]   # "
    assert all([calc == pytest.approx(exp, abs=0.1)
                for (calc, exp) in zip(pab_longitudes, expected_pab_longitudes)])
    assert all([calc == pytest.approx(exp, abs=0.1)
                for (calc, exp) in zip(pab_latitudes, expected_pab_latitudes)])


_____TEST_CLASS_SKYFIELD_ENGINE_CONSTRUCTORS_________________________________ = 0


def test_main_constructor():
    sfe = almanac.SkyfieldEngine(site_latitude=33, site_longitude=-105,
                                 site_elevation=2020, sun_altitude_dark=-9)
    assert type(sfe) == almanac.SkyfieldEngine
    assert (sfe.site_latitude, sfe.site_longitude, sfe.site_elevation) == \
           (33, -105, 2020)
    assert sfe.sun_alt_dark == -9


def test_constructor_from_site():
    site = make_new_site()
    sfe = almanac.SkyfieldEngine.from_site(site=site, sun_altitude_dark=-10)
    assert type(sfe) == almanac.SkyfieldEngine
    assert sfe.site_latitude == site.latitude
    assert sfe.site_longitude == site.longitude
    assert sfe.site_elevation == site.elevation
    assert sfe.sun_alt_dark == -10
    assert type(sfe.sf_obs) == skyfield.toposlib.GeographicPosition
    assert sfe.sf_obs.model.name == 'WGS84'


_____TEST_CLASS_SKYFIELD_ENGINE_PRIVATE_METHODS______________________________ = 0


# noinspection DuplicatedCode
def test__target_azalt():
    sfe = make_new_skyfield_engine()

    # Case: target is skyfield Star object:
    rigel = Star(ra_hours=5.2422978748, dec_degrees=-8.20164055)
    t = sfe.timescale.utc(2022, 5, 28, 7, 20, 0)
    az, alt = sfe._target_azalt(rigel, t)
    expected_az = 340.2435
    expected_alt = -64.0229
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)

    # Case: target is sun:
    sun = sfe.sun_eph
    t = sfe.timescale.utc(2022, 5, 28, 7, 20, 0)
    az, alt = sfe._target_azalt(sun, t)
    expected_az = 5.8897
    expected_alt = -35.4046
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)

    t = sfe.timescale.utc(2022, 2, 28, 7, 20, 0)
    az, alt = sfe._target_azalt(sun, t)
    expected_az = 3.1495
    expected_alt = -65.0093
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)

    t = sfe.timescale.utc(2022, 6, 28, 0, 20, 0)
    az, alt = sfe._target_azalt(sun, t)
    expected_az = 284.5591
    expected_alt = 21.4810
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)

    # Case: target is moon:
    moon = sfe.moon_eph
    t = sfe.timescale.utc(2022, 5, 28, 7, 20, 0)
    az, alt = sfe._target_azalt(moon, t)
    expected_az = 35.4763
    expected_alt = -35.2373
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)

    t = sfe.timescale.utc(2022, 6, 28, 0, 20, 0)
    az, alt = sfe._target_azalt(moon, t)
    expected_az = 293.2324
    expected_alt = 11.4594
    assert az == pytest.approx(expected_az, abs=0.1)
    assert alt == pytest.approx(expected_alt, abs=0.1)


# noinspection DuplicatedCode
def test__horizon_crossings():
    sfe = make_new_skyfield_engine()

    # Case: target is skyfield Star object:
    rigel = Star(ra_hours=5.2422978748, dec_degrees=-8.20164055)  # Rigel
    t_start = sfe.timescale.utc(2022, 5, 28, 7, 20, 0)
    t_end = sfe.timescale.utc(2022, 5, 30, 7, 20, 0)
    crossings = sfe._horizon_crossings(rigel, t_start, t_end, 30)
    assert isinstance(crossings, list)
    assert len(crossings) == 4
    assert isinstance(crossings[0], tuple)
    assert len(crossings[0]) == 2
    expected_times = [sfe.timescale.utc(2022, 5, 28, 16, 48, 51),
                      sfe.timescale.utc(2022, 5, 28, 22, 55, 52),
                      sfe.timescale.utc(2022, 5, 29, 16, 44, 56),
                      sfe.timescale.utc(2022, 5, 29, 22, 51, 57)]
    expected_types = [almanac.SkyfieldEngine._HorizonCrossing.RISING,
                      almanac.SkyfieldEngine._HorizonCrossing.SETTING,
                      almanac.SkyfieldEngine._HorizonCrossing.RISING,
                      almanac.SkyfieldEngine._HorizonCrossing.SETTING]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 1
        assert c[1] == etype

    # Case: target is sun:
    sun = sfe.sun_eph
    t_start = sfe.timescale.utc(2022, 5, 28, 17, 20, 0)
    t_end = sfe.timescale.utc(2022, 5, 30, 17, 20, 0)
    crossings = sfe._horizon_crossings(sun, t_start, t_end, -10)
    expected_times = [sfe.timescale.utc(2022, 5, 29, 2, 53, 52),
                      sfe.timescale.utc(2022, 5, 29, 11, 4, 58),
                      sfe.timescale.utc(2022, 5, 30, 2, 54, 35),
                      sfe.timescale.utc(2022, 5, 30, 11, 4, 32)]
    expected_types = [almanac.SkyfieldEngine._HorizonCrossing.SETTING,
                      almanac.SkyfieldEngine._HorizonCrossing.RISING,
                      almanac.SkyfieldEngine._HorizonCrossing.SETTING,
                      almanac.SkyfieldEngine._HorizonCrossing.RISING]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 1
        assert c[1] == etype

    # Case: target is moon:
    moon = sfe.moon_eph
    t_start = sfe.timescale.utc(2022, 1, 11, 17, 20, 0)
    t_end = sfe.timescale.utc(2022, 1, 13, 17, 20, 0)
    crossings = sfe._horizon_crossings(moon, t_start, t_end, 1)
    expected_times = [sfe.timescale.utc(2022, 1, 11, 19, 57, 11),
                      sfe.timescale.utc(2022, 1, 12, 9, 27, 46),
                      sfe.timescale.utc(2022, 1, 12, 20, 28, 57),
                      sfe.timescale.utc(2022, 1, 13, 10, 24, 24)]
    expected_types = [almanac.SkyfieldEngine._HorizonCrossing.RISING,
                      almanac.SkyfieldEngine._HorizonCrossing.SETTING,
                      almanac.SkyfieldEngine._HorizonCrossing.RISING,
                      almanac.SkyfieldEngine._HorizonCrossing.SETTING]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 10  # moon not nearly as good as more distant objects.
        assert c[1] == etype


# noinspection DuplicatedCode
def test__meridian_crossings():
    sfe = make_new_skyfield_engine()

    # Case: target is skyfield Star object:
    rigel = Star(ra_hours=5.2422978748, dec_degrees=-8.20164055)  # Rigel
    t_start = sfe.timescale.utc(2022, 5, 11, 17, 20, 0)
    t_end = sfe.timescale.utc(2022, 5, 13, 17, 20, 0)
    crossings = sfe._meridian_crossings(rigel, t_start, t_end)
    assert isinstance(crossings, list)
    assert len(crossings) == 4
    assert isinstance(crossings[0], tuple)
    assert len(crossings[0]) == 2
    expected_times = [sfe.timescale.utc(2022, 5, 11, 20, 59, 12),
                      sfe.timescale.utc(2022, 5, 12, 8, 57, 14),
                      sfe.timescale.utc(2022, 5, 12, 20, 55, 17),
                      sfe.timescale.utc(2022, 5, 13, 8, 53, 19)]
    expected_types = [almanac.SkyfieldEngine._MeridianCrossing.TRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.TRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 1
        assert c[1] == etype

    # Case: target is sun:
    sun = sfe.sun_eph
    t_start = sfe.timescale.utc(2022, 5, 11, 1, 20, 0)
    t_end = sfe.timescale.utc(2022, 5, 13, 1, 20, 0)
    crossings = sfe._meridian_crossings(sun, t_start, t_end)
    expected_times = [sfe.timescale.utc(2022, 5, 11, 6, 58, 31),
                      sfe.timescale.utc(2022, 5, 11, 18, 58, 30),
                      sfe.timescale.utc(2022, 5, 12, 6, 58, 30),
                      sfe.timescale.utc(2022, 5, 12, 18, 58, 29)]
    expected_types = [almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.TRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.TRANSIT]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 1
        assert c[1] == etype

    # Case: target is moon:
    moon = sfe.moon_eph
    t_start = sfe.timescale.utc(2022, 5, 11, 1, 20, 0)
    t_end = sfe.timescale.utc(2022, 5, 13, 1, 20, 0)
    crossings = sfe._meridian_crossings(moon, t_start, t_end)
    expected_times = [sfe.timescale.utc(2022, 5, 11, 2, 58, 31),
                      sfe.timescale.utc(2022, 5, 11, 15, 20, 59),
                      sfe.timescale.utc(2022, 5, 12, 3, 43, 30),
                      sfe.timescale.utc(2022, 5, 12, 16, 6, 16)]
    expected_types = [almanac.SkyfieldEngine._MeridianCrossing.TRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.TRANSIT,
                      almanac.SkyfieldEngine._MeridianCrossing.ANTITRANSIT]
    for c, etime, etype in zip(crossings, expected_times, expected_types):
        diff_seconds = abs(c[0] - etime) * 24 * 3600
        assert diff_seconds < 10  # moon not nearly as good as more distant objects.
        assert c[1] == etype


def test__closest_time():
    fn = almanac.SkyfieldEngine._nearest_skyfield_time
    timescale = load.timescale()
    times = [timescale.utc(2022, 9, 30, hour, 0, 0) for hour in range(6)]
    assert fn(timescale.utc(2021, 1, 1, 1, 1, 1), times) == times[0]
    assert fn(timescale.utc(2023, 5, 5, 5, 5, 5), times) == times[5]
    assert fn(timescale.utc(2022, 9, 30, 0, 1, 1), times) == times[0]
    assert fn(timescale.utc(2022, 9, 30, 3, 1, 1), times) == times[3]
    assert fn(timescale.utc(2022, 9, 30, 4, 30, 19), times) == times[5]
    assert fn(timescale.utc(2022, 9, 30, 1, 30, 1), times) == times[2]


def test__skycoord_to_skyfield_star():
    sc = SkyCoord(10.0, 20.0, unit="deg")
    sf_star = almanac.SkyfieldEngine._skycoord_to_skyfield_star(sc)
    assert isinstance(sf_star, Star)
    assert sf_star.ra.hours * 15.0 == pytest.approx(10.0)
    assert sf_star.dec.degrees == pytest.approx(20.0)
    del sc, sf_star


def test__astropy_time_to_skyfield():
    ts = load.timescale()
    apy_time = Time('2019-09-30 23:05:06', scale='utc')
    sf_time = almanac.SkyfieldEngine._astropy_time_to_skyfield(apy_time)
    assert isinstance(sf_time, skyfield.timelib.Time)
    assert sf_time.utc == (2019, 9, 30, 23, 5, 6.0)
    assert almanac.SkyfieldEngine._astropy_time_to_skyfield(apy_time) == \
           almanac.SkyfieldEngine._astropy_time_to_skyfield(apy_time, sf_timescale=ts)


def test__skyfield_time_to_astropy():
    ts = load.timescale()                            # type skyfield.timelib.Timescale
    sf_time = ts.utc(2019, 9, 30, 23, 5, 6)          # type skyfield.timelib.Time
    apy_time = almanac.SkyfieldEngine._skyfield_time_to_astropy(sf_time)
    assert isinstance(apy_time, Time)  # type astropy.time.Time
    assert apy_time.scale == 'utc'
    assert apy_time.format == 'iso'
    expected_apy_time = Time('2019-09-30 23:05:06')
    assert (apy_time - expected_apy_time).sec == pytest.approx(0, abs=1e-9)



_____TEST_CLASS_SKYFIELD_ENGINE_PUBLIC_METHODS_______________________________ = 0


def test_sun_azalt():
    sfe = make_new_skyfield_engine()
    az, alt = sfe.sun_azalt(Time('2022-02-28 07:20:00'))
    assert az == pytest.approx(3.1398, abs=0.02)
    assert alt == pytest.approx(-65.0095, abs=0.02)

    az, alt = sfe.sun_azalt(Time('2022-06-28 00:20:00'))
    assert az == pytest.approx(284.5572, abs=0.02)
    assert alt == pytest.approx(21.4844, abs=0.02)


def test_moon_azalt():
    sfe = make_new_skyfield_engine()
    az, alt = sfe.moon_azalt(Time('2022-02-28 07:20:00'))
    assert az ==  pytest.approx(78.1383, abs=0.1)  # skyfield not so accurate for moon.
    assert alt == pytest.approx(-61.3053, abs=0.1)   # "


def test_skycoord_azalt():
    sfe = make_new_skyfield_engine()

    # Case: scalar SkyCoord object passed in:
    rigel = SkyCoord(15 * 5.2422978748, -8.20164055, unit='deg')  # Rigel
    t = Time('2022-05-28 07:20:00')
    az, alt = sfe.skycoord_azalt(rigel, t)
    expected_az = 340.2345
    expected_alt = -64.0217
    assert az == pytest.approx(expected_az, abs=0.01)
    assert alt == pytest.approx(expected_alt, abs=0.01)


# noinspection DuplicatedCode
def test_prev_sunset_utc():
    sfe = make_new_skyfield_engine()

    # Case: near site's solar midnight.
    t = Time('2022-05-28 07:01:00')
    sunset_utc = sfe.prev_sunset_utc(t)
    expected_utc = Time('2022-05-28 02:02:38')
    diff_seconds = abs(sunset_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: near site's solar noon.
    t = Time('2022-05-26 19:01:00')
    sunset_utc = sfe.prev_sunset_utc(t)
    expected_utc = Time('2022-05-26 02:01:21')
    diff_seconds = abs(sunset_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: no prev sunset within allowable timespan:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site_barrow = ini.Site(site_fullpath)
    sfe_barrow = almanac.SkyfieldEngine.from_site(site_barrow, -10)
    t_ref = Time('2022-06-22 11:00:00')
    sunset_utc = sfe_barrow.prev_sunset_utc(ref_time_utc=t_ref, horizon_degrees=0)
    assert sunset_utc is None


# noinspection DuplicatedCode
def test_next_sunrise_utc():
    sfe = make_new_skyfield_engine()

    # Case: near site's solar midnight.
    t = Time('2022-05-28 07:01:00')
    sunrise_utc = sfe.next_sunrise_utc(t)
    expected_utc = Time('2022-05-28 11:55:56')
    diff_seconds = abs(sunrise_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: near site's solar noon.
    t = Time('2022-05-26 19:01:00')
    sunrise_utc = sfe.next_sunrise_utc(t)
    expected_utc = Time('2022-05-27 11:56:19')
    diff_seconds = abs(sunrise_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: no next sunrise within allowable timespan:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site_barrow = ini.Site(site_fullpath)
    sfe_barrow = almanac.SkyfieldEngine.from_site(site_barrow, -10)
    t_ref = Time('2022-06-22 11:00:00')
    sunrise_utc = sfe_barrow.next_sunrise_utc(ref_time_utc=t_ref)
    assert sunrise_utc is None


# noinspection DuplicatedCode
def test_prev_dark_start_utc():
    sfe = make_new_skyfield_engine()

    # Case: User provides sun altitude defining dark time:
    t = Time('2022-05-28 07:01:00')  # near site's solar midnight
    dark_start_utc = sfe.prev_dark_start_utc(ref_time_utc=t, sun_alt_dark=-11.75)
    expected_utc = Time('2022-05-28 03:03:16')
    diff_seconds = abs(dark_start_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: SkyfieldEngine object provides default sun altitude defining dark time:
    t = Time('2022-05-28 07:01:00')  # near site's solar midnight
    assert sfe.sun_alt_dark == -10
    dark_start_utc = sfe.prev_dark_start_utc(ref_time_utc=t)
    expected_utc = Time('2022-05-28 02:53:08')
    diff_seconds = abs(dark_start_utc - expected_utc).sec
    assert diff_seconds < 1


# noinspection DuplicatedCode
def test_next_dark_end_utc():
    sfe = make_new_skyfield_engine()

    # Case: User provides sun altitude defining dark time:
    t = Time('2022-05-28 07:01:00')  # near site's solar midnight
    dark_end_utc = sfe.next_dark_end_utc(ref_time_utc=t, sun_alt_dark=-11.75)
    expected_utc = Time('2022-05-28 10:55:18')
    diff_seconds = abs(dark_end_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: SkyfieldEngine object provides default sun altitude defining dark time:
    t = Time('2022-05-28 07:01:00')  # near site's solar midnight
    assert sfe.sun_alt_dark == -10
    dark_end_utc = sfe.next_dark_end_utc(ref_time_utc=t)
    expected_utc = Time('2022-05-28 11:05:26')
    diff_seconds = abs(dark_end_utc - expected_utc).sec
    assert diff_seconds < 1


# noinspection DuplicatedCode
def test_sun_antitransit_utc():
    sfe = make_new_skyfield_engine()

    # Case: near site's solar midnight.
    t = Time('2022-05-28 07:01:00')
    antitransit_utc = sfe.sun_antitransit_utc(t)
    expected_utc = Time('2022-05-28 06:59:23')
    diff_seconds = abs(antitransit_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: near site's solar noon.
    t = Time('2022-05-26 19:01:00')
    antitransit_utc = sfe.sun_antitransit_utc(t)
    expected_utc = Time('2022-05-27 06:59:16')
    diff_seconds = abs(antitransit_utc - expected_utc).sec
    assert diff_seconds < 1


def test_target_transit_utc():
    sfe = make_new_skyfield_engine()
    rigel = SkyCoord(15 * 5.2422978748, -8.20164055, unit='deg')  # Rigel
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-02-08 07:20:00'))
    expected_utc = Time('2022-02-08 03:04:53')
    diff_seconds = abs(transit_utc - expected_utc).sec
    assert diff_seconds < 1
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-07-22 23:20:00'))
    expected_utc = Time('2022-07-22 16:16:07')
    diff_seconds = abs(transit_utc - expected_utc).sec
    assert diff_seconds < 1


# noinspection DuplicatedCode
def test_prev_target_rise_utc():
    # Normal cases:
    sfe = make_new_skyfield_engine()
    rigel = SkyCoord(15 * 5.2422978748, -8.20164055, unit='deg')  # Rigel
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-02-08 07:20:00'))
    prev_rise_utc = sfe.prev_target_rise_utc(rigel, transit_utc, horizon_degrees=30)
    expected_utc = Time('2022-02-08 00:01:23')
    diff_seconds = abs(prev_rise_utc - expected_utc).sec
    assert diff_seconds < 1
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-07-22 23:20:00'))
    prev_rise_utc = sfe.prev_target_rise_utc(rigel, transit_utc, horizon_degrees=30)
    expected_utc = Time('2022-07-22 13:12:36')
    diff_seconds = abs(prev_rise_utc - expected_utc).sec
    assert diff_seconds < 1
    prev_rise_utc = sfe.prev_target_rise_utc(rigel, transit_utc, horizon_degrees=15)
    expected_utc = Time('2022-07-22 11:52:15')
    diff_seconds = abs(prev_rise_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: target always below horizon, never rises:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site_barrow = ini.Site(site_fullpath)
    sfe_barrow = almanac.SkyfieldEngine.from_site(site_barrow, -10)
    antares = SkyCoord(15 * 16.4901280308, -26.43200249, unit='deg')  # Antares
    prev_rise_utc = sfe_barrow.prev_target_rise_utc(antares, transit_utc,
                                                    horizon_degrees=0)
    assert prev_rise_utc is None


# noinspection DuplicatedCode
def test_next_target_set_utc():
    sfe = make_new_skyfield_engine()
    rigel = SkyCoord(15 * 5.2422978748, -8.20164055, unit='deg')  # Rigel
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-02-08 07:20:00'))
    next_set_utc = sfe.next_target_set_utc(rigel, transit_utc, horizon_degrees=30)
    expected_utc = Time('2022-02-08 06:08:23')
    diff_seconds = abs(next_set_utc - expected_utc).sec
    assert diff_seconds < 1
    transit_utc = sfe.target_transit_utc(rigel, Time('2022-07-22 23:20:00'))
    next_set_utc = sfe.next_target_set_utc(rigel, transit_utc, horizon_degrees=30)
    expected_utc = Time('2022-07-22 19:19:39')
    diff_seconds = abs(next_set_utc - expected_utc).sec
    assert diff_seconds < 1
    next_set_utc = sfe.next_target_set_utc(rigel, transit_utc, horizon_degrees=15)
    expected_utc = Time('2022-07-22 20:40:00')
    diff_seconds = abs(next_set_utc - expected_utc).sec
    assert diff_seconds < 1

    # Case: target always below horizon, never rises:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site_barrow = ini.Site(site_fullpath)
    sfe_barrow = almanac.SkyfieldEngine.from_site(site_barrow, -10)
    antares = SkyCoord(15 * 16.4901280308, -26.43200249, unit='deg')  # Antares
    next_set_utc = sfe_barrow.next_target_set_utc(antares, transit_utc,
                                                  horizon_degrees=0)
    assert next_set_utc is None


def test_moon_skycoord():
    sfe = make_new_skyfield_engine()
    sc = sfe.moon_skycoord(Time('2022-02-08 07:20:00'))
    assert sc.ra.deg == pytest.approx(15 * 2.8840125, abs=0.01)
    assert sc.dec.deg == pytest.approx(15.21992, abs=0.01)

    sc = sfe.moon_skycoord(Time('2022-02-28 07:20:00'))
    assert sc.ra.deg == pytest.approx(15 * 20.7681235, abs=0.01)
    assert sc.dec.deg == pytest.approx(-23.25534, abs=0.01)


def test_moon_transit_time():
    sfe = make_new_skyfield_engine()

    # Case: ref time near moon's transit:
    transit_utc = sfe.moon_transit_time(Time('2022-04-28 17:20:00'))
    expected_transit_utc = Time('2022-04-28 17:30:59')
    diff_seconds = (transit_utc - expected_transit_utc).sec
    assert diff_seconds < 5  # skyfield not as accurate for moon as for other targets.

    # Case: ref time near moon's antitransit:
    transit_utc = sfe.moon_transit_time(Time('2022-04-29 05:31:00'))
    expected_transit_utc = Time('2022-04-28 17:30:59')
    diff_seconds = (transit_utc - expected_transit_utc).sec
    assert diff_seconds < 5  # skyfield not as accurate for moon as for other targets.


def test_moonsets_and_moonrises():
    sfe = make_new_skyfield_engine()
    rising = almanac.AlmanacEngine.Event.RISING
    setting = almanac.AlmanacEngine.Event.SETTING

    # Get absence of crossings (empty list) in timespan between crossings:
    crossings = sfe.moonsets_and_moonrises(Time('2022-02-09 09:16:00'),
                                           Time('2022-02-09 17:16:00'))
    assert isinstance(crossings, list)
    assert len(crossings) == 0

    # Get crossings within 48-hour timespan:
    crossings = sfe.moonsets_and_moonrises(Time('2022-02-09 07:16:00'),
                                           Time('2022-02-11 07:16:00'))
    assert isinstance(crossings, list)
    assert len(crossings) == 4
    expected_crossings = [(Time('2022-02-09 08:21:09'), setting),
                          (Time('2022-02-09 18:56:42'), rising),
                          (Time('2022-02-10 09:18:12'), setting),
                          (Time('2022-02-10 19:35:17'), rising)]
    for c, ec in zip(crossings, expected_crossings):
        diff_seconds = abs(c[0] - ec[0]).sec  # compare times
        assert diff_seconds < 5
        assert c[1] == ec[1]  # compare crossing types (setting vs rising)


def test_calc_moon_illumination():
    sfe = make_new_skyfield_engine()

    mi = sfe.moon_illumination(Time('2022-02-09 07:16:00'))
    assert mi == pytest.approx(0.561, abs=0.01)

    mi = sfe.moon_illumination(Time('2022-02-20 04:16:00'))
    assert mi == pytest.approx(0.865, abs=0.01)

    mi = sfe.moon_illumination(Time('2022-02-27 04:16:00'))
    assert mi == pytest.approx(0.161, abs=0.01)


def test_local_sidereal_time():
    sfe = make_new_skyfield_engine()

    lst = sfe.local_sidereal_time(Time('2022-02-08 07:20:00'))
    expected_lst = 9 + 31/60 + 24.821/3600
    assert lst == pytest.approx(expected_lst, abs=0.001)

    lst = sfe.local_sidereal_time(Time('2022-07-11 20:20:11'))
    expected_lst = 8 + 36/60 + 56.967/3600
    assert lst == pytest.approx(expected_lst, abs=0.001)


_____TEST_CLASS_ASTRONIGHT_CONSTRUCTOR_________________________________________ = 0


def test_astronight_attributes_basic():
    # Normal case:
    an = make_new_astronight()
    assert an.site.name == an.site_name
    assert an.site_name == 'New Mexico Skies (Dome)'
    assert an.an_date.an_str == '20220208'
    assert an.sun_altitude_dark == -9.0
    assert isinstance(an.engine, almanac.AlmanacEngine)
    assert abs((an.sun_antitransit_utc - Time('2022-02-09 07:16:17')).sec) < 2
    assert abs((an.sunset_utc - Time('2022-02-09 00:41:13')).sec) < 2
    assert abs((an.sunrise_utc - Time('2022-02-09 13:50:51')).sec) < 2
    assert an.timespan_no_sun.seconds == pytest.approx(47378.631, abs=3)
    assert abs((an.dark_start_utc - Time('2022-02-09 01:21:25')).sec) < 2
    assert abs((an.dark_end_utc - Time('2022-02-09 13:10:43')).sec) < 2
    assert an.timespan_dark.seconds == pytest.approx(42558.126, abs=3)
    assert an.moon_illumination == pytest.approx(0.57, abs=0.01)
    assert an.moon_skycoord.ra.degree == pytest.approx(55.0, abs=0.2)
    assert an.moon_skycoord.dec.degree == pytest.approx(19.3, abs=0.2)
    assert an.moon_skycoord.isscalar == True
    assert (an.moon_transit_utc - Time('2022-02-09 01:18:23')).sec == \
           pytest.approx(0, abs=5)  # skyfield not so accurate re: moon.
    assert abs((an.timespan_dark_no_moon.start - Time('2022-02-09 08:21:09')).sec) < 5
    assert an.timespan_dark_no_moon.end == an.dark_end_utc

    # Case: missing sun_altitude_dark:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    an_1 = almanac.Astronight(site, 20220208)
    assert an_1.sun_altitude_dark == site.sun_altitude_dark


def test_astronight_exceptions():
    # Site ini file missing:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'INVALID.ini')
    with pytest.raises(FileNotFoundError):
        site = ini.Site(site_fullpath)
        _ = almanac.Astronight(site, 20220208, -9.0)
    # Sun always up:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site = ini.Site(site_fullpath)
    with pytest.raises(almanac.SunAlwaysUpError):
        _ = almanac.Astronight(site, 20220622, -9.0)
    # No dark time:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'PointBarrow.ini')
    site = ini.Site(site_fullpath)
    with pytest.raises(almanac.NoDarkTimeError):
        _ = almanac.Astronight(site, 20220815, -9.0)


_____TEST_ASTRONIGHT_METHODS___________________________________________________ = 0


def test_astronight_moon_distance():
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    an = almanac.Astronight(site, 20220402)
    betelgeuse = SkyCoord('05h 55m 10.30536s +07d 24m 25.4304s')
    assert an.moon_distance(betelgeuse) == pytest.approx(53.16, abs=0.1)


def test_astronight_target_transit():
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    an = almanac.Astronight(site, 20220402)
    betelgeuse = SkyCoord('05h 55m 10.30536s +07d 24m 25.4304s')
    transit_time = an.target_transit_utc(betelgeuse)
    assert isinstance(transit_time, Time)
    assert (transit_time - Time('2022-04-03 00:13:12.8')).to_value(u.second) == \
           pytest.approx(0, abs=2)


def test_astronight_target_observable():
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    an = almanac.Astronight(site, 20220402)
    betelgeuse = SkyCoord('05h 55m 10.30536s +07d 24m 25.4304s')

    # Case: moon distant enough:
    ts_obs = an.target_observable(betelgeuse, min_alt=30, min_moon_dist=45)
    assert ts_obs.start == an.dark_start_utc
    assert abs((ts_obs.end - Time('2022-04-03 04:08:12')).sec) < 2

    # Case: moon too close (min_moon_dist changes):
    ts_obs = an.target_observable(betelgeuse, min_alt=30, min_moon_dist=90)
    assert ts_obs.seconds == 0


def test_astronight_time_from_hhmm():
    an = make_new_astronight(an_date=20220208)  # sun antitransit ~ 07:16 UTC.
    assert an.time_from_hhmm('00:00') == Time('2022-02-09 00:00:00')
    assert an.time_from_hhmm('0001') == Time('2022-02-09 00:01:00')   # colon is absent.
    assert an.time_from_hhmm('12:00') == Time('2022-02-09 12:00:00')  # after night.
    assert an.time_from_hhmm('22:34') == Time('2022-02-08 22:34:00')  # before night.
    assert an.time_from_hhmm('19:10') == Time('2022-02-09 19:10:00')  # after night.
    assert an.time_from_hhmm('20:02') == Time('2022-02-08 20:02:00')  # before night.
    assert an.time_from_hhmm('07:16') == Time('2022-02-09 07:16:00')  # ~ sun antitr.


def test_astronight_times_at_sun_alt():
    an = make_new_astronight(an_date=20221210)
    prev, next = an.times_at_sun_alt(sun_alt=-9)
    assert abs(prev - Time('2022-12-11 00:39')).sec < 60
    assert abs(next - Time('2022-12-11 13:11')).sec < 60
    assert an.times_at_sun_alt(sun_alt=-89) == (None, None)


def test_astronight_acp_header_string():
    an = make_new_astronight(an_date=20230115)
    # an = make_new_astronight(an_date=20221213)
    header = an.acp_header_string
    assert len(header) == 3
    # Otherwise, test is largely by inspection.
    print('\n*************')
    for line in header:
        print(line)


def test_astronight__calc_dark_no_moon():
    an = make_new_astronight(an_date=20220208)
    assert (an.timespan_dark_no_moon.start - Time('2022-02-09 08:21:12')).sec == \
        pytest.approx(0, abs=5)
    assert an.timespan_dark_no_moon.end == an.dark_end_utc

    an = make_new_astronight(an_date=20220222)
    assert an.timespan_dark_no_moon.start == an.dark_start_utc
    assert (an.timespan_dark_no_moon.end - Time('2022-02-23 07:23:43')).sec == \
        pytest.approx(0, abs=5)

    an = make_new_astronight(an_date=20211218)  # moon up all night.
    assert an.timespan_dark_no_moon is None

    an = make_new_astronight(an_date=20211203)  # moon down all night.
    assert an.timespan_dark_no_moon == an.timespan_dark

    an = make_new_astronight(an_date=20220628)
    assert an.timespan_dark_no_moon == an.timespan_dark

    an = make_new_astronight(an_date=20220613)
    assert an.timespan_dark_no_moon is None