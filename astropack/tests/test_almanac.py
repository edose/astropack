""" test_almanac.py
    All test functions and Astronight class appear to pass all tests as of 2022-02-09.
"""

__author__ = "Eric Dose, Albuquerque"


# Python core:
import os

# External packages:
import pytest
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord
from skyfield.api import load, Star

# Author's packages:
from astropak import ini
from astropak.util import Timespan

# Test target:
from astropak import almanac

THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FOR_TEST_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')


__________TEST_FUNCTIONS__________________________________________________________ = 0


def test__skycoords_to_skyfield_star_list():
    # ***** OK 2022-02-08.
    # Case: one scalar SkyCoord passed in:
    sc = SkyCoord(10.0, 20.0, unit="deg")
    sl = almanac._skycoords_to_skyfield_star_list(sc)
    assert isinstance(sl, list)
    assert len(sl) == 1
    assert sl[0].ra.hours * 15.0 == pytest.approx(10.0)
    del sc, sl
    # Case: one array SkyCoord passed in:
    sc = SkyCoord([10.0, 11.0, 13.0], [20.0, 21.0, 55.0], unit="deg")
    sl = almanac._skycoords_to_skyfield_star_list(sc)
    assert isinstance(sl, list)
    assert len(sl) == 3
    assert sl[2].dec.degrees == pytest.approx(55.0)
    del sc, sl
    # Case: list of two scalar SkyCoords passed in:
    sc = [SkyCoord(10.0, 20.0, unit="deg"), SkyCoord(15.0, 24.0, unit="deg")]
    sl = almanac._skycoords_to_skyfield_star_list(sc)
    assert isinstance(sl, list)
    assert len(sl) == 2
    assert sl[1].dec.degrees == pytest.approx(24.0)
    # Case: list of both scalar and array Skycoords passed in:
    sc = [SkyCoord(10.0, 20.0, unit="deg"),
          SkyCoord([10.0, 11.0, 13.0], [20.0, 21.0, 55.0], unit="deg"),
          SkyCoord(111.0, 88.0, unit="deg"),
          SkyCoord([17.0, 18.0], [27.0, 25.0], unit="deg")]
    sl = almanac._skycoords_to_skyfield_star_list(sc)
    assert isinstance(sl, list)
    assert len(sl) == 7
    assert sl[5].dec.degrees == pytest.approx(27.0)


def test__skyfield_time_to_astropy():
    # ***** OK 2022-02-08.
    ts = load.timescale()
    t_sf = ts.utc(2019, 9, 30, 23, 5, 6)
    t_apy = almanac._skyfield_time_to_astropy(t_sf)
    assert isinstance(t_apy, Time)
    assert t_apy.scale == 'utc'
    assert t_apy.format == 'iso'
    assert (t_apy - Time('2019-09-30 23:05:06')).sec == pytest.approx(0, abs=1e-9)


def test_moon_illumination_pct():
    # ***** OK 2022-02-08.
    an = make_new_astronight()
    mi = almanac.moon_illumination_pct(an.master_eph, an.timescale, Time('2022-02-09 07:16:00'))
    assert mi == pytest.approx(56.9, abs=0.2)


def test_local_sidereal_time():
    # ***** OK 2022-02-08.
    timescale = load.timescale()
    lst = almanac.local_sidereal_time(-105.5, timescale, Time('2022-02-08 07:20:00'))
    assert lst == pytest.approx(9.5, abs=0.1)


def test__closest_transit():
    # ***** OK 2022-02-08.
    timescale = load.timescale()
    times = [timescale.utc(2019, 9, 30, hour, 5, 8) for hour in range(6)]
    # print()
    #  for time in times:
    #      print(str(time))
    values = [0, 1, 0, 1, 0, 1]
    target_time = timescale.utc(2019, 9, 30, 2, 44, 44)
    # print('target_time =', str(target_time))
    ct = almanac._closest_transit(times, values, target_time)
    # print('ct =', str(ct))
    assert ct == times[3]


def test_target_transit_time():
    # ***** OK 2022-02-09.
    an = make_new_astronight()
    # Case: one scalar SkyCoord object passed in:
    sc = SkyCoord(15 * 7.58, 25.15, unit='deg')  # MP 830
    target_time = Time('2022-02-09 07:00:00')
    transit = almanac.target_transit_time(an.obs, an.master_eph, an.timescale, sc, target_time)
    assert isinstance(transit, Time)  # scalar astropy Time result.
    assert (transit - Time('2022-02-09 05:21:07')).sec == pytest.approx(0, abs=5)
    # Case: one array SkyCoord object passed in:
    sc = SkyCoord([15 * 7.58, 15 * 12.462], [25.15, 6.556], unit='deg')  # MPs 830, 2802
    target_time = Time('2022-02-09 07:00:00')
    transit = almanac.target_transit_time(an.obs, an.master_eph, an.timescale, sc, target_time)
    assert isinstance(transit, list)  # list of astropy Time results.
    assert len(transit) == 2
    assert (transit[1] - Time('2022-02-09 10:13:01')).sec == pytest.approx(0, abs=5)


def test_moon_transit_time():
    # ***** OK 2022-02-09.
    an = make_new_astronight()
    mtt = almanac.moon_transit_time(an.obs, an.master_eph, an.timescale, an.local_middark_utc)
    assert (mtt - Time('2022-02-09 01:18:25')).sec == pytest.approx(0, abs=5)
    mtt = almanac.moon_transit_time(an.obs, an.master_eph, an.timescale,
                                    an.local_middark_utc + TimeDelta(24 * 3600, format='sec'))
    # print(str(mtt))
    assert (mtt - Time('2022-02-10 02:04:37')).sec == pytest.approx(0, abs=5)


def test_moon_ra_dec():
    # ***** OK 2022-02-09.
    an = make_new_astronight()
    ra, dec = almanac.moon_ra_dec(an.obs, an.master_eph, an.timescale, an.local_middark_utc)
    assert ra == pytest.approx(55.1, abs=0.2)
    assert dec == pytest.approx(19.3, abs=0.2)
    ra, dec = almanac.moon_ra_dec(an.obs, an.master_eph, an.timescale,
                                  an.local_middark_utc + TimeDelta(24 * 3600, format='sec'))
    assert ra == pytest.approx(67.3, abs=0.2)
    assert dec == pytest.approx(22.6, abs=0.2)


def test_find_target_up_down():
    # ***** OK 2022-02-09.
    an = make_new_astronight()
    with pytest.raises(ValueError):
        _ = almanac.find_target_up_down(an.obs, an.master_eph, an.master_eph['moon'], an.timescale,
                                        an.timespan_dark, up_down='INVALID')
    test_timespan = Timespan(Time('2022-02-09 21:00:00'), Time('2022-02-10 21:00:00'))
    # print()

    # ===== Moon, up and down:
    # print('\nMoon:')
    moon_set_actual = Time('2022-02-10 09:22:55')   # within a second.
    moon_rise_actual = Time('2022-02-10 19:30:41')  # "
    moon_up = almanac.find_target_up_down(an.obs, an.master_eph, an.master_eph['moon'], an.timescale,
                                          test_timespan, 'up', horizon=almanac.HORIZON_USNO)
    # print(str(moon_up))
    assert isinstance(moon_up, list)
    assert len(moon_up) == 2
    assert moon_up[0].start == test_timespan.start
    assert (moon_up[0].end - moon_set_actual).sec == pytest.approx(0, abs=2)
    assert (moon_up[1].start - moon_rise_actual).sec == pytest.approx(0, abs=2)
    assert moon_up[1].end == test_timespan.end
    moon_down = almanac.find_target_up_down(an.obs, an.master_eph, an.master_eph['moon'], an.timescale,
                                            test_timespan, 'down', horizon=almanac.HORIZON_USNO)
    # print(str(moon_down))
    assert isinstance(moon_down, list)
    assert len(moon_down) == 1
    assert (moon_down[0].start - moon_set_actual).sec == pytest.approx(0, abs=2)
    assert (moon_down[0].end - moon_rise_actual).sec == pytest.approx(0, abs=2)

    # ===== Sun, up and down:
    # print('\nSun:')
    sun_set_actual = Time('2022-02-10 00:42:08')
    sun_rise_actual = Time('2022-02-10 13:49:59')
    sun_up = almanac.find_target_up_down(an.obs, an.master_eph, an.master_eph['sun'], an.timescale,
                                         test_timespan, 'up', horizon=almanac.HORIZON_USNO)
    # print(str(sun_up))
    assert isinstance(sun_up, list)
    assert len(sun_up) == 2
    assert sun_up[0].start == test_timespan.start
    assert (sun_up[0].end - sun_set_actual).sec == pytest.approx(0, abs=2)
    assert (sun_up[1].start - sun_rise_actual).sec == pytest.approx(0, abs=2)
    assert sun_up[1].end == test_timespan.end
    sun_down = almanac.find_target_up_down(an.obs, an.master_eph, an.master_eph['sun'], an.timescale,
                                           test_timespan, 'down', horizon=almanac.HORIZON_USNO)
    # print(str(sun_down))
    assert isinstance(sun_down, list)
    assert len(sun_down) == 1
    assert (sun_down[0].start - sun_set_actual).sec == pytest.approx(0, abs=2)
    assert (sun_down[0].end - sun_rise_actual).sec == pytest.approx(0, abs=2)

    # ===== Star Procyon, rises and sets at this observatory:
    # print('\nProcyon:')
    procyon_rise_actual = Time('2022-02-09 23:05:07')
    procyon_set_actual = Time('2022-02-10 11:37:54')
    procyon = Star(ra_hours=(7, 39, 18), dec_degrees=(5, 13, 30))
    procyon_up = almanac.find_target_up_down(an.obs, an.master_eph, procyon, an.timescale,
                                             test_timespan, 'up', horizon=almanac.HORIZON_USNO)
    # print(str(procyon_up))
    assert isinstance(procyon_up, list)
    assert len(procyon_up) == 1
    assert (procyon_up[0].start - procyon_rise_actual).sec == pytest.approx(0, abs=2)
    assert (procyon_up[0].end - procyon_set_actual).sec == pytest.approx(0, abs=2)
    procyon_down = almanac.find_target_up_down(an.obs, an.master_eph, procyon, an.timescale,
                                               test_timespan, 'down', horizon=almanac.HORIZON_USNO)
    # print(str(procyon_down))
    assert isinstance(procyon_down, list)
    assert len(procyon_down) == 2
    assert procyon_down[0].start == test_timespan.start
    assert (procyon_down[0].end - procyon_rise_actual).sec == pytest.approx(0, abs=2)
    assert (procyon_down[1].start - procyon_set_actual).sec == pytest.approx(0, abs=2)
    assert procyon_down[1].end == test_timespan.end

    # ===== Star Kochab, always up at this observatory:
    # print('\nKochab:')
    kochab = Star(ra_hours=(14, 50, 42), dec_degrees=(74, 9, 20))
    kochab_up = almanac.find_target_up_down(an.obs, an.master_eph, kochab, an.timescale,
                                            test_timespan, 'up', horizon=almanac.HORIZON_USNO)
    # print(str(kochab_up))
    assert kochab_up == [test_timespan]
    kochab_down = almanac.find_target_up_down(an.obs, an.master_eph, kochab, an.timescale,
                                              test_timespan, 'down', horizon=almanac.HORIZON_USNO)
    # print(str(kochab_down))
    assert kochab_down is None

    # ===== Star Gatria, always down at this observatory:
    # print('\nGatria:')
    gatria = Star(ra_hours=(15, 18, 55), dec_degrees=(-68, 40, 46))
    gatria_up = almanac.find_target_up_down(an.obs, an.master_eph, gatria, an.timescale,
                                            test_timespan, 'up', horizon=almanac.HORIZON_USNO)
    # print(str(gatria_up))
    assert gatria_up is None
    gatria_down = almanac.find_target_up_down(an.obs, an.master_eph, gatria, an.timescale,
                                              test_timespan, 'down', horizon=almanac.HORIZON_USNO)
    # print(str(gatria_down))
    assert gatria_down == [test_timespan]


def test_find_no_sun():
    # ***** OK 2022-02-09.
    an = make_new_astronight()
    fn = almanac.calc_timespan_no_sun(an.obs, an.master_eph, an.timescale, an.local_middark_utc)
    assert isinstance(fn, Timespan)
    # print(str(fn))
    assert (fn.start - Time('2022-02-09 00:41:13')).sec == pytest.approx(0, abs=2)
    assert (fn.end - Time('2022-02-09 13:50:52')).sec == pytest.approx(0, abs=2)


def test_make_skyfield_observatory_from_site():
    # ***** OK 2022-02-09.
    from skyfield.toposlib import GeographicPosition
    an = make_new_astronight()
    assert an.site_name == 'New Mexico Skies (Dome)'
    obs = almanac.make_skyfield_observatory_from_site(an.site)
    assert isinstance(obs, GeographicPosition)
    assert obs.model.name == 'WGS84'
    assert obs.latitude.degrees == pytest.approx(an.site.latitude)
    assert obs.longitude.degrees == pytest.approx(an.site.longitude)
    assert obs.elevation.m == pytest.approx(an.site.elevation)


def test_ha_dec_from_az_alt():
    fn = almanac.ha_dec_from_az_alt
    result = fn(latitude=32.9, az_alt=(45, 45))
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert fn(latitude=32.9, az_alt=(45, 45)) == pytest.approx((57.2093, 53.5034), abs=0.0002)
    assert fn(latitude=32.9, az_alt=(315, 45)) == pytest.approx((-57.2093, 53.5034), abs=0.0002)
    result = fn(latitude=32.9, az_alt=[(45, 45), (0, 60)])
    assert isinstance(result, list)
    assert len(result) == 2
    assert all([isinstance(x, tuple) for x in result])
    assert result[1] == pytest.approx((0, 62.9), abs=0.0002)


__________TEST_CLASS_ASTRONIGHT__________________________________________________________ = 0


def test_astronight_constructor():
    # ***** OK 2022-02-08.
    # Normal case:
    an = make_new_astronight()
    assert an.site_name == 'New Mexico Skies (Dome)'
    assert an.an_date_string == '20220208'
    assert an.sun_altitude_dark == -9.0
    assert an.obs.elevation.m == an.site.elevation
    assert an.master_eph.filename == 'de440s.bsp'
    assert str(type(an.timescale)) == '<class \'skyfield.timelib.Timescale\'>'
    assert an.timespan_no_sun.seconds == pytest.approx(47378.6, abs=2)
    assert (an.timespan_no_sun.midpoint - Time('2022-02-09 07:16:02.284')).sec == pytest.approx(0, abs=2)
    assert an.timespan_dark.seconds == pytest.approx(42558.1, abs=2)
    assert (an.timespan_dark.midpoint - Time('2022-02-09 07:16:04.796')).sec == pytest.approx(0, abs=2)
    assert (an.local_middark_utc - Time('2022-02-09 07:16:04.796')).sec == pytest.approx(0, abs=2)
    assert an.local_middark_lst_hour_string == '09:31:26.25'
    assert an.moon_illumination == pytest.approx(57.0, abs=1.0)
    assert an.moon_ra == pytest.approx(55.0, abs=0.2)  # degrees
    assert an.moon_dec == pytest.approx(19.3, abs=0.2)  # degrees
    assert an.moon_skycoord.isscalar == True
    assert an.moon_skycoord.ra.degree == pytest.approx(55.0, abs=2.0)  # degrees
    assert (an.moon_transit - Time('2022-02-09 01:18:25')).sec == pytest.approx(0, abs=2)
    assert str(type(an.moon_up_timespans[0])) == '<class \'astropak.util.Timespan\'>'
    assert an.moon_up_timespans[0].seconds == pytest.approx(25456, abs=2)
    assert an.timespan_dark_no_moon.seconds == pytest.approx(17102, abs=2)
    # Case: missing sun_altitude_dark:
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    an_1 = almanac.Astronight(site, 20220208)
    assert an_1.sun_altitude_dark == site.sun_altitude_dark


def test_astronight_exceptions():
    # ***** OK 2022-02-08.
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
    # .timespan_observable() absent min_alt or min_moon_dist:
    an = make_new_astronight()
    sc = SkyCoord(10, 20, unit="deg")
    with pytest.raises(ValueError):
        _ = an.timespan_observable(sc)


__________TEST_SETUP_and_SUPPORT_________________________________________________ = 0


def make_new_astronight():
    """Returns a typical Astronight object."""
    site_fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_Dome.ini')
    site = ini.Site(site_fullpath)
    return almanac.Astronight(site, 20220208, -9.0)
