"""test_web.py"""

__author__ = "Eric Dose, Albuquerque"


# Python core packages:
import os
from datetime import datetime, timedelta, timezone

# External packages:
import pytest
import pandas as pd
from astropy.time import Time
from astroquery.exceptions import InvalidQueryError

# Author's other modules:
from astropack.ini import Site

# Test target:
import astropack.web as web


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_TOP_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, "tests")


def test_get_mp_ephem():
    site_fullpath = os.path.join(TEST_TOP_DIRECTORY, '$data_for_test', 'NMS_dome.ini')
    site = Site(site_fullpath)
    utc_start = Time('2022-03-01 00:11:22')
    hours_spacing = 2
    n_entries = 6

    # Case: MP id is integer:
    mp_id = 333
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert isinstance(df, pd.DataFrame)
    assert list(df.loc[:, 'Altitude'].values) == [-12, 13, 38, 59, 64, 45]
    assert list(df.columns) == ['Date', 'RA', 'Dec', 'Delta', 'r', 'Elongation',
                                'Phase', 'V', 'Proper motion', 'Direction',
                                'Azimuth', 'Altitude', 'Sun altitude',
                                'Moon phase', 'Moon distance', 'Moon altitude']

    # Case: MP id is string representing integer:
    mp_id = '123'
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert isinstance(df, pd.DataFrame)
    assert str(df['Date'].iloc[0]) == '2022-03-01 00:00:00'
    assert list(df.loc[:, 'Azimuth'].values) == [265, 286, 67, 93, 106, 120]

    # Case: MP id is string representing Name:
    mp_id = 'Badenia'  # (MP number 333)
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert list(df.loc[:, 'Altitude'].values) == [-12, 13, 38, 59, 64, 45]

    # Case: MP id is string representing Designation:
    mp_id = '2004 AA'
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert list(df.loc[:, 'Altitude'].values) == [12, -13, -38, -60, -64, -45]

    # Case: utc_start is string interpretable by astropy Time:
    mp_id = 333
    utc_start = '2022-03-01T00:11:22'
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert str(df['Date'].iloc[0]) == '2022-03-01 00:00:00'
    assert list(df.loc[:, 'Altitude'].values) == [-12, 13, 38, 59, 64, 45]

    # Case: utc_start is py datetime with non-UTC timezone:
    mp_id = 333
    utc_start = datetime(2022, 3, 12, 12, 45, 0).\
        replace(tzinfo=timezone(offset=timedelta(hours=3)))
    df = web.get_mp_ephem(mp_id=mp_id, utc_start=utc_start, step_hours=hours_spacing,
                          num_entries=n_entries, site=site)
    assert str(df['Date'].iloc[0]) == '2022-03-12 12:00:00'
    assert list(df.loc[:, 'Altitude'].values) == [10, -14, -36, -48, -43, -24]

    # Error case: utc_start is wrong type:
    with pytest.raises(TypeError):
        _ = web.get_mp_ephem(mp_id=mp_id, utc_start=345.67, step_hours=hours_spacing,
                             num_entries=n_entries, site=site)

    # Error case: MP id is not matched to a real MP:
    with pytest.raises(InvalidQueryError):
        _ = web.get_mp_ephem(mp_id='Not a real MP ID.',
                             utc_start='2022-03-01T00:11:22', step_hours=hours_spacing,
                             num_entries=n_entries, site=site)


def test_get_mp_info():
    # Case: MP number given:
    d = web.get_mp_info(mp_number=1)
    assert isinstance(d, dict)
    assert len(d) == 6
    assert set(d.keys()) == set(['name', 'number', 'designation', 'period', 'H', 'G'])
    assert d['name'] == 'Ceres'
    assert d['number'] == 1
    assert d['designation'] is None
    assert float(d['period']) == pytest.approx(4.6, abs=0.2)
    assert float(d['H']) == pytest.approx(3.52, abs=0.2)
    assert d['G'] == '0.15'

    # Case: MP name given:
    d = web.get_mp_info(mp_name='Ceres')
    assert set(d.keys()) == set(['name', 'number', 'designation', 'period', 'H', 'G'])
    assert d['number'] == 1
    assert float(d['H']) == pytest.approx(3.52, abs=0.2)

    # Case: Return None for non-existent MP:
    d = web.get_mp_info(mp_name='NotMP')
    assert d is None

    # Error cases:
    with pytest.raises(TypeError):
        _ = web.get_mp_info(mp_number='hahaha')
        _ = web.get_mp_info(mp_name=333)
    with pytest.raises(ValueError):
        _ = web.get_mp_info()
