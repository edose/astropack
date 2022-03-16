__author__ = "Eric Dose, Albuquerque"

""" test_ini.py """

# Python core:
import os
from datetime import datetime, timezone

# External packages:
import pytest

# Test target:
from astropak import ini


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_FOR_TEST_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'test', '$data_for_test')

__________TEST_INI_UTILITIES___________________________________ = 0


def test_parse_multiline():
    """OK 2022-03-10."""
    # Normal case:
    lines = """First CCC 0.18
    Second XZXC haha"""
    d = ini.parse_multiline(lines, min_words_per_line=3, max_words_per_line=3)
    assert len(d) == 2
    assert d['First'] == ('CCC', '0.18')
    d = ini.parse_multiline(lines, min_words_per_line=2, max_words_per_line=5)
    assert len(d) == 2
    assert d['Second'] == ('XZXC', 'haha')
    # Error: invalid number of words (substrings) per line:
    with pytest.raises(ini.MultilineParseError) as e:
        d = ini.parse_multiline(lines, min_words_per_line=1, max_words_per_line=2)
    with pytest.raises(ini.MultilineParseError) as e:
        d = ini.parse_multiline(lines, min_words_per_line=4, max_words_per_line=100)


def test_dict_to_floats():
    """OK 2022-03-10."""
    # Normal case:
    d = {'DSS': '0.18', 'HA': '-5.4'}
    f = ini._dict_to_floats(d)
    assert len(f) == 2
    assert f['HA'] == float('-5.4')
    # Error: value doesn't represent float:
    d = {'DSS': '0.18', 'HA': 'XXX'}
    with pytest.raises(ValueError) as e:
        f = ini._dict_to_floats(d)


def test_string_to_boolean():
    """OK 2022-03-10."""
    # Normal case:
    assert ini.string_to_boolean('True') == True
    assert ini.string_to_boolean('False') == False
    assert ini.string_to_boolean('Yes') == True
    assert ini.string_to_boolean('No') == False
    # Non-boolean string case:
    assert ini.string_to_boolean('hahaha', True) == True
    assert ini.string_to_boolean('hahaha', False) == False
    assert ini.string_to_boolean('hahaha') is None


def test_class_site():
    """OK 2022-03-10."""

    fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'NMS_dome.ini')
    site = ini.Site(fullpath)
    assert site.fullpath == fullpath
    assert site.filename == os.path.basename(fullpath)
    assert site.name == 'New Mexico Skies (Dome)'
    assert site.longitude == float('-105.528978')
    assert site.latitude == float('+32.903156')
    assert site.elevation == 2180.0
    assert site.utc_offset == -7.0
    assert site.coldest_date == (1, 25)
    assert site.summer_midnight_temperature == 20.0
    assert site.winter_midnight_temperature == -5.0
    assert site.midnight_temperature_for_date(datetime(2022, 4, 25).replace(tzinfo=timezone.utc)) == \
           pytest.approx(7.21779671198)
    assert site.summer_midnight_humidity == 40.0
    assert site.winter_midnight_humidity == 60.0
    assert site.midnight_humidity_for_date(datetime(2022, 5, 25).replace(tzinfo=timezone.utc)) == \
           pytest.approx(45.2629364962)
    assert site.extinction['Clear'] == (float('0.18'), float('0.14'))
    assert site.extinction_for_date(datetime(2022, 2, 15).replace(tzinfo=timezone.utc), 'Clear') == \
           pytest.approx(0.14129089)
    assert site.dome_present == True
    assert site.dome_slew_rate == float('2.85')

    # Case: if date passed in has no timezone info, ensure that it's given timezone UTC:
    assert site.midnight_temperature_for_date(datetime(2022, 4, 25)) == \
           site.midnight_temperature_for_date(datetime(2022, 4, 25).replace(tzinfo=timezone.utc))
    assert site.midnight_humidity_for_date(datetime(2022, 5, 25).replace(tzinfo=timezone.utc)) == \
           site.midnight_humidity_for_date(datetime(2022, 5, 25))
    assert site.extinction_for_date(datetime(2022, 2, 15), 'Clear') == \
           site.extinction_for_date(datetime(2022, 2, 15).replace(tzinfo=timezone.utc), 'Clear')


def test_class_instrument():
    fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'BoreaC14.ini')
    # Test static method _get_transforms():
    _, _, i = ini.get_ini_data(fullpath)
    transform_text = i.get('Filters', 'Transforms')
    transforms = ini.Instrument._get_transforms(transform_text)
    assert len(transforms) == 2
    assert transforms[('Clear', 'SR', 'SR', 'SI')] == (float('+0.4'), float('-0.6'))
    # Test class:
    fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'BoreaC14.ini')
    inst = ini.Instrument(fullpath)
    assert inst.fullpath == fullpath
    assert inst.filename == os.path.basename(fullpath)
    assert inst.mount_model == 'PlaneWave L-500'
    assert inst.nominal_slew_time == 10.0
    assert inst.ota_model == 'Celestron C14 Edge'
    assert inst.ota_aperture == float('0.35')
    assert inst.focal_length == 2710.0
    assert inst.camera_model == 'SBIG STXL-6303E'
    assert inst.x_pixels == 3072
    assert inst.y_pixels == 2047
    assert inst.pixel_size == 9.0
    assert inst.ccd_gain == float('1.57')
    assert inst.saturation_adu == 54000.0
    assert inst.vignetting_pct_at_corner == 38.0
    assert inst.nominal_cooling_time == 360.0
    assert inst.pinpoint_pixel_scale_multipler == float(0.99388)
    assert inst.filters_available == ('Clear', 'BB', 'SG', 'SR', 'SI')
    assert len(inst.transforms) == 2
    assert inst.transforms[('Clear', 'SR', 'SR', 'SI')] == (float('+0.4'), float('-0.6'))
    assert inst.transforms[('BB', 'SR', 'SR', 'SI')] == (float('-0.131'), )
    assert inst.min_fwhm_pixels == float('1.5')
    assert inst.max_fwhm_pixels == 14.0
    assert inst.nominal_fwhm_pixels == 7.0
    assert inst.exposure_overhead == 20.0
    assert inst.max_exposure_no_guiding == 119.0


def test_class_humanobserver():
    fullpath = os.path.join(DATA_FOR_TEST_DIRECTORY, 'EVD.ini')
    obs = ini.HumanObserver(fullpath)
    assert obs.fullpath == fullpath
    assert obs.filename == os.path.basename(fullpath)
    assert obs.name == 'Eric Dose'
    assert obs.alcdef_contact_name == 'Eric V. Dose'
    assert obs.alcdef_contact_info == 'mp@ericdose.com'
    assert obs.alcdef_observers == 'Dose, E.V.'
