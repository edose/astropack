__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os

# External packages:

# Author's packages:


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')

# This function is not needed, as:
# (1) it yields apparent Dec only, and
# (2) skyfield and astropy can do this for J2000.0 etc.
#
# def test_ha_dec_from_az_alt():
#     fn = almanac.ha_dec_from_az_alt
#
#     print()
#     ha, dec = fn(latitude=DSW_LAT, az_alt=(90, 89))
#     ha_expected, dec_expected = (-0.0794 * 15, 32.8975)
#     print(ha, ha_expected)
#     print(dec, dec_expected)
#     assert ha == pytest.approx(ha_expected, abs=0.001)
#     assert dec == pytest.approx(dec_expected, abs=0.001)
#
#     print()
#     ha, dec = fn(latitude=DSW_LAT, az_alt=(40, 55))
#     ha_expected = 15 * (-2.6257)
#     dec_expected = 54.57777
#     print(ha, ha_expected)
#     print(dec, dec_expected)
#     # assert ha == pytest.approx(ha_expected, abs=0.001)
#     # assert dec == pytest.approx(dec_expected, abs=0.001)
#
#     print()
#     ha, dec = fn(latitude=DSW_LAT, az_alt=(90, 10))
#     ha_expected = 15 * (-5.4386)
#     dec_expected = 5.4127
#     print(ha, ha_expected)
#     print(dec, dec_expected)
#     assert ha == pytest.approx(ha_expected, abs=0.001)
#     assert dec == pytest.approx(dec_expected, abs=0.001)
#
#     print()
#     ha, dec = fn(latitude=DSW_LAT, az_alt=(300, 77))
#     ha_expected = 15 * 0.9622
#     dec_expected = 38.5891
#     print(ha, ha_expected)
#     print(dec, dec_expected)
#     assert ha == pytest.approx(ha_expected, abs=0.001)
#     assert dec == pytest.approx(dec_expected, abs=0.001)