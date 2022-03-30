"""
    Constants: numerical, astronomical, datetime formats, and special characters.
"""

__author__ = "Eric Dose, Albuquerque"

import math


# NUMERICAL CONSTANTS:
RADIANS_PER_DEGREE = math.pi / 180.0
"""float, about 1/57.396"""
DEGREES_PER_RADIAN = 1 / RADIANS_PER_DEGREE          # ca. 57.396
"""float, about 57.396"""
ARCSECONDS_PER_RADIAN = 3600 * DEGREES_PER_RADIAN    # ca. 206265
"""float, about 206265"""
FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2))  # ca. 2.35482
"""float, about 2.35482"""
DAYS_PER_YEAR_NOMINAL = 365.25
"""float, equals 365.25"""


# ASTRONOMICAL and PHYSICS QUANTITIES:
SECONDS_PER_SIDEREAL_DAY = 86164.0905
"""float, equals 86164.0905"""


# MISCELLANEOUS constants and conventions:
FORMAT_ISO_8601 = '%Y-%m-%dT%H:%M:%S.%f'     # strict (scientific & automation)
"""string, strict ISO 8601, same as 'FITS' format, 
yielding e.g., '2022-04-02T12:34:42.432'"""
FORMAT_ISO_8601_HUMAN = '%Y-%m-%d %H:%M:%S'  # human-facing
"""string, common ISO 8601, slightly more human-readable format, 
yielding e.g., '2022-04-02 12:34:42.432'"""


# SPECIAL CHARACTERS (unicode), esp. for weather and astronomy:
CH_DEGREE = '\u00B0'
"""\u00B0  Unicode degree sign"""
CH_PLUS_MINUS = '\u00B1'
"""\u00B1  Unicode plus-minus sign"""
CH_MICRO = '\u00B5'
"""\u00B5  Unicode micro (mu) sign"""
CH_MOON = '\u263D'  # open to left. Open to right = '\u263E'
"""\u263D  Unicode moon character, works with matplotlib graphics package 
(unlike \\U0001F319 = \U0001F319)"""
# MOON_CHARACTER = '\U0001F319'  # Matplotlib says 'glyph missing from current font'.
CH_COMET = '\u2604'
"""\u2604  Unicode comet character"""
CH_STAR_WHITE = '\u2606'
"""\u2606  Unicode white (empty) star character"""
CH_STAR_BLACK = '\u2605'
"""\u2605  Unicode black (filled) star character"""
CH_SUN_BLACK = '\u2600'
"""\u2600  Unicode sun character"""
CH_CLOUD_BLACK = '\u2601'
"""\u2601  Unicode cloud character"""
CH_BULLET_POINT = '\u2022'
"""\u2022  Unicode bullet point (unordered list) character"""
CH_INFINITY = '\u221E'
"""\u221E  Unicode infinity character"""
CH_COPYRIGHT = '\u00A9'
"""\u00A9  Copyright character"""
CH_MUCH_LESS_THAN = '\u00AB'
"""\u00AB  Much-less-than sign"""
CH_MUCH_GREATER_THAN = '\u00BB'
"""\u00BB  Much-greater-than sign"""
CH_ALMOST_EQUAL = '\u2248'
"""\u2248  Almost-equal sign"""
CH_APPROX_EQUAL = '\u2245'
"""\u2245  Approximately-equal sign"""
CH_NOT_EQUAL = '\u2260'
"""\u2260  Not-equal sign"""
CH_IDENTICAL_TO = '\u2261'
"""\u2261  Identical-to sign (also used to mean 'is defined as')"""
