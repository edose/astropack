""" Module astropak.ini
    Several classes (e.g., Site, Instrument) that read from .ini file and serve them on request.
"""

__author__ = "Eric Dose, Albuquerque"
# Python core:
import os
import os.path
import configparser
from datetime import datetime, timezone
from math import pi, cos


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


__________SITE_INI_____________________________________________________ = 0

# Site file example "NMS_dome.ini":
"""
[Site]
Name = New Mexico Skies (Dome)
MPC Code = N/A

[Location]
# Longitude, Latitude: decimal degrees; Elevation: meters.
Longitude = -105.528978
Latitude = +32.903156
Elevation = 2180
# For correct timezone in standard (not daylight savings/summer time).
UTC Offset = -7
# User discretion, for sufficiently dark to observe.
Sun Altitude Dark = -9

[Climate]
# Coldest date of each year in mm-dd.
Coldest Date = 01-25
# Nominal midnight temperatures (deg C): summer winter
Midnight Temperatures = 20 -3
# Nominal midnight percent humidity: summer winter
Midnight Humidities = 40 60
# Each line: filter summer_extinction winter_extinction
# Approximate values, but that's ok.
Extinctions = Clear 0.18 0.14,
              V     0.20 0.15,
              R     0.16 0.12,
              I     0.11 0.08

[Dome]
Present = True
# Slew rate in (degrees az)/second.
Slew Rate = 2.85
"""


class SiteParseError(Exception):
    """When a line in Site ini file cannot be meaningfuly parsed."""
    pass


class SiteValueError(Exception):
    """When a line can be read, but the extracted value is clearly wrong."""
    pass


class Site:
    """ Holds one telescope site's information. Immutable.
        Gets info from one .ini file.
    """
    def __init__(self, fullpath):
        self.fullpath, self.filename, i = get_ini_data(fullpath)
        self.name = i.get('Site', 'Name')
        self.mpc_code = i.get('Site', 'MPC Code')
        self.longitude = float(i.get('Location', 'Longitude'))  # degrees
        self.latitude = float(i.get('Location', 'Latitude'))    # degrees
        self.elevation = float(i.get('Location', 'Elevation'))  # meters ASL
        self.utc_offset = float(i.get('Location', 'UTC Offset'))  # hours
        self.sun_altitude_dark = float(i.get('Location', 'Sun Altitude Dark'))  # degrees
        coldest_date_strings = i.get('Climate', 'Coldest Date').split('-')    # mm-dd
        try:
            self.coldest_date = int(coldest_date_strings[0]), int(coldest_date_strings[1])
        except ValueError:
            raise SiteParseError(i.get('Climate', 'Coldest Date'))

        midnight_temp_strings = i.get('Climate', 'Midnight Temperatures').split()
        try:
            self.summer_midnight_temperature = int(midnight_temp_strings[0])
            self.winter_midnight_temperature = int(midnight_temp_strings[1])
        except ValueError:
            raise SiteParseError('text = >' + i.get('Climate', 'Midnight Temperatures') + '<')

        midnight_humidity_strings = i.get('Climate', 'Midnight Humidities').split()
        try:
            self.summer_midnight_humidity = int(midnight_humidity_strings[0])
            self.winter_midnight_humidity = int(midnight_humidity_strings[1])
        except ValueError:
            raise SiteParseError('text = >' + i.get('Climate', 'Midnight Humidities') + '<')

        extinction_dict = parse_multiline(i.get('Climate', 'Extinctions'), 3)
        self.extinction = _dict_to_floats(extinction_dict)
        self.dome_present = string_to_boolean(i.get('Dome', 'Present'))
        if self.dome_present:
            self.dome_slew_rate = float(i.get('Dome', 'Slew Rate'))
        else:
            self.dome_present = None

        # Validate values:
        if self.longitude < -180 or self.longitude > 360:
            raise SiteValueError('Longitude must be in range [-180,360]')
        if self.latitude < -90 or self.latitude > 90:
            raise SiteValueError('Latitude must be in range [-90, 90]')
        if self.elevation < 0 or self.elevation > 5000:
            raise SiteValueError('Elevation (meters) must be in range [0, 5000]')
        # TODO: make the following conditional a try/except block, when converting to datetime:
        if self.coldest_date[0] < 1 or self.coldest_date[0] > 12 or \
           self.coldest_date[1] < 1 or self.coldest_date[1] > 31:
            raise SiteValueError('Coldest date must represent valid month-day date of year.')

    def _get_date_phase(self, date):
        """For date [py datetime object], return annual date phase between 0.0 and 1.0."""
        coldest_datetime = datetime(year=date.year,
                                    month=self.coldest_date[0],
                                    day=self.coldest_date[1]).replace(tzinfo=timezone.utc)
        phase = (date - coldest_datetime).days / 365.25
        return phase

    def _interpolate_for_date(self, date, summer_quantity, winter_quantity):
        """Private utility to interpolate quantities between summer and winter."""
        amplitude = summer_quantity - winter_quantity
        mean = (summer_quantity + winter_quantity) / 2
        phase_in_radians = self._get_date_phase(date) * 2 * pi
        quantity_for_date = mean - (amplitude / 2) * cos(phase_in_radians)
        return quantity_for_date

    def midnight_temperature_for_date(self, date):
        """ Return interpolated nominal midnight temperature for date, based on summer and winter values.
        :param: date: date in question, including year. UTC will be assumed if no timezone info
             is given. [py datetime object]
        :return: nominal midnight temperature for date, in deg C. [float]
        """
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        temp_for_date = self._interpolate_for_date(date, self.summer_midnight_temperature,
                                                   self.winter_midnight_temperature)
        return temp_for_date

    def midnight_humidity_for_date(self, date):
        """ Return interpolated nominal midnight humidity for date, based on summer and winter values.
        :param: date: date in question, including year. UTC will be assumed if no timezone info
             is given. [py datetime object]
        :return: nominal midnight relative humidity for date, in percent. [float]
        """
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        humidity_for_date = self._interpolate_for_date(date, self.summer_midnight_humidity,
                                                       self.winter_midnight_humidity)
        return humidity_for_date

    def extinction_for_date(self, date, filter):
        """ Return interpolated extinction for date, based on summer and winter values.
        :param: date: date in question, including year. UTC will be assumed if no timezone info
             is given.[py datetime object]
        :param: filter: filter for which extinction is to be computed. [string]
        :return: extinction [float]."""
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        summer_extinction, winter_extinction = self.extinction[filter][0], self.extinction[filter][1]
        extinction_for_date = self._interpolate_for_date(date, summer_extinction, winter_extinction)
        return extinction_for_date

    def __str__(self):
        return 'Site object: ' + self.name


__________INSTRUMENT_INI_____________________________________________________ = 0

# Instrument file example "BoreaC14.ini":
"""
[Mount]
Model = PlaneWave L-500
Nominal Slew Time = 10

[OTA]
Model = Celestron C14 Edge
Aperture = 0.35
Focal Length = 2710

[Camera]
Model = SBIG STXL-6303E
X Pixels = 3072
Y Pixels = 2047
Pixel Size = 9
# Gain is in electrons/ADU.
CCD Gain = 1.57
Saturation ADU = 54000
Vignetting Pct At Corner = 38
Nominal Cooling Time = 360

[Plate Solution]
Pinpoint Pixel Scale Multiplier = 0.99388

[Filters]
Available = Clear BB SG SR SI
V14 Time To SN 100 = Clear 20,
                     V     50,
                     BB    25
# Transforms = Filter Passband CI_pb1 CI_pb2 1st-order_tr [2nd-order tr] # one only per line
Transforms = Clear SR SR SI   +0.4  -0.6,
             BB    SR SR SI   -0.131

[Scale]
Min FWHM Pixels = 1.5
Max FWHM Pixels = 14
Nominal FWHM Pixels = 7

[Timing]
Exposure Overhead = 20
Max Exposure No Guiding = 119
"""


class TransformParseError(Exception):
    """When transform line cannot be properly parsed."""
    pass


class Instrument:
    """ Holds one instrument's information.
        Gets info from one .ini file.
    """
    def __init__(self, fullpath):
        self.fullpath, self.filename, i = get_ini_data(fullpath)
        self.mount_model = i.get('Mount', 'Model')
        self.nominal_slew_time = float(i.get('Mount', 'Nominal Slew Time'))
        self.ota_model = i.get('OTA', 'Model')
        self.ota_aperture = float(i.get('OTA', 'Aperture'))
        self.focal_length = float(i.get('OTA', 'Focal Length'))
        self.camera_model = i.get('Camera', 'Model')
        self.x_pixels = int(i.get('Camera', 'X Pixels'))
        self.y_pixels = int(i.get('Camera', 'Y Pixels'))
        self.pixel_size = float(i.get('Camera', 'Pixel Size'))
        self.ccd_gain = float(i.get('Camera', 'CCD Gain'))
        self.saturation_adu = float(i.get('Camera', 'Saturation ADU'))
        self.vignetting_pct_at_corner = float(i.get('Camera', 'Vignetting Pct At Corner'))
        self.nominal_cooling_time = float(i.get('Camera', 'Nominal Cooling Time'))
        self.pinpoint_pixel_scale_multipler = float(i.get('Plate Solution',
                                                          'Pinpoint Pixel Scale Multiplier'))
        self.filters_available = tuple(i.get('Filters', 'Available').replace(',', ' ').split())
        self.v14_time_to_sn100 = _dict_to_floats(parse_multiline(i.get('Filters', 'V14 Time To SN 100'),
                                                                 2, 2))
        self.transforms = self._get_transforms(i.get('Filters', 'Transforms'))
        self.min_fwhm_pixels = float(i.get('Scale', 'Min FWHM Pixels'))
        self.max_fwhm_pixels = float(i.get('Scale', 'Max FWHM Pixels'))
        self.nominal_fwhm_pixels = float(i.get('Scale', 'Nominal FWHM Pixels'))
        self.exposure_overhead = float(i.get('Timing', 'Exposure Overhead'))
        self.max_exposure_no_guiding = float(i.get('Timing', 'Max Exposure No Guiding'))

    @staticmethod
    def _get_transforms(multiline_string):
        """ Parse transform dictionary from the Transform multiline string.
        :param multiline_string: transform string from Instrument .ini file
        :return: transform dictionary. Dictionary keys are tuples in form
        (filter, passband, color passband 1, color passband 2), and dictionary values are
        a tuple of float(s) giving transform value.
        """
        min_values_per_transform = 5
        max_values_per_transform = 6
        # transforms_dict = parse_multiline(multiline_string, min_words_per_line=5, max_words_per_line=6)
        strings = multiline_string.splitlines()
        transforms_dict = {}
        for s in strings:
            values = s.replace(',', ' ').split()
            if len(values) < min_values_per_transform or len(values) > max_values_per_transform:
                raise TransformParseError('>' + s + '<')
            key = tuple(values[:4])
            value = tuple([float(v) for v in values[4:]])
            transforms_dict[key] = value
        return transforms_dict


__________HUMANOBSERVER_INI_____________________________________________________________ = 0

# HumanObserver file example "EVD.ini":
"""
[Identity]
Name = Eric Dose

[ALCDEF]
Contact Name = Eric V. Dose
Contact Info = mp@ericdose.com
Observers = Dose, E.V.
"""


class HumanObserver:
    """ Holds one human observer's information.
        Gets info from one .ini file.
    """
    def __init__(self, fullpath):
        self.fullpath, self.filename, i = get_ini_data(fullpath)
        self.name = i.get('Identity', 'Name')
        self.alcdef_contact_name = i.get('ALCDEF', 'Contact Name')
        self.alcdef_contact_info = i.get('ALCDEF', 'Contact Info')
        self.alcdef_observers = i.get('ALCDEF', 'Observers')

    def __str__(self):
        return 'HumanObserver object: ' + self.name


__________UTILITIES____________________________________________________________ = 0


class MultilineParseError(Exception):
    """When multi-line ini field cannot be properly parsed."""
    pass


def get_ini_data(fullpath):
    """ Read .ini file into ConfigParser object. Used at beginning of all astropak .ini file readers."""
    if not (os.path.exists(fullpath) and os.path.isfile(fullpath)):
        raise FileNotFoundError('fullpath >' + fullpath + '<')
    filename = os.path.basename(fullpath)
    ini_data = configparser.ConfigParser()
    ini_data.read(fullpath)
    return fullpath, filename, ini_data


def parse_multiline(multiline_string, min_words_per_line=None, max_words_per_line=None):
    """Parse multiline ini value to dict, where each line's first substring (item) is made the key and
       tuple of remaining substrings (items) is made the value.
    :param multiline_string: multiline string to be parsed. [string]
    :param min_words_per_line: minimum number of items (white-space delimited substrings per line [int],
        or None to allow any number.
    :param max_words_per_line: minimum number of items (white-space delimited substrings per line [int],
        or None to allow any number.
    :return: dictionary where first item of each line is the key, and remaining items of that line
        make up the value as a tuple. [python dict of string: tuple of strings]
    """
    lines = multiline_string.splitlines()
    this_dict = {}
    for line in lines:
        words = line.replace(',', ' ').split()
        if min_words_per_line is not None:
            if len(words) < min_words_per_line:
                raise MultilineParseError('\'' + line + '\'')
        if max_words_per_line is not None:
            if len(words) > max_words_per_line:
                raise MultilineParseError('\'' + line + '\'')
        this_dict[words[0]] = tuple(words[1:])
    return this_dict


def _dict_to_floats(d):
    """Takes dict (e.g., a parsed ini value), converts each item's value to float or to tuple of floats.
    :param d: input dictionary with each item's value being a float, or a string representing a float,
         or a tuple of strings representing floats. [python dict]
    :return: dictionary with same keys as input, and with values having been converted to floats or to
        tuples of floats. [python dict]
    """
    for k, v in d.items():
        if isinstance(d[k], tuple):
            d[k] = tuple([float(value) for value in v])
        else:
            d[k] = float(d[k])
    return d


def string_to_boolean(bool_string, default_value=None):
    """Attempt to extract Boolean meaning from string s, return default_value if failed.
    :param bool_string: string to be interpreted as boolean (interpretable values = 'true', 'yes, 'y',
        'false', 'no', 'n', all case-insensitive). [string]
    :param default_value: value to be returned if strings cannot be interpreted as boolean. [any obj]
    :return: True, False [bool], or default_value.
    """
    sl = bool_string.lower().strip()
    if sl in ['true', 'yes', 'y']:
        return True
    if sl in ['false', 'no', 'n']:
        return False
    return default_value
