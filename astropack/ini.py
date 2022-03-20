""" Module astropack.ini
    Several classes (e.g., Site, Instrument) that read from .ini file, and
    serve them on request.
"""

__author__ = "Eric Dose, Albuquerque"
# Python core:
import os
import os.path
import configparser
from datetime import datetime, timezone
from math import pi, cos


THIS_PACKAGE_ROOT_DIRECTORY = \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

__all__ = ['SiteParseError', 'SiteValueError', 'Site',
           'TransformParseError', 'Instrument',
           'HumanObserver',
           'MultilineParseError',
           'get_ini_data',
           'string_to_boolean']


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
    """Raised if a line in Site .ini file cannot be meaningfuly parsed."""
    pass


class SiteValueError(Exception):
    """Raised if a line in Site .ini has been read,
    but the extracted value is clearly wrong."""
    pass


class Site:
    """ Holds and makes available one telescope/observer site's information.
    Loads info from one .ini file.

    Parameters
    ----------
    fullpath : str
        Full path to site's .ini file.

    Attributes
    ----------
    fullpath : str
        Value of input parameter ``fullpath``.
    name : str
        Long name of site, e.g., 'New Mexico Skies (Dome)'
    mpc_code : str
        Minor Planet Center (IAU) location code for this site, e.g., 'H12'.
    longitude : float
        Site's longitude, in degrees.
    latitude : float
        Site's latitude, in degrees.
    elevation : float
        Site's elevation above mean sea level, in meters.
    utc_offset : float
        Site's local standard time offset from UTC time.

        The most likely use is
        in estimating approximate midnight at a given longitude,
        and is most important to distinguish the sign of large offsets from UTC,
        that is, for longitudes near the International Date Line.
        At those longitudes, date shift from UTC cannot be guessed from longitude alone.
        Thus, the value need be only approximate, but special attention must be
        paid for longitudes  approaching -180 or +180 degrees.
    sun_altitude_dark : float
        The maximum (least negative) sun altitude, in degrees, giving a sky
        sufficiently dark to allow observations.

        The value is at user discretion, and typically in the range [-8, -12].
    coldest_date : tuple of 2 int
        Tuple (month, day) representing the average date with the coldest midnight
        of each year.

        In the Northern Hemisphere, this is usually in the range (1, 23) to (2, 1).
        ``coldest_date`` is part of this class's estimation of midnight temperature,
        humidity, and extinctions for any date of the year.
    summer_midnight_temperature : float
        Estimated mean summer midnight temperature, in degrees C.
    winter_midnight_temperature : float
        Estimated mean winter midnight temperature, in degrees C.
    summer_midnight_humidity : float
        Estimated mean summer midnight relative humidity percent, in range [0, 100].
    winter_midnight_humidity : float
        Estimated mean winter midnight relative humidity percent, in range [0, 100].
    extinction : dict (str: tuple of float)
        Summer and winter estimated mean extinction values for site, by filter.
        If ``extinction``['V'] = (0.18, 0.14), the site's extinction values in 'V'
        filter are expected to be close to 0.18 in summer and 0.14 in winter.
    dome_present : bool
        If True, the site's telescope is housed in a dome.
    dome_slew_rate : float
        The dome's rotation rate when slewing, in degrees per second.

        Value is frequently near 3 for dome's of diameter 3-5 meters.
    """

    def __init__(self, fullpath):
        self.fullpath, self.filename, i = get_ini_data(fullpath)
        self.name = i.get('Site', 'Name')
        self.mpc_code = i.get('Site', 'MPC Code')
        self.longitude = float(i.get('Location', 'Longitude'))  # degrees
        self.latitude = float(i.get('Location', 'Latitude'))    # degrees
        self.elevation = float(i.get('Location', 'Elevation'))  # meters ASL
        self.utc_offset = float(i.get('Location', 'UTC Offset'))  # hours
        self.sun_altitude_dark = float(i.get('Location',
                                             'Sun Altitude Dark'))  # degrees
        coldest_date_strings = i.get('Climate', 'Coldest Date').split('-')    # mm-dd
        try:
            self.coldest_date = int(coldest_date_strings[0]), \
                                int(coldest_date_strings[1])
        except ValueError:
            raise SiteParseError(i.get('Climate', 'Coldest Date'))

        midnight_temp_strings = i.get('Climate', 'Midnight Temperatures').split()
        try:
            self.summer_midnight_temperature = int(midnight_temp_strings[0])
            self.winter_midnight_temperature = int(midnight_temp_strings[1])
        except ValueError:
            raise SiteParseError('text = >' + i.get('Climate',
                                                    'Midnight Temperatures') + '<')

        midnight_humidity_strings = i.get('Climate', 'Midnight Humidities').split()
        try:
            self.summer_midnight_humidity = int(midnight_humidity_strings[0])
            self.winter_midnight_humidity = int(midnight_humidity_strings[1])
        except ValueError:
            raise SiteParseError('text = >' + i.get('Climate',
                                                    'Midnight Humidities') + '<')

        extinction_dict = _parse_multiline(i.get('Climate', 'Extinctions'), 3)
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
        # TODO: wrap this in a try/except block, when converting to datetime:
        if self.coldest_date[0] < 1 or self.coldest_date[0] > 12 or \
           self.coldest_date[1] < 1 or self.coldest_date[1] > 31:
            raise SiteValueError('Coldest date must represent valid '
                                 'month-day date of year.')

    def _get_date_phase(self, date):
        """For ``date`` [py datetime], return annual date phase."""
        coldest_datetime = datetime(year=date.year,
                                    month=self.coldest_date[0],
                                    day=self.coldest_date[1]).\
            replace(tzinfo=timezone.utc)
        phase = (date - coldest_datetime).days / 365.25
        return phase

    def _interpolate_for_date(self, date, summer_quantity, winter_quantity):
        """Private utility to interpolate quantities between
        summer and winter nominal values."""
        amplitude = summer_quantity - winter_quantity
        mean = (summer_quantity + winter_quantity) / 2
        phase_in_radians = self._get_date_phase(date) * 2 * pi
        quantity_for_date = mean - (amplitude / 2) * cos(phase_in_radians)
        return quantity_for_date

    def midnight_temperature_for_date(self, date):
        """Return interpolated nominal midnight temperature for site and date,
        based on summer and winter values.

        Parameters
        ----------
        date : datetime
            Date, including year, for which midnight temperature estimate is needed.
            UTC will be assumed if no timezone info is given.

        Returns
        -------
        temp_for_date : float
            Midnight temperature estimate for date, in degrees C.
        """
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        temp_for_date = self._interpolate_for_date(date,
                                                   self.summer_midnight_temperature,
                                                   self.winter_midnight_temperature)
        return temp_for_date

    def midnight_humidity_for_date(self, date):
        """ Return interpolated nominal midnight humidity for date,
        based on summer and winter values.

        Parameters
        ----------
        date : datetime
            Date, including year, for which midnight humidity estimate is needed.
            UTC will be assumed if no timezone info is given.

        Returns
        -------
        humidity_for_date : float
            Midnight relative humidity estimate for date, in range [0, 100].
        """
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        humidity_for_date = self._interpolate_for_date(date,
                                                       self.summer_midnight_humidity,
                                                       self.winter_midnight_humidity)
        return humidity_for_date

    def extinction_for_date(self, date, filter):
        """ Return interpolated extinction for date, based on summer and winter values.

        Parameters
        ----------
        date : datetime
            Date, including year, for which atmospheric extinction estimate is wanted.
            UTC will be assumed if no timezone info is given.

        filter : str
            Filter for which extinction is wanted.

        Returns
        -------
        extinction_for_date : float
            Midnight atmospheric extinction estimate for date and filter,
            always positive.
        """
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        summer_extinction, winter_extinction = self.extinction[filter][0],\
                                               self.extinction[filter][1]
        extinction_for_date = self._interpolate_for_date(date, summer_extinction,
                                                         winter_extinction)
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
# Transforms = Filter Passband CI_pb1 CI_pb2 1st-order_tr [2nd-order tr] 
# One only per line.
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
    """Raised when transform line cannot be properly parsed."""
    pass


class Instrument:
    """Holds and makes available one telescope instrument's information.
    Includes information about mount performance, telescope, and camera.
    Loads info from one .ini file.

    Parameters
    ----------
    fullpath : str
        Full path to site's .ini file.

    Attributes
    ----------
    fullpath : str
        Full path to site's .ini file, from input parameter ``fullpath``.
    mount_model : str
        Name of telescope mount model.
    nominal_slew_time : float
        An estimate of the mean time required for mount to slew and settle, in seconds.
    ota_model : str
        Name of the Optical Tube Assembly (optical telescope proper) model.
    ota_aperture : float
        Diameter of the OTA's optical (sky-side) aperture, in meters.
    focal_length : float
        Instrument's focal length at the camera optical plane, in millimeters.

        This includes effects of any modifiers, such as focal reducers, within
        the optical path.
    camera_model : str
        Name of the camera's model.
    x_pixels : int
        Number of pixels in the camera's x-direction. X is taken as horizontal,
        as in the FITS convention.

        Usually this corresponds to Right Ascension in sky coordinates.
    y_pixels : int
        Number of pixels in the camera's y-direction. Y is taken as horizontal,
        as in the FITS convention.

        Usually this corresponds to Declination in sky coordinates.
    pixel_size : float
        Physical spacing between pixel centers, in microns.

        Assumes square pixels so that pixels sizes in the X and Y directions are equal.
    ccd_gain : float
        Camera chip gain in electrons per ADU. (Applies to CMOS sensors as well.)
    saturation_adu : float
        Highest ADU per pixel for which the user will trust his camera to remain linear.

        This need not equal full saturation ADU of the pixels. For non-anti-blooming
        scientific CCDs, this is often 70-90% of actual physical well depth.
    vignetting_pct_at_corner : float
        Estimate of the percent ADU drop at image corners relative to image center,
        for an evenly illuminated OTA front, as for panel flat images.

        Typically, values are in range (5, 40).
    nominal_cooling_time : float
        Estimate of the time required for camera cooling to complete and stabilize
        before imaging, in seconds.
    pinpoint_pixel_scale_multipler : float
        If images are plate-solved by the PinPoint solver, this is a correction to
        the calculated plate scale.

        PinPoint uses proprietary distortion constants that tend to make its WCS plate
        scale differ slightly from truly linear WCS plate scales from other solvers.
        This value is depends upon the optical details of a given optical train, e.g.,
        'pincushioning' and thus is typically constant for a given optical instrument.
        Determining the difference in plate scales between PinPoint and a
        linear-WCS solver will help find sky objects near image corners.
        Values will be very nearly, but not exactly one.
    filters_available : list of str
        Names of filters available in the instrument.
    v14_time_to_sn100 : dict (str: float)
        Estimate, by filter, of the time required for the instrument to achieve
        signal-to-noise ratio of 100 for a typical star of V magnitude 14.

        This value may be approximate, and is used when estimating appropriate
        extinction times. Keys are filter names, values are exposure times in seconds.
    transforms : dict (keys tuple of 4 str: values tuple of 1 or 2 float)
        Optical transform from filter data to standard passband, for example, from
        'V' filter data to 'V' standard passband.

        Keys are of form
        (f, pb, CI_pb1, CI_pb2), where f is the filter name, pb is the passband name,
        and CI_pb1 - CI_pb2 defines the light source's color index.
        Values are of the form (T1, ) where T1 is the transform value defined by the
        key, or rarely, (T1, T2) where T2 is the second-order transform value.
    min_fwhm_pixels : float
        Minimum full-width at half-maximum (FWHM) allowed for an image signal to be
        considered a legitimate point source (star, minor planet, etc.). In pixels.
    max_fwhm_pixels : float
        Maximum FWHM allowed for photometric processing, in pixels.
    nominal_fwhm_pixels : float
        Normal FWHM value expected for an image of satisfactory quality. In pixels.
    exposure_overhead : float
        Normal time expected from the end of one exposure to the beginning of the
        next exposure.

        Includes time for image download, plate-solving, FITS generation,
        saving to disk, and any CPU processing to initiate the next image. In seconds.
    max_exposure_no_guiding : float
        Maximum exposure length allowed without invoking guiding, in seconds.
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
        self.vignetting_pct_at_corner = float(i.get('Camera',
                                                    'Vignetting Pct At Corner'))
        self.nominal_cooling_time = float(i.get('Camera', 'Nominal Cooling Time'))
        self.pinpoint_pixel_scale_multipler = \
            float(i.get('Plate Solution', 'Pinpoint Pixel Scale Multiplier'))
        self.filters_available = \
            tuple(i.get('Filters', 'Available').replace(',', ' ').split())
        self.v14_time_to_sn100 = \
            _dict_to_floats(_parse_multiline(i.get('Filters',
                                                   'V14 Time To SN 100'), 2, 2))
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
        (filter, passband, color passband 1, color passband 2), and dictionary
        values are a tuple of float(s) giving transform value.
        """
        min_values_per_transform = 5
        max_values_per_transform = 6
        strings = multiline_string.splitlines()
        transforms_dict = {}
        for s in strings:
            values = s.replace(',', ' ').split()
            if len(values) < min_values_per_transform or \
                len(values) > max_values_per_transform:
                raise TransformParseError('>' + s + '<')
            key = tuple(values[:4])
            value = tuple([float(v) for v in values[4:]])
            transforms_dict[key] = value
        return transforms_dict


__________HUMANOBSERVER_INI_______________________________________________________ = 0

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
    """ Holds one human observer's information. Gets info from one .ini file.

    Parameters
    ----------
    fullpath : str
        Full path to site's .ini file.

    Attributes
    ----------
    fullpath : str
        Full path to site's .ini file, from input parameter ``fullpath``.
    name : str
        Observer's name.
    alcdef_contact_name : str
        Name as should be attached to data submitted to the ALCDEF (Asteroid
        Lightcurve Data Exchange Format database).
    alcdef_contact_info : str
        Typically the observers e-mail address.
    alcdef_observers : str
        Comma-separated list of observers' names.
        For one observer, may be the same as ``name`` above.
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
    """Raised when multi-line .ini field cannot be properly parsed."""
    pass


def get_ini_data(fullpath):
    """Read .ini file into ConfigParser object. Used at beginning of all
    astropak .ini file readers.

    Parameters
    ----------
    fullpath : str
        Full path to site's .ini file.

    Returns
    -------
    fullpath, filename, ini_data : tuple of str, str, dict
        Where:

        **fullpath** is the full path to site's .ini file,
        from input parameter ``fullpath``.

        **filename** is the base name of ``fullpath``.

        **ini_data** is a dict containing raw, unparsed data read from the
        .ini file at ``fullpath`'.
    """

    if not (os.path.exists(fullpath) and os.path.isfile(fullpath)):
        raise FileNotFoundError('fullpath >' + fullpath + '<')
    filename = os.path.basename(fullpath)
    ini_data = configparser.ConfigParser()
    ini_data.read(fullpath)
    return fullpath, filename, ini_data


def _parse_multiline(multiline_string,
                     min_words_per_line=None, max_words_per_line=None):
    """Parse one multiline ini value to a new dict, where each line's
    first substring (item) is made the key, and a tuple of remaining substrings
    (items) is made the value.

    Parameters
    ----------
    multiline_string : str
        The multiline string, extracted as an item from an ini file, typically by
        :fun:`.ini.get_ini_data()``. The multiline string to be parsed.

    min_words_per_line : int, optional
        The minimum number of items (white-space delimited substrings) per line
        that will be parsed per line, or None to allow one or more items per line.
        Default is None.

    max_words_per_line : int
        The maximum number of items (white-space delimited substrings) per line
        that will be parsed per line, or None to allow any number of items per line.
        Default is None.

    Returns
    -------
    result_dict: dict
        Results of parsing the multiline string.

        A python dictionary in which each
        key comprises the first white space delimited item of each line,
        and the value comprises a tuple of the remaining items of that line.
    """
    lines = multiline_string.splitlines()
    result_dict = {}
    for line in lines:
        words = line.replace(',', ' ').split()
        if min_words_per_line is not None:
            if len(words) < min_words_per_line:
                raise MultilineParseError('\'' + line + '\'')
        if max_words_per_line is not None:
            if len(words) > max_words_per_line:
                raise MultilineParseError('\'' + line + '\'')
        result_dict[words[0]] = tuple(words[1:])
    return result_dict


def _dict_to_floats(d):
    """Takes dict (e.g., a parsed ini value), and attempts to convert each item's
    value to a float or to a tuple of float.

    Change is made in-place to the original dict;
    if the original dict will be needed later, user should first make a (deep) copy.

    Parameters
    ----------
    d : dict
        Input dictionary with each item's value being a float, or a string
        representing a float, or a tuple of strings representing floats.

    Returns
    -------
    d : dict
        Modified dictionary with each item's value having been converted to float
        or to a tuple of float.
    """
    for k, v in d.items():
        if isinstance(d[k], tuple):
            d[k] = tuple([float(value) for value in v])
        else:
            d[k] = float(d[k])
    return d


def string_to_boolean(bool_string, default_value=None):
    """Attempt to extract Boolean meaning from string, return default_value if failed.

    Parameters
    ----------
    bool_string : str
        A string representing a boolean value. Interpretable values are:
        'true', 'yes, 'y', 'false', 'no', 'n', all case-insensitive.

    default_value : any object
        Value to be returned if ``bool_string`` is not interpreted as boolean.

    Returns
    -------
    bool_value : bool, or None
        Boolean value as interpreted from ``bool_string``, or None.
    """
    sl = bool_string.lower().strip()
    if sl in ['true', 'yes', 'y']:
        return True
    if sl in ['false', 'no', 'n']:
        return False
    return default_value
