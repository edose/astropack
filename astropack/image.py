""" Module astropack.image.
    Image, FITS, aperture handling.
    Adapted from author's photrix package, image module.
"""

__author__ = "Eric Dose :: Albuquerque"


# Python core packages:
from math import ceil, sqrt

# External packages:
import numpy as np
import astropy.io.fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
from astropy.time import Time, TimeDelta
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel
from photutils.morphology import data_properties
from photutils.segmentation import make_source_mask

# From elsewhere in this package:
from .reference import FWHM_PER_SIGMA, ARCSECONDS_PER_RADIAN
from .geometry import XY, DXY, Circle_in_2D, Rectangle_in_2D

FITS_EXTENSIONS = ['fts', 'fit', 'fits']  # allowed filename extensions
ISO_8601_FORMAT = '%Y-%m-%dT%H:%M:%S'


__all__ = ['FITS',
           'MaskError',
           'Ap',
           'PointSourceAp',
           'MovingSourceAp',
           'make_circular_mask',
           'make_pill_mask',
           'calc_background_value']


class MaskError(Exception):
    """Raised when encountering any mask error,
    most often when masks differ in shape from image cutout or from each other. """
    pass


__________CLASSES______________________________________________________________ = 0


class FITS:
    """Loads and makes available data from one FITS file.

    Intended to be immutable.
    Internally, all image data & coordinates are held as zero-based (y,x) arrays
    (python, first coordinate is y, second is x), and NOT as FITS which are (x,y),
    origin=1. But most users will want to use attribute ``.image_xy`` which puts x
    (horizontal) pixel axis as first index and y (vertical) axis as the second index.

    Parameters
    ----------
    fullpath : str
        Full path to the FITS file.

    pinpoint_pixel_scale_multiplier : float, optional
        Factor (prob. close to 1) by which to multiply pixel scale whenever
        PinPoint plate solution is detected. Corrects plate scale to linear WCS
        plate scale. Default is 1.

    Attributes
    ----------
    fullpath : str
        Full path to the FITS file, from input parameter ``fullpath``.

    pinpoint_pixel_scale_multiplier : float, optional
        Factor (prob. close to 1) by which to multiply pixel scale whenever
        PinPoint plate solution is detected. From input parameter
        ``pinpoint_pixel_scale_multiplier``.

    is_valid : bool
        True whenever all validation checks pass, else False.

    header : dict
        FITS header information.

    image_fits : 2-dimensional |ndarray| of float
        Image array (ADU per pixel) in FITS indexing convention.

    image_xy : 2-dimensional |ndarray| of float
        Image array (ADU per pixel) in MaxIm/Astrometric convention, indices are
        (horizontal x, vertical y) with pixel (0, 0) at top left. Transpose of
        ``.image_fits``.

    object : str
        Name of the image's target, as read from FITS header. E.g., 'NGC 3535'.

    focal_length : float
        Focal length, in millimeters, as read from FITS header, or as calculated
        from plate scale and pixel size, both as read from FITS header.

    exposure : float
        Exposure duration, in seconds, as read from FITS header.

    temperature : float
        CCD temperature, in degrees C, as read from FITS header.

    utc_start : |py.datetime|
        Image exposure start time, UTC, as read from FITS header.

    utc_mid : |py.datetime|
        Image exposure mid-time, UTC, as calculated from ``utc_start`` and ``exposure``.

    filter : str
        Filter through which exposure was made, as read from FITS header.

    airmass : float
        Airmass, as read from FITS header.

    guide_exposure : float
        Guider exposure duration, in seconds, as read from FITS header.

    fwhm : float
        Full width at half-maximum (FWHM) of stars detected in the image, in pixels,
        as read from FITS header.

    is_calibrated : bool
        True if image appears to be calibrated for flat and dark images,
        as determined from FITS header items. Currently, detects only those
        calibrations performed by MaxIm versions 5 or 6.

    wcs_fits : |WCS|
        World Coordinate System representation of a plate solution of the
        image, as read from FITS header. Linear solution only, no distortion
        parameters included.

    is_plate_solved : bool
        True if image appears to have been plate solved, as determined from FITS
        header items, else False.

    is_plate_solved_by_pinpoint : bool
        True if image appears to have been plate solved by PinPoint, as determined
        from presence or absence of PinPoint's proprietary distortion parameters,
        else False.

    wcs_corrected : |WCS|
        Linear WCS plate solution, as corrected by application of
        ``pinpoint_pixel_scale_multiplier``, or a copy of ``wcs_fits`` if
        ``pinpoint_pixel_scale_multiplier`` is absent or equals 1.

    ra_deg : float
        Right Ascension of image center, in degrees, from ``wcs_corrected``.

    dec_deg : float
        Declination of image center, in degrees, from ``wcs_corrected``.
    """

    """NB: ``pinpoint_pixel_scale_multiplier`` is a value (prob. near 1) by which to
    multiply pixel scale, iff PinPoint plate solution is detected.
    It's the best solution I can devise for the PinPoint's pixel scale deviating from
    those of linear WCS plate solvers.
    This arises because (sigh) Pinpoint plate solutions include "private"
    distortion parameters, so that their WCS values are not what they would be &
    should be for a linear-WCS-only solution.
    That is, PinPoint zero-order solution cannot be used as one, nor even correctable
    given its "private" distortion parameters.
    """

    def __init__(self, fullpath, pinpoint_pixel_scale_multiplier=1):
        self.fullpath = fullpath
        self.pinpoint_pixel_scale_multiplier = pinpoint_pixel_scale_multiplier
        try:
            hdulist = astropy.io.fits.open(self.fullpath)
        except IOError:
            self.is_valid = False
            return
        self.header = hdulist[0].header.copy()

        # FITS convention = (vert/Y, horiz/X), pixel (1,1) bottom left -- rarely used.
        self.image_fits = hdulist[0].data.astype(float)
        hdulist.close()

        # MaxIm/Astrometrica convention = (horiz/X, vert/Y) pixel (0,0) top left).
        # NB: self.image_fits, self.image_xy, and self.image are different views
        # of the SAME array. They are meant to be read-only--changing any one of
        # them *will* change the others.
        self.image_xy = np.transpose(self.image_fits)

        # Extract data from header:
        self.object = self.header_value('OBJECT')  # name of target, e.g., 'MP_845'.
        # self._is_calibrated = self._is_calibrated()  # True iff apparently calibrated.
        self.focal_length = self._get_focal_length()  # in mm.
        self.exposure = self.header_value(['EXPTIME', 'EXPOSURE'])  # in seconds
        self.temperature = self.header_value(['SET-TEMP', 'CCD-TEMP'])
        self.utc_start = self._get_utc_start()  # py datetime.
        self.utc_mid = self.utc_start + TimeDelta(self.exposure / 2.0, format='sec')
        self.filter = self.header_value('FILTER')
        self.airmass = self.header_value('AIRMASS')
        self.guide_exposure = self.header_value('TRAKTIME')  # in seconds
        self.fwhm = self.header_value('FWHM')
        self.is_calibrated = self._is_calibrated()

        # Make 2 WCS objects: (1) internal from FITS as-is,
        # (2) corrected pixel scale (normally used):
        self.wcs_fits = WCS(self.header)
        self.is_plate_solved = self._is_plate_solved()
        self.is_plate_solved_by_pinpoint = self._detect_pinpoint_plate_solution()
        self.wcs_corrected = self._make_corrected_wcs()

        # Determine image center (RA, Dec), whatever WCS chose as center pixels:
        xy_center = XY(self.image_xy.shape[0] / 2.0, self.image_xy.shape[1] / 2.0)
        skycoord_center = self.xy_to_skycoords(xy_center)
        self.ra_deg, self.dec_deg = skycoord_center.ra.degree,\
                                    skycoord_center.dec.degree

        self.is_valid = True  # if it got through all that initialization.

    def header_value(self, key):
        """Return value associated with given FITS header key, or None if key not found.

        Parameters
        ----------
        key : str, or list of str
            FITS header key, or list of keys to try.

        Returns
        -------
        value : str, or None
            If ``key`` is str, return value of the FITS header item under ``key``.
            If `key` is a list, return value of the FITS header item under the
            first valid key in ``key``.
            If ``key`` not found in header, return None.
        """
        if isinstance(key, str):
            return self.header.get(key, None)
        for k in key:
            value = self.header.get(k, None)
            if value is not None:
                return value
        return None

    def _detect_pinpoint_plate_solution(self):
        """True if proprietary PinPoint plate solution was detected, else False."""
        return all([key in self.header.keys()
                    for key in ['TR1_0', 'TR2_1', 'TR1_6', 'TR2_5']])

    def _make_corrected_wcs(self):
        """Make and return corrected WCS object.
        If `pinpoint_pixel_scale_multipler`` is one, return a copy of ``wcs_fits``."""
        corrected_header = self.header.copy()
        if self._detect_pinpoint_plate_solution():
            for key in (['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2']):
                if key in corrected_header:
                    corrected_header[key] = float(corrected_header[key]) *\
                                            self.pinpoint_pixel_scale_multiplier
        return WCS(corrected_header)

    def xy_to_skycoords(self, xy):
        """ Convert the image's (x,y) coordinates to (RA, Dec) sky coordinates.

        Wrapper for :meth:`astropy.wcs.WCS.pixel_to_world()`, using ``wcs_corrected``.

        Parameters
        ----------
        xy : |XY|, or tuple of 2 float, or a list of either
            Image (x,y) coordinates, or a list of such coordinates.
            If list, all elements must be of same type (all |XY|, or all tuple).

        Returns
        -------
        ra_dec : |SkyCoord|, scalar or array
            Sky coordinates corresponding to ``xy`` image coordinates, as derived
            using ``wcs_corrected``. |SkyCoord| is scalar type if ``xy`` is one tuple,
            or array type if ``xy`` is a list.
        """
        if isinstance(xy, XY):
            x, y = xy.x, xy.y
        elif isinstance(xy, tuple):
            x, y = xy
        elif isinstance(xy, list):
            if isinstance(xy[0], XY):
                xy_tuple_list = [(xy_element.x, xy_element.y) for xy_element in xy]
            elif isinstance(xy[0], tuple):
                xy_tuple_list = xy
            else:
                raise ValueError('FITS.xy_to_skycoords() requires a tuple'
                                 ' (x,y) or a list of such tuples.')
            x, y = tuple(zip(*xy_tuple_list))
        else:
            raise ValueError('FITS.xy_to_skycoords() requires a tuple'
                             ' (x,y) or a list of such tuples.')
        return pixel_to_skycoord(x, y, self.wcs_corrected, origin=0, mode='wcs')

    def skycoords_to_xy(self, skycoords):
        """Convert (RA, Dec) sky coordinates to image (x,y) coordinates.

        Wrapper for :meth:`astropy.wcs.WCS.world_to_pixel()`, using ``wcs_corrected``.

        Parameters
        ----------
        skycoords : |SkyCoord|, scalar or array
            RA, Dec sky coordinates.

        Returns
        -------
        xy : tuple of 2 float, or a list of such tuple
            Image (x,y) coordinates corresponding to sky coordinates, as
            derived using ``wcs_corrected``.
            ``xy`` is one tuple if ``skycoords`` is a |SkyCoord| of scalar type (one
            (RA, Dec) sky coordinate, or a list of tuples if ``skycoords``
            is array type.
        """
        # TODO: consider returning astropack.geometry.XY (namedtuple) or list of them.
        xy = skycoord_to_pixel(skycoords, self.wcs_corrected, origin=0, mode='wcs')
        if xy[0].size == 1:
            xy = float(xy[0]), float(xy[1])  # fix python weirdness with return types.
        else:
            xy = list(zip(xy[0], xy[1]))
        return xy

    def corner_skycoords(self):
        """Return sky coordinates for the four corners of this image.

        Returns
        -------
        skycoords : |SkyCoord|, array-type of length 4
            Sky coordinates for the 4 corners of ``image_xy``. Often used to make
            bounding box for catalog lookups
            (in combination with corner sky coordinates of other images).
        """
        x_size, y_size = self.image_xy.shape
        xy_list = [(0, 0), (x_size - 1, 0), (0, y_size - 1), (x_size - 1, y_size - 1)]
        skycoords = self.xy_to_skycoords(xy_list)
        return skycoords

    def _is_calibrated_by_maxim_5_6(self):
        hval = self.header_value('CALSTAT')
        if hval is not None:
            if hval.strip().upper() == 'BDF':  # calib. by MaxIm DL v. 5 or 6
                return True
        return False

    def _is_calibrated(self):
        """Returns True if this FITS file was photometrically calibrated
        (by Maxim 5 or 6), else False."""
        calib_fn_list = [self._is_calibrated_by_maxim_5_6()]
        return any([is_c for is_c in calib_fn_list])

    def _is_plate_solved(self):
        return self.wcs_fits is not None

    def _get_focal_length(self):
        # If FL available, return it. Else, compute FL from plate solution.
        value = self.header_value('FOCALLEN')
        if value is not None:
            return value  # mm
        x_pixel = self.header_value('XPIXSZ')
        y_pixel = self.header_value('YPIXSZ')
        x_scale = self.header_value('CDELT1')
        y_scale = self.header_value('CDELT2')
        if any([val is None for val in [x_pixel, y_pixel, x_scale, y_scale]]):
            return None
        fl_x = x_pixel / abs(x_scale) * (ARCSECONDS_PER_RADIAN / (3600 * 1800))
        fl_y = y_pixel / abs(y_scale) * (ARCSECONDS_PER_RADIAN / (3600 * 1800))
        return (fl_x + fl_y) / 2.0

    def _get_utc_start(self):
        utc_string = self.header_value('DATE-OBS')
        utc_start = Time(utc_string)
        utc_start.format = 'iso'
        return utc_start


class Ap:
    """PARENT CLASS of all apertures for aperture photometry.

     Not to be instantiated directly.
     Each aperture shape requires a specific subclass.

     Parameters
     ----------
     image_xy : 2-dimensional |ndarray| of float
        Image array, with (x,y) indexing.
        Most conveniently arranged via |FITS| instance and passing its ``image_xy``.

     xy_center : |XY| instance
        Pixel position (x,y) of light source within parent image.
        This should be the best prior estimate of the light source's centroid
        at mid-exposure.

     dxy_offset : |DXY| instance
        lowest (x,y) index of cutout (upper-left corner of image), that is,
        the offset of cutout origin from parent image's origin.

     foreground_mask : 2-dmensional |ndarray| of bool
        Mask array for pixels to be counted in flux, centroid, etc.
        Required. Uses (x,y) index convention and numpy mask convention
        (True = pixel is masked out and unused).
        Mask shape defines shape of cutout to be used.

     background_mask : 2-dimensional |ndarray| of bool, or int, or float, optional
        Specifies Mask array for pixels to be counted in background flux, centroid, etc.
        If |ndarray|: boolean array in same shape as foreground_mask, giving mask
        directly.
        If int or float: this mask will be set to that value per pixel. If zero,
        background mask is not used.
        If unspecified: this mask will be set to the inverse of `foreground_mask`.
        Numpy mask convention (True -> pixel masked out, unused).

    source_id : str, optional
        String identifying the source (e.g., comp star ID, MP number) that this
        aperture is intended to measure. Default is empty string.

    obs_id : str, optional
        String identifying the specific observation to which this aperture applies,
        typically unique among all observations (aperture instances) in one image.
        Default is empty string.

    Attributes
    ----------
    image_xy : |ndarray| of float
        Full image array with (x,y) indices, from input parameter ``image_xy``.

    xy_center : |XY| instance
        Pixel position (x,y) of light source within parent image, from input parameter
        ``xy_center``.

    dxy_offset : |DXY| instance
        lowest (x,y) index of cutout (upper-left corner of image), from input
        parameter ``xy_offset``.

    input_foreground_mask : |ndarray| of bool
        Mask array for pixels to be counted in flux, centroid, etc., from input
        parameter ``input_foreground_mask``.

    input_background_mask : |ndarray| of bool, or int, or float
        Mask array for pixels to be counted in background flux, centroid, etc.

    source_id : str
        String identifying the source, from input parameter ``source_id``.

    obs_id : str
        String identifying the specific observation to which this aperture applies,
        from input parameter ``obs_id``.

    is_valid : bool
        True if this aperture appears to be valid for use in photometry, else False.

    cutout : 2-dimensional |ndarray| of float
        A slice of the image that is located and is just sufficiently large to
        encompass the foreground and background masks.

    mask_overlap_pixel_count : int
        Number of pixels set to False (pixel used) in both foreground and
        background masks. Ideally zero when masks are defined explicitly as |ndarray|.

    foreground_mask : |ndarray| of float
        Mask array for pixels to be counted in flux, centroid, etc., from input
        parameter ``input_foreground_mask``.

    background_mask : |ndarray| of float
        Mask array for pixels to be counted in background flux, centroid, etc.,
        as given explicitly by or as derived from input parameter
        ``input_background_mask``. Same shape and image location as ``foreground_mask``.

    x_low, x_high : float
        Lowest and highest x pixel index defining location and shape of foreground mask,
        background mask, and image cutout.

    y_low, y_high : float
        Lowest and highest y pixel index defining location and shape of foreground mask,
        background mask, and image cutout.

    foreground_pixel_count : int
        Number of foreground mask pixels set to False, that is, number of image pixels
        used to compute raw flux.

    background_pixel_count : int
        Number of background mask pixels set to False, that is, number of image pixels
        used to compute background level.

    background_level : float
        Best estimate of the background level (ADU) for this aperture.

    background_std : float
        Best estimate of the per-pixel standard deviation of the background level
        (ADU) for this aperture.

    foreground_min : float
        Minimum ADU for any image pixel within the foreground aperture.

    foreground_max : float
        Minimum ADU for any image pixel within the foreground aperture.

    raw_flux : float
        Best estimate of the total ADU flux within the foreground aperture,
        without any background correction.

    net_flux : float
        Best estimate of the total ADU flux within the foreground aperture,
        corrected for background level.

    stats : `~photutils.segmentation.SourceCatalog
        Summary of statistics for `cutout`.

    xy_centroid : tuple of 2 float
        Best (current) estimate of the background-corrected centroid (x,y)
        for the light source in ``cutout``.

    sigma : float
        One-sigma width, in pixels, of the background-corrected flux in ``cutout``.

    fwhm : float
        Estimated full width at half-maximum, in pixels, of the background-corrected
        flux in ``cutout``.

    elongation : float
        Elongation of the background-corrected flux in ``cutout``.
    """

    def __init__(self, image_xy, xy_center, dxy_offset, foreground_mask,
                 background_mask=None, source_id='', obs_id=''):
        # Save inputs:
        self.image_xy = image_xy
        # Starting (x,y), relative to image (0,0):
        self.xy_center = xy_center
        # Cutout offset (dx,dy), relative to image (0,0).
        self.dxy_offset = dxy_offset
        input_foreground_mask = foreground_mask
        input_background_mask = background_mask
        self.source_id = str(source_id)
        self.obs_id = str(obs_id)

        # Default values:
        self.is_valid = None
        self.mask_overlap_pixel_count = None

        # Ensure background mask is boolean array of same shape as foreground mask:
        # If no background mask given, make one from inverse of foreground mask:
        if input_background_mask is None:
            self.background_mask = np.logical_not(input_foreground_mask)
        # If background mask is given but int or float array (rare),
        # convert it to a boolean numpy array:
        elif type(input_background_mask) in (int, float) and \
            input_background_mask == 0:
            self.background_mask = np.full_like(input_foreground_mask,
                                                fill_value=True, dtype=bool)
        # If background mask is valid array, use it directly
        elif isinstance(input_background_mask, np.ndarray):
            if input_background_mask.shape != input_foreground_mask.shape:
                raise MaskError('Foreground and background masks differ in shape.')
            self.background_mask = input_background_mask
        else:
            raise MaskError('Background mask type ' +
                            str(type(input_background_mask)) + ' is not valid.')

        # Determine x and y boundaries of cutout containing foreground and background:
        self.x_low = self.dxy_offset.dx
        self.y_low = self.dxy_offset.dy
        self.x_high = self.dxy_offset.dx + input_foreground_mask.shape[0] - 1
        self.y_high = self.dxy_offset.dy + input_foreground_mask.shape[1] - 1
        if not self._cutout_is_wholly_inside_image():
            self.is_valid = False
            return

        # Make the final cutout array and masks, all with indices in x,y order:
        self.cutout = image_xy[self.x_low: self.x_high + 1, self.y_low: self.y_high + 1]
        self.foreground_mask = input_foreground_mask.copy()  # indices x,y order.
        self.background_mask = self.background_mask.copy()        # "

        # Compute pixels in use by both masks (should always be zero):
        self.mask_overlap_pixel_count = \
            np.sum(np.logical_and((self.foreground_mask == False),
                                  (self.background_mask == False)))

        # Compute background mask statistics:
        self.foreground_pixel_count = np.sum(self.foreground_mask == False)
        self.background_pixel_count = np.sum(self.background_mask == False)
        if self.background_pixel_count >= 1:
            self.background_level, self.background_std = \
                calc_background_value(self.cutout, self.background_mask)
        else:
            self.background_level, self.background_std = 0.0, 0.0

        # Compute aperture statistics via numpy masked arrays (ma):
        foreground_ma = np.ma.array(data=self.cutout, mask=self.foreground_mask)
        self.foreground_max = np.ma.max(foreground_ma)
        self.foreground_min = np.ma.min(foreground_ma)
        self.raw_flux = np.ma.sum(foreground_ma)

        cutout_net_ma = np.ma.array(self.cutout - self.background_level,
                                    mask=self.foreground_mask)
        self.net_flux = np.sum(cutout_net_ma)
        self.stats = data_properties(
            data=cutout_net_ma.data,
            mask=self.foreground_mask,
            background=np.full_like(cutout_net_ma.data,
                                    fill_value=self.background_level))
        # NB: self.stats yields numpy (y,x) order.
        self.xy_centroid = (self.stats.ycentroid + self.dxy_offset.dx,
                            self.stats.xcentroid + self.dxy_offset.dy)

        # Sigma and FWHM represent spreading *other than* motion
        # (i.e., real optical dispersion).
        self.sigma = self.stats.semimajor_sigma.value
        self.fwhm = self.stats.fwhm.value  # given by photutils v.1.1.0
        self.elongation = self.stats.elongation.value
        self.is_valid = True

    def __str__(self):
        raise NotImplementedError

    def _cutout_is_wholly_inside_image(self):
        """Returns True iff cutout lies entirely within parent image."""
        return (self.x_low > 0) and (self.y_low > 0) and \
               (self.x_high < self.image_xy.shape[1]) and \
               (self.y_high < self.image_xy.shape[0])

    def flux_stddev(self, gain=1):
        """Returns estimated standard deviation of background-corrected flux.

        This was made a method so that ``gain`` can be passed in separately.

        Subclasses should inherit this method unchanged.

        Parameters
        ----------
        gain : float
            CCD-like gain in e-/ADU. Property of the specific camera (model).
            Needed only to improve accuracy of uncertainty estimations.

        Returns
        -------
        flux_stddev : float
            Standard deviation in ADU arising only from Poisson ('shot') noise.
        """
        # Poisson noise ~ flux in electrons.
        flux_variance_from_poisson_noise = self.raw_flux / gain
        flux_variance_from_background = self.foreground_pixel_count * \
            ((self.background_std ** 2) / self.background_pixel_count)
        flux_variance = flux_variance_from_poisson_noise + flux_variance_from_background
        flux_stddev = sqrt(flux_variance)
        return flux_stddev

    def _make_new_object(self, new_xy_center):
        """ Make new object from same image using new xy_center, with same mask shape.
            Each subclass of Ap must implement ``_make_new_object()``.
            Used mostly by ``recenter()``.
        """
        raise NotImplementedError("Every subclass of Ap must implement "
                                  "._make_new_object().")

    def recenter(self, max_adjustment=None, max_iterations=3):
        """Iteratively refine estimate of light source's center, and adjust the image
        cutout's location on the image if needed.

        Subclasses should inherit this method unchanged.

        Parameters
        ----------
        max_adjustment : float, optional
            The maximum adjustment of the center position that allows for
            convergence and exit from the refinement loop. If unspecified,
            refinement will continue until `max_iterations` refinement passes
            have been completed.

        max_iterations : int, optional
            The maximum number of refinement iterations to allow. Default is 3.

        Returns
        -------
        next_ap : |MovingSourceAp|
            A new `~.image.Ap` subclass instance with center closer to light
            source's centroid.
        """
        previous_ap, next_ap = self, self
        for i in range(max_iterations):
            previous_centroid = previous_ap.xy_centroid
            new_xy_center = previous_centroid
            next_ap = self._make_new_object(new_xy_center)
            new_centroid = next_ap.xy_centroid
            adjustment_distance = sqrt((new_centroid[0] - previous_centroid[0])**2 +
                                       (new_centroid[1] - previous_centroid[1])**2)
            if (max_adjustment is not None) and (adjustment_distance < max_adjustment):
                return next_ap
            previous_ap = next_ap
        return next_ap


class PointSourceAp(Ap):
    """Standard circular photometric aperture for stationary point source of light,
    esp. for a star.

    Subclass of `~.image.Ap`.
    Makes a circular foreground mask and an annular background mask, concentric,
    both centered on the given image coordinates of the point source.

    Parameters
    ----------

    image_xy : 2-dimensional |ndarray| of float
        Image array, with (x,y) indexing.
        Most conveniently arranged via |FITS| instance and passing its ``image_xy``.

    xy_center : |XY|, or tuple of 2 float
        Pixel position (x,y) of light source within parent image.
        This should be the best prior estimate of the light source's centroid
        at mid-exposure.

    foreground_radius : float, positive
        Radial size of foreground around point source centroid, in pixels.

    gap : float, non-negative
        Width of gap, that is, the difference between radius of foreground and
        inside radius of background annulus, in pixels.

    background_width : float, positive
        Radial width of annulus, that is, the difference between inside and
        outside radii of background annulus, in pixels.

    source_id : str, optional
        String identifying the source (e.g., comp star ID, MP number) that this
        aperture is intended to measure. Default is empty string.

    obs_id : str, optional
        String identifying the specific observation to which this aperture applies,
        typically unique among all observations (aperture instances) in one image.
        Default is empty string.

    Attributes
    ----------

    image_xy : 2-dimensional |ndarray| of float
        Image array, with (x,y) indexing.
        Most conveniently arranged via |FITS| instance and passing its ``image_xy``.

    xy_center : |XY| instance
        Pixel position (x,y) of light source within parent image.
        This should be the best prior estimate of the light source's centroid
        at mid-exposure.

    dxy_offset : |DXY| instance
        lowest (x,y) index of cutout (upper-left corner of image), that is,
        the offset of cutout origin from parent image's origin.

    foreground_mask : 2-dmensional |ndarray| of bool
        Mask array for pixels to be counted in flux, centroid, etc.
        Required. Uses (x,y) index convention and numpy mask convention
        (True = pixel is masked out and unused).
        Mask shape defines shape of cutout to be used.

    background_mask : 2-dimensional |ndarray| of bool, or int, or float, optional
        Specifies Mask array for pixels to be counted in background flux, centroid, etc.
        If |ndarray|: boolean array in same shape as foreground_mask, giving mask
        directly.
        If int or float: this mask will be set to that value per pixel. If zero,
        background mask is not used.
        If unspecified: this mask will be set to the inverse of `foreground_mask`.
        Numpy mask convention (True -> pixel masked out, unused).

    source_id : str, optional
        String identifying the source (e.g., comp star ID, MP number) that this
        aperture is intended to measure. Default is empty string.

    obs_id : str, optional
        String identifying the specific observation to which this aperture applies,
        typically unique among all observations (aperture instances) in one image.
        Default is empty string.

    foreground_radius : float
        Radial size of foreground around point source centroid, in pixels, from
        input parameter ``foreground_radius``.

    gap : float
        Width of gap, that is, the difference between radius of foreground and
        inside radius of background annulus, in pixels, from input parameter ``gap``.

    background_width : float
        Radial width of annulus, that is, the difference between inside and
        outside radii of background annulus, in pixels, from input parameter
        ``background_width``.

    annulus_inner_radius : float
        Radius of inner edge of background annulus, in pixels.

    annulus_outer_radius : float
        Radius of outer edge of background annulus, in pixels.

    is_valid : Bool
        True if aperture appears valid, else False.

    cutout : |ndarray|
        Small subset of image emcompassing the aperture and masks.

    mask_overlap_pixel_count : int
        Number of image pixels in both the foreground mask and background mask.
        Ideally zero.

    foreground_pixel_count : int
        Number of pixels for which foreground flux is summed.

    background_pixel_count : int
        Number of pixels for which background flux is summed.

    foreground_max : float
        Maximum ADU of any pixel in the foreground.

    background_min : float
        Maximum ADU of any pixel in the foreground.

    raw_flux : float
        Total flux in foreground, before subtracting background level.

    net_flux : float
        Background-corrected foreground flux in foreground.
        Used in computing the instrumental magnitude.

    stats : photutils :class:~photutils.segmentation.SourceCatalog instance
        A bundle of statistics derived from this aperture.

        See https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.
        SourceCatalog.html#photutils.segmentation.SourceCatalog

    xy_centroid : tuple of 2 float
        Pixel (x,y) position for centroid of background-corrected flux,
        for this target

    sigma : float
        RMS distance of pixels from centroid, weighted by background-corrected flux.
        A measure of flux spreading.

    fwhm : float
        Full-width at half-maximum of background-corrected flux. A measure of flux
        spreading, and estimated from ``sigma``.

    elongation : float
        A measure of non-circularity of background-corrected flux.
        Circular flux profiles will give zero.
        Specifically, elongation is the ratio of the lengths of semi-major and
        semi-minor axes of the corrected flux, considered as having best-fit
        elliptical shape.
    """

    def __init__(self, image_xy, xy_center, foreground_radius, gap, background_width,
                 source_id='', obs_id=''):
        if isinstance(xy_center, tuple):
            xy_center = XY(xy_center[0], xy_center[1])  # ensure is XY object.
        self.foreground_radius = foreground_radius
        self.gap = gap
        self.background_width = background_width
        self.annulus_inner_radius = self.foreground_radius + self.gap
        self.annulus_outer_radius = self.annulus_inner_radius + self.background_width
        cutout_size = int(ceil(2 * self.annulus_outer_radius)) + 4
        dxy_origin = DXY(int(round(xy_center.x - cutout_size / 2)),
                         int(round(xy_center.y - cutout_size / 2)))
        xy_center_in_cutout = xy_center - dxy_origin
        foreground_mask = make_circular_mask(mask_size=cutout_size,
                                             xy_origin=xy_center_in_cutout,
                                             radius=self.foreground_radius)
        background_mask_center_disc = \
            np.logical_not(make_circular_mask(cutout_size, xy_center_in_cutout,
                                              self.annulus_inner_radius))
        background_mask_outer_disc = make_circular_mask(cutout_size,
                                                        xy_center_in_cutout,
                                                        self.annulus_outer_radius)
        background_mask = np.logical_or(background_mask_center_disc,
                                        background_mask_outer_disc)
        super().__init__(image_xy, xy_center, dxy_origin,
                         foreground_mask, background_mask, source_id, obs_id)

    def _make_new_object(self, new_xy_center):
        """ Make new object using new xy_center. Overrides parent-class method, as
        required. Masks will be recreated by the constructor, using new xy_center.
        """
        return PointSourceAp(self.image_xy, new_xy_center,
                             self.foreground_radius, self.gap, self.background_width,
                             self.source_id, self.obs_id)

    def __str__(self):
        return 'PointSourceAp at x,y = ' + str(self.xy_center.x) + ', ' + \
               str(self.xy_center.y)


class MovingSourceAp(Ap):
    """ Elongated 'pill-shaped' photometric aperture for moving point source of light,
    esp. for a minor planet/asteroid.

    Subclass of `~.image.Ap`.
    Makes 'pill-shaped' foreground and background masks, concentric,
    both centered on the given image coordinates of the point source.

    A 'pill-shaped' mask comprises the union of all active pixels in three submasks:
    a circle mask centered at `xy_start`, a circle mask centered at `xy_end`, and
    all the pixels between the two circles, most conveniently derived as a rectangle
    with two opposite edges centered at ``xy_start`` and ``xy_end`` and having width
    of twice the circles' radius.

    Parameters
    ----------

    image_xy : 2-dimensional |ndarray| of float
        Full image array with (x,y) indices, from input parameter ``image_xy``.

    xy_start :|XY|, or tuple of 2 float
        Pixel position (x,y) of light source within parent image, at the beginning
        time of exposure.

    xy_end : |XY|, or tuple of 2 float
        Pixel position (x,y) of light source within parent image, at the end
        time of exposure.

    foreground_radius : float, positive
        Radial size of foreground around centroid, in pixels.

    gap : float, non-negative
        Width of gap, that is, the difference between radius of foreground and
        inside radius of background annulus, in pixels.

    background_width : float, positive
        Radial width of annulus, that is, the difference between inside and
        outside radii of background annulus, in pixels.

    source_id : str, optional
        String identifying the source (e.g., comp star ID, MP number) that this
        aperture is intended to measure. Default is empty string.

    obs_id= : str, optional
        String identifying the specific observation to which this aperture applies,
        typically unique among all observations (aperture instances) in one image.
        Default is empty string.

    Attributes
    ----------

    image_xy : 2-dimensional |ndarray| of float
        Image array, with (x,y) indexing.
        Most conveniently arranged via |FITS| instance and passing its ``image_xy``.

    xy_start :|XY| instance
        Pixel position (x,y) of light source within parent image, at the beginning
        time of exposure.

    xy_end : |XY| instance
        Pixel position (x,y) of light source within parent image, at the end
        time of exposure.

    dxy_offset : |DXY| instance
        lowest (x,y) index of cutout (upper-left corner of image), that is,
        the offset of cutout origin from parent image's origin.

    foreground_mask : 2-dmensional |ndarray| of bool
        Mask array for pixels to be counted in flux, centroid, etc.
        Required. Uses (x,y) index convention and numpy mask convention
        (True = pixel is masked out and unused).
        Mask shape defines shape of cutout to be used.

    background_mask : 2-dimensional |ndarray| of bool, or int, or float, optional
        Specifies Mask array for pixels to be counted in background flux, centroid, etc.
        If |ndarray|: boolean array in same shape as foreground_mask, giving mask
        directly.
        If int or float: this mask will be set to that value per pixel. If zero,
        background mask is not used.
        If unspecified: this mask will be set to the inverse of `foreground_mask`.
        Numpy mask convention (True -> pixel masked out, unused).

    source_id : str, optional
        String identifying the source (e.g., comp star ID, MP number) that this
        aperture is intended to measure. Default is empty string.

    obs_id : str, optional
        String identifying the specific observation to which this aperture applies,
        typically unique among all observations (aperture instances) in one image.
        Default is empty string.

    foreground_radius : float
        Radial size of foreground around point source centroid, in pixels, from
        input parameter ``foreground_radius``.

    gap : float
        Width of gap, that is, the difference between radius of foreground and
        inside radius of background annulus, in pixels, from input parameter ``gap``.

    background_width : float
        Radial width of annulus, that is, the difference between inside and
        outside radii of background annulus, in pixels, from input parameter
        ``background_width``.

    annulus_inner_radius : float
        Radius of inner edge of background annulus, in pixels.

    annulus_outer_radius : float
        Radius of outer edge of background annulus, in pixels.

    is_valid : Bool
        True if aperture appears valid, else False.

    cutout : |ndarray|
        Small subset of image emcompassing the aperture and masks.

    mask_overlap_pixel_count : int
        Number of image pixels in both the foreground mask and background mask.
        Ideally zero.

    foreground_pixel_count : int
        Number of pixels for which foreground flux is summed.

    background_pixel_count : int
        Number of pixels for which background flux is summed.

    foreground_max : float
        Maximum ADU of any pixel in the foreground.

    background_min : float
        Maximum ADU of any pixel in the foreground.

    raw_flux : float
        Total flux in foreground, before subtracting background level.

    net_flux : float
        Background-corrected foreground flux in foreground.
        Used in computing the instrumental magnitude.

    stats : photutils :class:~photutils.segmentation.SourceCatalog instance
        A bundle of statistics derived from this aperture.

        See https://photutils.readthedocs.io/en/stable/api/photutils.segmentation.
        SourceCatalog.html#photutils.segmentation.SourceCatalog

    xy_centroid : tuple of 2 float
        Pixel (x,y) position for centroid of background-corrected flux,
        for this target

    sigma : float
        RMS distance of pixels from centroid, weighted by background-corrected flux.
        A measure of flux spreading.

    fwhm : float
        Full-width at half-maximum of background-corrected flux. A measure of flux
        spreading, and estimated from ``sigma``.

    elongation : float
        A measure of non-circularity of background-corrected flux.
        Circular flux profiles will give zero.
        Specifically, elongation is the ratio of the lengths of semi-major and
        semi-minor axes of the corrected flux, considered as having best-fit
        elliptical shape.        Radius of outer edge of background annulus, in pixels.
    """

    def __init__(self, image_xy, xy_start, xy_end,
                 foreground_radius, gap, background_width,
                 source_id='', obs_id=''):
        self.xy_start = xy_start if isinstance(xy_start, XY) \
            else XY.from_tuple(xy_start)
        self.xy_end = xy_end if isinstance(xy_end, XY) \
            else XY.from_tuple(xy_end)
        self.foreground_radius = foreground_radius
        self.gap = gap
        self.background_width = background_width
        self.background_inner_radius = self.foreground_radius + self.gap
        self.background_outer_radius = self.foreground_radius + \
                                       self.gap + self.background_width
        xy_center = self.xy_start + (self.xy_end - self.xy_start) / 2
        corner_x_values = (self.xy_start.x + self.background_outer_radius,
                           self.xy_start.x - self.background_outer_radius,
                           self.xy_end.x + self.background_outer_radius,
                           self.xy_end.x - self.background_outer_radius)
        x_min = min(corner_x_values)
        x_max = max(corner_x_values)
        corner_y_values = (self.xy_start.y + self.background_outer_radius,
                           self.xy_start.y - self.background_outer_radius,
                           self.xy_end.y + self.background_outer_radius,
                           self.xy_end.y - self.background_outer_radius)
        y_min = min(corner_y_values)
        y_max = max(corner_y_values)
        dxy_cutout_size = DXY(int(round(x_max - x_min + 4)),
                              int(round(y_max - y_min + 4)))
        dxy_offset = DXY(int(round(xy_center.x) - dxy_cutout_size.dx / 2.0),
                         int(round(xy_center.y) - dxy_cutout_size.dy / 2.0))
        xy_start_cutout = self.xy_start - dxy_offset
        xy_end_cutout = self.xy_end - dxy_offset
        foreground_mask = make_pill_mask(dxy_cutout_size,
                                         xy_start_cutout, xy_end_cutout,
                                         self.foreground_radius)
        background_inner_mask = make_pill_mask(dxy_cutout_size,
                                               xy_start_cutout, xy_end_cutout,
                                               self.background_inner_radius)
        background_outer_mask = make_pill_mask(dxy_cutout_size,
                                               xy_start_cutout, xy_end_cutout,
                                               self.background_outer_radius)
        background_mask = np.logical_or(background_outer_mask,
                                        np.logical_not(background_inner_mask))
        super().__init__(image_xy, xy_center, dxy_offset,
                         foreground_mask, background_mask,
                         source_id, obs_id)
        self.sigma = self.stats.semiminor_sigma.value
        self.fwhm = self.sigma * FWHM_PER_SIGMA

    def _make_new_object(self, new_xy_center):
        """ Make new object using new xy_center. Overrides parent-class method,
        as required.
        Masks will be recreated by the constructor, using new xy_center.
        """
        if isinstance(new_xy_center, tuple):
            new_xy_center = XY.from_tuple(new_xy_center)
        current_xy_center = self.xy_start + (self.xy_end - self.xy_start) / 2
        dxy_shift = new_xy_center - current_xy_center
        new_xy_start = self.xy_start + dxy_shift
        new_xy_end = self.xy_end + dxy_shift
        return MovingSourceAp(self.image_xy, new_xy_start, new_xy_end,
                              self.foreground_radius, self.gap, self.background_width,
                              self.source_id, self.obs_id)

    def __str__(self):
        return 'MovingSourceAp at x,y = ' + str(self.xy_center.x) + ', ' + \
               str(self.xy_center.y)


_____IMAGE_and_GEOMETRY_SUPPORT____________________________________ = 0


def make_circular_mask(mask_size, xy_origin, radius):
    """Construct a traditional circular mask array for small, stationary object,
    esp. for a star.
    Unmask only those pixels *within* ``radius`` pixels of a given point.

    Parameters
    ----------
    mask_size : int
        Edge size of new mask array, which will be square.

    xy_origin : |XY|, tuple of 2 float
        Pixel (x, y) coordinates of circle's origin, relative to mask's (0, 0) origin.

    radius : float
        Radius of circle, that is, of unmasked pixels within mask.

    Returns
    -------
    mask : 2-dimensional |ndarray| of bool
        Circular mask array in which True denotes 'masked out' and not
        used in calculations.
    """
    if isinstance(xy_origin, tuple):
        xy_origin = XY.from_tuple(xy_origin)  # ensure is XY object.
    circle = Circle_in_2D(xy_origin=xy_origin, radius=radius)
    is_inside = circle.contains_points_unitgrid(0, mask_size - 1, 0, mask_size - 1,
                                                include_edges=True)
    mask = np.transpose(np.logical_not(is_inside))
    return mask


def make_pill_mask(mask_shape_xy, xya, xyb, radius):
    """Construct a mask array for light source in motion (e.g., minor planet)
     Unmask only those pixels within `radius` pixels of the line segment
     from ``xya`` to ``xyb``.

    Parameters
    ----------
    mask_shape_xy : tuple of 2 float
        Pixel size (x,y) of mask array to generate.

    xya : |XY|, or tuple of 2 float
        Pixel (xa, ya) coordinates of light source centroid at beginning of motion.

    xyb : |XY|, tuple of 2 float
        Pixel (xa, ya) coordinates of light source centroid at end of motion.

    radius : float
        Radius, in pixels,  of end arcs and half-width of center region.

    Returns
    -------
    mask : 2-dimensional |ndarray| of bool
        'Pill-shaped' mask array in which True denotes 'masked out'
        and not used in calculations.
    """
    if isinstance(xya, tuple):
        xya = XY.from_tuple(xya)  # ensure is XY object.
    if isinstance(xyb, tuple):
        xyb = XY.from_tuple(xyb)  # ensure is XY object.
    if xya == xyb:
        return make_circular_mask(max(mask_shape_xy), xya, radius)

    # Make circle and rectangle objects:
    circle_a = Circle_in_2D(xya, radius)
    circle_b = Circle_in_2D(xyb, radius)
    dxy_ab = xyb - xya
    length_ab = dxy_ab.length
    dxy_a_corner1 = (radius / length_ab) * DXY(dxy_ab.dy, -dxy_ab.dx)
    dxy_a_corner2 = (radius / length_ab) * DXY(-dxy_ab.dy, dxy_ab.dx)
    xy_corner1 = xya + dxy_a_corner1
    xy_corner2 = xya + dxy_a_corner2
    xy_corner3 = xyb + dxy_a_corner2
    rectangle = Rectangle_in_2D(xy_corner1, xy_corner2, xy_corner3)

    # Make mask, including edges so no gaps can appear at rectangle corners:
    circle_a_contains = \
        circle_a.contains_points_unitgrid(0, mask_shape_xy.as_tuple[0] - 1,
                                          0, mask_shape_xy.as_tuple[1] - 1, True)
    circle_b_contains = \
        circle_b.contains_points_unitgrid(0, mask_shape_xy.as_tuple[0] - 1,
                                          0, mask_shape_xy.as_tuple[1] - 1, True)
    rectangle_contains = \
        rectangle.contains_points_unitgrid(0, mask_shape_xy.as_tuple[0] - 1,
                                           0, mask_shape_xy.as_tuple[1] - 1, True)
    # Render each in numpy mask-boolean but (x,y) index conventions:
    circle_a_mask = np.logical_not(circle_a_contains)
    circle_b_mask = np.logical_not(circle_b_contains)
    rectangle_mask = np.logical_not(rectangle_contains)
    mask = np.logical_and(np.logical_and(circle_a_mask, circle_b_mask), rectangle_mask)
    return mask


def calc_background_value(data, mask=None, dilate_size=3):
    """Calculate the best estimate of background value.

    Iterative sigma-clipped median is the best known algorithm. Author's testing
    with synthetic 16-bit image data frequently gives accuracy to 1 ADU.

    Parameters
    ----------
    data : 2-dimensional |ndarray| of float
        Array of pixels to be used (subject to masking) in calculating
        the background value. Typically, an image cutout array.

    mask : |ndarray| of bool, or None
        If |ndarray|, the mask array defining which pixels in ``data`` are used in
        calculating the background value.
        If None (default, but **not recommended**), use all the pixels in ``data``.

    dilate_size : float
        For a detected outlier pixel that is masked out, all pixels within a distance
        of ``dilate_size`` pixels will also be masked out.

        If image is undersampled (e.g., FWHM < 3 pixels), ``dilate_size``
        should probably be set to 1-2 pixels.
        If image is oversampled (e.g., FWHM > 10-15 pixels), ``dilate_size`` might
        best be set to 25-50% of FWHM.

    Returns
    -------
    median, std : tuple of 2 float
        Best estimate of background level (flux per pixel) and
        standard deviation, both within active background pixels.
    """
    if mask is None:  # use all pixels.
        this_mask = np.full_like(data, False, dtype=bool)
    elif mask.shape != data.shape:  # bad mask shape.
        return None
    # nb: despite python warnings, to write ...(mask is False) would fail here.
    # Must be mask == False.
    # elif np.sum(mask is False) == 0:  # <-- i.e., this will not work.
    elif np.sum(mask == False) == 0:  # no valid pixels.
        return 0.0, 0.0
    else:
        this_mask = mask.copy()
    # Notes on make_source_mask():
    #    this_mask: user-supplied mask of pixels simply not to consider at all.
    #    npixels = 2, which will mask out practically all cosmic ray hits.
    #    kernel(fwhm=2), a bit of smoothing,
    #    dilate_size=3, small enough to preserve most background pixels,
    #    but user may override.
    #    source_mask: masks out detected light-source pixels.
    #    MAY INCLUDE pixels masked out by this_mask!
    #    stats_mask: the intersection of valid pixels from this_mask and source_mask.
    kernel = Gaussian2DKernel(2.0 / FWHM_PER_SIGMA)
    source_mask = make_source_mask(data, mask=this_mask, nsigma=2, npixels=5,
                                   kernel=kernel, dilate_size=int(dilate_size))
    stats_mask = np.logical_or(this_mask, source_mask)
    _, median, std = sigma_clipped_stats(data, sigma=3.0, mask=stats_mask)
    return median, std


# _____FITS_FILE_HANDLING_______________________________________________ = 0
#
# TODO: write all_fits_files().
