""" Module astropak.image.
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
from photutils.morphology import data_properties
from photutils.segmentation import make_source_mask

# From elsewhere in this package:
from .reference import FWHM_PER_SIGMA, ARCSECONDS_PER_RADIAN
from .geometry import XY, DXY, Circle_in_2D, Rectangle_in_2D

FITS_EXTENSIONS = ['fts', 'fit', 'fits']  # allowed filename extensions
ISO_8601_FORMAT = '%Y-%m-%dT%H:%M:%S'


class MaskError(Exception):
    """ Raised on any mask error, usually when masks differ in shape. """
    pass


__________CLASSES______________________________________________________________ = 0


class FITS:
    """ Holds data from a FITS file. Immutable. Used mostly by an Image object (class Image).
        Internally, ALL image data & coordinates are held as zero-based (y,x) arrays (python, first
            coordinate is y, second is x), and NOT as FITS which are (x,y), origin=1
        TESTS OK 2022-02-26.
    Usage: obj = FITS(fullpath='C:/Astro/Borea/xxx.fts')
    """
    def __init__(self, fullpath, pinpoint_pixel_scale_multiplier=1):
        """
        :param fullpath:
        :param pinpoint_pixel_scale_multiplier: value (prob. near 1) by which to multiply pixel scale
            IFF pinpoint plate solution is detected. The best solution I can devise for the PinPoint mess.
            Required because (sigh) Pinpoint plate solutions include "private" distortion parameters
            so that their WCS values are not what they would be & should be for a WCS-only solution.
            That is, the zero-order solution is *called* a proper WCS but IS NOT, and it cannot be
            used as one, nor even correctable given its "private" (user-hostile) distortion parameters.
        """
        self.fullpath = fullpath
        self.pinpoint_pixel_scale_multiplier = pinpoint_pixel_scale_multiplier
        try:
            hdulist = astropy.io.fits.open(self.fullpath)
        except IOError:
            self.is_valid = False
            return
        self.header = hdulist[0].header.copy()

        # FITS convention = (vert/Y, horiz/X), pixel (1,1) at bottom left -- typically never used.
        self.image_fits = hdulist[0].data.astype(np.float64)
        hdulist.close()

        # MaxIm/Astrometrica convention = (horiz/X, vert/Y) pixel (0,0 at top left). USE THIS.
        # NB: self.image_fits, self.image_xy, and self.image are different views of the SAME array.
        #     They are meant to be read-only--changing any one of them *will* change the others.
        self.image_xy = np.transpose(self.image_fits)  # x and y axes as expected (not like FITS).

        # Extract data from header:
        self.object = self.header_value('OBJECT')  # name of target, e.g., 'MP_845'.
        # self._is_calibrated = self._is_calibrated()  # True iff apparently calibrated.
        self.focal_length = self._get_focal_length()  # in mm.
        self.exposure = self.header_value(['EXPTIME', 'EXPOSURE'])  # in seconds
        self.temperature = self.header_value(['SET-TEMP', 'CCD-TEMP'])  # deg C --> rename?
        self.utc_start = self._get_utc_start()  # py datetime.
        # self.utc_mid = self.utc_start + timedelta(seconds=self.exposure / 2.0)  # py datetime
        self.utc_mid = self.utc_start + TimeDelta(self.exposure / 2.0, format='sec')
        self.filter = self.header_value('FILTER')
        self.airmass = self.header_value('AIRMASS')
        self.guide_exposure = self.header_value('TRAKTIME')  # in seconds
        self.fwhm = self.header_value('FWHM')  # in pixels; independent of plate solution.
        self.is_calibrated = self._is_calibrated()

        # Make 2 WCS objects: (1) internal from FITS as-is, (2) corrected pixel scale (normally used):
        self.wcs_fits = WCS(self.header)
        self.is_plate_solved = self._is_plate_solved()
        self.is_plate_solved_by_pinpoint = self._detect_pinpoint_plate_solution()
        self.wcs_corrected = self._make_corrected_wcs()

        # Determine image center (RA, Dec), whatever WCS chose as center pixels:
        xy_center = XY(self.image_xy.shape[0] / 2.0, self.image_xy.shape[1] / 2.0)
        skycoord_center = self.xy_to_skycoords(xy_center)
        self.ra_deg, self.dec_deg = skycoord_center.ra.degree, skycoord_center.dec.degree

        self.is_valid = True  # if it got through all that initialization.

    def header_value(self, key):
        """Return value associated with given FITS header key.
        :param key: FITS header key [string] or list of keys to try [list of strings]
        :return: value of FITS header entry, typically [float] if possible, else [string].
            None if key not found in header.
        """
        if isinstance(key, str):
            return self.header.get(key, None)
        for k in key:
            value = self.header.get(k, None)
            if value is not None:
                return value
        return None

    def _detect_pinpoint_plate_solution(self):
        """Detect and return boolean, True iff proprietary PinPoint plate solution was detected."""
        return all([key in self.header.keys()
                    for key in ['TR1_0', 'TR2_1', 'TR1_6', 'TR2_5']])

    def _make_corrected_wcs(self):
        """Make and return corrected WCS object. If pixel_scale_multipler is one, corrected WCS
            object will be identical to self.wcs_fits."""
        corrected_header = self.header.copy()
        if self._detect_pinpoint_plate_solution():
            for key in (['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2']):
                if key in corrected_header:
                    corrected_header[key] = float(corrected_header[key]) *\
                                            self.pinpoint_pixel_scale_multiplier
        return WCS(corrected_header)

    def xy_to_skycoords(self, xy):
        """Wrapper for astropy.wcs.pixel_to_world(), using corrected plate solution.
           Accept a scalar tuple (x,y) or a list of such tuples,
           return one corresponding SkyCoord object (which is scalar or array).
        :param xy: scalar (x,y) tuple, or list of (x,y) tuples.
        :return: astropy SkyCoord object (scalar or array type, matching parameter pixels).
        """
        if isinstance(xy, tuple):
            x, y = xy
        elif isinstance(xy, list):
            x, y = tuple(zip(*xy))
        else:
            raise ValueError('FITS.xy_to_skycoords() requires a tuple (x,y) or a list of such tuples.')
        return pixel_to_skycoord(x, y, self.wcs_corrected, origin=0, mode='wcs')

    def skycoords_to_xy(self, skycoords):
        """Wrapper for astropy.wcs.world_to_pixel(), using corrected plate solution.
           Accept a SkyCoord object (scalar or array),
           return a scalar (x,y) tuple for a single coordinate, else a list of (x,y) tuples.
        :param skycoords: astropy Skycoord object (scalar or array type).
        :return: pixel position as scalar (x,y) tuple, or list of (x,y) tuples.
        """
        # TODO: consider returning astropak.geometry.XY (which is a namedtuple) or list of them.
        xy = skycoord_to_pixel(skycoords, self.wcs_corrected, origin=0, mode='wcs')
        if xy[0].size == 1:
            xy = float(xy[0]), float(xy[1])  # due to python weirdness.
        else:
            xy = list(zip(xy[0], xy[1]))
        return xy

    def corner_skycoords(self,):
        """Return sky coordinates for 4 corners of this image.
           Used to make bounding box (of multiple images, probably).
        :return: RA,Dec values for image corners. [SkyCoord object of length 4]
        """
        x_size, y_size = self.image_xy.shape
        xy_list = [(0, 0), (x_size - 1, 0), (0, y_size - 1), (x_size - 1, y_size - 1)]
        sc = self.xy_to_skycoords(xy_list)
        return sc

    def _is_calibrated_by_maxim_5_6(self):
        hval = self.header_value('CALSTAT')
        if hval is not None:
            if hval.strip().upper() == 'BDF':  # calib. by MaxIm DL v. 5 or 6
                return True
        return False

    def _is_calibrated(self):
        """Returns True iff this FITS file was photometrically calibrated (by Maxim 5 or 6)."""
        calib_fn_list = [self._is_calibrated_by_maxim_5_6()]  # may add more fns when available.
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
    """ PARENT CLASS of all apertures for aperture photometry, each shape being a specific subclass.
        Defines a slice of a full image based on parameters passed in, computes and holds properties.
        Ap-class objects are IMMUTABLE once constructed; re-centering generates and returns a new object.
        Masks are boolean arrays, must match data in shape, and follow numpy's mask convention:
            mask pixel=True means corresponding data pixel is invalid ("masked away") and not used.
            mask pixel=False means data pixel is valid and used.
        Masks are foreground (light source) and background (surrounding sky).
        [New in astropak2022] Cutout (=bounding box of background and foreground masks) must lie entirely
            inside the parent image, else this object is invalidated.
    """
    def __init__(self, image_xy, xy_center, xy_offset, foreground_mask, background_mask=None,
                 source_id='', obs_id=''):
        """ General constructor, from explicitly passed-in parent data array and 2 mask arrays.
        :param image_xy: the parent image array with indices in x,y order (transpose of FITS order)
                   [numpy ndarray; to pass in CCDData or numpy masked array,
                   please see separate, specific constructors, below].
        :param xy_center: center pixel position in parent. This should be the best prior estimate
                   of the light source's centroid at mid-exposure, as (x,y) (not as numpy [y, x] array).
                   [XY object, 2-tuple, 2-list, or 2-array of floats]
        :param xy_offset: lowest index of cutout (upper-left corner of image), that is, the offset of
                   cutout indices from parent image's indices.
                   [XY object, 2-tuple, 2-list, or 2-array of floats]
        :param foreground_mask: mask array for pixels to be counted in flux, centroid, etc.
                   Required boolean array. Numpy mask convention (True -> pixel masked out, unused).
                   Mask shape defines shape of cutout to be used.
                   Indices must be in x,y order. [numpy ndarray of booleans]
        :param background_mask: mask array for pixels to be counted in background flux, centroid, etc.
                   If None, will be set to inverse of foreground_mask.
                   If zero [int or float], will be set to zero (that is, background is not used).
                   Otherwise (normal case), must be boolean array in same shape as foreground_mask.
                   Indices must be in x,y order.
                   Numpy mask convention (True -> pixel masked out, unused).
        :param source_id: optional string describing the source (e.g., comp star ID, MP number) that this
                   aperture is intended to measure. [string]
        :param obs_id: optional observation ID string, will typically be unique among all observations
                   (aperture objects) in one image. [string]
        """
        # Save inputs:
        self.image_xy = image_xy
        self.xy_center = XY.from_tuple(xy_center)  # starting x,y, relative to image (0,0).
        self.xy_offset = XY.from_tuple(xy_offset)  # cutout offset x,y, relative to image (0,0).
        self.input_foreground_mask = foreground_mask
        self.input_background_mask = background_mask
        self.source_id = str(source_id)
        self.obs_id = str(obs_id)

        # Default values:
        self.is_valid = None
        self.mask_overlap_pixel_count = None

        # Ensure background mask is boolean array of same shape as foreground mask, whatever was input:
        # If no background mask given, make one from inverse of foreground mask:
        if self.input_background_mask is None:
            self.background_mask = np.logical_not(self.input_foreground_mask)
        # If background mask is given but int or float array (rare), convert it to a boolean numpy array:
        elif type(self.input_background_mask) in (int, float) and self.input_background_mask == 0:
            self.background_mask = np.full_like(self.input_foreground_mask, fill_value=True, dtype=np.bool)
        # If background mask is valid array, use it directly
        elif isinstance(self.input_background_mask, np.ndarray):
            if self.input_background_mask.shape != self.input_foreground_mask.shape:
                raise MaskError('Foreground and background masks differ in shape.')
            self.background_mask = self.input_background_mask
        else:
            raise MaskError('Background mask type ' + str(type(self.input_background_mask)) +
                            ' is not valid.')

        # Determine x and y boundaries of cutout containing foreground and background:
        self.x_low = self.xy_offset.x
        self.y_low = self.xy_offset.y
        self.x_high = self.xy_offset.x + self.input_foreground_mask.shape[0] - 1
        self.y_high = self.xy_offset.y + self.input_foreground_mask.shape[1] - 1
        if not self._cutout_is_wholly_inside_image():
            self.is_valid = False
            return

        # Make the final cutout array and masks, all with indices in x,y order:
        self.cutout = image_xy[self.x_low: self.x_high + 1, self.y_low: self.y_high + 1]
        self.foreground_mask = self.input_foreground_mask.copy()  # indices in x,y order.
        self.background_mask = self.background_mask.copy()        # "

        # Compute pixels in use by both masks (should always be zero):
        self.mask_overlap_pixel_count = np.sum(np.logical_and((self.foreground_mask == False),
                                                              (self.background_mask == False)))

        # Compute background mask statistics:
        self.foreground_pixel_count = np.sum(self.foreground_mask == False)
        self.background_pixel_count = np.sum(self.background_mask == False)
        if self.background_pixel_count >= 1:
            self.background_level, self.background_std = calc_background_value(self.cutout,
                                                                               self.background_mask)
        else:
            self.background_level, self.background_std = 0.0, 0.0

        # Compute aperture statistics via numpy masked arrays (ma):
        foreground_ma = np.ma.array(data=self.cutout, mask=self.foreground_mask)
        self.foreground_max = np.ma.max(foreground_ma)
        self.foreground_min = np.ma.min(foreground_ma)
        self.raw_flux = np.ma.sum(foreground_ma)

        cutout_net_ma = np.ma.array(self.cutout - self.background_level, mask=self.foreground_mask)
        self.net_flux = np.sum(cutout_net_ma)
        self.stats = data_properties(data=cutout_net_ma.data,
                                     mask=self.foreground_mask,
                                     background=np.full_like(cutout_net_ma.data,
                                                             fill_value=self.background_level))
        # NB: self.stats yields numpy (y,x) order.
        self.xy_centroid = (self.stats.ycentroid + self.xy_offset.x,
                            self.stats.xcentroid + self.xy_offset.y)

        # Sigma and FWHM represent spreading *other than* motion (i.e., real optical dispersion).
        self.sigma = self.stats.semimajor_sigma.value
        self.fwhm = self.stats.fwhm.value  # given by photutils v.1.1.0
        self.elongation = self.stats.elongation.value
        self.is_valid = True

    def __str__(self):
        raise NotImplementedError

    def _cutout_is_wholly_inside_image(self):
        """Returns True iff cutout lies entirely within parent image."""
        return (self.x_low > 0) and (self.y_low > 0) and \
               (self.x_high < self.image_xy.shape[1]) and (self.y_high < self.image_xy.shape[0])

    def flux_stddev(self, gain=1):
        """ Returns tuple of stats related to net flux of foreground pixels.
            Made a property so that gain can be passed in separately from Ap construction.
        :param gain: CCD-like gain in e-/ADU. Property of the specific camera (model).
                     Needed only for accurate uncertainty estimation. [float]
        :return: flux standard deviation. [float]
        """
        flux_variance_from_poisson_noise = self.raw_flux / gain  # from var(e-) = flux in e-.
        flux_variance_from_background = self.foreground_pixel_count * \
            ((self.background_std ** 2) / self.background_pixel_count)
        flux_variance = flux_variance_from_poisson_noise + flux_variance_from_background
        flux_stddev = sqrt(flux_variance)
        return flux_stddev

    def _make_new_object(self, new_xy_center):
        """ Make new object from same image using new xy_center, with same mask shape.
            Each subclass of Ap must implement ._make_new_object(). Used mostly by .recenter().
        """
        raise NotImplementedError("Every subclass of Ap must implement ._make_new_object().")

    def recenter(self, max_adjustment=None, max_iterations=3):
        """ Subclasses should inherit this unchanged.
            (By contrast, ._make_new_object() must be written specially for every subclass).
        :param max_adjustment: max distance, in pixels, that allows for convergence (return immediately
             without further adjustment iterations. [float]
        :param max_iterations: absolute max number of adjustment cycles to permit. [int]
        :return: newly recentered object. [Ap subclass object]
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
    """ Standard photometric aperture for stationary point source of light, esp. for a star.
            Always makes a circular foreground mask and an annular background mask, both centered on the
            given image coordinates of the point source.
        (If we will need a background mask that bleeds to the cutout's edges rather than being restricted
            to an annulus--for example, to work close to the parent image's edges, then that will
            definitely require a new sibling class to this one so that recentering retains the mask shapes.)
    """
    def __init__(self, image_xy, xy_center, foreground_radius, gap, background_width,
                 source_id=None, obs_id=None):
        """ Main and probably sole constructor for this class.
        :param image_xy: the parent image array, indices in x,y order (transpose of FITS image)
               [numpy ndarray].
        :param xy_center: center pixel position in parent. This should be the best prior estimate
               of the light source's centroid at mid-exposure, as (x,y) (not as numpy [y, x] array).
               [XY object, 2-tuple, 2-list, or 2-array of floats]
        :param foreground_radius: radial size of foreground around point source, in pixels. [float]
        :param gap: width of gap, difference between radius of foreground and inside radius of
               background annulus, in pixels. [float]
        :param background_width: width of annulus, difference between inside and outside radii
               of background annulus, in pixels.
        :param source_id: string describing the source (e.g., comp star ID, MP number) that this
               aperture is intended to measure. A tag only, not used in calculations. [string]
        :param obs_id: observation ID string, will typically be unique among all observations
               (aperture objects) in one image. A tag only, not used in calculations. [string]
        """
        xy_center = XY.from_tuple(xy_center)  # ensure is XY object.
        self.foreground_radius = foreground_radius
        self.gap = gap
        self.background_width = background_width
        self.annulus_inner_radius = self.foreground_radius + self.gap
        self.annulus_outer_radius = self.annulus_inner_radius + self.background_width
        cutout_size = int(ceil(2 * self.annulus_outer_radius)) + 4  # generous, for safety.
        dxy_origin = DXY(int(round(xy_center.x - cutout_size / 2)),
                         int(round(xy_center.y - cutout_size / 2)))
        xy_center_in_cutout = xy_center - dxy_origin
        foreground_mask = make_circular_mask(mask_size=cutout_size, xy_origin=xy_center_in_cutout,
                                             radius=self.foreground_radius)
        background_mask_center_disc = np.logical_not(make_circular_mask(cutout_size,
                                                                        xy_center_in_cutout,
                                                                        self.annulus_inner_radius))
        background_mask_outer_disc = make_circular_mask(cutout_size,
                                                        xy_center_in_cutout,
                                                        self.annulus_outer_radius)
        background_mask = np.logical_or(background_mask_center_disc, background_mask_outer_disc)
        super().__init__(image_xy, xy_center, dxy_origin,
                         foreground_mask, background_mask, source_id, obs_id)

    def _make_new_object(self, new_xy_center):
        """ Make new object using new xy_center. Overrides parent-class method.
            Masks will be recreated by the constructor, using new xy_center.
         """
        return PointSourceAp(self.image_xy, new_xy_center,
                             self.foreground_radius, self.gap, self.background_width,
                             self.source_id, self.obs_id)

    def __str__(self):
        return 'PointSourceAp at x,y = ' + str(self.xy_center.x) + ', ' + str(self.xy_center.y)


class MovingSourceAp(Ap):
    """ Elongated 'pill-shaped' photometric aperture for moving point source of light,
            esp. for a minor planet/asteroid.
        (If we will need a background mask that bleeds to the cutout's edges rather than being restricted
            to a (pill-shaped) annulus--for example, to work close to the parent image's edges, that will
            definitely require a new sibling class to this one so that recentering retains the mask shapes.)
        """
    def __init__(self, image_xy, xy_start, xy_end, foreground_radius, gap, background_width,
                 source_id=None, obs_id=None):
        """ Main and probably sole constructor for this class.
        :param image_xy: the parent image array [numpy ndarray].
        :param xy_start: x,y pixel position in parent image of the beginning of the MP's motion.
               [XY object, 2-tuple, 2-list, or 2-array of floats]
        :param xy_end:x,y pixel position in parent image of the beginning of the MP's motion.
               [XY object, 2-tuple, 2-list, or 2-array of floats]
        :param foreground_radius: radius of the aperture source end-caps, in pixels.
               Does not include effect of MP motion. [float]
        :param gap: Gap in pixels between foreground mask and background mask. [float]
        :param background_width: Width in pixels of background mask. [float]
        :param source_id: optional string describing the source (e.g., comp star ID, MP number) that this
               aperture is intended to measure. A tag only, not used in calculations. [string]
        :param obs_id: optional observation ID string, will typically be unique among all observations
               (aperture objects) in one image. A tag only, not used in calculations. [string]
        """
        self.xy_start = XY.from_tuple(xy_start)
        self.xy_end = XY.from_tuple(xy_end)
        self.foreground_radius = foreground_radius
        self.gap = gap
        self.background_width = background_width
        self.background_inner_radius = self.foreground_radius + self.gap
        self.background_outer_radius = self.foreground_radius + self.gap + self.background_width
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
        dxy_cutout_size = DXY(int(round(x_max - x_min + 4)), int(round(y_max - y_min + 4)))
        dxy_offset = DXY(int(round(xy_center.x) - dxy_cutout_size.dx / 2.0),
                         int(round(xy_center.y) - dxy_cutout_size.dy / 2.0))
        xy_start_cutout = self.xy_start - dxy_offset
        xy_end_cutout = self.xy_end - dxy_offset
        foreground_mask = make_pill_mask(dxy_cutout_size, xy_start_cutout, xy_end_cutout,
                                         self.foreground_radius)
        background_inner_mask = make_pill_mask(dxy_cutout_size, xy_start_cutout, xy_end_cutout,
                                               self.background_inner_radius)
        background_outer_mask = make_pill_mask(dxy_cutout_size,  xy_start_cutout, xy_end_cutout,
                                               self.background_outer_radius)
        background_mask = np.logical_or(background_outer_mask,
                                        np.logical_not(background_inner_mask))
        super().__init__(image_xy, xy_center, dxy_offset, foreground_mask, background_mask,
                         source_id, obs_id)
        self.sigma = self.stats.semiminor_sigma.value
        self.fwhm = self.sigma * FWHM_PER_SIGMA

    def _make_new_object(self, new_xy_center):
        """ Make new object using new xy_center. Overrides parent-class method.
            Masks will be recreated by the constructor, using new xy_center.
        """
        if not isinstance(new_xy_center, XY):
            new_xy_center = XY.from_tuple(new_xy_center)
        current_xy_center = self.xy_start + (self.xy_end - self.xy_start) / 2
        dxy_shift = new_xy_center - current_xy_center
        new_xy_start = self.xy_start + dxy_shift
        new_xy_end = self.xy_end + dxy_shift
        return MovingSourceAp(self.image_xy, new_xy_start, new_xy_end,
                              self.foreground_radius, self.gap, self.background_width,
                              self.source_id, self.obs_id)

    def __str__(self):
        return 'MovingSourceAp at x,y = ' + str(self.xy_center.x) + ', ' + str(self.xy_center.y)


_____IMAGE_and_GEOMETRY_SUPPORT____________________________________ = 0


def make_circular_mask(mask_size, xy_origin, radius):
    """ Construct a traditional mask array for small, stationary object, esp. for a star.
        Unmask only those pixels *within* radius pixels of a given point. Invert the mask separately to
            mask the interior.
        Numpy mask bookean convention: pixel True -> MASKED OUT.
        Numpy indexing convention: (y,x)
    :param mask_size: edge size of new mask array, which will be a square. [int]
    :param xy_origin: (x, y) pixel coordinates of circle origin, rel. to mask origin. [2-tuple of floats]
    :param radius: radius of ends and half-width of center region. [float]
    :return: mask array, True -> MASKED out/invalid (numpy boolean convention),
             and indexed as (y,x) (numpy indexing convention). [np.ndarray of booleans]
    """
    xy_origin = XY.from_tuple(xy_origin)  # ensure is XY object.
    circle = Circle_in_2D(xy_origin=xy_origin, radius=radius)
    is_inside = circle.contains_points_unitgrid(0, mask_size - 1, 0, mask_size - 1, include_edges=True)
    mask = np.transpose(np.logical_not(is_inside))  # render in numpy mask-boolean and index conventions.
    return mask


def make_pill_mask(mask_shape_xy, xya, xyb, radius):
    """ Construct a mask array for MP in motion: unmask only those pixels within radius pixels of
        any point in line segment from xya to xyb. Convention: pixel True -> VALID (opposite of numpy).
    :param mask_shape_xy: (x,y) size of array to generate. [2-tuple of ints]
    :param xya: (xa, ya) pixel coordinates of start-motion point. [XY object, or 2-tuple of floats]
    :param xyb: (xb, yb) pixel coordinates of end-motion point. [XY object, or 2-tuple of floats]
    :param radius: radius of ends and half-width of center region. [float]
    :return: mask array, True -> MASKED out/invalid (numpy boolean convention),
             and indexed as (x,y) (indexing convention is x,y-image, not numpy). [np.ndarray of booleans]
    """
    xya = XY.from_tuple(xya)  # ensure is XY object.
    xyb = XY.from_tuple(xyb)  # "
    if xya == xyb:
        return make_circular_mask(max(mask_shape_xy), xya, radius)

    # Make circle and rectangle objects:
    circle_a = Circle_in_2D(xya, radius)
    circle_b = Circle_in_2D(xyb, radius)
    dxy_ab = xyb - xya
    length_ab = dxy_ab.length
    dxy_a_corner1 = (radius / length_ab) * DXY(dxy_ab.dy, -dxy_ab.dx)  # perpendicular to ab vector.
    dxy_a_corner2 = (radius / length_ab) * DXY(-dxy_ab.dy, dxy_ab.dx)  # "
    xy_corner1 = xya + dxy_a_corner1
    xy_corner2 = xya + dxy_a_corner2
    xy_corner3 = xyb + dxy_a_corner2
    rectangle = Rectangle_in_2D(xy_corner1, xy_corner2, xy_corner3)

    # Make mask, including edges so no gaps can appear at rectangle corners:
    circle_a_contains = circle_a.contains_points_unitgrid(0, mask_shape_xy[0] - 1,
                                                          0, mask_shape_xy[1] - 1, True)
    circle_b_contains = circle_b.contains_points_unitgrid(0, mask_shape_xy[0] - 1,
                                                          0, mask_shape_xy[1] - 1, True)
    rectangle_contains = rectangle.contains_points_unitgrid(0, mask_shape_xy[0] - 1,
                                                            0, mask_shape_xy[1] - 1, True)
    # Render each in numpy mask-boolean but (x,y) index conventions:
    circle_a_mask = np.logical_not(circle_a_contains)
    circle_b_mask = np.logical_not(circle_b_contains)
    rectangle_mask = np.logical_not(rectangle_contains)
    mask = np.logical_and(np.logical_and(circle_a_mask, circle_b_mask), rectangle_mask)
    # any_contains = np.logical_or(np.logical_or(circle_a_contains, circle_b_contains), rectangle_contains)
    # mask = np.transpose(np.logical_not(any_contains))
    return mask


def calc_background_value(data, mask=None, dilate_size=3):
    """ Calculate the sigma-clipped median value of a (possibly masked) data array.
    :param data: array of pixels. [2-D ndarray of floats]
    :param mask: mask array, True=masked, i.e., use only False pixels. [None, or 2-D nadarray of bool]
    :param dilate_size: the number of pixels out from a detected outlier pixels to also mask.
        If sample is heavily undersampled (e.g., FWHM < 3 pixels), one might try 1-2 pixels.
        If image is heavily oversampled (e.g., FWHM > 10-15 pixels), one might increase this
        to perhaps 0.25-0.5 FWHM. [int]
    :return: tuple of background adu level (flux per pixel), standard deviation within used pixels.
                 [2-tuple of floats]
    """
    if mask is None:  # use all pixels.
        this_mask = np.full_like(data, False, dtype=np.bool)
    elif mask.shape != data.shape:  # bad mask shape.
        return None
    # nb: despite python warnings, to write ...(mask is False) would fail here. Must be mask == False.
    # elif np.sum(mask is False) == 0:  # <-- i.e., this will not work.
    elif np.sum(mask == False) == 0:  # no valid pixels.
        return 0.0, 0.0
    else:
        this_mask = mask.copy()
    # Notes on make_source_mask():
    #    this_mask: user-supplied mask of pixels simply not to consider at all.
    #    npixels = 2, which will mask out practically all cosmic ray hits.
    #    filter_fwhm=2, small enough to still capture cosmic ray hits, etc.
    #        (make_source_mask() will probably fail if fwhm < 2.)
    #    dilate_size=3, small enough to preserve most background pixels, but user may override.
    #    source_mask: masks out detected light-source pixels. MAY INCLUDE pixels masked out by this_mask!
    #    stats_mask: the intersection of valid pixels from this_mask and source_mask.
    source_mask = make_source_mask(data, mask=this_mask, nsigma=2, npixels=5,
                                   filter_fwhm=2, dilate_size=int(dilate_size))
    stats_mask = np.logical_or(this_mask, source_mask)
    _, median, std = sigma_clipped_stats(data, sigma=3.0, mask=stats_mask)
    return median, std


# _____FITS_FILE_HANDLING_______________________________________________ = 0
#
#
# def all_fits_filenames(top_directory, rel_directory, validate_fits=False):
#     """  Return list of all FITS filenames (name.extension; no directory info) in given directory_path.
#          (Code for this exists already, somewhere.)
#     :param top_directory:
#     :param rel_directory:
#     :param validate_fits: If True, open FITS files and include only if valid.
#         If False, include filename if it appears valid without opening the FITS file.
#     :return: List of all FITS filenames in given directory_path [list of strings]
#     """
#     # TODO: write all_fits_files().
#     pass

