__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

# Python core packages:
import os

# External packages:
import pytest
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

# Test target:
from astropak import image

THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_TOP_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, "test")


__________TEST_CLASS_FITS____________________________________________________ = 0


def test_class_fits_constructor():
    fits = get_test_fits_file()
    # print('fullpath >', fullpath, '<')
    # fits = image.FITS(fullpath, pinpoint_pixel_scale_multiplier=0.99)
    print('exposure (630)', str(fits.exposure))
    assert fits.airmass == pytest.approx(1.01368, abs=0.00002)
    assert fits.exposure == 630.0
    assert fits.filter == 'BB'
    assert fits.focal_length == pytest.approx(2704.7)
    assert fits.guide_exposure == pytest.approx(0.5)
    assert fits.fullpath == os.path.join(TEST_TOP_DIRECTORY, '$data_for_test', 'MP_1300-0004-BB.fts')
    assert fits.image_fits[0][1] == 1047  # index order same as given in FITS header.
    assert fits.image_xy[0][1] == 982  # index order x,y, like image where 0,0 is top left, y is vertical.
    assert fits.is_calibrated == True
    assert fits.is_valid == True
    assert fits.is_plate_solved == True
    assert fits.is_plate_solved_by_pinpoint == True
    assert fits.object == 'MP_1300'
    assert fits.pinpoint_pixel_scale_multiplier == 0.99
    assert fits.temperature == -35.0
    assert fits.utc_mid == Time('2022-02-10 05:15:38.790', scale='utc')
    assert fits.utc_start == Time('2022-02-10 05:10:23.790', scale='utc')
    assert isinstance(fits.wcs_fits, WCS)
    assert all(fits.wcs_fits.wcs.crpix == [1536.0, 1023.5])
    assert fits.wcs_fits.wcs.crval[0] == pytest.approx(123.83159869)
    assert fits.wcs_fits.wcs.cd[0][0] == pytest.approx(-0.000190406559045)
    assert isinstance(fits.wcs_corrected, WCS)
    assert all(fits.wcs_corrected.wcs.crpix == [1536.0, 1023.5])
    assert fits.wcs_corrected.wcs.crval[0] == pytest.approx(123.83159869)
    assert fits.wcs_corrected.wcs.cd[0][0] == pytest.approx(-0.00018850249345455)
    assert fits.ra_deg == pytest.approx(123.83138202)
    assert fits.dec_deg == pytest.approx(30.12364316)


def test_class_fits_header_value():
    fits = get_test_fits_file()
    assert fits.header_value('OBJECT') == fits.object
    assert fits.header_value('BAD HEADER KEY') is None


def test_class_fits__detect_pinpoint_plate_solution():
    fits = get_test_fits_file()
    assert fits._detect_pinpoint_plate_solution() == True
    del fits.header['TR1_0']
    assert fits._detect_pinpoint_plate_solution() == False


def test_class_fits__make_corrected_wcs():
    pass  # already tested in constructor test.


def test_class_fits_xy_to_skycoords():
    fits = get_test_fits_file()
    xy = (0, 800)
    sc_tuple = fits.xy_to_skycoords(xy)
    assert sc_tuple.size == 1
    assert sc_tuple.ra.degree == pytest.approx(124.16599266442681)
    assert sc_tuple.dec.degree == pytest.approx(30.16703856)

    xy = [(0, 800)]
    sc_list1 = fits.xy_to_skycoords(xy)
    assert sc_list1 == sc_tuple

    xy = [(0, 800), (100, 810)]
    sc_list2 = fits.xy_to_skycoords(xy)
    assert sc_list2.size == 2
    assert sc_list2[0] == sc_list1
    assert sc_list2[1].ra.degree == pytest.approx(124.14419646)
    assert sc_list2[1].dec.degree == pytest.approx(30.16509713)

    with pytest.raises(ValueError):
        sc = fits.xy_to_skycoords(np.array([1, 2]))


def test_class_fits_skycoords_to_xy():
    fits = get_test_fits_file()
    sc1 = SkyCoord(ra=124.00 * u.degree, dec=30.10 * u.degree, frame='icrs')
    xy1 = fits.skycoords_to_xy(sc1)
    assert sc1.size == 1
    assert xy1[0] == pytest.approx(762.8360547833112)
    assert xy1[1] == pytest.approx(1152.8791205364719)

    sc3 = SkyCoord([124.00, 124.05, 124.10], [30.10, 30.12, 30.14], frame="icrs", unit="deg")
    xy3 = fits.skycoords_to_xy(sc3)
    assert sc3.size == 3
    assert xy3[0] == xy1
    assert xy3[2][1] == pytest.approx(942.4950312316655)


def test_class_fits_corner_skycoords():
    fits = get_test_fits_file()
    corners = fits.corner_skycoords()
    assert isinstance(corners, SkyCoord)
    assert len(corners) == 4
    assert corners[0].ra.degree == pytest.approx(124.1654873)
    assert corners[0].dec.degree == pytest.approx(30.31785206)
    assert corners[3].ra.degree == pytest.approx(123.4988068)
    assert corners[3].dec.degree == pytest.approx(29.92878204)


def test_class_fits__is_calibrated():
    """Test internal function."""
    fits = get_test_fits_file()
    assert fits._is_calibrated() == True
    fits.header['CALSTAT'] = None
    assert fits._is_calibrated() == False
    del fits.header['CALSTAT']
    assert fits._is_calibrated() == False
    fits.header['CALSTAT'] = 'BDF'
    assert fits._is_calibrated() == True


__________TEST_CLASS_POINTSOURCEAP_____________________________________ = 0

# Test subclass of class Ap (class Ap is never instantiated directly).


def test_class_pointsourceap_constructor():
    fits = get_test_fits_file()
    ap = image.PointSourceAp(image_xy=fits.image_xy, xy_center=(1489, 955),
                             foreground_radius=12, gap=5, background_width=8,
                             source_id='some star', obs_id='whatever')
    assert ap.foreground_radius == 12
    assert ap.gap == 5
    assert ap.background_width == 8
    assert ap.annulus_inner_radius == 17
    assert ap.annulus_outer_radius == 25
    assert ap.image_xy.shape == (3072, 2047)
    assert ap.xy_center == (1489, 955)
    assert ap.xy_offset == (1489 - 27, 955 - 27)
    assert ap.input_foreground_mask.shape == (54, 54)
    assert ap.input_background_mask.shape == ap.input_foreground_mask.shape
    assert ap.source_id == 'some star'
    assert ap.obs_id == 'whatever'
    assert ap.is_valid == True
    assert (ap.background_mask == ap.input_background_mask).all()
    assert (ap.x_low, ap.y_low) == ap.xy_offset
    assert ap.x_high == ap.x_low + 53
    assert ap.y_high == ap.y_low + 53
    assert ap.cutout.shape == ap.input_foreground_mask.shape
    assert ap.mask_overlap_pixel_count == 0
    assert ap.foreground_pixel_count == 441
    assert ap.background_pixel_count == 1060
    assert ap.background_level == pytest.approx(1078, abs=1)
    assert ap.background_std == pytest.approx(26.4, abs=1)
    assert ap.foreground_max == 23277
    assert ap.foreground_min == 1097
    assert ap.raw_flux == pytest.approx(1621816, abs=100)
    assert ap.net_flux == pytest.approx(1146418, abs=100)
    assert ap.xy_centroid == pytest.approx((1488.3887587250026, 955.2717900451668), abs=0.01)
    assert ap.sigma == pytest.approx(3.400, abs=0.002)
    assert ap.fwhm == pytest.approx(7.439, abs=0.002)
    assert ap.elongation == pytest.approx(1.173, abs=0.002)
    assert ap.flux_stddev(gain=1) == pytest.approx(1274, abs=2)
    assert ap.flux_stddev(gain=1.5) == pytest.approx(1040, abs=2)


def test_class_pointsourceap_recenter():
    fits = get_test_fits_file()
    ap = image.PointSourceAp(image_xy=fits.image_xy, xy_center=(1482, 952),
                             foreground_radius=12, gap=5, background_width=8,
                             source_id='some star', obs_id='whatever')
    assert ap.xy_center == pytest.approx((1482, 952))
    assert ap.xy_centroid == pytest.approx((1487.8586617491405, 955.097767706891), abs=0.01)

    ap_1 = ap.recenter(max_iterations=1)
    assert ap_1.xy_center == ap.xy_centroid
    assert ap_1.xy_centroid == pytest.approx((1488.350973107373, 955.2710427996668), abs=0.01)

    ap_2 = ap.recenter(max_iterations=2)
    assert ap_2.xy_centroid == pytest.approx((1488.3554435661015, 955.2871555171624), abs=0.01)

    ap_3 = ap.recenter(max_iterations=3)
    assert ap_3.xy_centroid == ap_2.xy_centroid


__________TEST_CLASS_MOVINGSOURCEAP_____________________________________ = 0

# Test subclass of class Ap (class Ap is never instantiated directly).


def test_class_movingsourceap_constructor():
    fits = get_test_fits_file()
    ap = image.MovingSourceAp(image_xy=fits.image_xy, xy_start=(1223, 972.8), xy_end=(1226.4, 973.8),
                              foreground_radius=12, gap=5, background_width=8,
                              source_id='some MP', obs_id='rock observation')
    assert ap.xy_start == pytest.approx((1223, 972.8))
    assert ap.xy_end == pytest.approx((1226.4, 973.8))
    assert ap.foreground_radius == 12
    assert ap.gap == 5
    assert ap.background_width == 8
    assert ap.background_inner_radius == 12 + 5
    assert ap.background_outer_radius == 12 + 5 + 8
    assert ap.source_id == 'some MP'
    assert ap.obs_id == 'rock observation'
    assert ap.is_valid == True
    assert ap.xy_center == pytest.approx(ap.xy_start + (ap.xy_end - ap.xy_start) / 2.0)
    assert ap.xy_offset == (1196, 945)
    assert ap.input_background_mask.shape == (57, 55)  # nb: x,y index order.
    assert ap.input_background_mask.shape == ap.input_foreground_mask.shape
    assert ap.cutout.shape == ap.input_foreground_mask.shape
    assert ap.cutout[34, 27] == ap.image_xy[1230, 972] == 2810.0                   # all x,y index order.
    assert ap.cutout[20, 29] == ap.image_xy[20 + ap.xy_offset.x, 29 + ap.xy_offset.y] == 1235.0  # "
    assert ap.mask_overlap_pixel_count == 0
    assert ap.foreground_pixel_count == 537
    assert ap.background_pixel_count == 1113
    assert ap.background_level == 1067.0
    assert ap.background_std == pytest.approx(25.938, abs=0.001)
    assert ap.foreground_max == 3822.0
    assert ap.foreground_min == 1008.0
    assert ap.raw_flux == pytest.approx(762126.0)
    assert ap.net_flux == pytest.approx(189147.0, abs=1)
    assert ap.xy_centroid == pytest.approx((1225.622, 973.242), abs=0.002)
    assert ap.sigma == pytest.approx(3.002, abs=0.002)
    assert ap.fwhm == pytest.approx(7.069, abs=0.002)  # seems suspiciously high.
    assert ap.elongation == pytest.approx(1.368, abs=0.002)  # reasonably high given MP motion.
    assert ap.flux_stddev(gain=1) == pytest.approx(873.184, abs=0.002)
    assert ap.flux_stddev(gain=1.5) == pytest.approx(713.028, abs=0.002)


def test_class_movingsourceap_recenter():
    fits = get_test_fits_file()
    ap = image.MovingSourceAp(image_xy=fits.image_xy, xy_start=(1223, 972.8), xy_end=(1226.4, 973.8),
                              foreground_radius=12, gap=5, background_width=8,
                              source_id='some MP', obs_id='rock observation')
    assert ap.xy_centroid == pytest.approx((1225.622, 973.242), abs=0.002)

    ap_1 = ap.recenter(max_iterations=1)
    assert ap_1.xy_center == ap.xy_centroid
    assert ap_1.xy_centroid == pytest.approx((1225.671140374866, 973.2309884931738), abs=0.002)

    ap_2 = ap.recenter(max_iterations=2)
    assert ap_2.xy_centroid == pytest.approx((1225.6740937846748, 973.2286770544649), abs=0.002)

    ap_3 = ap.recenter(max_iterations=3)
    assert ap_3.xy_centroid == ap_2.xy_centroid


__________HELPER_FUNCTIONS____________________________________________ = 0


def get_test_fits_file():
    """ Helper function to get FITS object."""
    fullpath = os.path.join(TEST_TOP_DIRECTORY, '$data_for_test', 'MP_1300-0004-BB.fts')
    # print('fullpath >', fullpath, '<')
    fits = image.FITS(fullpath, pinpoint_pixel_scale_multiplier=0.99)
    return fits


_____IMAGE_and_GEOMETRY_SUPPORT____________________________________ = 0


def test_calc_background_value():
    # Test simple case, no excluded pixels, no mask:
    im, _ = np.meshgrid(np.arange(10), np.arange(12))
    im = im.copy()
    median, std = image.calc_background_value(im)
    assert median == 4.5
    assert std == pytest.approx(2.872, abs=0.005)

    # Test simple case, no mask:
    im[5, 3:7] = 40  # outlier pixels.
    median, std = image.calc_background_value(im)
    assert median == 4.0
    assert std == pytest.approx(3.136, abs=0.005)

    # Test with added mask:
    mask = np.array(im <= 1.0)
    median, std = image.calc_background_value(im, mask)
    assert median == 6.0
    assert std == pytest.approx(2.445, abs=0.005)

    # Test uniform pixel (problematic) case:
    im = np.full(shape=(10, 20), fill_value=7, dtype=np.float).copy()
    median, std = image.calc_background_value(im)
    assert median == 7.0
    assert std == 0.0

    # Ensure dilate_size default, and that it gets used as an integer:
    assert image.calc_background_value(im) == \
           image.calc_background_value(im, dilate_size=3) == \
           image.calc_background_value(im, dilate_size=3.3)
