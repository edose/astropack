###################################################################
FITS reading, Aperture Photometry (`astropack.image`)
###################################################################

Class |FITS| for reading and then serving image-array and FITS-header data.

Class |PointSourceAp| for high-accuracy, low-bias aperture photometry of sky-stationary
point light sources, esp. target stars and comparison stars.

Classes |MovingSourceAp| for high-accuracy, low-bias aperture photometry of
moving light sources, especially minor planets/asteroids. Constructs and applies
a 'pill-shaped' aperture and background-annulus set that automatically include
the effects of target motion.

*********************************
Getting Started with class |FITS|
*********************************

|FITS| is another *server class* that reads data and organizes it, then serves
parsed and derived data on demand.

    >>> from astropack.image import FITS
    >>> fullpath = os.path.join(IMAGE_DIRECTORY, 'MP_1300-0004-BB.fts')
    >>> f = FITS(fullpath)
    >>> f.object, f.exposure, f.filter
    ('MP_1300', 630.0, 'BB')
    >>> f.header_value('airmass')
    1.01368456711
    >>> f.header_value('invalid header key') is None
    True
    #
    >>> f.xy_to_skycoords(xy=(500, 874))
    <SkyCoord (FK5: equinox=2000.0): (ra, dec) in deg
        (124.0593176, 30.15306119)>
    >>> f.corner_skycoords()
    <SkyCoord (FK5: equinox=2000.0): (ra, dec) in deg
        [(124.16886656, 30.31980748), (123.49149578, 30.31639601),
         (124.17015564, 29.93020564), (123.49545203, 29.92680762)]>
    #
    >>> from astropy.coordinates import SkyCoord
    >>> import astropy.units as u
    sc = SkyCoord(124.059, 30.153, unit=u.deg)
    >>> f.skycoords_to_xy(sc)
    (501.41738211369625, 874.3474219930808)

*****************************************
Getting Started with the Aperture classes
*****************************************

The Aperture classes |PointSourceAp| and |MovingSourceAp|
are designed for aperture photometry.

|PointSourceAp| applies to effectively stationary targets of approximately known
pixel location. It is suitable to stars in an image, including to comparison stars.
It generates and applies foreground and background pixel masks of circular shape.
|PointSourceAp| works most conveniently with instances of the class |FITS|,
but such are not required.

    >>> from astropack.image import FITS, PointSourceAp
    >>> fullpath = os.path.join(IMAGE_DIRECTORY, 'MP_1300-0004-BB.fts')
    >>> f = FITS(fullpath)
    >>> ap = PointSourceAp(image_xy=fits.image_xy, xy_center=(1489, 955),
                           foreground_radius=12, gap=5, background_width=8,
                           source_id='some star', obs_id='')
    >>> ap.is_valid, ap.foreground_pixel_count, ap.foreground_max
    (True, 441, 23277.0)
    >>> ap.xy_centroid, ap.sigma, ap.fwhm
    ((1488.3887587250026, 955.2717900451668), 3.400290617513252, 7.439113620441175)
    >>> ap.background_pixel_count, ap.background_level, ap.background_std
    (1060, 1078.0, 25.96445561050696)
    >>> ap.raw_flux, ap.net_flux, ap.flux_stddev(gain=1.5)
    (1621816.0, 1146418.0, 1039.947662017538)

|MovingSourceAp| applies to targets that move significantly during an image's
exposure time, especially for minor planets/asteroids.
It generates and applies foreground and background pixel masks of 'pill' or
'racetrack' shape, that is, shapes constructed of the union of 3 elementary shapes:
a circle centered over the target's position at exposure start,
a circle centered over the target's position at exposure end,
and a rectangle comprising all the pixels between the two circles.
|MovingSourceAp| works most conveniently with instances of the class |FITS|,
but such are not required.

|PointSourceAp| and |MovingSourceAp| differ only in the constructor and in the
instance attributes describing target pixel position (a single point for
|PointSourceAp| vs. start and end points for |MovingSourceAp|).
The APIs (pixel counts, fluxes, etc) are otherwise identical.

***************
Reference/API
***************

.. automodapi:: astropack.image
   :no-inheritance-diagram:
   :inherited-members:

.. Note:: ``pinpoint_pixel_scale_multiplier`` is a value (prob. near 1) by which to
    multiply pixel scale, iff PinPoint plate solution is detected.
    It's the best solution I can devise for the PinPoint's pixel scale deviating from
    those of linear WCS plate solvers.
    This arises because (sigh) Pinpoint plate solutions include "private"
    distortion parameters, so that their WCS values are not what they would be &
    should be for a linear-WCS-only solution.
    That is, PinPoint zero-order solution cannot be used as one, nor even correctable
    given its "private" distortion parameters.
