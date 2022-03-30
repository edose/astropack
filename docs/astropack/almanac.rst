###########################################
Observation Planning (`astropack.almanac`)
###########################################

The |Astronight| class forms the centerpiece of astropack's observation planning.
The user forms an Astronight instance by passing in a |Site| instance and
an "astronight date", which is practically always the *local* date when the sun sets.

Then, once the |Astronight| instance is in place, most almanac data
(e.g., dark time span,
moon's illumination percent, transit time, and time span above horizon) are available
as instance attributes, and observable time spans and transit times are available for
any sky object with known Right Ascension and Declination.

A number of almanac functions are also available in this module, but most users
would better construct an |Astronight| instance and refer to the attributes.

**************************************
Getting Started with class Astronight
**************************************

In user code for observation planning, constructing the |Astronight| instance is the first
step:

      >>> from astropack.almanac import Astronight
      >>> from astropack.ini import Site
      >>> site_ini_path = os.path.join(INI_DIRECTORY, 'my_site.ini')
      >>> my_site = Site(site_ini_path)
      >>> an_date = 20220402  # int, or as string '20220402'
      >>> an = Astronight(my_site, an_date)

That's all there is to making the instance.
The user may set the sun (negative) altitude denoting observable darkness either
explicitly through |Astronight|'s ``sun_altitude_dark`` parameter, or simply use
the default value in the Site .ini file.

Now that the instance is constructed,
many useful attributes are available for all kinds of observation planning tasks:

    >>> an.timespan_no_sun
    Timespan(2022-04-03 01:23:07.224, 2022-04-03 12:47:19.178)
    >>> an.timespan_no_sun.start
    <Time object: scale='utc' format='iso' value=2022-04-03 01:23:07.224>
    >>> an.timespan_no_sun.end
    <Time object: scale='utc' format='iso' value=2022-04-03 12:47:19.178>

|Timespan| is an astropack representation of a span of time having a definite ``start``
and ``end``. |Time| is astropy's standard Time object, representing a single instant of
time, which may be expressed in any of numerous timescales and formats.

    >>> an.timespan_dark
    Timespan(2022-04-03 02:02:36.490, 2022-04-03 12:07:51.845)
    >>> an.timespan_dark_no_moon
    Timespan(2022-04-03 03:01:28.676, 2022-04-03 12:07:51.845)

``.timespan_no_sun`` is the span of time (for this night, at this site)
during which the sun is below the horizon.
``.timespan_dark`` is the span of time during which the sun is far enough below the
horizon to be considered observable (if there is no moon).
``.timespan_dark_no_moon`` is the span of time during which the sun is far below the horizon
and the moon is down.

|Astronight| offers yet more planning data:

    >>> an.utc_local_middark
    <Time object: scale='utc' format='iso' value=2022-04-03 07:05:14.167>
    >>> an.local_middark_lst_hour_string
    '12:49:31.26'

``.utc_local_middark`` is the midpoint UTC |Time| of the ``.timespan_dark`` Timespan.
Astronight uses "mid-dark" rather than local clock midnight to remove dependence on
time zones, etc. (Also, "mid-dark" time is not quite the same as the
sun's anti-transit time, its time of farthest extent below the horizon.)

``.local_middark_lst_hour_string`` is the local sidereal time (the Right Ascension
directly overhead) at mid-dark time.

The moon is important to observing, so these attributes are available, again, for this
site during the specified observing night:

    >>> an.moon_illumination  # in percent, at middark
    4.3253057700347375
    >>> an.moon_transit
    <Time object: scale='utc' format='iso' value=2022-04-02 20:18:47.519>
    >>> an.moon_ra, an.moon_dec  # in degrees, at middark
    (34.98185287647763, 11.957508621623298)
    >>> an.moon_skycoord
    <SkyCoord (ICRS): (ra, dec) in deg
        (34.98185288, 11.95750862)>
    >>> an.moon_up_timespans
    [Timespan(2022-04-03 02:02:36.490, 2022-04-03 03:01:28.676)]

``.moon_transit`` is the transit time closest to the night's mid-dark time.
``.moon_skycoord`` is the astropy |SkyCoord| representation of the same data
in ``.moon_ra`` and ``.moon_dec``

``.moon_up_timespans`` is always given as a list, because it is possible, usually in
winter near full moon, for the moon to be up at sunset and sunrise, but down between,
giving two moon-lit time spans in the same night. The attribute name ends in a
plural to remind the user. But usually there is only one timespan, and the
|Timespan| itself is obtained as ``an.moon_up_timespans[0]``.

Finally, user-specified targets may have their observation availability characterized:

    >>> from astropy.coordinates import SkyCoord
    >>> betelgeuse = SkyCoord('05h 55m 10.30536s +07d 24m 25.4304s')
    >>> an.timespan_observable(betelgeuse, min_alt=30, min_moon_dist=45)
    Timespan(2022-04-03 02:02:36.490, 2022-04-03 04:08:12.314)
    >>> an.timespan_observable(betelgeuse, min_alt=10, min_moon_dist=45)
    Timespan(2022-04-03 02:02:36.490, 2022-04-03 05:43:43.096)
    >>> an.transit(betelgeuse)
    <Time object: scale='utc' format='iso' value=2022-04-03 00:13:12.843>

Note that targets' observable timespans are always limited to sky darkness.
Targets are not considered observable while (1) the moon is nearer than
``min_moon_dist`` in degrees, and (2) the moon is above the horizon.

*********************************************************
Getting Started with `astropack.almanac` module functions
*********************************************************

A number of free-standing almanac functions are also available in this module,
but most require familiarity with the skyfield package
(https://rhodesmill.org/skyfield/); in most cases,
it is better and much easier simply to construct an |Astronight| instance
and then call on the attributes as needed.
See the Reference/API section below if you're interested in usage details,
or if you're interested in the implementing code.

***************
Reference/API
***************

.. automodapi:: astropack.almanac
   :no-inheritance-diagram:
