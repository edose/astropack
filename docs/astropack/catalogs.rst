###########################################
Catalog retrieval (`astropack.catalogs`)
###########################################

ATLAS refcat2 catalog handling, via class |AtlasRefcat2|.

**********************************
Getting Started with AtlasRefcat2
**********************************

The wonderful ATLAS refcat2 astronomical catalog is the only catalog that most advanced
amateurs will need. It combines Gaia astrometry for outstandingly accurate star
positions with ATLAS photometry, mostly in Sloan-like filters and passbands.
It is essentially all-sky complete to magnitude 19.

In this `astropack.catalogs` module, the |AtlasRefcat2| class represents the catalog
data for all catalog stars within a user-specified range of Right Ascensions and
Declinations.


.. Note:: Before first use of |AtlasRefcat2|,
   the user must (1) download the catalog "chunk" subdirectories, and
   (2) make a small index file pointing to the local subdirectories' locations.
   This index file is very small and simple, e.g.,

   >>> with open('D:/Astro/Catalogs/ATLAS-refcat2/index.txt') as f:
   >>> print(f.read())
   16 : mag-0-16
   17 : mag-16-17
   18 : mag-17-18

   The catalog has one magnitude "chunk" per subdirectory.
   The user-created index file has one line per subdirectory, each line containing only:
   the maximum magnitude of the chunk, a colon, and the subdirectory name.
   The index file must reside in the same directory as the catalog subdirectories,
   in the above case, *D:/Astro/Catalogs/ATLAS-refcat2/*.

   Once this index file is created and the subdirectories are downloaded and copied in,
   they need never be touched again. The |AtlasRefcat2| class will merge star entries
   from the subdirectories as needed, automatically.

To get catalog data, instantiate the AtlasRefcat2 instance with the catalog
location and ranges of RA and Dec.
The ``target_epoch`` may be and usually should be specified as the approximate date
for which data are needed, for example, the date when images were taken for which
catalog data are needed.
This ensures that proper motions are accounted for in star positions.

    >>> from astropack.catalogs import AtlasRefcat2
    >>> from datetime import datetime, timezone
    >>> cat = AtlasRefcat2('D:/Astro/Catalogs/ATLAS-refcat2/', 'index.txt', \
    >>>   ra_deg_range=(3, 4.5), dec_deg_range=(34, 34.7), \
    >>>   target_epoch=datetime(2022, 4, 3).replace(tzinfo=timezone.utc))
    >>> len(cat.df_all)
    2178
    >>> min(cat.df_all.RA_deg), max(cat.df_all.RA_deg)
    3.000373075624063, 4.499971973935839

Then apply the class's methods to select the catalog stars by
magnitude ranges, etc, as needed:

    >>> len(cat.df_all), len(cat.df_selected)
    2178, 2178
    >>> cat.select_on_r_mag(11, 14.5)
    >>> cat.select_on_ri_color(0.2, 0.45)
    >>> len(cat.df_all), len(cat.df_selected)
    2178, 67

Users still working in legacy Johnson-Cousins UBVRI filters will have to apply
transforms to convert ATLAS refcat2 Sloan magnitudes to UBVRI magnitudes, and vice versa.
All new astronomical catalogs known to the author, including ATLAS refcat2, provide
Sloan magnitudes natively.
(Astropack's author gave up on Johnson-Cousins photometry years ago.)

***************
Reference/API
***************

.. automodapi:: astropack.catalogs
   :no-inheritance-diagram:
