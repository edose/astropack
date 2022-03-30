###################################################################
Web Retrieval for Minor Planets/Asteroids (`astropack.web`)
###################################################################

Functions that retrieve minor planet data from the
`Minor Planet Center Web Service <https://minorplanetcenter.net/web_service/>`__
and `Minor Planet Ephemeris Service
<https://www.minorplanetcenter.net/iau/MPEph/MPEph.html>`__
via package
`astroquery <https://astroquery.readthedocs.io/en/latest/index.html>`__.

Each function returns results as a convenient pandas |DataFrame|.

**************************************
Getting Started with the web functions
**************************************


>>> from astropack.web import get_mp_info
>>> df_info = get_mp_info(333)
>>> df_info
{'name': 'Badenia',
 'number': 333,
 'designation': None,
 'period': '5.54',
 'H': '9.49',
 'G': '0.15'}

>>> from astropack.web import get_mp_ephem
>>> from astropack.ini import Site
>>> import os
>>> my_site = Site(os.path.join(INI_DIRECTORY, 'MySite.ini'))
>>> df_ephem = get_mp_ephem(333, '2022-04-03 12:34:56', step_hours=2,
               num_entries=6, site=my_site)  # site=None returns geocentric data.
>>> df_ephem
                 Date          RA  ...  Moon distance  Moon altitude
0 2022-04-03 12:00:00  157.998750  ...            116            -24
1 2022-04-03 14:00:00  157.989583  ...            115             -2
2 2022-04-03 16:00:00  157.980417  ...            114             22
3 2022-04-03 18:00:00  157.971667  ...            113             47
4 2022-04-03 20:00:00  157.962917  ...            112             68
5 2022-04-03 22:00:00  157.953750  ...            111             68
[6 rows x 16 columns]
>>> df_ephem.columns
Index(['Date', 'RA', 'Dec', 'Delta', 'r', 'Elongation', 'Phase', 'V',
       'Proper motion', 'Direction', 'Azimuth', 'Altitude', 'Sun altitude',
       'Moon phase', 'Moon distance', 'Moon altitude'],
      dtype='object')

***************
Reference/API
***************

.. automodapi:: astropack.web
   :no-inheritance-diagram:
