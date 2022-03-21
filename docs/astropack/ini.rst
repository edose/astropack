##########################################################
Site, Instrument and Observer .ini files (`astropack.ini`)
##########################################################

classes +  code. (rst file)

***************
Introduction
***************

Example :class:`~.ini.Site` .ini file
------------------------------------------
::

   [Site]
   Name = New Mexico Skies (Dome)
   MPC Code = N/A

   [Location]
   # Longitude, Latitude: decimal degrees; Elevation: meters.
   Longitude = -105.53
   Latitude = +32.90
   Elevation = 2200
   # For correct timezone in standard (not daylight savings/summer) time.
   # Needed only to determine this site's side of International Date Line.
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
   # Dome slew rate in (degrees azimuth)/second.
   Slew Rate = 2.85

***************
Reference/API
***************

.. automodapi:: astropack.ini
   :no-inheritance-diagram:
