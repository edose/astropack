##########################################################
Site, Instrument and Observer .ini files (`astropack.ini`)
##########################################################

Module `astropack.ini` provides .ini file infrastructure, as well as three
*server classes* for .ini file formats useful to astronomy:

* Class |Site|, reading and serving data about an observing site,
* Class |Instrument|, reading and serving data about a
  mount-telescope-camera astronomical instrument, and
* Class |HumanObserver|, which contains a small amount of data about the
  person taking the astronomical images, especially in support of automatic
  generation of ALCDEF-format astroid lightcurve data submission files
  (where ALCDEF refers to the Asteroid Lightcurve Data Exchange Format,
  as well as to a database harboring and serving such data). See
  https://alcdef.org/index.php

The idea behind these simple .ini files and respective classes
is to separate code from data, that is, to put specific user data
in safe, separate places that can be updated freely without having
to modify or even to understand the python code at all.

*************************************
Getting Started with the |Site| class
*************************************

The |Site| class simply reads and parses a Site-format .ini file,
and returns data on demand.

>>> from astropack.ini import Site
>>> fullpath = os.path.join(INI_DIRECTORY, 'my_dome.ini')
>>> s = Site(fullpath)
>>> s.name
'My Observatory Site (Dome)'
>>> s.longitude, s.latitude
(-105.528978, 32.903156)
>>> site.extinction['Clear']  # for winter, summer
(0.18, 0.14)

There are many more attributes available to the user.
An example |Site| .ini file and all class methods and attributes
are given in |Site|'s API page.

*******************************************
Getting Started with the |Instrument| class
*******************************************

The |Instrument| class simply reads and parses an Instrument-format .ini file,
and returns data on demand.

>>> from astropack.ini import Instrument
>>> s = Instrument(os.path.join(INI_DIRECTORY, 'my_C14.ini')
>>> s.mount_model, s.camera_model, s.ota_aperture  # aperture diameter in meters
('PlaneWave L-500', 'SBIG STXL-6303E', 0.35)
>>> s.filters_available
('Clear', 'BB', 'SG', 'SR', 'SI')
# Transform keys are (filter, target passband, color passbands 1 and 2).
# Values are first- and optionally second-order color coefficients.
>>> s.transforms
{('Clear', 'SR', 'SR', 'SI'): (0.4, -0.6),
 ('BB', 'SR', 'SR', 'SI'): (-0.131,)}

There are many more attributes available to the user.
An example |Instrument| .ini file and all class methods and attributes
are given in |Instrument|'s API page.

***********************************************
Getting Started with the |HumanObserver| class
***********************************************

The |HumanObserver| class simply reads and parses a very small
HumanObserver-format .ini file, and returns data on demand.
Used by the author mostly in code that generates ALCDEF submissions.

>>> from astropack.ini import HumanObserver
>>> fullpath = (INI_DIRECTORY, 'myself.ini')
>>> obs = HumanObserver(fullpath)
>>> obs.alcdef_contact_name
>>> obs.alcdef_observers
'Dose, E.V.'
obs.alcdef_contact_info
'ahaha@noway.com'

An example |HumanObserver| .ini file and class attributes
are given in |HumanObserver|'s API page.

***************
Reference/API
***************

.. automodapi:: astropack.ini
   :no-inheritance-diagram:
