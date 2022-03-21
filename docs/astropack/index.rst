.. the "raw" directive below is used to hide the title in favor of
.. just the logo's being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 { display:none; }
    </style>

########################
Astropack Documentation
########################

.. raw:: html

   <img src="../_static/backpack-star-logo.png" width="500" style="margin-bottom:30px;"/>

This is the intro paragraph before TOC.

.. Important:: I suppose we should put something clever here.
    Maybe about |py.datetime| or |py.timedelta|, what?

******************
User Documentation
******************

.. toctree::
   :maxdepth: 1

   almanac
   catalogs
   geometry
   image
   ini
   reference
   stats
   util
   web

And this is after the TOC.

***************
Index
***************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
