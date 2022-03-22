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

   <img src="../_static/backpack-star-logo2.png" width="500" style="margin-bottom:30px;"/>

**Astropack** is `the author's <https://github.com/edose>`_
first attempt at a package that is tested
and documented to an `Astropy <http://www.astropy.org/about.html>`_-like level of detail.

.. Note:: Sheesh. I'm going to dispense with the "the author" references.
          I'll just write "I". I'm the author, after all. OK.

**Astropack** began as the author's...as my personal utilities repo
in support of my independent variable-star and asteroid/minor planet observing
programs. Astropack continues to grow.

.. Important:: 1. What's here is good to go.
   It's tested and documented. Indeed it's tested *against* the documentation.

   2. So you can install it, or take code from it, or adapt code to taste,
   But let's make a deal--if you take or use 10 lines of code or more, you'll please
   make an attribution back to me and Astropack.
   Fair enough? (That's pretty much what the license says, anyway.)

   3. **Astropack** is a work in progress (April 2022). It's marked as 0.1 beta
   version, but it's more in what the author calls **Deep Beta**, really more like a Delta version.
   So long as it's in version zero-point-anything, I may remove items, I may change
   the API, I may even make breaking changes to the API--I'll try not to, but
   in any version 0.x, all's fair. But whatever is in the repo, that code works
   as documented. The test files are also in the repo and they've all passed, always.

You'll see a pattern arising in much of the code, and that pattern is:
server classes. That is, classes and their objects designed to (1) get data, (2) organize
it automatically, and serve that data or derivative data on demand.
You'll see this pattern in Astropack classes |Astronight|, |AtlasRefcat2|, and
|Site|, for example. A class constructs an object and then delivers a clean,
reliable and extremely well-defined bit of data whenever the user (the user's code,
really) needs to pull it off the shelf. Like an agile and talented clerical assistant.

Maybe an quick example, from class |Astronight| (which represents one observing site
during on observing night):

      >>> from astropack.almanac import Astronight
      >>> from astropack.ini import Site
      >>> this_site = Site('My Dome')
      >>> an = Astronight(this_site, 20220402)
      >>>
      >>> an.timespan_dark
      Timespan(2022-04-03 02:02:36.490, 2022-04-03 12:07:51.845)
      >>> an.moon_skycoord
      <SkyCoord (ICRS): (ra, dec) in deg
          (34.98185288, 11.95750862)>
      >>>
      >>> betelgeuse = SkyCoord('05h 55m 10.30536s +07d 24m 25.4304s')
      >>> an.timespan_observable(betelgeuse, min_alt=30, min_moon_dist=45)



In other words: as I refer to it, a server class represents a **bundle of data that speaks for
itself**. A reference book that opens itself to just the right page whenever you ask,
and stays completely out of the way when you don't.
It's how much of Astropack is coded. It's ideal for technical coding.

Or so I say. Check it for yourself. Right here, in fact:

**************************
Documentation, by Module
**************************

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
Index & Search
***************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
