astropack: a pack of mostly astronomical routines
---------------------------------------------------

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge


This repo
--------------

**Astropack** began as `my <https://github.com/edose>`_
personal utilities repo in support of my independent variable-star and
asteroid/minor planet observing programs.
Consider this a stub code library. Astropack will grow to cover all manner
of needs not addressed by |astropy| and other libraries.

A pattern arising in much of the code is what I call **server classes**.
That is: classes and their objects that I've designed to (1) get data, (2) organize
and extend the data automatically, and (3) serve it to the user's code on demand.
Like a particularly talented tecnical assistant.
A class constructs an object and then delivers a clean,
reliable and extremely well-defined bit of data whenever the user (the user's code,
really) needs it. That bit of data, and not more.

In other words: what I mean by a server class is a **bundle of data that speaks for
itself**. A reference book that opens itself to just the right page when you need
it to, and hides away when you don't.

Upcoming plans
---------------

We're now at version 0.1 beta.
**Beta** version may be a stretch--maybe closer to sigma version.
Early days.
Even so, what's here is correct and
tested, but is a mere shell (maybe 10-15 %) of what it will be.

Next things to do:
   * Clean up the documentation. This is mostly done.
   * Install documentation onto `Read The Docs <https://readthedocs.org/>`_,
     you know, just as the grownups do. By mid-April 2022.
   * Devise better interfaces or wrappers (or later, an entire replacement)
     for the external package `skyfield <https://rhodesmill.org/skyfield/>`_,
     which has wonderful, wonderful math but a very difficult API as well as,
     um, documentation that I find
     impenetrable, opaque, obfuscated, befogged, and with a toy surprise:
     dead end or false lead at every turn.
     Oh, can I do better? (they ask)

     Um. Yeah.

     And by June 2022 if just refactoring will do.
     But more like mid-late 2022 if I decide I can't stand the pain and
     do away with skyfield, just as this very repo did away with
     its heiroglyphical predecessor
     package `ephem <https://rhodesmill.org/pyephem/>`_.

   * More clearly standardize units to SI (mostly done already), and
     set up a strongly preferred
     class for each kind of quantity, especially for time,
     which has become garbled across three different representations
     (just for UTC, forget about JD and TT): python's |datetime|,
     astropy's |Time|, and skyfield's Time object (they just *had*
     to reinvent the wheel and then *use the same name* for it, sigh).

License
-------

This project is Copyright (c) Eric Dose and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause license. See the licenses folder for
more information.

If you'd like to contribute
----------------------------

I'd love to have your contributions! Maybe not just yet--the code isn't quite
stable enough.
I'll signal when I'm ready to take pull requests. Soon.
Until then, please do post your thoughts to this github repo's Issues page.

Look, this thing is still in Deep Beta. Will be for a few months, no doubt.

*BUT!* Later, we'll be able to post this boilerplate from |astropy| and say it to you,
for real.

Even shout. Come on, I'm talking to you:

.. Important::
    We love contributions! astropack is open source,
    built on open source, and we'd love to have you hang out in our community.

    **Imposter syndrome disclaimer:** We want your help. No, really.

    There may be a little voice inside your head that is telling you that you're not
    ready to be an open source contributor; that your skills aren't nearly good
    enough to contribute. What could you possibly offer a project like this one?

    We assure you - the little voice in your head is wrong. If you can write code at
    all, you can contribute code to open source. Contributing to open source
    projects is a fantastic way to advance one's coding skills. Writing perfect code
    isn't the measure of a good developer (that would disqualify all of us!); it's
    trying to create something, making mistakes, and learning from those
    mistakes. That's how we all improve, and we are happy to help others learn.

    Being an open source contributor doesn't just mean writing code, either. You can
    help out by writing documentation, tests, or even giving feedback about the
    project (and yes - that includes giving feedback about the contribution
    process). Some of these contributions may be the most valuable to the project as
    a whole, because you're coming to the project with fresh eyes, so you can see
    the errors and assumptions that seasoned contributors have glossed over.

    Note: This disclaimer was originally written by
    `Adrienne Lowe <https://github.com/adriennefriend>`_ for a
    `PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_, and was adapted by
    astropack based on its use in the README file for the
    `MetPy project <https://github.com/Unidata/MetPy>`_.*

*Maybe in May or thereabouts. Hope so. But who knows?*