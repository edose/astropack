###############################################################
Point and Vectors classes for Photometry (`astropack.geometry`)
###############################################################

Classes |XY| (point) and |DXY| (vector) are lightweight classes
for managing locations in a 2-dimensional plane, as for photometric apertures in an image.
Instances are subclasses of namedtuple, and so are themselves both tuples and namedtuples.

Classes |Rectangle_in_2D| and |Circle_in_2D| for defining regions in a 2-dimensional plane.
They are normally used for constructing aperture masks for background-corrected
aperture photometry.

*******************************************
Getting Started with |XY| and |DXY| classes
*******************************************

First, let's use class |XY|:

    >>> from astropack.geometry import XY, DXY
    >>> xy1 = XY(1, 5)
    >>> xy1.x, xy1.y
    (1, 5)
    >>> xy1
    XY(x=1, y=5)
    >>> tuple2 = (5, 9)
    >>> xy2 = XY.from_tuple(tuple2)
    >>> xy2
    XY(x=2, y=9)
    >>> xy1 == xy2
    False
    >>> xy1 != xy2
    True
    >>> xy2 - xy1  # vector from xy1 to xy2
    DXY(dx=1, dy=4)
    >>> xy1.vector_to(xy2)  # synonym for subtraction
    DXY(dx=1, dy=4)

Now, let's add in class |DXY|:

    >>> # ...continuing the above:
    >>> dxya = DXY(3,4)
    >>> dxyb = DXY.from_tuple((5, 7))  # tuple inside parentheses
    >>> dxya == dxyb
    False
    >>> dxya != dxyb
    True
    >>> dxya.length, dxya.length2
    (5.0, 25)
    >>> dxyb + dxya  # vectors DXY may be added (points XY cannot)
    DXY(dx=8, dy=11)
    >>> dxyb - dxya
    DXY(dx=2, dy=3)
    >>> dxya * 5
    DXY(dx=15, dy=20)
    >>> dxya * 5 == 5 * dxya
    True
    >>> dxya / 2
    DXY(dx=1.5, dy=2.0)
    >>> xy1 + dxya  # point xy1 displaced by dxya
    XY(x=4, y=9)
    >>> xy1 + 15 * dxya - 8 * dxyb
    XY(x=6, y=9)
    >>> dxya.direction         # result in radians
    0.9272952180016122
    >>> dxya.angle_with(dxyb)  # result in radians
    0.023251622810463002
    >>> dxya.dot(dxyb)
    43


***************
Reference/API
***************

.. automodapi:: astropack.geometry
   :no-inheritance-diagram:
