"""The astropack.geometry package provides classes and functions to create and
manipulate elementary geometric objects, especially useful for building
apertures and aperture masks.
"""

__author__ = "Eric Dose, Albuquerque"

# Python base:
from collections import namedtuple
import numbers
from math import sqrt, atan2, pi
from typing import TypeAlias, Union, Tuple

# External packages:
import numpy as np

__all__ = ['XY',
           'DXY',
           'Rectangle_in_2D',
           'Circle_in_2D',
           'distance_to_line',
           'make_golden_spiral']


XY_type: TypeAlias = 'XY'
DXY_type: TypeAlias = 'DXY'


__________CLASSES_____________________________________________________________ = 0

# import astropack.geometry


class XY:
    """Represents a Cartesian point (x,y) in a plane.

    |XY| instances are capable of a few operations not listed in the methods.

    Parameters
    ----------
    x : float
        X point location.
    y : float
        Y point location.


    Attributes
    ----------
    x : float
        X point location.
    y : float
        Y point location.
    """

    def __init__(self, x: float, y: float):
        self.x, self.y = x, y

    @classmethod
    def from_tuple(cls, xy_tuple: tuple[float, float]) -> XY_type:
        """Alternate constructor, from a tuple of 2 floats.

        Parameters
        ----------
        xy_tuple : tuple of 2 float
            Tuple (x,y) to be converted to new |XY| instance.

        Returns
        -------
        new_xy : |XY|
            new |XY| instance from ``xy_tuple``.
        """
        return XY(xy_tuple[0], xy_tuple[1])

    @property
    def as_tuple(self) -> Tuple[float, float]:
        """Return (x,y) as a tuple of floats.

        Returns
        -------
        xy_tuple : tuple of 2 float
            (x,y) as a tuple of 2 float numbers
        """
        return self.x, self.y

    def __eq__(self, other: XY_type) -> bool:
        """True if this XY instance is numerically equivalent to ``other``,
        else False."""
        if not isinstance(other, XY):
            return False
        return (self.x == other.x) and (self.y == other.y)

    def __ne__(self, other: XY_type) -> bool:
        """True if this XY instance is numerically unequal to ``other``,
        else False."""
        if not isinstance(other, XY):
            return False
        return not self == other

    def __add__(self, dxy: DXY_type) -> XY_type:
        """Add a displacement to this point.

        Parameters
        ----------
        dxy : |DXY|
            Distance by which to displace this point to yield new |XY| instance.

        Returns
        -------
        new_xy : |XY|
            new |XY| (point) which is displaced from this point by ``other``.
        """
        if isinstance(dxy, DXY):
            return XY(self.x + dxy.dx, self.y + dxy.dy)
        raise TypeError('XY.__add__() requires type DXY as operand.')

    def __sub__(self, other: Union[XY_type, DXY_type]) -> Union[XY_type, DXY_type]:
        """Subtract a displacement to give a new point, or
        subtract another point to give the displacement from this point to new point.

        Parameters
        ----------
        other : |DXY| or |XY|
            If |DXY| instance, the inverse of distance by which to displace
            this point to yield new point.
            If |XY| instance, the starting point of a vector ending in this point.

        Returns
        -------
        new_xy : |XY| or |DXY|
            If ``other`` is |DXY|, a new |XY| (point) which is
            displaced from this point by ``other``.
            If ``other`` is |XY|, a new |DXY| (vector) which
            is the displacement ``other`` to this point.
        """
        if isinstance(other, XY):
            return DXY(self.x - other.x, self.y - other.y)
        elif isinstance(other, DXY):
            return XY(self.x - other.dx, self.y - other.dy)
        raise TypeError('XY.__sub__() requires type XY or DXY as operand.')

    def vector_to(self, other: XY_type) -> DXY_type:
        """Returns |DXY| vector extending from this point to another point.

        Parameters
        ----------
        other : |XY|
            Another point, defining the end of the result vector.

        Returns
        -------
        new_vector : |DXY|
            Vector from this point to ``other``.
        """
        if isinstance(other, XY):
            return other - self
        raise TypeError('XY.vector_to() requires type XY as operand.')

    def __str__(self) -> str:
        return f'XY: {", ".join([str(self.x), str(self.y)])}'


class DXY:
    """Represents a vector (dx, dy) in a plane. A two-dimensional displacement.

    |DXY| instances are capable of a few operations not listed in the methods.

    Parameters
    ----------
    dx : float
        Vector magnitude in X direction.
    dy : float
        Vector magnitude in Y direction.

    Attributes
    ----------
    dx : float
        Vector magnitude in X direction.
    dy : float
        Vector magnitude in Y direction.
    """

    def __init__(self, dx: float, dy: float):
        self.dx, self.dy = dx, dy

    @classmethod
    def from_tuple(cls, dxy_tuple: tuple[float, float]) -> DXY_type:
        """Alternate constructor, from a tuple of 2 floats.

        Parameters
        ----------
        dxy_tuple : tuple of 2 float
            Tuple (dx,dy) to be converted to new |DXY| instance.

        Returns
        -------
        new_dxy : |DXY|
            new |DXY| instance from ``dxy_tuple``.
        """
        return DXY(dxy_tuple[0], dxy_tuple[1])

    @property
    def as_tuple(self) -> tuple[float, ...]:
        """Return (dx,dy) as a tuple of floats.

        Returns
        -------
        dxy_tuple : tuple of 2 float
            (dx,dy) as a tuple of 2 float numbers
        """
        return tuple([self.dx, self.dy])

    def __eq__(self, other: DXY_type):
        """True if this DXY instance is numerically equivalent to ``other``,
        else False."""
        if not isinstance(other, DXY):
            return False
        return self.dx == other.dx and self.dy == other.dy

    def __add__(self, other: Union[XY_type, DXY_type]) -> Union[XY_type, DXY_type]:
        """Adds two vectors to make a sum vector,
        or adds this vector to a point to make a new, displaced point.

        Parameters
        ----------
        other : |DXY| or |XY|
            Vector to add to this vector, or
            a point to which this vector shall be added.

        Returns
        -------
        result : |DXY| or |XY|
            Sum of vectors, or a new point displaced by this vector from point 'other'.
        """
        if isinstance(other, DXY):
            return DXY(self.dx + other.dx, self.dy + other.dy)
        elif isinstance(other, XY):
            return XY(other.x + self.dx, other.y + self.dy)
        raise TypeError('DXY.__add__() requires type DXY or XY as operand.')

    def __mul__(self, factor: float) -> DXY_type:
        """Multiplies this vector by a scalar factor.

        Parameters
        ----------
        factor : float
            Factor by which to multiply this vector. May be negative or zero.
            Handles syntax case ``dxy * factor`` (in that order).

        Returns
        -------
        new_displacement : |DXY|
            This vector multiplied by factor ``other``.
        """
        if isinstance(factor, numbers.Real):
            return DXY(factor * self.dx, factor * self.dy)
        raise TypeError('DXY.__mul__() requires float scalar as operand.')

    def __rmul__(self, factor: float) -> DXY_type:
        """Alias of __mul__(), to handle case of ``dxy * factor`` (in that order).

        Parameters
        ----------
        factor : float
            Factor by which to multiply this vector. May be negative or zero.
            Handles syntax case ``dxy * factor`` (in that order).

        Returns
        -------
        new_vector : |DXY|
            This vector multiplied by factor `other`.
        """
        if isinstance(factor, numbers.Real):
            return DXY(factor * self.dx, factor * self.dy)
        raise TypeError('DXY.__rmul__() requires float scalar as operand.')

    def __sub__(self, other: DXY_type) -> DXY_type:
        """Subtracts vector ``other`` from this one, yielding new vector.

        Parameters
        ----------
        other : |DXY|
            Vector to subtract from this vector.

        Returns
        -------
        new_vector : |DXY|
            New vector after subtracting ``other``.
        """
        if isinstance(other, DXY):
            return DXY(self.dx - other.dx, self.dy - other.dy)
        raise TypeError('DXY.__sub__() requires type DXY as operand.')

    def __truediv__(self, other: float) -> DXY_type:
        """Divides this vector by a scalar factor to yield a new vector.

        Parameters
        ----------
        other : float
            Scalar number by which to divide this vector.

        Returns
        -------
        new_displacement : |DXY|
            Vector divided by scalar factor ``other``.
        """
        if isinstance(other, numbers.Real):
            if other != 0:
                return DXY(self.dx / other, self.dy / other)
            raise ZeroDivisionError
        raise TypeError('DXY.__div__() requires scalar float as operand.')

    def angle_with(self, other: DXY_type) -> float:
        """Returns angle with other |DXY| vector, in radians, within range [0, 2 * pi].

        Parameters
        ----------
        other : |DXY|
            Another vector.

        Returns
        -------
        angle : float
            Angle between this vector and ``other``.
        """
        if isinstance(other, DXY):
            angle = other.direction - self.direction
            if angle < 0:
                angle += 2 * pi
            return angle
        raise TypeError('DXY.angle_with() requires type DXY as operand.')

    def dot(self, other: DXY_type) -> float:
        """Returns dot product of this vector with another |DXY| vector.

        Parameters
        ----------
        other : |DXY|
            Vector with which to determine dot product with this vector.

        Returns
        -------
        dot_procduct : float
            Dot product of this vector with ``other``.
        """
        if isinstance(other, DXY):
            return self.dx * other.dx + self.dy * other.dy
        raise TypeError('DXY.dot () requires type DXY as operand.')

    @property
    def length2(self) -> float:
        """Return square of vector length.

        Returns
        -------
        l2 : float
            Square of this vector's length.
        """
        return self.dx ** 2 + self.dy ** 2

    @property
    def length(self) -> float:
        """Return vector length.

        Returns
        -------
        l : float
            This vector's length.
        """
        return sqrt(self.length2)

    @property
    def direction(self) -> float:
        """Return angle of this vector relative to positive x-axis
        (e.g., +y yields +pi/2), in radians. Returns zero if vector is zero-length.

        Returns
        -------
        direction : float
            Angle of this vector relative to positive x-axis, in radians.
            In range [0, 2 * pi].
        """
        if self.length2 > 0.0:
            return atan2(self.dy, self.dx)
        return 0.0

    def __str__(self) -> str:
        return f'DXY: {", ".join([str(self.dx), str(self.dy)])}'


class Rectangle_in_2D:
    """Represents a rectangle of any orientation within a 2-D Cartesian plane,
    able to test whether (x,y) point(s) are inside the rectangle, and with facility
    to construct a 2-D mask array in the shape of a rectangle of arbitrary location,
    shape, size, and orientation.

    Parameters
    ----------
    xy_a, xy_b, xy_c: each |XY|
        Each a vertex of the new rectangle.
        Vertices must be adjacent, clockwise or counter-clockwise.

    Attributes
    ----------
    area : float
        Area of this rectangle.
    """

    def __init__(self, xy_a, xy_b, xy_c):
        self.a, self.b, self.c = (xy_a, xy_b, xy_c)
        self.ab, self.bc = (self.b - self.a, self.c - self.b)
        if abs(self.ab.dot(self.bc)) > 1e-10 * min(self.ab.length, self.bc.length):
            raise ValueError('Rectangle_in_2D edges are not perpendicular.')
        self.area = self.ab.length * self.bc.length

    def contains_point(self, xy: XY_type, include_edges: bool = True) -> bool:
        """Returns True if this rectangle contains point ``xy``, else return False.

        Parameters
        ----------
        xy : |XY|
            Point to determine whether inside or outside rectangle.

        include_edges : bool, optional
            If and only if ``include_edges`` is True, a point on the rectangle's
            boundaries will be considered contained. Default is True.

        Returns
        -------
        is_contained : bool
            True if point ``xy`` is contained within this rectangle.
        """
        dot_ab = self.ab.length2  # reflexive dot product.
        dot_bc = self.bc.length2  # "
        dot_ab_pt = self.ab.dot(xy - self.a)
        dot_bc_pt = self.bc.dot(xy - self.b)
        if include_edges:
            return (0 <= dot_ab_pt <= dot_ab) and (0 <= dot_bc_pt <= dot_bc)
        return (0 < dot_ab_pt < dot_ab) and (0 < dot_bc_pt < dot_bc)

    def contains_points(self, xy_array: list[XY_type] | np.ndarray[XY_type],
                        include_edges: bool = True) -> bool | list[bool]:
        """Returns True for each corresponding point in ``xy_array`` if this rectangle
         contains that point, else return False. Array version of `contains_point()`.

        Parameters
        ----------
        xy_array : array or list of |XY|
            Points to determine whether inside or outside rectangle.

        include_edges : bool, optional
            If and only if ``include_edges`` is True, points on the rectangle's
            boundaries will be considered contained. Default is True.

        Returns
        -------
        are_contained : array or list of bool
            Boolean values, each True if corresponding point in ``xy_array``
            is contained within this rectangle.
        """
        dot_ab = self.ab.length2  # reflexive dot product, compute only once for array.
        dot_bc = self.bc.length2  # "
        dot_ab_array = [self.ab.dot(pt - self.a) for pt in xy_array]
        dot_bc_array = [self.bc.dot(pt - self.b) for pt in xy_array]
        if include_edges:
            return[(0 <= dot_a_pt <= dot_ab) and (0 <= dot_b_pt <= dot_bc)
                   for (dot_a_pt, dot_b_pt) in zip(dot_ab_array, dot_bc_array)]
        return[(0 < dot_a_pt < dot_ab) and (0 < dot_b_pt < dot_bc)
               for (dot_a_pt, dot_b_pt) in zip(dot_ab_array, dot_bc_array)]

    def contains_points_unitgrid(self, x_min: int, x_max: int,
                                 y_min: int, y_max: int,
                                 include_edges: bool = True) -> np.ndarray:
        """Constructs a unit grid in X and Y, then sets each point in the grid to
        True if and only if that point's coordinates fall within this rectangle,
        else sets the point to False.
        All 4 coordinate input values must be integers.
        The implementation is numpy-optimized, but
        the return array indices are as (x, y), NOT in numpy index convention (y, x).

        Parameters
        ----------
        x_min, x_max: each int
            Minimum and maximum X values in unit grid.
        y_min, y_max: each int
            Minimum and maximum Y values in unit grid.
        include_edges : bool, optional
            If and only if ``include_edges`` is True, grid points lying exactly on the
            rectangle's boundaries will be considered contained. Default is True.

        Returns
        -------
        grid : 2-dimensional |ndarray|
            Unit grid in X and Y, with each point in the grid set to
            True if and only if that point's coordinates fall within this rectangle,
            else False.
        """
        if any(not isinstance(val, int) for val in (x_min, x_max, y_min, y_max)):
            raise TypeError('All 4 grid edges must be integers but are not.')
        if x_max < x_min or y_max < y_min:
            raise ValueError('Grid point specifications must be in ascending order.')
        dot_ab = self.ab.length2  # reflexive dot product, precompute once for array.
        dot_bc = self.bc.length2  # "
        # NB: indexing='ij' actually gives (x,y) indices:
        x, y = np.meshgrid(range(x_min, x_max + 1), range(y_min, y_max + 1),
                           indexing='ij')
        # Grids (numpy arrays) of dot products:
        dot_ab_xy = self.ab.dx * (x - self.a.x) + self.ab.dy * (y - self.a.y)
        dot_bc_xy = self.bc.dx * (x - self.b.x) + self.bc.dy * (y - self.b.y)
        if include_edges:
            return (0 <= dot_ab_xy) & (dot_ab_xy <= dot_ab) & \
                   (0 <= dot_bc_xy) & (dot_bc_xy <= dot_bc)
        return (0 < dot_ab_xy) & (dot_ab_xy < dot_ab) & \
               (0 < dot_bc_xy) & (dot_bc_xy < dot_bc)


class Circle_in_2D:
    """Represents a circle in a 2-D Cartesian plane,
    able to test whether (x,y) point(s) are inside the circle, and with facility
    to construct a 2-D mask array in the shape of a circle of arbitrary location and
    size.

    Parameters
    ----------

    xy_origin : |XY|
        Location of circle's center.

    radius : float
        Radius of circle.

    Attributes
    ----------

    area : float
        Area of circle.
    """

    def __init__(self, xy_origin: XY, radius: float):
        self.origin = xy_origin
        self.radius = radius
        if not isinstance(self.origin, XY):
            raise TypeError('Circle origin must be a XY object')
        self.x, self.y = (self.origin.x, self.origin.y)
        self.area = pi * (self.radius ** 2)

    def contains_point(self, xy: XY_type, include_edges: bool = True) -> bool:
        """Returns True if this circle contains point ``xy``, else return False.

        Parameters
        ----------
        xy : |XY|
            Point to determine whether inside or outside this circle.

        include_edges : bool, optional
            If and only if ``include_edges`` is True, a point on the circle's
            boundaries will be considered contained. Default is True.

        Returns
        -------
        is_contained : bool
            True if point ``xy`` is contained within this circle.
        """
        distance2 = (xy - self.origin).length2
        if include_edges:
            return distance2 <= self.radius ** 2
        return distance2 < self.radius ** 2

    def contains_points(self, xy_array: np.ndarray[XY_type] | list[XY_type],
                        include_edges=True) -> np.ndarray[bool] | list[bool]:
        """Returns True for each corresponding point in ``xy_array`` if this circle
         contains that point, else return False. Array version of ``contains_point()``.

        Parameters
        ----------
        xy_array : array or list of |XY|
            Points to determine whether inside or outside this circle.

        include_edges : bool, optional
            If and only if ``include_edges`` is True, points on this circle's
            boundaries will be considered contained. Default is True.

        Returns
        -------
        are_contained : array or list of bool
            Boolean values, each True if corresponding point in ``xy_array`` is
            contained within this circle.
        """
        distances2 = [(xy - self.origin).length2 for xy in xy_array]
        if include_edges:
            return [distance2 <= self.radius ** 2 for distance2 in distances2]
        return [distance2 < self.radius ** 2 for distance2 in distances2]

    def contains_points_unitgrid(self, x_min: int, x_max: int,
                                 y_min: int, y_max: int,
                                 include_edges=True) -> np.ndarray:
        """Constructs a unit grid in X and Y, then sets each point in the grid to
        True if and only if that point's coordinates fall within this circle,
        else sets the point to False.
        All 4 coordinate input values must be integers.
        The implementation is numpy-optimized, but
        the return array indices are as (x, y), NOT in numpy index convention (y, x).

        Parameters
        ----------
        x_min, x_max: each int
            Minimum and maximum X values in unit grid.
        y_min, y_max: each int
            Minimum and maximum Y values in unit grid.
        include_edges : bool, optional
            If and only if ``include_edges`` is True, grid points lying exactly
            on the circle's boundaries will be considered contained. Default is True.

        Returns
        -------
        grid : 2-dimensional |ndarray|
            Unit grid in X and Y, with each point in the grid set to
            True if and only if that point's coordinates fall within this circle,
            else False.
        """
        if any(not isinstance(val, int) for val in (x_min, x_max, y_min, y_max)):
            raise TypeError('All 4 grid edges must be integers but are not.')
        if x_max < x_min or y_max < y_min:
            raise ValueError('Grid point specifications must be in ascending order.')
        dx2_values = [(x - self.origin.x) ** 2 for x in range(x_min, x_max + 1)]
        dy2_values = [(y - self.origin.y) ** 2 for y in range(y_min, y_max + 1)]
        dx2, dy2 = np.meshgrid(dx2_values, dy2_values, indexing='ij')
        if include_edges:
            return (dx2 + dy2) <= self.radius ** 2
        return (dx2 + dy2) < self.radius ** 2


__________FUNCTIONS_____________________________________________________________ = 0


def distance_to_line(xy_pt: Union[XY_type, tuple[float, float]],
                     xy_1: Union[XY_type, tuple[float, float]],
                     xy_2: Union[XY_type, tuple[float, float]],
                     dist_12: float | None = None) -> float:
    """Shortest distance from a point to a line defined by two non-coincident points.
    Line is consdered to extend infinitely in both directions (that is, it is a
    formal line, not a line *segment*).

    Parameters
    ----------
    xy_pt : |XY|, or tuple of 2 float
        The point whose distance to measure.

    xy_1, xy_2 : each  |XY|, or tuple of 2 float
        Two distinct points lying on the line.

    dist_12 : float or None, optional
        Square of distance between ``xy_1`` and ``xy_2`` if already computed
        and available. Otherwise, return None (the commonest case), which indicates
        that this function needs to compute ``dist_12``. Default is None

    Returns
    -------
    distance : float
        Shortest distance between ``xy_pt`` and the line defined by points
        ``xy_1`` and ``xy_2``.
    """
    if any([isinstance(xy, (XY, tuple)) is False for xy in [xy_pt, xy_1, xy_2]]):
        raise TypeError('Parameters \'xy_pt\', \'xy_1\', and \'xy_2\' '
                        'must each be 2-tuple or XY instance.')
    xpt, ypt = xy_pt if isinstance(xy_pt, tuple) else xy_pt.as_tuple
    x1, y1 = xy_1 if isinstance(xy_1, tuple) else xy_1.as_tuple
    x2, y2 = xy_2 if isinstance(xy_2, tuple) else xy_2.as_tuple

    if (x1, y1) == (x2, y2):
        return sqrt((xpt - x1)**2 + (ypt - y1)**2)
    if dist_12 is None:
        dist_12 = sqrt((y2 - y1)**2 + (x2 - x1)**2)
    distance = abs((x2 - x1) * (y1 - ypt) - (x1 - xpt) * (y2 - y1)) / dist_12
    return distance


def make_golden_spiral(n_points: int) -> list[namedtuple]:
    """Return points evenly spaced on a sphere.

.. note:: This function may be removed in favor of recommending astropy's
          function
          `golden_spiral_grid() <https://docs.astropy.org/en/stable/
          api/astropy.coordinates.golden_spiral_grid.html>`_

    Parameters
    ----------
    n_points : int
        Number of points wanted, which points will be distributed
        over the entire sphere.

    Returns
    -------
    az_alt_list : list of namedtuple 'azalt', each of 2 float
        Each namedtuple represents one (azimuth, altitude) sky position and has
        members ``.az`` and ``.alt``.
    """
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)  # altitude in radians
    theta = pi * (1 + 5**0.5) * indices  # azimuth in radians
    altitudes = (pi / 2.0 - phi) * (180.0 / pi)
    azimuths = (theta * 180.0 / pi) % 360.0
    azalt = namedtuple('azalt', ['az', 'alt'])
    azalt_list = [azalt(az, alt) for (az, alt) in zip(azimuths, altitudes)]
    return azalt_list
