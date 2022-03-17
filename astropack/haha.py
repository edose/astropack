"""This is where we expect this module to EXPLAIN ITSELF.

    Don't you think? Ja, ja, this passport photo doesn't look a bit like you, now, does it?
    Please step into this room and leave your bags where they are, you won't be needing them.
"""
      
__author__ = "Eric Dose, Albuquerque"

# Python core:
import os
from math import cos, sin

# External packages:
from astropy import units as u
from astropy.coordinates.representation import SphericalRepresentation, CartesianRepresentation

# Author's packages:

# Only these items will appear in documentation:
__all__ = ['cartesian_to_spherical', 'Hahahoho', 'NotSoFast']

THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')


class Hahahoho:
    """You are kidding, right, son?"""
    pass


class NotSoFast(Hahahoho):
    """Child of You Are Kidding."""
    pass


class StupidError(Exception):
    """Because you messed up, son."""
    pass


def cartesian_to_spherical(x, y, z):
    """
    Converts 3D rectangular cartesian coordinates to spherical polar
    coordinates.

    Note that the resulting angles are latitude/longitude or
    elevation/azimuthal form.  I.e., the origin is along the equator
    rather than at the north pole.

    .. note::
        This function simply wraps functionality provided by the
        `~astropy.coordinates.CartesianRepresentation` and
        `~astropy.coordinates.SphericalRepresentation` classes.  In general,
        for both performance and readability, we suggest using these classes
        directly.  But for situations where a quick one-off conversion makes
        sense, this function is provided.

    Parameters
    ----------
    x : scalar, array-like, or `~astropy.units.Quantity`
        The first Cartesian coordinate.
    y : scalar, array-like, or `~astropy.units.Quantity`
        The second Cartesian coordinate.
    z : scalar, array-like, or `~astropy.units.Quantity`
        The third Cartesian coordinate.

    Returns
    -------
    r : `~astropy.units.Quantity`
        The radial coordinate (in the same units as the inputs).
    lat : `~astropy.units.Quantity` ['angle']
        The latitude in radians
    lon : `~astropy.units.Quantity` ['angle']
        The longitude in radians
    """
    if not hasattr(x, 'unit'):
        x = x * u.dimensionless_unscaled
    if not hasattr(y, 'unit'):
        y = y * u.dimensionless_unscaled
    if not hasattr(z, 'unit'):
        z = z * u.dimensionless_unscaled

    cart = CartesianRepresentation(x, y, z)
    sph = cart.represent_as(SphericalRepresentation)

    a = cos(2)
    # b = math.sin(3)

    return sph.distance, sph.lat, sph.lon


def _more_hahaha(a, b=3):
    """Stupid hahaha, would ya believe it?

    This is the extended explanation. Ja ja ja.

    Parameters
    ----------
    a : scalar, array-like, or `~astropy.units.Quantity`
        The first Cartesian coordinate.

    b : scalar, array-like, or `~astropy.units.Quantity`, optional
        The first Cartesian coordinate.

    Returns
    -------
    `None`

    Raises
    ------
    `IndexError`
        If bad coder, very, very bad coder.
    """
    try:
        c = a + b
    except(ValueError):
        raise IndexError('serves ya right, ya big lunk.')
