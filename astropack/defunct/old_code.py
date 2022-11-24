__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os

# External packages:

# Author's packages:


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')

# This function is not needed, as:
# (1) it yields apparent Dec only, and
# (2) skyfield and astropy can do this for J2000.0 etc.
#
# def ha_dec_from_az_alt(latitude: float,
#                        az_alt: [(float, float), [(float, float)]]) \
#                        -> [tuple[float, float], list[tuple[float, float]]]:
#     """Return hour angle and declination (APPARENT, =current epoch! not J2000.0 etc.)
#     for a given azimuth and altitude, at a site's given latitude.
#
#      Parameters
#      ----------
#      latitude : float
#         Observing site latitude, in degrees.
#
#      az_alt : tuple of 2 float, or list of tuple
#         Tuple of (azimuth, altitude), in degrees.
#
#     Returns
#     -------
#     ha_dec : tuple of 2 float, or list of tuple
#        Tuple of (hourangle, declination), in degrees.
#     """
#     az_alt_is_list = isinstance(az_alt, list)
#     if not az_alt_is_list:
#         az_alt = [az_alt]
#
#     # Formulae (1) and (2) from http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm:
#     # phi=latitude, delta=declination, H=hourangle, A=target azimuth, a=target altitude.
#     cos_phi, sin_phi = cos(latitude / DEGREES_PER_RADIAN), \
#         sin(latitude / DEGREES_PER_RADIAN)
#     results = []
#     for (az, alt) in az_alt:
#         cos_alt, sin_alt = cos(alt / DEGREES_PER_RADIAN), sin(alt / DEGREES_PER_RADIAN)
#         cos_az = cos(az / DEGREES_PER_RADIAN)
#
#         # (1) sin(δ) = sin(a) sin(φ) + cos(a) cos(φ) cos(A)
#         sin_delta = (sin_alt * sin_phi) + (cos_alt * cos_phi * cos_az)
#         cos_delta = sqrt(1.0 - sin_delta ** 2)  # cosine of declination >= 0 always.
#         # Clip to correct any floating-point errors driving abs(sin(delta)) > 1:
#         sin_delta = clip(sin_delta, -1.0, 1.0)
#         delta = asin(sin_delta) * DEGREES_PER_RADIAN
#
#         # (2) cos(H) = { sin(a) - sin(δ) sin(φ) } / { cos(δ) cos(φ) }
#         cos_h = (sin_alt - sin_delta * sin_phi) / (cos_delta * cos_phi)
#         # Clip to correct any floating-point errors driving abs(cos(H)) > 1:
#         cos_h = clip(cos_h, -1.0, 1.0)
#         # print(f'cos_h={cos_h}')
#         sign_h = -1.0 if 0.0 <= az <= 180.0 else 1.0  # because sign of H is ambiguous.
#         h = sign_h * acos(cos_h) * DEGREES_PER_RADIAN
#         # print(f'sign_h={sign_h}, acos={acos(cos_h)}')
#         results.append((h, delta))
#     if not az_alt_is_list:
#         results = results[0]
#     return results
