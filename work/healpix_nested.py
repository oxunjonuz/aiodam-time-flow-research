#!/usr/bin/env python3
"""Self-contained HEALPix pixel -> (theta, phi), NESTED and RING, no healpy.

Written because the project's scripts run under the system python3, which has
neither healpy nor astropy_healpix.  Both geometries are validated against
healpy (in the local venv) by `verify_h3_joint.py`.

The FITS table of the GW150914 LALInference skymap stores PROB as a 2D array;
flattening it row-major and treating the result as NESTED gives a map whose
maximum sits at the WRONG sky position.  Which (ordering, geometry) pair is
correct is decided by PHYSICS in h3_joint_fit.py: only the correct one puts the
probability mass on the published 6.9 ms H1-L1 delay annulus.
"""
import math

import numpy as np


def nside2npix(nside):
    return 12 * nside * nside


def nest2ang(nside, pix):
    """Vectorised NESTED pixel index -> (theta, phi) in radians."""
    pix = np.asarray(pix, dtype=np.int64)
    n = int(round(math.log2(nside)))
    x = np.zeros_like(pix)
    y = np.zeros_like(pix)
    for i in range(n):
        x |= ((pix >> (2 * i)) & 1) << i
        y |= ((pix >> (2 * i + 1)) & 1) << i
    face = (pix >> (2 * n)) & 0b11
    jr = float(nside)
    z = np.empty(len(pix), dtype=float)
    phi = np.empty(len(pix), dtype=float)
    for f in (0, 1, 2, 3):
        m = face == f
        if not m.any():
            continue
        xx = x[m].astype(float)
        yy = y[m].astype(float)
        if f == 0:
            z[m] = 1.0 - (xx + 1.0) ** 2 / (3.0 * jr * jr)
            phi[m] = np.arctan2(yy + 0.5, xx + 0.5)
        elif f == 1:
            z[m] = -1.0 + (xx + 1.0) ** 2 / (3.0 * jr * jr)
            phi[m] = np.arctan2(yy + 0.5, xx + 0.5)
        elif f == 2:
            z[m] = (2.0 * yy + 1.0) / (3.0 * jr) - 1.0
            phi[m] = np.arctan2(xx + 0.5, yy + 0.5 - jr)
        else:
            z[m] = -(2.0 * yy + 1.0) / (3.0 * jr) + 1.0
            phi[m] = np.arctan2(xx + 0.5, yy + 0.5 - jr)
    f2 = face == 2
    f3 = face == 3
    phi[f2] = (phi[f2] + math.pi / 4.0) % (2.0 * math.pi)
    phi[f3] = (phi[f3] - math.pi / 4.0) % (2.0 * math.pi)
    theta = np.arccos(np.clip(z, -1.0, 1.0))
    return theta, phi


def ring2ang(nside, ipix):
    """Vectorised RING pixel index -> (theta, phi) in radians."""
    ipix = np.asarray(ipix, dtype=np.int64)
    npix = nside2npix(nside)
    ncap = 2 * nside * (nside - 1)
    z = np.empty(len(ipix), dtype=float)
    phi = np.empty(len(ipix), dtype=float)

    m = ipix < ncap
    if m.any():
        p = ipix[m].astype(float)
        iring = np.floor((1.0 + np.sqrt(1.0 + 2.0 * p)) / 2.0)
        iphi = p + 1.0 - 2.0 * iring * (iring - 1.0)
        z[m] = 1.0 - iring ** 2 / (3.0 * nside * nside)
        phi[m] = (iphi - 0.5) * math.pi / (2.0 * iring)

    m = (ipix >= ncap) & (ipix < npix - ncap)
    if m.any():
        p = ipix[m] - ncap
        iring = p // (4 * nside) + nside
        iphi = p % (4 * nside) + 1
        fodd = 0.5 * (1.0 + (iring + nside) % 2)
        z[m] = (2.0 * nside - iring) * 2.0 / (3.0 * nside)
        phi[m] = (iphi - fodd) * math.pi / (2.0 * nside)

    m = ipix >= npix - ncap
    if m.any():
        p = (npix - ipix[m]).astype(float)
        iring = np.floor((1.0 + np.sqrt(2.0 * p - 1.0)) / 2.0)
        iphi = 4.0 * iring + 1.0 - (p - 2.0 * iring * (iring - 1.0))
        z[m] = -(1.0 - iring ** 2 / (3.0 * nside * nside))
        phi[m] = (iphi - 0.5) * math.pi / (2.0 * iring)

    theta = np.arccos(np.clip(z, -1.0, 1.0))
    return theta, phi


def read_map_flat(path):
    """Read the PROB column of a HEALPix FITS map, flattened row-major."""
    from astropy.io import fits
    with fits.open(path) as h:
        data = np.asarray(h[1].data["PROB"]).astype(float)
    return data.reshape(-1)
