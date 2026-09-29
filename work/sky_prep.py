#!/usr/bin/env python3
"""Prepare the GW150914 sky samples ONCE, and decide the correct reading by PHYSICS.

WHY THIS FILE EXISTS.

The project's pipeline scripts run under the system python3, which has neither
healpy nor astropy_healpix.  My first two attempts to hand-roll a NESTED pixel
-> (theta, phi) reader were BOTH WRONG (checked against healpy: pixel 0 came out
at theta = 0.09 deg instead of 89.93 deg).  Re-implementing HEALPix is exactly
the class of error this project keeps catching, so it is not done: the sky is
read here with validated libraries and frozen to a .npz the fit reads.

THE TRAP, AND HOW IT WAS SETTLED.

The FITS table stores PROB as a 2D array (3072 x 1024).  Reading it with
`data['PROB'].reshape(-1)` (row-major) and `healpy.read_map` gives TWO DIFFERENT
arrays, and they put the probability mass on DIFFERENT parts of the sky.  Neither
library is trusted on its word: the reading is chosen by a PHYSICAL test.

The published constraint (arXiv:1602.03840) is that the source lies on the
annulus of constant H1-L1 arrival-time difference 6.9 (+0.5/-0.4) ms.  A correct
sky map must therefore concentrate its probability on that annulus.  Measured,
with a shuffled-map NULL as the control:

    reading                    peak concentration   at delay
    row-major + NESTED              0.1635          -6.90 ms   <-- accepted
    healpy.read_map + NESTED        0.0233          +3.70 ms
    shuffled map (20 draws)         0.0061 mean, 0.0064 max

The accepted reading is 25x the null and lands on the published 6.9 ms (the sign
is consistent: t_H1 - t_L1 = -(r_H1 - r_L1).n/c = +6.9 ms when the projection is
-6.9 ms, i.e. the signal reached L1 first).  Two INDEPENDENT geometry
implementations -- healpy.pix2ang(nest=True) and astropy_healpix.HEALPix -- must
agree on the sky position of every pixel, and they do to 1e-15.

WHAT IS WRITTEN
  artifacts/sky_samples.npz: ra, dec (radians) of the 90% region subsampled,
  their weights, and the scalar summary (ml position, areas, delay peak).
"""
import json
import math
import os
import sys

import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
SKYMAP = os.path.join(HERE, "sky", "LALInference_skymap.fits")
OUT = os.path.join(HERE, "artifacts", "sky_samples.npz")
NSIDE = 512

CL = 2.99792458e8
VERTEX = {
    "H1": np.array([-2.16141492636e6, -3.83469517889e6, 4.60035022664e6]),
    "L1": np.array([-7.42760447238e4, -5.49628371971e6, 3.22425701744e6]),
}
PUBLISHED_DELAY_MS = 6.9
PUBLISHED_AREA_90 = 610.0


def read_flat(path):
    with fits.open(path) as h:
        return np.asarray(h[1].data["PROB"]).astype(float).reshape(-1)


def geometry_healpy(nside):
    import healpy as hp
    th, ph = hp.pix2ang(nside, np.arange(12 * nside * nside), nest=True)
    return ph, np.pi / 2.0 - th


def geometry_astropy(nside):
    from astropy_healpix import HEALPix
    import astropy.units as u
    hpix = HEALPix(nside=nside, order="nested", frame="icrs")
    lon, lat = hpix.healpix_to_lonlat(np.arange(12 * nside * nside))
    return lon.to_value(u.rad), lat.to_value(u.rad)


def gmst_rad(gps):
    """Greenwich mean sidereal time (IAU-1982 family, as LAL uses).

    The skymap is in CELESTIAL (equatorial) coordinates; the detector vertices
    are in the EARTH-FIXED frame.  Projecting one onto the other without this
    rotation puts the probability on the wrong annulus -- measured: without it
    the peak sits at -3.4 ms with concentration 0.026, with it at -6.9 ms with
    concentration 0.164.
    """
    unix = 315964800.0 + gps - 18.0
    jd = unix / 86400.0 + 2440587.5
    T = (jd - 2451545.0) / 36525.0
    g = (280.46061837 + 360.98564736629 * (jd - 2451545.0)
         + 0.000387933 * T * T - T * T * T / 38710000.0)
    return math.radians(g % 360.0)


def delay_ms(ra, dec, gps=1126259462.4):
    n = np.stack([np.cos(dec) * np.cos(ra),
                  np.cos(dec) * np.sin(ra),
                  np.sin(dec)], axis=1)
    g = gmst_rad(gps)
    c, s = math.cos(g), math.sin(g)
    R = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    n_earth = (R @ n.T).T
    return (VERTEX["H1"] - VERTEX["L1"]) @ n_earth.T / CL * 1e3


def concentration(prob, proj):
    edges = np.arange(-10.05, 10.05, 0.1)
    hist, _ = np.histogram(proj, bins=edges, weights=prob)
    hh = hist / hist.sum()
    i = int(np.argmax(hh))
    return float(hh[i]), float(edges[i] + 0.05)


def main():
    out = {"skymap": os.path.relpath(SKYMAP, HERE), "nside": NSIDE}
    prob = read_flat(SKYMAP)
    nside = NSIDE
    npix = 12 * nside * nside
    assert len(prob) == npix, (len(prob), npix)
    prob = prob / prob.sum()

    ra_h, dec_h = geometry_healpy(nside)
    ra_a, dec_a = geometry_astropy(nside)
    geo_agree = float(max(np.abs(ra_h - ra_a).max(), np.abs(dec_h - dec_a).max()))
    out["geometry_agreement_rad"] = geo_agree

    proj = delay_ms(ra_h, dec_h)
    conc, delay_peak = concentration(prob, proj)

    # shuffled-map null
    rng = np.random.default_rng(12345)
    null = [concentration(prob[rng.permutation(npix)], proj)[0] for _ in range(20)]
    null_mean, null_max = float(np.mean(null)), float(np.max(null))

    k = int(np.argmax(prob))
    order = np.argsort(prob)[::-1]
    cum = np.cumsum(prob[order])
    pixarea = 4.0 * math.pi / npix * (180.0 / math.pi) ** 2
    area90 = float((int(np.searchsorted(cum, 0.9)) + 1) * pixarea)
    area50 = float((int(np.searchsorted(cum, 0.5)) + 1) * pixarea)

    out.update({
        "ml_pixel": k,
        "ml_ra_deg": math.degrees(float(ra_h[k])),
        "ml_dec_deg": math.degrees(float(dec_h[k])),
        "area90_deg2": area90, "area50_deg2": area50,
        "delay_annulus_peak_ms": delay_peak,
        "delay_concentration": conc,
        "null_concentration_mean": null_mean,
        "null_concentration_max": null_max,
        "concentration_over_null": conc / null_max,
        "published_delay_ms": PUBLISHED_DELAY_MS,
        "published_area90_deg2": PUBLISHED_AREA_90,
        "sign_note": "proj = (r_H1-r_L1).n/c is -6.9 ms; the arrival-time "
                     "difference t_H1-t_L1 = -proj = +6.9 ms, i.e. L1 first, "
                     "which is the published statement.",
    })

    ok = (geo_agree < 1e-9
          and abs(delay_peak + PUBLISHED_DELAY_MS) < 0.6
          and conc > 5.0 * null_max
          and abs(area90 - PUBLISHED_AREA_90) / PUBLISHED_AREA_90 < 0.05)
    out["accepted"] = bool(ok)

    n90 = int(np.searchsorted(cum, 0.9)) + 1
    sel = order[:n90]
    step = max(1, len(sel) // 400)
    sel = sel[::step]
    w = prob[sel] / prob[sel].sum()
    out["n_sky_samples"] = int(len(sel))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, ra=ra_h[sel], dec=dec_h[sel], weight=w,
             ml_ra=np.array([ra_h[k]]), ml_dec=np.array([dec_h[k]]),
             area90_deg2=np.array([area90]), area50_deg2=np.array([area50]),
             delay_peak_ms=np.array([delay_peak]), concentration=np.array([conc]))
    out["written"] = os.path.relpath(OUT, HERE)
    print(json.dumps(out, indent=2, sort_keys=True, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
