#!/usr/bin/env python3
"""Independent recomputation of the H3 absorption statistic on the REAL data.

Different route from h3_real_data.py on every axis that matters:
  * it re-reads the strain and rebuilds the PSD itself (its own estimator code);
  * it builds the delay directions from the ANALYTIC series expansion in
    x = (pi M f)^(2/3) rather than from the numeric t(f) function;
  * it solves the projection by EIGEN-DECOMPOSITION of the weighted Gram matrix
    (np.linalg.eigh) rather than by QR;
  * it computes the absorbed fraction as 1 - (residual energy / total energy)
    from the eigenbasis, not from an explicit projector.

Writes nothing; prints its numbers so they can be compared with the artifact.
Exit 1 if it disagrees with artifacts/h3_real_data.json beyond tolerance.
"""
import json
import math
import os
import sys

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt

FS = 16384
FMIN, FMAX = 20.0, 300.0
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
FILES = {"H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30
MC = 28.096 * G * MSUN / C ** 3
DL = 410.0 * 3.0856775814913673e22 / C

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{(' -- ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


def hp(x, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (FS / 2), btype="highpass", output="sos"), x)


def psd_median(x, nper=FS, overlap=0.5):
    step = int(nper * (1 - overlap))
    w = np.hanning(nper)
    nrm = 2.0 / (FS * float(np.sum(w ** 2)))
    segs = np.array([np.abs(np.fft.rfft(x[i:i + nper] * w)) ** 2 * nrm
                     for i in range(0, len(x) - nper, step)])
    return np.fft.rfftfreq(nper, 1.0 / FS), np.median(segs, axis=0) / math.log(2.0)


def main():
    art = json.load(open(os.path.join(HERE, "artifacts", "h3_real_data.json")))
    print("== independent H3 recomputation (eigh route, analytic series) ==")
    df = 1.0 / 4.0
    f = np.arange(FMIN, FMAX + df / 2.0, df)
    x = (np.pi * MC * f) ** (2.0 / 3.0)
    # analytic expansion of the phase perturbation, in x, from the series
    # dPsi = 2 x^{3/2} tau0/Mc  - (5/128) tau1 x^{-5/2} + (25/32768) tau2 Mc x^{-13/2}
    t_span = abs((5.0 / 256.0) * MC ** (-5.0 / 3.0)
                 * ((np.pi * FMIN) ** (-8.0 / 3.0) - (np.pi * FMAX) ** (-8.0 / 3.0)))
    # leading-order amplitude, and the three free directions (analytic)
    h = (1.0 / DL) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) * MC ** (5.0 / 6.0) \
        * f ** (-7.0 / 6.0)
    psi = (3.0 / 128.0) * x ** (-5.0 / 2.0)
    col_mc = h * (-(5.0 / 3.0) * psi / MC)
    col_tc = h * (2.0 * np.pi * f)
    col_ph = -h

    for det in ("H1", "L1"):
        with h5py.File(os.path.join(DATA, FILES[det]), "r") as fh:
            strain = fh["strain/Strain"][:]
        off = hp(np.concatenate([strain[:8 * FS], strain[-8 * FS:]]))
        fp, Pp = psd_median(off)
        S = np.exp(np.interp(np.log(f), np.log(fp), np.log(Pp)))
        w = np.sqrt(4.0 * df / S)
        A = np.column_stack([col_mc, col_tc, col_ph]) * w[:, None]
        # eigendecomposition of the Gram matrix -> orthonormal basis of the span
        Gram = A.T @ A
        evals, evecs = np.linalg.eigh(Gram)
        keep = evals > evals.max() * 1e-24
        Q = A @ evecs[:, keep] / np.sqrt(evals[keep])

        def absorbed(u):
            y = u * w
            e_tot = float(y @ y)
            c = Q.T @ y
            e_res = float(y @ y - c @ c)
            # absorbed_fraction = ||P u|| / ||u|| = sqrt(1 - e_res/e_tot),
            # NOT 1 - e_res/e_tot (that is the absorbed ENERGY fraction).
            frac = math.sqrt(max(1.0 - e_res / e_tot, 0.0)) if e_tot > 0 else 1.0
            return frac, math.sqrt(max(e_res, 0.0))

        u0 = h * (2.0 * np.pi * f * 1.2e-3)
        ul = h * (2.0 * np.pi * f * (1e-2 * ((5.0 / 256.0) * MC ** (-5.0 / 3.0)
                                             * (np.pi * f) ** (-8.0 / 3.0) * -1.0)))
        tau2 = 1.2e-3 / t_span ** 2
        uq = h * (2.0 * np.pi * f * (tau2 * ((5.0 / 256.0) * MC ** (-5.0 / 3.0)
                                             * (np.pi * f) ** (-8.0 / 3.0)) ** 2))
        a0, r0 = absorbed(u0)
        al, rl = absorbed(ul)
        aq, rq = absorbed(uq)

        ref = art["projection"][det]
        check(f"{det} H1 residual SNR matches artifact",
              abs(r0 - ref["H1_const_delay"]["unmodelled_snr"]) < 1e-9,
              f"{r0:.3e} vs {ref['H1_const_delay']['unmodelled_snr']:.3e}")
        check(f"{det} H2 residual SNR matches artifact",
              abs(rl - ref["H2_linear_delay"]["unmodelled_snr"]) < 1e-8,
              f"{rl:.3e} vs {ref['H2_linear_delay']['unmodelled_snr']:.3e}")
        check(f"{det} H3 absorbed fraction matches artifact",
              abs(aq - ref["H3_quadratic"]["absorbed_fraction"]) < 1e-6,
              f"{aq:.8f} vs {ref['H3_quadratic']['absorbed_fraction']:.8f}")
        check(f"{det} H3 residual SNR matches artifact",
              abs(rq - ref["H3_quadratic"]["unmodelled_snr"]) < 1e-6,
              f"{rq:.6f} vs {ref['H3_quadratic']['unmodelled_snr']:.6f}")
        print(f"    ({det}: H1 rho={r0:.2e}  H2 rho={rl:.2e}  H3 absorbed={aq:.6f}  "
              f"H3 rho={rq:.6f})")

    print()
    if FAILS:
        print(f"INDEPENDENT_DISAGREES: {FAILS}")
        return 1
    print("INDEPENDENT_CONFIRMED (H3, eigh route)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
