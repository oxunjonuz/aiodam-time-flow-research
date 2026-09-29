#!/usr/bin/env python3
"""Decisive PSD consistency check, FFT-binned (no analog-filter leakage).

For a segment x of length N, Parseval gives
    var_band = (2/N^2) * sum_{k in band} |X_k|^2
and with E|X_k|^2 = (fs*N/2) S_k this equals sum_{k in band} S_k * df.
So the FFT-binned band variance and the PSD integral MUST agree when both are
computed from the same data with the same estimator.  Any disagreement is either
estimator bias (median vs mean) or genuine non-Gaussian excess in the data.

The same statistic is computed on synthetic noise generated from the estimated
PSD, which isolates estimator bias from data excess.
"""
import math

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt

FS = 16384
FMIN, FMAX = 20.0, 300.0
PATHS = {"H1": "/work/time_flow_research_20260927/data/H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "/work/time_flow_research_20260927/data/L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}


def hp(x, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (FS / 2), btype="highpass", output="sos"), x)


def welch_all(x, nper, overlap=0.5):
    step = int(nper * (1 - overlap))
    w = np.hanning(nper)
    nrm = 2.0 / (FS * float(np.sum(w ** 2)))
    segs = np.array([np.abs(np.fft.rfft(x[i:i + nper] * w)) ** 2 * nrm
                     for i in range(0, len(x) - nper, step)])
    f = np.fft.rfftfreq(nper, 1.0 / FS)
    return f, segs


for name, path in PATHS.items():
    with h5py.File(path, "r") as h:
        strain = h["strain/Strain"][:]
    off = hp(np.concatenate([strain[:8 * FS], strain[-8 * FS:]]))
    N = len(off)
    f = np.fft.rfftfreq(N, 1.0 / FS)
    df = f[1] - f[0]
    band = (f >= FMIN) & (f <= FMAX)
    X = np.fft.rfft(off)
    var_fft = 2.0 / N ** 2 * float(np.sum(np.abs(X[band]) ** 2))

    fp, segs = welch_all(off, FS, 0.5)          # 1 s segments, 50% overlap
    bp = (fp >= FMIN) & (fp <= FMAX)
    Pmed = np.median(segs, axis=0) / math.log(2.0)
    Pmean = segs.mean(axis=0)
    int_med = float(np.sum(Pmed[bp]) * (fp[1] - fp[0]))
    int_mean = float(np.sum(Pmean[bp]) * (fp[1] - fp[0]))

    # synthetic control from the median PSD
    rng = np.random.default_rng(20260927)
    Sg = np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                          np.log(np.clip(fp, 1e-3, None)), np.log(Pmed)))
    Sg[0] = Sg[1]
    W = np.fft.rfft(rng.normal(0.0, 1.0, N)) * np.sqrt(Sg * FS / 2.0)
    syn = hp(np.fft.irfft(W, N))
    Xs = np.fft.rfft(syn)
    var_syn = 2.0 / N ** 2 * float(np.sum(np.abs(Xs[band]) ** 2))

    print(f"--- {name}  (nseg={len(segs)}, 1 s Hann 50%)")
    print(f"  band variance, FFT-binned, REAL      : {var_fft:.4e}")
    print(f"  band variance, FFT-binned, SYNTHETIC : {var_syn:.4e}")
    print(f"  PSD integral, median-Welch           : {int_med:.4e}")
    print(f"  PSD integral, mean-Welch             : {int_mean:.4e}")
    print(f"  K_real  = var_real / int_median      : {var_fft/int_med:.3f}")
    print(f"  K_synth = var_synth / int_median     : {var_syn/int_med:.3f}")
    print(f"  excess  = K_real / K_synth           : {var_fft/var_syn:.3f}")
    print(f"  mean/median PSD power ratio at 100 Hz: "
          f"{Pmean[np.argmin(abs(fp-100))]/Pmed[np.argmin(abs(fp-100))]:.1f}")
    print(f"  ASD median at 100 Hz: {math.sqrt(Pmed[np.argmin(abs(fp-100))]):.3e} /rtHz")
