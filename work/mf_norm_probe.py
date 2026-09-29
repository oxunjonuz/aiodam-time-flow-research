#!/usr/bin/env python3
"""Pin the matched-filter amplitude convention and the edge treatment.

Facts to establish, each by measurement:
  A. which template convention makes sigma equal the standard optimal SNR
     4*df*sum(|h~|^2/S)  (analytic value printed for reference)
  B. with the data windowed, does the real-data peak land on the merger time?
  C. what does the off-source control give (same filter, quiet segment)?
  D. does an SNR-20 injection come back as 20 at the injected time?
"""
import math

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

FS = 16384
SEG_N = 4 * FS
FMIN, FMAX = 20.0, 300.0
GPS0 = 1126259447
T_EVENT = 1126259462.4 - GPS0
PATH = "/work/time_flow_research_20260927/data/H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5"


def highpass(x, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (FS / 2), btype="highpass", output="sos"), x)


def median_welch(x, nper):
    w = np.hanning(nper)
    nrm = 2.0 / (FS * float(np.sum(w ** 2)))
    n = len(x) // nper
    acc = np.empty((n, nper // 2 + 1))
    for i in range(n):
        acc[i] = np.abs(np.fft.rfft(x[i * nper:(i + 1) * nper] * w)) ** 2 * nrm
    return np.fft.rfftfreq(nper, 1.0 / FS), np.median(acc, axis=0) / math.log(2.0)


with h5py.File(PATH, "r") as h:
    strain = h["strain/Strain"][:]
f_psd, P_psd = median_welch(highpass(np.concatenate([strain[:8 * FS], strain[-8 * FS:]])), 4 * FS)

f = np.fft.rfftfreq(SEG_N, 1.0 / FS)
df = f[1] - f[0]
m = (f >= FMIN) & (f <= FMAX)
fm = f[m]
Sn = np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                      np.log(np.clip(f_psd, 1e-3, None)), np.log(P_psd)))

G = 6.67430e-11
CL = 2.99792458e8
MSUN = 1.98892e30
MC = 28.096 * G * MSUN / CL ** 3
DL = 410.0 * 3.0856775814913673e22 / CL
H_TILDE = (1.0 / DL) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) * MC ** (5.0 / 6.0) \
    * fm ** (-7.0 / 6.0)
PHI = (3.0 / 128.0) * (np.pi * MC * fm) ** (-5.0 / 3.0) - np.pi / 4.0
SIG_ANALYTIC = math.sqrt(4.0 * df * float(np.sum(H_TILDE ** 2 / Sn[m])))
print(f"analytic optimal SNR of the template = {SIG_ANALYTIC:.4f}")

# A. conventions
for label, H, C in (("H=fs*h~, C=4df/fs^2", FS * H_TILDE, 4 * df / FS ** 2),
                    ("H=h~,    C=4df", H_TILDE, 4 * df),
                    ("H=h~/fs, C=4df/fs^2", H_TILDE / FS, 4 * df / FS ** 2)):
    sig = math.sqrt(C * float(np.sum(H ** 2 / Sn[m])))
    print(f"A  {label:24s} -> sigma={sig:.6f}")

# chosen convention: H = fs*h~, C = 4 df/fs^2
C = 4 * df / FS ** 2
H0 = FS * H_TILDE
SIG = math.sqrt(C * float(np.sum(H0 ** 2 / Sn[m])))
WIN = tukey(SEG_N, 0.25)


def rho_curve(data, tc_grid, window=True):
    x = data * WIN if window else data
    X = np.fft.rfft(x)[m]
    base = C * X * np.conj(H0) / Sn[m]
    ph = np.exp(2j * np.pi * np.outer(tc_grid, fm))
    return np.abs(ph @ base) / SIG


tc_grid = np.arange(0.0, SEG_N / FS, 1.0 / FS)

# B. real data
seg_t0 = T_EVENT - 2.0
i0 = int(round(seg_t0 * FS))
seg = highpass(strain[i0:i0 + SEG_N])
for window in (True, False):
    rho = rho_curve(seg, tc_grid, window)
    k = int(np.argmax(rho))
    print(f"B  window={window!s:5s} peak rho={rho[k]:8.2f} at seg t={k/FS:.4f} "
          f"(GPS {GPS0 + seg_t0 + k/FS:.2f}; merger at {GPS0 + T_EVENT:.2f})")

# C. off-source control
j0 = int(round(4.0 * FS))
seg_off = highpass(strain[j0:j0 + SEG_N])
rho_off = rho_curve(seg_off, tc_grid, True)
print(f"C  off-source (4-8 s) peak rho={rho_off.max():.2f} at t={rho_off.argmax()/FS:.3f}")

# D. injection recovery
rng = np.random.default_rng(20260927)
W = np.fft.rfft(rng.normal(0.0, 1.0, SEG_N)) * np.sqrt(Sn * FS / 2.0)
noise = np.fft.irfft(W, SEG_N)
Href = np.zeros(len(f), dtype=complex)
Href[m] = H0 * np.exp(-2j * np.pi * fm * 2.0)
inj = np.fft.irfft(Href, SEG_N)
sig_inj = math.sqrt(C * float(np.sum(np.abs(Href[m]) ** 2 / Sn[m])))
inj = inj * (20.0 / sig_inj)
for label, data in (("noise", noise), ("noise+inj20", noise + inj)):
    rho = rho_curve(data, tc_grid, True)
    k = int(np.argmax(rho))
    print(f"D  {label:12s} peak rho={rho[k]:7.3f} at t={k/FS:.4f}")
