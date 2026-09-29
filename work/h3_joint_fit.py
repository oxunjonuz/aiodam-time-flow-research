#!/usr/bin/env python3
"""H3 (owner's msg_00227): JOINT H1+L1 fit with ONE tau2 and free (M_c, t_c, phi_c),
with the ANTENNA RESPONSES F+, Fx at the GW150914 sky position -- instead of the
single-detector inclination scan of msg224/msg226.

WHY THIS REPLACES THE INCLINATION SCAN.

msg226 answered "does the published single-detector SNR get reached?" by scanning
inclination iota on ONE detector.  That is a question about one detector's
optimal SNR, not about the source.  The source has a SKY POSITION; once it is
known the two detectors are not independent -- the antenna beam patterns F+,d and
Fx,d at that position tie their amplitude AND phase together.

THE TWO STRUCTURAL FACTS THAT MAKE THIS TRACTABLE (both measured, not assumed).

(1) phenomxpy's IMRPhenomT returns the polarizations with the convention

        hp(iota, phi_ref) = A(iota) * exp(+2 i phi_ref) * hp0
        hc(iota, phi_ref) = B(iota) * exp(+2 i phi_ref) * hc0
        hc0 = -i * hp0            (measured: max deviation 3.3e-4)

    with A = (1+cos^2 iota)/2, B = cos iota.  So the strain projected onto a
    detector is

        F+ hp + Fx hc = exp(2 i phi_ref) * hp0 * ( F+ A - i Fx B )

    i.e. the ENTIRE detector response collapses into ONE COMPLEX COEFFICIENT

        C_d = F+,d * A(iota) - i * Fx,d * B(iota).

    Consequence: phi_ref multiplies every detector's C_d by the SAME phase, so it
    is exactly degenerate with an overall phase and is NOT identifiable in a
    coherent network fit with a fixed sky position.  It is profiled out, and that
    is reported as a finding rather than hidden.

(2) Because h_d = C_d * h for a common h, the coherent network statistic is

        rho_coh^2 = | sum_d conj(C_d) zeta_d |^2 / sum_d |C_d|^2 s_d^2

    with zeta_d the complex matched-filter output of the (1,0) waveform on
    detector d and s_d^2 its norm.  zeta_d does NOT depend on C_d, so it is
    computed ONCE per (M_c, tau2) per detector and the whole (iota, psi) grid is
    then free.  The fit is genuinely coherent: the two detectors' amplitude and
    phase ratio is PREDICTED by the sky position, not a free knob.

WHAT IS FITTED.  For a tau2 on a grid, the delay phase
    dPsi(f) = 2 pi f * tau2 * (t_kernel(f) - t_ref)^2
is applied to the IMR template; the network SNR is maximised over the coalescence
time t_c (dense scan) and over (M_c, iota, psi) on grids.  A genuine profile
likelihood over tau2 with everything else profiled.

TWO DELAY KERNELS, both reported (they differ by 155x at 250 Hz, msg226):
  * "pn"  : the inspiral emission-time formula (theory-side reading)
  * "imr" : the IMR waveform's OWN measured group delay (what the waveform does)

CONTROLS (a failed CONTROL means the instrument is broken and the run is void):
  NC1  antenna-off (incoherent sum) must reproduce the RECORDED numbers of
       h3_direct_fit_imr.json: profile at tau2=0 of 7.273249726912924 (H1) /
       5.577882654444806 (L1) and rho_opt 31.659512701612417 / 28.73504152105233.
  NC2  the frozen sky samples must carry the published delay annulus: their
       probability-weighted delay must peak at 6.9 ms (arXiv:1602.03840).
  NC3  the frozen sky samples' 90% area must be near the published 610 deg2.
  NC4  F+^2 + Fx^2 must be independent of psi -- checked, not assumed.
  NC5  the antenna-predicted optimal-SNR ratio rho_opt(H1)/rho_opt(L1) must land
       near the published single-detector ratio 19.5/13.3 = 1.466.
  NC6  noiseless round trip: a template injected at a target network SNR must
       come back at that SNR.
  NC7  power control: the same fit must RECOVER an injected tau2 of 20 ms and
       must NOT recover the paper's 1.2 ms.
  NC8  the coherent statistic must be <= the incoherent one for every template
       (Cauchy-Schwarz on the network sum).

HONEST LIMITS, repeated in the artifact: spins = 0; the delay is PHASE-ONLY;
eta is held at the published mass ratio 29/36; the sky position is READ from the
published LALInference skymap (frozen by sky_prep.py), not re-measured from the
strain; tau(t) is still not DERIVED from the theory.
"""
import json
import math
import os
import sys

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

from phenomxpy import IMRPhenomT

FS = 16384
GPS0 = 1126259447
GPS_EVENT = 1126259462.4
T_EVENT = GPS_EVENT - GPS0
FMIN, FMAX = 20.0, 300.0
MC_MSUN = 28.096
DL_MPC = 410.0
Q = 29.0 / 36.0
ETA = Q / (1.0 + Q) ** 2
F_MAX_FD = 1024.0
DF_FD = 0.25
SEED = 20260928
N_INJ = 24

G = 6.67430e-11
CL = 2.99792458e8
MSUN = 1.98892e30
TSUN = G * MSUN / CL ** 3
MPC = 3.0856775814913673e22
TMPC = MPC / CL
MC = MC_MSUN * TSUN
DL = DL_MPC * TMPC
T_SPAN_S = 0.8444822952069388
PAPER_AMPLITUDE_MS = 1.2

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
FILES = {"H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}
SKY_NPZ = os.path.join(HERE, "artifacts", "sky_samples.npz")

REC_LO_PROFILE_AT_ZERO = {"H1": 7.273249726912924, "L1": 5.577882654444806}
REC_IMR_PEAK = {"H1": 16.74195690868695, "L1": 12.844119129032624}
REC_IMR_RHO_OPT = {"H1": 31.659512701612417, "L1": 28.73504152105233}
PUBLISHED_DELAY_MS = 6.9
PUBLISHED_AREA_90 = 610.0
PUBLISHED_RHO_HAT = {"H1": 19.5, "L1": 13.3}

# official LAL detector geometry: lal/lib/tools/LALDetectors.h (src_750b9573aef0)
ARM = {
    "H1": (np.array([-0.22389266154, 0.79983062746, 0.55690487831]),
           np.array([-0.91397818574, 0.02609403989, -0.40492342125])),
    "L1": (np.array([-0.95457412153, -0.14158077340, -0.26218911324]),
           np.array([0.29774156894, -0.48791033647, -0.82054461286])),
}
VERTEX = {
    "H1": np.array([-2.16141492636e6, -3.83469517889e6, 4.60035022664e6]),
    "L1": np.array([-7.42760447238e4, -5.49628371971e6, 3.22425701744e6]),
}


def gmst_rad(gps):
    unix = 315964800.0 + gps - 18.0
    jd = unix / 86400.0 + 2440587.5
    T = (jd - 2451545.0) / 36525.0
    g = (280.46061837 + 360.98564736629 * (jd - 2451545.0)
         + 0.000387933 * T * T - T * T * T / 38710000.0)
    return math.radians(g % 360.0)


def antenna(det, ra, dec, psi, gps):
    """F+, Fx for a detector at a sky position (LAL/pycbc convention)."""
    gha = gmst_rad(gps) - ra
    cgha, sgha = math.cos(gha), math.sin(gha)
    cdec, sdec = math.cos(dec), math.sin(dec)
    cpsi, spsi = math.cos(psi), math.sin(psi)
    xa, ya = ARM[det]
    D = 0.5 * (np.outer(ya, ya) - np.outer(xa, xa))
    v = np.array([-cpsi * sgha - spsi * cgha * sdec,
                  -cpsi * cgha + spsi * sgha * sdec,
                  spsi * cdec])
    w = np.array([spsi * sgha - cpsi * cgha * sdec,
                  spsi * cgha + cpsi * sgha * sdec,
                  cpsi * cdec])
    return float(v @ D @ v - w @ D @ w), float(v @ D @ w + w @ D @ v)


def load(det):
    with h5py.File(os.path.join(DATA, FILES[det]), "r") as h:
        return h["strain/Strain"][:], int(h["meta/GPSstart"][()])


def highpass(x, fs, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (fs / 2.0), btype="highpass", output="sos"), x)


def median_welch(x, fs, nper, overlap=0.5):
    step = int(nper * (1.0 - overlap))
    w = np.hanning(nper)
    norm = 2.0 / (fs * float(np.sum(w ** 2)))
    segs = np.array([np.abs(np.fft.rfft(x[i:i + nper] * w)) ** 2 * norm
                     for i in range(0, len(x) - nper, step)])
    return np.fft.rfftfreq(nper, 1.0 / fs), np.median(segs, axis=0) / math.log(2.0)


def psd_on_grid(f, f_psd, P_psd):
    return np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                            np.log(np.clip(f_psd, 1e-3, None)), np.log(P_psd)))


def t_of_f_pn(f, Mc=MC):
    return -(5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)


def amp_cont(f, Mc=MC, DL_=DL):
    return (1.0 / DL_) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)


def imr_polarizations(mc_msun, freqs):
    mtot = mc_msun / ETA ** 0.6
    wf = IMRPhenomT(eta=ETA, s1=[0, 0, 0], s2=[0, 0, 0], f_min=FMIN, f_ref=FMIN,
                    total_mass=mtot, distance=DL_MPC, inclination=0.0,
                    phi_ref=0.0, f_max=F_MAX_FD, delta_t=0.5 / F_MAX_FD,
                    delta_f=DF_FD, condition=True)
    fd = wf.compute_fd_polarizations()
    hp = np.asarray(fd[0], dtype=complex)
    hc = np.asarray(fd[1], dtype=complex)
    f_model = np.arange(len(hp)) * DF_FD
    Hp = np.interp(freqs, f_model, hp.real) + 1j * np.interp(freqs, f_model, hp.imag)
    Hc = np.interp(freqs, f_model, hc.real) + 1j * np.interp(freqs, f_model, hc.imag)
    return Hp, Hc


def imr_group_delay(hp, f_model):
    ph = np.unwrap(np.angle(hp))
    return -(np.gradient(ph, DF_FD) / (2.0 * np.pi))


class DetectorFilter:
    """One detector's matched filter for a (1,0) waveform.

    zeta(t) = C4 * ifft( rfft(d w) conj(H) / S ) on the t = n/fs grid, where
    H = fs * base * exp(i dPsi) * exp(+2 pi i f tau_d).

    The light-travel phase exp(+2 pi i f tau_d), with tau_d = r_d . n / c, is the
    arrival-time offset of this detector relative to the geocenter: the wave
    reaches L1 6.9 ms before H1, and without this term the two detectors' peaks
    do NOT align and the coherent statistic is destroyed.  MEASURED over the
    (iota, psi) grid at tau2 = 0: coherent 21.08 with this sign, 16.62 with the
    opposite one, against 21.10 incoherent -- a correct phase relation must bring
    the coherent statistic up to the incoherent one, and only this sign does.

    The template's antenna coefficient C is NOT in zeta: z = conj(C) zeta and
    sigma = |C| s, so the coherent network statistic needs zeta and s only --
    that is what makes the (iota, psi) grid free.
    """

    def __init__(self, seg_n, fs, fmin, fmax, Sn, base, kernel_t=None, t_ref=0.0,
                 tau_d=0.0):
        self.n = seg_n
        self.fs = fs
        self.f = np.fft.rfftfreq(seg_n, 1.0 / fs)
        self.df = self.f[1] - self.f[0]
        self.m = (self.f >= fmin) & (self.f <= fmax)
        self.fm = self.f[self.m]
        self.Sn = Sn[self.m]
        self.C4 = 4.0 * self.df / fs ** 2
        self.base = base[self.m]
        self.kernel_t = None if kernel_t is None else kernel_t[self.m]
        self.t_ref = t_ref
        self.tau_d = tau_d
        self.win = tukey(seg_n, 0.25)
        self._X = None

    def set_data(self, data):
        self._X = np.fft.rfft(data * self.win)[self.m]

    def _psi(self, tau2):
        # arrival offset: a detector at r sees the wave at t - (r.n)/c, so the
        # template carries exp(+2 pi i f tau_d) with tau_d = r.n/c.
        # MEASURED, over the (iota, psi) grid at tau2 = 0: this sign gives a
        # coherent network SNR of 21.08 (against 21.10 incoherent -- they must
        # nearly coincide when the phase relation is right), while the opposite
        # sign gives 16.62.  The sign was settled by that measurement.
        psi = 2.0 * np.pi * self.fm * self.tau_d
        if tau2 != 0.0 and self.kernel_t is not None:
            psi = psi + 2.0 * np.pi * self.fm * tau2 * (self.kernel_t - self.t_ref) ** 2
        return psi

    def H1(self, tau2=0.0):
        return self.fs * self.base * np.exp(1j * self._psi(tau2))

    def s(self, tau2=0.0):
        return math.sqrt(self.C4 * float(np.sum(np.abs(self.H1(tau2)) ** 2 / self.Sn)))

    def zeta(self, tau2=0.0, stride=1):
        b = self.C4 * self._X * np.conj(self.H1(tau2)) / self.Sn
        Br = np.zeros(self.n // 2 + 1, dtype=complex)
        Br[self.m] = b
        B = np.zeros(self.n, dtype=complex)
        B[:self.n // 2 + 1] = Br
        return self.n * np.fft.ifft(B)[::stride]

    def rho_single(self, C, tau2=0.0, stride=1):
        return np.abs(self.zeta(tau2, stride)) * abs(C) / (abs(C) * self.s(tau2))

    def template_time_series(self, tc, C, tau2=0.0):
        H = np.zeros(len(self.f), dtype=complex)
        H[self.m] = self.fs * (C * self.base) * np.exp(1j * self._psi(tau2)) \
            * np.exp(-2j * np.pi * self.fm * tc)
        return np.fft.irfft(H, self.n)


def main():
    out = {"setup": {}, "sky": {}, "antenna": {}, "kernels": {}, "controls": {},
           "control_values": {}, "findings": {}, "thresholds": {}, "answer": {},
           "injection_recovery": {}, "per_grid": {}, "notes": {}}
    seg_n = 4 * FS
    i0 = int(round((T_EVENT - 2.0) * FS))
    freqs = np.fft.rfftfreq(seg_n, 1.0 / FS)
    assert abs((freqs[1] - freqs[0]) - DF_FD) < 1e-12

    psds, segs = {}, {}
    for det in ("H1", "L1"):
        x, gps0 = load(det)
        assert gps0 == GPS0
        off = np.concatenate([x[:8 * FS], x[-8 * FS:]])
        psds[det] = median_welch(highpass(off, FS), FS, FS, 0.5)
        segs[det] = highpass(x[i0:i0 + seg_n], FS)
    Sn = {d: psd_on_grid(freqs, *psds[d]) for d in ("H1", "L1")}

    # ---------------- frozen sky samples (produced by sky_prep.py)
    sky = np.load(SKY_NPZ)
    sky_ra, sky_dec, sky_w = sky["ra"], sky["dec"], sky["weight"]
    ra_ml = float(sky["ml_ra"][0])
    dec_ml = float(sky["ml_dec"][0])
    area90 = float(sky["area90_deg2"][0])
    delay_peak_ms = float(sky["delay_peak_ms"][0])
    # probability-weighted delay over the frozen samples (NC2)
    n_s = np.stack([np.cos(sky_dec) * np.cos(sky_ra),
                    np.cos(sky_dec) * np.sin(sky_ra),
                    np.sin(sky_dec)], axis=1)
    g = gmst_rad(GPS_EVENT)
    c, s = math.cos(g), math.sin(g)
    R = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])
    proj_s = (VERTEX["H1"] - VERTEX["L1"]) @ (R @ n_s.T) / CL * 1e3
    delay_weighted_ms = float((proj_s * sky_w).sum())
    out["sky"] = {
        "samples_file": os.path.relpath(SKY_NPZ, HERE),
        "n_samples": int(len(sky_ra)),
        "ml_ra_deg": math.degrees(ra_ml), "ml_dec_deg": math.degrees(dec_ml),
        "area90_deg2": area90,
        "delay_annulus_peak_ms": delay_peak_ms,
        "delay_weighted_over_samples_ms": delay_weighted_ms,
        "published_delay_ms": PUBLISHED_DELAY_MS,
        "published_area90_deg2": PUBLISHED_AREA_90,
        "sign": "proj = (r_H1-r_L1).n/c = -6.9 ms; the arrival difference "
                "t_H1-t_L1 = -proj = +6.9 ms (L1 first) -- the published statement.",
        "why": "the sky position is READ from the published LALInference skymap, "
               "not re-measured from the strain; it is what turns the two "
               "detectors from independent into coherent.",
    }

    ant = {}
    for det in ("H1", "L1"):
        row = {}
        for psi_deg in (0.0, 30.0, 60.0, 90.0, 120.0, 150.0):
            Fp, Fc = antenna(det, ra_ml, dec_ml, math.radians(psi_deg), GPS_EVENT)
            row[f"{psi_deg:.1f}"] = {"Fp": Fp, "Fx": Fc, "Fp2_Fx2": Fp * Fp + Fc * Fc}
        ant[det] = row
    out["antenna"] = {
        "at_ml_position": ant,
        "psi_independence_of_Fp2_Fx2": {
            d: float(np.ptp([v["Fp2_Fx2"] for v in ant[d].values()])) for d in ant},
        "convention": "F+ = x.D.x - y.D.y, Fx = x.D.y + y.D.x, D = (y y^T - x x^T)/2",
    }

    Hp, Hc = imr_polarizations(MC_MSUN, freqs)
    base = Hp
    f_model = np.arange(len(Hp)) * DF_FD
    t_imr_grid = np.interp(freqs, f_model, imr_group_delay(Hp, DF_FD))
    t_pn_grid = t_of_f_pn(np.clip(freqs, 1e-3, None))
    t_ref = float(t_pn_grid[np.argmin(abs(freqs - FMIN))])
    # leading-order inspiral template, on the same grid (for NC1)
    lo_base = amp_cont(np.clip(freqs, 1e-3, None)) * np.exp(
        1j * ((3.0 / 128.0) * (np.pi * MC * np.clip(freqs, 1e-3, None)) ** (-5.0 / 3.0)
              - np.pi / 4.0))
    # per-detector light-travel time tau_d = r_d . n / c, with n the source
    # direction in the EARTH-FIXED frame (GMST rotation applied)
    g_gmst = gmst_rad(GPS_EVENT)
    cg, sg = math.cos(g_gmst), math.sin(g_gmst)
    Rg = np.array([[cg, sg, 0.0], [-sg, cg, 0.0], [0.0, 0.0, 1.0]])
    n_ml = Rg @ np.array([math.cos(dec_ml) * math.cos(ra_ml),
                          math.cos(dec_ml) * math.sin(ra_ml),
                          math.sin(dec_ml)])
    tau_det = {d: float(VERTEX[d] @ n_ml) / CL for d in ("H1", "L1")}
    out["kernels"]["tau_det_s"] = tau_det
    out["kernels"]["tau_H1_minus_tau_L1_ms"] = (tau_det["H1"] - tau_det["L1"]) * 1e3
    out["kernels"] = {
        "t_pn_at_20_Hz_s": float(t_pn_grid[np.argmin(abs(freqs - 20.0))]),
        "t_imr_at_20_Hz_s": float(t_imr_grid[np.argmin(abs(freqs - 20.0))]),
        "t_pn_at_250_Hz_s": float(t_pn_grid[np.argmin(abs(freqs - 250.0))]),
        "t_imr_at_250_Hz_s": float(t_imr_grid[np.argmin(abs(freqs - 250.0))]),
        "t_ref_used_s": t_ref,
        "tau_det_s": tau_det,
        "tau_H1_minus_tau_L1_ms": (tau_det["H1"] - tau_det["L1"]) * 1e3,
        "hc0_over_hp0": "measured -i (max deviation 3.3e-4)",
    }

    sel_f = (freqs >= FMIN) & (freqs <= FMAX) & (np.abs(Hp) > 1e-25)
    ratio = Hc[sel_f] / Hp[sel_f]
    out["setup"] = {
        "eta": ETA, "chirp_mass_Msun": MC_MSUN, "distance_Mpc": DL_MPC,
        "df_Hz": DF_FD, "seg_n": seg_n, "phenomxpy": "IMRPhenomT",
        "spins": "zero", "delay": "phase-only", "profiled": "t_c (dense scan)",
        "scanned": "tau2, M_c, iota, psi",
        "separable_projection": "F+ hp + Fx hc = exp(2 i phi_ref) hp0 (F+ A - i Fx B)",
        "hc0_over_hp0_mean": [float(ratio.real.mean()), float(ratio.imag.mean())],
        "hc0_over_hp0_max_dev": float(np.abs(ratio - ratio.mean()).max()),
        "coherent_statistic": "rho_coh^2 = |sum_d conj(C_d) zeta_d|^2 / sum_d |C_d|^2 s_d^2",
    }

    tau2_paper = (PAPER_AMPLITUDE_MS * 1e-3) / T_SPAN_S ** 2
    grid_rec = np.linspace(-8.0, 8.0, 33) * tau2_paper
    grid_exp = np.linspace(-32.0, 32.0, 129) * tau2_paper
    amp = lambda x: float(x * T_SPAN_S ** 2 * 1e3)          # noqa: E731
    out["grids"] = {"recorded_half_span_ms": 9.6, "expanded_half_span_ms": 38.4,
                    "tau2_paper_per_s2": float(tau2_paper)}

    # ---------------- NC1: reproduces the recorded numbers, template by template
    # The recorded profile at tau2 = 0 (7.273249726912924 / 5.577882654444806)
    # came from the LEADING-ORDER INSPIRAL template; the recorded rho_opt
    # (31.659512701612417 / 28.73504152105233) and peak (16.74195690868695 /
    # 12.844119129032624) came from the FULL IMR template.  Comparing the IMR
    # template against the LO number is a category error -- my first version did
    # exactly that and NC1 failed for the wrong reason.
    #
    # The recorded pipeline had NO light-travel phase (each detector was fitted
    # independently), so the reproduction uses tau_d = 0; the joint fit adds it.
    # That difference is measured and reported rather than hidden: with tau_d = 0
    # the numbers match the recorded ones to <1e-6, and with tau_d the single
    # detector peak moves by ~1e-4 (the 6.9 ms shift is a fraction of the 1/16384 s
    # sample spacing times the number of samples, so it barely changes the peak).
    nc1 = {}
    for det in ("H1", "L1"):
        flt_lo = DetectorFilter(seg_n, FS, FMIN, FMAX, Sn[det], lo_base, tau_d=0.0)
        flt_lo.set_data(segs[det])
        lo_zero = float(np.max(flt_lo.rho_single(1.0, stride=2)))
        flt_imr = DetectorFilter(seg_n, FS, FMIN, FMAX, Sn[det], base, tau_d=0.0)
        flt_imr.set_data(segs[det])
        imr_zero = float(np.max(flt_imr.rho_single(1.0, stride=2)))
        nc1[det] = {
            "LO_profile_at_zero_coarse": lo_zero,
            "recorded_LO_profile_at_zero": REC_LO_PROFILE_AT_ZERO[det],
            "IMR_peak_coarse": imr_zero,
            "recorded_IMR_peak": REC_IMR_PEAK[det],
            "rho_opt": flt_imr.s(0.0),
            "recorded_rho_opt": REC_IMR_RHO_OPT[det],
        }
    out["control_values"]["NC1_recorded_numbers"] = nc1

    # ---------------- the joint fit
    mc_grid = MC_MSUN * np.array([0.95, 0.975, 1.0, 1.025, 1.05])
    iota_grid = np.radians(np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0,
                                     120.0, 140.0, 160.0, 180.0]))
    psi_grid = np.radians(np.array([0.0, 22.5, 45.0, 67.5, 90.0, 112.5, 135.0]))
    combos = []
    for i_iota, iota in enumerate(iota_grid):
        A = (1.0 + math.cos(iota) ** 2) / 2.0
        B = math.cos(iota)
        for i_psi, psi in enumerate(psi_grid):
            Cd = []
            for d in ("H1", "L1"):
                Fp, Fc = antenna(d, ra_ml, dec_ml, psi, GPS_EVENT)
                Cd.append(Fp * A - 1j * Fc * B)
            combos.append({"iota_deg": math.degrees(iota),
                           "psi_deg": math.degrees(psi), "C": np.array(Cd)})
    Cmat = np.array([c["C"] for c in combos])            # (ncombo, 2)

    def joint_profile(grid, kernel_t, stride=2):
        prof, inc_prof, best = [], [], []
        for t2 in grid:
            best_val, best_cfg, best_inc = -1.0, None, -1.0
            for mc in mc_grid:
                Hp_i, _ = imr_polarizations(float(mc), freqs)
                Z, s2 = [], []
                for d in ("H1", "L1"):
                    flt = DetectorFilter(seg_n, FS, FMIN, FMAX, Sn[d], Hp_i, tau_d=tau_det[d],
                                         kernel_t=kernel_t, t_ref=t_ref)
                    flt.set_data(segs[d])
                    Z.append(flt.zeta(float(t2), stride))
                    s2.append(flt.s(float(t2)) ** 2)
                Z = np.array(Z)                          # (2, L)
                s2 = np.array(s2)                        # (2,)
                comb = Cmat.conj() @ Z                   # (ncombo, L)
                num = np.abs(comb) ** 2
                den = (np.abs(Cmat) ** 2) @ s2           # (ncombo,)
                r = num.max(axis=1) / den
                j = int(np.argmax(r))
                if r[j] > best_val:
                    best_val = float(r[j])
                    best_cfg = {"mc": float(mc), "iota_deg": combos[j]["iota_deg"],
                                "psi_deg": combos[j]["psi_deg"]}
                inc = float(np.max(np.abs(Z[0]) ** 2)) / s2[0] \
                    + float(np.max(np.abs(Z[1]) ** 2)) / s2[1]
                best_inc = max(best_inc, inc)
            prof.append(math.sqrt(best_val))
            inc_prof.append(math.sqrt(best_inc))
            best.append(best_cfg)
        return np.array(prof), np.array(inc_prof), best

    res = {}
    for kname, kt in (("pn", t_pn_grid), ("imr", t_imr_grid)):
        for gname, grid in (("recorded", grid_rec), ("expanded", grid_exp)):
            prof, inc_prof, best = joint_profile(grid, kt)
            k = int(np.argmax(prof))
            target = prof[k] - 1.0
            left = grid[:k][prof[:k] < target]
            right = grid[k:][prof[k:] < target]
            lo = float(left[-1]) if len(left) else None
            hi = float(right[0]) if len(right) else None
            res[f"joint_{kname}_{gname}"] = {
                "delay_kernel": kname, "grid": gname,
                "profile_peak_net_snr": float(prof[k]),
                "profile_at_zero": float(prof[np.argmin(np.abs(grid))]),
                "profile_at_paper": float(prof[np.argmin(np.abs(grid - tau2_paper))]),
                "profile_drop_within_grid": float(prof[k] - float(np.min(prof))),
                "tau2_best_fit_amplitude_ms": amp(grid[k]),
                "best_fit_config": best[k],
                "best_fit_at_grid_edge": bool(k == 0 or k == len(grid) - 1),
                "one_sigma_low_amplitude_ms": None if lo is None else amp(lo),
                "one_sigma_high_amplitude_ms": None if hi is None else amp(hi),
                "unconstrained_on_grid": bool(lo is None and hi is None),
                "incoherent_peak_net_snr": float(np.max(inc_prof)),
                "coherent_le_incoherent": bool(
                    float(np.max(prof)) <= float(np.max(inc_prof)) * (1.0 + 1e-9)),
                "grid_half_span_amplitude_ms": amp(grid[-1]),
                "grid": grid.tolist(), "profile": prof.tolist(),
                "incoherent_profile": inc_prof.tolist(),
            }
    out["per_grid"] = res

    # ---------------- NC6 / NC7: power control on the JOINT statistic
    rng = np.random.default_rng(SEED)
    grid_wide = np.linspace(-64.0, 64.0, 65) * tau2_paper
    flt = {}
    for d in ("H1", "L1"):
        flt[d] = DetectorFilter(seg_n, FS, FMIN, FMAX, Sn[d], base, tau_d=tau_det[d],
                                kernel_t=t_imr_grid, t_ref=t_ref)
    iota0 = math.radians(30.0)
    A0 = (1.0 + math.cos(iota0) ** 2) / 2.0
    B0 = math.cos(iota0)
    C0 = {}
    for d in ("H1", "L1"):
        Fp, Fc = antenna(d, ra_ml, dec_ml, math.radians(0.0), GPS_EVENT)
        C0[d] = Fp * A0 - 1j * Fc * B0
    tc_inj = 2.0

    def net_snr_of(sigs, tau2):
        Z, s2 = [], []
        for d in ("H1", "L1"):
            flt[d].set_data(sigs[d])
            Z.append(flt[d].zeta(tau2, 1))
            s2.append(flt[d].s(tau2) ** 2)
        Z = np.array(Z)
        Cs = np.array([C0["H1"], C0["L1"]])
        comb = Cs.conj() @ Z
        return math.sqrt(float(np.abs(comb).max() ** 2 / (np.abs(Cs) ** 2 @ np.array(s2))))

    rt = {}
    for label, amp_ms in (("zero", 0.0), ("paper", PAPER_AMPLITUDE_MS), ("large", 20.0)):
        t2 = (amp_ms * 1e-3) / T_SPAN_S ** 2
        tot2 = 0.0
        sigs = {}
        for d in ("H1", "L1"):
            sig = flt[d].template_time_series(tc_inj, C0[d], tau2=t2)
            sigs[d] = sig
            tot2 += abs(C0[d]) ** 2 * flt[d].s(t2) ** 2
        scale = 25.0 / math.sqrt(tot2)
        got = net_snr_of({d: sigs[d] * scale for d in sigs}, t2)
        rt[label] = {"target_net_snr": 25.0, "recovered_net_snr": got,
                     "rel_err": abs(got - 25.0) / 25.0}

    found = {}
    for label, amp_ms, grid in (("zero", 0.0, grid_rec),
                                ("paper", PAPER_AMPLITUDE_MS, grid_rec),
                                ("large", 20.0, grid_wide)):
        t2 = (amp_ms * 1e-3) / T_SPAN_S ** 2
        sigs, tot2 = {}, 0.0
        for d in ("H1", "L1"):
            sig = flt[d].template_time_series(tc_inj, C0[d], tau2=t2)
            sigs[d] = sig
            tot2 += abs(C0[d]) ** 2 * flt[d].s(t2) ** 2
        scale = 25.0 / math.sqrt(tot2)
        Cs = np.array([C0["H1"], C0["L1"]])
        best = []
        for _ in range(N_INJ):
            noise = {}
            for d in ("H1", "L1"):
                W = np.fft.rfft(rng.normal(0.0, 1.0, seg_n)) * np.sqrt(Sn[d] * FS / 2.0)
                noise[d] = highpass(np.fft.irfft(W, seg_n), FS)
            # the SAME coherent statistic the fit uses
            zs, s2s = [], []
            for d in ("H1", "L1"):
                flt[d].set_data(noise[d] + sigs[d] * scale)
            pr = []
            for gg in grid:
                Zg, s2g = [], []
                for d in ("H1", "L1"):
                    Zg.append(flt[d].zeta(float(gg), 2))
                    s2g.append(flt[d].s(float(gg)) ** 2)
                comb = Cs.conj() @ np.array(Zg)
                pr.append(float(np.abs(comb).max() ** 2 / (np.abs(Cs) ** 2 @ np.array(s2g))))
            best.append(amp(float(grid[int(np.argmax(pr))])))
        best = np.array(best)
        found[label] = {
            "injected_amplitude_ms": amp_ms,
            "recovered_mean_ms": float(best.mean()),
            "recovered_std_ms": float(best.std()),
            "bias_ms": float(best.mean() - amp_ms),
            "scan_half_span_ms": amp(grid[-1]),
            "recovered_at_scan_edge": bool(
                abs(abs(best.mean()) - abs(amp(grid[-1]))) < 1e-6),
        }
    out["injection_recovery"] = {"noiseless_round_trip": rt, "recovery": found}

    # ---------------- NC5: antenna-predicted optimal-SNR ratio vs published
    ratio_rows = {}
    for iota_deg in (0.0, 30.0, 60.0, 90.0):
        A = (1.0 + math.cos(math.radians(iota_deg)) ** 2) / 2.0
        B = math.cos(math.radians(iota_deg))
        row = {}
        for d in ("H1", "L1"):
            Fp, Fc = antenna(d, ra_ml, dec_ml, math.radians(0.0), GPS_EVENT)
            C = Fp * A - 1j * Fc * B
            row[d] = abs(C) * DetectorFilter(seg_n, FS, FMIN, FMAX, Sn[d], base, tau_d=tau_det[d]).s(0.0)
        ratio_rows[f"iota_{iota_deg:.0f}"] = {
            "rho_opt_H1": row["H1"], "rho_opt_L1": row["L1"],
            "ratio": row["H1"] / row["L1"]}
    out["antenna"]["optimal_snr_ratio"] = {
        "rows": ratio_rows,
        "published_single_detector_ratio": PUBLISHED_RHO_HAT["H1"] / PUBLISHED_RHO_HAT["L1"],
        "why": "with the sky position fixed the two detectors are tied; the ratio "
               "of their optimal SNRs is a PREDICTION, not a free knob.",
    }

    # ---------------- controls
    c = {}
    c["NC1_recorded_numbers_reproduced"] = all(
        abs(nc1[d]["LO_profile_at_zero_coarse"] - nc1[d]["recorded_LO_profile_at_zero"]) < 1e-6
        and abs(nc1[d]["IMR_peak_coarse"] - nc1[d]["recorded_IMR_peak"]) < 1e-6
        and abs(nc1[d]["rho_opt"] - nc1[d]["recorded_rho_opt"]) < 1e-6
        for d in ("H1", "L1"))
    c["NC2_sky_samples_carry_published_delay"] = abs(delay_peak_ms + PUBLISHED_DELAY_MS) < 0.6
    c["NC3_sky_samples_area_matches_published"] = (
        abs(area90 - PUBLISHED_AREA_90) / PUBLISHED_AREA_90 < 0.05)
    c["NC4_Fp2_Fx2_independent_of_psi"] = all(
        v < 1e-9 for v in out["antenna"]["psi_independence_of_Fp2_Fx2"].values())
    r0 = ratio_rows["iota_0"]["ratio"]
    c["NC5_antenna_ratio_matches_published"] = abs(
        r0 - PUBLISHED_RHO_HAT["H1"] / PUBLISHED_RHO_HAT["L1"]) < 0.10
    c["NC6_noiseless_round_trip"] = all(v["rel_err"] < 1e-3 for v in rt.values())
    c["NC7_large_injected_delay_recovered"] = (
        abs(found["large"]["recovered_mean_ms"] - 20.0) < 6.0
        and not found["large"]["recovered_at_scan_edge"])
    c["NC8_paper_injection_not_resolved"] = (
        found["paper"]["recovered_std_ms"] > PAPER_AMPLITUDE_MS)
    c["NC9_coherent_le_incoherent"] = all(
        res[k]["coherent_le_incoherent"] for k in res)
    out["controls"] = {k: bool(v) for k, v in c.items()}
    out["control_values"].update({
        "delay_annulus_peak_ms": delay_peak_ms,
        "delay_weighted_over_samples_ms": delay_weighted_ms,
        "area90_deg2": area90,
        "antenna_ratio_at_iota0": r0,
        "published_ratio": PUBLISHED_RHO_HAT["H1"] / PUBLISHED_RHO_HAT["L1"],
        "noiseless_round_trip": rt,
        "recovery_ms": {k: v["recovered_mean_ms"] for k, v in found.items()},
    })

    T = dict(c)
    T["T_joint_one_sided_bound_with_imr_kernel"] = any(
        (res[f"joint_{k}_expanded"]["one_sigma_high_amplitude_ms"] is not None
         or res[f"joint_{k}_expanded"]["one_sigma_low_amplitude_ms"] is not None)
        for k in ("pn", "imr"))
    T["T_joint_unconstrained_on_recorded_grid"] = all(
        res[f"joint_{k}_recorded"]["unconstrained_on_grid"] for k in ("pn", "imr"))
    T["T_joint_recovers_large_injection"] = bool(
        abs(found["large"]["recovered_mean_ms"] - 20.0) < 6.0)
    # The paper's 1.2 ms is NOT RESOLVED: the recovered scatter (std) is several
    # times the injected amplitude itself.  My first version of this threshold
    # asked whether the recovered MEAN differed from 1.2 ms by more than 1 ms --
    # a bad criterion, because a mean can land near 1.2 ms by chance while the
    # uncertainty is 4.6 ms, which is the opposite of resolving it.  The honest
    # form is the one NC8 uses: the scatter exceeds the signal.
    T["T_joint_cannot_resolve_paper_injection"] = bool(
        found["paper"]["recovered_std_ms"] > PAPER_AMPLITUDE_MS)
    out["thresholds"] = {k: bool(v) for k, v in T.items()}
    out["findings"] = {k: bool(v) for k, v in T.items() if k not in c}
    out["controls_all_pass"] = bool(all(c.values()))

    ans = {}
    for k in ("pn", "imr"):
        fe = res[f"joint_{k}_expanded"]
        f9 = res[f"joint_{k}_recorded"]
        ans[k] = {
            "kernel": k,
            "net_snr_peak": fe["profile_peak_net_snr"],
            "drop_recorded_grid": f9["profile_drop_within_grid"],
            "drop_expanded_grid": fe["profile_drop_within_grid"],
            "one_sigma_low_ms": fe["one_sigma_low_amplitude_ms"],
            "one_sigma_high_ms": fe["one_sigma_high_amplitude_ms"],
            "unconstrained_recorded": f9["unconstrained_on_grid"],
            "best_fit_ms": fe["tau2_best_fit_amplitude_ms"],
            "best_fit_at_edge": fe["best_fit_at_grid_edge"],
        }
    out["answer"]["joint_fit"] = ans
    out["answer"]["statement"] = (
        "The joint H1+L1 fit with ONE tau2 and free (M_c, t_c, phi_c), with the "
        "antenna responses F+,Fx read from the published GW150914 sky position, "
        "is reported BEFORE any bound is stated.  Two structural facts came out "
        "of building it.  (1) phi_c is NOT identifiable: with the sky position "
        "fixed it multiplies every detector's antenna coefficient by the same "
        "phase and is exactly degenerate with an overall phase, so it is "
        "profiled out rather than quoted.  (2) The entire detector projection "
        "collapses to ONE complex coefficient C_d = F+ A(iota) - i Fx B(iota) "
        "(measured identity, hc0 = -i hp0), so fixing the sky position makes the "
        "two detectors COHERENT: their amplitude and phase ratio is predicted, "
        "not free."
    )
    out["notes"] = {
        "power_control_is_the_point": "without NC7, 'the data do not constrain "
                                      "tau2' is indistinguishable from 'the fit "
                                      "cannot see tau2 at all'.",
        "sky_is_read_not_measured": "the position comes from the published "
                                    "LALInference skymap (sky_prep.py); it is "
                                    "not re-derived from the strain.",
        "still_open": "tau(t) is still not DERIVED from the theory; the paper "
                      "does not specify it.",
        "phase_only": "the delay is a phase perturbation only.",
    }

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h3_joint_fit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("sky", "antenna", "controls",
                                          "control_values", "findings",
                                          "thresholds", "answer",
                                          "injection_recovery")},
                     indent=2, sort_keys=True, default=str))
    return 0 if all(c.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
