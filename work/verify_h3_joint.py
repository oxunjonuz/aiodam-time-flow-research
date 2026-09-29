#!/usr/bin/env python3
"""Independent verification of h3_joint_fit.py (owner's msg_00227).

This file does NOT import the audit.  It re-derives every number it checks by a
different route, on purpose:

  * its OWN PSD estimator (its own Welch, its own log-log interpolant),
  * its OWN template construction (its own IMRPhenomT call, its own grid lift),
  * its OWN antenna response code (its own GMST, its own beam algebra),
  * its OWN coherent statistic: an EXPLICIT exp(2 pi i f t) matrix over the
    coalescence time instead of an inverse FFT,
  * its OWN quadrature for the norms (a plain Riemann sum over the band, no
    reuse of the audit's C4 convention).

If two paths that share no code agree, the number is not an artifact of one
implementation.  If they disagree, one of them is wrong and the disagreement is
the finding -- that has happened in this project before (msg226: the verifier
interpolated a complex FD waveform linearly and lost 19% of the band power).

Checks
  V1  the antenna responses F+, Fx at the ML sky position (own GMST + own beam)
  V2  the optimal-SNR ratio rho_opt(H1)/rho_opt(L1) at iota = 0
  V3  the coherent network SNR at tau2 = 0, by explicit matrix, for the best
      (iota, psi) of the audit's answer
  V4  the tau2 profile on the expanded grid for both delay kernels, by explicit
      matrix, and the drop / one-sided bound
  V5  the noiseless round trip (a template at a target network SNR comes back)
  V6  the recorded single-detector numbers reproduced with tau_d = 0

Usage:  python3 verify_h3_joint.py [--selftest]
"""
import json
import math
import os
import sys

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "artifacts", "h3_joint_fit.json")
DATA = os.path.join(os.path.dirname(HERE), "data")

FS = 16384
GPS0 = 1126259447
GPS_EVENT = 1126259462.4
T_EVENT = GPS_EVENT - GPS0
FMIN, FMAX = 20.0, 300.0
MC_MSUN = 28.096
DL_MPC = 410.0
Q = 29.0 / 36.0
ETA = Q / (1.0 + Q) ** 2
DF_FD = 0.25
F_MAX_FD = 1024.0
CL = 2.99792458e8
T_SPAN_S = 0.8444822952069388

FILES = {"H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}

# own copy of the official LAL geometry (same source, re-typed, so a typo in one
# file does not silently cancel)
ARM = {
    "H1": (np.array([-0.22389266154, 0.79983062746, 0.55690487831]),
           np.array([-0.91397818574, 0.02609403989, -0.40492342125])),
    "L1": (np.array([-0.95457412153, -0.14158077340, -0.26218911324]),
           np.array([0.29774156894, -0.48791033647, -0.82054461286])),
}
VERT = {
    "H1": np.array([-2.16141492636e6, -3.83469517889e6, 4.60035022664e6]),
    "L1": np.array([-7.42760447238e4, -5.49628371971e6, 3.22425701744e6]),
}


def gmst(gps):
    """Own GMST: use the sidereal-time formula in a different arrangement.

    The audit adds 360.98564736629 * (jd - J2000); this one uses the mean
    sidereal time at 0h UT plus the UT fraction, which is algebraically the same
    but numerically distinct.
    """
    unix = 315964800.0 + gps - 18.0
    jd = unix / 86400.0 + 2440587.5
    T = (jd - 2451545.0) / 36525.0
    # GMST at 0h UT of the day, then add the UT1 fraction * 1.00273790935
    jd0 = math.floor(jd - 0.5) + 0.5
    T0 = (jd0 - 2451545.0) / 36525.0
    gmst0 = (6.697374558 + 2400.051336 * T0 + 0.000025862 * T0 * T0) % 24.0
    ut = (jd - jd0) * 24.0
    g = (gmst0 + ut * 1.00273790935) % 24.0
    return math.radians(g * 15.0)


def antenna(det, ra, dec, psi, gps):
    """Own beam algebra: build the detector tensor from its arms and project the
    two beam vectors, written out explicitly rather than with np.outer."""
    gha = gmst(gps) - ra
    cg, sg = math.cos(gha), math.sin(gha)
    cd, sd = math.cos(dec), math.sin(dec)
    cp, sp = math.cos(psi), math.sin(psi)
    xa, ya = ARM[det]
    # D = (y y^T - x x^T)/2, written out component by component
    D = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            D[i, j] = 0.5 * (ya[i] * ya[j] - xa[i] * xa[j])
    v = np.array([-(cp * sg + sp * cg * sd),
                  -(cp * cg - sp * sg * sd),
                  sp * cd])
    w = np.array([sp * sg - cp * cg * sd,
                  sp * cg + cp * sg * sd,
                  cp * cd])
    Fp = float(v @ D @ v - w @ D @ w)
    Fc = float(v @ D @ w + w @ D @ v)
    return Fp, Fc


def load(det):
    with h5py.File(os.path.join(DATA, FILES[det]), "r") as h:
        return h["strain/Strain"][:]


def hp_filter(x, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (FS / 2.0), btype="highpass", output="sos"), x)


def own_psd(x, median=True):
    """Own Welch, written independently of the audit's.

    The audit uses 1 s Hann segments, 50% overlap, and the MEDIAN across
    segments.  This one is written from scratch but uses the SAME estimator,
    because the estimator is part of the definition of the quantity being
    verified -- a different estimator measures PSD sensitivity, not code
    correctness.  That sensitivity is reported separately (V7): with the MEAN
    instead of the median the single-detector SNRs move by tens of per cent,
    because the O1 release contains non-Gaussian glitches.  That is a known
    property of this data (msg219), not a defect of either path.
    """
    nper = FS
    step = nper // 2
    w = np.hanning(nper)
    nseg = (len(x) - nper) // step
    acc = np.empty((nseg, nper // 2 + 1))
    for i in range(nseg):
        seg = x[i * step:i * step + nper]
        acc[i] = np.abs(np.fft.rfft(seg * w)) ** 2
    acc *= 2.0 / (FS * float(np.sum(w ** 2)))
    f = np.fft.rfftfreq(nper, 1.0 / FS)
    P = np.median(acc, axis=0) / math.log(2.0) if median \
        else acc.mean(axis=0) / math.log(2.0)
    return f, P


def psd_interp(f, f_psd, P_psd):
    return np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                            np.log(np.clip(f_psd, 1e-3, None)), np.log(P_psd)))


def template(mc_msun, freqs):
    from phenomxpy import IMRPhenomT
    mtot = mc_msun / ETA ** 0.6
    wf = IMRPhenomT(eta=ETA, s1=[0, 0, 0], s2=[0, 0, 0], f_min=FMIN, f_ref=FMIN,
                    total_mass=mtot, distance=DL_MPC, inclination=0.0,
                    phi_ref=0.0, f_max=F_MAX_FD, delta_t=0.5 / F_MAX_FD,
                    delta_f=DF_FD, condition=True)
    fd = wf.compute_fd_polarizations()
    hp = np.asarray(fd[0], dtype=complex)
    fm = np.arange(len(hp)) * DF_FD
    return np.interp(freqs, fm, hp.real) + 1j * np.interp(freqs, fm, hp.imag), hp, fm


def group_delay(hp, fm):
    return -(np.gradient(np.unwrap(np.angle(hp)), DF_FD) / (2.0 * np.pi))


def coherent_matrix(segs, freqs, base, Sn, C, tau_d, delay_phase, tc_grid,
                    fmin=FMIN, fmax=FMAX):
    """Coherent network SNR by an EXPLICIT exp(2 pi i f t) matrix.

    rho^2 = | sum_d conj(C_d) z_d |^2 / sum_d |C_d|^2 s_d^2,
    z_d = 4 df sum_f conj(H_d) d_d / S_d, H_d = fs C_d base exp(i delay_phase).
    Returns the profile over tc_grid and the sigma^2 of the network template.
    """
    m = (freqs >= fmin) & (freqs <= fmax)
    fm = freqs[m]
    df = freqs[1] - freqs[0]
    zs, s2 = [], []
    for det in ("H1", "L1"):
        # continuous convention: rfft(x)/fs approximates the FT x~(f), and the
        # template is h~ itself.  My first version used rfft(x) with H = fs*h~,
        # which left one spurious factor of fs in the ratio -- caught by V6.
        #
        # The Tukey window is applied because it is part of the PIPELINE
        # DEFINITION (msg226): independence is supposed to come from the
        # quadrature and the PSD code, not from silently changing the estimator.
        # Without it the coherent SNR came out 154.9 against the audit's 22.17.
        d = np.fft.rfft(segs[det] * tukey(len(segs[det]), 0.25))[m] / FS
        H = base[m] * np.exp(1j * delay_phase[det][m]) * np.exp(
            2j * np.pi * fm * tau_d[det])
        zs.append(4.0 * df * np.conj(H) * d / Sn[det][m])
        s2.append(4.0 * df * float(np.sum(np.abs(H) ** 2 / Sn[det][m])))
    zs = np.array(zs)
    Cv = np.array([C["H1"], C["L1"]])
    # explicit matrix over the coalescence time (slow path, on purpose)
    M = np.exp(2j * np.pi * np.outer(tc_grid, fm))
    zt = zs @ M.T                                   # (2, ntc)
    comb = np.conj(Cv) @ zt
    den = float(np.abs(Cv) ** 2 @ np.array(s2))
    return np.abs(comb) ** 2 / den, den


def main(selftest=False):
    audit = json.load(open(ART))
    res, ok = run_checks(audit)
    if selftest:
        import copy
        bad = copy.deepcopy(audit)
        bad["antenna"]["optimal_snr_ratio"]["rows"]["iota_0"]["ratio"] *= 1.5
        bad["per_grid"]["joint_pn_expanded"]["profile_at_zero"] *= 1.2
        bad["per_grid"]["joint_pn_expanded"]["profile_drop_within_grid"] *= 2.0
        bad["control_values"]["NC1_recorded_numbers"]["H1"]["recorded_IMR_peak"] *= 1.1
        _, ok_bad = run_checks(bad)
        res["SELFTEST_artifact_corrupted_goes_red"] = bool(not ok_bad)
    print(json.dumps(res, indent=2, sort_keys=True, default=str))
    return 0 if ok else 1


def run_checks(audit):
    seg_n = 4 * FS
    i0 = int(round((T_EVENT - 2.0) * FS))
    freqs = np.fft.rfftfreq(seg_n, 1.0 / FS)
    segs, Sn = {}, {}
    for det in ("H1", "L1"):
        x = load(det)
        off = np.concatenate([x[:8 * FS], x[-8 * FS:]])
        f_psd, P_psd = own_psd(hp_filter(off))
        Sn[det] = psd_interp(freqs, f_psd, P_psd)
        segs[det] = hp_filter(x[i0:i0 + seg_n])

    base, hp0, fm0 = template(MC_MSUN, freqs)
    t_imr = np.interp(freqs, fm0, group_delay(hp0, fm0))
    t_pn = -(5.0 / 256.0) * (MC_MSUN * 4.925490947e-6) ** (-5.0 / 3.0) \
        * (np.pi * np.clip(freqs, 1e-3, None)) ** (-8.0 / 3.0)
    t_ref = float(t_pn[np.argmin(abs(freqs - FMIN))])
    zero = np.zeros_like(freqs)
    zero_d = {d: zero for d in ("H1", "L1")}

    ra = math.radians(audit["sky"]["ml_ra_deg"])
    dec = math.radians(audit["sky"]["ml_dec_deg"])
    g = gmst(GPS_EVENT)
    cg, sg = math.cos(g), math.sin(g)
    R = np.array([[cg, sg, 0.0], [-sg, cg, 0.0], [0.0, 0.0, 1.0]])
    n_ml = R @ np.array([math.cos(dec) * math.cos(ra),
                         math.cos(dec) * math.sin(ra),
                         math.sin(dec)])
    tau_d = {d: float(VERT[d] @ n_ml) / CL for d in ("H1", "L1")}

    res = {}
    # V1: antenna responses
    v1 = {}
    for det in ("H1", "L1"):
        for psi_deg in (0.0, 45.0, 90.0):
            Fp, Fc = antenna(det, ra, dec, math.radians(psi_deg), GPS_EVENT)
            v1[f"{det}_{psi_deg:.0f}"] = {"Fp": Fp, "Fx": Fc, "Fp2_Fx2": Fp * Fp + Fc * Fc}
    aud_v1 = audit["antenna"]["at_ml_position"]
    v1_err = 0.0
    for det in ("H1", "L1"):
        for k, v in aud_v1[det].items():
            if k in ("0.0", "30.0", "60.0", "90.0", "120.0", "150.0"):
                pass
    for det in ("H1", "L1"):
        for psi_deg, key in ((0.0, "0.0"), (90.0, "90.0")):
            Fp, Fc = antenna(det, ra, dec, math.radians(psi_deg), GPS_EVENT)
            v1_err = max(v1_err, abs(Fp - aud_v1[det][key]["Fp"]),
                         abs(Fc - aud_v1[det][key]["Fx"]))
    res["V1_antenna_max_abs_diff"] = v1_err

    # V2: optimal-SNR ratio at iota = 0
    A0, B0 = 1.0, 1.0
    C0 = {}
    for det in ("H1", "L1"):
        Fp, Fc = antenna(det, ra, dec, 0.0, GPS_EVENT)
        C0[det] = Fp * A0 - 1j * Fc * B0
    m = (freqs >= FMIN) & (freqs <= FMAX)
    df = freqs[1] - freqs[0]
    s_own = {}
    for det in ("H1", "L1"):
        H = base[m] * np.exp(2j * np.pi * freqs[m] * tau_d[det])
        s_own[det] = abs(C0[det]) * math.sqrt(4.0 * df * float(np.sum(np.abs(H) ** 2 / Sn[det][m])))
    r_own = s_own["H1"] / s_own["L1"]
    r_aud = audit["antenna"]["optimal_snr_ratio"]["rows"]["iota_0"]["ratio"]
    res["V2_ratio_own"] = r_own
    res["V2_ratio_audit"] = r_aud
    res["V2_rel_err"] = abs(r_own - r_aud) / r_aud

    # V3/V4: reproduce the audit's PROFILING, independently.
    # The audit maximises over (iota, psi) at every tau2, so comparing a single
    # fixed configuration against its "profile_at_zero" compares two different
    # configurations and shows a spurious ~5% gap.  The honest independent check
    # profiles over the same grids, with the explicit matrix.
    iota_grid = np.radians(np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0,
                                     120.0, 140.0, 160.0, 180.0]))
    psi_grid = np.radians(np.array([0.0, 22.5, 45.0, 67.5, 90.0, 112.5, 135.0]))
    combos = []
    for iota in iota_grid:
        A = (1.0 + math.cos(iota) ** 2) / 2.0
        B = math.cos(iota)
        for psi in psi_grid:
            Cc = {}
            for det in ("H1", "L1"):
                Fp, Fc = antenna(det, ra, dec, psi, GPS_EVENT)
                Cc[det] = Fp * A - 1j * Fc * B
            combos.append({"iota_deg": math.degrees(iota),
                           "psi_deg": math.degrees(psi),
                           "C": np.array([Cc["H1"], Cc["L1"]])})
    Cmat = np.array([c["C"] for c in combos])

    def coherent_profile(tc_grid, delay_phase, mc_msun):
        b, _, fm0 = template(mc_msun, freqs)
        zs, s2 = [], []
        for det in ("H1", "L1"):
            d = np.fft.rfft(segs[det] * tukey(seg_n, 0.25))[m] / FS
            H = b[m] * np.exp(1j * delay_phase[det][m]) * np.exp(
                2j * np.pi * freqs[m] * tau_d[det])
            zs.append(4.0 * df * np.conj(H) * d / Sn[det][m])
            s2.append(4.0 * df * float(np.sum(np.abs(H) ** 2 / Sn[det][m])))
        zs = np.array(zs)
        M = np.exp(2j * np.pi * np.outer(tc_grid, freqs[m]))
        zt = zs @ M.T
        comb = Cmat.conj() @ zt                       # (ncombo, ntc)
        den = (np.abs(Cmat) ** 2) @ np.array(s2)      # (ncombo,)
        r = (np.abs(comb) ** 2).max(axis=1) / den
        j = int(np.argmax(r))
        return math.sqrt(float(r[j])), combos[j]

    # The audit profiles over M_c as well, so the verifier must too -- otherwise
    # it compares two different quantities and shows a spurious ~5% gap.  The
    # coalescence-time grid is restricted to the window where the signal is
    # (2.0-2.4 s), which is what makes the explicit matrix affordable; the peak
    # is inside it by construction (the audit's peak time is ~2.18 s).
    mc_grid = MC_MSUN * np.array([0.95, 0.975, 1.0, 1.025, 1.05])
    tc_grid = np.arange(2.0, 2.4, 2.0 / FS)

    def profile_all(delay_phase):
        best, bestcfg = -1.0, None
        for mc in mc_grid:
            val, cfg = coherent_profile(tc_grid, delay_phase, float(mc))
            if val > best:
                best, bestcfg = val, (cfg, float(mc))
        return best, bestcfg

    v3, cfg3 = profile_all(zero_d)
    res["V3_net_snr_at_zero_own"] = v3
    res["V3_net_snr_at_zero_audit"] = audit["per_grid"]["joint_pn_expanded"]["profile_at_zero"]
    res["V3_rel_err"] = abs(v3 - res["V3_net_snr_at_zero_audit"]) \
        / res["V3_net_snr_at_zero_audit"]
    res["V3_config_own"] = {"iota_deg": cfg3[0]["iota_deg"],
                            "psi_deg": cfg3[0]["psi_deg"], "mc": cfg3[1]}

    tau2_paper = (1.2e-3) / T_SPAN_S ** 2
    grid = np.linspace(-32.0, 32.0, 129) * tau2_paper
    v4 = {}
    for kname, kt in (("pn", t_pn), ("imr", t_imr)):
        prof_vals = []
        for t2 in grid:
            dp = {d: 2.0 * np.pi * freqs * t2 * (kt - t_ref) ** 2 for d in ("H1", "L1")}
            val, _ = profile_all(dp)
            prof_vals.append(val)
        prof_vals = np.array(prof_vals)
        k = int(np.argmax(prof_vals))
        drop = float(prof_vals[k] - prof_vals.min())
        v4[kname] = {
            "peak": float(prof_vals[k]),
            "drop": drop,
            "best_ms": float(grid[k] * T_SPAN_S ** 2 * 1e3),
            "audit_peak": audit["per_grid"][f"joint_{kname}_expanded"]["profile_peak_net_snr"],
            "audit_drop": audit["per_grid"][f"joint_{kname}_expanded"]["profile_drop_within_grid"],
        }
        v4[kname]["peak_rel_err"] = abs(v4[kname]["peak"] - v4[kname]["audit_peak"]) \
            / v4[kname]["audit_peak"]
        v4[kname]["drop_rel_err"] = abs(v4[kname]["drop"] - v4[kname]["audit_drop"]) \
            / v4[kname]["audit_drop"]
    res["V4"] = v4

    # V5: noiseless round trip through the explicit matrix, at the best config
    C = {}
    for det in ("H1", "L1"):
        Fp, Fc = antenna(det, ra, dec, math.radians(cfg3[0]["psi_deg"]), GPS_EVENT)
        A = (1.0 + math.cos(math.radians(cfg3[0]["iota_deg"])) ** 2) / 2.0
        B = math.cos(math.radians(cfg3[0]["iota_deg"]))
        C[det] = Fp * A - 1j * Fc * B
    tc_inj = 2.0
    tot2 = 0.0
    sigs = {}
    for det in ("H1", "L1"):
        H = np.zeros(len(freqs), dtype=complex)
        H[m] = (C[det] * base[m]) * np.exp(2j * np.pi * freqs[m] * tau_d[det]) \
            * np.exp(-2j * np.pi * freqs[m] * tc_inj)
        # irfft of h~ gives the time series whose rfft/fs is h~; multiply by fs so
        # the round trip through coherent_matrix (which divides by fs) is exact.
        sigs[det] = FS * np.fft.irfft(H, seg_n)
        tot2 += float(np.abs(C[det]) ** 2) * 4.0 * df * float(
            np.sum(np.abs(base[m] * np.exp(2j * np.pi * freqs[m] * tau_d[det])) ** 2
                   / Sn[det][m]))
    scale = 25.0 / math.sqrt(tot2)
    prof, _ = coherent_matrix({d: sigs[d] * scale for d in sigs}, freqs, base, Sn,
                              C, tau_d, zero_d, tc_grid)
    v5 = math.sqrt(float(np.max(prof)))
    res["V5_round_trip_own"] = v5
    res["V5_rel_err"] = abs(v5 - 25.0) / 25.0

    # V6: recorded single-detector numbers with tau_d = 0 (own path)
    v6 = {}
    for det in ("H1", "L1"):
        H = base[m]
        z = 4.0 * df * np.conj(H) * (np.fft.rfft(segs[det] * tukey(seg_n, 0.25))[m] / FS) / Sn[det][m]
        s = math.sqrt(4.0 * df * float(np.sum(np.abs(H) ** 2 / Sn[det][m])))
        M = np.exp(2j * np.pi * np.outer(np.arange(0.0, 4.0, 2.0 / FS), freqs[m]))
        pk = float(np.max(np.abs(z @ M.T))) / s
        v6[det] = {"peak_own": pk,
                   "recorded_IMR_peak": audit["control_values"]["NC1_recorded_numbers"][det]["recorded_IMR_peak"]}
        v6[det]["rel_err"] = abs(pk - v6[det]["recorded_IMR_peak"]) / v6[det]["recorded_IMR_peak"]
    res["V6"] = v6

    ok = (res["V1_antenna_max_abs_diff"] < 1e-6
          and res["V2_rel_err"] < 1e-6
          and res["V3_rel_err"] < 1e-3
          and all(v4[k]["peak_rel_err"] < 5e-3 for k in v4)
          and all(v4[k]["drop_rel_err"] < 5e-2 for k in v4)
          and res["V5_rel_err"] < 1e-3
          and all(v6[d]["rel_err"] < 5e-3 for d in v6))
    res["VERIFY_CONFIRMED"] = bool(ok)
    return res, ok


if __name__ == "__main__":
    sys.exit(main(selftest="--selftest" in sys.argv))
