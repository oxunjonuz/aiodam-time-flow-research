# Artifacts / Артефакты

Every measurement in this work is frozen in a JSON (or NPZ) artifact under
`work/artifacts/`, written by the script named in the "produced by" column. The
SHA-256 below was computed on the frozen bytes in this repository.

*Каждое измерение заморожено в артефакте в `work/artifacts/`; скрипт-производитель
указан в колонке. Хеши ниже посчитаны по замороженным байтам в этом репозитории.*

| artifact | produced by | SHA-256 |
|---|---|---|
| `analysis_results.json` | `analysis.py` | `ac2594b538b05f0aa1f8229c34b6b188bdd2109ff3aff8eff78b7bb39aec0709` |
| `h3_real_data.json` | `h3_real_data.py` | `00d338e464f8c5e27d4bf785cfbad1639e41d4f201b2ea3cc221a1ddeaf9372d` |
| `h3_joint_fit.json` | `h3_joint_fit.py` | `1daa84078cabac6f146d4586b2c0703738570d6c3126119012e7034537eb84c8` |
| `sky_samples.npz` | `sky_prep.py` | `0afea4a2de3e69c9e37ddd70d7d38400f0f43eeb6ea8059060abbc766f183e74` |
| `h3_imr_check.json` | `h3_imr_check.py` | `36ed90946a221691c26f87bfd827e6c85b20034b0d2a1b270bf91c15a0d6e755` |
| `h3_imr_spins.json` | `h3_imr_spins.py` | `d967c59e9bc839d62434ba8cc20484b4918c3b1abe103925148ee445c1888396` |
| `h3_direct_fit_imr.json` | `h3_direct_fit_imr.py` | `53b40ebcfcaf259101ddbeead3b3ebc6998c2d64cb3f0866171de4730c94578f` |
| `h3_direct_fit.json` | `h3_direct_fit.py` | `bae71da10862ce9300a3f6f0015739e25fd636d381dbddd3f1e16cdbe1c26290` |
| `h3_rho_opt_bands.json` | `h3_rho_opt_bands.py` | `332ee17c425443e7593819901a2e586d7cf9524928843d4053a5dd12a7661273` |
| `h3_grid_bounds.json` | `h3_grid_bounds.py` | `e12bab1f56212246c6f90584ba1af23a3de8e286749da8bdbafb969066b578a0` |
| `pn_orders.json` | `pn_orders.py` | `60589fa39eaa4316e954a41d87d29bf0d4f6b354c6c6ecd87280e4aaf59d9c4d` |
| `et_ce_scaling.json` | `et_ce_scaling.py` | `b01835cd209246004829437a577bbd6fcaad5f8b166a78108c9fc5beab955509` |
| `h4_paper_arithmetic.json` | `h4_paper_arithmetic.py` | `545ce4a6483cc13476073bf3d652b177e56327b7956a887a1823162ff26e5faa` |
| `h4_paper_text_audit.json` | `h4_paper_text_audit.py` | `bce8865bced7098232702567ead1eefbdf611ecc8d9616fa0eab59dfcaf0db17` |
| `h4_weight_audit.json` | `h4_weight_audit.py` | `c4b021601976c42ea45195b81f0867e7bb91cdd6053eb83c10a233ef59edbf0a` |
| `prereg_t8.json` | `prereg_t8.py` | `c3cd7c8a50f112d40713a3caa741e3f25ce0e80e7927413e129fc0df045c36a9` |

Formalization sources / Исходники формализации:

| file | SHA-256 |
|---|---|
| `work/lean/Identifiability.lean` | `01a4827356ea06364e28936edb88b4b9ae18d162204cd2e51a26c5c02690d0e0` |
| `work/lean/Instances.lean` | `5d31382b315c8a0ab5e3c93d4edde003901700e8d91f16d631aa4da8d7363768` |
| `work/lean/RelativityOfNow.lean` | `db150578e647daeff85bd1d7eee2f1efdc8fe3020ded8935282fbb5b43f4991b` |
| `work/lean/PastHypothesis.lean` | `5fc2f2556f8dc482a225032420fd195a47948c1a49f7ac4ddaf8d068bf08cccd` |
| `work/lean/HomogeneousUniverse.lean` | `236cf1d62682920dd50accdc1d0a79e6630b69f1075050e934b61ff094b990c0` |

## What each artifact contains / Что содержит каждый артефакт

| artifact | key numbers / ключевые числа |
|---|---|
| `analysis_results.json` | H1/H2/H3 on the analytic (calibrated) PSD: absorbed fractions `1 − 1.4·10⁻¹⁵` / `1 − 1.7·10⁻¹⁵`, quadratic SNR `0.087`, detectable amplitude `68.8 ms`, controls K1–K9 |
| `h3_real_data.json` | H1/L1 on the measured PSD: ASD(100 Hz) `1.03·10⁻²³` / `9.85·10⁻²⁴`, peaks `7.27` / `5.58`, H1−L1 `0.37 ms`, quadratic SNR `0.174` / `0.127`; 15/15 thresholds |
| `h3_joint_fit.json` | coherent H1+L1 fit: network SNR peak `22.93` (pn) / `23.07` (imr), `22.17` at τ₂ = 0; drops `3.496σ` / `2.352σ` on ±38.4 ms; one-sided crossings `+3.6 ms` / `+2.4 ms`; power control 20.3 ± 5.1 ms |
| `sky_samples.npz` | 402 frozen sky samples from the published LALInference map (RA, Dec, probability) |
| `h3_imr_check.json` | IMR vs inspiral optimal SNR: `31.6595` / `28.7350`; amplitude scale `0.5288` / `0.4470`; the 37.3↔7.27 gap decomposed |
| `h3_imr_spins.json` | spin-envelope scan: ρ_opt range `30.02–32.49` (H1), `27.16–29.54` (L1); inclination multiplier identity |
| `h3_direct_fit_imr.json` | τ₂ fitted to the strain with full IMR: drops `2.925` / `2.144` (pn kernel, ±38.4 ms); power control |
| `h3_direct_fit.json` | τ₂ fitted with the leading inspiral: unconstrained on every grid |
| `h3_rho_opt_bands.json` | band decomposition: `37.3118` (20–300 Hz), `20.7098` (to ISCO); the B2-vs-B3 test |
| `h3_grid_bounds.json` | expanded τ₂ grid ±38.4 ms; the physical phase bound `1.667 ms` (the paper’s 1.2 ms uses 72 % of it) |
| `pn_orders.json` | the identities `∂Ψ/∂t_c = 2πf` and `f·∂Ψ_chirp/∂f = M_c·∂Ψ_chirp/∂M_c` at PN orders 0…3.5, with the negative control |
| `et_ce_scaling.json` | the ET/CE dichotomy: `0.087 → 0.0047` (fixed amplitude) vs `0.087 → 7.625` (fixed τ₂) |
| `h4_paper_arithmetic.json` | the geometric factor `2.0083`; the root scan and its negative control |
| `h4_paper_text_audit.json` | the volume-weighted mean `1.25034`; the refutation of the author’s earlier claim, quoted next to it |
| `h4_weight_audit.json` | the three weights: `1.250339` / `1.285956` / `1.241729` |
| `prereg_t8.json` | the frozen T8 executed as written: 0 false alarms in 100 noise realizations; signal `z = 0.0872` |
