#!/usr/bin/env python3
"""Independent check of the PUBLICATION package against the FROZEN artifacts.

This script does NOT import any audit module. It reads
  * the frozen measurement artifacts in work/artifacts/, BY KEY;
  * the published documents in publication/;
and asserts that every headline number quoted in the paper and in the English
translations is the number actually recorded in the artifact. It also checks the
structural facts a reader would check first: author names present, the four-category
claim ledger present, the scope caveat present, and the metadata files parse.

Its job is to catch a specific class of error that this project has made before:
a document asserting a number that no measurement produced, or a translation that
drifts from the original.

Exit 0 => PUBLICATION_CONSISTENT. Exit 1 => a mismatch, printed.
"""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ART = os.path.join(ROOT, "work", "artifacts")
PUB = os.path.join(ROOT, "publication")

bad = []
ok = []


def check(name, cond, detail=""):
    if cond:
        ok.append(name)
    else:
        bad.append(f"{name}: {detail}")


def load(fn):
    with open(os.path.join(ART, fn), encoding="utf-8") as f:
        return json.load(f)


def read(fn):
    with open(os.path.join(PUB, fn), encoding="utf-8") as f:
        return f.read()


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def ru_num(s):
    """Russian documents use a decimal comma: '0.087' -> also accept '0,087'."""
    return s.replace(".", ",")


def dig(d, *path):
    """Read a nested key defensively; return None if any level is missing."""
    cur = d
    for p in path:
        if isinstance(cur, dict) and p in cur:
            cur = cur[p]
        else:
            return None
    return cur


# ---------------------------------------------------------------- artifacts ---
analysis = load("analysis_results.json")
real = load("h3_real_data.json")
joint = load("h3_joint_fit.json")
imr = load("h3_imr_check.json")
spins = load("h3_imr_spins.json")
dfit = load("h3_direct_fit_imr.json")
bands = load("h3_rho_opt_bands.json")
grid = load("h3_grid_bounds.json")
h4t = load("h4_paper_text_audit.json")
h4w = load("h4_weight_audit.json")
t8 = load("prereg_t8.json")

# --- the numbers the paper quotes, taken from the artifacts BY KEY ------------
n = {}
n["h1_absorbed"] = analysis["H1"]["one_minus_absorbed"]
n["h2_absorbed"] = analysis["H2"]["one_minus_absorbed"]
n["h3_quad_snr_analytic"] = analysis["H3"]["unmodelled_snr_calibrated"]
n["h3_detectable_ms"] = analysis["H3"]["detectable_amplitude_ms_at_snr5"]
n["geo_factor"] = 2.0083  # asserted below against h4

n["asd100_H1"] = real["inputs"]["H1_psd"]["asd_at_100Hz"]
n["asd100_L1"] = real["inputs"]["L1_psd"]["asd_at_100Hz"]
n["peak_snr_H1"] = real["matched_filter"]["H1"]["peak_snr"]
n["peak_snr_L1"] = real["matched_filter"]["L1"]["peak_snr"]
n["h1_minus_l1"] = real["matched_filter"]["H1_minus_L1_ms"]
n["proj_quad_snr_H1"] = real["projection"]["H1"]["H3_quadratic"]["unmodelled_snr"]
n["proj_quad_snr_L1"] = real["projection"]["L1"]["H3_quadratic"]["unmodelled_snr"]
n["proj_det_ms_H1"] = real["projection"]["H1"]["H3_quadratic"]["detectable_amplitude_ms_at_snr5"]
n["proj_det_ms_L1"] = real["projection"]["L1"]["H3_quadratic"]["detectable_amplitude_ms_at_snr5"]

jg = joint["answer"]["joint_fit"]
jl = joint["sky"]
n["net_snr_pn"] = jg["pn"]["net_snr_peak"]
n["net_snr_imr"] = jg["imr"]["net_snr_peak"]
n["drop_pn_expanded"] = jg["pn"]["drop_expanded_grid"]
n["drop_imr_expanded"] = jg["imr"]["drop_expanded_grid"]
n["cross_pn"] = jg["pn"]["one_sigma_high_ms"]
n["cross_imr"] = jg["imr"]["one_sigma_high_ms"]
n["sky_ra"] = jl["ml_ra_deg"]
n["sky_dec"] = jl["ml_dec_deg"]
n["sky_area90"] = jl["area90_deg2"]
n["tau_h1l1"] = dig(joint, "kernels", "tau_H1_minus_tau_L1_ms")
n["imr_h1"] = dig(imr, "bands", "H1", "rho_opt_imr_20_300") or dig(imr, "rho_opt", "H1", "imr")
n["bound_ms"] = dig(grid, "physical_bounds", "phase_sweep_ceiling", "amplitude_ms")
n["geo_factor"] = h4t["H4a"]["geometric_factor"]
n["imr_h1"] = dig(imr, "bands", "H1", "bands", "B_full_20_300", "rho_opt_IMR")
n["imr_l1"] = dig(imr, "bands", "L1", "bands", "B_full_20_300", "rho_opt_IMR")
n["w_full"] = h4w["weights"]["W1_full_proper_volume"]["mean_redshift"]
n["w_created"] = h4w["weights"]["W2_created_volume"]["mean_redshift"]
n["w_flat"] = h4w["weights"]["W3_flat_volume"]["mean_redshift"]
n["t8_false_alarms"] = t8["prereg_T8"]["n_false_alarms"]

# every referenced key must actually have been provided by an artifact
for k, v in n.items():
    check(f"artifact provided {k}", v is not None, "key missing")

# Each entry: (label, token in the ENGLISH doc, token in the RUSSIAN doc,
# automated value, quoted value). The Russian document legitimately uses a decimal
# comma; `ru_num` generates that form automatically when the two tokens are equal.
EXPECT = [
    ("quadratic SNR 0.087 (analytic PSD)", "0.087", "0.087", n["h3_quad_snr_analytic"], 0.087),
    ("detectable amplitude 69 ms (leading order)", "69 ms", "69 мс", n["h3_detectable_ms"], 68.8),
    ("ASD(100 Hz) H1 1.03e-23", "1.03·10⁻²³", "1.03·10⁻²³", n["asd100_H1"], 1.0275194697647844e-23),
    ("ASD(100 Hz) L1 9.85e-24", "9.85·10⁻²⁴", "9.85·10⁻²⁴", n["asd100_L1"], 9.848940548987532e-24),
    ("match-filter peak H1 7.27", "7.27", "7.27", n["peak_snr_H1"], 7.273928125661218),
    ("match-filter peak L1 5.58", "5.58", "5.58", n["peak_snr_L1"], 5.579860999865262),
    ("H1-L1 offset 0.37 ms", "0.37", "0.37", n["h1_minus_l1"], 0.3662109375),
    ("quadratic SNR 0.174 (measured PSD, H1)", "0.174", "0.174", n["proj_quad_snr_H1"], 0.17399876847706267),
    ("quadratic SNR 0.127 (measured PSD, L1)", "0.127", "0.127", n["proj_quad_snr_L1"], 0.12691700655029176),
    ("detectable 34.5 ms (H1)", "34.5", "34.5", n["proj_det_ms_H1"], 34.483002681659485),
    ("detectable 47.3 ms (L1)", "47.3", "47.3", n["proj_det_ms_L1"], 47.27498830207957),
    ("network SNR peak 22.93 (pn)", "22.93", "22.93", n["net_snr_pn"], 22.928036745737085),
    ("network SNR peak 23.07 (imr)", "23.07", "23.07", n["net_snr_imr"], 23.069285482122204),
    ("profile drop 3.496 (pn, expanded)", "3.496", "3.496", n["drop_pn_expanded"], 3.496088347885145),
    ("profile drop 2.352 (imr, expanded)", "2.352", "2.352", n["drop_imr_expanded"], 2.3522000877029754),
    ("sky area 616.4 deg2", "616.4", "616.4", n["sky_area90"], 616.4218405179827),
    ("geometric factor 2.0083", "2.0083", "2.0083", n["geo_factor"], 2.0083004031215705),
    ("volume mean 1.25034", "1.25034", "1.25034", n["w_full"], 1.2503387806021058),
    ("created-volume mean 1.28596", "1.28596", "1.28596", n["w_created"], 1.2859563125961027),
    ("flat-volume mean 1.24173", "1.24173", "1.24173", n["w_flat"], 1.2417289905474065),
    ("physical bound 1.667 ms", "1.667", "1,667", n["bound_ms"], 1.6666666666666665),
    ("IMR optimal SNR 31.66 (H1)", "31.66", "31,66", n["imr_h1"], 31.659512701612417),
]

EN = read("paper/paper_en.md")
RU = read("paper/paper_ru.md")
REP = read("docs_en/REPORT_en.md")
CLM = read("docs_en/CLAIMS_en.md")
NTS = read("docs_en/NOTES_en.md")
RND = read("docs_en/ROUNDS_en.md")
RDME = read("README.md")
ARTF = read("ARTIFACTS.md")
docs = {"paper_en": EN, "paper_ru": RU, "REPORT_en": REP, "CLAIMS_en": CLM,
        "NOTES_en": NTS, "ROUNDS_en": RND, "README": RDME, "ARTIFACTS": ARTF}

# The paper (EN and RU) must carry each headline number, and the artifact's own
# value must be consistent with the rounded figure the document quotes. Each
# language gets its own token because the Russian document legitimately uses a
# decimal comma and Russian unit words; that is a language convention, not a
# different number.
for label, en_tok, ru_tok, val, target in EXPECT:
    check(f"paper_en has {label}", en_tok in EN, f"token {en_tok!r} not found")
    found_ru = (ru_tok in RU) or (ru_num(ru_tok) in RU)
    check(f"paper_ru has {label}", found_ru,
          f"token {ru_tok!r} not found (nor its comma form {ru_num(ru_tok)!r})")
    check(f"artifact {label} agrees with the quoted figure",
          rel(val, target) < 5e-3, f"artifact {val!r} vs quoted {target!r}")

# --------------------------------------------------------------- structural ---
for dname in ("paper_en", "paper_ru"):
    check(f"{dname} names Aiodam", "Aiodam" in docs[dname], "missing")
    check(f"{dname} names Ubaydullayev", "Ubaydullayev" in docs[dname], "missing")
    # The disclosure must be a STATEMENT, not a bare phrase. An earlier version of
    # this check looked for "autonomous research agent" anywhere in the document —
    # but that string occurs TWICE in paper_en.md (once in the YAML author list,
    # once in the disclosure paragraph), so deleting the disclosure left the check
    # green. Measured, not assumed: a probe that removed the paragraph stayed green.
    # The check now requires the full sentence, so there is nothing left to hide
    # behind a duplicate.
    if dname == "paper_en":
        check(f"{dname} states the authorship disclosure",
              "carried out **by an autonomous research agent, Aiodam**, acting as" in docs[dname],
              "missing the full agent-authorship sentence")
    else:
        check(f"{dname} states the authorship disclosure",
              "Работа выполнена **автономным исследовательским агентом Aiodam**" in docs[dname],
              "missing the full agent-authorship sentence")
    # and the error-count disclosure, which is now part of the authorship claim
    if dname == "paper_en":
        check(f"{dname} states who found the errors",
              "found **eleven** of them" in docs[dname] and "count_errors.py" in docs[dname],
              "missing the error-attribution sentence")
    else:
        check(f"{dname} states who found the errors",
              "нашёл **одиннадцать**" in docs[dname] and "count_errors.py" in docs[dname],
              "missing the error-attribution sentence")

# the scope caveat must appear in the paper (EN and RU) and in REPORT_en.
# Match on NORMALIZED text (markdown asterisks and line wraps removed) so the
# check tests the presence of the statement, not its formatting.
def norm(s):
    return " ".join(s.replace("*", " ").replace("`", " ").split())


for dname in ("paper_en", "paper_ru", "REPORT_en"):
    t = norm(docs[dname])
    has = ("not derived" in t and "does not specify" in t) or \
          ("не выведена" in t and "не задаёт" in t) or \
          ("не выведен" in t and "не задаёт" in t)
    check(f"{dname} carries the scope caveat", has, "missing 'tau(t) not derived'")

# four-category ledger present in CLAIMS_en
for cat in ("PROVEN", "DERIVED symbolically", "MEASURED", "HYPOTHESIS"):
    check(f"CLAIMS_en has category {cat}", cat in CLM, "missing")

# the EN report must be a real translation of comparable size to the RU original
ru_rep = open(os.path.join(ROOT, "work", "REPORT.md"), encoding="utf-8").read()
check("REPORT_en is non-trivial", len(REP) > 15000, f"len={len(REP)}")
check("REPORT_en size comparable to RU original",
      0.4 < len(REP) / len(ru_rep) < 2.2, f"{len(REP)}/{len(ru_rep)}")

# every English translation must be non-empty
for k, v in docs.items():
    check(f"{k} non-empty", len(v) > 500, f"len={len(v)}")

# ------------------------------------------------- self-contradiction checks ---
# The owner found four places where the paper contradicted itself in its own
# headline caveat. Each is now a machine check, so it cannot come back by editing.
# (1) the figure caption must not call the tau2 profile a constraint while the
#     text beside it calls it an upper bound;
# (2) the inclination must not be called measured while it is a fit;
# (3) the 69 ms vs 0.61 ms comparison must not be called three orders of magnitude
#     (it is 113x = two orders), and the SNR numbers must be attributed;
# (4) §9 must not carry the same limitation twice.

for dname, lang in (("paper_en", "en"), ("paper_ru", "ru")):
    t = norm(docs[dname])
    # (1) the caption may state the drops, but not claim a measured constraint
    if lang == "en":
        check("caption does not overclaim the tau2 bound",
              "first constraint on" not in t,
              "caption still says 'first constraint'")
        check("caption states the upper-bound caveat",
              "upper bound on the strength of the constraint" in t,
              "caption lacks the upper-bound caveat")
        # (2) inclination
        check("inclination is not called measured",
              "orientation of the source, measured at" not in t,
              "text still says the orientation was measured")
        check("inclination is called fitted",
              "fitted inclination" in t or "fitted" in t,
              "text does not say the inclination is fitted")
        # (3) orders of magnitude and SNR attribution
        check("69 ms vs 0.61 ms is two orders, not three",
              "two orders of magnitude above" in t and "three orders of magnitude" not in t,
              "order-of-magnitude claim wrong")
        check("the three SNR numbers are attributed",
              "Which signal-to-noise ratio is which" in t,
              "no SNR provenance paragraph")
        check("the fitted inclination is a listed limitation",
              "is fitted, not measured" in t,
              "limitation 9 missing")
    else:
        check("caption does not overclaim the tau2 bound",
              "первое ограничение" not in t,
              "caption still says 'первое ограничение'")
        check("caption states the upper-bound caveat",
              "оценка силы ограничения сверху" in t,
              "caption lacks the upper-bound caveat")
        check("inclination is not called measured",
              "ориентация источника, измеренная при" not in t,
              "text still says the orientation was measured")
        check("inclination is called fitted",
              "подогнанное наклонение" in t,
              "text does not say the inclination is fitted")
        check("69 ms vs 0.61 ms is two orders, not three",
              "на два порядка выше" in t and "на три порядка выше" not in t,
              "order-of-magnitude claim wrong")
        check("the three SNR numbers are attributed",
              "Какое С/Ш какое" in t,
              "no SNR provenance paragraph")
        check("the fitted inclination is a listed limitation",
              "подогнанное, а не измеренное" in t,
              "limitation 9 missing")

# (4) no limitation may appear twice in §9 — count the numbered items and require
#     the limitation texts to be unique
for dname in ("paper_en", "paper_ru"):
    src = docs[dname]
    m = re.search(r"# 9\..*?(?=\n# 10\.)", src, re.S)
    check(f"{dname} has a §9 limitations block", m is not None, "no §9")
    if not m:
        continue
    body = m.group(0)
    items = re.findall(r"^\s*(\d+)\.\s", body, re.M)
    check(f"{dname} §9 items are numbered without repeats",
          len(items) == len(set(items)), f"numbers={items}")
    # each item's first 40 normalised characters must be unique
    starts = [norm(x)[:40] for x in re.split(r"^\s*\d+\.\s", body, flags=re.M)[1:]]
    check(f"{dname} §9 has no duplicated limitation",
          len(starts) == len(set(starts)), f"duplicates in {starts}")

# (5) §2 must not restrict the headline to register 1. That was a self-
#     contradiction: §2 said "only register 1 may enter the headline" while §8
#     said the headline rests on the formal results A1–A3 AND the symbolic
#     identities A4–A6 (register 2). The ledger (CLAIMS.md §F) settles it: the
#     headline rests on registers 1 and 2, and the excluded category is
#     [HYPOTHESIS]. The two statements must agree.
for dname, lang in (("paper_en", "en"), ("paper_ru", "ru")):
    t = norm(docs[dname])
    if lang == "en":
        check("§2 does not restrict the headline to register 1",
              "only allowed into the headline if it is in register 1" not in t,
              "§2 still says only register 1 may enter the headline")
        check("§2 states the headline rests on registers 1 and 2",
              "rests on registers 1 and 2" in t,
              "§2 does not say the headline rests on registers 1 and 2")
        check("§2 excludes the hypothesis category from the headline",
              "hypothesis category is never allowed into the headline" in t,
              "§2 does not name the excluded category")
    else:
        check("§2 does not restrict the headline to register 1",
              "в заголовок допускается только утверждение из регистра 1" not in t,
              "§2 still says only register 1 may enter the headline")
        check("§2 states the headline rests on registers 1 and 2",
              "держится только на регистрах 1 и 2" in t,
              "§2 does not say the headline rests on registers 1 and 2")
        check("§2 excludes the hypothesis category from the headline",
              "не допускается ни одно утверждение категории «гипотеза»" in t,
              "§2 does not name the excluded category")

# (6) §3.2 announces a count of formal results and then lists them. The number
#     word must match the number of bullets — "Two further" in front of three
#     bullets is the same class of error as a caption stronger than its figure.
_COUNT_WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}
_COUNT_WORD_RU = {1: "один", 2: "два", 3: "три", 4: "четыре", 5: "пять"}
for dname, lang in (("paper_en", "en"), ("paper_ru", "ru")):
    src = docs[dname]
    m = re.search(r"## 3\.2.*?(?=\n# 4\.)", src, re.S)
    check(f"{dname} has a §3.2 block", m is not None, "no §3.2")
    if not m:
        continue
    body = m.group(0)
    n_bullets = len(re.findall(r"^\s*\*\s", body, re.M))
    if lang == "en":
        lead = re.search(r"^(\w+) further formal results", body, re.M)
        word = lead.group(1).lower() if lead else None
        expected = _COUNT_WORD.get(n_bullets)
    else:
        lead = re.search(r"^(\w+) дальнейших формальных результата", body, re.M)
        word = lead.group(1).lower() if lead else None
        expected = _COUNT_WORD_RU.get(n_bullets)
    check(f"{dname} §3.2 count word matches its {n_bullets} bullets",
          word is not None and word == expected,
          f"word={word!r} bullets={n_bullets} expected={expected!r}")

# (8) the "four places" clause must NOT appear in any reader-facing document. The
#     owner asked (round 232) that it be removed: it changes nothing for the
#     reader, and a reader who sees the paper contradicting itself may simply
#     close the PDF. The count of four errors stays in the internal record
#     (NOTES.md, count_errors.py) — that is the honest log — but the wording that
#     advertises a self-contradiction does not belong in the paper, the README or
#     the deposit metadata. This is a machine check so the clause cannot return by
#     an edit that no test looks at.
_REMOVED_CLAUSE = [
    "contradicting its own",
    "противоречит собственной",
    "in four places",
    "в четырёх местах",
]
for _name in ("paper_en", "paper_ru", "README", "zenodo"):
    _txt = read(".zenodo.json") if _name == "zenodo" else docs[_name]
    _hit = [c for c in _REMOVED_CLAUSE if c in _txt]
    check(f"{_name} does not advertise the removed self-contradiction clause",
          not _hit, f"still present: {_hit}")

# (9) NO reader-facing document may say that errors were found while BUILDING or
#     PROOFREADING the PDF, or count the reviews the errors came from. The owner
#     asked for this twice — round 232 ("it changes nothing for the reader") and
#     round 233 ("there is no need at all to write about errors made and found
#     while the PDF was being built"). The number of errors and who found them
#     stays in the internal record (NOTES.md, count_errors.py); the reader gets the
#     split and nothing about how the proofreading went. A machine check, because
#     the same sentence has now been removed by hand twice.
_PDF_STAGE_CLAUSE = [
    "review of the built PDF",
    "built PDF and found",
    "reviews the errors",
    "in three reviews",
    "вычитке собранного PDF",
    "собранного PDF",
    "в трёх вычитках",
    "дефект вёрстки",
    "typesetting",
]
for _name in ("paper_en", "paper_ru", "README", "zenodo"):
    _txt = read(".zenodo.json") if _name == "zenodo" else docs[_name]
    _hit = [c for c in _PDF_STAGE_CLAUSE if c in _txt]
    check(f"{_name} says nothing about errors found while building the PDF",
          not _hit, f"still present: {_hit}")

# (7) §0 quotes an error split (co-author / agent, and a number of reviews). That
#     split must be the one count_errors.py produces from NOTES.md — otherwise the
#     paragraph is a remembered number, which is exactly what it claims not to be.
#     The counter is run as a subprocess so this check does not import its logic.
import subprocess as _sp

_ce = os.path.join(PUB, "count_errors.py")
if os.path.exists(_ce):
    _p = _sp.run([sys.executable, _ce], capture_output=True, text=True)
    _out = _p.stdout
    _m_owner = re.search(r"found by the co-author:\s*(\d+)", _out)
    _m_agent = re.search(r"found by the agent:\s*(\d+)", _out)
    _m_total = re.search(r"TOTAL recorded errors:\s*(\d+)", _out)
    check("count_errors.py ran and is self-consistent",
          _p.returncode == 0 and _m_owner and _m_agent and _m_total,
          f"rc={_p.returncode} out={_out[-200:]}")
    if _m_owner and _m_agent and _m_total:
        _ow, _ag, _tt = int(_m_owner.group(1)), int(_m_agent.group(1)), int(_m_total.group(1))
        _words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
                  7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
                  12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen"}
        _words_ru = {1: "один", 2: "два", 3: "три", 4: "четыре", 5: "пять",
                     6: "шесть", 7: "семь", 8: "восемь", 9: "девять", 10: "десять",
                     11: "одиннадцать", 12: "двенадцать"}
        _en, _ru = norm(EN), norm(RU)
        check("paper_en §0 co-author count equals count_errors.py",
              f"found {_words[_ow]} of them" in _en, f"expected {_words[_ow]}")
        check("paper_en §0 agent count equals count_errors.py",
              f"remaining {_ag} were found by the agent" in _en, f"expected {_ag}")
        check("paper_ru §0 co-author count equals count_errors.py",
              f"нашёл {_words_ru[_ow]} из них" in _ru, f"expected {_words_ru[_ow]}")
        check("paper_ru §0 agent count equals count_errors.py",
              f"Остальные {_ag} нашли" in _ru, f"expected {_ag}")
        check("paper_en §0 total equals count_errors.py",
              f"records {_tt} errors made on the way" in _en, f"expected {_tt}")
        check("paper_ru §0 total equals count_errors.py",
              f"записано {_tt} ошибок" in _ru, f"expected {_tt}")
        # the typesetting defect is PDF assembly, not a research error, and it is
        # not among the 47 — but that exclusion must NOT be argued in the reader
        # documents. The owner (round 233) asked that the disclosure of errors
        # found while building and proofreading the PDF be kept out of the paper
        # entirely: it changes nothing for the reader and only invites mockery.
        # So the §0 paragraph now states the split and stops. This negative check
        # replaces the earlier positive one that required the clause to be there.
        _typ = ["defect of PDF assembly", "not among the 47",
                "дефект вёрстки PDF", "в число 47 она не входит"]
        check("paper_en §0 does not argue the typesetting exclusion",
              not any(c in _en for c in _typ[:2]),
              "typesetting-exclusion clause is back in §0")
        check("paper_ru §0 does not argue the typesetting exclusion",
              not any(c in _ru for c in _typ[2:]),
              "typesetting-exclusion clause is back in §0")
        # the deposit metadata carries the same split; it must not drift from it
        _zen = open(os.path.join(PUB, ".zenodo.json"), encoding="utf-8").read()
        check("zenodo description carries the same co-author count",
              f"<strong>{_words[_ow]}</strong> of the {_tt} recorded errors" in _zen,
              f"expected {_words[_ow]} of {_tt}")
        check("zenodo description carries the same agent count",
              f"remaining <strong>{_ag}</strong>" in _zen, f"expected {_ag}")
        check("zenodo description does not argue the typesetting exclusion",
              "not among the 47" not in _zen,
              "typesetting-exclusion clause is back in the deposit metadata")

# ------------------------------------------------------------------ metadata ---
try:
    json.load(open(os.path.join(PUB, ".zenodo.json"), encoding="utf-8"))
    check(".zenodo.json parses", True)
except Exception as e:
    check(".zenodo.json parses", False, str(e))

try:
    import yaml  # optional
    yaml.safe_load(open(os.path.join(PUB, "CITATION.cff"), encoding="utf-8"))
    check("CITATION.cff parses", True)
except ImportError:
    # no yaml available: do a light structural check instead
    cff = open(os.path.join(PUB, "CITATION.cff"), encoding="utf-8").read()
    check("CITATION.cff looks structural (no yaml module)",
          cff.startswith("cff-version") and "authors:" in cff and "title:" in cff)
except Exception as e:
    check("CITATION.cff parses", False, str(e))

# LICENSE files present and non-empty
for fn in ("LICENSE", "LICENSE-CODE"):
    p = os.path.join(PUB, fn)
    check(f"{fn} present and non-empty",
          os.path.exists(p) and os.path.getsize(p) > 200,
          "missing or tiny")

# ARTIFACTS.md must quote the sha256 of EVERY artifact actually on disk (not just
# one): the index is a claim about the frozen bytes, so it is checked against them.
import hashlib
import re as _re

listed = set(_re.findall(r"`([0-9a-f]{64})`", ARTF))
disk_hashes = {}
for _name in sorted(os.listdir(ART)):
    with open(os.path.join(ART, _name), "rb") as f:
        disk_hashes[_name] = hashlib.sha256(f.read()).hexdigest()

missing_from_index = {k: v for k, v in disk_hashes.items() if v not in listed}
check("ARTIFACTS.md lists every artifact hash on disk",
      not missing_from_index,
      f"missing {sorted(missing_from_index)[:5]}")

# and every hash it DOES list must be a real artifact hash (no stale entries)
stale = [h for h in listed if h not in set(disk_hashes.values())]
# the Lean-source hashes in ARTIFACTS.md are also 64-hex, so allow those
lean_hashes = set()
lean_dir = os.path.join(ROOT, "work", "lean")
if os.path.isdir(lean_dir):
    for _name in os.listdir(lean_dir):
        if _name.endswith(".lean"):
            with open(os.path.join(lean_dir, _name), "rb") as f:
                lean_hashes.add(hashlib.sha256(f.read()).hexdigest())
stale = [h for h in stale if h not in lean_hashes]
check("ARTIFACTS.md has no stale hashes", not stale, f"stale={stale[:3]}")

# the specific joint-fit hash must be present (kept as an explicit named check)
_h = hashlib.sha256(open(os.path.join(ART, "h3_joint_fit.json"), "rb").read()).hexdigest()
check("ARTIFACTS.md quotes the real h3_joint_fit sha256", _h in ARTF, f"disk={_h}")

# MANIFEST.md is a claim about the published bytes, so it is checked against them:
# every file it lists must hash to what it says, and every published file must be
# listed. Without this, adding a file (e.g. the figures) would silently leave the
# manifest incomplete.
MAN = read("MANIFEST.md")
man_rows = _re.findall(r"\|\s*`([^`]+)`\s*\|\s*(\d+)\s*\|\s*`([0-9a-f]{64})`\s*\|", MAN)
check("MANIFEST.md has rows", len(man_rows) > 5, f"rows={len(man_rows)}")
man_bad = []
for _rel, _size, _h in man_rows:
    p = os.path.join(PUB, _rel)
    if not os.path.exists(p):
        man_bad.append(f"{_rel}: missing")
        continue
    with open(p, "rb") as f:
        b = f.read()
    if hashlib.sha256(b).hexdigest() != _h:
        man_bad.append(f"{_rel}: hash mismatch")
    elif len(b) != int(_size):
        man_bad.append(f"{_rel}: size {len(b)} != {_size}")
check("MANIFEST.md matches the bytes on disk", not man_bad, f"{man_bad[:4]}")

listed = {r for r, _, _ in man_rows}
on_disk = set()
for d in ("", "paper", "docs_en", "figures"):
    folder = os.path.join(PUB, d) if d else PUB
    if not os.path.isdir(folder):
        continue
    for name in os.listdir(folder):
        if not os.path.isfile(os.path.join(folder, name)):
            continue
        _relp = os.path.relpath(os.path.join(folder, name), PUB).replace(os.sep, "/")
        if _relp in ("MANIFEST.md", "make_manifest.py"):
            continue
        on_disk.add(_relp)
missing_from_manifest = on_disk - listed
check("MANIFEST.md lists every published file", not missing_from_manifest,
      f"missing {sorted(missing_from_manifest)[:5]}")
check("MANIFEST.md lists the figures", any(r.startswith("figures/") for r in listed),
      "no figures in the manifest")

# ------------------------------------------------------- Lean traceability ---
# The paper names specific Lean theorems. If a name is wrong, the paper claims a
# proof that does not exist — the exact "headline stronger than content" failure.
# So each named theorem is checked to be DEFINED in work/lean/*.lean.
LEAN_THEOREMS = [
    "identifiable_iff",
    "not_identifiable_of_range_le",
    "demo_identifiable",
    "demo_not_identifiable",
    "exists_boost_reversing_time_order",
    "arrow_sum_zero",
    "demo_closure_is_the_line",
    "paramArrow_ne_zero",
    "rev_isTraj",
    "arrow_vanishes_on_symmetric",
    "rev_allowed_of_invariant",
    "Bfuture_not_invariant",
    "demo_boundary_breaks_reversal",
    "creation_not_identifiable_unpinned",
    "creation_identifiable_pinned",
    "baseline_unique_when_pinned",
    "demo_unpinned_pair",
    "demo_pinned_identifiable",
]
lean_dir = os.path.join(ROOT, "work", "lean")
lean_text = ""
if os.path.isdir(lean_dir):
    for _n in os.listdir(lean_dir):
        if _n.endswith(".lean"):
            with open(os.path.join(lean_dir, _n), encoding="utf-8") as f:
                lean_text += f.read()

defined = set(re.findall(r"^(?:theorem|lemma|def)\s+([A-Za-z0-9_']+)",
                         lean_text, re.M))

ALL_DOCS = "\n".join(docs.values())

# (1) every theorem the PAPER names must be DEFINED in Lean
PAPER_NAMED = [
    "identifiable_iff",
    "not_identifiable_of_range_le",
    "demo_identifiable",
    "demo_not_identifiable",
    "exists_boost_reversing_time_order",
    "arrow_sum_zero",
    "demo_closure_is_the_line",
]
for thm in PAPER_NAMED:
    check(f"theorem named in the paper exists in Lean: {thm}",
          (thm in defined) and (thm in EN),
          f"defined={thm in defined} named_in_paper={thm in EN}")

# (2) the remaining theorems are cited in the companion documents, and all must
#     really exist (a cited name that is not defined would be a phantom proof)
COMPANION = [t for t in LEAN_THEOREMS if t not in PAPER_NAMED]
for thm in COMPANION:
    check(f"theorem {thm} exists in Lean and is cited somewhere in the package",
          (thm in defined) and (thm in ALL_DOCS),
          f"defined={thm in defined} cited={thm in ALL_DOCS}")

# (3) no theorem may be cited anywhere with a name that Lean does not define —
#     i.e. every name in our list must be defined (already covered above), and
#     the claims file must NOT imply a Lean proof for the symbolic claims
check("CLAIMS_en states A4-A6 have no Lean file",
      "no lean file" in CLM.lower(),
      "missing the 'no Lean file' statement")
check("CLAIMS_en states the narrowed claim 'in Lean A1-A3'",
      ("A1–A3" in CLM) or ("A1-A3" in CLM), "missing the narrowed claim")

# the reported physical bound must be in the artifact AND quoted in the paper
check("artifact provided physical bound (1.667 ms)",
      n["bound_ms"] is not None and rel(n["bound_ms"], 1.6666666666666665) < 5e-3,
      f"bound={n['bound_ms']}")
check("1.667 ms bound quoted in paper_en", "1.667" in EN)
check("1.667 ms bound quoted in paper_ru", "1,667" in RU)
check("1.667 ms bound present in ROUNDS_en", "1.667" in RND)

# The bound depends on WHICH delay convention is used, by a factor of 15. Both
# must exist in the shipment: the artifact value (the paper's accumulating-delay
# convention) AND the disclosure that the opposite convention gives 24.96 ms.
# Re-derive both here from the physics, independently of the artifact, and
# require the documents to carry the caveat. (This check exists because an
# independent re-derivation once returned 24.96 ms and the artefact was right.)
import math as _math

_MC = 28.096 * (6.67430e-11 * 1.98892e30 / 2.99792458e8 ** 3)


def _t_of_f(f):
    return -(5.0 / 256.0) * _MC ** (-5.0 / 3.0) * (_math.pi * f) ** (-8.0 / 3.0)


_fmin, _fmax = 20.0, 300.0
_tmin, _tmax = _t_of_f(_fmin), _t_of_f(_fmax)
_span = abs(_tmax - _tmin)
# convention A: delay vanishes at the band start, grows to merger (the paper's)
_kern_a = max(abs(2 * _math.pi * _fmin * (_tmin - _tmin) ** 2),
              abs(2 * _math.pi * _fmax * (_tmax - _tmin) ** 2))
_amp_a_ms = (_math.pi / _kern_a) * _span ** 2 * 1e3
# convention B: vertex at coalescence (delay vanishes AT merger)
_kern_b = max(abs(2 * _math.pi * _fmin * _tmin ** 2),
              abs(2 * _math.pi * _fmax * _tmax ** 2))
_amp_b_ms = (_math.pi / _kern_b) * _span ** 2 * 1e3

check("independent re-derivation reproduces the artifact's 1.667 ms",
      rel(_amp_a_ms, 1.6666666666666665) < 5e-3, f"re-derived {_amp_a_ms:.4f} ms")
check("independent re-derivation reproduces the 24.96 ms opposite-convention bound",
      rel(_amp_b_ms, 24.9635) < 5e-3, f"re-derived {_amp_b_ms:.4f} ms")
check("both conventions differ by ~15x (the caveat is not cosmetic)",
      13.0 < _amp_b_ms / _amp_a_ms < 17.0, f"ratio {_amp_b_ms / _amp_a_ms:.2f}")
check("paper_en discloses the convention dependence",
      "24.96" in EN and "convention" in EN.lower(),
      "missing the convention caveat")
check("paper_ru discloses the convention dependence",
      ("24,96" in RU) and ("конвенц" in RU.lower()),
      "missing the convention caveat")

# ------------------------------------------------------------------- report ---
print("PUBLICATION CHECK")
print("-" * 60)
for o in ok:
    print("  ok    ", o)
for b in bad:
    print("  FAIL  ", b)
print("-" * 60)
print(f"{len(ok)} ok, {len(bad)} failed")

if bad:
    print("VERDICT: PUBLICATION_INCONSISTENT")
    sys.exit(1)
print("VERDICT: PUBLICATION_CONSISTENT")
sys.exit(0)
