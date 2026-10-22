"""
Electrical design study for a 415 V utilities MCC fed by an 11/0.433 kV transformer.

    11 kV bus --[R3, VCB]-- 11/0.433 kV Dyn11 transformer --[R2, ACB]-- 415 V MCC bus --[R1]-- cooling-tower sub-MCC

Steps (each one is a few lines of textbook formula):
  1. Load list -> running kW, kVAR, kVA and plant power factor
  2. Power-factor correction -> capacitor bank kVAR, standard step size
  3. Transformer sizing before / after correction
  4. Three-phase fault level at the 415 V bus (per-unit method)
  5. IDMT overcurrent relay settings and grading (R1 -> R2 -> R3)

Run:  python src/study.py        (writes results/study_report.md and 3 plots)
"""
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

# ---------------- Design inputs (change these and re-run) ----------------
CFG = {
    "v_lv": 415.0,                  # V, MCC bus line voltage
    "v_lv_noload": 433.0,           # V, transformer secondary no-load voltage
    "v_hv": 11000.0,                # V
    "target_pf": 0.95,
    "cap_step_kvar": 25,            # APFC panel step size
    "future_margin": 0.20,          # 20 % spare capacity for future load
    "std_tx_kva": [315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500],
    "source_fault_mva": 250.0,      # 11 kV bus fault level (assumed)
    "motor_contribution_x": 4.0,    # running motors feed ~4 x FLC into a fault (first cycles)
    "std_breaking_ka": [25, 36, 50, 65, 80, 100],
    "tx_copper_loss_pct": 1.0,      # full-load copper loss ~1 % of rating (typical distribution transformer)
    "hours_per_year": 7200,         # continuous process plant
    "tariff_kva_month": 475.0,      # Rs/kVA/month demand charge   (illustrative)
    "tariff_kvah": 7.5,             # Rs/kVAh energy charge        (illustrative, kVAh billing)
    "grading_margin_s": 0.30,       # breaker time + relay overshoot + safety
    "psm_cap": 20.0,                # IDMT curves are defined up to PSM 20
    "r1_feeder_tag": "F-201",       # load fed by the outgoing feeder R1
    "motor_start_x": 6.0,           # DOL starting current = 6 x FLC
    "motor_start_time_s": 8.0,      # fan run-up time (high inertia)
}

# IEC 60255 / IEC 60255-151 standard curves:  t = TMS * k / (PSM^a - 1)
CURVES = {"SI": (0.14, 0.02), "VI": (13.5, 1.0), "EI": (80.0, 2.0)}


def idmt_time(curve, tms, current, pickup):
    k, a = CURVES[curve]
    psm = np.minimum(np.asarray(current, dtype=float) / pickup, CFG["psm_cap"])
    with np.errstate(divide="ignore", invalid="ignore"):
        t = tms * k / (psm ** a - 1.0)
    return np.where(psm > 1.0, t, np.inf)


def round_up(x, step):
    return math.ceil(x / step - 1e-9) * step


def amps(kva, volts):
    return kva * 1000.0 / (math.sqrt(3) * volts)


# ---------------- Step 1: load list ----------------
def load_summary(df):
    df = df.copy()
    df["kw"] = df.rated_kw * df.load_factor / df.efficiency * df.qty_running
    df["kvar"] = df.kw * np.tan(np.arccos(df.power_factor))
    df["kva"] = np.hypot(df.kw, df.kvar)
    df["flc_a"] = df.rated_kw * 1000 / (math.sqrt(3) * CFG["v_lv"] * df.efficiency * df.power_factor)
    return df


# ---------------- Step 2: power-factor correction ----------------
def pf_correction(p_kw, q_kvar):
    pf1 = p_kw / math.hypot(p_kw, q_kvar)
    qc_needed = p_kw * (q_kvar / p_kw - math.tan(math.acos(CFG["target_pf"])))
    qc = round_up(qc_needed, CFG["cap_step_kvar"])
    q2 = q_kvar - qc
    s2 = math.hypot(p_kw, q2)
    return {"pf_before": pf1, "qc_needed": qc_needed, "qc_selected": qc,
            "steps": int(qc / CFG["cap_step_kvar"]), "q_after": q2, "s_after": s2, "pf_after": p_kw / s2}


def pick_transformer(kva):
    need = kva * (1 + CFG["future_margin"])
    return need, next(r for r in CFG["std_tx_kva"] if r >= need)


def tx_impedance_pct(kva):
    # Typical %Z for distribution transformers (IS 1180 / IS 2026 practice)
    return 4.5 if kva <= 630 else 5.0 if kva <= 1250 else 6.25


# ---------------- Step 4: fault level ----------------
def fault_level(tx_kva, motor_flc_total):
    base_mva = tx_kva / 1000.0
    z_src = base_mva / CFG["source_fault_mva"]          # pu on transformer base
    z_tx = tx_impedance_pct(tx_kva) / 100.0
    z_tot = z_src + z_tx
    i_base = amps(tx_kva, CFG["v_lv_noload"])
    i_grid = i_base / z_tot
    i_motor = CFG["motor_contribution_x"] * motor_flc_total
    return {"z_src": z_src, "z_tx": z_tx, "z_tot": z_tot, "i_base": i_base, "fault_mva": base_mva / z_tot,
            "i_fault_grid": i_grid, "i_motor": i_motor, "i_fault_total": i_grid + i_motor}


# ---------------- Step 5: relay coordination ----------------
def plug_setting(load_a, ct_primary, min_pickup_a=0.0):
    """Smallest plug setting (0.5-2.0 x CT rating, 0.05 steps as on a numerical relay) that is above
    1.25 x load current AND above min_pickup_a (an upstream relay must pick up above the downstream one)."""
    need = max(1.25 * load_a, min_pickup_a) / ct_primary
    return min(max(round_up(need, 0.05), 0.5), 2.0)


def min_tms_for_margin(curve, pickup, downstream, i_max):
    """Smallest TMS (0.01 steps) that stays >= margin slower than `downstream` at every current up to i_max."""
    i_grid = np.geomspace(pickup * 1.05, i_max, 400)
    t_down = downstream(i_grid)
    ok = np.isfinite(t_down)
    unit = idmt_time(curve, 1.0, i_grid[ok], pickup)
    need = np.max((t_down[ok] + CFG["grading_margin_s"]) / unit)
    return round_up(max(need, 0.05), 0.01)


def coordinate(loads, tx_kva, fault):
    ratio = CFG["v_lv_noload"] / CFG["v_hv"]                # LV amps -> HV amps
    i_fault = fault["i_fault_grid"]                         # coordinate on the grid in-feed (worst case for R3)

    feeder = loads[loads.tag == CFG["r1_feeder_tag"]].iloc[0]
    r1_load = amps(feeder.kva, CFG["v_lv"])
    r1_ct = 300
    r1_pickup = plug_setting(r1_load, r1_ct) * r1_ct
    # Motor-start check: one fan starts DOL while the others run
    one_flc = feeder.flc_a
    i_start = CFG["motor_start_x"] * one_flc + (feeder.qty_running - 1) * r1_load / feeder.qty_running
    unit = idmt_time("EI", 1.0, i_start, r1_pickup)
    r1_tms = round_up(max(1.25 * CFG["motor_start_time_s"] / unit, 0.05), 0.01)

    r2_ct = 1250
    r2_load = amps(tx_kva, CFG["v_lv"])
    r2_pickup = plug_setting(r2_load, r2_ct) * r2_ct
    r1 = lambda i: idmt_time("EI", r1_tms, i, r1_pickup)
    r2_tms = min_tms_for_margin("SI", r2_pickup, r1, i_fault)

    r3_ct = 75
    r3_load_hv = amps(tx_kva, CFG["v_hv"])
    r3_pickup_hv = plug_setting(r3_load_hv, r3_ct, min_pickup_a=1.1 * r2_pickup * ratio) * r3_ct
    r3_pickup_lv = r3_pickup_hv / ratio                     # same relay, seen from the LV side
    r2 = lambda i: idmt_time("SI", r2_tms, i, r2_pickup)
    r3_tms = min_tms_for_margin("SI", r3_pickup_lv, r2, i_fault)
    r3_inst_hv = 1.3 * i_fault * ratio                      # high-set: above max LV through-fault

    relays = [
        {"name": "R1", "where": "415 V outgoing feeder to cooling-tower MCC", "curve": "EI", "ct": f"{r1_ct}/1",
         "pickup_a": r1_pickup, "ps": r1_pickup / r1_ct, "tms": r1_tms, "pickup_lv": r1_pickup},
        {"name": "R2", "where": "415 V transformer incomer (ACB)", "curve": "SI", "ct": f"{r2_ct}/1",
         "pickup_a": r2_pickup, "ps": r2_pickup / r2_ct, "tms": r2_tms, "pickup_lv": r2_pickup},
        {"name": "R3", "where": "11 kV transformer feeder (VCB)", "curve": "SI", "ct": f"{r3_ct}/1",
         "pickup_a": r3_pickup_hv, "ps": r3_pickup_hv / r3_ct, "tms": r3_tms, "pickup_lv": r3_pickup_lv},
    ]
    for r in relays:
        r["t_at_fault"] = float(idmt_time(r["curve"], r["tms"], i_fault, r["pickup_lv"]))
    start = {"i_start": i_start, "t_r1_at_start": float(idmt_time("EI", r1_tms, i_start, r1_pickup))}
    return relays, start, r3_inst_hv, i_fault


# ---------------- Plots ----------------
def plot_loads(df):
    fig, ax = plt.subplots(figsize=(8, 4))
    y = np.arange(len(df))
    ax.barh(y - 0.2, df.kw, 0.4, label="kW (active)")
    ax.barh(y + 0.2, df.kvar, 0.4, label="kVAR (reactive)")
    ax.set_yticks(y, [f"{t}  {d}" for t, d in zip(df.tag, df.description)], fontsize=8)
    ax.invert_yaxis()
    ax.set_title("Running load by consumer")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "load_breakdown.png", dpi=130)
    plt.close(fig)


def plot_triangle(p, q1, pfc):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for q, lab, c in [(q1, f"Before: PF {pfc['pf_before']:.3f}", "tab:red"),
                      (pfc["q_after"], f"After: PF {pfc['pf_after']:.3f}", "tab:green")]:
        ax.plot([0, p, p, 0], [0, 0, q, 0], color=c, lw=2, label=lab)
        ax.text(p * 0.45, q * 0.55 + 8, f"{math.hypot(p, q):.0f} kVA", color=c)
    ax.annotate("", xy=(p + 15, pfc["q_after"]), xytext=(p + 15, q1),
                arrowprops={"arrowstyle": "->", "color": "tab:blue"})
    ax.text(p + 22, (q1 + pfc["q_after"]) / 2, f"-{pfc['qc_selected']} kVAR\ncapacitor", color="tab:blue", fontsize=9)
    ax.set_xlabel("kW")
    ax.set_ylabel("kVAR")
    ax.set_xlim(0, p * 1.35)
    ax.set_title("Power triangle before and after correction")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(RESULTS / "power_triangle.png", dpi=130)
    plt.close(fig)


def plot_tcc(relays, start, i_fault):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    i = np.geomspace(200, 40000, 600)
    for r, c in zip(relays, ["tab:green", "tab:blue", "tab:red"]):
        t = idmt_time(r["curve"], r["tms"], i, r["pickup_lv"])
        ax.plot(i, t, color=c, lw=2, label=f"{r['name']} {r['curve']}  pickup {r['pickup_lv']:.0f} A (LV), TMS {r['tms']:.2f}")
        ax.plot(i_fault, r["t_at_fault"], "o", color=c)
        ax.text(i_fault * 1.07, r["t_at_fault"], f"{r['t_at_fault']:.2f} s", color=c, va="center", fontsize=9)
    ax.axvline(i_fault, color="grey", ls="--")
    ax.text(i_fault * 0.98, 60, f"max 3-ph fault\n{i_fault/1000:.1f} kA", ha="right", fontsize=9, color="grey")
    ax.plot(start["i_start"], CFG["motor_start_time_s"], "k^")
    ax.text(start["i_start"] * 1.08, CFG["motor_start_time_s"], "fan DOL start\n(must NOT trip)", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(200, 40000)
    ax.set_ylim(0.03, 100)
    ax.set_xlabel("Current referred to 415 V (A)")
    ax.set_ylabel("Operating time (s)")
    ax.set_title("Time-current curves: R1 -> R2 -> R3 graded by >= 0.3 s")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(RESULTS / "tcc_relay_coordination.png", dpi=130)
    plt.close(fig)


# ---------------- Main ----------------
def main():
    loads = load_summary(pd.read_csv(ROOT / "data" / "mcc_loads.csv"))
    p, q = loads.kw.sum(), loads.kvar.sum()
    s1 = math.hypot(p, q)
    pfc = pf_correction(p, q)

    need1, tx1 = pick_transformer(s1)
    need2, tx2 = pick_transformer(pfc["s_after"])
    i1, i2 = amps(s1, CFG["v_lv"]), amps(pfc["s_after"], CFG["v_lv"])
    cu_full = CFG["tx_copper_loss_pct"] / 100 * tx2
    cu1, cu2 = cu_full * (s1 / tx2) ** 2, cu_full * (pfc["s_after"] / tx2) ** 2
    kwh = p * CFG["hours_per_year"]
    save_kvah = (kwh / pfc["pf_before"] - kwh / pfc["pf_after"]) * CFG["tariff_kvah"]
    save_demand = (s1 - pfc["s_after"]) * CFG["tariff_kva_month"] * 12

    motors = loads[loads.type == "motor"]
    fault = fault_level(tx2, (motors.flc_a * motors.qty_running).sum())
    breaking = next(k for k in CFG["std_breaking_ka"] if k >= 1.1 * fault["i_fault_total"] / 1000)
    relays, start, r3_inst, i_fault = coordinate(loads, tx2, fault)

    plot_loads(loads)
    plot_triangle(p, q, pfc)
    plot_tcc(relays, start, i_fault)

    tbl = loads[["tag", "description", "qty_running", "rated_kw", "kw", "kvar", "kva"]].round(1)
    lines = [
        "# Electrical design study: 415 V utilities MCC (auto-generated by src/study.py)\n",
        "## 1. Load summary\n", tbl.to_markdown(index=False), "",
        f"**Total running load: {p:.1f} kW, {q:.1f} kVAR, {s1:.1f} kVA -> power factor {pfc['pf_before']:.3f}**\n",
        "## 2. Power-factor correction\n",
        f"- Qc = P (tan phi1 - tan phi2) = {p:.1f} x ({q/p:.4f} - {math.tan(math.acos(CFG['target_pf'])):.4f}) "
        f"= **{pfc['qc_needed']:.1f} kVAR** -> select **{pfc['qc_selected']} kVAR APFC panel ({pfc['steps']} x {CFG['cap_step_kvar']} kVAR steps)**",
        f"- After correction: {pfc['s_after']:.1f} kVA, PF **{pfc['pf_after']:.3f}**",
        f"- Line current at 415 V: {i1:.0f} A -> {i2:.0f} A (**{(1 - i2/i1)*100:.1f} % less current, "
        f"{(1 - (i2/i1)**2)*100:.1f} % less I2R loss** in cables)",
        f"- Transformer copper loss: {cu1:.2f} kW -> {cu2:.2f} kW = {(cu1-cu2)*CFG['hours_per_year']:.0f} kWh/year saved",
        f"- Illustrative bill impact: demand charge Rs {save_demand/1e5:.2f} lakh/yr; "
        f"if billed in kVAh, Rs {save_kvah/1e5:.1f} lakh/yr (tariffs in CFG are assumptions)\n",
        "## 3. Transformer sizing (with 20 % future margin)\n",
        "| Case | Load kVA | Required kVA | Standard rating |", "|---|---|---|---|",
        f"| Without correction | {s1:.0f} | {need1:.0f} | **{tx1} kVA** |",
        f"| With correction | {pfc['s_after']:.0f} | {need2:.0f} | **{tx2} kVA** |", "",
        f"Selected: **{tx2} kVA, 11/0.433 kV, Dyn11, %Z = {tx_impedance_pct(tx2)}**. "
        f"Loading {pfc['s_after']/tx2*100:.0f} % today.\n",
        "## 4. Fault level at the 415 V bus (per-unit, base = transformer rating)\n",
        f"- Z_source = {tx2/1000:.2f} MVA / {CFG['source_fault_mva']:.0f} MVA = {fault['z_src']:.4f} pu; "
        f"Z_transformer = {fault['z_tx']:.4f} pu; total {fault['z_tot']:.4f} pu",
        f"- I_base = {fault['i_base']:.0f} A -> grid in-feed fault current **{fault['i_fault_grid']/1000:.1f} kA** "
        f"({fault['fault_mva']:.1f} MVA)",
        f"- Motor contribution ~{CFG['motor_contribution_x']:.0f} x FLC = {fault['i_motor']/1000:.1f} kA -> "
        f"total **{fault['i_fault_total']/1000:.1f} kA** -> specify switchgear for **{breaking} kA**\n",
        "## 5. Relay settings and grading (IEC 60255 curves)\n",
        "| Relay | Location | Curve | CT | Plug setting | Pickup (A, primary) | TMS | Time at max fault |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in relays:
        lines.append(f"| {r['name']} | {r['where']} | {r['curve']} | {r['ct']} | {r['ps']*100:.0f} % | "
                     f"{r['pickup_a']:.1f} | {r['tms']:.2f} | {r['t_at_fault']:.3f} s |")
    lines += [
        "",
        f"- Grading margins at {i_fault/1000:.1f} kA: R2-R1 = {relays[1]['t_at_fault']-relays[0]['t_at_fault']:.2f} s, "
        f"R3-R2 = {relays[2]['t_at_fault']-relays[1]['t_at_fault']:.2f} s (required >= {CFG['grading_margin_s']} s, checked at every current)",
        f"- Motor-start check: R1 sees {start['i_start']:.0f} A during a fan DOL start and would need "
        f"{start['t_r1_at_start']:.1f} s to trip > {CFG['motor_start_time_s']:.0f} s run-up -> **no nuisance trip** "
        "(why R1 uses the Extremely Inverse curve)",
        f"- R3 instantaneous high-set: {r3_inst:.0f} A on 11 kV (1.3 x max LV through-fault), clears 11 kV-side faults fast",
    ]
    (RESULTS / "study_report.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
