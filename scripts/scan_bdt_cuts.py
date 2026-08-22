#!/usr/bin/env python3
"""
BDT Score Cut Optimization Study for ATobb Combine Analysis.
Scans multiple CR/SR boundary thresholds on bdt_score_avg (e.g. 0.15, 0.20, 0.25, 0.30, 0.35),
runs simultaneous CR+SR fits, compares expected median limit sensitivity across mass points,
and produces a summary comparison plot.

Usage:
  python3 scripts/scan_bdt_cuts.py [--cuts 0.15,0.20,0.25,0.30,0.35] [--rebin 2] [--outdir DIR]
"""

import argparse
import csv
import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    HAS_MPLHEP = True
except ImportError:
    HAS_MPLHEP = False

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import (
    COM_TEV,
    DEFAULT_REBIN,
    LUMI_FB,
    MASS_POINTS,
    STUDIES_SCORE_CUTS_DIR,
)
from modules.shape_extractor import extract_shapes_for_categories
from modules.datacard_generator import write_datacard
from modules.workspace_builder import build_workspace
from scripts.run_limits import run_asymptotic_limit


def main():
    parser = argparse.ArgumentParser(description="BDT Score Cut Optimization Scan")
    parser.add_argument("--cuts", default="0.15,0.20,0.25,0.30,0.35,0.40", help="Comma-separated list of cut thresholds to scan")
    parser.add_argument("--rebin", type=int, default=DEFAULT_REBIN, help="Rebin factor (default: 2)")
    parser.add_argument("--sr-var", choices=["avg", "mass_bdt"], default="avg", help="SR discriminator (default: avg)")
    parser.add_argument("--outdir", default=STUDIES_SCORE_CUTS_DIR, help=f"Study output directory (default: {STUDIES_SCORE_CUTS_DIR})")
    args = parser.parse_args()

    cut_values = [float(c.strip()) for c in args.cuts.split(",") if c.strip()]
    os.makedirs(args.outdir, exist_ok=True)

    print("\n" + "=" * 90)
    print(f" BDT Score Cut Scan Study: Thresholds = {cut_values}")
    print(f" Output Directory: {args.outdir}")
    print("=" * 90)

    scan_results = {}  # cut_val -> { ma -> limit_dict }

    for cut in cut_values:
        tag = f"cut_{int(round(cut*100)):02d}"
        tag_dir = os.path.join(args.outdir, tag)
        shapes_dir = os.path.join(tag_dir, "shapes")
        cards_dir = os.path.join(tag_dir, "datacards")
        ws_dir = os.path.join(tag_dir, "workspaces")
        for d in [shapes_dir, cards_dir, ws_dir]:
            os.makedirs(d, exist_ok=True)

        print(f"\n---> Evaluating CR cut <= {cut:.2f} (SR > {cut:.2f}) [workspace: {ws_dir}]")
        scan_results[cut] = {}

        for ma in MASS_POINTS:
            sr_var = f"bdt_score_{ma}" if args.sr_var == "mass_bdt" else "bdt_score_avg"
            categories = [
                {"cat_name": "Muon_CR", "channel": "Muon", "var_name": "bdt_score_avg", "x_min": 0.0, "x_max": cut},
                {"cat_name": "Electron_CR", "channel": "Electron", "var_name": "bdt_score_avg", "x_min": 0.0, "x_max": cut},
                {"cat_name": "Muon_SR", "channel": "Muon", "var_name": sr_var, "x_min": cut, "x_max": 1.0},
                {"cat_name": "Electron_SR", "channel": "Electron", "var_name": sr_var, "x_min": cut, "x_max": 1.0},
            ]

            shape_file = os.path.join(shapes_dir, f"shapes_{ma}.root")
            card_file = os.path.join(cards_dir, f"datacard_{ma}.txt")
            ws_file = os.path.join(ws_dir, f"workspace_{ma}.root")

            yields = extract_shapes_for_categories(categories, shape_file, ma_point=ma, rebin_factor=args.rebin)
            write_datacard(
                categories=["Muon_CR", "Electron_CR", "Muon_SR", "Electron_SR"],
                shapes_root_path=shape_file,
                yields=yields,
                output_card_path=card_file,
                ma_point=ma,
                title_comment=f"CR<={cut:.2f}, SR>{cut:.2f}",
            )
            build_workspace(card_file, ws_file, channel_masks=True)
            lim = run_asymptotic_limit(ws_file, ws_dir, ma, unblind=False)
            if lim:
                scan_results[cut][ma] = lim
                print(f"      [{ma}] Expected Median BR < {lim['exp_median']:.6f}  ([-1s, +1s]: [{lim['exp_minus1s']:.6f}, {lim['exp_plus1s']:.6f}])")

    # Export comparison CSV
    csv_file = os.path.join(args.outdir, "bdt_cut_scan_comparison.csv")
    fieldnames = ["mass_point", "mass_GeV"] + [f"median_cut_{int(round(c*100)):02d}" for c in cut_values]
    with open(csv_file, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fieldnames)
        for ma in MASS_POINTS:
            m_val = float(ma.replace("MA", ""))
            row = [ma, f"{m_val:.1f}"]
            for c in cut_values:
                val = scan_results.get(c, {}).get(ma, {}).get("exp_median", "")
                row.append(f"{val:.6f}" if val != "" else "")
            writer.writerow(row)
    print(f"\n[+] Saved scan comparison CSV -> {csv_file}")

    # Plot comparison curves
    fig, ax = plt.subplots(figsize=(10, 8))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"]
    markers = ["o", "s", "^", "v", "D", "P", "X"]

    masses_arr = np.array([float(ma.replace("MA", "")) for ma in MASS_POINTS])
    for idx, c in enumerate(cut_values):
        med_vals = [scan_results.get(c, {}).get(ma, {}).get("exp_median", np.nan) for ma in MASS_POINTS]
        ax.plot(
            masses_arr,
            med_vals,
            label=f"CR $\\leq {c:.2f}$ (SR $> {c:.2f}$)",
            color=colors[idx % len(colors)],
            marker=markers[idx % len(markers)],
            linewidth=2.2,
        )

    if HAS_MPLHEP:
        hep.cms.label(llabel="Preliminary", data=False, lumi=int(LUMI_FB), com=COM_TEV, loc=0, ax=ax)
    else:
        ax.set_title(f"CMS Preliminary   {int(LUMI_FB)} fb$^{{-1}}$ ({COM_TEV} TeV, 2024)", loc="right", fontsize=16)

    ax.text(0.04, 0.90, r"$t \to H^{+}b,\, H^{+} \to W^{+}A\,(\to b\bar{b})$", transform=ax.transAxes, fontsize=20, verticalalignment="top")
    ax.text(0.04, 0.83, r"$m_{H^{+}} = 130\,\mathrm{GeV}$", transform=ax.transAxes, fontsize=18, verticalalignment="top")
    ax.text(0.04, 0.77, r"BDT Cut Optimization Scan", transform=ax.transAxes, fontsize=16, verticalalignment="top", color="darkblue")

    ax.set_xlabel(r"$m_{A}\,\,[\mathrm{GeV}]$", fontsize=22)
    ax.set_ylabel(r"$95\%\,\mathrm{CL}\,\mathrm{upper\,limit\,on}\,\mathcal{B}(t \to H^{+}b)$", fontsize=20)
    ax.set_xlim(min(masses_arr) - 5, max(masses_arr) + 5)

    all_meds = [lim["exp_median"] for c in cut_values for lim in scan_results.get(c, {}).values()]
    max_val = max(all_meds) if all_meds else 5.0
    ax.set_ylim(0.0, max_val * 1.35)

    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9, fontsize=14)

    plt.tight_layout()
    out_png = os.path.join(args.outdir, "bdt_cut_scan_comparison.png")
    out_pdf = os.path.join(args.outdir, "bdt_cut_scan_comparison.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()

    print(f"[+] Saved scan comparison plot -> {out_png}")


if __name__ == "__main__":
    main()
