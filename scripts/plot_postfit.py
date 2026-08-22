#!/usr/bin/env python3
"""
Pre-fit and Post-fit Distribution Plotter for ATobb Combine Analysis.
Reads shapes from fitDiagnostics.root (produced by combine -M FitDiagnostics --saveShapes)
and generates Data-to-MC comparison plots with ratio panels.

Usage:
  python3 scripts/plot_postfit.py [--fit-root fitDiagnostics_CR.root] [--output-dir DIR]
"""

import argparse
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

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import (
    BACKGROUND_PROCESSES,
    COM_TEV,
    DEFAULT_OUTPUT_DIR,
    LUMI_FB,
)

# Colors for background processes in stacked plots
PROCESS_COLORS = {
    "TT_bb": "#b30000",      # Dark red
    "TT_bj": "#e6550d",      # Orange red
    "TT_cc": "#31a354",      # Dark green
    "TT_cj": "#74c476",      # Light green
    "TT_jj": "#fec44f",      # Yellow
    "ST": "#756bb1",         # Purple
    "DY": "#2b8cbe",         # Blue
    "signal": "#000000",     # Black
}


def plot_category_distribution(cat_dir, cat_name, output_path, fit_type="fit_b", is_cr=True):
    """Plot pre-fit or post-fit data vs stacked background with ratio panel."""
    if not cat_dir:
        return

    # Extract processes
    bkg_procs = ["DY", "ST", "TT_jj", "TT_cj", "TT_cc", "TT_bj", "TT_bb"]
    h_bkgs = []
    bkg_labels = []
    bkg_colors = []

    for p in bkg_procs:
        h = cat_dir.Get(p)
        if h and h.Integral() > 0:
            h_bkgs.append(h)
            bkg_labels.append(p)
            bkg_colors.append(PROCESS_COLORS.get(p, "#999999"))

    h_total = cat_dir.Get("total_background") or cat_dir.Get("total")
    h_data = cat_dir.Get("data")

    if not h_total or not h_bkgs:
        return

    nbins = h_total.GetNbinsX()
    bin_edges = [h_total.GetBinLowEdge(i) for i in range(1, nbins + 2)]
    bin_centers = [h_total.GetBinCenter(i) for i in range(1, nbins + 1)]
    bin_widths = [h_total.GetBinWidth(i) for i in range(1, nbins + 1)]

    # Prepare stacked array
    stack_data = []
    for h in h_bkgs:
        contents = [h.GetBinContent(i) for i in range(1, nbins + 1)]
        stack_data.append(contents)

    total_vals = np.array([h_total.GetBinContent(i) for i in range(1, nbins + 1)])
    total_errs = np.array([h_total.GetBinError(i) for i in range(1, nbins + 1)])

    # Prepare data points
    has_data = (h_data is not None and h_data.GetN() > 0)
    data_x, data_y, data_ey = [], [], []
    if has_data:
        import ctypes
        for i in range(h_data.GetN()):
            x = ctypes.c_double(0.0)
            y = ctypes.c_double(0.0)
            h_data.GetPoint(i, x, y)
            ey = h_data.GetErrorY(i)
            data_x.append(x.value)
            data_y.append(y.value)
            data_ey.append(ey)
        data_x = np.array(data_x)
        data_y = np.array(data_y)
        data_ey = np.array(data_ey)

    # Setup 2-panel figure
    fig, (ax_main, ax_ratio) = plt.subplots(
        2, 1, figsize=(10, 9),
        gridspec_kw={"height_ratios": [3.2, 1], "hspace": 0.08},
        sharex=True,
    )

    # 1. Main Panel (Stacked MC)
    ax_main.hist(
        [bin_centers] * len(stack_data),
        bins=bin_edges,
        weights=stack_data,
        stacked=True,
        color=bkg_colors,
        label=bkg_labels,
        edgecolor="none",
        zorder=1,
    )

    # Total MC uncertainty band
    for b_idx in range(len(bin_centers)):
        x_low = bin_edges[b_idx]
        x_high = bin_edges[b_idx + 1]
        y_low = max(0, total_vals[b_idx] - total_errs[b_idx])
        y_high = total_vals[b_idx] + total_errs[b_idx]
        ax_main.fill_between([x_low, x_high], [y_low, y_low], [y_high, y_high], color="gray", alpha=0.35, hatch="//", zorder=2)

    # Observed Data Points
    if has_data and len(data_x) > 0:
        ax_main.errorbar(
            data_x, data_y, yerr=data_ey,
            fmt="o", color="black", markersize=5, linewidth=1.5,
            label="Data" if not is_cr else "Data (CR)",
            zorder=5,
        )

    # Header and Annotations
    if HAS_MPLHEP:
        hep.cms.label(llabel="Preliminary", data=has_data, lumi=int(LUMI_FB), com=COM_TEV, loc=0, ax=ax_main)
    else:
        ax_main.set_title(f"CMS Preliminary   {int(LUMI_FB)} fb$^{{-1}}$ ({COM_TEV} TeV, 2024)", loc="right", fontsize=16)

    mode_label = "Post-fit" if "fit" in fit_type else "Pre-fit"
    ax_main.text(0.04, 0.90, f"{cat_name} ({mode_label})", transform=ax_main.transAxes, fontsize=18, fontweight="bold", verticalalignment="top")
    ax_main.set_ylabel("Events / bin", fontsize=20)
    ax_main.set_ylim(0.0, max(total_vals.max(), (data_y.max() if has_data and len(data_y) > 0 else 0)) * 1.35)
    ax_main.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9, fontsize=13, ncol=2)

    # 2. Ratio Panel (Data / MC)
    ax_ratio.axhline(1.0, color="gray", linestyle="--", linewidth=1.2)

    # Relative MC uncertainty band in ratio panel
    ratio_mc_err = np.where(total_vals > 0, total_errs / total_vals, 0.0)
    for b_idx in range(len(bin_centers)):
        x_low = bin_edges[b_idx]
        x_high = bin_edges[b_idx + 1]
        err_val = ratio_mc_err[b_idx]
        ax_ratio.fill_between([x_low, x_high], [1 - err_val, 1 - err_val], [1 + err_val, 1 + err_val], color="gray", alpha=0.35, hatch="//")

    if has_data and len(data_x) > 0:
        ratio_y = np.where(total_vals > 0, data_y / total_vals, np.nan)
        ratio_ey = np.where(total_vals > 0, data_ey / total_vals, np.nan)
        ax_ratio.errorbar(data_x, ratio_y, yerr=ratio_ey, fmt="o", color="black", markersize=5, linewidth=1.5)

    ax_ratio.set_xlabel(r"$\mathrm{bdt\_score\_avg}$", fontsize=20)
    ax_ratio.set_ylabel("Data / MC", fontsize=16)
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    pdf_path = os.path.splitext(output_path)[0] + ".pdf"
    plt.savefig(pdf_path)
    plt.close()

    print(f"  [+] Saved {mode_label} plot for {cat_name} -> {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot Pre-fit and Post-fit distributions from fitDiagnostics.root")
    parser.add_argument(
        "--fit-root",
        "-f",
        default=os.path.join(DEFAULT_OUTPUT_DIR, "outputs/cr_only/fitDiagnostics_CR.root"),
        help="Path to fitDiagnostics.root",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default=os.path.join(DEFAULT_OUTPUT_DIR, "outputs/cr_only/plots"),
        help="Output directory for plots",
    )
    args = parser.parse_args()

    if not os.path.exists(args.fit_root):
        print(f"Error: {args.fit_root} not found. Run run_cr_fit.py first.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    f = ROOT.TFile.Open(args.fit_root)
    if not f or f.IsZombie():
        print(f"Error opening {args.fit_root}")
        sys.exit(1)

    # Check available shape directories (e.g. shapes_prefit, shapes_fit_b, shapes_fit_s)
    shape_dirs = ["shapes_prefit", "shapes_fit_b"]
    for sdir_name in shape_dirs:
        sdir = f.Get(sdir_name)
        if not sdir:
            continue
        for key in sdir.GetListOfKeys():
            cat_name = key.GetName()
            cat_obj = sdir.Get(cat_name)
            out_png = os.path.join(args.output_dir, f"{cat_name}_{sdir_name}.png")
            plot_category_distribution(cat_obj, cat_name, out_png, fit_type=sdir_name, is_cr=("CR" in cat_name))

    f.Close()
    print(f"\n[+] All distribution plots saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
