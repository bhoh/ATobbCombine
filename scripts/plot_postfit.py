#!/usr/bin/env python3
"""
Pre-fit and Post-fit Distribution Plotter with MC Uncertainty Bands for ATobb Combine Analysis.
Reads shapes from fitDiagnostics.root (and nominal shapes from shapes_MA18.root if available)
and generates publication-quality Data-to-MC comparison plots with stacked backgrounds,
error bands, and Data/MC ratio panels.

Usage:
  python3 scripts/plot_postfit.py [--fit-root fitDiagnostics_CR.root] [--shapes-root shapes_MA18.root] [--output-dir DIR]
"""

import argparse
import os
import sys
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    BASELINE_DIR,
    LUMI_FB,
)

# Colors and clean labels for background processes in stacked plots
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

PROCESS_LABELS = {
    "TT_bb": r"$\mathrm{t\bar{t}+b\bar{b}}$",
    "TT_bj": r"$\mathrm{t\bar{t}+bj}$",
    "TT_cc": r"$\mathrm{t\bar{t}+c\bar{c}}$",
    "TT_cj": r"$\mathrm{t\bar{t}+cj}$",
    "TT_jj": r"$\mathrm{t\bar{t}+jj}$",
    "ST": "Single Top",
    "DY": r"Drell-Yan",
    "signal": r"$\mathrm{H^+ \to cb}$",
}


def compute_mc_uncertainties(cat_dir, cat_name, shapes_file=None):
    """
    Compute bin-by-bin MC total and statistical uncertainties.
    Combines ROOT histogram errors and nominal shapes file errors when available.
    """
    h_total = cat_dir.Get("total_background") or cat_dir.Get("total")
    if not h_total:
        return None, None, None

    nbins = h_total.GetNbinsX()
    total_vals = np.array([h_total.GetBinContent(i) for i in range(1, nbins + 1)])
    combine_errs = np.array([h_total.GetBinError(i) for i in range(1, nbins + 1)])

    # Calculate process-by-process statistical uncertainty from shapes file if available
    procs = ["DY", "ST", "TT_jj", "TT_cj", "TT_cc", "TT_bj", "TT_bb"]
    stat_err2 = np.zeros(nbins)

    if shapes_file and not shapes_file.IsZombie():
        for p in procs:
            h_fit = cat_dir.Get(p)
            h_nom = shapes_file.Get(f"histo_{cat_name}_{p}")
            if h_fit and h_nom:
                for b in range(1, nbins + 1):
                    c_nom = h_nom.GetBinContent(b)
                    e_nom = h_nom.GetBinError(b)
                    c_fit = h_fit.GetBinContent(b)
                    scale = (c_fit / c_nom) if c_nom > 0 else 1.0
                    stat_err2[b - 1] += (scale * e_nom) ** 2
    else:
        # Fallback to individual TH1 errors in cat_dir
        for p in procs:
            h_fit = cat_dir.Get(p)
            if h_fit:
                for b in range(1, nbins + 1):
                    stat_err2[b - 1] += (h_fit.GetBinError(b)) ** 2

    stat_errs = np.sqrt(stat_err2)

    # Use maximum of Combine Hesse error or MC stat error
    total_errs = np.maximum(combine_errs, stat_errs)
    return total_vals, total_errs, stat_errs


def plot_category_distribution(cat_dir, cat_name, output_path, fit_type="fit_b", is_cr=True, shapes_file=None):
    """Plot pre-fit or post-fit data vs stacked background with MC error band and ratio panel."""
    if not cat_dir:
        return

    # Extract background processes in stacked order
    bkg_procs = ["DY", "ST", "TT_jj", "TT_cj", "TT_cc", "TT_bj", "TT_bb"]
    h_bkgs = []
    bkg_labels = []
    bkg_colors = []

    for p in bkg_procs:
        h = cat_dir.Get(p)
        if h and h.Integral() > 0:
            h_bkgs.append(h)
            bkg_labels.append(PROCESS_LABELS.get(p, p))
            bkg_colors.append(PROCESS_COLORS.get(p, "#999999"))

    h_total = cat_dir.Get("total_background") or cat_dir.Get("total")
    h_data = cat_dir.Get("data")

    if not h_total or not h_bkgs:
        return

    nbins = h_total.GetNbinsX()
    bin_edges = np.array([h_total.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    bin_centers = np.array([h_total.GetBinCenter(i) for i in range(1, nbins + 1)])
    bin_widths = np.array([h_total.GetBinWidth(i) for i in range(1, nbins + 1)])

    # Compute stacked MC data
    stack_data = []
    for h in h_bkgs:
        contents = [h.GetBinContent(i) for i in range(1, nbins + 1)]
        stack_data.append(contents)

    # Calculate uncertainties
    total_vals, total_errs, stat_errs = compute_mc_uncertainties(cat_dir, cat_name, shapes_file)

    # Observed Data Points (only displayed in CR or when explicitly unblinded)
    has_data = (h_data is not None and h_data.GetN() > 0 and is_cr)
    data_x, data_y, data_ey_lo, data_ey_hi = [], [], [], []
    if has_data:
        import ctypes
        for i in range(h_data.GetN()):
            x = ctypes.c_double(0.0)
            y = ctypes.c_double(0.0)
            h_data.GetPoint(i, x, y)
            ey_lo = h_data.GetErrorYlow(i)
            ey_hi = h_data.GetErrorYhigh(i)
            data_x.append(x.value)
            data_y.append(y.value)
            data_ey_lo.append(ey_lo)
            data_ey_hi.append(ey_hi)
        data_x = np.array(data_x)
        data_y = np.array(data_y)
        data_ey = [np.array(data_ey_lo), np.array(data_ey_hi)]

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

    # Total MC uncertainty band (hatched error band)
    unc_boxes = []
    for b_idx in range(nbins):
        x_low = bin_edges[b_idx]
        x_high = bin_edges[b_idx + 1]
        y_low = max(0.0, total_vals[b_idx] - total_errs[b_idx])
        y_high = total_vals[b_idx] + total_errs[b_idx]
        ax_main.fill_between(
            [x_low, x_high],
            [y_low, y_low],
            [y_high, y_high],
            facecolor="none",
            edgecolor="#444444",
            hatch="//////",
            linewidth=0.0,
            zorder=3,
        )

    # Proxy patch for MC uncertainty band in legend
    unc_patch = mpatches.Patch(
        facecolor="none",
        edgecolor="#444444",
        hatch="//////",
        label="Bkg. Stat. Unc." if "prefit" in fit_type else "Bkg. Uncertainty",
    )

    # Observed Data Points
    if has_data and len(data_x) > 0:
        ax_main.errorbar(
            data_x,
            data_y,
            yerr=data_ey,
            fmt="o",
            color="black",
            markersize=5,
            linewidth=1.5,
            capsize=0,
            label="Data (CR)",
            zorder=5,
        )

    # Header and Annotations
    if HAS_MPLHEP:
        hep.cms.label(
            llabel="Preliminary",
            data=has_data,
            lumi=int(LUMI_FB),
            com=COM_TEV,
            loc=0,
            ax=ax_main,
        )
    else:
        ax_main.set_title(
            f"CMS Preliminary   {int(LUMI_FB)} fb$^{{-1}}$ ({COM_TEV} TeV, 2024)",
            loc="right",
            fontsize=16,
        )

    mode_label = "Post-fit" if "fit" in fit_type else "Pre-fit"
    cat_display = cat_name.replace("_", " ")
    ax_main.text(
        0.04,
        0.90,
        f"{cat_display} ({mode_label})",
        transform=ax_main.transAxes,
        fontsize=18,
        fontweight="bold",
        verticalalignment="top",
    )
    ax_main.set_ylabel("Events / bin", fontsize=20)

    # Y-axis scaling
    max_y = max(total_vals.max(), (data_y.max() if has_data and len(data_y) > 0 else 0))
    ax_main.set_ylim(0.0, max_y * 1.42)

    # Build handles and labels with uncertainty patch included
    handles, labels = ax_main.get_legend_handles_labels()
    handles.append(unc_patch)
    labels.append(unc_patch.get_label())

    # Reverse legend order to match stack order top-to-bottom
    ax_main.legend(
        handles[::-1],
        labels[::-1],
        loc="upper right",
        frameon=True,
        facecolor="white",
        framealpha=0.9,
        fontsize=12,
        ncol=2,
    )

    # 2. Ratio Panel (Data / MC)
    ax_ratio.axhline(1.0, color="#555555", linestyle="--", linewidth=1.2, zorder=1)

    # Relative MC uncertainty band centered at 1.0
    ratio_mc_err = np.where(total_vals > 0, total_errs / total_vals, 0.0)
    for b_idx in range(nbins):
        x_low = bin_edges[b_idx]
        x_high = bin_edges[b_idx + 1]
        err_val = ratio_mc_err[b_idx]
        ax_ratio.fill_between(
            [x_low, x_high],
            [1.0 - err_val, 1.0 - err_val],
            [1.0 + err_val, 1.0 + err_val],
            facecolor="none",
            edgecolor="#444444",
            hatch="//////",
            linewidth=0.0,
            zorder=2,
        )

    if has_data and len(data_x) > 0:
        ratio_y = np.where(total_vals > 0, data_y / total_vals, np.nan)
        ratio_ey_lo = np.where(total_vals > 0, data_ey[0] / total_vals, np.nan)
        ratio_ey_hi = np.where(total_vals > 0, data_ey[1] / total_vals, np.nan)
        ax_ratio.errorbar(
            data_x,
            ratio_y,
            yerr=[ratio_ey_lo, ratio_ey_hi],
            fmt="o",
            color="black",
            markersize=5,
            linewidth=1.5,
            capsize=0,
            zorder=5,
        )

    # Determine score cut domain label
    if is_cr:
        x_title = r"$\mathrm{bdt\_score\_avg} \leq 0.3$"
    else:
        x_title = r"$\mathrm{bdt\_score\_avg} > 0.3$"

    ax_ratio.set_xlabel(x_title, fontsize=20)
    ax_ratio.set_ylabel("Data / MC", fontsize=16)
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.set_yticks([0.6, 0.8, 1.0, 1.2, 1.4])
    ax_ratio.set_xlim(bin_edges[0], bin_edges[-1])

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    pdf_path = os.path.splitext(output_path)[0] + ".pdf"
    plt.savefig(pdf_path)
    plt.close()

    print(f"  [+] Saved {mode_label} plot for {cat_name} -> {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot Pre-fit and Post-fit distributions from fitDiagnostics.root with MC error bands")
    parser.add_argument(
        "--fit-root",
        "-f",
        default=os.path.join(BASELINE_DIR, "cr_fit/fitDiagnostics_CR.root"),
        help="Path to fitDiagnostics.root",
    )
    parser.add_argument(
        "--shapes-root",
        "-s",
        default=os.path.join(BASELINE_DIR, "shapes/shapes_MA18.root"),
        help="Path to shapes.root containing nominal MC statistical errors",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default=os.path.join(BASELINE_DIR, "cr_fit/plots"),
        help="Output directory for plots",
    )
    args = parser.parse_args()

    if not os.path.exists(args.fit_root):
        print(f"Error: {args.fit_root} not found. Run run_cr_fit.py first.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    f_fit = ROOT.TFile.Open(args.fit_root)
    if not f_fit or f_fit.IsZombie():
        print(f"Error opening {args.fit_root}")
        sys.exit(1)

    f_shape = None
    if os.path.exists(args.shapes_root):
        f_shape = ROOT.TFile.Open(args.shapes_root)

    # Check available shape directories (e.g. shapes_prefit, shapes_fit_b, shapes_fit_s)
    shape_dirs = ["shapes_prefit", "shapes_fit_b"]
    for sdir_name in shape_dirs:
        sdir = f_fit.Get(sdir_name)
        if not sdir:
            continue
        for key in sdir.GetListOfKeys():
            cat_name = key.GetName()
            cat_obj = sdir.Get(cat_name)
            out_png = os.path.join(args.output_dir, f"{cat_name}_{sdir_name}.png")
            plot_category_distribution(
                cat_obj,
                cat_name,
                out_png,
                fit_type=sdir_name,
                is_cr=("CR" in cat_name),
                shapes_file=f_shape,
            )

    if f_shape and not f_shape.IsZombie():
        f_shape.Close()
    f_fit.Close()
    print(f"\n[+] All distribution plots with MC error bands saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
