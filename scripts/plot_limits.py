#!/usr/bin/env python3
"""
Plot CMS 95% CL CLs Upper Limits (Brazilian Plot) using mplhep and ROOT.
Reads directly from higgsCombine*.AsymptoticLimits*.root files or limits CSV files.

Usage:
  python3 scripts/plot_limits.py [--dir DIR] [--csv CSV] [--output OUTPUT.png] [--lumi LUMI] [--unblind] [--logy]
"""

import argparse
import csv
import glob
import os
import re
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


def parse_mass_from_filename(filename):
    """Extract mass tag (e.g. MA18) and numerical value from filename."""
    base = os.path.basename(filename)
    match = re.search(r"(MA\d+)", base, re.IGNORECASE)
    if match:
        ma_str = match.group(1).upper()
        ma_val = float(ma_str.replace("MA", ""))
        return ma_str, ma_val
    match_mh = re.search(r"mH(\d+)", base)
    if match_mh:
        val = float(match_mh.group(1))
        return f"M{int(val)}", val
    return base, 0.0


def extract_limits_from_root_files(input_dir, pattern="higgsCombine*.AsymptoticLimits*.root"):
    """Read all limit trees from ROOT files in a directory."""
    search_pattern = os.path.join(input_dir, pattern)
    file_list = glob.glob(search_pattern)

    valid_files = [f for f in file_list if not os.path.basename(f).startswith("higgsCombineTest.")]
    if not valid_files:
        valid_files = file_list

    results = []
    for fpath in valid_files:
        ma_str, ma_val = parse_mass_from_filename(fpath)
        f = ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            continue
        tree = f.Get("limit")
        if not tree or tree.GetEntries() == 0:
            f.Close()
            continue

        lim = {"obs": None, "exp_2p5": None, "exp_16": None, "exp_50": None, "exp_84": None, "exp_97p5": None}
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            q = round(tree.quantileExpected, 4)
            val = float(tree.limit)
            if q < 0:
                lim["obs"] = val
            elif abs(q - 0.025) < 0.005:
                lim["exp_2p5"] = val
            elif abs(q - 0.160) < 0.005:
                lim["exp_16"] = val
            elif abs(q - 0.500) < 0.005:
                lim["exp_50"] = val
            elif abs(q - 0.840) < 0.005:
                lim["exp_84"] = val
            elif abs(q - 0.975) < 0.005:
                lim["exp_97p5"] = val
        f.Close()

        if lim["exp_50"] is not None:
            results.append({
                "mass_point": ma_str,
                "mass_GeV": ma_val,
                "exp_minus2s": lim["exp_2p5"],
                "exp_minus1s": lim["exp_16"],
                "exp_median": lim["exp_50"],
                "exp_plus1s": lim["exp_84"],
                "exp_plus2s": lim["exp_97p5"],
                "observed": lim["obs"],
            })

    results.sort(key=lambda x: (x["mass_GeV"], x["mass_point"]))
    return results


def extract_limits_from_csv(csv_path):
    """Read limits from a CSV file."""
    results = []
    with open(csv_path, mode="r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            obs_val = None
            if "observed" in row and row["observed"] != "" and row["observed"] != "BLINDED":
                try:
                    obs_val = float(row["observed"])
                except ValueError:
                    pass
            results.append({
                "mass_point": row["mass_point"],
                "mass_GeV": float(row["mass_GeV"]),
                "exp_minus2s": float(row["exp_minus2s"]) if row.get("exp_minus2s") else None,
                "exp_minus1s": float(row["exp_minus1s"]) if row.get("exp_minus1s") else None,
                "exp_median": float(row["exp_median"]) if row.get("exp_median") else None,
                "exp_plus1s": float(row["exp_plus1s"]) if row.get("exp_plus1s") else None,
                "exp_plus2s": float(row["exp_plus2s"]) if row.get("exp_plus2s") else None,
                "observed": obs_val,
            })
    results.sort(key=lambda x: x["mass_GeV"])
    return results


def plot_brazilian_limit(
    results,
    output_path="limit_cls.png",
    lumi=110.0,
    com=13.6,
    year="2024",
    status="Preliminary",
    show_observed=False,
    scale_factor=1.0,
    y_label=None,
    y_min=None,
    y_max=None,
    title_extra="",
    logy=False,
    is_percent=False,
):
    if not results:
        print("Error: No limit data to plot.")
        return

    eff_scale = 100.0 * scale_factor if is_percent else scale_factor
    if y_label is None:
        y_label = (
            r"$95\%\,\mathrm{CL}\,\mathrm{upper\,limit\,on}\,\mathcal{B}(t \to H^{+}b)\,\,[\%]$"
            if is_percent
            else r"$95\%\,\mathrm{CL}\,\mathrm{upper\,limit\,on}\,\mathcal{B}(t \to H^{+}b)$"
        )

    masses = np.array([r["mass_GeV"] for r in results])
    exp_m2s = np.array([r["exp_minus2s"] * eff_scale for r in results])
    exp_m1s = np.array([r["exp_minus1s"] * eff_scale for r in results])
    exp_med = np.array([r["exp_median"] * eff_scale for r in results])
    exp_p1s = np.array([r["exp_plus1s"] * eff_scale for r in results])
    exp_p2s = np.array([r["exp_plus2s"] * eff_scale for r in results])

    has_obs = show_observed and any(r["observed"] is not None for r in results)
    if has_obs:
        obs_vals = np.array([r["observed"] * eff_scale if r["observed"] is not None else np.nan for r in results])

    fig, ax = plt.subplots(figsize=(10, 8))

    # Standard CMS Brazilian flag colors
    green = "#607641"   # 68% expected (1 sigma)
    yellow = "#F5BB54"  # 95% expected (2 sigma)

    # Fill 95% (2 sigma) and 68% (1 sigma) bands
    ax.fill_between(masses, exp_m2s, exp_p2s, color=yellow, label=r"$95\%$ expected", zorder=1)
    ax.fill_between(masses, exp_m1s, exp_p1s, color=green, label=r"$68\%$ expected", zorder=2)

    # Plot Median expected limit
    ax.plot(masses, exp_med, color="black", linestyle="--", linewidth=2.5, label="Median expected", zorder=3)

    # Plot Observed limit if available
    if has_obs:
        ax.plot(masses, obs_vals, color="black", marker="o", markersize=6, linewidth=2.0, label="Observed", zorder=4)

    # CMS header
    lumi_val = int(lumi) if float(lumi).is_integer() else lumi
    if HAS_MPLHEP:
        hep.cms.label(
            llabel=status,
            data=has_obs,
            lumi=lumi_val,
            com=com,
            loc=0,
            ax=ax,
        )
    else:
        lumi_str = f"{int(lumi)}" if float(lumi).is_integer() else f"{lumi:.1f}"
        ax.set_title(f"CMS {status}   {lumi_str} fb$^{{-1}}$ ({com} TeV, {year})", loc="right", fontsize=16)

    # Physics process annotations
    ax.text(
        0.04, 0.90,
        r"$t \to H^{+}b,\, H^{+} \to W^{+}A\,(\to b\bar{b})$",
        transform=ax.transAxes,
        fontsize=20,
        verticalalignment="top",
    )
    ax.text(
        0.04, 0.83,
        r"$m_{H^{+}} = 130\,\mathrm{GeV}$",
        transform=ax.transAxes,
        fontsize=18,
        verticalalignment="top",
    )
    if title_extra:
        ax.text(
            0.04, 0.77,
            title_extra,
            transform=ax.transAxes,
            fontsize=16,
            verticalalignment="top",
            color="darkblue",
        )

    # Formatting axes
    ax.set_xlabel(r"$m_{A}\,\,[\mathrm{GeV}]$", fontsize=22)
    ax.set_ylabel(y_label, fontsize=20)
    ax.set_xlim(min(masses) - 5, max(masses) + 5)

    # Set y scale and limits
    if logy:
        ax.set_yscale("log")
        if y_min is not None and y_max is not None:
            ax.set_ylim(y_min, y_max)
        else:
            min_val = min(exp_m2s.min(), 0.0001)
            max_val = max(exp_p2s.max() * 2.5, 0.01)
            ax.set_ylim(min_val * 0.5, max_val)
    else:
        if y_min is not None and y_max is not None:
            ax.set_ylim(y_min, y_max)
        else:
            max_val = exp_p2s.max() * 1.35
            ax.set_ylim(0.0, max_val)

    # Legend
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9, fontsize=16)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    pdf_path = os.path.splitext(output_path)[0] + ".pdf"
    plt.savefig(pdf_path)
    plt.close()

    print(f"[+] Saved limit plot to:\n    - PNG: {output_path}\n    - PDF: {pdf_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot CMS 95% CL CLs Upper Limits (Brazilian Plot)")
    parser.add_argument("--dir", "-d", default=".", help="Directory containing higgsCombine*.root files (default: .)")
    parser.add_argument("--csv", "-c", default=None, help="Input CSV file (optional, overrides --dir)")
    parser.add_argument("--output", "-o", default="limit_cls.png", help="Output plot filename (default: limit_cls.png)")
    parser.add_argument("--lumi", type=float, default=110.0, help="Integrated luminosity in fb^-1 (default: 110.0)")
    parser.add_argument("--com", type=float, default=13.6, help="Center of mass energy in TeV (default: 13.6)")
    parser.add_argument("--status", default="Preliminary", help="CMS label status (default: Preliminary)")
    parser.add_argument("--unblind", action="store_true", help="Plot observed limit if present")
    parser.add_argument("--title-extra", default="", help="Extra annotation text")
    parser.add_argument("--ylabel", "-y", default=None, help="Custom Y axis label")
    parser.add_argument("--percent", action="store_true", help="Display BR in percentage [%]")
    parser.add_argument("--logy", action="store_true", help="Use log scale on y axis (default: linear)")
    args = parser.parse_args()

    if args.csv and os.path.exists(args.csv):
        print(f"Reading limits from CSV: {args.csv}")
        results = extract_limits_from_csv(args.csv)
    else:
        print(f"Reading limits from ROOT files in: {os.path.abspath(args.dir)}")
        results = extract_limits_from_root_files(args.dir)

    if not results:
        print("Error: No limit data found.")
        sys.exit(1)

    plot_brazilian_limit(
        results,
        output_path=args.output,
        lumi=args.lumi,
        com=args.com,
        status=args.status,
        show_observed=args.unblind,
        title_extra=args.title_extra,
        y_label=args.ylabel,
        logy=args.logy,
        is_percent=args.percent,
    )


if __name__ == "__main__":
    main()
