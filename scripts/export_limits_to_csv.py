#!/usr/bin/env python3
"""
Extract AsymptoticLimits from higgsCombine*.AsymptoticLimits*.root files
and export the results to CSV format.

Usage:
  python3 scripts/export_limits_to_csv.py [--dir DIR] [--pattern PATTERN] [--output OUTPUT.csv] [--recursive]
"""

import argparse
import csv
import glob
import os
import re
import sys
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


def extract_limit_from_root(file_path):
    """Read limit tree from a single higgsCombine ROOT file."""
    if not os.path.exists(file_path):
        return None

    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        return None

    tree = f.Get("limit")
    if not tree or tree.GetEntries() == 0:
        f.Close()
        return None

    limits = {
        "obs": None,
        "exp_2p5": None,
        "exp_16": None,
        "exp_50": None,
        "exp_84": None,
        "exp_97p5": None,
    }

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        q = round(tree.quantileExpected, 4)
        val = float(tree.limit)

        if q < 0:
            limits["obs"] = val
        elif abs(q - 0.025) < 0.005:
            limits["exp_2p5"] = val
        elif abs(q - 0.160) < 0.005:
            limits["exp_16"] = val
        elif abs(q - 0.500) < 0.005:
            limits["exp_50"] = val
        elif abs(q - 0.840) < 0.005:
            limits["exp_84"] = val
        elif abs(q - 0.975) < 0.005:
            limits["exp_97p5"] = val

    f.Close()
    return limits


def process_directory(input_dir, pattern="higgsCombine*.AsymptoticLimits*.root", recursive=False):
    """Find and process all matching ROOT files."""
    if recursive:
        search_pattern = os.path.join(input_dir, "**", pattern)
        file_list = glob.glob(search_pattern, recursive=True)
    else:
        search_pattern = os.path.join(input_dir, pattern)
        file_list = glob.glob(search_pattern)

    valid_files = [f for f in file_list if not os.path.basename(f).startswith("higgsCombineTest.")]
    if not valid_files:
        valid_files = file_list

    results = []
    for fpath in valid_files:
        ma_str, ma_val = parse_mass_from_filename(fpath)
        lim = extract_limit_from_root(fpath)
        if lim and lim["exp_50"] is not None:
            results.append({
                "file": fpath,
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


def write_csv(results, output_csv_path):
    """Write parsed limits to a CSV file."""
    fieldnames = [
        "mass_point",
        "mass_GeV",
        "exp_minus2s",
        "exp_minus1s",
        "exp_median",
        "exp_plus1s",
        "exp_plus2s",
        "observed",
        "file_path",
    ]

    with open(output_csv_path, mode="w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({
                "mass_point": r["mass_point"],
                "mass_GeV": f"{r['mass_GeV']:.1f}",
                "exp_minus2s": f"{r['exp_minus2s']:.6f}" if r["exp_minus2s"] is not None else "",
                "exp_minus1s": f"{r['exp_minus1s']:.6f}" if r["exp_minus1s"] is not None else "",
                "exp_median": f"{r['exp_median']:.6f}" if r["exp_median"] is not None else "",
                "exp_plus1s": f"{r['exp_plus1s']:.6f}" if r["exp_plus1s"] is not None else "",
                "exp_plus2s": f"{r['exp_plus2s']:.6f}" if r["exp_plus2s"] is not None else "",
                "observed": f"{r['observed']:.6f}" if r["observed"] is not None else "BLINDED",
                "file_path": r["file"],
            })

    print(f"[+] Successfully saved {len(results)} limit results to CSV: {output_csv_path}")


def print_table(results):
    """Print ASCII table in terminal."""
    if not results:
        print("No limit results found.")
        return

    print("\n" + "=" * 105)
    print(f"{'Mass Point':<12} | {'Mass [GeV]':<10} | {'-2 sigma':<12} | {'-1 sigma':<12} | {'Median [pb]':<14} | {'+1 sigma':<12} | {'+2 sigma':<12} | {'Observed':<10}")
    print("=" * 105)
    for r in results:
        obs_str = f"{r['observed']:.4f}" if r["observed"] is not None else "BLINDED"
        m2 = f"{r['exp_minus2s']:.4f}" if r["exp_minus2s"] is not None else "-"
        m1 = f"{r['exp_minus1s']:.4f}" if r["exp_minus1s"] is not None else "-"
        med = f"{r['exp_median']:.4f}" if r["exp_median"] is not None else "-"
        p1 = f"{r['exp_plus1s']:.4f}" if r["exp_plus1s"] is not None else "-"
        p2 = f"{r['exp_plus2s']:.4f}" if r["exp_plus2s"] is not None else "-"
        print(f"{r['mass_point']:<12} | {r['mass_GeV']:<10.1f} | {m2:<12} | {m1:<12} | {med:<14} | {p1:<12} | {p2:<12} | {obs_str:<10}")
    print("=" * 105)


def main():
    parser = argparse.ArgumentParser(description="Export Combine AsymptoticLimits ROOT files to CSV")
    parser.add_argument("--dir", "-d", default=".", help="Directory containing higgsCombine*.root files (default: .)")
    parser.add_argument("--pattern", "-p", default="higgsCombine*.AsymptoticLimits*.root", help="Glob pattern for root files")
    parser.add_argument("--output", "-o", default="limits_summary.csv", help="Output CSV file name (default: limits_summary.csv)")
    parser.add_argument("--recursive", "-r", action="store_true", help="Recursively search subdirectories")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.dir)
    output_csv = os.path.abspath(args.output) if not os.path.isabs(args.output) else args.output

    results = process_directory(input_dir, args.pattern, args.recursive)
    if not results:
        print(f"Warning: No matching higgsCombine AsymptoticLimits ROOT files found in {input_dir}")
        sys.exit(1)

    print_table(results)
    write_csv(results, output_csv)


if __name__ == "__main__":
    main()
