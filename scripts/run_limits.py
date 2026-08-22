#!/usr/bin/env python3
"""
Calculate 95% CL Asymptotic Limits for ATobb Combine Analysis.
Supports:
  - Combined CR+SR Fit (All 4 channels unmasked)
  - SR-Only Limits (CR channels masked: mask_Muon_CR=1, mask_Electron_CR=1)
  - Reuses the baseline master workspaces generated with --channel-masks.

Usage:
  python3 scripts/run_limits.py [--mode {combined,sr_only,both}] [--ws-dir DIR] [--outdir DIR] [--mass all] [--unblind]
"""

import argparse
import csv
import os
import re
import subprocess
import sys

# Add repository root to python path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import BASELINE_DIR, MASS_POINTS


def run_asymptotic_limit(workspace_path, work_dir, ma_point, masked_channels=None, unblind=False):
    """Run Combine AsymptoticLimits on a single workspace with optional channel masking."""
    abs_ws = os.path.abspath(workspace_path)
    abs_work = os.path.abspath(work_dir)
    mask_flags = ""
    if masked_channels:
        mask_params = ",".join(f"{ch}=1" for ch in masked_channels)
        freeze_params = ",".join(masked_channels)
        mask_flags = f"--setParameters {mask_params} --freezeParameters {freeze_params}"

    cmd = (
        f"combine -M AsymptoticLimits {abs_ws} "
        f"{'--run blind' if not unblind else ''} "
        f"{mask_flags} "
        f"--rAbsAcc 0 --rRelAcc 0.0005 -n _{ma_point}"
    )
    res = subprocess.run(cmd, shell=True, cwd=abs_work, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[-] Error running Combine for {ma_point}:\n{res.stderr}")
        return None

    exp_2p5 = re.search(r"Expected\s+2\.5%:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)
    exp_16 = re.search(r"Expected\s+16\.0%:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)
    exp_50 = re.search(r"Expected\s+50\.0%:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)
    exp_84 = re.search(r"Expected\s+84\.0%:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)
    exp_97p5 = re.search(r"Expected\s+97\.5%:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)
    obs = re.search(r"Observed Limit:\s+(?:r|BR)\s*<\s*([\d\.]+)", res.stdout)

    if exp_50:
        return {
            "mass_point": ma_point,
            "mass_GeV": float(ma_point.replace("MA", "")),
            "exp_minus2s": float(exp_2p5.group(1)) if exp_2p5 else None,
            "exp_minus1s": float(exp_16.group(1)) if exp_16 else None,
            "exp_median": float(exp_50.group(1)),
            "exp_plus1s": float(exp_84.group(1)) if exp_84 else None,
            "exp_plus2s": float(exp_97p5.group(1)) if exp_97p5 else None,
            "observed": float(obs.group(1)) if (obs and unblind) else None,
        }
    return None


def run_limits_for_directory(
    ws_dir,
    work_dir,
    label="Combined CR+SR",
    mass_list=MASS_POINTS,
    masked_channels=None,
    unblind=False,
):
    """Run limits across all mass points in a workspace directory."""
    os.makedirs(work_dir, exist_ok=True)
    print("\n" + "=" * 90)
    print(f" Computing 95% CL Asymptotic Limits: {label} ({'UNBLINDED' if unblind else 'BLINDED'})")
    if masked_channels:
        print(f" Masked Channels: {masked_channels}")
    print("=" * 90)

    results = []
    for ma in mass_list:
        ws_file = os.path.join(ws_dir, f"workspace_{ma}.root")
        if not os.path.exists(ws_file):
            print(f"[-] Warning: Workspace {ws_file} not found. Skipping {ma}.")
            continue

        lim = run_asymptotic_limit(
            ws_file,
            work_dir,
            ma,
            masked_channels=masked_channels,
            unblind=unblind,
        )
        if lim:
            results.append(lim)
            obs_str = f"Obs: {lim['observed']:.6f}" if lim["observed"] is not None else "BLINDED"
            print(f"  {ma:<8}: Expected Median = {lim['exp_median']:.6f} ([-1s, +1s]: [{lim['exp_minus1s']:.6f}, {lim['exp_plus1s']:.6f}]) | {obs_str}")

    results.sort(key=lambda x: x["mass_GeV"])
    return results


def print_summary_table(results, title="95% CL Limits Summary"):
    """Display an ASCII summary table in the terminal."""
    if not results:
        return
    print("\n" + "=" * 115)
    print(f" {title}")
    print("=" * 115)
    print(f"{'Mass Point':<12} | {'Mass [GeV]':<10} | {'-2 sigma':<12} | {'-1 sigma':<12} | {'Median BR':<14} | {'+1 sigma':<12} | {'+2 sigma':<12} | {'Observed':<10}")
    print("-" * 115)
    for r in results:
        obs_str = f"{r['observed']:.6f}" if r["observed"] is not None else "BLINDED"
        m2 = f"{r['exp_minus2s']:.6f}" if r["exp_minus2s"] is not None else "-"
        m1 = f"{r['exp_minus1s']:.6f}" if r["exp_minus1s"] is not None else "-"
        med = f"{r['exp_median']:.6f}" if r["exp_median"] is not None else "-"
        p1 = f"{r['exp_plus1s']:.6f}" if r["exp_plus1s"] is not None else "-"
        p2 = f"{r['exp_plus2s']:.6f}" if r["exp_plus2s"] is not None else "-"
        print(f"{r['mass_point']:<12} | {r['mass_GeV']:<10.1f} | {m2:<12} | {m1:<12} | {med:<14} | {p1:<12} | {p2:<12} | {obs_str:<10}")
    print("=" * 115 + "\n")


def write_csv(results, output_csv_path):
    """Write results dictionary to a clean CSV file."""
    if not results:
        return
    fieldnames = ["mass_point", "mass_GeV", "exp_minus2s", "exp_minus1s", "exp_median", "exp_plus1s", "exp_plus2s", "observed"]
    with open(output_csv_path, mode="w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
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
            })
    print(f"[+] Saved limit results to CSV -> {output_csv_path}")


def main():
    parser = argparse.ArgumentParser(description="Calculate Asymptotic Limits for ATobb Combine Analysis")
    parser.add_argument("--mode", choices=["combined", "sr_only", "both"], default="both", help="Analysis mode (default: both)")
    parser.add_argument("--ws-dir", default=os.path.join(BASELINE_DIR, "workspaces"), help="Directory with master workspaces")
    parser.add_argument("--outdir", default=BASELINE_DIR, help="Base output directory")
    parser.add_argument("--unblind", action="store_true", help="Run unblinded calculation on data")
    parser.add_argument("--mass", default="all", help="Target mass point (MA18..MA125 or all, default: all)")
    args = parser.parse_args()

    mass_list = MASS_POINTS if args.mass == "all" else [args.mass]

    if args.mode in ["combined", "both"]:
        comb_out_dir = os.path.join(args.outdir, "limits_combined")
        res_comb = run_limits_for_directory(
            args.ws_dir,
            comb_out_dir,
            label="Combined CR+SR",
            mass_list=mass_list,
            masked_channels=None,
            unblind=args.unblind,
        )
        print_summary_table(res_comb, title="Combined CR+SR 95% CL Upper Limits on BR(t -> H+ b)")
        write_csv(res_comb, os.path.join(comb_out_dir, "limits_summary.csv"))

    if args.mode in ["sr_only", "both"]:
        sr_out_dir = os.path.join(args.outdir, "limits_sr_only")
        res_sr = run_limits_for_directory(
            args.ws_dir,
            sr_out_dir,
            label="SR-Only (Masked CR)",
            mass_list=mass_list,
            masked_channels=["mask_Muon_CR", "mask_Electron_CR"],
            unblind=args.unblind,
        )
        print_summary_table(res_sr, title="SR-Only (Masked CR) 95% CL Upper Limits on BR(t -> H+ b)")
        write_csv(res_sr, os.path.join(sr_out_dir, "limits_summary.csv"))


if __name__ == "__main__":
    main()
