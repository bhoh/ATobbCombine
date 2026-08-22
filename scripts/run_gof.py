#!/usr/bin/env python3
"""
Goodness-of-Fit (GoF) Test Module for ATobb Combine Analysis.
Runs Saturated GoF test with Signal Region (SR) masked:
  --setParameters mask_Muon_SR=1,mask_Electron_SR=1
  --freezeParameters mask_Muon_SR,mask_Electron_SR
Calculates p-value from toys and produces publication-quality CMS GoF plots.

Usage:
  python3 scripts/run_gof.py [--workspace WORKSPACE.root] [--toys 100] [--algo saturated] [--outdir DIR] [--seed 1234]
"""

import argparse
import os
import subprocess
import sys
import numpy as np
import ROOT

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

from configs.analysis_config import BASELINE_DIR, COM_TEV, DEFAULT_OUTPUT_DIR, LUMI_FB


def extract_q_from_file(root_file):
    """Extract limit/test-statistic values from a GoodnessOfFit output ROOT file."""
    if not os.path.exists(root_file):
        return []
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        return []
    tree = f.Get("limit")
    if not tree:
        f.Close()
        return []
    values = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        values.append(tree.limit)
    f.Close()
    return values


def run_gof(
    workspace_path,
    output_dir,
    n_toys=100,
    algo="saturated",
    masked_channels=("mask_Muon_SR", "mask_Electron_SR"),
    seed=123456,
    frequentist=False,
    unblind=False,
):
    """
    Execute Goodness-of-Fit test with masked SR categories.
    """
    os.makedirs(output_dir, exist_ok=True)
    ws_abs = os.path.abspath(workspace_path)
    work_dir = os.path.abspath(output_dir)

    mask_params = ",".join(f"{ch}=1" for ch in masked_channels)
    freeze_params = ",".join(masked_channels)

    print("\n" + "=" * 80)
    print(f" Goodness-of-Fit Test: Algorithm = {algo.upper()} (Masked SR)")
    print(f" Workspace: {ws_abs}")
    print(f" Masked Channels: {masked_channels}")
    print(f" Toys: {n_toys}")
    print("=" * 80)

    # 1. Run Data / Observed GoF
    data_tag = "data"
    data_cmd = (
        f"combine -M GoodnessOfFit {ws_abs} --algo {algo} "
        f"--setParameters {mask_params} --freezeParameters {freeze_params} "
        f"--fixedSignalStrength 0 -n _{data_tag}"
    )
    print(f"\n[+] Running Observed GoF test:\n    {data_cmd}")
    res_data = subprocess.run(data_cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
    if res_data.returncode != 0:
        print(f"[-] Error running data GoF:\n{res_data.stderr}")
        return None

    data_file = os.path.join(work_dir, f"higgsCombine_{data_tag}.GoodnessOfFit.mH120.root")
    data_q_vals = extract_q_from_file(data_file)
    if not data_q_vals:
        print("[-] Could not extract observed test statistic.")
        return None
    q_obs = data_q_vals[0]
    print(f"  [+] Observed Test Statistic q_obs = {q_obs:.4f}")

    # 2. Run Toys GoF
    toys_tag = "toys"
    freq_flag = "--toysFrequentist" if frequentist else ""
    toys_cmd = (
        f"combine -M GoodnessOfFit {ws_abs} --algo {algo} -t {n_toys} "
        f"{freq_flag} "
        f"--setParameters {mask_params} --freezeParameters {freeze_params} "
        f"--fixedSignalStrength 0 -s {seed} -n _{toys_tag}"
    )
    print(f"\n[+] Running {n_toys} GoF Pseudo-experiments (Toys):\n    {toys_cmd}")
    res_toys = subprocess.run(toys_cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
    if res_toys.returncode != 0:
        print(f"[-] Error running toys GoF:\n{res_toys.stderr}")
        return None

    toys_file = os.path.join(work_dir, f"higgsCombine_{toys_tag}.GoodnessOfFit.mH120.{seed}.root")
    toy_q_vals = np.array(extract_q_from_file(toys_file))
    if len(toy_q_vals) == 0:
        print("[-] Could not extract toy test statistics.")
        return None

    # 3. Calculate p-value
    n_extreme = np.sum(toy_q_vals >= q_obs)
    p_val = n_extreme / len(toy_q_vals)

    print("\n" + "-" * 60)
    print(" Goodness-of-Fit Summary")
    print("-" * 60)
    print(f"  Algorithm:            {algo.upper()}")
    print(f"  Observed q_data:      {q_obs:.4f}")
    print(f"  Toys Mean +/- Std:    {np.mean(toy_q_vals):.4f} +/- {np.std(toy_q_vals):.4f}")
    print(f"  Toys Median:          {np.median(toy_q_vals):.4f}")
    print(f"  Extreme Toys (>=q):   {n_extreme} / {len(toy_q_vals)}")
    print(f"  p-value:              {p_val:.4f} ({p_val * 100:.1f}%)")
    print("-" * 60)

    # 4. Plot distribution
    plot_gof_distribution(
        q_obs=q_obs,
        toy_q_vals=toy_q_vals,
        p_val=p_val,
        algo=algo,
        output_dir=output_dir,
    )

    return {
        "algo": algo,
        "q_obs": q_obs,
        "toy_mean": float(np.mean(toy_q_vals)),
        "toy_std": float(np.std(toy_q_vals)),
        "toy_median": float(np.median(toy_q_vals)),
        "p_value": p_val,
        "n_toys": len(toy_q_vals),
    }


def plot_gof_distribution(q_obs, toy_q_vals, p_val, algo, output_dir):
    """Generate CMS Brazilian-style Goodness-of-Fit distribution plot."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Determine binning
    min_x = min(np.min(toy_q_vals), q_obs) * 0.8
    max_x = max(np.max(toy_q_vals), q_obs) * 1.2
    bins = np.linspace(min_x, max_x, 30)

    # Plot toys histogram
    counts, edges, _ = ax.hist(
        toy_q_vals,
        bins=bins,
        histtype="stepfilled",
        color="#708090",
        alpha=0.4,
        edgecolor="#2F4F4F",
        linewidth=1.8,
        label=f"Expected ({len(toy_q_vals)} toys)",
    )

    # Vertical line for observed
    ax.axvline(
        q_obs,
        color="#D9531E",
        linestyle="--",
        linewidth=2.8,
        label=f"Observed ($q_{{\\mathrm{{data}}}} = {q_obs:.2f}$)",
        zorder=5,
    )

    # CMS header
    if HAS_MPLHEP:
        hep.cms.label(
            llabel="Preliminary",
            data=True,
            lumi=int(LUMI_FB),
            com=COM_TEV,
            loc=0,
            ax=ax,
        )
    else:
        ax.set_title(
            f"CMS Preliminary   {int(LUMI_FB)} fb$^{{-1}}$ ({COM_TEV} TeV, 2024)",
            loc="right",
            fontsize=16,
        )

    # Physics process & GoF annotations
    ax.text(
        0.04,
        0.90,
        r"$t \to H^{+}b,\, H^{+} \to W^{+}A\,(\to b\bar{b})$",
        transform=ax.transAxes,
        fontsize=20,
        verticalalignment="top",
    )
    ax.text(
        0.04,
        0.83,
        f"Goodness-of-Fit ({algo.capitalize()} Test)",
        transform=ax.transAxes,
        fontsize=18,
        verticalalignment="top",
    )
    ax.text(
        0.04,
        0.77,
        "Signal Region Masked (CR-only GoF)",
        transform=ax.transAxes,
        fontsize=16,
        verticalalignment="top",
        color="darkblue",
    )

    # Results box
    res_text = (
        f"$q_{{\\mathrm{{data}}}} = {q_obs:.2f}$\n"
        f"$\\mu_{{\\mathrm{{toys}}}} = {np.mean(toy_q_vals):.2f} \\pm {np.std(toy_q_vals):.2f}$\n"
        f"$p\\text{{-value}} = {p_val:.3f}$"
    )
    ax.text(
        0.65,
        0.65,
        res_text,
        transform=ax.transAxes,
        fontsize=17,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="gray", alpha=0.9),
    )

    ax.set_xlabel(f"{algo.capitalize()} Test Statistic $q$", fontsize=22)
    ax.set_ylabel("Pseudo-experiments / Bin", fontsize=22)
    ax.set_ylim(0, max(counts) * 1.35)
    ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.9, fontsize=16)

    plt.tight_layout()
    out_png = os.path.join(output_dir, f"gof_{algo}_maskedSR.png")
    out_pdf = os.path.join(output_dir, f"gof_{algo}_maskedSR.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()

    print(f"\n[+] Saved GoF distribution plot to:\n    - PNG: {out_png}\n    - PDF: {out_pdf}")


def main():
    parser = argparse.ArgumentParser(description="Run Goodness-of-Fit test with masked SR")
    parser.add_argument(
        "--workspace",
        "-w",
        default=os.path.join(BASELINE_DIR, "workspaces/workspace_MA18.root"),
        help="Input Combined workspace with channel masks",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default=os.path.join(BASELINE_DIR, "gof"),
        help="Output directory for GoF results",
    )
    parser.add_argument("--toys", "-t", type=int, default=100, help="Number of toy pseudo-experiments (default: 100)")
    parser.add_argument("--algo", default="saturated", choices=["saturated", "KS", "AD"], help="GoF algorithm (default: saturated)")
    parser.add_argument("--seed", "-s", type=int, default=123456, help="Random seed for toys (default: 123456)")
    parser.add_argument("--frequentist", action="store_true", help="Use --toysFrequentist (fits nuisance parameters per toy)")
    parser.add_argument("--unblind", action="store_true", help="Run on real observed data")
    args = parser.parse_args()

    run_gof(
        workspace_path=args.workspace,
        output_dir=args.outdir,
        n_toys=args.toys,
        algo=args.algo,
        seed=args.seed,
        frequentist=args.frequentist,
        unblind=args.unblind,
    )


if __name__ == "__main__":
    main()
