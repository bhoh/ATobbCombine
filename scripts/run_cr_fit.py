#!/usr/bin/env python3
"""
CR-Only Fit and Nuisance Pull Analysis for ATobb Combine.
Runs combine FitDiagnostics masking the Signal Region (SR) on the baseline master workspace:
  --setParameters mask_Muon_SR=1,mask_Electron_SR=1
  --freezeParameters mask_Muon_SR,mask_Electron_SR
Measures background rateParams (ttbbXsec, ttccXsec) and extracts post-fit parameters.

Usage:
  python3 scripts/run_cr_fit.py [--workspace WORKSPACE] [--unblind] [--output-dir DIR]
"""

import argparse
import os
import re
import subprocess
import sys
import ROOT

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import BASELINE_DIR

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def run_fit_diagnostics(
    workspace_path,
    output_dir,
    masked_channels=("mask_Muon_SR", "mask_Electron_SR"),
    unblind=False,
    robust=False,
):
    """Execute Combine FitDiagnostics with masked SR categories."""
    os.makedirs(output_dir, exist_ok=True)
    abs_ws = os.path.abspath(workspace_path)

    mask_flags = ""
    if masked_channels:
        mask_params = ",".join(f"{ch}=1" for ch in masked_channels)
        freeze_params = ",".join(masked_channels)
        mask_flags = f"--setParameters {mask_params} --freezeParameters {freeze_params}"

    robust_flag = "--robustFit=1" if robust else "--robustFit=0"
    cmd = (
        f"combine -M FitDiagnostics {abs_ws} "
        f"{'' if unblind else '-t -1 --expectSignal=0'} "
        f"{mask_flags} "
        f"--saveShapes --saveNormalizations "
        f"{robust_flag} --cminDefaultMinimizerStrategy=0 "
        f"-n _CR"
    )

    print(f"\n[+] Running Combine FitDiagnostics on CR (Masked SR):")
    print(f"    Command: {cmd}")
    print(f"    Mode: {'UNBLINDED (Real Data)' if unblind else 'BLINDED (Asimov Toy)'}")

    res = subprocess.run(cmd, shell=True, cwd=output_dir, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[-] FitDiagnostics failed:\n{res.stderr}")
        return False

    fit_root = os.path.join(output_dir, "fitDiagnostics_CR.root")
    if not os.path.exists(fit_root):
        alt_root = os.path.join(output_dir, "fitDiagnostics.root")
        if os.path.exists(alt_root):
            os.rename(alt_root, fit_root)

    print(f"[+] FitDiagnostics completed successfully! Output: {fit_root}")
    return fit_root


def extract_rate_params(fit_root_path):
    """Extract post-fit values for rate parameters from fitDiagnostics output."""
    if not os.path.exists(fit_root_path):
        return {}

    f = ROOT.TFile.Open(fit_root_path)
    if not f or f.IsZombie():
        return {}

    fit_b = f.Get("fit_b")
    if not fit_b:
        f.Close()
        return {}

    params = {}
    pars = fit_b.floatParsFinal()
    for i in range(pars.getSize()):
        p = pars.at(i)
        name = p.GetName()
        if "Xsec" in name or "ratio_" in name:
            params[name] = {
                "val": p.getVal(),
                "err": p.getError(),
                "err_lo": p.getErrorLo(),
                "err_hi": p.getErrorHi(),
            }
    f.Close()
    return params


def print_fit_summary(rate_params):
    """Print clean summary table of post-fit rate parameters."""
    print("\n" + "=" * 60)
    print(" Post-fit Rate Parameters (CR-Only Fit)")
    print("=" * 60)
    print(f"{'Parameter':<24} | {'Pre-fit':<10} | {'Post-fit Value +/- Error':<20}")
    print("-" * 60)
    for name, d in rate_params.items():
        print(f"{name:<24} | {'1.000':<10} | {d['val']:.4f} +/- {d['err']:.4f}")
    print("=" * 60 + "\n")


def main():
    parser = argparse.ArgumentParser(description="CR-Only Background Fit and Nuisance Pull Analysis")
    parser.add_argument(
        "--workspace",
        "-w",
        default=os.path.join(BASELINE_DIR, "workspaces/workspace_MA18.root"),
        help="Input workspace (default: baseline workspace_MA18.root)",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default=os.path.join(BASELINE_DIR, "cr_fit"),
        help="Directory to save fitDiagnostics outputs",
    )
    parser.add_argument("--unblind", action="store_true", help="Run unblinded on real data (default: Asimov toy)")
    parser.add_argument("--robust", action="store_true", help="Use robustFit=1 (slower, for exact errors)")
    args = parser.parse_args()

    fit_root = run_fit_diagnostics(
        workspace_path=args.workspace,
        output_dir=args.output_dir,
        unblind=args.unblind,
        robust=args.robust,
    )

    if fit_root:
        rps = extract_rate_params(fit_root)
        if rps:
            print_fit_summary(rps)


if __name__ == "__main__":
    main()
