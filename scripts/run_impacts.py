#!/usr/bin/env python3
"""
Nuisance Parameter Impacts Evaluation & Plotting Module for ATobb Combine Analysis.

Workflows:
  1. CR-Only Mode (Real Data Fit):
     - Stage 1 (Bkg Snapshot): Fit CR with BR=0 frozen on real data to obtain post-fit rateParams (ttbbXsec, ttccXsec)
       and nuisance pulls, saving the 'MultiDimFit' snapshot.
     - Stage 2 (Impacts on BR): Load snapshot, set POI (BR), and run parameter scans using -t -1 --toysFrequentist
       with symmetric POI range (rMin=-0.05, rMax=0.05) to evaluate impacts on BR without distorting background fit.
  2. Combined CR+SR Mode (Pure Blinded Asimov Toy):
     - Directly evaluates impacts on POI (BR) using -t -1 on the master workspace without any snapshot or real data.
  3. SR-Only Mode (Pure Blinded Asimov Toy):
     - Directly evaluates impacts on POI (BR) using -t -1 on the masked master workspace without any snapshot.

Outputs:
  - Saved in dedicated per-mass subdirectories:
    * outputs/baseline/impacts_cr_only/MA<X>/
    * outputs/baseline/impacts_combined/MA<X>/

Usage:
  python3 scripts/run_impacts.py [--mass all] [--mode {cr_only,combined,sr_only}] [--poi BR] [--parallel 8]
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import BASELINE_DIR, MASS_POINTS


def run_impacts_for_mass(
    workspace_path,
    output_dir,
    ma_point="MA18",
    mode="cr_only",
    poi="BR",
    parallel=8,
    blind=False,
    br_inj=0.001,
    r_min=None,
    r_max=None,
    extra_freeze="",
):
    """Run Impacts pipeline with symmetric POI bounds, organized per-mass subdirectory."""
    mass_outdir = os.path.join(output_dir, ma_point)
    os.makedirs(mass_outdir, exist_ok=True)
    abs_ws = os.path.abspath(workspace_path)
    abs_out = os.path.abspath(mass_outdir)

    mass_val = int(ma_point.replace("MA", ""))

    # Determine POI if not provided
    if poi is None:
        poi = "BR"

    # Default POI ranges allowing negative values for symmetric parabolic impacts
    if poi == "BR":
        if r_min is None:
            r_min = -0.05
        if r_max is None:
            r_max = 0.05
        range_flags = f"--setParameterRanges BR={r_min},{r_max} --rMin {r_min} --rMax {r_max}"
    else:
        if r_min is None:
            r_min = 0.0
        if r_max is None:
            r_max = 5.0
        range_flags = f"--setParameterRanges {poi}={r_min},{r_max} --rMin {r_min} --rMax {r_max}"

    name_tag = f"{ma_point}_{mode}_{poi}"

    # Channel masking and workflow selection
    if mode == "cr_only":
        freeze_str = "mask_Muon_SR,mask_Electron_SR"
        if extra_freeze:
            freeze_str += f",{extra_freeze}"
        mask_flags = f"--setParameters mask_Muon_SR=1,mask_Electron_SR=1 --freezeParameters {freeze_str}"
        stage1_mask = f"--setParameters mask_Muon_SR=1,mask_Electron_SR=1,BR=0 --freezeParameters {freeze_str},BR"
        data_label = "BLINDED (Asimov Toy)" if blind else "UNBLINDED (Real CR Data with BR=0 frozen snapshot)"
        use_snapshot = True
        snapshot_ws_name = f"higgsCombine_bkgSnapshot_{name_tag}.MultiDimFit.mH{mass_val}.root"
        snapshot_ws_path = os.path.join(abs_out, snapshot_ws_name)
    elif mode == "sr_only":
        freeze_str = "mask_Muon_CR,mask_Electron_CR"
        if extra_freeze:
            freeze_str += f",{extra_freeze}"
        mask_flags = f"--setParameters mask_Muon_CR=1,mask_Electron_CR=1 --freezeParameters {freeze_str}"
        data_label = f"BLINDED (Pure Asimov Toy -t -1, injected BR={br_inj})"
        use_snapshot = False
    else:  # combined
        freeze_str = extra_freeze if extra_freeze else ""
        mask_flags = f"--freezeParameters {freeze_str}" if freeze_str else ""
        data_label = f"BLINDED (Pure Asimov Toy -t -1, injected BR={br_inj})"
        use_snapshot = False

    print("\n" + "=" * 95)
    print(f" Computing Nuisance Parameter Impacts for {ma_point} (Mode: {mode.upper()}, POI: {poi})")
    if use_snapshot:
        print(f" Workflow: 2-Stage Snapshot (Stage 1: CR Real Data Fit -> Stage 2: Symmetric Impact on {poi} with -t -1)")
    else:
        print(f" Workflow: Direct Blinded Asimov Impact (-t -1 on master workspace, NO snapshot)")
    print(f" POI Range: [{r_min}, {r_max}] (Symmetric 2-sided impacts enabled)")
    print(f" Data Mode: {data_label}")
    print(f" Workspace: {abs_ws}")
    print(f" Target Output Directory: {abs_out}")
    print(f" Parallel Cores: {parallel}")
    print("=" * 95)

    if use_snapshot:
        # Step 1: Stage 1 Background Fit on Real CR Data with BR=0 frozen, saving snapshot 'MultiDimFit'
        print(f"\n[+] Step 1: Running Stage 1 Background Fit (BR=0 frozen) & Saving Snapshot 'MultiDimFit'...")
        bkg_cmd = (
            f"combine -M MultiDimFit -d {abs_ws} -m {mass_val} "
            f"{stage1_mask} "
            f"--robustFit 1 --cminDefaultMinimizerStrategy 0 "
            f"--algo none --saveWorkspace -n _bkgSnapshot_{name_tag}"
        )
        print(f"    Command: {bkg_cmd}")
        res_bkg = subprocess.run(bkg_cmd, shell=True, cwd=abs_out, capture_output=True, text=True)
        if res_bkg.returncode != 0 or not os.path.exists(snapshot_ws_path):
            print(f"[-] Stage 1 background fit failed:\n{res_bkg.stderr}\n{res_bkg.stdout}")
            return False
        print(f"    [+] Background snapshot saved successfully: {snapshot_ws_path}")
        ws_for_scans = snapshot_ws_name
        snapshot_opt = "--snapshotName MultiDimFit --toysFrequentist"
    else:
        ws_for_scans = abs_ws
        snapshot_opt = ""

    # Step 2a: Initial Fit for POI (BR)
    print(f"\n[+] Step 2a: Running Initial Fit for POI={poi} with -t -1 ({range_flags})...")
    init_cmd = (
        f"combineTool.py -M Impacts -d {ws_for_scans} -m {mass_val} "
        f"{snapshot_opt} "
        f"--redefineSignalPOIs {poi} "
        f"{range_flags} {mask_flags} "
        f"-t -1 "
        f"--doInitialFit --robustFit 1 --cminDefaultMinimizerStrategy 0 "
        f"-n {name_tag}"
    )
    print(f"    Command: {init_cmd}")
    res_init = subprocess.run(init_cmd, shell=True, cwd=abs_out, capture_output=True, text=True)
    if res_init.returncode != 0:
        print(f"[-] Initial POI fit failed:\n{res_init.stderr}")
        return False

    # Step 2b: Fit Individual Nuisance Parameters
    print(f"\n[+] Step 2b: Running parallel parameter scans with -t -1 ({parallel} threads)...")
    fits_cmd = (
        f"combineTool.py -M Impacts -d {ws_for_scans} -m {mass_val} "
        f"{snapshot_opt} "
        f"--redefineSignalPOIs {poi} "
        f"{range_flags} {mask_flags} "
        f"-t -1 "
        f"--doFits --parallel {parallel} --robustFit 1 --cminDefaultMinimizerStrategy 0 "
        f"-n {name_tag}"
    )
    print(f"    Command: {fits_cmd}")
    res_fits = subprocess.run(fits_cmd, shell=True, cwd=abs_out, capture_output=True, text=True)
    if res_fits.returncode != 0:
        print(f"[-] Nuisance parameter fits failed:\n{res_fits.stderr}")
        return False

    # Step 3: Collect JSON
    print(f"\n[+] Step 3: Aggregating impacts into JSON...")
    json_name = f"impacts_{name_tag}.json"
    json_cmd = (
        f"combineTool.py -M Impacts -d {ws_for_scans} -m {mass_val} "
        f"{snapshot_opt} "
        f"--redefineSignalPOIs {poi} "
        f"-n {name_tag} "
        f"-o {json_name}"
    )
    res_json = subprocess.run(json_cmd, shell=True, cwd=abs_out, capture_output=True, text=True)
    if res_json.returncode != 0:
        print(f"[-] Failed to collect JSON:\n{res_json.stderr}")
        return False

    # Step 4: Plot Impacts PDF and PNG
    print(f"\n[+] Step 4: Rendering CMS Impact Plot...")
    plot_name = f"impacts_{name_tag}"
    plot_cmd = f"plotImpacts.py -i {json_name} -o {plot_name}"
    res_plot = subprocess.run(plot_cmd, shell=True, cwd=abs_out, capture_output=True, text=True)
    if res_plot.returncode != 0:
        print(f"[-] plotImpacts.py failed:\n{res_plot.stderr}")
        return False

    # Remove any existing PNGs before conversion to guarantee clean overwrite
    for old_png in glob.glob(os.path.join(abs_out, f"{plot_name}*.png")):
        try:
            os.remove(old_png)
        except OSError:
            pass

    # Convert PDF pages to PNG
    pdf_path = os.path.join(abs_out, f"{plot_name}.pdf")
    png_path = os.path.join(abs_out, f"{plot_name}.png")
    subprocess.run(
        f"convert -density 150 {plot_name}.pdf {plot_name}.png 2>/dev/null || true",
        shell=True,
        cwd=abs_out,
    )
    # Ensure plot_name.png points to the main overview page (page 0)
    cand_files = sorted(glob.glob(os.path.join(abs_out, f"{plot_name}-*.png")))
    if cand_files:
        shutil.copyfile(cand_files[0], png_path)

    # Cleanup temporary parameter fit ROOT files to keep directory tidy
    for tmp_f in glob.glob(os.path.join(abs_out, f"higgsCombine*Fit*_{name_tag}*.root")):
        try:
            os.remove(tmp_f)
        except OSError:
            pass

    print(f"\n[+] Successfully generated symmetric impacts for {ma_point} ({mode.upper()}, POI={poi})!")
    print(f"    - JSON:                {os.path.join(abs_out, json_name)}")
    print(f"    - PDF:                 {pdf_path}")
    if os.path.exists(png_path):
        print(f"    - PNG (Overview):      {png_path}")

    return True


def main():
    parser = argparse.ArgumentParser(description="Evaluate and plot nuisance parameter impacts using 2-Stage Snapshot / Direct Asimov workflow")
    parser.add_argument("--mass", default="all", help="Target mass point (e.g. MA18, MA25, MA125, or all, default: all)")
    parser.add_argument("--mode", choices=["cr_only", "combined", "sr_only"], default="cr_only", help="Channel mode (default: cr_only)")
    parser.add_argument("--poi", default="BR", help="Parameter of Interest (default: BR)")
    parser.add_argument("--parallel", "-j", type=int, default=8, help="Number of parallel fitting threads (default: 8)")
    parser.add_argument("--ws-dir", default=os.path.join(BASELINE_DIR, "workspaces"), help="Path to master workspaces")
    parser.add_argument("--outdir", default=None, help="Output directory for impacts (default: outputs/baseline/impacts_<mode>)")
    parser.add_argument("--blind", action="store_true", help="Use Asimov toy for initial CR fit (default: False, real data)")
    parser.add_argument("--br-inj", type=float, default=0.001, help="Injected BR value for Asimov toy (default: 0.001)")
    parser.add_argument("--rMin", type=float, default=None, help="Minimum POI value (e.g. -0.05 or -1.0, default: -0.05 for BR)")
    parser.add_argument("--rMax", type=float, default=None, help="Maximum POI value (e.g. 0.05 or 1.0, default: 0.05 for BR)")
    parser.add_argument("--freeze", default="", help="Extra parameters to freeze (comma-separated)")
    args = parser.parse_args()

    out_dir = os.path.join(BASELINE_DIR, f"impacts_{args.mode}") if args.outdir is None else args.outdir
    mass_list = MASS_POINTS if args.mass == "all" else [args.mass]

    for ma in mass_list:
        ws_file = os.path.join(args.ws_dir, f"workspace_{ma}.root")
        if not os.path.exists(ws_file):
            print(f"[-] Workspace not found: {ws_file}")
            continue

        run_impacts_for_mass(
            workspace_path=ws_file,
            output_dir=out_dir,
            ma_point=ma,
            mode=args.mode,
            poi=args.poi,
            parallel=args.parallel,
            blind=args.blind,
            br_inj=args.br_inj,
            r_min=args.rMin,
            r_max=args.rMax,
            extra_freeze=args.freeze,
        )


if __name__ == "__main__":
    main()
