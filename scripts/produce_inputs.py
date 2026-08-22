#!/usr/bin/env python3
"""
Main Input Producer for ATobb Combine Analysis.
Extracts shapes, builds datacards, and compiles RooWorkspaces (with channel masks):
  - Baseline Master Workspaces: 4 categories (Muon_CR, Electron_CR, Muon_SR, Electron_SR)
    with --channel-masks enabled for reusability in:
      * Combined Simultaneous Fit
      * SR-Only Limits (via masking CR channels)
      * CR-Only Fit & GoF (via masking SR channels)

Usage:
  python3 scripts/produce_inputs.py [--outdir DIR] [--cr-cut 0.3] [--rebin 2] [--mass all]
"""

import argparse
import os
import sys

# Add repository root to python path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from configs.analysis_config import (
    BASELINE_DIR,
    CR_CUT_DEFAULT,
    DEFAULT_REBIN,
    MASS_POINTS,
)
from modules.shape_extractor import extract_shapes_for_categories
from modules.datacard_generator import write_datacard
from modules.workspace_builder import build_workspace


def produce_baseline_inputs(
    output_dir=BASELINE_DIR,
    mass_list=MASS_POINTS,
    cr_cut=CR_CUT_DEFAULT,
    rebin_factor=DEFAULT_REBIN,
    sr_var_mode="avg",
):
    """
    Generate baseline 4-category simultaneous inputs (Muon_CR, Electron_CR, Muon_SR, Electron_SR)
    with channel masks enabled for maximal reusability.
    """
    print("\n" + "=" * 80)
    print(f" [Baseline] Generating Master Workspaces (CR <= {cr_cut}, SR > {cr_cut})")
    print(f" Output Directory: {output_dir}")
    print("=" * 80)

    shapes_dir = os.path.join(output_dir, "shapes")
    cards_dir = os.path.join(output_dir, "datacards")
    ws_dir = os.path.join(output_dir, "workspaces")
    for d in [shapes_dir, cards_dir, ws_dir]:
        os.makedirs(d, exist_ok=True)

    for ma in mass_list:
        sr_var = f"bdt_score_{ma}" if sr_var_mode == "mass_bdt" else "bdt_score_avg"
        categories = [
            {"cat_name": "Muon_CR", "channel": "Muon", "var_name": "bdt_score_avg", "x_min": 0.0, "x_max": cr_cut},
            {"cat_name": "Electron_CR", "channel": "Electron", "var_name": "bdt_score_avg", "x_min": 0.0, "x_max": cr_cut},
            {"cat_name": "Muon_SR", "channel": "Muon", "var_name": sr_var, "x_min": cr_cut, "x_max": 1.0},
            {"cat_name": "Electron_SR", "channel": "Electron", "var_name": sr_var, "x_min": cr_cut, "x_max": 1.0},
        ]

        shape_file = os.path.join(shapes_dir, f"shapes_{ma}.root")
        card_file = os.path.join(cards_dir, f"datacard_{ma}.txt")
        ws_file = os.path.join(ws_dir, f"workspace_{ma}.root")

        # 1. Extract shapes
        yields = extract_shapes_for_categories(categories, shape_file, ma_point=ma, rebin_factor=rebin_factor)

        # 2. Write 4-channel datacard
        write_datacard(
            categories=["Muon_CR", "Electron_CR", "Muon_SR", "Electron_SR"],
            shapes_root_path=shape_file,
            yields=yields,
            output_card_path=card_file,
            ma_point=ma,
            title_comment=f"Baseline Master Datacard (CR<={cr_cut}, SR>{cr_cut})",
        )

        # 3. Build workspace with --channel-masks and ChargedHiggsBR
        build_workspace(card_file, ws_file, channel_masks=True)
        print(f"  [+] [{ma}] Successfully built workspace -> {ws_file}")

    print("\n[+] All baseline inputs generated successfully!")


def main():
    parser = argparse.ArgumentParser(description="Produce Master Combine inputs for ATobb Analysis")
    parser.add_argument("--outdir", default=BASELINE_DIR, help=f"Output directory (default: {BASELINE_DIR})")
    parser.add_argument("--cr-cut", type=float, default=CR_CUT_DEFAULT, help=f"CR/SR boundary cut on bdt_score_avg (default: {CR_CUT_DEFAULT})")
    parser.add_argument("--rebin", type=int, default=DEFAULT_REBIN, help=f"Rebin factor (default: {DEFAULT_REBIN})")
    parser.add_argument("--mass", default="all", help="Target mass point (MA18..MA125 or all, default: all)")
    parser.add_argument("--sr-var", choices=["avg", "mass_bdt"], default="avg", help="SR discriminator: avg or mass_bdt (default: avg)")
    args = parser.parse_args()

    mass_list = MASS_POINTS if args.mass == "all" else [args.mass]
    produce_baseline_inputs(
        output_dir=args.outdir,
        mass_list=mass_list,
        cr_cut=args.cr_cut,
        rebin_factor=args.rebin,
        sr_var_mode=args.sr_var,
    )


if __name__ == "__main__":
    main()
