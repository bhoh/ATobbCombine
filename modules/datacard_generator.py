"""
Datacard Generator Module for ATobb Combine Analysis.
Constructs Combine datacards for CR-only, SR-only, and combined CR+SR categories.
"""

import os
from configs.analysis_config import (
    ALL_PROCESSES,
    BACKGROUND_PROCESSES,
    CHANNEL_RESTRICTED_SYST,
    LNN_UNCERTAINTIES,
    RATE_PARAMETERS,
    SYSTEMATIC_VARIATIONS,
    TTBAR_ONLY_SYST,
)


def write_datacard(
    categories,
    shapes_root_path,
    yields,
    output_card_path,
    ma_point="MA18",
    title_comment="ATobb Search",
):
    """
    Generate and write a Combine datacard for specified categories.
    categories: list of category names (e.g. ['Muon_CR', 'Electron_CR', 'Muon_SR', 'Electron_SR'])
    """
    procs = ["signal"] + BACKGROUND_PROCESSES
    proc_ids = list(range(len(procs)))  # signal has id 0, backgrounds have ids 1..7

    n_channels = len(categories)
    n_backgrounds = len(BACKGROUND_PROCESSES)

    lines = []
    lines.append(f"## Combine Datacard for {title_comment} (Target: {ma_point})")
    lines.append(f"imax {n_channels} number of channels")
    lines.append("jmax * number of backgrounds")
    lines.append("kmax * number of nuisance parameters")
    lines.append("-" * 100)

    # Shape mapping
    lines.append(f"shapes * * {shapes_root_path} histo_$CHANNEL_$PROCESS histo_$CHANNEL_$PROCESS_$SYSTEMATIC")
    lines.append(f"shapes data_obs * {shapes_root_path} histo_$CHANNEL_data_obs")
    lines.append("-" * 100)

    # Observation header
    ch_hdr = "bin          " + "  ".join(f"{cat:<18}" for cat in categories)
    obs_hdr = "observation  " + "  ".join(f"{int(round(yields.get(cat, {}).get('data_obs', 0.0))):<18}" for cat in categories)
    lines.append(ch_hdr)
    lines.append(obs_hdr)
    lines.append("-" * 100)

    # Process rate table
    bin_cols, proc_name_cols, proc_id_cols, rate_cols = [], [], [], []
    for cat in categories:
        for p, pid in zip(procs, proc_ids):
            rate_val = yields.get(cat, {}).get(p, 0.0)
            bin_cols.append(f"{cat:<18}")
            proc_name_cols.append(f"{p:<18}")
            proc_id_cols.append(f"{pid:<18}")
            rate_cols.append(f"{rate_val:<18.4f}")

    lines.append("bin          " + "  ".join(bin_cols))
    lines.append("process      " + "  ".join(proc_name_cols))
    lines.append("process      " + "  ".join(proc_id_cols))
    lines.append("rate         " + "  ".join(rate_cols))
    lines.append("-" * 100)

    # Normalization (lnN) uncertainties
    for name, info in LNN_UNCERTAINTIES.items():
        val = info["value"]
        target_procs = info["processes"]
        cols = []
        for cat in categories:
            for p in procs:
                if p in target_procs:
                    cols.append(f"{val:<18}")
                else:
                    cols.append(f"{'-':<18}")
        lines.append(f"{name:<24} lnN  " + "  ".join(cols))

    lines.append("-" * 100)

    # Shape systematics
    all_shape_systs = list(SYSTEMATIC_VARIATIONS.keys()) + ["CMS_PDF"]
    for syst_name in all_shape_systs:
        cols = []
        for cat in categories:
            # Determine if category is Muon or Electron based
            is_muon = "Muon" in cat
            is_ele = "Electron" in cat or "EGamma" in cat
            for p in procs:
                if syst_name in CHANNEL_RESTRICTED_SYST:
                    allowed_ch = CHANNEL_RESTRICTED_SYST[syst_name]
                    if ("Muon" in allowed_ch and not is_muon) or ("Electron" in allowed_ch and not is_ele):
                        cols.append(f"{'-':<18}")
                        continue
                if syst_name in TTBAR_ONLY_SYST and not p.startswith("TT_"):
                    cols.append(f"{'-':<18}")
                    continue
                cols.append(f"{'1.0':<18}")
        lines.append(f"{syst_name:<24} shape " + "  ".join(cols))

    lines.append("-" * 100)
    # Barlow-Beeston-lite autoMCStats
    lines.append("* autoMCStats 5 0")

    lines.append("-" * 100)
    # Rate parameters for ttbb/ttbj and ttcc/ttcj
    for rp in RATE_PARAMETERS:
        rp_line = f"{rp['name']:<24} rateParam           {rp['channel_pattern']:<18}  {rp['process']:<18}  {rp['init']:<18.1f}  [{rp['min']:.0f},{rp['max']:.0f}]"
        lines.append(rp_line)

    output_dir = os.path.dirname(output_card_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_card_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    return output_card_path
