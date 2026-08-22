"""
Shape Extraction Module for ATobb Combine Analysis.
Loads ROOT histograms, applies rebinning and score range slicing ([x_min, x_max]),
calculates PDF RMS envelopes, and writes sanitized TH1F shapes into output ROOT files.
"""

import math
import os
import uuid
import ROOT

from configs.analysis_config import (
    BASE_DIR,
    CHANNELS,
    DATA_FILES_BY_CHANNEL,
    DY_FILES_BY_CHANNEL,
    MIN_BIN_CONTENT,
    ST_FILES,
    SYSTEMATIC_VARIATIONS,
    CHANNEL_RESTRICTED_SYST,
    TTBAR_FILES,
    TTBAR_ONLY_SYST,
)

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def get_hist_from_file(file_path, hist_path):
    """Safely open a ROOT file and return a cloned detached TH1."""
    if not os.path.exists(file_path):
        return None
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        return None
    h = f.Get(hist_path)
    if not h:
        f.Close()
        return None
    h_clone = h.Clone(f"h_{uuid.uuid4().hex[:8]}")
    h_clone.SetDirectory(0)
    f.Close()
    return h_clone


def slice_and_sanitize_histogram(hist, x_min=0.0, x_max=1.0, is_data=False):
    """
    Slice a histogram to [x_min, x_max] and sanitize empty/negative bins for MC.
    Returns a new detached TH1F spanning the sub-range.
    """
    if not hist:
        return None

    x_min_val = 0.0 if x_min is None else max(0.0, float(x_min))
    x_max_val = 1.0 if x_max is None else min(1.0, float(x_max))

    # If full range [0, 1] requested, just sanitize in place
    if math.isclose(x_min_val, 0.0, abs_tol=1e-5) and math.isclose(x_max_val, 1.0, abs_tol=1e-5):
        nbins = hist.GetNbinsX()
        for b in range(1, nbins + 1):
            if not is_data and hist.GetBinContent(b) <= 0:
                hist.SetBinContent(b, MIN_BIN_CONTENT)
                hist.SetBinError(b, MIN_BIN_CONTENT)
        return hist

    # Determine bin indices
    bin_start = hist.FindBin(x_min_val + 1e-5)
    bin_end = hist.FindBin(x_max_val - 1e-5)
    bin_start = max(1, bin_start)
    bin_end = min(hist.GetNbinsX(), max(bin_start, bin_end))

    n_bins_new = bin_end - bin_start + 1
    new_x_low = hist.GetXaxis().GetBinLowEdge(bin_start)
    new_x_high = hist.GetXaxis().GetBinUpEdge(bin_end)

    h_sliced = ROOT.TH1F(f"h_slice_{uuid.uuid4().hex[:8]}", hist.GetTitle(), n_bins_new, new_x_low, new_x_high)
    h_sliced.SetDirectory(0)

    for new_b, old_b in enumerate(range(bin_start, bin_end + 1), start=1):
        content = hist.GetBinContent(old_b)
        error = hist.GetBinError(old_b)
        if not is_data and content <= 0:
            content = MIN_BIN_CONTENT
            error = MIN_BIN_CONTENT
        h_sliced.SetBinContent(new_b, content)
        h_sliced.SetBinError(new_b, error)

    return h_sliced


def load_process_hist(
    base_dir,
    channel,
    var_name,
    process_key,
    syst_suffix="Central",
    ma_point="MA18",
    rebin_factor=2,
    x_min=0.0,
    x_max=1.0,
):
    """
    Load, sum, rebin, and slice histograms for a given process and systematic variation.
    """
    suffix = "" if syst_suffix == "Central" else f"_{syst_suffix}"
    h_raw = None

    if process_key == "data_obs":
        for df in DATA_FILES_BY_CHANNEL[channel]:
            fpath = f"{base_dir}/{df}.root"
            hpath = f"{channel}/{var_name}/histo_Data"
            h = get_hist_from_file(fpath, hpath)
            if h:
                if h_raw is None:
                    h_raw = h.Clone(f"h_raw_{uuid.uuid4().hex[:8]}")
                    h_raw.SetDirectory(0)
                else:
                    h_raw.Add(h)

    elif process_key == "signal":
        sig_file = f"TTToHcToWAToBB-MHc130_{ma_point}_SingleLepFilter"
        fpath = f"{base_dir}/{sig_file}.root"
        hpath = f"{channel}/{var_name}/histo_{sig_file}{suffix}"
        h_raw = get_hist_from_file(fpath, hpath)
        if not h_raw and suffix != "":
            hpath_c = f"{channel}/{var_name}/histo_{sig_file}"
            h_raw = get_hist_from_file(fpath, hpath_c)

    elif process_key.startswith("TT_"):
        cat = process_key.replace("TT_", "")
        for tf in TTBAR_FILES:
            fpath = f"{base_dir}/{tf}.root"
            hpath = f"{channel}/{var_name}/histo_{tf}_{cat}{suffix}"
            h = get_hist_from_file(fpath, hpath)
            if not h and suffix != "":
                hpath_c = f"{channel}/{var_name}/histo_{tf}_{cat}"
                h = get_hist_from_file(fpath, hpath_c)
            if h:
                if h_raw is None:
                    h_raw = h.Clone(f"h_raw_{uuid.uuid4().hex[:8]}")
                    h_raw.SetDirectory(0)
                else:
                    h_raw.Add(h)

    elif process_key == "ST":
        for sf in ST_FILES:
            fpath = f"{base_dir}/{sf}.root"
            hpath = f"{channel}/{var_name}/histo_{sf}{suffix}"
            h = get_hist_from_file(fpath, hpath)
            if not h and suffix != "":
                hpath_c = f"{channel}/{var_name}/histo_{sf}"
                h = get_hist_from_file(fpath, hpath_c)
            if h:
                if h_raw is None:
                    h_raw = h.Clone(f"h_raw_{uuid.uuid4().hex[:8]}")
                    h_raw.SetDirectory(0)
                else:
                    h_raw.Add(h)

    elif process_key == "DY":
        for df in DY_FILES_BY_CHANNEL[channel]:
            fpath = f"{base_dir}/{df}.root"
            hpath = f"{channel}/{var_name}/histo_{df}{suffix}"
            h = get_hist_from_file(fpath, hpath)
            if not h and suffix != "":
                hpath_c = f"{channel}/{var_name}/histo_{df}"
                h = get_hist_from_file(fpath, hpath_c)
            if h:
                if h_raw is None:
                    h_raw = h.Clone(f"h_raw_{uuid.uuid4().hex[:8]}")
                    h_raw.SetDirectory(0)
                else:
                    h_raw.Add(h)

    if not h_raw:
        return None

    # Apply rebin factor
    if rebin_factor and rebin_factor > 1:
        h_raw.Rebin(rebin_factor)

    # Apply score range slicing and sanitization
    is_data = (process_key == "data_obs")
    return slice_and_sanitize_histogram(h_raw, x_min=x_min, x_max=x_max, is_data=is_data)


def calculate_pdf_envelope(
    base_dir,
    channel,
    var_name,
    process_key,
    h_nominal,
    ma_point="MA18",
    rebin_factor=2,
    x_min=0.0,
    x_max=1.0,
):
    """
    Compute Hessian PDF RMS variations (CMS_PDFUp and CMS_PDFDown) across 100 replicas.
    """
    if not h_nominal:
        return None, None
    nbins = h_nominal.GetNbinsX()
    pdf_hists = []

    for i in range(1, 101):
        var_suffix = f"PDF_{i:03d}"
        h_pdf = load_process_hist(
            base_dir,
            channel,
            var_name,
            process_key,
            var_suffix,
            ma_point=ma_point,
            rebin_factor=rebin_factor,
            x_min=x_min,
            x_max=x_max,
        )
        if h_pdf:
            pdf_hists.append(h_pdf)

    if not pdf_hists:
        return None, None

    h_up = h_nominal.Clone(f"h_pdf_up_{uuid.uuid4().hex[:8]}")
    h_down = h_nominal.Clone(f"h_pdf_down_{uuid.uuid4().hex[:8]}")
    h_up.SetDirectory(0)
    h_down.SetDirectory(0)

    n_replicas = len(pdf_hists)
    for b in range(1, nbins + 1):
        nom = h_nominal.GetBinContent(b)
        variance = sum((h.GetBinContent(b) - nom) ** 2 for h in pdf_hists) / max(1, (n_replicas - 1))
        sigma = math.sqrt(variance)
        h_up.SetBinContent(b, nom + sigma)
        h_down.SetBinContent(b, max(MIN_BIN_CONTENT, nom - sigma))

    return h_up, h_down


def extract_shapes_for_categories(
    categories,
    output_root_path,
    ma_point="MA18",
    rebin_factor=2,
):
    """
    Extract and write shape histograms for multiple analysis categories (e.g. Muon_CR, Electron_CR, Muon_SR, Electron_SR).
    categories: list of dicts with:
      - 'cat_name': e.g. 'Muon_CR' or 'Muon'
      - 'channel': 'Muon' or 'Electron'
      - 'var_name': 'bdt_score_avg' or f'bdt_score_{ma_point}'
      - 'x_min': float (e.g. 0.0 or 0.3)
      - 'x_max': float (e.g. 0.3 or 1.0)
    """
    out_f = ROOT.TFile.Open(output_root_path, "RECREATE")
    procs = ["signal", "TT_bb", "TT_bj", "TT_cc", "TT_cj", "TT_jj", "ST", "DY"]
    yields = {}

    for cat in categories:
        cat_name = cat["cat_name"]
        ch = cat["channel"]
        var_name = cat["var_name"]
        x_min = cat.get("x_min", 0.0)
        x_max = cat.get("x_max", 1.0)

        yields[cat_name] = {}

        # 1. Data observation
        h_data = load_process_hist(
            BASE_DIR,
            ch,
            var_name,
            "data_obs",
            "Central",
            ma_point=ma_point,
            rebin_factor=rebin_factor,
            x_min=x_min,
            x_max=x_max,
        )
        if h_data:
            out_f.cd()
            h_data.Write(f"histo_{cat_name}_data_obs", ROOT.TObject.kOverwrite)
            yields[cat_name]["data_obs"] = h_data.Integral()
        else:
            yields[cat_name]["data_obs"] = 0.0

        # 2. MC Processes
        for proc in procs:
            h_nom = load_process_hist(
                BASE_DIR,
                ch,
                var_name,
                proc,
                "Central",
                ma_point=ma_point,
                rebin_factor=rebin_factor,
                x_min=x_min,
                x_max=x_max,
            )
            if not h_nom:
                yields[cat_name][proc] = 0.0
                continue

            out_f.cd()
            h_nom.Write(f"histo_{cat_name}_{proc}", ROOT.TObject.kOverwrite)
            yields[cat_name][proc] = h_nom.Integral()

            # Shape systematics
            for syst_name, (up_var, down_var) in SYSTEMATIC_VARIATIONS.items():
                if syst_name in CHANNEL_RESTRICTED_SYST and ch not in CHANNEL_RESTRICTED_SYST[syst_name]:
                    continue
                if syst_name in TTBAR_ONLY_SYST and not proc.startswith("TT_"):
                    continue

                h_up = load_process_hist(
                    BASE_DIR,
                    ch,
                    var_name,
                    proc,
                    up_var,
                    ma_point=ma_point,
                    rebin_factor=rebin_factor,
                    x_min=x_min,
                    x_max=x_max,
                )
                h_down = load_process_hist(
                    BASE_DIR,
                    ch,
                    var_name,
                    proc,
                    down_var,
                    ma_point=ma_point,
                    rebin_factor=rebin_factor,
                    x_min=x_min,
                    x_max=x_max,
                )

                out_f.cd()
                if h_up:
                    h_up.Write(f"histo_{cat_name}_{proc}_{syst_name}Up", ROOT.TObject.kOverwrite)
                else:
                    h_nom.Write(f"histo_{cat_name}_{proc}_{syst_name}Up", ROOT.TObject.kOverwrite)

                if h_down:
                    h_down.Write(f"histo_{cat_name}_{proc}_{syst_name}Down", ROOT.TObject.kOverwrite)
                else:
                    h_nom.Write(f"histo_{cat_name}_{proc}_{syst_name}Down", ROOT.TObject.kOverwrite)

            # PDF RMS envelope
            h_pdf_up, h_pdf_down = calculate_pdf_envelope(
                BASE_DIR,
                ch,
                var_name,
                proc,
                h_nom,
                ma_point=ma_point,
                rebin_factor=rebin_factor,
                x_min=x_min,
                x_max=x_max,
            )
            if h_pdf_up and h_pdf_down:
                out_f.cd()
                h_pdf_up.Write(f"histo_{cat_name}_{proc}_CMS_PDFUp", ROOT.TObject.kOverwrite)
                h_pdf_down.Write(f"histo_{cat_name}_{proc}_CMS_PDFDown", ROOT.TObject.kOverwrite)

    out_f.Close()
    return yields
