"""
Global Analysis Configuration for ATobb Combine Analysis.
Defines channels, regions (CR/SR), processes, systematics, and rate parameters.
"""

import os

# Base directory for raw analyzer ROOT inputs
BASE_DIR = "/data9/Users/bhoh/SKNanoOutput/AtobbInference/2024"
DEFAULT_OUTPUT_DIR = "/home/bhoh/ATobbCombine"
BASELINE_DIR = os.path.join(DEFAULT_OUTPUT_DIR, "outputs/baseline")
STUDIES_SCORE_CUTS_DIR = os.path.join(DEFAULT_OUTPUT_DIR, "outputs/studies_score_cuts")

# Data-taking eras and channels
PERIODS = ["C", "D", "E", "F", "G", "H", "I"]
CHANNELS = ["Muon", "Electron"]
MASS_POINTS = ["MA18", "MA25", "MA36", "MA50", "MA70", "MA90", "MA110", "MA125"]

# Luminosity and collision energy
LUMI_FB = 110.0
COM_TEV = 13.6
YEAR = "2024"

# Default region definitions
# Control Region (CR): bdt_score_avg <= CR_CUT_DEFAULT (background dominated)
# Signal Region (SR): bdt_score_avg > CR_CUT_DEFAULT (signal enriched)
CR_CUT_DEFAULT = 0.3
DEFAULT_REBIN = 2  # 100 bins -> 50 bins (bin width 0.02)
MIN_BIN_CONTENT = 1e-6

# Sample process definitions
TTBAR_FILES = ["TTLJ_powheg", "TTLL_powheg"]
TTBAR_CATEGORIES = ["bb", "bj", "cc", "cj", "jj"]

BACKGROUND_PROCESSES = ["TT_bb", "TT_bj", "TT_cc", "TT_cj", "TT_jj", "ST", "DY"]
ALL_PROCESSES = ["signal"] + BACKGROUND_PROCESSES

DATA_FILES_BY_CHANNEL = {
    "Muon": [f"Muon{stream}_{period}" for period in PERIODS for stream in (0, 1)],
    "Electron": [f"EGamma{stream}_{period}" for period in PERIODS for stream in (0, 1)],
}

DY_FILES_BY_CHANNEL = {
    "Muon": ["DYto2Mu_MLL50to120"],
    "Electron": ["DYto2E_MLL50to120"],
}

ST_FILES = [
    "ST_tW_antitop_Lep",
    "ST_tW_top_Lep",
    "ST_tW_top_Semilep",
    "ST_tch_antitop_lep",
    "ST_tW_antitop_Semilep",
    "ST_tch_top_lep",
]

# Systematic variations (Shape)
SYSTEMATIC_VARIATIONS = {
    "CMS_JESTotal": ("JESTotal_Up", "JESTotal_Down"),
    "CMS_JER": ("JER_Up", "JER_Down"),
    "CMS_BTagSF": ("BTagSF_Up", "BTagSF_Down"),
    "CMS_Pileup": ("Pileup_Up", "Pileup_Down"),
    "CMS_BFragmentation": ("BFragmentation_Up", "BFragmentation_Down"),
    "CMS_UnclusteredEnergy": ("UnclusteredEnergy_Up", "UnclusteredEnergy_Down"),
    "CMS_Top_Pt_Reweight": ("Top_Pt_Reweight_Up", "Top_Pt_Reweight_Down"),
    "CMS_ISR": ("ISR_Up", "ISR_Down"),
    "CMS_FSR": ("FSR_Up", "FSR_Down"),
    "CMS_hdamp": ("hdamp_Up", "hdamp_Down"),
    "CMS_PDFAlphaS": ("AlphaS_Up", "AlphaS_Down"),
    "CMS_LHEScale_MuF": ("LHEScale_MuF_Up", "LHEScale_MuF_Down"),
    "CMS_LHEScale_MuR": ("LHEScale_MuR_Up", "LHEScale_MuR_Down"),
    "CMS_LHEScale_MuFMuR": ("LHEScale_MuFMuR_Up", "LHEScale_MuFMuR_Down"),
    "CMS_eff_muon_id": ("Muon_ID_Up", "Muon_ID_Down"),
    "CMS_eff_muon_iso": ("Muon_Iso_Up", "Muon_Iso_Down"),
    "CMS_eff_muon_trig": ("Muon_Trig_Up", "Muon_Trig_Down"),
    "CMS_scale_m": ("MuonEn_Up", "MuonEn_Down"),
    "CMS_eff_e_id": ("Electron_ID_Up", "Electron_ID_Down"),
    "CMS_eff_e_reco": ("Electron_Reco_Up", "Electron_Reco_Down"),
    "CMS_eff_e_trig": ("Electron_Trig_Up", "Electron_Trig_Down"),
    "CMS_scale_e": ("ElectronEn_Up", "ElectronEn_Down"),
    "CMS_res_e": ("ElectronRes_Up", "ElectronRes_Down"),
}

# Channel-specific systematics
CHANNEL_RESTRICTED_SYST = {
    "CMS_eff_muon_id": ["Muon"],
    "CMS_eff_muon_iso": ["Muon"],
    "CMS_eff_muon_trig": ["Muon"],
    "CMS_scale_m": ["Muon"],
    "CMS_eff_e_id": ["Electron"],
    "CMS_eff_e_reco": ["Electron"],
    "CMS_eff_e_trig": ["Electron"],
    "CMS_scale_e": ["Electron"],
    "CMS_res_e": ["Electron"],
}

# Process-specific systematics
TTBAR_ONLY_SYST = ["CMS_Top_Pt_Reweight", "CMS_hdamp"]

# Normalization uncertainties (lnN)
LNN_UNCERTAINTIES = {
    "lumi_2024": {
        "value": 1.025,
        "processes": ALL_PROCESSES,
        "description": "2024 Luminosity uncertainty (2.5%)",
    },
    "ttbar_xsec": {
        "value": 1.05,
        "processes": ["TT_jj"],
        "description": "TT_jj normalization uncertainty (5.0%, ttbb/ttbj and ttcc/ttcj floated via rateParams)",
    },
    "ST_xsec": {
        "value": 1.15,
        "processes": ["ST"],
        "description": "Single top cross section uncertainty (15.0%)",
    },
    "DY_xsec": {
        "value": 1.20,
        "processes": ["DY"],
        "description": "Drell-Yan cross section uncertainty (20.0%)",
    },
    "ratio_ttbb": {
        "value": 1.50,
        "processes": ["TT_bb"],
        "description": "Relative ratio uncertainty between ttbb and ttbj (50.0%)",
    },
    "ratio_ttcc": {
        "value": 1.50,
        "processes": ["TT_cc"],
        "description": "Relative ratio uncertainty between ttcc and ttcj (50.0%)",
    },
}

# Floating rate parameters (rateParam)
RATE_PARAMETERS = [
    {
        "name": "ttbbXsec",
        "channel_pattern": "*",
        "process": "TT_bb",
        "init": 1.0,
        "min": 0.0,
        "max": 5.0,
    },
    {
        "name": "ttbbXsec",
        "channel_pattern": "*",
        "process": "TT_bj",
        "init": 1.0,
        "min": 0.0,
        "max": 5.0,
    },
    {
        "name": "ttccXsec",
        "channel_pattern": "*",
        "process": "TT_cc",
        "init": 1.0,
        "min": 0.0,
        "max": 5.0,
    },
    {
        "name": "ttccXsec",
        "channel_pattern": "*",
        "process": "TT_cj",
        "init": 1.0,
        "min": 0.0,
        "max": 5.0,
    },
]
