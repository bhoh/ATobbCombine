# ATobbCombine: Statistical Inference Framework for $t \to H^{+}b,\, H^{+} \to W^{+}A\,(\to b\bar{b})$

This repository provides a modular, flexible statistical inference framework built on top of [CMS Higgs Combine](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) and [CombineHarvester](https://github.com/cms-analysis/CombineHarvester) for the $t \to H^{+}b,\, H^{+} \to W^{+}A\,(\to b\bar{b})$ search using 2024 Run 3 CMS data ($110\text{ fb}^{-1}$ at $\sqrt{s} = 13.6\text{ TeV}$).

---

## 1. Physics Model & POI Parameterization

The Parameter of Interest (POI) is the branching ratio:
$$\mathcal{B} \equiv \mathcal{B}(t \to H^{+}b)$$

Assuming $\mathcal{B}(t \to W b) = 1 - \mathcal{B}$, the process scalings in $pp \to t\bar{t}$ are parameterized as:

1. **Standard Model $t\bar{t}$ Background**:
   $$t\bar{t} \to W^{+}b\,W^{-}\bar{b} \quad \implies \quad \mu_{t\bar{t}} = (1 - \mathcal{B})^2$$

2. **Single Charged Higgs Signal**:
   $$t\bar{t} \to (H^{+}b)(W^{-}\bar{b}) + (W^{+}b)(H^{-}\bar{b}), \quad H^{+} \to W^{+}A\,(\to b\bar{b})$$
   $$\sigma_{\text{signal}}^{\text{eff}} = 2 \cdot \sigma_{tt} \cdot \left[ 2 \cdot (3 \cdot \mathcal{B}_{W\to\ell\nu}) \cdot \mathcal{B}_{W\to q\bar{q}} + (3 \cdot \mathcal{B}_{W\to\ell\nu})^2 \right] \cdot (1 - \mathcal{B}) \cdot \mathcal{B}$$
   where $\sigma_{tt} = 923.6\text{ pb}$, $\mathcal{B}_{W\to\ell\nu} = 0.1086$, $\mathcal{B}_{W\to q\bar{q}} = 0.6741$.
   $$\implies \text{signal\_scale} \approx 1007.293 \times (1 - \mathcal{B}) \cdot \mathcal{B}$$

This non-linear scaling is implemented in the custom Combine physics model [`modules/charged_higgs_br.py`](modules/charged_higgs_br.py).

---

## 2. Analysis Regions & Blinding Policy

- **Control Region (CR)**: $\text{bdt\_score\_avg} \le 0.3$ (background-dominated phase space).
  - **Real Data Fit Allowed**: Used to measure $t\bar{t}+bb$ and $t\bar{t}+cc$ rate parameters (`ttbbXsec`, `ttccXsec`) directly from data.
- **Signal Region (SR)**: $\text{bdt\_score\_avg} > 0.3$ (signal-enriched region for 8 mass points: $m_A \in \{18, 25, 36, 50, 70, 90, 110, 125\}\text{ GeV}$).
  - **Blinded Analysis**: Any measurement or fit involving the Signal Region (Combined CR+SR or SR-only) **strictly uses `-t -1` Asimov toys** to protect against unblinding bias.
- **Master 4-Category Workspaces with Channel Masks**:
  - Master workspaces compile all 4 joint categories: `Muon_CR`, `Electron_CR`, `Muon_SR`, `Electron_SR`.
  - Built with `--channel-masks` so that **Combined Fits**, **SR-Only Limits** (`mask_Muon_CR=1,mask_Electron_CR=1`), and **CR-Only Fits / GoF** (`mask_Muon_SR=1,mask_Electron_SR=1`) all reuse the single master workspace without redundant recreation.

---

## 3. Repository & Directory Structure

```
~/ATobbCombine/
├── configs/
│   ├── __init__.py
│   └── analysis_config.py      # Central physics config (channels, systematics, rateParams, paths)
│
├── modules/
│   ├── __init__.py
│   ├── charged_higgs_br.py     # Custom Combine PhysicsModel for BR(t -> H+ b)
│   ├── shape_extractor.py      # ROOT histogram extraction, slicing, & PDF envelope engine
│   ├── datacard_generator.py   # Datacard writer (4-channel simultaneous fit) with autoMCStats
│   └── workspace_builder.py    # text2workspace.py wrapper with --channel-masks support
│
├── scripts/
│   ├── produce_inputs.py       # [Main] Master workspace generator (shapes, datacards, workspaces)
│   ├── run_cr_fit.py           # CR background fit on real data & pull extraction (masked SR)
│   ├── run_gof.py              # Saturated Goodness-of-Fit test (masked SR)
│   ├── run_limits.py           # 95% CL limits on BR (Combined and SR-only via channel masks)
│   ├── run_impacts.py          # Nuisance parameter pulls & impact ranking plots
│   ├── plot_limits.py          # Brazilian limit plotter (Linear scale, decimal / percent)
│   ├── plot_postfit.py         # Pre-fit & Post-fit Data-to-MC distribution plotter
│   ├── scan_bdt_cuts.py        # Automated BDT boundary cut scan into isolated study workspaces
│   └── export_limits_to_csv.py # ROOT limit extractor to structured CSV
│
├── outputs/
│   ├── baseline/               # Master baseline analysis (CR <= 0.3, SR > 0.3)
│   │   ├── shapes/             # Input shape ROOT files
│   │   ├── datacards/          # Master 4-channel datacards
│   │   ├── workspaces/         # Master RooWorkspaces with channel masks
│   │   ├── cr_fit/             # CR-only real data fit results & pre/post-fit distribution plots
│   │   ├── gof/                # Saturated Goodness-of-Fit test results & plots
│   │   ├── impacts/            # Nuisance parameter impacts JSON, PDF, and PNG plots
│   │   ├── limits_combined/    # Combined CR+SR limits summary CSV & Brazilian plots (-t -1)
│   │   └── limits_sr_only/     # SR-only limits summary CSV & Brazilian plots (masked CR, -t -1)
│   │
│   └── studies_score_cuts/     # Isolated workspace directories per cut (0.15 ~ 0.40)
│       ├── cut_15/             # Isolated shapes, datacards, workspaces for cut <= 0.15
│       ├── cut_20/             # Isolated shapes, datacards, workspaces for cut <= 0.20
│       ├── ...
│       ├── bdt_cut_scan_comparison.csv
│       └── bdt_cut_scan_comparison.png
│
└── README.md                   # Analysis documentation and CLI guide
```

---

## 4. Environment Setup

Initialize CMSSW 16_0_0 and set up the Python environment:
```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ~/CMSSW_16_0_0/src && cmsenv
cd ~/ATobbCombine
export PYTHONPATH=/home/bhoh/ATobbCombine:$PYTHONPATH
chmod +x scripts/*.py
```

---

## 5. Workflow Guide

### Step 1: Generate Master Workspaces (Baseline: CR $\le 0.3$, SR $> 0.3$)
Compile the 4-channel master RooWorkspaces with `--channel-masks` and `ChargedHiggsBR`:
```bash
python3 scripts/produce_inputs.py
```

---

### Step 2: Control Region Fit & Post-fit Diagnostic Plots (Real Data)
Perform background profiling in the CR on real data (with masked SR) and inspect pre/post-fit Data/MC distributions:
```bash
# Run CR FitDiagnostics on master workspace with SR masked (Unblinded on CR real data)
python3 scripts/run_cr_fit.py --workspace outputs/baseline/workspaces/workspace_MA18.root --output-dir outputs/baseline/cr_fit

# Generate Pre-fit and Post-fit stacked distribution plots with Data/MC ratios
python3 scripts/plot_postfit.py --fit-root outputs/baseline/cr_fit/fitDiagnostics_CR.root --output-dir outputs/baseline/cr_fit/plots
```

---

### Step 3: Goodness-of-Fit (GoF) Test with Masked Signal Region
Perform the Saturated Goodness-of-Fit test on the background model by masking both Signal Region categories:
```bash
python3 scripts/run_gof.py --workspace outputs/baseline/workspaces/workspace_MA18.root --toys 100 --outdir outputs/baseline/gof
```

---

### Step 4: Compute $95\%\text{ CL}$ Upper Limits on $\mathcal{B}(t \to H^{+}b)$ (Blinded Asimov)
Calculate AsymptoticLimits for both Combined CR+SR and SR-only models:
```bash
# Compute limits across all 8 mass points (Combined and SR-only via channel masks)
python3 scripts/run_limits.py --mode both --ws-dir outputs/baseline/workspaces --outdir outputs/baseline

# Generate Brazilian limit plots (decimal format)
python3 scripts/plot_limits.py --csv outputs/baseline/limits_combined/limits_summary.csv --output outputs/baseline/limits_combined/limit_cls.png
python3 scripts/plot_limits.py --csv outputs/baseline/limits_sr_only/limits_summary.csv --output outputs/baseline/limits_sr_only/limit_cls.png

# Generate Brazilian limit plots (percentage [%] format)
python3 scripts/plot_limits.py --csv outputs/baseline/limits_combined/limits_summary.csv --output outputs/baseline/limits_combined/limit_cls_percent.png --percent
python3 scripts/plot_limits.py --csv outputs/baseline/limits_sr_only/limits_summary.csv --output outputs/baseline/limits_sr_only/limit_cls_percent.png --percent
```

---

### Step 5: Nuisance Parameter Impact Plots (CombineHarvester)
Evaluate nuisance parameter pulls and their ranking impacts on $\mathcal{B}(t \to H^+ b)$:
```bash
# Run Impacts for Combined CR+SR model (Strictly uses -t -1 Asimov toy)
python3 scripts/run_impacts.py --mass MA18 --mode combined --parallel 8 --outdir outputs/baseline/impacts

# Run Impacts for CR-only model (Real data fit)
python3 scripts/run_impacts.py --mass MA18 --mode cr_only --parallel 8 --outdir outputs/baseline/impacts
```
Outputs:
* `outputs/baseline/impacts/impacts_MA18_combined.json`
* `outputs/baseline/impacts/impacts_MA18_combined.pdf`
* `outputs/baseline/impacts/impacts_MA18_combined.png`

---

### Step 6: BDT Score Cut Optimization Scan Study
To evaluate the optimal CR/SR boundary threshold across a scan of cut values ($0.15, 0.20, 0.25, 0.30, 0.35, 0.40$), separate dedicated workspaces are generated:
```bash
python3 scripts/scan_bdt_cuts.py --cuts 0.15,0.20,0.25,0.30,0.35,0.40 --rebin 2
```

---

## 6. Expected Limits on $\mathcal{B}(t \to H^{+}b)$ (2024 CMS Data, $110\text{ fb}^{-1}$)

$95\%\text{ CL}$ Expected Upper Limits on $\mathcal{B}(t \to H^{+}b)$ including Barlow-Beeston-lite MC statistical uncertainties (`autoMCStats`):

| Mass Point ($m_A$) | Mass [GeV] | SR-Only Expected Median | **Combined CR+SR Expected Median** | **Combined Expected $\mathcal{B}\,\,[\%]$** | Sensitivity Gain |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **MA18** | 18 | $0.0033$ | **$0.0030$** | **$0.30\%$** | **$+9.1\%$** |
| **MA25** | 25 | $0.0018$ | **$0.0017$** | **$0.17\%$** | **$+5.6\%$** |
| **MA36** | 36 | $0.0015$ | **$0.0014$** | **$0.14\%$** | **$+6.7\%$** |
| **MA50** | 50 | $0.0010$ | **$0.0009$** | **$0.09\%$** | **$+10.0\%$** |
| **MA70** | 70 | $0.0006$ | **$0.0006$** | **$0.06\%$** | $-$ |
| **MA90** | 90 | $0.0004$ | **$0.0003$** | **$0.03\%$** | **$+25.0\%$** |
| **MA110** | 110 | $0.0013$ | **$0.0011$** | **$0.11\%$** | **$+15.4\%$** |
| **MA125** | 125 | $0.0061$ | **$0.0050$** | **$0.50\%$** | **$+18.0\%$** |
