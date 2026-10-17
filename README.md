# Vehicle-Mounted Radar Clutter Suppression — Simulation & Processing Code

MATLAB code for the manuscript:

> **A Vehicle-Mounted Radar Clutter Suppression Method on the Move Integrating IMU-Assisted Measurement and Dynamic Time-Frequency Ridge Tracking**, *Drones* (MDPI).

This repository contains the simulation, processing and evaluation code that produces the results, tables and figures of the paper (six-method comparison, module ablation, range-cell sensitivity and CNR sweep).

---

## 1. Requirements

- **MATLAB R2022b** (or newer). No additional toolboxes are strictly required; the code uses base MATLAB plus standard signal-processing functions that ship with MATLAB (e.g. `fft`, `awgn` from the Communications Toolbox — if `awgn` is unavailable, replace it with an equivalent complex-Gaussian noise addition).
- No GPU is needed; a desktop CPU is sufficient (see Table 8 of the paper for the per-CPI complexity).

---

## 2. Files to include in the repository

| File | Role | Needed? |
|---|---|---|
| **`ob_moco_v3_tf_ridge_EF.m`** | **Main script** — self-contained (all helper functions are local): builds the scene, runs the six methods (A MTI / B wide notch / C RD-soft notch / **D proposed** / E RPCA / F no-IMU), and produces all metrics, figures and the module-ablation, range-cell (N_c) and CNR arms. | **Required** |
| `run_single_experiment.m` | Single Monte-Carlo run of the four baseline schemes (A/B/C/D); returns the performance metrics. | Recommended |
| `run_expert_required_experiment.m` | Recorded-IMU replay with PCA and no-IMU baselines (supports §3.3 / Appendix A discussion). | Recommended |
| `run_expert_multitarget_sweep.m` | Multi-target SNR and target-speed sweep (supports the R2-10 CNR/velocity study). | Recommended |
| `ob_moco_v4_ablation_experiment.m`, `ob_moco_v4_ablation.m` | Module-ablation experiments (support Table 5). | Recommended |
| `run_v4_experiments.m` | Driver that runs the v4 experiment sets. | Recommended |
| `plot_expert_required_results.m` | Plots the metrics and RD comparisons. | Recommended |
| `prepare_real_imu.m` | Converts recorded IMU CSV files into the replay CSV used by `run_expert_required_experiment.m`. | Recommended |
| `ob_moco_v3_tf_ridge.m`, `ob_moco_v3_tf_ridge_1.m`, `ob_moco_v3_tf_ridge_2.m` | Earlier full versions of the pipeline (archived for reference). | Optional |
| `ob_moco_v2_dynamic_clutter.m` | Earlier clutter-model version (archived for reference). | Optional |
| `s_clutter_full.mat` (≈143 MB) | Pre-computed dynamic-clutter matrix. **Too large for a normal GitHub file (>100 MB).** Either omit it (the main script regenerates the clutter via `sim_dynamic_clutter` when the file is absent) or store it with **Git LFS**. | Optional |

> **Self-containedness:** `ob_moco_v3_tf_ridge_EF.m` defines every helper it needs (`my_go_cfar`, `compute_scr`, `tf_ridge_filter`, `viterbi_ridge`, `sim_dynamic_clutter`, `find_target_peak_prior`, `detect_candidate_targets`, …) as **local functions**, so no extra `.m` file is required to run it.

---

## 3. Quick start

```matlab
% add the repository to the path, then run the main script
addpath(genpath(pwd));
ob_moco_v3_tf_ridge_EF;
```

The script prints the metrics to the command window and writes, for every run, a timestamped folder

```
results_<yyyymmdd_HHMMSS>/
├── data/      % metrics_SCR.mat, metrics_detection.mat, metrics_clutter.mat,
│              % params_system.mat, params_target.mat, mid_frame_diagnostic.mat,
│              % metrics_summary.txt
├── figs/      % all figures (.png + .fig)
└── index.html % an index of the generated figures
```

### Reproducing the CNR sweep (Table 7)

The CNR sweep re-runs the pipeline at different additive-noise levels by writing the desired level to `snr_override.txt` before each run:

```matlab
for s = [-25 -20 -15 -10 -5]
    fid = fopen('snr_override.txt','w'); fprintf(fid,'%g',s); fclose(fid);
    ob_moco_v3_tf_ridge_EF;   % one run per noise level
end
```

Each run appends one line `SNR_db, SCR, Pd_slow, Pd_fast` to `snr_sweep.csv`.

> Note: in the paper the additive-noise level is denoted **N_L** (dB); in the code it is the `SNR_db` variable (the value is read from `snr_override.txt` when that file exists, otherwise the default `-15` dB is used).

---

## 4. Mapping to the paper

| Paper item | Produced by |
|---|---|
| Fig. 3 (time-frequency spectrum / mask) | `ob_moco_v3_tf_ridge_EF.m` (`fig1a`/`fig1b`) |
| Fig. 4 (RD maps for the six methods) | `ob_moco_v3_tf_ridge_EF.m` (`fig2a–fig2e`, RPCA/no-IMU panels) |
| Fig. 5 / Table 4 (SCR, FAR ratio, TERI) | `ob_moco_v3_tf_ridge_EF.m` (summary block) |
| Fig. 6 / Fig. 7 (Pd_slow/Pd_fast; Pd–Pfa, FAR, ROC) | `ob_moco_v3_tf_ridge_EF.m` (detection block) |
| Table 5 (module ablation) | `ob_moco_v3_tf_ridge_EF.m` (D vs D0/Dfx/Dfn arms) |
| Table 6 (N_c sensitivity) | `ob_moco_v3_tf_ridge_EF.m` (`SCR_Nc`, `Pd_slow_Nc`) |
| Table 7 (CNR sweep) | driver in §3 above (`snr_sweep.csv`) |
| Table 8 (computational cost) | complexity analysis in §4 of the paper (per-stage FLOPs) |
| Table 9 / §4.3 (IMU-robustness) | `run_expert_required_experiment.m` (IMU-error injection) |
| Appendix A (RSS / noise / Viterbi validation) | `run_expert_multitarget_sweep.m` and the derivation scripts |

All parameters used in the paper are listed in **Appendix B (Table A2)** of the manuscript.

---

## 5. Data

All data in the paper are **simulated** using the models described in §3 of the manuscript. The pre-computed clutter matrix `s_clutter_full.mat` is optional: when it is absent, the main script synthesises the dynamic clutter with `sim_dynamic_clutter` (matching the parameters of Table 2), so the code runs out of the box.

---

## 6. License & contact

- **License:** (to be specified — e.g. MIT License).
- **Contact:** the corresponding author of the paper (see the manuscript) for questions about the code or the simulation settings.

---

## 7. Citation

If you use this code, please cite the paper:

> *(full citation to be inserted upon publication)*