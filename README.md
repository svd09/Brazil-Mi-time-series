# Brazil Temperature Extremes & MI Mortality — Analysis Pipeline


## 1. What this project does

This is the analysis pipeline for a paper studying the association between
**extreme temperature events (heat and cold) and myocardial infarction (MI /
STEMI) mortality in Brazil**. The pipeline:

1. Defines heat and cold "extreme temperature episode" (ETE) indicators from
   daily mean/heat-index temperature, both nationally and per region (since
   temperature percentiles are location-specific).
2. Fits **distributed lag non-linear models (DLNM)** with a linear exposure–response
   and a spline lag–response, inside a **time-stratified case-crossover design**
   implemented as a **conditional Poisson model** (`gnm` with `eliminate = stratum`,
   quasi-Poisson family), adjusted for air pollutants (SO₂, CO, NO₂, O₃, PM2.5).
3. Extracts odds ratios (OR) for each exposure definition, for:
   - Brazil as a whole
   - Each region (North, Northeast, Southeast, South)
   - Demographic subgroups (male/female, young/old)
4. Computes **attributable fraction (AF)** and **attributable number (AN)** of
   deaths due to each exposure, using parametric simulation of the model
   coefficients (multivariate normal draws from `coef`/`vcov`), both:
   - over **all days**, and
   - restricted to **event days only** (to correct for the fact that
     multi-day/consecutive-day events are rarer, so raw totals understate
     their per-day impact).
5. Produces forest plots and subgroup comparison plots summarizing the OR/AF/AN
   results for the paper.

## 2. Exposure definitions

Extreme temperature episodes are derived from `tmean` (heat index/mean
temperature), computed **within each location** (`loc6digit`) using local
percentile thresholds:

| Code | Definition |
|---|---|
| `EHE_90` | Day ≥ 90th percentile of local `tmean` (moderate heat day) |
| `EHE_95` | Day ≥ 95th percentile (defined but not consistently modeled) |
| `EHE_99` | Day ≥ 99th percentile (extreme heat day) |
| `EHE_10` | Day ≤ 10th percentile (moderate cold day) |
| `EHE_5`  | Day ≤ 5th percentile (defined but not consistently modeled) |
| `EHE_01` | Day ≤ 1st percentile (extreme cold day) |

Each of `EHE_90`, `EHE_99`, `EHE_10`, `EHE_01` also has **`_2`** and **`_3`**
suffixed versions marking the day as an event only if the threshold was also
met on the **following 1 or 2 days** (i.e., 2‑day and 3‑day consecutive
events / "heatwave" or "coldwave" style indicators), e.g. `EHE_90_2`,
`EHE_90_3`.

> Despite the `EHE_*` (Extreme **H**eat Episode) naming convention, these same
> variables are used for cold thresholds too (`EHE_10`, `EHE_01`, ...). The
> plot scripts refer to the heat set as **EHE** and the cold set as **ECE**
> (Extreme **C**old Episode) — the naming is inconsistent between scripts,
> not a difference in the underlying variables.

## 3. Repository structure (inferred)

The scripts assume the following folders exist relative to the working
directory (paths are used as-is in the code, both relative and — in a couple
of places — hardcoded absolute HPC paths):

```
data/                # .qs / .rds analytic datasets (whole country + per region)
model objects/        # saved cb / model / crosspred objects (.rds / .qs)
results/              # csv outputs: OR tables, AF/AN tables, forest-plot tables
graphs/                # .tiff figures (forest plots, subgroup plots)
codes/                 # shared function scripts, sourced by the paper scripts
```

## 4. Script inventory & suggested run order

The scripts are not currently numbered, and several read data/objects saved
by earlier scripts. Based on `source()`, `read*()`, and `write*()` calls, the
inferred dependency order is:

| Order | Script | Purpose |
|---|---|---|
| 0 | `AN_ON_EVENT_DAYS_FUNCTION.R` | Defines `af_an_from_beta_event_days()` — the core AF/AN function (see §5). Save this under `codes/` so other scripts can `source()` it. |
| 1 | `ETE_AND_ECE_MODELS_FOR_WHOLE_BRAZIL_.R` | Fits national heat (`EHE_90...`) and cold (`EHE_10...`) DLNM/conditional-Poisson models on the whole dataset; writes `results/brazil_heat_res.csv` and `results/brazil_cold_res.csv`. |
| 2 | `PAPER_SCRIPT_MODEL_ETE_AND_AF_AND_AN_.R` | Re-fits whole-Brazil models for `EHE_01`, `EHE_10` (and their `_2`/`_3` and `_90`/`_99` counterparts), extracts cumulative OR by lag, and computes AF/AN via coefficient simulation. Saves model objects (`model objects/results_ehe_*.rds`). |
| 3 | `PAPER_SCRIPT_GETTING_AF___AN_PER_EVENT_DAY.R` | Revised AF/AN calculation restricted to event days only, using `af_an_from_beta_event_days()`. Produces `results/ehe_01_af_an_event_days_per_day.csv`, `results/ehe_10_af_an_event_days_per_day.csv`, etc. |
| 4 | `PAPER_SCRIPT_MODEL_ETE_FOR_REGIONS.R` | Builds region-specific cold exposure variables (percentiles computed **within region**) for North, Northeast, Southeast, South; fits cold models per region; saves `model objects/<region>_cold_results.{rds,qs}` and `results/<region>_cold_res.csv`. |
| 5 | `PAPER_SCRIPT_ETE_MODELS_HEAT_REGIONS.R` | Same as above but for heat exposure variables/models, per region; saves `model objects/<region>_heat_results.qs` and `results/<region>_heat_results.csv`. |
| 6 | `PAPER_SCRIPT_AF_AND_AN_FOR_REGIONAL_MODELS.R` | Loads the saved regional model objects (currently wired up for **North** and **Southeast** only — see §7) and computes AF/AN per exposure variable per region; writes `results/<region>_cold_af_an_results.csv`. |
| 7 | `PAPER_SCRIPT_MODELS_BY_AGE_AND_SEX_HPC_PART_III.R` | Subcohort analysis: splits data into males/females/young/old, builds exposure variables, fits heat and cold models for each subgroup, writes `<group>_or.csv` / `<group>_res_cold.csv`. Designed to run on the HPC cluster (absolute paths, `qs`/`dtplyr` for memory efficiency). |
| 8 | `MAKING_FORESTPLOTS.R` | Reads the region-level result tables (`results/ehe_by_location_grouped.csv`, `results/ehe_estimates_long.csv`, `results/ece_rr_table.csv`) and builds forest plots with `ckbplotr::forest_plot()`, combined with `cowplot`. Outputs `graphs/forest_plot_ehe.tiff` and `graphs/forest_plot_ece.tiff`. |
| 9 | `EHE_PLOT_SUBCOHORT_HEAT.R` / `ECE_PLOT_SUBCOHORT_COLD.R` | Base-R plots comparing OR (with 95% CI) across the four subcohorts (Males/Females/Young/Old) for the heat and cold event definitions respectively. Output: `graphs/EHE_subgroup_plot.tiff`, `graphs/ECE_effects_by_subgroup.tiff`. **Currently use inline/hardcoded data** rather than reading from `results/` — see §7. |

## 5. Core function: `af_an_from_beta_event_days()`

Defined in `AN_ON_EVENT_DAYS_FUNCTION.R`.

**Purpose:** given a simulated lag-coefficient vector, compute the daily
attributable number (AN) of deaths via the lag convolution of exposure `x`
with the coefficients, then summarize AF/AN both over all days and restricted
to event days.

**Signature:**
```r
af_an_from_beta_event_days(
  beta_l,            # numeric vector of lag coefficients (length L+1)
  x,                 # exposure vector (e.g., EHE_10)
  y,                 # outcome vector (daily event/death counts)
  trunc_excess = TRUE,   # TRUE: zero out "protective" (RR<1) days; FALSE: allow negative AN
  use_mu = FALSE,        # use fitted means instead of raw y as the denominator for AF
  mu = NULL,
  event_mask = NULL      # logical vector marking "event days"; if NULL and x is 0/1, uses x==1
)
```

**Returns** a one-row `data.frame` with:
`AF_all`, `AN_tot_all`, `n_days_all`, `n_event_days`, `AN_sum_event_days`,
`AN_per_event_day`, `AF_on_event_days`.

**Used by:** `PAPER_SCRIPT_GETTING_AF___AN_PER_EVENT_DAY.R`, typically called
inside a `pblapply()` loop over Monte Carlo draws of the coefficient vector
(`MASS::mvrnorm(nsim, mu = b, Sigma = V)`), then row-bound to get a simulation-based
distribution of AF/AN.

## 6. Dependencies

Installed as a single call via `easypackages::libraries(...)` at the top of
most scripts:

```
tidyverse, tidylog, vroom, lubridate, stringr, purrr, tidyr, sf, geobr,
readxl, ggrepel, weathermetrics, zoo, survival, dlnm, splines, mvmeta,
metafor, gnm, broom, ggh4x, rlist, mixmeta, ckbplotr, janitor,
bannerCommenter, qs, data.table, pbapply, MASS, Matrix, dtplyr, cowplot,
progress
```

Key modeling packages:
- **`dlnm`** — crossbasis / crosspred / crossreduce for distributed-lag non-linear models
- **`gnm`** — conditional (quasi-)Poisson regression with `eliminate=` for the case-crossover stratum
- **`MASS`** — `mvrnorm()` for coefficient simulation (AF/AN uncertainty)
- **`ckbplotr`** + **`cowplot`** — forest plots
- **`qs`** — fast serialization of large model objects

## 7. Known issues / cleanup suggestions

These are observations from reading the code, not changes that have been made:

1. **Missing function file:** `af_an_from_beta()` is called but never defined
   in the uploaded scripts (see §5). Locate `codes/AF and AN for conditional
   poisson models function.R` and add it to the repo, or confirm it's just
   `af_an_from_beta_event_days()` under an older name.
2. **Inconsistent I/O paths:** mix of relative (`results/...`) and absolute,
   user-specific paths (`/home/svd14/...`, `/mnt/vstor/...`). Recommend
   centralizing via a project-root helper (e.g. `here::here()`) or a `config.R`.
3. **Regional AF/AN script incomplete:** `PAPER_SCRIPT_AF_AND_AN_FOR_REGIONAL_MODELS.R`
   only has working blocks for **North** and **Southeast**; Northeast and South
   are not yet wired up (their model objects are produced by
   `PAPER_SCRIPT_MODEL_ETE_FOR_REGIONS.R` / `PAPER_SCRIPT_ETE_MODELS_HEAT_REGIONS.R`
   but never consumed here).
4. **`PAPER_SCRIPT_GETTING_AF___AN_PER_EVENT_DAY.R`** has leftover/duplicated
   code near the end (a second, unused `EHE_01`/`model_01_2` block, and a
   `rm(list = setdiff(ls(), c("results_01_2", ...)))` referencing objects that
   are never created — `models_ehe`, `cb_ehe`). Worth trimming before final use.
5. **`PAPER_SCRIPT_MODELS_BY_AGE_AND_SEX_HPC_PART_III.R`** has a large block of
   "junk code" explicitly marked as such at the bottom (duplicate crossbasis/
   model loops, a `progress_bar` demo snippet) — flagged in-code by the
   original author as safe to delete.
6. **Plot scripts use hardcoded data, not pipeline output:**
   `EHE_PLOT_SUBCOHORT_HEAT.R` and `ECE_PLOT_SUBCOHORT_COLD.R` currently build
   their `df` from literal numbers typed into the script rather than reading
   from the `results/*_or.csv` / `results/*_res_cold.csv` files written by
   `PAPER_SCRIPT_MODELS_BY_AGE_AND_SEX_HPC_PART_III.R`. If those numbers were
   copy-pasted from that script's console output, consider having the plot
   scripts read the CSVs directly so the figures update automatically when
   models are re-run.
7. **Naming/typo:** in `ECE_PLOT_SUBCOHORT_COLD.R`, one event level is
   `"ECE_02_3"` — almost certainly meant to be `"ECE_01_3"` (i.e. `EHE_01_3`),
   consistent with the `EHE_01`/`EHE_01_2`/`EHE_01_3` naming used everywhere else.
8. **`gnm(..., eliminate = factor(stratum))` re-fit repeatedly per exposure
   variable** — each region/subgroup script loops over 6–7 exposure variables,
   refitting a full model each time. This is fine for correctness but is the
   main runtime/memory cost of the pipeline (hence the aggressive `rm()`/`gc()`
   calls and HPC use for the subgroup script).

## 8. Output summary (what each script ultimately produces)

| Output | Produced by | Description |
|---|---|---|
| `results/brazil_heat_res.csv`, `results/brazil_cold_res.csv` | `ETE_AND_ECE_MODELS_FOR_WHOLE_BRAZIL_.R` | National OR per exposure definition |
| `results/ehe_01_overall.csv`, `results/ehe_10_overall.csv` | `PAPER_SCRIPT_MODEL_ETE_AND_AF_AND_AN_.R` | Cumulative OR by lag day (0–7) |
| `model objects/results_ehe_*.rds` | `PAPER_SCRIPT_MODEL_ETE_AND_AF_AND_AN_.R` | Saved `cb`/`model`/`crosspred` objects per exposure |
| `results/ehe_*_af_an_event_days_per_day.csv` | `PAPER_SCRIPT_GETTING_AF___AN_PER_EVENT_DAY.R` | Simulation-based AF/AN, all-days & event-days-only |
| `results/<region>_cold_res.csv`, `results/<region>_heat_results.csv` | `PAPER_SCRIPT_MODEL_ETE_FOR_REGIONS.R`, `PAPER_SCRIPT_ETE_MODELS_HEAT_REGIONS.R` | Regional OR per exposure |
| `model objects/<region>_cold_results.*`, `model objects/<region>_heat_results.qs` | same | Regional model objects for downstream AF/AN |
| `results/<region>_cold_af_an_results.csv` | `PAPER_SCRIPT_AF_AND_AN_FOR_REGIONAL_MODELS.R` | Regional AF/AN (North, Southeast only currently) |
| `results/<group>_or.csv`, `results/<group>_res_cold.csv` | `PAPER_SCRIPT_MODELS_BY_AGE_AND_SEX_HPC_PART_III.R` | Subgroup (sex/age) OR for heat & cold |
| `graphs/forest_plot_ehe.tiff`, `graphs/forest_plot_ece.tiff` | `MAKING_FORESTPLOTS.R` | Forest plots of regional OR, heat & cold |
| `graphs/EHE_subgroup_plot.tiff`, `graphs/ECE_effects_by_subgroup.tiff` | `EHE_PLOT_SUBCOHORT_HEAT.R`, `ECE_PLOT_SUBCOHORT_COLD.R` | OR comparison across Males/Females/Young/Old |




