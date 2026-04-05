# Pairwise U-Statistics for Heavy-Tailed Inference

This repository contains code for simulation and numerical illustration of pairwise U-statistic methods for variance and fourth-order inference under heavy-tailed sampling.

The code supports:
- Pairwise variance estimation via `SPDV_n`
- Pairwise quartic dispersion estimation via `SPFD_n`
- Confidence intervals for `sigma^2`, `theta4`, and pairwise `mu4`
- Normal-theory, percentile bootstrap, and studentized bootstrap procedures
- Numerical illustrations of sharp moment thresholds and projection dominance

## Overview

This repository accompanies a paper on projection-driven sharp moment thresholds for pairwise U-statistics. The main theoretical message is:

- Variance inference based on the quadratic pairwise statistic has a sharp fourth-moment threshold for asymptotic normality.
- Fourth-order inference based on the quartic pairwise statistic has a sharp eighth-moment threshold for asymptotic normality.
- Studentized bootstrap calibration for the quartic statistic requires stronger moment conditions, with a twelfth-moment threshold appearing in heavy-tailed settings.

The code is organized so that simulation scripts call a central helper file, `functions_pairwise.R`, which implements estimators, bootstrap routines, truth functions, and Monte Carlo wrappers.

## Repository structure

```text
.
├── functions_pairwise.R
├── sim_sigma2_grid.R
├── sim_theta4_grid.R
├── fig_coverage_vs_df.R
├── fig_coverage_vs_df_theta4.R
├── fig_projection_dominance.R
├── results/
└── README.md
```

You may also add other scripts for supplementary tables, additional figures, or exploratory simulations.

## Main files

### `functions_pairwise.R`

Core implementation file. It contains:

- Pairwise estimators:
  - `spdv()`
  - `spfd()`
  - `mu4pairwise()`

- Projection and variance-related helpers:
  - `projectionspdv()`
  - `projectionspfdempirical()`
  - `projectionspfdplugin()`
  - `sespdvasymptotic()`
  - `sespdvprojection()`
  - `sespfdprojection()`

- Population truth functions:
  - `truthsigma2()`
  - `truththeta4()`
  - `truthmu4()`

- Data generators:
  - `simulateonesample()`
  - `rstdt()`
  - `rnorm1()`
  - `rcontamnorm()`

- Confidence interval routines:
  - `cisigma2normal()`
  - `cisigma2percentile()`
  - `cisigma2studentized()`
  - `citheta4normal()`
  - `citheta4percentile()`
  - `citheta4studentized()`
  - `cimu4pairwisepercentile()`

- Simulation wrappers:
  - `runonerepsigma2()`
  - `runonereptheta4()`
  - `runonedesign()`

### `sim_sigma2_grid.R`

Main simulation driver for variance inference. It runs Monte Carlo experiments for confidence intervals targeting `sigma^2` using the pairwise variance statistic `SPDV_n` across a grid of sample sizes and distributions.

Typical outputs include:
- empirical coverage for normal, percentile, and studentized intervals
- average interval widths
- failure counts or instability diagnostics
- summary tables across sample sizes and distribution families

### `sim_theta4_grid.R`

Main simulation driver for fourth-order inference. It runs Monte Carlo experiments for confidence intervals targeting `theta4` using the quartic pairwise statistic `SPFD_n` across a grid of sample sizes and distributions.

Typical outputs include:
- empirical coverage for `theta4` intervals
- comparisons among normal-theory, percentile bootstrap, and studentized bootstrap procedures
- behavior across heavy-tail regimes
- numerical evidence for the eighth- and twelfth-moment threshold phenomena

### `fig_coverage_vs_df.R`

Builds a coverage-versus-degrees-of-freedom figure for variance inference, usually at fixed sample size.

### `fig_coverage_vs_df_theta4.R`

Builds a coverage-versus-degrees-of-freedom figure for fourth-order inference, highlighting the sharper threshold behavior for the quartic statistic.

### `fig_projection_dominance.R`

Builds the projection-dominance diagnostic for the quartic statistic. The main quantity is

\[
R_n = \frac{\operatorname{Var}(\text{linear projection term})}{\operatorname{Var}(\text{full centered statistic})}.
\]

Values near 1 indicate that the linear Hoeffding projection dominates the finite-sample variability.

## Requirements

This project uses R.

### Recommended R version

- R 4.2 or newer

### Required packages

- `ggplot2`

Install the main package dependency with:

```r
install.packages("ggplot2")
```

If you plan to extend plotting or table generation, you may also want:
- `dplyr`
- `tidyr`
- `readr`

though they are not required by the core scripts listed above.

## Supported distributions

The helper code supports the following data-generating distributions:

- `"normal"`
- `"contam"` or `"contamnorm"`
- `"t<df>"`, such as:
  - `"t3"`
  - `"t4"`
  - `"t8"`
  - `"t12"`
  - `"t20"`

For Student-\(t\) designs, the generator uses a standardized \(t\) variable with variance 1 whenever that scaling is well-defined.

## Running the scripts

From an R session:

```r
source("functions_pairwise.R")
source("sim_sigma2_grid.R")
source("sim_theta4_grid.R")
source("fig_coverage_vs_df.R")
source("fig_coverage_vs_df_theta4.R")
source("fig_projection_dominance.R")
```

Or from the shell:

```bash
Rscript sim_sigma2_grid.R
Rscript sim_theta4_grid.R
Rscript fig_coverage_vs_df.R
Rscript fig_coverage_vs_df_theta4.R
Rscript fig_projection_dominance.R
```

## Output

Most scripts save outputs under a `results/` directory. Typical outputs include:

- `.csv` files with raw Monte Carlo results
- `.csv` files with summarized coverage or variance ratios
- `.pdf` figures
- `.png` figures

Example:

```text
results/
├── fig_coverage_vs_df_theta4_raw.csv
├── fig_coverage_vs_df_theta4_summary.csv
├── fig_coverage_vs_df_theta4.pdf
├── fig_coverage_vs_df_theta4.png
├── fig_projection_dominance_raw.csv
├── fig_projection_dominance_summary.csv
├── fig_projection_dominance.pdf
└── fig_projection_dominance.png
```

Depending on your simulation scripts, the `results/` folder may also contain outputs from:
- `sim_sigma2_grid.R`
- `sim_theta4_grid.R`

## Reproducibility

Each script sets a random seed at the top. For publication-quality runs, increase the Monte Carlo and bootstrap replication counts in the scripts.

Typical development settings:
- `M <- 1000` or `2000`
- `B <- 499`

Typical final settings:
- `M <- 5000`
- `B <- 1999`

These larger runs may take substantial time depending on sample size and bootstrap settings.

## Naming conventions

The helper file uses base-style function names without snake_case conversion. For example:

- `simulateonesample()` instead of `simulate_one_sample()`
- `truththeta4()` instead of `truth_theta4()`
- `projectionspfdempirical()` instead of `projection_spfd_empirical()`

If you write new scripts, make sure you match the exact function names defined in `functions_pairwise.R`.

## Example workflow

A typical workflow is:

1. Source `functions_pairwise.R`.
2. Run `sim_sigma2_grid.R` or `sim_theta4_grid.R` for the main Monte Carlo study.
3. Generate figure scripts such as `fig_coverage_vs_df.R`, `fig_coverage_vs_df_theta4.R`, or `fig_projection_dominance.R`.
4. Save numerical summaries and figures to `results/`.
5. Use the generated outputs in the manuscript or supplement.

## Minimal example

```r
source("functions_pairwise.R")

set.seed(1)

x <- simulateonesample(n = 100, distname = "t8")

spdv(x)
spfd(x)
mu4pairwise(x)

citheta4normal(x)
citheta4percentile(x, B = 499)
citheta4studentized(x, B = 499)
```

## Notes

- The quartic statistic is more sensitive to heavy tails than the variance statistic.
- Studentized bootstrap intervals can become unstable in regimes where asymptotic normality still appears reasonable.
- The projection-dominance diagnostic is especially useful for understanding when the linear term adequately explains finite-sample behavior.

## Suggested citation

If you use this repository, please cite the associated paper or manuscript.

Example BibTeX entry:

```bibtex
@article{srivastav_pairwise_heavytails,
  title   = {Projection-Driven Sharp Moment Thresholds for Pairwise U-Statistics},
  author  = {Srivastav, Sudesh K. and Srivastav, Apurv},
  journal = {Working paper / manuscript},
  year    = {2026}
}
```

Update the entry with the final journal information once available.

## License

Add your preferred license here, for example:

- MIT License
- GPL-3.0
- CC BY 4.0 for documentation

Example placeholder:

```text
This project is licensed under the MIT License - see the LICENSE file for details.
```

## Contact

For questions, issues, or collaboration, please open an issue on GitHub or contact the repository author.