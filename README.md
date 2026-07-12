# easyRasch2jmv

<!-- badges: start -->
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<a href="https://doi.org/10.5281/zenodo.20136864" target="_blank"><img src="https://zenodo.org/badge/1199494308.svg" alt="DOI"></a>
<!-- badges: end -->

A [Jamovi](https://www.jamovi.org/) module for Rasch Measurement Theory
analysis. All analyses delegate their computations to the
[easyRasch2](https://github.com/pgmj/easyRasch2) R package, so results
are numerically identical to easyRasch2 with the same seeds and
iterations.

## Analyses

### Item fit

- **Conditional Item Infit** — conditional infit MSQ via the
  `easyRasch2` package (`iarm::out_infit()` on a CML fit with WLE person
  estimates; Müller, 2020), with optional simulation-based cutoffs
  (Johansson, 2025) and a dot-plot of observed vs simulated
  distributions. Numerically identical to `easyRasch2::RMitemInfit()` /
  `RMitemInfitCutoff()` with the same seed and iterations.
- **Conditional Item Infit (Multiple Imputation)** — pooled infit MSQ
  via Rubin's rules across `m` `mice` imputations (with optional
  auxiliary variables in the imputation model), pooled by
  `easyRasch2::RMitemInfitMI()`.
- **Bootstrap Item-Restscore** — non-parametric bootstrap of
  `iarm::item_restscore()` for use with large samples where the
  asymptotic test over-rejects. The bootstrap runs via
  `easyRasch2::RMitemRestscoreBoot()` (numerically identical draws and
  classifications with the same seed).
- **Item-Restscore Correlations** — observed vs model-expected
  item-restscore correlations via Goodman-Kruskal's gamma (Kreiner, 2011).
  Supports dichotomous and polytomous (Partial Credit Model) data.
  Delegates to `easyRasch2::RMitemRestscore()` (CML item locations via
  `psychotools`, WLE person-mean reference).

### Local dependence

- **Q3 Residual Correlation Matrix** — Yen's Q3 via the `easyRasch2`
  package (CML item estimation with WLE person estimates), with
  optional simulation-based cutoffs (Christensen et al., 2017;
  [Johansson, 2024](https://pgmj.github.io/simcutoffs.html)), a Q3
  heatmap, and a per-pair distribution plot. Numerically identical to
  `easyRasch2::RMlocdepQ3()` / `RMlocdepQ3Cutoff()` / `RMlocdepQ3Plot()`
  with the same seed and iterations.
- **Partial Gamma Local Dependence** — partial-gamma coefficients per
  item pair in both rest-score directions, with optional
  simulation-based per-pair expected ranges and a per-pair distribution
  plot; via `easyRasch2::RMlocdepGamma()` / `RMlocdepGammaCutoff()` /
  `RMlocdepGammaPlot()`.

### Dimensionality / unidimensionality

- **PCA of Standardized Residuals** — eigenvalues from PCA on the
  CML/WLE standardized Rasch residuals, with optional simulation-based
  cutoff for the first contrast (Chou & Wang, 2010) and a
  PC1-loading-vs-item-location plot; computed by
  `easyRasch2::RMdimResidualPCA()` / `RMdimResidualPCACutoff()`.
- **Dynamic/adaptive CFA Fit-Index Cutoffs** — observed one-factor
  categorical-CFA fit indices (CFI, RMSEA, SRMR) compared to a
  parametric-bootstrap null distribution simulated under the fitted
  PCM/RM, via `lavaan` WLSMV / ULSMV (avoids rule-of-thumb cutoffs),
  plus each item's standardized loading against its simulated expected
  range; computed by `easyRasch2::RMdimCFACutoff()` / `RMdimCFA()` with
  figures from `RMdimCFAPlot()`.

### Differential item functioning

- **Andersen LR-test DIF** — Andersen's likelihood-ratio test via
  `easyRasch2::RMdifLR()` (which remains based on `eRm::LRtest()`), with
  per-group item or threshold locations, MaxDiff flagging, a faceted
  figure of group-by-item locations, and a response-distribution
  tileplot.
- **Partial Gamma DIF** — partial-gamma coefficients for categorical
  DIF variables, with optional simulation-based expected ranges and a
  response-distribution tileplot; via `easyRasch2::RMdifGamma()` /
  `RMdifGammaCutoff()` / `RMplotTile()`.

### Reliability, targeting, score conversion

- **Reliability** — Cronbach's α, the WLE-based PSI, marginal
  reliability (Green, 1984), and RMU from plausible values (Bignardi,
  Kievit & Bürkner, 2025), computed by `easyRasch2::RMreliability()`;
  optional bootstrap HDCIs for α, PSI, and Marginal.
- **Targeting Plot** — Wright-map style person-item targeting with
  back-to-back histograms of person and item threshold locations, plus
  a threshold-location dot plot, drawn by `easyRasch2::RMtargeting()`
  (CML item parameters via `psychotools`, WLE person locations; MML
  fallback under sparse categories). Threshold table with Wald CIs from
  `easyRasch2::RMitemParameters()`.
- **Sum Score to Logit Transformation** — raw-score → person-location
  lookup, with WLE (CML via `psychotools`) or EAP (MML via `mirt`),
  computed by `easyRasch2::RMscoreSE()`.

### Visualization

- **Item Probability Curves** — model-implied category probability
  curves (polytomous, faceted per item) or joint item characteristic
  curves (dichotomous), computed by `easyRasch2::RMitemCatProb()`.
- **Item Characteristic Curves (CICC)** — expected item score curves
  with observed class-interval averages overlaid (Kreiner-style
  graphical item fit), via `easyRasch2::RMitemICCPlot()`; optional DIF
  variable with per-group observed averages and partial-gamma
  annotations.

## Requirements

All R-package dependencies are bundled as Imports and are installed
alongside the jamovi module: `easyRasch2`, `eRm`, `iarm`, `lavaan`,
`mirt`, `psychotools`, `mice`, `ggplot2`, `ggdist`, `patchwork`,
`scales`.

## Installation

In Jamovi, click on **Modules** (far top right) and choose **"Jamovi Library"**. Search for *easyRasch* and install from there.

For development version (see NEWS.md for changelog):

1. Download the latest `.jmo` file from the [Releases](https://github.com/pgmj/easyRasch2jmv/releases) page.
2. In Jamovi, go to the **Modules** menu (⊞) → **Sideload**.
3. Select the downloaded `.jmo` file.

## Sample Data

Two sample datasets are bundled with the module, both from R package `eRm`:

- **pcmdat2** — Polytomous dataset (50 persons × 5 items, scored 0–3)
- **raschdat3** — Dichotomous dataset (50 persons × 8 items, scored 0–1)

These can be loaded from Jamovi's **Open** → **Data Library** after installing
the module.


## References

- Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
  estimating reliability using Bayesian measurement uncertainty. PsyArXiv.
  <https://doi.org/10.31234/osf.io/h54k8>
- Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item-response
  models with principal component analysis on standardized residuals.
  *Educational and Psychological Measurement, 70*(5), 717–731.
  <https://doi.org/10.1177/0013164410379322>
- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for
  Yen's Q3: Identification of local dependence in the Rasch model using
  residual correlations. *Applied Psychological Measurement, 41*(3), 178–194.
  <https://doi.org/10.1177/0146621616677520>
- Johansson, M. (2024). Simulation-based cutoff values for Rasch item fit and
  residual correlations. <https://pgmj.github.io/simcutoffs.html>
- Johansson, M. (2025). Detecting item misfit in Rasch models.
  *Educational Methods & Psychometrics, 3*(18).
  <https://doi.org/10.61186/emp.2025.5>
- Kreiner, S. (2011). A note on item-restscore association in Rasch models.
  *Applied Psychological Measurement, 35*(7), 557–561.
  <https://doi.org/10.1177/0146621611410227>
- Mair, P., & Hatzinger, R. (2007). Extended Rasch modeling: The eRm package
  for the application of IRT models in R. *Journal of Statistical Software,
  20*(9). <https://doi.org/10.18637/jss.v020.i09>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust them?
  *Journal of Statistical Distributions and Applications, 7*(1), 5.
  <https://doi.org/10.1186/s40488-020-00108-7>
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
  *Journal of Statistical Software, 48*(2), 1–36.
  <https://doi.org/10.18637/jss.v048.i02>


## Credits

This is largely based on my `easyRasch` package, and I am using Claude Opus 
to "transfer" functions to a more properly formatted package - `easyRasch2` - 
which is the foundation for all code in the Jamovi module. 

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed 
psychologist with a PhD in behavior analysis. He works as a research specialist 
at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), 
Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
