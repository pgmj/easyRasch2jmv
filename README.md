# easyRasch2jmv

<!-- badges: start -->
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<!-- badges: end -->

A [Jamovi](https://www.jamovi.org/) module wrapping analysis logic from the
[easyRasch2](https://github.com/pgmj/easyRasch2) R package for Rasch
Measurement Theory analysis.

## Analyses

### Item fit

- **Conditional Item Infit** — conditional infit MSQ via `iarm`
  (Müller, 2020), with optional simulation-based cutoffs (Johansson,
  2025) and a dot-plot of observed vs simulated distributions.
- **Conditional Item Infit (Multiple Imputation)** — pooled infit MSQ
  via Rubin's rules across `m` `mice` imputations.
- **Bootstrap Item-Restscore** — non-parametric bootstrap of
  `iarm::item_restscore()` for use with large samples where the
  asymptotic test over-rejects.
- **Item-Restscore Correlations** — observed vs model-expected
  item-restscore correlations via Goodman-Kruskal's gamma (Kreiner, 2011).
  Supports dichotomous and polytomous (Partial Credit Model) data. CML
  item locations from `eRm`.

### Local dependence

- **Q3 Residual Correlation Matrix** — Yen's Q3 via `mirt` (MML), with
  optional simulation-based cutoffs (Christensen et al., 2017;
  [Johansson, 2024](https://pgmj.github.io/simcutoffs.html)).

### Dimensionality / unidimensionality

- **PCA of Standardized Residuals** — eigenvalues from PCA on Rasch
  residuals, with optional simulation-based cutoff for the first
  contrast (Chou & Wang, 2010) and a PC1-loading-vs-item-location plot.
- **Dynamic/adaptive CFA Fit-Index Cutoffs** — observed one-factor
  categorical-CFA fit indices (CFI, RMSEA, SRMR) compared to a
  parametric-bootstrap null distribution simulated under the fitted
  PCM/RM, via `lavaan` WLSMV / ULSMV. Avoids rule-of-thumb cutoffs.

### Differential item functioning

- **Andersen LR-test DIF** — Andersen's likelihood-ratio test via
  `eRm::LRtest()`, with per-group item or threshold locations and a
  faceted figure of group-by-item locations.
- **Partial Gamma DIF** — partial-gamma coefficients via `iarm` for
  categorical DIF variables, with optional simulation-based cutoffs.

### Reliability, targeting, score conversion

- **Reliability** — Cronbach's α, PSI (`eRm::SepRel()`), empirical
  reliability (`mirt::empirical_rxx()`), and RMU from plausible values
  (Bignardi, Kievit & Bürkner, 2025).
- **Targeting Plot** — Wright-map style person-item targeting with
  back-to-back histograms of person and item threshold locations, plus
  a threshold-location dot plot. CML item parameters from `eRm`, MLE
  person parameters.
- **Sum Score to Logit Transformation** — raw-score → person-location
  lookup, with WLE (CML via `eRm`) or EAP (MML via `mirt`).

### Visualization

- **Item Probability Curves** — category probability curves for PCM
  items via `eRm::plotICC()`.

## Requirements

All R-package dependencies are bundled as Imports and are installed
alongside the jamovi module: `eRm`, `iarm`, `lavaan`, `mirt`,
`psychotools`, `mice`, `ggplot2`, `ggdist`, `patchwork`, `scales`.

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
