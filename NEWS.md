# easyRasch2jmv 0.3.0 

New functions:

- Four reliability metrics: Person Separation Index (PSI) via
  eRm::SepRel(); empirical reliability via mirt::empirical_rxx(); Cronbach's 
  alpha (including bootstrap confidence intervals); and
  Relative Measurement Uncertainty (RMU) computed from mirt plausible
  values using the Bignardi, Kievit & Bürkner (2025) split-half
  correlation method.
- Transformation table from ordinal sum score to logit score, with either WLE or
  EAP scores.
- Non-parametric bootstrap of item-restscore associations,
  recommended for use with large sample sizes (n > 800). 
  See <https://pgmj.github.io/rasch_itemfit/> for more details.
- Conditional item infit with missing data, using multiple imputation. Since conditional 
  item infit needs complete data, incomplete responses are discarded. This function
  uses `mice` to impute n datasets with complete responses, then calculates infit from
  all datasets and pools the results using Rubin's rules. Optionally also uses
  simulation to determine appropriate cutoff for interpreting infit.
- New function for transformation of ordinal sum scores to interval logit scores.

Fixes/modifications:

- Partial gamma DIF bug fix: no longer allows numeric type variable as DIF variable.
- New data validation to better handle labelled ordinal data from SPSS, etc.
- Conditional item infit now clearly states it only uses complete response data,
  and also notes the complete response sample size in the table/figure footnote/caption.

# easyRasch2jmv 0.2.4

- Targeting (Wright map) figure added.

# easyRasch2jmv 0.2.3.2

- Updated description of module
- Fix for clearWith showLegend in item threshold probability figures.
- Added references

# easyRasch2jmv 0.2.3

- Now available from Jamovi Library
- Five functions implemented. Beyond the three listed previously:
  - Partial gamma DIF
  - Item threshold probability figures

# easyRasch2jmv 0.1.0

- Initial release with 3 functions
  - Item-restscore
  - Q3 residual correlations
  - Conditional item infit
