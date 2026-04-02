# easyRasch2jmv

A [Jamovi](https://www.jamovi.org/) module wrapping analysis logic from the
[easyRasch2](https://github.com/pgmj/easyRasch2) R package for Rasch
Measurement Theory analysis.

## Analyses

- **Conditional item infit MSQ** from R package `iarm` (Müller, 2020), with 
  simulation-based cutoff values (Johansson, 2025).

- **Yen's Q3 residuals** for assessing local dependencies, with simulation-based
  cutoff values (Christensen, et al., 2017; [Johansson, 2024](https://pgmj.github.io/simcutoffs.html)).

- **Item-Restscore Correlations** — Computes observed and model-expected
  item-restscore correlations (Goodman-Kruskal's gamma; Kreiner, 2011) via the `iarm`
  package. Supports dichotomous (Rasch model) and polytomous (Partial Credit
  Model) data. Reports observed vs. expected correlations, absolute
  differences, adjusted p-values, item locations, and item locations relative
  to the sample mean person location.

## Installation

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

- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical Values for Yen’s Q3: 
Identification of Local Dependence in the Rasch Model Using Residual Correlations. 
*Applied Psychological Measurement, 41*(3), 178–194. <https://doi.org/10.1177/0146621616677520>
- Johansson, M. (2025). Detecting Item Misfit in Rasch Models. 
*Educational Methods & Psychometrics, 3*(18). <https://doi.org/10.61186/emp.2025.5>
- Kreiner, S. (2011). A Note on Item–Restscore Association in Rasch Models. 
*Applied Psychological Measurement, 35*(7), 557–561. <https://doi.org/10.1177/0146621611410227>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust them? 
*Journal of Statistical Distributions and Applications, 7*(1), 5. 
<https://doi.org/10.1186/s40488-020-00108-7>


## Credits

As stated, this is based on my `easyRasch` package, and I am using Claude Opus 4.6 
to "transfer" functions to this more properly formatted package. While it uses my
earlier code, most of the code in this package is produced by the LLM and bug fixed by me.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed 
psychologist with a PhD in behavior analysis. He works as a research specialist 
at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), 
Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
