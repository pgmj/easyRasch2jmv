# easyRasch2jmv

A [Jamovi](https://www.jamovi.org/) module wrapping analysis logic from the
[easyRasch2](https://github.com/pgmj/easyRasch2) R package for Rasch
Measurement Theory analysis.

## Analyses

- **Item-Restscore Correlations** — Computes observed and model-expected
  item-restscore correlations (Goodman-Kruskal's gamma) via the `iarm`
  package. Supports dichotomous (Rasch model) and polytomous (Partial Credit
  Model) data. Reports observed vs. expected correlations, absolute
  differences, adjusted p-values, item locations, and item locations relative
  to the sample mean person location.

## Installation

1. Download the latest `.jmo` file from the [Releases](https://github.com/pgmj/easyRasch2jmv/releases) page.
2. In Jamovi, go to the **Modules** menu (⊞) → **Sideload**.
3. Select the downloaded `.jmo` file.

Alternatively, the module can be installed from the Jamovi library once published.

## Sample Data

Two sample datasets are bundled with the module:

- **pcmdat2** — Polytomous dataset (50 persons × 5 items, scored 0–3)
- **raschdat3** — Dichotomous dataset (50 persons × 8 items, scored 0–1)

These can be loaded from Jamovi's **Open** → **Data Library** after installing
the module.

## Dependencies

This module requires the following R packages (automatically handled by Jamovi):

- `eRm` — Rasch model fitting
- `iarm` — Item-restscore statistics
- `psychotools` (>= 0.7-3)
- `jmvcore` (>= 0.8.5)

## Author

Magnus Johansson <pgmj@pm.me>

## License

GPL (>= 3)
