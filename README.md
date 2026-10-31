# Monsaraz em Mudança — Superfície de Água do Alqueva

A reproducible geospatial pipeline and web application that tracks the change in the
**Alqueva reservoir's water surface within the parish (*freguesia*) of Monsaraz**, Portugal,
from Copernicus Sentinel-2 imagery (2017–2025).

**Live demo:** https://gustavosanta-star.github.io/monsaraz-agua-alqueva/

![The Monsaraz em Mudança web application](app.png)

---

## Overview

The application lets anyone compare the **Q3 (July–September) water surface** of the parish between
any two years: the delineated water extent on a satellite basemap, the area of each date, the
absolute and relative difference, a quarterly time series, and a documented methodology.

The engineering goal was a pipeline where **adding a new date is a single function call** — not a
manual GIS session — with reproducible analysis, documented indicators and quality control.

Between the lowest and highest Q3 values in the series, the parish's water surface grew by
~276 ha (+27%); the change map below shows **where** it grew (green = water gained, 2020 → 2025):

![Water gained between 2020 and 2025 (green)](mudanca_2020_2025.png)

## Study area

Freguesia de **Monsaraz** (municipality of Reguengos de Monsaraz, Évora district), on the right
bank of the Guadiana / Alqueva reservoir. Official boundary from **CAOP (DGT)**. Area ≈ 8,829 ha.
The large body of the Alqueva lies mostly outside the parish; the indicator measures only the water
**inside the administrative boundary**.

## Data & method

| Step | Choice |
|------|--------|
| Imagery | Copernicus **Sentinel-2 L2A** (10 m), via Google Earth Engine |
| Temporal composite | Per-quarter **median**, with cloud/shadow masking (Sentinel-2 **SCL** band) |
| Water detection | **MNDWI** = (B3 − B11) / (B3 + B11), threshold **0** |
| Post-processing | Vectorization + 10 m simplification |
| Areas | Computed in **ETRS89/PT-TM06** (EPSG:3763); served in WGS84 |

## Consistency check — JRC Global Surface Water

Comparison with **JRC Global Surface Water** (Pekel et al., 2016), the EU reference dataset for
surface water. JRC GSW is itself derived from **optical imagery (Landsat)**, so agreement shows
**cross-sensor consistency, not independent validation**: two optical methods can share the same
errors. An independent test (area vs. reservoir water level, EDIA/SNIRH) is planned.

- **Spatial agreement (IoU): 92.7%** — July 2020, compared only over JRC-observed pixels, at 30 m.
- **Temporal correlation (Pearson r): 0.80** — 41 months (2017–2021) with >80% JRC coverage, areas
  compared over the same observed pixels, at 30 m.

## Uncertainty

Two terms are reported. Neither is a statistical confidence interval.

1. **Threshold sensitivity:** area change for a ±0.05 shift of the MNDWI threshold — mean
   ≈ ±15 ha (~1.3%).
2. **Shoreline (mixed-pixel) term:** MNDWI uses B11, native **20 m**. Shoreline length inside the
   parish, counted on the 20 m grid, is **46–54 km**. Worst case — every edge pixel biased by half
   a pixel in the same direction — gives **±46–54 ha (~4%)**; with independent errors it is
   **<1 ha**. The threshold term falls between the two.

A systematic edge bias affects both dates in the same direction and largely cancels in a
difference. Even under the most conservative assumption (opposite biases, added linearly), the
2020 → 2025 change (+276 ha) exceeds the combined ±100 ha; small differences such as
2022 → 2023 (+35 ha) do not, and the app flags them as within the uncertainty.

## Quality checks

- **Non-reservoir water** (farm ponds, small dams): isolated using a reservoir envelope built from
  the JRC GSW maximum water extent (largest connected body, ≥ ~92 ha, dilated 30 m). It accounts for
  **0.2–1.0 ha in Q3 (<0.1%)** — negligible, so it is reported rather than removed. Many ponds are
  not detected at all in Q3 (dry, turbid or vegetated).
- **Cloud gaps:** ~0.1% of the parish without valid data in every quarterly composite.

## Limitations

- **Fixed MNDWI threshold (0)** — not scene-adaptive (turbidity, sun angle).
- The **quarterly median smooths** within-quarter water-level changes: each value is a *typical*
  extent, not a single instant.
- The dominant area uncertainty is at the **shoreline**: MNDWI uses B11 (native **20 m**), so mixed
  edge pixels along this very irregular shoreline may weigh as much as the threshold choice.
- The uncertainty terms cover the **threshold and shoreline mixed pixels** only; geolocation,
  radiometric and temporal-mixing errors are not quantified.
- **No independent validation yet:** the JRC comparison is optical-vs-optical. Its scope is also
  limited: the IoU is a single date; *r* covers 2017–2021 only (the JRC monthly series ends in
  2021); the >80%-coverage filter biases toward clearer months; and part of *r* reflects the
  **shared seasonal cycle** of both series.
- **Administrative boundary:** from a WFS-derived CAOP GeoJSON; the exact CAOP edition is
  unconfirmed. The parish geometry is stable since the 2013 reorganization.
- The water area is a **proxy for the system state**, driven by both rainfall and reservoir
  management — the comparison describes differences, it does **not** attribute their causes.

## Reproduce

Requires a free (non-commercial) **Google Earth Engine** account and Cloud project.

```bash
pip install -r requirements.txt
earthengine authenticate
```

Set your Earth Engine project id in the scripts, then run in order:

| Script | Does |
|--------|------|
| `01_explorar.py` | GEE connectivity + first Sentinel-2 preview |
| `02_agua.py` | Parish boundary + MNDWI water + area |
| `03_serie.py` | Annual summer series |
| `04_trimestral.py` | Quarterly series with cloud masking + coverage flag |
| `05_comparar.py` | Per-year water polygons (GeoJSON) + change map |
| `06_validacao.py` | Consistency check vs JRC GSW (IoU, Pearson r) |
| `07_incerteza.py` | Threshold-sensitivity uncertainty |
| `08_margem.py` | Shoreline length and mixed-pixel uncertainty (20 m grid) |
| `09_envolvente.py` | Quality check: water outside the reservoir envelope (ponds) |

Serve the web app locally (it reads the generated data files, so it needs a server — not a
double-click):

```bash
python -m http.server 8000
# then open http://localhost:8000
```

## Data sources & credits

- **Copernicus Sentinel-2** (ESA), via Google Earth Engine.
- **JRC Global Surface Water** — Pekel, J.-F., Cottam, A., Gorelick, N. & Belward, A. S. (2016).
  *High-resolution mapping of global surface water and its long-term changes.* Nature 540, 418–422.
- **CAOP** — Carta Administrativa Oficial de Portugal, Direção-Geral do Território.
- Basemap tiles: Esri, Maxar, Earthstar Geographics.

## Author

**Gustavo Santa** · ORCID [0009-0008-4266-0732](https://orcid.org/0009-0008-4266-0732)

Personal portfolio project. _License: MIT (suggested — add a LICENSE file if you agree)._
