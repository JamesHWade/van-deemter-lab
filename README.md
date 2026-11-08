# Van Deemter Lab

An undergraduate teaching app built with Shiny for Python. Explore how velocity affects plate height, theoretical plate count and resolution, while controlling length, retention and selectivity.

## Run

Requires Python 3.10 or newer. From this directory:

```sh
python -m venv .venv
source .venv/bin/activate
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
shiny run --reload app.py
```

Open the local address printed by Shiny. For Posit Connect, deploy this directory as a Shiny for Python application with `app.py` as the entry point.

## Included

- Linked van Deemter, chromatogram, resolution and retention-time plots.
- Adjustable A, B, C, linear velocity, column length, retention and selectivity.
- Optimum-velocity button; saved comparison; CSV export of settings, metrics and peaks.
- Normalized time view to distinguish relative overlap from absolute peak width.
- Equation explanations, model assumptions and five guided experiments with instructor notes.

## Scientific model

All scenarios are illustrative, not fitted to a particular column. H is in micrometres, L in millimetres and u in millimetres per second. N = 1000 L/H. Retention factors are k1 and k2 = alpha k1. Times are in seconds internally. Both solutes have equal plate count and unit integrated peak area. Peak standard deviations are tR/sqrt(N); exact model resolution uses tangent-baseline widths of 4 sigma. The familiar later-peak-width approximation is reported separately.

Assumes isocratic elution, constant coefficients, Gaussian peaks and no extra-column broadening, overload, pressure limitation or temperature effect. Optimum velocity minimizes H, not run time. Volumetric flow cannot be inferred without column geometry and porosity.

Background: [Agilent, The LC Handbook](https://www.agilent.com/cs/library/primers/public/LC-Handbook-Complete-2.pdf). Framework: [Shiny for Python](https://shiny.posit.co/py/).

## Validation

```sh
python -m unittest -v test_model.py
```

Seven automated checks cover units, optimum, length scaling, coelution, peak integration, resolution and boundary inputs. App import and server startup were checked. Full browser interaction testing was attempted but the execution environment blocked Chromium process sockets; browser-level behavior remains unverified.
