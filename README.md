# Copper Price Forecasting — ARMA & ARMA‑GARCH

Time‑series modelling of **daily copper prices** (USD/tonne, Dec 2000 – May 2025) in R.
The project starts from the classic Box–Jenkins workflow — identify a model from the
ACF/PACF, select the order by BIC — lands on an **AR(2)** model for the conditional mean,
then asks an honest follow‑up question: *is that model actually any good?* It isn't, quite —
so the analysis extends it to an **ARMA(2,0)–GARCH(1,1)** model that captures what the AR(2)
misses.

> **The one‑sentence finding:** copper prices are close to a random walk *in the mean*
> (you can't reliably forecast the direction of tomorrow's price), but their **volatility is
> highly predictable** — and that is the thing worth modelling.

📄 **Full write‑up:** [`report/copper-forecasting-report.pdf`](report/copper-forecasting-report.pdf)
🔗 **Portfolio:** [choolwecheelo.com](https://www.choolwecheelo.com)

---

## Why copper?

Copper is the backbone of the Zambian economy — it consistently accounts for the majority of
the country's export earnings, so the behaviour of its price feeds directly into the national
budget, the exchange rate, and mining‑sector investment. Understanding how predictable that
price is (and isn't) is therefore of real practical interest.

## The data

| | |
|---|---|
| File | `data/Daily copper with dates.xlsx` |
| Observations | 6,298 daily prices (after cleaning) |
| Period | 1 Dec 2000 → 30 May 2025 |
| Units | US$ / tonne (LME‑style cash settlement) |
| Range | ≈ US$1,319 → US$10,889 |

**Cleaning applied** (all reproduced in code):
- **36 year typos fixed** — the entire Jan–Feb **2022** block had been keyed as **1922**;
  any pre‑2000 date is shifted forward 100 years.
- **1 duplicate date** removed; series sorted chronologically.

## Method

1. **Stationarity** — the price *level* is non‑stationary (ADF fails to reject a unit root;
   KPSS rejects stationarity), so the analysis works with percent **log‑returns**
   `r_t = 100·ln(P_t / P_{t‑1})`, which are stationary.
2. **Identification** — ACF and PACF of returns both spike at lags 1–2 then cut off.
3. **Order selection** — a BIC grid search over `p,q ∈ {0..5}` (cross‑checked with
   `auto.arima`) selects **ARMA(2,0) = AR(2)**.
4. **Diagnostics** — Ljung–Box, ARCH‑LM, Jarque–Bera.
5. **Out‑of‑sample backtest** — one‑step forecasts over the last 252 trading days vs a naive
   random walk (RMSE / MAE / directional accuracy).
6. **Improvement** — **ARMA(2,0)–GARCH(1,1)** with Student‑*t* innovations.

## Key results

**The AR(2) is valid in the mean but inadequate overall.** Its residuals are free of
autocorrelation (Ljung–Box p ≈ 0.10) but fail decisively on volatility clustering
(ARCH‑LM p ≈ 0) and normality, and out of sample it does **not** beat a random walk:

| Model (1‑step, last 252 days) | RMSE | MAE |
|---|---:|---:|
| ARMA(2,0) | 118.3 | 86.2 |
| Random walk | **117.7** | **84.9** |

Directional accuracy of the AR(2): **48.8%** — no better than a coin toss.

**The GARCH extension fixes the diagnostics and wins decisively on fit:**

| Model | k | logLik | AIC | BIC | ARCH‑LM p |
|---|--:|--:|--:|--:|--:|
| ARMA(2,0) | 4 | −12,354 | 24,715 | 24,742 | ~0 |
| ARMA(2,0)–GARCH(1,1)‑t | 7 | **−10,906** | **21,825** | **21,872** | **1.0** |

- Volatility **persistence** `α₁+β₁ = 0.974` (shocks decay slowly).
- Student‑*t* **degrees of freedom ≈ 4.6** (heavy tails, strongly significant).
- The standardised residuals show **no remaining ARCH** — the clustering is fully captured.

**Why the GARCH model is the better model:** it gives valid standard errors, *time‑varying*
(honest) forecast intervals, and accurate tail‑risk / value‑at‑risk — even though it doesn't
make the price more predictable in the mean. For copper, the risk is the predictable part.

## Reproduce it

Requires **R ≥ 4.x** and the packages `readxl`, `tseries`, `forecast`, `lmtest`, `FinTS`,
`rugarch`:

```r
install.packages(c("readxl","tseries","forecast","lmtest","FinTS","rugarch"))
```

From the repository root:

```bash
Rscript R/copper_analysis.R
```

This prints the full analysis to the console and writes the figures to `figures/`.
The PDF report is built from the literate‑programming source `report/copper_report.Rnw`
(knitr + LaTeX), which embeds the live R output.

## Repository layout

```
copper-forecasting/
├── data/
│   └── Daily copper with dates.xlsx     # raw daily price series
├── R/
│   └── copper_analysis.R                # full, commented analysis (ARMA + backtest + GARCH)
├── report/
│   ├── copper_report.Rnw                # literate source (LaTeX + R)
│   └── copper-forecasting-report.pdf    # compiled report
├── figures/                             # generated plots (PNG)
└── README.md
```

## Caveats

The model is **univariate** — it uses only the price's own history, with no macro covariates
(USD index, Chinese demand, inventories). The point‑forecast conclusion (random walk is
hard to beat) is specific to short horizons. Natural extensions: asymmetric volatility
(EGARCH/GJR), ARIMAX with exogenous drivers, and regime‑switching for the boom/crash phases.

---

*Choolwe Cheelo · BSc Actuarial Science · analysis in R.*
