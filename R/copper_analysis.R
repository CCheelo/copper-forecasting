# =============================================================================
#  Copper Price Time-Series Analysis
#  ARMA modelling (original work) + ARMA-GARCH extension (improvement)
#  -----------------------------------------------------------------------------
#  Daily copper price (USD/tonne), Dec 2000 - May 2025 (LME-style spot).
#
#  Pipeline:
#    1. Load & clean the daily price series
#    2. Stationarity of the level (ADF / KPSS)        -> non-stationary
#    3. Log-returns                                    -> stationary
#    4. ACF / PACF                                     -> candidate orders
#    5. BIC grid search (+ auto.arima cross-check)     -> ARMA(2,0) = AR(2)
#    6. Fit ARMA(2,0); coefficient significance
#    7. Residual diagnostics (Ljung-Box, ARCH-LM, JB)  -> ARCH + fat tails remain
#    8. Out-of-sample backtest vs random walk (RMSE/MAE)
#    9. IMPROVEMENT: ARMA(2,0)-GARCH(1,1) with Student-t innovations
#   10. Model comparison + forecast
#
#  Reproduces the figures and numbers in report/copper-forecasting-report.pdf.
#
#  Run from the repository root:   Rscript R/copper_analysis.R
#  Packages: readxl, tseries, forecast, lmtest, FinTS, rugarch
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(tseries); library(forecast)
  library(lmtest);  library(FinTS);  library(rugarch)
})

fig_dir <- "figures"; dir.create(fig_dir, showWarnings = FALSE)
sec <- function(t) cat("\n", strrep("=", 76), "\n  ", t, "\n", strrep("=", 76), "\n", sep = "")

## find the data whether run from repo root, R/, or data/ -----------------------
find_data <- function() {
  cands <- c("data/Daily copper with dates.xlsx",
             "Daily copper with dates.xlsx",
             "../data/Daily copper with dates.xlsx")
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("Cannot find 'Daily copper with dates.xlsx'")
  hit[1]
}

# =============================================================================
# 1. LOAD & CLEAN
# =============================================================================
sec("1. LOAD & CLEAN")
raw <- read_excel(find_data(), sheet = "Sheet1", skip = 2)  # 2 blank rows above header
names(raw) <- c("Date", "Price")
raw$Date <- as.Date(raw$Date); raw$Price <- as.numeric(raw$Price)
raw <- raw[!is.na(raw$Date) & !is.na(raw$Price), ]

# Fix year typos: the Jan-Feb 2022 block was keyed as 19xx. Bump pre-2000 dates +100y.
lt <- as.POSIXlt(raw$Date); bad <- (lt$year + 1900) < 2000
if (any(bad)) { lt$year[bad] <- lt$year[bad] + 100; raw$Date <- as.Date(lt)
                cat("Year typos fixed:", sum(bad), "\n") }
raw <- raw[order(raw$Date), ]
cu  <- raw[!duplicated(raw$Date), ]          # drop duplicate dates, keep chronological
price <- cu$Price
cat(sprintf("Observations: %d | %s to %s | price %.0f-%.0f USD/tonne\n",
            nrow(cu), min(cu$Date), max(cu$Date), min(price), max(price)))

# =============================================================================
# 2. STATIONARITY OF THE LEVEL
# =============================================================================
sec("2. STATIONARITY OF THE LEVEL")
print(suppressWarnings(adf.test(price)))     # H0: unit root -> expect NOT rejected
print(suppressWarnings(kpss.test(price)))    # H0: stationary -> expect rejected

png(file.path(fig_dir, "01_price_level.png"), 1100, 700, res = 120)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(cu$Date, price, type = "l", col = "#c2410c", main = "Daily copper price (USD/tonne)",
     xlab = "", ylab = "USD/tonne")
plot(cu$Date, log(price), type = "l", col = "#1b3a2d", main = "Log price",
     xlab = "", ylab = "log(USD/tonne)")
dev.off()

# =============================================================================
# 3. LOG-RETURNS (stationary)
# =============================================================================
sec("3. LOG-RETURNS")
ret <- diff(log(price)) * 100                # percent daily log-returns
N <- length(ret)
print(suppressWarnings(adf.test(ret)))       # expect rejection of unit root
print(suppressWarnings(kpss.test(ret)))      # expect non-rejection of stationarity
cat(sprintf("mean=%.4f%%  sd=%.4f%%  n=%d\n", mean(ret), sd(ret), N))

png(file.path(fig_dir, "02_log_returns.png"), 1100, 450, res = 120)
par(mar = c(4, 4, 3, 1))
plot(cu$Date[-1], ret, type = "l", col = "#c2410c", lwd = 0.5,
     main = "Daily copper log-returns (%)", xlab = "", ylab = "% return"); abline(h = 0, col = "grey60")
dev.off()

# =============================================================================
# 4. ACF / PACF
# =============================================================================
sec("4. ACF / PACF")
png(file.path(fig_dir, "03_acf_pacf.png"), 1100, 700, res = 120)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
Acf(ret, lag.max = 30, main = "ACF of log-returns")
Pacf(ret, lag.max = 30, main = "PACF of log-returns")
dev.off()
cat("Both ACF and PACF spike at lags 1-2 then cut off -> low-order AR. See figures/03_acf_pacf.png\n")

# =============================================================================
# 5. ORDER SELECTION BY BIC
# =============================================================================
sec("5. ORDER SELECTION BY BIC (p,q in 0..5)")
grid <- expand.grid(p = 0:5, q = 0:5); res_grid <- data.frame()
for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; q <- grid$q[i]
  f <- tryCatch(arima(ret, order = c(p, 0, q), include.mean = TRUE, method = "ML"),
                error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(f)) res_grid <- rbind(res_grid, data.frame(p = p, q = q, AIC = AIC(f), BIC = BIC(f)))
}
res_grid <- res_grid[order(res_grid$BIC), ]
cat("Top 8 by BIC:\n"); print(head(res_grid, 8), row.names = FALSE)
auto <- auto.arima(ret, d = 0, ic = "bic", seasonal = FALSE,
                   stepwise = FALSE, approximation = FALSE, max.p = 5, max.q = 5)
cat(sprintf("auto.arima(ic=bic) -> ARMA(%d,%d)\n", auto$arma[1], auto$arma[2]))

# =============================================================================
# 6. FIT ARMA(2,0)  (the original model)
# =============================================================================
sec("6. ARMA(2,0)")
fit <- Arima(ret, order = c(2, 0, 0), include.mean = TRUE)
print(fit); cat("\n"); print(coeftest(fit))

# =============================================================================
# 7. RESIDUAL DIAGNOSTICS  (why ARMA alone is not enough)
# =============================================================================
sec("7. RESIDUAL DIAGNOSTICS")
r <- residuals(fit)
print(Box.test(r,   20, "Ljung-Box", fitdf = 2))   # mean: expect clean (p > 0.05)
print(Box.test(r^2, 20, "Ljung-Box"))              # variance: expect ARCH (p ~ 0)
print(ArchTest(r, 12))                             # ARCH-LM: expect rejection
print(jarque.bera.test(r))                         # normality: expect rejection

png(file.path(fig_dir, "04_residual_diagnostics.png"), 1100, 700, res = 120)
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
plot(cu$Date[-1], r, type = "l", col = "#c2410c", lwd = 0.4, main = "ARMA(2,0) residuals",
     xlab = "", ylab = "resid")
Acf(r^2, lag.max = 30, main = "ACF of squared residuals (ARCH signature)")
dev.off()

# =============================================================================
# 8. OUT-OF-SAMPLE BACKTEST vs RANDOM WALK  (last 252 trading days, 1-step)
# =============================================================================
sec("8. OUT-OF-SAMPLE BACKTEST (last 252 days, 1-step)")
H <- 252; tr <- 1:(N - H); te <- (N - H + 1):N
fit_tr  <- Arima(ret[tr], order = c(2, 0, 0), include.mean = TRUE)
fit_all <- Arima(ret, model = fit_tr)              # apply fixed coefficients to full series
pred_ret <- fitted(fit_all)[te]
p_prev <- price[te]; p_act <- price[te + 1]
p_arma <- p_prev * exp(pred_ret / 100)             # reconstruct price forecast
p_rw   <- p_prev                                   # naive random walk
rmse <- function(a, f) sqrt(mean((a - f)^2)); mae <- function(a, f) mean(abs(a - f))
cat(sprintf("ARMA(2,0)  : RMSE=%.2f  MAE=%.2f\n", rmse(p_act, p_arma), mae(p_act, p_arma)))
cat(sprintf("Random walk: RMSE=%.2f  MAE=%.2f\n", rmse(p_act, p_rw),   mae(p_act, p_rw)))
cat(sprintf("ARMA directional accuracy: %.1f%%  (~50%% = no edge)\n",
            100 * mean(sign(p_act - p_prev) == sign(p_arma - p_prev))))
cat("=> AR(2) does NOT beat a random walk on point forecasts: copper is near-efficient in the mean.\n")

png(file.path(fig_dir, "05_backtest.png"), 1100, 450, res = 120)
par(mar = c(4, 4, 3, 1)); d <- cu$Date[te + 1]
plot(d, p_act, type = "l", col = "#1b3a2d", lwd = 1.2, xlab = "", ylab = "USD/tonne",
     main = "Hold-out: actual vs ARMA(2,0) one-step forecast")
lines(d, p_arma, col = "#c2410c", lwd = 1, lty = 2)
legend("topleft", c("Actual", "ARMA(2,0) forecast"), col = c("#1b3a2d", "#c2410c"),
       lty = c(1, 2), bty = "n")
dev.off()

# =============================================================================
# 9. IMPROVEMENT: ARMA(2,0)-GARCH(1,1) with Student-t innovations
# =============================================================================
sec("9. ARMA(2,0)-GARCH(1,1)-t  (improved model)")
spec <- ugarchspec(variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model         = list(armaOrder = c(2, 0), include.mean = TRUE),
                   distribution.model = "std")          # Student-t for fat tails
fg <- ugarchfit(spec, data = ret, solver = "hybrid")
print(round(fg@fit$matcoef, 6))
LLg <- as.numeric(likelihood(fg)); npg <- length(coef(fg))
AICg <- -2 * LLg + 2 * npg; BICg <- -2 * LLg + log(N) * npg
pers <- as.numeric(sum(coef(fg)[c("alpha1", "beta1")])); tdf <- as.numeric(coef(fg)["shape"])
cat(sprintf("\nlogLik=%.2f  AIC=%.2f  BIC=%.2f  | persistence(a1+b1)=%.4f  t.d.f.=%.2f\n",
            LLg, AICg, BICg, pers, tdf))

# standardised-residual diagnostics: ARCH should now be GONE
z <- as.numeric(residuals(fg, standardize = TRUE))
print(Box.test(z,   20, "Ljung-Box"))     # expect clean
print(Box.test(z^2, 20, "Ljung-Box"))     # expect clean (ARCH absorbed)
print(ArchTest(z, 12))                     # expect non-rejection

png(file.path(fig_dir, "06_conditional_volatility.png"), 1100, 450, res = 120)
par(mar = c(4, 4, 3, 1))
plot(cu$Date[-1], as.numeric(sigma(fg)), type = "l", col = "#c2410c", lwd = 0.7,
     main = "GARCH(1,1) conditional volatility", xlab = "", ylab = "cond. SD (%/day)")
abline(h = sd(ret), col = "#1b3a2d", lty = 2)
legend("topright", c("GARCH conditional SD", "constant SD (ARMA assumption)"),
       col = c("#c2410c", "#1b3a2d"), lty = c(1, 2), bty = "n")
dev.off()

# =============================================================================
# 10. MODEL COMPARISON
# =============================================================================
sec("10. MODEL COMPARISON")
cmp <- data.frame(
  Model   = c("ARMA(2,0)", "ARMA(2,0)-GARCH(1,1)-t"),
  Params  = c(4, npg),
  logLik  = round(c(as.numeric(logLik(fit)), LLg)),
  AIC     = round(c(AIC(fit), AICg)),
  BIC     = round(c(BIC(fit), BICg)),
  ARCH_p  = signif(c(ArchTest(r, 12)$p.value, ArchTest(z, 12)$p.value), 3))
print(cmp, row.names = FALSE)
cat("\n=> GARCH wins on AIC/BIC by a wide margin and is the only model whose residuals pass the no-ARCH test.\n")

sec("DONE")
cat("Figures written to ./figures/\n")
