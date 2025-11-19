# Libraries 
library(quantmod)    
library(tseries)     
library(forecast)    

#Download data from Yahoo
getSymbols("^GSPC", 
           src = "yahoo", 
           from = "1990-01-01", 
           to   = "2025-10-24",  
           auto.assign = TRUE)

# Use Adjusted Close for return calculations
px_d <- Ad(GSPC)  #

# Aggregate to weekly 
wk   <- to.weekly(px_d, indexAt = "endof", drop.time = TRUE)
px_w <- Cl(wk)                 # weekly adjusted close

# Box–Cox Transformation Check
bc <- BoxCox.lambda(px_w)
bc

# Weekly log returns (stationary)
r_w  <- na.omit(diff(log(px_w)))

# Plots: price & returns 
par(mfrow = c(2,1))
plot(px_w, main = "S&P 500 — Weekly Adjusted Close", ylab = "Index level")
plot(r_w,  main = "S&P 500 — Weekly Log Returns", ylab = "log return")
par(mfrow = c(1,1))

# ACF/PACF on log returns 
par(mfrow = c(1,2))
acf(r_w,  main = "ACF: Weekly returns")
pacf(r_w, main = "PACF: Weekly returns")
par(mfrow = c(1,1))

# Stationarity (returns) ADF test 
adf_ret <- adf.test(as.numeric(r_w)); adf_ret

# Rolling mean/variance (descriptive; variance clustering expected) 
roll_mean <- rollapply(r_w, width = 52, FUN = mean, by = 1, align = "right", fill = NA)
roll_var  <- rollapply(r_w, width = 52, FUN = var,  by = 1, align = "right", fill = NA)
par(mfrow=c(2,1))
plot(roll_mean, main = "52-week rolling mean of returns", ylab = "mean"); abline(h = 0, lty = 2)
plot(roll_var,  main = "52-week rolling variance of returns", ylab = "variance")
par(mfrow=c(1,1))

# test if mean is statistically non-zero 
mean_ret <- mean(r_w)
t_mean   <- t.test(as.numeric(r_w), mu = 0)
mean_ret; t_mean$p.value
inc_mean <- TRUE

# Candidate grid + fits 
cands <- list(
  c(1,0,0),  # AR(1)
  c(2,0,0),  # AR(2)
  c(0,0,1),  # MA(1)
  c(0,0,2),  # MA(2)
  c(1,0,1),  # ARMA(1,1)
  c(2,0,1),  # ARMA(2,1)
  c(1,0,2)   # ARMA(1,2)
)
cand_names <- c("AR(1)","AR(2)","MA(1)","MA(2)","ARMA(1,1)","ARMA(2,1)","ARMA(1,2)")

fit_one_base <- function(ord, x, inc_mean=TRUE){
  suppressWarnings( arima(x, order = ord, include.mean = inc_mean, method = "ML") )
}

# AICc 
AICc_manual <- function(model) {
  n <- length(model$residuals)
  k <- length(coef(model)) + 1   # +1 for sigma^2
  AIC(model) + (2*k*(k+1)) / (n - k - 1)
}

fits_base <- lapply(cands, \(o) fit_one_base(o, r_w, inc_mean))
names(fits_base) <- cand_names

IC <- data.frame(
  Model = names(fits_base),
  p = sapply(cands, \(z) z[1]),
  d = sapply(cands, \(z) z[2]),
  q = sapply(cands, \(z) z[3]),
  AIC  = sapply(fits_base, AIC),
  AICc = sapply(fits_base, AICc_manual),
  BIC  = sapply(fits_base, BIC),
  row.names = NULL
)
IC <- IC[order(IC$AICc), ]
IC


#Top Models (ARMA(1,1) & AR(2))
top <- IC[1:2, c("Model","p","d","q")]
top

fit_top1 <- fits_base[[ top$Model[1] ]]
fit_top2 <- fits_base[[ top$Model[2] ]]

#Summary
print(summary(fit_top1))
print(summary(fit_top2))



#Residual Dianostics
library(astsa)
astsa::sarima(r_w, top$p[1], top$d[1], top$q[1],
              no.constant = !inc_mean, details = TRUE)
astsa::sarima(r_w, top$p[2], top$d[2], top$q[2],
              no.constant = !inc_mean, details = TRUE)


# 52w holdout

# Split data
y_xts <- r_w
y     <- as.numeric(y_xts)
ti    <- index(y_xts)

h <- 52
n <- length(y)
train_idx <- 1:(n - h)
test_idx  <- (n - h + 1):n

y_train  <- y[train_idx]
y_test   <- y[test_idx]
ti_train <- ti[train_idx]
ti_test  <- ti[test_idx]

# Top-2 orders 
ord1  <- c(top$p[1], top$d[1], top$q[1]); name1 <- top$Model[1]
ord2  <- c(top$p[2], top$d[2], top$q[2]); name2 <- top$Model[2]

# Walk-forward helper 
wf_forecast <- function(y_all, train_end_index, horizon_idx, order_vec,
                        include.mean = TRUE, band = 0.95) {
  k <- length(horizon_idx)
  fmean <- fse <- err <- rep(NA_real_, k)
  L <- U <- rep(NA_real_, k)
  z <- qnorm((1 + band)/2)
  for (i in seq_len(k)) {
    fit_end <- train_end_index + i - 1
    fit_y   <- y_all[1:fit_end]
    fit     <- arima(fit_y, order = order_vec, include.mean = include.mean, method = "ML")
    pr      <- predict(fit, n.ahead = 1)
    fmean[i] <- pr$pred[1]
    fse[i]   <- pr$se[1]
    L[i]     <- fmean[i] - z * fse[i]
    U[i]     <- fmean[i] + z * fse[i]
    y_act    <- y_all[horizon_idx[i]]
    err[i]   <- y_act - fmean[i]
  }
  list(mean = fmean, se = fse, L = L, U = U, error = err)
}

# 52-week walk forward evaluation
res1 <- wf_forecast(y, max(train_idx), test_idx, ord1, include.mean = inc_mean, band = 0.95)
res2 <- wf_forecast(y, max(train_idx), test_idx, ord2, include.mean = inc_mean, band = 0.95)

RMSE1 <- sqrt(mean(res1$error^2)); MAE1 <- mean(abs(res1$error))
RMSE2 <- sqrt(mean(res2$error^2)); MAE2 <- mean(abs(res2$error))

cat("\n--- Holdout (52w) 1-step-ahead errors ---\n")
cat(sprintf("%-12s RMSE = %.6f | MAE = %.6f | order=c(%d,%d,%d)\n",
            name1, RMSE1, MAE1, ord1[1], ord1[2], ord1[3]))
cat(sprintf("%-12s RMSE = %.6f | MAE = %.6f | order=c(%d,%d,%d)\n",
            name2, RMSE2, MAE2, ord2[1], ord2[2], ord2[3]))

# 52w plots
years_back <- 156  # ~3 years; change to 104 for ~2 years
start_zoom <- max(1, (n - h) - years_back + 1)

plot_zoom_panel <- function(tag, res_obj) {
  ylim <- range(c(y[start_zoom:n], res_obj$L, res_obj$U), na.rm = TRUE)
  plot(ti[start_zoom:n], y[start_zoom:n],
       type = "l", col = "black", lwd = 1.2,
       xlab = "", ylab = "Weekly log return",
       main = sprintf("%s — 52w Walk-forward", tag),
       ylim = ylim, xaxt = "n", yaxt = "s", bty = "l")
  axis.Date(1, at = pretty(ti[start_zoom:n]), format = "%Y-%m")
  abline(v = ti_test[1], col = "green")                 # holdout start
  lines(ti_test, y_test, col = "black", lwd = 1.5)      # actual in holdout
  lines(ti_test, res_obj$mean, col = "red", type = "o") # forecast mean
  lines(ti_test, res_obj$U, col = "blue", lty = "dashed")
  lines(ti_test, res_obj$L, col = "blue", lty = "dashed")
}

op <- par(no.readonly = TRUE)


par(mfrow = c(2,1), mar = c(4,4,2,1), mgp = c(2.2, .8, 0))

plot_zoom_panel(name1, res1)  # ARMA(1,1)
plot_zoom_panel(name2, res2)  # AR(2)

par(op)

# Diebold–Mariano test 
suppressWarnings(suppressMessages(require(forecast)))
dm_out <- dm.test(e1 = res1$error, e2 = res2$error, alternative = "two.sided", h = 1, power = 2)
print(dm_out)



# RMSE / MAE comparison table
results <- data.frame(
  Model = c(name1, name2),
  p = c(ord1[1], ord2[1]),
  d = c(ord1[2], ord2[2]),
  q = c(ord1[3], ord2[3]),
  RMSE = c(RMSE1, RMSE2),
  MAE  = c(MAE1,  MAE2)
)

print(results, row.names = FALSE)

winner <- list(name = name1, order = ord1)

# Compute t-values and p-values 
est <- fit_top1$coef
se  <- sqrt(diag(fit_top1$var.coef))
tval <- est / se
pval <- 2 * (1 - pnorm(abs(tval)))

data.frame(
  Parameter = names(est),
  Estimate  = est,
  Std.Error = se,
  t.value   = tval,
  p.value   = pval
)


# 26-week RETURNS plot
fit_full <- arima(y, order = ord1, include.mean = inc_mean, method = "ML")
fc26     <- predict(fit_full, n.ahead = 26)

z95    <- qnorm(0.975)
L95ret <- fc26$pred - z95 * fc26$se
U95ret <- fc26$pred + z95 * fc26$se

ti_future <- seq(from = ti[n], by = 7, length.out = 27)[-1]

pre_weeks  <- 52
post_weeks <- 26
x_left  <- as.Date(ti[n]) - pre_weeks*7
x_right <- as.Date(ti[n]) + post_weeks*7

start_idx  <- max(1, n - pre_weeks + 1)
ti_obs_win <- ti[start_idx:n]
y_obs_win  <- y[start_idx:n]
ylim_zoom  <- range(c(y_obs_win, L95ret, U95ret), na.rm = TRUE)

op <- par(no.readonly = TRUE)
par(mfrow = c(1,1), mar = c(4,4,2,1), mgp = c(2.2, .8, 0))

plot(ti_obs_win, y_obs_win, type = "l", col = "black", lwd = 1.2,
     xlab = "Week", ylab = "Weekly log return",
     main = sprintf("26-week forecast — %s (red) with 95%% bands", winner$name),
     xlim = c(x_left, x_right), ylim = ylim_zoom,
     xaxt = "n", yaxt = "s", bty = "l")
axis.Date(1, at = seq(x_left, x_right, by = "1 month"), format = "%Y-%m")
abline(v = ti[n], col = "green")
lines(ti_future, fc26$pred, col = "red", type = "o")
lines(ti_future, U95ret,   col = "blue", lty = "dashed")
lines(ti_future, L95ret,   col = "blue", lty = "dashed")

legend("topleft",
       c("Actual (last 52w)", "Forecast mean (26w)", "95% band", "Last observed"),
       col = c("black","red","blue","green"), lty = c(1,1,2,1), pch = c(NA,1,NA,NA), bty = "n")
par(op)


# 26-week PRICE plot
P_last <- as.numeric(last(px_w))
mu_ret <- as.numeric(fc26$pred)
se_ret <- as.numeric(fc26$se)

cum_mu <- cumsum(mu_ret)
P_mean <- P_last * exp(cum_mu)

set.seed(123)
B <- 5000
R_sim <- matrix(rnorm(B * length(mu_ret), mean = rep(mu_ret, each=B), sd = rep(se_ret, each=B)),
                nrow = B)
Cum_sim <- t(apply(R_sim, 1, cumsum))
P_sim   <- P_last * exp(Cum_sim)
P_L95   <- apply(P_sim, 2, quantile, probs = 0.025, na.rm = TRUE)
P_U95   <- apply(P_sim, 2, quantile, probs = 0.975, na.rm = TRUE)

pre_weeks_price <- 52
x_leftP  <- as.Date(ti[n]) - pre_weeks_price*7
x_rightP <- as.Date(ti[n]) + 26*7
px_obs   <- px_w[paste0(as.Date(x_leftP), "/", as.Date(ti[n]))]
ylimP    <- range(c(as.numeric(px_obs), P_L95, P_U95), na.rm = TRUE)

op <- par(no.readonly = TRUE)
par(mfrow = c(1,1), mar = c(4,4,2,1), mgp = c(2.2, .8, 0))

plot(index(px_obs), as.numeric(px_obs),
     type = "l", col = "black", lwd = 1.2,
     xlab = "Week", ylab = "Price (Adj Close)",
     main = sprintf("Price-level 26-week forecast — %s", winner$name),
     xlim = c(x_leftP, x_rightP), ylim = ylimP, xaxt = "n", bty = "l")
axis.Date(1, at = seq(x_leftP, x_rightP, by = "1 month"), format = "%Y-%m")
abline(v = ti[n], col = "green")
lines(ti_future, P_mean, col = "red", type = "o")
lines(ti_future, P_U95,  col = "blue", lty = "dashed")
lines(ti_future, P_L95,  col = "blue", lty = "dashed")

legend("topleft",
       c("Actual price (last 52w)", "Forecast mean (price)", "95% band (sim.)", "Last observed"),
       col = c("black","red","blue","green"), lty = c(1,1,2,1), pch = c(NA,1,NA,NA), bty = "n")
par(op)

