# sp500-arima-model

# S&P 500 Time Series Forecasting: An ARMA Modeling Approach

## Overview

This project implements a comprehensive time series analysis workflow for forecasting S&P 500 weekly log returns using ARMA (Autoregressive Moving Average) models. The analysis spans 35 years of market data (January 1990 - October 2025) and demonstrates the complete process from exploratory data analysis to model validation and forecasting.

## Project Highlights

- **Data**: 1,868 weeks of S&P 500 adjusted closing prices from Yahoo Finance
- **Final Model**: ARMA(1,1) selected through rigorous comparison using AICc, BIC, and out-of-sample validation
- **Validation**: 52-week walk-forward cross-validation with Diebold-Mariano testing
- **Forecasts**: 26-week ahead predictions for both returns and price levels with 95% confidence intervals

## Key Features

### Statistical Analysis
- Stationarity testing using Augmented Dickey-Fuller (ADF) test
- ACF/PACF analysis for model identification
- Rolling window analysis (52-week mean and variance)
- Comprehensive residual diagnostics (Ljung-Box tests, Q-Q plots)

### Model Selection
- Systematic comparison of 7 candidate models (AR, MA, and ARMA specifications)
- Information criteria comparison (AIC, AICc, BIC)
- Out-of-sample validation using walk-forward methodology
- Statistical comparison via Diebold-Mariano test

### Forecasting
- Return-level forecasts with constant uncertainty bands
- Price-level forecasts using Monte Carlo simulation (5,000 paths)
- Visualization of both historical data and future projections

## Methodology

1. **Data Preparation**: Convert daily prices to weekly log returns to achieve stationarity
2. **Model Identification**: Use ACF/PACF patterns to identify candidate models
3. **Model Estimation**: Fit models using maximum likelihood estimation
4. **Model Selection**: Compare models using information criteria and residual diagnostics
5. **Validation**: Perform 52-week walk-forward cross-validation
6. **Forecasting**: Generate 26-week forecasts with confidence intervals

## Requirements
```r
library(quantmod)   # Data retrieval
library(tseries)    # Time series tests
library(forecast)   # Forecasting tools
library(astsa)      # Diagnostic plots
```

## Usage

1. **Clone the repository**
```bash
git clone https://github.com/yourusername/sp500-time-series-forecast.git
cd sp500-time-series-forecast
```

2. **Run the R script**
```r
source("Time_Series_Project_1.R")
```

The script will:
- Download S&P 500 data automatically
- Perform all analyses and generate diagnostic plots
- Output model comparison tables
- Generate forecast plots for both returns and prices

## Key Results

### Selected Model: ARMA(1,1)
```
r_t = μ + φ₁r_{t-1} + θ₁ε_{t-1} + ε_t

φ₁ = -0.5425 (p = 0.0036)
θ₁ = 0.4601 (p = 0.019)
μ = 0.0016 (p = 0.0014)
```

### Out-of-Sample Performance
- **RMSE**: 0.0241
- **MAE**: 0.0177
- No significant difference from AR(2) competitor (Diebold-Mariano p = 0.217)

### Key Findings
- Weak negative autocorrelation at lag 1 suggests mild mean reversion
- Forecasts rapidly converge to unconditional mean (~0.16% weekly return)
- Price-level forecast uncertainty widens dramatically over 26-week horizon
- Results consistent with efficient market hypothesis

## Visualizations

The project generates several key visualizations:
- Price levels and log returns time series
- ACF/PACF plots
- 52-week rolling statistics
- Residual diagnostic plots
- Walk-forward forecast plots
- 26-week return and price forecasts with confidence bands

## Limitations & Future Work

### Current Limitations
- Assumes constant variance (ignores volatility clustering)
- Uses only historical returns (no external variables)
- Weekly aggregation may mask higher-frequency dynamics

### Planned Extensions
- **GARCH Models**: Capture time-varying volatility
- **ARIMAX/VAR Models**: Incorporate macroeconomic variables
- **Higher Frequency**: Explore daily data patterns
- **Alternative Models**: Compare with machine learning approaches

## Project Structure
```
├── S&P 500 ARIMA.R          # Main R script
├── S&P 500 ARIMA.pdf        # Full project report
└── README.md                        # This file
```

## References

Shumway, R. H., & Stoffer, D. S. (2017). *Time Series Analysis and Its Applications: With R Examples* (4th ed.). Springer.

## Author

Matthew Dikman  
CUNY Baruch - Statistics Graduate Program  
October 2025

## License

This project is available for educational and research purposes.

---

**Note**: This analysis is for educational purposes only and should not be used as financial advice. Past performance does not guarantee future results.
