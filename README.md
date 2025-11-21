# M-quantile Spectral Discriminant Analysis (`mqper` and related functions)

A collection of R functions for performing robust discriminant analysis in the frequency domain using M-quantile regression methods.

---

## Description
This suite of functions implements M-quantile spectral analysis and discriminant classification for time series data. The main components are:

1. **`mqper`**: Computes M-quantile periodograms for robust spectral estimation
2. **`lperd.mtm`**: Smooths and logs the periodogram
3. **`cep.mtm`**: Computes cepstral coefficients from the log-periodogram
4. **`cep.lda`**: Main function performing discriminant analysis using cepstral features

Key features:
- Robust spectral estimation using M-quantile regression
- Cepstral coefficient extraction for feature representation
- Linear discriminant analysis (LDA) for classification
- Automatic selection of optimal number of cepstral coefficients
- Cross-validation support

---

## Installation
Ensure the following packages are installed:
```R
install.packages(c("MASS", "stats"))

source('mquantile.R') # to estimate the M-quantile regression periodogram
source("mquantile.rlda.R")

# Generate sample time series data
set.seed(123)
x <- rnorm(256)  # 256-point time series

# Compute M-quantile periodogram (median)
spec_result <- mqper(x, tau = 0.5)

# Plot spectrum
plot(spec_result$freq, spec_result$spec, type = "l", 
     main = "M-Quantile Periodogram", xlab = "Frequency", ylab = "Spectrum")

# Compute cepstral coefficients
cep_result <- cep.mtm(x, tau = 0.5)

# Plot cepstral coefficients
plot(cep_result$quef, cep_result$cep, type = "h",
     main = "Cepstral Coefficients", xlab = "Quefrency", ylab = "Amplitude")

# Create sample data: 3 classes, 50 time series each
y <- rep(1:3, each = 50)
x <- matrix(rnorm(256*150), nrow = 256, ncol = 150)

# Run discriminant analysis
result <- cep.lda(y, x, tau = 0.5)

# View confusion matrix
table(y, result$predict$class)

# Get optimal number of cepstral coefficients
print(result$Lopt)
