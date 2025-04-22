# Clear all existing objects from the workspace
rm(list = ls())

# Load necessary libraries
library(gamlss)
library(gamlss.dist)
library(KSgeneral)
library(tseries)

# Load the data from a text file and convert the date column to Date type
data <- read.table("rain_2015.txt", header = TRUE)
data$date <- as.Date(data$date)

# Calculate the number of rows to use (half of the dataset)
n <- floor(nrow(data)/2)

# Extract the first half of the data
Z <- data$rain[1:n]

# Convert the vector Z into a time series object with a frequency of 365 (daily data)
Z_ts <- ts(Z, frequency = 365)

# Calculate the autocorrelation function (ACF) of the time series without plotting it
acf_result <- acf(Z_ts, plot = FALSE, lag.max = 19)

# Set the file for output and specify plot parameters
postscript(file = "acf_app.eps", family = "Times", horizontal=TRUE, paper="letter")
par(mar=c(5.5, 5.5, 2.5, 2.5), mgp=c(3, 1, 0))
plot(acf_result, ci = 0.99, main = "", cex.axis = 1.7, cex.lab = 1.7)
dev.off()

# Summary statistics and standard deviation of the data
summary(Z)
sd(Z)

# Count and proportion of zero values in the dataset
zeros <- sum(Z == 0)
prop_zeros <- zeros / length(Z)
cat("Proportion of zeros:", round(prop_zeros * 100, 2), "%\n")

# Fit the Zero-Adjusted Gamma (ZAGA) and Zero-Adjusted Inverse Gaussian (ZAIG) models
fit1 <- gamlss(Z ~ 1, family = ZAGA)
fit2 <- gamlss(Z ~ 1, family = ZAIG)

# Extract the estimated parameters from the fitted ZAGA model
mu_ZAGA <- fit1$mu.fv[1]; sigma_ZAGA <- fit1$sigma.fv[1]; nu_ZAGA <- fit1$nu.fv[1]

# Extract the estimated parameters from the fitted ZAIG model
mu_ZAIG <- fit2$mu.fv[1]; sigma_ZAIG <- fit2$sigma.fv[1]; nu_ZAIG <- fit2$nu.fv[1]

# Define the cumulative distribution function (CDF) for the ZAGA distribution
ZAGA_cdf <- function(x, mu = mu_ZAGA, sigma = sigma_ZAGA, nu = nu_ZAGA) {
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi) || xi < 0) {
      out[i] <- 0
    } else if (xi == 0) {
      out[i] <- nu
    } else {
      shape <- 1 / sigma^2
      scale <- mu * sigma^2
      if (is.finite(shape) && is.finite(scale) && shape > 0 && scale > 0) {
        out[i] <- nu + (1 - nu) * pgamma(xi, shape = shape, scale = scale)
      } else {
        out[i] <- NA
      }
    }
  }
  return(out)
}

# Define the cumulative distribution function (CDF) for the ZAIG distribution
ZAIG_cdf <- function(x, mu = mu_ZAIG, sigma = sigma_ZAIG, nu = nu_ZAIG) {
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    xi <- x[i]
    if (is.na(xi) || xi < 0) {
      out[i] <- 0
    } else if (xi == 0) {
      out[i] <- nu
    } else {
      denom <- sigma * sqrt(xi)
      if (is.finite(denom) && denom > 0) {
        z1 <- ((xi / mu) - 1) / denom
        z2 <- -((xi / mu) + 1) / denom
        cdf1 <- pnorm(z1)
        cdf2 <- exp(2 / (mu * sigma^2)) * pnorm(z2)
        out[i] <- nu + (1 - nu) * (cdf1 + cdf2)
      } else {
        out[i] <- NA
      }
    }
  }
  return(out)
}

# Calculate Akaike Information Criterion (AIC) and Schwarz Bayesian Criterion (SBC)
AIC_zaga <- fit1$aic; SBC_zaga <- fit1$sbc
AIC_zaig <- fit2$aic; SBC_zaig <- fit2$sbc
cat("AIC ZAGA:", AIC_zaga, " | SBC ZAGA:", SBC_zaga, "\n")
cat("AIC ZAIG:", AIC_zaig, " | SBC ZAIG:", SBC_zaig, "\n")

# Performing goodness-of-fit tests using the cumulative distribution functions
Z <- Z[!is.na(Z) & is.finite(Z)]
ks_zaga <- KSgeneral::mixed_ks_test(Z, range(Z), Mixed_dist = ZAGA_cdf)
ks_zaig <- KSgeneral::mixed_ks_test(Z, range(Z), Mixed_dist = function(x) ZAIG_cdf(x))
print(ks_zaga)
print(ks_zaig)

# Setting the probability of false alarm
alpha <- 1/370

# For ZAGA: Lower limit is 0, calculate the upper control limit
LCL_zaga <- 0
UCL_zaga <- qZAGA(1 - alpha, mu_ZAGA, sigma_ZAGA, nu_ZAGA)

# For ZAIG: Lower limit is 0, calculate the upper control limit
LCL_zaig <- 0
UCL_zaig <- qZAIG(1 - alpha, mu_ZAIG, sigma_ZAIG, nu_ZAIG)

# Outliers ZAGA e ZAIG
out_zaga <- which(data$rain > UCL_zaga)
outliers_zaga <- data$rain[out_zaga]

# Identify outliers that exceed the ZAGA upper control limit
out_zaig <- which(data$rain > UCL_zaig)
outliers_zaig <- data$rain[out_zaig]

cat("ZAGA:", length(out_zaga), "observations out of control\n")
cat("ZAIG:", length(out_zaig), "observations out of control\n")

# Create a control chart for the rainfall data

# Set the file for output and specify plot parameters
postscript(file = "cc_app.eps", family = "Times", horizontal=TRUE, paper="letter")
par(mar=c(5.5, 5.5, 2.5, 2.5), mgp=c(3, 1, 0))
plot(data$rain, xaxs = "r", pch = "+", ylab = "Rainfall", 
     xlab = "Observations", ylim = c(0, max(data$rain)+10), 
     cex.axis = 1.7, cex.lab = 1.7)
points(out_zaga, outliers_zaga, pch = "+", col = "red")
points(out_zaig, outliers_zaig, pch = "+", col = "blue")
abline(h = UCL_zaga, lty = 1)
abline(h = UCL_zaig, lty = 2)
abline(h = 0, col = "black")
abline(v = n, lty = 3, lwd = 2 ,col = "black")
text(70, max(data$rain)+5, labels="Phase I", cex=1.7)
text(n + 60, max(data$rain)+5, labels="Phase II", cex=1.7)
dev.off()

# ECDF + CDF ZAGA
postscript(file = "ecdf_ZAGA_app.eps", family = "Times", horizontal=TRUE, paper="letter")
par(mar=c(5.5, 5.5, 2.5, 2.5), mgp=c(3, 1, 0))

ecdfZ <- ecdf(Z)
curve(ecdfZ, from = min(Z), to = max(Z), xlab = "Rainfall", ylab = "Cumulative distribution",
      ylim = c(0, 1), lty = 1, cex.lab = 1.7, cex.axis = 1.7)
curve(ZAGA_cdf(x), from = min(Z), to = max(Z), lty = 2, add = TRUE)
dev.off()

# Theoretical ARL 
deltas <- seq(1.0, 1.5, 0.1)
ARL_zaga <- ARL_zaig <- SDRL_zaga <- SDRL_zaig <- MRL_zaga <- MRL_zaig <- numeric(length(deltas))

for (i in seq_along(deltas)) {
  mu_shifted <- mu_ZAGA * deltas[i]
  alpha_zaga <- 1 - pZAGA(UCL_zaga, mu = mu_shifted, sigma = sigma_ZAGA, nu = nu_ZAGA)
  alpha_zaig <- 1 - pZAIG(UCL_zaig, mu = mu_shifted, sigma = sigma_ZAIG, nu = nu_ZAIG)
  
  ARL_zaga[i] <- 1 / alpha_zaga
  SDRL_zaga[i] <- sqrt((1 - alpha_zaga) / alpha_zaga^2)
  MRL_zaga[i] <- log(0.5) / log(1 - alpha_zaga)
  
  ARL_zaig[i] <- 1 / alpha_zaig
  SDRL_zaig[i] <- sqrt((1 - alpha_zaig) / alpha_zaig^2)
  MRL_zaig[i] <- log(0.5) / log(1 - alpha_zaig)
}

table_arl <- data.frame(
  delta = deltas,
  ARL_ZAGA = round(ARL_zaga, 2),
  SDRL_ZAGA = round(SDRL_zaga, 2),
  MRL_ZAGA = round(MRL_zaga, 2),
  ARL_ZAIG = round(ARL_zaig, 2),
  SDRL_ZAIG = round(SDRL_zaig, 2),
  MRL_ZAIG = round(MRL_zaig, 2)
)

print(table_arl)
write.table(table_arl, file = "table_app_ARL.txt", sep = "\t", row.names = FALSE)

