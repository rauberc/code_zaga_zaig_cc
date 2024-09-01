# Clear all existing objects from the workspace
rm(list = ls())

# Load necessary libraries
library(gamlss)      # For generalized additive models for location, scale, and shape
library(gamlss.dist) # For distributions used in gamlss
library(KSgeneral)   # For Kolmogorov-Smirnov tests

# Load the data from a text file and convert the date column to Date type
data <- read.table("rain_2015.txt", header = TRUE)
data$date <- as.Date(data$date)

# Convert the vector Z into a time series object with a frequency of 12 (monthly data)
Z_ts <- ts(data$rain, frequency = 365)

# Calculate the autocorrelation function (ACF) of the time series without plotting it
acf_result <- acf(Z_ts, plot = FALSE, lag.max = 19)  # Analyze autocorrelation up to lag 19
acf_result$acf  # Display the calculated ACF values

# Set the file for output and specify plot parameters
postscript(file = "acf_app1.eps", family = "Times", horizontal=T, paper="letter")
{
par(mar=c(5.5, 5.5, 2.5, 2.5)) # Margins: c(bottom, left, top, right)
par(mgp=c(3, 1, 0))  # Axis label positions: c(axis title, axis labels, axis line)

# Plot the autocorrelation function (ACF) with 99% confidence intervals 
plot(acf_result, ci = 0.99, main = "", cex.axis = 1.7, cex.lab = 1.7)
}
# Close the graphics device 
dev.off()

# Calculate the number of rows to use (half of the filtered dataset)
n = floor(nrow(data)/2)

# Extract the first half of the rainfall data
Z <- data$rain[1:n]

# Plot the extracted data
plot(Z)

# Summary statistics and standard deviation of the extracted data
summary(Z)
sd(Z)

# Count and proportion of zero values in the dataset
zeros <- length(Z[Z == 0])
zeros # Display the number of zeros
prop_zeros <- zeros/nrow(data)
prop_zeros # Display the proportion of zeros

# Function to handle errors while fitting models
myoptim <- function(...) tryCatch(expr = gamlss(...), error = function(e) NA)

# Fit the Zero-Adjusted Gamma (ZAGA) and Zero-Adjusted Inverse Gaussian (ZAIG) models
fit1 <- myoptim(Z~1, family=ZAGA, trace = F)
fit2 <- myoptim(Z~1, family=ZAIG, trace = F)

# Extract the estimated parameters from the fitted ZAGA model
mu_ZAGA <- fit1$mu.fv[1]
sigma_ZAGA <- fit1$sigma.fv[1]
nu_ZAGA <- fit1$nu.fv[1]

# Extract the estimated parameters from the fitted ZAIG model
mu_ZAIG <- fit2$mu.fv[1]
sigma_ZAIG <- fit2$sigma.fv[1]
nu_ZAIG <- fit2$nu.fv[1]

# Define the cumulative distribution function (CDF) for the ZAGA distribution
ZAGA_cdf <- function(x)
{
  result <- 0
  if (x < 0){
    result <- 0
  }
  else if (x == 0){
    result <- nu_ZAGA
  }
  else{
    result <- nu_ZAGA + (1 - nu_ZAGA) * pgamma(x, shape = 1/sigma_ZAGA^2, scale = mu_ZAGA * sigma_ZAGA^2)
  }
  return(result)
}

# Define the cumulative distribution function (CDF) for the ZAIG distribution
ZAIG_cdf <- function(x)
{
  result <- 0
  if (x < 0){
    result <- 0
  }
  else if (x == 0){
    result <- nu_ZAIG
  }
  else{
    cdf1 <- pnorm(((x/mu_ZAIG) - 1)/(sigma_ZAIG * sqrt(x)))
    cdf2 <- exp(2/(mu_ZAIG * sigma_ZAIG^2)) * pnorm((-((x/mu_ZAIG) + 1))/(sigma_ZAIG * sqrt(x)))
    cdf <- cdf1 + cdf2
    result <- nu_ZAIG + (1 - nu_ZAIG) * cdf
  }
  return(result)
}

# Define the range for variable of interest
points <- c(min(Z), max(Z))

# Calculate Akaike Information Criterion (AIC) and Schwarz Bayesian Criterion (SBC)
AIC_zaga <- fit1$aic 
AIC_zaga
SBC_zaga <- fit1$sbc 
SBC_zaga

AIC_zaig <- fit2$aic 
AIC_zaig
SBC_zaig <- fit2$sbc 
SBC_zaig

# Performing goodness-of-fit tests using the cumulative distribution functions
KSgeneral::mixed_ks_test(Z, points, Mixed_dist = ZAGA_cdf) 
KSgeneral::mixed_ks_test(Z, points, Mixed_dist = ZAIG_cdf)

# Setting the probability of false alarm
alpha <- 1/370

# For ZAGA: Lower limit is 0, calculate the upper control limit
LCL_zaga <- 0
UCL_zaga <- qZAGA(1-alpha, mu_ZAGA, sigma_ZAGA, nu_ZAGA)
round(UCL_zaga,2)

# For ZAIG: Lower limit is 0, calculate the upper control limit
LCL_zaig <- 0
UCL_zaig <- qZAIG(1-alpha, mu_ZAIG, sigma_ZAIG, nu_ZAIG)
round(UCL_zaig,2)

# Identify outliers that exceed the ZAGA upper control limit
out <- which(data$rain > UCL_zaga)
out_zaga <- data$rain[out]

# Create a control chart for the rainfall data

# Set the file for output and specify plot parameters
postscript(file = "cc_app1.eps", family = "Times", horizontal=T, paper="letter")
{
par(mar=c(5.5, 5.5, 2.5, 2.5)) # Margins: c(bottom, left, top, right)
par(mgp=c(3, 1, 0))

# Plot the rainfall data excluding outliers, with specified labels and axis limits
plot(data$rain[-out], xaxs = "r", pch = "+", ylab = expression("Rainfall"), 
     xlab = expression("Observations"), ylim = c(0, UCL_zaig+50), cex.axis = 1.7, cex.lab = 1.7)

# Add the outliers to the plot in red
points(out, out_zaga, pch = "+", col = "red")

# Add horizontal lines representing the upper control limits for ZAGA and ZAIG control charts
abline(h = UCL_zaga, lty = 1) # Solid line for ZAGA
abline(h = UCL_zaig, lty = 2) # Dashed line for ZAIG
abline(h = 0, col = "black")

# Add a vertical line to divide the plot into Phase I and Phase II
abline(v=n, lty=3, lwd=2 ,col="black") 

# Add labels for Phase I and Phase II
text(70, UCL_zaig+20, labels="Phase I", cex=1.7, col=1)
text(n+100, UCL_zaig+20, labels="Phase II", cex=1.7, col=1)
}
# Close the graphics device
dev.off()

# Plot the empirical cumulative distribution function (ECDF) for the ZAGA distribution
ZAGA_cdf <- function(x, nu_ZAGA, mu_ZAGA, sigma_ZAGA) 
{
  result=ifelse(x==0, nu_ZAGA, 
                nu_ZAGA + (1 - nu_ZAGA) * pgamma(x, shape = 1/sigma_ZAGA^2, 
                          scale = mu_ZAGA * sigma_ZAGA^2))
  return(result)
}

# Set the file for output and specify plot parameters
postscript(file = "ecdf_ZAGA_app1.eps", family = "Times", horizontal=T, paper="letter")
{
# Set up plot margins and axis label positions
par(mar=c(5.5, 5.5, 2.5, 2.5)) # Margins: c(bottom, left, top, right)
par(mgp=c(3, 1, 0)) # Axis label positions: c(axis title, axis labels, axis line)

# Calculate the empirical cumulative distribution function (ECDF) for the data
ecdfZ = stats::ecdf(Z)

# Define the cumulative distribution function (CDF) for the fitted ZAGA distribution
cdf = function(x) { 
  ZAGA_cdf(x, nu_ZAGA = nu_ZAGA, mu_ZAGA = mu_ZAGA, sigma_ZAGA = sigma_ZAGA) 
}

# Plot the ECDF (empirical cumulative distribution) of the data
curve(ecdfZ, from = min(Z), to = max(Z), xlim = c(min(Z), max(Z)), ylim = c(0,1), 
      xlab = "Rainfall", ylab = "Cumulative distribution", 
      lty = 1, # Solid line for the empirical CDF
      main = "", cex.lab = 1.7, cex.axis = 1.7) # Labels and axis font size

# Add the fitted ZAGA CDF to the plot with a dashed line
curve(cdf, from = min(Z), to = max(Z), lty = 2, cex.lab = 1.7, cex.axis = 1.7, add = TRUE)
}
# Close the graphics device 
dev.off()

# ARL performance

# Perturbation vector with values ranging from 1.0 to 1.5 in increments of 0.1
perturbation_vector <- c(seq(1.0, 1.5, 0.1)) # 1 means in control
n_perturbation_vector <- length(perturbation_vector)

# Initialize vectors to store the theoretical alpha values for ZAGA and ZAIG
alpha_theoretical_zaga = alpha_theoretical_zaig <- vector()

# Calculate the theoretical alpha values 
for(i in 1:n_perturbation_vector){
  mu <- mu_ZAGA * (perturbation_vector[i]) # Multiply the mean by the perturbation factor

  # Calculate the probability for lower and upper limits in ZAGA
  a <- pZAGA(LCL_zaga, mu, sigma_ZAGA, nu_ZAGA) 
  b <- pZAGA(UCL_zaga, mu, sigma_ZAGA, nu_ZAGA) 
  
  # Append the unilateral alpha value (left side) to the alpha theoretical vector for ZAGA
  alpha_theoretical_zaga <- c(alpha_theoretical_zaga, (1-(b))) 
  
  # Repeat the calculations for ZAIG 
  a1 <- pZAGA(LCL_zaig, mu, sigma_ZAGA, nu_ZAGA) 
  b1 <- pZAGA(UCL_zaig, mu, sigma_ZAGA, nu_ZAGA)
  
  alpha_theoretical_zaig <- c(alpha_theoretical_zaig, (1-(b1)))
}

# Combine the control limits for ZAGA and ZAIG into a matrix for comparison
results <- matrix(rbind(LCL_zaga, UCL_zaga, LCL_zaig, UCL_zaig), 
                           nrow=2, ncol=2)
colnames(results) <- c("Zaga","Zaig")
rownames(results) <- c("Lower Limit","Upper Limit")

# Output the results matrix to a text file
write.table(round(results, digits = 4), 
            file="results_control_limits_app1.txt", sep=" ")

# Calculate various performance metrics based on the theoretical alpha values
resulting_alpha_theoretical <- matrix(rbind(alpha_theoretical_zaga, 
                                            1/alpha_theoretical_zaga, 
                                            sqrt((1-alpha_theoretical_zaga)/alpha_theoretical_zaga^2), 
                                            log(0.5)/log(1-alpha_theoretical_zaga), 
                                            alpha_theoretical_zaig, 
                                            1/alpha_theoretical_zaig, 
                                            sqrt((1-alpha_theoretical_zaig)/alpha_theoretical_zaig^2), 
                                            log(0.5)/log(1-alpha_theoretical_zaig)),
                                      nrow=8, ncol=n_perturbation_vector)

# Name the columns by the perturbation vector values
colnames(resulting_alpha_theoretical) <- c(perturbation_vector)

# Name the rows to identify the different metrics for ZAGA and ZAIG
rownames(resulting_alpha_theoretical) <- c("alpha Zaga","ARL1 Zaga", "SDRL1 Zaga", "MRL1 Zaga", 
                                           "alpha Zaig","ARL1 Zaig", "SDRL1 Zaig", "MRL1 Zaig")

# Print the resulting theoretical alpha values rounded to 4 decimal places
print(round(resulting_alpha_theoretical, digits = 2))

# Output the resulting theoretical alpha values matrix to a text file
write.table(round(resulting_alpha_theoretical, digits = 2), 
            file="results_ARL_app1.txt", sep=" ")


