# =========================================================
# 0. CLEAR ENVIRONMENT AND LOAD REQUIRED PACKAGES
# =========================================================

rm()
graphics.off()

# Install required packages
install.packages("AER")
install.packages("plm")
install.packages("lmtest")
install.packages("sandwich")
install.packages("stargazer")
install.packages("AER")   # Install AER package
install.packages("gam")                     # Install GAM package

# Load required libraries
library(splines)                            # Load spline functions
library(gam)                                # Load GAM functions
library(AER) 
library(stargazer) 
library(plm)
library(ggplot2)
library(quantreg)
library(lmtest)
library(sandwich)


# =========================================================
# 1. DEFINE QUANTILE GRID
# =========================================================

tau_min  <- 0.05
tau_max  <- 0.95
tau_step <- 0.05
taus <- seq(tau_min, tau_max, by = tau_step)


# =========================================================
# 2. SET WORKING DIRECTORY AND LOAD DATA
# =========================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the dataset
df <- read.csv("data1graph.csv")

# Display variable names
names(df)


# =========================================================
# 3. ROBUSTNESS SAMPLE BASED ON MATURITY
# =========================================================

# Define lower and upper maturity thresholds
x <- 0.3
y <- 0.7

# Compute maturity quantile cutoffs
cut_x <- quantile(df$years_to_maturity, probs = x, na.rm = TRUE)
cut_y <- quantile(df$years_to_maturity, probs = 1 - y, na.rm = TRUE)

# Create low- and high-maturity subsamples
df_low_x  <- df[df$years_to_maturity <= cut_x, ]
df_high_y <- df[df$years_to_maturity >= cut_y, ]

# Optional inspection of maturity subsamples
#df_low_x
#df_high_y

# Select the high-maturity subsample for subsequent analysis
df<-df_high_y


# =========================================================
# 4. OLS AND QUANTILE REGRESSION PLOT
# =========================================================

# Convert auction date to Date format
df$auction_date <- as.Date(df$auction_date, format = "%m/%d/%Y")

# Display variable names
names (df)

# Estimate OLS and quantile regression models
mod_ols <- lm(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, data = df)
mod_q25 <- rq(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, tau = 0.25, data = df)
mod_q50 <- rq(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, tau = 0.50, data = df)
mod_q75 <- rq(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, tau = 0.75, data = df)

# Set plot margins and use a single plotting panel
par(mar = c(4,4,2,1))
par(mfrow = c(1,1))

# Plot auction premia against issuance shocks
plot(df$shock_log, df$auction_premium_bp,
     xlab = "shock_log",
     ylab = "auction_premium_bp",
     main = "OLS and quantile regression",
     pch = 1)

# Create a prediction grid for fitted regression lines
x_grid <- seq(min(df$shock_log, na.rm = TRUE),
              max(df$shock_log, na.rm = TRUE),
              length.out = 200)

# Refine the prediction grid for smoother fitted lines
#######################
x_grid <- seq(min(df$shock_log, na.rm = TRUE),
              max(df$shock_log, na.rm = TRUE),
              length.out = 300)

# Construct prediction data while holding controls at their sample means
new_df <- data.frame(
  shock_log = x_grid,
  Vix_euro = mean(df$Vix_euro, na.rm = TRUE),
  outstand_notional = mean(df$outstand_notional, na.rm = TRUE),
  auction_date = rep(mean(df$auction_date, na.rm = TRUE), length(x_grid))
)

# Generate fitted values from quantile and OLS models
y_q25 <- predict(mod_q25, newdata = new_df)
y_q75 <- predict(mod_q75, newdata = new_df)
y_ols <- predict(mod_ols, newdata = new_df)

# Add fitted quantile and OLS regression lines
lines(x_grid, y_q25, lwd = 2, lty = 2, col = "blue")
lines(x_grid, y_q75, lwd = 2, lty = 2, col = "red")
lines(x_grid, y_ols, lwd = 2, lty = 1, col = "black")

# Add plot legend
legend("topleft",
       legend = c("ols", "q25", "q75"),
       col = c("black", "blue", "red"),
       lwd = 2,
       lty = c(1, 2, 2),
       bty = "n")


# =========================================================
# 5. RESIDUAL DIAGNOSTIC PLOT
# =========================================================

# Extract OLS residuals and fitted values
res_lin <- residuals(mod_ols)
fit_lin <- fitted(mod_ols)

# Plot residuals against fitted values
plot(fit_lin, res_lin,
     xlab = "fitted values",
     ylab = "residuals",
     main = "residuals vs fitted - linear model",
     pch = 19)

# Add horizontal zero-reference line
abline(h = 0, lty = 2)


# =========================================================
# 6. HETEROSKEDASTICITY TEST AND ROBUST INFERENCE
# =========================================================

# Perform the Breusch-Pagan test
bptest(mod_ols)

# Perform the studentized Breusch-Pagan test
bptest(mod_ols, studentize = TRUE)

# Compute heteroskedasticity-robust standard errors
coeftest(mod_ols, vcov = vcovHC(mod_ols, type = "HC1"))


# =========================================================
# 7. LINEAR, QUADRATIC, OLS, AND QUANTILE REGRESSIONS
# =========================================================

# Estimate linear quantile regression models
mod_q25 <- rq(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, tau = 0.25, data = df)
mod_q75 <- rq(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, tau = 0.75, data = df)

# Estimate OLS specifications with and without a time control
mod_ols <- lm(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional, data = df)
mod_olsd <- lm(auction_premium_bp ~ shock_log+Vix_euro+outstand_notional+auction_date, data = df)

# Estimate quadratic OLS and quantile regression specifications
mod_olsd2 <- lm(auction_premium_bp ~ shock_log+I(shock_log^2)+Vix_euro+outstand_notional+auction_date, data = df)
mod2_q25 <- rq(auction_premium_bp ~ shock_log+I(shock_log^2)+Vix_euro+outstand_notional+auction_date, tau = 0.25, data = df)
mod2_q75 <- rq(auction_premium_bp ~ shock_log+I(shock_log^2)+Vix_euro+outstand_notional+auction_date, tau = 0.75, data = df)

# Compute HC1 robust standard errors for OLS specifications
mod_olsr<-coeftest(mod_ols, vcov = vcovHC(mod_ols, type = "HC1"))
mod_olsdr<-coeftest(mod_olsd, vcov = vcovHC(mod_olsd, type = "HC1"))
mod_olsd2r<-coeftest(mod_olsd2, vcov = vcovHC(mod_olsd2, type = "HC1"))

# Display regression results
stargazer(mod_q25, mod_q75, mod_olsr, mod_olsdr, mod_olsd2r,  type = "text")

# Display additional model results
stargazer(mod_amount, type="text")

# Extract R-squared and adjusted R-squared statistics
summary(mod_ols)$r.squared
summary(mod_ols)$adj.r.squared
summary(mod_olsd)$r.squared
summary(mod_olsd)$adj.r.squared
summary(mod_olsd2)$r.squared
summary(mod_olsd2)$adj.r.squared


# =========================================================
# 8. SELECT THE BEST POLYNOMIAL DEGREE USING LOOCV
# =========================================================

# Load cross-validation package
library(boot)

# Create a complete-case dataset for model selection
# =========================
df2 <- df[, c("auction_premium_bp", "shock_log", "Vix_euro", "outstand_notional", "auction_date")]
df2 <- na.omit(df2)

# Set seed for reproducibility
set.seed(14)

# Use one observation per fold for leave-one-out cross-validation
k <- nrow(df2)

# Set the maximum polynomial degree to evaluate
max_deg <- 9  # Change this value to test a different maximum degree

# Initialize vector to store cross-validation MSE values
mse_cv <- numeric(max_deg)

# Display variable names
names(df)

# Estimate polynomial models and compute LOOCV errors
#####################
for(d in 1:max_deg){
  
  # Construct polynomial regression formula for degree d
  form <- as.formula(
    paste0("auction_premium_bp ~ poly(shock_log, ", d, ", raw=TRUE) + Vix_euro + outstand_notional+auction_date")
  )
  
  # Estimate the polynomial model
  mod <- glm(form, data = df2)
  
  # Perform leave-one-out cross-validation
  cv_res <- cv.glm(df2, mod, K = k)
  
  # Store the corresponding cross-validation MSE
  mse_cv[d] <- cv_res$delta[1]
}

# Identify the polynomial degree with the lowest MSE
best_deg <- which.min(mse_cv)

# Print optimal polynomial degree and minimum MSE
cat("meilleur degré =", best_deg, "\n")
cat("mse minimum =", mse_cv[best_deg], "\n")

# Plot LOOCV MSE across polynomial degrees
plot(1:max_deg, mse_cv,
     type = "b",
     pch = 19,
     xlab = "Polynomial degree of shock_log",
     ylab = "LOOCV MSE",
     main = "Leave-one-out CV for polynomial")

# Mark the selected polynomial degree
abline(v = best_deg, lty = 2)


# =========================================================
# 9. QUADRATIC QUANTILE AND OLS SPECIFICATIONS
# =========================================================

# Estimate quadratic quantile regressions at the 25th and 75th quantiles
############## POLYNONIAL 2

poly25 <- rq(auction_premium_bp ~ poly(shock_log,  2, raw=TRUE) 
             + Vix_euro + outstand_notional+auction_date, tau = 0.25, data = df2)
poly75 <- rq(auction_premium_bp ~ poly(shock_log,  2, raw=TRUE) 
             + Vix_euro + outstand_notional+auction_date, tau = 0.75, data = df2)

# Estimate the corresponding quadratic OLS model
mod_ols <- lm(auction_premium_bp ~ poly(shock_log,  2, raw=TRUE) +Vix_euro+outstand_notional+auction_date, data = df)

# Compute HC1 robust standard errors
mod_olsr<-coeftest(mod_ols, vcov = vcovHC(mod_ols, type = "HC1"))

# Display quadratic regression results
stargazer(poly25, poly75, mod_olsr, type = "text")


# =========================================================
# 10. COEFFICIENT EVOLUTION ACROSS QUANTILES
# =========================================================

# Initialize vectors for coefficient estimates and uncertainty measures
####### PLOT IMAGE DISTRIBUTION IMAGE
beta_vec <- numeric(length(taus))
low_vec  <- numeric(length(taus))
high_vec <- numeric(length(taus))

# Initialize storage vectors for the quantile loop
eta_vec <- numeric(length(taus))
low_vec  <- numeric(length(taus))
high_vec <- numeric(length(taus))

# Ensure auction dates use Date format
df$auction_date <- as.Date(df$auction_date, format = "%m/%d/%Y")

# Estimate the quadratic quantile model across the full quantile grid
for (i in seq_along(taus)) {
  
  # Select the current quantile
  tau_i <- taus[i]
  
  # Estimate the quantile regression model
  model_i <- rq(auction_premium_bp ~ shock_log+I(shock_log^2)+Vix_euro+outstand_notional+auction_date, tau = tau_i, data = df)
  
  # Extract model summary
  sum_i   <- summary(model_i)
  
  # Extract coefficient matrix
  coef_mat <- sum_i$coefficients
  
  # Locate the linear shock coefficient
  row_id   <- which(rownames(coef_mat) == "shock_log")
  
  # Store coefficient estimates and associated statistics
  beta_vec[i] <- coef_mat[row_id, 1]
  low_vec[i]  <- coef_mat[row_id, 2]
  high_vec[i] <- coef_mat[row_id, 3]
}

# Combine quantile-specific estimates into a results table
results_qr <- data.frame(
  tau   = taus,
  beta  = beta_vec,
  lower = low_vec,
  upper = high_vec
)

# Display the final quantile results table
print(results_qr)

# Plot the evolution of the linear shock coefficient across quantiles
plot(results_qr$tau, results_qr$beta,
     type = "b",
     pch = 19,
     lwd = 2,
     xlab = "quantile",
     ylab = "coefficient on linear term")