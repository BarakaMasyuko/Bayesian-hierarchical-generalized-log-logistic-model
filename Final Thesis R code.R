# Load Required Libraries
library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)

# Set parallel computing for faster sampling
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Load the cleaned dataset
child_survival_data_final <- read.csv("C:/Users/BarakaNgumbauMasyuko/OneDrive - Girl Child Network/Desktop/final_cleaned_data.csv")

# Step 1: Check and Fix Survival Time Variable (Ensure No Zeros)
child_survival_data_final <- child_survival_data_final %>%
  mutate(
    survival_time = log(1 + ifelse(Current.age.of.child.in.months..months.since.birth.for.dead.children. == 0, 
                                   0.1, 
                                   Current.age.of.child.in.months..months.since.birth.for.dead.children.)),
    event_indicator = as.integer(Child.is.alive == "Yes"),  # Convert survival status to binary (1 = alive, 0 = dead)
    region = as.numeric(as.factor(Region)),  
    maternal_education = as.numeric(Highest.educational.level),  
    wealth_index = as.numeric(Wealth.index.combined),  
    iptp_use = as.numeric(During.pregnancy.took..SP.fansidar.for.malaria),  
    malaria_endemicity = as.numeric(Malaria.endemicity.zone),  
    antenatal_visits = as.numeric(Number.of.antenatal.visits.during.pregnancy),  
    timing_first_antenatal = as.numeric(Timing.of.1st.antenatal.check..months.),  
    mosquito_net_use = as.numeric(Children.under.5.slept.under.mosquito.bed.net.last.night..household.questionnaire.),  
    anemia_level = as.numeric(Anemia.level),  
    malaria_history = as.numeric(Told.child.had.malaria)  # Convert malaria history to numeric (1 = Yes, 0 = No)
  )



####################################################################Kaplan Meier
# Bin wealth index into tertiles for clearer visualization
child_survival_data_final <- child_survival_data_final %>%
  mutate(wealth_tertile = ntile(wealth_index, 3))

km_wealth <- survfit(
  Surv(Current.age.of.child.in.months..months.since.birth.for.dead.children., 
       1 - event_indicator) ~ wealth_tertile, 
  data = child_survival_data_final
)

ggsurvplot(
  km_wealth,
  data = child_survival_data_final,
  conf.int = TRUE,
  pval = TRUE,
  legend.title = "Wealth Index",
  legend.labs = c("Low (T1)", "Middle (T2)", "High (T3)"),
  palette = "lancet",  # Lancet journal colors
  title = "Survival by Wealth Tertile"
)

############################################################################
# Recode education levels if needed
child_survival_data_final <- child_survival_data_final %>%
  mutate(edu_cat = case_when(
    maternal_education <= 2 ~ "Primary or less",
    maternal_education == 3 ~ "Secondary",
    maternal_education >= 4 ~ "Higher"
  ))

km_edu <- survfit(
  Surv(Current.age.of.child.in.months..months.since.birth.for.dead.children.,
       1 - event_indicator) ~ edu_cat,
  data = child_survival_data_final
)

ggsurvplot(
  km_edu,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "Maternal Education",
  palette = "aaas",  # Science journal colors
  title = "Survival by Maternal Education Level",
  font.legend = 10  # Adjust legend font size
)
##################################################################################
# Load the survival package
library(survival)

# Create survival object
surv_obj <- Surv(time = child_survival_data_final$Current.age.of.child.in.months..months.since.birth.for.dead.children.,
                 event = 1 - child_survival_data_final$event_indicator)  # Note: event=1 for death, 0 for censored

# Fit Kaplan-Meier estimator
km_fit <- survfit(surv_obj ~ 1)

# Plot Kaplan-Meier curve
plot(km_fit, 
     main = "Kaplan-Meier Survival Curve",
     xlab = "Time (months)", 
     ylab = "Survival Probability",
     col = "blue",
     conf.int = TRUE)

# Add grid lines for better readability
grid()

# Add legend
legend("topright", 
       legend = c("KM Estimate", "95% CI"),
       col = c("blue", "blue"),
       lty = c(1, 2),
       bty = "n")


#############################################################################
# Example: KM curve by malaria history
surv_obj_group <- Surv(time = child_survival_data_final$Current.age.of.child.in.months..months.since.birth.for.dead.children.,
                       event = 1 - child_survival_data_final$event_indicator)

km_fit_group <- survfit(surv_obj_group ~ child_survival_data_final$malaria_history)

ggsurvplot(km_fit_group,
           data = child_survival_data_final,
           risk.table = TRUE,
           pval = TRUE,  # Add log-rank test p-value
           conf.int = TRUE,
           xlab = "Time (months)",
           ylab = "Survival Probability",
           title = "KM Curve by Malaria History",
           legend.title = "Malaria History",
           legend.labs = c("No", "Yes"),
           ggtheme = theme_minimal(),
           palette = c("#E7B800", "#2E9FDF"))
#############################################################################

# Then re-run KM analysis
# KM by Region (using the 'Region' column)
km_region <- survfit(
  Surv(Current.age.of.child.in.months..months.since.birth.for.dead.children., 
       1 - event_indicator) ~ Region, 
  data = child_survival_data_final
)

# Plot with survminer
ggsurvplot(
  km_region,
  data = child_survival_data_final,
  risk.table = TRUE,
  pval = TRUE,               # Log-rank test for group differences
  conf.int = TRUE,           # Show confidence intervals
  xlab = "Time (months)",
  ylab = "Survival Probability",
  title = "Child Survival by Region",
  legend.title = "Region",
  legend.labs = levels(factor(child_survival_data_final$Region)),  # Auto-generate labels
  ggtheme = theme_minimal(),
  palette = "jco",           # Journal of Clinical Oncology color palette
  risk.table.height = 0.25   # Adjust risk table size
)


#########################################################################
# KM by Endemicity (using 'Malaria.endemicity.zone')
km_endemic <- survfit(
  Surv(Current.age.of.child.in.months..months.since.birth.for.dead.children., 
       1 - event_indicator) ~ Malaria.endemicity.zone, 
  data = child_survival_data_final
)

# Plot with customizations
ggsurvplot(
  km_endemic,
  data = child_survival_data_final,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Survival Probability",
  title = "Child Survival by Malaria Endemicity Zone",
  legend.title = "Endemicity Zone",
  legend.labs = levels(factor(child_survival_data_final$Malaria.endemicity.zone)),
  ggtheme = theme_bw(),       # Clean black-and-white theme
  palette = "npg",           # Nature Publishing Group colors
  break.time.by = 12,        # X-axis breaks every 12 months
  risk.table.y.text = FALSE   # Cleaner risk table
)



############################################################################
# Load libraries
library(survival)
library(ggplot2)

# Fit a null Cox model (no predictors)
cox_fit <- coxph(Surv(Current.age.of.child.in.months..months.since.birth.for.dead.children., 
                      1 - event_indicator) ~ 1, 
                 data = child_survival_data_final)

# Now get baseline hazard
hazard <- basehaz(cox_fit, centered = FALSE)

# Proceed with your log-transform and plotting
hazard$logH <- log(hazard$hazard)
hazard$logT <- log(hazard$time)

ggplot(hazard, aes(x = logT, y = logH)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Log-Cumulative Hazard Plot",
       x = "log(Time in months)",
       y = "log(Cumulative Hazard)") +
  theme_minimal()

#############################################################################
# Prepare Data for Stan
stan_data_list_new <- list(
  N = nrow(child_survival_data_final),
  K = 9,  # Now 9 covariates (added malaria history)
  R = length(unique(child_survival_data_final$region)),  
  survival_time = child_survival_data_final$survival_time,  
  event_indicator = child_survival_data_final$event_indicator,  
  region = child_survival_data_final$region,  
  maternal_education = child_survival_data_final$maternal_education,  
  wealth_index = child_survival_data_final$wealth_index,  
  iptp_use = child_survival_data_final$iptp_use,  
  malaria_endemicity = child_survival_data_final$malaria_endemicity,  
  antenatal_visits = child_survival_data_final$antenatal_visits,  
  timing_first_antenatal = child_survival_data_final$timing_first_antenatal,  
  mosquito_net_use = child_survival_data_final$mosquito_net_use,  
  anemia_level = child_survival_data_final$anemia_level,  
  malaria_history = child_survival_data_final$malaria_history  # New variable
)

# Define the Updated Stan Model with Malaria History
stan_model_code_new <- "
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> R;
  vector[N] survival_time;
  int<lower=0, upper=1> event_indicator[N];
  int<lower=1, upper=R> region[N];
  vector[N] maternal_education;
  vector[N] wealth_index;
  vector[N] iptp_use;
  vector[N] malaria_endemicity;
  vector[N] antenatal_visits;
  vector[N] timing_first_antenatal;
  vector[N] mosquito_net_use;
  vector[N] anemia_level;
  vector[N] malaria_history;  // New predictor
}

parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real mu_region;
  real<lower=0> sigma_region;
  vector[R] region_effects;
  vector[K] theta;
}

model {
  // Updated Priors
  alpha ~ normal(1, 0.5);
  beta ~ gamma(2, 1);  // More flexible prior

  mu_region ~ normal(0, 1);
  sigma_region ~ normal(0, 1);

  // Updated Informative Priors for Covariates
  theta[1] ~ normal(0.5, 0.5);  // Maternal education (positive effect expected)
  theta[2] ~ normal(-0.5, 0.5);  // Wealth index (negative effect expected)
  theta[3] ~ normal(-0.2, 0.3);  // IPTp use (negative effect expected)
  theta[4] ~ normal(-0.3, 0.3);  // Malaria endemicity (negative effect expected)
  theta[5] ~ normal(-0.1, 0.3);  // Number of antenatal visits (negative effect expected)
  theta[6] ~ normal(-0.2, 0.3);  // Timing of 1st antenatal check (negative effect expected)
  theta[7] ~ normal(-0.1, 0.3);  // Mosquito net use (negative effect expected)
  theta[8] ~ normal(0.2, 0.3);   // Anemia level (positive effect expected)
  theta[9] ~ normal(-0.3, 0.3);  // Malaria history (expected to increase mortality risk)

  region_effects ~ normal(mu_region, sigma_region);

  for (n in 1:N) {
    real lambda = exp(region_effects[region[n]] + 
                   theta[1] * maternal_education[n] + 
                   theta[2] * wealth_index[n] + 
                   theta[3] * iptp_use[n] + 
                   theta[4] * malaria_endemicity[n] +
                   theta[5] * antenatal_visits[n] +
                   theta[6] * timing_first_antenatal[n] +
                   theta[7] * mosquito_net_use[n] +
                   theta[8] * anemia_level[n] +
                   theta[9] * malaria_history[n]);  // Added malaria history
    real log_hazard = log(alpha) + log(beta) + (beta - 1) * log(survival_time[n]) -
                      log(1 + (alpha * survival_time[n])^beta);
    real log_survival = -log(1 + (alpha * survival_time[n])^beta);
    
    target += event_indicator[n] * log_hazard + (1 - event_indicator[n]) * log_survival;
  }
}

generated quantities {
  vector[N] y_rep;
  for (n in 1:N) {
    y_rep[n] = survival_time[n] + normal_rng(0, 0.1);
  }
}
"

# Save & Compile the Stan Model
writeLines(stan_model_code_new, "gll_model_hierarchical_updated.stan")
stan_model_updated <- stan_model("gll_model_hierarchical_updated.stan")

# Run MCMC Sampling with Optimized Settings
fit_new <- sampling(
  stan_model_updated,
  data = stan_data_list_new,
  iter = 4000, 
  warmup = 2000,  
  chains = 4,  
  seed = 123,  
  control = list(adapt_delta = 0.99, max_treedepth = 15)  
)

# Check Convergence Diagnostics
print(fit_new)  
check_hmc_diagnostics(fit_new)  
traceplot(fit_new, pars = c("alpha", "beta", "mu_region", "sigma_region", "theta"))  

# Posterior Predictive Checks (PPCs)
# Extract posterior predictive samples using as.array()
y_rep_new <- as.array(fit_new, pars = "y_rep")

# Alternatively, use rstan::extract() explicitly
y_rep_new <- rstan::extract(fit_new, pars = "y_rep")$y_rep


# Posterior Distributions for Parameters
mcmc_areas(as.array(fit_new), pars = c("alpha", "beta", "mu_region", "sigma_region", "theta[1]", "theta[2]", "theta[3]", "theta[4]", "theta[5]", "theta[6]", "theta[7]", "theta[8]","theta[9]"))

ggplot() +
  geom_density(aes(x = stan_data_list_new$survival_time, color = "Observed"), size = 1) + 
  geom_density(aes(x = apply(y_rep_new, 2, mean), color = "Posterior Predictive"), size = 1, linetype = "dashed") +
  scale_color_manual(name = "Legend", values = c("Observed" = "black", "Posterior Predictive" = "red")) + 
  labs(title = "Posterior Predictive Check: Observed vs. Posterior Predictive Distribution", 
       x = "Survival Time", 
       y = "Density") +
  theme_minimal()

# Posterior Distributions for Parameters
mcmc_areas(as.array(fit_new), pars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]", "theta[5]", "theta[6]", "theta[7]", "theta[8]","theta[9]"))




#################################################################

# Observed statistics
obs_median <- median(child_survival_data_final$Current.age.of.child.in.months..months.since.birth.for.dead.children.)
obs_event_rate <- mean(1 - child_survival_data_final$event_indicator)  # Death rate



# For each posterior draw, compute test statistics
sim_median <- apply(y_rep_new, 1, median)
sim_event_rate <- apply(y_rep_new, 1, mean)

# Two-tailed p-values
p_median <- mean(sim_median > obs_median) * 2  # Multiply by 2 for two-tailed
p_event <- mean(sim_event_rate > obs_event_rate) * 2



# Assuming 'y_rep' is a matrix from rstanarm::posterior_predict()
# Dimensions: [n_observations x n_posterior_draws]

# Correct: Compute statistics per posterior draw (column-wise)
sim_median <- apply(y_rep_new, 2, median)  # Column-wise median
sim_event_rate <- apply(y_rep_new, 2, mean) # Column-wise mean

# Two-tailed p-values (adjust if obs is outside simulated distribution)
p_median <- min(
  mean(sim_median <= obs_median) * 2,
  mean(sim_median >= obs_median) * 2
) %>% pmin(1.0)  # Cap at 1.0

p_event <- min(
  mean(sim_event_rate <= obs_event_rate) * 2,
  mean(sim_event_rate >= obs_event_rate) * 2
) %>% pmin(1.0)

# Print results
cat(sprintf(
  "Posterior predictive p-values:\n- Median survival: p = %.3f\n- Event rate: p = %.3f",
  p_median, p_event
))


# Check observed vs. predicted median survival times
observed_median <- quantile(child_survival_data_final$survival_time, 0.5) 
predicted_median_range <- quantile(apply(y_rep_new, 2, median)) 

# Print results
print(observed_median)
print(predicted_median_range)

# Ensure yrep is correctly passed
bayesplot::ppc_dens_overlay(
  y = child_survival_data_final$survival_time,  # Observed data
  yrep = y_rep_new[1:100, ]  # Predicted survival times (first 100 simulations)
)
################################################################
# Compute observed statistics
mean_obs <- mean(child_survival_data_final$survival_time)
var_obs <- var(child_survival_data_final$survival_time)
quant_obs <- quantile(child_survival_data_final$survival_time, probs = c(0.25, 0.5, 0.75))

# Compute posterior predictive statistics
mean_pred <- apply(y_rep_new, 2, mean)
var_pred <- apply(y_rep_new, 2, var)
quant_pred <- apply(y_rep_new, 2, function(x) quantile(x, probs = c(0.25, 0.5, 0.75)))

# Compute posterior predictive p-values
ppp_mean <- mean(mean_pred >= mean_obs)
ppp_var <- mean(var_pred >= var_obs)
ppp_quant <- apply(quant_pred, 1, function(q) mean(q >= quant_obs))

# Print results
list(
  Observed = list(mean = mean_obs, variance = var_obs, quantiles = quant_obs),
  Posterior_Predictive_CI = list(
    mean = quantile(mean_pred, probs = c(0.025, 0.5, 0.975)),
    variance = quantile(var_pred, probs = c(0.025, 0.5, 0.975)),
    quantiles = apply(quant_pred, 1, function(q) quantile(q, probs = c(0.025, 0.5, 0.975)))
  ),
  Posterior_Predictive_p_values = list(mean = ppp_mean, variance = ppp_var, quantiles = ppp_quant)
)

#################################################################



















