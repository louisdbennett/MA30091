library(dplyr)
library(ggplot2)

load("Datasets/hsesub.Rdata")

# Assume drinkYN_19 is Bernoulli distributed with unknown probability p and calculate confidence intervals this way
# Can do similar for cigarette/e-cig consumption using cigsta3_19 and NDPNow_19 variables

calculate_binomial_ci <- function(dat_vec, weights = rep(1, length.out = length(dat_vec)), alpha = 0.05) {
  dat_vec <- as.numeric(dat_vec)
  
  n <- length(dat_vec)
  
  # Need to check if you can weight like this!
  p_hat <- sum(dat_vec * weights) / sum(weights)
  
  var_est <- p_hat * (1 - p_hat) / sum(weights)
  
  ci <- p_hat + c(-1, 1) * qnorm(1 - alpha / 2) * sqrt(var_est)
  
  list(n = n, p_hat = p_hat, alpha = alpha, ci = ci)
}

# drinking
drinking_na <- is.na(subdat$drinkYN_19)
overall_drinking <- calculate_binomial_ci(subdat$drinkYN_19[!drinking_na] == 2, subdat$wt_int[!drinking_na])

100 * overall_drinking$ci

# smoking
smoking_na <- is.na(subdat$cigsta3_19)
overall_smoking <- calculate_binomial_ci(subdat$cigsta3_19[!smoking_na] == 1, subdat$wt_int[!smoking_na])

100 * overall_smoking$ci

# e-cigs
ecig_na <- is.na(subdat$NDPNow_19)
overall_ecigs <- calculate_binomial_ci(subdat$NDPNow_19[!ecig_na] %in% c(1, 3), subdat$wt_int[!ecig_na])

100 * overall_ecigs$ci

# Can also split this over age groups (or UK location/ethnic groups)

age_groups <- attr(subdat$ag16g10, 'labels')

age_groups_drinking <- lapply(age_groups, \(x) {
  dat <- subdat %>% dplyr::filter(ag16g10 == x, !is.na(drinkYN_19)) %>% select(wt_int, drinkYN_19)
  
  dat_vec <- dat$drinkYN_19 == 2
  
  summary <- calculate_binomial_ci(dat_vec, dat$wt_int)
  
  tibble::tibble(p_hat = summary$p_hat, lower = summary$ci[1], upper = summary$ci[2])
}) %>% 
  bind_rows(.id = 'age_group') %>% 
  filter(!age_group %in% names(age_groups[1:3]))

ggplot(data = age_groups_drinking, aes(x = age_group, y = p_hat, ymax = upper, ymin = lower)) +
  geom_bar(stat = 'identity', alpha = 0.8, fill = '#b41313') +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  xlab('Age Group') +
  ylab('% Drank Alcohol Last 12 Months') +
  theme_minimal()

# Current smokers

age_groups_smoking <- lapply(age_groups, \(x) {
  dat <- subdat %>% dplyr::filter(ag16g10 == x, !is.na(cigsta3_19)) %>% select(wt_int, cigsta3_19)
  
  dat_vec <- dat$cigsta3_19 == 1
  
  summary <- calculate_binomial_ci(dat_vec, dat$wt_int)
  
  tibble::tibble(p_hat = summary$p_hat, lower = summary$ci[1], upper = summary$ci[2])
}) %>% 
  bind_rows(.id = 'age_group') %>% 
  filter(!age_group %in% names(age_groups[1:3]))

ggplot(data = age_groups_smoking, aes(x = age_group, y = p_hat, ymax = upper, ymin = lower)) +
  geom_bar(stat = 'identity', alpha = 0.8, fill = '#b41313') +
  geom_errorbar(width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  xlab('Age Group') +
  ylab('% Current Smokers') +
  theme_minimal()


# Assume a poisson distribution for the number of cigarettes smoked a week
dat_vec <- na.omit(subdat$cigdyal_19[subdat$cigdyal_19 != 0]) * 7
sum(dat_vec) / length(dat_vec)

## Question two model

# Cleaning data poisson assumption

relevant_factors <- c('Sex', 'topqual2', 'marstatD', 'urban14b', 'origin2', 'GOR1', 'Age35g', 'qimd19')

subdat_clean <- subdat %>% 
  filter(!ag16g10 %in% c(1:3), !is.na(cigdyal_19)) %>% 
  mutate(
    across(all_of(relevant_factors), as.factor), cigweek_19 = log(7 * cigdyal_19 + 0.000001)
  ) %>% 
  select(
    cigweek_19, all_of(relevant_factors)
  )

full_model <- glm(formula = cigweek_19 ~ Age35g:qimd19 + Age35g:Sex + ., family = gaussian, data = subdat_clean)
null_model <- glm(formula = cigweek_19 ~ 1, family = gaussian, data = subdat_clean)

fitted_model <- step(full_model, scope = list(lower = null_model, upper = full_model), direction = 'backward', trace = 0)

summary(fitted_model)

formula <- cigweek_19 ~ Sex + topqual2 + marstatD + qimd19 + origin2 + Age35g

mod <- glm(formula = formula, family = gaussian, data = subdat_clean)

summary(mod)

plot(mod, 2)

# Try with binomial model

subdat_clean <- subdat %>% 
  filter(!ag16g10 %in% c(1:3), !is.na(cigsta3_19)) %>% 
  mutate(
    across(all_of(relevant_factors), as.factor), smokes = as.numeric(cigsta3_19 == 1)
  ) %>% 
  select(
    smokes, all_of(relevant_factors)
  )

full_model <- glm(formula = smokes ~ Age35g:qimd19 + Age35g:Sex + ., family = binomial, data = subdat_clean)
null_model <- glm(formula = smokes ~ 1, family = binomial, data = subdat_clean)

fitted_model <- step(full_model, scope = list(lower = null_model, upper = full_model), direction = 'backward', trace = 0)

summary(fitted_model)

plot(fitted_model, 2)
