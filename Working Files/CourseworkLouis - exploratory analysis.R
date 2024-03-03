library(dplyr)
library(ggplot2)

load("Datasets/hsesub.Rdata")

# Assume drinkYN_19 is binomally distributed with unknown probability p and calculate confidence intervals this way
# Can do similar for cigarette/e-cig consumption using cigsta3_19 and NDPNow_19 variables

calculate_binomial_ci <- function(dat_vec, alpha = 0.05) {
  dat_vec <- as.numeric(dat_vec)
  
  n <- length(dat_vec)
  
  p_hat <- mean(dat_vec)
  
  var_est <- p_hat * (1 - p_hat) / n
  
  ci <- p_hat + c(-1, 1) * qnorm(1 - alpha / 2) * sqrt(var_est)
  
  list(n = n, p_hat = p_hat, alpha = alpha, ci = ci)
}

# drinking
overall_drinking <- calculate_binomial_ci(na.omit(subdat$drinkYN_19) == 2)

100 * overall_drinking$ci

# smoking
overall_smoking <- calculate_binomial_ci(na.omit(subdat$cigsta3_19) == 1)

100 * overall_smoking$ci

# e-cigs
overall_ecigs <- calculate_binomial_ci(na.omit(subdat$NDPNow_19) %in% c(1, 3))

100 * overall_ecigs$ci

# Can also split this over age groups (or UK location/ethnic groups)

age_groups <- attr(subdat$ag16g10, 'labels')

age_groups_drinking <- lapply(age_groups, \(x) {
  dat_vec <- subdat %>% dplyr::filter(ag16g10 == x) %>% pull(drinkYN_19)
  
  dat_vec <- na.omit(dat_vec) == 2
  
  summary <- calculate_binomial_ci(dat_vec)
  
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
  dat_vec <- subdat %>% dplyr::filter(ag16g10 == x) %>% pull(cigsta3_19)
  
  dat_vec <- na.omit(dat_vec) == 1
  
  summary <- calculate_binomial_ci(dat_vec)
  
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

