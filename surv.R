# from
# https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops4/Buros/weibull-survival-model.html


library(biostan)
library(dplyr)
library(ggplot2)
library(survival)
library(survMisc)
library(rstan)
library(purrr)
library(GGally)
library(ggfortify)


stan_file <- system.file('stan', 'weibull_survival_null_model.stan',
                         package='biostan')

stan_file <- '~/projectes/survival_example/weibull_survival_null_model.stan'
stan_file2 <- '~/projectes/survival_example/weibull_survival_null_model2.stan'

print_stan_file(stan_file)
print_stan_file(stan_file, section='model')


# Data simulation (wrong!)

sim_data <- function(alpha, mu, Nobs, Ncen){
  observed_data <- data.frame(os_status= rep_len('DECEASED', Nobs),
                              os_months= rweibull(Nobs, alpha, exp(-(mu/alpha))),
                              stringsAsFactors = FALSE)

  censored_data <- data.frame(os_status= rep_len('LIVING', Ncen),
                              os_months= runif(Ncen) * rweibull(Ncen, alpha,
                                                  exp(-(mu/alpha))),
                              stringsAsFactors = FALSE)

  return(observed_data %>% bind_rows(censored_data))
}


# Data simulation (correct)

sim_data <- function(alpha, mu, Nobs, Ncen) {
  data <- data.frame(
    surv_months = rweibull(n = Nobs + Ncen, alpha, exp(-(mu)/alpha)),
    censor_months = rexp(n = Nobs + Ncen,
                         rate = 1/100),
    stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      os_status = ifelse( surv_months < censor_months, 'DECEASED', 'LIVING'),
      os_months = ifelse(surv_months < censor_months,
                         surv_months,
                         censor_months))
    
    return(data)
}

test_alpha <- 0.8
test_mu <- -3 
test_Nobs <- 200
test_Ncen <- 300

# sim1 <- sim_data(test_alpha, test_mu, test_Nobs, test_Ncen)
sim1 <- sim_data(test_alpha, test_mu, test_Nobs, test_Ncen)
summary(sim1$os_months)
summary(as.factor(sim1$os_status))
head(sim1)

sim1 %>%
  dplyr::mutate(os_deceased= os_status == 'DECEASED') -> sim1_censored

head(sim1_censored)

p <- autoplot(survival::survfit(Surv(os_months, os_deceased) ~ os_status,
                 data= sim1_censored), conf.int= FALSE)

p2 <- autoplot(survival::survfit(Surv(os_months, os_deceased) ~ 1,
                 data= sim1_censored), conf.int= FALSE)

p
p2
#################

# Modelling simulated data

generate_stan_data <- function(data){
  data %>%
    filter(os_status == 'DECEASED') -> observed_data

  data %>%
    filter(os_status != 'DECEASED') -> censored_data


  stan_data <- list(
                    Nobs= nrow(observed_data),
                    Ncen= nrow(censored_data),
                    yobs= observed_data$os_months,
                    ycen= censored_data$os_months)
}

stan_data <- generate_stan_data(sim1)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

recover_simulated <- 
  stan(stan_file2,
       data= stan_data,
       chains= 4,
       iter = 1000,
       seed = 1111)

print(recover_simulated)

print_stan_file(stan_file, section = 'model')


gen_inits <- function(){
  list(
       alpha_raw= 0.01 * rnorm(1),
       mu= rnorm(1)
       )
}

recover_simulated2 <- 
  stan(stan_file,
       data= stan_data,
       chains= 4,
       iter= 1000,
       init= gen_inits)

print(recover_simulated2)

traceplot(recover_simulated2, 'lp__')
traceplot(recover_simulated2, c('alpha', 'mu'), ncol=1)


# Just non censored
recover_simulated_obs <- 
  stan(stan_file,
       data= generate_stan_data(sim1 %>% filter(os_status == 'DECEASED')),
       chains= 4,
       iter= 1000,
       init= gen_inits)

print(recover_simulated_obs)


# Just censored
recover_simulated_cens <- 
  stan(stan_file,
       data= generate_stan_data(sim1 %>% filter(os_status != 'DECEASED')),
       chains= 4,
       iter= 1000,
       init= gen_inits)

print(recover_simulated_cens)


pp_alpha <- extract(recover_simulated, 'alpha')$alpha
pp_mu <- extract(recover_simulated, 'mu')$mu

pp_newdata <- 
  map2(.x= pp_alpha,
       .y= pp_mu,
       .f= ~ sim_data(alpha= .x,
                      mu= .y,
                      Nobs= test_Nobs,
                      Ncen= test_Ncen))

ggplot(pp_newdata %>%
         dplyr::bind_rows() %>%
         dplyr::mutate(type = 'posterior predicted values') %>%
         bind_rows(sim1 %>% dplyr::mutate(type = 'actual data')),
  aes(x = os_months,
      group = os_status,
      colour = os_status,
      fill = os_status)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~type, ncol = 1)


pp_newdata %>%
  purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
  purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
  purrr::map(fortify) -> pp_survdata

pp_survdata %>%
  purrr::map(~ dplyr::mutate(., time_group = floor(time))) %>%
  dplyr::bind_rows() %>%
  dplyr::group_by(time_group) %>%
  dplyr::summarize(surv_mean = mean(surv),
                   surv_p50 = median(surv),
                   surv_lower = quantile(surv, probs = 0.025),
                   surv_upper = quantile(surv, probs = 0.975)) %>%
  dplyr::ungroup() -> pp_survdata_agg


## km-curve for test data 

sim1 %>% 
  dplyr::mutate(os_deceased = os_status == 'DECEASED') -> sim_dec

test_data_kmcurve <- 
  fortify(
      survfit(Surv(os_months, os_deceased) ~ 1,
             data= sim_dec))

pp_newdata %>%
  purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
  purrr::map(~ survfit(Surv(os_months, os_deceased) ~ 1, data = .)) ->
    new_data

ggplot(pp_survdata_agg %>%
       dplyr::mutate(type = 'posterior predicted values') %>%
       dplyr::rename(surv = surv_p50, lower = surv_lower, upper = surv_upper,
                     time = time_group) %>%
       bind_rows(test_data_kmcurve %>%
       dplyr::mutate(type = 'actual data')),
       aes(x = time,
           group = type,
           linetype = type)) + 
  geom_line(aes(y = surv,
                colour = type)) +
  geom_ribbon(aes(ymin = lower,
                  ymax = upper),
              alpha = 0.2) +
  xlim(c(0, 200)) +
  theme_bw()
