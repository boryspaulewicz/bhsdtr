## -*- coding: utf-8 -*-

library(rstan)
library(ggplot2)
library(Hmisc) # rMultinom for simulations
devtools::document('~/cs/code/r/bhsdtr')
devtools::install('~/cs/code/r/bhsdtr')
library(bhsdtr)
data(gabor)
## Standard Stan optimizations
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(plyr)

## Combined response has to be calculated for this dataset
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)

## Aggregation without information loss
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))

## Model specification
fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
model = make_stan_model(random)

## Stan data structure
sdata = make_stan_data(adata, fixed, random)

## Main model fit
fit = stan(model_code = model,
           data = sdata,
           pars = c('delta_fixed', 'gamma_fixed',
                    'delta_sd_1', 'gamma_sd_1',
                    'delta_random_1', 'gamma_random_1',
                    'Corr_delta_1', 'Corr_gamma_1',
                    ## we want counts_new for plotting
                    'counts_new'),
           iter = 8000,
           chains = 4)
## save(fit, file = '~/temp/fit')
load('~/temp/fit')

print(fit, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed',
               'delta_sd_1', 'gamma_sd_1',
               'Corr_delta_1', 'Corr_gamma_1'))

library(xtable)
smr = as.data.frame(round(summary(fit)$summary[,c(1, 2, 3, 4, 8, 9, 10)], 2))
smr$n_eff = as.integer(round(smr$n_eff))
smr = smr[-grep('random', rownames(smr)),]
smr = smr[-grep('counts_new', rownames(smr)),]
print(xtable(smr[-nrow(smr),]), file = 'fit_table.tex')

## Wykresy dopasowania
(p1 = plot_sdt_fit(fit, adata, c('order', 'duration')))
ggsave(p1, file = 'roc_fit.pdf')
(p2 = plot_sdt_fit(fit, adata, c('order', 'duration'), type = ''))
ggsave(p2, file = 'response_fit.pdf')

## Testujemy izomorficzność mapowania gamma-c
criteria = c(-2, -1, 0, 1, 2)
cumprobs = c(pnorm(criteria, sd = 2), 1)
areas = c(cumprobs[1], cumprobs[-1] - cumprobs[-length(cumprobs)])
gamma = log(areas / areas[length(areas)])
g_to_c = function(x, s = 2){
    x = exp(x)
    s * qnorm(cumsum(x/sum(x))[-length(x)])
}
rbind(g_to_c(gamma), criteria)

## Oglądamy kryteria w pierwszym warunku
crit = gamma_to_c(as.data.frame(fit))
round(apply(crit, 2, mean), 2)
HPDinterval(as.mcmc(crit))

######################################################################
## Dopasowanie do symulacji

data = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))

## To będą nasze realistyczne oszacowania parametrów
s = apply(as.data.frame(fit), 2, mean)
## Używamy tych samych danych do symulacji
sdata = make_stan_data(adata, fixed, random)
gamma_fixed = matrix(nrow = ncol(sdata$X_gamma), ncol = sdata$K - 1)
for(r in 1:nrow(gamma_fixed))
  gamma_fixed[r,] = s[paste('gamma_fixed[', 1:(sdata$K - 1), ',', r, ']', sep = '')]
gamma_random = list()
for(r in 1:sdata$G_1)
  ##! Nie wiem, czy nie trzeba dać byrow = T
  gamma_random[[r]] = matrix(s[paste('gamma_random_1[', r, ',', 1:(sdata$K - 1), ',1]', sep = '')],
                             nrow = sdata$K - 1, ncol = sdata$Z_gamma_ncol)
delta_fixed = s[grep('delta_fixed', names(s))]
delta_random = matrix(nrow = sdata$G_1, ncol = sdata$Z_delta_ncol_1)
for(r in 1:nrow(delta_random))
  for(i in 1:ncol(delta_random))
    delta_random[r, i] = s[sprintf('delta_random_1[%d,%d]', r, i)]
multinomial_p = matrix(ncol = sdata$K, nrow = sdata$N)
thr = gamma = matrix(ncol = sdata$K - 1, nrow = sdata$N)
dprim = rep(NA, sdata$N)
for(n in 1:sdata$N){
  dprim[n] = exp(sdata$X_delta[n,] %*% delta_fixed + sdata$Z_delta_1[n,] %*% delta_random[sdata$group_1[n],])
  gamma[n,] = t(gamma_fixed) %*% sdata$X_gamma[n,] + gamma_random[[sdata$group_1[n]]] %*% sdata$Z_gamma_1[n,]
}
criteria = t(apply(gamma, 1, function(x){
  exp_x_0 = exp(c(x, 0))
  sdata$criteria_scale * qnorm(cumsum(exp_x_0 / sum(exp_x_0))[-sdata$K])
}))
distr_shift = -sdata$stim_sign * dprim / 2
multinomial_p[, 1] = pnorm(criteria[, 1] + distr_shift)
for(k in 2:(sdata$K - 1)){
  multinomial_p[, k] = pnorm(criteria[, k] + distr_shift) - pnorm(criteria[, k - 1] + distr_shift)
}
multinomial_p[, sdata$K] = pnorm(-(criteria[, sdata$K - 1] + distr_shift))
## Możemy już symulować
data_sim = data
for(r in 1:nrow(data$counts))
  data_sim$counts[r,] = table(c(rMultinom(t(multinomial_p[r,]), sum(data$counts[r,])), 1:sdata$K)) - 1

fit.sim = stan(model_code = make_stan_model(random),
               pars = c('delta_fixed', 'gamma_fixed',
                          'delta_sd_1', 'gamma_sd_1',
                        'delta_random_1', 'gamma_random_1',
                        'counts_new'),
               data = make_stan_data(data_sim, list(delta = ~ -1 + duration:order, gamma = ~ order),
                                     list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))),
               chains = 4,
               iter = 8000)
## save(fit.sim, file = '~/temp/fit.sim')
## save(data_sim, file = '~/temp/data_sim')

load('~/temp/fit.sim')

print(fit.sim, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed', 'delta_sd_1', 'gamma_sd_1'))

######################################################################
## Wykresy dopasowania

load('~/temp/fit.sim')
load('~/temp/data_sim')

(p1 = plot_sdt_fit(fit.sim, data_sim, c('duration', 'order')))
ggsave('roc_sim_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.sim, data_sim, c('duration', 'order'), type = ''))
ggsave('response_sim_fit.pdf', p2)

######################################################################
## Porównanie oszacowań z wartościami zakładanymi w symulacji

s2 = as.data.frame(fit.sim)

res = as.data.frame(cbind(t(apply(s2[,grep('delta_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])), delta_fixed))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## 1

res = as.data.frame(cbind(t(apply(s2[,grep('gamma_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])),
                          as.vector(t(gamma_fixed))))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## 1

res = as.data.frame(cbind(t(apply(s2[,grep('delta_random', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])),
                          as.vector(delta_random)))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## .99

gamma_random_mat = matrix(ncol = length(gamma_random[[1]]), nrow = length(gamma_random))
for(r in 1:nrow(gamma_random_mat))
  for(col in 1:ncol(gamma_random_mat))
    gamma_random_mat[r, col] = gamma_random[[r]][col]
res = as.data.frame(cbind(t(apply(s2[,grep('gamma_random', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])),
                          as.vector(gamma_random_mat)))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## 0.98
c(sum(!res$in_ci), nrow(res))
## 6 / 329 = 0.18

######################################################################
## Fitting the model to simulated data aggregated over subjects

load('~/temp/data_sim')
df = data_sim$data
df$stimulus = data_sim$stimulus
df$i = 1:nrow(df)
df = ddply(df, c('order', 'duration', 'stimulus'), function(x)apply(data_sim$counts[x$i,], 2, sum))
data_sim_2 = list(data = df[,c('order', 'duration')], stimulus = df$stimulus, counts = df[, paste(1:8)])

fit.aggr = stan(model_code = make_stan_model(),
                pars = c('delta_fixed', 'gamma_fixed', 'counts_new'),
                data = make_stan_data(data_sim_2, list(delta = ~ -1 + duration:order, gamma = ~ 1 + order)),
                chains = 4,
                iter = 8000)

## save(fit.aggr, file = '~/temp/fit_aggr')

load('~/temp/fit_aggr')

print(fit.aggr, pars = c('delta_fixed', 'gamma_fixed'),
      probs = c(.025, .975))
## Próbkowanie bdb

## Model fit plots
(p1 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration')))
ggsave('roc_sim_aggr_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration'), type = ''))
ggsave('response_sim_aggr_fit.pdf', p2)

## Sprawdzamy, czy wartości prawdziwe są zawarte w CI

s2 = as.data.frame(fit.aggr)

res = as.data.frame(cbind(t(apply(s2[,grep('delta_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])), delta_fixed))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## .25

## (p = ggplot(res, aes(true, fit)) + geom_point() + geom_abline(slope = 1, intercept = 0) +
##   geom_errorbar(aes(ymin = lo, ymax = hi)) + labs(x = 'True value', y = 'Estimate'))
## ggsave('aggregated_fit_vs_true_delta_fixed.pdf', p)

res = as.data.frame(cbind(t(apply(s2[,grep('gamma_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])),
                          as.vector(t(gamma_fixed))))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## .36
cbind(sum(res$in_ci), nrow(res))
## 5 14
res$Condition = rep(1:(ncol(model.matrix(fixed$gamma, data_sim_2$data))), each = ncol(data_sim_2$counts_1) - 1)
## p = ggplot(res, aes(true, fit)) + geom_point() + geom_abline(slope = 1, intercept = 0) + geom_errorbar(aes(ymin = lo, ymax = hi)) +
##   facet_wrap(~Condition) + labs(x = 'True value', y = 'Estimate')
## ggsave('aggregated_fit_vs_true_gamma_fixed.pdf', p)

## Porównujemy posterior sd
aggr.posterior.sd = apply(as.data.frame(fit.aggr), 2, sd)
true.posterior.sd = apply(as.data.frame(fit.sim), 2, sd)
sd.ratios = true.posterior.sd[grep('fixed', names(true.posterior.sd))] /
  aggr.posterior.sd[grep('fixed', names(aggr.posterior.sd))]
mean(sd.ratios)
## 3.126772

######################################################################
## Rozkłady efektów losowych

(p = ggplot(data.frame(delta = as.vector(t(delta_random)), i = 1:ncol(delta_random)), aes(sample = delta)) +
     stat_qq() +
     facet_wrap(~i) +
     ylab('delta sample quantile') +
     xlab('Theoretical normal quantile'))
ggsave('delta_qq_plot.pdf', p)

shapiro.test(delta_random[,1])
## W = 0.98841, p-value = 0.9181
shapiro.test(delta_random[,2])
## W = 0.96939, p-value = 0.2512
gamma_normality = data.frame(k = 1:7, statistic = NA, p.value = NA)
for(i in 1:ncol(gamma_random_mat)){
  res = shapiro.test(gamma_random_mat[,i])
  gamma_normality$statistic[i] = res$statistic
  gamma_normality$p.value[i] = res$p.value
}
round(gamma_normality, 3)
##   k statistic p.value
## 1 1     0.974   0.372
## 2 2     0.992   0.983
## 3 3     0.958   0.092
## 4 4     0.979   0.552
## 5 5     0.978   0.527
## 6 6     0.945   0.028
## 7 7     0.985   0.792

gamma_rnd = gamma_random_mat[,1]
for(i in 2:ncol(gamma_random_mat))
  gamma_rnd = c(gamma_rnd, gamma_random_mat[,i])
k = rep(1:ncol(gamma_random_mat), each = nrow(gamma_random_mat))
(p = ggplot(data.frame(gamma = gamma_rnd, k = k), aes(sample = gamma)) + stat_qq() + facet_wrap(~k) +
     ylab('gamma sample quantile') +
     xlab('Theoretical normal quantile'))
ggsave('gamma_qq_plot.pdf', p)
