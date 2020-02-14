ggs = function(f)ggsave(f, width = 10, height = 10, dpi = 200)

######################################################################
## KOD 1

library(ggplot2)
library(bhsdtr)
library(Hmisc)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## r = 1 ~ absolutely clear image left, r = 4 ~ no experience left, r = 5 ~ no experience right, ... r = 8 ~  absolutely clear image right
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
d1 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
d1$r = combined_response()
## Only the left-tilted gabor patches
d1 = d1[d1$stim == 0,]

######################################################################
## KOD 2

## This order-preserving link function represents the main criterion
## directly as a free paramter and all the other criteria as log of
## distance between the criterion and the criterion adjacent to it
## which is closer to the main criterion

adata = aggregate_responses(d1, response = 'r', variables = c('id'))
fixed = list(eta = ~ 1, gamma = ~ 1)
random = list(list(group = ~ id, gamma = ~ 1))
model = make_stan_model(random, model = 'ordinal')
sdata = make_stan_data(adata, fixed, random, model = 'ordinal')

######################################################################
## KOD 3

o = optimizing(stan_model(model_code = model), sdata, init = 0)
if(o$return_code)warning('not ok')

######################################################################
## Wykres

sdt_ml_plot(model, adata, sdata)
ggs('plot1.jpg')

######################################################################
## Wykres

gamma_fixed = o$par[grep('gamma_fixed\\[', names(o$par))]
gamma_random = matrix(o$par[grep('gamma_random_1', names(o$par))], nrow = sdata$group_1_size)
for(i in 1:nrow(gamma_random))
    gamma_random[i,] = gamma_random[i,] + gamma_fixed
criteria = gamma_to_crit(gamma_random, beta_index = NULL, gamma_link = link)
df = expand.grid(criterion = 1:7, id = 1:nrow(criteria))
df$estimate = as.vector(t(criteria))
ggplot(df, aes(id, estimate, group = id)) + geom_line() + geom_point()

######################################################################
## Model nasycony a efektem bod¼ca

fixed = list(eta = ~ 1, gamma = ~ stim)
random = list(list(group = ~ id, gamma = ~ stim))
d2 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
## Note that because for an ordinal model the stimulus ('stim') qvariable is not
## special we have to include it in the additional variables vector 
adata = aggregate_responses(d2, response = 'r', variables = c('id', 'stim'))
sdata = make_stan_data(adata, fixed, random, model = 'ordinal', gamma_link = link)
model = make_stan_model(random, model = 'ordinal', gamma_link = link)

fit0 = stan(model_code = model, data = sdata,
            pars = c('gamma_fixed', 'gamma_random_1', 'gamma_sd_1'),
            init_r = .5,
            iter = 5000,
            chains = 4)
## saveRDS(fit0, '~/windows/temp/ordinal_fit_0.rds')
print(fit0, probs = c(.025, .957), pars = c('gamma_fixed', 'gamma_sd_1', 'gamma_random_1'))

## _sd_1 ma byæ macierz±, ¿eby by³o wiadomo o co chodzi

######################################################################
## KOD 4

fixed = list(delta = ~ 1, gamma = ~ 1)
random = list(list(group = ~ id, delta = ~ 1, gamma = ~ 1))
d2 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
adata = aggregate_responses(d2, 'stim', 'r', 'id')
sdata = make_stan_data(adata, fixed, random, model = 'sdt')
model = make_stan_model(random, model = 'sdt')

fit1 = stan(model_code = model, data = sdata,
           pars = c('delta_fixed', 'gamma_fixed', 'delta_random_1', 'gamma_random_1'),
           init_r = .5,
           iter = 5000,
           chains = 4)
## saveRDS(fit1, '~/windows/temp/ordinal_fit_1.rds')
print(fit1, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))

## Now we will simulate from the known model
adata_sim = simulate_from_fit(fit1, adata, fixed, random)
sdata_sim = make_stan_data(adata_sim, fixed, random, model = 'sdt')
fit1.1 = stan(model_code = model, data = sdata_sim,
              pars = c('delta_fixed', 'gamma_fixed', 'delta_random_1', 'gamma_random_1'),
              init_r = .5,
              iter = 5000,
              chains = 4)
## saveRDS(fit1.1, '~/windows/temp/ordinal_fit_1.1.rds')
print(fit1.1, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))

cor(apply(as.data.frame(fit1)[,1:10], 2, mean), apply(as.data.frame(fit1.1)[,1:10], 2, mean))
## .97

s1 = as.data.frame(fit1.1)
## First we have to add fixed and random effects to obtain participant
## specific delta estimates and than we have to translate between the
## delta and the d' parameters
dprim1 = exp(s1[['delta_fixed[1,1]']] + s1[, grep('delta_random', names(s1))])

######################################################################
## KOD 5

fixed = list(delta = ~ 1, gamma = ~ 1)
## There are no participant effects in the thresholds in this model
random = list(list(group = ~ id, delta = ~ 1))
model = make_stan_model(random, model = 'sdt')

fit2 = stan(model_code = model, data = sdata_sim,
            pars = c('delta_fixed', 'gamma_fixed', 'delta_random_1'),
            init_r = .5,
            iter = 8000,
            chains = 4)
## saveRDS(fit2, '~/windows/temp/ordinal_fit_2.rds')
print(fit2, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))

s2 = as.data.frame(fit2)
dprim2 = exp(s2[['delta_fixed[1,1]']] + s2[, grep('delta_random', names(s2))])

df = data.frame(dprim2 = apply(dprim2, 2, mean), dprim1 = apply(dprim1, 2, mean))
ggplot(df, aes(dprim1, dprim2)) + geom_point() + labs(title = sprintf('r = %.2f', cor(df$dprim1, df$dprim2))) +
    xlab('d\' estimates from the true model') + ylab('d\' estimates from the simplified model')
ggs('plot2.jpg')
## This does not look good for the simplified model

## What about the interval estimates?
fun = function(x)quantile(x, c(.025, .975))
dprim1_ci = as.data.frame(t(apply(dprim1, 2, fun)))
dprim2_ci = as.data.frame(t(apply(dprim2, 2, fun)))
dprim1_ci$model = 'true'
dprim1_ci$id = 1:ncol(dprim1) - .1
dprim2_ci$model = 'simplified'
dprim2_ci$id = 1:ncol(dprim1) + .1
df = rbind(dprim1_ci, dprim2_ci)
names(df)[1:2] = c('lo', 'hi')
ggplot(df, aes(id, group = model, color = model)) + geom_errorbar(aes(ymin = lo, ymax = hi)) + ylab('d\'')
ggs('plot3.jpg')

dprim1_ci$width = dprim1_ci[,2] - dprim1_ci[,1]
dprim2_ci$width = dprim2_ci[,2] - dprim2_ci[,1]
df = data.frame(ci_width_ratio = dprim2_ci$width / dprim1_ci$width)
df$i = 1:nrow(df)
ggplot(df, aes(i, ci_width_ratio)) + geom_point() + geom_abline(intercept = 1, slope = 0) +
    ylab('CI width ratio: simplified model / true model')
ggs('plot4.jpg')

plot(df$simplified_ci / df$true_ci, ylab = 'CI width ratio: simplified model / true model') + abline(1, 0)

######################################################################
## Saturated model with stimulus effects

link = 'log_distance'
fixed = list(delta = ~ 1, gamma = ~ 1)
random = list(list(group = ~ id, delta = ~ 1, gamma = ~ 1))
d2 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
adata = aggregate_responses(d2, 'stim', 'r', 'id')
sdata = make_stan_data(adata, fixed, random, gamma_link = link, model = 'sdt')
model = make_stan_model(random, gamma_link = link, model = 'sdt')

fit1 = stan(model_code = model, data = sdata,
           pars = c('delta_fixed', 'gamma_fixed', 'delta_random_1', 'gamma_random_1'),
           init_r = .5,
           iter = 5000,
           chains = 4)
## saveRDS(fit1, '~/windows/temp/ordinal_fit_1.rds')
print(fit1, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))
