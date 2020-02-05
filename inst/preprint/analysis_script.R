## -*- coding: utf-8 -*-

## TODO Include the analysis for the nested model (comparing the true
## hierarchical and non-hierarchical estimates)

######################################################################
## IMPORTANT GLOBAL VARIABLES
##
## Change this to FALSE once all the stanfit objects are stored to
## save time and memory
fresh_start = FALSE
## Where do you want to store all the (large) stanfit objects? (I use
## Linux so '~/path' makes sense to me)
temp_path = '~/temp'

## WARNING we are not setting the random seed here because 1) R and
## Stan have separate seeds and 2) the user of this script may use a
## different number of mcmc threads.

library(rstan)
library(ggplot2)
library(Hmisc) # rMultinom for simulations
library(bhsdtr)
library(plyr)
library(gridExtra)
## Standard Stan optimizations
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data(gabor)

######################################################################
## Simulating individual data for sdt models

sim_sdt = function(n = 1, dprim = 1.5, criteria = c(-2.1, -1.4, -.7, 0, .7, 1.4, 2.1), model = 'sdt', sd_ratio = 1){
    which_bin = function(x, thr)min(which(x <= c(-Inf, thr, Inf)) - 1)
    d = data.frame(stim = rep(1:2, each = n), e = rnorm(n * 2), r = NA)
    for(i in 1:nrow(d))
        d$r[i] = which_bin(d$e[i] + .5 * dprim * c(-1, 1)[d$stim[i]], criteria)
    attr(d, 'dprim') = dprim
    attr(d, 'criteria') = criteria
    attr(d, 'sd_ratio') = sd_ratio
    d
}

######################################################################
## Quick test of all the possible types of non-hierarchical models
## fitted to simulated data

d = sim_sdt(n = 10000)
table(d[, c('stim', 'r')])
adata = aggregate_responses(d, 'stim', 'r')

res = expand.grid(par = c('dprim', paste('c', 1:length(attr(d, 'criteria')), sep = '')),
                  model = c('sdt', 'uvsdt', 'metad'),
                  link = c('softmax', 'log_distance', 'log_ratio'), true = NA, est = NA)

for(link in levels(res$link)){
    for(model in levels(res$model)){
        print(sprintf('Fitting model %s with link %s', model, link))
        fixed = list(delta = ~ 1, gamma = ~ 1)
        if(model == 'uvsdt')fixed$theta = ~ 1
        fit = stan(model_code = make_stan_model(gamma_link = link, model = model),
                   data = make_stan_data(adata, fixed, gamma_link = link, model = model),
                   pars = c('delta_fixed', 'gamma_fixed',
                            ## we need counts_new for plotting
                            'counts_new'),
                   init_r = .5,
                   iter = 4000,
                   chains = 4)
        s = as.data.frame(fit)
        res[(res$link == link) & (res$model == model), 'est'] =
            c(mean(exp(s[, grep('delta_fixed\\[1', names(s))])), apply(gamma_to_crit(s, gamma_link = link), 2, mean))
        res[(res$link == link) & (res$model == model), 'true'] = c(dprim = attr(d, 'dprim'), attr(d, 'criteria'))
        rm(s)
        rm(fit)
        gc()
    }
}

round(cor(res$est, res$true, use = 'pairwise.complete.obs'), 2)
## 1

## Comparing point estimates to known true values
res$par_type = c('criteria', 'dprim')[(as.character(res$par) == 'dprim') + 1]
ggplot(res, aes(est, true, group = par_type, color = par_type)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(size = 5) +
    facet_grid(model ~ link)
## All is good

######################################################################
## Fitting the model to real study data

## We will be using the gabor dataset provided with the package
data(gabor)
?gabor

## Combined response has to be calculated for this dataset
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)

## Aggregation without information loss
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))

## Model specification
fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
model = make_stan_model(random)
## The model code is human-readable
cat(model)

## Required Stan data structure
sdata = make_stan_data(adata, fixed, random)

## Main model fit
if(fresh_start){
    fit = stan(model_code = model,
               data = sdata,
               pars = c('delta_fixed', 'gamma_fixed',
                        'delta_sd_1', 'gamma_sd_1',
                        'delta_random_1', 'gamma_random_1',
                        'Corr_delta_1', 'Corr_gamma_1',
                        ## we need counts_new for plotting
                        'counts_new'),
               iter = 8000,
               chains = 4)
    save(fit, file = paste(temp_path, 'fit', sep = '/'))
}else{
    load(paste(temp_path, 'fit', sep = '/'))
}

## if(fresh_start){
##     fixed = list(delta = ~ -1 + duration:order, gamma = ~ -1 + order)
##     random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
##     model = make_stan_model(random)
##     sdata = make_stan_data(adata, fixed, random)
##     fit_nested = stan(model_code = model,
##                       data = sdata,
##                       pars = c('delta_fixed', 'gamma_fixed',
##                                'delta_sd_1', 'gamma_sd_1',
##                                'delta_random_1', 'gamma_random_1',
##                                'Corr_delta_1', 'Corr_gamma_1',
##                                ## we need counts_new for plotting
##                                'counts_new'),
##                       iter = 8000,
##                       chains = 4)
##     save(fit_nested, file = paste(temp_path, 'fit_nested', sep = '/'))
## }else{
##     load(paste(temp_path, 'fit_nested', sep = '/'))
## }

## Model fit summary does not indicate any convergence issues
print(fit, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed',
               'delta_sd_1', 'gamma_sd_1',
               'Corr_delta_1', 'Corr_gamma_1'))

## Model fit summary table for the paper
library(xtable)
smr = as.data.frame(round(summary(fit)$summary[,c(1, 2, 3, 4, 8, 9, 10)], 2))
smr$n_eff = as.integer(round(smr$n_eff))
smr = smr[-grep('random', rownames(smr)),]
smr = smr[-grep('counts_new', rownames(smr)),]
print(xtable(smr[-nrow(smr),]), file = 'fit_table.tex')

## Some plots
(p1 = plot_sdt_fit(fit, adata, c('order', 'duration')))
ggsave(p1, file = 'roc_fit.pdf')
ggsave(p1, file = 'roc_fit.png')
(p2 = plot_sdt_fit(fit, adata, c('order', 'duration'), type = 'response'))
ggsave(p2, file = 'response_fit.pdf')
ggsave(p2, file = 'response_fit.png')

## The criteria are placed fairly symmetrically
crit = gamma_to_crit(as.data.frame(fit))
round(apply(crit, 2, mean), 2)

######################################################################
## Fitting the single criterion hierarchical SDT model

fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
## this time we ignore the ratings
gabor$r = combined_response(gabor$stim, accuracy = gabor$acc)
if(fresh_start){
    fit.s = stan(model_code = make_stan_model(random),
                 data = make_stan_data(aggregate_responses(gabor, 'stim', 'r', c('duration', 'order', 'id')),
                                       fixed, random),
                 pars = c('delta_fixed', 'gamma_fixed',
                          'delta_sd_1', 'gamma_sd_1',
                          'delta_random_1', 'gamma_random_1',
                          'Corr_delta_1', 'Corr_gamma_1',
                          'counts_new'),
                 iter = 8000,
                 chains = 4)
    save(fit.s, file = paste(temp_path, 'fit.s', sep = '/'))
}else{
    load(paste(temp_path, 'fit.s', sep = '/'))
}

print(fit.s, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed',
               'delta_sd_1', 'gamma_sd_1',
               'Corr_delta_1', 'Corr_gamma_1'))
## All is fine

## It does not make sense to use plot_sdt_fit here, because the
## single-criterion SDT model is untestable

## Compare delta estimates
smr1 = as.data.frame(summary(fit)$summary[,c(1, 2, 3, 4, 8, 9, 10)])
smr2 = as.data.frame(summary(fit.s)$summary[,c(1, 2, 3, 4, 8, 9, 10)])
smr = rbind(smr1[grep('delta', rownames(smr1)),],
            smr2[grep('delta', rownames(smr2)),])
smr$model = as.factor(rep(c('Multiple criteria', 'Single criterion'), each = nrow(smr) / 2))
names(smr)[4:5] = c('ci.lo', 'ci.hi')
smr$par = rep(rownames(smr)[1:(nrow(smr) / 2)], 2)
(p = ggplot(smr, aes(par, mean, color = model)) +
     geom_point() +
     geom_errorbar(aes(ymin = ci.lo, ymax = ci.hi)) +
     theme(axis.text.x = element_text(angle = 90)) +
     labs(color = 'Model', y = 'Posterior point and interval estimates'))
ggsave('single_vs_multiple_c.pdf', p)
## Similar but not exactly the same

## Let's compare the main criteria in the base condition
c1 = gamma_to_crit(as.data.frame(fit))
c2 = gamma_to_crit(as.data.frame(fit.s))
round(rbind(quantile(c1[,4], c(.025, .5, .957)),
            quantile(c2, c(.025, .5, .957))), 2)
## Similar but not exactly the same
##
##       2.5%  50% 95.7%
## [1,] -0.05 0.11  0.25
## [2,] -0.13 0.06  0.21

smr2$n_eff = as.integer(round(smr2$n_eff))
smr2 = smr2[-grep('random', rownames(smr2)),]
smr2 = smr2[-grep('counts_new', rownames(smr2)),]
print(xtable(round(smr2[-nrow(smr2),], 2)), file = 'fit_table_single_c.tex')

######################################################################
## Fitting the model to simulated data to test if the model recovers
## known realistic parameter values

data = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))

simulate_from_fit = function(fit, adata, fixed, random){
    ## Realistic known parameter values
    s = apply(as.data.frame(fit), 2, mean)
    sdata = make_stan_data(adata, fixed, random)
    gamma_fixed = matrix(nrow = ncol(sdata$X_gamma), ncol = sdata$K - 1)
    for(r in 1:nrow(gamma_fixed))
        gamma_fixed[r,] = s[paste('gamma_fixed[', 1:(sdata$K - 1), ',', r, ']', sep = '')]
    gamma_random = list()
    for(r in 1:sdata$G_1)
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
    data_sim = adata
    ## Now we can simulate
    for(r in 1:nrow(adata$counts))
        data_sim$counts[r,] = table(c(rMultinom(t(multinomial_p[r,]), sum(adata$counts[r,])), 1:sdata$K)) - 1
    data_sim$dprim = dprim
    data_sim$gamma = gamma
    data_sim$criteria = criteria
    data_sim
}

if(fresh_start){
    adata_sim = simulate_from_fit(fit, adata, fixed, random)
    save(adata_sim, file = paste(temp_path, 'data_sim', sep = '/'))
}else{
    load(paste(temp_path, 'data_sim', sep = '/'))
}

## if(fresh_start){
##     fixed = list(delta = ~ -1 + duration:order, gamma = ~ -1 + order)
##     random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
##     adata_sim_nested = simulate_from_fit(fit_nested, adata, fixed, random)
##     save(adata_sim, file = paste(temp_path, 'adata_sim_nested', sep = '/'))
## }else{
##     load(paste(temp_path, 'adata_sim_nested', sep = '/'))
## }

if(fresh_start){
    fit.sim = stan(model_code = make_stan_model(random),
                   pars = c('delta_fixed', 'gamma_fixed',
                            'delta_sd_1', 'gamma_sd_1',
                            'delta_random_1', 'gamma_random_1',
                            'counts_new'),
                   data = make_stan_data(adata_sim, list(delta = ~ -1 + duration:order, gamma = ~ order),
                                         list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))),
                   chains = 4,
                   iter = 8000)
    save(fit.sim, file = paste(temp_path, 'fit.sim', sep = '/'))
}else{
    load(paste(temp_path, 'fit.sim', sep = '/'))
}

## if(fresh_start){
##     fit.sim_nested = stan(model_code = make_stan_model(random),
##                           pars = c('delta_fixed', 'gamma_fixed',
##                                    'delta_sd_1', 'gamma_sd_1',
##                                    'delta_random_1', 'gamma_random_1',
##                                    'counts_new'),
##                           data = make_stan_data(adata_sim_nested, list(delta = ~ -1 + duration:order, gamma = ~ -1 + order),
##                                                 list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))),
##                           chains = 4,
##                           iter = 8000)
##     save(fit.sim_nested, file = paste(temp_path, 'fit.sim', sep = '/'))
## }else{
##     load(paste(temp_path, 'fit.sim_nested', sep = '/'))
## }

print(fit.sim, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed', 'delta_sd_1', 'gamma_sd_1'))
## All is fine

(p1 = plot_sdt_fit(fit.sim, adata_sim, c('duration', 'order')))
ggsave('roc_sim_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.sim, adata_sim, c('duration', 'order'), type = 'response'))
ggsave('response_sim_fit.pdf', p2)

######################################################################
## Comparing true values with the estimates based on the true model

s2 = as.data.frame(fit.sim)
s = apply(as.data.frame(fit), 2, mean)
sdata = make_stan_data(adata, fixed, random)
gamma_fixed = matrix(nrow = ncol(sdata$X_gamma), ncol = sdata$K - 1)
for(r in 1:nrow(gamma_fixed))
    gamma_fixed[r,] = s[paste('gamma_fixed[', 1:(sdata$K - 1), ',', r, ']', sep = '')]
gamma_random = list()
for(r in 1:sdata$G_1)
    gamma_random[[r]] = matrix(s[paste('gamma_random_1[', r, ',', 1:(sdata$K - 1), ',1]', sep = '')],
                               nrow = sdata$K - 1, ncol = sdata$Z_gamma_ncol)
delta_fixed = s[grep('delta_fixed', names(s))]
delta_random = matrix(nrow = sdata$G_1, ncol = sdata$Z_delta_ncol_1)
for(r in 1:nrow(delta_random))
    for(i in 1:ncol(delta_random))
        delta_random[r, i] = s[sprintf('delta_random_1[%d,%d]', r, i)]

res = as.data.frame(cbind(t(apply(s2[,grep('delta_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])), delta_fixed))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## 1 (all the true values are within the 95% credible intervals)

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
## .98

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
## 5 / 329 = 0.01

######################################################################
## Fitting the simplified model to simulated data aggregated over
## participants. This is how SDT models with additional criteria are
## often (mis)used.

load(paste(temp_path, 'data_sim', sep = '/'))
df = adata_sim$data
df$stimulus = adata_sim$stimulus
df$i = 1:nrow(df)
df = ddply(df, c('order', 'duration', 'stimulus'), function(x)apply(adata_sim$counts[x$i,], 2, sum))
data_sim_2 = list(data = df[,c('order', 'duration')], stimulus = df$stimulus, counts = df[, paste(1:8)])

if(fresh_start){
    fit.aggr = stan(model_code = make_stan_model(),
                    pars = c('delta_fixed', 'gamma_fixed', 'counts_new'),
                    data = make_stan_data(data_sim_2, list(delta = ~ -1 + duration:order, gamma = ~ order)),
                    chains = 4,
                    iter = 8000)
    save(fit.aggr, file = paste(temp_path, 'fit_aggr', sep = '/'))
}else{
    load(paste(temp_path, 'fit_aggr', sep = '/'))
}

## if(fresh_start){
##     fit.aggr_nested = stan(model_code = make_stan_model(),
##                     pars = c('delta_fixed', 'gamma_fixed', 'counts_new'),
##                     data = make_stan_data(data_sim_2, list(delta = ~ -1 + duration:order, gamma = ~ -1 + order)),
##                     chains = 4,
##                     iter = 8000)
##     save(fit.aggr_nested, file = paste(temp_path, 'fit_aggr_nested', sep = '/'))
## }else{
##     load(paste(temp_path, 'fit_aggr_nested', sep = '/'))
## }

## Excellent sampling
print(fit.aggr, pars = c('delta_fixed', 'gamma_fixed'),
      probs = c(.025, .975))

## The simplified model seems to fit the data very well, but...
(p1 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration')))
ggsave('roc_sim_aggr_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration'), type = 'response'))
ggsave('response_sim_aggr_fit.pdf', p2)

######################################################################
## A model with fixed participant effects

sdata.fixed = make_stan_data(adata_sim, fixed = list(delta = ~ -1 + as.factor(id):duration:order,
                                                    gamma = ~ -1 + as.factor(id) + order))
if(fresh_start){
    fit.fixed = stan(model_code = make_stan_model(),
                     data = sdata.fixed,
                     pars = c('delta_fixed', 'gamma_fixed'),
                     iter = 8000,
                     chains = 4)
    save(fit.fixed, file = paste(temp_path, 'fit.fixed', sep = '/'))
}else{
    load(paste(temp_path, 'fit.fixed', sep = '/'))
}

print(fit.fixed, pars = c('delta_fixed', 'gamma_fixed'),
      probs = c(.025, .975))
## All is good

######################################################################
## Comparison between the models fitted to the simulated data

ss = list('True hierarchical model' = as.data.frame(fit.sim),
          'Non-hierarchical model' = as.data.frame(fit.aggr))
## Here we are calculating the posterior samples for gamma for the
## second condition by adding the effect of the condition on gamma and
## than using the appropriate transformation.
ss_ = ss
for(i in 1:2){
    ss_[[i]][grep('delta_fixed', names(ss_[[i]]))] = apply(ss_[[i]][grep('delta_fixed', names(ss_[[i]]))], 2, exp)
    ss_[[i]][grep('gamma_fixed', names(ss_[[i]]))][,8:14] = ss_[[i]][grep('gamma_fixed', names(ss_[[i]]))][,1:7] +
        ss_[[i]][grep('gamma_fixed', names(ss_[[i]]))][,8:14]
    ss_[[i]][grep('gamma_fixed', names(ss_[[i]]))][,1:7] = gamma_to_crit(ss_[[i]])
    ss_[[i]][grep('gamma_fixed', names(ss_[[i]]))][,8:14] = gamma_to_crit(ss_[[i]], 2)
    names(ss_[[i]])= gsub('delta_fixed', 'dprim_fixed', names(ss_[[i]]))
    names(ss_[[i]])= gsub('gamma_fixed', 'criteria_fixed', names(ss_[[i]]))
}
round(apply(ss_[[1]][, grep('criteria_fixed', names(ss_[[1]]))], 2, mean), 2)
round(apply(ss_[[2]][, grep('criteria_fixed', names(ss_[[2]]))], 2, mean), 2)
## Ok

## ## This will add the fixed effects averaged over the participants
## sf = as.data.frame(fit.fixed)
## ss[['Fixed participants effects model']] = ss[[2]]
## sdata = make_stan_data(data_sim_2, list(delta = ~ -1 + duration:order, gamma = ~ 1 + order))
## for(i in 1:sdata$X_delta_ncol){
##     ss[[3]][,sprintf('delta_fixed[%d]', i)] =
##         apply(sf[, paste('delta_fixed[',
##                          grep(colnames(sdata$X_delta)[i], colnames(sdata.fixed$X_delta)),
##                          ']', sep = '')], 1, mean)
## }
## colnames(sdata.fixed$X_gamma)[-grep('order', colnames(sdata.fixed$X_gamma))] = '(Intercept)'
## for(i in 1:sdata$X_gamma_ncol){
##     for(k in 1:7)
##     ss[[3]][,sprintf('gamma_fixed[%d,%d]', k, i)] =
##         apply(sf[, paste('gamma_fixed[', k, ',',
##                          grep(colnames(sdata$X_gamma)[i], colnames(sdata.fixed$X_gamma)),
##                          ']', sep = '')], 1, mean)
## }

bias_plot = function(true_values, ss){
    df = data.frame(model = rep(names(ss), each = length(grep('fixed', names(ss[[1]])))),
                    par = c(names(ss[[1]])[grep('fixed', names(ss[[1]]))],
                            names(ss[[2]])[grep('fixed', names(ss[[2]]))]))
    if(length(grep('delta', levels(df$par))) > 0){
        df$par = factor(df$par, c(grep('delta', levels(df$par), value = T),
                                  grep('gamma', levels(df$par), value = T)))
    }else{
        df$par = factor(df$par, c(grep('dprim', levels(df$par), value = T),
                                  grep('criteria', levels(df$par), value = T)))
    }
    for(m in levels(df$model)){
        df$point.est[df$model == m] = apply(ss[[m]][,as.character(df$par[df$model == m])], 2, mean)
        df$ci.lo[df$model == m] = apply(ss[[m]][,as.character(df$par[df$model == m])], 2, function(x)quantile(x, .025))
        df$ci.hi[df$model == m] = apply(ss[[m]][,as.character(df$par[df$model == m])], 2, function(x)quantile(x, .975))
    }
    df$true = true_values[as.character(df$par)]
##    df$par.type = NA
##    df$par.type[grep('delta', as.character(df$par))] = 'delta'
##    df$par.type[grep('gamma', as.character(df$par))] = 'gamma'
    ggplot(df, aes(par, point.est - true)) +
        geom_abline(intercept = 0, slope = 0, lty = 2, alpha = .5) +
        geom_point() +
        geom_errorbar(aes(ymin = ci.lo - true, ymax = ci.hi - true)) +
        facet_grid(~ model) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(y = 'Point and interval estimates centered on true values',
             x = 'Model parameter') +
        theme_minimalist()
}

s = as.data.frame(fit)
true_values = apply(s[, grep('fixed', names(s))], 2, mean)
(p1 = bias_plot(true_values, ss))
ggsave('true_vs_nonhier.pdf')
## The estimates based on overly aggregated data are useless

true_values_ = apply(apply(s[, grep('delta_fixed', names(s))], 2, exp), 2, mean)
names(true_values_) = gsub('delta', 'dprim', names(true_values_))
## Translating between the posterior gamma and the criteria for the
## true model that was used to simulate the data
s_ = s
s_[grep('gamma_fixed', names(s_))][,8:14] = s_[grep('gamma_fixed', names(s_))][,1:7] +
    s_[grep('gamma_fixed', names(s_))][,8:14]
true_values_ = c(true_values_, apply(gamma_to_crit(s_), 2, mean))
true_values_ = c(true_values_, apply(gamma_to_crit(s_, 2), 2, mean))
(p2 = bias_plot(true_values_, ss_))
ggsave('true_vs_nonhier_sdt.pdf')

p3 = grid.arrange(p1, p2, ncol = 2)
ggsave('true_vs_nonhier_both.pdf', p3, width = 13, height = 7)
ggsave('true_vs_nonhier_both.png')
ggsave('true_vs_nonhier_both.tiff')

res = as.data.frame(cbind(t(apply(ss[[2]][,grep('delta_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])), delta_fixed))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## .5

res = as.data.frame(cbind(t(apply(ss[[2]][,grep('gamma_fixed', names(s2))], 2, function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)])),
                          as.vector(t(gamma_fixed))))
names(res) = c('lo', 'fit', 'hi', 'true')
res$in_ci = (res$true >= res$lo) & (res$true <= res$hi)
res
mean(res$in_ci)
## .28
cbind(sum(res$in_ci), nrow(res))
## 4 14

##!?
##
##res$Condition = rep(1:(ncol(model.matrix(fixed$gamma, data_sim_2$data))), each = ncol(data_sim_2$counts_1) - 1)

## Let's compare the posterior SDs
aggr.posterior.sd = apply(as.data.frame(fit.aggr), 2, sd)
true.posterior.sd = apply(as.data.frame(fit.sim), 2, sd)
sd.ratios = true.posterior.sd[grep('fixed', names(true.posterior.sd))] /
  aggr.posterior.sd[grep('fixed', names(aggr.posterior.sd))]
mean(sd.ratios)
## 3.122806

######################################################################
## Approximate normality of random effects distributions

(p = ggplot(data.frame(delta = as.vector(t(delta_random)), i = c('32 ms', '64 ms')), aes(sample = delta)) +
     stat_qq() +
     facet_wrap(~i) +
     ylab('delta sample quantile') +
     xlab('Theoretical normal quantile') + theme_minimalist())
ggsave('delta_qq_plot.pdf', p)
## this is perfectly acceptable

shapiro.test(delta_random[,1])
## W = 0.9884, p-value = 0.918
shapiro.test(delta_random[,2])
## W = 0.96955, p-value = 0.2512
gamma_normality = data.frame(k = 1:7, statistic = NA, p.value = NA)
for(i in 1:ncol(gamma_random_mat)){
  res = shapiro.test(gamma_random_mat[,i])
  gamma_normality$statistic[i] = res$statistic
  gamma_normality$p.value[i] = res$p.value
}
round(gamma_normality, 3)
##   k statistic p.value
## 1 1     0.975   0.411
## 2 2     0.992   0.984
## 3 3     0.958   0.091
## 4 4     0.979   0.562
## 5 5     0.978   0.519
## 6 6     0.945   0.028
## 7 7     0.985   0.784

gamma_rnd = gamma_random_mat[,1]
for(i in 2:ncol(gamma_random_mat))
  gamma_rnd = c(gamma_rnd, gamma_random_mat[,i])
k = rep(1:ncol(gamma_random_mat), each = nrow(gamma_random_mat))
(p = ggplot(data.frame(gamma = gamma_rnd, k = k), aes(sample = gamma)) + stat_qq() + facet_wrap(~k) +
     ylab('gamma sample quantile') +
     xlab('Theoretical normal quantile') +
     theme_minimalist())
ggsave('gamma_qq_plot.pdf', p)
## this is also perfectly acceptable

######################################################################
## Fitting the meta-d' model to simulated data

simulate_metad_data = function(span = 2, d1 = 2, d2 = 1, nof_crit = 10 + 1, bias = 0, crit_scale = 1){
    ## outermost criteria spread = d' + 2 * span
    if((nof_crit / 2) == floor(nof_crit / 2)){
        print("nof_crit must be ann odd number, increasing by 1")
        nof_crit = nof_crit + 1
    }
    nof_r = nof_crit + 1
    ## the main criterion
    k = floor(nof_crit / 2) + 1
    crit = (seq(-d1/2 - span, d1/2 + span, length.out = nof_crit) * crit_scale) + bias
    crit_ex = c(-Inf, crit, Inf)
    d = expand.grid(r = 1:(nof_crit + 1), s = c(-0.5, .5))
    ## p1 contains combined response probabilities for the metad' = d'
    ## case, p2 contains combined response probabilities for the metad' !=
    ## d' case.
    d$p1 = d$p2 = NA
    for(i in 1:nrow(d)){
        d$p1[i] = pnorm(crit_ex[d$r[i]+1] - d$s[i] * d1) - pnorm(crit_ex[d$r[i]] - d$s[i] * d1)
        if(d$r[i] <= nof_r / 2){
            d$p2[i] = (pnorm(crit_ex[d$r[i]+1] - d$s[i] * d2) - pnorm(crit_ex[d$r[i]] - d$s[i] * d2)) /
                pnorm(crit[k] - d$s[i] * d2) * pnorm(crit[k] - d$s[i] * d1)
        }else{
            d$p2[i] = (pnorm(crit_ex[d$r[i]+1] - d$s[i] * d2) - pnorm(crit_ex[d$r[i]] - d$s[i] * d2)) /
                pnorm(-(crit[k] - d$s[i] * d2)) * pnorm(-(crit[k] - d$s[i] * d1))
        }
    }
    d$stim = as.factor(as.numeric(as.factor(d$s)))
    d
}

d = simulate_metad_data(nof_crit = 100)
## save(d, file = '~/temp/d')
## load('~/temp/d')

## Let's see how the combined response distributions differ
ggplot(d[d$r > 1 & d$r < max(d$r),], aes(r, p1, group = stim, color = stim)) + geom_line() +
    xlab("Combined response distributions for the meta-d' = d' case") + ylab('Combined response probability')
ggplot(d[d$r > 1 & d$r < max(d$r),], aes(r, p2, group = stim, color = stim)) + geom_line() +
    xlab("Combined response distributions for the meta-d' != d' case") + ylab('Combined response probability')
## It looks weird, but this is exactly what the meta-d' predicts in this case

## Let's simulate some data
d = simulate_metad_data()
N = 1000
ds = expand.grid(stim = 1:2, i = 1:(N/2))
ds$r = NA
for(i in 1:N)
    ds$r[i] = rMultinom(t(d$p2[d$stim == ds$stim[i]]), 1)
## save(ds, file = '~/temp/ds')
## load('~/temp/ds')
ds$dec = as.numeric(ds$r > 6) + 1
ds$rating = ds$r
ds$rating[ds$r > 6] = ds$rating[ds$r > 6] - 6
ds$rating[ds$r <= 6] = 7 - ds$rating[ds$r <= 6]

adata = aggregate_responses(ds, 'stim', 'r')
fixed = list(delta = ~ 1, gamma = ~ 1)

if(fresh_start){
    fit.metad = stan(model_code = make_stan_model(metad = TRUE),
                     data = make_stan_data(adata, fixed, metad = TRUE),
                     pars = c('delta_fixed', 'gamma_fixed',
                              ## we need counts_new for plotting
                              'counts_new'),
                     iter = 4000,
                     chains = 3)
    print(fit.metad, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))
    save(fit.metad, file = paste(temp_path, 'fit.metad', sep = '/'))
}else{
    load(paste(temp_path, 'fit.metad', sep = '/'))
}

if(fresh_start){
    fit.sdt = stan(model_code = make_stan_model(),
                   data = make_stan_data(adata, fixed),
                   pars = c('delta_fixed', 'gamma_fixed',
                            ## we need counts_new for plotting
                            'counts_new'),
                   iter = 4000,
                   chains = 3)
    print(fit.sdt, probs = c(.025, .957), pars = c('delta_fixed', 'gamma_fixed'))
    save(fit.sdt, file = paste(temp_path, 'fit.sdt', sep = '/'))
}else{
    load(paste(temp_path, 'fit.sdt', sep = '/'))
}

## Let's see how it performs
s = as.data.frame(fit.metad)
round(apply(exp(s[,1:2]), 2, mean), 2)
## The dprim and meta-d' parameters are correctly recovered
round(HPDinterval(as.mcmc(exp(s[,1]) - exp(s[,2]))), 2)
res = as.data.frame(cbind(round(apply(gamma_to_crit(s), 2, mean), 1), crit))
names(res) = c('estimated', 'true')
ggplot(res, aes(true, estimated)) + geom_point() + geom_abline(intercept = 0, slope = 1) +
    xlab('True criteria') + ylab('Estimated criteria')
## The criteria are correctly recovered

dm = as.matrix(read.csv('maniscalco.csv'))[,-3]
res = rbind(dm[c(1:2,4:8,3,9:13)],
            c(apply(s[,1:2], 2, function(x)mean(exp(x))),
              apply(gamma_to_crit(s), 2, mean)))
colnames(res) = c('dprim', 'meta-dprim', paste('c', 1:11, sep = ''))
round(res[,1:2], 2)
round(res[,3:13], 2)
crit_t = as.data.frame(t(res[,3:13]))
names(crit_t) = c('Maniscalco', 'bhsdtr')
ggplot(crit_t, aes(bhsdtr, Maniscalco)) + geom_point() + geom_line() + labs(title = 'Criteria estimates')
ggsave('criteria_maniscalco_vs_moi.pdf')

plot_sdt_fit(fit.metad, adata)
ggsave('metad_roc_fit.pdf')
## Fits nicely
plot_sdt_fit(fit.sdt, adata)
ggsave('metad_roc_sdt_fit.pdf')
## Also seems to fit nicely
plot_sdt_fit(fit.metad, adata, type = F)
ggsave('metad_resp_fit.pdf')
plot_sdt_fit(fit.sdt, adata, type = F)
ggsave('metad_resp_sdt_fit.pdf')
## The lack of fit is very apparent here

######################################################################
## Here we demonstrate the consequences of using Fleming's criteria
## parametrization (HMeta-d')

gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
adata = aggregate_responses(gabor[gabor$duration == '64' & gabor$order == 'DR',], 'stim', 'r', c('id', 'order'))
fixed = list(delta = ~ 1, gamma = ~ 1)
random = list(list(group = ~ id, delta = ~ 1, gamma = ~ 1))
model = make_stan_model(random)
sdata = make_stan_data(adata, fixed, random)

if(fresh_start){
    fit.2 = stan(model_code = model,
                 data = sdata,
                 pars = c('delta_fixed', 'gamma_fixed',
                          'delta_sd_1', 'gamma_sd_1',
                          'delta_random_1', 'gamma_random_1',
                          'Corr_delta_1', 'Corr_gamma_1',
                          'counts_new'),
                 iter = 8000,
                 chains = 4)
    fit.2
    save(fit.2, file = paste(temp_path, 'fit.2', sep = '/'))
}else{
    load(paste(temp_path, 'fit.2', sep = '/'))
}

fname = paste(temp_path, 'adata_sim', sep = '/')
if(fresh_start){
    adata_sim = simulate_from_fit(fit.2, adata, fixed, random)
    save(adata_sim, file = fname)
}else{
    load(fname)
}

## We will fit two versions of the model to the same dataset. The
## models differ in the way the criteria are parametrized. We will
## make it as easy as possible for the incorrect model to fit the data
## by using strongly informative priors (based on the true model
## parameters) for the criteria.
sdata_c = sdata

##! WARNING You have to download the functions from the HMeta-d'
##! github repository yourself.
source('Function_metad_group.R')
nR_S1 = as.data.frame(t(adata_sim$counts[adata_sim$stimulus == 1,]))
nR_S2 = as.data.frame(t(adata_sim$counts[adata_sim$stimulus == 2,]))
fit.hm1 = metad_group(list(nR_S1), list(nR_S2))
eff = effectiveSize(fit.hm1)
smr = summary(fit.hm1)$statistics
round(smr[-grep('raw', rownames(smr)),], 2)
