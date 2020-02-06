## -*- coding: utf-8 -*-

######################################################################
## WARNING The interface has changed since the Behavior Research
## Methods paper was accepted for publication. In the current version
## of the package every parameter (delta, gamma, theta - in uv sdt
## models) is represented by a matrix even when it is unidimensional
## in a given model. So, for example, in an SDT model delta fixed
## effects are represented by parameters delta_fixed[1,1], ...,
## delta_fixed[1,n], where n is the number of delta fixed effects
## (i.e., the number of columns in the fixed effects model matrix for
## delta).

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

## We will try evey kind of model with every link function
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

## What is the correlation between the estimates and the true values
## of the parameters?
round(cor(res$est, res$true, use = 'pairwise.complete.obs'), 2)
## 1

## A plot comparing point estimates to known true values
res$par_type = c('criteria', 'dprim')[(as.character(res$par) == 'dprim') + 1]
ggplot(res, aes(est, true, group = par_type, color = par_type)) +
    geom_abline(slope = 1, intercept = 0) +
    geom_point(size = 5) +
    facet_grid(model ~ link)
## All is good

######################################################################
## Fitting the model described in the BRM paper to real study data

## We will be using the gabor dataset provided with the package
data(gabor)
?gabor

## Combined response has to be calculated for this dataset
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)

## Aggregation without information loss - we keep keep the variability
## due to the participants and all the experimental conditions.
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'order', 'id'))

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
                        ## we need counts_new if we want to use the
                        ## plot_sdt_fit function
                        'counts_new'),
               ## iter = .5 helps a lot with rejected initial samples,
               ## the default range of initial samples in stan is (-2,
               ## 2), which is too wide for the parameters in our
               ## models.
               iter = .5,
               iter = 8000,
               chains = 4)
    save(fit, file = paste(temp_path, 'fit', sep = '/'))
}else{
    load(paste(temp_path, 'fit', sep = '/'))
}

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

## According to the fitted model, the criteria are placed fairly
## symmetrically
crit = gamma_to_crit(as.data.frame(fit))
round(apply(crit, 2, mean), 2)

######################################################################
## Fitting the single criterion hierarchical SDT model (just for fun)

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
## single-criterion SDT model is untestable (i.e., it is saturated)

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
## Fairly similar
##
##       2.5%  50% 95.7%
## [1,] -0.05 0.11  0.25
## [2,] -0.13 0.06  0.21

## I am not sure why this is here...
smr2$n_eff = as.integer(round(smr2$n_eff))
smr2 = smr2[-grep('random', rownames(smr2)),]
smr2 = smr2[-grep('counts_new', rownames(smr2)),]
print(xtable(round(smr2[-nrow(smr2),], 2)), file = 'fit_table_single_c.tex')

######################################################################
## Fitting the model to data simulated from itself to test if the
## model recovers known realistic parameter values

data = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))

## This function takes a fitted hierarchical model and uses the point
## estimates to simulate new responses. The sample size is the same as
## in the aggregated data object.
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

print(fit.sim, probs = c(.025, .957),
      pars = c('delta_fixed', 'gamma_fixed', 'delta_sd_1', 'gamma_sd_1'))
## All is fine

(p1 = plot_sdt_fit(fit.sim, adata_sim, c('duration', 'order')))
ggsave('roc_sim_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.sim, adata_sim, c('duration', 'order'), type = 'response'))
ggsave('response_sim_fit.pdf', p2)

######################################################################
## Comparing true values with the point estimates based on the true
## model

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
## The proportion of true delta_fixed values which are within the 95%
## credible intervals
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
## Fitting the simplified, non-hierarchical (i.e., false in this case)
## model to simulated data aggregated over participants (i.e., too
## much aggregation). This is how SDT models are often (mis)used.

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

## Excellent sampling
print(fit.aggr, pars = c('delta_fixed', 'gamma_fixed'),
      probs = c(.025, .975))

## Based on the ROC curve plots the false simplified model seems to
## fit the data very
(p1 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration')))
ggsave('roc_sim_aggr_fit.pdf', p1)
(p2 = plot_sdt_fit(fit.aggr, data_sim_2, c('order', 'duration'), type = 'response'))
ggsave('response_sim_aggr_fit.pdf', p2)

######################################################################
## Comparison between the true and the simplified models fitted to the
## simulated data

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
## The estimates based on overly aggregated data are COMPLETELY
## USELESS

true_values_ = apply(apply(s[, grep('delta_fixed', names(s))], 2, exp), 2, mean)
names(true_values_) = gsub('delta', 'dprim', names(true_values_))
## Translating between gamma and criteria for the true model that was
## used to simulate the data
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

## Let's compare the posterior SDs
aggr.posterior.sd = apply(as.data.frame(fit.aggr), 2, sd)
true.posterior.sd = apply(as.data.frame(fit.sim), 2, sd)
sd.ratios = true.posterior.sd[grep('fixed', names(true.posterior.sd))] /
  aggr.posterior.sd[grep('fixed', names(aggr.posterior.sd))]
mean(sd.ratios)
## 3.122806 is the ratio of the correct posterior SD to posterior SD
## based on the simplified model

######################################################################
## Demonstration of approximately normal distribution of delta and
## gamma random effects

(p = ggplot(data.frame(delta = as.vector(t(delta_random)), i = c('32 ms', '64 ms')), aes(sample = delta)) +
     stat_qq() +
     facet_wrap(~i) +
     ylab('delta sample quantile') +
     xlab('Theoretical normal quantile') + theme_minimalist())
ggsave('delta_qq_plot.pdf', p)
## this is perfectly acceptable

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
