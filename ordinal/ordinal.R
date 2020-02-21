ggs = function(f)ggsave(f, width = 10, height = 10, dpi = 200)

######################################################################
## KOD 1

library(ggplot2)
library(bhsdtr)
library(Hmisc)
library(coda)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
## source('../tests/code.R')

######################################################################
## Ordinal model for responses 'noise'

## r = 1 ~ absolutely clear image left, r = 4 ~ no experience left, r
## = 5 ~ no experience right, ... r = 8 ~ absolutely clear image right
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
d1 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
## Only the left-tilted gabor patches
d1 = d1[d1$stim == 0,]

## This order-preserving link function represents the main criterion
## directly as a free paramter and all the other criteria as log of
## distance between the criterion and the criterion adjacent to it
## which is closer to the main criterion
link = 'log_distance'
adata0 = aggregate_responses(d1, response = 'r', variables = c('id'))
fixed0 = list(eta = ~ 1, gamma = ~ 1)
## Each participant has a unique patterh of the thresholds
random0 = list(list(group = ~ id, gamma = ~ 1))
model0 = make_stan_model(random, model = 'ordinal', gamma_link = link)
sdata0 = make_stan_data(adata0, fixed0, random0, model = 'ordinal', gamma_link = link)

## The model is saturated, so we know that it will fit perfectly.
o = optimizing(stan_model(model_code = model), sdata, init = 0)
if(o$return_code)warning('not ok')

sdt_ml_plot(model0, adata0, sdata0)
ggs('ordinal_fit_plot.jpg')

## The pattern of the thresholds seems to differ between the participants
gamma_fixed = o$par[grep('gamma_fixed\\[', names(o$par))]
gamma_random = matrix(o$par[grep('gamma_random_1', names(o$par))], nrow = sdata$group_1_size)
for(i in 1:nrow(gamma_random))
    gamma_random[i,] = gamma_random[i,] + gamma_fixed
criteria = gamma_to_crit(gamma_random, beta_index = NULL, gamma_link = link)
df = expand.grid(criterion = 1:7, id = 1:nrow(criteria))
df$estimate = as.vector(t(criteria))
ggplot(df, aes(id, estimate, group = id)) + geom_line() + geom_point()

## Now we will fit a model which is similar to an SDT model with two
## important differences. We will not assume that the answer "signal"
## cannot be less likely for the stimulus "signa" than for the
## stimulus "noise", i.e., we will not assume that d' is non-negative,
## and we will not assume that the pattern of the thresholds is the
## same for "noise" and "signal" stimuli. In this model each
## participant may have two different sets of thresholds, one for
## "noise" and one for "signal" stimuli. Since the thresholds are
## represented as random effects, we will account for the fact that
## thresholds are possibly correlated; For example, if some threshold
## is higher then the threshold above (below) it may also be higher
## (lower). This is an important assumption, without it the important
## assumption that the data points are independent given the model
## will be seriously violated. If this assumption is seriously
## violated, for example by using a fixed anova model to analyze
## repeated measures data, the uncertainty in the data becomes
## strongly underestimated and the whole inference thing fails
## miserably.

## Here we add the effect of stimulus on the thresholds/gamma
## vector. We will use separate intercepts and slopes parametrization,
## i.e., the thresholds will be estimated for ech stimulus class
## instead of estimating them for "noise" stimuli and estimating the
## difference in the thresholds between "signal" and "noise".
fixed0 = list(eta = ~ 1, gamma = ~ -1 + stim)
## Here we assume that each participant has two sets of thresholds for
## the "noise" and "signal" stimuli.
random0 = list(list(group = ~ id, gamma = ~ -1 + stim))
d2 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
## We have to make sure that stim is a factor for the separate
## intercepts and slopes parametrization to work as intended.
d2$stim = as.factor(d2$stim)
## Note that because for an ordinal model the stimulus ('stim')
## variable is not special we have to include it in the additional
## variables vector
adata0 = aggregate_responses(d2, response = 'r', variables = c('id', 'stim'))
sdata0 = make_stan_data(adata0, fixed0, random0, model = 'ordinal', gamma_link = link)
model0 = make_stan_model(random0, model = 'ordinal', gamma_link = link)

fit0 = stan(model_code = model0, data = sdata0,
            pars = c('gamma_fixed', 'gamma_random_1', 'gamma_sd_1'),
            init_r = .5,
            iter = 5000,
            chains = 4)
## saveRDS(fit0, '~/windows/temp/ordinal_fit_0.rds')

## We will not look at the eta parameter since it was fixed at 0
print(fit0, probs = c(.025, .975), pars = c('gamma_fixed', 'gamma_sd_1'))
## neff and Rhat look good
s0 = as.data.frame(fit0)

## Here we translate between the gamma and the criteria vectors (fixed
## effects = sample average) for each stimulus
crit1 = gamma_to_crit(s0, 1, gamma_link = link)
crit2 = gamma_to_crit(s0, 1, beta_index = 2, gamma_link = link)

sumfun = function(x)c(quantile(x, c(.025, .975)), mean(x))[c(1,3,2)]
crit1mean = apply(crit1, 2, sumfun)
round(t(crit1mean), 2)
##       2.5%       97.5%
## [1,] -1.82 -1.30 -0.89
## [2,] -1.37 -0.95 -0.59
## [3,] -0.66 -0.34 -0.05
## [4,]  0.49  0.73  0.97
## [5,]  0.93  1.21  1.54
## [6,]  1.16  1.48  1.87
## [7,]  1.83  2.31  3.16

crit2mean = apply(crit2, 2, sumfun)
round(t(crit2mean), 2)

## Note that the average criterion position is the same (in terms of
## the mean average and the credible intervals) as the fourth
## criterion which is the same as the main decision criterion which is
## the same as the fourth element of the gamma vector. That is because
## the log_distance link function represents this main decision
## criterion as an unconstrained parameter, whereas all the other
## thresholds are represented in the gamma parameter vector as log of
## distance between the respective threshold and the threshold
## adjacent to it, e.g., gamma_5 = log(criterion5 - criterion 4),
## gamma_3 = log(criterion4 - criterion3), etc.

## We can recover an estimate of the average d' by calculating the
## distance between the main (or average) criteria for "signal" and
## "noise" stimuli":

## Note that shifting the criteria to the right is equivalent to
## shifting the underlying distribution of "latent values" (or
## evidence) to the left, so we have to subtract the average position
## of the criteria for the "signal" stimuli from the average position
## of the criteria for the "noise" stimuli
dprim0 = sumfun(s0[, 'gamma_fixed[4,1]'] - s0[, 'gamma_fixed[4,2]'])
round(dprim0, 2)
## 2.5%       97.5% 
## 0.96  1.31  1.68 

## In an SDT model we assume that the stimulus class affects only the
## d' parameter, and that both criteria (gamma) and d' may differ
## between the participants (or items, etc.). In bhsdtr d' is modelled
## as log(d') = delta, which forces it to be non-negative. This is an
## important assumption of Signal Detection Theory - if this
## assumption is violated, for example because of the response
## reversal in some participants or conditions, than an SDT model is
## not valid and a non-trivial generalization of SDT has to be
## used. Otherwise this is the same general hierarchical ordinal
## regression model.
fixed = list(delta = ~ 1, gamma = ~ 1)
random = list(list(group = ~ id, delta = ~ 1, gamma = ~ 1))
## Because of the importance of the stimulus class ("noise" and
## "signal") in SDT models we have to indicate which variable
## represents this information
adata = aggregate_responses(d2, stimulus = 'stim', response = 'r', variables = c('id'))
sdata = make_stan_data(adata, fixed, random, model = 'sdt', gamma_link = link)
model = make_stan_model(random, model = 'sdt', gamma_link = link)

fit1 = stan(model_code = model, data = sdata,
            pars = c('delta_fixed', 'delta_sd_1', 'gamma_fixed', 'gamma_sd_1', 'delta_random_1', 'gamma_random_1',
                     ## We will need this for the plot_sdt_fit function
                     'counts_new'),
            init_r = .5,
            iter = 5000,
            chains = 4)
## saveRDS(fit1, '~/windows/temp/ordinal_fit_1.rds')

print(fit1, probs = c(.025, .975), pars = c('delta_fixed', 'gamma_fixed', 'delta_sd_1', 'gamma_sd_1'))
## neff and Rhat look good

s1 = as.data.frame(fit1)

## So how well do the two d' estimates match?
##
## Here we translate between the delta and the d' parameters by using
## exponentiation
dprim1 = sumfun(exp(s1[,'delta_fixed[1,1]']))

round(cbind(dprim0, dprim1), 2)
## They estimates match poorly. In fact, the 2.5 percentile of the
## posterior distribution of d' is above the 97.5 percentile of the
## posterior distribution of the average shift in the criteria.

## What about the individual d' estimates? To calculate the individual
## shifts in the thresholds/criteria we have to add the random effects
## to the fixed effect:
dprim0i = (s0[['gamma_fixed[4,1]']] + s0[, grep('gamma_random_1\\[.+,4,1]', names(s0))]) -
    (s0[['gamma_fixed[4,2]']] + s0[, grep('gamma_random_1\\[.+,4,2]', names(s0))])
dprim0i = as.data.frame(t(apply(dprim0i, 2, sumfun)))
## Note that we have to apply the exponential function first before
## computing the posterior means and quantiles (non-linearity)
dprim1i = as.data.frame(t(apply(exp(s1[['delta_fixed[1,1]']] + s1[, grep('delta_random_1', names(s1))]), 2, sumfun)))
names(dprim0i) = names(dprim1i) = c('lo', 'mean', 'hi')

dprim0i$model = 'ordinal'
dprim0i$i = 1:nrow(dprim0i) - .1
dprim1i$model = 'sdt'
dprim1i$i = 1:nrow(dprim1i) + .1
ggplot(rbind(dprim0i, dprim1i), aes(i, mean, group = model, color = model)) +
    geom_point() + geom_errorbar(aes(ymin = lo, ymax = hi)) +
    geom_abline(intercept = 0, slope = 0) +
    labs(title = sprintf('Correlation between individual d\' estimates: r = %.2f', cor(dprim0i$mean, dprim1i$mean)))    
ggs('individual_dprim_ordinal_vs_sdt.jpg')
## The point and interval estimates clearly match quite well. However,
## note that according to the ordinal model some credible intervals
## around individual d' point estimates include 0, which is
## contradicts one assumption of the SDT model. We will address this
## issue later.

## Maybe the SDT model does not fit the data well?
plot_sdt_fit(fit1, adata, c('id'), type = 'response', bw = F)
## The SDT model seems to fit very well in a sense that the observed
## distributions of responses are within the 95% posterior predictive
## intervals
ggs('fit1_plot.jpg')

## The discrepancy between the sample average d' estimates from both
## models is not a consequence of allowing for negative d' values in
## the ordinal model. This discrepancy is a consequence of the assumed
## distribution of individual shifts / d' values. In the ordinal model
## we assumed that the shifts (= d') are normally distributed whereas
## in the SDT model we assumed that log(d') are normally
## distributed. We can see if the assumption that log(d') are normally
## distributed is approximately correct in our dataset. First let's
## take a look at estimates based on the SDT model:

ggplot(data.frame(delta = log(dprim1i$mean)), aes(sample = delta)) +
    stat_qq() +
    ylab('Sample quantile of log(d\') = delta random effects') +
    xlab('Theoretical normal quantile') +
    labs(title = 'SDT model')
ggs('qq_delta_sdt.jpg')
## The individual log(d') estimates based on the SDT model seem to be
## approximately normally distributed. We cannot plot log(d')
## estimates based on the ordinal model, because according to this
## model one participant has a negative d' value. However, we can see
## if the d' point estimates based on the ordinal model are
## approximately normally distributed:

ggplot(data.frame(delta = dprim1i$mean), aes(sample = delta)) +
    stat_qq() +
    ylab('Sample quantile of d\' random effects') +
    xlab('Theoretical normal quantile') +
    labs('Ordinal model')
ggs('qq_dprim_ordinal.jpg')
## They are clearly not. This was to be expected; In typical binary
## classification tasks below-chance performance is either a
## consequence of response reversal in a small number of participants
## or, in case of low but non-negative true d' values, of sampling
## error.

## The average / fixed effect d' in each model is estimated as the
## mean of different distributions, but in the ordinal model shifts
## are (incorrectly) assumed to normally distributed and in the SDT
## model delta = log(d') random effects are assumed to be normally
## distributed. We can easily account for this:

round(cbind(sumfun(exp(s1[['delta_sd_1[1,1]']]^2 / 2 + s1[['delta_fixed[1,1]']])),
            dprim0), 2)
## 2.5%  0.90   0.96
##       1.27   1.31
## 97.5% 1.89   1.68

## Now we see that the two models give very similar point and interval
## estimates of sample-average d'.

## To summarize, the SDT model is much simpler (it has only one set of
## thresholds for each participants), it fits the data well at the
## individual level, and it correctly assumes that individual d' /
## shift values are log-normally distributed. SDT is a clear winner,
## right? Wrong.

## The SDT model is not only a statistical model, it has causal
## assumptions too. Because the stimulus class was chosen at random on
## every trial every estimate of the statistical effect of the
## stimulus class on any parameter of the model is an unbiased
## estimate of the total causal effect of the stimulus class on this
## parameter. Moreover, in such designs every statistical effect of
## the participant factor is an estimate of the property of the
## participants. This means that both the individual d' parameters and
## the individual thresholds are participant-, not item-,
## parameters. So far so good.

## However, in the SDT model we have also assumed that the stimulus
## class can *only* affect the distance between the evidence
## distribution means (= d'), whereas in the ordinal model we allowed
## for the possibility of the stimulus class having an effect on the
## thresholds. It may seem odd to you to allow for this, but there is
## nothing in the design of the study that justifies the selectivity
## assumption. The estimates of the thresholds are based on the data
## that were generated after the stimulus was presented, and so,
## theoretically, they could be affected by the stimulus class.

## The selectivity assumption is central to Signal Detection Theory,
## and I do not recall even a single paper or a book where this
## assumption was tested. It does not matter if an SDT model with
## selectivity assumption fits the data well, because if the effect of
## the stimulus class is not selective, the model is simply false. In
## particular, when the selectivity assumption is false an SDT model
## cannot deconfound sensitivity / d' from bias / threshold effects.

## We can at least try to test this assumption. We will not alter the
## ordinal model code to account for the non-negativity of d' (and the
## log-normal distribution of random d' effects), because this would
## be complicated. As we already know, even without this correction
## the estimates of participants' shifts are similar to the d'
## estimates, and we only want to give an example of how this problem
## could be addressed.

## Note that we do not use the separate intercepts and slopes
## parametrization (i.e., ~ -1 + stim) here, because we are very much
## be interested in the *effect* of the stimulus class on the pattern
## of the thresholds.
fixed2 = list(eta = ~ 1, gamma = ~ stim)
random2 = list(list(group = ~ id, gamma = ~ stim))
adata2 = aggregate_responses(d2, response = 'r', variables = c('id', 'stim'))
sdata2 = make_stan_data(adata2, fixed2, random2, model = 'ordinal', gamma_link = link)
model2 = make_stan_model(random0, model = 'ordinal', gamma_link = link)

fit2 = stan(model_code = model2, data = sdata2,
            pars = c('gamma_fixed', 'gamma_random_1', 'gamma_sd_1'),
            init_r = .5,
            iter = 5000,
            chains = 4)
## saveRDS(fit2, '~/windows/temp/ordinal_fit_2.rds')
print(fit2, probs = c(.025, .975), pars = c('gamma_fixed', 'gamma_sd_1'))

s2 = as.data.frame(fit2)

## Here we see the average shift / d' estimated as a mean of normally
## distributed d' (not log(d')) random effects again
round(cbind(sumfun(-s2[['gamma_fixed[4,2]']]), dprim0), 2)
##            dprim0
## 2.5%  0.97   0.96
##       1.29   1.31
## 97.5% 1.62   1.68

## Ok, if we did simplify things significantly, but it still makes
## sense to take a look at the standard deviations of the other gamma
## random effects, because these random effects directly correspond to
## the *pattern* of the thresholds.

round(HPDinterval(as.mcmc(s2[, grep('gamma_fixed', names(s2))])), 2)
##                  lower upper
## gamma_fixed[1,1] -1.80 -0.45
## gamma_fixed[2,1] -0.88 -0.22
## gamma_fixed[3,1] -0.20  0.33
## gamma_fixed[4,1]  0.50  0.92
## gamma_fixed[5,1] -1.07 -0.37
## gamma_fixed[6,1] -1.76 -0.91
## gamma_fixed[7,1] -0.77 -0.06
## gamma_fixed[1,2] -0.74  1.67
## gamma_fixed[2,2] -1.37  0.22
## gamma_fixed[3,2] -0.37  0.12
## gamma_fixed[4,2] -1.62 -0.98
## gamma_fixed[5,2]  0.04  0.62
## gamma_fixed[6,2]  0.01  1.00
## gamma_fixed[7,2] -0.40  0.41

## Among the parameters that represent the stimulus induced
## differences between the gamma parameters (i.e., gamma_fixed[.,2])
## the gamma_fixed[4,2] parameter does not seem to be the only one
## that is non-zero. This indicates that the *pattern* of the
## thresholds is affected by the stimulus class, not only the constant
## shift. Interestingly, the random effects' standard deviations tell
## a different story:

round(HPDinterval(as.mcmc(s2[, grep('gamma_sd_1', names(s2))])), 2)
##                 lower upper
## gamma_sd_1[1,1]  0.45  1.87
## gamma_sd_1[2,1]  0.31  0.98
## gamma_sd_1[3,1]  0.44  0.85
## gamma_sd_1[4,1]  0.33  0.67
## gamma_sd_1[5,1]  0.48  1.00
## gamma_sd_1[6,1]  0.02  0.79
## gamma_sd_1[7,1]  0.00  0.60
## gamma_sd_1[1,2]  0.00  2.55
## gamma_sd_1[2,2]  0.00  1.84
## gamma_sd_1[3,2]  0.00  0.53
## gamma_sd_1[4,2]  0.53  1.02
## gamma_sd_1[5,2]  0.00  0.51
## gamma_sd_1[6,2]  0.02  0.97
## gamma_sd_1[7,2]  0.00  0.66

## Even though the average pattern of the thresholds seems to be
## different for the two stimulus classes, this difference in pattern
## seems to be fairly constant accross the participants. This is
## *really* weird, since the pattern of the thresholds for the "noise"
## stimuli clearly varies with participants - none of the 95% HPD
## intervals for the gamma[.,1] random effects' standard deviations
## includes zero. In terms of standard deviations of posterior
## distributions of standard deviations (phew) most are in fact quite
## far from zero.

## If only we could do a better job and use the model wich respects
## the non-negativity assumption and accounts for the log-normal
## distribution of d' random effects. Is it so complicated? Let's
## index the rows of our dataset by i. Because of the way the model
## was parametrized...

head(sdata2$X_gamma)

## when stim = "noise" (= 0) gamma[4]_i = gamma_fixed[4,1] +
## gamma_random_1[g][4,1], where g_i is the participant index in row
## i. We can leave this part of the model intact, since this is the
## main criterion position, which is an unconstrained parameter, but
## also, as we can see below, the distribution of the main criterion
## random effects is approximately normal:

ggplot(data.frame(gamma_4 = apply(s2[, grep('gamma_random_1\\[.+,4,1]', names(s2))], 2, mean)),
       aes(sample = gamma_4)) +
    stat_qq() +
    ylab('Sample quantile of main criterion random effects') +
    xlab('Theoretical normal quantile') +
    labs('Ordinal model')
ggs('qq_gamma_4_ordinal.jpg')

## However, when stim = "signal" (= 1), gamma[4]_i = gamma_fixed[4,1]
## + gamma_random_1[g][4,1] + gamma_fixed[4,2] +
## gamma_random_1[g][4,2]. We want this other part to represent log of
## average distance between the thresholds = shift = d'.

model3_ = make_stan_model(random0, model = 'ordinal', gamma_link = link)
model3 = gsub('criteria\\[Kb2] = gamma\\[Kb2];',
              'if(X_gamma[n,2] == 0){
                  criteria[Kb2] = gamma[Kb2];
               }else{ 
                  criteria[Kb2] = gamma[Kb2] - (gamma_fixed[Kb2,2] + gamma_random_1[group_1[n]][Kb2,2]) - exp(gamma_fixed[Kb2,2] + gamma_random_1[group_1[n]][Kb2,2]);
               }',
              model3_)

## This way the gamma[4,1] parameter represents the position of the
## main criterion relative to the midpoint between the evidence
## distribution means, and gamma[4,2] = delta = log(d').

fit3 = stan(model_code = model3, data = sdata2,
            pars = c('gamma_fixed', 'gamma_random_1', 'gamma_sd_1'),
            init_r = .5,
            iter = 5000,
            chains = 4)
## saveRDS(fit3, '~/windows/temp/ordinal_fit_2.rds')
print(fit3, probs = c(.025, .975), pars = c('gamma_fixed', 'gamma_sd_1'))

s3 = as.data.frame(fit3)

## This is close. Not the same, but the model is different.
round(cbind(sumfun(exp(s3[['gamma_fixed[4,2]']])), dprim1), 2)
## 2.5%  0.83   0.63
##       1.12   0.92
## 97.5% 1.43   1.26

## What about the bias = main criterion relative the the midpoint
## between the main criteria?
round(cbind(sumfun(s3[['gamma_fixed[4,1]']] - exp(s3[['gamma_fixed[4,2]']]) / 2),
      sumfun(s1[['gamma_fixed[4,1]']])), 2)
##        [,1]  [,2]
## 2.5%  -0.04 -0.07
##        0.15  0.08
## 97.5%  0.35  0.24

## Also close.

## What we witness here is perhaps a historic moment. We have fitted a
## hierarchical SDT model wich respects the non-negativity assumption,
## the ordering of the thresholds, accounts for the log-normal
## distribution of the d' random effects, allows for random effects in
## both the d' and the individual thresholds while still respecting
## the assumptions (thanks to an order-preserving gamma link function
## and the delta log link function), and accounts for the correlations
## between the random effects. First such model was described in the
## bhsdtr preprint. This model however does not assume that the
## stimulus class affects only the d' parameter. We can now again see
## if there is any evidence of non-selectivity:

round(HPDinterval(as.mcmc(s3[, grep('gamma_fixed', names(s3))])), 2)
##                  lower upper
## gamma_fixed[1,1] -1.78 -0.44
## gamma_fixed[2,1] -0.91 -0.24
## gamma_fixed[3,1] -0.20  0.32
## gamma_fixed[4,1]  0.50  0.92
## gamma_fixed[5,1] -1.09 -0.40
## gamma_fixed[6,1] -1.77 -0.94
## gamma_fixed[7,1] -0.77 -0.07
## gamma_fixed[1,2] -0.75  1.74
## gamma_fixed[2,2] -1.35  0.19
## gamma_fixed[3,2] -0.38  0.09
## gamma_fixed[4,2] -0.18  0.36
## gamma_fixed[5,2]  0.05  0.61
## gamma_fixed[6,2]  0.04  1.02
## gamma_fixed[7,2] -0.37  0.44

round(HPDinterval(as.mcmc(s3[, grep('gamma_sd_1', names(s3))])), 2)
##                 lower upper
## gamma_sd_1[1,1]  0.45  1.82
## gamma_sd_1[2,1]  0.29  0.97
## gamma_sd_1[3,1]  0.44  0.84
## gamma_sd_1[4,1]  0.32  0.66
## gamma_sd_1[5,1]  0.49  1.01
## gamma_sd_1[6,1]  0.05  0.81
## gamma_sd_1[7,1]  0.00  0.60
## gamma_sd_1[1,2]  0.00  2.54
## gamma_sd_1[2,2]  0.00  1.78
## gamma_sd_1[3,2]  0.00  0.48
## gamma_sd_1[4,2]  0.40  0.86
## gamma_sd_1[5,2]  0.00  0.53
## gamma_sd_1[6,2]  0.00  0.96
## gamma_sd_1[7,2]  0.00  0.65

## The pattern of results is very similar with one notable
## exception. This time every gamma_sd_1[4,.] parameter except for
## gamma_sd_1[4,2] = log(d') may be zero. We still see evidence of the
## effect of the stimulus on the average position of the thresholds,
## but the effects for all the log-distances seem to be more or less
## the same for each participant.


## Interaction is a strange beast. For example, if you have an
## additive effect of two factors on the original scale of the
## dependent variable...

df = expand.grid(f1 = 1:2, f2 = 1:2)
df$y = df$f1 + df$f2
df$f1 = as.factor(df$f1)
df$f2 = as.factor(df$f2)
ggplot(df, aes(f1, y, group = f2, color = f2)) + geom_line()

## ...then after applying a non-linear transformation to the dependent
## variable you may obtain a significant interacive effect:

## So, let's have a look at the average criteria:

crit1 = gamma_to_crit(s3, gamma_link = link)
## Notice that we have to add the effect on gamma to gamma for
## stimulus = "noise" to obtain gamma for stimulus = "signal"
crit2 = gamma_to_crit(s3[, grep('gamma_fixed\\[.+,1]', names(s3))] + s3[, grep('gamma_fixed\\[.+,2]', names(s3))], gamma_link = link)

## Now we can directly test if the pattern may be different:
round(HPDinterval(as.mcmc((crit2 - apply(crit2, 1, mean)) - (crit1 - apply(crit1, 1, mean)))), 2)
## When we centered the criteria for the two classes of stimuli we no
## longer see any evidence for the differences in their relative
## positions. The only thing that seems to differ, and not by much, is
## the main criterion. Selectivity holds! And this is how you *begin*
## to test a hierarchical SDT model.

## An SDT model has other problematic assumptions. One is that the
## time-accuracy tradeof can be ignored. If particpants or conditions
## differ in the response time, than the differences in d' do not mean
## that sensitivity is different, because participants may simply take
## more time to answer in certain experimental conditions. You need
## something like the diffusion model to deal with this problem. If
## you combine the diffusion model with ratings, you also need
## something like a hierarchical diffusion model with ratings and an
## order-preserving link function. I plan to define and implement such
## models in the future.

## The other problems that need to be addressed have to do with
## non-linearity and aggregation. If you fit an SDT model, or any
## other non-linear model, to data aggregated over grouping factors,
## such as participants or items, because of non-linearity (which does
## not play well with aggregation or averaging) the point and interval
## estimates will be biased and the inference will be invalid. You can
## see how horribly biased the point *and* interval estimates of every
## SDT parameter become when you overly aggregate your data. I have
## already tested this on a bunch of datasets from the excellent
## Confidence database created by Doby Rahnev and the results of
## over-aggregation were dramatic in every single case.

## This means that every SDT analysis to this day that was done using
## data aggregated over participants, and almost all of them were, is
## invalid. In particular, to my knowledge there is no valid evidence
## that Unequal Variance SDT model fits the data better then the Equal
## Variance SDT model. I am working on a paper that makes this exact
## point rather clear. Note that here we were aggregating over trials,
## which means that, if the SDT parameters may differ between the
## trials, and they certainly do differ, the only question is by how
## much, then all the SDT estimates that I have presented to you are
## biased as well. There is a lot to be done before we can say that
## the humble SDT model deconfounds sensitivity from bias in any
## particular case.

## Now we that we did all this stuff with SDT models we can generate a
## bunch of interesting theoretical results about ordinal models and
## Item Response Theory.

sumfun(apply(crit2, 1, mean) - apply(crit1, 1, mean))


lm(apply(crit1, 2, mean) ~ apply(crit2, 2, mean))

## Unequal Variance SDT? We cannot introduce any new parameters to
## this model, because it is already saturated. What we can do is fit
## an Unequal Variance SDT model, estimate the ratio of the variances,
## fix the ratio and introduce it to the saturated ordinal model. If
## there is still evidence of non-selectivity, we can keep trying with
## different modifications. If non-selectivity disappears, we will
## have clear evidence that the Unequal Variance SDT model is better
## than the Equal Variance SDT model. We will not do that now. We have
## more interesting things to say about ordinal models and Item
## Response Theory.

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
