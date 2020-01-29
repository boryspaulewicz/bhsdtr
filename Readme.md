# NEWS ALERT! Two new additional link functions for the criteria (experimental)

In the current version of the package there are three link functions
for the SDT criteria to choose from. One is the link function
described in the preprint - this is now called 'softmax'. This link
function (softmax followed by inverse normal CDF) is computationally
intensive and makes the task of specifying the priors for the gamma
vector difficult.

The unconstrained gamma vector can be mapped to the ordered criteria
vector in other ways. Note that the main criterion (the K / 2
threshold) considered in isolation is an uncostrained parameter. The
rest of the criteria can be represented as log-distances between
adjacent criteria or log-ratios of distances between adjacent
criteria. For example, the K/2 + 1 criterion can be represented as
log(c_<sub>K+1</sub> - c<sub>K/2</sub>). This general idea leads
to several simple solutions. One is:

the main criterion is unconstrained:

c_<sub>K/2</sub> = &gamma<sub>K/2</sub>

the criteria above the main criterion are represented as
log-distances, e.g.:

c_<sub>K/2+3</sub> = c_<sub>K/2+2</sub> + exp(&gamma<sub>K/2+3</sub>)

and similarly for the criteria below the main criterion, e.g.:

c_<sub>K/2-3</sub> = c_<sub>K/2-2</sub> - exp(&gamma<sub>K/2-3</sub>)

This is now the 'log-distance' gamma link function. The prior for
&&gamma<sub>K/2</sub> is now easy to specify, because this element of
the &gamma vector represents the position of the main criterion
relative to the midpoint between the evidence distribution means,
i.e., the value of 0 corresponds to no bias and the positive
(negative) values correspond to the tendency to respond 'noise'
(signal). The priors for all the other elements of the &gamma vector
are almost as easy to specify. For example, the assumption that the
average distance between the criteria is probably .5 can be
represented by setting the means of the priors for the &gamma vector
(except for &gamma<sub>K/2</sub>) at log(.5).

The other link function is now called 'log-ratio'. The K/2 element
again represents the main criterion, the &gamma<sub>K/2+1</sub>
element represents log(c<sub>K/2+1</sub> - c<sub>K/2</sub>), which I
like to call the 'spread' parameter. The &gamma<sub>K/2-1</sub>
element represents the log-ratio of distances (hence the name of the
link function), i.e., log((c<sub>K/2</sub> - c<sub>K/2-1</sub>) /
(c<sub>K/2+1</sub> - c<sub>K/2</sub>)) - this is the assymetry between
the lower and the upper spread of the criteria which are next to the
main criterion. The &gamma<sub>K/2+i</sub> where i > 1 represent the
ratios of distances, i.e., &gamma<sub>K/2+i</sub> =
log((c<sub>K/2+i</sub> - c<sub>K/2+i-1</sub>) / (c<sub>K/2+1</sub> -
c<sub>K/2</sub>)), and I like to call them 'upper consistency'
parameters. The &gamma<sub>K/2-i</sub> elements, where i > 1 are
'lower consistency' parameters, i.e., &gamma<sub>K/2-i</sub> =
log((c<sub>K/2-i+1</sub> - c<sub>K/2-1</sub>) / (c<sub>K/2</sub> -
c<sub>K/2-1</sub>)). The reasonable prior for the log-ratio parameters
has mean = log(1) = 0.

A broad class of meaningfully simplified models can be obtained just
by restricting the values of the elements of the &gamma vector. For
example, by fixing all the ratios at log(1) = 0 we get, as a special
case, the parsimonious SDT model as described in this great
[paper](https://link.springer.com/article/10.3758/s13428-019-01231-3)
by Selker, van den Bergh, Criss, and Wagenmakers. The gamma vector can
be constrained in other ways and the constrains can be soft (i.e.,
priors with small SDs) which means that a continuum of more and more
simplified models can be obtained.

The two new simple link functions are *experimental*, in particular
the default priors are not well calibrated, meaning that some tweaking
may be necessary to achieve efficient sampling.

In order to use the new link functions the appropriate name has to be
specified when calling the make_stan_data, make_stan_model, and
gamma_to_crit functions, as described in the documentation.

# bhsdtr

The bhsdtr (short for Bayesian Hierarchical Signal Detection Theory
with Ratings) package implements a novel method of Bayesian inference
for hierarchical or non-hierarchical equal variance normal Signal
Detection Theory models with one or more criteria. It uses the
state-of-the-art platform [Stan](http://mc-stan.org/) for sampling
from posterior distributions. Our method can accommodate binary
responses as well as additional ratings and an arbitrary number of
nested or crossed random grouping factors. SDT parameters can be
regressed on additional predictors within the same model via
intermediate unconstrained parameters, and the model can be extended
by using automatically generated human-readable Stan code as a
template.

## Background

The equal-variance SDT with one criterion is almost (d' is not constrained to be non-negative) quivalent to probit
regression (see [this
paper](http://www.columbia.edu/~ld208/psymeth98.pdf) by DeCarlo) which
means that any software capable of fitting hierarchical generalized
linear models can be used to fit the hierarchical version of
equal-variance SDT *with one criterion and possibly negative d'*. However, the single-criterion
SDT model is untestable, because the data and the model have the same
dimensionality. The main reason for using SDT is to deconfound
sensitivity and bias. This can only be achieved if an SDT model is
approximately true, but there is no way to test it in the
single-criterion case. An SDT model becomes testable (e.g., by
comparing the theoretical and the observed ROC curves) when it is
generalized - by introducing additional criteria - to the version that
accomodates ratings (e.g., "I am almost certain that this item is
new").

SDT is a *non-linear* model. An immediate consequence of non-linearity
is that inference based on data aggregated over ''random'' grouping
factors (such as subjects or items) is invalid, because the resulting
estimates are biased (see [this
paper](http://rouder.psyc.missouri.edu/sites/default/files/morey-jmp-zROC-2008_0.pdf)
by Morey, Pratte, and Rouder for a demonstration, or see our [preprint](http://dx.doi.org/10.23668/psycharchives.2725)
for an even more striking demonstration). The only way to
avoid this problem is to model the (possibly
correlated) effects of all the relevant random grouping factors.

A subset of hierarchical SDT models with ratings can be fitted using
hierarchical ordered regression models, such as the cumulative model
in the excellent
[brms](https://cran.r-project.org/web/packages/brms/index.html)
package. As we explain [in the preprint](http://dx.doi.org/10.23668/psycharchives.2725), the d' parameter is
non-negative by definition and ignoring this assumption may lead to
problems if a bayesian SDT model is used. Without some modifications (which can be done in the brms package) hierarchical ordinal regression models do not restrict
the d' to be non-negative, becaues in such models d' is just the
unconstrained linear regression slope that represents the effect of
the stimulus class (''noise'' or ''signal''). Moreover, in typical
situations it does not make much sense to assume that the d' random
effects normally distributed. Finally, in the cumulative
model the parameters that correspond to the criteria in an SDT model
cannot be affected differently by the same grouping factor (i.e., the
effects are constant across categories), because the criteria in this
model are simply additional effects in the linear part. This means
that the model assumes that the pattern of the criteria is the same
for every participant (or item, etc.). Participants differ in their criteria
placement patterns and so the data from a typical rating experiment
cannot be independent given an SDT model with ratings represented as
a cumulative model.

In the bhsdtr package the generalized SDT model is supplemented with a
hierarchical linear regression structure (normally distributed
correlated random effects) thanks to a novel parametrization described
in [this preprint](http://dx.doi.org/10.23668/psycharchives.2725)
(which is now under review), and (more concisely) in the package
documentation. [Here](https://github.com/boryspaulewicz/bhsdtr/tree/master/inst/preprint/analysis_script.R)
is the annotated R script that performs all the analyses and produces
all the tables and some of the figures in the paper.

## Features

The bhsdtr package can be used to:

- fit generalized (more than one criterion), [meta-d'](http://www.columbia.edu/~bsm2105/type2sdt/), or basic (one criterion) equal variance SDT models
- fit hierarchical or non-hierarchical (e.g., single participant) models
- assess the fit using publication-ready ROC curve and combined response distribution plots with predictive intervals calculated for the chosen alpha level
- model the dependence of the SDT parameters on additional variables (e.g., task difficulty) using separate linear regression structures for the delta (d', meta-d') and gamma (criteria) parameters

### Prerequisites

A fairly up-to-date version of [R](https://www.r-project.org/) with
[the devtools
package](https://cran.r-project.org/web/packages/devtools/index.html)
already installed.

## Installing

The bhsdtr package, together will all of its dependencies, can be
installed directly from this github repository using the devtools
package:

```
devtools::install_git('git://github.com/boryspaulewicz/bhsdtr')
```

The installed package can be loaded using:

```
library(bhsdtr)
```

## Usage example

The package contains the gabor dataset


```
data(gabor)
head(gabor)
?gabor
```

To fit a hierarchical SDT model to this data we need to create some
data structures required by the stan function. This is how you can
create the combined response variable that encodes both the binary
classification decision and rating:

```
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
```

The combined responses have to be aggregated to make the sampling more
efficient. This is done using the aggregate_resoponses function, which
requires the names of the stimulus and the combined response
variables. If the data have a hierarchical structure, this structure
has to be preserved by listing the variables that cannot be collapsed
by the aggregation step. Here we list three variables that have to be
preserved: duration, id and order.

```
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
```

Finally, the fixed and random effects structure has to be specified
using lists of R model formulae. Here we assume that d' (= exp(delta))
depends on duration (a within-subject variable) and order (a
between-subject variable), but gamma (from which the criteria
parameter vector is derived) depends only on order. There is only one
random grouping factor - id - which represents the subjects. Note that
the random effects specification is a *list of lists* of model
formulate. That's because there can be more than one random grouping
factor.

```
fixed = list(delta = ~ -1 + duration:order, gamma = ~ -1 + order)
random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
```

Now we can start sampling:

```
fit = stan(model_code = make_stan_model(random),
    data = make_stan_data(adata, fixed, random),
    pars = c('delta_fixed', 'gamma_fixed',
        'delta_sd_1', 'gamma_sd_1',
        'delta_random_1', 'gamma_random_1',
        'Corr_delta_1', 'Corr_gamma_1',
        ## we need counts_new for plotting
        'counts_new'),
    iter = 8000,
    chains = 4)
```

When the make_stan_model and make_stan_data functions are called with
the optional metad=TRUE argument the meta-d' model is fitted. There
are two delta (d') parameters in the meta-d' model and so the
delta_fixed regression coefficients form a two-row matrix: the first
row represents the fixed effects for the d' parameter and the
second row represents the fixed effects for the meta-d' parameter.

The model fit can be assessed using the plot_sdt_fit function, which
produces ROC curve plots ...


```
plot_sdt_fit(fit, adata, c('order', 'duration')))
```

![ROC curve](inst/preprint/roc_fit.png)

... or combined response distribution plots:

```
plot_sdt_fit(fit, adata, c('order', 'duration'), type = 'response')
```

![Combined response distributions](inst/preprint/response_fit.png)
