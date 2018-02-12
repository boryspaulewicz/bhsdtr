# bhsdtr

The bhsdtr R package implements a novel method of Bayesian inference
for hierarchical or non-hierarchical equal variance normal Signal
Detection Theory models with one or more criteria. It uses the
state-of-the-art platform Stan for sampling from posterior
distributions. Our method can accommodate binary responses as well as
additional ratings and an arbitrary number of nested or crossed random
grouping factors. The SDT parameters can be regressed on additional
predictors within the same model via intermediate unconstrained
parameters, and the model can be extended by using automatically
generated human-readable Stan code as a template.

The equal-variance SDT with one criterion is trivially equivalent to
probit regression (see
[this](http://www.columbia.edu/~ld208/psymeth98.pdf) paper by DeCarlo)
which means that any software capable of fitting hierarhical
generalized linear models can be used to fit the hierarchical version
of equal-variance SDT *with one criterion*. However, the
single-criterion SDT model is untestable, because the data and the
model have the same dimensionality (=2). The SDT model becomes
testable (e.g., by comparing the theoretical and the observed ROC
curves) when it is generalized - by introducing additional criteria -
to the version that accomodates ratings.

In the bhsdtr package the generalized SDT model is supplemented with a
general hierarchical linear regression structure thanks to a novel
parametrization which is described in this non peer-reviewed
[preliminary
paper](https://github.com/boryspaulewicz/bhsdtr/tree/master/inst/preprint/paper.pdf). This
parametrization requires some getting used to but it is the necessary
price to pay for the correctness of the implementation of the general
hierarchical linear regression structure.

### Prerequisites

All you need is a fairly up-to-date version of
[R](https://www.r-project.org/).

## Installing

The bhsdtr package, together will all its dependencies, can be
installed directly from this github repository using the devtools
package:

```
devtools::install_git('boryspaulewicz/bhsdtr')
```

The installed package can be loaded using:

```
library(bhsdtr)
```

## Usage example

The bhsdtr package contains the gabor dataset


```
data(gabor)
head(gabor)
```

We need to create some data structures required by the stan
function. This is how you can create the combined response variable
that encodes both the binary classification decision and rating:

```
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
```

The combined responses have to be aggregated to make the sampling more
efficient. This is done using the aggregate_resoponses function, which
requires the names of the stimulus variable and the combined response
variable. If the data have a hierarchical structure, this structure
has to be preserved by listing the variables that cannot be collapsed
by the aggregation step. Here we list three variables that have to be
preserved: duration, id and orded.

```
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
```

Finally, the fixed and random effects structure has to be specified
using lists of R model formulae. Here we assume that d' (= exp(delta))
depends on duration (a within-subject variable) and order (a
between-subject variable), but gamma (from which the criteria
parameter vector is derived) depends only on order. There is only one
random grouping factor - id - which represents the subjects.

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
