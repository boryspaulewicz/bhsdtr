# bhsdtr

The bhsdtr R package implements a novel method of Bayesian inference
for hierarchical or non-hierarchical equal variance normal Signal
Detection Theory models with one or more criteria. It the
state-of-the-art platform Stan for sampling from posterior
distributions. Our method can accommodate binary responses as well as
additional ratings and an arbitrary number of nested or crossed random
grouping factors. The SDT parameters can be regressed on additional
predictors within the same model via intermediate unconstrained
parameters, and the model can be extended by using automatically
generated human-readable Stan code as a template.

### Prerequisites

All you need is a fairly up-to-date version of R.

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
by the aggregation step:

```
adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
```