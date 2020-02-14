

It is now possible to fit hierarchical ordinal polytomous regression
models in bhsdtr by using order-preserving link functions for the
thresholds. This way not only the constant shift in the thresholds,
but also the differences in the pattern of the thresholds can be
accounted for in certain situations. To my knowledge this is a new
kind of hierarchical ordinal model that can be used to estimate
ordinal polytomous item-parameter variability in a way that until now
was not possible. If you have a bunch of Likert-type items, or
confidence ratings, or anything of this sort in your data, than you
may find this model useful.

Because by definition a polytomous item has a multinomial distribution
and the multinomial distribution does not have a separate parameter
for the residual variance, when item-parameter variability is not
accounted for by the model it can only be misplaced
(overdispersion). The correct way to account for item-parameter
variability in ordinal polytomous items is to allow for individual
threshold effects while respecting the ordering of the
thresholds. This can only be achieved by using the order-preserving
link functions.

Assume that the ordered responses ranging from 1 to K come from the
same general kind of process as the responses in an SDT model when the
stimulus is 1 ("noise"). For example, when John is faced with the
Likert-type item "How often did you cry during the last two weeks" an
unknown process produces some "internal experience" value. This latent
value manifests itself in the response "never" (= 1) if it is low
enough, "a few times" (= 2) if it is a bit higher than some latent
threshold, or in the response "frequently" (= 3) if it is higher than
some other, still higher, threshold. The thresholds are assumed to be
ordered, because otherwise the item would not have the intended
meaning and thus would not be a valid instrument of measuring
depression.

Now consider a participant who repeatedly responds to the same item
and who has the memory of this event wiped out after each
response. Despite the wiping out of memory we can still expect some
variability in the responses due to the instability of the complex
process which produces the internal value, as well as due to the
instability of the decision criteria. This variability is certainly
there, the question is if it is non-negligible, and if so, what can be
done about it. Unfortunately, without additional assumptions it is
impossible to model both kinds of variability at the same time even in
the hypothetical repeated measurement with memory wipeout scenario,
because every shift in the internal value is indistiguishible from a
uniform effect in the thresholds, i.e., a constant shift of all the
thresholds. As I will explain shortly, the required additional
assumptions that make such models identifiable under certain
conditions are necessarily *causal*.

In many situations this is both a statistical and a conceptual
problem. In particular, when the origin of the internal value is
unknown because the design is purely observational (e.g.,
questionnaire measurement), the way of describing the process of
responding to a Likert-type item in terms of internal values and
thresholds is simply redundant. As far as the design of the study is
concerned, the value of the internal evidence "comes from nowhere" -
all we can say about it is that it is relative to the thresholds - and
the position of each threshold is relative to the positions of other
thresholds and to the position of the internal value. In this case to
say that John thinks that he was not crying at all during the last two
weeks and to say that *John's* threshold for the response "a few
times" is above the amount of crying that John thinks he experienced,
or decided to reveal, is to say the same thing.

It follows that *in such situations we do not need the internal value
parameter at all*, unless there is some factor about which we can
perhaps assume that it can *only* produce a uniform effect in the
thresholds, such as the stimulus class ("noise" or "signal") in a
binary classification task, which can *perhaps* be assumed to induce a
change (*d'*) in the internal value only; In such cases the internal
value parameter may be justified and useful, because it does constrain
the model: in the repeated measurement case the model is no longer
saturated.

But what about the promise of Item Response Theory? Isn't it true that
we can separate the internal values (so called "person parameters")
from the thresholds (so called "item parameters") and place them on a
common scale? Once the participants leave the room the items retain
their estimated "difficulty levels", right? The item difficulty
(=thresholds) is not defined relative to some sample of participants
who filled the questionnaire. How is this possible if the two kinds of
effects are completely confounded? In particular, how does Item
Response Theory achieve when each person provides only a single
response to each item?

It is not an *achievment* of Item Response Theory, but it's
(unrealistic) *assumption*. The illusory nature of this achievement
becomes obvious as soon as we realize that without additional
(causal!)  assumptions it is impossible to distinguish between genuine
internal value effects and (dis)simulation (i.e., faking), which
conteptually corresponds to the change in the *participant*-specific,
not item-specific, decision thresholds, also known as
"item-parameters". There is nothing in the *desing* of the process of
measurement by a questionnaire that makes it possible to distinguish
genuine differences in internal values from the differences in the way
of responding. This "achievement" of Item Response Theory is a
consequence of treating the "item-parameters" as constant accross
participants when fitting an IRT model.

Because of random assignment of the stimulus class above chance
performance in a binary classification task cannot be faked, and so
the *d'* paramater definately captures at least part of the effect of
the stimulus on the internal value. On the other hand, every possible
pattern of responses to a questionnaire can be faked. The only way to
distinguish the "item" and "person" effects is to *assume selective
influence*, i.e., that there is some factor f (such as the stimulus
class in a binary classification task, or, in case of a depression
questionnaire, participants trait, or a therapeutic intervention)
which can *only* affect the internal value, thus producing a
*constant* shift in the thresholds. The assumption of selective
influence makes the model with the two kinds of parameters
identifiable if there are repeated measurements for at least one
participant and one item or a large number of participants provides
single responses to all of the items. In particular, in typical
applications of IRT it is simply *assumed* in the model that the
participants differ in the internal values, and the item parameters
are constant accross the sample of participants (selectivity).

In every book on Item Response Theory that I know of the issue of
measurement invariance due to the variability in the "item-parameters"
is reduced to the problem of Differential Item Functioning
(DIF). However, the variability in "item-parameters" that I talk about
here has nothing to do with DIF, because DIF is a *population* level
effect; An estimate of DIF is obtained when the estimates of
"item-parameters" obtained for two different (large) samples of
participants are compared. However, these estimates are obtained by
fixing the "item-parameters" within each sample. It does not take much
thought to see that the assumption of constant "item-parameters" at
the group/population/sample level is utterly unrealistic. It follows
that the achievement of separating participant and item effects in IRT
models is not an achievement at all, it is just wishful thinking
combined with aggressive marketing.

Some of the problems related to item-parameter variability can be
addressed *when additional causal assumptions of selective influence
are introduced and the design of the study makes the model
identifiable*. For example, we can use the data from a binary
classification task and model only the responses for one kind of
stimulus, e.g., "noise".

Here is a simple example which illustrates how we can fit a
hierarchical ordinal polytomous regression model using the bhsdtr
package. We will be using the data from one condition in the gabor
dataset:


```r
library(ggplot2)
library(bhsdtr)
library(Hmisc)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## r = 1 ~ absolutely clear image left, r = 4 ~ no experience left, r = 5 ~ no experience right, ... r = 8 ~  absolutely clear image right
gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
d1 = gabor[gabor$duration == '32 ms' & gabor$order == 'DECISION-RATING',]
## Only the left-tilted gabor patches
d1 = d1[d1$stim == 0,]
```
        
To fit a hierachical model we have to aggregate the data (this makes
the process of sampling much faster), create the fixed- and
random-effects lists, and create the data structures and the model
code required by stan. Note that we have introduced a new parameter -
eta - which represents the internal value. As we already know, without
the selectivity assumption the difference between item and person
parameters is arbitrary and so in bhsdtr the first element of the
fixed effects vector for the eta parameter is automatically fixed at 0
eto make the model identifiable (this bahavior may change in the near
future):

KOD 2

To simplify the task we will use maximum likelihood optimization,
since the model is saturated anyway:

KOD 3

Unsuprisingly, the model fits perfectly:

WYKRES 1

Now we will add the estimated gamma/criteria random effects to the
fixed effects to obtain participant-specific gamma/criteria estimates:

KOD 4

We can now translate between the gamma and the criteria vectors:

KOD 5

This is how the estimated participant-specific thresholds look like:

What we see here is *a mixture of d' and threshold effects*. The
estimates are somewhat similar to what we would have obtained if we
fitted an SDT model, but they clearly do not have the intended
interpretation; They do not deconfound the two things that an SDT
model is supposed to deconfound, namely sensitivity and bias. Contrary
to what seems to be commonly assumed, if we included the stimulus
class we would not solve the problem - for that we have to also
introduce the strong causal assumption of selective influence, which,
like any other strong assumption, can easily be false.

The justification of the selective influence assumption comes from a
psychological theory of the task, not from the design of the
study. All that is guaranteed by randomly assigning the stimulus class
(i.e., by choosing the stimulus class at random on every trial) is
that any statistical effect of the stimulus class is an unbiased
estimate of the total causal effect of the stimulus class. In order to
empirically justify the assumption of *selective* influence of the
stimulus class on the internal evidence value we would have to *test*
if the stimulus class does not also affect the thresholds.

We can do exactly this. We will use the 'log_distance' link function,
because when this link function is used one of the elements of the
gamma vector directly represents the position of the main decision
criterion, and all the other elements represent the pattern of the
distances between the thresholds. If the stimulus does not affect the
pattern of the thresholds we should observe an effect only for the
main criterion parameter:

KOD

There at least two important advantages that an SDT model fitted to
binary classification task data has over an IRT model fitted to
questionnaire data. One is that the causal assumption of selective
influence of the stimulus class on the internal values is much more
plausible than the - implicit in most applications of IRT - causal
assumption of selective influence of participants on the
"person-parameters". The other is that in a typical binary
classification task case the process of responding is repeated many
times, the task is relatively simple ("was the gabor patch tilted left
or right" vs "do you feel that other people are happier than you") and
the processes responsible for the internal values and for the
placement of the decision thresholds are likely to be much more
stable.

If we campared the estimates of the main decision criterion (i.e.,
bias) in SDT models fitted to two different large samples of
participants the difference would likely be negligible with very
narrow confidence/credible intervals. By the same token, instability
of item parameters in IRT models cannot be estimated by assuming
selective influence and comparing estimates of "item-parameters"
obtained in two different large samples, because large-sample average
estimates of item parameters are the best way to *hide or misplace*
(it becomes part of the variability in the "person-parameters") the
intra- or inter-individual variability in item parameters; By assuming
that internal values are participant-specific and item parameters are
item-specific we force the estimates of item parameters to become
estimates of average person-specific thresholds. This not only begs
the question, but also, because the model is non-linear, makes the
estimates of item-parameters *asymptotically biased*.

Back to our simulation. The variability due to *d'* is still here, we
just cannot deconfound it from response bias. This means that, just
like in an IRT model that does not fully account for true person- and
and item-parameter effects, there is *no way to correctly estimate
reliability of the estimates of person-parameters/internal values*. We
will now demonstrate this fact.

First, we will fit a hierarchical SDT model with multiple criteria and
an order-preserving link function. This will give us realistic
estimates of "person-" (*d'*) and "item"-parameters (the thresholds).

KOD 4

What will happen if we emulate a typical IRT analysis and fit a
hierarchical ordinal model in which the thresholds are assumed to be
constant accross participants? After all, we can just invoke the
magical power of wishful thinking and assume that the thresholds are
item-parameters, right? How well will the internal value estimates
based on multiple responses to a simple, non-personal question match
the estimates based on a model that accounts for participant specific
threshold patterns?

KOD 5

The answer is - the two kinds of estimates will match poorly. Here is
the comparison of point estimates:

and here is the comparison of interval estimates:

WYKRES 3

Let's see more directly how the lengths of the credible intervals
differ. If the interval estimate widths were comparable they should
all be scattered closely to the horizontal line at 1. Unsurprisingly,
this is now what we see here:

WYKRES

The credible intervals directly represent reliability. For the
interval estimates based on a simplified model to contain the great
majority of internal values *or* account for the variability in the
thresholds they would have to be wider. We do not need to look at
"item-parameter" estimates in the simplified model, because it is
already obvious that the point and interval estimates of the
thresholds based on this model have to be severily biased.

An ordinal model for ordered polytomous responses that does not
correctly account for the variability in the thresholds produces
biased point and interval estimates in *every* parameter, and the
estimates seem more reliable than they truly are. That's because it
assumes that the thresholds, which are determined both by the items
and by the properties of the participants, are constant accross the
participants. This assumption, which is universally accepted in Item
Response Theory, is so ridiculous it should be laughet at on the
street.
