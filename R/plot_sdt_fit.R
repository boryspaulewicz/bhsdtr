## -*- coding: utf-8 -*-

#' Creates ROC or response distribution model fit plots with posterior
#' predictive intervals.
#'
#' @param fit a stanfit object.
#' @param adata an aggregated data object produced by the
#'   \code{aggregate_responses} function. It must be the same aggregated data
#'   object that was used to produce the stanfit object.
#' @param variables an optional vector of names of variables. These variables,
#'   if specified, define the subsets of data that will be represented by
#'   separate plots.
#' @param type a type of plot required. If type = 'roc' then ROC curve plots are
#'   created, if type != 'roc' then response distribution plots are created.
#' @param alpha an alpha level for the posterior predictive intervals, e.g., if
#'   alpha = .05 then 100\% - alpha = 95\% posterior predictive intervals will
#'   be calculated.
#' @return a plot object.
#' @export
plot_sdt_fit = function(fit, adata, variables = NULL, type = 'roc', alpha = .05){
    s = as.data.frame(fit)
    rm(fit)
    cnt_new = t(s[,grep('counts_new', names(s))])
    rm(s)
    ## Tworzymy zbiór danych rozwinięty pionowo
    df = cbind(adata$data, adata$stimulus)
    names(df)[ncol(df)] = 'stimulus'
    stim = ncol(df)
    K = ncol(adata$counts)
    df = cbind(plyr::ddply(df, names(df), function(x)data.frame(1:K)), as.vector(t(adata$counts)))
    df = cbind(df, rep(1:nrow(adata$data), each = ncol(adata$counts)))
    ## Indeksy ważnych kolumn
    resp = ncol(df) - 2
    cnt = ncol(df) - 1
    obs = ncol(df)
    ## Kolejność taka, jak w próbkach ze Stan-a
    df = df[order(df[[resp]], df[[obs]]),]
    ## Agregacja obserwacji i próbek dla dowolnego wykresu
    if(length(variables) > 0){
        f = as.factor(df[[variables[1]]])
        if(length(variables) > 1)
            for(v in variables[-1])
                f = f:as.factor(df[[v]])
    }else{
        f = as.factor(rep('', nrow(df)))
    }
    dfa = stats::aggregate(df[[cnt]] ~ df[[resp]] + df[[stim]] + f, FUN = sum)
    cnt_new_a = matrix(nrow = ncol(cnt_new), ncol = nrow(dfa))
    pb = utils::txtProgressBar(min = 1, max = ncol(cnt_new), style = 3)
    print('Aggregating posterior samples...')
    for(r in 1:ncol(cnt_new)){
        cnt_new_a[r,] = stats::aggregate(cnt_new[,r] ~ df[[resp]] + df[[stim]] + f, FUN = sum)[[4]]
        utils::setTxtProgressBar(pb, r)
    }
    close(pb)
    rm(cnt_new)
    names(dfa) = c('response', 'stimulus', 'f', 'count')
    dfa$stimulus = as.factor(dfa$stimulus)
    dfa$n = plyr::ddply(dfa, c('f', 'stimulus'), function(x)data.frame(n = rep(sum(x$count), nrow(x))))[[3]]
    dfa$i = 1:nrow(dfa)
    if(type == 'roc'){
        ## Wykres ROC: wyliczamy p(Hit) ~ p(FA)
        print('Calculating ROC curves...')
        dfroc = plyr::ddply(dfa, c('stimulus', 'f'),
                      function(x){
                          cumfr_new = apply(t(cnt_new_a[,x$i]), 2, function(v)rev(cumsum(rev(v / x$n))))
                          data.frame(response = x$response,
                                     cumfr = rev(cumsum(rev(x$count / x$n))),
                                     cumfr.fit = apply(cumfr_new, 1, mean),
                                     cumfr.lo = apply(cumfr_new, 1, function(x)stats::quantile(x, alpha / 2)),
                                     cumfr.hi = apply(cumfr_new, 1, function(x)stats::quantile(x, 1 - alpha / 2)))
                      }, .progress = 'text')
        rm(dfa)
        ## dfroc$in.pi = as.numeric((dfroc$cumfr >= dfroc$cumfr.lo) && (dfroc$cumfr <= dfroc$cumfr.hi))
        ## dfroc$in.pi[dfroc$in.pi == 0] = .5
        dfrocs = dfroc[dfroc$stimulus == '1',]
        dfrocs$stim2 = dfroc[dfroc$stimulus == '2',]
        rm(dfroc)
        ggplot(dfrocs, aes(cumfr, stim2$cumfr)) +
            geom_line(aes(x = cumfr.fit, y = stim2$cumfr.fit), lty = 2) +
            geom_errorbar(aes(ymin = stim2$cumfr.lo, ymax = stim2$cumfr.hi, x = cumfr.fit), width = 0.02) +
            geom_errorbarh(aes(xmin = cumfr.lo, xmax = cumfr.hi, y = stim2$cumfr.fit), height = 0.02) +
            geom_point() +
            labs(x = 'p(Hit)', y = 'p(FA)') +
            coord_fixed() +
            facet_wrap(~f)
    }else{
        ## Wykres rozkładów odpowiedzi
        dfa[,c('count.lo', 'count.hi', 'count.fit')] = cbind(apply(cnt_new_a, 2, function(x)stats::quantile(x, alpha / 2)),
                                                             apply(cnt_new_a, 2, function(x)stats::quantile(x, 1 - alpha / 2)),
                                                             apply(cnt_new_a, 2, mean))
        ggplot(dfa, aes(response, count / n, color = stimulus), group = stimulus) +
            geom_errorbar(aes(ymin = count.lo / n, ymax = count.hi / n, lty = stimulus), width = 0.2) +
            geom_line(aes(y = count.fit / n, lty = stimulus)) +
            geom_point(aes(pch = stimulus)) +
            labs(x = 'Response', y = 'Frequency', color = 'Stimulus', pch = 'Stimulus', lty = 'Stimulus') +
            facet_wrap(~f)
    }
}
## ok
