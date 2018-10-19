################################################################################
###
###  TITLE:             proposal_variances.R
###
###  AUTHOR:            Mark Wheldon
###
###  DATE CREATED:      2018-10-19
###
###  DESCRIPTION:       Proposal variances for Bayesian reconstruction.
###
###  DISCLAIMER:        The views and opinions expressed in this work
###                     are those of the authors and do not
###                     necessarily represent those of the United
###                     Nations. This work has not been formally
###                     edited and cleared by the United Nations.
###
###-----------------------------------------------------------------------------
###
###  SYNOPSIS
###
###  Tuning diagnostics: conditional variances to choose Metropolis proposal
###  variances.
###
###  MCMC convergence diagnostic.
###
###-----------------------------------------------------------------------------
###
################################################################################


##' Plot acceptance proportions for a Bayesian popualtion reconstruction
##'
##' Metropolis acceptance proportions should be between about 0.2 and
##' 0.5 for 'good' mixing. This function plots them for a completed
##' run of \code{\link{pop.est.sampler}}. They can be accessed
##' directly via
##' \code{\var{results.recon}$alg.stats$acceptance.proportions}.
##'
##' @param results.recon Output from \code{\link{pop.est.sampler}}.
##' @param ylim Limits for y axes for plots.
##' @return Nothing. Creates a plot of acceptance proportions.
##' @author Mark C. Wheldon
##' @export plot.acceptance.props
plot.acceptance.props <- function(results.recon, ylim = c(0, 1)) {

    ## Algorithm stats
    results.recon$alg.stats$acceptance.proportions


    ## Plot acceptance ratios
    par(mfrow = c(3,2))
    for(i in 1:3) {
        z <- results.recon$alg.stats$acceptance.proportions[[i]]
        zn <- names(results.recon$alg.stats$acceptance.proportions)[i]

        matplot(x=colnames(z), y = t(z), type = "b", pch = 1, main = zn
               ,ylim = ylim, xlab = "year", ylab = "acc. prop.")

        rect(xleft = as.numeric(min(colnames(z))) - 5
            ,ybottom = 0.2
            ,xright = as.numeric(max(colnames(z))) + 5
            ,ytop = 0.5
            ,col = "grey"
            ,density = 10
             )
    }
    z <- results.recon$alg.stats$acceptance.proportions[["baseline.count"]]
    rownames(z) <- sapply(strsplit(rownames(z), "[^0-9]"), "[[", 1)
    plot(x =rownames(z)
        ,y = z
        ,type = "b", main = "baseline.count"
        ,ylim = ylim, xlab = "age", ylab = "acc. prop."
         )
    rect(xleft = as.numeric(min(rownames(z))) - 5
        ,ybottom = 0.2
        ,xright = as.numeric(max(rownames(z))) + 5
        ,ytop = 0.5
        ,col = "grey"
        ,density = 10
         )
    z <- results.recon$alg.stats$acceptance.proportions[["srb"]]
    plot(x =colnames(z)
        ,y = z
        ,type = "b", main = "srb"
        ,ylim = ylim, xlab = "age", ylab = "acc. prop."
         )
    rect(xleft = as.numeric(min(colnames(z))) - 5
        ,ybottom = 0.2
        ,xright = as.numeric(max(colnames(z))) + 5
        ,ytop = 0.5
        ,col = "grey"
        ,density = 10
         )
    with(results.recon$alg.stats$acceptance.proportions, {
        barplot(height = c(fert = sigmasq.f, surv = sigmasq.s, mig = sigmasq.g
                         , pop = sigmasq.n, srb = sigmasq.srb)
               ,ylim = c(0,1)
               ,main = "variances"
                )
        abline(h = c(0.1, 0.5, 0.9), lty = 3)
    }
    )

}

##' Calculate conditional variances from MCMC chains of a population reconstruction
##'
##' To aid tuning of the Metropolis-Hastings algorithm, conditional
##' variances can be computed and inspected.
##'
##' @param results.recon Output from \code{\link{pop.est.sampler}}.
##' @param plot Logical: plot variances?
##' @param return.res Logical: return calculated conditional variances
##'     and suggested proposal variances (\cite{Roberts & Rosenthal, 2001}) in a list?
##' @return If \code{isTRUE{return.res}}, a list with conditional
##'     variances for each parameter, otherwise nothing. A plot is
##'     generated if \code{isTRUE{plot}}.
##' @author Mark C. Wheldon
##' @references Roberts, G. O., and Rosenthal, J. S. (2001),
##'     "Optimal Scaling for Various Metropolis-Hastings Algorithms",
##'     Statistical Science, 16, 351-367.
##' @export
conditional.variances <- function(results.recon, plot = TRUE, return.res = TRUE) {

    ## Check

    if(ncol(results.recon$surv.prop.mcmc[[1]]) > nrow(results.recon$surv.prop.mcmc[[1]])) {
        stop("Too few iterations. Save at least ", 2 * ncol(results.recon$surv.prop.mcmc[[1]]), ".")
        }

    ## Calculate conditional variances

    vitalCondVars <- list()

    vitalCondVars$fert.rate <-
        chain.cond.vars.feb08(log(results.recon$fert.rate.mcmc))

    vitalCondVars$surv.prop <-
        (chain.cond.vars.feb08(logit(results.recon$surv.prop.mcmc[["female"]])) +
         chain.cond.vars.feb08(logit(results.recon$surv.prop.mcmc[["male"]]))) / 2

    vitalCondVars$mig <-
        (chain.cond.vars.feb08(results.recon$mig.prop.mcmc[["female"]]) +
         chain.cond.vars.feb08(results.recon$mig.prop.mcmc[["male"]])) / 2

    vitalCondVars$population.count <-
        (chain.cond.vars.feb08(log(results.recon$baseline.count.mcmc[["female"]])) +
         chain.cond.vars.feb08(log(results.recon$baseline.count.mcmc[["male"]]))) / 2

    vitalCondVars$srb <-
        chain.cond.vars.feb08(log(results.recon$srb.mcmc))

    vitalCondVars$variances <-
        chain.cond.vars.feb08(log(results.recon$variances.mcmc))

    new.prop.vars <- lapply(vitalCondVars, "*", 2.3^2)


    if(plot) {

        par(mfrow = c(3,2))
        plot(c(0,0), xlim = c(0, max(results.recon$alg.params$prop.varcovar$fert.rate
                                    ,new.prop.vars$fert.rate)) * 1.01
            ,ylim = c(0, max(results.recon$alg.params$prop.varcovar$fert.rate
                            ,new.prop.vars$fert.rate)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "fert.rate")
        points(results.recon$alg.params$prop.varcovar$fert.rate, new.prop.vars$fert.rate)
        abline(a = 0, b = 1, col = "blue")
        plot(c(0,0), xlim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$surv.prop, "[", 1, 1)
                                    ,new.prop.vars$surv.prop)) * 1.01
            ,ylim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$surv.prop, "[", 1, 1)
                            ,new.prop.vars$surv.prop)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "surv.prop")
        points(sapply(results.recon$alg.params$prop.varcovar$surv.prop, "[", 1, 1), new.prop.vars$surv.prop)
        abline(a = 0, b = 1, col = "blue")
        plot(c(0,0), xlim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$mig, "[", 1, 1)
                                    ,new.prop.vars$mig)) * 1.01
            ,ylim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$mig, "[", 1, 1)
                            ,new.prop.vars$mig)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "mig")
        points(sapply(results.recon$alg.params$prop.varcovar$mig, "[", 1, 1), new.prop.vars$mig)
        abline(a = 0, b = 1, col = "blue")
        plot(c(0,0), xlim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
                                    ,new.prop.vars$baseline.pop.count)) * 1.01
            ,ylim = c(0, max(sapply(results.recon$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
                            ,new.prop.vars$baseline.pop.count)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "baseline.pop.count")
        points(sapply(results.recon$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
             , new.prop.vars$baseline.pop.count)
        abline(a = 0, b = 1, col = "blue")
        plot(c(0,0), xlim = c(0, max(unlist(results.recon$alg.params$prop.varcovar$srb)
                                    ,new.prop.vars$srb)) * 1.01
            ,ylim = c(0, max(unlist(results.recon$alg.params$prop.varcovar$srb)
                            ,new.prop.vars$srb)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "srb")
        points(unlist(results.recon$alg.params$prop.varcovar$srb), new.prop.vars$srb)
        abline(a = 0, b = 1, col = "blue")
        plot(c(0,0), xlim = c(0, max(unlist(results.recon$alg.params$prop.varcovar$variances)
                                    ,new.prop.vars$variances)) * 1.01
            ,ylim = c(0, max(unlist(results.recon$alg.params$prop.varcovar$variances)
                            ,new.prop.vars$variances)) * 1.01
            ,type = "n", xlab = "old", ylab = "new", main = "variances")
        points(unlist(results.recon$alg.params$prop.varcovar$variances), new.prop.vars$variances)
        abline(a = 0, b = 1, col = "blue")
    }

    ## RETURN
    if(return.res) {
        return(list(conditional.variances = vitalCondVars,
                    suggested.proposal.variances = new.prop.vars))
    } else {
        return(invisible())
        }

}



chain.cond.vars.feb08 <-
    function(chain, tryC = TRUE, ..., try.finally = NULL)

    #-- Calculate conditional variances in mcmc output --#
    #
    # DATE         : 8 November 2010
    # AUTHOR       : Mark Wheldon
    #
    # ARGUMENTS
    #   chain      : mcmc object
    #   tryC       : logical. Should lm be wrapped by tryCatch()?
    #   try.finally: finally function for tryCatch
    #   ...        : passed to tryCatch()
    #
    # DESCRIPTION
    #   'Chain' is the output from an MCMC sampler consisting of
    #   N iterations, K variables, x_1, ..., x_K.
    #
    #   This function calculates Var(x_k | x_{-k}) for all variables
    #   in 'chain' using linear regression of x_k on x_{-k}.
    #

{

    #.. Coerce 'chain' to a data frame
    x <- as.data.frame(chain)

    #.. List of formulae for all regressions with each
    #   variable as response, all others as predictors
    form.list <- list()
    for(i in 1:ncol(x)) {
      form.list[[i]] <- formula(x[1,c(i, (1:ncol(x))[-i])])
    }
    names(form.list) <- colnames(x)

    #.. Perform regressions, take residual *std deviations*
    resid.std.devs <-
      sapply(form.list, FUN = function(z)
             {
                 if(tryC) {
                     tryCatch(summary(lm(z, data = x, model = FALSE))$sigma
                              ,...
                              ,finally = try.finally
                              )
                 } else {
                     summary(lm(z, data = x, model = FALSE))$sigma
                 }
             }
             )

    #.. Convert to variances
    return(resid.std.devs^2)

  }
