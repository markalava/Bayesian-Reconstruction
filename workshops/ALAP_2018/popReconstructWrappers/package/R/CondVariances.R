################################ ///// ################################
#####
#####                    Population Estimation
#####
#####                  RESULTS SUMMARY FUNCTIONS
#####
#####               CALCULATE CONDITIONAL VARIANCES
#####
################################ ///// ################################

chain.cond.vars <-
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
