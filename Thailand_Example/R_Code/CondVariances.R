################################ ///// ################################
#####
#####                    Population Estimation
#####
#####                  RESULTS SUMMARY FUNCTIONS
#####
#####               CALCULATE CONDITIONAL VARIANCES
###
### LICENCE:            Released under the GNU General Public License v3. See the
###                     accompanying file 'LICENSE.txt' or visit
###                     <http://www.gnu.org/licenses/>. This code is distributed
###                     in the hope that it will be useful, but WITHOUT ANY
###                     WARRANTY; without even the implied warranty of
###                     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
###                     See the GNU General Public License for more details.
###
### DISCLAIMER:         The views and opinions expressed in this work
###                     are those of the authors and do not
###                     necessarily represent those of the United
###                     Nations. This work has not been formally
###                     edited and cleared by the United Nations.
#####
################################ ///// ################################

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
    require(coda)

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
