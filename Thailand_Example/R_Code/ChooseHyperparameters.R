################################################################################
###
###  DESC:          Choose hyperparameters
###
###  DATE CREATED:  26th January 2012
###
###  AUTHOR:        Mark C Wheldon
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
###
################################################################################

### ############################################################################
### * SYNOPSIS

## Choose values for alpha and beta. See also 'WheldonReport20111103.Rnw'.


### ############################################################################
### * FUNCTIONS

################################################################################
### ** Inverse Gamma Functions

rinvGamma <- function(n, shape, scale)
  {
    return(1/rgamma(n, shape = shape, rate = scale))
  }

dinvGamma <- function(x, shape, scale, log = FALSE)
  {
      if(log) {
          return(log(scale^shape / gamma(shape) * (1/x)^(shape + 1) *
           exp(-scale/x)))
      } else {
          return(scale^shape / gamma(shape) * (1/x)^(shape + 1) *
           exp(-scale/x))
      }
  }

qinvGamma <- function(p, shape, scale, lower.tail = TRUE)
  {
      return(1/qgamma(1-p, shape = shape, rate = scale
                      ,lower.tail = lower.tail)
             )
  }

pinvGamma <- function(q, shape, scale, lower.tail = TRUE)
{
    pgamma(1/q, shape = shape, rate = scale, lower.tail = !lower.tail)
}


################################################################################
### ** t-distribution Functions

mae.t.ab.f <- function(a,b)
{
    sqrt(2 * b / pi) * gamma(a-0.5) / gamma(a)
}


################################################################################
### ** Hyperparameter functions

################################################################################
### *** Set by MAD / quantiles of the marginal

absDev.beta <- function(absDev, prob = 0.5, alpha = 0.5)
    {
        ##
        ## Calculates beta such that Pr( |X-mu| <= absDev ) = prob
        ## for X a location-scale t distribution with 2*alpha df
        ##

        return(alpha * (absDev / qt(p = (prob + 1)/2, df = 2 * alpha))^2)
    }

make.hyper.params <-
    function(absDev = list(fert = 0.1, surv = 0.1, mig = 0.2
             ,pop = 0.1)
             ,prob = list(fert = 0.9, surv = 0.9, mig = 0.9, pop = 0.9)
             ,alpha = list(fert = 0.5, surv = 0.5, mig = 0.5, pop = 0.5)
             ,s.star
             )
{
    ##
    ## ARGUMENTS:
    ##
    ## absDev           List of the 'p's in the elicitation statement: "With
    ##                  probability 'prob', the true values are within 'p'
    ##                  percent of the initial estimates".
    ##
    ## prob             List giving probabilities with which the elicitation
    ##                  statement holds.
    ##
    ## alpha            List of inverse gamma alpha parameters.
    ##
    ## s.star           Initial estimates of survival proportions.
    ##

    ## -------* Functions

    logit <- function(p) log(p / (1-p))
    inv.logit <- function(x) exp(x) / (1 + exp(x))


    ## -------* Fert, Mig, Pop

    beta.fert.rate <-
        absDev.beta(absDev = log(1 - absDev$fert)
                          ,prob = prob$fert, alpha = alpha$fert
                          )

    beta.mig.prop <-
        absDev.beta(absDev = absDev$mig
                          ,prob = prob$mig, alpha = alpha$mig
                          )

    beta.population.count <-
        absDev.beta(absDev = log(1 - absDev$pop)
                          ,prob = prob$pop, alpha = alpha$pop
                          )

    beta.srb <-
        absDev.beta(absDev = log(1 - absDev$srb)
                          ,prob = prob$srb, alpha = alpha$srb
                          )


    ## -------* Surv

    ## The absDev is changed so that the min absDev on the (1-surv) scale is at
    ## least the input absDev.

    ## -------** Lower and upper quantiles of logit s

    logit.death.LQ <- function(beta, logit.death.prob)
    {
        logit.death.prob - qt(p = (1 + prob$surv) / 2, df = 2 * alpha$surv) *
            sqrt(beta / alpha$surv)
    }

    logit.death.UQ <- function(beta, logit.death.prob)
    {
        logit.death.prob + qt(p = (1 + prob$surv) / 2, df = 2 * alpha$surv) *
            sqrt(beta / alpha$surv)
    }


    ## -------** Target function for surv rel error

    surv.target <- function(aD, max.death.prop)
    {

        be <-
            absDev.beta(absDev = aD, prob = prob$surv, alpha = alpha$surv)

        L <- logit.death.LQ(be, logit(max.death.prop))
        U <- logit.death.UQ(be, logit(max.death.prop))

        if(U <= L) return(Inf)

        if(max.death.prop <= 0.5)
            x <- (max.death.prop - inv.logit(L)) / max.death.prop
        else
            x <- (inv.logit(U) - max.death.prop) / max.death.prop

        return((x - absDev$surv)^2)
    }


    ## -------** Calculate beta surv

    absDev.death.prop <-
        optimize(surv.target, interval = c(log(1 + absDev$surv), log(2))
                 ,max.death.prop = 1 - min(s.star)
                 ,tol = 1E-7
                 )$minimum

    beta.surv.prop <-
        absDev.beta(absDev = absDev.death.prop
                          ,prob = prob$surv, alpha = alpha$surv
                          )


    ## -------* Output

    list(al.f = alpha$fert
         ,be.f = beta.fert.rate
         ,al.s = alpha$surv
         ,be.s = beta.surv.prop
         ,al.g = alpha$mig
         ,be.g =  beta.mig.prop
         ,al.n = alpha$pop
        ,be.n = beta.population.count
        ,al.srb = alpha$srb
         ,be.srb = beta.srb
         )

}
