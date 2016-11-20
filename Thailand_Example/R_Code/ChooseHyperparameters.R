################################################################################
###
###  DESC:          Choose hyperparameters
###
###  DATE CREATED:  26th January 2012
###
###  AUTHOR:        Mark C Wheldon
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

hyper.params.f <-
    function(absDev.f, absDev.s, absDev.g, absDev.n, al.f = 0.5, al.s = 0.5
             ,al.g = 0.5, al.n = 0.5)
{
    list(alpha.fert.rate = al.f
         ,beta.fert.rate = absDev.beta(absDev.f, prob = 0.9, alpha = 0.5)
         ,alpha.surv.prop = al.s
         ,beta.surv.prop = absDev.beta(absDev.s, prob = 0.9, alpha = 0.5)
         ,alpha.mig.prop = al.g
         ,beta.mig.prop = absDev.beta(absDev.g, prob = 0.9, alpha = 0.5)
         ,alpha.population.count = al.n
         ,beta.population.count = absDev.beta(absDev.n, prob = 0.9, alpha = 0.5)
         )
}
