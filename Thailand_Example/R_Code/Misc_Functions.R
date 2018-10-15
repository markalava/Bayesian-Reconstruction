################################################################################
###
###  DESC :            Miscellaneous functions for Bayesian estimation model
###
###  DATE CREATED :    9 Jan 2012
###
###  AUTHOR :          Mark Wheldon
###
################################################################################

### ############################################################################
### * Misc


#--- Odds and inverse ---#

odds <- function(p) p/(1-p)
inv.odds <- function(y) y/(1+y)


#--- Logit and inverse ---#

logit <- function(p) log(p / (1 - p))
invlogit <- function(x) exp(x) / (1 + exp(x))


#--- Half o logit and inverse ---#

hlogit <- function(x) 0.5*log(x/(1-x))
inv.hlogit <- function(x) exp(2*x)/(1+exp(2*x))


#--- Timestamp a Directory for 'make' ---#

dir.tstamp <- function(x = getwd()) {
    y <- round(runif(1) * 1E6)
    path <- file.path(x, paste("tstamp", y, sep = "_"))
    file.create(path)
    file.remove(path)
}


#--- Add a suffix to end of object names ---#

add.suff <-
    function(..., list = character()
             ,suff, separator = "."
             ,value = "character"
             ,verbose = getOption("verbose")
             )
{
    if(!is.character(list) || is.recursive(list))
        stop("'list' must be a character vector")
    names <- as.character(substitute(list(...)))[-1L]
    list <- c(list, names)
    if(identical(value, "character")) {
        if(verbose) message("value is 'character'")
        return(paste(list, suff, sep = separator))
    } else if(identical(value, "objects")) {
       if(verbose) message("value is 'objects'")
    parEnv <- parent.frame()
    for(y in list) {
        if(y %in% ls(envir = parEnv)) {
        nm <- paste0(c(y, suff), collapse = separator)
        assign(x = nm, value = get(y), envir = parEnv)
        if(verbose) message("'", nm, "' placed in calling environment")
    }
    }
}
}


#--- matrix.along ---#

matrix.along <- function(x, m, ...)
{
    matrix(x, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m), ...)
}


#--- mat.2.v ---#

mat.2.v <- function(m, ...)
{
    ## UPDATED 2014-02-24
    ## ~ Made 'names' incorrectly.

    for(i in 1:length(dim(m))) {
        if(is.null(dimnames(m)[[i]]))
            dimnames(m)[[i]] <- 1:(dim(m)[i])
    }
    x <- as.vector(m)
    names.gr <- expand.grid(dimnames(m)[[1]], dimnames(m)[[2]])
    names(x) <- paste(names.gr[,2], names.gr[,1], sep = ".")
    return(x)
}


#--- proposal variances ---#

make.prop.varcovar <-
    function(var, corr, age = 0:100, year = 1960:2000)
{
    ## UPDATED 2014-02-24
    ## ~ Only appeared to use the first element of var.

    ag.yr <- expand.grid(age, year)
    ag.yr <- paste(ag.yr[,2], ag.yr[,1], sep = ".")

    pv <- replicate(length(ag.yr), matrix(0, nrow = 2, ncol = 2)
                    ,simplify = FALSE
                    )
    for(x in 1:length(var)) {
        pv[[x]] <- diag(var[x] * (1 - corr), nrow = 2) +
                           matrix(var[x] * corr, nrow = 2, ncol = 2)
    }
    names(pv) <- ag.yr
    return(pv)
}


### ############################################################################
### * Inverse gamma functions

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


#--- Choose inverse gamma distn based on quantiles ---#

IG.chooseByQ <- function(probs, val, init.vals = c(1,2))
    # probs: the quantiles to match against
    # val: values of quantiles to match
  {
    #.. Function palatable to optim; simply minimizes absolute
    #   difference of quantiles
    qIG.paramsfirst <- function(theta)
      {
        a <- theta[1]
        b <- theta[2]
        q <- qinvGamma(p = probs, shape = a, scale = b)
        return(max(abs(q - val)))
      }

    params <- optim(init.vals, qIG.paramsfirst
                    ,control = list(trace = 0))$par
    names(params) <- c("shape", "scale")

    #.. calculate true quantiles to check
    true.q <- qinvGamma(probs, params[1], params[2])

    return(list(params = params
                ,true.q = true.q))
  }


IG.chooseByQ.fixScale <-
    function(probs, vals, init.scale = 2)
    # probs: the quantiles to match against
    # vals: values of quantiles to match
  {
    #.. Function palatable to optimize; simply minimizes squared
    #   difference of quantiles
    qIG.scaleOptim <- function(scale, shape, probs, vals)
      {
        q <- qinvGamma(p = probs, shape = shape, scale = scale)
        return(max((q - vals)^2))
      }
    qIG.shapeOptim <- function(shape, scale, probs, vals)
      {
        q <- qinvGamma(p = probs, shape = shape, scale = scale)
        return(max((q - vals)^2))
      }

    new.shape <- optimize(qIG.shapeOptim, interval = c(1E-4, 100)
             ,probs = probs, scale = init.scale, vals = vals
             )$minimum
    new.scale <- optimize(qIG.scaleOptim, interval = c(0, 100)
                       ,probs = probs, shape = new.shape, vals = vals
                       )$minimum

    #.. calculate true quantiles to check
    true.q <- qinvGamma(probs, new.shape, new.scale)

    return(list(params = c(shape = new.shape, scale = new.scale)
                ,true.q = true.q))
  }

IG.chooseByQ.fixShape <-
    function(probs, vals, fixed.shape = 2
             ,optim.interval = c(1E-10, 1E10))
    # probs: the quantiles to match against
    # vals: values of quantiles to match
  {
    #.. Function palatable to optimize; simply minimizes squared
    #   difference of quantiles
    qIG.scaleOptim <- function(scale, shape, probs, vals)
      {
        q <- qinvGamma(p = probs, shape = shape, scale = scale)
        return(max((q - vals)^2))
      }

    IG.scale = optimise(qIG.scaleOptim, interval = optim.interval
             ,shape = fixed.shape, probs = probs, vals = vals)

    #.. calculate true quantiles to check
    true.q <- qinvGamma(probs, fixed.shape, IG.scale$minimum)

    #.. calculate true probabilities
    true.p <- pinvGamma(vals, shape = fixed.shape
                        ,scale = IG.scale$minimum
                        )

    return(list(params = c(shape = fixed.shape, scale = IG.scale$minimum)
                ,true.q = true.q, true.p = true.p
                , optimise.objective = IG.scale$objective))
  }


IG.chooseByQv2 <- function(probs, vals, init.vals)
    # probs: the quantiles to match against
    # vals: values of quantiles to match
{
    # probs and val should both be length 2

    if(length(vals) != 2 || length(probs) != 2) {
        stop("'probs' and 'vals' must be of length 2")
    }

    #.. Function palatable to optim; simply minimizes absolute
    #   difference of quantiles

    qIG.paramsfirst <- function(theta, vals)
    {
        a <- theta[1]
        b <- theta[2]
        q <- qinvGamma(p = probs, shape = a, scale = b)
        return(max(abs(q - vals)))
    }


    #.. Setting init.vals if missing: gamma(k, l) is sum of
    #   k indpt expo(l) variables, so approximate as 1/k the
    #   quantile of an expo(l) variable

    if(missing(init.vals)) {
        p1 <- probs[1]
        p2 <- probs[2]
        v1 <- vals[1]
        v2 <- vals[2]

        # gamma rate. Start with shape = 1.
        rate.Gamma <- -1*log(p1)/(1/v1)
        shape.Gamma <- rate.Gamma*(1/v2)/(-log(p2))
        rate.Gamma <- -shape.Gamma*log(p1)/(1/v1)
        shape.Gamma <- rate.Gamma*(1/v2)/(-log(p2))

        init.vals <- c(shape.Gamma, 1/rate.Gamma)

    }


    #.. Run optim

    params <- optim(init.vals, qIG.paramsfirst, vals = vals
                    ,control = list(trace = 0)
                    )$par
    names(params) <- c("shape", "scale")

    #.. calculate true quantiles to check
    true.q <- qinvGamma(probs, params[1], params[2])

    return(list(params = params
                ,true.q = true.q))

}


##
##--- Function to generate inverse gamma parameters

invGam.paramGen.jan28 <-
    function(#.. *Variance* parameters
             fert.var, surv.var, mig.var, baseline.var, census.var

             #.. Functions required
             ,InvGamma.fun = "IG.chooseByQ")
{

    #--- Match functions ---#

    IG.fun.matched <- match.fun(InvGamma.fun)


    #--- Variance hyperparameters ---#

    asFert.IGparams <- suppressWarnings(
        IG.fun.matched(probs = c(0.5, 0.975)
                     ,val = c(fert.var, 2*fert.var))$params)
    asSurv.IGparams <- suppressWarnings(
        IG.fun.matched(probs = c(0.5, 0.975)
                     ,val = c(surv.var, 2*surv.var))$params)
    asMig.IGparams <- suppressWarnings(
        IG.fun.matched(probs = c(0.5, 0.975)
                     ,val = c(mig.var, 2*mig.var)
                     ,init.vals = c(1, mig.var)
                     )$params)
    baseline.IGparams <- suppressWarnings(
        IG.fun.matched(probs = c(0.5, 0.975)
                     ,val = c(baseline.var, 2*baseline.var))$params)
    census.IGparams <- suppressWarnings(
        IG.fun.matched(probs = c(0.5, 0.975)
                     ,val = c(census.var, 2*census.var))$params)


    #--- Output ---#

    out <-
        list(al.f = asFert.IGparams[1], be.f = asFert.IGparams[2]
             ,al.s = asSurv.IGparams[1], be.s = asSurv.IGparams[2]
             ,al.g = asMig.IGparams[1], be.g = asMig.IGparams[2]
             ,al.b = baseline.IGparams[1], be.b = baseline.IGparams[2]
             ,al.lhood = census.IGparams[1]
             ,be.lhood = census.IGparams[2]
             )

    return(out)

}


#--- Choose normal variance based on asymmetric confidence interval ---#

norm.chooseByACI <-
    function(l, u, prob, sigma.interval = c(0.0001, 1000)) {

    if(l >= u) stop("'u' must be greater than 'l'")

    #.. function for optimize
    fun <- function(sigma, ll, uu, pr) {
        normlP <- pnorm(q = uu, mean = 0, sd = sigma) -
            pnorm(q = ll, mean = 0, sd = sigma)
        return(abs(normlP - pr))
    }
    op <- optimize(fun, interval = sigma.interval, ll = l, uu = u
             ,pr = prob)
    s <- op$minimum
    return(list(variance = s^2
                ,true.p = pnorm(u, sd = s) - pnorm(l, sd = s)
                ))
}


#--- Choose gamma distn based on quantiles ---#

gamma.choose.by.quant <- function(probs, vals, init.vals, ...)
    ## probs: the quantiles to match against
    ## vals: valsues of quantiles to match
  {

      ## Initial values
      if(missing(init.vals)) {
          init.shape <- abs(vals[2]-vals[1])/2 + 1
          init.vals <- c(init.shape, 1)
      }

    ## Function palatable to optim; simply minimizes absolute difference of
    ## quantiles
    qGamm.paramsfirst <- function(theta)
      {
        a <- theta[1]
        b <- theta[2]
        p <- pgamma(q = vals, shape = a, scale = b)
        return(sum((p - probs)^2))
      }

    params <- optim(par = init.vals, fn = qGamm.paramsfirst
                    #,control = list(trace = 0)
                    ,...
                    )$par
    names(params) <- c("shape", "scale")

    #.. calculate true quantiles to check
    true.p <- pgamma(vals, params[1], params[2])

    return(list(params = params
                ,true.p = true.p))
  }

gamma.choose.by.mode <- function(mode, upper, upper.p, init.vals = c(1, 1), ...)
    ## probs: the quantiles to match against
    ## vals: valsues of quantiles to match
  {

     if(init.vals[1] < 1) stop("initial value for mode < 1")

    ## Function palatable to optim; simply minimizes absolute difference of
    ## quantiles
    qGamm.paramsfirst <- function(theta)
      {
        a <- theta[1]
        b <- theta[2]
        if(a >= 1) optim.mode <- (a - 1) * b
        else optim.mode <- 0
        optim.upper <- qgamma(p = upper.p, shape = a, scale = b)
        return((optim.mode - mode)^2 + (optim.upper - upper)^2)
      }

    params <- optim(par = init.vals, fn = qGamm.paramsfirst
                    #,control = list(trace = 0)
                    ,...
                    )$par
    names(params) <- c("shape", "scale")

    #.. calculate true quantiles to check
    achieved.vals <-
        c(mode = as.numeric((params[1] - 1) * params[2])
          ,upper = qgamma(upper.p, shape = params[1], scale = params[2])
          )

    return(list(params = params
                ,achieved.vals = achieved.vals))
  }


################################################################################
### * Re-start previous chain

start.vals.nov11 <-
    function(prev.chain, mid.run = FALSE, init.values, assign.global = TRUE)
{
    ##
    ## ARGUMENTS:
    ##
    ##   prev.chain:
    ##   Character vector giving file path of '.RData' file of chain to
    ##   re-start.
    ##
    ##   mid.run:
    ##   Logical; is the saved chain a 'mid-run' save?
    ##
    ##   init.values:
    ##   Named list of initial value matrices or lists for vital rate
    ##   parameters.
    ##
    ##   assign.global:
    ##   Logical; assign start values to global environment or return as an
    ##   object?
    ##

    ## -------* Set up

    ## Save workspace
    si.path <- file.path(tempdir(), "start_vals_fun_save.RData")
    save.image(si.path)

    ## Load back on exit
    on.exit(load(si.path))


    ## -------* Load previous chain

    ## Load previous chain and catch object name
    prev.chain.name <- withVisible(load(prev.chain))$value

    prev.saved.iters <- nrow(get(prev.chain.name)$fert.rate.mcmc)

    if(mid.run || is.null(get(prev.chain.name)$alg.params)) {
        get(prev.chain.name)$alg.params$non.zero.fert.rows <-
            as.logical(apply(asFertINDI.mat == 0L, 1, function(z) !all(z)))
    }


    ## -------* Make start values

    ## Make a new matrix like this one
    matrix.along <- function(x, m, ...)
    {
        matrix(x, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m), ...)
    }

    ## Start values for VR parameters
    fert.start <- matrix.along(0, init.values$fert)
    fert.start[get(prev.chain.name)$alg.params$non.zero.fert.rows,] <-
        get(prev.chain.name)$fert.rate.mcmc[prev.saved.iters,]

    surv.start <-
        list(female = matrix.along(get(prev.chain.name)$surv.prop.mcmc$female[prev.saved.iters,]
             ,init.values$surv$female, byrow = FALSE)
             ,male = matrix.along(get(prev.chain.name)$surv.prop.mcmc$male[prev.saved.iters,]
              ,init.values$surv$male, byrow = FALSE)
             )

    mig.start <-
        list(female = matrix.along(get(prev.chain.name)$mig.prop.mcmc$female[prev.saved.iters,]
             ,init.values$mig$female, byrow = FALSE)
             ,male = matrix.along(get(prev.chain.name)$mig.prop.mcmc$male[prev.saved.iters,]
              ,init.values$mig$male, byrow = FALSE)
             )

    baseline.start <-
        list(female = matrix.along(get(prev.chain.name)$baseline.count.mcmc$female[prev.saved.iters,]
             ,init.values$baseline$female, byrow = FALSE)
             ,male = matrix.along(get(prev.chain.name)$baseline.count.mcmc$male[prev.saved.iters,]
              ,init.values$baseline$male, byrow = FALSE)
             )

    srb.start <-
        matrix.along(get(prev.chain.name)$srb.mcmc[prev.saved.iters,]
                     ,init.values$srb, byrow = FALSE)


    ## Make start values for variance parameters

    fert.var.start <-
        get(prev.chain.name)$variances.mcmc[,"fert.rate.var"][prev.saved.iters]

    surv.var.start <-
        get(prev.chain.name)$variances.mcmc[,"surv.prop.var"][prev.saved.iters]

    mig.var.start <-
        get(prev.chain.name)$variances.mcmc[,"mig.var"][prev.saved.iters]

    pop.var.start <-
        get(prev.chain.name)$variances.mcmc[,"population.count.var"][prev.saved.iters]

    srb.var.start <-
        get(prev.chain.name)$variances.mcmc[,"srb.var"][prev.saved.iters]


    ## -------** Assign start values to global environment

    if(assign.global) {
        for(z in c("fert.start", "surv.start", "mig.start", "baseline.start"
                   ,"fert.var.start", "surv.var.start", "mig.var.start"
                   ,"pop.var.start", "srb.var.start"
                   )) {
            assign(x = z, value = get(z), envir = .GlobalEnv)
            message(z, " assigned to .GlobalEnv")
        }
    } else {
        return(list(fert.start = fert.start, surv.start = surv.start
                       ,mig.start = mig.start, baseline.start = baseline.start
                       ,fert.var.start = fert.var.start
                       ,surv.var.start = surv.var.start
                       ,mig.var.start = mig.var.start
                       ,pop.var.start = pop.var.start
                       ,srb.var.start = srb.var.start
                       ))
    }
}


### ############################################################################
### * MCMC Diagnostics

#--- Conditional variances ---#

mcmc.cond.var <- function(mcmc.chain)
  {
    ## Caculates conditional variances from a sample of a set of
    #    variables via linear regression.

    # mcmc.chain:    An object coercible to a data frame with chains
    #                  as columns.

    require(coda)

    #.. Create data frame
    v.df <- as.data.frame(mcmc.chain)

    #.. List of formulae for all regressions with each
    #   variable as response, all others as predictors
    form.list <- list()
    for(i in 1:ncol(v.df)) {
      form.list[[i]] <- formula(v.df[1,c(i, (1:ncol(v.df))[-i])])
    }
    names(form.list) <- colnames(v.df)

    #.. Perform regressions, take residual *std deviations*
    resid.std.devs <-
      sapply(form.list, FUN = function(z)
             {
               summary(lm(z, data = v.df, model = FALSE))$sigma
             }
             )

    #.. Convert to variances
    return(resid.std.devs^2)

  }


### ############################################################################
### * Save Female, Male ff List

## 'proj.to.counts' and other functions output lists with ff pointers as
## elements. This function extracts the elements and saves them as a regular
## list in an '.RData' file.

save.MF.fflist <- function(x, f, tstamp.dir = NULL)
{
    ##
    ## Arguments:
    ##
    ## x                The list of ff objects
    ## f                The file name to save to
    ## tstamp.dir       timestamp directory with this function?
    ##

    ## name of list of ff objects
    nm <- deparse(substitute(x))

    ## Handle single-sex parameters (fert, srb)
    if(!inherits(x, "list")) {
        if(!inherits(x, "ff")) stop("'x' must be a list or inherit from 'ff'")
        assign(nm, x[])
        save(list = nm, file = f)
    } else {
        ## Two-sex parameters
        assign(nm, list(female = x$female[], male = x$male[]))
        save(list = nm, file = f)
    }

    ## time stamp
    if(!is.null(tstamp.dir)) {
        fs <- unlist(strsplit(f, split = .Platform$file.sep))
        if(length(fs) > 1) {
            d <- paste0(head(fs, -1), collapse = .Platform$file.sep)
            tstamp.dir(d)
        }
    }
}


### ############################################################################
### * Post-process sampler output

### ############################################################################
### ** Convert outputs to counts, etc.

## These are now obsolete. Use programs in the 'ResultsSummary' folder instead.

####
####
####-----------------------------------------------------------------
#### MAKE LESLIE MATRIX
####-----------------------------------------------------------------
####
leslie.oct27 <-
    function(pop, surv, fert, srb = 1.05, age.int = 5, label.dims = FALSE)

    #-- Make the leslie matrix for CCMPP --#
    #
    #   pop     :  population count at baseline
    #   fert    :  matrix of age specific fertility rates NOT yet
    #                mulitplied by age.int
    #   srb     :  sex ratio at birth matrix
    #   surv    :  Survivorship probabilities: the probability of
    #                reaching the age at the start of the interval.
    #              The first row should be nL0/(n*l0).
    #              The last row is survival for age.int years in the open
    #                interval
    #   proj.steps
    #   age.int :  needed for correct interpretation of survival
    #                and fertility rates
    #   label.dims
    #           :  should output have dimnames set?

{
    #
    #.. Constants
    #
    n.age.grps <- length(pop)
    n.surv <- length(surv)

    #
    #.. Make Leslie matrix
    #
    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)

    #.. first row = fert and birth survival
    k <- 1/(1+srb) * surv[1] * 0.5
    dbl.fert <- age.int*fert + c(age.int*fert[-1], 0) * surv[-1]
    lesM[1,] <- k * dbl.fert

    #.. rows 2:(n.age.grps) = survival ratios
    lesM[2:n.age.grps,1:(n.age.grps-1)] <- diag(surv[-c(1,n.surv)])
    lesM[n.age.grps,n.age.grps] <- surv[n.surv]

    if(label.dims) {
        age.labs <- seq(from = 0, by = 5, length = n.age.grps)
        dimnames(lesM) <- list(age.labs, age.labs)
    }

    #
    #.. return
    #
    return(lesM)
}

####
####
####-----------------------------------------------------------------
#### Calculate net number of migrants
####-----------------------------------------------------------------
####
netMig.27oct <- function(n1, n2, L)
{
    #-- Find net number of migrants in a CCMPP projection --#
    #
    # ARGUMENTS
    #
    #   n1      :  Population count vector at time t
    #   n2      :  Population count vector at time t + delta
    #   L       :  Leslie matrix used to get population at t + delta
    #
    #
    # METHOD
    #
    # Invert n2 = L(n1 + 0.5 mig) + (0.5)*mig
    # Can get proportions by pre-multiplying output by 'solve(diag(n1))'

    #
    #.. Make sure inputs are of correct form
    #
    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    L <- as.matrix(L)

    return(2 * solve(L + diag(nrow(L))) %*% (n2 - L %*% n1))
}


### ############################################################################
### ** Calculate life table columns from MCMC survival ratio output

## These are now obsolete. Use the lifeTable class instead.

nsx.2.nLx <- function(nsx, age.scale = 5, radix = 100000)
{
    #.. nsx is vector of survival ratios with first
    #   element ns0

    #.. checks
    #if(sum(nsx > 1) > 0 || sum(nsx < 0) > 0) stop("nsx must be a vector of survival ratios: i.e., 0 <= nsx <= 1")
    #if(age.scale < 0) stop("age.scale must be positive")
    #if(radix < 1) stop("radix must be >= 1")

    nLx <- vector("numeric", length = length(nsx))
    nLx[1] <- 5*radix*nsx[1]
    for(i in 2:length(nLx)) nLx[i] <- nsx[i] * nLx[i-1]
    return(nLx)
}

nLx.2.lx <- function(nLx, nax, age.scale = 5, radix = 100000
                      ,check.nax = TRUE)
{
    # here is the nax schedule from the UN General model LT for
    # females, leb = 46:
    # c(2.5, 2.5, 2.649, 2.578, 2.536, 2.528, 2.818, 2.533, 2.566, 2.587, 2.588, 2.573, 2.546, 2.498, 2.409, 2.312, 4.084)

    # need a vector of nax values same length as nLX
    if(check.nax && !identical(length(nLx), length(nax)))
        stop("nax and nLx must be of same length")

    lx <- vector("numeric", length = length(nLx))
    lx[1] <- radix
    for(i in 2:length(lx))
        lx[i] <- (nLx[i-1] - lx[i-1]*nax[i-1])/(age.scale - nax[i-1])
    return(lx)
}


####
####
####-----------------------------------------------------------------
#### CALCULATE nmx FROM SURV PROPS ASSUMING STATIONARY POP
####-----------------------------------------------------------------
####
nsx.to.nmx.nov01 <-
    function(nsx, l0 = 1000, n = 5, fix.names = TRUE)

    #-- Calculate nmx from surv props assuming stationary pop --#
    #
    #   nsx     :  survival proportions
    #   l0      :  radix to use for life table
    #   n       :  size of age interval
    #   fix.names
    #           :  sort out names

{
    #
    #.. Transform last nsx since it is T_{x+5} / T_x
    #
    k <- c(head(nsx, -1), tail(nsx,1) / (1 - tail(nsx,1)))

    #
    #.. Use stationary assumption to get nLx and infTx
    #
    nLx <- 5 * l0 * cumprod(nsx)

    #
    #.. life table
    #
    lx <- c(l0, head(nLx, -1) / n)
    ndx <- c(head(lx, -1) - tail(lx, -1), tail(lx, 1))

    #
    #.. return
    #
    out <- ndx / nLx
    if(fix.names) names(out) <- names(nsx) else names(nsx) <- NULL

    return(out)
}


nsx.to.nmx.jan07 <-
    function(nsx, l0 = 1000, n = 5, fix.names = TRUE)

    #-- Calculate nmx from surv props assuming stationary pop --#
    #
    #   nsx     :  survival proportions
    #   l0      :  radix to use for life table
    #   n       :  size of age interval
    #   fix.names
    #           :  sort out names

{
    #
    #.. Transform last nsx since it is T_{x+5} / T_x
    #
    k <- c(head(nsx, -1), tail(nsx,1) / (1 - tail(nsx,1)))

    #
    #.. Use stationary assumption to get nLx and infTx
    #
    nLx <- n * l0 * cumprod(nsx)

    #
    #.. life table
    #
    lx <- c(l0, head(nLx, -1) / n)
    ndx <- c(head(lx, -1) - tail(lx, -1), tail(lx, 1))

    #
    #.. return
    #
    out <- ndx / nLx
    if(fix.names) names(out) <- names(nsx) else names(nsx) <- NULL

    return(out)
}


### ############################################################################
### * Summarization functions

#--- Plot function ---#

pop.compare <- function(bline.mat, proj.mat, cen.mat)
  {
    require(reshape)
    require(lattice)
    proj.df <- melt(proj.mat[,-1])
    proj.df$src <- "projection"
    cen.df <- melt(cen.mat)
    cen.df$src <- "census"
    bline.df <- melt(bline.mat)
    bline.df$src <- "baseline"
    bline.df$X1 <- (cen.df$X1)[1:nrow(bline.df)]
    pl.df <- rbind(proj.df, cen.df, bline.df)
    return(pl.df)
  }


#--- Plot posterior and prior densitites for inverse gamma ---#

postPrior.plot <- function(x, prior, main = "", xlab = NULL, ...)
  {
    # x       data vector to plot
    # prior   prior density function: must take x-axis vector as
    #          first argument
    # ...     passed to prior()

    #.. x vector for prior
    xvec <- seq(from = 0, to = max(x), length = 5000)

    #.. Determine y-axis
    dx <- density(x)
    d.ymax <- max(dx$y)
    px <- prior(xvec, ...)
    p.ymax <- max(px, na.rm = TRUE)

    #.. Plot
    plot(dx, main = main, xlab = xlab
         ,ylim = c(0, max(d.ymax, p.ymax)))
    lines(xvec, px, col = "blue", lty = 2)
    legend("topright", lty = c(1, 2), col = c("black", "blue")
           ,legend = c("Posterior", "Prior")
           )
  }

postPrior.plot.sd <- function(x, prior, main = "", xlab = NULL, ...)
  {
    # x       data vector to plot
    # prior   prior density function: must take x-axis vector as
    #          first argument
    # ...     passed to prior()

    #.. x vector for prior
    sqrt.xvec <- seq(from = 0, to = 1.1*max(sqrt(x)), length = 5000)

    #.. Determine y-axis
    dx <- density(sqrt(x))
    d.ymax <- max(dx$y)
    px <- prior(sqrt.xvec^2, ...) *2*sqrt.xvec
    p.ymax <- max(px, na.rm = TRUE)

    #.. Plot
    plot(dx, main = main, xlab = xlab, xlim = c(0, 1.1*max(sqrt(x)))
         ,ylim = c(0, max(d.ymax, p.ymax)))
    lines(sqrt.xvec, px, col = "blue", lty = 2)
    legend("topright", lty = c(1, 2), col = c("black", "blue")
           ,legend = c("Posterior", "Prior")
           )
  }


#--- Mode ---#

statmod <- function(x) {
    z <- table(as.vector(x))
    names(z)[z == max(z)]
 }


#--- Calculate variance of log vital rate parameter given prop error ---#

varLogVR <- function(prop.err, prob)
{
    if(prob < 0 | prob > 1) stop("'prob' must be a valid probability")
    return((log(1-prop.err)/qnorm(p = prob, mean = 0, sd = 1))^2)
}


#--- Posterior predictive interval plots ---#

#plot.bigmcmcQuants <- function(m, q = c(0.2, 0.5, 0.8)
#                               ,t = NULL, t.name = "truth"
#                               ,t.col
#                               ,meas, m.name = "observed"
#                               ,m.col
#                               ,q.col
#                               ,ages.years, ...)
#  {
#    #--- Plot distribution of vital rate trajectories ---#

#    # m          mcmc object
#    # q          vector of quantiles to plot
#    # t          the 'true' value (if exists) subset to make same
#    #            dims as vital
#    # t.col      colour to plot "truth"
#    # meas       measured values (e.g., truth plus noise)
#    # m.col      colour to plot measured
#    # q.col      colour
#    # ages.years either "cols" or list with $ages and $years as
#    #            character vectors
#    # ...        passed to xyplot()

#    require(lattice)
#    require(reshape)

#    if(ages.years = "cols") {

#      #.. Make ages and years from column names

#      if(is.null(colnames(m))) stop("must specify ages and years")

#      colspl <- strsplit(colnames(m), ".", fixed = TRUE)
#      years <- sapply(colspl, FUN = function(z) z[1])
#      ages <- sapply(colspl, FUN = function(z) z[2])
#      ages.numeric <- gsub("[^0-9]", "", ages)

#    } else {
#      if(!is.list(ages.years)) stop("'ages.years' must be a list")

#      years <- ages.years$years
#      ages <- ages.years$ages
#      ages.numeric <- gsub("[^0-9]", "", ages)

#    }


#    #.. Set colours

#    if(missing(t.col)) t.col <- trellis.par.get("superpose.line")$col[1]
#    if(missing(m.col)) m.col <- trellis.par.get("superpose.line")$col[2]
#    if(missing(q.col)) q.col <- trellis.par.get("superpose.line")$col[3]


#    #.. Calculate quantiles

#    q.vital <- matrix(nrow = length(q), ncol = ncol(m))
#    dimnames(q.vital) <- list(as.character(q), colnames(m))

#    for(j in 1:ncol(m)) q.vital[,j] <- quantile(m[,j], probs = q)


#    #.. Prepare data sets

#    mqvit <- melt(q.vital)
#    mqvit <- rename(mqvit, c(X1 = "quant"))
#    mqvit.col <- cbind(mqvit
#                       ,year = rep(years, rep(length(q), length(years)))
#                       ,age = rep(rep(1:length(unique(ages.numeric))
#                          ,length(unique(years)))
#                          ,rep(length(q), length(ages.numeric))
#                          )
#                       )

#    if(!missing(meas)) {
#      mmeas <- melt(meas)
#      mmeas <- rename(mmeas, c(X2 = "year"))
#      mmeas.col <-
#        cbind(mmeas
#              ,year = mmeas$year
#              ,age = sapply(strsplit(as.character(mmeas$X1), "\\.")
#                 ,FUN = function(z) z[[1]])
#              ,quant = 99
#              )
#    }

#    plot.df <- rbind(mqvit.col[,c("value", "year", "age", "quant")]
#                     ,mmeas.col[,c("value", "year", "age", "quant")]
#                     )



#    #.. Plot quantiles

#    print(xyplot(value ~ age | ordered(year)
#                 ,data = plot.df
#                 ,groups = quant
#                 ,type = "b"
#                 ,xlab = "age"
#                 ,ylab = "vital"
#                 ,col = c(q.col, q.col, q.col, m.col)
#                 ,lty = c(2,1,2,1)
#                 ,key = list(text = list(c("median", "80% PI", m.name))
#                    ,lines = list(lty = c(1, 2, 1))
#                    ,col = c(q.col, q.col, m.col)
#                    ,type = c("l", "l", "l"))
#                 ,scales = list(x = list(at = 1:length(unique(ages))
#                                  ,labels = unique(ages)
#                                  ,cex = 0.5
#                                  ,rot = -20)
#                    )
#                 #,...
#                 )
#          )

#  }



####
####
####
####=================================================================
#### INVERSE GAMMA PARAMETERS FROM mae OF A T DISTN
####=================================================================
####
####
###
###
###-----------------------------------------------------------------
### MAE of a loc-scale t
###-----------------------------------------------------------------
###
MAE.t <- function(s, n)
{
    #--- MAE of a loc-scale t ---#
    #
    # DESCRIPTION:
    #   The density of the loc-scale t is:
    #      f(x) = gamma((n+1)/2)/gamma(n/2) * 1/sqrt(n * pi) *
    #               (1 + (log(x) - m)^2/(s * n))^(-(n+1)/2)
    #
    # ARGUMENTS:
    #   s:       scale parameter
    #   n:       degrees of freedom
    #

    return(1/(beta(n/2, 1/2)) * (1/(n - 1)) * 2 * sqrt(n) * s)

}

###
###
###-----------------------------------------------------------------
### Find inverse gamma parameters
###-----------------------------------------------------------------
###
IGparam.t <- function(e, s, n)
{
    #--- Find inverse gamma parameters from MAE of a loc-sale t ---#
    #
    # DESCRIPTION:
    #   Given the median and upper 0.975 quantile of
    #
    #   The density of the loc-scale t is:
    #      f(x) = gamma((n+1)/2)/gamma(n/2) * 1/sqrt(n * pi) *
    #               (1 + (log(x) - m)^2/(s * n))^(-(n+1)/2)
    #
    # ARGUMENTS:
    #   e:       mean absolute error
    #   s:       scale parameter
    #   n:       shape parameter

}






## ... Miscellaneous ...
#.....................................................................#

#.. Search for objects of particular class

lsClass <- function(...)
  {
    cl <- list(...)
    l <- lapply(cl, FUN = function(x)
                {
                  ls(envir = .GlobalEnv)[
                       sapply(as.list(ls(envir = .GlobalEnv))
                              ,FUN = function(z)
                              {
                                x %in% class(get(z))
                              }
                              )]
                }
                )
    names(l) <- as.character(cl)
    l
  }







## ... OBSOLETE functions ...
#.....................................................................#

#--- Vital rate plots USES OBSOLETE OUTPUT FORMAT ---#

#.. Many ways to improve, e.g.:
#   - take ages and years automatically
#   - cannot handle d not an array

plot.vitalQuants <- function(d, v, q = c(0.2, 0.5, 0.8)
                             ,t = NULL, t.name = "truth"
                             ,t.col = NULL
                             ,meas = NULL, m.name = "observed"
                             ,m.col = NULL
                             ,q.col = NULL
                             ,ages, years, ...)
  {
    #--- Plot distribution of vital rate trajectories ---#

    # d       the output (list) from sampler function
    # v       Vital rate to plot (must be named as in output)
    # q       vector of quantiles to plot
    # t       the 'true' value (if exists) subset to make same
    #         dims as vital
    # meas    measured values (e.g., truth plus noise)
    # ages    vector of age-group labels
    # years   vector of years
    # ...     passed to xyplot()

    require(lattice)
    require(reshape)

    if(is.null(ages)) ages <- 1:(dim(d[[v]])[1])
    if(is.null(years)) years <- 1:(dim(d[[v]])[2])

    if(is.null(t.col)) t.col <- trellis.par.get("superpose.line")$col[1]
    if(is.null(m.col)) m.col <- trellis.par.get("superpose.line")$col[2]
    if(is.null(q.col)) q.col <- trellis.par.get("superpose.line")$col[3]

    vital <- d[[v,drop=F]][,-(dim(d[[v,drop=F]])[2]), ,drop = F]
    dimnames(vital) <- list(ages, years, dimnames(vital)[3])

    q.vital <- apply(vital, c(1,2), FUN = quantile,
                     probs = q)
    dimnames(q.vital) <- list(dimnames(q.vital)[[1]]
                             ,ages
                              ,years
                              )

    if(!is.null(t) && !is.null(meas)) {
      mt <- melt(t)
      true <- cbind(rep(t.name, nrow(mt)), mt)
      colnames(true) <- c("X1", "X2", "X3", "value")
      mmeas <- melt(meas)
      measured <- cbind(rep(m.name, nrow(mmeas)), mmeas)
      colnames(measured) <- colnames(true)

      pldf.v <- rbind(melt(q.vital), true, measured)

      print(xyplot(value ~ X2 | ordered(X3)
                   ,data = pldf.v
                   ,groups = X1
                   ,type = "b"
                   ,xlab = "age"
                   ,ylab = "vital"
                   ,col = c(q.col, q.col, q.col, t.col, m.col)
                   ,lty = c(2, 1, 2, 1, 1)
                   ,key = list(text = list(c("median", "80% PI", t.name
                                 ,m.name))
                      ,lines = list(lty = c(1, 2, 1, 1))
                      ,col = c(q.col, q.col, t.col, m.col)
                      ,type = c("l", "l", "l", "l"))
                   ,...
                   )
            )
      #.. attempt to free some memory
      rm(true)
      rm(meas)

    } else if(!is.null(t)) {
      mt <- melt(t)
      true <- cbind(rep(t.name, nrow(mt)), mt)
      colnames(true) <- c("X1", "X2", "X3", "value")

      pldf.v <- rbind(melt(q.vital), true)

      print(xyplot(value ~ X2 | ordered(X3)
                   ,data = pldf.v
                   ,groups = X1
                   ,type = "b"
                   ,xlab = "age"
                   ,ylab = "vital"
                   ,col = c(q.col, q.col, q.col, m.col)
                   ,lty = c(2, 1, 2, 1)
                   ,key = list(text = list(c("median", "80% PI", m.name))
                      ,lines = list(lty = c(1, 2, 1))
                      ,col = c(q.col, q.col, m.col)
                      ,type = c("l", "l", "l"))
                   ,...
                   )
            )
      #.. attempt to free some memory
      rm(true)

    } else {
      pldf.v <- melt(q.vital)

      print(xyplot(value ~ X2 | ordered(X3)
                   ,data = pldf.v
                   ,groups = X1
                   ,type = "b"
                   ,xlab = "age"
                   ,col = c(q.col, q.col, q.col)
                   ,lty = c(2, 1, 2)
                   ,xlab = "age"
                   ,ylab = "vital"
                   ,key = list(text = list(c("median", "80% PI"))
                      ,lines = list(lty = c(1, 2))
                      ,col = c(q.col, q.col)
                      ,type = c("l", "l"))
                   ,...
                   )
            )
    }

    #.. Attempt to free memory
    rm(vital, q.vital, mt, pldf.v)
    garb <- gc(verbose = FALSE)
  }


#--- Histogram plots USES OBSOLETE OUTPUT FORMAT ---#

plot.vitalHists <- function(d, v, ages, years, thin = NULL, trim = NULL
                            ,parStripText.list = NULL)
  {
    #--- Plot distribution of vital rate trajectories ---#

    # d       the output (list) from sampler function
    # t       the 'true' value (if exists) subset to make same
    #         dims as vital
    # v       Vital rate to plot (must be named as in output)
    # ages    vector of age-group labels
    # years   vector of years

    require(lattice)
    require(reshape)

    if(is.null(ages)) ages <- 1:(dim(d[[v]])[1])
    if(is.null(years)) years <- 1:(dim(d[[v]])[2])
    if(is.null(thin)) {
      thin.v <-
        round(seq(from = 1, to = dim(d[[v]])[3], by = dim(d[[v]])[3]*0.1))
    } else {
      thin.v <- seq(from = 1, to = dim(d[[v]])[3], by = thin)
    }

    vital <- d[[v,drop=F]][,-(dim(d[[v,drop=F]])[2]), ,drop = F]
    dimnames(vital) <- list(ages, years, 1:dim(vital)[3])

    vital.thinned <- vital[,,thin.v]

    mv <- melt(vital.thinned)
    if(is.null(trim)) trim <- max(mv$value)

    print(histogram(~ value | ordered(X1) + ordered(X2)
                    ,data = mv
                    ,subset = value <= trim
                    ,xlab = "vital"
                    ,par.strip.text = parStripText.list
                    )
          )

    #.. Attempt to free memory
    rm(vital, vital.thinned, mv)
    garb <- gc(verbose = FALSE)
  }
