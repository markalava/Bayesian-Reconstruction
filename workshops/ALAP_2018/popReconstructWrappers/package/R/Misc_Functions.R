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

invGam.paramGen <-
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


