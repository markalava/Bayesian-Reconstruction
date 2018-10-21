################################################################################
###
### DATE OF ORIGINAL:   8th August 2012
###
### AUTHOR:             Mark Wheldon
###
### DESC:               Two-sex version of Bayesian population reconstruction
###                     with SRB also estimated.
###
###-----------------------------------------------------------------------------
###
################################################################################

### * HELPER FUNCTIONS
################################################################################

ls.before <- ls()


### ** Misc Functions
################################################################################

logit <- function(p) log(p / (1 - p))

invlogit <- function(x)
{
    if(any(is.infinite(exp(x)))) {
        y <- x
        y[is.infinite(exp(x))] <- 1
        y[!is.infinite(exp(x))] <-
            invlogit(y[!is.infinite(exp(x))])
        return(y)
    }
    else return(exp(x) / (1 + exp(x)))
}


##--- Generates random draws from inverse gamma ---##

rinvGamma <- function(n, shape, scale)
{
    return(1/rgamma(n, shape = shape, rate = scale))
}


##--- Returns value of inverse gamma pdf ---##

dinvGamma <- function(x, shape, scale, log = FALSE)
{
    if(log) d <-
        shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
    else d <- scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
    return(d)
}


##--- Creates column names for mcmc objects ---##

makeColNames <- function(m)
{
    ## m:    matrix of input values

    e <- expand.grid(rownames(m), colnames(m))
    apply(e, 1, FUN = function(z) paste(z[2], z[1], sep = "."))

}


### ** Projection for census years
################################################################################

proj.cen.yrs <-
    function(full.proj, bline.yr, vr.yrs, cen.yrs, proj.yrs
             ,labels = FALSE)

{

    ##-------* HOUSKEEPING

    bline.yr <- as.numeric(bline.yr)
    vr.yrs <- as.numeric(vr.yrs)
    cen.yrs <- as.numeric(cen.yrs)
    proj.yrs <- as.numeric(proj.yrs)

    interp.counts <-
        list(female = matrix(NA, nrow = nrow(full.proj[["female"]])
             ,ncol = length(cen.yrs))
             ,male = matrix(NA, nrow = nrow(full.proj[["male"]])
              ,ncol = length(cen.yrs))
             )

    if(labels) {
        if(!is.null(rownames(full.proj[["female"]]))) {
            rownames(interp.counts[["female"]]) <- rownames(full.proj[["female"]])
        } else {
            rownames(interp.counts[["female"]]) <- seq(from = 0, by = 5
                                       ,length.out = nrow(full.proj[["female"]])
                                       )
        }
        if(!is.null(rownames(full.proj[["male"]]))) {
            rownames(interp.counts[["male"]]) <- rownames(full.proj[["male"]])
        } else {
            rownames(interp.counts[["male"]]) <- seq(from = 0, by = 5
                                       ,length.out = nrow(full.proj[["male"]])
                                       )
        }
        colnames(interp.counts[["female"]]) <- cen.yrs
        colnames(interp.counts[["male"]]) <- cen.yrs
    }

    cen.in.proj <- cen.yrs %in% proj.yrs
    cen.notin.proj <- !cen.in.proj
    proj.in.cen <- proj.yrs %in% cen.yrs


    ##------- * INTERPOLATE IF NECESSARY

    ## Check to see if interpolation necessary

    if(all(cen.in.proj)) {
        interp.counts <- lapply(full.proj, function(z) z[,proj.in.cen,drop=FALSE])
        if(labels) {
            dimnames(interp.counts[["female"]]) <-
                list(rownames(full.proj[["female"]]), cen.yrs)
            dimnames(interp.counts[["male"]]) <-
                list(rownames(full.proj[["male"]]), cen.yrs)
        }
    } else {

        ##-------** Female

        ## Growth rates
        ## (Does this for all, including those not needing interpolation)

        ratio.mat <- matrix(NA,
                            ,nrow = nrow(full.proj[["female"]])
                            ,ncol = ncol(full.proj[["female"]]) - 1)
        for(j in 1:ncol(ratio.mat)) {
            ratio.mat[,j] <-
                full.proj[["female"]][,(j+1)] / full.proj[["female"]][,j]
        }

        yr.diffs <- diff(proj.yrs)

        ## if(any(is.nan(log(ratio.mat)))) {
        ##     cat("\nratio.mat=\n");  print(ratio.mat)
        ##     stop()
        ## }
        r.mat <- log(ratio.mat) / yr.diffs


        ## Duplicate last column of r.mat so that "interpolation"
        ## can be done for years beyond the end of proj.yrs,
        ## assuming that the growth rate remains constant
        ## (although don't know why we'd do this)

        r.mat <- cbind(r.mat, r.mat[,ncol(r.mat)])


        ## Interpolate

        for(j in which(cen.notin.proj)) {
            base.yr.ind <-
                max(which(proj.yrs <= cen.yrs[j]))
            time.gap <- cen.yrs[j] - proj.yrs[base.yr.ind]
            growth <- r.mat[,base.yr.ind]
            interp.counts[["female"]][,j] <-
                full.proj[["female"]][,base.yr.ind] * exp(growth * time.gap)
        }
        interp.counts[["female"]][,cen.in.proj] <-
            full.proj[["female"]][,proj.in.cen]

        if(labels) {
            dimnames(interp.counts[["female"]]) <-
                list(rownames(full.proj[["female"]]), cen.yrs)
        }


        ##-------** Male

        ## Growth rates
        ## (Does this for all, including those not needing interpolation)

        ratio.mat <- matrix(NA,
                            ,nrow = nrow(full.proj[["male"]])
                            ,ncol = ncol(full.proj[["male"]]) - 1)
        for(j in 1:ncol(ratio.mat)) {
            ratio.mat[,j] <-
                full.proj[["male"]][,(j+1)] / full.proj[["male"]][,j]
        }

        yr.diffs <- diff(proj.yrs)

        ## if(any(is.nan(log(ratio.mat)))) {
        ##     cat("\nratio.mat=\n");  print(ratio.mat)
        ##     stop()
        ## }
        r.mat <- log(ratio.mat) / yr.diffs


        ## Duplicate last column of r.mat so that "interpolation"
        ## can be done for years beyond the end of proj.yrs,
        ## assuming that the growth rate remains constant
        ## (although don't know why we'd do this)

        r.mat <- cbind(r.mat, r.mat[,ncol(r.mat)])


        ## Interpolate

        for(j in which(cen.notin.proj)) {
            base.yr.ind <-
                max(which(proj.yrs <= cen.yrs[j]))
            time.gap <- cen.yrs[j] - proj.yrs[base.yr.ind]
            growth <- r.mat[,base.yr.ind]
            interp.counts[["male"]][,j] <-
                full.proj[["male"]][,base.yr.ind] * exp(growth * time.gap)
        }
        interp.counts[["male"]][,cen.in.proj] <-
            full.proj[["male"]][,proj.in.cen]

        if(labels) {
            dimnames(interp.counts[["male"]]) <-
                list(rownames(full.proj[["male"]]), cen.yrs)
        }

    }

    return(interp.counts)

}


### ** Likelihood
################################################################################

propRecon.log.lhood <-
    function(log.n.census, log.n.hat, ll.var)
{
    ##.. log.n.census and log.n.hat should already be logged

    ##.. log.n.hat should be log projected counts for census years
    ##   interpolated if necessary


    ##-- value of log likelihoods --##

    density <- unlist(mapply(FUN = "dnorm", x = log.n.census, mean = log.n.hat
                             ,MoreArgs = list(sd = sqrt(ll.var), log = TRUE)
                             ,SIMPLIFY = FALSE
                             )
                      ,use.names = FALSE
                      )

    ##-- joint log likelihood --##

    return(sum(density))
}


### ** Prior
################################################################################

log.prior <-
    function(## estimated vitals (transformed)
             f, s, g, baseline.n, srb
             ## fixed prior means on vitals (transformed)
             ,prior.mean.f, prior.mean.s
             ,prior.mean.g, prior.mean.b
             ,prior.mean.srb
             ## fixed prior parameters on variance distns
             ,alpha.f, beta.f, alpha.s, beta.s
             ,alpha.g, beta.g
             ,alpha.n, beta.n
             ,alpha.srb, beta.srb
             ## updated variances on prior distns
             ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n, sigmasq.srb
             ## non zero rows of fertility matrix
             ,non.zero.fert
             )
{


    ##-- Values of prior densities for vitals --##

    ##.. Note that log densities are calculated for numerical stability.
    ##     f, baseline.n, prior.mean.f, prior.mean.b are logged coming
    ##     in, s, prior.mean.s is logit transformed coming in, g and
    ##     prior.mean.g are not transformed coming in.
    log.f.prior <- dnorm(as.vector(f[non.zero.fert,])
                         ,mean = as.vector(prior.mean.f[non.zero.fert,])
                         ,sd = sqrt(sigmasq.f)
                         ,log = TRUE)
    log.s.prior <- unlist(mapply(FUN = "dnorm", x = s, mean = prior.mean.s
                                 ,MoreArgs = list(sd = sqrt(sigmasq.s), log = TRUE)
                                 ,SIMPLIFY = FALSE
                                 )
                          ,use.names = FALSE
                          )
    log.g.prior <- unlist(mapply(FUN = "dnorm", x = g, mean = prior.mean.g
                                 ,MoreArgs = list(sd = sqrt(sigmasq.g), log = TRUE)
                                 ,SIMPLIFY = FALSE
                                 )
                          ,use.names = FALSE
                          )
    log.b.prior <- unlist(mapply(FUN = "dnorm", x = baseline.n, mean = prior.mean.b
                                 ,MoreArgs = list(sd = sqrt(sigmasq.n), log = TRUE)
                                 ,SIMPLIFY = FALSE
                                 )
                          ,use.names = FALSE
                          )
    log.srb.prior <-
        dnorm(drop(srb), mean = drop(prior.mean.srb), sd = sqrt(sigmasq.srb), log = TRUE)


    ##-- Values of prior densities for variances --##

    log.sigmasq.f.prior <-
        dinvGamma(sigmasq.f, alpha.f, beta.f, log = TRUE)
    log.sigmasq.s.prior <-
        dinvGamma(sigmasq.s, alpha.s, beta.s, log = TRUE)
    log.sigmasq.g.prior <-
        dinvGamma(sigmasq.g, alpha.g, beta.g, log = TRUE)
    log.sigmasq.n.prior <-
       dinvGamma(sigmasq.n, alpha.n, beta.n, log = TRUE)
    log.sigmasq.srb.prior <-
       dinvGamma(sigmasq.srb, alpha.srb, beta.srb, log = TRUE)


    ##-- The log posterior is the SUM of these --##

    return(sum(log.f.prior, log.s.prior, log.g.prior, log.b.prior, log.srb.prior
               ,log.sigmasq.f.prior
               ,log.sigmasq.s.prior
               ,log.sigmasq.g.prior
               ,log.sigmasq.n.prior
               ,log.sigmasq.srb.prior
               ))

}


### ** Posterior
################################################################################

log.post <-
    function(## estimated vitals
             f, s, g, baseline.n, srb
             ## fixed prior means on vitals
             ,prior.mean.f, prior.mean.s
             ,prior.mean.g, prior.mean.b
             ,prior.mean.srb
             ## fixed prior parameters on variance distns
             ,alpha.f, beta.f, alpha.s, beta.s
             ,alpha.g, beta.g
             ,alpha.n, beta.n, alpha.srb, beta.srb
             ## updated variances on prior distns
             ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n, sigmasq.srb
             ## value of the log likelihood
             ,log.like
             ## non zero rows of fertility matrix
             ,non.zero.fert
             )
{

    ##-- The log posterior is the SUM of log.prior and log.like --##

    return(sum(log.prior(f, s, g, baseline.n, srb
                                        ,prior.mean.f, prior.mean.s, prior.mean.g
                                        ,prior.mean.b, prior.mean.srb
                                        ,alpha.f, beta.f, alpha.s, beta.s, alpha.g, beta.g
                                        ,alpha.n, beta.n, alpha.srb, beta.srb
                                        ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n
                                        ,sigmasq.srb
                               ,non.zero.fert
                               )
               ,log.like))

}


### ** Acceptance Ratio
################################################################################

acc.ra <- function(log.prop, log.current)
{
    min(1, exp(log.prop - log.current))
}

acc.ra.var <-
    function(log.prop.post, log.curr.post, log.prop.var, log.curr.var)
{
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post
               ))
}

acc.ra.srb <-
    function(log.prop.post, log.curr.post, log.prop.srb, log.curr.srb)
{
    min(1, exp(log.prop.post + log.curr.srb - log.curr.post - log.prop.srb))
}


### * SAMPLER
################################################################################

##' Generate MCMC samples for two-sex Bayesian reconstruction
##'
##' The main function that does the estimation of parameters for the reconstruction.
##'
##' @param n.iter,burn.in Number of saved iterations and burn-in iterations.
##' @param al.f,be.f,al.s,be.s,al.g,be.g,al.n,be.n,al.srb,be.srb Fixed hyperparameters (same for male and female). See \code{\link{make.hyper.params}}.
##' @param mean.f,mean.s,mean.g,mean.b,mean.srb Fixed prior means. Lists with \code{female} and \code{male} components.
##' @param start.f,start.s, start.g, start.b, start.srb Inital values for vitals. Not transformed coming in. Lists with \code{female} and \code{male} components.
##' @param start.sigmasq.f,start.sigmasq.s,start.sigmasq.g,start.sigmasq.n,start.sigmasq.srb Inital values for variances (same for male and female).
##' @param pop.data Census data. Not transformed coming in. List with \code{female} and \code{male} components.
##' @param fert.rows Rows of \code{mean.f} and \code{start.f} that are non-zero.
##' @param proj.periods Number of periods to project forward over (e.g., number of five-year steps).
##' @param age.size Age group width.
##' @param verb,progress.step print algorithm progress every 'progress.step' iterations
##' @param prop.varcovar Variances for proposal distributions used in M-H steps which update vital rates. Lists with \dQuote{female} and \dQuote{male} components.
##' @param save.mid.run.every,save.mid.run.name How often should the chain be saved, and with what name?
##' @return List with components for CCMPP and other algorithm parameters.
##' @author Mark C. Wheldon
##' @export
pop.est.sampler <- function(n.iter, burn.in = 0,
             al.f, be.f, al.s, be.s, al.g, be.g, al.n, be.n, al.srb, be.srb,
             mean.f, mean.s, mean.g, mean.b, mean.srb,
             start.f, start.s, start.g, start.b, start.srb,
             start.sigmasq.f, start.sigmasq.s, start.sigmasq.g,
             start.sigmasq.n, start.sigmasq.srb,
             pop.data,
             fert.rows = as.logical(apply(mean.f == 0L, 1, function(z) !all(z))),
             proj.periods = ncol(mean.f),

             ## age group width
             age.size = 5,

             ## print algorithm progress every 'progress.step' iterations
             verb = FALSE, progress.step = 1000,

             ## **variances** for proposal distributions used in M-H steps which
             ## update vital rates.
             ## (lists with "female" and "male" components)
             prop.varcovar,

             ## save mid run
             save.mid.run.every = progress.step,
             save.mid.run.name = NULL
             ) {

    ## -------* HOUSEKEEPING

    ## Begin timing
    ptm <- proc.time()

    ## CCMP function
    ccmp.function <- ccmp.femDom
    mig.string <- switch(EXPR = formals(ccmp.function)$mig.type[[2]]
                         ,prop.prev.pop = "prop"
                        ,net.count = "count")

    ## Other variables

    ## should chains be kept as filebacked matrices
    ## using package ff ?
    filebacked.chains <- FALSE
    file.suff <- ""

    ## tolerance defining allowable survival probabilities
    s.tol <-  10^(-10)

    ## variance proposal type
    v.prop.type <- "MH.scaled"

    ## -------** Check inputs

    ## DIMENSIONS
    ##

    ## Female
    input.dims <-
        sapply(list(start.s = start.s[["female"]][1:(nrow(start.s[["female"]])-1),]
                              , start.g = start.g[["female"]]
                              ,mean.f = mean.f
                              , mean.s = mean.s[["female"]][1:(nrow(mean.s[["female"]])-1),]
                              ,mean.g = mean.g[["female"]]
                              )
                         ,"dim")
    mismatch.dims <- apply(input.dims, 2, "identical", dim(start.f))
    if(!all(mismatch.dims))
        stop("Dims of all female inputs do not match 'dim(start.f)'", "\n"
             ,paste(names(mismatch.dims)[!mismatch.dims], collapse = "  "))

    ## Male
    input.dims <- sapply(list(start.s = start.s[["male"]][1:(nrow(start.s[["male"]])-1),]
                              , start.g = start.g[["male"]]
                              ,mean.s = mean.s[["male"]][1:(nrow(mean.s[["male"]])-1),]
                              ,mean.g = mean.g[["male"]]
                              )
                         ,"dim")
    mismatch.dims <- apply(input.dims, 2, "identical", dim(start.f))
    if(!all(mismatch.dims))
        stop("Dims of all male inputs do not match 'dim(start.f)'", "\n"
             ,paste(names(mismatch.dims)[!mismatch.dims], collapse = "  "))


    ## YEARS
    ##

    ## Female
    all.vr.years <-
        list(start.s = colnames(start.s[["female"]])
             , start.g = colnames(start.g[["female"]])
             ,mean.f = colnames(mean.f), mean.s = colnames(mean.s[["female"]])
             ,mean.g = colnames(mean.g[["female"]]), srb = colnames(start.srb)
             )
    mismatch.yrs <- sapply(all.vr.years, "identical", colnames(start.f))
    if(!all(mismatch.dims))
        stop("Years of these inputs do not match years of 'start.f'", "\n"
             ,paste(names(mismatch.yrs)[!mismatch.yrs], collapse = "  "))

    all.vr.years.eq <-
        sapply(all.vr.years, FUN = function(z) all.equal(colnames(start.f), z))
    if(!all(all.vr.years.eq))
        stop(paste("colnames(", names(all.vr.years)[min(which(!all.vr.years.eq))], ")"
                   ," != colnames(start.f). There may be more...", sep = "")
             )

    ## Male
    all.vr.years <-
        list(start.s = colnames(start.s[["male"]]), start.g = colnames(start.g[["male"]])
             ,mean.s = colnames(mean.s[["male"]])
             ,mean.g = colnames(mean.g[["male"]])
             )
    mismatch.yrs <- sapply(all.vr.years, "identical", colnames(start.f))
    if(!all(mismatch.dims))
        stop("Years of these inputs do not match years of 'start.f'", "\n"
             ,paste(names(mismatch.yrs)[!mismatch.yrs], collapse = "  "))

    all.vr.years.eq <-
        sapply(all.vr.years, FUN = function(z) all.equal(colnames(start.f), z))
    if(!all(all.vr.years.eq))
        stop(paste("colnames(", names(all.vr.years)[min(which(!all.vr.years.eq))], ")"
                   ," != colnames('start.f'). There may be more...", sep = "")
             )


    ## DETERMINE VARIOUS SETS OF YEARS
    ##

    proj.years <-
        seq(from = as.numeric(colnames(start.f)[1]), by = age.size
            ,length = proj.periods + 1)

    if(!all.equal(proj.years[1:ncol(start.f)]
                  ,as.numeric(colnames(start.f))))
        stop("colnames(start.f) !=  seq(from = as.numeric(colnames(start.f)[1]), by = age.size, length = ncol(start.f))")

    vr.years <- as.numeric(colnames(start.f))
    baseline.year <- as.numeric(colnames(start.b[["female"]]))
    census.years <- as.numeric(colnames(pop.data[["female"]]))


    ## -------* STORAGE

    ## -------** Objects

    ## MCMC objects for posterior samples
    ## Samples are stored as 2D arrays for compatibility with coda's
    ## mcmc format with iterations as rows, year*age.group as columns.
    ## Age.group cycles fastest across columns, e.g.,
    ## _____________________________________________________
    ##   1960  | 1960  | 1960  | ... | 1965  | 1965  | ...
    ##   15.19 | 20.24 | 25.29 | ... | 15.19 | 20.24 | ...
    ## 1  --   |  --   |  --   | ... |  --   |  --   | ...
    ## 2  --   |  --   |  --   | ... |  --   |  --   | ...
    ##   etc.
    ## _____________________________________________________

    if(filebacked.chains) {

        warning("FILEBACKED.CHAINS NOT WELL TESTED")

        ## -------*** File backed

        ## Store chains as 'ff' objects. Turn into ff.mcmc at end.

        fert.rate.mcmc <-
            ff(dim = c(n.iter, sum(fert.rows) * ncol(start.f))
               ,dimnames = list(NULL
                ,makeColNames(start.f[fert.rows,]))
               ,vmode = "double"
               ,filename = paste("fert_mcmc", file.suff, ".bin"
                ,sep = "")
               ,overwrite = TRUE
               ,finalizer = "close"
               )

        surv.prop.mcmc <-
            list(female =
                 ff(dim = c(n.iter, nrow(start.s[["female"]]) * ncol(start.s[["female"]]))
                    ,dimnames = list(NULL
                     ,makeColNames(start.s[["female"]]))
                    ,vmode = "double"
                    ,filename = paste("surv_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 ,male =
                 ff(dim = c(n.iter, nrow(start.s[["male"]]) * ncol(start.s[["male"]]))
                    ,dimnames = list(NULL
                     ,makeColNames(start.s[["male"]]))
                    ,vmode = "double"
                    ,filename = paste("surv_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 )

        lx.mcmc <-
            list(female =
                 ff(dim = c(n.iter, nrow(pop.data[["female"]]) * (proj.periods))
                    ,dimnames = list(NULL
                     ,makeColNames(matrix(0, nrow = nrow(start.b[["female"]])
                                                         ,ncol = proj.periods
                                                         ,dimnames = list(rownames(start.b[["female"]])
                                                          ,seq(from = as.numeric(colnames(start.b[["female"]])[1]) +
                                                               age.size, by = age.size, length = proj.periods))
                                                         )
                                                  )
                     )
                    ,vmode = "double"
                    ,filename = paste("lx_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 ,male =
                 ff(dim = c(n.iter, nrow(pop.data[["male"]]) * (proj.periods))
                    ,dimnames = list(NULL
                     ,makeColNames(matrix(0, nrow = nrow(start.b[["male"]])
                                                         ,ncol = proj.periods
                                                         ,dimnames = list(rownames(start.b[["male"]])
                                                          ,seq(from = as.numeric(colnames(start.b[["male"]])[1]) +
                                                               age.size, by = age.size, length = proj.periods))
                                                         )
                                                  )
                     )
                    ,vmode = "double"
                    ,filename = paste("lx_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 )

        mig.mcmc <-
            list(female =
                 ff(dim = c(n.iter, nrow(start.g[["female"]]) * ncol(start.g[["female"]]))
                    ,dimnames = list(NULL
                     ,makeColNames(start.g[["female"]]))
                    ,vmode = "double"
                    ,filename = paste("migProp_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 ,male =
                 ff(dim = c(n.iter, nrow(start.g[["male"]]) * ncol(start.g[["male"]]))
                    ,dimnames = list(NULL
                     ,makeColNames(start.g[["male"]]))
                    ,vmode = "double"
                    ,filename = paste("migProp_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 )

        baseline.count.mcmc <-
            list(female =
                 ff(dim = c(n.iter, nrow(start.b[["female"]]))
                    ,dimnames = list(NULL, makeColNames(start.b[["female"]]))
                    ,vmode = "double"
                    ,filename = paste("baseline_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 ,male =
                 ff(dim = c(n.iter, nrow(start.b[["male"]]))
                    ,dimnames = list(NULL, makeColNames(start.b[["male"]]))
                    ,vmode = "double"
                    ,filename = paste("baseline_mcmc", file.suff, ".bin"
                     ,sep = "")
                    ,overwrite = TRUE
                    ,finalizer = "close"
                    )
                 )

        variances.mcmc <-
            ff(dim = c(n.iter, 5)
               ,dimnames = list(NULL, c("fert.rate.var", "surv.prob.var"
                ,"mig.var", "population.count.var", "srb.var"))
               ,vmode = "double"
               ,filename = paste("variances_mcmc", file.suff, ".bin"
                ,sep = "")
               ,overwrite = TRUE
               ,finalizer = "close"
               )

        srb.mcmc <-
            ff(dim = c(n.iter, ncol(start.srb))
               ,dimnames = list(NULL, colnames(start.srb))
               ,vmode = "double"
               ,filename = paste("srb_mcmc", file.suff, ".bin"
                ,sep = "")
               ,overwrite = TRUE
               ,finalizer = "close"
               )

    } else {

        ## -------*** Store in memory

        ## Fertility
        fert.rate.mcmc <-
            coda::mcmc(matrix(nrow = n.iter
                        ,ncol = sum(fert.rows) * ncol(start.f)
                        ,dimnames = list(NULL
                         ,makeColNames(start.f[fert.rows,]))
                        ))

        ## Survival proportions
        surv.prop.mcmc <-
            list(female = coda::mcmc(matrix(nrow = n.iter
                 ,ncol = nrow(start.s[["female"]]) * ncol(start.s[["female"]])
                 ,dimnames = list(NULL
                  ,makeColNames(start.s[["female"]])
                  )))
                 ,male = coda::mcmc(matrix(nrow = n.iter
                  ,ncol = nrow(start.s[["female"]]) * ncol(start.s[["female"]])
                  ,dimnames = list(NULL
                   ,makeColNames(start.s[["female"]])
                   )))
                 )

        ## lx
        lx.mcmc <-
            list(female =
                 coda::mcmc(matrix(nrow = n.iter
                             ,ncol = nrow(start.b[["female"]]) * (proj.periods)
                             ,dimnames = list(NULL
                              ,makeColNames(matrix(0, nrow = nrow(start.b[["female"]])
                                                                  ,ncol = proj.periods
                                                                  ,dimnames = list(rownames(start.b[["female"]])
                                                                   ,seq(from = as.numeric(colnames(start.b[["female"]])[1]) +
                                                                        age.size, by = age.size, length = proj.periods))
                                                                  )
                                                           )
                              )))
                 ,male =
                 coda::mcmc(matrix(nrow = n.iter
                             ,ncol = nrow(start.b[["male"]]) * (proj.periods)
                             ,dimnames = list(NULL
                              ,makeColNames(matrix(0, nrow = nrow(start.b[["male"]])
                                                                  ,ncol = proj.periods
                                                                  ,dimnames = list(rownames(start.b[["male"]])
                                                                   ,seq(from = as.numeric(colnames(start.b[["male"]])[1]) +
                                                                        age.size, by = age.size, length = proj.periods))
                                                                  )
                                                           )
                              )))
                 )

        ## migration proportions
        mig.mcmc <-
            list(female =
                 coda::mcmc(matrix(nrow = n.iter
                             ,ncol = nrow(start.g[["female"]]) * ncol(start.g[["female"]])
                             ,dimnames = list(NULL
                              ,makeColNames(start.g[["female"]]))
                             ))
                 ,male =
                 coda::mcmc(matrix(nrow = n.iter
                             ,ncol = nrow(start.g[["female"]]) * ncol(start.g[["female"]])
                             ,dimnames = list(NULL
                              ,makeColNames(start.g[["female"]]))
                             ))
                 )

        ## baseline counts
        baseline.count.mcmc <-
            list(female = coda::mcmc(matrix(nrow = n.iter, ncol = nrow(start.b[["female"]])
                 ,dimnames = list(NULL
                  ,makeColNames(start.b[["female"]]))
                 ))
                 ,male = coda::mcmc(matrix(nrow = n.iter, ncol = nrow(start.b[["male"]])
                  ,dimnames = list(NULL
                   ,makeColNames(start.b[["male"]]))
                  ))
                 )

        ## srb
        srb.mcmc <-
            coda::mcmc(matrix(nrow = n.iter, ncol = ncol(start.srb)
                        ,dimnames = list(NULL, colnames(start.srb))
                        ))

        ## variances
        variances.mcmc <-
            coda::mcmc(matrix(nrow = n.iter, ncol = 5))
        colnames(variances.mcmc) <-
            c("fert.rate.var", "surv.prop.var", "mig.var"
              ,"population.count.var", "srb.var")

    }


    ## -------*** Monitors

    ## Record acceptance rate
    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s[["female"]])
              ,ncol = ncol(mean.s[["female"]])
              ,dimnames = dimnames(mean.s[["female"]])
              )
             ,mig = matrix(0, nrow = nrow(mean.g[["female"]])
              , ncol = ncol(mean.g[["female"]])
              ,dimnames = dimnames(mean.g[["female"]])
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b[["female"]])
              ,dimnames = dimnames(mean.b[["female"]])
              )
             ,srb = matrix(0, ncol = ncol(start.srb)
              ,dimnames = list(NULL, colnames(start.srb)))
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.g = 0
             ,sigmasq.n = 0
             ,sigmasq.srb = 0
             )


    ## Count how often acceptance ratio missing or na
    ar.na <- acc.count

    ## Count how often projection gives negative population
    pop.negative <- acc.count


    ## Count how often surv probs are outside tolerance
    s.out.tol <- matrix(0, nrow = nrow(mean.s[["female"]])
                        , ncol = ncol(mean.s[["female"]])
                        ,dimnames = dimnames(mean.s[["female"]]))


    ## -------** Initialize

    ## Set current vitals and variances to inital values.
    ## Take logs/logits here where required.
    log.curr.f <- log(start.f) #<-- log(0) stored as "-Inf". Gets
    logit.curr.s <-
        list(female = logit(start.s[["female"]])
             ,male = logit(start.s[["male"]])
             )
    curr.g <- list(female = start.g[["female"]], male = start.g[["male"]])
    log.curr.b <- list(female = log(start.b[["female"]]), male = log(start.b[["male"]]))
    log.curr.srb <- log(start.srb)

    curr.sigmasq.f <- start.sigmasq.f
    curr.sigmasq.s <- start.sigmasq.s
    curr.sigmasq.g <- start.sigmasq.g
    curr.sigmasq.n <- start.sigmasq.n
    curr.sigmasq.srb <- start.sigmasq.srb


    ## Fixed means for vitals and baseline
    ## Set these to inputs, take logs where required.
    log.mean.f <- log(mean.f)
    logit.mean.s <-
        list(female = logit(mean.s[["female"]])
             ,male = logit(mean.s[["male"]])
             )
    mean.g <- list(female = mean.g[["female"]], male = mean.g[["male"]])
    log.mean.b <- list(female = log(mean.b[["female"]]), male = log(mean.b[["male"]]))
    log.mean.srb <- log(mean.srb)


    ## Fixed census data
    ## Take logs here
    log.census.mat <-
        list(female = log(pop.data[["female"]]), male = log(pop.data[["male"]]))


    ## Set current projection: based on initial values
    log.curr.proj <-
        lapply(proj.cen.yrs(full.proj =
                                           ccmp.function(pop = lapply(log.curr.b, "exp")
                                                         ,surv = lapply(logit.curr.s
                                                          ,"invlogit")
                                                         ,fert = exp(log.curr.f)
                                                         ,mig = curr.g
                                                         ,proj.steps = proj.periods
                                                         ,age.int = age.size
                                                         ,srb = exp(log.curr.srb)
                                                         )
                                           ,bline.yr = baseline.year
                                           ,vr.yrs = vr.years
                                           ,cen.yrs = census.years, proj.yrs = proj.years
                                           ,labels = FALSE
                                           )
               ,"log")


    ## Current log posterior
    log.curr.posterior <-
        log.post(f = log.curr.f
                                ,s = logit.curr.s
                                ,g = curr.g
                                ,baseline.n = log.curr.b
                                ,srb = log.curr.srb
                                ,prior.mean.f = log.mean.f
                                ,prior.mean.s = logit.mean.s
                                ,prior.mean.g = mean.g
                                ,prior.mean.b = log.mean.b
                                ,prior.mean.srb = log.mean.srb
                                ,alpha.f = al.f, beta.f = be.f
                                ,alpha.s = al.s, beta.s = be.s
                                ,alpha.g = al.g, beta.g = be.g
                                ,alpha.n = al.n, beta.n = be.n
                                ,alpha.srb = al.srb, beta.srb = be.srb
                                ,sigmasq.f = curr.sigmasq.f
                                ,sigmasq.s = curr.sigmasq.s
                                ,sigmasq.g = curr.sigmasq.g
                                ,sigmasq.n = curr.sigmasq.n
                                ,sigmasq.srb = curr.sigmasq.srb
                                ,non.zero.fert = fert.rows
                                ,log.like = propRecon.log.lhood(
                                 log.n.census = log.census.mat
                                 ,log.n.hat = log.curr.proj
                                 ,ll.var = curr.sigmasq.n)
                                )


    ## -------* BEGIN LOOP

    if(verb) {
        cat("\n\ntotal iterations = ", n.iter+burn.in
            ,"\nburn in = ", burn.in
            ,", stored = ", n.iter, sep = "")
        cat("\n\nfert.rows = ", which(fert.rows)
            ,"\ncensus years = ", census.years
            ,"\nvital rate years = ", vr.years
            ,"\nprojection years = ", proj.years
            )
        cat("\n\n"
            ,"iter ", " quantity\n", "---- ", " --------", "\n"
            ,sep = "")
    }

    for(i in 1:(n.iter + burn.in)) {

        ## Burn.in is not stored so create index for assignment later
        k <- i - burn.in


        ## -------** Vital Rate M-H Steps

        ## -------*** Fertility

        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " Fertility")

        ## Proposal

        ## cycle through components
        for(j in 1:length(log.curr.f[fert.rows,])) {

            ## make a matrix conformable w fertility rate matrix
            log.prop.f.mat <-
                matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
            log.prop.f.mat[fert.rows,][j] <-
                rnorm(1, 0, sqrt(prop.varcovar$fert.rate[[j]][1]))

            ## make proposal
            log.prop.f <- log.curr.f + log.prop.f.mat

            ## Run CCMP (project on the original scale)
            ## ** Don't allow negative population
            full.proj <- ccmp.function(pop = lapply(log.curr.b, "exp")
                                       ,fert = exp(log.prop.f) #<-- use proposal
                                       ,surv = lapply(logit.curr.s, "invlogit")
                                       ,mig = curr.g
                                       ,srb = exp(log.curr.srb)
                                       ,proj.steps = proj.periods
                                       ,age.int = age.size
                                       ,label.dims = FALSE
                                       ,base.year = "1985"
                                       )

            full.proj.unlist <- unlist(full.proj, use.names = FALSE)

            if(sum(full.proj.unlist < 0) > 0 || is.na(sum(full.proj.unlist)) ||
               is.nan(sum(full.proj.unlist))) {
                if(i > burn.in) {
                    pop.negative[["female"]]$fert.rate[j] <-
                        pop.negative[["female"]]$fert.rate[j] + 1/n.iter
                }
            } else {
                prop.proj <-
                    proj.cen.yrs(full.proj = full.proj
                                                ,bline.yr = baseline.year
                                                ,vr.yrs = vr.years
                                                ,cen.yrs = census.years, proj.yrs = proj.years
                                                ,labels = FALSE
                                                )
                log.prop.proj <- lapply(prop.proj, "log")

                ## Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
                    log.post(f = log.prop.f #<-- use proposal
                                            ,s = logit.curr.s
                                            ,g = curr.g
                                            ,baseline.n = log.curr.b
                                            ,srb = log.curr.srb
                                            ,prior.mean.f = log.mean.f
                                            ,prior.mean.s = logit.mean.s
                                            ,prior.mean.g = mean.g
                                            ,prior.mean.b = log.mean.b
                                            ,prior.mean.srb = log.mean.srb
                                            ,alpha.f = al.f, beta.f = be.f
                                            ,alpha.s = al.s, beta.s = be.s
                                            ,alpha.g = al.g, beta.g = be.g
                                            ,alpha.n = al.n, beta.n = be.n
                                            ,alpha.srb = al.srb, beta.srb = be.srb
                                            ,sigmasq.f = curr.sigmasq.f
                                            ,sigmasq.s = curr.sigmasq.s
                                            ,sigmasq.g = curr.sigmasq.g
                                            ,sigmasq.n = curr.sigmasq.n
                                            ,sigmasq.srb = curr.sigmasq.srb
                                            ,non.zero.fert = fert.rows
                                            ,log.like = propRecon.log.lhood(
                                             log.n.census = log.census.mat
                                             ,log.n.hat = log.prop.proj #<-- use proposal
                                             ,ll.var = curr.sigmasq.n)
                                            )

                ## Acceptance ratio
                ar <- acc.ra(log.prop = log.prop.posterior,
                                            log.current = log.curr.posterior)

                ## Move or stay
                ## stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                    if(i > burn.in) ar.na$fert.rate[j] <-
                        ar.na$fert.rate[j] + 1/n.iter
                } else {
                    ## if accept, update current fert rates, store proposed
                    ## rate, update current projection and count acceptance
                    if(runif(1) <= ar) {
                        if(i > burn.in) acc.count$fert.rate[j] <-
                            acc.count$fert.rate[j] + 1/n.iter
                        log.curr.f <- log.prop.f
                        log.curr.proj <- log.prop.proj
                        log.curr.posterior <- log.prop.posterior

                    }
                    ## if reject, leave current fert rates and projections
                    ## alone

                } # close else after checking for ar=0, missing, inf

            } # close else after checking neg or zero population

        } # close loop over all age-spec fertility rates

        ## Store proposed fertility rate matrix
        if(i > burn.in) fert.rate.mcmc[k,] <-
            as.vector(exp(log.curr.f[fert.rows,]))


        ## -------*** Survival

        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " Survival")

        ## Proposal (multivariate normal proposals; propose male and female simultaneously)

        ## cycle through components
        for(j in 1:length(logit.curr.s[["female"]])) {

            logit.fm.prop <- mvtnorm::rmvnorm(1, sigma = prop.varcovar[["surv.prop"]][[j]])
                                        # first is female

            ## make proposal
            logit.prop.s <- logit.curr.s
            for(l in 1:2) logit.prop.s[[l]][j] <- logit.prop.s[[l]][j] + logit.fm.prop[l]

            ## If proposal resulted in back-transformed s = 0 or 1, do
            ##   nothing
            logit.prop.s.unlist <- unlist(logit.prop.s, use.names = FALSE)
            if(invlogit(logit.prop.s.unlist[j]) > 1 - s.tol ||
               invlogit(logit.prop.s.unlist[j]) < s.tol) {
                ## leave current surv rates and projections
                ## alone (simply do not propose
                ## extreme survival probabilities)
                s.out.tol[j] <- s.out.tol[j] + 1/n.iter
            } else {

                ## Run CCMP (project on the original scale)
                ## ** Don't allow negative population; again, simply treat
                ##    this as if the proposal were never made
                full.proj <-
                    ccmp.function(pop = lapply(log.curr.b, "exp")
                                  ,fert = exp(log.curr.f)
                                  ,surv = lapply(logit.prop.s, "invlogit") #<-- use prop
                                  ,srb = exp(log.curr.srb)
                                  ,mig = curr.g
                                  ,proj.steps = proj.periods
                                  ,age.int = age.size)

                full.proj.unlist <- unlist(full.proj, use.names = FALSE)

                if(sum(full.proj.unlist < 0) > 0 || is.na(sum(full.proj.unlist))
                   || is.nan(sum(full.proj.unlist))) {
                    if(i > burn.in) {
                        pop.negative$surv.prop[j] <-
                            pop.negative$surv.prop[j] + 1/n.iter
                    }
                } else {
                    prop.proj <-
                        proj.cen.yrs(full.proj = full.proj
                                                    ,bline.yr = baseline.year
                                                    ,vr.yrs = vr.years
                                                    ,cen.yrs = census.years, proj.yrs = proj.years
                                                    )
                    log.prop.proj <- lapply(prop.proj, "log")

                    ## Calculate log posterior of proposed vital under projection
                    log.prop.posterior <-
                        log.post(f = log.curr.f
                                                ,s = logit.prop.s #<-- use proposal
                                                ,g = curr.g
                                                ,baseline.n = log.curr.b
                                                ,srb = log.curr.srb
                                                ,prior.mean.f = log.mean.f
                                                ,prior.mean.s = logit.mean.s
                                                ,prior.mean.g = mean.g
                                                ,prior.mean.b = log.mean.b
                                                ,prior.mean.srb = log.mean.srb
                                                ,alpha.f = al.f, beta.f = be.f
                                                ,alpha.s = al.s, beta.s = be.s
                                                ,alpha.g = al.g, beta.g = be.g
                                                ,alpha.n = al.n, beta.n = be.n
                                                ,alpha.srb = al.srb, beta.srb = be.srb
                                                ,sigmasq.f = curr.sigmasq.f
                                                ,sigmasq.s = curr.sigmasq.s
                                                ,sigmasq.g = curr.sigmasq.g
                                                ,sigmasq.n = curr.sigmasq.n
                                                ,sigmasq.srb = curr.sigmasq.srb
                                                ,non.zero.fert = fert.rows
                                                ,log.like = propRecon.log.lhood(
                                                 log.n.census = log.census.mat
                                                 ,log.n.hat = log.prop.proj #<-- use proposal
                                                 ,ll.var = curr.sigmasq.n)
                                                )

                    ## Acceptance ratio
                    ar <- acc.ra(log.prop = log.prop.posterior,
                                                log.current = log.curr.posterior)

                    ## Move or stay
                    ## stay if acceptance ratio 0, missing, infinity, etc.
                    if(is.na(ar) || is.nan(ar) || ar < 0) {
                        if(i > burn.in) ar.na$surv.prop[j] <-
                            ar.na$surv.prop[j] + 1/n.iter
                    } else {
                        ## if accept, update current surv rates, update current projection
                        ## and count acceptance
                        if(runif(1) <= ar) {
                            if(i > burn.in) acc.count$surv.prop[j] <-
                                acc.count$surv.prop[j] + 1/n.iter
                            logit.curr.s <- logit.prop.s
                            log.curr.proj <- log.prop.proj
                            log.curr.posterior <- log.prop.posterior
                        }
                        ## if reject, leave current surv rates and projections alone

                    } # close else{ after checking for undefined ar

                } # close else{ after checking for negative pop

            } # close else{ after checking for s outside tol

        } # close loop over all age-spec survival probabilities

        ## Store proposed survival probability matrix
        if(i > burn.in) for(l in 1:2)
            surv.prop.mcmc[[l]][k,] <- as.vector(invlogit(logit.curr.s[[l]]))


        ## -------*** Migration

        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " Migration")

        ## Proposal (multivariate normal proposals; propose male and female simultaneously)

        ## cycle through components
        for(j in 1:length(curr.g[["female"]])) {

            fm.prop <- mvtnorm::rmvnorm(1, sigma = prop.varcovar[["mig"]][[j]])
                                        # first is female

            ## make proposal
            prop.g <- curr.g
            for(l in 1:2) prop.g[[l]][j] <- prop.g[[l]][j] + fm.prop[l]

            ## Run CCMP (project on the original scale)
            ## ** Don't allow negative population; again, simply treat
            ##    this as if the proposal were never made
            full.proj <- ccmp.function(pop = lapply(log.curr.b, "exp")
                                       ,fert = exp(log.curr.f)
                                       ,surv = lapply(logit.curr.s, "invlogit")
                                       ,mig = prop.g #<-- use prop
                                       ,proj.steps = proj.periods
                                       ,age.int = age.size
                                       ,srb = exp(log.curr.srb)
                                       )

            full.proj.unlist <- unlist(full.proj, use.names = FALSE)

            if(sum(full.proj.unlist < 0) > 0 || is.na(sum(full.proj.unlist))
               || is.nan(sum(full.proj.unlist))) {
                if(i > burn.in) {
                    pop.negative$mig[j] <-
                        pop.negative$mig[j] + 1/n.iter
                }
            } else {
                prop.proj <-
                    proj.cen.yrs(full.proj = full.proj
                                                ,bline.yr = baseline.year
                                                ,vr.yrs = vr.years
                                                ,cen.yrs = census.years, proj.yrs = proj.years
                                                )
                log.prop.proj <- lapply(prop.proj, "log")

                ## Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
                    log.post(f = log.curr.f
                                            ,s = logit.curr.s
                                            ,g = prop.g #<-- use proposal
                                            ,baseline.n = log.curr.b
                                            ,srb = log.curr.srb
                                            ,prior.mean.f = log.mean.f
                                            ,prior.mean.s = logit.mean.s
                                            ,prior.mean.g = mean.g
                                            ,prior.mean.b = log.mean.b
                                            ,prior.mean.srb = log.mean.srb
                                            ,alpha.f = al.f, beta.f = be.f
                                            ,alpha.s = al.s, beta.s = be.s
                                            ,alpha.g = al.g, beta.g = be.g
                                            ,alpha.n = al.n, beta.n = be.n
                                            ,alpha.srb = al.srb, beta.srb = be.srb
                                            ,sigmasq.f = curr.sigmasq.f
                                            ,sigmasq.s = curr.sigmasq.s
                                            ,sigmasq.g = curr.sigmasq.g
                                            ,sigmasq.n = curr.sigmasq.n
                                            ,sigmasq.srb = curr.sigmasq.srb
                                            ,non.zero.fert = fert.rows
                                            ,log.like = propRecon.log.lhood(
                                             log.n.census = log.census.mat
                                             ,log.n.hat = log.prop.proj #<-- use proposal
                                             ,ll.var = curr.sigmasq.n)
                                            )

                ## Acceptance ratio
                ar <- acc.ra(log.prop = log.prop.posterior,
                                            log.current = log.curr.posterior)

                ## Move or stay
                ## stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                    if(i > burn.in) ar.na$mig[j] <-
                        ar.na$mig[j] + 1/n.iter
                } else {
                    ## if accept, update current surv rates,
                    ##  update current projection and count acceptance
                    if(runif(1) <= ar) {
                        if(i > burn.in) acc.count$mig[j] <-
                            acc.count$mig[j] + 1/n.iter
                        curr.g <- prop.g
                        log.curr.proj <- log.prop.proj
                        log.curr.posterior <- log.prop.posterior
                    }
                    ## if reject, leave current surv rates and projections alone

                } # close else{ after checking for undefined ar

            } # close else{ after checking for negative pop

        } # close loop over all age-spec migration

        ## Store proposed survival probability matrix
        if(i > burn.in) for(l in 1:2)
            mig.mcmc[[l]][k,] <- as.vector(curr.g[[l]])


        ## -------*** Baseline population

        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " Baseline")


        ## Proposal

        ## cycle through components
        for(j in 1:length(log.curr.b[["female"]])) {

            fm.prop <- mvtnorm::rmvnorm(1, sigma = prop.varcovar[["baseline.pop.count"]][[j]])
                                        # first is female

            ## make proposal
            log.prop.b <- log.curr.b
            for(l in 1:2) log.prop.b[[l]][j] <- log.prop.b[[l]][j] + fm.prop[l]

            ## Run CCMP (project on the original scale)
            ## ** Don't allow negative population
            full.proj <- ccmp.function(pop = lapply(log.prop.b, "exp"), #<-- use proposal
                                       fert = exp(log.curr.f),
                                       surv = lapply(logit.curr.s, "invlogit"),
                                       mig = curr.g,
                                       proj.steps = proj.periods,
                                       age.int = age.size,
                                       srb = exp(log.curr.srb))

            full.proj.unlist <- unlist(full.proj, use.names = FALSE)

            if(sum(full.proj.unlist < 0) > 0 || is.na(sum(full.proj.unlist))
               || is.nan(sum(full.proj.unlist))) {
                if(i > burn.in) {
                    pop.negative$baseline.count[j] <-
                        pop.negative$baseline.count[j] + 1/n.iter
                }
            } else {
                prop.proj <-
                    proj.cen.yrs(full.proj = full.proj
                                                ,bline.yr = baseline.year
                                                ,vr.yrs = vr.years
                                                ,cen.yrs = census.years
                                                , proj.yrs = proj.years
                                                )
                log.prop.proj <- lapply(prop.proj, "log")

                ## Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
                    log.post(f = log.curr.f
                                            ,s = logit.curr.s
                                            ,g = curr.g
                                            ,baseline.n = log.prop.b #<-- use proposal
                                            ,srb = log.curr.srb
                                            ,prior.mean.f = log.mean.f
                                            ,prior.mean.s = logit.mean.s
                                            ,prior.mean.g = mean.g
                                            ,prior.mean.b = log.mean.b
                                            ,prior.mean.srb = log.mean.srb
                                            ,alpha.f = al.f, beta.f = be.f
                                            ,alpha.s = al.s, beta.s = be.s
                                            ,alpha.g = al.g, beta.g = be.g
                                            ,alpha.n = al.n, beta.n = be.n
                                            ,alpha.srb = al.srb, beta.srb = be.srb
                                            ,sigmasq.f = curr.sigmasq.f
                                            ,sigmasq.s = curr.sigmasq.s
                                            ,sigmasq.g = curr.sigmasq.g
                                            ,sigmasq.n = curr.sigmasq.n
                                            ,sigmasq.srb = curr.sigmasq.srb
                                            ,non.zero.fert = fert.rows
                                            ,log.like = propRecon.log.lhood(
                                             log.n.census = log.census.mat
                                             ,log.n.hat = log.prop.proj #<-- use proposal
                                             ,ll.var = curr.sigmasq.n)
                                            )

                ## Acceptance ratio
                ar <- acc.ra(log.prop = log.prop.posterior,
                                            log.current = log.curr.posterior)

                ## Move or stay
                ## stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                    if(i > burn.in) ar.na$baseline.count[j] <-
                        ar.na$baseline.count[j] + 1/n.iter
                } else {
                    ## if accept, update current mig rates, store proposed
                    ## rate, update current projection and count acceptance
                    if(runif(1) <= ar) {
                        if(i > burn.in) acc.count$baseline.count[j] <-
                            acc.count$baseline.count[j] + 1/n.iter
                        log.curr.b <- log.prop.b
                        log.curr.proj <- log.prop.proj
                        log.curr.posterior <- log.prop.posterior
                    }
                    ## if reject, leave current fert rates and projections
                    ## alone, store current rate

                } # close else after checking for ar=na, nan, zero

            } # close else after checking for negative population

        } # close loop over all age-specific baseline counts

        ## Store proposed baseline count matrix
        if(i > burn.in) for(l in 1:2)
            baseline.count.mcmc[[l]][k,] <- as.vector(exp(log.curr.b[[l]]))


        ## -------*** SRB

        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " SRB")

        ## Proposal

        ## cycle through components
        for(j in 1:length(log.curr.srb)) {

            ## make proposal
            log.prop.srb <- log.curr.srb
            log.prop.srb[j] <-
                log.curr.srb[j] + rnorm(1, 0, sqrt(prop.varcovar$srb[[j]]))

            ## Run CCMP (project on the original scale)
            ## ** Don't allow negative population
            full.proj <- ccmp.function(pop = lapply(log.curr.b, "exp")
                                       ,fert = exp(log.curr.f)
                                       ,surv = lapply(logit.curr.s, "invlogit")
                                       ,mig = curr.g
                                       ,proj.steps = proj.periods
                                       ,age.int = age.size
                                       ,srb = exp(log.prop.srb) #<-- use proposal
                                       )

            full.proj.unlist <- unlist(full.proj, use.names = FALSE)

            if(sum(full.proj.unlist < 0) > 0 || is.na(sum(full.proj.unlist))
               || is.nan(sum(full.proj.unlist))) {
                if(i > burn.in) {
                    pop.negative$srb[j] <-
                        pop.negative$srb[j] + 1/n.iter
                }
            } else {
                prop.proj <-
                    proj.cen.yrs(full.proj = full.proj
                                                ,bline.yr = baseline.year
                                                ,vr.yrs = vr.years
                                                ,cen.yrs = census.years, proj.yrs = proj.years
                                                )
                log.prop.proj <- lapply(prop.proj, "log")

                ## Calculate log posterior of proposed vital under projection
                log.prop.posterior <-
                    log.post(f = log.curr.f
                                            ,s = logit.curr.s
                                            ,g = curr.g
                                            ,baseline.n = log.curr.b
                                            ,srb = log.prop.srb #<-- use prop
                                            ,prior.mean.f = log.mean.f
                                            ,prior.mean.s = logit.mean.s
                                            ,prior.mean.g = mean.g
                                            ,prior.mean.b = log.mean.b
                                            ,prior.mean.srb = log.mean.srb
                                            ,alpha.f = al.f, beta.f = be.f
                                            ,alpha.s = al.s, beta.s = be.s
                                            ,alpha.g = al.g, beta.g = be.g
                                            ,alpha.n = al.n, beta.n = be.n
                                            ,alpha.srb = al.srb, beta.srb = be.srb
                                            ,sigmasq.f = curr.sigmasq.f
                                            ,sigmasq.s = curr.sigmasq.s
                                            ,sigmasq.g = curr.sigmasq.g
                                            ,sigmasq.n = curr.sigmasq.n
                                            ,sigmasq.srb = curr.sigmasq.srb
                                            ,non.zero.fert = fert.rows
                                            ,log.like = propRecon.log.lhood(
                                             log.n.census = log.census.mat
                                             ,log.n.hat = log.prop.proj #<-- use proposal
                                             ,ll.var = curr.sigmasq.n)
                                            )

                ## Acceptance ratio
                ar <- acc.ra(log.prop = log.prop.posterior,
                                            log.current = log.curr.posterior)

                ## Move or stay
                ## stay if acceptance ratio 0, missing, infinity, etc.
                if(is.na(ar) || is.nan(ar) || ar < 0) {
                    if(i > burn.in) ar.na$srb[j] <-
                        ar.na$srb[j] + 1/n.iter
                } else {
                    ## if accept, update current mig rates, store proposed
                    ## rate, update current projection and count acceptance
                    if(runif(1) <= ar) {
                        if(i > burn.in) acc.count$srb[j] <-
                            acc.count$srb[j] + 1/n.iter
                        log.curr.srb <- log.prop.srb
                        log.curr.proj <- log.prop.proj
                        log.curr.posterior <- log.prop.posterior
                    }
                    ## if reject, leave current fert rates and projections
                    ## alone, store current rate

                } # close else after checking for ar=na, nan, zero

            } # close else after checking for negative population

        } # close loop over all age-specific baseline counts

        ## Store proposed baseline count matrix
        if(i > burn.in) srb.mcmc[k,] <- exp(log.curr.srb)


        ## -------** Variance Updates


        if(verb && identical(i%%progress.step, 0)) cat("\n", i, " Variances")

        ## For M-H Updates use asymmetric proposal:
        ## x' = x*exp(lambda * U - 0.5)
        ## where x' is the proposed move, x is the current value,
        ## U ~ Unif(0,1) and lambda is a tuning parameter to
        ## control the acceptance probabilities.
        ## The proposal density is then
        ## q(x'|x) = 1/x' => q(x|x') = 1/x != q(x'|x)


        ## -------*** Fertility rate

        if(v.prop.type == "inverse.gamma") {
            prop.sigmasq.f <-
                rinvGamma(1, al.f +
                                         length(mean.f[fert.rows,])/2,
                                         (be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                                          log.mean.f[fert.rows,])^2)) *
                                         prop.varcovar$variances["fert.rate.var"]
                                         )
        } else if (v.prop.type == "MH.scaled") {
            prop.sigmasq.f <-
                curr.sigmasq.f *
                    exp(prop.varcovar$variances["fert.rate.var"]*(runif(1,0,1)-0.5))
        }

        ## Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
            log.post(f = log.curr.f
                                    ,s = logit.curr.s
                                    ,g = curr.g
                                    ,baseline.n = log.curr.b
                                    ,srb = log.curr.srb
                                    ,prior.mean.f = log.mean.f
                                    ,prior.mean.s = logit.mean.s
                                    ,prior.mean.g = mean.g
                                    ,prior.mean.b = log.mean.b
                                    ,prior.mean.srb = log.mean.srb
                                    ,alpha.f = al.f, beta.f = be.f
                                    ,alpha.s = al.s, beta.s = be.s
                                    ,alpha.g = al.g, beta.g = be.g
                                    ,alpha.n = al.n, beta.n = be.n
                                    ,alpha.srb = al.srb, beta.srb = be.srb
                                    ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                                    ,sigmasq.s = curr.sigmasq.s
                                    ,sigmasq.g = curr.sigmasq.g
                                    ,sigmasq.n = curr.sigmasq.n
                                    ,sigmasq.srb = curr.sigmasq.srb
                                    ,non.zero.fert = fert.rows
                                    ,log.like = propRecon.log.lhood(
                                     log.n.census = log.census.mat
                                     ,log.n.hat = log.curr.proj
                                     ,ll.var = curr.sigmasq.n)
                                    )

        ## Acceptance ratio
        if(v.prop.type == "inverse.gamma") {
            ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                            ,log.curr.post = log.curr.posterior
                                            ,log.prop.var = dinvGamma(prop.sigmasq.f
                                             ,al.f + length(mean.f[fert.rows,])/2
                                             ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                                              log.mean.f[fert.rows,])^2)
                                             ,log = TRUE)
                                            ,log.curr.var = dinvGamma(curr.sigmasq.f
                                             ,al.f + length(mean.f[fert.rows,])/2
                                             ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                                              log.mean.f[fert.rows,])^2)
                                             ,log = TRUE)
                                            )
        } else if (v.prop.type == "MH.scaled") {
            ar <-
                acc.ra.var(log.prop.post = log.prop.posterior
                                          ,log.curr.post = log.curr.posterior
                                          ,log.prop.var = -log(prop.varcovar$variances["fert.rate.var"]
                                                               * prop.sigmasq.f)
                                          ,log.curr.var = -log(prop.varcovar$variances["fert.rate.var"]
                                                               * curr.sigmasq.f)
                                          )
        }

        ## Move or stay
        ## stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.f <-
                ar.na$sigmasq.f + 1/n.iter
        } else {
            ## if accept, update current, store proposed and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.f <-
                    acc.count$sigmasq.f + 1/n.iter
                curr.sigmasq.f <- prop.sigmasq.f
                log.curr.posterior <- log.prop.posterior
            }
            ## if reject, leave current and posterior

        } # close else after checking for ar=na, nan, zero

        if(i > burn.in) variances.mcmc[k,"fert.rate.var"] <- curr.sigmasq.f


        ## -------*** Survival Proportion

        if(v.prop.type == "inverse.gamma") {

            ## Constants
            hlf.l.mean <- length(unlist(mean.s, use.names = FALSE)) / 2
            hlf.sum.diff.sqrd <-
                0.5 * sum(unlist(mapply(FUN = "-"
                                        ,logit.curr.s, logit.mean.s
                                        ,SIMPLIFY = FALSE
                                        ), use.names = FALSE)^2)

            ## Proposal from inv gamma distn
            prop.sigmasq.s <-
                rinvGamma(1, al.s + hlf.l.mean,
                                         (be.s + hlf.sum.diff.sqrd)  *
                                         prop.varcovar$variances["surv.prop.var"]
                                         )

        } else if(v.prop.type == "MH.scaled") {

            prop.sigmasq.s <-
                curr.sigmasq.s *
                    exp(prop.varcovar$variances["surv.prop.var"]*(runif(1,0,1)-0.5))

        }

        ## Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
            log.post(f = log.curr.f
                                    ,s = logit.curr.s
                                    ,g = curr.g
                                    ,baseline.n = log.curr.b
                                    ,srb = log.curr.srb
                                    ,prior.mean.f = log.mean.f
                                    ,prior.mean.s = logit.mean.s
                                    ,prior.mean.g = mean.g
                                    ,prior.mean.b = log.mean.b
                                    ,prior.mean.srb = log.mean.srb
                                    ,alpha.f = al.f, beta.f = be.f
                                    ,alpha.s = al.s, beta.s = be.s
                                    ,alpha.g = al.g, beta.g = be.g
                                    ,alpha.n = al.n, beta.n = be.n
                                    ,alpha.srb = al.srb, beta.srb = be.srb
                                    ,sigmasq.f = curr.sigmasq.f
                                    ,sigmasq.s = prop.sigmasq.s #<-- use proposal
                                    ,sigmasq.g = curr.sigmasq.g
                                    ,sigmasq.n = curr.sigmasq.n
                                    ,sigmasq.srb = curr.sigmasq.srb
                                    ,non.zero.fert = fert.rows
                                    ,log.like = propRecon.log.lhood(
                                     log.n.census = log.census.mat
                                     ,log.n.hat = log.curr.proj
                                     ,ll.var = curr.sigmasq.n)
                                    )

        ## Acceptance ratio
        if(v.prop.type == "inverse.gamma") {
            ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                            ,log.curr.post = log.curr.posterior
                                            ,log.prop.var = dinvGamma(prop.sigmasq.s
                                             ,al.s + hlf.l.mean
                                             ,be.s + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            ,log.curr.var = dinvGamma(curr.sigmasq.s
                                             ,al.s + hlf.l.mean
                                             ,be.s + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            )
        } else if(v.prop.type == "MH.scaled") {
            ar <-
                acc.ra.var(log.prop.post = log.prop.posterior
                                          ,log.curr.post = log.curr.posterior
                                          ,log.prop.var = -log(prop.varcovar$variances["surv.prop.var"]
                                                               * prop.sigmasq.s)
                                          ,log.curr.var = -log(prop.varcovar$variances["surv.prop.var"]
                                                               * curr.sigmasq.s)
                                          )
        }

        ## Move or stay
        ## stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.s <-
                ar.na$sigmasq.s + 1/n.iter
        } else {
            ## if accept, update current, store proposed
            ## and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.s <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.s <- prop.sigmasq.s
                log.curr.posterior <- log.prop.posterior

            }
            ## if reject, leave current and posterior

        } # close else after checking for ar=na, nan, zero

        if(i > burn.in) variances.mcmc[k,"surv.prop.var"] <- curr.sigmasq.s


        ## -------*** Migration Proportion

        if(v.prop.type == "inverse.gamma") {

            ## Constants
            hlf.l.mean <- length(unlist(mean.g, use.names = FALSE)) / 2
            hlf.sum.diff.sqrd <-
                0.5 * sum(unlist(mapply(FUN = "-"
                                        ,curr.g, mean.g
                                        ,SIMPLIFY = FALSE
                                        ), use.names = FALSE)^2)

            ## Proposal
            prop.sigmasq.g <-
                rinvGamma(1, al.g + hlf.l.mean,
                                         (be.g + hlf.sum.diff.sqrd) * prop.varcovar$variances["mig.var"]
                                         )

        } else if(v.prop.type == "MH.scaled") {

            prop.sigmasq.g <-
                curr.sigmasq.g *
                    exp(prop.varcovar$variances["mig.var"]*(runif(1,0,1)-0.5))
        }

        ## Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
            log.post(f = log.curr.f
                                    ,s = logit.curr.s
                                    ,g = curr.g
                                    ,baseline.n = log.curr.b
                                    ,srb = log.curr.srb
                                    ,prior.mean.f = log.mean.f
                                    ,prior.mean.s = logit.mean.s
                                    ,prior.mean.g = mean.g
                                    ,prior.mean.b = log.mean.b
                                    ,prior.mean.srb = log.mean.srb
                                    ,alpha.f = al.f, beta.f = be.f
                                    ,alpha.s = al.s, beta.s = be.s
                                    ,alpha.g = al.g, beta.g = be.g
                                    ,alpha.n = al.n, beta.n = be.n
                                    ,alpha.srb = al.srb, beta.srb = be.srb
                                    ,sigmasq.f = curr.sigmasq.f
                                    ,sigmasq.s = curr.sigmasq.s
                                    ,sigmasq.g = prop.sigmasq.g #<-- use proposal
                                    ,sigmasq.n = curr.sigmasq.n
                                    ,sigmasq.srb = curr.sigmasq.srb
                                    ,non.zero.fert = fert.rows
                                    ,log.like = propRecon.log.lhood(
                                     log.n.census = log.census.mat
                                     ,log.n.hat = log.curr.proj
                                     ,ll.var = curr.sigmasq.n)
                                    )

        ## Acceptance ratio
        if(v.prop.type == "inverse.gamma") {

            ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                            ,log.curr.post = log.curr.posterior
                                            ,log.prop.var = dinvGamma(prop.sigmasq.g
                                             ,al.g + hlf.l.mean
                                             ,be.g + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            ,log.curr.var = dinvGamma(curr.sigmasq.g
                                             ,al.g + hlf.l.mean
                                             ,be.g + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            )

        } else if(v.prop.type == "MH.scaled") {

            ar <-
                acc.ra.var(log.prop.post = log.prop.posterior
                                          ,log.curr.post = log.curr.posterior
                                          ,log.prop.var = -log(prop.varcovar$variances["mig.var"]
                                                               * prop.sigmasq.g)
                                          ,log.curr.var = -log(prop.varcovar$variances["mig.var"]
                                                               * curr.sigmasq.g)
                                          )
        }

        ## Move or stay
        ## stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.g <-
                ar.na$sigmasq.g + 1/n.iter
        } else {
            ## if accept, update current, store proposed
            ## and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.g <-
                    acc.count$sigmasq.g + 1/n.iter
                curr.sigmasq.g <- prop.sigmasq.g
                log.curr.posterior <- log.prop.posterior

            }
            ## if reject, leave current and posterior

        } # close else after checking for ar=na, nan, zero

        if(i > burn.in) variances.mcmc[k,"mig.var"] <- curr.sigmasq.g


        ## -------*** Population Count

        if(v.prop.type == "inverse.gamma") {

            ## Constants
            hlf.l.mean <-
                (length(unlist(mean.b, use.names = FALSE)) +
                 length(unlist(log.census.mat, use.names = FALSE))
                 ) / 2
            hlf.sum.diff.sqrd <-
                0.5 * (sum(unlist(mapply(FUN = "-"
                                         ,log.curr.b, log.mean.b
                                         ,SIMPLIFY = FALSE
                                         ), use.names = FALSE)^2) +
                       sum(unlist(mapply(FUN = "-"
                                         ,log.census.mat, log.curr.proj
                                         ,SIMPLIFY = FALSE
                                         ), use.names = FALSE)^2)
                       )

            ## Proposal
            prop.sigmasq.n <-
                rinvGamma(1, al.n + hlf.l.mean / 2,
                                         (be.n + hlf.sum.diff.sqrd) * prop.varcovar$variances["population.count.var"]
                                         )

        } else if(v.prop.type == "MH.scaled") {

            prop.sigmasq.n <-
                curr.sigmasq.n *
                    exp(prop.varcovar$variances["population.count.var"]*(runif(1,0,1)-0.5))
        }

        ## Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
            log.post(f = log.curr.f
                                    ,s = logit.curr.s
                                    ,g = curr.g
                                    ,baseline.n = log.curr.b
                                    ,srb = log.curr.srb
                                    ,prior.mean.f = log.mean.f
                                    ,prior.mean.s = logit.mean.s
                                    ,prior.mean.g = mean.g
                                    ,prior.mean.b = log.mean.b
                                    ,prior.mean.srb = log.mean.srb
                                    ,alpha.f = al.f, beta.f = be.f
                                    ,alpha.s = al.s, beta.s = be.s
                                    ,alpha.g = al.g, beta.g = be.g
                                    ,alpha.n = al.n, beta.n = be.n
                                    ,alpha.srb = al.srb, beta.srb = be.srb
                                    ,sigmasq.f = curr.sigmasq.f
                                    ,sigmasq.s = curr.sigmasq.s
                                    ,sigmasq.g = curr.sigmasq.g
                                    ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                                    ,sigmasq.srb = curr.sigmasq.srb
                                    ,non.zero.fert = fert.rows
                                    ,log.like = propRecon.log.lhood(
                                     log.n.census = log.census.mat
                                     ,log.n.hat = log.curr.proj
                                     ,ll.var = prop.sigmasq.n)
                                    )

        ## Acceptance ratio
        if(v.prop.type == "inverse.gamma") {

            ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                            ,log.curr.post = log.curr.posterior
                                            ,log.prop.var = dinvGamma(prop.sigmasq.n
                                             ,al.n + hlf.l.mean
                                             ,be.n + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            ,log.curr.var = dinvGamma(curr.sigmasq.n
                                             ,al.n + hlf.l.mean
                                             ,be.n + hlf.sum.diff.sqrd
                                             ,log = TRUE)
                                            )

        } else if(v.prop.type == "MH.scaled") {

            ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                            ,log.curr.post = log.curr.posterior
                                            ,log.prop.var = -log(prop.varcovar$variances["population.count.var"]
                                                                 * prop.sigmasq.n)
                                            ,log.curr.var = -log(prop.varcovar$variances["population.count.var"]
                                                                 * curr.sigmasq.n)
                                            )
        }

        ## Move or stay
        ## stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.n <-
                ar.na$sigmasq.n + 1/n.iter
        } else {
            ## if accept, update current, store proposed
            ## and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.n <-
                    acc.count$sigmasq.n + 1/n.iter
                curr.sigmasq.n <- prop.sigmasq.n
                log.curr.posterior <- log.prop.posterior

            }
            ## if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

        if(i > burn.in) variances.mcmc[k,"population.count.var"] <- curr.sigmasq.n


        ## -------*** Sex-ratio-at-birth

        prop.sigmasq.srb <-
            curr.sigmasq.srb *
                exp(prop.varcovar$variances["srb.var"]*(runif(1,0,1)-0.5))

        ## Calculate log posterior
        log.prop.posterior <-
            log.post(f = log.curr.f
                                    ,s = logit.curr.s
                                    ,g = curr.g
                                    ,baseline.n = log.curr.b
                                    ,srb = log.curr.srb
                                    ,prior.mean.f = log.mean.f
                                    ,prior.mean.s = logit.mean.s
                                    ,prior.mean.g = mean.g
                                    ,prior.mean.b = log.mean.b
                                    ,prior.mean.srb = log.mean.srb
                                    ,alpha.f = al.f, beta.f = be.f
                                    ,alpha.s = al.s, beta.s = be.s
                                    ,alpha.g = al.g, beta.g = be.g
                                    ,alpha.n = al.n, beta.n = be.n
                                    ,alpha.srb = al.srb, beta.srb = be.srb
                                    ,sigmasq.f = curr.sigmasq.f
                                    ,sigmasq.s = curr.sigmasq.s
                                    ,sigmasq.g = curr.sigmasq.g
                                    ,sigmasq.n = curr.sigmasq.n
                                    ,sigmasq.srb = prop.sigmasq.srb #<-- use proposal
                                    ,non.zero.fert = fert.rows
                                    ,log.like = propRecon.log.lhood(
                                     log.n.census = log.census.mat
                                     ,log.n.hat = log.curr.proj
                                     ,ll.var = curr.sigmasq.n)
                                    )

        ## Acceptance ratio
        ar <- acc.ra.var(log.prop.post = log.prop.posterior
                                        ,log.curr.post = log.curr.posterior
                                        ,log.prop.var = -log(prop.varcovar$variances["srb.var"]
                                                             * prop.sigmasq.srb)
                                        ,log.curr.var = -log(prop.varcovar$variances["srb.var"]
                                                             * curr.sigmasq.srb)
                                        )

        ## Move or stay
        ## stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.srb <-
                ar.na$sigmasq.srb + 1/n.iter
        } else {
            ## if accept, update current, store proposed
            ## and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.srb <-
                    acc.count$sigmasq.srb + 1/n.iter
                curr.sigmasq.srb <- prop.sigmasq.srb
                log.curr.posterior <- log.prop.posterior
            }
        } # close else after checking for ar=na, nan, zero

        ## Store proposed srb
        if(i > burn.in) variances.mcmc[k,"srb.var"] <- curr.sigmasq.srb

        if(verb && identical(i%%progress.step, 0)) cat("\n\n")


        ## -------** Store current population

        if(i > burn.in) {
            full.curr.proj <-
                ccmp.function(pop = lapply(log.curr.b, "exp"),
                              surv = lapply(logit.curr.s, "invlogit"),
                              fert = exp(log.curr.f),
                              mig = curr.g,
                              proj.steps = proj.periods,
                              age.int = age.size
                              ,label.dims = TRUE
                              ,base.year = "1985"
                              )

            for(l in 1:2) lx.mcmc[[l]][k,] <- full.curr.proj[[l]][,-1]
        }


        ## -------** Save mid-run

        if(!is.null(save.mid.run.every) && !is.null(save.mid.run.name) &&
           i > burn.in) {
            if(identical(k %% save.mid.run.every, 0)) {
                mid.save.list <-
                    list(fert.rate.mcmc = fert.rate.mcmc[1:k,]
                         ,surv.prop.mcmc = lapply(surv.prop.mcmc
                          ,function(z) z[1:k,])
                         ,mig.XXX.mcmc = lapply(mig.mcmc
                          ,function(z) z[1:k,])
                         ,baseline.count.mcmc = lapply(baseline.count.mcmc
                          ,function(z) z[1:k,])
                         ,srb.mcmc = srb.mcmc[1:k,]
                         ,lx.mcmc = lapply(lx.mcmc
                          ,function(z) z[1:k,])
                         ,variances.mcmc = variances.mcmc[1:k,]
                         ,alg.params = list(non.zero.fert.rows = fert.rows
                          ,iters = n.iter, burn.in = burn.in
                         )
                         )
                names(mid.save.list) <- sub("XXX", mig.string, names(mid.save.list))
                fn <- paste(save.mid.run.name, "_OUTPUT_mid_run.Rdata", sep = "")
                if(verb) message("saving 'mid.save.list' to file ", fn)
                save(mid.save.list, file = fn)
                rm(mid.save.list); gc()
            }
        }

    } # Ends outer-most loop

    ## ......... End Loop ........ ##
    ##.............................##


    ## -------* OUTPUT

    ## Initial values
    start.vals <- list(fert.rate = start.f
                       ,surv.prop = start.s
                       ,mig.XXX = start.g
                       ,baseline.count = start.b
                       ,srb = start.srb
                       ,start.sigmasq.f = start.sigmasq.f
                       ,start.sigmasq.s = start.sigmasq.s
                       ,start.sigmasq.g = start.sigmasq.g
                       ,start.sigmasq.n = start.sigmasq.n
                       ,start.sigmasq.srb = start.sigmasq.srb
                       ,pop.data = pop.data
                       )

    ## fixed parameters
    fixed.params <- list(alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.mig.XXX = al.g
                         ,beta.mig.XXX = be.g
                         ,alpha.population.count = al.n
                         ,beta.population.count = be.n
                         ,alpha.srb = al.srb
                         ,beta.srb = be.srb
                         ,mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.mig.XXX = mean.g
                         ,mean.baseline.count = mean.b
                         ,mean.srb = mean.srb
                         ,mean.pop.data = pop.data
                         )

    ## migration type
    new.names <-
        lapply(list(start.vals, fixed.params), FUN = function(z) {
            sub("XXX", mig.string, names(z))
        })
    names(start.vals) <- new.names[[1]]
    names(fixed.params) <- new.names[[2]]

    ## algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    ## algorithm parameters
    alg.params <- list(prop.varcovar = prop.varcovar
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit", mig.XXX = "I"
                        ,baseline.count = "log"
                        ,population.count = "log"
                        ,srb = "log")
                       ,projection.periods = proj.periods
                       ,age.gp.size = age.size
                       ,non.zero.fert.rows = fert.rows
                       ,surv.tolerance = s.tol
                       ,burn.in = burn.in
                       ,iters = n.iter
                       ,years = list(vital.rate.years = vr.years
                        ,baseline.year = baseline.year
                        ,census.years = census.years
                        ,projection.years = proj.years
                        )
                       )
    names(alg.params) <- sub("XXX", mig.string, names(alg.params))

    ## results
    if(filebacked.chains) {
        fert.rate.mcmc <-
            ff.mcmc(fert.rate.mcmc, start = 1, end = n.iter
                    ,thin = 1)
        surv.prop.mcmc <-
            ff.mcmc(surv.prop.mcmc, start = 1, end = n.iter
                    ,thin = 1)
        lx.mcmc <-
            ff.mcmc(lx.mcmc, start = 1, end = n.iter
                    ,thin = 1)
        mig.mcmc <-
            ff.mcmc(mig.mcmc, start = 1, end = n.iter, thin = 1)
        baseline.count.mcmc <-
            ff.mcmc(baseline.count.mcmc, start = 1, end = n.iter
                    ,thin = 1)
        variances.mcmc <-
            ff.mcmc(variances.mcmc, start = 1, end = n.iter
                    ,thin = 1)
        srb.mcmc <-
            ff.mcmc(srb.mcmc, start = 1, end = n.iter
                    ,thin = 1)
    }

    ret.list <- list(fert.rate.mcmc = fert.rate.mcmc
                     ,surv.prop.mcmc = surv.prop.mcmc
                     ,mig.XXX.mcmc = mig.mcmc
                     ,baseline.count.mcmc = baseline.count.mcmc
                     ,lx.mcmc = lx.mcmc
                     ,variances.mcmc = variances.mcmc
                     ,srb.mcmc = srb.mcmc
                     ,alg.stats = alg.stats
                     ,fixed.params = fixed.params
                     ,start.vals = start.vals
                     ,alg.params = alg.params
                     )

    names(ret.list) <- sub("XXX", mig.string, names(ret.list))

    return(ret.list)

}
