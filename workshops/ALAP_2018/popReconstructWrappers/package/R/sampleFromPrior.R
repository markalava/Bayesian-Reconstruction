
##' Draw a random sample from the joint prior distribution
##'
##' Sample from the prior to determine distributions of transformed
##' parameters such as mortality rate and migration counts.
##'
##' @param n.iter Size of sample to draw.
##' @param batch.size Size of batches in which to make draws.
##' @param al.f,be.f,al.s,be.s,al.g,be.g,al.n,be.n,al.srb,be.srb,mean.f,mean.s,mean.g,mean.b,mean.srb,start.f,start.s,start.g,start.b,start.srb,start.sigmasq.f,start.sigmasq.s,start.sigmasq.g,start.sigmasq.n,start.sigmasq.srb,fert.rows,proj.periods,age.size See \code{\link{pop.est.sampler}}.
##' @param message.every Control progress reporting; print a message every \code{message.every} iterations.
##' @param parallelize Run in parallel?
##' @param cores Number of parallel processes (cores) to use if \code{parallelize} is true.
##' @param pvm.hostfile pvm hosts for parallelization via SNOW on *NIX.
##' @param pos.sample.max.tries
##' @param max.elapsed.time
##' @param pvm.hostfile
##' @return
##' @author Mark C. Wheldon
##' @export
##' @import ff
sample.from.prior <- function(n.iter = 1E3,
                              batch.size = ceiling(n.iter/4),
                              al.f, be.f, al.s, be.s, al.g, be.g, al.n, be.n, al.srb, be.srb,
                              mean.f, mean.s, mean.g, mean.b, mean.srb,
                              fert.rows = as.logical(apply(mean.f == 0L, 1, function(z) !all(z))),
                              proj.periods = ncol(mean.f),
                              name.pref = "", name.suff = "",
                              age.size = 5,
                              message.every = 1000,
                              neg.pop.tol = 0,
                              parallelize = TRUE,
                              cores = NULL,
                              pvm.hostfile = "~/.pvm_hosts"
                              ) {

    message("\nSAMPLING FROM PRIOR")

    ##-------* Inputs

    baseline <- mean.b
    fert <- mean.f
    surv <- mean.s
    mig <- mean.g
    srb <- mean.srb
    hyper.params <- list(alpha.fert.rate = al.f, beta.fert.rate = be.f,
                         alpha.surv.prop = al.s, beta.surv.prop = be.s,
                         alpha.mig.prop = al.g, beta.mig.prop = be.g,
                         alpha.population.count = al.n, beta.population.count = be.n,
                         alpha.srb = al.srb, beta.srb = be.srb)
    age.int <- age.size
    proj.steps <- proj.periods
    base.year <- colnames(baseline)[1]
    ccmp.f <- ccmp.femDom

    label.dims <- TRUE

    pos.sample.max.tries <- 1E6
    max.elapsed.time <- 3600*24*2 # 2 days

    name.pref <- ""
    name.suf <- ""
    output.dir <- "."
    return.list <- TRUE

    platform.type <- .Platform$OS.type

    ## -------* Do input checks

    ## Batchsize should be <= n.iter
    if(batch.size > n.iter) {
        batch.size <- n.iter
        message("batch.size > n.iter; batch.size now ", batch.size)
    }

    ## CCMPP function
    ccmp.f <- match.fun(ccmp.f)

    ## Basic check of female/male lists

    ## Female and male components must have same dim
    check.FM.dims <- function(obj) {
        if(!is.null(dim(obj[["female"]])) && !is.null(dim(obj[["male"]]))) {
            if(!identical(dim(obj[["female"]]), dim(obj[["male"]])))
                stop("Dimensions of female and male components of "
                    ,deparse(substitute(obj), " not equal"))
            if(!identical(dimnames(obj[["female"]]), dimnames(obj[["male"]])))
                warning("Dimnames of female and male components of "
                       ,deparse(substitute(obj)
                               ," not equal. Dimnames of female component will be used."
                                ))
        } else {
            if(!identical(length(obj[["female"]]), length(obj[["male"]])))
                stop("Lengths of female and male components of "
                    ,deparse(substitute(obj), " not equal"))
            if(!identical(names(obj[["female"]]), names(obj[["male"]])))
                warning("Names of female and male components of "
                       ,deparse(substitute(obj)
                               ," not equal. Names of female component will be used."
                                ))
        }
    }
    check.FM.dims(surv)
    check.FM.dims(mig)
    check.FM.dims(baseline)

    ## All inputs must be finite
    all.finite <- function(x) {
        fcall <- as.list(match.call())
        if(sum(!is.finite(x)) > 0)
            stop("'", deparse(fcall$x), "' has non-finite values")
    }
    all.finite(fert)
    lapply(surv, "all.finite")
    lapply(mig, "all.finite")
    lapply(baseline, "all.finite")


    ## -------* Set up parallel backends

    if(parallelize) {

        cores <- min(parallel::detectCores() - 1, cores)

            if(requireNamespace("doParallel", quietly = TRUE)) {

                ## Parllel setup
                w <- parallel::makeCluster(cores)
                doParallel::registerDoParallel(w)

                message("Parallel Workers: ", foreach::getDoParWorkers())

                on.exit(parallel::stopCluster(w))

            } else {
                stop("You must install 'doParallel' package.")
            }

        }


    ## -------* Functions

    ## -------** Basic

    `%dopar%` <- foreach::`%dopar%`

    logit <- function(p) log(p / (1-p))
    invlogit <- function(x) exp(x) / (1 + exp(x))

    ## make a matrix with same dims and names as this one
    matrix.along <- function(x, m, byrow = FALSE) {
        matrix(x, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m)
              ,byrow = byrow
               )
    }

    ## make a file name from a string possibly containing illegal and /or
    ## undesirable characters
    make.filename <-
        function(s
                ,forbidden= "\\||\\\\|\\?|\\*|<|:|>|\\[|\\]|/|\\."
                ,replacement = rep("_", length(forbidden))
                ,fixed = FALSE
                 )
        {
            ## AGRUGMENTS
            ## s:            a character string to make safe
            ## forb:    character vector of forbidden characters
            ## replacement: character vector same length as 'forb'
            ##               which gives replacement to be used
            ## fixed:        logical. Indicates whether forbidden is a string
            ##               to be matched as is. Otherwise it is a regexp.

            for(i in 1:length(forbidden)) {
                s <- gsub(pattern = forbidden[i]
                         ,x = s
                         ,replacement = replacement[i]
                         ,fixed = fixed
                          )
            }
            s
        }


    ## -------** Location scale t functions

    ## Random sample
    loc.scale.t.rs <- function(mu = 0, alpha, beta, n = 1) {
        rt(n = n, df = 2 * alpha) * sqrt(beta/alpha) + mu
    }

    ## Quantile
    loc.scale.t.q <- function(p, mu = 0, alpha, beta) {
        mu + qt(p = p, df = 2 * alpha) * sqrt(beta/alpha)
    }

    ## Density
    ## loc.scale.t.d <- function(x, mu = 0, alpha, beta) {
    ##     s <- sqrt(beta/alpha)
    ##     dt((x-mu)/s, df = 2*alpha) / s
    ## }

    ## Truncated, quantile
    trunc.loc.scale.t.q <- function(p, mu = 0, alpha, beta, llim = -1) {
        nu <- 2 * alpha
        sc <- sqrt(beta/alpha)
        k <- pt(q = (llim - mu)/sc, df = nu)
        return(qt(p = p * (1 - k) + k, df = nu) * sc + mu)
    }

    ## Truncated, random sample (inverse CDF method)
    trunc.loc.scale.t.rs <- function(n = 1, mu = 0, alpha, beta, llim = -1) {
        u <- runif(n = n, 0,1)
        return(trunc.loc.scale.t.q(u, mu = mu
                                  ,alpha = alpha, beta = beta, llim = llim
                                   ))
    }


    ## -------** Sampling functions

    ## -------*** Sample vital rates

    ## loc.scale.t not resulting in NaN
    sampleVR <- function(init.est, a, b) {
        out <- NaN
        out <-
            loc.scale.t.rs(mu = init.est, alpha = a, beta = b
                          ,n = length(init.est)
                           )
        out
    }

    ## trunc.loc.scale.t not resultg in in NaN
    sampleMig <- function(init.est, a, b, llim = -1) {
        out <- NaN
        while(any(!is.finite(out))) {
            out <-
                trunc.loc.scale.t.rs(mu = init.est, alpha = a, beta = b
                                    ,n = length(init.est), llim = llim
                                     )
        }
        out
    }

    ## Sample non-migration vital rates
    drawVRSample <- function() {
        out <- list(fert.sample = NaN
                   ,surv.sample = NaN
                   ,mig.sample = NaN
                   ,baseline.sample = NaN
                   ,srb.sample = NaN
                    )
        while(any(unlist(rapply(out, function(z) any(!is.finite(z))
                               ,how = "replace"
                                ), recursive = TRUE))) {
                                    out <- list(fert.sample =
                                                    exp(sampleVR(log(fert[fert.rows,])
                                                                ,hyper.params$alpha.fert.rate
                                                                ,hyper.params$beta.fert.rate
                                                                 ))
                                               ,surv.sample = lapply(surv, function(z) {
                                                    invlogit(sampleVR(logit(z)
                                                                     ,hyper.params$alpha.surv.prop
                                                                     ,hyper.params$beta.surv.prop
                                                                      ))
                                                })
                                               ,mig.sample = lapply(mig, function(z) {
                                                    sampleMig(z
                                                             ,hyper.params$alpha.mig.prop
                                                             ,hyper.params$beta.mig.prop
                                                              )
                                                })
                                               ,baseline.sample = lapply(baseline, function(z) {
                                                    exp(sampleVR(log(z)
                                                                ,hyper.params$alpha.population.count
                                                                ,hyper.params$beta.population.count
                                                                 ))
                                                })
                                               ,srb.sample =
                                                    exp(sampleVR(log(srb)
                                                                ,hyper.params$alpha.srb
                                                                ,hyper.params$beta.srb
                                                                 ))
                                                )
                                }
        out
    }


    ## -------*** Project sample

    projectSample <- function(input.list) {
        fertMatrix <-
            matrix(0, nrow = nrow(fert), ncol = ncol(fert))
        fertMatrix[fert.rows,] <- input.list$fert.sample
        ccmp.f(fert = fertMatrix
              ,surv = input.list$surv.sample
              ,mig = input.list$mig.sample
              ,pop = input.list$baseline.sample
              ,srb = input.list$srb.sample
              ,base.year = colnames(input.list$baseline.sample[[1]])
              ,age.int = age.int
               )
    }


    ## -------*** Draw sample resulting in positive population

    positivePopSample <- function() {
        projPop <- numeric()
        neg.pop <- neg.pop.tol + 1

        no.tries <- 0

        while(neg.pop > neg.pop.tol && no.tries < pos.sample.max.tries) {
            VRSample <- drawVRSample()
            projPop <- projectSample(VRSample)
            neg.pop <- sum(unlist(lapply(projPop, function(z) sum(z < 0))))
            if(!is.finite(neg.pop)) neg.pop <- neg.pop.tol + 1
            no.tries <- no.tries + 1
        }

        if(neg.pop > neg.pop.tol)
            stop("No non-negative pop counts found in under 'pos.sample.max.tries'")

        return(c(VRSample, list(pop.proj = projPop)))
    }


    ## -------** Storage functions

    ## -------*** Copy batch to ff

    ## Copy batch output to the ff object
    cpBatch2ff <- function(batch) {
        start.at <- batch.size * (batches.done) + 1 # save is done before
                                # incrementing
                                # 'batches.done'
        for(l in 1:min(n.iter - start.at + 1, length(batch))) {
            fert.rate.prior[start.at + l - 1,] <-
                as.vector(batch[[l]]$fert.sample)
            srb.prior[start.at + l - 1,] <-
                as.vector(batch[[l]]$srb.sample)
            for(s in c("female", "male")) {
                surv.prop.prior[[s]][start.at + l - 1,] <-
                    as.vector(batch[[l]]$surv.sample[[s]])
                mig.prop.prior[[s]][start.at + l - 1,] <-
                    as.vector(batch[[l]]$mig.sample[[s]])
                baseline.count.prior[[s]][start.at + l - 1,] <-
                    as.vector(batch[[l]]$baseline.sample[[s]])
                pop.count.prior[[s]][start.at + l - 1,] <-
                    as.vector(batch[[l]]$pop.proj[[s]][,-1])
            }
        }
    }


    ## -------*** Save ff objects

    saveff <- function(end.point, filename.tag, save.dir) {

        ## Make sure don't try to save beyond end of array
        end.point <- min(n.iter, end.point)

        ## Change to save.dir for duration of function
        oldwd <- getwd()
        on.exit(setwd(oldwd))

        if(!missing(save.dir)) setwd(file.path(getwd(), save.dir))

        ## Fert
        assign(paste(name.pref, "fert.rate.prior", name.suf, sep = "")
              ,fert.rate.prior[1:end.point,]
               )
        save(list = paste(name.pref, "fert.rate.prior", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "fert_rate_prior"
                                             ,filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## SRB
        assign(paste(name.pref, "srb.prior", name.suf, sep = "")
              ,srb.prior[1:end.point,])
        save(list = paste(name.pref, "srb.prior", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "srb_prior"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Surv female
        assign(paste(name.pref, "surv.prop.prior.female", name.suf, sep = "")
              ,surv.prop.prior[["female"]][1:end.point,]
               )
        save(list = paste(name.pref, "surv.prop.prior.female", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "surv_prop_prior_female"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Mig female
        assign(paste(name.pref, "mig.prop.prior.female", name.suf, sep = "")
              ,mig.prop.prior[["female"]][1:end.point,]
               )
        save(list = paste(name.pref, "mig.prop.prior.female", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "mig_prop_prior_female"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Baseline female
        assign(paste(name.pref, "baseline.count.prior.female", name.suf, sep = "")
              ,baseline.count.prior[["female"]][1:end.point,]
               )
        save(list = paste(name.pref, "baseline.count.prior.female", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "baseline_count_prior_female"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Pop count female
        assign(paste(name.pref, "pop.count.prior.female", name.suf, sep = "")
              ,pop.count.prior[["female"]][1:end.point,]
               )
        save(list = paste(name.pref, "pop.count.prior.female", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "pop_count_prior_female"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Surv male
        assign(paste(name.pref, "surv.prop.prior.male", name.suf, sep = "")
              ,surv.prop.prior[["male"]][1:end.point,]
               )
        save(list = paste(name.pref, "surv.prop.prior.male", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "surv_prop_prior_male"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Mig male
        assign(paste(name.pref, "mig.prop.prior.male", name.suf, sep = "")
              ,mig.prop.prior[["male"]][1:end.point,]
               )
        save(list = paste(name.pref, "mig.prop.prior.male", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "mig_prop_prior_male"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Baseline male
        assign(paste(name.pref, "baseline.count.prior.male", name.suf, sep = "")
              ,baseline.count.prior[["male"]][1:end.point,]
               )
        save(list = paste(name.pref, "baseline.count.prior.male", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "baseline_count_prior_male"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )

        ## Pop count male
        assign(paste(name.pref, "pop.count.prior.male", name.suf, sep = "")
              ,pop.count.prior[["male"]][1:end.point,]
               )
        save(list = paste(name.pref, "pop.count.prior.male", name.suf, sep = "")
            ,file = paste(make.filename(paste(name.pref, "pop_count_prior_male"
                                            , filename.tag, name.suf, sep = ""))
                         ,".RData", sep = "")
             )
    }


    ## -------* Constants

    ## -------** 'Measure' inputs

    fert.years <- unlist(strsplit(colnames(fert), "[^0-9]"))
    surv.years <- unlist(strsplit(colnames(surv[["female"]]), "[^0-9]"))
    mig.years <- unlist(strsplit(colnames(mig[["female"]]), "[^0-9]"))
    baseline.years <- unlist(strsplit(colnames(baseline[["female"]]), "[^0-9]"))
    srb.years <- unlist(strsplit(colnames(srb), "[^0-9]"))

    fert.ages <- unlist(strsplit(rownames(fert), "[^0-9]"))
    surv.ages <- unlist(strsplit(rownames(surv[["female"]]), "[^0-9]"))
    mig.ages <- unlist(strsplit(rownames(mig[["female"]]), "[^0-9]"))
    baseline.ages <-
        sapply(strsplit(rownames(baseline[["female"]]), "[^0-9]"), "[[", 1)

    len.mig.fem <- length(mig[["female"]])

    ## baseline dims
    if(length(baseline.years) != 1) stop("baseline must be a single column")

    ## srb dims
    if(dim(srb)[1] !=1) stop("srb must be a single row")

    ## years had better be the same
    if(!all(c(identical(fert.years, surv.years), identical(fert.years, mig.years)
             ,baseline.years %in% fert.years, srb.years %in% fert.years))
       )
        stop("colnames of inputs (years) should be consistent")
    VR.years <- fert.years

    ## ages had better be the same
    if(!all(c(fert.ages %in% surv.ages, mig.ages %in% surv.ages
             ,baseline.ages %in% surv.ages)
            ))
        stop("rownames of inputs (ages) should be consistent")


    ## -------** Projection years

    proj.years <-
        seq(from = as.numeric(baseline.years) + age.int, by = age.int
           ,length = ncol(fert))

    message("Projection years: ", paste(proj.years, collapse = ", "))


    ## -------* Storage for output

    ## Use 'ff' objects to prevent running out of memory

    ## Fert rates
    fert.rate.prior <-
        ff(0, dim = c(n.iter, nrow(fert[fert.rows,]) * ncol(fert)))
    colnames(fert.rate.prior) <-
        paste(rep(fert.years, each = length(fert.ages[fert.rows]))
             ,rep(fert.ages[fert.rows], length(fert.years))
             ,sep = ".")

    ## Surv props
    surv.prop.prior <-
        list(female = ff(0
                       , dim = c(n.iter, nrow(surv[["female"]]) *
                                         ncol(surv[["female"]])))
            ,male = ff(0
                     , dim = c(n.iter, nrow(surv[["male"]]) * ncol(surv[["male"]])))
             )
    for(i in 1:2) colnames(surv.prop.prior[[i]]) <-
                      paste(rep(surv.years, each = length(surv.ages))
                           ,rep(surv.ages, length(surv.years))
                           ,sep = ".")

    ## Mig props
    mig.prop.prior <-
        list(female = ff(0
                       , dim = c(n.iter, nrow(mig[["female"]]) * ncol(mig[["female"]])))
            ,male = ff(0
                     , dim = c(n.iter, nrow(mig[["male"]]) * ncol(mig[["male"]])))
             )
    for(i in 1:2) colnames(mig.prop.prior[[i]]) <-
                      paste(rep(mig.years, each = length(mig.ages))
                           ,rep(mig.ages, length(mig.years))
                           ,sep = ".")

    ## Baseline counts
    baseline.count.prior <-
        list(female = ff(0, dim = c(n.iter
                                   ,nrow(baseline[["female"]]) * ncol(baseline[["female"]])))
            ,male = ff(0, dim = c(n.iter
                                , nrow(baseline[["male"]]) * ncol(baseline[["male"]])))
             )
    for(i in 1:2) colnames(baseline.count.prior[[i]]) <-
                      paste(rep(baseline.years, each = length(baseline.ages))
                           ,rep(baseline.ages, length(baseline.years))
                           ,sep = ".")

    ## SRB
    srb.prior <- ff(0, dim = c(n.iter, nrow(srb) * ncol(srb)))
    colnames(srb.prior) <- srb.years

    ## Projected population counts
    pop.count.prior <-
        list(female = ff(0
                       , dim = c(n.iter, nrow(mig[["female"]]) * ncol(mig[["female"]])))
            ,male = ff(0
                     , dim = c(n.iter, nrow(mig[["male"]]) * ncol(mig[["male"]])))
             )
    for(i in 1:2) colnames(pop.count.prior[[i]]) <-
                      paste(rep(proj.years, each = length(mig.ages))
                           ,rep(mig.ages, length(proj.years))
                           ,sep = ".")


    ## -------* PERFORM SAMPLING

    ## Control counters
    size.achieved <- 0
    elap.time <- 0
    batches.done <- 0

    ## Run 'batch.size' parallel samplers until they find a non-negative pop count
    while(elap.time < max.elapsed.time && size.achieved < n.iter) {

        bsNA <- 1

        ## Guard against any 'NA' creeping in due to faults with nodes
        while(bsNA > 0) {

            st <- system.time({
                batchSample <- foreach::foreach(iterators::icount(batch.size)) %dopar% positivePopSample()
            })

            bsNA <- sum(sapply(batchSample, function(z) {
                unlist(sapply(z, function(y) {
                    if(is.recursive(y)) sapply(y, "is.na")
                    else is.na(y)
                }))
            }))
            if(bsNA > 0) message(" --- NAs detected")
        }

        ## Copy batch output to the ff object
        cpBatch2ff(batchSample)

        ## Save temporary copy of ff objects
        saveff(end.point = (batches.done + 1) * batch.size # save is done before
                                # incrementing
                                # 'batches.done'
              ,filename.tag = "_TEMP"
              ,save.dir = output.dir
               )

        ## Update counters
        batches.done <- batches.done + 1
        size.achieved <- size.achieved + batch.size
        elap.time <- elap.time + st[3]

        if(batches.done %% message.every == 0)
            message("\nelapsed time = ", elap.time, "\nbatches.done = "
                   ,batches.done, "\nsize.achieved = "
                   ,size.achieved)

    }


    ## -------** How many achieved?

    if(size.achieved < n.iter) {

        warning("size.achieved < n.iter: size.achieved = ", size.achieved)


        ## -------*** Trim ff objects

        new.clone <- function(...) suppressWarnings(clone(...))

        ## Fert
        fert.rate.prior <-
            vt(new.clone(vt(fert.rate.prior)
                        ,length = end.point * ncol(fert.rate.prior)
                        ,dim = c(ncol(fert.rate.prior), end.point)
                         ))
        colnames(fert.rate.prior) <-
            paste(rep(fert.years, each = length(fert.ages[fert.rows]))
                 ,rep(fert.ages[fert.rows], length(fert.years))
                 ,sep = ".")

        ## Surv props
        surv.prop.prior <-
            list(female = vt(new.clone(vt(surv.prop.prior[["female"]])
                                      ,length = end.point * ncol(surv.prop.prior[["female"]])
                                      ,dim = c(ncol(surv.prop.prior[["female"]]), end.point)
                                       ))
                ,male = vt(new.clone(vt(surv.prop.prior[["male"]])
                                    ,length = end.point * ncol(surv.prop.prior[["male"]])
                                    ,dim = c(ncol(surv.prop.prior[["male"]]), end.point)
                                     ))
                 )
        for(i in 1:2) colnames(surv.prop.prior[[i]]) <-
                          paste(rep(surv.years, each = length(surv.ages))
                               ,rep(surv.ages, length(surv.years))
                               ,sep = ".")

        ## Mig props
        mig.prop.prior <-
            list(female = vt(new.clone(vt(mig.prop.prior[["female"]])
                                      ,length = end.point * ncol(mig.prop.prior[["female"]])
                                      ,dim = c(ncol(mig.prop.prior[["female"]]), end.point)
                                       ))
                ,male = vt(new.clone(vt(mig.prop.prior[["male"]])
                                    ,length = end.point * ncol(mig.prop.prior[["male"]])
                                    ,dim = c(ncol(mig.prop.prior[["male"]]), end.point)
                                     ))
                 )
        for(i in 1:2) colnames(mig.prop.prior[[i]]) <-
                          paste(rep(mig.years, each = length(mig.ages))
                               ,rep(mig.ages, length(mig.years))
                               ,sep = ".")

        ## Baseline counts
        baseline.count.prior <-
            list(female = vt(new.clone(vt(baseline.count.prior[["female"]])
                                      ,length = end.point * ncol(baseline.count.prior[["female"]])
                                      ,dim = c(ncol(baseline.count.prior[["female"]]), end.point)
                                       ))
                ,male = vt(new.clone(vt(baseline.count.prior[["male"]])
                                    ,length = end.point * ncol(baseline.count.prior[["male"]])
                                    ,dim = c(ncol(baseline.count.prior[["male"]]), end.point)
                                     ))
                 )
        for(i in 1:2) colnames(baseline.count.prior[[i]]) <-
                          paste(rep(baseline.years, each = length(baseline.ages))
                               ,rep(baseline.ages, length(baseline.years))
                               ,sep = ".")

        ## SRB
        srb.prior <-
            vt(new.clone(vt(srb.prior)
                        ,length = end.point * ncol(srb.prior)
                        ,dim = c(ncol(srb.prior), end.point)
                         ))
        colnames(srb.prior) <- srb.years

        ## Pop counts
        pop.count.prior <-
            list(female = vt(new.clone(vt(pop.count.prior[["female"]])
                                      ,length = end.point * ncol(pop.count.prior[["female"]])
                                      ,dim = c(ncol(pop.count.prior[["female"]]), end.point)
                                       ))
                ,male = vt(new.clone(vt(pop.count.prior[["male"]])
                                    ,length = end.point * ncol(pop.count.prior[["male"]])
                                    ,dim = c(ncol(pop.count.prior[["male"]]), end.point)
                                     ))
                 )
        for(i in 1:2) colnames(pop.count.prior[[i]]) <-
                          paste(rep(proj.years, each = length(mig.ages))
                               ,rep(mig.ages, length(proj.years))
                               ,sep = ".")
    }


    ## -------* Save and exit

    if(return.list) {
        out <- list(fert.rate.prior = fert.rate.prior[],
                   surv.prop.prior = list(female = surv.prop.prior$female[],
                                           male = surv.prop.prior$male[]),
                    mig.prop.prior = list(female = mig.prop.prior$female[],
                                          male = mig.prop.prior$male[]),
                   baseline.count.prior = list(female = baseline.count.prior$female[],
                                               male = baseline.count.prior$male[]),
                   srb.prior = srb.prior[],
                   pop.count.prior = list(female = pop.count.prior$female[],
                                               male = pop.count.prior$male[])
                   )
        names(out) <- paste(name.pref, names(out), name.suf, sep = "")
        return(out)

    } else {

        message("Saving all objects to directory ", output.dir)
        saveff(end.point = (batches.done + 1) * batch.size
                                # save is done before
                                # incrementing
                                # 'batches.done'
              ,filename.tag = ""
              ,save.dir = output.dir
               )
        for(y in c("fert.rate.prior", "surv.prop.prior", "mig.prop.prior"
                  ,"baseline.count.prior", "srb.prior", "pop.count.prior")) {
            assign(paste(name.pref, y, name.suf, sep = "")
                  ,get(y)
                  ,envir = .GlobalEnv
                   )
            message(paste(name.pref, y, name.suf, sep = "")
                   ," placed in .GlobalEnv")
            rm(list = y)
            gc(verbose = FALSE)
        }
    }

    ## remove '_TEMP' objects
    file.remove(paste(output.dir
                     ,grep(x = dir(output.dir), pattern = "_TEMP", value = TRUE)
                     ,sep = "/"
                      ))

}
