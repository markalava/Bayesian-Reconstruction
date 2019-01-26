################################################################################
###
### TITLE:              pop_reconstruction_functions.R
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Functions to do population reconstruction developed for
###                     "Reconstructing Past Populations with Uncertainty from
###                     Fragmentary Data" submitted to Journal of the American
###                     Statistical Assocation".
###
### REFERENCE:          Wheldon, M. C., Raftery, A. E., Clark, S. J.,
###                     & Gerland, P. (2013). Reconstructing Past
###                     Populations with Uncertainty from Fragmentary
###                     Data. Journal of the American Statistical
###                     Association, 108(501),
###                     96–110. http://doi.org/10.1080/01621459.2012.737729
###
### LICENCE:            Released under the Creative Commons BY-NC-SA Licence
###                     (https://creativecommons.org).
###
### DISCLAIMER:         The views and opinions expressed in this work
###                     are those of the authors and do not
###                     necessarily represent those of the United
###                     Nations. This work has not been formally
###                     edited and cleared by the United Nations.
###
###-----------------------------------------------------------------------------
###
### SYNOPSIS:
###
### Three sets of functions used in the simulation, reconstruction, and model
### checking process are defined here, as well as Metropolis proposal variances.
###
###
### _Auxiliary Functions_
###
### Auxiliary functions used by the main functions for the simulation study and
### the Burkina Faso reconstruction. Function names are prefixed 'popReconAux'.
###
###
### _ Simulation study functions_
###
### This function produces a single replicate of the simulation study. It is
### designed to be used with snowFT() and performParallel() but can be used
### without either, in which case the simulation is done in serial. Function
### names are prefixed 'simStudy'.
###
###
### _Reconstruction Functions_
###
### These are the main functions used to reconstruct the female population of
### Burkina Faso as described in Section 5. There are separate functions for the
### 'normal likelihood' and 't2 likelihood' models, as described in Section 5.
### Function names are prefixed 'popRecon' and 'popReconModCheck'.
###
###
### _Metropolis proposal variances_
###
### Proposal variances used for the MCMC samplers, set to ensure acceptance
### proportions lie between 0.1 and 0.5.
###
################################################################################


###
### * Auxiliary Functions
###
################################################################################

###
### Logit and inverse
###

popReconAux.logit <- function(p) log(p / (1 - p))
popReconAux.invlogit <- function(x)
    {
     if(any(is.infinite(exp(x)))) {
         y <- x
         y[is.infinite(exp(x))] <- 1
         y[!is.infinite(exp(x))] <-
             estMod.invlogit.nov17(y[!is.infinite(exp(x))])
         return(y)
     }
     else return(exp(x) / (1 + exp(x)))
    }


###
### Inverse gamma
###

popReconAux.rinvGamma <- function(n, shape, scale)
  {
    return(1/rgamma(n, shape = shape, rate = scale))
  }

popReconAux.dinvGamma <- function(x, shape, scale, log = FALSE)
  {
    if(log) d <- shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
    else d <- scale^shape / gamma(shape) * (1/x)^(shape + 1) * exp(-scale/x)
    return(d)
  }

popReconAux.qinvGamma <- function(p, shape, scale, lower.tail = TRUE)
  {
      return(1/qgamma(1-p, shape = shape, rate = scale
                      ,lower.tail = lower.tail)
             )
  }

popReconAux.pinvGamma <- function(q, shape, scale, lower.tail = TRUE)
{
    pgamma(1/q, shape = shape, rate = scale, lower.tail = !lower.tail)
}


###
### Calculate conditional variances in an MCMC chain.
###

popReconAux.cond.vars <-
    function(chain, tryC = TRUE, ..., try.finally = NULL)

    ##------------------------------------------------------------
    ##
    ## Calculate conditional variances in an MCMC chain. Output is used to tune
    ## Metropolis proposals.
    ##
    ## ARGUMENTS
    ##   chain      : mcmc object
    ##   tryC       : logical. Should lm be wrapped by tryCatch()?
    ##   try.finally: finally function for tryCatch
    ##   ...        : passed to tryCatch()
    ##
    ## DESCRIPTION
    ##   'Chain' is the output from an MCMC sampler consisting of
    ##   N iterations, K variables, x_1, ..., x_K.
    ##
    ##   This function calculates Var(x_k | x_{-k}) for all variables
    ##   in 'chain' using linear regression of x_k on x_{-k}.
    ##
    ## CREATOR
    ##   Mark C. Wheldon
    ##
    ## REFERENCE
    ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
    ##   P. (2013). Reconstructing Past Populations with Uncertainty
    ##   from Fragmentary Data. Journal of the American Statistical
    ##   Association, 108(501), 96–110.
    ##   http://doi.org/10.1080/01621459.2012.737729
    ##
    ## LICENCE:
    ##   Released under the Creative Commons BY-NC-SA Licence
    ##   (https://creativecommons.org).
    ##
    ##------------------------------------------------------------

{
    require(coda)

    ## Coerce 'chain' to a data frame
    x <- as.data.frame(chain)

    ## List of formulae for all regressions with each
    ## variable as response, all others as predictors
    form.list <- list()
    for(i in 1:ncol(x)) {
        form.list[[i]] <- formula(x[1,c(i, (1:ncol(x))[-i])])
    }
    names(form.list) <- colnames(x)

    ## Perform regressions, take residual *std deviations*
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

    ## Convert to variances
    return(resid.std.devs^2)

}


###
### * Simulation Study Functions
###
################################################################################

simStudy.estimation.once <-
    function(## the array determining number of times to run
             rep.ar.count

             ## sim stats
             ,runs.per.node, desired.overall.reps, overall.reps
             ,cluster.size

             ## coverage level
             ,alpha

             ## algorithm parameters
             ,start.iter, start.burn.in, prop.vars
             ,max.tune.reruns = 5
             ,max.iter = 15000, max.burn.in = 5000
             ,min.iter = 9000, min.burn.in = 600
             ,runRaftLew = TRUE
             ,rlPctl = 0.9 # take rlPctl of Raftery-Lewis 'N' within each
                           # parameter type
             ,checkAR = TRUE # check acceptance proportions
             ,condVarMaxRows = 1E4 # maximum number of rows to use for cond vars
             ,ar.lower = 0.1, ar.upper = 0.5

             ## est model arguments
             ,ccmp.f, age.size, census.columns, fert.rows
             ,s.tol, verb

             ## hyper parameters
             ,hyper.params

             ## source functions
             ,estMod.f.source
             ,aux.fs.source = NULL

             ## function names
             ,estMod.f

             ## true values
             ,true.values

             ## profile?
             ,Rprof.file.path = FALSE

             ## sink?
             ,sink.file.path = FALSE

             ## dump?
             ,dump.file.path = FALSE

             ## root path
             ,working.directory = getwd()

             ## save output?
             ,save.file = NULL
             )

    ##
    ## CREATOR
    ##   Mark C. Wheldon
    ##
    ## REFERENCE
    ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
    ##   P. (2013). Reconstructing Past Populations with Uncertainty
    ##   from Fragmentary Data. Journal of the American Statistical
    ##   Association, 108(501), 96–110.
    ##   http://doi.org/10.1080/01621459.2012.737729
    ##
    ## LICENCE
    ##   Released under the Creative Commons BY-NC-SA Licence
    ##   (https://creativecommons.org).
    ##
{

    ##-------* HOUSEKEEPING

    setwd(working.directory)


    ##-------** Check inputs (not exhaustive)

    if(runRaftLew) {
        rd.nmin <- qnorm((0.95+1)/2)^2 * 0.025*(1-0.025)/0.0125^2
        if(min.iter < rd.nmin) {
            warning(paste("'min.iter' too low for Raftery-Lewis diagnostic; 'min.iter' changed to", ceiling(rd.nmin)))
            min.iter <- ceiling(rd.nmin)
        }
    }
    if(start.iter < min.iter) {
        warning("'start.iter' < 'min.iter'; 'start.iter' changed to ", min.iter)
        start.iter <- min.iter
    }
    if(max.iter < start.iter) {
        warning("'max.iter' < 'start.iter'; 'max.iter' changed to ", start.iter)
        max.iter <- start.iter
    }
    if(start.burn.in < min.burn.in) {
        warning("'start.burn.in' < 'min.burn.in'; 'start.burn.in' changed to ", min.burn.in)
        start.burn.in <- min.burn.in
    }
    if(max.burn.in < start.burn.in) {
        warning("'max.burn.in' < 'start.burn.in'; 'max.burn.in' changed to ", start.burn.in)
    }
    if(verb) cat("\n\n *** Iterations ***\nmin.iter = ", min.iter, "\nstart.iter = ", start.iter, "\nmax.iter = ", max.iter)


    ##-------* Source or define helper functions

    ##-------** Aux functions

    if(!missing(aux.fs.source)) source(aux.fs.source)


    ##-------** Estimation model function

    source(estMod.f.source)


    ##-------** Generate initial estimates function

    ## This function performs step (2) in Section 4.2 of the manuscript.

    generate.new.initial.estimates <-
        function(## True vitals
                 fertRate.true, survProp.true, migProp.true
                 ,baselineCount.true, census.true

                 ## *Variance* parameters
                 ,fert.var, surv.var, mig.var, popCount.var
                 )
        {
            asFertMEAS.mat <-
                exp(log(fertRate.true) + rnorm(length(fertRate.true)
                                               ,mean = 0, sd = sqrt(fert.var)
                                               )
                    )
            asSurvMEAS.mat <-
                popReconAux.invlogit(popReconAux.logit(survProp.true) +
                                     rnorm(length(survProp.true)
                                           ,mean = 0, sd = sqrt(surv.var))
                                     )
            asMigMEAS.mat <-
                migProp.true + rnorm(length(migProp.true)
                                     ,mean = 0, sd = sqrt(mig.var)
                                     )
            baselineMEAS.mat <-
                exp(log(baselineCount.true) + rnorm(length(baselineCount.true)
                                                    ,mean = 0, sd = sqrt(popCount.var)
                                                    )
                    )
            censusMEAS.mat <-
                exp(log(census.true) + rnorm(length(census.true)
                                             ,mean = 0
                                             ,sd = sqrt(popCount.var))
                    )

            out <-
                list(asFertMEAS.mat = asFertMEAS.mat, sigmasq.f = fert.var
                     ,asSurvMEAS.mat = asSurvMEAS.mat, sigmasq.s = surv.var
                     ,asMigMEAS.mat = asMigMEAS.mat, sigmasq.g = mig.var
                     ,baselineMEAS.mat = baselineMEAS.mat
                     ,censusMEAS.mat = censusMEAS.mat
                     ,sigmasq.n = popCount.var
                     )
            return(out)
        }


    ##-------** Results summary function

    ## This function calculates the coverage of the posterior intervals.

    posterior.interval.coverage <-
        function(chain  # chain from mcmc sampler
                 ,truth # true values in matrix form
                 ,alpha # 1-alpha = desired coverage. can be a vector.
                 )
        {
            require(reshape)

            ##--- Prepare truth matrix ---#
            true.m <- melt(truth)
            colnames(true.m) <- c("age", "years", "truth")

            ##--- Calculate quantiles of chain ---#
            prob.vec <-
                c(0.5, as.numeric(sapply(alpha, function(z) c(z/2, 1-z/2))))
            q.vital1 <-
                apply(chain, 2
                      ,function(z) quantile(z, probs = prob.vec)
                      )

            ##--- Make list, one component for each alpha ---#
            q.vital.li <- as.list(rep(0, length(alpha)))
            names(q.vital.li) <- alpha
            for(i in 1:length(q.vital.li))
                q.vital.li[[i]] <- q.vital1[c(2*i, 1, 2*i+1),]

            ##--- Check coverage ---#
            out <- lapply(q.vital.li, function(z, t) {
                rownames(z) <- c("lower", "median", "upper")
                q.vital.m <- as.data.frame.table(z)
                q.vital2 <-data.frame(quant = q.vital.m$Var1
                                      ,colsplit(q.vital.m$Var2, split = "\\."
                                                ,names = c("years", "age"))
                                      ,value = q.vital.m$Freq)
                q.vital3 <-
                    cast(q.vital2, years + age ~ quant)

                ##--- Merge and compare ---#
                merged <-
                    merge(t, q.vital3, sort = FALSE)
                merged$covered <-
                    apply(merged[,3:6], 1, FUN = function(z) {
                        z["truth"] >= z["lower"] & z["truth"] <= z["upper"]
                    })
                return(merged)
            }
                          ,t = true.m)

            ##--- output ---#
            return(out)
        }


    ##-------** Conditional variances function

    ## This function is used for 'auto-tuning' the MCMC. It estimates the
    ## conditional variances of each variable in the MCMC chain, given all the
    ## other variables.

    mcmc.chain.conditional.variances <-
        function(chain, tryC, ..., try.finally = NULL)

            ##-- Calculate conditional variances in mcmc output --#
            ##
            ## ARGUMENTS
            ##   chain      : mcmc object
            ##   tryC       : logical. Should lm be wrapped by tryCatch()?
            ##   try.finally: finally function for tryCatch
            ##   ...        : passed to tryCatch()
            ##
            ## DESCRIPTION
            ##   'Chain' is the output from an MCMC sampler consisting of
            ##   N iterations, K variables, x_1, ..., x_K.
            ##
            ##   This function calculates Var(x_k | x_{-k}) for all variables
            ##   in 'chain' using linear regression of x_k on x_{-k}.
            ##
            ## CREATOR
            ##   Mark C. Wheldon
            ##
            ## REFERENCE
            ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
            ##   P. (2013). Reconstructing Past Populations with Uncertainty
            ##   from Fragmentary Data. Journal of the American Statistical
            ##   Association, 108(501), 96–110.
            ##   http://doi.org/10.1080/01621459.2012.737729
            ##
            ## LICENCE
            ##   Released under the Creative Commons BY-NC-SA Licence
            ##   (https://creativecommons.org).
            ##

        {
            require(coda)

            ## Coerce 'chain' to a data frame
            x <- as.data.frame(chain)

            ## List of formulae for all regressions with each variable as
            ## response, all others as predictors
            form.list <- list()
            for(i in 1:ncol(x)) {
                form.list[[i]] <- formula(x[1,c(i, (1:ncol(x))[-i])])
            }
            names(form.list) <- colnames(x)

            ## Perform regressions, take residual *std deviations*
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

            ## Convert to variances
            return(resid.std.devs^2)
        }



    ##-------** Some error handling functions

    ## tryCatch() raftery.diag
    raftery.diag.tryC <-
        function(data, q = 0.025, r = 0.005, s = 0.95
                 ,converge.eps = 0.001
                 ,try.error = function(e) {
                     cat("\n\nError with raftery.diag:\n")
                     print(e)
                     cat("\n")
                 }
                 )
        {
            tryCatch(raftery.diag(data, q, r, s, converge.eps)
                     ,error = try.error)
        }

    ## Function to handle error with lm() when calculating conditional
    ## variances.
    condVar.tryError.fun <- function(e) {
        cat("\n\nError with lm() when calculating conditional variances:\n")
        print(e)
        cat("\nFull MCMC sample saved as\n", getwd(), "/lm_Error_"
            ,rep.ar.count, "_", i, "_"
            ,"FULLOUT_"
            ,format(Sys.time(), format = "%Y%m%d_%H%M")
            ,".RData\n"
            ,sep = "")
        save(simulationStudy.mcmc.samples
             ,file = paste(getwd(), "/lm_Error_"
              ,rep.ar.count, "_", i, "_"
              ,"FULLOUT_"
              ,format(Sys.time(), format = "%Y%m%d_%H%M")
              ,".RData", sep = "")
             )

    }


    ##-------* DO THE SIMULATION

    ## Storage

    mcmc.info <- list()

    ##-------** Loop over runs.per.node

    for(i in 1:runs.per.node) {


        ##-------*** Save Rprof and some output logs if filenames given

        ## Rprofile files
        if(is.character(Rprof.file.path)) {
            fp <- paste(Rprof.file.path
                        ,"RprofFile_node", rep.ar.count
                        ,"_iter", i, sep = "")
            Rprof(filename = fp)
        }

        ## Sink
        if(is.character(sink.file.path)) {
            fs <- paste(sink.file.path
                        ,"Sink_node", rep.ar.count
                        ,"_iter", i, ".Rout", sep = "")
            sink(file = fs)
        }


        ##-------*** Steps (1)--(3): Generate initial estimates

        ## Make sure no negative populations (Step (3))

        projMat <-
            matrix(-1, nrow = nrow(true.values$asFertTRUE.mat)
                   ,ncol = ncol(true.values$asFertTRUE.mat))

        while(sum(projMat < 1) > 0) {

            ## Draw variance parameters
            fert.varR <-
                popReconAux.rinvGamma(1, hyper.params$al.f, hyper.params$be.f)
            surv.varR <-
                popReconAux.rinvGamma(1, hyper.params$al.s, hyper.params$be.s)
            mig.varR <-
                popReconAux.rinvGamma(1, hyper.params$al.g, hyper.params$be.g)
            popCount.varR <-
                popReconAux.rinvGamma(1, hyper.params$al.g, hyper.params$be.g)

            ## Draw vital rates
            dataGenArgs <-
                list(fertRate.true = true.values$asFertTRUE.mat
                     ,survProp.true = true.values$asSurvTRUE.mat
                     ,migProp.true = true.values$asMigTRUE.mat
                     ,baselineCount.true = true.values$baselineTRUE.mat
                     ,census.true = true.values$censusTRUE.mat
                     ,fert.var = fert.varR, surv.var = surv.varR
                     ,mig.var = mig.varR, popCount.var = popCount.varR
                     )

            ## Make initial estimates for this replicate
            simData <-
                do.call(generate.new.initial.estimates, args = dataGenArgs)

            ## Project
            projMat <-
                do.call(match.fun(ccmp.f)
                        ,args = list(pop = simData$baselineMEAS.mat
                         ,fert = simData$asFertMEAS.mat
                         ,surv = simData$asSurvMEAS.mat
                         ,mig = simData$asMigMEAS.mat
                         ,proj.steps = ncol(simData$asFertMEAS.mat)
                         ,age.int = 5
                         ))
        }

        if(verb) {
            cat("\n\n\n****** Node ", rep.ar.count, "******\n****** this node iter ", i, " ******\n")
            cat("\nSimData\n")
            print(simData)
            cat("\n\nprojMat\n")
            print(projMat)
        }


        ##-------*** Step (4): Run the MCMC sampler

        ##-------**** Set-up

        ## Arguments for estimation model
        estModArgs <-
            c(## inverse gamma parameters
              hyper.params

              ## the rest
              ,list(n.iter = start.iter, burn.in = start.burn.in
                    ,ccmp.f = ccmp.f
                    ,init.f = simData$asFertMEAS.mat
                    ,init.s = simData$asSurvMEAS.mat
                    ,init.g = simData$asMigMEAS.mat
                    ,init.b = simData$baselineMEAS.mat
                    ,init.sigmasq.f = 5
                    ,init.sigmasq.s = 5
                    ,init.sigmasq.g = 5
                    ,init.sigmasq.n = 5
                    ,mean.f = simData$asFertMEAS.mat
                    ,mean.s = simData$asSurvMEAS.mat
                    ,mean.g = simData$asMigMEAS.mat
                    ,mean.b = simData$baselineMEAS.mat
                    ,prop.vars = prop.vars
                    ,pop.data = simData$censusMEAS.mat
                    ,proj.periods = ncol(simData$asFertMEAS.mat)
                    ,age.size = age.size
                    ,census.columns = census.columns
                    ,fert.rows = fert.rows
                    ,s.tol = s.tol
                    ,verb = verb
                    )
              )

        ## Dump objects if file path given (dump some objects into a file)
        if(is.character(dump.file.path)) {
            fd <- paste(dump.file.path
                        ,"Dump_node", rep.ar.count
                        ,"_iter", i, ".R", sep = "")
            DOTrandomDOTseed <- .Random.seed
            dump(list = c("simData", "projMat", "estModArgs"
                 ,"DOTrandomDOTseed")
                 ,file = fd)
        }


        ##-------**** Run sampler with some auto-tuning

        tuneRepeat <- TRUE
        reruns <- 0

        while(tuneRepeat && reruns < max.tune.reruns) {

            ## Increment
            reruns <- reruns + 1
            if(estModArgs$verb) {
                cat("\n\n*** RUN ", reruns, " ***\n\n")
                print(estModArgs)
            }

            ## Initialize indicators and other objects:

            ## vital rate conditional variances
            vitalCondVars <- NULL

            ## Metropolis acceptance rates OK? (logical; updated if 'checkAR' == TRUE)
            condn.acc.prop <- TRUE

            ## Iterations and burn-in sufficient? (logical; updated if 'runRaftLew' == TRUE)
            condn.max.iter <- TRUE
            condn.max.burn.in <- TRUE

            ## Run the MCMC
            simulationStudy.mcmc.samples <-
                do.call(estMod.f, args = estModArgs)

            ## Try to auto-tune by re-running if acceptance rates not within
            ## 'ar.lower'--'ar.upper' percent.
            if(checkAR) {
                condn.acc.prop <-
                    all(unlist(sapply(simulationStudy.mcmc.samples$alg.stats$acceptance.proportions[1:4]
                                      ,function(z) return(z >= ar.lower &
                                                          z <= ar.upper)
                                      )))

                if(estModArgs$verb) {
                    cat("\n\nAcc Propns\n")
                    print(simulationStudy.mcmc.samples$alg.stats$acceptance.proportions)
                }
            }

            ## Calculate Raftery-Lewis recommendations
            if(runRaftLew) {

                raftLew <- list()

                raftLew$var.q025.r0125.s95 <-
                    raftery.diag.tryC(simulationStudy.mcmc.samples$variances.mcmc
                                      ,q = 0.025, r = 0.0125, s = 0.8)
                raftLew$var.q975.r0125.s95 <-
                    raftery.diag.tryC(simulationStudy.mcmc.samples$variances.mcmc
                                      ,q = 0.975, r = 0.0125, s = 0.8)

                raftLew$fert.q025.r0125.s95 <-
                    raftery.diag.tryC(log(simulationStudy.mcmc.samples$fert.rate.mcmc)
                                      ,q = 0.025, r = 0.0125, s = 0.8)
                raftLew$fert.q975.r0125.s95 <-
                    raftery.diag.tryC(log(simulationStudy.mcmc.samples$fert.rate.mcmc)
                                      ,q = 0.975, r = 0.0125, s = 0.8)

                raftLew$surv.q025.r0125.s95 <-
                    raftery.diag.tryC(popReconAux.logit(simulationStudy.mcmc.samples$surv.prop.mcmc)
                                      ,q = 0.025, r = 0.0125, s = 0.8)
                raftLew$surv.q975.r0125.s95 <-
                    raftery.diag.tryC(popReconAux.logit(simulationStudy.mcmc.samples$surv.prop.mcmc)
                                      ,q = 0.975, r = 0.0125, s = 0.8)

                raftLew$mig.q025.r0125.s95 <-
                    raftery.diag.tryC(simulationStudy.mcmc.samples$mig.prop.mcmc
                                      ,q = 0.025, r = 0.0125, s = 0.8)
                raftLew$mig.q975.r0125.s95 <-
                    raftery.diag.tryC(simulationStudy.mcmc.samples$mig.prop.mcmc
                                      ,q = 0.975, r = 0.0125, s = 0.8)

                raftLew$baseline.q025.r0125.s95 <-
                    raftery.diag.tryC(log(simulationStudy.mcmc.samples$baseline.count.mcmc)
                                      ,q = 0.025, r = 0.0125, s = 0.8)
                raftLew$baseline.q975.r0125.s95 <-
                    raftery.diag.tryC(log(simulationStudy.mcmc.samples$baseline.count.mcmc)
                                      ,q = 0.975, r = 0.0125, s = 0.8)

                ## Upper rlPctl percentiles within each level of hierarchy (each
                ## parameter type vital rates and variances)
                raftLew.type <-
                    with(raftLew
                         ,list(variances = list(var.q025.r0125.s95
                               ,var.q975.r0125.s95)
                               ,fert.rate = list(fert.q025.r0125.s95
                                ,fert.q975.r0125.s95)
                               ,surv.prop = list(surv.q025.r0125.s95
                                ,surv.q975.r0125.s95)
                               ,mig.prop = list(mig.q025.r0125.s95
                                ,mig.q975.r0125.s95)
                               ,baseline.count = list(baseline.q025.r0125.s95
                                ,baseline.q975.r0125.s95)
                               )
                         )
                raftLew.pctl <-
                    sapply(raftLew.type[-1], function(z) {
                        quantile(sapply(z, function(y) {
                            y[["resmatrix"]][,"N"]
                        }), rlPctl)
                    })
                raftLew.pctlVars <-
                    sapply(raftLew.type[1], function(z) {
                        max(sapply(z, function(y) {
                            y[["resmatrix"]][,"N"]
                        }))
                    })


                ## Max rlPctl from each variable type
                raftLew.maxPctl <- max(raftLew.pctl, raftLew.pctlVars)

                ## Maximum of burn-in (over all parameters)
                raftLew.maxBI <-
                    max(unlist(sapply(raftLew
                                      ,FUN = function(z) z[["resmatrix"]][,"M"])))

                ## Compare R-L recommendations with actual iterations
                rlCheck <-
                    max(min.iter, min(raftLew.maxPctl, max.iter))
                condn.max.iter <-
                    estModArgs$n.iter >= rlCheck

                if(estModArgs$verb) cat("\n\nrl Iters ", rlCheck)

                ## burn in
                rlCheck <-
                    max(min.burn.in, min(raftLew.maxBI, max.burn.in))
                condn.max.burn.in <-
                    estModArgs$burn.in >= rlCheck

                if(estModArgs$verb) cat("\n\nrl Burn in ", rlCheck)

            }

            ## If conditional acceptance ratios within 'ar.lower'--'ar.upper'
            ## and chain lengths OK, step out of loop...

            if(condn.acc.prop && condn.max.iter && condn.max.burn.in) {
                tuneRepeat <- FALSE
            } else {
                ## ...otherwise set chain lengths to raftery lewis
                ## recommendations, set proposal variances and re-run.
                if(!condn.acc.prop) {
                    if(estModArgs$verb) cat(paste("\n\nAcceptance ratios not between"
                                                  ,ar.lower, "and", ar.upper, "\n"))

                    vitalCondVars <- list()

                    ## Compute conditional variances of each variable using
                    ## linear regression and use result to set new Metropolis
                    ## proposal variances. Use no more than condVarMaxRows rows
                    ## of data to ensure do not run out of memory.
                    if(estModArgs$n.iter > condVarMaxRows) {
                        condVarRowStart <-
                            estModArgs$n.iter - condVarMaxRows
                        condVarRowEnd <- estModArgs$n.iter
                    } else {
                        condVarRowStart <- 1
                        condVarRowEnd <- estModArgs$n.iter
                    }

                    vitalCondVars$fert.rate <-
                        matrix(mcmc.chain.conditional.variances(chain = log(simulationStudy.mcmc.samples$fert.rate.mcmc[condVarRowStart:condVarRowEnd,])
                                                                ,tryC = TRUE
                                                                ,error = condVar.tryError.fun
                                                                )
                               ,nrow = nrow(estModArgs$mean.f[estModArgs$fert.rows,])
                               ,ncol = ncol(estModArgs$mean.f[estModArgs$fert.rows,])
                               ,dimnames = dimnames(estModArgs$mean.f[estModArgs$fert.rows,])
                               )

                    ## Infininte values might be a problem if proportions are
                    ## close to 0 or 1. Remove any infinities from logit
                    ## transformed by deleting the entire row. This is
                    ## acceptable since the conditional variances are only used
                    ## for tuning.
                    logitSurv <-
                        popReconAux.logit(simulationStudy.mcmc.samples$surv.prop.mcmc)
                    logitSurv.infRows <-
                        which(is.infinite(logitSurv), arr.ind = TRUE)[,1]
                    if(length(logitSurv.infRows) > 0) {
                        logitSurv.notInf <-
                            as.data.frame(logitSurv[-logitSurv.infRows,])
                    } else logitSurv.notInf <- logitSurv
                    if(nrow(logitSurv.notInf) < min.iter)
                        warning("Number of non infinite rows in 'logit(simulationStudy.mcmc.samples$surv.prop.mcmc)' < 'min.iter'\nnrow(logitSurv.notInf) == "
                                ,nrow(logitSurv.notInf), "\n"
                                )
                    if(nrow(logitSurv.notInf) > condVarMaxRows) {
                        condVarRowStart.surv <-
                            nrow(logitSurv.notInf) - condVarMaxRows
                        condVarRowEnd.surv <- nrow(logitSurv.notInf)
                    } else {
                        condVarRowStart.surv <- 1
                        condVarRowEnd.surv <- nrow(logitSurv.notInf)
                    }
                    vitalCondVars$surv.prop <-
                        matrix(mcmc.chain.conditional.variances(chain = logitSurv.notInf[condVarRowStart.surv:condVarRowEnd.surv,]
                                                                ,tryC = TRUE
                                                                ,error = condVar.tryError.fun
                                                                )
                               ,nrow = nrow(estModArgs$mean.s)
                               ,ncol = ncol(estModArgs$mean.s)
                               ,dimnames = dimnames(estModArgs$mean.s)
                               )

                    vitalCondVars$mig <-
                        matrix(mcmc.chain.conditional.variances(chain = simulationStudy.mcmc.samples$mig.prop.mcmc[condVarRowStart:condVarRowEnd,]
                                                                ,tryC = TRUE
                                                                ,error = condVar.tryError.fun
                                                                )
                               ,nrow = nrow(estModArgs$mean.g)
                               ,ncol = ncol(estModArgs$mean.g)
                               ,dimnames = dimnames(estModArgs$mean.g)
                               )

                    vitalCondVars$population.count <-
                        matrix(mcmc.chain.conditional.variances(log(simulationStudy.mcmc.samples$baseline.count.mcmc[condVarRowStart:condVarRowEnd,])
                                                                ,tryC = TRUE
                                                                ,error = condVar.tryError.fun
                                                                )
                               ,nrow = nrow(estModArgs$mean.b)
                               ,ncol = ncol(estModArgs$mean.b)
                               ,dimnames = dimnames(estModArgs$mean.b)
                               )

                    if(any(unlist(lapply(vitalCondVars, "is.na")))) {
                        warning("Some conditional variances estimated to be zero; perhaps n.iter is too low. Conditional variances NOT updated.")
                    } else {
                        estModArgs$prop.vars <-
                            lapply(vitalCondVars, "*", 2.3^2)
                    }
                }

                if(!condn.max.iter) {
                    if(estModArgs$verb) cat("\nIterations not sufficient")

                    estModArgs$n.iter <-
                        round(as.numeric(max(min.iter
                                             ,min(max.iter
                                                  ,raftLew.maxPctl))
                                         ))
                }
                if(!condn.max.burn.in) {
                    if(estModArgs$verb) cat("\nBurn in not sufficient\n\n")
                    round(estModArgs$burn.in <-
                          as.numeric(max(min.burn.in
                                         ,min(max.burn.in
                                              ,raftLew.maxBI))
                                     ))
                }
            }
        }


        ##-------**** Save all output

        if(!is.null(save.file)) save(simulationStudy.mcmc.samples
                                     ,file = paste(save.file, sep = "")
                                     ,compress = TRUE
                                     )


        ##-------* SUMMARIZE

        ##-------** Coverage

        ## Fertility
        fert.rate.postQuants <-
            do.call(posterior.interval.coverage
                    ,args = list(chain = simulationStudy.mcmc.samples$fert.rate.mcmc
                     ,truth = true.values$asFertTRUE.mat, alpha = alpha)
                    )

        ## Survival
        surv.prop.postQuants <-
            do.call(posterior.interval.coverage
                    ,args = list(chain = simulationStudy.mcmc.samples$surv.prop.mcmc
                     ,truth = true.values$asSurvTRUE.mat, alpha = alpha)
                    )

        ## Migration
        mig.prop.postQuants <-
            do.call(posterior.interval.coverage
                    ,args = list(chain = simulationStudy.mcmc.samples$mig.prop.mcmc
                     ,truth = true.values$asMigTRUE.mat, alpha = alpha)
                    )

        ## Baseline
        baseline.count.postQuants <-
            do.call(posterior.interval.coverage
                    ,args =
                    list(chain = simulationStudy.mcmc.samples$baseline.count.mcmc
                         ,truth = true.values$baselineTRUE.mat, alpha = alpha)
                    )

        ## lx
        lx.postQuants <-
            do.call(posterior.interval.coverage
                    ,args =
                    list(chain = simulationStudy.mcmc.samples$lx.mcmc
                         ,truth = true.values$censusTRUE.mat, alpha = alpha)
                    )


        ##-------** List for output

        mcmc.info[[i]] <-
            list(init.vals = simulationStudy.mcmc.samples$init.vals
                 ,fixed.params = simulationStudy.mcmc.samples$fixed.params
                 ,alg.stats = simulationStudy.mcmc.samples$alg.stats
                 ,alg.params = simulationStudy.mcmc.samples$alg.params
                 ,raftLew = switch(runRaftLew, "TRUE" = raftLew
                  ,"FALSE" = NULL)
                 ,var.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$variances.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,fert.rate.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$fert.rate.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,fert.rate.postQuants = fert.rate.postQuants
                 ,surv.prop.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$surv.prop.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,surv.prop.postQuants = surv.prop.postQuants
                 ,mig.prop.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$mig.prop.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,mig.prop.postQuants = mig.prop.postQuants
                 ,baseline.count.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$baseline.count.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,baseline.count.postQuants =
                 baseline.count.postQuants
                 ,lx.summstat = switch(estModArgs$n.iter > 1E3
                  ,"TRUE" = summary(simulationStudy.mcmc.samples$lx.mcmc)
                  ,"FALSE" = "'n.iter' < 1E3")
                 ,lx.postQuants = lx.postQuants
                 ,cond.vars = switch(!is.null(vitalCondVars) && !any(is.na(unlist(vitalCondVars)))
                  ,"TRUE" = vitalCondVars
                  ,"FALSE" = NULL
                  )
                 ,sim.stats = list(cluster.size = cluster.size
                  ,desired.overall.reps = desired.overall.reps
                  ,runs.per.node = runs.per.node
                  ,overall.reps = overall.reps
                  ,reruns = reruns
                  ,alpha = alpha
                  )
                 ,tuning.indicators = list(chain.length.OK = condn.max.iter
                  ,burn.in.OK = condn.max.burn.in
                  ,acceptance.props.OK = condn.acc.prop
                  )
                 )
        print(traceback())
        Rprof(NULL)
        sink()
    }
    print(traceback())
    return(mcmc.info)

}


###
### * Reconstruction Functions
###
################################################################################

###
### ** 'Normal-likelihood' SAMPLER
###
################################################################################

###
### *** Cohort component population projection (CCMPP)
###
################################################################################

## This function implements the deterministic cohort component method of
## population projection. It is described in Chapter 6 of
##      Preston, S. H., Heuveline, P. & Guillot, M. "Demography: Measuring and
##      Modeling Population Processes."  Blackwell: Malden, Massachusetts. 2001.

popRecon.ccmp.female <-
    function(pop, surv, fert, srb = 1.05, mig
             ,proj.steps, age.int = 5
             ,label.dims = FALSE, base.year = "1960"
             )

    ##--------------------------------------------------
    ##
    ## ARGUMENTS:
    ##
    ##   pop     :  population count at baseline
    ##   fert    :  matrix of age specific fertility rates NOT yet
    ##                mulitplied by age.int
    ##   srb     :  sex ratio at birth matrix
    ##   surv    :  Survivorship probabilities: the probability of
    ##                reaching the age at the start of the interval.
    ##              The first row should be nL0/(n*l0).
    ##              The last row is survival for age.int years in the open
    ##                interval
    ##   mig     :  Migration. Net number of migrants as a
    ##                 _proportion_ of prev time period's population
    ##   proj.steps
    ##           :  number of time periods to project forward
    ##                if missing, set to ncol(fert)
    ##   age.int :  needed for correct interpretation of survival
    ##                and fertility rates
    ##   label.dims
    ##           :  should output have dimnames set?
    ##   base.year
    ##           :  start year for projections (aesthetic)
    ##
    ## CREATOR
    ##   Mark C. Wheldon
    ##
    ## REFERENCE
    ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
    ##   P. (2013). Reconstructing Past Populations with Uncertainty
    ##   from Fragmentary Data. Journal of the American Statistical
    ##   Association, 108(501), 96–110.
    ##   http://doi.org/10.1080/01621459.2012.737729
    ##
    ## LICENCE
    ##   Released under the Creative Commons BY-NC-SA Licence
    ##   (https://creativecommons.org).
    ##
    ##
    ##--------------------------------------------------

{

    ##-- Checks --##

    ## If proj.steps is greater than the number of columns in surv,
    ##   reduce or recycle fert, surv matrices


    ##-- Constants --##

    n.age.grps <- length(pop)
    n.surv <- nrow(surv)


    ##-- Derive proj.steps from ncol(fert) --##

    if(missing(proj.steps)) proj.steps <- ncol(fert)


    ##-- Make SRB into matrix if not already --##

    if(is.null(dim(srb))) srb <- matrix(rep(srb, proj.steps))


    ##-- Loop over number of required time periods --##

    ##-- Initialise pop.matrix --##

    pop.mat <- matrix(0, nrow = n.age.grps, ncol = 1 + proj.steps)
    pop.mat[,1] <- pop


    ##-- Initialize leslie matrix --##

    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)


    ##-- Project --##

    for(i in 1:proj.steps)
    {

        ##-- Make the leslie matrix --##

        ##-- Rows 2:(n.age.grps) = survival ratios --##

        lesM[2:n.age.grps,1:(n.age.grps-1)] <-
            diag(surv[-c(1,n.surv),i])
        lesM[n.age.grps,n.age.grps] <- surv[(n.surv),i]

        ##-- First row = fert and survival ratios --##

        k <- 1/(1+srb[i]) * surv[1,i] * 0.5

        dbl.fert <- age.int*fert[,i] + c(age.int*fert[-1,i], 0) *
            surv[-1,i]

        lesM[1,] <- k * dbl.fert

        ##-- Migrants --##

        net.numb.mig <- mig[,i] * pop.mat[,i]

        ##-- Project --##

        pop.mat[,i+1] <-
            lesM %*% (pop.mat[,i] + 0.5 * net.numb.mig) +
                0.5 * net.numb.mig

    }


    ##-- Add dim names --##

    if(label.dims) {
        ages <- seq(from = 0
                    ,to = age.int*(nrow(as.matrix(pop))-1)
                    ,by = age.int)
        yrs <- (0:proj.steps) * age.int + as.numeric(base.year)
        dimnames(pop.mat) <- list(ages, yrs)
    }


    ##-- Output --##

    return(pop.mat)

}


###
### *** Other secondary functions
###
################################################################################

## These are other functions used by the main population reconstruction function.

###
### Creates column names for mcmc objects
###

popRecon.makeColNames <- function(m)
  {
    # m:    matrix of input values

    e <- expand.grid(rownames(m), colnames(m))
    apply(e, 1, FUN = function(z) paste(z[2], z[1], sep = "."))

  }


###
### Log likelihood
###

popRecon.log.lhood <- function(log.n.census, log.n.hat, ll.var, cen.col)
  {
    ## log.n.census and log.n.hat should already be logged

    ## Choose years for which likelihood exists. Passed in from sampler function
    log.n.hat.select <- log.n.hat[,cen.col]

    ## value of log likelihoods
    density <- dnorm(log.n.census,
                     mean = log.n.hat.select,
                     sd = sqrt(ll.var),
                     log = TRUE
                     )

    ## joint log likelihood
    return(sum(density))
  }


###
### Log posterior
###

popRecon.log.post <-
    function(
             ## estimated vitals
             f, s, g, baseline.n
             ## fixed prior means on vitals
             ,prior.mean.f, prior.mean.s
             ,prior.mean.g, prior.mean.b
             ## fixed prior parameters on variance distns
             ,alpha.f, beta.f, alpha.s, beta.s
             ,alpha.g, beta.g
             ,alpha.n, beta.n
             ## updated variances on prior distns
             ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n
             ## value of the log likelihood
             ,log.like
             ## non zero rows of fertility matrix
             ,non.zero.fert
             )
{

    ## Values of prior densities for vitals

    ## Note that log densities are calculated for numerical stability.

    ## f, baseline.n, prior.mean.f, prior.mean.b are logged coming in, s,
    ## prior.mean.s is logit transformed coming in, g and prior.mean.g are not
    ## transformed coming in.
    log.f.prior <- dnorm(as.vector(f[non.zero.fert,])
                         ,mean = as.vector(prior.mean.f[non.zero.fert,])
                         ,sd = sqrt(sigmasq.f)
                         ,log = TRUE)
    log.s.prior <- dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s)
                         ,log = TRUE)
    log.g.prior <- dnorm(g, mean = prior.mean.g
                         ,sd = sqrt(sigmasq.g)
                         ,log = TRUE)
    log.b.prior <- dnorm(baseline.n, mean = prior.mean.b
                         ,sd = sqrt(sigmasq.n)
                         ,log = TRUE)


    ## Values of prior densities for variances

    log.sigmasq.f.prior <-
        log(popReconAux.dinvGamma(sigmasq.f, alpha.f, beta.f))
    log.sigmasq.s.prior <-
        log(popReconAux.dinvGamma(sigmasq.s, alpha.s, beta.s))
    log.sigmasq.g.prior <-
        log(popReconAux.dinvGamma(sigmasq.g, alpha.g, beta.g))
    log.sigmasq.n.prior <-
        log(popReconAux.dinvGamma(sigmasq.n, alpha.n, beta.n))


    ## The log posterior is the SUM of these with the log.like

    return(sum(log.f.prior, log.s.prior, log.g.prior, log.b.prior
               ,log.sigmasq.f.prior
               ,log.sigmasq.s.prior
               ,log.sigmasq.g.prior
               ,log.sigmasq.n.prior
               ,log.like))

}


###
### Acceptance Ratio
###

popRecon.acc.ra <- function(log.prop, log.current)
  {
    min(1, exp(log.prop - log.current))
  }

popRecon.acc.ra.var <-
    function(log.prop.post, log.curr.post, log.prop.var, log.curr.var)
{
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post
               ))
}


###
### *** SAMPLER
###
################################################################################

popRecon.sampler <-
    function(#.. number of iterations and burn-in (not saved)
             n.iter, burn.in = 0

             #.. ccmp function: list with components 'FUN' and
             #   'formals'. Use formals to set, eg, mig.type.
             ,ccmp.f = "popRecon.ccmp.female"

             #.. fixed variance hyper-parameters
             ,al.f, be.f, al.s, be.s, al.g, be.g, al.n, be.n

             #.. fixed prior means
             ,mean.f, mean.s, mean.g, mean.b

             #.. inital values for vitals and variances
             #   *vitals not transformed coming in*
             ,init.f, init.s, init.g, init.b
             ,init.sigmasq.f, init.sigmasq.s, init.sigmasq.g
             ,init.sigmasq.n

             #.. census data: use 'census.columns' to match this
             #   to columns in the projection returned by ccmp
             #   *not transformed coming in*
             ,pop.data

             #.. periods in projection range for which census data exists
             ,census.columns = c(4,6,8)

             #.. number of periods to project forward over (e.g.,
             #     number of five-year steps)
             ,proj.periods = 8

             #.. age group width
             ,age.size = 5

             #.. rows in fertility rate matrix corresponding to ages with
             #     non-zero fertility
             ,fert.rows = 4:10

             #.. print algorithm progress
             ,verb = FALSE

             #.. tolerance defining allowable survival probabilities
             ,s.tol = 10^(-10)

             #.. **variances** for proposal distributions used in M-H
             #   steps which update vital rates.
             ,prop.vars

             )

    ##
    ## CREATOR
    ##   Mark C. Wheldon
    ##
    ## REFERENCE
    ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
    ##   P. (2013). Reconstructing Past Populations with Uncertainty
    ##   from Fragmentary Data. Journal of the American Statistical
    ##   Association, 108(501), 96–110.
    ##   http://doi.org/10.1080/01621459.2012.737729
    ##
    ## LICENCE
    ##   Released under the Creative Commons BY-NC-SA Licence
    ##   (https://creativecommons.org).
    ##
{

    ## -------- Begin timing ------- ##

    ptm <- proc.time()


    ## --------- Libraries --------- ##

    require(coda)


    ## ------ Match functions ------ ##

    #.. CCMP function
    ccmp.function <- match.fun(ccmp.f)
    mig.string <- "prop"


    ## ---------- Storage ---------- ##

    #.. Matrices for current and proposed values

    log.curr.f <- matrix(NA, nrow = nrow(init.f), ncol = ncol(init.f))
    logit.curr.s <- matrix(NA, nrow = nrow(init.s), ncol = ncol(init.s))
    curr.g <- matrix(NA, nrow = nrow(init.g), ncol = ncol(init.g))
    log.curr.b <- vector("numeric", length = length(init.b))

    log.prop.f <- matrix(NA, nrow = nrow(init.f), ncol = ncol(init.f))
    logit.prop.s <- matrix(NA, nrow = nrow(init.s), ncol = ncol(init.s))
    prop.g <- matrix(NA, nrow = nrow(init.g), ncol = ncol(init.g))
    log.prop.b <- vector("numeric", length = length(init.b))

    curr.sigmasq.f <- NA
    curr.sigmasq.s <- NA
    curr.sigmasq.g <- NA
    curr.sigmasq.n <- NA

    prop.sigmasq.f <- NA
    prop.sigmasq.s <- NA
    prop.sigmasq.g <- NA
    prop.sigmasq.n <- NA

    # current version of ccmp outputs a matrix with baseline
    #   in first column, so add extra column to these matrices
    log.curr.proj <- matrix(NA, nrow = nrow(init.b),
                     ncol = proj.periods + 1)
    log.prop.proj <- matrix(NA, nrow = nrow(init.b),
                     ncol = proj.periods + 1)


    #.. Matrices to store increments used make the proposals

    log.prop.f.mat <-
        matrix(0, nrow = nrow(init.f), ncol = ncol(init.f))
    logit.prop.s.mat <-
        matrix(0, nrow = nrow(init.s), ncol = ncol(init.s))
    prop.g.mat <-
        matrix(0, nrow = nrow(init.g), ncol = ncol(init.g))
    log.prop.b.mat <- matrix(0, nrow = nrow(init.b), ncol = 1)


    #.. MCMC objects for posterior samples
    # Samples are stored as 2D arrays for compatibility with coda's
    # mcmc format with iterations as rows, year*age.group as columns.
    # Age.group cycles fastest across columns, e.g.,
    # _____________________________________________________
    #   1960  | 1960  | 1960  | ... | 1965  | 1965  | ...
    #   15.19 | 20.24 | 25.29 | ... | 15.19 | 20.24 | ...
    # 1  --   |  --   |  --   | ... |  --   |  --   | ...
    # 2  --   |  --   |  --   | ... |  --   |  --   | ...
    #   etc.
    # _____________________________________________________

      # Fertility
      fert.rate.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = length(fert.rows) * ncol(init.f)))
      colnames(fert.rate.mcmc) <-
          popRecon.makeColNames(init.f[fert.rows,])

      # Survival proportions
      surv.prop.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.s) * ncol(init.s)))
      colnames(surv.prop.mcmc) <-
          popRecon.makeColNames(init.s)

      # lx
      lx.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.b) * (proj.periods)))
      colnames(lx.mcmc) <-
          popRecon.makeColNames(matrix(0, nrow = nrow(init.b)
                                 ,ncol = proj.periods
          ,dimnames = list(rownames(init.b)
           ,seq(from = as.numeric(colnames(init.b)[1]) +
           age.size, by = age.size, length = proj.periods))
                                 )
                          )

      # migration proportions
      mig.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.g) * ncol(init.g)))
      colnames(mig.mcmc) <-
          popRecon.makeColNames(init.g)

      # migration counts
      mig.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.g) * ncol(init.g)))
      colnames(mig.mcmc) <-
          popRecon.makeColNames(init.g)

      # baseline counts
      baseline.count.mcmc <-
          mcmc(matrix(nrow = n.iter, ncol = nrow(init.b)))
      colnames(baseline.count.mcmc) <- popRecon.makeColNames(init.b)

      # variances
      variances.mcmc <-
          mcmc(matrix(nrow = n.iter, ncol = 4))
      colnames(variances.mcmc) <-
          c("fert.rate.var", "surv.prop.var", "mig.var"
            ,"population.count.var")


    #.. Record acceptance rate

    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.g = 0
             ,sigmasq.n = 0
             )


    #.. Count how often acceptance ratio missing or na

    ar.na <- acc.count


    #.. Count how often projection gives negative population

    pop.negative <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             )


    #.. Count how often surv probs are outside tolerance

    s.out.tol <- matrix(0, nrow = nrow(mean.s), ncol = ncol(mean.s)
                        ,dimnames = dimnames(mean.s))


    ## -------- Initialize -------- ##

    #.. Set current vitals and variances to inital values
    #   Take logs/logits here where required

    log.curr.f <- log(init.f) #<-- log(0) stored as "-Inf". Gets
    log.prop.f <- log(init.f) #    converted to 0 under exponentiation
    logit.curr.s <- popReconAux.logit(init.s)
    curr.g <- init.g
    log.curr.b <- log(init.b)

    curr.sigmasq.f <- init.sigmasq.f
    curr.sigmasq.s <- init.sigmasq.s
    curr.sigmasq.g <- init.sigmasq.g
    curr.sigmasq.n <- init.sigmasq.n


    #.. Fixed means for vitals and baseline
    #   Set these to inputs, take logs where required.

    log.mean.f <- log(mean.f)
    logit.mean.s <- popReconAux.logit(mean.s)
    mean.g <- mean.g
    log.mean.b <- log(mean.b)


    #.. Fixed census data
    #   Take logs here

    log.census.mat <- log(pop.data)


    #.. Set current projection: base on initial values

    log.curr.proj <-
      log(
          ccmp.function(pop = exp(log.curr.b),
                            surv = popReconAux.invlogit(logit.curr.s),
                            fert = exp(log.curr.f),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size)
          )


    #.. Current log posterior

    log.curr.posterior <-
        popRecon.log.post(f = log.curr.f
                       ,s = logit.curr.s
                       ,g = curr.g
                       ,baseline.n = log.curr.b
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.g = mean.g
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.g = al.g, beta.g = be.g
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.g = curr.sigmasq.g
                       ,sigmasq.n = curr.sigmasq.n
                       ,log.like = popRecon.log.lhood(
                        log.n.census = log.census.mat
                        ,log.n.hat = log.curr.proj
                        ,ll.var = curr.sigmasq.n
                        ,cen.col = census.columns
                        )
                       ,non.zero.fert = fert.rows
                       )


    ## -------- Begin loop ------- ##
    #...............................#

    if(verb) {
        cat("\n\ntotal iterations = ", n.iter+burn.in
            ,"\nburn in = ", burn.in
            ,", stored = ", n.iter, "\n\n"
            ,"iter ", " quantity\n", "---- ", " --------"
            ,sep = "")
          }

    for(i in 1:(n.iter + burn.in)) {

      # Burn.in is not stored so create index for assignment later
      k <- i - burn.in


      ## -------- Vital Rate M-H Steps ------- ##

      ##...... Fertility .....##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Fertility")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(log.curr.f[fert.rows,])) {

        #.. make a matrix conformable w fertitlity rate matrix
        log.prop.f.mat <-
            matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
        log.prop.f.mat[fert.rows,][j] <-
            rnorm(1, 0, sqrt(prop.vars$fert.rate[j]))

        #.. make proposal
        log.prop.f <- log.curr.f + log.prop.f.mat

        # - Run CCMP (project on the original scale)
        #   ** Don't allow negative population
        prop.proj <-
          ccmp.function(pop = exp(log.curr.b),
                            fert = exp(log.prop.f), #<-- use proposal
                            surv = popReconAux.invlogit(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size,
                            )

        if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
           || is.nan(sum(prop.proj))) {
          if(i > burn.in) {
            pop.negative$fert.rate[j] <-
                pop.negative$fert.rate[j] + 1/n.iter
          }
      } else {
          log.prop.proj <- log(prop.proj)

          # - Calculate log posterior of proposed vital under projection
          log.prop.posterior <-
              popRecon.log.post(f = log.prop.f #<-- use proposal
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

          #- Acceptance ratio
          ar <- popRecon.acc.ra(log.prop = log.prop.posterior,
                             log.current = log.curr.posterior)

          # - Move or stay
          #.. stay if acceptance ratio 0, missing, infinity, etc.
          if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$fert.rate[j] <-
                ar.na$fert.rate[j] + 1/n.iter
          } else {
            #.. if accept, update current fert rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
              if(i > burn.in) acc.count$fert.rate[j] <-
                  acc.count$fert.rate[j] + 1/n.iter
              log.curr.f <- log.prop.f
              log.curr.proj <- log.prop.proj
              log.curr.posterior <- log.prop.posterior
            }
            #.. if reject, leave current fert rates and projections
            #   alone

          } # close else after checking for ar=0, missing, inf

        } # close else after checking neg or zero population

      } # close loop over all age-spec fertility rates

      #.. Store proposed fertility rate matrix
      if(i > burn.in) fert.rate.mcmc[k,] <-
          as.vector(exp(log.curr.f[fert.rows,]))


      ##...... Survival ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Survival")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.s)) {

        #.. make a matrix conformable w rate matrix
        logit.prop.s.mat <-
            matrix(0, nrow = nrow(logit.curr.s)
                   ,ncol = ncol(logit.curr.s))
        logit.prop.s.mat[j] <- rnorm(1, 0, sqrt(prop.vars$surv.prop[j]))

        #.. make proposal
        logit.prop.s <- logit.curr.s + logit.prop.s.mat

        #.. If proposal resulted in back-transformed s = 0 or 1, do
        #   nothing
        if(popReconAux.invlogit(logit.prop.s[j]) > 1 - s.tol ||
           popReconAux.invlogit(logit.prop.s[j]) < s.tol) {
          #.. leave current surv rates and projections
          #   alone (simply do not propose
          #   extreme survival probabilities)
          s.out.tol[j] <- s.out.tol[j] + 1/n.iter
        } else {

          # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
          prop.proj <-
            ccmp.function(pop = exp(log.curr.b),
                              fert = exp(log.curr.f),
                              surv = popReconAux.invlogit(logit.prop.s), #<-- use prop
                              mig = curr.g,
                              proj.steps = proj.periods,
                              age.int = age.size,
                              )

          if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
             || is.nan(sum(prop.proj))) {
            if(i > burn.in) {
              pop.negative$surv.prop[j] <-
                  pop.negative$surv.prop[j] + 1/n.iter
            }
          } else {
            log.prop.proj <- log(prop.proj)

            # - Calculate log posterior of proposed vital under projection
            log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.prop.s #<-- use proposal
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

            #- Acceptance ratio
            ar <- popRecon.acc.ra(log.prop = log.prop.posterior,
                               log.current = log.curr.posterior)

            # - Move or stay
            #.. stay if acceptance ratio 0, missing, infinity, etc.
            if(is.na(ar) || is.nan(ar) || ar < 0) {
              if(i > burn.in) ar.na$surv.prop[j] <-
                  ar.na$surv.prop[j] + 1/n.iter
            } else {
              #.. if accept, update current surv rates,
              #   update current projection and count acceptance
              if(runif(1) <= ar) {
                if(i > burn.in) acc.count$surv.prop[j] <-
                    acc.count$surv.prop[j] + 1/n.iter
                logit.curr.s <- logit.prop.s
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
              } #.. if reject, leave current surv rates and projections
                #   alone

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(i > burn.in) surv.prop.mcmc[k,] <-
        as.vector(popReconAux.invlogit(logit.curr.s))


      ##...... Migration ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Migration")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(curr.g)) {

        #.. make a matrix conformable w rate matrix
        prop.g.mat <-
            matrix(0, nrow = nrow(curr.g), ncol = ncol(curr.g))
        prop.g.mat[j] <- rnorm(1, 0, sqrt(prop.vars$mig[j]))

        #.. make proposal
        prop.g <- curr.g + prop.g.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
        prop.proj <-
            ccmp.function(pop = exp(log.curr.b),
                              fert = exp(log.curr.f),
                              surv = popReconAux.invlogit(logit.curr.s),
                              mig = prop.g, #<-- use proposal
                              proj.steps = proj.periods,
                              age.int = age.size,
                              )

      if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
         || is.nan(sum(prop.proj))) {
        if(i > burn.in) {
          pop.negative$mig[j] <-
              pop.negative$mig[j] + 1/n.iter
        }
      } else {
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = prop.g #<-- use proposal
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- popRecon.acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$mig[j] <-
                ar.na$mig[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$mig[j] <-
                    acc.count$mig[j] + 1/n.iter
                curr.g <- prop.g
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

    } # close loop over all age-specific migration proportions

      #.. Store proposed migration proportion matrix
      if(i > burn.in) mig.mcmc[k,] <- as.vector(curr.g)


      ##...... Baseline population ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

      # - Proposal

      #.. cycle through components (never update last
      #   value as this affects years beyond the estimation period)
      for(j in 1:length(log.curr.b)) {

      #.. make a matrix conformable w rate matrix
      log.prop.b.mat <- matrix(0, nrow = nrow(log.curr.b), ncol = 1)
      log.prop.b.mat[j] <- rnorm(1, 0, sqrt(prop.vars$population.count[j]))

      #.. make proposal
      log.prop.b <- log.curr.b + log.prop.b.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
      prop.proj <-
          ccmp.function(pop = exp(log.prop.b), #<-- use proposal
                            fert = exp(log.curr.f),
                            surv = popReconAux.invlogit(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size,
                            )

      if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
         || is.nan(sum(prop.proj))) {
        if(i > burn.in) {
          pop.negative$baseline.count[j] <-
              pop.negative$baseline.count[j] + 1/n.iter
        }
      } else {
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.prop.b #<-- use proposal
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- popRecon.acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$baseline.count[j] <-
                ar.na$baseline.count[j] + 1/n.iter
        } else {
            #.. if accept, update current mig rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$baseline.count[j] <-
                    acc.count$baseline.count[j] + 1/n.iter
                log.curr.b <- log.prop.b
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

  } # close loop over all age-specific baseline counts

      #.. Store proposed baseline count matrix
      if(i > burn.in) baseline.count.mcmc[k,] <-
          as.vector(exp(log.curr.b))


      ## ------- Variance Updates ------- ##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Variances")

      ##...... Fertility rate ......##

      prop.sigmasq.f <-
        popReconAux.rinvGamma(1, al.f +
                         length(mean.f[fert.rows,])/2,
                  be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                  log.mean.f[fert.rows,])^2)
                  )

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popRecon.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.f <-
                ar.na$sigmasq.f + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.f <-
                    acc.count$sigmasq.f + 1/n.iter
                curr.sigmasq.f <- prop.sigmasq.f
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"fert.rate.var"] <- curr.sigmasq.f


      ##...... Survival Proportion ......##

      prop.sigmasq.s <-
        popReconAux.rinvGamma(1, al.s + length(mean.s)/2,
                  be.s +
                    0.5*sum((logit.curr.s - logit.mean.s)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = prop.sigmasq.s  #<-- use proposal
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popRecon.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.s <-
                ar.na$sigmasq.s + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.s <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.s <- prop.sigmasq.s
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"surv.prop.var"] <- curr.sigmasq.s


      ##...... Migration Proportion ......##

      prop.sigmasq.g <-
        popReconAux.rinvGamma(1, al.g + length(mean.g)/2,
                  be.g +
                    0.5*sum((curr.g - mean.g)^2))

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = prop.sigmasq.g #<-- use proposal
                             ,sigmasq.n = curr.sigmasq.n
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popRecon.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.g <-
                ar.na$sigmasq.g + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.g <-
                    acc.count$sigmasq.g + 1/n.iter
                curr.sigmasq.g <- prop.sigmasq.g
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"mig.var"] <- curr.sigmasq.g


      ##...... Population Count ......##

      prop.sigmasq.n <-
        popReconAux.rinvGamma(1, al.n + (length(mean.b) +
                                    length(log.census.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2) +
                  sum((log.census.mat - log.curr.proj[,census.columns])^2)
                  )
                         )

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popRecon.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                             ,log.like = popRecon.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = prop.sigmasq.n #<-- use proposal
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popRecon.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.census.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.census.mat - log.curr.proj[,census.columns])^2))
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.n
                              ,al.n + (length(mean.b) +
                                       length(log.census.mat))/2
                              ,be.n + 0.5 * (sum((log.curr.b - log.mean.b)^2) + sum((log.census.mat - log.curr.proj[,census.columns])^2))
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.n <-
                ar.na$sigmasq.n + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.n <-
                    acc.count$sigmasq.n + 1/n.iter
                curr.sigmasq.n <- prop.sigmasq.n
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) {
        variances.mcmc[k,"population.count.var"] <- curr.sigmasq.n
      }


      ## ------- Store current population and migration ------- ##

      lx.mcmc[k,] <-
          as.vector(exp(log.curr.proj))[-(1:ncol(baseline.count.mcmc))]


      if(verb && identical(i%%1000, 0)) cat("\n\n")

  } # Ends outer-most loop

    ## ......... End Loop ........ ##
    #...............................#


    ## ---------- Output --------- ##

    #cat("inital values", "\n\n")
    #.. initial values
    init.vals <- list(fert.rate = init.f
                      ,surv.prop = init.s
                      ,mig.XXX = init.g
                      ,baseline.count = init.b
                      ,init.sigmasq.f = init.sigmasq.f
                      ,init.sigmasq.s = init.sigmasq.s
                      ,init.sigmasq.g = init.sigmasq.g
                      ,init.sigmasq.n = init.sigmasq.n
                      ,pop.data = pop.data
                      )

    #.. fixed parameters
    fixed.params <- list(mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.mig.XXX = mean.g
                         ,mean.baseline.count = mean.b
                         ,pop.counts = pop.data
                         ,alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.mig.XXX = al.g
                         ,beta.mig.XXX = be.g
                         ,alpha.population.count = al.n
                         ,beta.population.count = be.n
                         )

    #.. migration type
    new.names <-
        lapply(list(init.vals, fixed.params), FUN = function(z) {
        sub("XXX", mig.string, names(z))
    })
    names(init.vals) <- new.names[[1]]
    names(fixed.params) <- new.names[[2]]


    #cat("algorithm statistics", "\n\n")
    #.. algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             ,acc.prop.adj4neg = mapply(FUN = function(a, b, n) {
                 (a * n) / (n - b)
             },
              acc.count[1:4], pop.negative, MoreArgs = list(n = n.iter)
              )
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    #cat("algorithm parameters", "\n\n")
    #.. algorithm parameters
    alg.params <- list(prop.vars = prop.vars
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit", mig.XXX = "I"
                        ,baseline.count = "log"
                        ,population.count = "log")
                       ,projection.periods = proj.periods
                       ,age.gp.size = age.size
                       ,cen.columns = census.columns
                       ,non.zero.fert.rows = fert.rows
                       ,surv.tolerance = s.tol
                       ,burn.in = burn.in
                       ,iters = n.iter
                       )
    names(alg.params) <- sub("XXX", mig.string, names(alg.params))

    #.. results

    ret.list <- list(fert.rate.mcmc = fert.rate.mcmc
                  ,surv.prop.mcmc = surv.prop.mcmc
                  ,mig.XXX.mcmc = mig.mcmc
                  ,baseline.count.mcmc = baseline.count.mcmc
                  ,lx.mcmc = lx.mcmc
                  ,variances.mcmc = variances.mcmc
                  ,alg.stats = alg.stats
                  ,fixed.params = fixed.params
                  ,init.vals = init.vals
                  ,alg.params = alg.params
                  )
    names(ret.list) <- sub("XXX", mig.string, names(ret.list))

    return(ret.list)

  }


###
### ** t2 likelihood SAMPLER
###
################################################################################

###
### *** Secondary Functions
###
################################################################################

###
### Likelihood
###

popReconModCheck.log.lhood <-
    function(log.n.census, log.n.hat, ll.var, ll.lambda, cen.col)
  {
    #.. log.n.census and log.n.hat should already be logged


    #-- Choose years for which likelihood exists --#
    #.. passed in from sampler function

    log.n.hat.select <- log.n.hat[,cen.col]


    #-- value of log likelihoods --#

    if(any(is.na(ll.lambda))) ll.lambda <- 1

    density <- dnorm(log.n.census,
                     mean = log.n.hat.select,
                     sd = sqrt(ll.var * ll.lambda),
                     log = TRUE
                     )

    #-- joint log likelihood --#

    return(sum(density))
  }


###
### Posterior
###

popReconModCheck.log.post <- function(# estimated vitals
                     f, s, g, baseline.n
                     # fixed prior means on vitals
                     ,prior.mean.f, prior.mean.s
                     ,prior.mean.g, prior.mean.b
                     # fixed prior parameters on variance distns
                     ,alpha.f, beta.f, alpha.s, beta.s
                     ,alpha.g, beta.g
                     ,alpha.n, beta.n
                     # updated variances on prior distns
                     ,sigmasq.f, sigmasq.s, sigmasq.g, sigmasq.n
                     # lambdas
                     ,lambda.f, lambda.s, lambda.g, lambda.n
                     # value of the log likelihood
                     ,log.like
                     # non zero rows of fertility matrix
                     ,non.zero.fert
                     )
  {

      #-- Values of prior densities for vitals --#

      #.. Note that log densities are calculated for numerical stability.
      #     f, baseline.n, prior.mean.f, prior.mean.b are logged coming
      #     in, s, prior.mean.s is logit transformed coming in, g and
      #     prior.mean.g are not transformed coming in.

      if(all(!is.na(lambda.f))) log.f.prior <-
          dnorm(as.vector(f[non.zero.fert,])
                           ,mean = as.vector(prior.mean.f[non.zero.fert,])
                           ,sd = sqrt(sigmasq.f * lambda.f)
                           ,log = TRUE)
      else log.f.prior <- dnorm(as.vector(f[non.zero.fert,])
                           ,mean = as.vector(prior.mean.f[non.zero.fert,])
                           ,sd = sqrt(sigmasq.f)
                           ,log = TRUE)

      if(all(!is.na(lambda.s))) log.s.prior <-
          dnorm(s, mean = prior.mean.s, sd = sqrt(sigmasq.s * lambda.s)
                           ,log = TRUE)
      else log.s.prior <- dnorm(s, mean = prior.mean.s
                           ,sd = sqrt(sigmasq.s)
                           ,log = TRUE)

      if(all(!is.na(lambda.g))) log.g.prior <-
          dnorm(g, mean = prior.mean.g
                 ,sd = sqrt(sigmasq.g * lambda.g)
                 ,log = TRUE)
      else log.g.prior <- dnorm(g, mean = prior.mean.g
                                 ,sd = sqrt(sigmasq.g)
                                 ,log = TRUE)

      if(all(!is.na(lambda.n))) log.b.prior <-
          dnorm(baseline.n, mean = prior.mean.b
                           ,sd = sqrt(sigmasq.n * lambda.n[,1])
                           ,log = TRUE)
      else log.b.prior <- dnorm(baseline.n, mean = prior.mean.b
                           ,sd = sqrt(sigmasq.n)
                           ,log = TRUE)

      #-- Values of prior densities for variances --#

      log.sigmasq.f.prior <-
          log(popReconAux.dinvGamma(sigmasq.f, alpha.f, beta.f))
      log.sigmasq.s.prior <-
          log(popReconAux.dinvGamma(sigmasq.s, alpha.s, beta.s))
      log.sigmasq.g.prior <-
          log(popReconAux.dinvGamma(sigmasq.g, alpha.g, beta.g))
      log.sigmasq.n.prior <-
          log(popReconAux.dinvGamma(sigmasq.n, alpha.n, beta.n))

      #-- Values of prior densities for lambdas --#

      if(all(!is.na(lambda.f))) log.lambda.f.prior <-
          log(popReconAux.dinvGamma(lambda.f, 1, 1))
      else log.lambda.f.prior <- 0
      if(all(!is.na(lambda.s))) log.lambda.s.prior <-
          log(popReconAux.dinvGamma(lambda.s, 1, 1))
      else log.lambda.s.prior <- 0
      if(all(!is.na(lambda.g))) log.lambda.g.prior <-
          log(popReconAux.dinvGamma(lambda.g, 1, 1))
      else log.lambda.g.prior <- 0
      if(all(!is.na(lambda.n))) log.lambda.n.prior <-
          log(popReconAux.dinvGamma(lambda.n, 1, 1))
      else log.lambda.n.prior <- 0


      #-- The log posterior is the SUM of these with the log.like --#

      return(sum(log.f.prior, log.s.prior, log.g.prior, log.b.prior
                 ,log.sigmasq.f.prior
                 ,log.sigmasq.s.prior
                 ,log.sigmasq.g.prior
                 ,log.sigmasq.n.prior
                 ,log.lambda.f.prior
                 ,log.lambda.s.prior
                 ,log.lambda.g.prior
                 ,log.lambda.n.prior
                 ,log.like))

  }


###
### Acceptance Ratio
###

popReconModCheck.acc.ra <- function(log.prop, log.current)
  {
    min(1, exp(log.prop - log.current))
  }

popReconModCheck.acc.ra.var <-
    function(log.prop.post, log.curr.post, log.prop.var, log.curr.var)
{
    min(1, exp(log.curr.var + log.prop.post - log.prop.var - log.curr.post
               ))
}

popReconModCheck.acc.ra.lambda <-
    function(log.prop.post, log.curr.post, log.prop.lambda
             ,log.curr.lambda)
{
    min(1, exp(log.curr.lambda + log.prop.post - log.prop.lambda -
               log.curr.post
               ))
}


###
### *** SAMPLER
###
################################################################################

popReconModCheck.sampler <-
    function(#.. number of iterations and burn-in (not saved)
             n.iter, burn.in = 0

             #.. ccmp function: list with components 'FUN' and
             #   'formals'. Use formals to set, eg, mig.type.
             ,ccmp.f

             #.. fixed variance hyper-parameters
             ,al.f, be.f, al.s, be.s, al.g, be.g, al.n, be.n

             #.. fixed prior means
             ,mean.f, mean.s, mean.g, mean.b

             #.. inital values for vitals and variances
             #   *vitals not transformed coming in*
             ,init.f, init.s, init.g, init.b
             ,init.sigmasq.f, init.sigmasq.s, init.sigmasq.g
             ,init.sigmasq.n

             #.. lambdas to include
             ,lambda.n = NA

             #.. census data: use 'census.columns' to match this
             #   to columns in the projection returned by ccmp
             #   *not transformed coming in*
             ,pop.data

             #.. periods in projection range for which census data exists
             ,census.columns = c(4,6,8)

             #.. number of periods to project forward over (e.g.,
             #     number of five-year steps)
             ,proj.periods = 8

             #.. age group width
             ,age.size = 5

             #.. rows in fertility rate matrix corresponding to ages with
             #     non-zero fertility
             ,fert.rows = 4:10

             #.. print algorithm progress
             ,verb = FALSE

             #.. tolerance defining allowable survival probabilities
             ,s.tol = 10^(-10)

             #.. **variances** for proposal distributions used in M-H
             #   steps which update vital rates.
             ,prop.vars

             #.. save traces of log likelihood and log posterior
             #   densities?

             ,save.log.lhood = FALSE
             ,save.log.post = FALSE

             )

    ##
    ## CREATOR
    ##   Mark C. Wheldon
    ##
    ## REFERENCE
    ##   Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
    ##   P. (2013). Reconstructing Past Populations with Uncertainty
    ##   from Fragmentary Data. Journal of the American Statistical
    ##   Association, 108(501), 96–110.
    ##   http://doi.org/10.1080/01621459.2012.737729
    ##
    ## LICENCE
    ##   Released under the Creative Commons BY-NC-SA Licence
    ##   (https://creativecommons.org).
    ##
{

    ## -------- Begin timing ------- ##

    ptm <- proc.time()


    ## --------- Libraries --------- ##
    require(coda)


    ## ------ Match functions ------ ##

    #.. CCMP function
    ccmp.function <- match.fun(ccmp.f)
    mig.string <- "prop"


    ## ---------- Storage ---------- ##

    #.. Matrices for current and proposed values
    #.. Set current vitals and variances to inital values
    #   Take logs/logits here where required

    log.curr.f <- log(init.f) #<-- log(0) stored as "-Inf". Gets
                              #     converted to 0 under exponentiation
    logit.curr.s <- popReconAux.logit(init.s)
    curr.g <- init.g
    log.curr.b <- log(init.b)

    log.prop.f <- log.curr.f
    log.prop.s <- logit.curr.s
    log.prop.g <- curr.g
    log.prop.b <- log.curr.b

    curr.sigmasq.f <- init.sigmasq.f
    curr.sigmasq.s <- init.sigmasq.s
    curr.sigmasq.g <- init.sigmasq.g
    curr.sigmasq.n <- init.sigmasq.n

    prop.sigmasq.f <- curr.sigmasq.f
    prop.sigmasq.s <- curr.sigmasq.s
    prop.sigmasq.g <- curr.sigmasq.g
    prop.sigmasq.n <- curr.sigmasq.n


    #.. Fixed means for vitals and baseline
    #   Set these to inputs, take logs where required.

    log.mean.f <- log(mean.f)
    logit.mean.s <- popReconAux.logit(mean.s)
    mean.g <- mean.g
    log.mean.b <- log(mean.b)


    #.. Fixed census data
    #   Take logs here

    log.census.mat <- log(pop.data)


    #.. Mixing parameter

    if(any(is.na(lambda.n))) {
        curr.lambda.n <- matrix(NA, ncol = 2) # Keep
                                              # 'curr.lambda.n[,-1]'
                                              # happy
    } else {
        curr.lambda.n <- lambda.n
        prop.lambda.n <- curr.lambda.n
    }

    ## other mixing not implemented
    lambda.f <- NA
    lambda.s <- NA
    lambda.g <- NA
    curr.lambda.f <- NA
    curr.lambda.s <- NA
    curr.lambda.g <- NA


    #.. Set current projection: base on initial values

    log.curr.proj <-
      log(
          ccmp.function(pop = exp(log.curr.b),
                            surv = popReconAux.invlogit(logit.curr.s),
                            fert = exp(log.curr.f),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size)
          )

    # For lambda.n

    log.mean.lx <- log(cbind(mean.b, pop.data))


    #.. Current log posterior

    log.curr.posterior <-
        popReconModCheck.log.post(f = log.curr.f
                       ,s = logit.curr.s
                       ,g = curr.g
                       ,baseline.n = log.curr.b
                       ,prior.mean.f = log.mean.f
                       ,prior.mean.s = logit.mean.s
                       ,prior.mean.g = mean.g
                       ,prior.mean.b = log.mean.b
                       ,alpha.f = al.f, beta.f = be.f
                       ,alpha.s = al.s, beta.s = be.s
                       ,alpha.g = al.g, beta.g = be.g
                       ,alpha.n = al.n, beta.n = be.n
                       ,sigmasq.f = curr.sigmasq.f
                       ,sigmasq.s = curr.sigmasq.s
                       ,sigmasq.g = curr.sigmasq.g
                       ,sigmasq.n = curr.sigmasq.n
                       ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                       ,log.like = popReconModCheck.log.lhood(
                        log.n.census = log.census.mat
                        ,log.n.hat = log.curr.proj
                        ,ll.var = curr.sigmasq.n
                        ,ll.lambda = curr.lambda.n[,-1]
                        ,cen.col = census.columns
                        )
                       ,non.zero.fert = fert.rows
                       )


    #.. Matrices to store increments used make the proposals

    log.prop.f.mat <-
        matrix(0, nrow = nrow(init.f), ncol = ncol(init.f))
    logit.prop.s.mat <-
        matrix(0, nrow = nrow(init.s), ncol = ncol(init.s))
    prop.g.mat <-
        matrix(0, nrow = nrow(init.g), ncol = ncol(init.g))
    log.prop.b.mat <- matrix(0, nrow = nrow(init.b), ncol = 1)


    #.. MCMC objects for posterior samples (OUTPUTS)

    # Samples are stored as 2D arrays for compatibility with coda's
    # mcmc format with iterations as rows, year*age.group as columns.
    # Age.group cycles fastest across columns, e.g.,
    # _____________________________________________________
    #   1960  | 1960  | 1960  | ... | 1965  | 1965  | ...
    #   15.19 | 20.24 | 25.29 | ... | 15.19 | 20.24 | ...
    # 1  --   |  --   |  --   | ... |  --   |  --   | ...
    # 2  --   |  --   |  --   | ... |  --   |  --   | ...
    #   etc.
    # _____________________________________________________

      #.. store in memory

      # Fertility
      fert.rate.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = length(fert.rows) * ncol(init.f)))
      colnames(fert.rate.mcmc) <-
          popRecon.makeColNames(init.f[fert.rows,])

      # Survival proportions
      surv.prop.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.s) * ncol(init.s)))
      colnames(surv.prop.mcmc) <-
          popRecon.makeColNames(init.s)

      # lx
      lx.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.b) * (proj.periods)))
      colnames(lx.mcmc) <-
          popRecon.makeColNames(matrix(0, nrow = nrow(init.b)
                                 ,ncol = proj.periods
          ,dimnames = list(rownames(init.b)
           ,seq(from = as.numeric(colnames(init.b)[1]) +
           age.size, by = age.size, length = proj.periods))
                                 )
                          )

      # migration proportions
      mig.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.g) * ncol(init.g)))
      colnames(mig.mcmc) <-
          popRecon.makeColNames(init.g)

      # migration counts
      mig.mcmc <-
          mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(init.g) * ncol(init.g)))
      colnames(mig.mcmc) <-
          popRecon.makeColNames(init.g)

      # baseline counts
      baseline.count.mcmc <-
          mcmc(matrix(nrow = n.iter, ncol = nrow(init.b)))
      colnames(baseline.count.mcmc) <- popRecon.makeColNames(init.b)

      # variances
      variances.mcmc <-
          mcmc(matrix(nrow = n.iter, ncol = 4))
      colnames(variances.mcmc) <-
          c("fert.rate.var", "surv.prop.var", "mig.var"
            ,"population.count.var")


      # mixing parameters

      if(all(!is.na(lambda.f)))
          lambda.f.mcmc <- fert.rate.mcmc

      if(all(!is.na(lambda.s)))
          lambda.s.mcmc <- surv.prop.mcmc

      if(all(!is.na(lambda.g)))
          lambda.g.mcmc <- mig.mcmc

      if(all(!is.na(lambda.n))) {
      lambda.n.mcmc <- mcmc(matrix(nrow = n.iter
                      ,ncol = nrow(pop.data) * (ncol(pop.data) + 1)
                                   ))
      colnames(lambda.n.mcmc) <-
          popRecon.makeColNames(cbind(init.b, pop.data))
  }

    #.. Record acceptance rate

    acc.count <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             ,sigmasq.f = 0
             ,sigmasq.s = 0
             ,sigmasq.g = 0
             ,sigmasq.n = 0
             )
    if(all(!is.na(lambda.f)))
        acc.count$lambda.f <- acc.count$fert.rate
    if(all(!is.na(lambda.s)))
        acc.count$lambda.s <- acc.count$surv.prop
    if(all(!is.na(lambda.g)))
        acc.count$lambda.g <- acc.count$mig
    if(all(!is.na(lambda.n)))
        acc.count$lambda.n <- matrix(0, nrow = nrow(mean.b)
             ,ncol = ncol(pop.data) + 1
             ,dimnames = list(rownames(pop.data)
              ,c(colnames(mean.b), colnames(pop.data))
             )
                                     )


    #.. Count how often acceptance ratio missing or na

    ar.na <- acc.count


    #.. Count how often projection gives negative population

    pop.negative <-
        list(fert.rate = matrix(0, nrow = nrow(mean.f[fert.rows,])
             ,ncol = ncol(mean.f[fert.rows,])
             ,dimnames = dimnames(mean.f[fert.rows,])
             )
             ,surv.prop = matrix(0, nrow = nrow(mean.s)
              ,ncol = ncol(mean.s)
              ,dimnames = dimnames(mean.s)
              )
             ,mig = matrix(0, nrow = nrow(mean.g), ncol = ncol(mean.g)
              ,dimnames = dimnames(mean.g)
              )
             ,baseline.count = matrix(0, nrow = nrow(mean.b)
              ,dimnames = dimnames(mean.b)
              )
             )


    #.. Count how often surv probs are outside tolerance

    s.out.tol <- matrix(0, nrow = nrow(mean.s), ncol = ncol(mean.s)
                        ,dimnames = dimnames(mean.s))


    #.. Save log likelihood

    if(save.log.lhood) log.lhood.mcmc <-
        mcmc(matrix(NA, nrow = n.iter, ncol = 1
                    ,dimnames = list(NULL, "log.lhood")
                    ))


    #.. Save log posterior density

    if(save.log.post) log.posterior.mcmc <-
        mcmc(matrix(NA, nrow = n.iter, ncol = 1
                    ,dimnames = list(NULL, "log.posterior")
                    ))



    ## -------- Begin loop ------- ##
    #...............................#

    if(verb) {
        cat("\n\ntotal iterations = ", n.iter+burn.in
            ,"\nburn in = ", burn.in
            ,", stored = ", n.iter, "\n\n"
            ,"iter ", " quantity\n", "---- ", " --------"
            ,sep = "")
          }

    for(i in 1:(n.iter + burn.in)) {

      # Burn.in is not stored so create index for assignment later
      k <- i - burn.in


      ## -------- Vital Rate M-H Steps ------- ##

      ##...... Fertility .....##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Fertility")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(log.curr.f[fert.rows,])) {

        #.. make a matrix conformable w fertitlity rate matrix
        log.prop.f.mat <-
            matrix(0, nrow = nrow(log.curr.f), ncol = ncol(log.curr.f))
        log.prop.f.mat[fert.rows,][j] <-
            rnorm(1, 0, sqrt(prop.vars$fert.rate[j]))

        #.. make proposal
        log.prop.f <- log.curr.f + log.prop.f.mat

        # - Run CCMP (project on the original scale)
        #   ** Don't allow negative population
        prop.proj <-
          ccmp.function(pop = exp(log.curr.b),
                            fert = exp(log.prop.f), #<-- use proposal
                            surv = popReconAux.invlogit(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size,
                            )

        if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
           || is.nan(sum(prop.proj))) {
          if(i > burn.in) {
            pop.negative$fert.rate[j] <-
                pop.negative$fert.rate[j] + 1/n.iter
          }
      } else {
          log.prop.proj <- log(prop.proj)

          # - Calculate log posterior of proposed vital under projection
          log.prop.posterior <-
              popReconModCheck.log.post(f = log.prop.f #<-- use proposal
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

          #- Acceptance ratio
          ar <- popReconModCheck.acc.ra(log.prop = log.prop.posterior,
                             log.current = log.curr.posterior)

          # - Move or stay
          #.. stay if acceptance ratio 0, missing, infinity, etc.
          if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$fert.rate[j] <-
                ar.na$fert.rate[j] + 1/n.iter
          } else {
            #.. if accept, update current fert rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
              if(i > burn.in) acc.count$fert.rate[j] <-
                  acc.count$fert.rate[j] + 1/n.iter
              log.curr.f <- log.prop.f
              log.curr.proj <- log.prop.proj
              log.curr.posterior <- log.prop.posterior
            }
            #.. if reject, leave current fert rates and projections
            #   alone

          } # close else after checking for ar=0, missing, inf

        } # close else after checking neg or zero population

      } # close loop over all age-spec fertility rates

      #.. Store proposed fertility rate matrix
      if(i > burn.in) fert.rate.mcmc[k,] <-
          as.vector(exp(log.curr.f[fert.rows,]))


      ##...... Survival ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Survival")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(logit.curr.s)) {

        #.. make a matrix conformable w rate matrix
        logit.prop.s.mat <-
            matrix(0, nrow = nrow(logit.curr.s)
                   ,ncol = ncol(logit.curr.s))
        logit.prop.s.mat[j] <- rnorm(1, 0, sqrt(prop.vars$surv.prop[j]))

        #.. make proposal
        logit.prop.s <- logit.curr.s + logit.prop.s.mat

        #.. If proposal resulted in back-transformed s = 0 or 1, do
        #   nothing
        if(popReconAux.invlogit(logit.prop.s[j]) > 1 - s.tol ||
           popReconAux.invlogit(logit.prop.s[j]) < s.tol) {
          #.. leave current surv rates and projections
          #   alone (simply do not propose
          #   extreme survival probabilities)
          s.out.tol[j] <- s.out.tol[j] + 1/n.iter
        } else {

          # - Run CCMP (project on the original scale)
          #   ** Don't allow negative population; again, simply treat
          #      this as if the proposal were never made
          prop.proj <-
            ccmp.function(pop = exp(log.curr.b),
                              fert = exp(log.curr.f),
                              surv = popReconAux.invlogit(logit.prop.s), #<-- use prop
                              mig = curr.g,
                              proj.steps = proj.periods,
                              age.int = age.size,
                              )

          if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
             || is.nan(sum(prop.proj))) {
            if(i > burn.in) {
              pop.negative$surv.prop[j] <-
                  pop.negative$surv.prop[j] + 1/n.iter
            }
          } else {
            log.prop.proj <- log(prop.proj)

            # - Calculate log posterior of proposed vital under projection
            log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.prop.s #<-- use proposal
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

            #- Acceptance ratio
            ar <- popReconModCheck.acc.ra(log.prop = log.prop.posterior,
                               log.current = log.curr.posterior)

            # - Move or stay
            #.. stay if acceptance ratio 0, missing, infinity, etc.
            if(is.na(ar) || is.nan(ar) || ar < 0) {
              if(i > burn.in) ar.na$surv.prop[j] <-
                  ar.na$surv.prop[j] + 1/n.iter
            } else {
              #.. if accept, update current surv rates,
              #   update current projection and count acceptance
              if(runif(1) <= ar) {
                if(i > burn.in) acc.count$surv.prop[j] <-
                    acc.count$surv.prop[j] + 1/n.iter
                logit.curr.s <- logit.prop.s
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
              } #.. if reject, leave current surv rates and projections
                #   alone

            } # close else{ after checking for undefined ar

          } # close else{ after checking for negative pop

        } # close else{ after checking for s outside tol

      } # close loop over all age-spec survival probabilities

      #.. Store proposed survival probability matrix
      if(i > burn.in) surv.prop.mcmc[k,] <-
        as.vector(popReconAux.invlogit(logit.curr.s))


      ##...... Migration ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Migration")

      # - Proposal

      #.. cycle through components
      for(j in 1:length(curr.g)) {

        #.. make a matrix conformable w rate matrix
        prop.g.mat <-
            matrix(0, nrow = nrow(curr.g), ncol = ncol(curr.g))
        prop.g.mat[j] <- rnorm(1, 0, sqrt(prop.vars$mig[j]))

        #.. make proposal
        prop.g <- curr.g + prop.g.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
        prop.proj <-
            ccmp.function(pop = exp(log.curr.b),
                              fert = exp(log.curr.f),
                              surv = popReconAux.invlogit(logit.curr.s),
                              mig = prop.g, #<-- use proposal
                              proj.steps = proj.periods,
                              age.int = age.size,
                              )

      if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
         || is.nan(sum(prop.proj))) {
        if(i > burn.in) {
          pop.negative$mig[j] <-
              pop.negative$mig[j] + 1/n.iter
        }
      } else {
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = prop.g #<-- use proposal
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- popReconModCheck.acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$mig[j] <-
                ar.na$mig[j] + 1/n.iter
        } else {
            #.. if accept, update current vital rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$mig[j] <-
                    acc.count$mig[j] + 1/n.iter
                curr.g <- prop.g
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

    } # close loop over all age-specific migration proportions

      #.. Store proposed migration proportion matrix
      if(i > burn.in) mig.mcmc[k,] <- as.vector(curr.g)


      ##...... Baseline population ......##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Baseline")

      # - Proposal

      #.. cycle through components (never update last
      #   value as this affects years beyond the estimation period)
      for(j in 1:length(log.curr.b)) {

      #.. make a matrix conformable w rate matrix
      log.prop.b.mat <- matrix(0, nrow = nrow(log.curr.b), ncol = 1)
      log.prop.b.mat[j] <- rnorm(1, 0, sqrt(prop.vars$population.count[j]))

      #.. make proposal
      log.prop.b <- log.curr.b + log.prop.b.mat

      # - Run CCMP (project on the original scale)
      #   ** Don't allow negative population
      prop.proj <-
          ccmp.function(pop = exp(log.prop.b), #<-- use proposal
                            fert = exp(log.curr.f),
                            surv = popReconAux.invlogit(logit.curr.s),
                            mig = curr.g,
                            proj.steps = proj.periods,
                            age.int = age.size,
                            )

      if(sum(prop.proj < 0) > 0 || is.na(sum(prop.proj))
         || is.nan(sum(prop.proj))) {
        if(i > burn.in) {
          pop.negative$baseline.count[j] <-
              pop.negative$baseline.count[j] + 1/n.iter
        }
      } else {
        log.prop.proj <- log(prop.proj)

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.prop.b #<-- use proposal
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.prop.proj #<-- use proposal
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

        #- Acceptance ratio
        ar <- popReconModCheck.acc.ra(log.prop = log.prop.posterior,
                           log.current = log.curr.posterior)

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$baseline.count[j] <-
                ar.na$baseline.count[j] + 1/n.iter
        } else {
            #.. if accept, update current mig rates, store proposed
            #   rate, update current projection and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$baseline.count[j] <-
                    acc.count$baseline.count[j] + 1/n.iter
                log.curr.b <- log.prop.b
                log.curr.proj <- log.prop.proj
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current fert rates and projections
            #   alone, store current rate

        } # close else after checking for ar=na, nan, zero

    } # close else after checking for negative population

  } # close loop over all age-specific baseline counts

      #.. Store proposed baseline count matrix
      if(i > burn.in) baseline.count.mcmc[k,] <-
          as.vector(exp(log.curr.b))


      ## ------- Variance Updates ------- ##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Variances")

      ##...... Fertility rate ......##

      #.. Make proposal

      if(all(!is.na(lambda.f))) {
          prop.sigmasq.f <-
        popReconAux.rinvGamma(1, al.f +
                         length(mean.f[fert.rows,])/2,
                  be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                  log.mean.f[fert.rows,])^2 /
                                 curr.lambda.f
                                 )
                  )
      } else {
          prop.sigmasq.f <-
        popReconAux.rinvGamma(1, al.f +
                         length(mean.f[fert.rows,])/2,
                  be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                  log.mean.f[fert.rows,])^2)
                  )
      }

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = prop.sigmasq.f #<-- use proposal
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popReconModCheck.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.f
                              ,al.f + length(mean.f[fert.rows,])/2
                              ,be.f + 0.5*sum((log.curr.f[fert.rows,] -
                                               log.mean.f[fert.rows,])^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.f <-
                ar.na$sigmasq.f + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.f <-
                    acc.count$sigmasq.f + 1/n.iter
                curr.sigmasq.f <- prop.sigmasq.f
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"fert.rate.var"] <- curr.sigmasq.f


      ##...... Survival Proportion ......##

      if(all(!is.na(lambda.s))) {
          prop.sigmasq.f <-
        popReconAux.rinvGamma(1, al.s + length(mean.s)/2,
                  be.s +
                    0.5*sum((logit.curr.s - logit.mean.s)^2) /
                                  curr.lambda.s
                                  )
      } else {
      prop.sigmasq.s <-
        popReconAux.rinvGamma(1, al.s + length(mean.s)/2,
                  be.s +
                    0.5*sum((logit.curr.s - logit.mean.s)^2))
  }

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = prop.sigmasq.s  #<-- use proposal
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popReconModCheck.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.s
                              ,al.s + length(mean.s)/2
                              ,be.s + 0.5*sum((logit.curr.s -
                                               logit.mean.s)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.s <-
                ar.na$sigmasq.s + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.s <-
                    acc.count$sigmasq.s + 1/n.iter
                curr.sigmasq.s <- prop.sigmasq.s
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"surv.prop.var"] <- curr.sigmasq.s


      ##...... Migration Proportion ......##

      if(all(!is.na(lambda.g))) {
      prop.sigmasq.g <-
        popReconAux.rinvGamma(1, al.g + length(mean.g)/2,
                  be.g +
                    0.5*sum((curr.g - mean.g)^2) /
                                  curr.lambda.g
                                  )

      } else {
      prop.sigmasq.g <-
        popReconAux.rinvGamma(1, al.g + length(mean.g)/2,
                  be.g +
                    0.5*sum((curr.g - mean.g)^2))
  }

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = prop.sigmasq.g #<-- use proposal
                             ,sigmasq.n = curr.sigmasq.n
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n #<-- use current
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      #- Acceptance ratio
      ar <- popReconModCheck.acc.ra.var(log.prop.post = log.prop.posterior
                             ,log.curr.post = log.curr.posterior
                             ,log.prop.var = popReconAux.dinvGamma(prop.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             ,log.curr.var = popReconAux.dinvGamma(curr.sigmasq.g
                              ,al.g + length(mean.g)/2
                              ,be.g + 0.5*sum((curr.g -
                                               mean.g)^2)
                              ,log = TRUE)
                             )

        # - Move or stay
        #.. stay if acceptance ratio 0, missing, infinity, etc.
        if(is.na(ar) || is.nan(ar) || ar < 0) {
            if(i > burn.in) ar.na$sigmasq.g <-
                ar.na$sigmasq.g + 1/n.iter
        } else {
            #.. if accept, update current, store proposed
            #   and count acceptance
            if(runif(1) <= ar) {
                if(i > burn.in) acc.count$sigmasq.g <-
                    acc.count$sigmasq.g + 1/n.iter
                curr.sigmasq.g <- prop.sigmasq.g
                log.curr.posterior <- log.prop.posterior
            } #.. if reject, leave current and posterior
        } # close else after checking for ar=na, nan, zero

      if(i > burn.in) variances.mcmc[k,"mig.var"] <- curr.sigmasq.g


      ##...... Population Count ......##

      if(all(!is.na(lambda.n))) {
          prop.sigmasq.n <-
        popReconAux.rinvGamma(1, al.n + (length(mean.b) +
                                    length(log.census.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2/curr.lambda.n[,1]) +
                  sum((log.census.mat - log.curr.proj[,census.columns])^2 /
                      curr.lambda.n[,-1])
                  )
                         )
      } else {
      prop.sigmasq.n <-
        popReconAux.rinvGamma(1, al.n + (length(mean.b) +
                                    length(log.census.mat))/2,
                be.n + 0.5 * (
                  sum((log.curr.b - log.mean.b)^2) +
                  sum((log.census.mat - log.curr.proj[,census.columns])^2)
                  )
                         )
  }

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = prop.sigmasq.n #<-- use proposal
                                 ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = curr.lambda.n
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = prop.sigmasq.n #<-- use proposal
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      if(i > burn.in) acc.count$sigmasq.n <-
                    acc.count$sigmasq.n + 1/n.iter
                curr.sigmasq.n <- prop.sigmasq.n
                log.curr.posterior <- log.prop.posterior

      if(i > burn.in) {
        variances.mcmc[k,"population.count.var"] <- curr.sigmasq.n
      }


      ## ------- Store current population --------------------- ##

      lx.mcmc[k,] <-
          as.vector(exp(log.curr.proj))[-(1:ncol(baseline.count.mcmc))]


      ## -------- Mixing parameters -------- ##

      if(verb && identical(i%%1000, 0)) cat("\n", i, " Mixing Parameters")


      ## ...... Pop Count ...... ##

      if(all(!is.na(lambda.n))) {

      ## For lambda.n

      log.curr.lx <- log.curr.proj[,c(1, census.columns)]

          # - Proposal

          #.. cycle through components
          for(j in 1:length(curr.lambda.n)) {

              prop.lambda.n[j] <-
                  popReconAux.rinvGamma(1, 3/2,
                    1 + 0.5/curr.sigmasq.n*(log.curr.lx[j] -
                                    log.mean.lx[j])^2
                                            )

        # - Calculate log posterior of proposed vital under projection
        log.prop.posterior <-
              popReconModCheck.log.post(f = log.curr.f
                             ,s = logit.curr.s
                             ,g = curr.g
                             ,baseline.n = log.curr.b
                             ,prior.mean.f = log.mean.f
                             ,prior.mean.s = logit.mean.s
                             ,prior.mean.g = mean.g
                             ,prior.mean.b = log.mean.b
                             ,alpha.f = al.f, beta.f = be.f
                             ,alpha.s = al.s, beta.s = be.s
                             ,alpha.g = al.g, beta.g = be.g
                             ,alpha.n = al.n, beta.n = be.n
                             ,sigmasq.f = curr.sigmasq.f
                             ,sigmasq.s = curr.sigmasq.s
                             ,sigmasq.g = curr.sigmasq.g
                             ,sigmasq.n = curr.sigmasq.n
                            ,lambda.f = curr.lambda.f
                       ,lambda.s = curr.lambda.s
                       ,lambda.g = curr.lambda.g
                       ,lambda.n = prop.lambda.n   #<-- use proposal
                             ,log.like = popReconModCheck.log.lhood(
                              log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = prop.lambda.n[,-1] #<-- use proposal
                              ,cen.col = census.columns
                              )
                             ,non.zero.fert = fert.rows
                             )

      if(i > burn.in) acc.count$lambda.n[j] <-
                    acc.count$lambda.n[j] + 1/n.iter
                curr.lambda.n[j] <- prop.lambda.n[j]
                log.curr.posterior <- log.prop.posterior

          } # ends: 'for(j in 1:length(curr.lambda.n))'

          if(i > burn.in) {
              lambda.n.mcmc[k,] <- as.vector(curr.lambda.n)
          }

      } # ends: 'if(!is.na(lambda.n) && !missing(lambda.n))'


      ##...... Log likelihood and log posterior

      if(save.log.lhood && i > burn.in) log.lhood.mcmc[k,] <-
          popReconModCheck.log.lhood(log.n.census = log.census.mat
                              ,log.n.hat = log.curr.proj
                              ,ll.var = curr.sigmasq.n
                              ,ll.lambda = curr.lambda.n[,-1]
                              ,cen.col = census.columns
                              )

      if(save.log.post && i > burn.in) log.posterior.mcmc[k,] <-
          log.curr.posterior


      if(verb && identical(i%%1000, 0)) cat("\n")

  } # Ends outer-most loop


    ## ......... End Loop ........ ##
    #...............................#


    ## ---------- Output --------- ##

    #cat("inital values", "\n\n")
    #.. initial values
    init.vals <- list(fert.rate = init.f
                      ,surv.prop = init.s
                      ,mig.XXX = init.g
                      ,baseline.count = init.b
                      ,init.sigmasq.f = init.sigmasq.f
                      ,init.sigmasq.s = init.sigmasq.s
                      ,init.sigmasq.g = init.sigmasq.g
                      ,init.sigmasq.n = init.sigmasq.n
                      ,pop.data = pop.data
                      )

    #.. fixed parameters
    fixed.params <- list(mean.fert.rate = mean.f
                         ,mean.surv.prop = mean.s
                         ,mean.mig.XXX = mean.g
                         ,mean.baseline.count = mean.b
                         ,pop.counts = pop.data
                         ,alpha.fert.rate = al.f
                         ,beta.fert.rate = be.f
                         ,alpha.surv.prop = al.s
                         ,beta.surv.prop = be.s
                         ,alpha.mig.XXX = al.g
                         ,beta.mig.XXX = be.g
                         ,alpha.population.count = al.n
                         ,beta.population.count = be.n
                         )

    #.. migration type
    new.names <-
        lapply(list(init.vals, fixed.params), FUN = function(z) {
        sub("XXX", mig.string, names(z))
    })
    names(init.vals) <- new.names[[1]]
    names(fixed.params) <- new.names[[2]]


    #cat("algorithm statistics", "\n\n")
    #.. algorithm statistics
    alg.stats <-
        list(acceptance.proportions = acc.count
             ,pop.went.neg = pop.negative
             ,acc.prop.adj4neg = mapply(FUN = function(a, b, n) {
                 (a * n) / (n - b)
             },
              acc.count[1:4], pop.negative, MoreArgs = list(n = n.iter)
              )
             ,acc.rat.na = ar.na
             ,surv.outside.tol = s.out.tol
             ,run.time = proc.time() - ptm
             )

    #cat("algorithm parameters", "\n\n")
    #.. algorithm parameters
    alg.params <- list(prop.vars = prop.vars
                       ,vital.transformations = list(fert.rate = "log"
                        ,surv.prob = "logit", mig.XXX = "I"
                        ,baseline.count = "log"
                        ,population.count = "log")
                       ,projection.periods = proj.periods
                       ,age.gp.size = age.size
                       ,cen.columns = census.columns
                       ,non.zero.fert.rows = fert.rows
                       ,surv.tolerance = s.tol
                       ,burn.in = burn.in
                       ,iters = n.iter
                       )
    names(alg.params) <- sub("XXX", mig.string, names(alg.params))

    ret.list <- list(fert.rate.mcmc = fert.rate.mcmc
                  ,surv.prop.mcmc = surv.prop.mcmc
                  ,mig.XXX.mcmc = mig.mcmc
                  ,baseline.count.mcmc = baseline.count.mcmc
                  ,lx.mcmc = lx.mcmc
                  ,variances.mcmc = variances.mcmc
                     ,log.lhood.mcmc = log.lhood.mcmc
                     ,log.posterior.mcmc = log.posterior.mcmc
                  ,alg.stats = alg.stats
                  ,fixed.params = fixed.params
                  ,init.vals = init.vals
                  ,alg.params = alg.params
                  )

    if(all(!is.na(lambda.f))) ret.list$lambda.f.mcmc = lambda.f.mcmc
    if(all(!is.na(lambda.s))) ret.list$lambda.s.mcmc = lambda.s.mcmc
    if(all(!is.na(lambda.g))) ret.list$lambda.g.mcmc = lambda.g.mcmc
    if(all(!is.na(lambda.n))) ret.list$lambda.n.mcmc = lambda.n.mcmc

    names(ret.list) <- sub("XXX", mig.string, names(ret.list))

    return(ret.list)

}


###
### * Metropolis proposal variances
###
################################################################################

###
### Simulation Study
###

sim.study.prop.vars <-
structure(list(fert.rate = structure(c(0.448568913121326, 0.678085805607733,
0.376340007056286, 0.446762140335542, 0.270932616397807, 0.429787179221642,
0.472826286250801, 0.265989972268545), .Dim = c(2L, 4L), .Dimnames = list(
    c("5", "10"), c("[1960, 1965)", "[1965, 1970)", "[1970, 1975)",
    "[1975, 1980)"))), surv.prop = structure(c(0.31961023033962,
0.232144042768031, 0.18106929310677, 0.304251566349687, 0.228382784230468,
0.204672022127278, 0.31757857318862, 0.201971113268454, 0.252914221932894,
0.174910693964099, 0.191639667748973, 0.184714418878724, 0.29203338625309,
0.264296899855242, 0.306949322757357, 0.227235773817815, 0.265817871771542,
0.219088345142732, 0.236019712394802, 0.218993812239991), .Dim = c(5L,
4L), .Dimnames = list(c("0", "5", "10", "15", "20+"), c("[1960, 1965)",
"[1965, 1970)", "[1970, 1975)", "[1975, 1980)"))), mig.prop = structure(c(0.29043064523009,
0.229320276154867, 0.282999839676462, 0.219184457148056, 0.141505473574661,
0.230219635648095, 0.211965862273689, 0.304040767787117, 0.287305839913259,
0.154829536476942, 0.247980244509981, 0.295432082070418, 0.221995569051617,
0.248493112069421, 0.211853343091724, 0.283736632321806), .Dim = c(4L,
4L), .Dimnames = list(c("0", "5", "10", "15"), c("[1960, 1965)",
"[1965, 1970)", "[1970, 1975)", "[1975, 1980)"))), population.count = structure(c(0.202344041554399,
0.215208028372875, 0.297643263171607, 0.418076005514405), .Dim = c(4L,
1L), .Dimnames = list(c("0", "5", "10", "15"), "1960"))), .Names = c("fert.rate",
"surv.prop", "mig.prop", "population.count"))


###
### Burkina Faso Reconstruction
###

bkfas.recon.prop.vars <-
structure(list(fert.rate = structure(c(0.15890025128247, 0.128435887523435,
0.158800127971759, 0.157539610301814, 0.172212860597094, 0.179320809119963,
0.173950687806731, 0.151053345520253, 0.126080751651391, 0.146912647038513,
0.176699221147778, 0.188621799120083, 0.178089815928788, 0.176672868398622,
0.138081226329162, 0.122907529478712, 0.138550870039475, 0.148735542340215,
0.150547371102324, 0.18655741310071, 0.177416246118741, 0.146469011635455,
0.135732232488195, 0.148409207346217, 0.147582819854948, 0.161958798745249,
0.183899134866859, 0.171518258512194, 0.133163074795245, 0.114175949998489,
0.136569208380193, 0.153080360452969, 0.178245079031769, 0.166741294173792,
0.204218894719422, 0.139423362521161, 0.120397439946617, 0.129656584640329,
0.160137696483466, 0.15810361397889, 0.179600749856407, 0.186884646864762,
0.147266225439485, 0.129861494316158, 0.130799987836118, 0.152810481465888,
0.166151534095606, 0.177543910276845, 0.164936682097412, 0.161887374162411,
0.14076660324777, 0.151803839414205, 0.162168935386418, 0.159076489265271,
0.178722493549484, 0.187032515067474, 0.152201615060823, 0.128270024978973,
0.148164971651057, 0.158534504056821, 0.174243121601706, 0.16559737991186,
0.179791489586071), .Names = c("1960.15", "1960.20", "1960.25",
"1960.30", "1960.35", "1960.40", "1960.45", "1965.15", "1965.20",
"1965.25", "1965.30", "1965.35", "1965.40", "1965.45", "1970.15",
"1970.20", "1970.25", "1970.30", "1970.35", "1970.40", "1970.45",
"1975.15", "1975.20", "1975.25", "1975.30", "1975.35", "1975.40",
"1975.45", "1980.15", "1980.20", "1980.25", "1980.30", "1980.35",
"1980.40", "1980.45", "1985.15", "1985.20", "1985.25", "1985.30",
"1985.35", "1985.40", "1985.45", "1990.15", "1990.20", "1990.25",
"1990.30", "1990.35", "1990.40", "1990.45", "1995.15", "1995.20",
"1995.25", "1995.30", "1995.35", "1995.40", "1995.45", "2000.15",
"2000.20", "2000.25", "2000.30", "2000.35", "2000.40", "2000.45"
)), surv.prop = structure(c(0.137260141511502, 0.159267739154683,
0.177502063818399, 0.172607787424726, 0.193741942730232, 0.160325571879268,
0.17192937342665, 0.193299274478368, 0.17456700740519, 0.16013406003179,
0.179958175892577, 0.154233212736466, 0.159087457263132, 0.163884233185136,
0.151762332096177, 0.18201231256642, 0.174736089141227, 0.175284444506721,
0.143722270904652, 0.158176825811976, 0.185453802702741, 0.165962733073859,
0.181420797249671, 0.175692520028528, 0.186523335963186, 0.163060179982782,
0.165088724468986, 0.171983910122981, 0.179283348055636, 0.151833740073798,
0.168026543723285, 0.158335601187387, 0.139798503455505, 0.147005901296327,
0.192905138877997, 0.172385190480876, 0.152694074345174, 0.168135155075832,
0.190012075687225, 0.173771372025207, 0.176614486793443, 0.182060972143073,
0.17876184307892, 0.174083017311653, 0.173884362216034, 0.16190160440031,
0.167508835263949, 0.160132013198273, 0.16489997500489, 0.136490057641145,
0.142710031220513, 0.122350789937005, 0.132974573588302, 0.179761827279125,
0.144860277566249, 0.156156077704994, 0.171170989743701, 0.178227072768334,
0.176524054724672, 0.164083456128903, 0.16323626968824, 0.173907332127376,
0.171075346139424, 0.157316954632342, 0.164263251168538, 0.178556714386141,
0.147516287443203, 0.13793951283755, 0.144659489448394, 0.133777458403318,
0.17134558108594, 0.178669958605215, 0.162615668763272, 0.164547587661282,
0.167899096533453, 0.179524775854362, 0.180025111632765, 0.162442160650646,
0.162926699869164, 0.182883941590439, 0.181874224305221, 0.17079398999539,
0.157924108212964, 0.156591895510267, 0.15487439771208, 0.134322388219807,
0.116410597700305, 0.114163112425803, 0.121037402863299, 0.156753198610677,
0.177503947216426, 0.167679410183366, 0.198426005221371, 0.169700283106636,
0.168888043035313, 0.170962024606019, 0.195876394577763, 0.182457281703608,
0.179258009749344, 0.171365004772588, 0.166310719613666, 0.158930079940361,
0.146353486998369, 0.157073548095607, 0.140275690241442, 0.135188261613274,
0.180434675558486, 0.188703908784251, 0.164271205279469, 0.175256951545873,
0.169718429380865, 0.190561476195709, 0.170888379233459, 0.16666009076583,
0.164625550750267, 0.174271388557425, 0.171031107567962, 0.156014506689061,
0.165769074595638, 0.155657812795294, 0.161673183173946, 0.152193138993943,
0.132763638402091, 0.105757264681696, 0.120209856651568, 0.162145728610074,
0.163001531937429, 0.184888839960852, 0.174065773566214, 0.170120970407304,
0.168890155088214, 0.17296104776491, 0.180586517843368, 0.152334839297501,
0.180267158628313, 0.177943762360936, 0.170192397708555, 0.171086998201058,
0.169990713469832, 0.151502058361146, 0.146629993888468, 0.132228926817977,
0.175570156817795, 0.196060278107832, 0.169285964560373, 0.172802971075143,
0.169873560429982, 0.182479446023037, 0.175435522194558, 0.190683628389136,
0.187710534126738, 0.175533179973077, 0.181780063216036, 0.16695871990962,
0.163437219463117, 0.168485561148718, 0.16533744291764, 0.152306901711214,
0.148290803675376, 0.117249368229747, 0.133280805825443, 0.161876865380265
), .Names = c("1960.0", "1960.5", "1960.10", "1960.15", "1960.20",
"1960.25", "1960.30", "1960.35", "1960.40", "1960.45", "1960.50",
"1960.55", "1960.60", "1960.65", "1960.70", "1960.75", "1960.80",
"1960.85", "1965.0", "1965.5", "1965.10", "1965.15", "1965.20",
"1965.25", "1965.30", "1965.35", "1965.40", "1965.45", "1965.50",
"1965.55", "1965.60", "1965.65", "1965.70", "1965.75", "1965.80",
"1965.85", "1970.0", "1970.5", "1970.10", "1970.15", "1970.20",
"1970.25", "1970.30", "1970.35", "1970.40", "1970.45", "1970.50",
"1970.55", "1970.60", "1970.65", "1970.70", "1970.75", "1970.80",
"1970.85", "1975.0", "1975.5", "1975.10", "1975.15", "1975.20",
"1975.25", "1975.30", "1975.35", "1975.40", "1975.45", "1975.50",
"1975.55", "1975.60", "1975.65", "1975.70", "1975.75", "1975.80",
"1975.85", "1980.0", "1980.5", "1980.10", "1980.15", "1980.20",
"1980.25", "1980.30", "1980.35", "1980.40", "1980.45", "1980.50",
"1980.55", "1980.60", "1980.65", "1980.70", "1980.75", "1980.80",
"1980.85", "1985.0", "1985.5", "1985.10", "1985.15", "1985.20",
"1985.25", "1985.30", "1985.35", "1985.40", "1985.45", "1985.50",
"1985.55", "1985.60", "1985.65", "1985.70", "1985.75", "1985.80",
"1985.85", "1990.0", "1990.5", "1990.10", "1990.15", "1990.20",
"1990.25", "1990.30", "1990.35", "1990.40", "1990.45", "1990.50",
"1990.55", "1990.60", "1990.65", "1990.70", "1990.75", "1990.80",
"1990.85", "1995.0", "1995.5", "1995.10", "1995.15", "1995.20",
"1995.25", "1995.30", "1995.35", "1995.40", "1995.45", "1995.50",
"1995.55", "1995.60", "1995.65", "1995.70", "1995.75", "1995.80",
"1995.85", "2000.0", "2000.5", "2000.10", "2000.15", "2000.20",
"2000.25", "2000.30", "2000.35", "2000.40", "2000.45", "2000.50",
"2000.55", "2000.60", "2000.65", "2000.70", "2000.75", "2000.80",
"2000.85")), mig.prop = structure(c(0.0367838682332873, 0.0334926424453817,
0.0287718492403349, 0.0286623236245622, 0.0272341222518445, 0.0300545300880446,
0.0328370306385146, 0.0321744384785064, 0.0367231600746433, 0.0458440191609768,
0.0511693268544198, 0.0596111431918158, 0.0766909542873833, 0.0732515778264404,
0.0948787183838824, 0.117752080066202, 0.123272862782605, 0.0398398021607279,
0.0359070266081482, 0.0318229889397586, 0.0302620839443824, 0.0312161590306619,
0.0330942484569774, 0.0353917391054951, 0.0344462691009954, 0.0347153471495791,
0.0414944963934467, 0.0429031338679918, 0.0565876095125219, 0.0592819148028749,
0.0576804298324959, 0.0768618083251519, 0.104942007595148, 0.124965523913564,
0.0365743057937032, 0.0419611611820298, 0.0378936882081968, 0.0357790314274174,
0.0308932721743071, 0.0305239167051393, 0.0316677278096958, 0.0328382112953607,
0.0303339667068596, 0.0348711240568865, 0.0434224498093322, 0.0484096979568276,
0.0484884480452959, 0.0549605777203043, 0.0672467430974446, 0.080162236591403,
0.0956099058248957, 0.0400118133389447, 0.0289433115402191, 0.0307470086843747,
0.026571661418974, 0.0235994977938011, 0.0262015868047237, 0.0258589304983994,
0.0278469378348158, 0.0269464002184104, 0.0306673407768787, 0.0353287698428431,
0.0371330154677031, 0.0452073014979137, 0.0623635024117327, 0.0846565104915526,
0.0917632546190314, 0.109887474239002, 0.0417698476196917, 0.0380117850354946,
0.0282597859190647, 0.0305245888819423, 0.032345935748905, 0.0249079161583914,
0.0272206720523319, 0.0261765997349203, 0.0296311250990949, 0.0272270246747285,
0.0326526487006855, 0.0374838711415061, 0.0406229928319442, 0.0495991414678239,
0.0600534974739972, 0.080972193288707, 0.0732398002622561, 0.0411508795503781,
0.0362668093564274, 0.0331859372818122, 0.0299670690936307, 0.0308715867223346,
0.030548371471404, 0.0282848006314123, 0.0300961300822108, 0.0293918755840028,
0.0297992817797593, 0.0330521472548004, 0.0376533359337197, 0.0453675323385973,
0.0510728521754444, 0.0709280937050579, 0.0791411835673709, 0.120994738309224,
0.0464035374164592, 0.0415857829431642, 0.0359681006826152, 0.0337690328455106,
0.0316070972035732, 0.0310107381373758, 0.0300393735157866, 0.0308944906648491,
0.0279083411074039, 0.0276735446536492, 0.0334027911599686, 0.0358183061712309,
0.0361892722791032, 0.048675829512846, 0.0571735320272543, 0.0758430373567911,
0.0708926777601453, 0.048244438638433, 0.0482588239214644, 0.0447463637693681,
0.039425636187294, 0.0406180963728261, 0.0386320550086539, 0.0410260220251613,
0.037969317184402, 0.0404084353768513, 0.0383613938986795, 0.04007874632952,
0.0456739013541493, 0.0492403080501581, 0.0523096284499509, 0.0736581565758044,
0.0835419115693889, 0.116587123354021, 0.0548912040972381, 0.0544971366508735,
0.0468130020189355, 0.0503422089151418, 0.046498746731505, 0.0483318848303559,
0.0464411478697621, 0.0473965821498648, 0.0460287198644804, 0.0482099160661111,
0.0416229080553777, 0.0392506707628495, 0.0494621433547527, 0.0465297761657499,
0.0540418873686472, 0.0756026451363889, 0.0710070878779001), .Names = c("1960.0",
"1960.5", "1960.10", "1960.15", "1960.20", "1960.25", "1960.30",
"1960.35", "1960.40", "1960.45", "1960.50", "1960.55", "1960.60",
"1960.65", "1960.70", "1960.75", "1960.80", "1965.0", "1965.5",
"1965.10", "1965.15", "1965.20", "1965.25", "1965.30", "1965.35",
"1965.40", "1965.45", "1965.50", "1965.55", "1965.60", "1965.65",
"1965.70", "1965.75", "1965.80", "1970.0", "1970.5", "1970.10",
"1970.15", "1970.20", "1970.25", "1970.30", "1970.35", "1970.40",
"1970.45", "1970.50", "1970.55", "1970.60", "1970.65", "1970.70",
"1970.75", "1970.80", "1975.0", "1975.5", "1975.10", "1975.15",
"1975.20", "1975.25", "1975.30", "1975.35", "1975.40", "1975.45",
"1975.50", "1975.55", "1975.60", "1975.65", "1975.70", "1975.75",
"1975.80", "1980.0", "1980.5", "1980.10", "1980.15", "1980.20",
"1980.25", "1980.30", "1980.35", "1980.40", "1980.45", "1980.50",
"1980.55", "1980.60", "1980.65", "1980.70", "1980.75", "1980.80",
"1985.0", "1985.5", "1985.10", "1985.15", "1985.20", "1985.25",
"1985.30", "1985.35", "1985.40", "1985.45", "1985.50", "1985.55",
"1985.60", "1985.65", "1985.70", "1985.75", "1985.80", "1990.0",
"1990.5", "1990.10", "1990.15", "1990.20", "1990.25", "1990.30",
"1990.35", "1990.40", "1990.45", "1990.50", "1990.55", "1990.60",
"1990.65", "1990.70", "1990.75", "1990.80", "1995.0", "1995.5",
"1995.10", "1995.15", "1995.20", "1995.25", "1995.30", "1995.35",
"1995.40", "1995.45", "1995.50", "1995.55", "1995.60", "1995.65",
"1995.70", "1995.75", "1995.80", "2000.0", "2000.5", "2000.10",
"2000.15", "2000.20", "2000.25", "2000.30", "2000.35", "2000.40",
"2000.45", "2000.50", "2000.55", "2000.60", "2000.65", "2000.70",
"2000.75", "2000.80")), population.count = structure(c(0.0180188163635564,
0.0165199993263465, 0.0179579947966866, 0.0156115741969201, 0.0178513102509598,
0.0162693262907368, 0.0140648016206893, 0.0181309005246593, 0.0188722416962788,
0.0193527869193709, 0.0198231112837106, 0.0204236548990628, 0.0223064462910951,
0.0234623698129105, 0.0222003552507938, 0.0237709222501987, 0.023660576047669
), .Names = c("1960.0", "1960.5", "1960.10", "1960.15", "1960.20",
"1960.25", "1960.30", "1960.35", "1960.40", "1960.45", "1960.50",
"1960.55", "1960.60", "1960.65", "1960.70", "1960.75", "1960.80"
))), .Names = c("fert.rate", "surv.prop", "mig.prop", "population.count"
))


###
### Model checking
###

mod.check.prop.vars <-
    structure(list(fert.rate = structure(c(0.0679657894611446, 0.0872058314869657,
0.0740033557933105, 0.0983934665242664, 0.113524244724523, 0.102564312075842,
0.101843476514153, 0.0835533788454877, 0.0715614432150811, 0.0658950166581847,
0.081631578382636, 0.100891404955363, 0.105822096632364, 0.115392975596451,
0.0883948288757947, 0.066886245737006, 0.0920079465665007, 0.0882533808311936,
0.0990115250156988, 0.105605030283607, 0.12093520050231, 0.0715330952383307,
0.0819217815928455, 0.0743844118335317, 0.0800943688822749, 0.0948175273211659,
0.122966116905382, 0.108824680899105, 0.0790456311958763, 0.070163339775771,
0.0712766962225196, 0.0874247857092805, 0.101308301787391, 0.106164430111148,
0.111269587789541, 0.0819196713306691, 0.068689628130059, 0.0810834068233205,
0.0963229633227391, 0.0980054151880907, 0.106172777021201, 0.11708353037486,
0.091362888985388, 0.0608036451743309, 0.0759288163744214, 0.0900554434009259,
0.102275820119393, 0.112074597301743, 0.113750157973617, 0.0914056169167423,
0.0756073724738213, 0.0929447598374884, 0.0885915902156074, 0.106234262598218,
0.110756682961373, 0.115498654318946, 0.0905451334321341, 0.0769193367822608,
0.0818065948408167, 0.087874225036128, 0.10427753293744, 0.107679419492818,
0.121960496692417), .Names = c("1960.15", "1960.20", "1960.25",
"1960.30", "1960.35", "1960.40", "1960.45", "1965.15", "1965.20",
"1965.25", "1965.30", "1965.35", "1965.40", "1965.45", "1970.15",
"1970.20", "1970.25", "1970.30", "1970.35", "1970.40", "1970.45",
"1975.15", "1975.20", "1975.25", "1975.30", "1975.35", "1975.40",
"1975.45", "1980.15", "1980.20", "1980.25", "1980.30", "1980.35",
"1980.40", "1980.45", "1985.15", "1985.20", "1985.25", "1985.30",
"1985.35", "1985.40", "1985.45", "1990.15", "1990.20", "1990.25",
"1990.30", "1990.35", "1990.40", "1990.45", "1995.15", "1995.20",
"1995.25", "1995.30", "1995.35", "1995.40", "1995.45", "2000.15",
"2000.20", "2000.25", "2000.30", "2000.35", "2000.40", "2000.45"
)), surv.prop = structure(c(0.0956530713737909, 0.119219973708042,
0.13363585944471, 0.138201430453909, 0.13580072658239, 0.135607102414324,
0.128671747089422, 0.126371534557923, 0.123378051469723, 0.117226812248627,
0.135170134407619, 0.109551664729342, 0.0981169430800515, 0.111728976445404,
0.11987741055006, 0.133443231190329, 0.12957103614845, 0.136149355809066,
0.0809620389701372, 0.115828952312828, 0.145284954099682, 0.134469487275181,
0.128612224321961, 0.123890061269458, 0.130747329288238, 0.144664901901044,
0.127967315143937, 0.120457289146385, 0.116751208183692, 0.110716806541191,
0.111083485246993, 0.0978788206303483, 0.0933856869357953, 0.103967509208338,
0.135965804940648, 0.149010211348832, 0.105237244032774, 0.112912681898563,
0.131229877454335, 0.14100354776511, 0.126323702907069, 0.141996773953816,
0.145206579478077, 0.126532426579024, 0.136934322130555, 0.127785892444828,
0.105811286931665, 0.106109708581669, 0.113857400562148, 0.112963689080217,
0.0855635261625545, 0.098203003831793, 0.108388898087797, 0.137648076696581,
0.120303946335625, 0.127397077881253, 0.132082186315459, 0.136489149244836,
0.137710432425407, 0.133613645359137, 0.127907222459612, 0.141394714821169,
0.132474613581119, 0.140963317112955, 0.129593522208004, 0.116164990023376,
0.108697727841022, 0.0979762705155692, 0.0964260163433794, 0.0950954233938129,
0.12895995750729, 0.141853762735825, 0.104141972782925, 0.127603620184087,
0.132352116249389, 0.130206196956693, 0.138206680581263, 0.125165513920202,
0.133775840842667, 0.12793297211632, 0.140069230879924, 0.134787023551576,
0.122756393941599, 0.117810230638655, 0.106182588305529, 0.112128843966124,
0.0934409957616603, 0.0756124155066622, 0.080196147502314, 0.120079250897558,
0.121599242529941, 0.130002994495422, 0.138997045784952, 0.13457622980335,
0.138422286525152, 0.12413823148119, 0.122432560546141, 0.134067833579873,
0.140812260134032, 0.121475282341184, 0.123687946444476, 0.114241932554314,
0.110696471699952, 0.106716891056951, 0.0882841003707008, 0.0874618607895278,
0.131864413101215, 0.135411121010439, 0.11266815820858, 0.126937061750499,
0.150793357034225, 0.14206511697302, 0.131022241788874, 0.131535118596494,
0.147419082437802, 0.130085019970602, 0.134557651582294, 0.146537816566472,
0.131095453917605, 0.129465921410239, 0.119908010450344, 0.111061314324475,
0.100689911847845, 0.0817012573918061, 0.0877559049390283, 0.110633811401435,
0.128769539102236, 0.126872147648235, 0.132221092069512, 0.125003153354043,
0.134451103825487, 0.145135296505107, 0.135412089017849, 0.12558585498981,
0.136648889588795, 0.141600190883658, 0.135992201354243, 0.12518250178506,
0.120047536477499, 0.101609708821878, 0.107292074682368, 0.107313824253305,
0.124338491156851, 0.140429916315076, 0.116784843563008, 0.13516637152537,
0.142255298563698, 0.140198789641248, 0.137778439328925, 0.129753585767575,
0.134525594472373, 0.129438177051389, 0.140173384882168, 0.149286743082801,
0.136462381701476, 0.1245567441415, 0.120653684613281, 0.107441125388483,
0.100115167487333, 0.0907245562362851, 0.0863110965153816, 0.118458770174707
), .Names = c("1960.0", "1960.5", "1960.10", "1960.15", "1960.20",
"1960.25", "1960.30", "1960.35", "1960.40", "1960.45", "1960.50",
"1960.55", "1960.60", "1960.65", "1960.70", "1960.75", "1960.80",
"1960.85", "1965.0", "1965.5", "1965.10", "1965.15", "1965.20",
"1965.25", "1965.30", "1965.35", "1965.40", "1965.45", "1965.50",
"1965.55", "1965.60", "1965.65", "1965.70", "1965.75", "1965.80",
"1965.85", "1970.0", "1970.5", "1970.10", "1970.15", "1970.20",
"1970.25", "1970.30", "1970.35", "1970.40", "1970.45", "1970.50",
"1970.55", "1970.60", "1970.65", "1970.70", "1970.75", "1970.80",
"1970.85", "1975.0", "1975.5", "1975.10", "1975.15", "1975.20",
"1975.25", "1975.30", "1975.35", "1975.40", "1975.45", "1975.50",
"1975.55", "1975.60", "1975.65", "1975.70", "1975.75", "1975.80",
"1975.85", "1980.0", "1980.5", "1980.10", "1980.15", "1980.20",
"1980.25", "1980.30", "1980.35", "1980.40", "1980.45", "1980.50",
"1980.55", "1980.60", "1980.65", "1980.70", "1980.75", "1980.80",
"1980.85", "1985.0", "1985.5", "1985.10", "1985.15", "1985.20",
"1985.25", "1985.30", "1985.35", "1985.40", "1985.45", "1985.50",
"1985.55", "1985.60", "1985.65", "1985.70", "1985.75", "1985.80",
"1985.85", "1990.0", "1990.5", "1990.10", "1990.15", "1990.20",
"1990.25", "1990.30", "1990.35", "1990.40", "1990.45", "1990.50",
"1990.55", "1990.60", "1990.65", "1990.70", "1990.75", "1990.80",
"1990.85", "1995.0", "1995.5", "1995.10", "1995.15", "1995.20",
"1995.25", "1995.30", "1995.35", "1995.40", "1995.45", "1995.50",
"1995.55", "1995.60", "1995.65", "1995.70", "1995.75", "1995.80",
"1995.85", "2000.0", "2000.5", "2000.10", "2000.15", "2000.20",
"2000.25", "2000.30", "2000.35", "2000.40", "2000.45", "2000.50",
"2000.55", "2000.60", "2000.65", "2000.70", "2000.75", "2000.80",
"2000.85")), mig.prop = structure(c(0.0121918502868898, 0.00999565637258205,
0.0107371654916488, 0.0062865492679386, 0.0105086515894261, 0.010357445967544,
0.0137781950441233, 0.0140316176613681, 0.0167443260512334, 0.0187694620080219,
0.0204114222878633, 0.0242912330504285, 0.0299623028046512, 0.0461195045071805,
0.0598230073300098, 0.0654935547608869, 0.0649277572912325, 0.0118651420121959,
0.00951024448619292, 0.0090990203248332, 0.00872726469887084,
0.00855568141101466, 0.0107874577426659, 0.0106809056800966,
0.0119147426423729, 0.0109884707877641, 0.0159414313546304, 0.021436832855021,
0.0224041126024035, 0.028725642229372, 0.0337855363325949, 0.0434414221705617,
0.0529184914181428, 0.0655746950333709, 0.0150862295823367, 0.0139204831346726,
0.0120476781507473, 0.0086484104627752, 0.0121935195099419, 0.00800848754854404,
0.0100859987191787, 0.0125958129643897, 0.0133900797050486, 0.0129330581822785,
0.0156259652854098, 0.0214886267635095, 0.0180946033568137, 0.0286856303513036,
0.0335244693998931, 0.0376777778989273, 0.050597984274937, 0.015436601333411,
0.0125559439054929, 0.0113550342463928, 0.0115402741186818, 0.00961629107266085,
0.00941846939116588, 0.00936162492549021, 0.0105200198356926,
0.0111142755819808, 0.0131134574203158, 0.0126340470305322, 0.0184645031485599,
0.0192521599172882, 0.0310364169274558, 0.0404796812875343, 0.0424260006685085,
0.0625141255685481, 0.0139567870975637, 0.0174099092817319, 0.0132230694660476,
0.0104156753602705, 0.0105612921389045, 0.00868404324024959,
0.00958406879182552, 0.00908895508979391, 0.0102372615910236,
0.0112696806849857, 0.0120767892782536, 0.0156574353283235, 0.0199558552762386,
0.0245825847501297, 0.0279706940295383, 0.039144310452476, 0.0416407856929598,
0.0163661012696416, 0.0141189224226879, 0.0132397606516459, 0.0118899761610524,
0.0121092124398829, 0.0105895550702899, 0.0118709374035985, 0.0126519751835992,
0.012589875742888, 0.0130427485262672, 0.0159760523830636, 0.0148629368409338,
0.018928242240839, 0.0272715338025926, 0.0377971435958206, 0.0511373220314246,
0.0626401007242425, 0.0208685990807146, 0.015973969245736, 0.0148775323289653,
0.0131864366189477, 0.0123732695427402, 0.0114208767911365, 0.0131195442819032,
0.0119399168764373, 0.0119869538345245, 0.0127315216945412, 0.0138866146986357,
0.0170270909314589, 0.0169737352767235, 0.0188441715219736, 0.0262835259834077,
0.0359253278620205, 0.043812035085478, 0.0241840967214855, 0.0199044454999602,
0.0173995492509994, 0.0173176591668357, 0.0167186289469382, 0.017564838238766,
0.0156364920104624, 0.0180543399422108, 0.0193010955585084, 0.0194926194195306,
0.0189455456913936, 0.0199059271858012, 0.0242989958850205, 0.0300546854504128,
0.0362136522459139, 0.0467670999690369, 0.0653581153143469, 0.0254481820548441,
0.0243081290327067, 0.0212844925178571, 0.0191665331079444, 0.0211396051488891,
0.0211513306193079, 0.0190923501603394, 0.0194579382652995, 0.0204166915420872,
0.0213113068892425, 0.0196498300649752, 0.0195065648747406, 0.0218003266450716,
0.0239445919521544, 0.0266076890431113, 0.0375126246250088, 0.0382198434005055
), .Names = c("1960.0", "1960.5", "1960.10", "1960.15", "1960.20",
"1960.25", "1960.30", "1960.35", "1960.40", "1960.45", "1960.50",
"1960.55", "1960.60", "1960.65", "1960.70", "1960.75", "1960.80",
"1965.0", "1965.5", "1965.10", "1965.15", "1965.20", "1965.25",
"1965.30", "1965.35", "1965.40", "1965.45", "1965.50", "1965.55",
"1965.60", "1965.65", "1965.70", "1965.75", "1965.80", "1970.0",
"1970.5", "1970.10", "1970.15", "1970.20", "1970.25", "1970.30",
"1970.35", "1970.40", "1970.45", "1970.50", "1970.55", "1970.60",
"1970.65", "1970.70", "1970.75", "1970.80", "1975.0", "1975.5",
"1975.10", "1975.15", "1975.20", "1975.25", "1975.30", "1975.35",
"1975.40", "1975.45", "1975.50", "1975.55", "1975.60", "1975.65",
"1975.70", "1975.75", "1975.80", "1980.0", "1980.5", "1980.10",
"1980.15", "1980.20", "1980.25", "1980.30", "1980.35", "1980.40",
"1980.45", "1980.50", "1980.55", "1980.60", "1980.65", "1980.70",
"1980.75", "1980.80", "1985.0", "1985.5", "1985.10", "1985.15",
"1985.20", "1985.25", "1985.30", "1985.35", "1985.40", "1985.45",
"1985.50", "1985.55", "1985.60", "1985.65", "1985.70", "1985.75",
"1985.80", "1990.0", "1990.5", "1990.10", "1990.15", "1990.20",
"1990.25", "1990.30", "1990.35", "1990.40", "1990.45", "1990.50",
"1990.55", "1990.60", "1990.65", "1990.70", "1990.75", "1990.80",
"1995.0", "1995.5", "1995.10", "1995.15", "1995.20", "1995.25",
"1995.30", "1995.35", "1995.40", "1995.45", "1995.50", "1995.55",
"1995.60", "1995.65", "1995.70", "1995.75", "1995.80", "2000.0",
"2000.5", "2000.10", "2000.15", "2000.20", "2000.25", "2000.30",
"2000.35", "2000.40", "2000.45", "2000.50", "2000.55", "2000.60",
"2000.65", "2000.70", "2000.75", "2000.80")), population.count = structure(c(0.00394836932472217,
0.003582345496587, 0.00379050683653265, 0.00256835606955148,
0.00324380011401299, 0.00431199230880955, 0.0032939923589788,
0.00360739727624942, 0.00494228528548945, 0.00556983322013338,
0.00697721719293354, 0.00660221996581503, 0.00892637263757422,
0.0407024175282445, 0.0173838025414991, 0.0116054746526456, 0.0608714728668632
), .Names = c("1960.0", "1960.5", "1960.10", "1960.15", "1960.20",
"1960.25", "1960.30", "1960.35", "1960.40", "1960.45", "1960.50",
"1960.55", "1960.60", "1960.65", "1960.70", "1960.75", "1960.80"
))), .Names = c("fert.rate", "surv.prop", "mig.prop", "population.count"
))

