################################################################################
###
###  TITLE:             1_estimation.R
###
###  AUTHOR:            Mark Wheldon
###
###  DATE CREATED:      2016-11-16
###
###  DESCRIPTION:       Two-sex population reconstruction for the population of
###                     Thailand as described in the paper "Bayesian
###                     reconstruction of two-sex populations by age: estimating
###                     sex ratios at birth and sex ratios of mortality".
###
###  REFERENCE:         Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland,
###                     P. (2015). Bayesian Reconstruction of Two-Sex
###                     Populations by Age: Estimating Sex Ratios at Birth and
###                     Sex Ratios of Mortality. Journal of the Royal
###                     Statistical Society Series A, 178(4), 977--1007.
###
###  LICENCE:           Released under the Creative Commons BY-NC-SA Licence
###                     (https://creativecommons.org).
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
###  Main population reconstruction script.
###
###-----------------------------------------------------------------------------
###
### THIS RUN
###
message("DATE-TIME STAMP: ", (dateTimeStamp <- format(Sys.time(), "%Y%m%d_%H-%M")))
###
###-----------------------------------------------------------------------------
###
### CONTROLS
###
### SHORT RUN as an example only.
###
n.iter <- 1E3
burn.in <- 1E3
###
###
################################################################################

### * SET UP
################################################################################

set.seed(1)

library(coda)

data.path <- "data"
prog.path <- "R_Code"
results.path <- "results"

if(!file.exists(results.path)) dir.create(results.path)

example(source, echo = FALSE)
sourceDir(prog.path)

### * Inputs
################################################################################

### ** Initial Estimates
################################################################################

load(file.path(data.path, "thai_initial_ests.RData"))

### ** Arguments
################################################################################

### Inverse Gamma parameters
###

invGam.params <-
    make.hyper.params(absDev = list(fert = 0.1, surv = 0.1, mig = 0.2
             ,pop = 0.1, srb = 0.1)
             ,prob = list(fert = 0.9, surv = 0.9, mig = 0.9, pop = 0.9, srb = 0.9)
             ,alpha = list(fert = 0.5, surv = 0.5, mig = 0.5, pop = 0.5, srb = 0.5)
             ,s.star = unlist(asSurvTHAI.mat)
             )

invGam.params

### ** Proposal variances
################################################################################

withVisible(load(file = file.path("data", "thai_propvars.RData")))$value

### ** Estimation arguments
################################################################################

estModArgs <-
    list(## inverse gamma parameters
        al.f = invGam.params$al.f
       ,be.f = invGam.params$be.f
       ,al.s = invGam.params$al.s
       ,be.s = invGam.params$be.s
       ,al.g = invGam.params$al.g
       ,be.g = invGam.params$be.g
       ,al.n = invGam.params$al.n
       ,be.n = invGam.params$be.n
       ,al.srb = invGam.params$al.srb
       ,be.srb = invGam.params$be.srb

        ## the rest
       ,n.iter = n.iter, burn.in = burn.in
       ,ccmp.f = list(FUN = "ccmp.femDom.jun12")
       ,start.f = asFertTHAI.mat
       ,start.s = asSurvTHAI.mat
       ,start.g = asMigTHAI.mat
       ,start.b = baselineTHAI.mat
       ,start.srb = srbTHAI.mat
       ,start.sigmasq.f = 5
       ,start.sigmasq.s = 5
       ,start.sigmasq.g = 5
       ,start.sigmasq.n = 5
       ,start.sigmasq.srb = 5
       ,mean.f = asFertTHAI.mat
       ,mean.s = asSurvTHAI.mat
       ,mean.g = asMigTHAI.mat
       ,mean.b = baselineTHAI.mat
       ,mean.srb = srbTHAI.mat
       ,fert.rows = as.logical(apply(asFertTHAI.mat == 0L, 1, function(z) !all(z)))
       ,prop.varcovar = thai.propvar
       ,pop.data = censusTHAI.mat
       ,proj.periods = ncol(asFertTHAI.mat)
       ,age.size = 5
       ,s.tol = 10^(-10)
       ,verb = TRUE, progress.step = 1E3
       ,filebacked.chains = FALSE
       ,file.suff = ""
       ,v.prop.type = "MH.scaled"
       ,save.mid.run.name = "Thai"
    )

### * Run Estimation
################################################################################

ThaiMcmc <- do.call(popRecon.pop.est.sampler.aug07, args = estModArgs)

save(ThaiMcmc, file = file.path(results.path, "thai_mcmc.RData"))
