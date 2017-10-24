################################################################################
###
###  TITLE:             3_post_process.R
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
###  Post-process population reconstruction: convert rates, proportions, to
###  counts.
###
###-----------------------------------------------------------------------------
###
### THIS RUN
###
message("DATE-TIME STAMP: ", (dateTimeStamp <- format(Sys.time(), "%Y%m%d_%H-%M")))
###
###
################################################################################

### ############################################################################
### * SET UP

set.seed(1)

library(ggplot2)
## theme_set(theme_bw())
library(ff)
library(lattice)
library(coda)
library(gdata)
library(reshape2); dmelt <- function(...) reshape2:::melt.data.frame(...)

data.path <- "data"
prog.path <- "R_Code"
results.path <- "results"

example(source, echo = FALSE)
sourceDir(prog.path)


### ############################################################################
### * INPUTS: Data, initial estimates, life table columns, etc.

### ############################################################################
### ** MCMC output

load(file = file.path(results.path, "thai_mcmc.RData"))


### ############################################################################
### ** Separation Factors

load(file = file.path(data.path, "thai_sep_factors.RData"))


### ############################################################################
### ** Constants

## Quantiles to calculate

quants.to.plot <- c(0.025, 0.1, 0.5, 0.9, 0.975)


### ############################################################################
### ** Calculate Counts from Posterior

nF <- ncol(ThaiMcmc$fixed.params$mean.mig.prop$female) *
    nrow(ThaiMcmc$fixed.params$mean.mig.prop$female)

proj.to.counts.aug22(fert.rate.mcmc = ThaiMcmc$fert.rate.mcmc
                     ,surv.prop.mcmc = list(female = ThaiMcmc$surv.prop.mcmc$female
                      ,male = ThaiMcmc$surv.prop.mcmc$male)
                     ,mig.prop.mcmc = list(female = ThaiMcmc$mig.prop.mcmc$female
                      ,male = ThaiMcmc$mig.prop.mcmc$male)
                     ,baseline.count.mcmc = list(female = ThaiMcmc$baseline.count.mcmc$female
                      ,male = ThaiMcmc$baseline.count.mcmc$male)
                     ,pop.count.mcmc = list(female = ThaiMcmc$lx.mcmc$female
                      ,male = ThaiMcmc$lx.mcmc$male)
                     ,sep.factors = list(female = thaiFemale.sf[1:nF]
                      ,male = thaiMale.sf[1:nF])
                     ,srb.mcmc = ThaiMcmc$srb.mcmc
                     ,name.pref = "Thai."
                     ,name.suf = ""
                     ,file.backed = TRUE
                     ,ccmp.f = "ccmp.femDom.jun12"
                     )

save.MF.fflist(Thai.birth.count.mcmc
               ,f = file.path(results.path,"Thai_birth_count_posterior.RData"))
save.MF.fflist(Thai.birth.count.mcmc
               ,f = file.path(results.path,"Thai_total_birth_count_posterior.RData"))
save.MF.fflist(Thai.total.birth.death.count.mcmc
               ,f = file.path(results.path,"Thai_total_birth_death_count_posterior.RData"))
save.MF.fflist(Thai.death.count.mcmc
               ,f = file.path(results.path,"Thai_death_count_posterior.RData"))
save.MF.fflist(Thai.mig.count.mcmc
               ,f = file.path(results.path, "Thai_mig_count_posterior.RData"))
save.MF.fflist(Thai.total.mig.count.mcmc
               ,f = file.path(results.path,"Thai_total_mig_count_posterior.RData"))
save.MF.fflist(Thai.mort.rate.mcmc
               ,f = file.path(results.path, "Thai_mort_rate_posterior.RData"))
save.MF.fflist(Thai.mig.rate.mcmc
               ,f = file.path(results.path,"Thai_mig_rate_posterior.RData"))
save.MF.fflist(Thai.person.years.mcmc
               ,f = file.path(results.path,"Thai_person_years_posterior.RData"))
save.MF.fflist(Thai.cohort.nq0.mcmc
               ,f = file.path(results.path,"Thai_cohort_nq0_posterior.RData"))
save.MF.fflist(Thai.period.nq0.mcmc
               ,f = file.path(results.path,"Thai_period_nq0_posterior.RData"))
save.MF.fflist(Thai.IMR.mcmc
               ,f = file.path(results.path,"Thai_IMR_posterior.RData"))
