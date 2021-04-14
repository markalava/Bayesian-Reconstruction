################################################################################
###
###  TITLE:             1_sample_from_prior.R
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
###-----------------------------------------------------------------------------
###
###  SYNOPSIS
###
###  Generates a sample from the joint prior distribution.
###
###  Run 1_estimation.R first.
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
### SHORT RUN AS EXAMPLE ONLY
###
n.iter <- 100
###
###
################################################################################

### ############################################################################
### * SET UP

set.seed(1)

library(ggplot2)
## theme_set(theme_bw())
library(lattice)
library(coda)
library(gdata)
library(reshape2); dmelt <- function(...) reshape2:::melt.data.frame(...)

data.path <- "data"
prog.path <- "R_Code"
results.path <- "results"

example(source, echo = FALSE)
sourceDir(prog.path)

#
#.. pairs functions
#
example(pairs, ask = FALSE)
graphics.off()


## Default for quantile() is to remove NA
quantile <- function(x, probs = seq(0, 1, 0.25), na.rm = TRUE,
              names = TRUE, type = 7, ...)
{
    stats::quantile(x, probs = probs, na.rm = na.rm,
              names = names, type = type, ...)
}

graphics.off()


### ############################################################################
### * INPUTS: Data, initial estimates, life table columns, etc.

### ############################################################################
### ** MCMC output: Only need the meta data

load(file = file.path(results.path, "thai_mcmc.RData"))


### ############################################################################
### ** Separation Factors

load(file = file.path(data.path, "thai_sep_factors.RData"))


### ############################################################################
### * Sample From Prior

###
### Sample
###

sample.from.prior.feb07(n.iter = n.iter
                        ,batch.size = 1E2
                        ,pos.sample.max.tries = 1E9
                        ,max.elapsed.time = 3600*24*14
                        ,baseline = ThaiMcmc$fixed.params$mean.baseline.count
                        ,fert = ThaiMcmc$fixed.params$mean.fert.rate
                        ,surv = ThaiMcmc$fixed.params$mean.surv.prop
                        ,mig = ThaiMcmc$fixed.params$mean.mig.prop
                        ,srb = ThaiMcmc$fixed.params$mean.srb
                        ,hyper.params = ThaiMcmc$fixed.params
                        ,age.int = 5
                        ,label.dims = FALSE
                        ,ccmp.f = ccmp.femDom.jun12
                        ,name.pref = "Thai."
                        ,name.suf = ""
                        ,return.list = FALSE
             ,parallelize = TRUE
              ,pvm.hostfile = "~/.pvm_hosts_fal"
               ,output.dir = results.path
                        )


###
### Convert to counts
###

nF <- ncol(ThaiMcmc$fixed.params$mean.mig.prop$female) *
    nrow(ThaiMcmc$fixed.params$mean.mig.prop$female)

proj.to.counts.aug22(Thai.fert.rate.prior
                     ,Thai.surv.prop.prior
                     ,Thai.mig.prop.prior
                     ,Thai.baseline.count.prior
                     ,Thai.pop.count.prior
                     ,srb.mcmc = Thai.srb.prior
                     ,name.pref = "Thai."
                     ,name.suf = ".prior"
                     ,file.backed = TRUE
                     ,ccmp.f = "ccmp.femDom.jun12"
                     ,sep.factors = list(female = thaiFemale.sf[1:nF]
                      ,male = thaiMale.sf[1:nF]
                     )
                     )


## Save
save.MF.fflist(x = Thai.birth.count.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_birth_count_prior.RData"))
save.MF.fflist(Thai.birth.count.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_total_birth_count_prior.RData"))
save.MF.fflist(Thai.total.birth.death.count.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_total_birth_death_count_prior.RData"))
save.MF.fflist(Thai.death.count.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_death_count_prior.RData"))
save.MF.fflist(Thai.mig.count.mcmc.prior
     ,f = file.path(results.path
      , "Thai_mig_count_prior.RData"))
save.MF.fflist(Thai.total.mig.count.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_total_mig_count_prior.RData"))
save.MF.fflist(Thai.mort.rate.mcmc.prior
     ,f = file.path(results.path
      , "Thai_mort_rate_prior.RData"))
save.MF.fflist(Thai.mig.rate.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_mig_rate_prior.RData"))
save.MF.fflist(Thai.person.years.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_person_years_prior.RData"))
save.MF.fflist(Thai.cohort.nq0.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_cohort_nq0_prior.RData"))
save.MF.fflist(Thai.period.nq0.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_period_nq0_prior.RData"))
save.MF.fflist(Thai.IMR.mcmc.prior
     ,f = file.path(results.path
      ,"Thai_IMR_prior.RData"))
