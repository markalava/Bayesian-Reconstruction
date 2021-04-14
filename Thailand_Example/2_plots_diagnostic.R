################################################################################
###
###  TITLE:             2_plots_diagnostic.R
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
###  Diagnostic plots for a population reconstruction.
###
###  Run 1_estimation.R first.
###
###-----------------------------------------------------------------------------
###
### THIS RUN
###
message("DATE-TIME STAMP: ", (dateTimeStamp <- format(Sys.time(), "%Y%m%d_%H-%M")))
###
###
###-----------------------------------------------------------------------------
###
### CONTROLS
###
chain.length <- "all"
chain.thin <- 1
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

## Directories
dir.create(file.path("results", "plots", "traceplots"), recursive = TRUE)

### ############################################################################
### * INPUTS: Data, initial estimates, life table columns, etc.

### ############################################################################
### ** Constants

## Quantiles to calculate
quants.to.plot <- c(0.025, 0.1, 0.5, 0.9, 0.975)


## Use file.backed data?
fb <- TRUE


### ############################################################################
### ** MCMC output

load(file = file.path(results.path, "thai_mcmc.RData"), verbose = TRUE)



attach(ThaiMcmc)



### ##################################################################
### * PLOT POSTERIOR AND PRIOR STD. DEVIATIONS

### Inverse Gamma parameters

pdf(height = 7, width = 7
    ,file = "results/plots/traceplots/Thai_sd_priorPost_density.pdf")
par(mfrow = c(3,2))
postPrior.plot.sd(variances.mcmc[,"fert.rate.var"]
          ,dinvGamma, shape = fixed.params$alpha.fert.rate
               ,scale = fixed.params$beta.fert.rate
               ,main = expression(sigma[f])
               )
postPrior.plot.sd(variances.mcmc[,"surv.prop.var"]
          ,dinvGamma, shape = fixed.params$alpha.surv.prop
               ,scale = fixed.params$beta.surv.prop
               ,main = expression(sigma[s])
               )
postPrior.plot.sd(variances.mcmc[,"mig.var"]
          ,dinvGamma, shape = fixed.params$alpha.mig.prop
               ,scale = fixed.params$beta.mig.prop
               ,main = expression(sigma[g])
               )
postPrior.plot.sd(variances.mcmc[,"population.count.var"]
          ,dinvGamma, shape = fixed.params$alpha.population.count
               ,scale = fixed.params$beta.population.count
               ,main = expression(sigma[n])
               )
postPrior.plot.sd(variances.mcmc[,"srb.var"]
          ,dinvGamma, shape = fixed.params$alpha.srb
               ,scale = fixed.params$beta.srb
               ,main = expression(sigma[srb])
               )
dev.off()


### ##################################################################
### * TRACE PLOTS

###
###
### Variances
###
###
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
     ,file = "results/plots/traceplots/Thai_var_traceplot_1.png"
    )
print(plot(variances.mcmc[,1:2]))
dev.off()

png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
     ,file = "results/plots/traceplots/Thai_var_traceplot_2.png")
plot(variances.mcmc[,3:4])
dev.off()

png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
     ,file = "results/plots/traceplots/Thai_var_traceplot_3.png")
print(plot(variances.mcmc[,5,drop=FALSE]))
dev.off()

###
###
### ##################################################################
### SRB
###
###
if(ncol(srb.mcmc) > 4) {
srbVars.seq <- seq(from = 1, to = ncol(srb.mcmc)
                    ,by = 4)
for(i in 1:(length(srbVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logsrb_traceplot20130612srb10_", i, ".png"
       ,sep = ""))
  print(plot(log(srb.mcmc[,srbVars.seq[i]:
                                            (srbVars.seq[i]+3),drop = FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logsrb_traceplot20130612srb10_",
       length(srbVars.seq), ".png", sep = ""))
print(plot(log(
        srb.mcmc[,srbVars.seq[length(srbVars.seq)]:
                                   ncol(srb.mcmc),drop = FALSE])))
dev.off()
} else {
    png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logsrb_traceplot20130612srb10.png"
       ,sep = ""))
  print(plot(log(srb.mcmc)))
  dev.off()
}

###
###
### ##################################################################
### Fertility
###
###
fertVars.seq <- seq(from = 1, to = ncol(fert.rate.mcmc)
                    ,by = 4)
for(i in 1:(length(fertVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logfert_traceplot20130612srb10_", i, ".png"
       ,sep = ""))
  print(plot(log(fert.rate.mcmc[,fertVars.seq[i]:
                                            (fertVars.seq[i]+3),drop=FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logfert_traceplot20130612srb10_",
       length(fertVars.seq), ".png", sep = ""))
print(plot(log(
        fert.rate.mcmc[,fertVars.seq[length(fertVars.seq)]:
                                   ncol(fert.rate.mcmc),drop = FALSE])))
dev.off()

###
###
### ##################################################################
### Survival
###
###
survVars.seq <- seq(from = 1, to = ncol(surv.prop.mcmc$female)
                    ,by = 4)

## FEMALES
for(i in 1:(length(survVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logitsurv_traceplot20130612srb10_", i
       ,"_female.png", sep = ""))
  print(plot(logit(surv.prop.mcmc$female[,survVars.seq[i]:
                                            (survVars.seq[i]+3),drop = FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logitsurv_traceplot20130612srb10_",
       length(survVars.seq), "_female.png", sep = ""))
print(plot(logit(
        surv.prop.mcmc$female[,survVars.seq[length(survVars.seq)]:
                                   ncol(surv.prop.mcmc$female),drop = FALSE])))
dev.off()

## MALE
for(i in 1:(length(survVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logitsurv_traceplot20130612srb10_", i, "_male.png"
       , sep = ""))
  print(plot(logit(surv.prop.mcmc$male[,survVars.seq[i]:
                                            (survVars.seq[i]+3),drop = FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logitsurv_traceplot20130612srb10_",
       length(survVars.seq), "_male.png", sep = ""))
print(plot(logit(
        surv.prop.mcmc$male[,survVars.seq[length(survVars.seq)]:
                                   ncol(surv.prop.mcmc$male),drop = FALSE])))
dev.off()

###
###
### ##################################################################
### Migration
###
###
migVars.seq <- seq(from = 1, to = ncol(mig.prop.mcmc$female)
                    ,by = 4)

## FEMALE
for(i in 1:(length(migVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_mig_traceplot20130612srb10_", i, "_female.png", sep = ""))
  print(plot(mig.prop.mcmc$female[,migVars.seq[i]:
                                            (migVars.seq[i]+3),drop = FALSE]))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_mig_traceplot20130612srb10_",
       length(migVars.seq), "_female.png", sep = ""))
print(plot(
        mig.prop.mcmc$female[,migVars.seq[length(migVars.seq)]:
                                   ncol(mig.prop.mcmc$female),drop = FALSE]))
dev.off()

## MALE
for(i in 1:(length(migVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_mig_traceplot20130612srb10_", i, "_male.png", sep = ""))
  print(plot(mig.prop.mcmc$male[,migVars.seq[i]:
                                            (migVars.seq[i]+3),drop = FALSE]))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_mig_traceplot20130612srb10_",
       length(migVars.seq), "_male.png", sep = ""))
print(plot(
        mig.prop.mcmc$male[,migVars.seq[length(migVars.seq)]:
                                   ncol(mig.prop.mcmc$male),drop = FALSE]))
dev.off()

###
###
### ##################################################################
### Baseline
###
###
baselineVars.seq <- seq(from = 1, to = ncol(baseline.count.mcmc$female)
                    ,by = 4)

## FEMALE
for(i in 1:(length(baselineVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logbaseline_traceplot20130612srb10_", i, "_female.png", sep = ""))
  print(plot(log(baseline.count.mcmc$female[,baselineVars.seq[i]:
                                            (baselineVars.seq[i]+3),drop = FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logbaseline_traceplot20130612srb10_",
       length(baselineVars.seq), "_female.png", sep = ""))
print(plot(
        log(baseline.count.mcmc$female[,baselineVars.seq[length(baselineVars.seq)]:
                                   ncol(baseline.count.mcmc$female),drop = FALSE])))
dev.off()

## MALE
for(i in 1:(length(baselineVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_logbaseline_traceplot20130612srb10_", i, "_male.png", sep = ""))
  print(plot(log(baseline.count.mcmc$male[,baselineVars.seq[i]:
                                            (baselineVars.seq[i]+3),drop = FALSE])))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_logbaseline_traceplot20130612srb10_",
       length(baselineVars.seq), "_male.png", sep = ""))
print(plot(
        log(baseline.count.mcmc$male[,baselineVars.seq[length(baselineVars.seq)]:
                                   ncol(baseline.count.mcmc$male),drop = FALSE])))
dev.off()

###
###
### ##################################################################
### lx
###
###
lxVars.seq <- seq(from = 1, to = ncol(lx.mcmc$female)
                    ,by = 4)

## FEMALE
for(i in 1:(length(lxVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_lx_traceplot20130612srb10_", i, "_female.png", sep = ""))
  print(plot(lx.mcmc$female[,lxVars.seq[i]:
                                            (lxVars.seq[i]+3),drop = FALSE]))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_lx_traceplot20130612srb10_",
       length(lxVars.seq), "_female.png", sep = ""))
print(plot(
        lx.mcmc$female[,lxVars.seq[length(lxVars.seq)]:
                                   ncol(lx.mcmc$female),drop = FALSE]))
dev.off()

## MALE
for(i in 1:(length(lxVars.seq)-1)) {
  png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
      ,file = paste("results/plots/traceplots/Thai_lx_traceplot20130612srb10_", i, "_male.png", sep = ""))
  print(plot(lx.mcmc$male[,lxVars.seq[i]:
                                            (lxVars.seq[i]+3),drop = FALSE]))
  dev.off()
}
png(width = 9, height = 7, units = "in", res = 360, pointsize = 10
    ,file = paste("results/plots/traceplots/Thai_lx_traceplot20130612srb10_",
       length(lxVars.seq), "_male.png", sep = ""))
print(plot(
        lx.mcmc$male[,lxVars.seq[length(lxVars.seq)]:
                                   ncol(lx.mcmc$male),drop = FALSE]))
dev.off()
