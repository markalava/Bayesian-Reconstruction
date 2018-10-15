################################################################################
###
###  TITLE:             4_plot_results.R
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
###  Plot prior and posterior confidence intervals for population reconstruction.
###
###  Run at least 1_estimation.R and 1_sample_from_prior.R first.
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
library(lattice)
library(coda)
library(gdata)
library(reshape2); dmelt <- function(...) reshape2:::melt.data.frame(...)

data.path <- "data"
prog.path <- "R_Code"
results.path <- "results"

if(!file.exists(file.path("results", "plots"))) {
    dir.create(file.path("results", "plots"), recursive = TRUE)
    }

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
### ** Constants

## Quantiles to calculate

quants.to.plot <- c(0.025, 0.05, 0.5, 0.95, 0.975)


## Use file.backed data?

fb <- TRUE


### ############################################################################
### ** MCMC output

for(x in dir(results.path, pattern = "\\.RData")) {
    load(file.path(results.path, x), verbose = TRUE)
    }


### ############################################################################
### ** Separation Factors

load(file = file.path(data.path, "thai_sep_factors.RData"))


### ############################################################################
### ** Sample From Prior

## Make one list for each parameter
Thai.surv.prop.prior <-
    list(female = Thai.surv.prop.prior.female
         ,male =  Thai.surv.prop.prior.male
         )
rm(Thai.surv.prop.prior.female
   ,Thai.surv.prop.prior.male)

Thai.mig.prop.prior <-
    list(female = Thai.mig.prop.prior.female
         ,male =  Thai.mig.prop.prior.male
         )
rm(Thai.mig.prop.prior.female
   ,Thai.mig.prop.prior.male)

Thai.baseline.count.prior <-
    list(female = Thai.baseline.count.prior.female
         ,male =  Thai.baseline.count.prior.male
         )
rm(Thai.baseline.count.prior.female
   ,Thai.baseline.count.prior.male)

Thai.pop.count.prior <-
    list(female = Thai.pop.count.prior.female
         ,male =  Thai.pop.count.prior.male
         )
rm(Thai.pop.count.prior.female
   ,Thai.pop.count.prior.male)


################################################################################
### * VITAL RATE AND POSTERIOR PREDICTIVE INTERVALS

################################################################################
### ** Sex ratio at birth

m <- ThaiMcmc$srb.mcmc

#
#.. Calculate posterior quantiles
#

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#

years <- colnames(m)

#
#.. Prepare data sets
#
alpha <- ThaiMcmc$fixed.params$alpha.srb
beta <- ThaiMcmc$fixed.params$beta.srb

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("srb.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

meas <- melt(ThaiMcmc$fixed.params$mean.srb)[,-1]
colnames(meas) <- c("year", "srb.50pctl")

upperQ95 <- exp(log(meas$srb.50pctl) + qt(p = 1-0.975, df = 2 * alpha
                                          , lower.tail = FALSE) *
                sqrt(beta/alpha)
                )
lowerQ95 <- exp(log(meas$srb.50pctl) - qt(p = 1-0.975, df = 2 * alpha
                                          , lower.tail = FALSE) *
                sqrt(beta/alpha)
                )
upperQ90 <- exp(log(meas$srb.50pctl) + qt(p = 1-0.9, df = 2 * alpha
                                          , lower.tail = FALSE) *
                sqrt(beta/alpha)
                )
lowerQ90 <- exp(log(meas$srb.50pctl) - qt(p = 1-0.9, df = 2 * alpha
                                          , lower.tail = FALSE) *
                sqrt(beta/alpha)
                )
mmeas.df <-
    data.frame(meas, srb.97.5pctl = melt(upperQ95)$value
               ,srb.95pctl = melt(upperQ90)$value
               ,srb.5pctl = melt(lowerQ90)$value
               ,srb.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est."
               )

plot.df <- rbind(mmeas.df, mqvit.df)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")


## Pr SRB above 1.06

Thai.srb.pr.gt.1.06 <-
    apply(m, 2, function(z, N) sum(z > 1.06) / N, N = nrow(m))

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_srb_priorpost_q95.pdf")
print(
      ggplot(data = plot.df, aes(x = year, y = srb.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = srb.2.5pctl
                      ,ymax = srb.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA) +
      ylab("srb")      )
dev.off()

#
#.. Save for table
#
Thai.srb.output.df <- plot.df
save(Thai.srb.output.df, Thai.srb.pr.gt.1.06
     ,file = file.path("results", "Thai_srb_output.RData")
     )
write.csv(Thai.srb.output.df
          ,file = file.path("results", "Thai_srb_output.csv")
          ,row.names = FALSE
     )


################################################################################
### *** Distribution of Change


################################################################################
### *** Distribution of OLS coefficient

## Look at distribution of change in SRB over the reconstruction period


###
### Posterior
###

srb.yrs <- as.numeric(colnames(ThaiMcmc$srb.mcmc))
srb.n <- ncol(ThaiMcmc$srb.mcmc)
X <- cbind(rep(1, srb.n), srb.yrs)
hat <- solve(t(X) %*% X) %*% t(X)
hat2 <- hat[2,]

m <- apply(ThaiMcmc$srb.mcmc, 1, function(z) {
    hat2 %*% matrix(z, ncol = 1)
})

q.vital <- quantile(m, probs = quants.to.plot)
Thai.srb.OLS.slope.q.vital <- q.vital
Thai.srb.OLS.slope.q.vital


## Probability > 0
Thai.srb.OLS.slope.pr.gt.zero <- sum(m > 0) / length(m)
Thai.srb.OLS.slope.pr.gt.zero


## Positive/negative m
thai.srb.pos.slope.0612srb10 <- ThaiMcmc$srb.mcmc[which(m > 0)[c(100,1000,5000)],]
thai.srb.neg.slope.0612srb10 <- ThaiMcmc$srb.mcmc[which(m < 0)[3000],]


## Plot


## Save for table
save(Thai.srb.OLS.slope.q.vital, Thai.srb.OLS.slope.pr.gt.zero
     ,file = file.path("results", "Thai_srb_OLS_slope.RData")
     )


################################################################################
### ** Fertility

m <- ThaiMcmc$fert.rate.mcmc
fert.rows <- ThaiMcmc$alg.params$non.zero.fert.rows

#
#.. Calculate posterior quantiles
#
q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
fert.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
fert.ages.numeric <- as.numeric(gsub("[^0-9]", "", fert.ages))

#
#.. Prepare data sets
#
alpha <- ThaiMcmc$fixed.params$alpha.fert.rate
beta <- ThaiMcmc$fixed.params$beta.fert.rate

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("fert.rate.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

meas <- melt(ThaiMcmc$fixed.params$mean.fert.rate[fert.rows,])
meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "fert.rate.50pctl"))

upperQ95 <- exp(log(meas$fert.rate.50pctl) + qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
              )
lowerQ95 <- exp(log(meas$fert.rate.50pctl) - qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
              )
upperQ90 <- exp(log(meas$fert.rate.50pctl) + qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
              )
lowerQ90 <- exp(log(meas$fert.rate.50pctl) - qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
              )
mmeas.df <-
    data.frame(meas, fert.rate.97.5pctl = melt(upperQ95)$value
               ,fert.rate.95pctl = melt(upperQ90)$value
               ,fert.rate.5pctl = melt(lowerQ90)$value
               ,fert.rate.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est."
               )

plot.df <- rbind(mmeas.df, mqvit.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_fert_rate_priorpost_q95.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = fert.rate.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = fert.rate.2.5pctl
                      ,ymax = fert.rate.97.5pctl, fill = legend), alpha = 0.15
                  ,color = NA) +
      ylab("fert. rate")
      )
dev.off()

#
#.. Save for table
#
Thai.fert.rate.output.df <- plot.df
save(Thai.fert.rate.output.df
     ,file = file.path("results", "Thai_fert_rate_output.RData")
     )
write.csv(Thai.fert.rate.output.df
          ,file = file.path("results", "Thai_fert_rate_output.csv")
          ,row.names = FALSE
     )


### ############################################################################
### ** Mortality

################################################################################
### *** Mortality rates

################################################################################
### **** Female

###
### Posterior

##
##
##--- Calculate quantiles
##
m <- Thai.mort.rate.mcmc[[1]][,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mort.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mort.ages.numeric <- as.numeric(gsub("[^0-9]", "", mort.ages))

#
#.. Posterior quantiles
#
mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("mort.rate.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )


###
### Prior

##.. Prior quantiles by sampling Not possible to do this without the age
##   structure because if you change the nLx (by re-sampling the nSx) the nax
##   from the initial life table may not necessarily yield a valid stationary
##   population. Must sample from the /whole/ prior; fertility, survival,
##   baseline and migration counts all.

m <- Thai.mort.rate.mcmc.prior[[1]][,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mort.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mort.ages.numeric <- as.numeric(gsub("[^0-9]", "", mort.ages))

#
#.. Posterior quantiles
#
mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("mort.rate.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )


###
### Plot

#
#.. Prepare datasets
#

plot.df <- rbind(mmeas.df, mqvit.df)

plot.df$age <- as.numeric(as.character(plot.df$age))
plot.df$year <- as.numeric(as.character(plot.df$year))
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
plot.df$sex <- "female"

## Express as per 1,000
plot.df[, grep("pctl", x = colnames(plot.df))] <-
    1000 * plot.df[, grep("pctl", x = colnames(plot.df))]

#
#.. Plot quantiles of under five mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mort_rate_u5_priorpost_q95_female.pdf")
print(
      ggplot(data = subset(plot.df, age < 5), aes(x = year
             , y = mort.rate.50pctl, color = legend
             ## , linetype = sex
             )) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = mort.rate.2.5pctl
                      ,ymax = mort.rate.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("mort. rate (per 1,000)")
      )
dev.off()


#
#.. Plot quantiles of log mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mort_rate_log_priorpost_q95_female.pdf")
print(
      ggplot(data = plot.df, aes(x = age
             , y = log(mort.rate.50pctl), color = legend
             ## , linetype = sex
             )) +
      facet_wrap(~ year) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = log(mort.rate.2.5pctl)
                      ,ymax = log(mort.rate.97.5pctl)
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("log mort. rate (per 1,000)")
      )
dev.off()

#
#.. Save for table
#
Thai.mort.rate.output.female.df <- plot.df
save(Thai.mort.rate.output.female.df
     ,file = file.path("results", "Thai_mort_rate_output_female.RData")
     )
write.csv(Thai.mort.rate.output.female.df
          ,file = file.path("results", "Thai_mort_rate_output_female.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Male

###
### Posterior

##
##
##--- Calculate quantiles
##
m <- Thai.mort.rate.mcmc[[1]][,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mort.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mort.ages.numeric <- as.numeric(gsub("[^0-9]", "", mort.ages))

#
#.. Posterior quantiles
#
mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("mort.rate.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )


###
### Prior

##.. Prior quantiles by sampling Not possible to do this without the age
##   structure because if you change the nLx (by re-sampling the nSx) the nax
##   from the initial life table may not necessarily yield a valid stationary
##   population. Must sample from the /whole/ prior; fertility, survival,
##   baseline and migration counts all.

m <- Thai.mort.rate.mcmc.prior[[1]][,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mort.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mort.ages.numeric <- as.numeric(gsub("[^0-9]", "", mort.ages))

#
#.. Posterior quantiles
#
mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("mort.rate.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )


###
### Plot

#
#.. Prepare datasets
#

plot.df <- rbind(mmeas.df, mqvit.df)

plot.df$age <- as.numeric(as.character(plot.df$age))
plot.df$year <- as.numeric(as.character(plot.df$year))
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")
plot.df$sex <- "male"

## Express as per 1,000
plot.df[, grep("pctl", x = colnames(plot.df))] <-
    1000 * plot.df[, grep("pctl", x = colnames(plot.df))]

#
#.. Plot quantiles of under five mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mort_rate_u5_priorpost_q95_male.pdf")
print(
      ggplot(data = subset(plot.df, age < 5), aes(x = year
             , y = mort.rate.50pctl, color = legend
             ## , linetype = sex
             )) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = mort.rate.2.5pctl
                      ,ymax = mort.rate.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("mort. rate (per 1,000)")
      )
dev.off()


#
#.. Plot quantiles of log mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mort_rate_log_priorpost_q95_male.pdf")
print(
      ggplot(data = plot.df, aes(x = age
             , y = log(mort.rate.50pctl), color = legend
             ## , linetype = sex
             )) +
      facet_wrap(~ year) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = log(mort.rate.2.5pctl)
                      ,ymax = log(mort.rate.97.5pctl)
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("log mort. rate (per 1,000)")
      )
dev.off()

#
#.. Save for table
#
Thai.mort.rate.output.male.df <- plot.df
save(Thai.mort.rate.output.male.df
     ,file = file.path("results", "Thai_mort_rate_output_male.RData")
     )
write.csv(Thai.mort.rate.output.male.df
          ,file = file.path("results", "Thai_mort_rate_output_male.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Both

Thai.mort.rate.output.both.df <-
    rbind(data.frame(Thai.mort.rate.output.female.df)
          ,data.frame(Thai.mort.rate.output.male.df)
          )

pdf(width = 9, height = 9, file = "results/plots/Thai_mort_rate_log_priorpost_q95_both.pdf")
print(
      ggplot(data = Thai.mort.rate.output.both.df
             ,aes(x = age, y = log(mort.rate.50pctl), color = legend, linetype = sex
                 )) +
      facet_wrap(~ year) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = log(mort.rate.2.5pctl), ymax = log(mort.rate.97.5pctl)
                      ,fill = legend), alpha = 0.15) +
      ylab("mort. rate (per 1,000)")
      )
dev.off()


### ############################################################################
### ** Migration

### ############################################################################
### *** Counts

### ############################################################################
### **** Females

## Counts generated from the posterior are for the whole five-year
## intervals, so divide by 5.

m <- Thai.mig.count.mcmc$female[,] / 5


###
### Calculate posterior quantiles

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))


## Make ages and years

colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))


## Prepare data sets

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("mig.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )
mqvit.df$age <- as.numeric(levels(mqvit.df$age)[mqvit.df$age])


###
### Priors

## DERIVATION:
##
## g_{a,t} ~ t(g_{a,t}^*, 2\alpha, \beta/\alpha) so
## g_{a,t}^* * n_{a,t}^* ~ t(g_{a,t}^*n_{a,t}^*, 2\alpha, (n_{a,t}^*)^2(\beta/\alpha))

alpha <- ThaiMcmc$fixed.params$alpha.mig.prop
beta <- ThaiMcmc$fixed.params$beta.mig.prop

proj.years <-
    colnames(ThaiMcmc$fixed.params$mean.fert.rate)
census.years <-
   c(colnames(ThaiMcmc$fixed.params$mean.pop.data$female)
      ,colnames(ThaiMcmc$fixed.params$mean.baseline.count$female)
      )
cyrs.proj <- census.years[census.years <= max(proj.years)]

meas <-
    melt(ThaiMcmc$fixed.params$mean.mig.prop$female[,colnames(ThaiMcmc$fixed.params$mean.mig.prop$female) %in% cyrs.proj] *
         cbind(ThaiMcmc$fixed.params$mean.baseline.count$female
               ,ThaiMcmc$fixed.params$mean.pop.data$female
               )[,cyrs.proj]
         )

meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "mig.count.50pctl"))
meas$age <- as.numeric(unlist(strsplit(as.character(meas$age), "[^0-9]")))

meas$n.star <-
    melt(cbind(ThaiMcmc$fixed.params$mean.baseline.count$female
               ,ThaiMcmc$fixed.params$mean.pop.data$female
               )[,cyrs.proj]
         )$value

upperQ95 <- meas$mig.count.50pctl + qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

lowerQ95 <- meas$mig.count.50pctl - qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

upperQ90 <- meas$mig.count.50pctl + qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

lowerQ90 <- meas$mig.count.50pctl - qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

mmeas.df <-
    data.frame(subset(meas, select = -n.star)
               ,mig.count.97.5pctl = melt(upperQ95)$value
               ,mig.count.95pctl = melt(upperQ90)$value
               ,mig.count.5pctl = melt(lowerQ90)$value
               ,mig.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est. (census fixed)"
               )


###
### Using sample from the prior

m <- Thai.mig.count.mcmc.prior$female[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

#
#.. Posterior quantiles
#
mmeas.df.check <- t(q.vital)
colnames(mmeas.df.check) <-
    paste("mig.count.", prettyNum(as.numeric(colnames(mmeas.df.check)) * 100)
          , "pctl", sep = "")
mmeas.df.check <-
    data.frame(mmeas.df.check
               ,age = sapply(strsplit(rownames(mmeas.df.check), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df.check), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )
mmeas.df.check$age <- as.numeric(levels(mmeas.df.check$age)[mmeas.df.check$age])


###
### Plot

plot.df <- rbind(mmeas.df.check, mmeas.df, mqvit.df)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est. (census fixed)")
plot.df$age <- as.numeric(plot.df$age)

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mig_count_priorpost_q95_female.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = mig.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = mig.count.2.5pctl, ymax = mig.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.2
                  ,color = NA) +
      ylab("mig count")
      )
dev.off()

#
#.. Save for table
#
Thai.mig.count.output.female.df <- plot.df
save(Thai.mig.count.output.female.df
     ,file = file.path("results", "Thai_mig_count_output_female.RData")
     )
write.csv(Thai.mig.count.output.female.df
          ,file = file.path("results", "Thai_mig_count_output_female.csv")
          ,row.names = FALSE
     )


### ############################################################################
### **** Males

## Counts generated from the posterior are for the whole five-year
## intervals, so divide by 5.

m <- Thai.mig.count.mcmc$male[,] / 5


###
### Calculate posterior quantiles

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))


## Make ages and years

colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))


## Prepare data sets

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("mig.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )
mqvit.df$age <- as.numeric(levels(mqvit.df$age)[mqvit.df$age])


###
### Priors

## DERIVATION:
##
## g_{a,t} ~ t(g_{a,t}^*, 2\alpha, \beta/\alpha) so
## g_{a,t}^* * n_{a,t}^* ~ t(g_{a,t}^*n_{a,t}^*, 2\alpha, (n_{a,t}^*)^2(\beta/\alpha))

alpha <- ThaiMcmc$fixed.params$alpha.mig.prop
beta <- ThaiMcmc$fixed.params$beta.mig.prop

proj.years <-
    colnames(ThaiMcmc$fixed.params$mean.fert.rate)
census.years <-
   c(colnames(ThaiMcmc$fixed.params$mean.pop.data$male)
      ,colnames(ThaiMcmc$fixed.params$mean.baseline.count$male)
      )
cyrs.proj <- census.years[census.years <= max(proj.years)]

meas <-
    melt(ThaiMcmc$fixed.params$mean.mig.prop$male[,colnames(ThaiMcmc$fixed.params$mean.mig.prop$male) %in% cyrs.proj] *
         cbind(ThaiMcmc$fixed.params$mean.baseline.count$male
               ,ThaiMcmc$fixed.params$mean.pop.data$male
               )[,cyrs.proj]
         )

meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "mig.count.50pctl"))
meas$age <- as.numeric(unlist(strsplit(as.character(meas$age), "[^0-9]")))

meas$n.star <-
    melt(cbind(ThaiMcmc$fixed.params$mean.baseline.count$male
               ,ThaiMcmc$fixed.params$mean.pop.data$male
               )[,cyrs.proj]
         )$value

upperQ95 <- meas$mig.count.50pctl + qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

lowerQ95 <- meas$mig.count.50pctl - qt(p = 1-0.975, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

upperQ90 <- meas$mig.count.50pctl + qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

lowerQ90 <- meas$mig.count.50pctl - qt(p = 1-0.9, df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha) * meas$n.star

mmeas.df <-
    data.frame(subset(meas, select = -n.star)
               ,mig.count.97.5pctl = melt(upperQ95)$value
               ,mig.count.95pctl = melt(upperQ90)$value
               ,mig.count.5pctl = melt(lowerQ90)$value
               ,mig.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est. (census fixed)"
               )


###
### Using sample from the prior

m <- Thai.mig.count.mcmc.prior$male[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

#
#.. Posterior quantiles
#
mmeas.df.check <- t(q.vital)
colnames(mmeas.df.check) <-
    paste("mig.count.", prettyNum(as.numeric(colnames(mmeas.df.check)) * 100)
          , "pctl", sep = "")
mmeas.df.check <-
    data.frame(mmeas.df.check
               ,age = sapply(strsplit(rownames(mmeas.df.check), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df.check), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )
mmeas.df.check$age <- as.numeric(levels(mmeas.df.check$age)[mmeas.df.check$age])


###
### Plot

plot.df <- rbind(mmeas.df.check, mmeas.df, mqvit.df)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est. (census fixed)")
plot.df$age <- as.numeric(plot.df$age)

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_mig_count_priorpost_q95_male.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = mig.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = mig.count.2.5pctl, ymax = mig.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.2
                  ,color = NA) +
      ylab("mig count")
      )
dev.off()

#
#.. Save for table
#
Thai.mig.count.output.male.df <- plot.df
save(Thai.mig.count.output.male.df
     ,file = file.path("results", "Thai_mig_count_output_male.RData")
     )
write.csv(Thai.mig.count.output.male.df
          ,file = file.path("results", "Thai_mig_count_output_male.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Both

Thai.mig.count.output.both.df <-
    rbind(data.frame(Thai.mig.count.output.female.df, sex = "female")
          ,data.frame(Thai.mig.count.output.male.df, sex = "male")
          )

pdf(width = 9, height = 9, file = "results/plots/Thai_mig_count_priorpost_q95_both.pdf")
print(
      ggplot(data = Thai.mig.count.output.both.df
             ,aes(x = age, y = mig.count.50pctl, color = legend, linetype = sex
                 )) +
      facet_wrap(~ year) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = mig.count.2.5pctl, ymax = mig.count.97.5pctl
                      ,fill = legend), alpha = 0.15) +
      ylab("mig. count")
      )
dev.off()


### ############################################################################
### ** Baseline

################################################################################
### *** Female

m <- ThaiMcmc$baseline.count.mcmc$female

#
#.. Calculate posterior quantiles
#
q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
baseline.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
baseline.ages.numeric <- as.numeric(gsub("[^0-9]", "", baseline.ages))

#
#.. Prepare data sets
#
alpha <- ThaiMcmc$fixed.params$alpha.population.count
beta <- ThaiMcmc$fixed.params$beta.population.count

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("baseline.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

meas <- melt(ThaiMcmc$fixed.params$mean.baseline.count$female)
meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "baseline.count.50pctl"))
meas$age <-
    sapply(strsplit(as.character(meas$age), split = "[^0-9]"), "[[", 1)

upperQ95 <- exp(log(meas$baseline.count.50pctl) + qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ95 <- exp(log(meas$baseline.count.50pctl) - qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
upperQ90 <- exp(log(meas$baseline.count.50pctl) + qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ90 <- exp(log(meas$baseline.count.50pctl) - qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
mmeas.df <-
    data.frame(meas, baseline.count.97.5pctl = melt(upperQ95)$value
               ,baseline.count.95pctl = melt(upperQ90)$value
               ,baseline.count.5pctl = melt(lowerQ90)$value
               ,baseline.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est."
               )

plot.df <- rbind(mmeas.df, mqvit.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_baseline_count_priorpost_q95_female.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = baseline.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = baseline.count.2.5pctl, ymax = baseline.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.15
                  ,color = NA) +
      ylab("fert. rate")
      )
dev.off()

#
#.. Save for table
#
Thai.baseline.count.output.female.df <- plot.df
save(Thai.baseline.count.output.female.df
     ,file = file.path("results", "Thai_baseline_count_output_female.RData")
     )
write.csv(Thai.baseline.count.output.female.df
          ,file = file.path("results", "Thai_baseline_count_output_female.csv")
          ,row.names = FALSE
     )


################################################################################
### *** Male

m <- ThaiMcmc$baseline.count.mcmc$male

#
#.. Calculate posterior quantiles
#
q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
baseline.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
baseline.ages.numeric <- as.numeric(gsub("[^0-9]", "", baseline.ages))

#
#.. Prepare data sets
#
alpha <- ThaiMcmc$fixed.params$alpha.population.count
beta <- ThaiMcmc$fixed.params$beta.population.count

mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("baseline.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )

meas <- melt(ThaiMcmc$fixed.params$mean.baseline.count$male)
meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "baseline.count.50pctl"))
meas$age <-
    sapply(strsplit(as.character(meas$age), split = "[^0-9]"), "[[", 1)

upperQ95 <- exp(log(meas$baseline.count.50pctl) + qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ95 <- exp(log(meas$baseline.count.50pctl) - qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
upperQ90 <- exp(log(meas$baseline.count.50pctl) + qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ90 <- exp(log(meas$baseline.count.50pctl) - qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
mmeas.df <-
    data.frame(meas, baseline.count.97.5pctl = melt(upperQ95)$value
               ,baseline.count.95pctl = melt(upperQ90)$value
               ,baseline.count.5pctl = melt(lowerQ90)$value
               ,baseline.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est."
               )

plot.df <- rbind(mmeas.df, mqvit.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

#
#.. Plot quantiles
#
pdf(width = 9, height = 9, file = "results/plots/Thai_baseline_count_priorpost_q95_male.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = baseline.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = baseline.count.2.5pctl, ymax = baseline.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.15
                  ,color = NA) +
      ylab("fert. rate")
      )
dev.off()

#
#.. Save for table
#
Thai.baseline.count.output.male.df <- plot.df
save(Thai.baseline.count.output.male.df
     ,file = file.path("results", "Thai_baseline_count_output_male.RData")
     )
write.csv(Thai.baseline.count.output.male.df
          ,file = file.path("results", "Thai_baseline_count_output_male.csv")
          ,row.names = FALSE
     )


### ############################################################################
### ** All Population counts

################################################################################
### *** Female

###
### Posterior

m <- ThaiMcmc$lx.mcmc$female

#
#.. Calculate posterior quantiles
#
q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
population.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
population.ages.numeric <- as.numeric(gsub("[^0-9]", "", population.ages))


mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("population.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )
if(is.factor(mqvit.df$age)) mqvit.df$age <- as.numeric(levels(mqvit.df$age)[mqvit.df$age])
if(is.factor(mqvit.df$year)) mqvit.df$year <- as.numeric(levels(mqvit.df$year)[mqvit.df$year])


###
### Prior

m <- Thai.pop.count.prior$female[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("population.count.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )
if(is.factor(mmeas.df$age)) mmeas.df$age <- as.numeric(levels(mmeas.df$age)[mmeas.df$age])
if(is.factor(mmeas.df$year)) mmeas.df$year <- as.numeric(levels(mmeas.df$year)[mmeas.df$year])


###
### Prior (Census fixed)

alpha <- ThaiMcmc$fixed.params$alpha.population.count
beta <- ThaiMcmc$fixed.params$beta.population.count

meas <- melt(ThaiMcmc$fixed.params$mean.pop.data$female)
meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "population.count.50pctl"))
meas$age <-
    sapply(strsplit(as.character(meas$age), split = "[^0-9]"), "[[", 1)

upperQ95 <- exp(log(meas$population.count.50pctl) + qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ95 <- exp(log(meas$population.count.50pctl) - qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
upperQ90 <- exp(log(meas$population.count.50pctl) + qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ90 <- exp(log(meas$population.count.50pctl) - qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
mmeas.cen.fixed.df <-
    data.frame(meas, population.count.97.5pctl = melt(upperQ95)$value
               ,population.count.95pctl = melt(upperQ90)$value
               ,population.count.5pctl = melt(lowerQ90)$value
               ,population.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est. (census fixed)"
               )
if(is.factor(mmeas.cen.fixed.df$age))
    mmeas.cen.fixed.df$age <- as.numeric(levels(mmeas.cen.fixed.df$age)[mmeas.cen.fixed.df$age])
if(is.factor(mmeas.cen.fixed.df$year))
    mmeas.cen.fixed.df$year <- as.numeric(levels(mmeas.cen.fixed.df$year)[mmeas.cen.fixed.df$year])


###
### Plot

plot.df <- rbind(mmeas.df, mmeas.cen.fixed.df, mqvit.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

#
#.. Plot quantiles
#
pdf(width = 9, height = 9
    , file = "results/plots/Thai_population_count_priorpost_q95_female.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = population.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = population.count.2.5pctl, ymax = population.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.15
                  ,color = NA) +
      ylab("population count")
      )
dev.off()

#
#.. Save for table
#
Thai.population.count.output.female.df <- plot.df
save(Thai.population.count.output.female.df
     ,file = file.path("results", "Thai_population_count_output_female.RData")
     )
write.csv(Thai.population.count.output.female.df
          ,file = file.path("results", "Thai_population_count_output_female.csv")
          ,row.names = FALSE
     )

################################################################################
### *** Male

###
### Posterior

m <- ThaiMcmc$lx.mcmc$male

#
#.. Calculate posterior quantiles
#
q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))
population.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
population.ages.numeric <- as.numeric(gsub("[^0-9]", "", population.ages))


mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("population.count.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )
if(is.factor(mqvit.df$age)) mqvit.df$age <-
    as.numeric(levels(mqvit.df$age)[mqvit.df$age])
if(is.factor(mqvit.df$year)) mqvit.df$year <-
    as.numeric(levels(mqvit.df$year)[mqvit.df$year])


###
### Prior

m <- Thai.pop.count.prior$male[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("population.count.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 2)
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )
if(is.factor(mmeas.df$age)) mmeas.df$age <- as.numeric(levels(mmeas.df$age)[mmeas.df$age])
if(is.factor(mmeas.df$year)) mmeas.df$year <- as.numeric(levels(mmeas.df$year)[mmeas.df$year])


###
### Prior (Census fixed)

alpha <- ThaiMcmc$fixed.params$alpha.population.count
beta <- ThaiMcmc$fixed.params$beta.population.count

meas <- melt(ThaiMcmc$fixed.params$mean.pop.data$male)
meas <-
    rename.vars(meas, from = colnames(meas)
                ,to = c("age", "year", "population.count.50pctl"))
meas$age <-
    sapply(strsplit(as.character(meas$age), split = "[^0-9]"), "[[", 1)

upperQ95 <- exp(log(meas$population.count.50pctl) + qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ95 <- exp(log(meas$population.count.50pctl) - qt(p = 1-0.975
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
upperQ90 <- exp(log(meas$population.count.50pctl) + qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
lowerQ90 <- exp(log(meas$population.count.50pctl) - qt(p = 1-0.9
                                                     , df = 2 * alpha
                                      , lower.tail = FALSE) *
                                          sqrt(beta/alpha)
                )
mmeas.cen.fixed.df <-
    data.frame(meas, population.count.97.5pctl = melt(upperQ95)$value
               ,population.count.95pctl = melt(upperQ90)$value
               ,population.count.5pctl = melt(lowerQ90)$value
               ,population.count.2.5pctl = melt(lowerQ95)$value
               ,legend = "init. est. (census fixed)"
               )
if(is.factor(mmeas.cen.fixed.df$age))
    mmeas.cen.fixed.df$age <- as.numeric(levels(mmeas.cen.fixed.df$age)[mmeas.cen.fixed.df$age])
if(is.factor(mmeas.cen.fixed.df$year))
    mmeas.cen.fixed.df$year <- as.numeric(levels(mmeas.cen.fixed.df$year)[mmeas.cen.fixed.df$year])


###
### Plot

plot.df <- rbind(mmeas.df, mmeas.cen.fixed.df, mqvit.df)
plot.df$age <- as.numeric(plot.df$age)
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

#
#.. Plot quantiles
#
pdf(width = 9, height = 9
    , file = "results/plots/Thai_population_count_priorpost_q95_male.pdf")
print(
      ggplot(data = plot.df, aes(x = age, y = population.count.50pctl, color = legend)) +
      facet_wrap(~ year) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = population.count.2.5pctl, ymax = population.count.97.5pctl
                      ,fill = legend)
                  ,alpha = 0.15
                  ,color = NA) +
      ylab("population count")
      )
dev.off()

#
#.. Save for table
#
Thai.population.count.output.male.df <- plot.df
save(Thai.population.count.output.male.df
     ,file = file.path("results", "Thai_population_count_output_male.RData")
     )
write.csv(Thai.population.count.output.male.df
          ,file = file.path("results", "Thai_population_count_output_male.csv")
          ,row.names = FALSE
          )


### ############################################################################
### * VITAL RATE SUMMARY MEASURE POSTERIOR PREDICTIVE INTERVALS

################################################################################
### ** Mortality

### ############################################################################
### *** Life Expectancy

################################################################################
### **** Female

## Posterior by converting posterior quantiles of survival and assuming
## stationary population relation holds

leb.f <- function(z) {
    ## z is a vector of age-specific survival proportions
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }

surv.prop.years <-
    sapply(strsplit(colnames(ThaiMcmc$surv.prop.mcmc$female), "\\."), "[[", 1)

Thai.leb.stationary.female.df <-
    apply(ThaiMcmc$surv.prop.mcmc$female[,], 1, function(z) {
        tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
      })

Thai.leb.stationary.Quantiles <-
    apply(Thai.leb.stationary.female.df, 1, "quantile", probs = quants.to.plot)

Thai.leb.stationary.Quantiles.df <-
    as.data.frame(t(Thai.leb.stationary.Quantiles))

colnames(Thai.leb.stationary.Quantiles.df) <-
    paste("leb."
          , strsplit(colnames(Thai.leb.stationary.Quantiles.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.leb.stationary.Quantiles.df$legend <- "posterior (stationary)"
Thai.leb.stationary.Quantiles.df$year <-
    as.numeric(rownames(Thai.leb.stationary.Quantiles.df))

## For female vs male
Thai.leb.stationary.Quantiles.female.df <-
    Thai.leb.stationary.Quantiles.df


## Prior by converting posterior quantiles of survival and assuming
## stationary population relation holds

Thai.lebPrior.stationary.females.df <-
    apply(Thai.surv.prop.prior$female[,], 1, function(z) {
        tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
      })

Thai.lebPrior.stationary.females.Quantiles <-
    apply(Thai.lebPrior.stationary.females.df, 1, "quantile", probs = quants.to.plot)

Thai.lebPrior.stationary.females.Quantiles.df <-
    as.data.frame(t(Thai.lebPrior.stationary.females.Quantiles))

colnames(Thai.lebPrior.stationary.females.Quantiles.df) <-
    paste("leb."
          , strsplit(colnames(Thai.lebPrior.stationary.females.Quantiles.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.lebPrior.stationary.females.Quantiles.df$legend <- "init. est. (stationary)"
Thai.lebPrior.stationary.females.Quantiles.df$year <-
    as.numeric(rownames(Thai.lebPrior.stationary.females.Quantiles.df))


###
### Plot

#
#.. Prepare datasets
#
Thai.lebTBL.df <-
    rbind(Thai.lebPrior.stationary.females.Quantiles.df, Thai.leb.stationary.Quantiles.df)

#
#.. Plot
#

plot.df <- Thai.lebTBL.df
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est. (stationary)")

pdf("results/plots/Thai_leb_priorpost_q95_female.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = leb.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = leb.2.5pctl
                                     ,ymax = leb.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("life expectancy at birth (years)")
      )
dev.off()

#
#.. Save for table
#
Thai.leb.output.female.df <- plot.df
save(Thai.leb.output.female.df
     ,file = file.path("results", "Thai_leb_output_female.RData")
     )
write.csv(Thai.leb.output.female.df
          ,file = file.path("results", "Thai_leb_output_female.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Male

## Posterior by converting posterior quantiles of survival and assuming
## stationary population relation holds

leb.f <- function(z) {
    ## z is a vector of age-specific survival proportions
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }

surv.prop.years <-
    sapply(strsplit(colnames(ThaiMcmc$surv.prop.mcmc$male), "\\."), "[[", 1)

Thai.leb.stationary.male.df <-
    apply(ThaiMcmc$surv.prop.mcmc$male[,], 1, function(z) {
        tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
      })

Thai.leb.stationary.Quantiles <-
    apply(Thai.leb.stationary.male.df, 1, "quantile", probs = quants.to.plot)

Thai.leb.stationary.Quantiles.df <-
    as.data.frame(t(Thai.leb.stationary.Quantiles))

colnames(Thai.leb.stationary.Quantiles.df) <-
    paste("leb."
          , strsplit(colnames(Thai.leb.stationary.Quantiles.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.leb.stationary.Quantiles.df$legend <- "posterior (stationary)"
Thai.leb.stationary.Quantiles.df$year <-
    as.numeric(rownames(Thai.leb.stationary.Quantiles.df))

## For female vs male
Thai.leb.stationary.Quantiles.male.df <-
    Thai.leb.stationary.Quantiles.df


## Prior by converting posterior quantiles of survival and assuming
## stationary population relation holds

Thai.lebPrior.stationary.males.df <-
    apply(Thai.surv.prop.prior$male[,], 1, function(z) {
        tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
      })

Thai.lebPrior.stationary.males.Quantiles <-
    apply(Thai.lebPrior.stationary.males.df, 1, "quantile", probs = quants.to.plot)

Thai.lebPrior.stationary.males.Quantiles.df <-
    as.data.frame(t(Thai.lebPrior.stationary.males.Quantiles))

colnames(Thai.lebPrior.stationary.males.Quantiles.df) <-
    paste("leb."
          , strsplit(colnames(Thai.lebPrior.stationary.males.Quantiles.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.lebPrior.stationary.males.Quantiles.df$legend <- "init. est. (stationary)"
Thai.lebPrior.stationary.males.Quantiles.df$year <-
    as.numeric(rownames(Thai.lebPrior.stationary.males.Quantiles.df))


###
### Plot

#
#.. Prepare datasets
#
Thai.lebTBL.df <-
    rbind(Thai.lebPrior.stationary.males.Quantiles.df, Thai.leb.stationary.Quantiles.df)

#
#.. Plot
#

plot.df <- Thai.lebTBL.df
plot.df$year <- as.numeric(plot.df$year)
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est. (stationary)")

pdf("results/plots/Thai_leb_priorpost_q95_male.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = leb.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = leb.2.5pctl
                                     ,ymax = leb.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("life expectancy at birth (years)")
      )
dev.off()

#
#.. Save for table
#
Thai.leb.output.male.df <- plot.df
save(Thai.leb.output.male.df
     ,file = file.path("results", "Thai_leb_output_male.RData")
     )
write.csv(Thai.leb.output.male.df
          ,file = file.path("results", "Thai_leb_output_male.csv")
          ,row.names = FALSE
     )


################################################################################
### *** IMR [0--5]

################################################################################
### **** Female

###
### Posterior

##
##
##--- Calculate quantiles
##
m <- Thai.IMR.mcmc$female[]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))


#
#.. Posterior quantiles
#
mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )


###
### Prior

m <- Thai.IMR.mcmc.prior$female[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))

#
#.. Posterior quantiles
#
mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )


###
### Plot

#
#.. Prepare datasets
#

plot.df <- rbind(mmeas.df, mqvit.df)

plot.df$year <- as.numeric(as.character(plot.df$year))
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")


#
#.. Plot quantiles of log mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_IMR_log_priorpost_q95_female.pdf")
print(
      ggplot(data = plot.df, aes(x = year
             , y = log(IMR.50pctl), color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = log(IMR.2.5pctl)
                      ,ymax = log(IMR.97.5pctl)
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("log IMR (0--5)")
      )
dev.off()

#
#.. Save for table
#
Thai.IMR.output.female.df <- plot.df
save(Thai.IMR.output.female.df
     ,file = file.path("results", "Thai_IMR_output_female.RData")
     )
write.csv(Thai.IMR.output.female.df
          ,file = file.path("results", "Thai_IMR_output_female.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Male

###
### Posterior

##
##
##--- Calculate quantiles
##
m <- Thai.IMR.mcmc$male[]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))


#
#.. Posterior quantiles
#
mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )


###
### Prior

m <- Thai.IMR.mcmc.prior$male[,]

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))

#
#.. Posterior quantiles
#
mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "init. est."
               )


###
### Plot

#
#.. Prepare datasets
#

plot.df <- rbind(mmeas.df, mqvit.df)

plot.df$year <- as.numeric(as.character(plot.df$year))
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")


#
#.. Plot quantiles of log mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_IMR_log_priorpost_q95_male.pdf")
print(
      ggplot(data = plot.df, aes(x = year
             , y = log(IMR.50pctl), color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = log(IMR.2.5pctl)
                      ,ymax = log(IMR.97.5pctl)
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("log IMR (0--5)")
      )
dev.off()

#
#.. Save for table
#
Thai.IMR.output.male.df <- plot.df
save(Thai.IMR.output.male.df
     ,file = file.path("results", "Thai_IMR_output_male.RData")
     )
write.csv(Thai.IMR.output.male.df
          ,file = file.path("results", "Thai_IMR_output_male.csv")
          ,row.names = FALSE
     )


################################################################################
### **** FemaleVSmale

### Posterior
###

Thai.IMR.stationary.output.femaleVSmale.df <-
    rbind(transform(subset(Thai.IMR.output.female.df
                            ,subset = legend == "posterior"
                            )
                     ,legend = "female"
                     )
          ,transform(subset(Thai.IMR.output.male.df
                            ,subset = legend == "posterior"
                            )
                     ,legend = "male"
                     )
          )


###
### Plot

#
#.. Plot
#

plot.df <- Thai.IMR.stationary.output.femaleVSmale.df

pdf("results/plots/Thai_IMR_priorpost_q95_femaleVSmale.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = IMR.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = IMR.2.5pctl
                                     ,ymax = IMR.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("life expectancy at birth (years)")
      )
dev.off()

#
#.. Save for table
#
Thai.IMR.output.femaleVSmale.df <- plot.df
save(Thai.IMR.output.femaleVSmale.df
     ,file = file.path("results", "Thai_IMR_output_femaleVSmale.RData")
     )
write.csv(Thai.IMR.output.femaleVSmale.df
          ,file = file.path("results", "Thai_IMR_output_femaleVSmale.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Male to Female Ratio

###
### Posterior

##
##
##--- Calculate quantiles

mf <- Thai.IMR.mcmc$female[]
mm <- Thai.IMR.mcmc$male[]

m <- mm / mf

rm(mf); rm(mm)

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Make ages and years
#
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
year <- unique(sapply(colspl, FUN = function(z) z[1]))


##
## Pr male U5MR < Female
##

## H0 is that male U5MR > female U5MR, i.e., male / fem > 1. So calculate
## probability that male/fem < 1.

Thai.IMR.sex.ratio.lt.1 <-
    apply(m, 2, function(z, N) sum(z < 1) / N, N = nrow(m))
Thai.IMR.sex.ratio.lt.1


#
#.. Posterior quantiles
#
mqvit.df <- t(q.vital)
colnames(mqvit.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
          , "pctl", sep = "")
mqvit.df <-
    data.frame(mqvit.df
               ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "posterior"
               )


###
### Prior

mf <- Thai.IMR.mcmc.prior$female[,]
mm <- Thai.IMR.mcmc.prior$male[,]

m <- mm / mf

rm(mm); rm(mf)

q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

#
#.. Posterior quantiles
#
mmeas.df <- t(q.vital)
colnames(mmeas.df) <-
    paste("IMR.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
          , "pctl", sep = "")
mmeas.df <-
    data.frame(mmeas.df
               ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                ,"[[", 1)
               ,legend = "prior"
               )


###
### Plot

#
#.. Prepare datasets
#

plot.df <- rbind(mmeas.df, mqvit.df)

plot.df$year <- as.numeric(as.character(plot.df$year))
plot.df$legend <- relevel(factor(plot.df$legend), ref = "prior")


#
#.. Plot quantiles of log mortality rate
#
pdf(width = 9, height = 9, file = "results/plots/Thai_IMR_log_priorpost_q95_male_to_female_ratio.pdf")
print(
      ggplot(data = plot.df, aes(x = year
             , y = log(IMR.50pctl), color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_ribbon(aes(ymin = log(IMR.2.5pctl)
                      ,ymax = log(IMR.97.5pctl)
                      ,fill = legend)
                  ,alpha = 0.25
                  ,color = NA) +
      ylab("log IMR (0--5)")
      )
dev.off()

#
#.. Save for table
#
Thai.IMR.output.male.to.female.ratio.df <- plot.df
save(Thai.IMR.output.male.to.female.ratio.df, Thai.IMR.sex.ratio.lt.1
     ,file = file.path("results", "Thai_IMR_output_male_to_female_ratio.RData")
     )
write.csv(Thai.IMR.output.male.to.female.ratio.df
          ,file = file.path("results", "Thai_IMR_output_male_to_female_ratio.csv")
          ,row.names = FALSE
     )


################################################################################
### ** Fertility

### ############################################################################
### *** Total Fertility Rate

###
### Posterior

dn <- list(NULL,
            unique(sapply(strsplit(colnames(ThaiMcmc$fert.rate.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
#
#.. Calculate
#
Thai.tfr <-
  matrix(0, nrow = nrow(ThaiMcmc$fert.rate.mcmc)
         ,ncol = length(dn[[2]])
         ,dimnames = dn
         )

fert.rate.mcmc.colYrs <-
  sapply(strsplit(colnames(ThaiMcmc$fert.rate.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])

for(i in 1:ncol(Thai.tfr)) {
  colYrs.index <- fert.rate.mcmc.colYrs == colnames(Thai.tfr)[i]
  Thai.tfr[,i] <-
    apply(ThaiMcmc$fert.rate.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) sum(z)
        )
}

#
#.. Posterior quantiles
#
Thai.tfrQuant <- apply(Thai.tfr, 2, FUN = function(z)
                        {
                          5 * quantile(z, probs = quants.to.plot)
                        })

Thai.tfrQuant.df <-
    as.data.frame(t(Thai.tfrQuant))

colnames(Thai.tfrQuant.df) <-
    paste("tfr.", strsplit(colnames(Thai.tfrQuant.df), split = "%")
          ,"pctl", sep = "")

Thai.tfrQuant.df$legend = "posterior"
Thai.tfrQuant.df$year = as.numeric(rownames(Thai.tfrQuant.df))


###
### Prior

#
#.. Calculate
#
Thai.tfr.prior <-
  matrix(0, nrow = nrow(Thai.fert.rate.prior)
         ,ncol = length(dn[[2]])
         ,dimnames = dn
         )

fert.rate.mcmc.colYrs <-
  sapply(strsplit(colnames(Thai.fert.rate.prior)
                         ,"\\."), FUN = function(z) z[[1]])

for(i in 1:ncol(Thai.tfr.prior)) {
  colYrs.index <- fert.rate.mcmc.colYrs == colnames(Thai.tfr.prior)[i]
  Thai.tfr.prior[,i] <-
    apply(Thai.fert.rate.prior[,colYrs.index]
        ,1
        ,FUN = function(z) sum(z)
        )
}

#
#.. Init. Est. quantiles
#
Thai.tfr.priorQuant <- apply(Thai.tfr.prior, 2, FUN = function(z)
                        {
                          5 * quantile(z, probs = quants.to.plot)
                        })

Thai.tfr.priorQuant.df <-
    as.data.frame(t(Thai.tfr.priorQuant))

colnames(Thai.tfr.priorQuant.df) <-
    paste("tfr.", strsplit(colnames(Thai.tfr.priorQuant.df), split = "%")
          ,"pctl", sep = "")

Thai.tfr.priorQuant.df$legend = "init. est."
Thai.tfr.priorQuant.df$year = as.numeric(rownames(Thai.tfr.priorQuant.df))


###
### Plot

#
#.. Prepare datasets
#
Thai.tfrTBL.df <- rbind(Thai.tfrQuant.df, Thai.tfr.priorQuant.df)

#
#.. Plot
#

plot.df <- Thai.tfrTBL.df
plot.df$legend <- relevel(factor(plot.df$legend), ref = "init. est.")

pdf("results/plots/Thai_tfr_priorpost_q95.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = tfr.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_point() + geom_ribbon(aes(ymin = tfr.2.5pctl
                                     ,ymax = tfr.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("total fertility rate")
      )
dev.off()

#
#.. Save for table
#
Thai.tfr.output.df <- plot.df
save(Thai.tfr.output.df
     ,file = file.path("results", "Thai_tfr_output.RData")
     )
write.csv(Thai.tfr.output.df
          ,file = file.path("results", "Thai_tfr_output.csv")
          ,row.names = FALSE
     )


################################################################################
### ** Migration

### ############################################################################
### *** Total Net Migration

################################################################################
### **** Female

###
### Posterior

Thai.total.mig.count.mcmcQuant <- apply(Thai.total.mig.count.mcmc$female[], 2
                                         , FUN = function(z)
                        {
                          quantile(z, quants.to.plot)
                        })

Thai.total.mig.count.mcmcQuant.df <-
    as.data.frame(t(Thai.total.mig.count.mcmcQuant))

colnames(Thai.total.mig.count.mcmcQuant.df) <-
    paste("total.mig.count.", strsplit(colnames(Thai.total.mig.count.mcmcQuant.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.total.mig.count.mcmcQuant.df$legend <- "posterior"
Thai.total.mig.count.mcmcQuant.df$year <-
    as.numeric(rownames(Thai.total.mig.count.mcmcQuant.df))


###
### Prior

Thai.total.mig.count.mcmc.input <-
    apply(Thai.total.mig.count.mcmc.prior$female[], 2
                                         , FUN = function(z)
                        {
                          quantile(z, quants.to.plot)
                        })

Thai.total.mig.count.mcmc.input.df <-
    as.data.frame(t(Thai.total.mig.count.mcmc.input))

colnames(Thai.total.mig.count.mcmc.input.df) <-
    paste("total.mig.count.", strsplit(colnames(Thai.total.mig.count.mcmc.input.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.total.mig.count.mcmc.input.df$legend <- "init. est."
Thai.total.mig.count.mcmc.input.df$year <-
    as.numeric(rownames(Thai.total.mig.count.mcmc.input.df))


###
### Plot

Thai.total.mig.count.mcmcTBL.df <-
    rbind(Thai.total.mig.count.mcmcQuant.df, Thai.total.mig.count.mcmc.input.df)

## Dataset
plot.df1 <- Thai.total.mig.count.mcmcTBL.df
plot.df1$legend <- relevel(factor(plot.df1$legend), ref = "init. est.")

## init. est. and posteriors are for the whole five year period. Divide by 5 to
## get average annual migration proportion.

plot.df <- plot.df1

plot.df[plot.df1$legend %in% c("init. est.", "posterior")
        ,grep(pattern = "pctl", x = colnames(plot.df))] <-
    plot.df[plot.df1$legend %in% c("init. est.", "posterior")
        ,grep(pattern = "pctl", x = colnames(plot.df))] / 5

## As 1000s
plot.df[, grep(pattern = "pctl", x = colnames(plot.df))] <-
    plot.df[, grep(pattern = "pctl", x = colnames(plot.df))] / 1E3


pdf("results/plots/Thai_total_mig_count_postonly_q95_female.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year
             , y = total.mig.count.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_point() + geom_ribbon(aes(ymin = total.mig.count.2.5pctl
                                     ,ymax = total.mig.count.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("net number of migrants (000s)")
      )
dev.off()

##
##.. Save for table
##
Thai.total.mig.count.output.female.df <- plot.df
save(Thai.total.mig.count.output.female.df
     ,file = file.path("results", "Thai_total_mig_count_output_female.RData")
     )
write.csv(Thai.total.mig.count.output.female.df
          ,file = file.path("results", "Thai_total_mig_count_output_female.csv")
          ,row.names = FALSE
     )


################################################################################
### **** Male

###
### Posterior

Thai.total.mig.count.mcmcQuant <- apply(Thai.total.mig.count.mcmc$male[], 2
                                         , FUN = function(z)
                        {
                          quantile(z, quants.to.plot)
                        })

Thai.total.mig.count.mcmcQuant.df <-
    as.data.frame(t(Thai.total.mig.count.mcmcQuant))

colnames(Thai.total.mig.count.mcmcQuant.df) <-
    paste("total.mig.count.", strsplit(colnames(Thai.total.mig.count.mcmcQuant.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.total.mig.count.mcmcQuant.df$legend <- "posterior"
Thai.total.mig.count.mcmcQuant.df$year <-
    as.numeric(rownames(Thai.total.mig.count.mcmcQuant.df))


###
### Prior

Thai.total.mig.count.mcmc.input <-
    apply(Thai.total.mig.count.mcmc.prior$male[], 2
                                         , FUN = function(z)
                        {
                          quantile(z, quants.to.plot)
                        })

Thai.total.mig.count.mcmc.input.df <-
    as.data.frame(t(Thai.total.mig.count.mcmc.input))

colnames(Thai.total.mig.count.mcmc.input.df) <-
    paste("total.mig.count.", strsplit(colnames(Thai.total.mig.count.mcmc.input.df)
                                            , split = "%")
          ,"pctl", sep = "")

Thai.total.mig.count.mcmc.input.df$legend <- "init. est."
Thai.total.mig.count.mcmc.input.df$year <-
    as.numeric(rownames(Thai.total.mig.count.mcmc.input.df))


###
### Plot

Thai.total.mig.count.mcmcTBL.df <-
    rbind(Thai.total.mig.count.mcmcQuant.df, Thai.total.mig.count.mcmc.input.df)

## Dataset
plot.df1 <- Thai.total.mig.count.mcmcTBL.df
plot.df1$legend <- relevel(factor(plot.df1$legend), ref = "init. est.")

## init. est. and posteriors are for the whole five year period. Divide by 5 to
## get average annual migration proportion.

plot.df <- plot.df1

plot.df[plot.df1$legend %in% c("init. est.", "posterior")
        ,grep(pattern = "pctl", x = colnames(plot.df))] <-
    plot.df[plot.df1$legend %in% c("init. est.", "posterior")
        ,grep(pattern = "pctl", x = colnames(plot.df))] / 5

## As 1000s
plot.df[, grep(pattern = "pctl", x = colnames(plot.df))] <-
    plot.df[, grep(pattern = "pctl", x = colnames(plot.df))] / 1E3


pdf("results/plots/Thai_total_mig_count_postonly_q95_male.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year
             , y = total.mig.count.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_point() + geom_ribbon(aes(ymin = total.mig.count.2.5pctl
                                     ,ymax = total.mig.count.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("net number of migrants (000s)")
      )
dev.off()

##
##.. Save for table
##
Thai.total.mig.count.output.male.df <- plot.df
save(Thai.total.mig.count.output.male.df
     ,file = file.path("results", "Thai_total_mig_count_output_male.RData")
     )
write.csv(Thai.total.mig.count.output.male.df
          ,file = file.path("results", "Thai_total_mig_count_output_male.csv")
          ,row.names = FALSE
     )


################################################################################
### ** Pop Sex Ratios

thai.years <- seq(from = 1960, to = 2000, by = 5)
thai.cen.years <- seq(from = 1960, to = 2000, by = 10)


################################################################################
### *** SRTP

###
### Posterior Quantiles
###

year.cols <-
    sapply(strsplit(colnames(ThaiMcmc$lx.mcmc$female), split = "\\.")
           ,"[[", 1)

uyc <- unique(year.cols)

fsums <-
    matrix(nrow = nrow(ThaiMcmc$lx.mcmc$female), ncol = length(uyc)
           ,dimnames = list(NULL, uyc)
           )
for(y in uyc) fsums[,y] <- rowSums(ThaiMcmc$lx.mcmc$female[,y == year.cols])
mf <- cbind(rowSums(ThaiMcmc$baseline.count.mcmc$female)
           ,fsums
            )

msums <-
    matrix(nrow = nrow(ThaiMcmc$lx.mcmc$male), ncol = length(uyc)
           ,dimnames = list(NULL, uyc)
           )
for(y in uyc) msums[,y] <- rowSums(ThaiMcmc$lx.mcmc$male[,y == year.cols])
mm <- cbind(rowSums(ThaiMcmc$baseline.count.mcmc$male)
           ,msums
            )

Thai.popRatio <- mm / mf

Thai.popRatio.Quantiles.df <-
    apply(Thai.popRatio, 2, "quantile", probs = quants.to.plot)
colnames(Thai.popRatio.Quantiles.df) <- thai.years

rm(mf, mm)


###
### Prior
###

year.cols <-
    sapply(strsplit(colnames(Thai.pop.count.prior[[1]]), split = "\\.")
           ,"[[", 1)

uyc <- unique(year.cols)

fsums <-
    matrix(nrow = nrow(Thai.pop.count.prior$female), ncol = length(uyc)
           ,dimnames = list(NULL, uyc)
           )
for(y in uyc) {
    ## message(y)
    ## message(y == year.cols)
    fsums[,y] <- rowSums(Thai.pop.count.prior$female[,y == year.cols])
}
mf <- cbind(rowSums(Thai.baseline.count.prior$female)
           ,fsums
            )

msums <-
    matrix(nrow = nrow(Thai.pop.count.prior$male), ncol = length(uyc)
           ,dimnames = list(NULL, uyc)
           )
for(y in uyc) msums[,y] <- rowSums(Thai.pop.count.prior$male[,y == year.cols])
mm <- cbind(rowSums(Thai.baseline.count.prior$male)
           ,msums
            )

m <- mm / mf

Thai.popRatio.prior.Quantiles.df <-
    apply(m, 2, "quantile", probs = quants.to.plot)
colnames(Thai.popRatio.prior.Quantiles.df) <- thai.years

rm(m, mf, mm)


###
### Census
###

Thai.cen.popRatio.df <-
    c(sum(ThaiMcmc$fixed.params$mean.baseline.count$male)
      ,colSums(ThaiMcmc$fixed.params$mean.pop.data$male)) /
    c(sum(ThaiMcmc$fixed.params$mean.baseline.count$female)
      ,colSums(ThaiMcmc$fixed.params$mean.pop.data$female))

names(Thai.cen.popRatio.df) <- thai.cen.years


###
### Plot
###

Thai.popRatio.df <-
    rbind(data.frame(t(Thai.popRatio.Quantiles.df), legend = "posterior"
                     ,check.names = FALSE)
          ,data.frame(t(Thai.popRatio.prior.Quantiles.df), legend = "prior"
                     ,check.names = FALSE)
          )

colnames(Thai.popRatio.df)[1:5] <-
    paste("srtp."
          ,gsub(pattern = "%", replacement = "", x = colnames(Thai.popRatio.df[1:5]))
          ,"pctl", sep = "")

Thai.popRatio.df <-
    rbind(Thai.popRatio.df
          ,data.frame(srtp.50pctl = Thai.cen.popRatio.df, legend = "census"
                      ,srtp.2.5pctl = NA, srtp.5pctl = NA, srtp.95pctl = NA
                      ,srtp.97.5pctl = NA
                      ,check.names = FALSE)
          )

Thai.popRatio.df$year <- as.numeric(substr(rownames(Thai.popRatio.df), 1, 4))

Thai.popRatio.df$legend <-
    relevel(factor(Thai.popRatio.df$legend), ref = "prior")

plot.df <- Thai.popRatio.df

#
#.. Plot
#


pdf("results/plots/Thai_srtp_priorpost_q95.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = srtp.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65
                ) +
      geom_point() +
      geom_ribbon(aes(ymin = srtp.2.5pctl
                      ,ymax = srtp.97.5pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("ratio (M/F)") +
      coord_cartesian(ylim = c(0.94, 1.05))
      )
dev.off()

pdf("results/plots/Thai_srtp_priorpost_q90.pdf", width = 7, height = 7)
print(ggplot(data = plot.df, aes(x = year, y = srtp.50pctl, color = legend)) +
      geom_line(size = 1.1, alpha = 0.65) +
      geom_point() +
      geom_ribbon(aes(ymin = srtp.5pctl
                                     ,ymax = srtp.95pctl, fill = legend)
                                 ,alpha = 0.15, color = NA) +
      ylab("ratio (M/F)")
      )
dev.off()

#
#.. Save for table
#
Thai.srtp.output.df <- plot.df
save(Thai.srtp.output.df
     ,file = file.path("results", "Thai_srtp_output.RData")
     )
write.csv(Thai.srtp.output.df
          ,file = file.path("results", "Thai_srtp_output.csv")
          ,row.names = FALSE
     )


################################################################################
