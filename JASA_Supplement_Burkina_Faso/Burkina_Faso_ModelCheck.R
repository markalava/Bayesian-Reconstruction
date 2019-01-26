################################################################################
###
### TITLE:              Burkina_Faso_ModelCheck.R
###
### DATE:               23 May 2012
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Reconstruction of the female population of Burkina Faso
###                     using the t-2 likelihood model described the paper
###                     "Reconstructing Past Populations with Uncertainty from
###                     Fragmentary Data" submitted to Journal of the American
###                     Statistical Association".
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
### This script performs the comparison between the t2- and normal-likelihood
### models described in the 'Model Checking and Sensitivity Analysis' section.
###
###-----------------------------------------------------------------------------
###
### INPUT FILES:
###
### This script expects to find the following files:
### - pop_reconstruction_functions.R
### - burkina_faso_initial_estimates_fertility_rates.csv
### - burkina_faso_initial_estimates_survival_proportions.csv
### - burkina_faso_initial_estimates_migration_proportions.csv
### - burkina_faso_initial_estimates_population_counts.csv
###
### The initial estimates used in the reconstruction are in the
### ..._initial_estimates_....csv files.
###
### The MCMC sampler is in pop_reconstruction_functions.R.
###
###-----------------------------------------------------------------------------
###
### OUTPUTS:
###
### This script creates the following directories:
### - outputs/
### - outputs/tables/
### - outputs/additional_plots/
###
### Tables that appear in the paper are placed in the directory outputs/tables/.
### These directories also contain .Rdata and .csv files containing the data
### values used in the plots.
###
################################################################################


################################################################################
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### * !!! THINGS YOU MUST SET !!!
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

## Change the working directory if needed (i.e., the directory this script is
## going to be run from).
if(identical(.Platform$OS.type, "unix")) {
    home.path <-  "~"
} else if(identical(.Platform$OS.type, "windows")) {
    home.path <- "T:/"
}
setwd(file.path(home.path
                ,"Documents", "PPGp_Working", "TEST_JASA_Scripts"
                ))

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

###
### * SET UP
###
################################################################################


## Libraries
library(coda)
library(reshape)
library(gdata)
library(lattice)

## Source R scripts
source("pop_reconstruction_functions.R")

## Create directories for output
dir.create("outputs")
dir.create("outputs/plots")
dir.create("outputs/additional_plots")
dir.create("outputs/tables")


###
### * Inputs
###
################################################################################

###
### Initial Estimates
###

## Fertility rates
asFertBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_fertility_rates.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )
                )

## Survival proportions
asSurvBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_survival_proportions.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )
                )

## Migration proportions
asMigBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_migration_proportions.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )
                )

## Baseline population
baselineBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_population_counts.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )[,1,drop=FALSE]
                )

## Rest of population counts
censusBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_population_counts.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )[,-1]
                )


###
### Alpha and beta hyper-parameters
###

invGam.params <-
    list(al.f = 1
         ,be.f = 0.0109
         ,al.s = 1
         ,be.s = 0.0109
         ,al.g = 1
         ,be.g = 0.0436
         ,al.n = 1
         ,be.n = 0.003637524
         )

###
### Arguments
###

## Fert rows
fert.rows <- 4:10


## Algorithm parameters
n.iter <- 5000#9E4
burn.in <- 50#1E3


## Start values of mixing parameters
lambda.n <-
    matrix(sqrt(5),
           nrow = nrow(cbind(baselineBKFEM.mat
           , censusBKFEM.mat))
           ,ncol = ncol(cbind(baselineBKFEM.mat
           , censusBKFEM.mat))
           , dimnames = dimnames(cbind(baselineBKFEM.mat
           , censusBKFEM.mat))
           )

## Estimation arguments
BKFem.ModCheck.Arguments <-
    list(## inverse gamma parameters
         al.f = invGam.params$al.f
         ,be.f = invGam.params$be.f
         ,al.s = invGam.params$al.s
         ,be.s = invGam.params$be.s
         ,al.g = invGam.params$al.g
         ,be.g = invGam.params$be.g
         ,al.n = invGam.params$al.n
         ,be.n = invGam.params$be.n

          ## start values of mixing parameters
         ,lambda.n = lambda.n

      ## the rest
      ,n.iter = n.iter, burn.in = burn.in
            ,ccmp.f = "popRecon.ccmp.female"
            ,init.f = asFertBKFEM.mat
            ,init.s = asSurvBKFEM.mat
            ,init.g = asMigBKFEM.mat
            ,init.b = baselineBKFEM.mat
            ,init.sigmasq.f = 5
            ,init.sigmasq.s = 5
            ,init.sigmasq.g = 5
            ,init.sigmasq.n = sqrt(5)
            ,mean.f = asFertBKFEM.mat
            ,mean.s = asSurvBKFEM.mat
            ,mean.g = asMigBKFEM.mat
            ,mean.b = baselineBKFEM.mat
            ,prop.vars = mod.check.prop.vars
            ,pop.data = censusBKFEM.mat
            ,proj.periods = ncol(asFertBKFEM.mat)
            ,age.size = 5
            ,census.columns = c(4,6,8,10)
            ,fert.rows = 4:10
            ,s.tol = 10^(-10)
            ,verb = TRUE
         ,save.log.lhood = TRUE
         ,save.log.post = TRUE
         )


###
### Run Estimation
###

set.seed(1)
BKFem.ModCheck.MCMC <-
    do.call(popReconModCheck.sampler, args = BKFem.ModCheck.Arguments)

save(BKFem.ModCheck.MCMC
     ,file = "outputs/Burkina_Faso_ModelCheck_RESULTS.Rdata")


###
### * Tuning
###
################################################################################

###
### Metropolis acceptance proportions
###

BKFem.ModCheck.MCMC$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

BKFem.raftLew.SMN <- list()

BKFem.raftLew.SMN$var.q025.r0125.s95 <-
    raftery.diag(BKFem.ModCheck.MCMC$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$var.q975.r0125.s95 <-
    raftery.diag(BKFem.ModCheck.MCMC$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.ModCheck.MCMC$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.ModCheck.MCMC$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.ModCheck.MCMC$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.ModCheck.MCMC$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

BKFem.raftLew.SMN$lambda.n.q025.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$lambda.n.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
BKFem.raftLew.SMN$lambda.n.q975.r0125.s95 <-
    raftery.diag(log(BKFem.ModCheck.MCMC$lambda.n.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(BKFem.raftLew.SMN, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(BKFem.raftLew.SMN$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(BKFem.raftLew.SMN$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
BKFem.raftLew.SMN$max.vars
BKFem.raftLew.SMN$max.vitals

## SAVE
save(BKFem.raftLew.SMN
     , file = "outputs/Burkina_Faso_ModelCheck_RAFTERYLEWIS.Rdata")


###
### * PLOTS AND TABLES
###
################################################################################

###
### Load normal likelihood chains
###

load(file = "outputs/Burkina_Faso_Recon_RESULTS.Rdata")


###
### ** Table 5
###
################################################################################

###
### Fertility Rate
###

old.post.q <-
    5 * apply(BKFem.Recon.MCMC$fert.rate.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))
new.post.q <-
    5 * apply(BKFem.ModCheck.MCMC$fert.rate.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))


### Absolute Relative Differences

fert.rate.q.abs.diff <- abs(new.post.q - old.post.q)
fert.rate.q.abs.rel.diff <- fert.rate.q.abs.diff / old.post.q
fert.rate.q.mean.abs.rel.diff <-
    apply(fert.rate.q.abs.rel.diff, 1, "mean")
fert.rate.q.max.abs.rel.diff <-
    apply(fert.rate.q.abs.rel.diff, 1, "max")


###
### Survival Proportion
###

old.post.q <-
    apply(BKFem.Recon.MCMC$surv.prop.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))
new.post.q <-
    apply(BKFem.ModCheck.MCMC$surv.prop.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))


### Absolute Relative Differences

surv.prop.q.abs.diff <- abs(new.post.q - old.post.q)
surv.prop.q.abs.rel.diff <- surv.prop.q.abs.diff / old.post.q
surv.prop.q.mean.abs.rel.diff <-
    apply(surv.prop.q.abs.rel.diff, 1, "mean")
surv.prop.q.max.abs.rel.diff <-
    apply(surv.prop.q.abs.rel.diff, 1, "max")


###
### Migration Proportion
###

old.post.q <-
    apply(BKFem.Recon.MCMC$mig.prop.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))

new.post.q <-
    apply(BKFem.ModCheck.MCMC$mig.prop.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))


### Absolute Differences

mig.prop.q.abs.diff <- abs(new.post.q - old.post.q)
mig.prop.q.mean.abs.diff <-
    apply(mig.prop.q.abs.diff, 1, "mean")
mig.prop.q.max.abs.diff <-
    apply(mig.prop.q.abs.diff, 1, "max")


### Absolute Relative Differences

mig.prop.q.abs.rel.diff <- mig.prop.q.abs.diff / abs(old.post.q)
mig.prop.q.mean.abs.rel.diff <-
    apply(mig.prop.q.abs.rel.diff, 1, "mean")
mig.prop.q.max.abs.rel.diff <-
    apply(mig.prop.q.abs.rel.diff, 1, "max")


###
### Baseline Count
###

old.post.q <-
    apply(BKFem.Recon.MCMC$baseline.count.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))
new.post.q <-
    apply(BKFem.ModCheck.MCMC$baseline.count.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))


### Absolute Relative Differences

baseline.count.q.abs.diff <- abs(new.post.q - old.post.q)
baseline.count.q.abs.rel.diff <- baseline.count.q.abs.diff / abs(old.post.q)
baseline.count.q.mean.abs.rel.diff <-
    apply(baseline.count.q.abs.rel.diff, 1, "mean")
baseline.count.q.max.abs.rel.diff <-
    apply(baseline.count.q.abs.rel.diff, 1, "max")


###
### Population counts
###

old.post.q <-
    apply(BKFem.Recon.MCMC$lx.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))
new.post.q <-
    apply(BKFem.ModCheck.MCMC$lx.mcmc, 2, FUN = "quantile"
          ,probs = c(0.025, 0.5, 0.975))


### Absolute Relative Differences

lx.q.abs.diff <- abs(new.post.q - old.post.q)
lx.q.abs.rel.diff <- lx.q.abs.diff / abs(old.post.q)
lx.q.mean.abs.rel.diff <-
    apply(lx.q.abs.rel.diff, 1, "mean")
lx.q.max.abs.rel.diff <-
    apply(lx.q.abs.rel.diff, 1, "max")


###
### Combine
###

q.mean.abs.rel.diff <-
    rbind(fert.rate.q.mean.abs.rel.diff
      ,surv.prop.q.mean.abs.rel.diff
      ,mig.prop.q.mean.abs.rel.diff
      ,baseline.count.q.mean.abs.rel.diff
      ,lx.q.mean.abs.rel.diff
      )

q.max.abs.rel.diff <-
    rbind(fert.rate.q.max.abs.rel.diff
      ,surv.prop.q.max.abs.rel.diff
      ,mig.prop.q.max.abs.rel.diff
      ,baseline.count.q.max.abs.rel.diff
      ,lx.q.max.abs.rel.diff
      )

BKFem.quantile.diff <-
    list(q.mean.abs.rel.diff = q.mean.abs.rel.diff
         ,q.max.abs.rel.diff = q.max.abs.rel.diff
         ,mig.prop.q.mean.abs.diff = mig.prop.q.mean.abs.diff
         ,mig.prop.q.max.abs.diff = mig.prop.q.max.abs.diff
         )


###
### TABLE 5
###

t2.cf.table1 <-
    as.data.frame(interleave(BKFem.quantile.diff[[1]]
      ,BKFem.quantile.diff[[2]]))

## use non-relative diffs for migration
t2.cf.table1[c("mig.prop.q.mean.abs.rel.diff", "mig.prop.q.max.abs.rel.diff"),] <-
    c(BKFem.quantile.diff[[3]], BKFem.quantile.diff[[4]])

t2.cf.table1$stat <- rep(c("mean", "max"), 5)

t2.cf.table1$name <-
    rep(c("fertility rate", "survival proportion"
            ,"migration proportion", "pop.\\ count, 1960"
            ,"pop.\\ count, 1965--2005"), each = 2)

t2.cf.table2 <- recast(t2.cf.table1, formula = name + variable ~ stat)

t2.cf.table3 <-
    within(t2.cf.table2, {
              variable <- as.numeric(sapply(strsplit(levels(variable)[variable], "%")
               ,"[[", 1)) / 100
              name[c(2, 3, 5, 6, 8, 9, 11, 12, 14, 15)] <- ""
              max <- 100 * max
              mean <- 100* mean
          })

colnames(t2.cf.table3) <-
    c("parameter", "quantile", "max ARD (\\%)", "mean ARD (\\%)")

BKFem.modCheck.table5 <- t2.cf.table3[c(1:3, 13:15, 4:12), c(1, 2, 4, 3)]

## SAVE
save(BKFem.modCheck.table5,
     file = "outputs/tables/Table_5.Rdata"
     )

write.csv(BKFem.modCheck.table5,
     file = "outputs/tables/Table_5.csv"
     )


###
### ** Lambda.n boxplots
###
################################################################################

## Save space
rm(BKFem.Recon.MCMC); gc(); warning("'BKFem.Recon.MCMC' removed from workspace")

plot.df <- melt(as.data.frame(BKFem.ModCheck.MCMC$lambda.n.mcmc)
                ,measure.vars = colnames(BKFem.ModCheck.MCMC$lambda.n.mcmc))
cnames <-
    strsplit(levels(plot.df[,"variable"])[plot.df[,"variable"]], split = "[^0-9]")
plot.df$ageF <- sapply(cnames, "[[", 2)
plot.df$yearF <- sapply(cnames, "[[", 1)

png(width = 7, height = 7, units = "in", res = 360, pointsize = 10
    ,file = file.path("outputs", "additional_plots", "BKFem_lambda_n_boxplot20110720_SMN.png")
    )
bwplot(log(value) ~ ageF | yearF, data = plot.df
       ,panel = function(x, y, ...) {
           panel.refline(h = log(1))
           panel.bwplot(x, y, ...)
       }
       ,as.table = TRUE
       ,ylab = expression(log(lambda[n]))
       ,xlab = "age"
       ,scales = list(x = list(rot = 45))
       )
dev.off()


###
### Largest median
###

lambda.quantiles <- apply(BKFem.ModCheck.MCMC$lambda.n.mcmc, 2, "quantile"
                          ,c(0.025, 0.25, 0.5, 0.75, 0.975))
Burkina.ModCheck.largestLambda <-
    lambda.quantiles[,which.max(lambda.quantiles["50%",]),drop=FALSE]

save(Burkina.ModCheck.largestLambda
     ,file = "outputs/Burkina_Faso_ModelCheck_LARGESTLAMBDA.Rdata"
     )
