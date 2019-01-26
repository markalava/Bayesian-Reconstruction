################################################################################
###
### TITLE:              Burkina_Faso_OOS_Pred_Check_results.R
###
### DATE:               23 May 2012
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Model checking for reconstruction of the female
###                     population of Burkina Faso without using the census data
###                     from the "Model Checking" section of the paper
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
### This script performs the out-of-sample predictive model check by re-running
### the Burkina Faso reconstruction four times, each time omitting one of the
### census data sets. Running this script as-is will perform the reconstructions
### in serial. Copying the relevant code into separate files and running them in
### parallel will be much faster.
###
###-----------------------------------------------------------------------------
###
### INPUT FILES:
###
### This script expects to find the following files:
### - pop_reconstruction_functions.R
### - burkina_faso_initial_estimates_population_counts.csv
###
###-----------------------------------------------------------------------------
###
### OUTPUTS:
###
### This script creates the directory outputs/, and subdirectories, into which
### outputs are saved. Figures and tables that appear in the paper are placed in
### the directories outputs/plots/ and outputs/tables/. These directories also
### contain .Rdata and .csv files containing the data values used in the plots.
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

memory.limit(4*1023)


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
                         )[,-(1:2)]
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
         ,be.n = 0.0109
         )


###
### Chain size
###

n.iter <- 5000#5E4
burn.in <- 50#500


###
### Other arguments
###

BKFem.recon.arguments <-
    list(#.. inverse gamma parameters
         al.f = invGam.params$al.f
         ,be.f = invGam.params$be.f
         ,al.s = invGam.params$al.s
         ,be.s = invGam.params$be.s
         ,al.g = invGam.params$al.g
         ,be.g = invGam.params$be.g
         ,al.n = invGam.params$al.n
         ,be.n = invGam.params$be.n

      #.. the rest
      ,n.iter = n.iter, burn.in = burn.in
            ,ccmp.f = "popRecon.ccmp.female"
            ,init.f = asFertBKFEM.mat
            ,init.s = asSurvBKFEM.mat
            ,init.g = asMigBKFEM.mat
            ,init.b = baselineBKFEM.mat
            ,init.sigmasq.f = 5
            ,init.sigmasq.s = 5
            ,init.sigmasq.g = 5
            ,init.sigmasq.n = 5
            ,mean.f = asFertBKFEM.mat
            ,mean.s = asSurvBKFEM.mat
            ,mean.g = asMigBKFEM.mat
            ,mean.b = baselineBKFEM.mat
            ,prop.vars = bkfas.recon.prop.vars
            ,pop.data = censusBKFEM.mat
            ,proj.periods = ncol(asFertBKFEM.mat)
            ,age.size = 5
            ,census.columns = c(6,8,10)
            ,fert.rows = 4:10
            ,s.tol = 10^(-10)
            ,verb = TRUE
         )


###
### * Run Reconstruction omitting 1975 census
###
################################################################################


## Run it *** TAKES A LONG TIME! ***
set.seed(1)
BKFem.Recon.MCMC.excl1975 <-
    do.call(popRecon.sampler, args = BKFem.recon.arguments)


## SAVE
save(BKFem.Recon.MCMC.excl1975
     ,file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1975_RESULTS.Rdata")


###
### ** Performance
###
################################################################################

###
### Metropolis acceptance ratios
###

BKFem.Recon.MCMC.excl1975$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

Burkina.Recon.RafteryLewis.excl1975 <- list()

Burkina.Recon.RafteryLewis.excl1975$var.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1975$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$var.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1975$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1975$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1975$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1975$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1975$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1975$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1975$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1975$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1975$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1975$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1975$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1975$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(Burkina.Recon.RafteryLewis.excl1975, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1975$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1975$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
Burkina.Recon.RafteryLewis.excl1975$max.vars
Burkina.Recon.RafteryLewis.excl1975$max.vitals

## SAVE
save(Burkina.Recon.RafteryLewis.excl1975, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1975_RAFTERYLEWIS.Rdata")


###
### ** Save space: remove MCMC objects for 1975 run
###
################################################################################

## This is probably a good idea if this script is run as-is.

rm(BKFem.Recon.MCMC.excl1975)


###
### * Run Reconstruction omitting 1985 census
###
################################################################################

## Run it *** TAKES A LONG TIME! ***
set.seed(1)
BKFem.Recon.MCMC.excl1985 <-
    do.call(popRecon.sampler, args = BKFem.recon.arguments)


## SAVE
save(BKFem.Recon.MCMC.excl1985, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1985_RESULTS.Rdata")


###
### ** Performance for 1985 run
###
################################################################################

###
### Metropolis acceptance ratios
###

BKFem.Recon.MCMC.excl1985$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

Burkina.Recon.RafteryLewis.excl1985 <- list()

Burkina.Recon.RafteryLewis.excl1985$var.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1985$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$var.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1985$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1985$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1985$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1985$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1985$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1985$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1985$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1985$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1985$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1985$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1985$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1985$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(Burkina.Recon.RafteryLewis.excl1985, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1985$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1985$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
Burkina.Recon.RafteryLewis.excl1985$max.vars
Burkina.Recon.RafteryLewis.excl1985$max.vitals

## SAVE
save(Burkina.Recon.RafteryLewis.excl1985, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1985_RAFTERYLEWIS.Rdata")


###
### ** Save space: remove MCMC objects for 1985 run
###
################################################################################

## This is probably a good idea if this script is run as-is.

rm(BKFem.Recon.MCMC.excl1985)


###
### * Run Reconstruction omitting 1995 census
###
################################################################################

## Run it *** TAKES A LONG TIME! ***
set.seed(1)
BKFem.Recon.MCMC.excl1995 <-
    do.call(popRecon.sampler, args = BKFem.recon.arguments)


## SAVE
save(BKFem.Recon.MCMC.excl1995, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1995_RESULTS.Rdata")


###
### ** Performance for 1985 run
###
################################################################################

###
### Metropolis acceptance ratios
###

BKFem.Recon.MCMC.excl1995$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

Burkina.Recon.RafteryLewis.excl1995 <- list()

Burkina.Recon.RafteryLewis.excl1995$var.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1995$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$var.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1995$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1995$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1995$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1995$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl1995$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1995$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1995$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl1995$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1995$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl1995$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl1995$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl1995$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(Burkina.Recon.RafteryLewis.excl1995, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1995$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl1995$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
Burkina.Recon.RafteryLewis.excl1995$max.vars
Burkina.Recon.RafteryLewis.excl1995$max.vitals

## SAVE
save(Burkina.Recon.RafteryLewis.excl1995, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl1995_RAFTERYLEWIS.Rdata")


###
### ** Save space: remove MCMC objects for 1995 run
###
################################################################################

## This is probably a good idea if this script is run as-is.

rm(BKFem.Recon.MCMC.excl1995)


###
### * Run Reconstruction omitting 2005 census
###
################################################################################

## Run it *** TAKES A LONG TIME! ***
set.seed(1)
BKFem.Recon.MCMC.excl2005 <-
    do.call(popRecon.sampler, args = BKFem.recon.arguments)


## SAVE
save(BKFem.Recon.MCMC.excl2005, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl2005_RESULTS.Rdata")


###
### ** Performance for 2005 run
###
################################################################################

###
### Metropolis acceptance ratios
###

BKFem.Recon.MCMC.excl2005$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

Burkina.Recon.RafteryLewis.excl2005 <- list()

Burkina.Recon.RafteryLewis.excl2005$var.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl2005$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$var.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl2005$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl2005$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl2005$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl2005$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC.excl2005$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl2005$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl2005$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC.excl2005$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl2005$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis.excl2005$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis.excl2005$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC.excl2005$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(Burkina.Recon.RafteryLewis.excl2005, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl2005$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis.excl2005$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
Burkina.Recon.RafteryLewis.excl2005$max.vars
Burkina.Recon.RafteryLewis.excl2005$max.vitals

## SAVE
save(Burkina.Recon.RafteryLewis.excl2005, file = "outputs/Burkina_Faso_OOS_Pred_Check_excl2005_RAFTERYLEWIS.Rdata")


###
### ** Save space: remove MCMC objects for 2005 run
###
################################################################################

## This is probably a good idea if this script is run as-is.

rm(BKFem.Recon.MCMC.excl2005)


###
### * Create results
###
################################################################################

###
### Load census counts
###

censusBKFEM.mat <-
    data.matrix(read.csv(file = "burkina_faso_initial_estimates_population_counts.csv"
                         ,header = TRUE, check.names = FALSE, row.names = 1
                         )[,-1]
                )

census.matrix <-
    matrix(censusBKFEM.mat, nrow = 1
           ,dimnames = list(NULL
            ,paste(expand.grid(rownames(censusBKFEM.mat), colnames(censusBKFEM.mat))[,2]
                   ,expand.grid(rownames(censusBKFEM.mat), colnames(censusBKFEM.mat))[,1]
                   ,sep = "."
                   )
           ))

###
### ** Re-Load MCMC Chains
###
################################################################################

###
### How many iterations should be kept?
###

## Use the Raftery-Lewis diagnostic to shorten the chains. All must be same
## length so choose maximum of recommended lengths.

max.M <- -1
max.N <- -1

for(yyyy in c("1975", "1985", "1995", "2005")) {
    load(paste("outputs/Burkina_Faso_OOS_Pred_Check_excl"
               ,yyyy
               ,"_RAFTERYLEWIS.Rdata"
               ,sep = ""
               ))
    max.M <- max(max.M
            ,get(paste("Burkina.Recon.RafteryLewis.excl"
                       ,yyyy, sep = ""))$max.vars["M",1]
            ,get(paste("Burkina.Recon.RafteryLewis.excl"
                       ,yyyy, sep = ""))$max.vitals["M",1]
            )
    max.N <- max(max.N
            ,get(paste("Burkina.Recon.RafteryLewis.excl"
                       ,yyyy, sep = ""))$max.vars["N",1]
            ,get(paste("Burkina.Recon.RafteryLewis.excl"
                       ,yyyy, sep = ""))$max.vitals["N",1]
            )
}


###
### 1975 Excluded
###

## Load and shorten MCMC chains
load("outputs/Burkina_Faso_OOS_Pred_Check_excl1975_RESULTS.Rdata")
rl.lx.mcmc.excl1975 <-
    window(BKFem.Recon.MCMC.excl1975$lx.mcmc
           ,start = max.M + 1
           ,end = max.M + max.N
           )
rm(BKFem.Recon.MCMC.excl1975)           # save space


###
### 1985 Excluded
###

## Load and shorten MCMC chains
load("outputs/Burkina_Faso_OOS_Pred_Check_excl1985_RESULTS.Rdata")
rl.lx.mcmc.excl1985 <-
    window(BKFem.Recon.MCMC.excl1985$lx.mcmc
           ,start = max.M + 1
           ,end = max.M + max.N
           )
rm(BKFem.Recon.MCMC.excl1985)           # save space


###
### 1995 Excluded
###

## Load and shorten MCMC chains
load("outputs/Burkina_Faso_OOS_Pred_Check_excl1995_RESULTS.Rdata")
rl.lx.mcmc.excl1995 <-
    window(BKFem.Recon.MCMC.excl1995$lx.mcmc
           ,start = max.M + 1
           ,end = max.M + max.N
           )
rm(BKFem.Recon.MCMC.excl1995)           # save space


###
### 2005 Excluded
###

## Load and shorten MCMC chains
load("outputs/Burkina_Faso_OOS_Pred_Check_excl2005_RESULTS.Rdata")
rl.lx.mcmc.excl2005 <-
    window(BKFem.Recon.MCMC.excl2005$lx.mcmc
           ,start = max.M + 1
           ,end = max.M + max.N
           )
rm(BKFem.Recon.MCMC.excl2005)           # save space


###
### Combine into a single mcmc object which contains the estimated counts for
### each census year from the run which excluded census data from that year
###

rl.lx.cen.deleted <- NULL

cat("\n\nMake rl.lx.cen.deleted\n----------------------\n\n")
for(yyyy in c("1975", "1985", "1995", "2005")) {
    cat(yyyy, "\n")
    rl.chain <-get(paste("rl.lx.mcmc.excl", yyyy, sep = ""))
    rl.lx.cen.deleted <-
        cbind(rl.lx.cen.deleted
          ,rl.chain[, (grep(yyyy, colnames(rl.chain)))]
          )
    rm(list = paste("rl.lx.mcmc.excl", yyyy, sep = "")) # save space
}

rl.lx.cen.deleted <- mcmc(rl.lx.cen.deleted)


###
### ** Mean absolute relative error
###
################################################################################

pred.medians <- apply(rl.lx.cen.deleted, 2, "median")
census.vec <- as.numeric(census.matrix)

ARE.medians <-
    matrix(abs((pred.medians - census.vec)/census.vec)
           ,nrow = nrow(censusBKFEM.mat), ncol = ncol(censusBKFEM.mat)
           ,dimnames = dimnames(censusBKFEM.mat)
           )

(BKFem.OOS.MARE <- mean(ARE.medians) * 100)

save(BKFem.OOS.MARE, file = "outputs/Burkina_Faso_OOS_MARE_RESULTS.Rdata")
