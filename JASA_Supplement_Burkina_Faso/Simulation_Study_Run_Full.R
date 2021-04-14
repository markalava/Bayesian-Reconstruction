################################################################################
###
### TITLE:              Simulation_Study_Run_Full.R
###
### DATE:               23 May 2012
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Run simulation study in "Reconstructing
###                     Past Populations with Uncertainty from Fragmentary Data"
###                     submitted to Journal of the American Statistical
###                     Association".
###
### REFERENCE:          Wheldon, M. C., Raftery, A. E., Clark, S. J.,
###                     & Gerland, P. (2013). Reconstructing Past
###                     Populations with Uncertainty from Fragmentary
###                     Data. Journal of the American Statistical
###                     Association, 108(501),
###                     96–110. http://doi.org/10.1080/01621459.2012.737729
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
### SYNOPSIS:
###
### This script runs the simulation study in Section 4 of the paper.
###
### It was designed to be run on a multi-processor UNIX-alike cluster using the
### R package 'snowFT()'. snowFT() will look for the file '.pvm_hosts' in your
### home directory. If you are not on UNIX-alike or do not wish to use snowFT(),
### the simulation can be run in serial; this will be /much/ slower. Change the
### variable 'use.snowFT' to control this behavior. Default is 'FALSE' but the
### results in the manuscript are from a run using snowFT(). Due to possible
### differences in random seeds, running the serial version may produce results
### slightly different from those in the manuscript but only within Monte Carlo
### error.
###
### If you use snowFT(), remember to set the home directory in the argument to
### the simulation function near line 367.
###
###-----------------------------------------------------------------------------
###
### INPUT FILES:
###
### This script expects to find the following files:
### - pop_reconstruction_functions.R
###
###-----------------------------------------------------------------------------
###
### OUTPUTS:
###
### Figures and tables that appear in the paper are placed in the directories
### outputs/plots/ and outputs/tables/. These directories also contain .Rdata and
### .csv files containing the data values used in the plots.
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

## Run in parallel using snowFT()? (Only available for *NIX systems)
## Assumes '.pvm_hosts' is in 'Sys.getenv("HOME")'
use.snowFT <- TRUE

## !!! YOU MUST ALSO SET THE VARIABLE 'working.directory' in the argument list
## !!! to performParallel() below. This should be the same as the working
## !!! directory specified above and it must be a string literal. It cannot be
## !!! an R symbol.

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

## Create directories for output
dir.create("outputs")
dir.create("outputs/plots")
dir.create("outputs/additional_plots")
dir.create("outputs/tables")
dir.create("perform_parallel_logs")
dir.create("perform_parallel_logs/DumpFiles")
dir.create("perform_parallel_logs/RprofFiles")
dir.create("perform_parallel_logs/SinkFiles")


###
### If on *NIX, start snowFT
###

if(use.snowFT && .Platform$OS.type == "unix") {
    library(snowFT)
    #library(rpvm)

  ##   hostfile <- paste(Sys.getenv("HOME"),"/.pvm_hosts",sep="")
  ## .PVM.start.pvmd(hostfile)

  ##   print(.PVM.config())

  ##   #.. how big is the cluster?
  ##   pvm.size <- nrow(.PVM.config())

    #.. run X times as many on each node
    ## cluster.size <- 2 * pvm.size
## } else {
    cluster.size <- 3
}
cluster.size


###
### Source files
###

## Population reconstruction functions
source("pop_reconstruction_functions.R")


###
### * CREATE NECESSARY OBJECTS
###
################################################################################

###
### ** 'True' vital rates (Table 1 in the paper)
###
################################################################################

asFertTRUE.mat <-
    matrix(c(0, 0.4, 0.3, 0, 0, 0.4, 0.3, 0, 0, 0.4, 0.3, 0, 0,
             0.4, 0.3, 0)
           ,nrow = 4
           ,dimnames = list(c("0", "5", "10", "15")
            ,c("[1960, 1965)", "[1965, 1970)", "[1970, 1975)", "[1975, 1980)"))
           )

asSurvTRUE.mat <-
    matrix(c(0.9, 0.95, 0.85, 0.8, 0.1, 0.9, 0.95, 0.85, 0.8,
             0.1, 0.9, 0.95, 0.85, 0.8, 0.1, 0.9, 0.95, 0.85, 0.8, 0.1)
           ,nrow = 5
           ,dimnames = list(c("0", "5", "10", "15", "20+"), c("[1960, 1965)",
            "[1965, 1970)", "[1970, 1975)", "[1975, 1980)"))
           )

asMigTRUE.mat <-
    matrix(c(-0.025, -0.05, -0.055, -0.005, -0.05, -0.1, -0.11,
             -0.01, 0.025, 0.05, 0.055, 0.005, 0.05, 0.1, 0.11, 0.01)
           ,nrow = 4
           ,dimnames = list(c("0", "5", "10", "15"), c("[1960, 1965)",
            "[1965, 1970)", "[1970, 1975)", "[1975, 1980)"))
           )

baselineTRUE.mat <-
    matrix(c(7500, 6000, 4000, 3000), nrow = 4
           ,dimnames = list(c("0", "5", "10", "15"), "1960")
           )

censusTRUE.mat <-
    matrix(c(8482, 6886, 4862, 3404, 9453, 7512, 5293, 3998, 11436
             ,9280, 6690, 4762, 14504, 11600, 8651, 6149)
           ,nrow = 4
           ,dimnames = list(c("0", "5", "10", "15")
            ,c("1965", "1970", "1975", "1980"))
           )

tvr <- list(asFertTRUE.mat = asFertTRUE.mat
            ,asSurvTRUE.mat = asSurvTRUE.mat
            ,asMigTRUE.mat = asMigTRUE.mat
            ,baselineTRUE.mat = baselineTRUE.mat
            ,censusTRUE.mat = censusTRUE.mat
            )

## Save as .Rdata and .csv
Table.1 <- list(asFertTRUE.mat = asFertTRUE.mat
            ,asSurvTRUE.mat = asSurvTRUE.mat
            ,asMigTRUE.mat = asMigTRUE.mat
            )
save(Table.1, file = "outputs/tables/Table_1.Rdata")
write.csv(asFertTRUE.mat, file = "outputs/tables/Table_1_fertility_rates.csv", row.names = FALSE)
write.csv(asSurvTRUE.mat, file = "outputs/tables/Table_1_survival_proportions.csv", row.names = FALSE)
write.csv(asMigTRUE.mat, file = "outputs/tables/Table_1_migration_proportions.csv", row.names = FALSE)

Table.2 <- cbind(baselineTRUE.mat, censusTRUE.mat)
save(Table.2, file = "outputs/tables/Table_2.Rdata")
write.csv(cbind(baselineTRUE.mat, censusTRUE.mat)
          ,file = "outputs/tables/Table_2.csv", row.names = FALSE)


###
### ** Level 4 hyperparameters (alpha and beta in Table 3 of the paper)
###
################################################################################

hyper.params <-
    list(al.f = 1
         ,be.f = 0.0109
         ,al.s = 1
         ,be.s = 0.0109
         ,al.g = 1
         ,be.g = 0.0436
         ,al.n = 1
         ,be.n = 0.0109
         )

## Save as .Rdata and .csv
Table.3 <-
    data.frame(v = c("f, s, n", "g")
               ,alpha = c(hyper.params$al.f, hyper.params$al.g)
               ,beta = c(hyper.params$be.f, hyper.params$be.g)
               ,quant.mae.0.025 = c(sqrt(popReconAux.qinvGamma(0.025
                ,hyper.params$al.f
                ,hyper.params$be.f
                )) * sqrt(2/pi)
                ,sqrt(popReconAux.qinvGamma(0.025
                ,hyper.params$al.g
                ,hyper.params$be.g
                )) * sqrt(2/pi)
                )
               ,quant.mae.0.25 = c(sqrt(popReconAux.qinvGamma(0.25
                ,hyper.params$al.f
                ,hyper.params$be.f
                )) * sqrt(2/pi)
                ,sqrt(popReconAux.qinvGamma(0.25
                ,hyper.params$al.g
                ,hyper.params$be.g
                )) * sqrt(2/pi)
                )
               ,quant.mae.0.5 = c(sqrt(popReconAux.qinvGamma(0.5
                ,hyper.params$al.f
                ,hyper.params$be.f
                )) * sqrt(2/pi)
                ,sqrt(popReconAux.qinvGamma(0.5
                ,hyper.params$al.g
                ,hyper.params$be.g
                )) * sqrt(2/pi)
                )
               ,quant.mae.0.5 = c(sqrt(popReconAux.qinvGamma(0.75
                ,hyper.params$al.f
                ,hyper.params$be.f
                )) * sqrt(2/pi)
                ,sqrt(popReconAux.qinvGamma(0.75
                ,hyper.params$al.g
                ,hyper.params$be.g
                )) * sqrt(2/pi)
                )
               ,quant.mae.0.5 = c(sqrt(popReconAux.qinvGamma(0.975
                ,hyper.params$al.f
                ,hyper.params$be.f
                )) * sqrt(2/pi)
                ,sqrt(popReconAux.qinvGamma(0.975
                ,hyper.params$al.g
                ,hyper.params$be.g
                )) * sqrt(2/pi)
                )
               )

save(Table.3, file = "outputs/tables/Table_3.Rdata")
write.csv(Table.3, file = "outputs/tables/Table_3.csv", row.names = FALSE)


###
### ** Control parameters for simulation and reconstruction
###
################################################################################

###
### Number of iterations to use for the reconstruction
###

## Don't set this to less than 5000 otherwise raftery.diag might fail.
n.iter <- 50#9000
burn.in <- 10


###
### Iterations in the 'outer' loop ('J' in the paper).
###

desired.overall.reps <- 3#200


###
### Runs per node (if using 'snowFT()').
###

## Depending on the cluster size, J might end up slightly larger than
## 'desired.overall.reps', but never smaller.

(runs.per.node <- ceiling(desired.overall.reps/cluster.size))
(overall.reps <- runs.per.node * cluster.size)


###
### * RUN THE SIMULATION
###
################################################################################

## !!! Remember to set the 'working.directory' argument !!!

## Seed
set.seed(1)

if(use.snowFT && .Platform$OS.type == "unix") {

    ##
    ## Use snowFT()
    ##

    printfun <- function(res, n, args = NULL)
    {
        cat(paste("\n", "replicate number ", n, "\n", sep = ""))
    }

    simStudy.pp.output <-
        performParallel(count = cluster.size, x = 1:cluster.size
                        ,fun = simStudy.estimation.once
                        ,printrepl = 1
                        ,printfun = "printfun"
                        ## ,initfun = function() {
                        ##     cat("\n\npvm.tid ", .PVM.mytid(), "\n"
                        ##         ,"pvm.parent ", .PVM.parent(), "\n"
                        ##         ,"pvm.siblings ", .PVM.siblings(), "\n"
                        ##         ,sep = ""
                        ##         )
                        ## }
                        ,exitfun = function() {
                            cat("\n\n WARNINGS\n")
                            warnings()
                            cat("\n\nTRACEBACK\n")
                            traceback()
                        }
                        ,seed = c(6,5,4,3,2,1)
                        ,ft_verbose = TRUE

                        ## ###########################################################
                        ## !! YOU MUST SET THIS !!
                        ## It should be the same as the working directory
                        ## specified above and it must be a string literal. It
                        ## cannot be an R symbol.
                        ,working.directory = "~/Documents/PPGp_Working/TEST_JASA_Scripts"
                        ## ###########################################################

                        ## snow arguments
                        ,runs.per.node = runs.per.node

                        ## alg control
                        ,start.iter = 5E3
                        ,start.burn.in = 100
                        ,prop.vars = sim.study.prop.vars

                        ## est model arguments
                        ,ccmp.f = "popRecon.ccmp.female"
                        ,age.size = 5
                        ,census.columns = 2:5
                        ,fert.rows = c(2, 3)
                        ,s.tol = 10^(-10)
                        ,verb = TRUE

                        ## hyper parameters
                        ,hyper.params = hyper.params

                        ## sources
                        ,estMod.f.source = "pop_reconstruction_functions.R"

                        ## true values
                        ,true.values = tvr

                        ## function names
                        ,estMod.f = "popRecon.sampler"

                        ## tuning control
                        ,max.tune.reruns = 5
                        ,max.iter = 5E4, max.burn.in = 2E3
                        ,min.iter = 5E3, min.burn.in = 600
                        ,runRaftLew = TRUE
                        ,checkAR = TRUE
                        ,ar.lower = 0.1, ar.upper = 0.5

                        ## sim stats
                        ,desired.overall.reps = desired.overall.reps
                        ,overall.reps = overall.reps
                        ,cluster.size = cluster.size
                        ,alpha = 0.05

                        ## sink file
                        ,sink.file.path = "perform_parallel_logs/SinkFiles/"

                        ## dump file
                        ,dump.file.path = "perform_parallel_logs/DumpFiles/"

                        ,Rprof.file.path = NULL

                        ## save full data file name
                        ,save.file = NULL

                        )

    ## Untangle all the lists
    simStudy.results <- list()
    k <- 1
    for(i in 1:(cluster.size)) {
        for(j in 1:runs.per.node) {
            simStudy.results[[k]] <- simStudy.pp.output[[i]][[j]]
            k <- k + 1
        }
    }

    length(simStudy.results)

    ## Keep only first 'desired.overall.reps' for reproducibility
    simStudy.results <- simStudy.results[1:desired.overall.reps]

} else {

    ##
    ## Run in serial
    ##

    ## Not optimal! Use snowFT() and performParallel() if possible

    simStudy.results <- list()
    for(j in 1:desired.overall.reps) {
        simStudy.results <-
            c(simStudy.results,
              simStudy.estimation.once(j
                                      ,runs.per.node = desired.overall.reps

                                      ## alg control
                                      ,start.iter = 5E3
                                      ,start.burn.in = 600
                                      ,prop.vars = sim.study.prop.vars

                                      ## est model arguments
                                      ,ccmp.f = "popRecon.ccmp.female"
                                      ,age.size = 5
                                      ,census.columns = 2:5
                                      ,fert.rows = c(2, 3)
                                      ,s.tol = 10^(-10)
                                      ,verb = TRUE

                                      ## hyper parameters
                                      ,hyper.params = hyper.params

                                      ## sources
                                      ,estMod.f.source = "pop_reconstruction_functions.R"

                                      ## true values
                                      ,true.values = tvr

                                      ## function names
                                      ,estMod.f = "popRecon.sampler"

                                      ## tuning control
                                      ,max.tune.reruns = 5
                                      ,max.iter = 5E4, max.burn.in = 2E3
                                      ,min.iter = 5E3, min.burn.in = 600
                                      ,runRaftLew = TRUE
                                      ,checkAR = TRUE
                                      ,ar.lower = 0.1, ar.upper = 0.5

                                      ## sim stats
                                      ,desired.overall.reps = desired.overall.reps
                                      ,overall.reps = overall.reps
                                      ,cluster.size = cluster.size
                                      ,alpha = 0.05

                                      ## sink file
                                      ,sink.file.path = NULL

                                      ## dump file
                                      ,dump.file.path = NULL

                                      ,Rprof.file.path = NULL

                                      ## save full data file name
                                      ,save.file = NULL
                                      )
              )

    }
}


###
### Save results
###

save(simStudy.results, file = "outputs/Simulation_Study_COVERAGE_RESULTS.Rdata")


###
### * GENERATE OUTPUTS
###
################################################################################

###
### Tuning
###

sapply(simStudy.results
       ,function(z) z$tuning.indicators
       )


###
### ** Table 4. Estimated Coverage Probabilities of 95 Percent (etc.)
###
################################################################################

fert.rate.coverage <-
    rowMeans(sapply(simStudy.results
                    ,function(z) z$fert.rate.postQuants[[1]][,"covered"]
                    )
             )
surv.prop.coverage <-
    rowMeans(sapply(simStudy.results
                    ,function(z) z$surv.prop.postQuants[[1]][,"covered"]
                    )
             )
mig.prop.coverage <-
    rowMeans(sapply(simStudy.results
                    ,function(z) z$mig.prop.postQuants[[1]][,"covered"]
                    )
             )
baseline.count.coverage <-
    rowMeans(sapply(simStudy.results
                    ,function(z) z$baseline.count.postQuants[[1]][,"covered"]
                    )
             )

t4.years.ages <-
    expand.grid(Ages = c("[0, 5)", "[5, 10)", "[10, 15)", "[15, 20)", "20+")
            ,Years = c("1960", "[1960, 1965)", "[1965, 1970)"
            ,"[1970, 1975)", "[1975, 1980)"
            )
            )

Table4 <-
    data.frame(t4.years.ages[,2:1]
               ,Population = NA, Fertility = NA, Survival = NA, Migration = NA
               )
Table4[1:length(baseline.count.coverage), "Population"] <- baseline.count.coverage
Table4[c(8,9,12,13,17,18,22,23), "Fertility"] <- fert.rate.coverage
Table4[6:nrow(Table4), "Survival"] <- surv.prop.coverage
Table4[c(6:9, 11:14, 16:19, 21:24), "Migration"] <- mig.prop.coverage

save(Table4, file = "outputs/tables/Table_4.Rdata")
write.csv(Table4, file = "outputs/tables/Table_4.csv")


################################################################################
### HALT PVM
###
### !!! note: THIS WILL TERMINATE THE R SESSION AS WELL

## if(use.snowFT) .PVM.halt()

###
################################################################################
