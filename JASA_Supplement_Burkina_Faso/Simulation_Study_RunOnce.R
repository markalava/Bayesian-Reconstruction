################################################################################
###
### TITLE:              Simulation_Study_RunOnce.R
###
### DATE:               23 May 2012
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Run one replicate of the simulation study in
###                     "Reconstructing Past Populations with Uncertainty from
###                     Fragmentary Data" submitted to Journal of the American
###                     Statistical Association" to produce output in
###                     Figure 2.
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
### This script runs the simulation study in Section 4 of the paper once to
### produce the plots in Figure 2. The full simulation study is in the file
### 'Simulation_Study_Run_Full.R'.
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

## Combine all in a list for input into simulation study function
tvr <- list(asFertTRUE.mat = asFertTRUE.mat
            ,asSurvTRUE.mat = asSurvTRUE.mat
            ,asMigTRUE.mat = asMigTRUE.mat
            ,baselineTRUE.mat = baselineTRUE.mat
            ,censusTRUE.mat = censusTRUE.mat
            )


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


###
### ** Control parameters for simulation and reconstruction
###
################################################################################

###
### Number of iterations to use for the reconstruction
###

## Don't set this to less than 5000 otherwise raftery.diag might fail.
## Chain is re-run if Metropolis acceptance proportions are outside [0.1,0.5] or
## Raftery-Lewis diagnostic suggests a longer chain.
n.iter <- 9000
burn.in <- 100


###
### * RUN THE SIMULATION
###
################################################################################

###
### ** Run Once (generates output for Figure 2)
###
################################################################################

## Seed
set.seed(1)

simulationStudy.output.run.once <-
    simStudy.estimation.once(#.. the array determining number of times to run
             rep.ar.count =1

             #.. sim stats
             ,runs.per.node=1, desired.overall.reps=1, overall.reps=1
             ,cluster.size=1

             #.. coverage level
             ,alpha = 0.05

             #.. algorithm parameters

                  ,start.iter = 5E3
                  ,start.burn.in = 600
                  ,prop.vars = sim.study.prop.vars

             ,max.tune.reruns = 5
             ,max.iter = 5E4, max.burn.in = 2E3
             ,min.iter = 5E3, min.burn.in = 600
             ,runRaftLew = TRUE
             ,checkAR = TRUE # check acceptance proportions
             ,ar.lower = 0.1, ar.upper = 0.5

             #.. est model arguments
             ,ccmp.f = "popRecon.ccmp.female"
                             , age.size=5, census.columns =2:5, fert.rows=c(2, 3)
                  ,s.tol = 10^(-10)
                  ,verb = TRUE

             #.. hyper parameters
             ,hyper.params = hyper.params

             #.. source functions
             ,estMod.f.source = "pop_reconstruction_functions.R"

             #.. function names
             ,estMod.f = "popRecon.sampler"

             #.. true values
             ,true.values = tvr

             #.. profile?
             ,Rprof.file.path = "perform_parallel_logs/RprofFiles/"

             #.. sink?
             ,sink.file.path = "perform_parallel_logs/SinkFiles/"

             #.. dump?
             ,dump.file.path =  "perform_parallel_logs/DumpFiles/"

             #.. save output?
             ,save.file = "outputs/Simulation_Study_run_once_FULL_RESULTS.Rdata"
             )[[1]]          # untangle list

save(simulationStudy.output.run.once
     ,file = "outputs/Simulation_Study_run_once_COVERAGE_RESULTS.Rdata"
     )


###
### ** Figure 2
###
################################################################################

load(file = "outputs/Simulation_Study_run_once_FULL_RESULTS.Rdata")

## Plotting parameters
t.col <- trellis.par.get("superpose.line")$col[1]
m.col <- trellis.par.get("superpose.line")$col[2]
qMed.col <- trellis.par.get("superpose.line")$col[3]
qLim.col <- qMed.col
t.pch <- 1
m.pch <- 3
qMed.pch <- 4
qLim.pch <- NA
t.lty <- 5
m.lty <- 2
qMed.lty <- 1
qLim.lty <- 4
m.cex <- 0.8
qMed.cex <- 0.8
qLim.cex <- 0.8
t.lwd <- 2
m.lwd <- 2
qMed.lwd <- 2
qLim.lwd <- 1
graphics.off()


###
### Figure 2 (a). Age Specific Fertility Rate
###

## Quantiles of posterior sample
emp.quantiles <-
     apply(simulationStudy.mcmc.samples$fert.rate.mcmc, 2
           ,function(z) quantile(z, probs = c(0.025, 0.5, 0.975))
           )
dimnames(emp.quantiles) <-
    list(as.character(c(0.025, 0.5, 0.975))
         ,colnames(simulationStudy.mcmc.samples$fert.rate.mcmc)
         )
emp.quantiles.melt <- melt(emp.quantiles)
X2.split <-
    strsplit(as.character(levels(emp.quantiles.melt$X2)[emp.quantiles.melt$X2])
             ,"\\."
             )
emp.quantiles.melt$years <- sapply(X2.split, "[[", 1)
emp.quantiles.melt$ages <- as.numeric(sapply(X2.split, "[[", 2))
emp.quantiles.df <-
    rename.vars(emp.quantiles.melt[,c("X1", "value", "years", "ages")]
                ,from = c("X1", "value"), to = c("quant", "fert.rate")
                )

## Initial estimates generated in this run
init.est.df <-
    rename.vars(melt(simulationStudy.mcmc.samples$fixed.params$mean.fert.rate[2:3,])
                ,from = c("X1", "X2", "value")
                ,to = c("ages", "years", "fert.rate")
                )
init.est.df$quant <- 77

## True values
fert.rate.true.df <-
    rename.vars(melt(asFertTRUE.mat[2:3,])
                ,from = c("X1", "X2", "value")
                ,to = c("ages", "years", "fert.rate")
                )
fert.rate.true.df$quant <- 99

## Bind all data
fert.rate.final.df <-
    rbind(fert.rate.true.df, init.est.df, emp.quantiles.df)

## Plot
pdf(file = "outputs/plots/simulation_study_Figure_2_a.pdf", width = 7, height = 7)
xyplot(fert.rate ~ ages | ordered(years)
             ,data = fert.rate.final.df
             ,groups = quant
             ##,subset = quant < 1
             ,panel = function(x, y, ...) {
                 panel.refline(h = c(1.5, 2))
                 panel.xyplot(x, y, ...)
             }
             ,type = "b"
             ,xlab = "age"
             ,ylab = "age-specific fertility rate"
             ,ylim = c(0, 0.8)
             ,col = c(qLim.col, qMed.col, qLim.col, m.col, t.col)
             ,pch = c(qLim.pch, qMed.pch, qLim.pch, m.pch, t.pch)
             ,lty = c(qLim.lty, qMed.lty, qLim.lty, m.lty, t.lty)
             ,lwd = c(qLim.lwd, qMed.lwd, qLim.lwd, m.lwd, t.lwd)
             ,key = list(text = list(c("Post. median", "95% Post. Int."
                         ,"initial est.", "truth"))
                ,lines = list(lty = c(qMed.lty, qLim.lty, m.lty, t.lty))
                ,col = c(qMed.col, qLim.col, m.col, t.col)
              ,pch = c(qMed.pch, qLim.pch, m.pch, t.pch)
               ,type = c("b", "b", "b", "b")
              ,columns = 3)
             ,as.table = TRUE
             ,las = 1
             ,par.settings = list(superpose.symbol = list(cex = 1))
             )
dev.off()


###
### Figure 2 (b). Total Fertility Rate
###

## Sum sample of age-specific fertility rates to get TFR
tfr.mcmc.mat <-
  matrix(0, nrow = nrow(simulationStudy.mcmc.samples$fert.rate.mcmc)
         ,ncol = 4
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(simulationStudy.mcmc.samples$fert.rate.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )
fert.rate.mcmc.colYrs <-
  sapply(strsplit(colnames(simulationStudy.mcmc.samples$fert.rate.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
for(i in 1:ncol(tfr.mcmc.mat)) {
  colYrs.index <- fert.rate.mcmc.colYrs == colnames(tfr.mcmc.mat)[i]
  tfr.mcmc.mat[,i] <-
    apply(simulationStudy.mcmc.samples$fert.rate.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) 5 * sum(z)
          )
}
tfr.quantiles.mat <-
    apply(tfr.mcmc.mat, 2, FUN = function(z)
                        {
                          quantile(z, probs = c(0.025, 0.5, 0.975))
                        })

## Initial estimates generated in this run
tfr.init.est.mat <-
    apply(simulationStudy.mcmc.samples$fixed.params$mean.fert.rate[2:3,], 2
          ,FUN = function(z) 5 * sum(z)
           )

## True TFR
tfr.true.mat <-
    apply(asFertTRUE.mat[2:3,], 2
          ,FUN = function(z) 5 * sum(z)
          )

## Plot
pdf(file = "outputs/plots/simulation_study_Figure_2_b.pdf", width = 7, height = 7)
plot(seq(from = 1960, to = 1975, by = 5), tfr.quantiles.mat[3,]
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
     ,lwd = qLim.lwd
     ,ylim = 5 * c(0.4, 1.1)
     ,ylab = "Total Fertility Rate"
     ,xlab = "year"
     ,las = 1)
lines(seq(from = 1960, to = 1975, by = 5), tfr.quantiles.mat[1,]
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
      ,lwd = qLim.lwd)
lines(seq(from = 1960, to = 1975, by = 5), tfr.quantiles.mat[2,]
     ,type = "b", lty = qMed.lty, col = qMed.col, pch = qMed.pch
      ,lwd = qMed.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,tfr.init.est.mat
      ,type = "b"
      ,col = m.col, pch = m.pch, lty = m.lty, lwd = m.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,tfr.true.mat
      ,type = "b"
      ,col = t.col, pch = t.pch, lty = t.lty, lwd = t.lwd)
legend("topright", lty = c(qMed.lty, qLim.lty, m.lty, t.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd, t.lwd)
       ,col = c(qMed.col, qMed.col, t.col, m.col)
       ,pch = c(qMed.pch, qLim.pch, t.pch, m.pch)
       ,legend = c("median", "95% PI", "truth", "initial est." )
       ,cex = 0.85
       )
dev.off()


###
### Figure 2 (c). Life expectancy at birth.
###

leb.mcmc.mat <-
  matrix(0, nrow = nrow(simulationStudy.mcmc.samples$surv.prop.mcmc)
         ,ncol = 4
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(simulationStudy.mcmc.samples$surv.prop.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )
surv.prop.mcmc.colYrs <-
    sapply(strsplit(colnames(simulationStudy.mcmc.samples$surv.prop.mcmc)
                    ,"\\."), FUN = function(z) z[[1]])
for(i in 1:ncol(leb.mcmc.mat)) {
  colYrs.index <- surv.prop.mcmc.colYrs == colnames(leb.mcmc.mat)[i]
  leb.mcmc.mat[,i] <-
    apply(simulationStudy.mcmc.samples$surv.prop.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) {
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }
        )
}
leb.quantiles.mat <- apply(leb.mcmc.mat, 2, FUN = function(z)
                        {
                          quantile(z, probs = c(0.025, 0.5, 0.975))
                        })

## Initial estimates generated in this run
leb.init.est.mat <-
  apply(simulationStudy.mcmc.samples$fixed.params$mean.surv.prop, 2
        ,FUN = function(z) {
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }
        )

## True leb
leb.true.mat <-
  apply(asSurvTRUE.mat, 2
        ,FUN = function(z) {
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }
        )


## Plot
pdf(file = "outputs/plots/simulation_study_Figure_2_c.pdf", width = 7, height = 7)
plot(seq(1960, 1975, by = 5), leb.quantiles.mat[3,]
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd
     ,ylim = c(14, 17.5)
     ,ylab = "Life expectancy at birth (years)"
     ,xlab = "year"
     ,las = 1)
lines(seq(1960, 1975, by = 5), leb.quantiles.mat[1,]
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd)
lines(seq(1960, 1975, by = 5), leb.quantiles.mat[2,]
     ,type = "b", lty = qMed.lty, col = qMed.col, lwd = qMed.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,leb.init.est.mat
      ,type = "b"
      ,col = m.col, lty = m.lty, pch = m.pch, lwd = m.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,leb.true.mat
      ,type = "b"
      ,col = t.col, lty = t.lty, pch = t.pch, lwd = t.lwd)
legend("topright", lty = c(qMed.lty, qLim.lty, m.lty, t.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd, t.lwd)
       ,col = c(qMed.col, qMed.col, t.col, m.col)
       ,pch = c(qMed.pch, qLim.pch, t.pch, m.pch)
       ,legend = c("median", "95% PI", "truth", "initial est.")
       ,cex = 0.85
       )
dev.off()


###
### Figure 2 (d). Net number of migrants.
###

## NB: Can't simply sum migration proportions because they are based on
##     different population totals. The two functions below are needed to
##     convert migration proportions into counts collapsed over age.

## Functions to create the Leslie matrix (Section 3.2.1 of the paper)
make.leslie.matrix <-
    function(pop, surv, fert, srb = 1.05, age.int = 5, label.dims = FALSE)
{
    ##-- Make the leslie matrix for CCMPP --##
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
    ##   proj.steps
    ##   age.int :  needed for correct interpretation of survival
    ##                and fertility rates
    ##   label.dims
    ##           :  should output have dimnames set?
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

    n.age.grps <- length(pop)
    n.surv <- length(surv)

    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)

    k <- 1/(1+srb) * surv[1] * 0.5
    dbl.fert <- age.int*fert + c(age.int*fert[-1], 0) * surv[-1]
    lesM[1,] <- k * dbl.fert

    lesM[2:n.age.grps,1:(n.age.grps-1)] <- diag(surv[-c(1,n.surv)])
    lesM[n.age.grps,n.age.grps] <- surv[n.surv]

    if(label.dims) {
        age.labs <- seq(from = 0, by = 5, length = n.age.grps)
        dimnames(lesM) <- list(age.labs, age.labs)
    }
    return(lesM)
}

## Function to calculate total migration count
total.mig.count <- function(n1, n2, L)
{
    ##-- Find net number of migrants in a CCMPP projection --##
    ##
    ## ARGUMENTS
    ##
    ##   n1      :  Population count vector at time t
    ##   n2      :  Population count vector at time t + delta
    ##   L       :  Leslie matrix used to get population at t + delta
    ##
    ##
    ## METHOD
    ##
    ## Invert n2 = L(n1 + 0.5 mig) + (0.5)*mig
    ## Can get proportions by pre-multiplying output by 'solve(diag(n1))'
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

    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    L <- as.matrix(L)

    return(2 * solve(L + diag(nrow(L))) %*% (n2 - L %*% n1))
}

## Function 'popRecon.ccmp.female' needed; sourced in from
## 'pop_reconstruction_functions.R' at top of this file.


## MCMC SAMPLE FOR TOTAL NET NUMBER OF MIGRANTS

## Prepare output matrix
net.mig.mcmc.mat <-
  matrix(0, nrow = nrow(simulationStudy.mcmc.samples$mig.prop.mcmc)
         ,ncol = 4
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(simulationStudy.mcmc.samples$mig.prop.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )
mig.prop.mcmc.colYrs <-
  sapply(strsplit(colnames(simulationStudy.mcmc.samples$mig.prop.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
mig.prop.mcmc.colYrsUniq <- unique(mig.prop.mcmc.colYrs)

## Combine population counts at baseline and in subsequent years
pop.mat <- cbind(simulationStudy.mcmc.samples$baseline.count.mcmc
                ,simulationStudy.mcmc.samples$lx.mcmc)

## Index for years
pop.mat.colYrs <- sapply(strsplit(colnames(pop.mat)
                         ,"\\."), FUN = function(z) z[[1]])
pop.mat.colYrsUniq <- unique(pop.mat.colYrs)

## cycle through mcmc sample
for(k in 1:nrow(simulationStudy.mcmc.samples$mig.prop.mcmc)) {

    ## cycle through years
    for(i in 1:ncol(net.mig.mcmc.mat)) {

        ## 5-year sub-intervals for indexing columns
        mig.colYrs.index <-
            colnames(net.mig.mcmc.mat) == mig.prop.mcmc.colYrsUniq[i]
        surv.colYrs.index <-
            surv.prop.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        fert.colYrs.index <-
            fert.rate.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        pop.colYrs.index1 <-
            pop.mat.colYrs == substr(mig.prop.mcmc.colYrsUniq[i]
            ,start = 2, stop = 5)
        pop.colYrs.index2 <-
            pop.mat.colYrs == as.numeric(substr(mig.prop.mcmc.colYrsUniq[i]
            ,start = 2, stop = 5)) + 5

        ## get vital rates and make leslie matrix
        sk <- simulationStudy.mcmc.samples$surv.prop.mcmc[k,surv.colYrs.index]
        fk <- rep(0, 4)
        fk[2:3] <- simulationStudy.mcmc.samples$fert.rate.mcmc[k,fert.colYrs.index]
        popk1 <- pop.mat[k,pop.colYrs.index1]
        popk2 <- pop.mat[k,pop.colYrs.index2]
        Lk <- make.leslie.matrix(pop = popk1, surv = sk, fert = fk, srb = 1.05
                           ,age.int = 5)

        ## calculate net number of migrants
        netMigk <- total.mig.count(n1 = popk1, n2 = popk2, L = Lk)

        ## store
        net.mig.mcmc.mat[k, mig.colYrs.index] <- sum(netMigk)
    }
}

## Posterior quantiles
net.mig.quantiles.mat <- apply(net.mig.mcmc.mat, 2, FUN = function(z)
                        {
                          quantile(z, probs = c(0.025, 0.5, 0.975))
                        })


## INITIAL ESTIMATE OF NET NUMBER OF MIGRANTS

## Prepare output matrix
net.mig.init.est.mat <- rep(0, 4)
names(net.mig.init.est.mat) <-
    colnames(simulationStudy.mcmc.samples$fixed.params$mean.mig.prop)

## Combine population counts at baseline and in subsequent years
pop.input.mat <-
    popRecon.ccmp.female(pop=simulationStudy.mcmc.samples$fixed.params$mean.baseline.count
                      ,surv=simulationStudy.mcmc.samples$fixed.params$mean.surv.prop
                      ,fert=simulationStudy.mcmc.samples$fixed.params$mean.fert.rate
                      ,mig=simulationStudy.mcmc.samples$fixed.params$mean.mig.prop
                      )

## Calculate input net migration
for(k in 1:(ncol(pop.input.mat)-1)) {
    Lk <- make.leslie.matrix(pop = pop.input.mat[,k]
                       ,surv = simulationStudy.mcmc.samples$fixed.params$mean.surv.prop[,k]
                       ,fert = simulationStudy.mcmc.samples$fixed.params$mean.fert.rate[,k]
                       ,srb = 1.05
                       ,age.int = 5)
        netMigk <- total.mig.count(n1 = pop.input.mat[,k]
                                ,n2 = pop.input.mat[,k+1]
                                ,L = Lk)
    net.mig.init.est.mat[k] <- sum(netMigk)
}


## TRUE NET NUMBER OF MIGRANTS

## Prepare output matrix
net.mig.true.mat <- rep(0, 4)
names(net.mig.true.mat) <- colnames(asMigTRUE.mat)

## True population counts
pop.true.mat <-
    popRecon.ccmp.female(pop=baselineTRUE.mat, surv=asSurvTRUE.mat
                      ,fert=asFertTRUE.mat, mig=asMigTRUE.mat)

## Calculate true net migration
for(k in 1:(ncol(pop.true.mat)-1)) {
    Lk <- make.leslie.matrix(pop = pop.true.mat[,k], surv = asSurvTRUE.mat[,k]
                       ,fert = asFertTRUE.mat[,k], srb = 1.05
                       ,age.int = 5)
        netMigk <- total.mig.count(n1 = pop.true.mat[,k]
                                ,n2 = pop.true.mat[,k+1]
                                ,L = Lk)
    net.mig.true.mat[k] <- sum(netMigk)
}


## PLOT

pdf(file = "outputs/plots/simulation_study_Figure_2_d.pdf", width = 7, height = 7)
plot(seq(from = 1960, to = 1975, by = 5), net.mig.quantiles.mat[3,]/1E3
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd
     ,ylim = c(-10, 10)
     ,ylab = "Total Net Number of Migrations (000s)"
     ,xlab = "year"
     ,las = 1)
lines(seq(from = 1960, to = 1975, by = 5), net.mig.quantiles.mat[1,]/1E3
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
     ,lwd = qLim.lwd)
lines(seq(from = 1960, to = 1975, by = 5), net.mig.quantiles.mat[2,]/1E3
     ,type = "b", lty = qMed.lty, col = qMed.col, pch = qMed.pch
      ,lwd = qMed.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,net.mig.init.est.mat/1E3
      ,type = "b"
      ,col = m.col, lty = m.lty, pch = m.pch, lwd = m.lwd)
lines(seq(from = 1960, to = 1975, by = 5)
      ,net.mig.true.mat/1E3
      ,type = "b"
      ,col = t.col, lty = t.lty, pch = t.pch, lwd = t.lwd)
abline(h = 0, col = "grey")
legend("topleft", lty = c(qMed.lty, qLim.lty, m.lty, t.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd, t.lwd)
       ,col = c(qMed.col, qMed.col, m.col, t.col)
       ,pch = c(qMed.pch, qMed.pch, m.pch, t.pch)
       ,legend = c("median", "95% PI", "initial est.", "true")
       ,cex = 0.85
       )
dev.off()
