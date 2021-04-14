################################################################################
###
### TITLE:              Burkina_Faso_Reconstruction.R
###
### DATE:               23 May 2012
###
### AUTHOR:             Mark C. Wheldon
###
### DESC:               Reconstruction of the female population of Burkina Faso
###                     from the paper "Reconstructing Past Populations with
###                     Uncertainty from Fragmentary Data" submitted to "Journal
###                     of the American Statistical Association".
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
### This script runs the reconstruction described in Section 5 of the
### paper. It sources functions and loads data from directories called 'source'
### and 'data' and outputs to a directory called 'outputs'.
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
         ,be.n = 0.0109
         )


###
### * Run Reconstruction
###
################################################################################

## Chain size
n.iter <- 5000#4E4
burn.in <- 100#500


## Other parameters
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
            ,census.columns = c(4,6,8,10)
            ,fert.rows = 4:10
            ,s.tol = 10^(-10)
            ,verb = TRUE
         )


## Run it *** TAKES A LONG TIME! ***
set.seed(1)
BKFem.Recon.MCMC <-
    do.call(popRecon.sampler, args = BKFem.recon.arguments)


## SAVE
save(BKFem.Recon.MCMC, file = "outputs/Burkina_Faso_Recon_RESULTS.Rdata")


###
### * Performance
###
################################################################################

###
### Metropolis acceptance ratios
###

BKFem.Recon.MCMC$alg.stats$acceptance.proportions[1:4]


###
### Raftery-Lewis Diagnostics
###

Burkina.Recon.RafteryLewis <- list()

Burkina.Recon.RafteryLewis$var.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$var.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis$fert.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$fert.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis$surv.q025.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC$surv.prop.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$surv.q975.r0125.s95 <-
    raftery.diag(popReconAux.logit(BKFem.Recon.MCMC$surv.prop.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis$mig.q025.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC$mig.prop.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$mig.q975.r0125.s95 <-
    raftery.diag(BKFem.Recon.MCMC$mig.prop.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis$baseline.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$baseline.count.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$baseline.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$baseline.count.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

Burkina.Recon.RafteryLewis$lx.q025.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$lx.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Burkina.Recon.RafteryLewis$lx.q975.r0125.s95 <-
    raftery.diag(log(BKFem.Recon.MCMC$lx.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    sapply(Burkina.Recon.RafteryLewis, FUN = function(z)
       {
           maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                               ,arr.ind = TRUE)
           z$resmatrix[maxVar.pos[1],]
       })

## Maximum N for variances
raftLewMax.var <-
    which(raftLewMax[,1:2] == max(raftLewMax[,1:2])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis$max.vars <-
 raftLewMax[1:2,raftLewMax.var[2], drop = FALSE])

## Maximum N for vital rates
raftLewMax.vitals <-
    which(raftLewMax[,3:10] == max(raftLewMax[,3:10])
          ,arr.ind = TRUE)
(Burkina.Recon.RafteryLewis$max.vitals <-
 raftLewMax[1:2, (2+raftLewMax.vitals[2]), drop = FALSE])

## PRINT
Burkina.Recon.RafteryLewis$max.vars
Burkina.Recon.RafteryLewis$max.vitals


###
### * Make Plots
###
################################################################################

###
### Graphical parameters
###

m.col <- trellis.par.get("superpose.line")$col[2]
q.col <- trellis.par.get("superpose.line")$col[3]
qMed.col <- q.col
qLim.col <- q.col

m.pch <- 3
qMed.pch <- 4
qLim.pch <- NA

m.lty <- 2
qMed.lty <- 1
qLim.lty <- 4

m.cex <- 0.8
qMed.cex <- 0.8
qLim.cex <- 0.8

m.lwd <- 2
qMed.lwd <- 2
qLim.lwd <- 1

m.name <- "initial est."

dev.off()


###
### ** Age-Specific Parameters
###
################################################################################

###
### Fertility
###

## Posterior quantiles
m <- BKFem.Recon.MCMC$fert.rate.mcmc
q = c(0.025, 0.5, 0.975)
q.vital <- apply(m, 2, function(z) quantile(z, probs = q))
dimnames(q.vital) <- list(as.character(q), colnames(m))

## Age, year labels
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
fert.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
fert.ages.numeric <- as.numeric(gsub("[^0-9]", "", fert.ages))

## Reshape data frame
mqvit <- melt(q.vital)
mqvit.col <- cbind(mqvit
                   ,expand.grid(quant = q, ages = fert.ages.numeric
                               ,years = years)
                   )

## Initial estimates
meas <- BKFem.Recon.MCMC$init.vals$fert.rate[BKFem.Recon.MCMC$alg.params$non.zero.fert.rows,]
mmeas.col <-
    cbind(value = melt(meas)$value
          ,expand.grid(ages = fert.ages.numeric
                       ,years = years, quant = 5) # use quant=5 for init.est
        )

## Plot
plot.df <-
    rbind(mqvit.col[,c("value", "years", "ages", "quant")]
          ,mmeas.col[,c("value", "years", "ages", "quant")]
          )
plot.df$quant <-
    factor(as.numeric(plot.df$quant)
           ,levels = unique(as.numeric(plot.df$quant))
           )
pdf(width = 7, height = 7, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_5.pdf")
print(xyplot(value ~ ages | ordered(years)
             ,data = plot.df
             ,groups = factor(quant)[,drop = TRUE]
             ,type = "b", cex = 1.1
             ,xlab = "age"
             ,ylab = "age-specific fertility rate"
             ,col = c(qLim.col, q.col, qLim.col, m.col)
             ,pch = c(qLim.pch, qMed.pch, qLim.pch, m.pch)
             ,lty = c(qLim.lty, qMed.lty, qLim.lty, m.lty)
             ,lwd = c(qLim.lwd, qMed.lwd, qLim.lwd, m.lwd)
             ,key = list(text = list(c("Post. median", "95% Post. Int."
                         ,m.name))
                ,lines = list(lty = c(qMed.lty, qLim.lty, m.lty))
                ,col = c(q.col, qLim.col, m.col)
              ,pch = c(qMed.pch, qLim.pch, m.pch)
               ,type = c("b", "b", "b")
              ,columns = 3)
             ,as.table = TRUE
             ,las = 1
             )
      )
dev.off()

## Save Values
asFert.post.quantiles.df <- plot.df
save(asFert.post.quantiles.df, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_5.Rdata")
write.csv(asFert.post.quantiles.df, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_5.csv"
          ,row.names = TRUE)

###
### Survival
###

## Calculate quantiles
m <- BKFem.Recon.MCMC$surv.prop.mcmc
q = c(0.025, 0.5, 0.975)
q.vital <- apply(m, 2, function(z) quantile(z, probs = q))
dimnames(q.vital) <- list(as.character(q), colnames(m))

## Age, year labels
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
surv.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
surv.ages.numeric <- as.numeric(gsub("[^0-9]", "", surv.ages))
ages.numeric <- seq(from = 0, to = 20, by = 5)

## Reshape data frame
mqvit <- melt(q.vital)
mqvit.col <- cbind(mqvit
                   ,expand.grid(quant = q, ages = surv.ages.numeric
                               ,years = years)
                   )

## Initial estimates
meas <- BKFem.Recon.MCMC$init.vals$surv.prop
mmeas.col <-
    cbind(value = melt(meas)$value
          ,expand.grid(ages = surv.ages.numeric
                       ,years = years, quant = 5)
        )

## Plot
plot.df <-
    rbind(mqvit.col[,c("value", "years", "ages", "quant")]
          ,mmeas.col[,c("value", "years", "ages", "quant")]
          )
plot.df$quant <-
    factor(as.numeric(plot.df$quant)
           ,levels = unique(as.numeric(plot.df$quant))
           )
pdf(width = 7, height = 7
    ,file = "outputs/additional_plots/Burkina_Faso_Reconstruction_asSurv_post_quantiles.pdf")
print(xyplot(value ~ ages | ordered(years)
             ,data = plot.df
             ,groups = quant
             ,type = "b", cex = 1.1
             ,xlab = "age"
             ,ylab = "age-specific survival proportion"
             ,col = c(qLim.col, q.col, qLim.col, m.col)
             ,pch = c(qLim.pch, qMed.pch, qLim.pch)
             ,lty = c(qLim.lty, qMed.lty, qLim.lty, m.lty, m.lty, m.lty)
             ,lwd = c(qLim.lwd, qMed.lwd, qLim.lwd, m.lwd, m.lwd, m.lwd)
             ,key = list(text = list(c("Post. median", "95% Post. Int."
                         ,m.name))
                ,lines = list(lty = c(qMed.lty, qLim.lty, m.lty))
                ,col = c(q.col, qLim.col, m.col)
              ,pch = c(qMed.pch, qLim.pch, m.pch)
               ,type = c("b", "b", "b")
              ,columns = 3)
             ,as.table = TRUE
             ,las = 1
             )
      )
dev.off()

## Save Values
asSurv.post.quantiles.df <- plot.df
save(asSurv.post.quantiles.df
     ,file = "outputs/additional_plots/Burkina_Faso_Reconstruction_asSurv_post_quantiles.Rdata")
write.csv(asSurv.post.quantiles.df
          ,file = "outputs/additional_plots/Burkina_Faso_Reconstruction_asSurv_post_quantiles.csv"
               ,row.names = TRUE)

###
### Just 1 - 5s0
###

## Plot
pdf("outputs/plots/Burkina_Faso_Reconstruction_Figure_6_c.pdf", width = 7, height = 7)
plot5q0.df <- subset(plot.df, subset = ages == 0)
plot5q0.df$s.5.0 <- 1 - plot5q0.df$value
plot(colnames(asSurvBKFEM.mat)
     ,subset(plot5q0.df, subset = quant == 0.975, select = "s.5.0"
             ,drop = TRUE)
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
     ,lwd = qLim.lwd
     ,ylim = c(1 - 0.9, 1 - 0.7)
     ,ylab = expression(1-paste(""[5], s[0], sep = ""))
     ,xlab = "year"
     ,las = 1)
lines(colnames(asSurvBKFEM.mat)
     ,subset(plot5q0.df, subset = quant == 0.025, select = "s.5.0"
             ,drop = TRUE)
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
      ,lwd = qLim.lwd)
lines(colnames(asSurvBKFEM.mat)
     ,subset(plot5q0.df, subset = quant == 0.5, select = "s.5.0"
             ,drop = TRUE)
     ,type = "b", lty = qMed.lty, col = qMed.col, pch = qMed.pch
      ,lwd = qMed.lwd)

lines(colnames(asSurvBKFEM.mat)
     ,subset(plot5q0.df, subset = quant == 5, select = "s.5.0"
             ,drop = TRUE)
      ,type = "b"
      ,col = m.col, pch = m.pch, lty = m.lty, lwd = m.lwd)

legend("topright", lty = c(qMed.lty, qLim.lty, m.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd)
       ,col = c(qMed.col, qLim.col, m.col)
       ,pch = c(qMed.pch, qLim.pch, m.pch)
       ,legend = c("median", "95% PI", "initial est.")
       ,cex = 0.85
       )
dev.off()

## Save Values
asSurv.5s0.post.quantiles.df <-
    subset(plot5q0.df, subset = quant %in% c("0.025", "0.5", "0.975"))
save(asSurv.5s0.post.quantiles.df
     ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_c.Rdata")
write.csv(asSurv.5s0.post.quantiles.df
          ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_c.csv"
          ,row.names = TRUE)


###
### Baseline
###

## Calculate quantiles
m <- BKFem.Recon.MCMC$baseline.count.mcmc
q = c(0.025, 0.5, 0.975)
q.vital <- apply(m, 2, function(z) quantile(z, probs = q))
dimnames(q.vital) <- list(as.character(q), colnames(m))

## Make ages and years labels
colspl <- strsplit(colnames(m), ".", fixed = TRUE)
years <- unique(sapply(colspl, FUN = function(z) z[1]))
baseline.ages <-
  unique(sapply(colspl, FUN = function(z) z[2]))[-length(surv.ages)]
baseline.ages.numeric <- as.numeric(gsub("[^0-9]", "", baseline.ages))

ages.numeric <- seq(from = 0, to = 15, by = 5)

## Reshape data frames
mqvit <- melt(q.vital)
mqvit.col <- cbind(mqvit
                   ,expand.grid(quant = q, ages = baseline.ages.numeric
                               ,years = years)
                   )

## Initial Estimates
meas <- BKFem.Recon.MCMC$init.vals$baseline.count
mmeas.col <-
    cbind(value = melt(meas)$value
          ,expand.grid(ages = baseline.ages.numeric
                       ,years = years, quant = 5)
        )
plot.df <-
    rbind(mqvit.col[,c("value", "years", "ages", "quant")]
          ,mmeas.col[,c("value", "years", "ages", "quant")]
          )
plot.df$quant <-
    factor(as.numeric(plot.df$quant)
           ,levels = unique(as.numeric(plot.df$quant))
           )

## Plot
pdf(width = 5, height = 5, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_4.pdf")
print(xyplot(value ~ ages | ordered(years)
             ,data = plot.df
             ,groups = quant
             ,type = "b", cex = 1.1
             ,xlab = "age"
             ,ylab = "age-specific baseline count"
             ,col = c(qLim.col, q.col, qLim.col, m.col)
             ,pch = c(qLim.pch, qMed.pch, qLim.pch, m.pch)
             ,lty = c(qLim.lty, qMed.lty, qLim.lty, m.lty)
             ,lwd = c(qLim.lwd, qMed.lwd, qLim.lwd, m.lwd)
             ,key = list(text = list(c("Post. median", "95% Post. Int."
                         ,m.name))
                ,lines = list(lty = c(qMed.lty, qLim.lty, m.lty))
                ,col = c(q.col, qLim.col, m.col)
              ,pch = c(qMed.pch, qLim.pch, m.pch)
               ,type = c("b", "b", "b")
              ,columns = 2)
             ,as.table = TRUE
             ,las = 1
             )
      )
dev.off()

## Save Values
baseline.count.post.quantiles.df <- plot.df
save(baseline.count.post.quantiles.df
     ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_4.Rdata")
write.csv(baseline.count.post.quantiles.df
          ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_4.csv"
          ,row.names = TRUE)


###
### ** Age-Summarized Measures
###
################################################################################

###
### Life Expectancy
###

## Calculate LEB from MCMC output
BKFem.Recon.LEB <-
  matrix(0, nrow = nrow(BKFem.Recon.MCMC$surv.prop.mcmc)
         ,ncol = ncol(asSurvBKFEM.mat)
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(BKFem.Recon.MCMC$surv.prop.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )
surv.prop.mcmc.colYrs <-
    sapply(strsplit(colnames(BKFem.Recon.MCMC$surv.prop.mcmc)
                    ,"\\."), FUN = function(z) z[[1]])
for(i in 1:ncol(BKFem.Recon.LEB)) {
  colYrs.index <- surv.prop.mcmc.colYrs == colnames(BKFem.Recon.LEB)[i]
  BKFem.Recon.LEB[,i] <-
    apply(BKFem.Recon.MCMC$surv.prop.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) {
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }
        )
}


## LEB from initial estimates
BKFem.leb.input <-
  apply(BKFem.Recon.MCMC$init.vals$surv.prop, 2
        ,FUN = function(z) {
            x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
            5 * sum(cumprod(x))
        }
        )

## Plot
pdf("outputs/plots/Burkina_Faso_Reconstruction_Figure_6_d.pdf", width = 7, height = 7)
plot.df <- apply(BKFem.Recon.LEB, 2, FUN = function(z)
                        {
                          quantile(z, probs = c(0.025, 0.5, 0.975))
                        })

plot(colnames(asSurvBKFEM.mat), plot.df[3,]
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd
     ,ylim = c(30, 60)
     ,ylab = "Life expectancy at birth (years)"
     ,xlab = "year"
     ,las = 1)
lines(colnames(asSurvBKFEM.mat), plot.df[1,]
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd)
lines(colnames(asSurvBKFEM.mat), plot.df[2,]
     ,type = "b", lty = qMed.lty, col = q.col, pch = qMed.pch
      ,lwd = qMed.lwd)

lines(colnames(asSurvBKFEM.mat)
      ,BKFem.leb.input
      ,type = "b"
      ,col = m.col, lty = m.lty, pch = m.pch, lwd = m.lwd)
legend("topleft", lty = c(qMed.lty, qLim.lty, m.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd)
       ,col = c(q.col, q.col, m.col)
       ,pch = c(qMed.pch, qLim.pch, m.pch)
       ,legend = c("median", "95% PI", "initial est.")
       ,cex = 0.85
       )
dev.off()


## Save Values
leb.post.quantiles.df <- plot.df
save(leb.post.quantiles.df, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_d.Rdata")
write.csv(leb.post.quantiles.df, file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_d.csv"
          ,row.names = TRUE)


###
### Total Fertility Rate
###

## Calculate from MCMC chains
BKFem.tfr.1004.loesFert <-
  matrix(0, nrow = nrow(BKFem.Recon.MCMC$fert.rate.mcmc)
         ,ncol = ncol(asFertBKFEM.mat)
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(BKFem.Recon.MCMC$fert.rate.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )
fert.rate.mcmc.colYrs <-
  sapply(strsplit(colnames(BKFem.Recon.MCMC$fert.rate.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
for(i in 1:ncol(BKFem.tfr.1004.loesFert)) {
  colYrs.index <- fert.rate.mcmc.colYrs == colnames(BKFem.tfr.1004.loesFert)[i]
  BKFem.tfr.1004.loesFert[,i] <-
    apply(BKFem.Recon.MCMC$fert.rate.mcmc[,colYrs.index]
        ,1
        ,FUN = function(z) sum(z)
        )
}

## TFR from input values
BKFem.tfr.input <-
  apply(BKFem.Recon.MCMC$init.vals$fert.rate[BKFem.Recon.MCMC$alg.params$non.zero.fert.rows,]
        ,2
        ,FUN = function(z) 5 * sum(z)
        )

## Plot
pdf("outputs/plots/Burkina_Faso_Reconstruction_Figure_6_a.pdf", width = 7, height = 7)
plot.df <- apply(BKFem.tfr.1004.loesFert, 2, FUN = function(z)
                        {
                          5*quantile(z, probs = c(0.025, 0.5, 0.975))
                        })
plot(colnames(asFertBKFEM.mat), plot.df[3,]
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
     ,lwd = qLim.lwd
     ,ylim = c(5, 9)
     ,ylab = "Total Fertility Rate"
     ,xlab = "year"
     ,las = 1)
lines(colnames(asFertBKFEM.mat), plot.df[1,]
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
      ,lwd = qLim.lwd)
lines(colnames(asFertBKFEM.mat), plot.df[2,]
     ,type = "b", lty = qMed.lty, col = qMed.col, pch = qMed.pch
      ,lwd = qMed.lwd)

lines(colnames(asFertBKFEM.mat)
      ,BKFem.tfr.input
      ,type = "b"
      ,col = m.col, pch = m.pch, lty = m.lty, lwd = m.lwd)
legend("topright", lty = c(qMed.lty, qLim.lty, m.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd)
       ,col = c(q.col, q.col, m.col)
       ,pch = c(qMed.pch, qLim.pch, m.pch)
       ,legend = c("median", "95% PI", "initial est.")
       ,cex = 0.85
       )
dev.off()

## Save Values
tfr.post.quantiles.df <- plot.df
save(tfr.post.quantiles.df
     ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_a.Rdata")
write.csv(tfr.post.quantiles.df
          ,file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_a.csv"
          ,row.names = TRUE)

###
### Total Net Migration
###

## NB: Can't simply sum migration proportions because they are based
##     on different population totals. Need to get net number of migrants
##     and convert back into proportions. Use Leslie matrix formula from
##     article draft.

## Make Leslie matrix function
make.leslie.matrix <-
    function(pop, surv, fert, srb = 1.05, age.int = 5, label.dims = FALSE)

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

{
    ## Constants
    n.age.grps <- length(pop)
    n.surv <- length(surv)

    ## Make Leslie matrix
    lesM <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)

    ## first row = fert and birth survival
    k <- 1/(1+srb) * surv[1] * 0.5
    dbl.fert <- age.int*fert + c(age.int*fert[-1], 0) * surv[-1]
    lesM[1,] <- k * dbl.fert

    ## rows 2:(n.age.grps) = survival ratios
    lesM[2:n.age.grps,1:(n.age.grps-1)] <- diag(surv[-c(1,n.surv)])
    lesM[n.age.grps,n.age.grps] <- surv[n.surv]

    if(label.dims) {
        age.labs <- seq(from = 0, by = 5, length = n.age.grps)
        dimnames(lesM) <- list(age.labs, age.labs)
    }

    ## return
    return(lesM)
}

## Calculate net number of migrants Function
netMig.calc <- function(n1, n2, L)
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

    ## Make sure inputs are of correct form
    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    L <- as.matrix(L)

    return(2 * solve(L + diag(nrow(L))) %*% (n2 - L %*% n1))
}

## Prepare output matrix
BKFem.Recon.netMig <-
  matrix(0, nrow = nrow(BKFem.Recon.MCMC$mig.prop.mcmc)
         ,ncol = ncol(asMigBKFEM.mat)
         ,dimnames = list(NULL,
            unique(sapply(strsplit(colnames(BKFem.Recon.MCMC$mig.prop.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])
                   )
            )
         )

## The 5-year sub-intervals to be used as an index into the columns of
## BKFem.Recon.netMig
mig.prop.mcmc.colYrs <-
  sapply(strsplit(colnames(BKFem.Recon.MCMC$mig.prop.mcmc)
                         ,"\\."), FUN = function(z) z[[1]])
mig.prop.mcmc.colYrsUniq <- unique(mig.prop.mcmc.colYrs)

## Concatenate baseline and lx to get a single matrix with population
## counts
pop.mat <- cbind(BKFem.Recon.MCMC$baseline.count.mcmc
                ,BKFem.Recon.MCMC$lx.mcmc)

## Index for population years
pop.mat.colYrs <- sapply(strsplit(colnames(pop.mat)
                         ,"\\."), FUN = function(z) z[[1]])
pop.mat.colYrsUniq <- unique(pop.mat.colYrs)

## cycle through mcmc sample trajectories
for(k in 1:nrow(BKFem.Recon.MCMC$mig.prop.mcmc)) {

    ## cycle through years
    for(i in 1:ncol(BKFem.Recon.netMig)) {

        ## 5-year sub-intervals for indexing columns
        mig.colYrs.index <-
            colnames(BKFem.Recon.netMig) == mig.prop.mcmc.colYrsUniq[i]
        surv.colYrs.index <-
            surv.prop.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        fert.colYrs.index <-
            fert.rate.mcmc.colYrs == mig.prop.mcmc.colYrsUniq[i]
        pop.colYrs.index1 <-
            pop.mat.colYrs == mig.prop.mcmc.colYrsUniq[i]
        pop.colYrs.index2 <-
            pop.mat.colYrs == as.numeric(mig.prop.mcmc.colYrsUniq[i]) + 5

        ## get vital rates and make leslie matrix
        sk <- BKFem.Recon.MCMC$surv.prop.mcmc[k,surv.colYrs.index]
        fk <- rep(0, nrow(asFertBKFEM.mat))
        fk[BKFem.Recon.MCMC$alg.params$non.zero.fert.rows] <- BKFem.Recon.MCMC$fert.rate.mcmc[k,fert.colYrs.index]
        popk1 <- pop.mat[k,pop.colYrs.index1]
        popk2 <- pop.mat[k,pop.colYrs.index2]
        Lk <- make.leslie.matrix(pop = popk1, surv = sk, fert = fk, srb = 1.05
                           ,age.int = 5)

        ## calculate net number of migrants
        netMigk <- netMig.calc(n1 = popk1, n2 = popk2, L = Lk)

        ## store
        BKFem.Recon.netMig[k, mig.colYrs.index] <- sum(netMigk)
    }
}

## NMIG from input values

## Prepare output matrix
BKFem.nmig.input <- rep(0, ncol(asMigBKFEM.mat))
names(BKFem.nmig.input) <-
    colnames(BKFem.Recon.MCMC$init.vals$mig.prop)

## Input population counts
pop.input.mat <-
    popRecon.ccmp.female(pop=BKFem.Recon.MCMC$init.vals$baseline.count
                      ,surv=BKFem.Recon.MCMC$init.vals$surv.prop
                      ,fert=BKFem.Recon.MCMC$init.vals$fert.rate
                      ,mig=BKFem.Recon.MCMC$init.vals$mig.prop
                      )

## Calculate input net migration
for(k in 1:(ncol(pop.input.mat)-1)) {
    Lk <- make.leslie.matrix(pop = pop.input.mat[,k]
                       ,surv = BKFem.Recon.MCMC$init.vals$surv.prop[,k]
                       ,fert = BKFem.Recon.MCMC$init.vals$fert.rate[,k]
                       ,srb = 1.05
                       ,age.int = 5)
        netMigk <- netMig.calc(n1 = pop.input.mat[,k]
                                ,n2 = pop.input.mat[,k+1]
                                ,L = Lk)
    BKFem.nmig.input[k] <- sum(netMigk)
}

## Plot
pdf("outputs/plots/Burkina_Faso_Reconstruction_Figure_6_b.pdf", width = 7, height = 7)
plot.df <- apply(BKFem.Recon.netMig, 2, FUN = function(z)
                        {
                          quantile(z, probs = c(0.025, 0.5, 0.975))
                        })
plot(colnames(asMigBKFEM.mat), plot.df[3,]/5E3
     ,type = "l", lty = qLim.lty, col = qLim.col, lwd = qLim.lwd
     ,ylim = c(-100, 100)
     ,ylab = "Total Net Number of Migrations (000s)"
     ,xlab = "year"
     ,las = 1)
lines(colnames(asMigBKFEM.mat), plot.df[1,]/5E3
     ,type = "l", lty = qLim.lty, col = qLim.col, pch = qLim.pch
     ,lwd = qLim.lwd)
lines(colnames(asMigBKFEM.mat), plot.df[2,]/5E3
     ,type = "b", lty = qMed.lty, col = qMed.col, pch = qMed.pch
      ,lwd = qMed.lwd)

lines(colnames(asMigBKFEM.mat)
      ,BKFem.nmig.input/5E3
      ,type = "b"
      ,col = m.col, lty = m.lty, pch = m.pch, lwd = m.lwd)
abline(h = 0, col = "grey")
legend("topleft", lty = c(qMed.lty, qLim.lty, m.lty)
       ,lwd = c(qMed.lwd, qLim.lwd, m.lwd)
       ,col = c(q.col, q.col, m.col)
       ,pch = c(qMed.pch, qMed.pch, m.pch)
       ,legend = c("median", "95% PI", "initial est.")
       ,cex = 0.85
       )
dev.off()

## Output table of quantiles for report
total.net.mig.post.quantiles <- plot.df/5E3
save(total.net.mig.post.quantiles
     , file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_b.Rdata")
write.csv(total.net.mig.post.quantiles
          , file = "outputs/plots/Burkina_Faso_Reconstruction_Figure_6_b.csv"
          ,row.names = TRUE)

###
### Total Population at Baseline
###

m <- rowSums(BKFem.Recon.MCMC$baseline.count.mcmc)
q = c(0.025, 0.5, 0.975)
q.vital <- quantile(m, probs = q)

