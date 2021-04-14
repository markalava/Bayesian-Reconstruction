################################################################################
###
###  TITLE:             2_tuning.R
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
###  Tuning diagnostics: conditional variances to choose Metropolis proposal
###  variances.
###
###  MCMC convergence diagnostic.
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
################################################################################

### * SET UP
################################################################################

set.seed(1)

library(coda)

data.path <- "data"
prog.path <- "R_Code"
results.path <- "results"

example(source, echo = FALSE)
sourceDir(prog.path)


### * Inputs
################################################################################

## MCMC Chain
load(file = file.path(results.path, "thai_mcmc.RData"), verbose = TRUE)


### * Tuning
################################################################################

### ** Check acceptance ratios, conditional vars
################################################################################

## Algorithm stats
ThaiMcmc$alg.stats$acceptance.proportions


### Plot acceptance ratios
###

pdf(width = 9, height = 9
    ,file = file.path(results.path, "Thai_20130612srb10_accProps.pdf")
    )
par(mfrow = c(3,2))
for(i in 1:3) {
    z <- ThaiMcmc$alg.stats$acceptance.proportions[[i]]
    zn <- names(ThaiMcmc$alg.stats$acceptance.proportions)[i]

    matplot(x=colnames(z), y = t(z), type = "b", pch = 1, main = zn
     ,ylim = c(0.1, 0.6), xlab = "year", ylab = "acc. prop.")

    rect(xleft = as.numeric(min(colnames(z))) - 5
         ,ybottom = 0.2
         ,xright = as.numeric(max(colnames(z))) + 5
         ,ytop = 0.5
         ,col = "grey"
         ,density = 10
         )
}
z <- ThaiMcmc$alg.stats$acceptance.proportions[["baseline.count"]]
rownames(z) <- sapply(strsplit(rownames(z), "[^0-9]"), "[[", 1)
plot(x =rownames(z)
     ,y = z
     ,type = "b", main = "baseline.count"
     ,ylim = c(0.1, 0.6), xlab = "age", ylab = "acc. prop."
     )
rect(xleft = as.numeric(min(rownames(z))) - 5
     ,ybottom = 0.2
     ,xright = as.numeric(max(rownames(z))) + 5
     ,ytop = 0.5
     ,col = "grey"
     ,density = 10
     )
z <- ThaiMcmc$alg.stats$acceptance.proportions[["srb"]]
plot(x =colnames(z)
     ,y = z
     ,type = "b", main = "srb"
     ,ylim = c(0.1, 0.6), xlab = "age", ylab = "acc. prop."
     )
rect(xleft = as.numeric(min(colnames(z))) - 5
     ,ybottom = 0.2
     ,xright = as.numeric(max(colnames(z))) + 5
     ,ytop = 0.5
     ,col = "grey"
     ,density = 10
     )
with(ThaiMcmc$alg.stats$acceptance.proportions, {
     barplot(height = c(fert = sigmasq.f, surv = sigmasq.s, mig = sigmasq.g
             , pop = sigmasq.n, srb = sigmasq.srb)
             ,ylim = c(0,1)
             ,main = "variances"
             )
     abline(h = c(0.1, 0.5, 0.9), lty = 3)
         }
     )
dev.off()


## Fixed params
ThaiMcmc$fixed.params

## init.vals
ThaiMcmc$init.vals

## alg.params
ThaiMcmc$alg.params


### ** Calculate conditional variances
################################################################################

vitalCondVars <- list()

vitalCondVars$fert.rate <-
  chain.cond.vars.feb08(log(ThaiMcmc$fert.rate.mcmc))

vitalCondVars$surv.prop <-
    (chain.cond.vars.feb08(logit(ThaiMcmc$surv.prop.mcmc[["female"]])) +
        chain.cond.vars.feb08(logit(ThaiMcmc$surv.prop.mcmc[["male"]]))) / 2

vitalCondVars$mig <-
  (chain.cond.vars.feb08(ThaiMcmc$mig.prop.mcmc[["female"]]) +
    chain.cond.vars.feb08(ThaiMcmc$mig.prop.mcmc[["male"]])) / 2

vitalCondVars$population.count <-
    (chain.cond.vars.feb08(log(ThaiMcmc$baseline.count.mcmc[["female"]])) +
     chain.cond.vars.feb08(log(ThaiMcmc$baseline.count.mcmc[["male"]]))) / 2

vitalCondVars$srb <-
    chain.cond.vars.feb08(log(ThaiMcmc$srb.mcmc))

vitalCondVars$variances <-
    chain.cond.vars.feb08(log(ThaiMcmc$variances.mcmc))

## PRINT
vitalCondVars

## Save
save(vitalCondVars, file = file.path(results.path, "vitalCondVars.Rdata"))

## Plot
new.prop.vars <- lapply(vitalCondVars, "*", 2.3^2)
pdf("results/propVars_compare_0612srb10.pdf", pointsize = 9)
par(mfrow = c(3,2))
plot(c(0,0), xlim = c(0, max(ThaiMcmc$alg.params$prop.varcovar$fert.rate
     ,new.prop.vars$fert.rate)) * 1.01
     ,ylim = c(0, max(ThaiMcmc$alg.params$prop.varcovar$fert.rate
     ,new.prop.vars$fert.rate)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "fert.rate")
points(ThaiMcmc$alg.params$prop.varcovar$fert.rate, new.prop.vars$fert.rate)
abline(a = 0, b = 1, col = "blue")
plot(c(0,0), xlim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$surv.prop, "[", 1, 1)
     ,new.prop.vars$surv.prop)) * 1.01
     ,ylim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$surv.prop, "[", 1, 1)
     ,new.prop.vars$surv.prop)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "surv.prop")
points(sapply(ThaiMcmc$alg.params$prop.varcovar$surv.prop, "[", 1, 1), new.prop.vars$surv.prop)
abline(a = 0, b = 1, col = "blue")
plot(c(0,0), xlim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$mig, "[", 1, 1)
     ,new.prop.vars$mig)) * 1.01
     ,ylim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$mig, "[", 1, 1)
     ,new.prop.vars$mig)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "mig")
points(sapply(ThaiMcmc$alg.params$prop.varcovar$mig, "[", 1, 1), new.prop.vars$mig)
abline(a = 0, b = 1, col = "blue")
plot(c(0,0), xlim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
     ,new.prop.vars$baseline.pop.count)) * 1.01
     ,ylim = c(0, max(sapply(ThaiMcmc$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
     ,new.prop.vars$baseline.pop.count)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "baseline.pop.count")
points(sapply(ThaiMcmc$alg.params$prop.varcovar$baseline.pop.count, "[", 1, 1)
       , new.prop.vars$baseline.pop.count)
abline(a = 0, b = 1, col = "blue")
plot(c(0,0), xlim = c(0, max(unlist(ThaiMcmc$alg.params$prop.varcovar$srb)
     ,new.prop.vars$srb)) * 1.01
     ,ylim = c(0, max(unlist(ThaiMcmc$alg.params$prop.varcovar$srb)
     ,new.prop.vars$srb)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "srb")
points(unlist(ThaiMcmc$alg.params$prop.varcovar$srb), new.prop.vars$srb)
abline(a = 0, b = 1, col = "blue")
plot(c(0,0), xlim = c(0, max(unlist(ThaiMcmc$alg.params$prop.varcovar$variances)
     ,new.prop.vars$variances)) * 1.01
     ,ylim = c(0, max(unlist(ThaiMcmc$alg.params$prop.varcovar$variances)
     ,new.prop.vars$variances)) * 1.01
     ,type = "n", xlab = "old", ylab = "new", main = "variances")
points(unlist(ThaiMcmc$alg.params$prop.varcovar$variances), new.prop.vars$variances)
abline(a = 0, b = 1, col = "blue")
dev.off()


### ** Raftery-Lewis Diagnostics
################################################################################

Thai.raftLew <- list()


## Fertility, variances, SRB
Thai.raftLew$var.q025.r0125.s95 <-
    raftery.diag(ThaiMcmc$variances.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$var.q975.r0125.s95 <-
    raftery.diag(ThaiMcmc$variances.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$fert.q025.r0125.s95 <-
    raftery.diag(log(ThaiMcmc$fert.rate.mcmc)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$fert.q975.r0125.s95 <-
    raftery.diag(log(ThaiMcmc$fert.rate.mcmc)
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$srb.q025.r0125.s95 <-
    raftery.diag(ThaiMcmc$srb.mcmc
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$srb.q975.r0125.s95 <-
    raftery.diag(ThaiMcmc$srb.mcmc
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Females
Thai.raftLew$surv.q025.r0125.s95$female <-
    raftery.diag(logit(ThaiMcmc$surv.prop.mcmc$female)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$surv.q975.r0125.s95$female <-
    raftery.diag(logit(ThaiMcmc$surv.prop.mcmc$female)
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$mig.q025.r0125.s95$female <-
    raftery.diag(ThaiMcmc$mig.prop.mcmc$female
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$mig.q975.r0125.s95$female <-
    raftery.diag(ThaiMcmc$mig.prop.mcmc$female
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$baseline.q025.r0125.s95$female <-
    raftery.diag(log(ThaiMcmc$baseline.count.mcmc$female)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$baseline.q975.r0125.s95$female <-
    raftery.diag(log(ThaiMcmc$baseline.count.mcmc$female)
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$lx.q025.r0125.s95$female <-
    raftery.diag(log(ThaiMcmc$lx.mcmc$female)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$lx.q975.r0125.s95$female <-
    raftery.diag(log(ThaiMcmc$lx.mcmc$female)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Males
Thai.raftLew$surv.q025.r0125.s95$male <-
    raftery.diag(logit(ThaiMcmc$surv.prop.mcmc$male)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$surv.q975.r0125.s95$male <-
    raftery.diag(logit(ThaiMcmc$surv.prop.mcmc$male)
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$mig.q025.r0125.s95$male <-
    raftery.diag(ThaiMcmc$mig.prop.mcmc$male
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$mig.q975.r0125.s95$male <-
    raftery.diag(ThaiMcmc$mig.prop.mcmc$male
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$baseline.q025.r0125.s95$male <-
    raftery.diag(log(ThaiMcmc$baseline.count.mcmc$male)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$baseline.q975.r0125.s95$male <-
    raftery.diag(log(ThaiMcmc$baseline.count.mcmc$male)
                 ,q = 0.975, r = 0.0125, s = 0.95)
Thai.raftLew$lx.q025.r0125.s95$male <-
    raftery.diag(log(ThaiMcmc$lx.mcmc$male)
                 ,q = 0.025, r = 0.0125, s = 0.95)
Thai.raftLew$lx.q975.r0125.s95$male <-
    raftery.diag(log(ThaiMcmc$lx.mcmc$male)
                 ,q = 0.975, r = 0.0125, s = 0.95)

## Maximum N for each variable
raftLewMax <-
    try(
    as.data.frame(lapply(Thai.raftLew, FUN = function(y) {
        mvp <- function(z) {
            maxVar.pos <- which(z$resmatrix == max(z$resmatrix[,2])
                                ,arr.ind = TRUE)
            z$resmatrix[maxVar.pos[1],]
        }
        if(identical(names(y), c("female", "male"))) lapply(y, "mvp")
        else mvp(y)
    })), FALSE)

if(!inherits(raftLewMax, "try-error")) {

    ## Maximum N for variances
    raftLewMax.vars <-
        which(raftLewMax[,1:2] == max(raftLewMax["N",1:2]), arr.ind = TRUE)
    Thai.raftLew$max.vars <-
        raftLewMax[1:2,raftLewMax.vars[2], drop = FALSE]

    ## Maximum N for vital rates, srb
    raftLewMax.vitals <-
        which(raftLewMax[,3:ncol(raftLewMax)] == max(raftLewMax["N",3:ncol(raftLewMax)])
              ,arr.ind = TRUE)
    Thai.raftLew$max.vitals <-
        raftLewMax[,-(1:2)][1:2,raftLewMax.vitals[2], drop = FALSE]

    ## PRINT
    Thai.raftLew$max.vars
    Thai.raftLew$max.vitals


### USE RAFTERY-LEWIS TO DETERMINE HOW MUCH OF CHAIN IS USED FROM
### NOW ON
###

    ## Use 0.9 quantile of Ns across variables and max burn in
    Thai.raftLew.M <-
        unlist(sapply(Thai.raftLew[1:10], FUN = function(z) {
            as.numeric(z[["resmatrix"]][,"M"])
        }))

    ## Use 0.9 quantile of Ns across variables and max burn in
    Thai.raftLew.N <-
        unlist(sapply(Thai.raftLew[1:10], FUN = function(z) {
            as.numeric(z[["resmatrix"]][,"N"])
        }))

    (Thai.raftLew.N.09Q <-
     quantile(Thai.raftLew.N, probs = 0.9))

    ## SAVE
    save(Thai.raftLew.M, Thai.raftLew.N
         ,Thai.raftLew.N.09Q
         ,file = file.path(results.path
          ,"Thai_raftLew_09Q_20130612srb10.Rdata"
          )
         )


### Which variables need more than n.iter?
###

    lapply(Thai.raftLew, function(y) {
        foo <- function(z) z$resmatrix[z$resmatrix[,"N"] > ThaiMcmc$alg.params$iters,,drop=FALSE]
        if(identical(names(y), c("female", "male"))) lapply(y, "foo") else foo(y)
    })


### Which variables need more than 90% ?
###

    lapply(Thai.raftLew, function(y) {
        foo <- function(z) {
            z$resmatrix[z$resmatrix[,"N"] > as.numeric(Thai.raftLew.N.09Q),,drop=FALSE]
        }
        if(identical(names(y), c("female", "male"))) lapply(y, "foo") else foo(y)
    })

}


### SAVE
###
save(Thai.raftLew
     ,file = file.path(results.path, "Thai_raftLew.Rdata")
     )
