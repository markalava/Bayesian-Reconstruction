################################################################################
###
###  DESC:          Plot survival proportions
###
###  DATE CREATED:  2017-09-09
###
###  AUTHOR:        Mark C Wheldon
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
### SYNOPSIS
###
### Calculate hyperparameters for survival proportions in two ways:
### i) Approximate, conservative way
### ii) More accurate way using numerical optimization
### then compare by plotting the 0.05 and 0.95 quantiles of the prior
### distribution for survival proportions, back-transformed to the proportion
### scale.
###
###-----------------------------------------------------------------------------
###
### REFERENCE
###
### Online appendix to:
###
### Wheldon, M. C., Raftery, A. E., Clark, S. J., & Gerland, P. (2015). Bayesian
### Reconstruction of Two-Sex Populations by Age: Estimating Sex Ratios at Birth
### and Sex Ratios of Mortality. Journal of the Royal Statistical Society:
### Series A (Statistics in Society), 178(4),
### 977–1007. https://doi.org/10.1111/rssa.12104
###
################################################################################

### * SET UP
################################################################################

set.seed(1)

data.path <- "data"
prog.path <- "R_Code"

example(source, echo = FALSE)
sourceDir(prog.path)

logit <- function(p) log(p / (1 - p))
invlogit <- function(x) exp(x) / (1 + exp(x))

### * Inputs
################################################################################

### ** Initial Estimates
################################################################################

load(file.path(data.path, "thai_initial_ests.RData"))

### * Compare Methods of Calculating Hyperparams
################################################################################

### ** Approximation ('absDev.s <- log(1.4)')
################################################################################

### Absolute deviation with 90 percent probability

absDev.f <- log(0.9)
absDev.s <- log(1.4)
absDev.g <- 0.2
absDev.n <- log(0.9)
absDev.srb <- log(0.9)

### Make the betas

invGam.params <-
    list(al.f = 0.5
         ,be.f = absDev.beta(absDev.f, prob = 0.9, alpha = 0.5)
         ,al.s = 0.5
         ,be.s = absDev.beta(absDev.s, prob = 0.9, alpha = 0.5)
         ,al.g = 0.5
         ,be.g = absDev.beta(absDev.g, prob = 0.9, alpha = 0.5)
         ,al.n = 0.5
         ,be.n = absDev.beta(absDev.n, prob = 0.9, alpha = 0.5)
         ,al.srb = 0.5
         ,be.srb = absDev.beta(absDev.srb, prob = 0.9, alpha = 0.5)
         )

### Plot U and L for each initial estimate

logit.death.init <-
    logit(1 - c(as.data.frame.table(asSurvTHAI.mat$female)$Freq
          ,as.data.frame.table(asSurvTHAI.mat$male)$Freq))
logit.death.5 <-
    logit.death.init -
    qt(p = 0.95, df = 2 * invGam.params$al.s) * sqrt(invGam.params$be.s/invGam.params$al.s)
logit.death.95 <-
    logit.death.init +
    qt(p = 0.95, df = 2 * invGam.params$al.s) * sqrt(invGam.params$be.s/invGam.params$al.s)
orig.lower.death.diff <-
    (invlogit(logit.death.init) - invlogit(logit.death.5)) / invlogit(logit.death.init)
orig.upper.death.diff <-
    (invlogit(logit.death.95) - invlogit(logit.death.init)) / invlogit(logit.death.init)
plot(rep(invlogit(logit.death.init), 2)
     ,c(orig.lower.death.diff, orig.upper.death.diff)
     ,col = c(rep("black", length(orig.lower.death.diff))
      ,rep("red", length(orig.upper.death.diff)))
     ,ylim = c(0, 0.6)
     ,xlab = "death prop"
     ,ylab = ""
     ,sub = paste("median = "
      ,prettyNum(median(c(orig.lower.death.diff, orig.upper.death.diff)))
      ,", range = ["
      ,prettyNum(range(c(orig.lower.death.diff, orig.upper.death.diff)))[1]
      ,", "
      ,prettyNum(range(c(orig.lower.death.diff, orig.upper.death.diff)))[2]
      ,"]"
      ,sep = ""
      )
     )
abline(h = 0.1, lty = 2, col = "blue")
legend("topright", pch = 1, col = c("black", "red")
       ,legend = c(expression((s^{"*"} - s^{0.05})/s^{"*"})
       ,expression((s^{0.95} - s^{"*"})/s^{"*"})
        )
       )

### ** Optimal (use 'make.hyper.params()')
################################################################################

invGam.params <-
    make.hyper.params(absDev = list(fert = 0.1, surv = 0.1, mig = 0.2
             ,pop = 0.1, srb = 0.1)
             ,prob = list(fert = 0.9, surv = 0.9, mig = 0.9, pop = 0.9, srb = 0.9)
             ,alpha = list(fert = 0.5, surv = 0.5, mig = 0.5, pop = 0.5, srb = 0.5)
             ,s.star = unlist(asSurvTHAI.mat)
             )

invGam.params

### Plot U and L for each initial estimate

logit.death.init <-
    logit(1 - c(as.data.frame.table(asSurvTHAI.mat$female)$Freq
          ,as.data.frame.table(asSurvTHAI.mat$male)$Freq))

logit.death.5 <-
    logit.death.init -
    qt(p = 0.95, df = 2 * invGam.params$al.s) * sqrt(invGam.params$be.s/invGam.params$al.s)

logit.death.95 <-
    logit.death.init +
    qt(p = 0.95, df = 2 * invGam.params$al.s) * sqrt(invGam.params$be.s/invGam.params$al.s)

orig.lower.death.diff <-
    (invlogit(logit.death.init) - invlogit(logit.death.5)) / invlogit(logit.death.init)
orig.upper.death.diff <-
    (invlogit(logit.death.95) - invlogit(logit.death.init)) / invlogit(logit.death.init)

plot(rep(invlogit(logit.death.init), 2)
     ,c(orig.lower.death.diff, orig.upper.death.diff)
     ,col = c(rep("black", length(orig.lower.death.diff))
      ,rep("red", length(orig.upper.death.diff)))
     ,ylim = c(0, 0.6)
     ,xlab = "death prop"
     ,ylab = ""
     ,sub = paste("median = "
      ,prettyNum(median(c(orig.lower.death.diff, orig.upper.death.diff)))
      ,", range = ["
      ,prettyNum(range(c(orig.lower.death.diff, orig.upper.death.diff)))[1]
      ,", "
      ,prettyNum(range(c(orig.lower.death.diff, orig.upper.death.diff)))[2]
      ,"]"
      ,sep = ""
      )
     )
abline(h = 0.1, lty = 2, col = "blue")
legend("topright", pch = 1, col = c("black", "red")
       ,legend = c(expression((s^{"*"} - s^{0.05})/s^{"*"})
       ,expression((s^{0.95} - s^{"*"})/s^{"*"})
        )
       )
