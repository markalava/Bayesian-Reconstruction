################################################################################
###
### TITLE:              ccmp_femDom_20120612.R
###
### AUTHOR:             Mark Wheldon
###
### DESC:               CCMPP, female dominant
###
###-----------------------------------------------------------------------------
###
### SYNOPSIS:
###
### Implements cohort component method of population projection as described in
### Preston et. al (2001) Ch. 6.
###
###-----------------------------------------------------------------------------
###
### CHANGES:
###
### (from ccmp_femDom_20120420.R)
### - do.checks default value now FALSE
###
################################################################################

ccmp.femDom.jun12 <-
    function(pop, surv, fert, srb = matrix(1.05, ncol = proj.steps)
             ,mig
             ,mig.type = c("prop.prev.pop", "net.count")
             ,proj.steps = ncol(fert), age.int = 5
             ,label.dims = FALSE, base.year, first.age = 0
             ,do.checks = FALSE, verb = FALSE, return.list = FALSE
             )

    ## Female dominant cohort-component method of population projection.
    ## Preston et al. (2001), Ch. 6

    ## ARGUMENTS
    ##
    ## - pop:
    ## Population count at baseline. List with up to two components: 'female' and
    ## 'male' each a column vector of age specific counts.
    ##
    ## - surv:
    ## List with up to two components: 'female' and 'male' each a matrix of
    ## survivorship proportions: the probability of reaching the age at the
    ## start of the interval for each projection interval. Years as columns,
    ## ages as rows. The first row should be nL0/(age.int*l0). The last row is
    ## survival for age.int years in the open interval.
    ##
    ## - fert:
    ## Matrix of (annual) age-specific fertility rates (NOT yet multiplied by
    ## age.int) for each projection interval. Years as columns, ages as rows.
    ##
    ## - srb:
    ## Matrix of sex ratios at birth. Should be of dim 1 x proj.steps.
    ##
    ## - mig.type:
    ## Indicates the type of migration data. Character, one of "prop.prev.pop",
    ## "net.count".
    ##
    ## - mig:
    ## List with up to two components: 'female' and 'male' each a matrix of
    ## age-specific proportions or counts for each projection interval,
    ## beginning with baseline. Years as columns, ages as rows.
    ##
    ## - proj.steps:
    ## The number of time periods to project forward.
    ##
    ## - age.int:
    ## Size of projection intervals (years).
    ##
    ## - label.dims:
    ## Should output have dimnames set? (cosmetic).
    ##
    ## - base.year:
    ## Label of baseline year (cosmetic).
    ##
    ## - first.age:
    ## Label of first age group (cosmetic).
    ##
    ## - do.checks:
    ## Do checks? Recommended.
    ##
    ## - verb:
    ## Logical. Controls messages. TRUE == messages, FALSE == no messages.
    ##
    ## - return.list:
    ## Logical. If female only projection, should a list be returned with one
    ## component called 'female'?
    ##
    ##
    ## OUTPUTS
    ##
    ## A list with up to two components: 'female' and 'male', each a matrix of
    ## age-specific population counts, years as columns, ages as rows.

{

    ## -------* Housekeeping

    ## Two-sexes or one?
    is2sex <- identical(length(pop), 2L)


    ## -------** Checks

    if(do.checks) {

        ## -------*** Check class of inputs

        ## Make sure pop, surv and mig are lists w correct components
        if(!is2sex) {
            if(!is.list(pop)) pop <- list(female = as.matrix(pop, ncol = 1))
            else if(is.null(dim(pop[["female"]])))
                pop <- list(female = as.matrix(pop[["female"]], ncol = 1))
            if(!is.list(surv)) surv <- list(female = as.matrix(surv))
            if(!is.list(mig)) mig <- list(female = as.matrix(mig))
        } else {
            ## Make sure pop is a matrix
            if(any(sapply(pop, function(z) is.null(dim(z)))))
                pop <- lapply(pop, function(z) as.matrix(z, ncol = 1))
        }

        if(is2sex) {
            ## Female and male components of input lists must be matrices with
            ## matching dimensions
            invisible(lapply(list(pop, surv, mig), FUN = function(z) {
                if(!identical(dim(z[["female"]]), dim(z[["male"]]))) stop("Male and female components must be of equal dimension within 'pop', 'surv' and 'mig'")
            }
                             )
                      )
        }

        ## Make sure srb is a matrix
        if(identical(length(srb), 1L)) srb <- matrix(srb, ncol = proj.steps)
        else if(!identical(ncol(srb), proj.steps)) stop("!identical(ncol(srb), proj.steps)")


        ## -------*** Check dimensions

        ## Number of years and ages in inputs must match
        input.dims <-
            lapply(list(fert = fert, surv = surv[["female"]], mig = mig[["female"]]
                        ,pop = pop[["female"]], srb = srb), FUN = "dim"
                   )

        input.ages <- sapply(input.dims[names(input.dims) != "srb"], "[[", 1)
        input.ages[["surv"]] <- input.ages[["surv"]] - 1
        if(!all.equal(input.ages - input.ages[[1]], rep.int(0, length(input.ages))
                      ,check.attributes = FALSE)) stop("Number of age groups in inputs not equal")

        input.years <-
            sapply(input.dims[names(input.dims) != "pop"], "[[", 2)
        if(!all.equal(input.years - input.years[[1]], rep.int(0, length(input.years))
                      ,check.attributes = FALSE)) stop("Number of years in inputs not equal")

    }


    ## Constants
    n.age.grps <- nrow(pop[["female"]])
    n.surv <- nrow(surv[["female"]])

    ## Derive proj.steps from ncol(fert)
    if(missing(proj.steps)) proj.steps <- ncol(fert)

    ## Halve the migrants
    half.mig <- lapply(mig, "*", 0.5)


    ## -------* Loop over number of required time periods

    ## Initialize pop.matrix and Leslie matrices
    if(is2sex) {
        pop.lis <-
            list(female = matrix(0, nrow = n.age.grps, ncol = 1 + proj.steps)
                 ,male = matrix(0, nrow = n.age.grps, ncol = 1 + proj.steps))
        pop.lis[["female"]][,1] <- pop[["female"]]
        pop.lis[["male"]][,1] <- pop[["male"]]
    lesM.female <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)
    lesM.male <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)
    } else {
        pop.lis <-
            list(female = matrix(0, nrow = n.age.grps, ncol = 1 + proj.steps))
        pop.lis[["female"]][,1] <- pop[["female"]]
    lesM.female <- matrix(0, nrow = n.age.grps, ncol = n.age.grps)
    }


    ## -------* Project in loop

    for(i in 1:proj.steps) {

        ## -------** Net migrants

        if(match.arg(mig.type) == "prop.prev.pop") {
            half.net.numb.mig <-
                mapply(FUN = function(a, b) a[,i] * b[,i] * age.int
                       ,half.mig, pop.lis, SIMPLIFY = FALSE
                       )
        } else if(match.arg(mig.type) == "net.count") {
            half.net.numb.mig <-
                lapply(X = half.mig, FUN = function(z) z[,i])
        }


        ## -------** Females

        ## Leslie Matrix for females: diagonal, then bottom right entry
        diag(lesM.female[-1,]) <- surv[["female"]][2:(n.age.grps),i]
        lesM.female[n.age.grps,n.age.grps] <- surv[["female"]][n.surv,i]

        ## First row
        k <- 1/(1+srb[i]) * surv[["female"]][1,i] * 0.5
        dbl.fert <- age.int * (fert[,i] + c(fert[-1,i], 0) * surv[["female"]][-1,i])
        lesM.female[1,] <- k * dbl.fert

        ## Project
        pop.lis[["female"]][,i+1] <-
            lesM.female %*% (pop.lis[["female"]][,i] +
                             half.net.numb.mig[["female"]]
                             ) + half.net.numb.mig[["female"]]


        ## -------** Males

        if(is2sex) {
            ## Leslie Matrix for males: diagonal, then bottom right entry.
            diag(lesM.male[-1,]) <- surv[["male"]][2:(n.age.grps),i]
            lesM.male[n.age.grps,n.age.grps] <- surv[["male"]][n.surv,i]
            pop.lis[["male"]][,i+1] <-
                lesM.male %*% (pop.lis[["male"]][,i] + half.net.numb.mig[["male"]]) +
                    half.net.numb.mig[["male"]]

            ## Add male births. Migration of mothers in this interval already
            ## taken care of in previous step.
            male.k <- srb[,i] / surv[["female"]][1,i] * surv[["male"]][1,i]
            pop.lis[["male"]][1,i+1] <- pop.lis[["male"]][1,i+1] +
                male.k * (pop.lis[["female"]][1,i+1] - half.net.numb.mig[["female"]][1])
        }

    }


    ## -------* OUTPUT

    ## Add dim names
    if(label.dims) {
        ages <- seq(from = first.age
                    ,to = age.int*(nrow(as.matrix(pop.lis[["female"]]))-1)
                    ,by = age.int)
        years <- (0:proj.steps) * age.int + as.numeric(base.year)
        for(i in 1:length(pop.lis))
            dimnames(pop.lis[[i]]) <- list(age = ages, year = years)
    }


    ## Return
    if(!is2sex && !return.list) return(pop.lis[[1]])
    else return(pop.lis)

}


## ###
## ### * EXAMPLES
## ###
## ################################################################################

## srb <- 1.05
## mig.type <- "prop.prev.pop"
## age.int <- 5
## label.dims <- TRUE
## base.year <- "1960"
## first.age <- 0
## do.checks <- TRUE
## verb <- FALSE

## ###
## ### No Migration, same baseline pop and survival for male and female, SRB = 1
## ###

## pop <- list(female = c(10, 12, 9, 4), male = c(10, 12, 9, 4))
## fert <- matrix(c(0, 3, 1.5, 0), nrow = 4, ncol = 12)
## proj.steps <- ncol(fert)
## surv0 <- list(female = matrix(rep(c(0.7, 0.95, 0.8, 0.75, 0.02), 12), ncol = 12)
##              ,male = matrix(rep(c(0.7, 0.95, 0.8, 0.75, 0.02), 12), ncol = 12)
##              )
## mig0 <- list(female = matrix(0, ncol = 12, nrow = 4)
##             ,male = matrix(0, ncol = 12, nrow = 4)
##             )

## eg.proj0 <- ccmp.femDom.apr20(pop, surv0, fert, srb = 1, mig0
##                              ,label.dims = TRUE, age.int = 1
##                              ,base.year = "1960"
##                              )

## all(eg.proj0[[1]] - eg.proj0[[2]] == 0)

## ## -> TRUE


## ###
## ### Change only SRB
## ###

## eg.proj1 <- ccmp.femDom.apr20(pop, surv0, fert, srb = 2, mig0
##                              ,label.dims = TRUE, age.int = 1
##                              ,base.year = "1960"
##                              )

## all((eg.proj1[[1]][1,] == 0.5 * eg.proj1[[2]][1,])[-1])

## ## -> TRUE


## ###
## ### Different everything
## ###

## pop <- list(female = c(10, 12, 9, 4), male = c(9, 8, 10, 2))
## fert <- cbind(matrix(c(0, 3, 1.5, 0, 0, 3, 1.5, 0), ncol = 2)
##               ,matrix(rep(c(0, 3, 1.5, 0), 10), ncol = 10)
##               )
## proj.steps <- ncol(fert)
## surv <- list(female = cbind(matrix(c(0.5, 0.9, 1, .9, .02, 0.7, 0.95
##              ,0.8, 0.75, 0.02), ncol = 2)
##              ,matrix(rep(c(0.7, 0.95, 0.8, 0.75, 0.02), 10), ncol = 10)
##              )
##              ,male = cbind(matrix(c(0.5, 0.9, 0.9, .7, .01, 0.7, 0.9
##               ,0.81, 0.7, 0.01), ncol = 2)
##               ,matrix(rep(c(0.7, 0.9, 0.81, 0.7, 0.01), 10), ncol = 10)
##               )
##              )
## mig <- list(female = cbind(matrix(c(1, 3, -2, 0, 1, 2, -1, 0), ncol = 2)
##             ,matrix(rep(c(1,2,-1,0), 10), ncol = 10)) / 100
##             ,male = cbind(matrix(c(2, 2, -3, 0, 2, 2, -2, 0), ncol = 2)
##              ,matrix(rep(c(2, 2, -2, 0), 10), ncol = 10)
##              ) / 100
##             )

## eg.proj <- ccmp.femDom.apr20(pop, surv, fert, srb = 1.05, mig
##                              ,label.dims = TRUE, age.int = 1
##                              ,base.year = "1960"
##                              )

## par(mfrow = c(1,2))
## for(i in 1:2) {
##     matplot(t(eg.proj[[i]]), type = "b", main = names(eg.proj)[[i]]
##             ,xlab = "period", ylab = "pop size"
##             )
## }

## par(mfrow = c(1,2))
## for(i in 1:2) {
##     matplot(eg.proj[[i]], type = "b", main = names(eg.proj)[[i]]
##             ,xlab = "age group", ylab = "pop size"
##             )
## }


## ###
## ### ** Compare female part of projection with projection using old female-only
## ###    function
## ###
## ################################################################################

## source("T:/PPGp/RPrograms/CCMPProjection/ccmp_female_20120310.R")

## fem.proj <-
##     ccmp.female.mar10(pop = pop[["female"]], surv = surv[["female"]]
##                       ,fert = fert, srb = 1.05, mig = mig[["female"]]
##                       ,age.int = 1, label.dims = TRUE, base.year = "1960")

## twoSex.proj.fem <-
##     ccmp.femDom.apr20(pop, surv, fert, srb = 1.05, mig
##                              ,label.dims = TRUE, age.int = 1
##                              ,base.year = "1960"
##                              )$female

## all(fem.proj == twoSex.proj.fem)

## ## -> TRUE


## ###
## ### ** Female only
## ###
## ################################################################################

## popF <- list(female = c(10, 12, 9, 4))
## fert <- cbind(matrix(c(0, 3, 1.5, 0, 0, 3, 1.5, 0), ncol = 2)
##               ,matrix(rep(c(0, 3, 1.5, 0), 10), ncol = 10)
##               )
## proj.steps <- ncol(fert)
## survF <- cbind(matrix(c(0.5, 0.9, 1, .9, .02, 0.7, 0.95
##              ,0.8, 0.75, 0.02), ncol = 2)
##              ,matrix(rep(c(0.7, 0.95, 0.8, 0.75, 0.02), 10), ncol = 10)
##              )
## migF <- cbind(matrix(c(1, 3, -2, 0, 1, 2, -1, 0), ncol = 2)
##             ,matrix(rep(c(1,2,-1,0), 10), ncol = 10)) / 100

## eg.projF <- ccmp.femDom.apr20(popF, survF, fert, srb = 1.05, migF
##                              ,label.dims = TRUE, age.int = 1
##                              ,base.year = "1960"
##                              )

## identical(eg.proj$female, eg.projF$female)
