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
################################################################################

ccmp.femDom <-
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
