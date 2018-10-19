################################################################################
###
###  DESC:          Transform simulation output (fert rates, surv props,
###                 mig props) to /counts/ of births, deaths and migrants.
###
###  DATE ORGINAL:  9th January 2012
###
###  AUTHOR:        Mark C Wheldon
###
################################################################################

### * SYNOPSIS
################################################################################

## These functions take a single set of age-specific fertility rates, survival
## proportions, migration proportions and baseline population counts used in a
## population projection and decomposes the projection into the number of births,
## deaths and net migrants.

### * FUNCTIONS
################################################################################

### ** Net number of migrants from n1, n2 and Leslie matrix
################################################################################

netMig.n1n2Leslie <- function(n1, n2, L) {
    ## ------------------------------------------------------------
    ##
    ## PURPOSE:
    ##
    ## Given two age structures and a leslie matrix defining a
    ## projection from the first to the second, back-out the net
    ## number of migrants
    ##
    ##
    ## ARGUMENTS:
    ##
    ##   n1      :  Population count vector at time t
    ##   n2      :  Population count vector at time t + delta
    ##   L       :  Leslie matrix used to get population at t + delta
    ##
    ##
    ## METHOD:
    ##
    ## Invert n2 = L(n1 + 0.5 mig) + (0.5)*mig
    ## Can get proportions by pre-multiplying output by 'solve(diag(n1))'
    ##
    ##------------------------------------------------------------

    ## Make sure inputs are of correct form
    n1 <- as.numeric(n1)
    n2 <- as.numeric(n2)
    L <- as.matrix(L)

    return(2 * solve(L + diag(nrow(L))) %*% (n2 - L %*% n1))
}


### ** Make leslie matrix from fert rates, surv props, pop counts
################################################################################

leslie <- function(pop, surv, fert, srb = 1.05, age.int = 5, label.dims = FALSE)

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


### ** Counts from Leslie Matrix and Population Vectors
################################################################################


##' Convert CCMPP parameters (rates, proportions, etc.) into counts
##'
##' Converts the quantities used in population projection (fert rates
##' suvival proportions, migration proportions, counts) into counts of
##' births, deaths and net number of migrants.
##'
##' @param results.recon List with components \code{fert.rate.mcmc},
##'     \code{surv.prop.mcmc}, \code{mig.prop.mcmc},
##'     \code{baseline.count.mcmc}, \code{pop.count.mcmc},
##'     \code{srb.mcmc}. All but the last are lists of mcmc objects with two
##'     components named 'female' and 'male'. Column names are
##'     yyyy.aa.aa (year.age-start.age-end). \code{mig.prop.mcmc} is
##'     assumed to be the average _annual_ proportion. It is
##'     multiplied by age.int to get migration for the entire
##'     interval. \code{pop.count.mcmc} does not contain
##'     baseline.count.mcmc. These should all have the same number of
##'     rows. The output will have this many rows too, otherwise the
##'     same number of rows as the smallest of these inputs.
##'
##' \code{srb.mcmc} is a numeric vector, length 1, or an mcmc/matrix
##' with one column per period (usually 5-year period).
##'
##' @param sep.factors List of separation factors used to separate
##'     cohort deaths into period deaths with components names
##'     'female' and 'male. These can be a scalar, in which case the
##'     same factor is used for all ages, year, or a _vector_ with
##'     names yyyy.aa indicating the 'y'ear and 'a'ge. Age goes from 0
##'     to the open ended interval in the population counts (not the
##'     survival proportions).
##'
##' For sep factor k_a, {delta}_D_a^c = k_a * D_a^c, i.e., k_a is the
##' fraction of cohort deaths for cohort aged 'a' that occurs in the
##' first triangle of the cohort parallelogram. sep factor for oldest
##' age group is currently re-set to 1.
##'
##' @param name.pref String to prepend to output objects.
##' @param name.suf String to append to output objects.
##' @param outputs Character vector listing the counts required as
##'     output. Counts are calculated regardless but are not output
##'     unless requested.
##' @return The requested outputs are matrices and are placed in the
##'     Global environment.
##' @author Mark C. Wheldon
##' @export post.process.recon
##' @import ff
post.process.recon <- function(results.recon,
             sep.factors = list(female = 0.5, male = 0.5),
             outputs = c("birth.count", "total.birth.count", "total.birth.death.count", "death.count", "mig.count", "total.mig.count", "mort.rate", "mig.rate", "person.years", "cohort.nq0", "period.nq0", "IMR"),
             name.pref = character(0), name.suf = character(0)
             ) {

    message("\nCONVERTING INPUT PARAMETERS TO COUNTS")


    ## -------* House keeping

    ## ------- ** Inputs

    results.params.post <-
        paste0(c("fert.rate", "surv.prop", "mig.prop",
                 "baseline.count", "srb",
                 "lx"),
               ".mcmc")
    results.params.prior <-
        paste0(c("fert.rate", "surv.prop", "mig.prop",
                 "baseline.count", "srb",
                 "pop.count"),
               ".prior")

    if(all(results.params.post %in% names(results.recon))) {

        results.are.post <- TRUE

        fert.rate.mcmc <- results.recon$fert.rate.mcmc
        surv.prop.mcmc <- results.recon$surv.prop.mcmc
        mig.prop.mcmc <- results.recon$mig.prop.mcmc
        baseline.count.mcmc <- results.recon$baseline.count.mcmc
        pop.count.mcmc <- results.recon$lx.mcmc
        srb.mcmc <- results.recon$srb.mcmc

    } else if(all(results.params.prior %in% names(results.recon))) {

        results.are.post <- FALSE

        fert.rate.mcmc <- results.recon$fert.rate.prior
        surv.prop.mcmc <- results.recon$surv.prop.prior
        mig.prop.mcmc <- results.recon$mig.prop.prior
        baseline.count.mcmc <- results.recon$baseline.count.prior
        pop.count.mcmc <- results.recon$pop.count.prior
        srb.mcmc <- results.recon$srb.prior
    }

    nF <- ncol(mig.prop.mcmc[[1]])

    sep.factors <- lapply(sep.factors, function(z, len = nF) {
        if(length(z) > 1) z[1:len]
        else z
        })


    ## ------- ** CCMPP

    ccmp.f <- ccmp.femDom


    ## -------** Basic check of female/male lists

    ## Female and male components must have same dim
    check.FM.dims <- function(obj) {
        if(!is.null(dim(obj[["female"]])) && !is.null(dim(obj[["male"]]))) {
            if(!identical(dim(obj[["female"]]), dim(obj[["male"]])))
                stop("Dimensions of female and male components of "
                     ,deparse(substitute(obj), " not equal"))
            if(!identical(dimnames(obj[["female"]]), dimnames(obj[["male"]])))
                warning("Dimnames of female and male components of "
                        ,deparse(substitute(obj)
                                 ," not equal. Dimnames of female component will be used."
                                 ))
        } else {
            if(!identical(length(obj[["female"]]), length(obj[["male"]])))
                stop("Lengths of female and male components of "
                     ,deparse(substitute(obj), " not equal"))
            if(!identical(names(obj[["female"]]), names(obj[["male"]])))
                warning("Names of female and male components of "
                        ,deparse(substitute(obj)
                                 ," not equal. Names of female component will be used."
                                 ))
        }
    }
    check.FM.dims(surv.prop.mcmc)
    check.FM.dims(mig.prop.mcmc)
    check.FM.dims(baseline.count.mcmc)
    check.FM.dims(pop.count.mcmc)
    check.FM.dims(sep.factors)


    ## -------** Make all mcmc inputs ff objects

    file.backed  <-  FALSE

    if(file.backed) {

        make.ff <- function(obj) {
            if(is.recursive(obj)) {
                if(!all(sapply(obj, "inherits", "ff"))) {
                    assign(deparse(substitute(obj))
                          ,value = list(female = ff(obj[["female"]]
                                                   ,dim = c(nrow(obj[["female"]]), ncol = ncol(obj[["female"]]))
                                                   ,dimnames = list(NULL, colnames(obj[["female"]]))
                                                    )
                                       ,male = ff(obj[["male"]]
                                                 ,dim = c(nrow(obj[["male"]]), ncol = ncol(obj[["male"]]))
                                                 ,dimnames = list(NULL, colnames(obj[["male"]]))
                                                  )
                                        )
                         , envir = .GlobalEnv
                           )
                }
            } else if(!inherits(obj, "ff")) {
                assign(deparse(substitute(obj))
                      ,value = ff(obj, dim = c(nrow(obj), ncol(obj))
                                 ,dimnames = list(NULL, colnames(obj))
                                  )
                     , envir = .GlobalEnv
                       )
            }
        }
        make.ff(fert.rate.mcmc)
        make.ff(surv.prop.mcmc)
        make.ff(mig.prop.mcmc)
        make.ff(baseline.count.mcmc)
        make.ff(pop.count.mcmc)
        make.ff(srb.mcmc)

    }


    ## -------** Strip any "+"s from column names

    colnames(fert.rate.mcmc) <-
        sapply(strsplit(colnames(fert.rate.mcmc), "\\+"), "[[", 1)
    for(s in c("female", "male")) {
        colnames(surv.prop.mcmc[[s]]) <-
            sapply(strsplit(colnames(surv.prop.mcmc[[s]]), "\\+"), "[[", 1)
        colnames(mig.prop.mcmc[[s]]) <-
            sapply(strsplit(colnames(mig.prop.mcmc[[s]]), "\\+"), "[[", 1)
        colnames(baseline.count.mcmc[[s]]) <-
            sapply(strsplit(colnames(baseline.count.mcmc[[s]]), "\\+"), "[[", 1)
        colnames(pop.count.mcmc[[s]]) <-
            sapply(strsplit(colnames(pop.count.mcmc[[s]]), "\\+"), "[[", 1)
    }


    ## -------** Set up required indices into the inputs

    ## surv years
    surv.prop.years <-
        sapply(strsplit(colnames(surv.prop.mcmc[["female"]])
                        ,"\\."), FUN = function(z) z[[1]])
    surv.prop.ages <-
        sapply(strsplit(colnames(surv.prop.mcmc[["female"]])
                        ,"[^0-9]"), FUN = function(z) z[[2]])

    ## fert years
    fert.rate.years <-
        sapply(strsplit(colnames(fert.rate.mcmc)
                        ,"\\."), FUN = function(z) z[[1]])
    fert.rate.ages <-
        sapply(strsplit(colnames(fert.rate.mcmc)
                        ,"[^0-9]"), FUN = function(z) z[[2]])

    ## mig years
    mig.prop.years <-
        sapply(strsplit(colnames(mig.prop.mcmc[["female"]])
                        ,"\\."), FUN = function(z) z[[1]])
    mig.prop.ages <-
        sapply(strsplit(colnames(mig.prop.mcmc[["female"]])
                        ,"[^0-9]"), FUN = function(z) z[[2]])

    ## baseline count
    baseline.count.years <-
        sapply(strsplit(colnames(baseline.count.mcmc[["female"]])
                        ,"\\."), FUN = function(z) z[[1]])
    baseline.count.ages <-
        sapply(strsplit(colnames(baseline.count.mcmc[["female"]])
                        ,"[^0-9]"), FUN = function(z) z[[2]])

    ## pop count
    pop.count.years <-
        sapply(strsplit(colnames(pop.count.mcmc[["female"]])
                       ,"\\."), FUN = function(z) z[[1]])
    pop.count.ages <-
        sapply(strsplit(colnames(pop.count.mcmc[["female"]])
                       ,"[^0-9]"), FUN = function(z) z[[2]])

    ## Concatenate baseline and pop count to get a single matrix with population
    ## counts
    if(!any(unique(baseline.count.years) %in% unique(pop.count.years))) {
        if(file.backed) {
            pc2 <- list(female = ff(rep(0, length(baseline.count.mcmc[["female"]]) +
                                           length(pop.count.mcmc[["female"]]))
                                   ,dim = c(nrow(baseline.count.mcmc[["female"]])
                                           ,ncol(pop.count.mcmc[["female"]]) +
                                            ncol(baseline.count.mcmc[["female"]])
                                            )
                                   ,dimnames = list(NULL, c(colnames(baseline.count.mcmc[["female"]])
                                                           ,colnames(pop.count.mcmc[["female"]])
                                                            ))
                                    )
                       ,male = ff(rep(0, length(baseline.count.mcmc[["male"]]) +
                                         length(pop.count.mcmc[["male"]]))
                                 ,dim = c(nrow(baseline.count.mcmc[["male"]])
                                         ,ncol(pop.count.mcmc[["male"]]) + ncol(baseline.count.mcmc[["male"]])
                                          )
                                 ,dimnames = list(NULL, c(colnames(baseline.count.mcmc[["male"]])
                                                         ,colnames(pop.count.mcmc[["male"]])
                                                          ))
                                  )
                        )
            for(s in c("female", "male")) {
                pc2[[s]][,1:ncol(baseline.count.mcmc[[s]])] <- baseline.count.mcmc[[s]][,]
                pc2[[s]][,(ncol(baseline.count.mcmc[[s]]) + 1):ncol(pc2[[s]])] <-
                    pop.count.mcmc[[s]][,]
            }
            pop.count.mcmc <- pc2
            rm(pc2)
            rm(baseline.count.mcmc)
        }  else pop.count.mcmc <-
                    list(female = cbind(baseline.count.mcmc[["female"]], pop.count.mcmc[["female"]])
                        ,male = cbind(baseline.count.mcmc[["male"]], pop.count.mcmc[["male"]])
                         )
    }

    ## Fix colnames for populations. Should just be yyyy.ageStart
    for(s in c("female", "male")) {
        colnames(pop.count.mcmc[[s]]) <-
            sapply(strsplit(colnames(pop.count.mcmc[[s]]), "[^0-9]"), function(z) {
                paste(z[1:2], collapse = ".")
            }
            )
    }

    ## Index for augmented population years
    pop.count.years <- sapply(strsplit(colnames(pop.count.mcmc[["female"]])
                                      ,"\\."), FUN = function(z) z[[1]])
    pop.years <- unique(pop.count.years)
    pop.count.ages <- sapply(strsplit(colnames(pop.count.mcmc[["female"]])
                                     ,"[^0-9]"), FUN = function(z) z[[2]])


    ## -------** Separation factors

    if(length(sep.factors[["female"]]) == 1) {
        message("Separation factors constant over age and time (and same for all iterations)")
        sep.factors <-
            list(female = rep(sep.factors[["female"]], ncol(mig.prop.mcmc[["female"]]))
                 ,male = rep(sep.factors[["male"]], ncol(mig.prop.mcmc[["male"]]))
                 )
        for(s in c("female", "male")) names(sep.factors[[s]]) <- colnames(mig.prop.mcmc[[s]])
    } else if(length(sep.factors[["female"]]) == length(unique(mig.prop.ages))) {
        message("Separation factors constant over time but differ by age (but same for all iterations)")
        sep.factors <-
            list(female = rep(sep.factors[["female"]]
                 ,length(unique(mig.prop.years[["female"]])))
                 ,male = rep(sep.factors[["male"]]
                 ,length(unique(mig.prop.years[["male"]])))
                 )
        for(s in c("female", "male")) names(sep.factors[[s]]) <- colnames(mig.prop.mcmc[[s]])
    } else if(length(sep.factors[["female"]]) == ncol(mig.prop.mcmc[["female"]])) {
        message("Separation factors differ by time and age (but same for all iterations)")
    } else stop("Do not know how to handle 'sep.factors' of this shape")

    sep.factors.years <-
        sapply(strsplit(names(sep.factors[["female"]])
                        ,"\\."), FUN = function(z) z[[1]])
    sep.factors.ages <-
        sapply(strsplit(names(sep.factors[["female"]])
                        ,"[^0-9]"), FUN = function(z) z[[2]])


    ## -------** Checks

    ## Years had better be the same
    if(!identical(unique(surv.prop.years)
                  ,unique(fert.rate.years)) ||
       !identical(unique(surv.prop.years)
                  ,unique(mig.prop.years)))
        stop("years must be the same for fert.rate, surv.prop and mig.prop")

    VR.years <- unique(surv.prop.years)

    ## Pop years had better contain vital rate years
    if(!all(unique(surv.prop.years) %in% pop.years))
        stop("vital rate years not in population years")

    ## Ages had better be the same
    if(!identical(head(unique(surv.prop.ages), -1)
                 ,unique(mig.prop.ages))) stop("ages not all the same")

    if(!identical(unique(mig.prop.ages)
                 ,unique(pop.count.ages)
                  )) stop("ages not all the same")

    if(!all(unique(fert.rate.ages) %in% unique(surv.prop.ages)))
        stop("fert.rate.mcmc ages not subset of surv.prop.mcmc ages")

    ## Separation factors
    if(!identical(sep.factors.years, mig.prop.years))
        stop("separation factor years must be the same as VR years")
    if(!identical(sep.factors.ages, mig.prop.ages))
        stop("separation factor ages must be same as mig ages")

    ## Calculate age/year intervals. Had better be consistent.
    age.diff <- sum(abs(diff(diff(as.numeric(unique(fert.rate.ages))))))
    year.diff <- sum(abs(diff(diff(as.numeric(unique(fert.rate.years))))))
    year.diff.pop <- sum(abs(diff(diff(as.numeric(unique(pop.count.years))))))
    if(!all(identical(age.diff, year.diff, year.diff.pop)))
        stop("age and year intervals not identical")

    age.int <- diff(as.numeric(unique(fert.rate.ages)))[1]

    ## Report years and ages
    message(paste(c("Vital rate years are: ", unique(VR.years)), collapse = " "))
    message(paste(c("Pop count years are: ", unique(pop.years)), collapse = " "))
    message(paste(c("Ages are: ", unique(surv.prop.ages)), collapse = " "))
    message(paste(c("Ages of non-zero fertility are: "
                    , unique(fert.rate.ages)), collapse = " "))
    message("Un-tested if these are not in ascending order")


    ## -------* CALCULATE

    ## -------** Some more useful objects

    ## The population counts at the start of the vital rate year periods

    pop.count.VRyears <-
            list(female = ff(pop.count.mcmc[["female"]][, pop.count.years %in% VR.years]
               ,dim = c(nrow(pop.count.mcmc[["female"]]), sum(pop.count.years %in% VR.years))
                 ,dimnames = list(NULL
                  ,colnames(pop.count.mcmc[["female"]])[pop.count.years %in% VR.years])
               )
                 ,male = ff(pop.count.mcmc[["male"]][, pop.count.years %in% VR.years]
               ,dim = c(nrow(pop.count.mcmc[["male"]]), sum(pop.count.years %in% VR.years))
                 ,dimnames = list(NULL
                  ,colnames(pop.count.mcmc[["male"]])[pop.count.years %in% VR.years])
               )
                 )


    ## -------** Counts

    ## -------*** Migration

    ## Migration count is just mig.prop * pop.count because it is input as a
    ## period measure. The sharing of migrants across the projection interval
    ## (half at the beginning, half at the end) is done in the projection
    ## itself. The inputs /are/ the age-period quantities.

    message("Calculating migration counts. Counts are for the whole
    'age.int' year long interval...")

    mig.count.mcmc <- list()

    mig.count.mcmc[["female"]] <-
            ffrowapply(mig.prop.mcmc[["female"]][i1:i2,] *
                       pop.count.VRyears[["female"]][i1:i2,] * age.int
                       ,X = mig.prop.mcmc$female
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
        mig.count.mcmc[["male"]] <-
            ffrowapply(mig.prop.mcmc[["male"]][i1:i2,] *
                       pop.count.VRyears[["male"]][i1:i2,] *
                       age.int
                       ,X = mig.prop.mcmc$male
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
    for(s in c("female", "male")) dimnames(mig.count.mcmc[[s]]) <-
        dimnames(mig.prop.mcmc[[s]])

    ## Total migration: sum over ages.
    message("Calculating total mig counts...")

    total.mig.count.mcmc <- list()

    dummy <- t(apply(mig.count.mcmc[["female"]][1,,drop=FALSE], 1, function(z) {
            tapply(z, mig.prop.years, "sum")}
                         ))           # nrow doesn't matter
        total.mig.count.mcmc[["female"]] <-
            ffrowapply(t(apply(mig.count.mcmc[["female"]][i1:i2,,drop=FALSE], 1
                               , function(z) {
                tapply(z, mig.prop.years, "sum")}
                               ))
                       ,X = ff(0, dim = c(nrow(mig.count.mcmc[["female"]]), ncol(dummy)))
                       ,RETURN = TRUE
                       )
        colnames(total.mig.count.mcmc[["female"]]) <- colnames(dummy)
        total.mig.count.mcmc[["male"]] <-
            ffrowapply(t(apply(mig.count.mcmc[["male"]][i1:i2,,drop=FALSE], 1, function(z) {
                tapply(z, mig.prop.years, "sum")}
                               ))
                       ,X = ff(0, dim = c(nrow(mig.count.mcmc[["male"]]), ncol(dummy)))
                       ,RETURN = TRUE
                       )
        colnames(total.mig.count.mcmc[["male"]]) <- colnames(dummy)


    ## -------*** Births

    ## Birth count is just population aged 0--5 next year, minus
    ## migrants, reverse survived. But don't calculate like this. Just
    ## re-project the 0--5 age group without mortality and don't add
    ## half of the migrants at the end.

    ## Don't add migration at beginning either because this is the count of
    ## _births_, not births and half of the migrants.

    ## birth.count is the
    ## number born in the interval (allowed to survive). Not the number left at
    ## the end of the interval (that would be _{0}{N}_{5}(t+5) which is in pop.count).

    message("Calculating birth counts...")

    ## Make a matrix with all ages and fertility rates filled in, zero elsewhere.
    fert.full <- ff(0.0, dim = c(nrow(fert.rate.mcmc), ncol(mig.prop.mcmc[["female"]]))
                    ## ,dimnames = dimnames(mig.prop.mcmc[["female"]])
                    )
    fert.full[, mig.prop.ages %in% fert.rate.ages] <- fert.rate.mcmc[,]


    ## Make the first row of the Leslie matrix (Preston et al. p.130)
    ## EXCLUDING survival ratio to get ALL births.

    ## Initialize
    fert.LeslieRow <- ff(0, dim = c(nrow(fert.rate.mcmc), length(mig.prop.years))
           ,dimnames = list(NULL, colnames(mig.prop.mcmc[["female"]]))
           )

    ## Need to do year by year
    for(y in VR.years) {
        ffrowapply(fert.LeslieRow[i1:i2,mig.prop.years==y] <-
                   (1/2) * age.int *
                   (fert.full[i1:i2,mig.prop.years==y,drop=FALSE] +
                    cbind(fert.full[i1:i2,mig.prop.years==y,drop=FALSE][,-1,drop=FALSE]
                          ,rep(0, i2 - i1 + 1)
                          ) *
                    cbind(surv.prop.mcmc[["female"]][i1:i2,surv.prop.years==y,drop=FALSE][,-1,drop=FALSE])
                    )
                   ,X = fert.LeslieRow[,mig.prop.years==y,drop=FALSE]
                   ,USE.NAMES = TRUE
                   ,RETURN = TRUE
                   )
    }

    ## Calculate births by projection: multiply _female_ population by
    ## fert.LeslieRow and appropriate sex ratio
    if(identical(length(srb.mcmc), 1L)) {
        birth.count.mcmc <-
            list(female = ffrowapply(fert.LeslieRow[i1:i2,,drop=FALSE] *
                 (1/(1+srb.mcmc[,])) *
                 (pop.count.VRyears[["female"]][i1:i2,,drop=FALSE] +
                  0.5 * mig.count.mcmc[["female"]][i1:i2,,drop=FALSE])
                 ,X = pop.count.VRyears[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(fert.LeslieRow[i1:i2,,drop=FALSE] *
                  (srb.mcmc[,]/(1+srb.mcmc[,])) *
                  (pop.count.VRyears[["female"]][i1:i2,,drop=FALSE] +
                   0.5 * mig.count.mcmc[["female"]][i1:i2,,drop=FALSE])
                  ,X = pop.count.VRyears[["female"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )
    } else {
        birth.count.mcmc <-
            list(female = ffrowapply(fert.LeslieRow[i1:i2,] * (1/(1+srb.mcmc[i1:i2,drop=FALSE])) *
                 (pop.count.VRyears[["female"]][i1:i2,,drop=FALSE] +
                  0.5 * mig.count.mcmc[["female"]][i1:i2,,drop=FALSE])
                 ,X = pop.count.VRyears[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(fert.LeslieRow[i1:i2,] *
                  (srb.mcmc[i1:i2,drop=FALSE]/(1+srb.mcmc[i1:i2,drop=FALSE])) *
                  (pop.count.VRyears[["female"]][i1:i2,,drop=FALSE] +
                   0.5 * mig.count.mcmc[["female"]][i1:i2,,drop=FALSE])
                  ,X = pop.count.VRyears[["female"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )
    }

    dimnames(birth.count.mcmc) <- dimnames(mig.prop.mcmc)


    ## -------*** Total Births

    ## Total births: sum over ages.

    message("Calculating total birth counts...")

    dummy <- t(apply(birth.count.mcmc[["female"]][1,,drop=FALSE], 1, function(z) {
            tapply(z, mig.prop.years, "sum")}
                         ))
        total.birth.count.mcmc <-
            list(female =
            ffrowapply(t(apply(birth.count.mcmc[["female"]][i1:i2,,drop=FALSE]
                               ,1, function(z) {
                                   tapply(z, mig.prop.years, "sum")}
                               ))
                       ,X = ff(0, dim = c(nrow(birth.count.mcmc[["female"]])
                                  ,ncol(dummy)
                                  )
                        )
                       ,RETURN = TRUE
                       )
                 ,male =
            ffrowapply(t(apply(birth.count.mcmc[["male"]][i1:i2,,drop=FALSE]
                               ,1, function(z) {
                                   tapply(z, mig.prop.years, "sum")}
                               ))
                       ,X = ff(0, dim = c(nrow(birth.count.mcmc[["male"]])
                                  ,ncol(dummy)
                                  )
                        )
                       ,RETURN = TRUE
                       )
                 )
        for(s in c("female", "male")) colnames(total.birth.count.mcmc[[s]]) <-
            colnames(dummy)


    ## -------*** Deaths

    message("Calculating death counts...")

    ## Calculate the deaths TO EACH COHORT (except births).
    ## 'death.count.cohorts.mcmc' yyyy.x gives the number aged [x,x+n) in
    ## yyyy who die before yyyy + n

    ## i.e., make 'death.count.cohorts.mcmc' which has first three columns:
    ## +---------+---------+---------+-
    ## | yyyy.0  | yyyy.5  | yyyy.10 | ...
    ## +---------+---------+---------+-
    ## | D_0^c   | D_5^c   | D_10^c  | ...
    ## +---------+---------+---------+-
    ##
    ## where D_0^c is the number of deaths to the cohort aged [0,5) at exact
    ## year yyyy during the projection interval [yyyy, yyyy+5)

    death.count.cohorts.mcmc <-
        list(female =
            ffrowapply((pop.count.VRyears[["female"]][i1:i2,] + (0.5) *
                        mig.count.mcmc[["female"]][i1:i2,]) *
                       (1 - surv.prop.mcmc[["female"]][i1:i2, surv.prop.ages != 0])
                       ,X = pop.count.VRyears[["female"]]
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
             ,male = ffrowapply((pop.count.VRyears[["male"]][i1:i2,] + (0.5) *
                        mig.count.mcmc[["male"]][i1:i2,]) *
                       (1 - surv.prop.mcmc[["male"]][i1:i2, surv.prop.ages != 0])
                       ,X = pop.count.VRyears[["male"]]
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
             )

    ## Number of births that died /in the projection interval/.
    ## 'total.birth.death.count.mcmc' yyyy gives the number of births that
    ## died between yyyy and yyyy + n.
   total.birth.death.count.mcmc <-
            list(female = ffrowapply(total.birth.count.mcmc[["female"]][i1:i2,] *
                       (1 - surv.prop.mcmc[["female"]][i1:i2, surv.prop.ages=="0"])
                       ,X = total.birth.count.mcmc[["female"]]
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
                 ,male = ffrowapply(total.birth.count.mcmc[["male"]][i1:i2,] *
                       (1 - surv.prop.mcmc[["male"]][i1:i2, surv.prop.ages=="0"])
                       ,X = total.birth.count.mcmc[["male"]]
                       ,RETURN = TRUE
                       ,USE.NAMES = TRUE
                       )
                 )


    ## Redistribute

    ## Need to re-distribute the deaths that occur to cohort aged [x,x+n) in
    ## year yyyy between yyyy and yyyy + n (which is a parallelogram on
    ## Lexis surface) to the age-time square bounded by ages x and x+n and
    ## years yyyy and yyyy + n.

    ## The parallelograms stretch over two age-time squares: the upper half
    ## of the square [x,x+n), [yyyy,yyyy+n) and the lower half of the square
    ## [x+n,x+2n), [yyyy, yyyy+n). These are denoted _{\delta}{D}_{x}^c and
    ## _{\alpha}{D}_{x+5}^c, respectively.
    ## _{n}{D}_{x}^c = _{\delta}{D}_{x}^c + _{\alpha}{D}_{x+5}^c.

    ## Use 'separation factors' to determine how many deaths to give
    ## to each square. The square [x+n,x+2n), [yyyy, yyyy+n) gets nax/n *
    ## _{n}{D}_{x}^c, the other square gets the remainder (i.e., ndx - nax/n * ndx =
    ## _{n}{D}_{x}^c*(1-nax/n)). The period death count, ndx is
    ## ndx = nax/n * _{\delta}{D}_{x}^c + (1-nax/n) * _{\alpha}{D}_{x}^c
    ## (note different alpha term to that in expression for _{n}{D}_{x}^c).

    ## All deaths to proj-interval births occur in that
    ## age-range (don't get redistributed) because the projection only goes
    ## to t+5.

    ## Some indices, constants
    death.count.cohorts.ages <-
        sapply(strsplit(colnames(death.count.cohorts.mcmc[["female"]]), "\\.")
               ,"[[", 2)
    death.count.cohorts.years <-
        sapply(strsplit(colnames(death.count.cohorts.mcmc[["female"]]), "\\.")
               ,"[[", 1)
    total.birth.death.count.years <-
        sapply(strsplit(colnames(total.birth.death.count.mcmc[["female"]]), "\\.")
               ,"[[", 1)
    nages <- length(unique(death.count.cohorts.ages))
    oldest.age <-
        as.character(max(as.numeric(unique(death.count.cohorts.ages))))

    ## Initialize a matrix to store the number of deaths in the /previous/
    ## age cohort, fractioned appropriately. Thus:
    ##
    ## 'death.agePrevious' will have first three columns:
    ## +---------+------------+-----------+-
    ## | yyyy.B  | yyyy.0     | yyyy.5    | ...
    ## +---------+------------+-----------+-
    ## | D_B^c   | k_0*D_0^c  | k_5*D_5^c | ...
    ## +---------+------------+-----------+-
    ##
    ## where D_B^c is the number of deaths to the births that occur in the
    ## projection interval and k_x are the fractions to share deaths into
    ## the periods (k_B = 1). 'death.agePrevious' is so called because its
    ## columns correspond to the age prior to the age in the same column in
    ## 'death.count.cohorts.mcmc'. The column headers are correct for
    ## 'death.agePrevious'.

    ## Initialize:
    death.agePrevious <-
        list(female = ff(initdata = 0
             ,dim = c(nrow(death.count.cohorts.mcmc[["female"]])
              , ncol(death.count.cohorts.mcmc[["female"]]))
             ,dimnames = NULL
             )
             ,male = ff(initdata = 0
              ,dim = c(nrow(death.count.cohorts.mcmc[["male"]])
               , ncol(death.count.cohorts.mcmc[["male"]]))
              ,dimnames = NULL
              )
             )

    sep.factors2.ff <-
        list(female = ff(initdata = 1 - sep.factors[["female"]]
             ,dim = c(nrow(death.count.cohorts.mcmc[["female"]])
              ,length(sep.factors[["female"]]))
             ,bydim = c(2,1)
             )
             ,male = ff(initdata = 1 - sep.factors[["male"]]
              ,dim = c(nrow(death.count.cohorts.mcmc[["male"]])
               ,length(sep.factors[["male"]]))
              ,bydim = c(2,1)
              )
             )

    ## Fill:
    ## Need to do year-by-year because the birth deaths need to be interleaved
    ## into the storage matrix.
    for(y in total.birth.death.count.years) {

        ffrowapply(death.agePrevious[["female"]][i1:i2
                                                 ,death.count.cohorts.years==y] <-
                   cbind(total.birth.death.count.mcmc[["female"]][i1:i2,
                                                                  total.birth.death.count.years==y
                                                                  ,drop=FALSE]
                         ,death.count.cohorts.mcmc[["female"]][i1:i2,
                                                               death.count.cohorts.years == y
                                                               ,drop=FALSE][,-nages,drop=FALSE] *
                         sep.factors2.ff[["female"]][i1:i2,death.count.cohorts.years == y
                                                     ,drop=FALSE][,-nages,drop=FALSE]
                         )
                   ,X = death.agePrevious[["female"]][,death.count.cohorts.years==y
                    ,drop = FALSE]
                   ,RETURN = TRUE
                   ,USE.NAMES = TRUE
                   )

        ffrowapply(death.agePrevious[["male"]][i1:i2
                                               ,death.count.cohorts.years==y] <-
                   cbind(total.birth.death.count.mcmc[["male"]][i1:i2,
                                                                total.birth.death.count.years==y
                                                                ,drop=FALSE]
                         ,death.count.cohorts.mcmc[["male"]][i1:i2,
                                                             death.count.cohorts.years == y
                                                             ,drop=FALSE][,-nages,drop=FALSE] *
                         sep.factors2.ff[["male"]][i1:i2,death.count.cohorts.years == y
                                                   ,drop=FALSE][,-nages,drop=FALSE]
                         )
                   ,X = death.agePrevious[["male"]][,death.count.cohorts.years==y
                    ,drop = FALSE]
                   ,RETURN = TRUE
                   ,USE.NAMES = TRUE
                   )
    }


    ## for debugging: colnames would be (but don't bother setting them):
    ## --COMMENTED CODE-->|
    ## temp <- expand.grid(c("B", head(unique(death.count.cohorts.ages), -1))
    ##                   ,total.birth.death.count.years
    ##                     )
    ## for(s in c("female", "male")) colnames(death.agePrevious[[s]]) <-
    ##     paste(temp[,2], temp[,1], sep = ".")
    ## |<--

    ## The approximate number of deaths in the age range for the fixed
    ## period is sep.factor[yyyy.aa] of the number of deaths to the cohort
    ## plus (1 - sep.factor[yyyy.aa]) times of the number of deaths in the
    ## previous age. The matrices 'death.count.cohorts.mcmc' and
    ## 'death.agePrevious' should align correctly to make this operation a
    ## simple addition of matrices (addition is done element-by-element).

    ## All deaths of births in the projection interval are kept
    ## (half of the cohort deaths from the 0--5 age group have
    ## been given to the next age group anyway).

    ## All deaths to the open-ended interval are kept. The
    ## operation below halves them so remember to add half back.

    sep.factors.ff <-
        list(female = ff(sep.factors[["female"]]
             ,dim = c(nrow(death.count.cohorts.mcmc[["female"]])
              ,length(sep.factors[["female"]])
              )
             ,bydim = c(2,1)
             )
             ,male = ff(sep.factors[["male"]]
              ,dim = c(nrow(death.count.cohorts.mcmc[["male"]])
               ,length(sep.factors[["male"]])
               )
              ,bydim = c(2,1)
              )
             )
    death.count.mcmc <-
        list(female = ffrowapply(death.count.cohorts.mcmc[["female"]][i1:i2,] *
             sep.factors.ff[["female"]][i1:i2,] +
             death.agePrevious[["female"]][i1:i2,]
             ,X = death.count.cohorts.mcmc[["female"]]
             ,RETURN = TRUE
             ,USE.NAMES = TRUE
             )
             ,male = ffrowapply(death.count.cohorts.mcmc[["male"]][i1:i2,] *
              sep.factors.ff[["male"]][i1:i2,] +
              death.agePrevious[["male"]][i1:i2,]
              ,X = death.count.cohorts.mcmc[["male"]]
              ,RETURN = TRUE
              ,USE.NAMES = TRUE
              )
             )


    ## CHECK: conservation of deaths.
    ##      total number of death.count.cohort + births.deaths should equal
    ##      total death.count.mcmc.
    ## --COMMENTED CODE -->|
    ## ae <- all.equal(rowSums(cbind(total.birth.death.count.mcmc[["female"]][]
    ##                               ,death.count.cohorts.mcmc[["female"]][]))
    ##                 ,rowSums(death.count.mcmc[["female"]][])
    ##                 )
    ## if(!isTRUE(ae)){
    ##     plot(death.count.mcmc[["female"]][,1:17],  type = "b")
    ##     lines(c(total.birth.death.count.mcmc[["female"]][1]
    ##             , death.count.cohorts.mcmc[["female"]][,1:17]),type = "b", col = "red")
    ##     stop("CHECK: female deaths not conserved")
    ## }

    ## ae <- all.equal(rowSums(cbind(total.birth.death.count.mcmc[["male"]][]
    ##                               ,death.count.cohorts.mcmc[["male"]][]))
    ##                 ,rowSums(death.count.mcmc[["male"]][])
    ##                 )
    ## if(!isTRUE(ae)){
    ##     plot(death.count.mcmc[["male"]][,1:17],  type = "b")
    ##     lines(c(total.birth.death.count.mcmc[["male"]][1]
    ##             , death.count.cohorts.mcmc[["male"]][,1:17]),type = "b", col = "red")
    ##     stop("CHECK: male deaths not conserved")
    ## }
    ## |<--


    ## -------** Infant/Child mortality

    ## Calculate cohort and period nq0 and IMR (for age 0--age.int).

    if(any(c("cohort.nq0", "period.nq0", "IMR") %in% outputs)) {

        message("Calculating infant/child mortality measures...")

        ## Years beginning projection intervals for which number of births are
        ## available
        total.birth.count.years <-
            sapply(strsplit(colnames(total.birth.count.mcmc[["female"]]), "\\."), "[[", 1)

        ## Get deaths occur to cohorts between ages 0--5. These are not the
        ## deaths to the births; they are the deaths to the cohort aged 0--5
        ## at the beginning of the projection interval. We need the 'delta'
        ## portion of these deaths.
        death.count.cohort.delta.0 <-
            list(female = ffrowapply(death.count.cohorts.mcmc[["female"]][i1:i2
                 ,death.count.cohorts.ages=="0"
                 ,drop=FALSE] *
                 sep.factors[["female"]][sep.factors.ages == "0"]
                 ,X = death.count.cohorts.mcmc[["female"]][,death.count.cohorts.ages=="0"
                  ,drop=FALSE]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(death.count.cohorts.mcmc[["male"]][i1:i2
                  ,death.count.cohorts.ages=="0"
                  ,drop=FALSE] *
                  sep.factors[["male"]][sep.factors.ages == "0"]
                  ,X = death.count.cohorts.mcmc[["male"]][,death.count.cohorts.ages=="0"
                   ,drop=FALSE]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

        death.count.delta.years <-
            sapply(strsplit(colnames(death.count.cohort.delta.0[["female"]]), "\\.")
                   ,"[[", 1)

    }


    ## Cohort.nq0:

    if("cohort.nq0" %in% outputs) {

        ## This is
        ## ( _{\alpha}{D}_{0}^c[t,t+5) + k_0 * _{\delta}{D}_{x}^c[t+5,t) ) / no. births
        ## Can only calculate up to penultimate year in VR.years because we
        ## don't have _{\delta}{D}_{x}^c[t+5,t) for t > the last year.

        cohort.nq0.numerator <-
            list(female = ffrowapply(total.birth.death.count.mcmc[["female"]][i1:i2,-ncol(total.birth.death.count.mcmc[["female"]])] +
                 death.count.cohort.delta.0[["female"]][i1:i2,-1]
                 ,X = total.birth.death.count.mcmc[["female"]][,-ncol(total.birth.death.count.mcmc[["female"]])
                  ,drop = FALSE]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(total.birth.death.count.mcmc[["male"]][i1:i2,-ncol(total.birth.death.count.mcmc[["male"]])] +
                  death.count.cohort.delta.0[["male"]][i1:i2,-1]
                  ,X = total.birth.death.count.mcmc[["male"]][,-ncol(total.birth.death.count.mcmc[["male"]])
                   ,drop = FALSE]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

        cohort.nq0.mcmc <-
            list(female = ffrowapply(cohort.nq0.numerator[["female"]][i1:i2,] /
                 total.birth.count.mcmc[["female"]][i1:i2,-ncol(total.birth.death.count.mcmc[["female"]])]
                 ,X = cohort.nq0.numerator[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(cohort.nq0.numerator[["male"]][i1:i2,] /
                  total.birth.count.mcmc[["male"]][i1:i2,-ncol(total.birth.death.count.mcmc[["male"]])]
                  ,X = cohort.nq0.numerator[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

    }

    if("period.nq0" %in% outputs) {

        ## This is a bit complicated since it mixes cohorts. Let B_0[t,t+5)
        ## be the number of births during interval [t,t+5). Then
        ## period.nq0[t,t+5) = _{\alpha}{D}_{0}^c[t,t+5) / B_0[t,t+5) +
        ##                       (_{\alpha}{D}_{0}^c[t,t+5) - B_0[t,t+5)) / B_0[t,t+5) *
        ##                       _{\delta}{D}_{0}^c[t,t+5) / (_{\alpha}{D}_{0}^c[t-5,t) - B_0[t-5,t))

        ## Cannot do for the first time interval because we don't have
        ## _{\alpha}{D}_{0}^c[t-5,t) or B_0[t-5,t) for this interval.

        ## Ratios of deaths in projection interval to births in that
        ## interval (except first time interval)
        projection.deaths.per.births <-
            list(female = ffrowapply(total.birth.death.count.mcmc[["female"]][i1:i2,] /
                 total.birth.count.mcmc[["female"]][i1:i2,]
                 ,X = total.birth.death.count.mcmc[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(total.birth.death.count.mcmc[["male"]][i1:i2,] /
                  total.birth.count.mcmc[["male"]][i1:i2,]
                  ,X = total.birth.death.count.mcmc[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )


        ## Differences of deaths in projection interval to births in
        ## that interval (all intervals needed for this one)
        projection.births.minus.deaths <-
            list(female = ffrowapply(total.birth.count.mcmc[["female"]][i1:i2,] -
                 total.birth.death.count.mcmc[["female"]][i1:i2,]
                 ,X = total.birth.count.mcmc[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(total.birth.count.mcmc[["male"]][i1:i2,] -
                  total.birth.death.count.mcmc[["male"]][i1:i2,]
                  ,X = total.birth.count.mcmc[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

        period.nq0.mcmc <-
            list(female = ffrowapply(projection.deaths.per.births[["female"]][i1:i2,-1,drop = FALSE] + # Red/Black
                 projection.births.minus.deaths[["female"]][i1:i2,-1,drop = FALSE] / # Red-Black
                 total.birth.count.mcmc[["female"]][i1:i2,-1,drop = FALSE] * # Black
                 death.count.cohort.delta.0[["female"]][i1:i2,-1,drop = FALSE] / # Shaded
                 projection.births.minus.deaths[["female"]][i1:i2,-(ncol(projection.births.minus.deaths[["female"]]))
                                                            ,drop = FALSE]
                 ,X = projection.deaths.per.births[["female"]][,-1,drop = FALSE]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(projection.deaths.per.births[["male"]][i1:i2,-1,drop = FALSE] + # Red/Black
                  projection.births.minus.deaths[["male"]][i1:i2,-1,drop = FALSE] / # Red-Black
                  total.birth.count.mcmc[["male"]][i1:i2,-1,drop = FALSE] * # Black
                  death.count.cohort.delta.0[["male"]][i1:i2,-1,drop = FALSE] / # Shaded
                  projection.births.minus.deaths[["male"]][i1:i2,-(ncol(projection.births.minus.deaths[["male"]]))
                                                           ,drop = FALSE]
                  ,X = projection.deaths.per.births[["male"]][,-1,drop = FALSE]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

    }

    if("IMR" %in% outputs) {

        ## This is the easy one (and the one which makes least sense as a
        ## measure).

        IMR.mcmc <-
            list(female = ffrowapply((total.birth.death.count.mcmc[["female"]][i1:i2,] +
                 death.count.cohort.delta.0[["female"]][i1:i2,]) /
                 total.birth.count.mcmc[["female"]][i1:i2,]
                 ,X = total.birth.death.count.mcmc[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply((total.birth.death.count.mcmc[["male"]][i1:i2,] +
                  death.count.cohort.delta.0[["male"]][i1:i2,]) /
                  total.birth.count.mcmc[["male"]][i1:i2,]
                  ,X = total.birth.death.count.mcmc[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )

    }


    ## -------** Rates

    ## Birth rates are an input; don't need to calculate those. Moreover, TFR is
    ## simply the sum of these, so don't need to calculate that either.


    ## -------*** Person Years

    ## Approximate person-years lived in the interval by same
    ## method as Preston et al. (2001), p. 122; i.e.,
    ## ({n}{N}{x}(t) + {n}{N}{x-n}(t){n}{S}{x}(t))/2.
    ## BUT {n}{N}{x-n}(t){n}{S}{x}(t) is {n}{N}{x}(t+n) in
    ## pop.data.mcmc.

    if(any(c("person.years", "mort.rate", "mig.rate") %in% outputs)) {

        message("Calculating person years...")

        pop.head <-
            list(female = ff(pop.count.mcmc[["female"]][, pop.count.years %in% VR.years]
                 ,dim = c(nrow(pop.count.mcmc[["female"]]), sum(pop.count.years %in% VR.years))
                 )
                 ,male = ff(pop.count.mcmc[["male"]][, pop.count.years %in% VR.years]
                  ,dim = c(nrow(pop.count.mcmc[["male"]]), sum(pop.count.years %in% VR.years))
                  )
                 )
        pop.tail <-
            list(female = ff(pop.count.mcmc[["female"]][, pop.count.years %in%
                 as.character(as.numeric(VR.years) + age.int)]
                 ,dim = c(nrow(pop.count.mcmc[["female"]])
                  ,sum(pop.count.years %in% as.character(as.numeric(VR.years) + age.int))
                  )
                 )
                 ,male = ff(pop.count.mcmc[["male"]][, pop.count.years %in%
                  as.character(as.numeric(VR.years) + age.int)]
                  ,dim = c(nrow(pop.count.mcmc[["male"]])
                   ,sum(pop.count.years %in% as.character(as.numeric(VR.years) + age.int))
                   )
                  )
                 )
        person.years.mcmc <-
            list(female = ffrowapply((5/2) * (pop.head[["female"]][i1:i2,] + pop.tail[["female"]][i1:i2,])
                 ,X = pop.head[["female"]]
                 ,RETURN = TRUE)
                 ,male = ffrowapply((5/2) * (pop.head[["male"]][i1:i2,] + pop.tail[["male"]][i1:i2,])
                  ,X = pop.head[["male"]]
                  ,RETURN = TRUE)
                 )
        for(s in c("female", "male")) colnames(person.years.mcmc[[s]]) <-
            colnames(pop.count.mcmc[[s]])[pop.count.years %in% VR.years]

    }


    ## Mortality rates

    if("mort.rate" %in% outputs) {

        message("Calculating mortality rates...")

        mort.rate.mcmc <-
            list(female = ffrowapply(death.count.mcmc[["female"]][i1:i2,] /
                 person.years.mcmc[["female"]][i1:i2,]
                 ,X = death.count.mcmc[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(death.count.mcmc[["male"]][i1:i2,] /
                  person.years.mcmc[["male"]][i1:i2,]
                  ,X = death.count.mcmc[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )
    }


    ## Migration rates

    if("mig.rate" %in% outputs) {

        message("Calculating migration rates...")

        mig.rate.mcmc <-
            list(female = ffrowapply(mig.count.mcmc[["female"]][i1:i2,] /
                 person.years.mcmc[["female"]][i1:i2,]
                 ,X = mig.count.mcmc[["female"]]
                 ,RETURN = TRUE
                 ,USE.NAMES = TRUE
                 )
                 ,male = ffrowapply(mig.count.mcmc[["male"]][i1:i2,] /
                  person.years.mcmc[["male"]][i1:i2,]
                  ,X = mig.count.mcmc[["male"]]
                  ,RETURN = TRUE
                  ,USE.NAMES = TRUE
                  )
                 )
    }


    ## -------* OUTPUT

    if(file.backed) {

        for(y in outputs) {
            assign(paste(name.pref, y, ".mcmc", name.suf, sep = "")
                  ,get(paste(y, ".mcmc", sep = ""))
                  ,envir = .GlobalEnv
                   )
            message(paste(name.pref, y, ".mcmc", name.suf, sep = "")
                  , " placed in .GlobalEnv")
        }
    } else {
        out <- list()
        z <- get(paste(outputs[1], ".mcmc", sep = ""))
        if(is.recursive(z)) z <- list(female = z$female[], male = z$male[])
        else z <- list(z[])
        out <- z
        if(length(outputs) > 1) {
            out <- list(out)
            for(y in outputs[-1]) {
                z <- get(paste(y, ".mcmc", sep = ""))
                if(is.recursive(z)) z <- list(female = z$female[], male = z$male[])
                else z <- list(z[])
                out <- c(out, list(z))
            }
        }
        names(out) <- outputs
        return(out)
    }

}
