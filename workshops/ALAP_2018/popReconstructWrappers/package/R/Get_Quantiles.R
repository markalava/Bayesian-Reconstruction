################################################################################
###-----------------------------------------------------------------------------
###
### SYNOPSIS:
###
### Plot results of population reconstruction
###
###-----------------------------------------------------------------------------
################################################################################

##' Get prior and posterior quantiles for results of Bayesian population reconstruction.
##'
##' Convert outputs of \code{\link{pop.est.sampler}} and
##' \code{\link{sample.from.prior}} and produce marginal quantiles fo
##' various parameters, returned in a conventient form for plotting.
##'
##' @param param Character vector of parameters for which quantiles are wanted.
##' @param results.recon Result of running \code{\link{pop.est.sampler}}.
##' @param results.post.process.recon Result of running \code{\link{post.process.recon}} on \code{results.recon}.
##' @param results.prior Result of running \code{\link{sample.from.prior}}.
##' @param results.post.process.prior Result of running \code{\link{post.process.recon}} on \code{results.prior}.
##' @return Dataframe with prior and posterior marginal quantiles for requested parameters.
##' @author Mark C. Wheldon
##' @export
get.quantiles.recon <- function(param = c("srb", "fert.rate", "mort.rate", "mig.count",
                                          "baseline.count", "tfr", "e0", "IMR"),
                                results.recon,
                                results.post.process.recon,
                                results.prior = NULL,
                                results.post.process.prior = NULL
                                ) {

    ## -------* Set-up

    param.df <-
        ##      name                       by.age   by.sex   by.year  post.proc.param needs.prior.sample  needs.transform
        rbind(c("srb",                     FALSE,   FALSE,   TRUE,    TRUE,           FALSE,              FALSE),
              c("mort.rate",               TRUE,    TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("surv.prop",               TRUE,    TRUE,    TRUE,    FALSE,          TRUE,               FALSE),
              c("mig.prop",                TRUE,    TRUE,    TRUE,    FALSE,          TRUE,               FALSE),
              c("mig.rate",                TRUE,    TRUE,    TRUE,    TRUE ,          TRUE,               FALSE),
              c("birth.count",             TRUE,    TRUE,    TRUE,    TRUE ,          TRUE,               FALSE),
              c("death.count",             TRUE,    TRUE,    TRUE,    TRUE ,          TRUE,               FALSE),
              c("mig.count",               TRUE,    TRUE,    TRUE,    TRUE ,          TRUE,               TRUE),
              c("baseline.count",          TRUE,    TRUE,    FALSE,   FALSE,          FALSE,              FALSE),
              c("total.birth.count",       FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("total.birth.death.count", FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("total.mig.count",         FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("cohort.nq0",              FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("period.nq0",              FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("IMR",                     FALSE,   TRUE,    TRUE,    TRUE,           TRUE,               FALSE),
              c("fert.rate",               TRUE,    FALSE,   TRUE,    FALSE,          TRUE,               FALSE),
              c("tfr",                     FALSE,   FALSE,   TRUE,    FALSE,          FALSE,              TRUE),
              c("e0",                      FALSE,   TRUE,    TRUE,    FALSE,          TRUE,               TRUE)
              )
    param.df <- as.data.frame(param.df, stringsAsFactors = FALSE)
    for(j in (1:ncol(param.df))[-1]) param.df[,j] <- as.logical(param.df[,j])
    colnames(param.df) <- c("param", "by.age", "by.sex", "by.year", "post.processed.param",
                            "needs.prior.sample", "needs.transform")

    y <- character(0)
    for(i in 1:nrow(param.df)) { y <- c(y, paste(as.numeric(param.df[i,-1]), collapse = ""))}
    y <- factor(y)
    levels(y) <- seq_along(levels(y))
    param.df$class <- factor(y)

    param.df <- param.df[order(param.df$class),]

    allowed.param <- param.df$param

    if(!all(param %in% allowed.param)) {
        stop("'param' must be one of '", paste(allowed.param, collapse = "', '"), "'.")
    }

    param.req.df <- param.df[param.df$param %in% param,]

    ## For now, make prior sample and post processed results a requirement.

    ## These are hard-coded for now
    quants.to.plot = c(0.025, 0.1, 0.5, 0.9, 0.975)

    ## -------** Functions

    ## Default for quantile() is to remove NA
    quantile <- function(x, probs = seq(0, 1, 0.25), na.rm = TRUE,
                         names = TRUE, type = 7, ...) {
        stats::quantile(x, probs = probs, na.rm = na.rm,
                        names = names, type = type, ...)
    }

    leb.f <- function(z) {
        ## z is a vector of age-specific survival proportions
        x <- c(head(z, -1), tail(z,1) / (1-tail(z,1)))
        5 * sum(cumprod(x))
    }

    ## -------* Quantiles

    out.df <- data.frame()

    for(this.param in param) {

        if(this.param == "srb") {

            m <- results.recon$srb.mcmc

            ## Calculate posterior quantiles
            q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
            dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

            ## Make ages and years
            years <- colnames(m)

            ## Prepare data sets
            alpha <- results.recon$fixed.params$alpha.srb
            beta <- results.recon$fixed.params$beta.srb
            mqvit.df <- t(q.vital)

            colnames(mqvit.df) <-
                paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                    , "pctl", sep = "")

            mqvit.df <-
                data.frame(mqvit.df
                          ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                        ,"[[", 1)
                          ,legend = "posterior"
                           )

            if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                meas <- try(reshape::melt(results.recon$fixed.params$mean.srb)[,-1],
                            silent = TRUE)
                if(class(meas) == "try-error") {
                    meas <- try(reshape::melt.array(results.recon$fixed.params$mean.srb)[,-1],
                                silent = TRUE)
                }
                colnames(meas) <- c("year", "param.50pctl")

                upperQ95 <- exp(log(meas$param.50pctl) + qt(p = 1-0.975, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                lowerQ95 <- exp(log(meas$param.50pctl) - qt(p = 1-0.975, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                upperQ90 <- exp(log(meas$param.50pctl) + qt(p = 1-0.9, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                lowerQ90 <- exp(log(meas$param.50pctl) - qt(p = 1-0.9, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )

                mmeas.df <-
                    data.frame(meas, param.97.5pctl = reshape::melt(upperQ95)$value
                              ,param.90pctl = reshape::melt(upperQ90)$value
                              ,param.10pctl = reshape::melt(lowerQ90)$value
                              ,param.2.5pctl = reshape::melt(lowerQ95)$value
                              ,legend = "init. est."
                               )
            } else {
                mmeas.df <- data.frame()
            }

            ##
            ## Combine
            ##

            out.df <- rbind(out.df,
                            data.frame(rbind(mmeas.df,
                                             mqvit.df
                                             ),
                                       param = this.param, sex = NA, age = NA)
                            )

        } else if(this.param == "fert.rate") {

            m <- results.recon$fert.rate.mcmc
            fert.rows <- results.recon$alg.params$non.zero.fert.rows

            ## Calculate posterior quantiles
            q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
            dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

            ## Make ages and years
            colspl <- strsplit(colnames(m), ".", fixed = TRUE)
            years <- unique(sapply(colspl, FUN = function(z) z[1]))
            fert.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
            fert.ages.numeric <- as.numeric(gsub("[^0-9]", "", fert.ages))

            ## Prepare data sets
            alpha <- results.recon$fixed.params$alpha.fert.rate
            beta <- results.recon$fixed.params$beta.fert.rate

            mqvit.df <- t(q.vital)
            colnames(mqvit.df) <-
                paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                    , "pctl", sep = "")
            mqvit.df <-
                data.frame(mqvit.df
                          ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                       ,"[[", 2)
                          ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                        ,"[[", 1)
                          ,legend = "posterior"
                           )

            if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                meas <- reshape::melt(results.recon$fixed.params$mean.fert.rate[fert.rows,])
                meas <-
                    gdata::rename.vars(meas, from = colnames(meas)
                                      ,to = c("age", "year", "param.50pctl"), info = FALSE)

                upperQ95 <- exp(log(meas$param.50pctl) + qt(p = 1-0.975, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                lowerQ95 <- exp(log(meas$param.50pctl) - qt(p = 1-0.975, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                upperQ90 <- exp(log(meas$param.50pctl) + qt(p = 1-0.9, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                lowerQ90 <- exp(log(meas$param.50pctl) - qt(p = 1-0.9, df = 2 * alpha
                                                          , lower.tail = FALSE) *
                                sqrt(beta/alpha)
                                )
                mmeas.df <-
                    data.frame(meas, param.97.5pctl = reshape::melt(upperQ95)$value
                              ,param.90pctl = reshape::melt(upperQ90)$value
                              ,param.10pctl = reshape::melt(lowerQ90)$value
                              ,param.2.5pctl = reshape::melt(lowerQ95)$value
                              ,legend = "init. est."
                               )
            } else {
                mmeas.df <- data.frame()
            }

            ##
            ## Combine
            ##

            out.df <- rbind(out.df,
                            data.frame(rbind(mmeas.df, mqvit.df),
                                       param = this.param, sex = "female")
                            )

        } else if(this.param == "baseline.count") {

            for(this.sex in c("female", "male")) {

                m <- results.recon$baseline.count.mcmc[[this.sex]]

                ## Calculate posterior quantiles
                q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
                dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                ## Make ages and years
                colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                year <- unique(sapply(colspl, FUN = function(z) z[1]))
                baseline.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
                baseline.ages.numeric <- as.numeric(gsub("[^0-9]", "", baseline.ages))

                ## Prepare data sets
                alpha <- results.recon$fixed.params$alpha.population.count
                beta <- results.recon$fixed.params$beta.population.count

                mqvit.df <- t(q.vital)
                colnames(mqvit.df) <-
                    paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                        , "pctl", sep = "")
                mqvit.df <-
                    data.frame(mqvit.df
                              ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                           ,"[[", 2)
                              ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                            ,"[[", 1)
                              ,legend = "posterior"
                              ,sex = this.sex
                               )

                if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                    meas <- reshape::melt(results.recon$fixed.params$mean.baseline.count[[this.sex]])
                    meas <-
                        gdata::rename.vars(meas, from = colnames(meas)
                                          ,to = c("age", "year", "param.50pctl"), info = FALSE)
                    meas$age <-
                        sapply(strsplit(as.character(meas$age), split = "[^0-9]"), "[[", 1)

                    upperQ95 <- exp(log(meas$param.50pctl) + qt(p = 1-0.975
                                                              , df = 2 * alpha
                                                              , lower.tail = FALSE) *
                                    sqrt(beta/alpha)
                                    )
                    lowerQ95 <- exp(log(meas$param.50pctl) - qt(p = 1-0.975
                                                              , df = 2 * alpha
                                                              , lower.tail = FALSE) *
                                    sqrt(beta/alpha)
                                    )
                    upperQ90 <- exp(log(meas$param.50pctl) + qt(p = 1-0.9
                                                              , df = 2 * alpha
                                                              , lower.tail = FALSE) *
                                    sqrt(beta/alpha)
                                    )
                    lowerQ90 <- exp(log(meas$param.50pctl) - qt(p = 1-0.9
                                                              , df = 2 * alpha
                                                              , lower.tail = FALSE) *
                                    sqrt(beta/alpha)
                                    )
                    mmeas.df <-
                        data.frame(meas, param.97.5pctl = reshape::melt(upperQ95)$value
                                  ,param.90pctl = reshape::melt(upperQ90)$value
                                  ,param.10pctl = reshape::melt(lowerQ90)$value
                                  ,param.2.5pctl = reshape::melt(lowerQ95)$value
                                  ,legend = "init. est."
                                  ,sex = this.sex
                                   )
                } else {
                    mmeas.df <- data.frame()
                }


                ##
                ## Combine
                ##

                out.df <- rbind(out.df,
                                data.frame(rbind(mmeas.df, mqvit.df),
                                           param = this.param)
                                )
            }

        } else if(this.param %in% param.df$param[with(param.df, by.age & by.sex & by.year &!post.processed.param)]) {

            for(this.sex in c("female", "male")) {

                m <- results.recon[[paste0(this.param, ".mcmc")]][[this.sex]][,]

                ## Calculate posterior quantiles
                q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
                dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                if(this.param == "mig.count") q.vital <- q.vital/5

                ## Make ages and years
                colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                year <- unique(sapply(colspl, FUN = function(z) z[1]))
                mig.ages <- unique(sapply(colspl, FUN = function(z) z[2])) #CHANGE NAME!
                mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

                ## Prepare data sets
                mqvit.df <- t(q.vital)
                colnames(mqvit.df) <-
                    paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                        , "pctl", sep = "")
                mqvit.df <-
                    data.frame(mqvit.df
                              ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                           ,"[[", 2)
                              ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                            ,"[[", 1)
                              ,legend = "posterior", sex = this.sex
                               )
                mqvit.df$age <- as.numeric(levels(mqvit.df$age)[mqvit.df$age])

                if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                    ## Priors
                    ## Using sample from the prior

                    m <- results.prior[[paste0(this.param, ".prior")]][[this.sex]][,]

                    q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
                    dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                    if(this.param == "mig.count") q.vital <- q.vital/5

                    ## Make ages and years
                    colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                    year <- unique(sapply(colspl, FUN = function(z) z[1]))
                    mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
                    mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

                        ## Prior quantiles
                        mmeas.df <- t(q.vital)
                        colnames(mmeas.df) <-
                            paste("param.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
                                , "pctl", sep = "")
                        mmeas.df <-
                            data.frame(mmeas.df
                                      ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                                                   ,"[[", 2)
                                      ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                                                    ,"[[", 1)
                                      ,legend = "init. est.", sex = this.sex
                                       )
                        mmeas.df$age <- as.numeric(levels(mmeas.df$age)[mmeas.df$age])

                    } else {
                        mmeas.df <- data.frame()
                    }

                    ##
                    ## Combine
                    ##

                    out.df <- rbind(out.df,
                                    data.frame(rbind(mmeas.df, mqvit.df),
                                               param = this.param)
                                    )
                }

            } else if(this.param %in% param.df$param[with(param.df, by.age & by.sex & by.year & post.processed.param & needs.prior.sample & !needs.transform)]) {

                for(this.sex in c("female", "male")) {

                    m <- results.post.process.recon[[this.param]][[this.sex]][,]

                    ## Calculate posterior quantiles
                    q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
                    dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                    ## Make ages and years
                    colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                    year <- unique(sapply(colspl, FUN = function(z) z[1]))
                    mig.ages <- unique(sapply(colspl, FUN = function(z) z[2])) #CHANGE NAME!
                    mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

                    ## Prepare data sets
                    mqvit.df <- t(q.vital)
                    colnames(mqvit.df) <-
                        paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                            , "pctl", sep = "")
                    mqvit.df <-
                        data.frame(mqvit.df
                                  ,age = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                               ,"[[", 2)
                                  ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                                ,"[[", 1)
                                  ,legend = "posterior", sex = this.sex
                                   )
                    mqvit.df$age <- as.numeric(levels(mqvit.df$age)[mqvit.df$age])

                    if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                        ## Priors
                        ## Using sample from the prior

                        m <- results.post.process.prior[[this.param]][[this.sex]][,]

                        q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
                        dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                        ## Make ages and years
                        colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                        year <- unique(sapply(colspl, FUN = function(z) z[1]))
                        mig.ages <- unique(sapply(colspl, FUN = function(z) z[2]))
                        mig.ages.numeric <- as.numeric(gsub("[^0-9]", "", mig.ages))

                        ## Posterior quantiles
                        mmeas.df <- t(q.vital)
                        colnames(mmeas.df) <-
                            paste("param.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
                                , "pctl", sep = "")
                        mmeas.df <-
                            data.frame(mmeas.df
                                      ,age = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                                                   ,"[[", 2)
                                      ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                                                    ,"[[", 1)
                                      ,legend = "init. est.", sex = this.sex
                                       )
                        mmeas.df$age <- as.numeric(levels(mmeas.df$age)[mmeas.df$age])

                    } else {
                        mmeas.df <- data.frame()
                    }

                    ##
                    ## Combine
                    ##

                    out.df <- rbind(out.df,
                                    data.frame(rbind(mmeas.df, mqvit.df),
                                               param = this.param)
                                    )
                }
            } else if(this.param %in% param.df$param[with(param.df, by.sex & !by.age & by.year & post.processed.param & needs.prior.sample & !needs.transform)]) {

                for(this.sex in c("female", "male")) {

                    ##
                    ## Female
                    ##

                    ## Posterior

                    ## Calculate quantiles
                    m <- results.post.process.recon[[this.param]][[this.sex]][,]

                    q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot))
                    dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                    ## Make ages and years
                    colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                    year <- unique(sapply(colspl, FUN = function(z) z[1]))

                    ## Posterior quantiles
                    mqvit.df <- t(q.vital)
                    colnames(mqvit.df) <-
                        paste("param.", prettyNum(as.numeric(colnames(mqvit.df)) * 100)
                            , "pctl", sep = "")
                    mqvit.df <-
                        data.frame(mqvit.df
                                  ,year = sapply(strsplit(rownames(mqvit.df), split = "[^0-9]")
                                                ,"[[", 1)
                                  ,legend = "posterior", sex = this.sex, age = NA
                                   )

                    if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                        ## Prior
                        m <- results.post.process.prior[[this.param]][[this.sex]][,]

                        q.vital <- apply(m, 2, function(z) quantile(z, probs = quants.to.plot, na.rm = TRUE))
                        dimnames(q.vital) <- list(as.character(quants.to.plot), colnames(m))

                        ## Make ages and years
                        colspl <- strsplit(colnames(m), ".", fixed = TRUE)
                        year <- unique(sapply(colspl, FUN = function(z) z[1]))

                        ## Prior quantiles
                        mmeas.df <- t(q.vital)
                        colnames(mmeas.df) <-
                            paste("param.", prettyNum(as.numeric(colnames(mmeas.df)) * 100)
                                , "pctl", sep = "")
                        mmeas.df <-
                            data.frame(mmeas.df
                                      ,year = sapply(strsplit(rownames(mmeas.df), split = "[^0-9]")
                                                    ,"[[", 1)
                                      ,legend = "init. est.", sex = this.sex, age = NA
                                       )

                    } else {
                        mmeas.df <- data.frame()
                    }

                    ##
                    ## Combine
                    ##

                    out.df <- rbind(out.df,
                                    data.frame(rbind(mmeas.df, mqvit.df),
                                               param = this.param)
                                    )

                }

            } else if(this.param == "e0") {

                for(this.sex in c("female", "male")) {

                    surv.prop.years <-
                        sapply(strsplit(colnames(results.recon$surv.prop.mcmc[[this.sex]]), "\\."), "[[", 1)

                    leb.stationary.df <-
                        apply(results.recon$surv.prop.mcmc[[this.sex]][,], 1, function(z) {
                            tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
                        })

                    leb.stationary.Quantiles <-
                        apply(leb.stationary.df, 1, "quantile", probs = quants.to.plot)

                    leb.stationary.Quantiles.df <-
                        as.data.frame(t(leb.stationary.Quantiles))

                    colnames(leb.stationary.Quantiles.df) <-
                        paste("param."
                            , strsplit(colnames(leb.stationary.Quantiles.df)
                                     , split = "%")
                             ,"pctl", sep = "")

                    leb.stationary.Quantiles.df$legend <- "posterior"
                    leb.stationary.Quantiles.df$year <-
                        as.numeric(rownames(leb.stationary.Quantiles.df))

                    ## For female vs male
                    leb.stationary.Quantiles.df <-
                        leb.stationary.Quantiles.df

                    if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                        ## Prior by converting posterior quantiles of survival and assuming
                        ## stationary population relation holds

                        lebPrior.stationary.df <-
                            apply(results.prior$surv.prop.prior[[this.sex]][,], 1, function(z) {
                                tapply(z, INDEX = surv.prop.years, FUN = "leb.f")
                            })

                        lebPrior.stationary.Quantiles <-
                            apply(lebPrior.stationary.df, 1, "quantile", probs = quants.to.plot)

                        lebPrior.stationary.Quantiles.df <-
                            as.data.frame(t(lebPrior.stationary.Quantiles))

                        colnames(lebPrior.stationary.Quantiles.df) <-
                            paste("param."
                                , strsplit(colnames(lebPrior.stationary.Quantiles.df)
                                         , split = "%")
                                 ,"pctl", sep = "")

                        lebPrior.stationary.Quantiles.df$legend <- "init. est."
                        lebPrior.stationary.Quantiles.df$year <-
                            as.numeric(rownames(lebPrior.stationary.Quantiles.df))

                    } else {
                        lebPrior.stationary.Quantiles.df <- data.frame()
                    }

                    ##
                    ## Combine
                    ##

                    out.df <- rbind(out.df,
                                    data.frame(rbind(leb.stationary.Quantiles.df,
                                                     lebPrior.stationary.Quantiles.df
                                                     ),
                                               param = this.param, age = NA, sex = this.sex)
                                    )
                }

            } else if(this.param == "tfr") {

                ##
                ## Posterior
                ##

                dn <- list(NULL,
                           unique(sapply(strsplit(colnames(results.recon$fert.rate.mcmc)
                                                 ,"\\."), FUN = function(z) z[[1]])
                                  )
                           )
                ## Calculate
                ctry.tfr <-
                    matrix(0, nrow = nrow(results.recon$fert.rate.mcmc)
                          ,ncol = length(dn[[2]])
                          ,dimnames = dn
                           )

                fert.rate.mcmc.colYrs <-
                    sapply(strsplit(colnames(results.recon$fert.rate.mcmc)
                                   ,"\\."), FUN = function(z) z[[1]])

                for(i in 1:ncol(ctry.tfr)) {
                    colYrs.index <- fert.rate.mcmc.colYrs == colnames(ctry.tfr)[i]
                    ctry.tfr[,i] <-
                        apply(results.recon$fert.rate.mcmc[,colYrs.index]
                             ,1
                             ,FUN = function(z) sum(z)
                              )
                }

                ## Posterior quantiles
                ctry.tfrQuant <- apply(ctry.tfr, 2, FUN = function(z)
                {
                    5 * quantile(z, probs = quants.to.plot)
                })

                ctry.tfrQuant.df <-
                    as.data.frame(t(ctry.tfrQuant))

                colnames(ctry.tfrQuant.df) <-
                    paste("param.", strsplit(colnames(ctry.tfrQuant.df), split = "%")
                         ,"pctl", sep = "")

                ctry.tfrQuant.df$legend = "posterior"
                ctry.tfrQuant.df$year = as.numeric(rownames(ctry.tfrQuant.df))

                if(!is.null(results.post.process.prior) && !is.null(results.prior)) {

                    ##
                    ## Prior
                    ##

                    ## Calculate
                    ctry.tfr.prior <-
                        matrix(0, nrow = nrow(results.prior$fert.rate.prior)
                              ,ncol = length(dn[[2]])
                              ,dimnames = dn
                               )

                    fert.rate.mcmc.colYrs <-
                        sapply(strsplit(colnames(results.prior$fert.rate.prior)
                                       ,"\\."), FUN = function(z) z[[1]])

                    for(i in 1:ncol(ctry.tfr.prior)) {
                        colYrs.index <- fert.rate.mcmc.colYrs == colnames(ctry.tfr.prior)[i]
                        ctry.tfr.prior[,i] <-
                            apply(results.prior$fert.rate.prior[,colYrs.index]
                                 ,1
                                 ,FUN = function(z) sum(z)
                                  )
                    }

                    ## Init. Est. quantiles
                    ctry.tfr.priorQuant <- apply(ctry.tfr.prior, 2, FUN = function(z)
                    {
                        5 * quantile(z, probs = quants.to.plot)
                    })

                    ctry.tfr.priorQuant.df <-
                        as.data.frame(t(ctry.tfr.priorQuant))

                    colnames(ctry.tfr.priorQuant.df) <-
                        paste("param.", strsplit(colnames(ctry.tfr.priorQuant.df), split = "%")
                             ,"pctl", sep = "")

                    ctry.tfr.priorQuant.df$legend = "init. est."
                    ctry.tfr.priorQuant.df$year = as.numeric(rownames(ctry.tfr.priorQuant.df))

                } else {
                    ctry.tfr.priorQuant.df <- data.frame()
                }


                ##
                ## Combine
                ##

                out.df <- rbind(out.df,
                                data.frame(rbind(ctry.tfrQuant.df,
                                                 ctry.tfr.priorQuant.df
                                                 ),
                                           param = this.param, age = NA, sex = NA)
                                )
            }
        }

        ## Kludge any pesky factors back to numeric
        out.df$age <- as.numeric(levels(as.factor(out.df$age))[as.factor(out.df$age)])
        out.df$year <- as.numeric(levels(as.factor(out.df$year))[as.factor(out.df$year)])

        ## Wipe out rownames
        rownames(out.df) <- NULL

        return(out.df)
    }



