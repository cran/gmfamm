#' Draft of family for traffic example
#'
#' Fix four distributional assumptions and supply derivatives. Use BCCGo for
#' speed data. Use negative binomial for count data.
#' @param ... Not used.
#' @returns A bamlss family object.
trafficfam <- function(...) {

  BCCGo <- utils::getFromNamespace("BCCGo", "gamlss.dist")

  links <- list("mu1" = make.link("log"),
                "mu2" = make.link("log"),
                "mu3" = make.link("log"),
                "mu4" = make.link("log"))

  f <- list(
    "names" = c("mu1", "sigma1", "mu2", "sigma2", "mu3", "sigma3", "nu3",
                "mu4", "sigma4", "nu4", "Lambda"),
    "links" = c("mu1" = "log", "sigma1" = "log",
                "mu2" = "log", "sigma2" = "log",
                "mu3" = "log", "sigma3" = "log", "nu3" = "identity",
                "mu4" = "log", "sigma4" = "log", "nu4" = "identity",
                "Lambda" = "identity"),
    "d" = function(y, par, log = FALSE, ...) {

      i <- y[, 2] == 1
      mu1 <- links$mu1$linkfun(par$mu1[i])
      mu1 <- links$mu1$linkinv(mu1 + par$Lambda[i])
      sigma1 <- par$sigma1[i]
      y1 <- y[i, 1]

      i <- y[, 2] == 2
      mu2 <- links$mu2$linkfun(par$mu2[i])
      mu2 <- links$mu2$linkinv(mu2 + par$Lambda[i])
      sigma2 <- par$sigma2[i]
      y2 <- y[i, 1]

      i <- y[, 2] == 3
      mu3 <- links$mu3$linkfun(par$mu3[i])
      mu3 <- links$mu3$linkinv(mu3 + par$Lambda[i])
      sigma3 <- par$sigma3[i]
      nu3 <- par$sigma[i]
      y3 <- y[i, 1]

      i <- y[, 2] == 4
      mu4 <- links$mu4$linkfun(par$mu4[i])
      mu4 <- links$mu3$linkinv(mu4 + par$Lambda[i])
      sigma4 <- par$sigma4[i]
      nu4 <- par$sigma[i]
      y4 <- y[i, 1]

      d <- c(
        bamlss.family("nbinom")$d(y1, par = list("mu" = mu1,
                                                 "theta" = sigma1),
                                  log = log),
        bamlss.family("nbinom")$d(y2, par = list("mu" = mu2, "theta" = sigma2),
                                  log = log),
        bamlss.family(BCCGo)$d(y3, par = list("mu" = mu3, "sigma" = sigma3,
                                              "nu" = nu3)),
        bamlss.family(BCCGo)$d(y4, par = list("mu" = mu4, "sigma" = sigma4,
                                              "nu" = nu4))
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("nbinom")$score$mu(y1, par = par1)
        y_score

      },
      sigma1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("nbinom")$score$theta(y1, par = par1)
        y_score

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("nbinom")$score$mu(y2, par = par2)
        y_score

      },
      sigma2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("nbinom")$score$theta(y2, par = par2)
        y_score

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$mu(y3, par = par3)
        y_score

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$sigma(y3, par = par3)
        y_score

      },
      nu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$nu(y3, par = par3)
        y_score

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$mu(y4, par = par4)
        y_score

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$sigma(y4, par = par4)
        y_score

      },
      nu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$nu(y4, par = par4)
        y_score

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(
            y[i, 2],
            bamlss.family("nbinom")$score$mu(
              y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
            bamlss.family("nbinom")$score$mu(
              y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
            bamlss.family(BCCGo)$score$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i],
                                  "nu" = par$nu3[i])),
            bamlss.family(BCCGo)$score$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i],
                                  "nu" = par$nu4[i])))
        }
        y_score
      }

    ),
    "hess" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_hess[i] <- bamlss.family("nbinom")$hess$mu(y1, par = par1)
        y_hess

      },
      sigma1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_hess[i] <- bamlss.family("nbinom")$hess$theta(y1, par = par1)
        y_hess

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_hess[i] <- bamlss.family("nbinom")$hess$mu(y2, par = par2)
        y_hess

      },
      sigma2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_hess[i] <- bamlss.family("nbinom")$hess$theta(y2, par = par2)
        y_hess

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$mu(y3, par = par3)
        y_hess

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$sigma(y3, par = par3)
        y_hess

      },
      nu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$nu(y3, par = par3)
        y_hess

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$mu(y4, par = par4)
        y_hess

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$sigma(y4, par = par4)
        y_hess

      },
      nu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$nu(y4, par = par4)
        y_hess

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_hess <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_hess[i] <- switch(
            y[i, 2],
            bamlss.family("nbinom")$hess$mu(
              y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
            bamlss.family("nbinom")$hess$mu(
              y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
            bamlss.family(BCCGo)$hess$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i],
                                  "nu" = par$nu3[i])),
            bamlss.family(BCCGo)$hess$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i],
                                  "nu" = par$nu4[i])))
        }
        y_hess
      }

    )
  )
  class(f) <- "family.bamlss"
  return(f)
}

#' Draft 2 of family for traffic example
#'
#' Fix four distributional assumptions and supply derivatives. Use BCCGo for
#' speed data. Use zero-truncated negative binomial for count data (no second
#' derivatives available).
#' @param ... Not used.
#' @returns A bamlss family object.
trafficfam2 <- function(...) {

  BCCGo <- utils::getFromNamespace("BCCGo", "gamlss.dist")

  links <- list("mu1" = make.link("log"),
                "mu2" = make.link("log"),
                "mu3" = make.link("log"),
                "mu4" = make.link("log"))

  f <- list(
    "names" = c("mu1", "sigma1", "mu2", "sigma2", "mu3", "sigma3", "nu3",
                "mu4", "sigma4", "nu4", "Lambda"),
    "links" = c("mu1" = "log", "sigma1" = "log",
                "mu2" = "log", "sigma2" = "log",
                "mu3" = "log", "sigma3" = "log", "nu3" = "identity",
                "mu4" = "log", "sigma4" = "log", "nu4" = "identity",
                "Lambda" = "identity"),
    "d" = function(y, par, log = FALSE, ...) {

      i <- y[, 2] == 1
      mu1 <- links$mu1$linkfun(par$mu1[i])
      mu1 <- links$mu1$linkinv(mu1 + par$Lambda[i])
      sigma1 <- par$sigma1[i]
      y1 <- y[i, 1]

      i <- y[, 2] == 2
      mu2 <- links$mu2$linkfun(par$mu2[i])
      mu2 <- links$mu2$linkinv(mu2 + par$Lambda[i])
      sigma2 <- par$sigma2[i]
      y2 <- y[i, 1]

      i <- y[, 2] == 3
      mu3 <- links$mu3$linkfun(par$mu3[i])
      mu3 <- links$mu3$linkinv(mu3 + par$Lambda[i])
      sigma3 <- par$sigma3[i]
      nu3 <- par$sigma[i]
      y3 <- y[i, 1]

      i <- y[, 2] == 4
      mu4 <- links$mu4$linkfun(par$mu4[i])
      mu4 <- links$mu3$linkinv(mu4 + par$Lambda[i])
      sigma4 <- par$sigma4[i]
      nu4 <- par$sigma[i]
      y4 <- y[i, 1]

      d <- c(
        bamlss.family("ztnbinom")$d(y1, par = list("mu" = mu1,
                                                   "theta" = sigma1),
                                    log = log),
        bamlss.family("ztnbinom")$d(y2, par = list("mu" = mu2,
                                                   "theta" = sigma2),
                                    log = log),
        bamlss.family(BCCGo)$d(y3, par = list("mu" = mu3, "sigma" = sigma3,
                                              "nu" = nu3)),
        bamlss.family(BCCGo)$d(y4, par = list("mu" = mu4, "sigma" = sigma4,
                                              "nu" = nu4))
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("ztnbinom")$score$mu(y1, par = par1)
        y_score

      },
      sigma1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("ztnbinom")$score$theta(y1, par = par1)
        y_score

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("ztnbinom")$score$mu(y2, par = par2)
        y_score

      },
      sigma2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("ztnbinom")$score$theta(y2, par = par2)
        y_score

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$mu(y3, par = par3)
        y_score

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$sigma(y3, par = par3)
        y_score

      },
      nu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_score[i] <- bamlss.family(BCCGo)$score$nu(y3, par = par3)
        y_score

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$mu(y4, par = par4)
        y_score

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$sigma(y4, par = par4)
        y_score

      },
      nu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_score[i] <- bamlss.family(BCCGo)$score$nu(y4, par = par4)
        y_score

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(
            y[i, 2],
            bamlss.family("ztnbinom")$score$mu(
              y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
            bamlss.family("ztnbinom")$score$mu(
              y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
            bamlss.family(BCCGo)$score$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i],
                                  "nu" = par$nu3[i])),
            bamlss.family(BCCGo)$score$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i],
                                  "nu" = par$nu4[i])))
        }
        y_score
      }

    ),
    "hess" = list(
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$mu(y3, par = par3)
        y_hess

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$sigma(y3, par = par3)
        y_hess

      },
      nu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        nu3 <- par$nu3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3, "nu" = nu3)

        y_hess[i] <- bamlss.family(BCCGo)$hess$nu(y3, par = par3)
        y_hess

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$mu(y4, par = par4)
        y_hess

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$sigma(y4, par = par4)
        y_hess

      },
      nu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        nu4 <- par$nu4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4, "nu" = nu4)

        y_hess[i] <- bamlss.family(BCCGo)$hess$nu(y4, par = par4)
        y_hess

      }

    )
  )
  class(f) <- "family.bamlss"
  return(f)
}

#' Draft of family for traffic example
#'
#' Fix four distributional assumptions and supply derivatives. Use gamma for
#' speed data. Use negative binomial for count data.
#' @param ... Not used.
#' @returns A bamlss family object.
#' @export
#' @examples
#' # Construct data
#' set.seed(123)
#' # Number of subjects
#' n <- 10
#'
#' # Number of observations
#' ni <- 3
#'
#' # Covariate vector
#' x <- rep(rnorm(n), each = ni)
#' t <- rep(c(0, 0.5, 1), times = n)
#'
#' # Additive predictor
#' eta_1 <- eta_2 <- eta_3 <- eta_4 <- t + 0.5*x
#'
#' # Outcomes
#' y1 <- rnbinom(n*ni, exp(eta_1), 0.3)
#' y2 <- rnbinom(n*ni, exp(eta_2), 0.4)
#' y3 <- rgamma(n*ni, shape = 0.3, scale = exp(eta_3) / 0.3)
#' y4 <- rgamma(n*ni, shape = 0.4, scale = exp(eta_4) / 0.4)
#'
#' # Data format
#' dat <- data.frame(
#'    id = factor(rep(seq_len(n), each = ni)),
#'    y = c(y1, y2, y3, y4),
#'    dim = factor(rep(c(1:4), each = n*ni)),
#'    t = t,
#'    x = x,
#'    fpc = 1
#' )
#'
#' # Specify formula
#' f <- list(
#'   gm(y, dim) ~ t + x,
#'   sigma1 ~ 1,
#'   mu2 ~ t + x,
#'   sigma2 ~ 1,
#'   mu3 ~ t + x,
#'   sigma3 ~ 1,
#'   mu4 ~ t + x,
#'   sigma4 ~ 1,
#'   Lambda ~ -1 + s(id, by = fpc, bs = "re")
#' )
#'
#' # Model
#' b <- bamlss(f, family = trafficfam3, n.iter = 20, burnin = 10,
#'             data = dat)
trafficfam3 <- function(...) {

  links <- list("mu1" = make.link("log"),
                "mu2" = make.link("log"),
                "mu3" = make.link("log"),
                "mu4" = make.link("log"))

  f <- list(
    "names" = c("mu1", "sigma1", "mu2", "sigma2", "mu3", "sigma3",
                "mu4", "sigma4", "Lambda"),
    "links" = c("mu1" = "log", "sigma1" = "log",
                "mu2" = "log", "sigma2" = "log",
                "mu3" = "log", "sigma3" = "log",
                "mu4" = "log", "sigma4" = "log",
                "Lambda" = "identity"),
    "d" = function(y, par, log = FALSE, ...) {

      i <- y[, 2] == 1
      mu1 <- links$mu1$linkfun(par$mu1[i])
      mu1 <- links$mu1$linkinv(mu1 + par$Lambda[i])
      sigma1 <- par$sigma1[i]
      y1 <- y[i, 1]

      i <- y[, 2] == 2
      mu2 <- links$mu2$linkfun(par$mu2[i])
      mu2 <- links$mu2$linkinv(mu2 + par$Lambda[i])
      sigma2 <- par$sigma2[i]
      y2 <- y[i, 1]

      i <- y[, 2] == 3
      mu3 <- links$mu3$linkfun(par$mu3[i])
      mu3 <- links$mu3$linkinv(mu3 + par$Lambda[i])
      sigma3 <- par$sigma3[i]
      nu3 <- par$sigma[i]
      y3 <- y[i, 1]

      i <- y[, 2] == 4
      mu4 <- links$mu4$linkfun(par$mu4[i])
      mu4 <- links$mu3$linkinv(mu4 + par$Lambda[i])
      sigma4 <- par$sigma4[i]
      nu4 <- par$sigma[i]
      y4 <- y[i, 1]

      d <- c(
        bamlss.family("nbinom")$d(y1, par = list("mu" = mu1,
                                                 "theta" = sigma1),
                                  log = log),
        bamlss.family("nbinom")$d(y2, par = list("mu" = mu2, "theta" = sigma2),
                                  log = log),
        bamlss.family("gamma")$d(y3, par = list("mu" = mu3, "sigma" = sigma3),
                                 log = log),
        bamlss.family("gamma")$d(y4, par = list("mu" = mu4, "sigma" = sigma4),
                                 log = log)
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("nbinom")$score$mu(y1, par = par1)
        y_score

      },
      sigma1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_score[i] <- bamlss.family("nbinom")$score$theta(y1, par = par1)
        y_score

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("nbinom")$score$mu(y2, par = par2)
        y_score

      },
      sigma2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_score[i] <- bamlss.family("nbinom")$score$theta(y2, par = par2)
        y_score

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_score[i] <- bamlss.family("gamma")$score$mu(y3, par = par3)
        y_score

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_score[i] <- bamlss.family("gamma")$score$sigma(y3, par = par3)
        y_score

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_score[i] <- bamlss.family("gamma")$score$mu(y4, par = par4)
        y_score

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_score[i] <- bamlss.family("gamma")$score$sigma(y4, par = par4)
        y_score

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(y[i, 2],
             bamlss.family("nbinom")$score$mu(
               y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
             bamlss.family("nbinom")$score$mu(
               y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
             bamlss.family("gamma")$score$mu(
               y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i])),
             bamlss.family("gamma")$score$mu(
               y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i])))
        }
        y_score
      }

    ),
    "hess" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_hess[i] <- bamlss.family("nbinom")$hess$mu(y1, par = par1)
        y_hess

      },
      sigma1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        sigma1 <- par$sigma1[i]
        par1 <- list("mu" = mu1, "theta" = sigma1)

        y_hess[i] <- bamlss.family("nbinom")$hess$theta(y1, par = par1)
        y_hess

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_hess[i] <- bamlss.family("nbinom")$hess$mu(y2, par = par2)
        y_hess

      },
      sigma2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("mu" = mu2, "theta" = sigma2)

        y_hess[i] <- bamlss.family("nbinom")$hess$theta(y2, par = par2)
        y_hess

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_hess[i] <- bamlss.family("gamma")$hess$mu(y3, par = par3)
        y_hess

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_hess[i] <- bamlss.family("gamma")$hess$sigma(y3, par = par3)
        y_hess

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_hess[i] <- bamlss.family("gamma")$hess$mu(y4, par = par4)
        y_hess

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_hess[i] <- bamlss.family("gamma")$hess$sigma(y4, par = par4)
        y_hess

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_hess <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_hess[i] <- switch(y[i, 2],
            bamlss.family("nbinom")$hess$mu(
              y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
            bamlss.family("nbinom")$hess$mu(
              y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
            bamlss.family("gamma")$hess$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i])),
            bamlss.family("gamma")$hess$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i])))
        }
        y_hess
      }

    )
  )
  class(f) <- "family.bamlss"
  return(f)
}

#' Draft of family for traffic example
#'
#' Fix four distributional assumptions and supply derivatives. Use lognormal for
#' speed data. Use Poisson for count data.
#' @param ... Not used.
#' @returns A bamlss family object.
trafficfam4 <- function(...) {

  links <- list("mu1" = make.link("log"),
                "mu2" = make.link("log"),
                "mu3" = make.link("identity"),
                "mu4" = make.link("identity"))

  f <- list(
    "names" = c("mu1", "mu2", "mu3", "sigma3", "mu4", "sigma4", "Lambda"),
    "links" = c("mu1" = "log", "mu2" = "log",
                "mu3" = "identity", "sigma3" = "log",
                "mu4" = "identity", "sigma4" = "log",
                "Lambda" = "identity"),
    "d" = function(y, par, log = FALSE, ...) {

      i <- y[, 2] == 1
      mu1 <- links$mu1$linkfun(par$mu1[i])
      mu1 <- links$mu1$linkinv(mu1 + par$Lambda[i])
      y1 <- y[i, 1]

      i <- y[, 2] == 2
      mu2 <- links$mu2$linkfun(par$mu2[i])
      mu2 <- links$mu2$linkinv(mu2 + par$Lambda[i])
      y2 <- y[i, 1]

      i <- y[, 2] == 3
      mu3 <- links$mu3$linkfun(par$mu3[i])
      mu3 <- links$mu3$linkinv(mu3 + par$Lambda[i])
      sigma3 <- par$sigma3[i]
      nu3 <- par$sigma[i]
      y3 <- y[i, 1]

      i <- y[, 2] == 4
      mu4 <- links$mu4$linkfun(par$mu4[i])
      mu4 <- links$mu3$linkinv(mu4 + par$Lambda[i])
      sigma4 <- par$sigma4[i]
      nu4 <- par$sigma[i]
      y4 <- y[i, 1]

      d <- c(
        bamlss.family("poisson")$d(y1, par = list("lambda" = mu1), log = log),
        bamlss.family("poisson")$d(y2, par = list("lambda" = mu2), log = log),
        bamlss.family("lognormal")$d(y3, par = list("mu" = mu3, "sigma" = sigma3),
                                     log = log),
        bamlss.family("lognormal")$d(y4, par = list("mu" = mu4, "sigma" = sigma4),
                                     log = log)
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        par1 <- list("lambda" = mu1)

        y_score[i] <- bamlss.family("poisson")$score$lambda(y1, par = par1)
        y_score

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        sigma2 <- par$sigma2[i]
        par2 <- list("lambda" = mu2)

        y_score[i] <- bamlss.family("poisson")$score$lambda(y2, par = par2)
        y_score

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_score[i] <- bamlss.family("lognormal")$score$mu(y3, par = par3)
        y_score

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_score[i] <- bamlss.family("lognormal")$score$sigma(y3, par = par3)
        y_score

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_score[i] <- bamlss.family("lognormal")$score$mu(y4, par = par4)
        y_score

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_score <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_score[i] <- bamlss.family("lognormal")$score$sigma(y4, par = par4)
        y_score

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(
            y[i, 2],
            bamlss.family("poisson")$score$lambda(
              y[i, 1], par = list("lambda" = mu1[i])),
            bamlss.family("poisson")$score$lambda(
              y[i, 1], par = list("lambda" = mu2[i])),
            bamlss.family("lognormal")$score$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i])),
            bamlss.family("lognormal")$score$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i]))
          )
        }
        y_score
      }

    ),
    "hess" = list(
      mu1 = function(y, par, ...) {

        i <- y[, 2] == 1
        y1 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta1 <- links$mu1$linkfun(par$mu1[i]) + par$Lambda[i]
        mu1 <- links$mu1$linkinv(eta1)
        par1 <- list("lambda" = mu1)

        y_hess[i] <- bamlss.family("poisson")$hess$lambda(y1, par = par1)
        y_hess

      },
      mu2 = function(y, par, ...) {

        i <- y[, 2] == 2
        y2 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta2 <- links$mu1$linkfun(par$mu2[i]) + par$Lambda[i]
        mu2 <- links$mu2$linkinv(eta2)
        par2 <- list("lambda" = mu2)

        y_hess[i] <- bamlss.family("poisson")$hess$lambda(y2, par = par2)
        y_hess

      },
      mu3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_hess[i] <- bamlss.family("lognormal")$hess$mu(y3, par = par3)
        y_hess

      },
      sigma3 = function(y, par, ...) {

        i <- y[, 2] == 3
        y3 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta3 <- links$mu3$linkfun(par$mu3[i]) + par$Lambda[i]
        mu3 <- links$mu3$linkinv(eta3)
        sigma3 <- par$sigma3[i]
        par3 <- list("mu" = mu3, "sigma" = sigma3)

        y_hess[i] <- bamlss.family("lognormal")$hess$sigma(y3, par = par3)
        y_hess

      },
      mu4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu3$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu3$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_hess[i] <- bamlss.family("lognormal")$hess$mu(y4, par = par4)
        y_hess

      },
      sigma4 = function(y, par, ...) {

        i <- y[, 2] == 4
        y4 <- y[i, 1]
        y_hess <- rep(0, nrow(y))

        eta4 <- links$mu4$linkfun(par$mu4[i]) + par$Lambda[i]
        mu4 <- links$mu4$linkinv(eta4)
        sigma4 <- par$sigma4[i]
        par4 <- list("mu" = mu4, "sigma" = sigma4)

        y_hess[i] <- bamlss.family("lognormal")$hess$sigma(y4, par = par4)
        y_hess

      },
      Lambda = function(y, par, ...) {

        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        eta3 <- links$mu3$linkfun(par$mu3) + par$Lambda
        mu3 <- links$mu3$linkinv(eta3)
        eta4 <- links$mu4$linkfun(par$mu4) + par$Lambda
        mu4 <- links$mu4$linkinv(eta4)

        y_hess <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_hess[i] <- switch(
            y[i, 2],
            bamlss.family("poisson")$hess$lambda(
              y[i, 1], par = list("lambda" = mu1[i])),
            bamlss.family("poisson")$hess$lambda(
              y[i, 1], par = list("lambda" = mu2[i])),
            bamlss.family("lognormal")$hess$mu(
              y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i])),
            bamlss.family("lognormal")$hess$mu(
              y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i]))
          )
        }
        y_hess
      }

    )
  )
  class(f) <- "family.bamlss"
  return(f)
}
