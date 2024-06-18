#' First draft of new family
#'
#' Fix three distributional assumptions and do not supply any derivatives.
#' @param ... Not used.
#' @returns A bamlss family object.
fam <- function(...) {

  links <- list("mu1" = make.link("logit"),
                "mu2" = make.link("log"),
                "mu3" = make.link("identity"))

  f <- list(
    "names" = c("mu1", "mu2", "mu3", "sigma3", "Lambda"),
    "links" = c("mu1" = "logit", "mu2" = "log", "mu3" = "identity",
                "sigma3" = "log",
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
      y3 <- y[i, 1]
      sigma3 <- par$sigma3[i]

      d <- c(
        dbinom(y1, size = 1, prob = mu1, log = log),
        dpois(y2, lambda = mu2, log = log),
        dnorm(y3, mean = mu3, sd = sigma3, log = log)
      )

      return(d)
    },
    "predict" = gmfamm_predict
  )
  class(f) <- "family.bamlss"
  return(f)
}

#' Next draft of new family
#'
#' Fix three distributional assumptions but supply derivatives.
#' @param ... Not used.
#' @returns A bamlss family object.
fam2 <- function(...) {

  links <- list("mu1" = make.link("logit"),
                "mu2" = make.link("log"),
                "mu3" = make.link("identity"))

  f <- list(
    "names" = c("mu1", "mu2", "mu3", "sigma3", "Lambda"),
    "links" = c("mu1" = "logit", "mu2" = "log", "mu3" = "identity",
                "sigma3" = "log",
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
      y3 <- y[i, 1]
      sigma3 <- par$sigma3[i]

      d <- c(
        dbinom(y1, size = 1, prob = mu1, log = log),
        dpois(y2, lambda = mu2, log = log),
        dnorm(y3, mean = mu3, sd = sigma3, log = log)
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        y_score <- y[, 1] - mu1
        y_score[y[, 2] != 1] <- 0
        y_score
      },
      mu2 = function(y, par, ...) {
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        y_score <- y[, 1] - mu2
        y_score[y[, 2] != 2] <- 0
        y_score
      },
      mu3 = function(y, par, ...) {
        y_score <- drop((y[, 1] - par$mu3 - par$Lambda)/(par$sigma3^2))
        y_score[y[, 2] != 3] <- 0
        y_score
      },
      sigma3 = function(y, par, ...) {
        y_score <- drop(-1 + (y[, 1] - par$mu3 - par$Lambda)^2/(par$sigma3^2))
        y_score[y[, 2] != 3] <- 0
        y_score
      },
      Lambda = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(y[i, 2],
                               y[i, 1] - mu1[i],
                               y[i, 1] - mu2[i],
                               drop((y[i, 1] - par$mu3[i] - par$Lambda[i]) /
                                      (par$sigma3[i]^2)))
        }
        y_score
      }
      ),
    "hess" = list(
      mu1 = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        hess <- mu1 * (1 - mu1)
        hess[y[, 2] != 1] <- 0
        hess
      },
      mu2 = function(y, par, ...) {
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        hess <- mu2
        hess[y[, 2] != 2] <- 0
        hess
      },
      mu3 = function(y, par, ...) {
        hess <- drop(1/(par$sigma3^2))
        hess[y[, 2] != 3] <- 0
        hess
      },
      sigma3 = function(y, par, ...) {
        hess <- rep(2, nrow(y))
        hess[y[, 2] != 3] <- 0
        hess
      },
      Lambda = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        hess <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          hess[i] <- switch(y[i, 2],
                            mu1[i] * (1 - mu1[i]),
                            mu2[i],
                            drop(1/(par$sigma3[i]^2)))
        }
        hess
      }
      )
  )
  class(f) <- "family.bamlss"
  return(f)
}

#' Draft of new family for gamlss2
#'
#' Fix three distributional assumptions but supply derivatives.
#' @param ... Not used.
#' @returns A gamlss2 family object.
famg <- function(...) {

  links <- list("mu1" = make.link("logit"),
                "mu2" = make.link("log"),
                "mu3" = make.link("identity"))

  f <- list(
    "names" = c("mu1", "mu2", "mu3", "sigma3", "Lambda"),
    "links" = c("mu1" = "logit", "mu2" = "log", "mu3" = "identity",
                "sigma3" = "log",
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
      y3 <- y[i, 1]
      sigma3 <- par$sigma3[i]

      d <- c(
        dbinom(y1, size = 1, prob = mu1, log = log),
        dpois(y2, lambda = mu2, log = log),
        dnorm(y3, mean = mu3, sd = sigma3, log = log)
      )

      return(d)
    },
    "predict" = gmfamm_predict,
    "score" = list(
      mu1 = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        y_score <- y[, 1] - mu1
        y_score[y[, 2] != 1] <- 0
        y_score
      },
      mu2 = function(y, par, ...) {
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        y_score <- y[, 1] - mu2
        y_score[y[, 2] != 2] <- 0
        y_score
      },
      mu3 = function(y, par, ...) {
        y_score <- drop((y[, 1] - par$mu3 - par$Lambda)/(par$sigma3^2))
        y_score[y[, 2] != 3] <- 0
        y_score
      },
      sigma3 = function(y, par, ...) {
        y_score <- drop(-1 + (y[, 1] - par$mu3 - par$Lambda)^2/(par$sigma3^2))
        y_score[y[, 2] != 3] <- 0
        y_score
      },
      Lambda = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        y_score <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          y_score[i] <- switch(y[i, 2],
                               y[i, 1] - mu1[i],
                               y[i, 1] - mu2[i],
                               drop((y[i, 1] - par$mu3[i] - par$Lambda[i]) /
                                      (par$sigma3[i]^2)))
        }
        y_score
      }
    ),
    "hess" = list(
      mu1 = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        hess <- mu1 * (1 - mu1)
        hess[y[, 2] != 1] <- 0
        hess
      },
      mu2 = function(y, par, ...) {
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        hess <- mu2
        hess[y[, 2] != 2] <- 0
        hess
      },
      mu3 = function(y, par, ...) {
        hess <- drop(1/(par$sigma3^2))
        hess[y[, 2] != 3] <- 0
        hess
      },
      sigma3 = function(y, par, ...) {
        hess <- rep(2, nrow(y))
        hess[y[, 2] != 3] <- 0
        hess
      },
      Lambda = function(y, par, ...) {
        eta1 <- links$mu1$linkfun(par$mu1) + par$Lambda
        mu1 <- links$mu1$linkinv(eta1)
        eta2 <- links$mu2$linkfun(par$mu2) + par$Lambda
        mu2 <- links$mu2$linkinv(eta2)
        hess <- rep(0, nrow(y))
        for (i in seq_len(nrow(y))) {
          hess[i] <- switch(y[i, 2],
                            mu1[i] * (1 - mu1[i]),
                            mu2[i],
                            drop(1/(par$sigma3[i]^2)))
        }
        hess
      }
    )
  )
  class(f) <- "gamlss2.family"
  return(f)
}

