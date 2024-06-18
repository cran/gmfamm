#' Generalized Multivariate Functional Additive Models
#'
#' This package does things.
#' _PACKAGE
#' @name varbinq
#' @import bamlss
#' @import stats
#' @importFrom mgcv Predict.matrix gam
#' @importFrom utils getFromNamespace
#' @importFrom MASS ginv
#' @importFrom splines spline.des
#' @importFrom Matrix crossprod
NULL


#' Family object for bamlss for Generalized Multivariate Functional Additive
#'  Mixed Models
#'
#' @param family Vector of bamlss family names to construct the full family.
#' @param ... Not used at the moment.
#'
#' @returns An object of class \code{family.bamlss}
#' @export
#'
#' @examples
#'
#' # Short example to see how a family can be specified.
#' gmfamm(family = c("binomial", "poisson", "gaussian"))
#'
#' # Long example to see how an analysis can be done.
#' \donttest{
#' library(tidyverse)
#' library(registr)
#' library(funData)
#' library(MFPCA)
#' library(MJMbamlss)
#' library(refund)
#'
#' # Take only three outcomes (normal, binary, poisson)
#' # Log-transformation of serBilir to get normal distribution
#' pbc <- pbc_gmfamm %>%
#'   filter(outcome %in% c("serBilir", "hepatomegaly", "platelets")) %>%
#'   droplevels() %>%
#'   mutate(y = case_when(outcome == "serBilir" ~ log(y),
#'                        outcome != "serBilir" ~ y),
#'          year = ifelse(year > 9.99, 9.99, year))
#'
#' pbc_list <- split(pbc, pbc$outcome) %>%
#'   lapply(function (dat) {
#'     dat <- dat %>%
#'       mutate(value = y, index = year) %>%
#'       select(id, value, index) %>%
#'       arrange(id, index)
#'   })
#'
#' # Fit separate univariate GPFCAs
#' # Two numbers (x, y) in npc criterion indicate x% total variance but each pc
#' # hast to contribute at least y%
#' gfpcs <- mapply(function (data, fams) {
#'   gfpca_twoStep(Y = data, family = fams, npc_criterion = c(0.99, 0.001),
#'                 verbose = FALSE)
#' }, data = pbc_list, fams = list("binomial", "poisson", "gaussian"),
#' SIMPLIFY = FALSE)
#'
#' # Convert fitted values to funData
#' mfdata <- multiFunData(lapply(gfpcs, function (x) {
#'   funData(argvals = x$t_vec,
#'           X = matrix(x$Yhat$value, ncol = length(x$t_vec), byrow = TRUE))
#' }))
#'
#' # Convert estimated eigenfunctions to funData
#' uniexpansions <- lapply(gfpcs, function (x) {
#'   list(type = "given",
#'        functions =  funData(argvals = x$t_vec, X = t(x$efunctions)))
#' })
#'
#' # Calculate the maximal number of MFPCs
#' m <- sum(sapply(gfpcs, "[[", "npc"))
#'
#' # Estimate the MFPCs with weights 1
#' mfpca <- MFPCA(mFData = mfdata, M = m, uniExpansions = uniexpansions)
#'
#' # Choose number of MFPCs based on threshold
#' nfpc <- min(which(cumsum(mfpca$values) / sum(mfpca$values) > 0.95))
#'
#' # Attach estimated MFPCs
#' pbc <- attach_wfpc(mfpca, pbc, n = nfpc, marker = "outcome", obstime = "year")
#'
#' # Specify formula
#' f <- list(
#'   gm(y, outcome) ~ year + drug + sex, # hepatomegaly
#'   mu2 ~ year, # platelets
#'   mu3 ~ year + age, # serBilir
#'   sigma3 ~ 1, # serBilir sd
#'   Lambda ~ -1 + s(id, fpc.1, bs = "pcre") +
#'     s(id, fpc.2, bs = "pcre") + s(id, fpc.3, bs = "pcre") +
#'     s(id, fpc.4, bs = "pcre")
#' )
#'
#' b <- bamlss(f,
#'             family = gmfamm(c("binomial", "poisson", "gaussian")),
#'             data = pbc)
#'
#' }
gmfamm <- function(family, ...) {

  # Let bamlss construct the families
  family <- lapply(family, bamlss.family)

  # Create parameter names, ignoring original parameter names for ease of
  # implementation
  nam_len <- lapply(family, function(x) length(x$names))
  if (any(unlist(nam_len) > 2)) {
    stop(paste0("Family only implemented for distribution functions with one ",
                "or two distributional parameters"))
  }
  nams <- sapply(nam_len, function (x) c("mu", "sigma")[seq_len(x)])
  nams <- unlist(Map(paste0, nams, seq_along(nams)))
  nams <- c(nams, "Lambda")

  # Extract the link functions from the family argument
  links <- c(unlist(lapply(family, "[[", "links")), "identity")
  names(links) <- nams

  # Link functions for mu parameters
  linkfuns <- lapply(links, function (x) make.link(x))

  # Create the density function
  d <- function(y, par, log = FALSE, ...) {

    out <- lapply(seq_along(family), function (dim) {

      # Dimwise helper objects
      i <- y[, 2] == dim
      mu0 <- paste0("mu", dim)
      sigma0 <- paste0("sigma", dim)
      dimnames <- family[[dim]]$names

      # Construct mu and sigma
      mu <- linkfuns[[mu0]]$linkfun(par[[mu0]][i])
      mu <- linkfuns[[mu0]]$linkinv(mu + par$Lambda[i])
      sigma <- if (length(family[[dim]]$names) == 2) {
        par[[sigma0]][i]
      } else NULL

      # Update parameter list for dimension
      dimpar <- list(mu, sigma)
      names(dimpar) <- dimnames

      family[[dim]]$d(y = y[i, 1], par = dimpar, log = log)
    })
    unlist(out)
  }

  # # Create the score functions
  # scores <- lapply(seq_along(linkfuns), function(pred) {
  #   browser()
  #
  #   # First function for mu or sigma predictors, second for lambda
  #   dim <- grep("[0-9]+", names(linkfuns)[pred])
  #   sig <- names(linkfuns)[pred] == paste0("sigma", dim)
  #
  #   if(length(dim)) {
  #     function(y, par, ...) {
  #
  #       # Dimwise helper objects
  #       i <- y[, 2] == dim
  #       mu0 <- paste0("mu", dim)
  #       sigma0 <- paste0("sigma", dim)
  #       dimnames <- family[[dim]]$names
  #       y_score <- rep(0, nrow(y))
  #
  #       # Construct mu and sigma
  #       eta <- linkfuns[[mu0]]$linkfun(par[[mu0]][i]) + par$Lambda[i]
  #       mu <- linkfuns[[mu0]]$linkinv(eta)
  #       sigma <- if (length(family[[dim]]$names) == 2) {
  #         par[[sigma0]][i]
  #       } else NULL
  #
  #       # Update parameter list for dimension
  #       dimpar <- list(mu, sigma)
  #       names(dimpar) <- dimnames
  #
  #       y_score[i] <- family[[dim]]$score[[if (sig) 2 else 1]](
  #         y[i, 1], par = dimpar)
  #       y_score
  #
  #     }
  #   } else {
  #     function(y, par, ...) {
  #
  #       ndim <- length(family)
  #       mupar <- grep("mu", names(par))
  #       lapply(par[])
  #
  #     }
  #   }
  #
  #
  # })
  # names(scores) <- nams

  # list(
  #   delta = function(y, par, ...) {
  #
  #     eta1 <- links$mu1$linkfun(par$mu1) + par$delta
  #     mu1 <- links$mu1$linkinv(eta1)
  #     eta2 <- links$mu2$linkfun(par$mu2) + par$delta
  #     mu2 <- links$mu2$linkinv(eta2)
  #     eta3 <- links$mu3$linkfun(par$mu3) + par$delta
  #     mu3 <- links$mu3$linkinv(eta3)
  #     eta4 <- links$mu4$linkfun(par$mu4) + par$delta
  #     mu4 <- links$mu4$linkinv(eta4)
  #
  #     y_score <- rep(0, nrow(y))
  #     for (i in seq_len(nrow(y))) {
  #       y_score[i] <- switch(
  #         y[i, 2],
  #         bamlss.family("nbinom")$score$mu(
  #           y[i, 1], par = list("mu" = mu1[i], "theta" = par$sigma1[i])),
  #         bamlss.family("nbinom")$score$mu(
  #           y[i, 1], par = list("mu" = mu2[i], "theta" = par$sigma2[i])),
  #         bamlss.family("gamma")$score$mu(
  #           y[i, 1], par = list("mu" = mu3[i], "sigma" = par$sigma3[i])),
  #         bamlss.family("gamma")$score$mu(
  #           y[i, 1], par = list("mu" = mu4[i], "sigma" = par$sigma4[i])))
  #     }
  #     y_score
  #   }
  # )

  # Create output family
  f <- list(
    "names" = nams,
    "links" = links,
    "d" = d,
    "predict" = gmfamm_predict
  )

  class(f) <- "family.bamlss"
  f

}


#' Indicate Generalized Multivariate Model
#'
#' This function is used in the formula call of a generalized multivariate
#' functional additive mixed model to supply the information of the outcome and
#' factor variables to bamlss.
#'
#' @param y Name of variable in data set which contains the values of the
#'   longitudinal outcome.
#' @param outcome Name of variable in data set which is the factor variable
#'   indicating which outcome the value is from. Note that only the ordering
#'   not the factor levels are used in the estimation process.
#' @param ... Additional arguments not used at the moment.
#' @returns Matrix combining y and outcomes of class 'matrix' and 'gm'.
#' @export
#' @examples
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
#' eta_1 <- t + 0.5*x
#' eta_2 <- t + 0.5*x
#'
#' # Outcomes
#' y1 <- rnorm(n*ni, eta_1, 0.3)
#' y2 <- rbinom(n*ni, 1, 1/(1 + exp(-eta_2)))
#'
#' # Data format
#' dat <- data.frame(
#'    id = factor(rep(seq_len(n), each = ni)),
#'    y = c(y1, y2),
#'    dim = factor(rep(c(1, 2), each = n*ni)),
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
#'   Lambda ~ -1 + s(id, by = fpc, bs = "re")
#' )
#'
gm <- function(y, outcome, ...) {
  rval <- cbind(y, outcome)
  class(rval) <- c("matrix", "gm")
  rval
}

