#' Prediction of Generalized Multivariate Functional Additive Mixed model
#'
#' Note: FPC basis has to be evaluated for newdata before the predict function.
#'
#' Functionality of some arguments are restricted.
#'
#' @param object bamlss-model object to be predicted.
#' @param newdata Dataset for which to create predictions. Not needed for
#'   conditional survival probabilities.
#' @param type Character string indicating which type of predictions to compute.
#'   \code{link} returns the predictors of the corresponding model., \code{"parameter"} returns the estimates for all
#'   pedictors, \code{"probabilities"} returns the survival probabilities
#'   conditional on the survival up to the last longitudinal measurement, and
#'   \code{"cumhaz"} return the cumulative hazard up to the survival time or for
#'   a time window after the last longitudinal measurement. If \code{type} is
#'   set to \code{"loglik"}, the log-likelihood of the joint model is returned.
#'   Note that types \code{"probabilities"} and \code{"cumhaz"} are not yet
#'   implemented.
#' @param compress TRUE if the
#' @param FUN A function that should be applied on the samples of predictors or
#'   parameters, depending on argument \code{type}.
#' @param cores Specifies the number of cores that should be used for
#'   prediction. Note that this functionality is based on the
#'   \code{\link[parallel]{parallel}} package.
#' @param chunks Should computations be split into \code{chunks}? Prediction is
#'   then processed sequentially.
#' @param verbose Print information during runtime of the algorithm.
#' @param ... Currently not used.
#' @inheritParams bamlss::predict.bamlss
gmfamm_predict <- function(object, newdata, model = NULL, term = NULL,
                           match.names = TRUE, intercept = TRUE,
                           type = c("link", "parameter"), compress = TRUE,
                           FUN = function(x) { mean(x, na.rm = TRUE) },
                           trans = NULL, what = c("samples", "parameters"),
                           nsamps = NULL, verbose = FALSE, drop = TRUE,
                           cores = NULL, chunks = 1, ...)
{

  # Load internal bamlss function
  predict.bamlss <- getFromNamespace("predict.bamlss", "bamlss")
  make.link2 <- getFromNamespace("make.link2", "bamlss")

  # Sort out what to predict
  what <- match.arg(what)
  if(is.null(object$samples) & what == "samples") {
    what <- "parameters"
  }
  if (what == "parameters") {
    if(is.null(object$parameters))
      stop("cannot find any parameters!")
    FUN <- function(x) {x}
  }

  # Get information from object and newdata
  y_id <- as.character(object$formula$mu1$formula[[2]][[2]])
  outcome_id <- as.character(object$formula$mu1$formula[[2]][[3]])
  if (is.null(newdata)) {
    newdata <- cbind(factor(object$y[[1]][, 2]), object$model.frame)
    colnames(newdata)[1] <- outcome_id
  }
  if(!outcome_id %in% colnames(newdata)) {
    stop("newdata needs to contain identifier for outcome variables")
  }

  # Create an identifier for which observations belong to which outcome
  outcome_ids <- model.matrix(~ -1 + newdata[, outcome_id])
  rownames(outcome_ids) <- NULL
  outcome_levels <- levels(newdata[, outcome_id])

  # Get info from formula on contained predictors
  preds <- names(object$formula)
  mu_names <- preds[grep("^mu[0-9]+", preds)]
  mus <- if (length(mu_names) > 0) {
    as.integer(gsub("mu", "", mu_names))
  } else NULL
  sigma_names <- preds[grep("^sigma[0-9]+", preds)]
  sigmas <- if (length(sigma_names) > 0) {
    as.integer(gsub("sigma", "", sigma_names))
  } else NULL

  # Remove predict family from object so that the original bamlss functions can
  # be applied
  family <- object$family
  object$family$predict <- NULL
  type <- match.arg(type)

  # If only specific terms are wanted, use the default predict function
  if (!is.null(term)) {

    if (type == "parameter") {
      warning("Check whether your code makes sense!")
    }
    pred_term <-
      predict.bamlss(object = object, newdata = newdata, model = model,
                     term = term, match.names = match.names,
                     intercept = intercept, type = type, FUN = FUN,
                     trans = trans, what = what, nsamps = nsamps,
                     verbose = verbose, drop = drop, cores = cores,
                     chunks = chunks, ...)
    pred_term <- incorporate_outcome(pred_list = pred_term, mus = mus,
                                     sigmas = sigmas, outcome_ids = outcome_ids,
                                     outcome_levels = outcome_levels)
    return(pred_term)
  }

  # If only specific model is wanted, check whether it is a Lambda model, then
  # for type = 'parameter' handle it differently
  if (!is.null(model)) {
    if (model == "Lambda") {
      pred_model <- predict.bamlss(
        object = object, newdata = newdata, model = model, term = term,
        match.names = match.names, intercept = intercept,
        type = "link",
        FUN = FUN, trans = NULL, what = what, nsamps = nsamps,
        verbose = verbose, drop = drop, cores = cores, chunks = chunks, ...)

      if (type == "parameter") {
        warning("Check whether your code makes sense!")

        pred_model <- apply_respfun_outcome(x = pred_model,
                                            outcome = newdata[, outcome_id],
                                            links = object$family$links)

      }

      return(pred_model)

    } else {
      pred_model <-
        predict.bamlss(object = object, newdata = newdata, model = model,
                       term = term, match.names = match.names,
                       intercept = intercept, type = type, FUN = FUN,
                       trans = trans, what = what, nsamps = nsamps,
                       verbose = verbose, drop = drop, cores = cores,
                       chunks = chunks, ...)
      outcome_level <- as.integer(sub("mu|sigma", "", model))
      pred_model <- pred_model[outcome_ids[, paste0("newdata[, outcome_id]",
                                                      outcome_level)] == 1]
      return(pred_model)
    }
  }

  # Simple prediction for parameter prediction
  if (what == "parameters") {

    # Use first the predict function for the additive predictor
    preds_link <- predict.bamlss(
      object = object, newdata = newdata, model = model, term = term,
      match.names = match.names, intercept = intercept,
      type = "link",
      FUN = FUN, trans = NULL, what = what, nsamps = nsamps, verbose = verbose,
      drop = drop, cores = cores, chunks = chunks, ...)

  } else {

    # Construct the full additive predictors from the samples
    preds_link <- predict.bamlss(
      object = object, newdata = newdata, model = model, term = term,
      match.names = match.names, intercept = intercept,
      type = "link",
      FUN = function(x) x, trans = NULL, what = what, nsamps = nsamps,
      verbose = verbose, drop = drop, cores = cores, chunks = chunks, ...)
    preds_link <- lapply(preds_link, as.matrix)

  }

  # Combine the additive predictors for mu with mfpc-Lambda
  for (mu in mus) {
    preds_link[[paste0("mu", mu)]] <-
      preds_link[[paste0("mu", mu)]] + preds_link[["Lambda"]]
  }
  preds_link$Lambda <- NULL

  # Summary function for samples has not yet been applied
  if (what == "samples") {
    preds_link <- lapply(preds_link, function (prd) {
      apply(prd, 1, FUN)
    })
  }

  # Apply response functions to linear predictors
  if (type == "parameter") {
    for (pred_name in names(preds_link)) {
      preds_link[[pred_name]] <-
        make.link2(object$family$links[[pred_name]])$linkinv(
          preds_link[[pred_name]]
        )
    }
  }

  # Combine information to single vectors for mu and sigma models
  preds_list <- compress_outcomes(pred_list = preds_link, mus = mu_names,
                                  sigmas = sigma_names,
                                  outcome = newdata[, outcome_id])
  preds_list

}

#' Incorporate outcome information into
#'
#' This is an internal function multiplying all outcome predictions with 0 if
#' the respective row is not part of the outcome.
#'
#' @param pred_list List of predictions for each outcome.
#' @param mus Integer vector for numbering available mus. Can be NULL but
#'   shouldn't.
#' @param sigmas Integer vector for numbering available sigmas. Can be NULL.
#' @param outcome_ids Numeric matrix resulting from model.matrix call containing
#'   the info about the outcomes. Column names are hard coded.
#' @param outcome_levels Character string containing the outcome names.
#' @returns List but now with 0 elements where the rows are not corresponding to
#'   outcomes.
incorporate_outcome <- function (pred_list, mus, sigmas, outcome_ids,
                                 outcome_levels) {

  if (!is.null(mus)) {
    for (mu in mus) {
      pred_list[[paste0("mu", mu)]] <-
        pred_list[[paste0("mu", mu)]] *
        outcome_ids[, paste0("newdata[, outcome_id]", outcome_levels[mu])]
    }
  }

  if (!is.null(sigmas)) {
    for (sigma in sigmas) {
      pred_list[[paste0("sigma", sigma)]] <-
        pred_list[[paste0("sigma", sigma)]] *
        outcome_ids[, paste0("newdata[, outcome_id]", outcome_levels[sigma])]
    }
  }

  pred_list

}

#' Compress the outcome list of predictions into single vectors
#'
#' This is an internal function combining all mu and sigma outcomes,
#' respectively, taking into account the outcome information.
#'
#' @param pred_list List of predictions for each outcome.
#' @param mus Character vector with names of included mu models.
#' @param sigmas Character vector with names of included sigma models.
#' @param outcome Factor vector containing the information of which row
#'   corresponds to which outcome.
#' @returns List with two elements containing predictions for mu and sigma model
#'   terms. If a some model parameters are missing (such as sigma for binomial
#'   distributional assumption) NA elements are contained.
compress_outcomes <- function (pred_list, mus, sigmas, outcome) {

  # Is the output a vector or a matrix?
  mat <- is.matrix(pred_list[[1]])
  if (mat) {
    mu_pred <- matrix(NA, nrow = length(outcome), ncol = nrow(pred_list[[1]]))
    sigma_pred <- matrix(NA, nrow = length(outcome),
                         ncol = nrow(pred_list[[1]]))
  } else {
    mu_pred <- rep(NA, length(outcome))
    sigma_pred <- rep(NA, length(outcome))
  }

  # Combine mu elements
  for(i in seq_along(outcome)) {
    mu_i <- paste0("mu", as.numeric(outcome[i]))
    if (mu_i %in% mus) {
      if (mat) {
        mu_pred[i, ] <- pred_list[[mu_i]][, i]
      } else {
        mu_pred[i] <- pred_list[[mu_i]][i]
      }
    }
  }

  # Combine sigma elements
  for(i in seq_along(outcome)) {
    sigma_i <- paste0("sigma", as.numeric(outcome[i]))
    if (sigma_i %in% sigmas) {
      if (mat) {
        sigma_pred[i, ] <- pred_list[[sigma_i]][, i]
      } else {
        sigma_pred[i] <- pred_list[[sigma_i]][i]
      }
    }
  }

  # Attach colnames if available
  if (mat) {
    colnames(mu_pred) <- rownames(pred_list[[1]])
    colnames(sigma_pred) <- rownames(pred_list[[1]])
  }

  preds <- list(
    "mu" = mu_pred,
    "sigma" = sigma_pred
  )
  preds

}

#' Apply link functions based on outcome information
#'
#' This is an internal function for the extreme case that a vector is plugged
#' into a response function depending on outcome information.
#'
#' @param x Vector of additive predictors.
#' @param outcome Factor vector containing information on the outcome of
#'   the corresponding element of vector x.
#' @param links Vector containing the names of the respective links for the mu
#'   outcomes.
#' @returns Vector of lenght x where different response functions have been
#'   applied.
apply_respfun_outcome <- function(x, outcome, links) {

  make.link2 <- getFromNamespace("make.link2", "bamlss")

  out <- x
  for (i in seq_along(as.numeric(outcome))) {
    out[i] <- make.link2(links[paste0("mu", outcome[i])])$linkinv(x[i])
  }

  out
}
