test_that("prediction for single term works", {
  b <- readRDS(test_path("fixtures", "fam_example.rds"))

  # Expectation for parameters
  m1 <- b$parameters$mu1$p["x1"]
  m2 <- b$parameters$mu2$p["x1"]
  m3 <- b$parameters$mu3$p["x1"]
  o1 <- c(rep(1, 500), rep(0, 1000))
  o2 <- c(rep(0, 500), rep(1, 500), rep(0, 500))
  o3 <- c(rep(0, 1000), rep(1, 500))
  expecto <- list("mu1" = b$model.frame[["x1"]]*m1*o1,
                  "mu2" = b$model.frame[["x1"]]*m2*o2,
                  "mu3" = b$model.frame[["x1"]]*m3*o3,
                  "sigma3" = rep(0, nrow(b$model.frame)),
                  "Lambda" = rep(0, nrow(b$model.frame)))

  out <- predict(b, term = "x1", type = "link", intercept = FALSE,
                 what = "parameters")
  expect_equal(out, expecto)

  expect_warning(predict(b, term = "x1", type = "parameter"),
                 "Check whether your code makes sense!")

  # Test if samples are correctly used
  s1 <- mean(b$samples[[1]][, "mu1.p.x1"])
  s2 <- mean(b$samples[[1]][, "mu2.p.x1"])
  s3 <- mean(b$samples[[1]][, "mu3.p.x1"])
  expecto_s <- list("mu1" = b$model.frame[["x1"]]*s1*o1,
                    "mu2" = b$model.frame[["x1"]]*s2*o2,
                    "mu3" = b$model.frame[["x1"]]*s3*o3,
                    "sigma3" = rep(0, nrow(b$model.frame)),
                    "Lambda" = rep(0, nrow(b$model.frame)))
  out_s <- predict(b, term = "x1", type = "link", intercept = FALSE,
                   what = "samples")
  expect_equal(out_s, expecto_s)

  # Test if sampling functions are correctly used
  sf1 <- matrix(c95(b$samples[[1]][, "mu1.p.x1"]), ncol = 3, byrow = TRUE,
                nrow = nrow(b$model.frame),
                dimnames = list(NULL, c("2.5%", "Mean", "97.5%")))
  sf2 <- matrix(c95(b$samples[[1]][, "mu2.p.x1"]), ncol = 3, byrow = TRUE,
                nrow = nrow(b$model.frame),
                dimnames = list(NULL, c("2.5%", "Mean", "97.5%")))
  sf3 <- matrix(c95(b$samples[[1]][, "mu3.p.x1"]), ncol = 3, byrow = TRUE,
                nrow = nrow(b$model.frame),
                dimnames = list(NULL, c("2.5%", "Mean", "97.5%")))
  nulldf <- data.frame(rep(0, 1500), rep(0, 1500), rep(0, 1500))
  colnames(nulldf) <- c("2.5%", "Mean", "97.5%")
  expecto_sf <- list("mu1" = as.data.frame(b$model.frame[["x1"]]*sf1*o1),
                     "mu2" = as.data.frame(b$model.frame[["x1"]]*sf2*o2),
                     "mu3" = as.data.frame(b$model.frame[["x1"]]*sf3*o3),
                     "sigma3" = nulldf,
                     "Lambda" = nulldf)
  expecto_sf <- lapply(expecto_sf, function (x) {
    attr(x, "row.names") <- as.character(attr(x, "row.names"))
    x
  })
  out_sf <- predict(b, term = "x1", type = "link", intercept = FALSE,
                    what = "samples", FUN = c95)
  expect_equal(out_sf, expecto_sf)

})

test_that("missing outcome variable throws error", {
  b <- readRDS(test_path("fixtures", "fam_example.rds"))

  newdata <- b$model.frame

  expect_error(predict(b, newdata = newdata, term = "x1"),
               "newdata needs to contain identifier for outcome variables")
  expect_error(predict(b, newdata = newdata),
               "newdata needs to contain identifier for outcome variables")
})

test_that("prediction for full model works", {
  b <- readRDS(test_path("fixtures", "fam_example.rds"))

  # Indicator for outcomes
  o1 <- c(rep(1, 500), rep(0, 1000))
  o2 <- c(rep(0, 500), rep(1, 500), rep(0, 500))
  o3 <- c(rep(0, 1000), rep(1, 500))
  outc <- factor(rep(1:3, each = 500))

  # Expectation for full prediction of type = link using optimized parameters
  x1 <- predict(b, term = "x1", type = "link", intercept = TRUE,
                what = "parameters")
  x2 <- predict(b, term = "x2", type = "link", intercept = FALSE,
                what = "parameters")
  fpc1 <- predict(b, term = "fpc.1", type = "link", intercept = FALSE,
                what = "parameters")
  fpc2 <- predict(b, term = "fpc.2", type = "link", intercept = FALSE,
                  what = "parameters")
  fpc3 <- predict(b, term = "fpc.3", type = "link", intercept = FALSE,
                  what = "parameters")
  pred_ls <- list("mu1" = x1$mu1 + x2$mu1 + fpc1$Lambda*o1 + fpc2$Lambda*o1 +
                    fpc3$Lambda*o1,
                  "mu2" = x1$mu2 + x2$mu2 + fpc1$Lambda*o2 + fpc2$Lambda*o2 +
                    fpc3$Lambda*o2,
                  "mu3" = x1$mu3 + x2$mu3 + fpc1$Lambda*o3 + fpc2$Lambda*o3 +
                    fpc3$Lambda*o3,
                  "sigma3" = b$parameters$sigma3$p * o3)
  expecto <- list("mu" = c(pred_ls$mu1[1:500], pred_ls$mu2[501:1000],
                           pred_ls$mu3[1001:1500]),
                  "sigma" = c(rep(NA, 1000), pred_ls$sigma3[1001:1500]))

  out <- predict(b, type = "link", what = "parameters")
  expect_equal(out, expecto)

  # Expectation for full prediction of type = link using samples
  p1 <- sapply(seq_len(nrow(b$model.frame)), function (s) {
    x <- b$model.frame[s, ]
    smp <- b$samples[[1]]
    lp <- smp[, "mu1.p.(Intercept)"] + x$x1 * smp[, "mu1.p.x1"] +
      x$x2 * smp[, "mu1.p.x2"] +
      smp[, paste0("Lambda.s.s(id,fpc.1).b", x$id)] * x$fpc.1 +
      smp[, paste0("Lambda.s.s(id,fpc.2).b", x$id)] * x$fpc.2 +
      smp[, paste0("Lambda.s.s(id,fpc.3).b", x$id)] * x$fpc.3
    mean(lp)
  })
  p2 <- sapply(seq_len(nrow(b$model.frame)), function (s) {
    x <- b$model.frame[s, ]
    smp <- b$samples[[1]]
    lp <- smp[, "mu2.p.(Intercept)"] + x$x1 * smp[, "mu2.p.x1"] +
      x$x2 * smp[, "mu2.p.x2"] +
      smp[, paste0("Lambda.s.s(id,fpc.1).b", x$id)] * x$fpc.1 +
      smp[, paste0("Lambda.s.s(id,fpc.2).b", x$id)] * x$fpc.2 +
      smp[, paste0("Lambda.s.s(id,fpc.3).b", x$id)] * x$fpc.3
    mean(lp)
  })
  p3 <- sapply(seq_len(nrow(b$model.frame)), function (s) {
    x <- b$model.frame[s, ]
    smp <- b$samples[[1]]
    lp <- smp[, "mu3.p.(Intercept)"] + x$x1 * smp[, "mu3.p.x1"] +
      x$x2 * smp[, "mu3.p.x2"] +
      smp[, paste0("Lambda.s.s(id,fpc.1).b", x$id)] * x$fpc.1 +
      smp[, paste0("Lambda.s.s(id,fpc.2).b", x$id)] * x$fpc.2 +
      smp[, paste0("Lambda.s.s(id,fpc.3).b", x$id)] * x$fpc.3
    mean(lp)
  })
  ps <- o3 * mean(b$samples[[1]][, "sigma3.p.(Intercept)"])
  expects <- list("mu" = c(p1[1:500], p2[501:1000], p3[1001:1500]),
                  "sigma" = c(rep(NA, 1000), ps[1001:1500]))
  outs <- predict(b, type = "link", what = "samples")
  expect_equal(outs, expects)

  # Exp


  # Test function compress_outcomes
  expect_equal(out, compress_outcomes(pred_list = pred_ls,
                                      mus = c("mu1", "mu2", "mu3"),
                                      sigmas = "sigma3",
                                      outcome = outc))


  # Expectation for full prediction of type = parameter and what = parameters
  resp1 <- function(x) 1/(1+exp(-x))
  resp2 <- function(x) exp(x)
  expecto2 <- list("mu" = c(resp1(expecto$mu[1:500]),
                            resp2(expecto$mu[501:1000]),
                            expecto$mu[1001:1500]),
                   "sigma" = c(rep(NA, 1000),
                               resp2(expecto$sigma[1001:1500])))

  out2 <- predict(b, type = "parameter", what = "parameters")
  expect_equal(out2, expecto2)

  # Expectation for full prediction of type = parameter and what = samples
  expects2 <- list("mu" = c(resp1(expects$mu[1:500]),
                            resp2(expects$mu[501:1000]),
                            expects$mu[1001:1500]),
                   "sigma" = c(rep(NA, 1000),
                               resp2(expects$sigma[1001:1500])))

  outs2 <- predict(b, type = "parameter", what = "samples")
  expect_equal(outs2, expects2)

})

test_that("prediction for single model terms works", {
  b <- readRDS(test_path("fixtures", "fam_example.rds"))

  # Expectation for predicting only model mu1
  x1 <- predict(b, term = "x1", type = "link", intercept = TRUE)
  x2 <- predict(b, term = "x2", type = "link", intercept = FALSE)
  expecto <- x1$mu1[1:500] + x2$mu1[1:500]
  out <- predict(b, model = "mu1")
  expect_equal(out, expecto)

  # Expectation for predicting only model lambda on additive predictor level
  fpc1 <- predict(b, term = "fpc.1", type = "link", intercept = FALSE)
  fpc2 <- predict(b, term = "fpc.2", type = "link", intercept = FALSE)
  fpc3 <- predict(b, term = "fpc.3", type = "link", intercept = FALSE)
  expecto2 <- fpc1$Lambda + fpc2$Lambda + fpc3$Lambda
  out2 <- predict(b, model = "Lambda")
  expect_equal(out2, expecto2)

  # Test apply_respfun_outcome
  resp1 <- function(x) 1/(1+exp(-x))
  resp2 <- function(x) exp(x)
  x <- rep(c(-2, 2), times = 3)
  expecto_fun <- c(-2, 2, resp2(-2), resp2(2), resp1(-2), resp1(2))
  outc <- factor(rep(1:3, each = 2))
  links <- c("mu1" = "identity", "mu2" = "log", "mu3" = "logit")
  out_fun <- apply_respfun_outcome(x = x, outcome = outc, links = links)
  expect_equal(out_fun, expecto_fun)

  # Expectation for predicting only model lambda on parameter level

  expecto3 <- c(resp1(expecto2[1:500]), resp2(expecto2[501:1000]),
                expecto2[1001:1500])
  expect_warning(out3 <- predict(b, model = "Lambda", type = "parameter"),
                 "Check whether your code makes sense!")
  expect_equal(out3, expecto3)

})

