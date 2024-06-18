

# Create model for testing ------------------------------------------------

old_dir <- getwd()
setwd(getSrcDirectory(function(){})[1])

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(MJMbamlss)
library(funData)
library(mvtnorm)
library(tidyverse)

# Number of individuals and other quantities
n <- 100
argvals <- seq(0, 1, by = 0.01)
x <- seq(0, 1, by = 0.1)

# Random covariance matrix
# Set the eigenvalues but the eigenvectors are random
set.seed(1705)
p <- 3
P <- qr.Q(qr(matrix(rnorm(p^2), ncol = p)))
evals <- c(4, 1, 0.5)
cov <- crossprod(P, P*(evals))

# Find spline functions
# Outcome1
m1sp1 <- splinefun(x, c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
# Outcome2
m2sp2 <- splinefun(x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
# Outcome3
m3sp2 <- splinefun(x, c(0, 0.2, 0.2, 0, -0.2, -0.2, 0, 0.2, 0.2, 0, -0.2))

m1 <- funData(argvals = argvals,
              X = matrix(m1sp1(argvals), nrow = 1, byrow = TRUE))
m2 <- funData(argvals = argvals,
              X = matrix(m2sp2(argvals), nrow = 1, byrow = TRUE))
m3 <- funData(argvals = argvals,
              X = matrix(m3sp2(argvals), nrow = 1, byrow = TRUE))

# True multivariate covariance structure
m <- MFPCA_cov(cov = cov, basis_funs = list(m1, m2, m3))

# Sample the covariates
set.seed(1450)
rho <- rmvnorm(n = n, sigma = diag(m$values))
x1 <- sample(c(0, 1), size = n, replace = TRUE)
x2 <- runif(n)

# Construct the data.frame
dat <- data.frame(id = factor(seq_len(n)),
                  x1 = x1,
                  x2 = x2,
                  rho = rho)
dat <- dat[rep(seq_len(n), each = length(argvals)), ]
dat$obstime <- rep(argvals, times = n)
dat <- rbind(dat, dat, dat)
dat$outcome <- factor(rep(c("hepatomegaly", "platelets", "serBilir"),
                          each = n*length(argvals)))
dat <- attach_wfpc(m, dat, n = 3, marker = "outcome")

# Construct the linear predictors and mu
pbc <- dat %>%
  mutate(
    eta = case_match(
      outcome,
      "hepatomegaly" ~ 1 - 0.5*x1 + x2^2 + rho.1*fpc.1 + rho.2*fpc.2 +
        rho.3*fpc.3,
      "platelets" ~ 1 + 0.5*x1 - x2^2 + rho.1*fpc.1 + rho.2*fpc.2 + rho.3*fpc.3,
      "serBilir" ~ 1 + 2*x1 + rho.1*fpc.1 + rho.2*fpc.2 + rho.3*fpc.3
    ),
    mu = case_match(
      outcome,
      "hepatomegaly" ~ 1 / (1 + exp(-eta)),
      "platelets" ~ exp(eta),
      "serBilir" ~ eta
    ))

# Draw observations
set.seed(1619)
pbc$y <- c(
  rbinom(n*length(argvals), size = 1,
         prob = pbc$mu[pbc$outcome == "hepatomegaly"]),
  rpois(n*length(argvals), lambda = pbc$mu[pbc$outcome == "platelets"]),
  rnorm(n*length(argvals), mean = pbc$mu[pbc$outcome == "serBilir"],
        sd = 0.01))

# Scale down the data set
pbc <- pbc[pbc$obstime %in% seq(0, 1, by = 0.2), ]

# Specify formula
f <- list(
  gm(y, outcome) ~ x1 + x2^2, # hepatomegaly
  mu2 ~ x1 + x2^2, # platelets
  mu3 ~ x1 + x2^2, # serBilir
  sigma3 ~ 1, # serBilir sd
  Lambda ~ -1 + s(id, fpc.1, bs = "mfpcre") +
    s(id, fpc.2, bs = "mfpcre") + s(id, fpc.3, bs = "mfpcre")
)
b <- bamlss(f, data = pbc, family = gmfamm:::fam2, optimizer = opt_bfit)
saveRDS(b, file = "tests/testthat/fixtures/fam_example.rds")

setwd(old_dir)
