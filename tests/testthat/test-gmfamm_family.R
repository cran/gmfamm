
test_that("family object is created correctly", {

  # Check whether family object can be created
  fam <- gmfamm(family = c("binomial", "poisson", "gaussian"))

  expect_equal(class(fam), "family.bamlss")
})


test_that("PBC example can be run", {

  pbc <- readRDS(test_path("fixtures", "pbc_example.rds"))

  # Specify formula
  form <- list(
    gm(y, outcome) ~ year + drug + sex, # hepatomegaly
    mu2 ~ year, # platelets
    mu3 ~ year + age, # serBilir
    sigma3 ~ 1, # serBilir sd
    Lambda ~ -1 + s(id, by = fpc.1, bs = "re") +
      s(id, by = fpc.2, bs = "re") + s(id, by = fpc.3, bs = "re") +
      s(id, by = fpc.4, bs = "re")
  )

  expect_warning({
    b <- bamlss(form,
                family = gmfamm(c("binomial", "poisson", "gaussian")),
                data = pbc, maxit = 10, n.iter = 20, burnin = 10)
  }, "the backfitting algorithm did not converge!")

  expect_equal(class(b), c("bamlss", "bamlss.frame", "list"))

})
