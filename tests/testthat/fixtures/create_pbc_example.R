

# Create PBC data for testing --------------------------------------------

old_dir <- getwd()
setwd(getSrcDirectory(function(){})[1])

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(tidyverse)
library(registr)
library(funData)
library(MFPCA)
library(MJMbamlss)
library(refund)

# Take only three outcomes (normal, binary, poisson)
# Log-transformation of serBilir to get normal distribution
pbc <- pbc_gmfamm %>%
  filter(outcome %in% c("serBilir", "hepatomegaly", "platelets")) %>%
  droplevels() %>%
  mutate(y = case_when(outcome == "serBilir" ~ log(y),
                       outcome != "serBilir" ~ y),
         year = ifelse(year > 9.99, 9.99, year))

pbc_list <- split(pbc, pbc$outcome) %>%
  lapply(function (dat) {
    dat <- dat %>%
      mutate(value = y, index = year) %>%
      select(id, value, index) %>%
      arrange(id, index)
  })

# Fit separate univariate GPFCAs
# Two numbers (x, y) in npc criterion indicate x% total variance but each pc
# hast to contribute at least y%
gfpcs <- mapply(function (data, fams) {
  gfpca_twoStep(Y = data, family = fams, npc_criterion = c(0.99, 0.001),
                verbose = FALSE)
}, data = pbc_list, fams = list("binomial", "poisson", "gaussian"),
SIMPLIFY = FALSE)

# Convert fitted values to funData
mfdata <- multiFunData(lapply(gfpcs, function (x) {
  funData(argvals = x$t_vec,
          X = matrix(x$Yhat$value, ncol = length(x$t_vec), byrow = TRUE))
}))

# Convert estimated eigenfunctions to funData
uniexpansions <- lapply(gfpcs, function (x) {
  list(type = "given",
       functions =  funData(argvals = x$t_vec, X = t(x$efunctions)))
})

# Calculate the maximal number of MFPCs
m <- sum(sapply(gfpcs, "[[", "npc"))

# Estimate the MFPCs with weights 1
mfpca <- MFPCA(mFData = mfdata, M = m, uniExpansions = uniexpansions)

# Choose number of MFPCs based on threshold
nfpc <- min(which(cumsum(mfpca$values) / sum(mfpca$values) > 0.95))

# Attach estimated MFPCs
pbc <- attach_wfpc(mfpca, pbc, n = nfpc, marker = "outcome", obstime = "year")
saveRDS(pbc, file = "tests/testthat/fixtures/pbc_example.rds")

setwd(old_dir)
