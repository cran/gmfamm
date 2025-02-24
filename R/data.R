#' Subset of PBC data set for GMFAMM
#'
#' A subset of data from the \link[JMbayes2]{pbc2} data set which is the Mayo
#' Clinic Primary Biliary Cirrhosis Data, where only patients who survived at
#' least 10 years since they entered the study and were alive and had not had a
#' transplant at the end of the 10th year.
#'
#' Additionally, subject 124 is excluded as it has only one longitudinal
#' measurement per outcome. Function \link[registr]{gfpca_twoStep}, however,
#' assumes at least two longitudinal observations per subject.
#'
#' @format ## `pbc_gmfamm`
#' A data frame with 5,943 rows and 10 columns:
#' \describe{
#'   \item{id}{patients identifier; in the subset, there are 50 patients
#'     included.}
#'   \item{years}{number of years in the study without event}
#'   \item{status}{a factor with levels \code{alive}, \code{transplanted}, and
#'     \code{dead}.}
#'   \item{drug}{a factor with levels \code{placebo} and \code{D-penicilin}.}
#'   \item{age}{at registration in years.}
#'   \item{sex}{a factor with levels \code{male} and \code{female}.}
#'   \item{year}{number of years between enrollment and this visit date.}
#'   \item{status2}{a numeric vector with value \code{1} denotign if the patient
#'     was dead, and \code{0} if the patient was alive or transplanted.}
#'   \item{outcome}{a factor with levels \code{albumin}, \code{alkaline},
#'     \code{ascites}, \code{edema}, \code{hepatomegaly}, \code{histologic},
#'     \code{platelets}, \code{prothrombin}, \code{serBilir}, \code{serChol},
#'     \code{SGOT}, \code{spiders}.}
#'   \item{y}{value of the corresponding outcome at the visit date.}
#' }
#' @source \link[JMbayes2]{pbc2}
#' @references Hall et al. (2008): Modelling sparse generalized
#'  longitudinal observations with latent gaussian processes. Journal of the
#'  Royal Statistical Society Series B: Statistical Methodology, 70(4), 703-723.
"pbc_gmfamm"
