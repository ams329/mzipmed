#### Dataset documentation: mzipmed_data ####

#' Data to be used in the mzipmed package examples
#'
#' @format A dataframe with 500 rows and 10 variables.
#'  \describe{
#'     \item{X}{Simulated binary exposure ~Bernoulli(0.5)}
#'     \item{C1}{Simulated covariate ~Normal(0,1)}
#'     \item{C2}{Simulated covariate ~Beta(2,2)}
#'     \item{ziM}{Zero-inflated count mediator based on X,C1,C2}
#'     \item{lmM}{Continuous mediator based on X,C1,C2 with error term ~Normal(0,4)}
#'     \item{binM}{Binary mediator based on X,C1,C2}
#'     \item{lmY}{Continuous outcome to be used for ziM}
#'     \item{binY}{Binary outcome to be used for ziM}
#'     \item{ziY1}{Zero-inflated count outcome to be used for lmM}
#'     \item{ziY2}{Zero-inflated count outcome to be used for binM}
#'  }
#'
#'  @source {Simulated to serve as an example}
#'
#'  @examples
#'  data(mzipmed_data)
"mzipmed_data"
