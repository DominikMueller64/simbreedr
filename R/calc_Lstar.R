#' Calculate adjusted chromosome length for obligate chiasma.
#'
#' Calculate the reduced chromosome length that will give the target
#' expected number of chiasmata when conditioning on there being at
#' least one chiasma on the four-strand bundle.
#'
#' @param L Length of chromosome (in cM); must be > 50
#' @param m Interference parameter for chi-square model
#' @param p Proportion of chiasmata coming from no-interference
#' process
#'
#' @return Adjusted length of chromosome
#'
#' @author Karl Broman, minor modifications by Dominik Mueller
#'
##' @keywords utilities
#'
#' @importFrom stats uniroot dpois
#'
#' @seealso \code{\link{cross}}, \code{\link{sim_meiosis}},
#' \code{\link{sim_crossovers}}
#'
#' @examples
#' calc_Lstar(100, 0, 0)
#' calc_Lstar(60, 10, 0.1)
#'
#' @export
calc_Lstar <- function(L, m = 0L, p = 0.0) {

    if (L <= 50)
      stop("Must have L > 50")

    if (m < 0)
      stop("Must have m >= 0")

    if (!is.integer(m) && !isTRUE(all.equal(m, round(m)))) {
        warning("m must be an non-negative integer; rounding")
        m <- round(m)
    }
    if (p < 0 || p > 1)
        stop("p must be in [0, 1]")

    if (isTRUE(all.equal(p, 1.0))) { # if p == 1, might as well take m=0, p=0
        m <- p <- 0
    }

    func_to_zero <- function(Lstar) {
        if (m == 0L)
            denom <- 1.0 - exp(-Lstar/50)
        else {
            lambda1 <- Lstar / 50 * (m + 1.0) * (1.0 - p)
            lambda2 <- Lstar / 50 * p
            s <- seq.int(0L, m, by = 1L)
            denom <- 1.0 - sum(dpois(s, lambda1) * (m + 1.0 - s) / (m + 1.0)) * exp(-lambda2)
        }
        2.0 * L - 2.0 * Lstar / denom
    }

    uniroot(func_to_zero, interval = c(1e-8, L), tol = sqrt(.Machine$double.eps))$root
}
