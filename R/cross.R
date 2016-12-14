#' @title Cross individuals, self individual, create doubled haploid.
#'
#' @param mother An individual.
#' @param father Another individual.
#' @param xo_params A list with crossover simulation parameters as produced by
#' \code{\link{create_xo_params}}.
#'
#' @return A list of the progeny's chromosomes. For \code{dh}, this is reduced to a single
#' chromatid to save memory space.
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @seealso \code{\link{create_xo_params}}, \code{\link{create_parents}}
#' @name cross
NULL

#' @rdname cross
#' @description Simulate the cross of two individuals.
#' @examples
#' L <- c(100, 200)
# 'parents <- create_parents(2, L)
#' xo_params <- create_xo_params(L)
#' cross(parents[[1]], parents[[2]], xo_params)
#' @export
cross <- function(mother, father, xo_params) {
  # .cross(mother, father, xo_params$L, xo_params$m, xo_params$p,
  #        xo_params$obligate_chiasma, xo_params$Lstar)
  .cross(mother = mother, father = father, params = xo_params)
}

#' @rdname cross
#' @description Self an individual.
#' @param parent An individual.
#' @examples
#' self(parents[[1]], xo_params)
#' @export
self <- function(parent, xo_params) {
  # .cross(parent, parent, xo_params$L, xo_params$m, xo_params$p,
         # xo_params$obligate_chiasma, xo_params$Lstar)
  .cross(mother = parent, father = parent, params = xo_params)
}

#' @rdname cross
#' @description Derive a doubled haploid from an indivdiual.
#' @param parent An individual.
#' @examples
#' dh(parents[[1]], xo_params)
#' @export
dh <- function(parent, xo_params) {
  ret <- .gamete(parent = parent, params = xo_params)
  attr(x = ret, which = 'homozygous') <- TRUE
  attr(x = ret, which = 'individual') <- TRUE
  ret
}

#' @rdname cross
#' @description Derive a gamete from an indivdiual.
#' @param parent An individual.
#' @examples
#' gamete(parents[[1]], xo_params)
#' @export
gamete <- function(parent, xo_params) {
  .gamete(parent = parent, params = xo_params)
}

#' @rdname cross
#' @description Simulate meiosis.
#' @param parent An individual.
#' @examples
#' sim_meiosis(parents[[1]], xo_params)
#' @export
sim_meiosis <- function(parent, xo_params) {
  .sim_meiosis(parent, L = xo_params$L, m = xo_params$m, p = xo_params$p,
               obligate_chiasma = xo_params$obligate_chiasma, Lstar = xo_params$Lstar)
}
