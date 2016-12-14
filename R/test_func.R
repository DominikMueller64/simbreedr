#' Simulate a chromatid.
#'
#' @description Simulate a single chromatid.
#'
#' @param nX0 A numeric. The number of crossover locaitons.
#' @param founder_alleles A numeric vector. The founder alleles used.
#' @param len A numeric. The length of the chromosome.
#'
#' @return An list with elements \code{alleles} and \code{locations}.
#'
#' @details This function is only for testing purposes.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' nXO <- 5L
#' len <- 100L
#' sim_xodat_chromatid(nXO, len, c(0L, 1L))
#'
#' @export
.sim_xodat_chromatid <-
  function(nXO, len = len, founder_alleles = c(0L, 1L)){
      locations <- c(sort(runif(n = nXO, min = 0.0, max = len)), len)
      alleles <- sample(x = founder_alleles, size = nXO + 1L, replace = TRUE)
      list('alleles' = alleles, 'locations' = locations)
  }

#' Simulate an individual.
#'
#' @description Simulate crossover data of an individual.
#'
#' @param nX0 A numeric vector. The number of crossover locations on each chromosome.
#' @param len A numeric vector. The length of the chromosomes.
#' @param founder_alleles A numeric vector. The founder alleles used.
#' @param homozygous A logical. Is the individual homozygous?
#'
#' @return An individual.
#'
#' @details This function is only for testing purposes.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' nXO <- c(1L, 3L, 6L)
#' len <- c(53, 103, 249)
#' founder_alleles <- c(0L, 1L)
#' .sim_xodat_individual(nXO, len, founder_alleles, homozygous = FALSE)
#' .sim_xodat_individual(nXO, len, founder_alleles, homozygous = TRUE)
#'
#' @export
.sim_xodat_individual <-
  function(nXO, len = rep(x = 100.0, times = length(nXO)), founder_alleles = c(0L, 1L),
           homozygous = FALSE) {
    tmp <- purrr::map2(nXO, len, function(nxo, l) {
      if (homozygous) {
        .sim_xodat_chromatid(nXO = nxo, founder_alleles = founder_alleles, len = len)
      } else {
        purrr::set_names(x = purrr::rerun(2L, {
          .sim_xodat_chromatid(nXO = nxo, founder_alleles = founder_alleles, len = l)
        }), nm = c('pat', 'mat'))
      }
    })
    attr(x = tmp, which = 'homozygous') <- homozygous
    attr(x = tmp, which = 'individual') <- TRUE
    tmp
  }


#' Simulate map
#'
#' @description Simulate a genetic map
#'
#' @param m A numeric vector. The number of markers on each chromosome.
#' @param len A numeric vector. The length of each chromosome.
#'
#' @return A list with the (incresingly sorted) map positions.
#'
#' @details This function is only for testing purposes.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' m <- c(135L, 534L, 834L)
#' len <- c(53, 103, 249)
#' sim_map(m, len)
#' sim_map(m)
#'
#' @export
.sim_map <- function(m, len = rep(x = 100.0, times = length(m))) {
  purrr::map2(m, len, ~ sort(stats::runif(n = .x, min = 0.0, max = .y)))
}
#' Simualte founder.
#'
#' @description Simulate founder genotypes.
#'
#' @param n An integer. The number of founders.
#' @param m A integer vector. The number of loci on each chromosome.
#'
#' @return A list of matrices with biallelic loci, coded with 0 and 1.
#'
#' @details This function is only for testing purposes.
#'
#' @author Dominik Mueller (\email{dominikmueller64@@yahoo.de})
#'
#' @examples
#' n <- 2L
#' m <- c(5L, 7L, 3L)
#' sim_founder(n ,m)
#'
#' @export
sim_founder <- function(n, m) {
  purrr::map(m, ~ matrix(data = sample(x = c(0, 1), size = n * .x, replace = TRUE), ncol = n))
}
