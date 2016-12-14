#' Convert to genotypes.
#'
#' @description Convert crossover data to genotypic data
#'
#' @param ind Crossover data for single indivdiuals, as produced by functions
#' \code{\link{create_parents}}, \code{\link{cross}}, \code{\link{self}} or \code{\link{dh}}.
#'
#' @param map A list with genetic positions of all loci on each chromosome, must be
#' increasingly sorted.
#'
#' @param founder A list with the genotypes of the founder alleles. Each element of the
#' list refers to one chromosome and is a numeric matrix where the number of columns
#' corresponds to the number of loci on the respective chromosomes. The rows of the matrix refer
#' to the genotypes of the founder alleles as assigned by \code{\link{create_parents}}.
#'
#' @param level The level of condensation of the output, one of 'chromatid', 'homolog' or
#' 'individual'. The level of aggregation increasing in this order.
#'
#' @details It is the responsibility of the user to make sure that matrices in
#' \code{founder} are correctly sorted (both rows and columns) and that genotypes of all
#' founder alleles present in \code{ind} are reported as rows in those matrices.
#'
#' @return The return value depends on the argument to \code{level}.
#' If \code{level = 'chromatid'}, chromatids are separated as a list of vectors
#' (same structure as the original \code{ind}).
#' If \code{level = 'homolog'}, chromatids  are row-bounded into a matrix.
#' If \code{level = 'individual'}, homologs are column-bounded into a matrix.
#' consecutive rows refer the one individual.
#'
#' @seealso \code{\link{create_parents}}
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @examples
#' set.seed(1234)
#' n_loci <- 3L
#' L <- c(80, 100, 70, 120)
#' xo_params <- create_xo_params(L)
#' parents <- create_parents(2, L)
#' map <- lapply(L, function(x) sort(runif(n_loci, 0, x)))
#' founder <- replicate(length(L), simplify = FALSE,
#'  matrix(sample(c(-1, 1), size = 2 * n_loci, replace = TRUE),  nrow = 2)
#' )
#' F1 <- cross(parents[[1]], parents[[2]], xo_params)
#' DH <- dh(F1, xo_params)
#' convert(F1, map, founder, level = 'individual')
#' convert(DH, map, founder, level = 'individual')
#'
#' @export
convert <- function(ind, map, founder,
                    level = c('chromatid', 'homolog', 'individual')) {

  level = match.arg(arg = level)

  if (length(map) != length(founder) || length(map) != length(ind))
    stop("'map' and 'founder' must be lists of equal length.")

  geno <- xo2geno(ind, map, founder)
  if (level == 'chromatid')
    return(geno)
  geno <- purrr::at_depth(geno, .depth = 1L, ~ do.call(what = rbind, args = .x))
  if (level == 'homolog')
    return(geno)
  geno <- purrr::at_depth(geno, .depth = 0L, ~ do.call(what = cbind, args = .x))
  if (level == 'individual')
    return(geno)
}
