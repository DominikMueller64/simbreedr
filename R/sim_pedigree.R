# sim_from_pedigree
#
#' Simulate genotypes for pedigree
#'
#' Simulate genotypes along one chromosome for a pedigree
#'
#' @param pedigree Matrix or data frame describing a pedigree, with first four
#' columns being individual ID, mom ID, dad ID, and sex (female as
#' \code{0}, male as \code{1}).
#' @param L Length of chromosome in cM (or a vector of chromosome lengths)
#' @param xchr If TRUE, simulate X chromosome. (If \code{L} is a
#' vector, this should be a vector of TRUE/FALSE values, of the same
#' length as \code{L}, or a character string with the name of the X
#' chromosome, in \code{L}.)
#' @param m Crossover interference parameter, for chi-square model
#' (m=0 corresponds to no interference).
#' @param p proportion of crossovers coming from no-interference process
#' @param obligate_chiasma If TRUE, require an obligate chiasma on the
#' 4-strand bundle at meiosis.
#'
#' @return A list with each component being the data for one
#' individual, as produced by the \code{\link{cross}} function.  Those
#' results are a list with two components, corresponding to the
#' maternal and paternal chromosomes. The chromosomes are represented
#' as lists with two components: an integer vector of alleles in
#' chromosome intervals, and a numeric vector of locations of the
#' right-endpoints of those intervals; these two vectors should have
#' the same length.
#'
#' If the input \code{L} is a vector, in order to simulate multiple
#' chromosomes at once, then the output will be a list with length
#' \code{length(L)}, each component being a chromosome and having the
#' form described above.
#'
#' @export
#' @keywords datagen
#' @seealso \code{\link{check_pedigree}},
#' \code{\link{sim_ril_pedigree}}, \code{\link{sim_ail_pedigree}},
#' \code{\link{sim_ril_pedigree}}
#'
#' @examples
#' # simulate AIL pedigree
#' tab <- sim_ail_pedigree(12, 30)
#' # simulate data from that pedigree
#' dat <- sim_from_pedigree(tab)
#' # simulate multiple chromosomes
#' dat <- sim_from_pedigree(tab, c("1"=100, "2"=75, "X"=100), xchr="X")

# id <- c(1, 2, 3, 4, 5, 6)
# mother <- c(NA, NA, 2, 2, 3, 4)
# father <- c(NA, NA, 1, 2, 3, 3)
# is_dh <- c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
# pedigree <- data.frame(id, mother, father, is_dh)
# xo_params <- create_xo_params(L = c(50, 100, 200))
#
# pop <- sim_pedigree(pedigree, xo_params)

sim_pedigree <- function(pedigree, xo_params) {
  list2env(xo_params, envir = environment())
  list2env(pedigree, envir = environment())

  nrow_ped <- nrow(pedigree)
  cur_allele <- 0L

  if(length(unique(id)) != nrow_ped)
    stop("IDs must be unique.")

  results <-  vector("list", nrow_ped)
  names(result) <- rownames(pedigree) <- as.character(id)

  for (i in seq_len(nrow_ped)) {
    homozygous <- is_dh[i]
    if (is.na(mother[i]) || is.na(father[i])) {
      if (homozygous) {
        alleles <- cur_allele
      } else {
        alleles <- cur_allele + c(0L, 1L)
      }
      parent <- create_parent(L, alleles, homozygous)
      results[[i]] <- parent
      cur_allele <- cur_allele + length(alleles)
    } else {
      mother_i <- match(x = mother[i], table = id)
      father_i <- match(x = father[i], table = id)

      if (is.na(mother_i) || is.na(father_i))
        stop("parents not found")
      if (mother_i >= i || father_i >= i)
        stop("Pedigree problem: parents follow individual")


      result[[i]] <- cross(results[[mother_i]], results[[father_i]],
                           xo_params)
    }
  }
  results
}
