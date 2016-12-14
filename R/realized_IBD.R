#' @title Compute realized IBD coefficients.
#'
#' @details The realized IBD coefficient between two individuals is here defined as
#' the probability that two alleles randomly drawn each from one of the individuals
#' are IBD, i.e., refer to the same founder allele. This is analogous to the definintion of
#' the coefficient of coancestry in classical quantitative genetics, but the difference is
#' that while the coefficient of coancestry is calculated conditional on the pedigree and
#' hence can be regarded as an expected IBD coefficient, the realized IBD coefficient
#' takes into account mendelian sampling.
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @name realized_IBD
NULL

#' @rdname realized_IBD
#' @description Compute realized IBD coefficients between two individuals.
#'
#' @param ind1 first individual
#' @param ind2 second individual
#' @examples
#' L <- c(100, 200)
#' parents <- create_parents(2, L)
#' xo_params <- create_xo_params(L)
#' F1 <- cross(parents[[1]], parents[[2]], xo_params)
#' realized_IBD(self(F1, xo_params), self(F1, xo_params))
#' NULL
#' @export
realized_IBD <- function(ind1, ind2) {
  .realized_IBD(ind1, ind2)
  # n_chrom = length(ind1)
  # len <- 0
  # tot <- 0
  # for (i in seq_len(n_chrom)) {
  #   tmp <- ind1[[i]]$mat$locations
  #   len <- len + tmp[length(tmp)]
  #   tot <- tot +
  #     .realized_IBD_chromatids(ind1[[i]]$mat, ind2[[i]]$mat) +
  #     .realized_IBD_chromatids(ind1[[i]]$mat, ind2[[i]]$pat) +
  #     .realized_IBD_chromatids(ind1[[i]]$pat, ind2[[i]]$mat) +
  #     .realized_IBD_chromatids(ind1[[i]]$pat, ind2[[i]]$pat)
  # }
  # # print(tot); print(len)
  # tot / (4.0 * len)
}

#' @rdname realized_IBD
#' @description Compute the realized inbreeding coefficient of an individual.
#'
#' @param ind individual
#' @examples
#' # After three generations of selfing.
#' realized_F(self(self(self(F1, xo_params), xo_params), xo_params))
#' NULL
#' @export
realized_F <- function(ind) {
  realized_IBD(ind, ind)
}


#' @rdname realized_IBD
#' @description Compute realized IBD coefficients within one or between two
#' populations.
#'
#' @param pop A population as a list of indivdiuals.
#' @param second_pop A second population.
#' @examples
#' pop <- replicate(4, self(F1, xo_params), simplify = FALSE)
#' realized_IBD_matrix(pop)
#' @export
realized_IBD_matrix <- function(pop, second_pop = NULL) {
  n <- length(pop)
  if (is.null(second_pop)) {
    mat <- matrix(data = NA_real_, nrow = n, ncol = n)
    for (i in 1L:(n - 1L)) {
      for (j in (i + 1L):n) {
        mat[i, j] <- mat[j, i] <- realized_IBD(pop[[i]], pop[[j]])
      }
      mat[i, i] <- realized_IBD(pop[[i]], pop[[i]])
    }
    mat[n, n] <- realized_IBD(pop[[n]], pop[[n]])
    colnames(mat) <- rownames(mat) <- names(pop)
    return(mat)
  }
  n2 <- length(second_pop)
  mat <- matrix(data = NA_real_, nrow = n, ncol = n2)
  for (i in seq_len(n)) {
    for (j in seq_len(n2)) {
      mat[i, j] <- realized_IBD(pop[[i]], second_pop[[j]])
    }
  }
  rownames(mat) <- names(pop)
  colnames(mat) <- names(second_pop)
  mat
}

# realized_IBD_chromatids <- function(chromatid1, chromatid2) {
#   tot <- 0.0
#   mat <- chromatid1
#   pat <- chromatid2
#   matloc <- mat$locations
#   patloc <- pat$locations
#   matalle <- mat$alleles
#   patalle <- pat$alleles
#   ipatmax <- length(patloc)
#   imatmax <- length(matloc)
#   imat = 1L
#   ipat = 1L
#   cur_alle <- as.integer(matloc[ipat] <= patloc[imat]) # 0: pat, 1: mat
#   left <- 0.0
#
#   while(TRUE) {
#
#     if (cur_alle) {
#       repeat {
#         right <- matloc[imat]
#         if (patalle[ipat] == matalle[imat]) {
#           tot <- tot + (right - left)
#         }
#         # cat(sprintf("%.3f \t %.3f \n", left, right))
#         imat <- imat + 1L
#         left <- right
#         if (patloc[ipat] < matloc[imat] || imat > imatmax)
#           break
#       }
#       cur_alle <- 1L - cur_alle
#     } else {
#       repeat {
#         right <- patloc[ipat]
#         if (patalle[ipat] == matalle[imat]) {
#           tot <- tot + (right - left)
#         }
#         # cat(sprintf("%.3f \t %.3f \n", left, right))
#         ipat <- ipat + 1L
#         left <- right
#         if (patloc[ipat] >= matloc[imat] || ipat > ipatmax)
#           break
#       }
#       cur_alle <- 1L - cur_alle
#     }
#     if (ipat > ipatmax || imat > imatmax)
#       break
#   }
#   return(tot)
# }

