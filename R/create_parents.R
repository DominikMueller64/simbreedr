#' Create parents.
#'
#' Create parents.
#'
#' @param n Number of parents to create.
#' @param len Vector of chromosome lengths in cM.
#' @param homozygous Should the parents be considered completely inbred?
#'
#' @return A nested list of lists where the first level refers
#' to the individuals and the second level to the different
#' chromosomes.
#' @seealso \code{\link{cross}}, \code{\link{sim_meiosis}}
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @examples
#' create_parents(n = 2, L = c(100, 200), inbred = FALSE)
#'
#' @export

create_parents <- function(n, len, homozygous = TRUE) {

  if (n < 1.0 || !isTRUE(all.equal(n, round(n))))
    stop("'n' must be a positive integer.")

  if (!is.atomic(len) || !is.numeric(len) || any(len <= 0.0))
    stop("'len' must be a numeric vector of (positive) chromosome lengths in CentiMorgan.")

  n <- round(n)
  len <- as.double(len)
  homozygous <- as.logical(homozygous)

  # pn <- paste('individual', seq_len(n), sep = ':')
  # cn <- paste('homolog', seq_along(L), sep = ':')

  ret <- replicate(n = n, expr = vector('list', length(len)), simplify = FALSE)

  # ret <- replicate(n = n,
  #           expr = setNames(vector('list', length(len)), cn),
  #           simplify = FALSE)
  # names(ret) <- pn

  cur_allele <- 0L
  for (i in seq_len(n)) {
    for (j in seq_along(len)) {
      l <- len[j]
      par <- list('alleles' = cur_allele, 'locations' = l)
      if (homozygous) {
        # tmp <- list('mat' = par, 'pat' = par)
        tmp <- par
      } else {
        pat <- list('alleles' = cur_allele + 1L, 'locations' = l)
        tmp <- list('mat' = par, 'pat' = pat)
      }
      # attr(tmp, 'homozygous') = homozygous
      # attr(tmp, 'L') = Lj
      ret[[i]][[j]] <- tmp
    } # for (j in seq_along(ret))
    attr(ret[[i]], 'homozygous') = homozygous
    attr(x = ret[[i]], which = 'individual') <- TRUE
    # class(ret[[i]]) <- append(class(ret[[i]]), values = 'individual')
    if (homozygous) cur_allele <- cur_allele + 1L
    else cur_allele <- cur_allele + 2L
  } # for (i in seq_along(L))
  ret
}

#' @export
create_parent <- function(len, alleles, homozygous = TRUE) {

  if (!is.atomic(len) || !is.numeric(len) || any(len <= 0.0))
    stop("'len' must be a numeric vector of (positive) chromosome lengths in CentiMorgan.")

  if (!is.integer(alleles))
    stop("'alleles must be an integer vector'")

  if (homozygous) {
    if (length(alleles) != 1L)
      stop("If 'homozygous = TRUE', 'alleles' must be of length 1.")
  } else {
    if (length(alleles) != 2L)
      stop("If 'homozygous = FALSE', 'alleles' must be of length 2.")
  }

  len <- as.double(len)
  homozygous <- as.logical(homozygous)
  # cn <- paste('homolog', seq_along(L), sep = ':')
  ret <- vector('list', length(len))

  for (j in seq_along(len)) {
    l <- len[j]
    par <- list('alleles' = alleles[1L], 'locations' = l)
    if (homozygous) {
      # tmp <- list('mat' = par, 'pat' = par)
      tmp <- par
    } else {
      pat <- list('alleles' = alleles[2L], 'locations' = l)
      tmp <- list('mat' = par, 'pat' = pat)
    }
    # attr(tmp, 'homozygous') <-  homozygous
    # attr(tmp, 'L') <-  Lj
    ret[[j]] <- tmp
  }
  attr(x = ret, which = 'homozygous') <- homozygous
  attr(x = ret, which = 'individual') <- TRUE
  ret
}
