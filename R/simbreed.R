#' A package for simulation of populations in plant breeding.
#'
#' This package provides flexible tools for simulating any kind of population encountered
#' in plant breeding, like for instance bi-parental populations, synthethic populations, back-
#' crosses, recombinant inbred lines etc.
#'
#'
#' @section Functions:
#' There are three kinds of central functions:
#'
#' \enumerate{
#' \item \code{\link{create_parents}}: Create parents as a starting point for the simulation of
#' populations.
#' \item \code{\link{cross}}, \code{\link{self}}, \code{\link{dh}}: Cross two individuals, self
#' a single individual or produce a doubled haploid from an individual.
#' \item \code{\link{convert}} Convert the crossover data into genotypic data.
#' }
#'
#' There are additional functions for the simulation of common pedigrees encountered in plant
#' breeding, the simulation of crossover data from these pedigrees and the analysis of crossover
#' data.
#'
#' @section Workflow:
#'
#' There are two kinds of basic workflows. The first workflow is centered around building the
#' population "manually" by using the function \item\code{\link{create_parents}} for obtaining
#' population parents and the functions \code{\link{cross}}, \code{\link{self}}, \code{\link{dh}}
#' for successively deriving progenys. The second workflow first constructs a pedigree and
#' then simulates the population accordingly. Note that only the second workflow allows for
#' computation of pedigree relationships between individuals.
#'
#' Please see the vignette accompanying this package for more details and examples.
#'
#'
#' @section Limitations:
#'
#' \enumerate{
#' \item This software is only able to deal with diploid or allopolyploid species, where
#' exactly two homologous chromosomes pair up in a four-strand bundle during meiosis.
#' \item There is currently no convenient possibility to simulate mutations without manually
#' editing output and input data after each round of recombination.
#' }
#'
#' @section Acknowledgement:
#'
#' This package is built upon previous work by \href{http://kbroman.org/}{Karl Broman} on the
#' R-package \href{https://github.com/kbroman/simcross}{simcross}, from which concepts,
#' data structures and some C++ routines were adapted.
#'
#'
#'
#'
#' @author Dominik Mueller (\email{dominikmueller64@yahoo.de})
#'
#' @docType package
#' @name simbreed
#'
#' @useDynLib simbreed
#' @importFrom Rcpp sourceCpp
NULL
