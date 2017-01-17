library('simbreed')

simbreed::gamete(parent = individual, xo_params)

n_chr <- 10L
m <- rep(100, n_chr) 
len <- rep(200, n_chr)
n_founder <- 2L

map <- .sim_map(m = m , len = len)

(founder <- sim_founder(n = n_founder, m = m))
parent <- create_parent(len = len, alleles = c(0L, 1L), homozygous = FALSE)
(eff <- purrr::map(m, rnorm))
(xo_params <- create_xo_params(L = len, m = 0, p = 1.0))
(xodat <- gamete(parent = parent, xo_params = xo_params))



gam <- simbreed::.xo2geno_gamete(xodat = xodat, map = map, founder = founder)
u <- cbind(do.call(rbind, founder),t(gam))
colMeans(u[u[,1]+u[,2] ==1,])

individual <- parent 
for (i in 1:5 ){
  individual <- self(individual, xo_params)
}
individual


simbreed::.bcgv(parent = individual, n_gam = 10, n_rep = 10, se_thresh = 0.05,
      params = xo_params, map = map, founder = founder, eff = eff)

library(simbreed)
devtools::unload()
devtools::uninstall()
devtools::install()
devtools::load_all()

bcgv <- function(ind, n_gam, n_rep, xo_params, map, founder, effects) {
  cgv <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    tmp <- numeric(n_gam)
    for (i in seq_len(n_gam)) {
      ## xodat <- gamete(parent = ind, xo_params = xo_params)
      ## tmp[i] <- .gamete_value(xodat = xodat, map = map, founder = founder, eff = eff)
      tmp[i] <- .gamete_value2(parent = ind, params = xo_params,
                               map = map, founder= founder, eff = eff)
    }
    cgv[r] <- max(tmp)
  }
  list('bcgv' = mean(cgv), 'se' = sd(cgv) / sqrt(n_rep))
}


microbenchmark::microbenchmark(times = 10,
bcgv(individual, 10, 100, xo_params, map, founder, eff),
.bcgv(individual, 10, 100, xo_params, map, founder, eff)
)

bcgv(individual, 10, 1000, xo_params, map, founder, eff)
.bcgv(individual, 10, 10, xo_params, map, founder, eff)

confi <- function(bcgv) {
  m <- bcgv$bcgv
  se <- bcgv$se
  c(m - 1.96*se, m + 1.96*se)
}

confi(.bcgv(individual, 10, 10, xo_params, map, founder, eff))

microbenchmark::microbenchmark(times = 1,
x <- purrr::flatten_dbl(purrr::rerun(1000, max(ohv(n = 10, individual, xo_params, eff)))),
x <- purrr::flatten_dbl(purrr::rerun(1000, max(ohv2(n = 10, individual, xo_params, eff))))
## x <- purrr::flatten_dbl(purrr::rerun(100, max(ohv(n = 1, individual, xo_params, eff))))
)

quantile(x)
(mu <- mean(x));
sdv <- sd(x) / sqrt(length(x))
pnorm(mu-0.01, mean = mu, sd = sdv, lower.tail = TRUE)
qnorm(0.01, mean = mu, sd = sdv, lower.tail = TRUE)
plot(density(x))



microbenchmark::microbenchmark(times = 100,
ohv1 = sum(purrr::pmap_dbl(.l = list(xodat, map, founder, eff),
                .f = function(x, m, f, e) {
                  .xo2geno_chromatid(xodat = x, map = m, founder = f) %*% e
                })),
ohv2 = sum(purrr::pmap_dbl(.l = list(xodat, map, founder, eff),
                .f = function(x, m, f, e) {
                  .chromatid_value(xodat = x, map = m, founder = f, eff = e)
                })),
ohv = .xo2geno_gamete(xodat = xodat, map = map, founder = founder) %*% unlist(eff),
ohv4 = .gamete_value(xodat, map, founder, eff)
)

individual <- parent


library('simbreed')
library('hypred')
library('magrittr')




nXO <- 0L
n_founder <- 2L
n_chr <- 10L
len <- runif(n = n_chr, min = 200, max = 200)
n_loci <- seq(from = log(5), to = log(1e5), length.out = 15L) %>% exp() %>% ceiling()
# Benchmarking an comparing the simulation of meiosis in hypred and simbreed.
ret <- purrr::map_df(n_loci, function(m) {
  map <- .sim_map(m = m, len = len)
  genome <- hypred::hypredGenome(num.chr = n_chr, len.chr = len / 100.0, num.snp.chr = m)
  # parent_xodat <- create_parent(len = len, alleles = seq.int(1L, n_founder) - 1L, homozygous = FALSE)
  parent_xodat <- .sim_xodat_individual(nXO, len, c(0L, 1L), FALSE)
  founder <- sim_founder(n = n_founder, m = rep(m, n_chr))

  parent <- .xo2geno_individual(parent_xodat, map, founder)

  bench_hypred <- microbenchmark::microbenchmark(times = 5L,
    hypred = hypred::hypredRecombine(genome, genomeA = parent[1L, ], genomeB = parent[2L, ],
                                     mutate = FALSE, block = FALSE),
    simbreed = for (j in seq_len(n_chr)) {
      tmp <- simbreed::.sim_meiosis(parent = parent_xodat[[j]], L = len[j], m = 0L, p = 1.0,
                                   obligate_chiasma = FALSE, Lstar = len[j])
      simbreed:::.xo2geno_chromatid(xodat = tmp, map = map[[j]], founder = founder[[j]])
    }
    # },
    #
    # simcross = for (j in seq_len(n_chr)) {
    #   tmp <- simcross::sim_meiosis(parent = parent_xodat[[j]], m = 0L, p = 1.0,
    #                                obligate_chiasma = FALSE, Lstar = len[j])
    #   simbreed:::.xo2geno_chromatid(xodat = tmp, map = map[[j]], founder = founder[[j]])
    # }
  )
  ret <- as.data.frame(bench_hypred)
  ret$m = m
  ret$time <- ret$time / 1e6
  ret
})


head(ret)

library('ggplot2')
plt <- ggplot(data = ret, mapping = aes(x = m, y = time, color = expr)) +
  stat_summary(geom = 'point', fun.y = 'mean') +
  scale_x_log10(breaks = n_loci) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  # scale_y_log10() +
  # scale_y_log10(breaks = scales::pretty_breaks(n = 10L)) +
  labs(y = 'execution time in milliseconds', x = 'number of loci per chromosome') #+
  # coord_cartesian(ylim = c(0.001, max(ret$time)))

ggsave('/tmp/mytmpfile.pdf')


# Test conversion
library('simbreed')
nXO <- 10L
m <- 10000L

xodat <- sim_xodat_chromatid(nXO)
map <- sim_map(m = m)
founder <- sim_founder(n = 2L, m = m)
t_founder <- t(founder)

microbenchmark::microbenchmark(times = 1e5L,
                               a = .xo2geno_chromatid(xodat, map, founder)
                               # b = xo2geno_chromatid_dev(xodat, map, t_founder),
                               # c = xo2geno_chromatid_dev2(xodat, map, t_founder)
)
a[1:5]
a = .xo2geno_chromatid(xodat, map, founder)
b = xo2geno_chromatid_dev(xodat, map, t_founder)
c = xo2geno_chromatid_dev2(xodat, map, t_founder)
sum((a-b)^2)
sum((a-c)^2)

# Convert chromosome
library('simbreed')
n_chr <- 10L
nXO <- 10L
m <- 10000L

map <- sim_map(m = m)
xodat <- purrr::rerun(2L, sim_xodat_chromatid(nXO))
founder <- t(sim_founder(n = 2L, m = m))

microbenchmark::microbenchmark(times = 1e4L,
.xo2geno_chromatid(xodat[[1]], map, founder),
.xo2geno_chromosome(xodat = xodat, map = map, founder = founder, FALSE),
.xo2geno_chromosome_vert(xodat = xodat, map = map, founder = founder, FALSE)
)

u<-.xo2geno_chromosome_vert(xodat = xodat, map = map, founder = founder, FALSE)
dim(u)


# Convert entire indivdiual/gamete
library('simbreed')
n_chr <- 10L
nXO <- 1L
m <- 100L

map <- purrr::rerun(n_chr, sim_map(m = m))
xodat <- purrr::rerun(n_chr, purrr::rerun(2L, sim_xodat_chromatid(nXO)))
founder <- purrr::rerun(n_chr, t(sim_founder(n = 2L, m = m)))
attr(xodat, 'homozygous') = FALSE
xodat_gam <- purrr::rerun(n_chr, sim_xodat_chromatid(nXO))
attr(xodat_gam, 'homozygous') = TRUE
population <- purrr::rerun(1000, xodat)

microbenchmark::microbenchmark(times = 1e3L,
.xo2geno_chromatid(xodat[[1]][[1]], map[[1]], founder[[1]]),
.xo2geno_chromosome(xodat[[1]], map[[1]], founder[[1]], FALSE),
.xo2geno_individual(xodat, map, founder),
.xo2geno_gamete(xodat_gam, map, founder),
.xo2geno_individual(xodat_gam, map, founder)
)

microbenchmark::microbenchmark(times = 1e2L,
.xo2geno_population(population, map, founder),
purrr:::rerun(length(population), .xo2geno_individual(xodat, map, founder)) %>%
  do.call(what=rbind)
)

microbenchmark::microbenchmark(times = 1e4L,
purrr:::pmap(list(xodat_gam, map, founder), .xo2geno_chromatid) %>%
  do.call(what=cbind)
)

.xo2geno_chromatid(xodat_gam[[1]], map[[1]], founder[[1]])



 u<-.xo2geno_individual(xodat, map, founder)


                               t(u)
)

.xo2geno_chromosome(xodat, map[[1]], founder[[1]], FALSE)


# Test simulation of meiosis
library('LDtools')
nXO <- 10000L
L <- 3000.0
mat <- sim_xodat_chromatid(nXO, max = L)
pat <- sim_xodat_chromatid(nXO, max = L)
parent <- list('mat' = mat, 'pat' = pat)

# Xlocations = .sim_crossovers(L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L)
# Xlocations <- seq(100, L - 100, length.out = 2L)
L <- 500
microbenchmark::microbenchmark(times = 1e4L,
# sim_meiosis_xo(parent, .sim_crossovers(L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L)),
# sim_meiosis_xo_test(parent, .sim_crossovers(L, m = 0L, p = 1.0, obligate_chiasma = FALSE, Lstar = L)),
# sim_meiosis_bro(parent = parent, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
# sim_meiosis_reliable(parent = parent, L = L, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
# sim_meiosis(parent = parent, L = L, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
# sim_meiosis_while(parent = parent, L = L, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
.sim_crossovers(L, m = 10L, p = 0.0, obligate_chiasma = FALSE, Lstar = L),
simcross::sim_crossovers(L, m = 10L, p = 0.0, obligate_chiasma = FALSE, Lstar = L)
)


.sim_crossovers_arma(L, m = 10L, p = 0.5, obligate_chiasma = FALSE, Lstar = L)


ct <- 0L
for (i in 1:1000) {
a <- sim_meiosis_xo(parent, Xlocations)
b <- sim_meiosis_xo_test(parent, Xlocations)
ct <- ct + identical(a, b)
}
ct/1000


head(a$locations)
head(b$locations)


microbenchmark::microbenchmark(times = 1e4L,
a =  sim_meiosis_reliable(parent = parent, L = L, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
b =  sim_meiosis(parent = parent, L = L, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L),
c =  sim_meiosis_bro(parent = parent, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = L)
)






str(parent)
str(F1[[1]])

sim_meiosis(parent = parent, L = 100, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = 100)
sim_meiosis(parent = F1[[1]], L = 100, m = 0, p = 1, obligate_chiasma = FALSE, Lstar = 100)


# Test speed of simulation of families/populations.
library('simbreed')
n_chr <- 10L
n_snp <- 10000L
L <- rep(200, n_chr)
xo_params <- create_xo_params(L, m = 0L, p = 1.0)
parents <- create_parents(n = 2, L = L , homozygous = T)
pryr::object_size(parents)

map <- replicate(n_chr, simplify = FALSE,
                 expr = sort(runif(n_snp, 0.0, 100.0)))

founder <- replicate(n = length(L), simplify = FALSE,
                     expr = matrix(sample(c(0, 1),
                                          size = 2* n_snp, TRUE), ncol = n_snp))
size_fam <- 1000L
n_gen <- 30L

microbenchmark::microbenchmark(times = 1, {
f1 = cross(parents[[1]], parents[[2]], xo_params)
fam = replicate(size_fam, simplify = FALSE, self(f1, xo_params))
out = purrr::map(fam, ~xo2geno(xodat = .x, map = map, founder = founder))
})

pryr::object_size(out)

sample_one <- function(x) sample(x, size = 1L)[[1L]]

# L <- params$L
# m <- params$m
# p <- params$p
# obligate_chiasma <- params$obligate_chiasma
# Lstar <- params$Lstar

microbenchmark::microbenchmark(times = 5, j={
  f1 <- cross(parents[[1]], parents[[2]], xo_params)
  fn <- replicate(size_fam, simplify = FALSE, self(f1, xo_params))
  for (i in seq_len(n_gen)) {
    ftmp <- fn
    for (fam in seq_len(size_fam)) {
      fn[[i]] <- cross(sample_one(ftmp), sample_one(ftmp), xo_params)
      # fn[[i]] <- .cross(sample_one(ftmp), sample_one(ftmp), L, m, p, obligate_chiasma, Lstar)
    }
  }
  # for (i in seq_len(n_gen))
  #   fn <- lapply(seq_len(size_fam), function(fam) cross(sample_one(fn), sample_one(fn), params))
    # fn <- purrr::map(seq_len(size_fam), ~ cross(sample_one(fn), sample_one(fn), params))
})

fn[[1]]




attributes(fam[[1]])

microbenchmark::microbenchmark(times = 10, bla = {
  purrr::map(fam, ~xo2geno(xodat = .x, map = map, founder = founder))
})



microbenchmark::microbenchmark(times = 1e4,
.cross(parents[[1]], parents[[2]], params),
dh(parents[[1]], params),
self(parents[[1]], params),
.gamete(parents[[1]], params)
)

