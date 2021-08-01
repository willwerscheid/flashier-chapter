if (exists("test") && test) {
  ntrials <- 2
  ns <- c(100, 500, 1000)
  k <- 2
} else {
  ntrials <- 10
  ns <- c(100, 200, 500, 1000, 2000, 5000, 10000)
  k <- 10
}
cat("\nRunning flashier with misspecified noise models:", ntrials, "trials...\n\n")


# simulation hyperparameters

pif_prior_a <- 0.5
pif_prior_b <- 0.5
pil_prior_a <- 0.5
pil_prior_b <- 0.5

s2f_prior_shape <- 1
s2f_prior_rate <- 2
s2l_prior_shape <- 1
s2l_prior_rate <- 2

s2_resid_sd <- 0.5


# set up results arrays

nfactors <- matrix(nrow = ntrials, ncol = length(ns))
nfactors <- rep(list(nfactors), 4)


# run simulations

for (trial in 1:ntrials) {
  cat("  Trial:", trial, "\n")
  set.seed(trial)

  pif <- rbeta(k, pif_prior_a, pif_prior_b)
  pil <- rbeta(k, pil_prior_a, pil_prior_b)
  s2f <- rgamma(k, s2f_prior_shape, s2f_prior_rate)
  s2l <- rgamma(k, s2l_prior_shape, s2l_prior_rate)

  for (n in ns) {
    cat("    n:", n, "\n")
    p <- n

    LL <- matrix(nrow = n, ncol = k)
    FF <- matrix(nrow = p, ncol = k)

    # simulate point-normal factors and loadings
    for (i in 1:k) {
      LL[, i] <- rnorm(n, sd = sqrt(s2l[i])) * rbinom(n, 1, 1 - pil[i])
      FF[, i] <- rnorm(n, sd = sqrt(s2f[i])) * rbinom(n, 1, 1 - pif[i])
    }

    Y_norm <- tcrossprod(LL, FF) + s2_resid_sd * rnorm(n * p)
    fl_norm <- flash.init(Y_norm) %>%
      flash.add.greedy(Kmax = 100, extrapolate = FALSE, verbose.lvl = 0)
    nfactors[[1]][trial, which(ns == n)] <- fl_norm$n.factors

    Y_t5 <- tcrossprod(LL, FF) + s2_resid_sd * rt(n * p, df = 5)
    fl_t5 <- flash.init(Y_t5) %>%
      flash.add.greedy(Kmax = 100, extrapolate = FALSE, verbose.lvl = 0)
    nfactors[[2]][trial, which(ns == n)] <- fl_t5$n.factors

    Y_pois <- matrix(log1p(rpois(n * p, exp(tcrossprod(LL, FF)))), nrow = n, ncol = p)
    fl_pois <- flash.init(Y_pois) %>%
      flash.add.greedy(Kmax = 100, extrapolate = FALSE, verbose.lvl = 0)
    nfactors[[3]][trial, which(ns == n)] <- fl_pois$n.factors

    Y_poisln <- matrix(log1p(rpois(n * p, exp(Y_norm))), nrow = n, ncol = p)
    fl_poisln <- flash.init(Y_poisln) %>%
      flash.add.greedy(Kmax = 100, extrapolate = FALSE, verbose.lvl = 0)
    nfactors[[4]][trial, which(ns == n)] <- fl_poisln$n.factors
  }
}


# process results

res <- lapply(nfactors, function(mat) {
  gather(data.frame(mat), key = "n", value = "nfactors")
})

distns <- c("normal", "t5", "poisson", "poislognorm")
for (i in 1:length(distns)) {
  res[[i]]$distn <- distns[i]
}

res <- bind_rows(res)

res$n <- factor(res$n)
levels(res$n) <- ns

res$distn <- factor(res$distn,  levels = distns)

saveRDS(res, "../../output/misspec_noise.rds")
