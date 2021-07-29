cat("Running pseudocount experiments...\n\n")
if (exists("test") && test) {
  Kmax <- 1
} else {
  Kmax <- 20
}

# Assumes rows are cells and cols are genes:
preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- rowSums(dat)
  size.factors <- size.factors / mean(size.factors)
  dat <- dat / size.factors
  gene_cts <- colSums(dat > 0)
  dat <- dat[, gene_cts >= min.nzcts]
  return(list(dat = dat, size.factors = size.factors, excluded.genes = gene_cts < min.nzcts))
}

adj.elbo <- function(dat, fl, pseudocount) {
  return(fl$elbo - prod(dim(dat)) * log(pseudocount) - sum(log1p(dat / pseudocount)))
}

llik.idd <- function(dat, fl, size.factors, pseudocount) {
  upper.cdf <- log1p((dat / pseudocount) / (1 + 1 / (2 * pseudocount * size.factors)))
  upper.cdf <- upper.cdf + log1p(1 / (2 * pseudocount * size.factors))
  upper.cdf <- upper.cdf - fitted(fl)
  upper.cdf <- t(t(upper.cdf) / fl$residuals.sd)

  zeros.llik <- pnorm(upper.cdf[which(dat == 0)], log.p = TRUE)

  nz.upper.cdf <- upper.cdf[which(dat > 0)]
  nz.upper.cdf.sign <- sign(nz.upper.cdf)
  nz.upper.cdf.log <- pnorm(as.matrix(-abs(nz.upper.cdf)), log.p = TRUE)

  nz.lower.cdf <- log(dat@x / pseudocount - 0.5 / (pseudocount * size.factors[dat@i + 1]) + 1)
  nz.lower.cdf <- nz.lower.cdf - fitted(fl)[which(dat > 0)]
  nz.lower.cdf <- nz.lower.cdf / fl$residuals.sd[findInterval(seq(dat@x) - 1, dat@p[-1]) + 1]
  nz.lower.cdf.sign <- sign(nz.lower.cdf)
  nz.lower.cdf.log <- pnorm(-abs(nz.lower.cdf), log.p = TRUE)

  diff.signs <- (nz.upper.cdf.sign > 0 & nz.lower.cdf.sign < 0)
  diff.signs.llik <- log(1 - exp(nz.upper.cdf.log[diff.signs]) - exp(nz.lower.cdf.log[diff.signs]))

  neg.signs <- (nz.upper.cdf.sign < 0 & nz.lower.cdf.sign < 0)
  neg.signs.llik <- nz.upper.cdf.log[neg.signs] +
    log(1 + exp(nz.lower.cdf.log[neg.signs] - nz.upper.cdf.log[neg.signs]))

  pos.signs <- (nz.upper.cdf.sign > 0 & nz.lower.cdf.sign > 0)
  pos.signs.llik <- nz.lower.cdf.log[pos.signs] +
    log(1 + exp(nz.upper.cdf.log[pos.signs] - nz.lower.cdf.log[pos.signs]))

  all.llik <- sum(zeros.llik) + sum(diff.signs.llik) + sum(neg.signs.llik) + sum(pos.signs.llik)
}

do.fit <- function(dat, size.factors, pseudocount, Kmax = 30) {
  cat("  Fitting Pseudocount:", as.character(pseudocount), "\n")

  fl <- flash.init(log1p(dat / pseudocount), var.type = 2) %>%
    flash.add.greedy(Kmax = Kmax, verbose.lvl = 1, extrapolate = FALSE) # %>%
    # flash.backfit(verbose.lvl = 3)

  fl$adj.elbo <- adj.elbo(dat, fl, pseudocount)

  cat("  Calculating IDD log likelihood...\n")

  fl$llik.idd <- llik.idd(dat, fl, size.factors, pseudocount)
  fl$elbo.idd <- fl$llik.idd + sum(fl$flash.fit$KL[[1]]) + sum(fl$flash.fit$KL[[2]])

  fl$flash.fit <- NULL
  fl$sampler <- NULL

  return(fl)
}


pseudocounts <- 2^seq(-4, 6, by = 0.5)


pbmc <- t(readRDS("../../data/pbmc.rds"))
pbmc.pp <- preprocess(pbmc)
pbmc.res <- list()
for (pc in pseudocounts) {
  pc.str <- as.character(pc)
  pbmc.res[[pc.str]] <- do.fit(pbmc.pp$dat, pbmc.pp$size.factors, pc, Kmax = Kmax)
}

rm(pbmc)
rm(pbmc.pp)

saveRDS(pbmc.res, "../../output/pseudocount_pbmc_greedy.rds")


trachea <- t(Matrix(readRDS("../../data/trachea.rds")))
trachea.pp <- preprocess(trachea)
trachea.res <- list()
for (pc in pseudocounts) {
  pc.str <- as.character(pc)
  trachea.res[[pc.str]] <- do.fit(trachea.pp$dat, trachea.pp$size.factors, pc, Kmax = Kmax)
}

rm(trachea)
rm(trachea.pp)

saveRDS(trachea.res, "../../output/pseudocount_trachea_greedy.rds")
