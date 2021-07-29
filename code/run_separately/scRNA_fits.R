library(tidyverse)
library(Matrix)
library(flashier)
library(glmpca)
library(fastTopics)

# Assumes rows are genes and cols are cells:
preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  gene_cts <- rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]

  lunpc <- max(1 / min(size.factors) - 1 / max(size.factors), 1)
  fl.dat <- log1p(t(t(dat) / size.factors) / lunpc)

  return(list(
    dat = dat,
    fl.dat = fl.dat,
    size.factors = size.factors,
    excluded.genes = gene_cts < min.nzcts)
  )
}

init.snn.LL <- function(fl) {
  LL <- fl$flash.fit$EF[[1]]
  FF <- fl$flash.fit$EF[[2]]
  LL <- cbind(LL, -LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))

  to.remove <- (colSums(FF) < .Machine$double.eps)
  LL <- LL[, -to.remove, drop = FALSE]
  FF <- FF[, -to.remove, drop = FALSE]

  return(list(LL, FF))
}

do.fits <- function(pp, K, do.snn = TRUE) {
  dat <- pp$dat
  fl.dat <- pp$fl.dat
  sf <- pp$size.factors

  set.seed(666)

  all.t <- system.time({
    fit.flashier.pn <- flash.init(
      fl.dat,
      var.type = 1,
      var.reg.fn = ebpm_exponential
    ) %>% flash.add.greedy(
        Kmax = K,
        prior.family = prior.point.normal(),
        extrapolate = FALSE
      ) %>%
      flash.backfit(verbose.lvl = 3)
  })[3]
  LL <- fit.flashier.pn$loadings.pm[[2]]
  colnames(LL) <- paste0("k=", formatC(1:ncol(LL), width = 2, flag = "0"))
  res <- as_tibble(LL) %>% add_column(Method = "flashier.pn")

  if (do.snn) {
    all.t <- c(all.t, system.time({
      snn.init <- init.snn.LL(fit.flashier.pn)
      fit.flashier.snn <- flash.init(
        fl.dat,
        var.type = 1,
        var.reg.fn = ebpm_exponential
      ) %>% flash.init.factors(
        EF = snn.init,
        prior.family = c(prior.normal.scale.mix(), prior.nonnegative())
      ) %>% flash.backfit(
        verbose.lvl = 3,
        maxiter = 5
      )
      kset <- (length(fit.flashier.snn$pve) - rank(fit.flashier.snn$pve) < K)
      fit.flashier.snn <- flash.init(
        fl.dat,
        var.type = 1,
        var.reg.fn = ebpm_exponential
      ) %>% flash.init.factors(
        EF = lapply(fit.flashier.snn$flash.fit$EF, function(x) x[, kset]),
        EF2 = lapply(fit.flashier.snn$flash.fit$EF2, function(x) x[, kset]),
        prior.family = c(prior.normal.scale.mix(), prior.nonnegative())
      ) %>% flash.backfit(
        verbose.lvl = 3
      )
    })[3])
    LL <- fit.flashier.snn$loadings.pm[[2]]
    colnames(LL) <- paste0("k=", formatC(1:ncol(LL), width = 2, flag = "0"))
    res <- res %>%
      bind_rows(as_tibble(LL) %>% add_column(Method = "flashier.snn"))
  }

  all.t <- c(all.t, system.time({
    fit.fasttopics <- fit_topic_model(
      t(dat),
      k = K,
      verbose = "detailed"
    )
  })[3])
  LL <- fit.fasttopics$L
  colnames(LL) <- paste0("k=", formatC(1:ncol(LL), width = 2, flag = "0"))
  res <- res %>%
    bind_rows(as_tibble(LL) %>% add_column(Method = "topic.model"))

  all.t <- c(all.t, system.time({
    fit.glmpca.pois <- glmpca(
      dat,
      L = K,
      fam = "poi",
      sz = sf,
      ctl = list(verbose = TRUE)
    )
  })[3])
  LL <- fit.glmpca.pois$factors
  colnames(LL) <- paste0("k=", formatC(1:ncol(LL), width = 2, flag = "0"))
  res <- res %>%
    bind_rows(as_tibble(LL) %>% add_column(Method = "glmpca.pois"))

  all.t <- c(all.t, system.time({
    fit.glmpca.nb <- glmpca(
      dat,
      L = K,
      fam = "nb2",
      sz = sf,
      ctl = list(verbose = TRUE)
    )
  })[3])
  LL <- fit.glmpca.nb$factors
  colnames(LL) <- paste0("k=", formatC(1:ncol(LL), width = 2, flag = "0"))
  res <- res %>%
    bind_rows(as_tibble(LL) %>% add_column(Method = "glmpca.nb"))

  fit.flashier.snn$sampler <- NULL
  fit.flashier.snn$flash.fit <- NULL

  names(all.t) <- unique(res$Method)

  if (do.snn) {
    return(list(t = all.t, res = res, flashier.fit = fit.flashier.snn))
  } else {
    return(list(t = all.t, res = res, flashier.fit = fit.flashier.pn))
  }
}

pbmc <- readRDS("./data/pbmc.rds")
pbmc <- pbmc[, colSums(pbmc) < 15000]
pbmc.res <- do.fits(preprocess(pbmc), K = 12)
saveRDS(pbmc.res, "./output/pbmc_res.rds")

rm(pbmc)

montoro <- readRDS("./data/sp_trachea.rds")
montoro <- expm1(montoro)
montoro <- montoro[, colSums(montoro) < 30000]
montoro.res <- do.fits(preprocess(montoro), K = 20)
saveRDS(montoro.res, "./output/montoro_res.rds")
