cat("\nBenchmarking backfits...\n\n")

if (exists("test") && test) {
  K <- 2
} else {
  K <- 10
}

parse.timing.res <- function(is.flashr, elapsed.t) {
  if (is.flashr) {
    tib <- read_table("Output.out", skip = 2, col_names = c(
      "Iteration", "Objective", "Obj Diff"
    ), col_types = "ddc") %>%
      mutate(Timestamp = elapsed.t * Iteration / max(Iteration)) %>%
      rename(Obj = Objective) %>%
      select(Iteration, Obj, Timestamp) %>%
      mutate(Timestamp = Timestamp - min(Timestamp))
  } else {
    tib <- read_tsv("Output.out", col_names = c(
      "Type", "Factor", "Iteration", "Obj", "Timestamp"
    ), show_col_types = FALSE)

    tib <- tib %>%
      select(Iteration, Obj, Timestamp) %>%
      mutate(Timestamp = Timestamp - min(Timestamp))
  }

  return(tib)
}

do.timing <- function(fl.expr, dat, K, method.name, is.flashr = FALSE) {
  maxiter <- 100
  tol <- sqrt(.Machine$double.eps) * prod(dim(dat))
  svd.init <- irlba::irlba(dat, nv = K, nu = K)

  flashr.init <- function(Y, K) {
    return(svd.init)
  }

  fns <- c(function(new, old, k) {new$obj},
           function(new, old, k) {Sys.time()})
  colnames <- c("ELBO", "Timestamp")
  colwidths <- c(18, 36)

  zz <- file("Output.out", open = "wt")
  sink(zz)
  sink(zz, type = "message")

  elapsed.t <- system.time(eval(fl.expr))[3]

  sink(NULL, type = "message")
  sink(NULL)
  close(zz)

  tib <- parse.timing.res(is.flashr = is.flashr, elapsed.t = elapsed.t)
  tib <- tib %>% mutate(Method = method.name)

  zz <- file.remove("Output.out")

  return(tib)
}

do.all.timings <- function(dat, K) {
  res <- do.timing(
    quote(
      fl <- flashr::flash_add_factors_from_data(
        dat,
        K = K,
        var_type = "constant",
        backfit = TRUE,
        nullcheck = FALSE,
        tol = tol,
        maxiter = maxiter
      )
    ),
    dat = dat,
    K = K,
    method.name = "flashr",
    is.flashr = TRUE
  )

  res <- res %>% bind_rows(do.timing(
    substitute(
      fl <- flash.init(dat) %>%
        flash.set.verbose(-1, fns = fns, colnames = colnames, colwidths = colwidths) %>%
        flash.init.factors(svd.init) %>%
        flash.backfit(tol = tol / K, maxiter = maxiter, method = "sequential")
    ),
    dat = dat,
    K = K,
    method.name = "efficient_ops"
  ))

  res <- res %>% bind_rows(do.timing(
    substitute(
      fl <- flash.init(dat) %>%
        flash.set.verbose(-1, fns = fns, colnames = colnames, colwidths = colwidths) %>%
        flash.init.factors(svd.init) %>%
        flash.backfit(tol = tol, maxiter = maxiter, method = "extrapolate")
    ),
    dat = dat,
    K = K,
    method.name = "extrapolation"
  ))

  if (mean(dat == 0) > 0.5) {
    dat <- Matrix(dat)

    res <- res %>% bind_rows(do.timing(
      substitute(
        fl <- flash.init(dat) %>%
          flash.set.verbose(-1, fns = fns, colnames = colnames, colwidths = colwidths) %>%
          flash.init.factors(svd.init) %>%
          flash.backfit(tol = tol, maxiter = maxiter, method = "extrapolate")
      ),
      dat = dat,
      K = K,
      method.name = "sparse"
    ))
  }

  return(res)
}


cat("  GTEx data...\n")
dat <- readRDS("../../data/gtex.rds")
gtex <- do.all.timings(
  dat = dat,
  K = K
)

all.res <- gtex %>% add_column(Dataset = "GTEx")


if (!exists("test") || !test) {
  cat("  PBMCs data...\n")
  dat <- readRDS("../../data/pbmc.rds")
  dat <- as.matrix(log1p(dat))
  pbmc <- do.all.timings(
    dat = dat,
    K = K
  )

  cat("  Montoro data...\n")
  dat <- readRDS("../../data/trachea.rds")
  dat <- log1p(dat)
  trachea <- do.all.timings(
    dat = dat,
    K = K
  )

  all.res <- all.res %>%
    bind_rows(
      pbmc %>% add_column(Dataset = "PBMC")
    ) %>%
    bind_rows(
      trachea %>% add_column(Dataset = "Montoro")
    )
}

rm(dat)


saveRDS(all.res, "../../output/res_bf.rds")
