cat("\nBenchmarking greedy fits...\n\n")

if (exists("test") && test) {
  Kmax <- 2
} else {
  Kmax <- 10
}

parse.timing.res <- function(interval, is.flashr) {
  if (is.flashr) {
    init.str <- "add_factors_from_data"
    wrapup.str <- "construct_flash_object"
  } else {
    init.str <- "init.factor"
    wrapup.str <- "wrapup.flash"
  }

  Rprof.res <- read_lines("Rprof.out")
  Rprof.res <- Rprof.res[-1]

  is.wrapup <- str_detect(Rprof.res, wrapup.str)
  Rprof.res <- Rprof.res[!is.wrapup]

  is.init <- str_detect(Rprof.res, init.str)
  init.chg <- which(is.init[1:(length(is.init) - 1)] != is.init[2:length(is.init)])

  all.chg <- c(init.chg, length(Rprof.res)) * interval

  if (is.flashr) {
    fl.res <- read_lines("Output.out")

    tbl.begins <- which(str_detect(fl.res, "Iteration")) + 1
    tbl.ends <- c(tbl.begins[2:length(tbl.begins)] - 3, length(fl.res))

    tib <- tibble()
    for (k in 1:length(tbl.begins)) {
      next.tbl <- fl.res[tbl.begins[k]:tbl.ends[k]]
      next.tbl <- lapply(next.tbl, str_trim)
      next.tbl <- sapply(next.tbl, str_split, " +")
      next.tbl <- do.call(rbind, next.tbl)
      next.tbl <- next.tbl[, -3]
      colnames(next.tbl) <- c("Iteration", "Obj")
      next.tbl <- as_tibble(next.tbl)
      next.tbl <- next.tbl %>%
        mutate_all(as.numeric) %>%
        mutate(Factor = k)
      first.t <- all.chg[2 * k]
      last.t <- all.chg[2 * k + 1]
      next.tbl <- next.tbl %>%
        mutate(Timestamp = seq(first.t, last.t, length.out = nrow(next.tbl)))

      tib <- tib %>%
        bind_rows(next.tbl)
    }
  } else {
    tib <- read_tsv("Output.out", show_col_types = FALSE)

    tib <- tib %>%
      select(Iter, Obj, Factor) %>%
      rename(Iteration = Iter)

    iters <- tib %>% group_by(Factor) %>% summarize(n = n()) %>% pull(n)
    all.chg <- all.chg[-1]
    ts <- numeric(0)
    for (k in 1:((length(all.chg)/2))) {
      ts <- c(ts, seq(all.chg[2 * k - 1], all.chg[2 * k], length.out = iters[k]))
    }
    tib <- tib %>%
      mutate(Timestamp = ts)
  }

  return(tib)
}

do.timing <- function(fl.expr, dat, interval, Kmax, method.name, is.flashr = FALSE) {
  tol <- sqrt(.Machine$double.eps) * prod(dim(dat))

  zz <- file("Output.out", open = "wt")
  sink(zz)
  sink(zz, type = "message")

  # flashier can provide timestamps without profiling, but since profiling
  #   itself adds time, I need to profile both to do a fair comparison.
  Rprof("Rprof.out", interval = interval)
  eval(fl.expr)
  Rprof(NULL)

  sink(NULL, type = "message")
  sink(NULL)
  close(zz)

  tib <- parse.timing.res(interval, is.flashr = is.flashr)
  tib <- tib %>% mutate(Method = method.name)

  zz <- file.remove("Rprof.out")
  zz <- file.remove("Output.out")

  return(tib)
}

do.all.timings <- function(dat, interval, Kmax) {
  res <- do.timing(
    quote(
      fl <- flashr::flash(
        dat,
        Kmax = Kmax,
        var_type = "constant",
        backfit = FALSE,
        nullcheck = FALSE,
        tol = tol
      )
    ),
    dat = dat,
    interval = interval,
    Kmax = Kmax,
    method.name = "flashr",
    is.flashr = TRUE
  )

  res <- res %>% bind_rows(do.timing(
    substitute(
      fl <- flash.init(dat) %>%
        flash.set.verbose(-1) %>%
        flash.add.greedy(Kmax = Kmax, init.fn = init.fn.softImpute, tol = tol)
    ),
    dat = dat,
    interval = interval,
    Kmax = Kmax,
    method.name = "efficient_ops"
  ))

  res <- res %>% bind_rows(do.timing(
    substitute(
      fl <- flash.init(dat) %>%
        flash.set.verbose(-1) %>%
        flash.add.greedy(Kmax = Kmax, tol = tol)
    ),
    dat = dat,
    interval = interval,
    Kmax = Kmax,
    method.name = "fast_init"
  ))

  if (mean(dat == 0) > 0.5) {
    dat <- Matrix(dat)

    res <- res %>% bind_rows(do.timing(
      substitute(
        fl <- flash.init(dat) %>%
          flash.set.verbose(-1) %>%
          flash.add.greedy(Kmax = Kmax, tol = tol)
      ),
      dat = dat,
      interval = interval,
      Kmax = Kmax,
      method.name = "sparse"
    ))
  }

  return(res)
}


cat("  GTEx data...\n")
dat <- readRDS("../../data/gtex.rds")
gtex <- do.all.timings(
  dat = dat,
  interval = .002,
  Kmax = Kmax
)

all.res <- gtex %>% add_column(Dataset = "GTEx")


if (!exists("test") || !test) {
  cat("  PBMCs data...\n")
  dat <- readRDS("../../data/pbmc.rds")
  dat <- as.matrix(log1p(dat))
  pbmc <- do.all.timings(
    dat = dat,
    interval = .01,
    Kmax = Kmax
  )

  cat("  Montoro data...\n")
  dat <- readRDS("../../data/trachea.rds")
  dat <- log1p(dat)
  trachea <- do.all.timings(
    dat = dat,
    interval = .01,
    Kmax = Kmax
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


saveRDS(all.res, "../../output/res_greedy.rds")
