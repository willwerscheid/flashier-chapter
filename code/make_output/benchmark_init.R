if (exists("test") && test) {
  Kmax <- 2
} else {
  Kmax <- 10
}

init.fn.irlba <- function(flash, mode.signs) {
  R <- flash$Y - flashier:::lowrank.expand(flash$EF)
  res <- irlba::irlba(R, nv = 1, nu = 1, maxit = 100)
  return(list(as.vector(res$u * sqrt(res$d)), as.vector(res$v * sqrt(res$d))))
}

parse.timing.res <- function(interval) {
  init.str <- "init.factor"
  wrapup.str <- "wrapup.flash"

  Rprof.res <- read_lines("Rprof.out")
  Rprof.res <- Rprof.res[-1]

  is.wrapup <- str_detect(Rprof.res, wrapup.str)
  Rprof.res <- Rprof.res[!is.wrapup]

  is.init <- str_detect(Rprof.res, init.str)
  init.chg <- which(is.init[1:(length(is.init) - 1)] != is.init[2:length(is.init)])
  all.chg <- c(init.chg, length(Rprof.res)) * interval
  all.t <- all.chg[2:length(all.chg)] - all.chg[1:(length(all.chg) - 1)]

  K <- length(all.t) / 2
  tib <- tibble(Factor = 1:K,
                InitTime = all.t[2 * (1:K) - 1],
                RefineTime = all.t[2 * (1:K)]) %>%
    mutate(TotalTime = InitTime + RefineTime)

  return(tib)
}

do.timing <- function(init.fn, dat, interval, Kmax, method.name, any.na = FALSE) {
  cat(method.name, "\n")

  tol <- sqrt(.Machine$double.eps) * prod(dim(dat))
  maxiter <- 100

  Rprof("Rprof.out", interval = interval)
  fl <- flash.init(dat) %>%
    flash.set.verbose(0) %>%
    flash.add.greedy(Kmax = Kmax, extrapolate = FALSE, init.fn = init.fn)
  Rprof(NULL)

  tib <- parse.timing.res(interval)
  tib <- tib %>%
    mutate(Method = method.name,
           AnyNA = any.na)

  zz <- file.remove("Rprof.out")

  return(tib)
}

do.all.timings <- function(dat, interval, Kmax) {
  res <- do.timing(init.fn.default, dat, interval, Kmax, "default") %>%
    bind_rows(do.timing(init.fn.softImpute, dat, interval, Kmax, "si")) %>%
    bind_rows(do.timing(init.fn.irlba, dat, interval, Kmax, "irlba"))

  dat <- Matrix(dat)
  res <- res %>%
    bind_rows(do.timing(init.fn.default, dat, interval, Kmax, "sparse"))

  dat <- as.matrix(dat)
  dat <- dat[rowSums(dat > 0) > 5, ]

  set.seed(666)
  is.missing <- sample(1:length(dat), ceiling(.2 * length(dat)))

  dat <- as.matrix(dat)
  dat[is.missing] <- NA

  res <- res %>%
    bind_rows(do.timing(init.fn.default, dat, interval, Kmax, "default", TRUE)) %>%
    bind_rows(do.timing(init.fn.softImpute, dat, interval, Kmax, "si", TRUE))

  dat <- Matrix(dat)
  res <- res %>%
    bind_rows(do.timing(init.fn.default, dat, interval, Kmax, "sparse", TRUE))

  return(res)
}


dat <- readRDS("../../data/pbmc.rds")
dat <- as.matrix(log1p(dat))
pbmc <- do.all.timings(
  dat = dat,
  interval = .025,
  Kmax = Kmax
)

dat <- readRDS("../../data/trachea.rds")
dat <- log1p(dat)
trachea <- do.all.timings(
  dat = dat,
  interval = .025,
  Kmax = Kmax
)

rm(dat)


all.res <- pbmc %>% add_column(Dataset = "PBMC") %>%
  bind_rows(
    trachea %>% add_column(Dataset = "Montoro")
  )

saveRDS(all.res, "../../output/res_init.rds")
