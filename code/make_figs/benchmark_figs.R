plot.fits <- function(greedy.res, bf.res) {
  greedy.tib <- greedy.res %>%
    filter(Method != "extrapolation") %>%
    group_by(Dataset) %>%
    mutate(Obj = Obj - max(Obj)) %>%
    group_by(Dataset, Method, Factor) %>%
    summarize(Obj = max(Obj), Timestamp = max(Timestamp)) %>%
    ungroup() %>%
    select(-Factor)

  bf.tib <- bf.res %>%
    group_by(Dataset, Iteration, Method) %>%
    filter(Timestamp == max(Timestamp)) %>%
    group_by(Dataset) %>%
    mutate(Obj = Obj - max(Obj)) %>%
    ungroup() %>%
    filter(Timestamp > 0) %>%
    select(-Iteration)

  tib <- greedy.tib %>% add_column(Algorithm = "Greedy") %>%
    bind_rows(
      bf.tib %>% add_column(Algorithm = "Backfit")
    ) %>%
    mutate(Algorithm = factor(Algorithm, levels = c("Greedy", "Backfit")),
           Dataset = factor(Dataset, levels = c("GTEx", "PBMC", "Montoro")),
           Method = factor(Method, levels = c(
             "flashr", "efficient_ops", "fast_init", "extrapolation", "sparse"
           )))

  plt <- ggplot(tib, aes(x = Timestamp, y = Obj, col = Method)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Elapsed Time (s)", y = "ELBO (difference from optimal)") +
    scale_x_log10() +
    facet_wrap(~ Algorithm + Dataset, nrow = 3, ncol = 2, scales = "free", dir = "v")

  plot(plt)

  return(plt)
}

plot.init <- function(res) {
  tib <- res %>%
    mutate(Factor = factor(Factor)) %>%
    select(-TotalTime) %>%
    pivot_longer(c("InitTime", "RefineTime"), names_to = "Process", values_to = "Time") %>%
    mutate(Process = str_remove(Process, "Time"),
           Method = factor(Method, levels = c("default", "si", "irlba", "sparse")),
           AnyNA = factor(AnyNA, levels = c("FALSE", "TRUE")),
           Process = factor(Process, levels = c("Refine", "Init")),
           Method = recode(Method, "si" = "softImpute"),
           AnyNA = recode(AnyNA, "FALSE" = "Complete", "TRUE" = "Incomplete"))

  plt <- ggplot(tib, aes(x = Factor, y = Time)) +
    geom_bar(aes(fill = Process), stat = "identity") +
    theme_minimal() +
    facet_grid(rows = vars(Dataset, AnyNA), cols = vars(Method), scales = "free_y") +
    labs(y = "Elapsed Time (s)")

  plot(plt)

  return(plt)
}

greedy.res <- readRDS("../../output/res_greedy.rds")
bf.res <- readRDS("../../output/res_bf.rds")
plot.fits(greedy.res, bf.res)
ggsave("../../figs/benchmark_fits.png", height = 8, width = 6.5)

init.res <- readRDS("../../output/res_init.rds")
plot.init(init.res)
ggsave("../../figs/benchmark_init.png", height = 8, width = 6.5)
