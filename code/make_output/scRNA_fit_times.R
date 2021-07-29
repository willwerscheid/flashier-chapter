pbmc.res <- readRDS("../../output/pbmc_res.rds")
montoro.res <- readRDS("../../output/montoro_res.rds")

pbmc.t <- pbmc.res$t
pbmc.t["flashier.snn"] <- pbmc.t["flashier.pn"] + pbmc.t["flashier.snn"]
montoro.t <- montoro.res$t

tib <- tibble(
  Method = c(names(pbmc.t), names(montoro.t)),
  Time = c(pbmc.t, montoro.t),
  Dataset = c(rep("PBMC-3k", length(pbmc.t)), rep("Montoro-droplet", length(montoro.t)))
)

sink("../../output/final_fit_times.txt")
tib %>%
  mutate(Time = Time / 3600) %>%
  pivot_wider(names_from = Method, values_from = Time) %>%
  column_to_rownames("Dataset") %>%
  xtable(digits = 2, align = "||rrrrrr||") %>%
  print(include.rownames = TRUE)
sink(NULL)
