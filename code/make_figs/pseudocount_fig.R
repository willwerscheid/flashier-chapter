dat_to_tib <- function(dat) {
  pc <- as.numeric(names(dat))
  elbo <- sapply(dat, `[[`, "elbo")
  adj.elbo <- sapply(dat, `[[`, "adj.elbo")
  elbo.idd <- sapply(dat, `[[`, "elbo.idd")
  tib <- tibble(Pseudocount = rep(pc, 2),
                Val = c(adj.elbo, elbo.idd),
                Type = rep(c("Adjusted ELBO", "IDD ELBO"), each = length(pc)))
  tib <- tib %>%
    group_by(Type) %>%
    mutate(isMax = (Val == max(Val))) %>%
    ungroup()

  return(tib)
}

pbmc <- readRDS("../../output/pseudocount_pbmc_greedy.rds")
pbmc_tib <- dat_to_tib(pbmc)

montoro <- readRDS("../../output/pseudocount_trachea_greedy.rds")
montoro_tib <- dat_to_tib(montoro)

tib <- pbmc_tib %>% add_column(Dataset = "PBMC-3k") %>%
  bind_rows(montoro_tib %>% add_column(Dataset = "Montoro-droplet")) %>%
  mutate(Dataset = factor(Dataset, levels = c("PBMC-3k", "Montoro-droplet")))

plt <- ggplot(tib, aes(x = Pseudocount, y = Val, col = isMax)) +
  geom_point() +
  scale_x_log10() +
  facet_wrap(~Dataset + Type, scales = "free") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(y = "")

ggsave("../../figs/pseudocounts.png", height = 4, width = 6)
