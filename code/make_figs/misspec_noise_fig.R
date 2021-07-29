misspec_res <- readRDS("../../output/misspec_noise.rds")

misspec_res <- misspec_res %>%
  mutate(distn = recode(distn,
                        "normal" = "Normal",
                        "poisson" = "Poisson",
                        "poislognorm" = "Poisson-lognormal")) %>%
  mutate(distn = factor(distn, levels = c("Normal", "t5", "Poisson", "Poisson-lognormal")))

ggplot(misspec_res, aes(x = n, y = nfactors, col = n)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 10), linetype = "dashed") +
  facet_wrap(~distn) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(y = "Number of Flash Factors Added")

ggsave("../../figs/misspec_noise.png", height = 4.5, width = 6)
