do.volcano.plot <- function(tib, fname) {
  ggplot(tib, aes(x = pm, y = z, color = exprmean, label = SYMBOL)) +
    geom_point() +
    scale_color_gradient2(low = "deepskyblue", mid = "gold", high = "orangered",
                          na.value = "gainsboro",
                          midpoint = mean(range(tib$exprmean))) +
    scale_y_sqrt() +
    geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                    segment.color = "darkgray", segment.size = 0.25,
                    min.segment.length = 0, na.rm = TRUE) +
    theme_minimal() +
    labs(
      x = "Factor Loading (posterior mean)",
      y = "|z-score|",
      color = "Mean Expression (log10)"
    ) +
    theme(legend.position = "bottom")

  ggsave(paste0("../../figs/", fname, ".png"), height = 6, width = 6)
}


pbmc <- readRDS("../../data/pbmc.rds")
pbmc <- pbmc[, colSums(pbmc) < 15000]

pbmc.celltype <- sapply(strsplit(colnames(pbmc), "_"), `[[`, 2)
dendritic <- which(pbmc.celltype == "Dendritic")

pbmc.res <- readRDS("../../output/pbmc_res.rds")
fl <- pbmc.res$flashier.fit
dendritic.k <- which.max(abs(colMeans(fl$loadings.pm[[2]][dendritic, ])))
dendritic.pm <- fl$loadings.pm[[1]][, dendritic.k]
dendritic.z <- abs(fl$loadings.pm[[1]][, dendritic.k] / fl$loadings.psd[[1]][, dendritic.k])

pbmc.genemap <- readRDS("../../output/pbmc_genemap.rds")

pbmc.tib <- tibble(
  ENSEMBL = names(dendritic.pm),
  pm = dendritic.pm,
  z = dendritic.z
) %>%
  left_join(pbmc.genemap) %>%
  mutate(exprmean = log10(exprmean))

pbmc.tib <- pbmc.tib %>%
  mutate(SYMBOL = ifelse(z < 20 & pm < 0.1, "", SYMBOL))

do.volcano.plot(pbmc.tib, "dendritic")


montoro <- readRDS("../../data/sp_trachea.rds")
montoro <- expm1(montoro)
montoro <- montoro[, colSums(montoro) < 30000]

montoro.celltype <- sapply(strsplit(colnames(montoro), "_"), `[[`, 3)
ionocyte <- which(montoro.celltype == "Ionocyte")

montoro.res <- readRDS("../../output/montoro_res.rds")
fl <- montoro.res$flashier.fit
ionocyte.k <- which.max(abs(colMeans(fl$loadings.pm[[2]][ionocyte, ])))
ionocyte.pm <- fl$loadings.pm[[1]][, ionocyte.k]
ionocyte.z <- abs(fl$loadings.pm[[1]][, ionocyte.k] / fl$loadings.psd[[1]][, ionocyte.k])

montoro.genemap <- readRDS("../../output/montoro_genemap.rds")

montoro.tib <- tibble(
  SYMBOL = names(ionocyte.pm),
  pm = ionocyte.pm,
  z = ionocyte.z
) %>%
  left_join(montoro.genemap) %>%
  mutate(exprmean = log10(exprmean))

montoro.tib <- montoro.tib %>%
  mutate(SYMBOL = ifelse(z < 35 & pm < 0.1, "", SYMBOL))

do.volcano.plot(montoro.tib, "ionocyte")
