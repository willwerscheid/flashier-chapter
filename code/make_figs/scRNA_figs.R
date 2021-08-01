plot.factors <- function(res,
                         dataset_name,
                         cell.types,
                         genemap,
                         max.pt.size = 2) {
  tib <- res$res %>%
    mutate(Cell.type = rep(cell.types, length.out = nrow(res$res)))
  tib <- tib %>%
    pivot_longer(
      -c(Method, Cell.type),
      names_to = "Factor",
      values_to = "Loading",
      values_drop_na = TRUE
    )

  # Re-normalize loadings so that factors are equally spread out; to make
  #   it easier to compare factors, flip them to make the largest loadings
  #   positive.
  tib <- tib %>%
    group_by(Method, Factor) %>%
    mutate(Flip = 2L * (max(abs(Loading)) == max(Loading)) - 1) %>%
    mutate(Loading = Flip * Loading / max(abs(Loading))) %>%
    select(-Flip) %>%
    ungroup()

  # Order factors by mean(abs(x)) as a proxy for sparsity.
  factor.tib <- tib %>%
    group_by(Method, Factor) %>%
    summarize(Loading_var = mean(abs(Loading))) %>%
    mutate(Factor_order = n() - rank(Loading_var) + 1) %>%
    select(-Loading_var) %>%
    ungroup()
  tib <- tib %>%
    left_join(factor.tib, by = c("Method", "Factor"))

  # Make the size of the point depend on how many of that type there are.
  celltype.tib <- tib %>%
    filter(Method == tib$Method[1], Factor == tib$Factor[1]) %>%
    group_by(Cell.type) %>%
    summarize(n = n())
  celltype.tib <- celltype.tib %>%
    mutate(Pt.size = max.pt.size / sqrt(n / min(n)))
  tib <- tib %>%
    left_join(celltype.tib %>% select(Cell.type, Pt.size), by = "Cell.type")

  tib <- tib %>%
    mutate(Method = factor(Method, levels = unique(res$res$Method)))

  plt <- ggplot(tib, aes(x = Factor_order, y = Loading, color = Cell.type)) +
    geom_jitter(position = position_jitter(width = 0.4, height = 0),
                size = tib$Pt.size) +
    scale_x_continuous(breaks = seq(from = 2, to = max(tib$Factor_order), by = 2)) +
    scale_color_brewer(palette = "Set3") +
    labs(x = "factor index", y = "factor loading", color = "Cell Type") +
    theme_minimal() +
    facet_wrap(~Method, ncol = 1, scales = "free_y")

  ggsave(paste0("../../figs/", dataset_name, "_factors.png"), plot = plt, width = 8, height = 10.5)

  z.scores <- abs(res$flashier.fit$loadings.pm[[1]] / res$flashier.fit$loadings.psd[[1]])
  top.genes <- apply(z.scores, 2, function(k) {
    rownames(res$flashier.fit$loadings.pm[[1]])[order(k, decreasing = TRUE)[1:10]]
  })
  colnames(top.genes) <- paste0("k=", formatC(1:ncol(top.genes), width = 2, flag = "0"))
  if (dataset_name == "pbmc") {
    genes.tib <- as_tibble(top.genes) %>%
      pivot_longer(cols = everything(), names_to = "Factor", values_to = "ENSEMBL") %>%
      left_join(genemap) %>%
      filter(!is.na(SYMBOL)) %>%
      select(Factor, SYMBOL)
  } else {
    genes.tib <- as_tibble(top.genes) %>%
      pivot_longer(cols = everything(), names_to = "Factor", values_to = "SYMBOL")
  }
  genes.tib <- genes.tib %>%
    group_by(Factor) %>%
    summarize(TopGenes = paste(SYMBOL, collapse = ", ")) %>%
    left_join(factor.tib %>% filter(Method == "flashier.snn")) %>%
    select(Factor_order, TopGenes) %>%
    arrange(Factor_order)

  tbl <- genes.tib %>%
    gt() %>%
    cols_label(Factor_order = "Factor", TopGenes = "Top Genes (by |z|-score)") %>%
    cols_align("left", columns = Factor_order) %>%
    cols_align("left", columns = TopGenes) %>%
    cols_width(
      Factor_order ~ pct(10),
      TopGenes ~ pct(90)
    ) %>%
    opt_row_striping()

  gtsave(tbl, paste0("../../figs/topgenes_", dataset_name, ".png"), vwidth = 800, vheight = 50 * nrow(genes.tib))

  return(NULL)
}

plot.t <- function(res,
                   dataset_name) {
  elapsed_t <- res$t
  elapsed_t["flashier.snn"] <- elapsed_t["flashier.snn"] + elapsed_t["flashier.pn"]

  tib <- tibble(Method = names(elapsed_t),
                Elapsed_min = elapsed_t / 60)
}


set.seed(666)

pbmc <- readRDS("../../data/pbmc.rds")
pbmc <- pbmc[, colSums(pbmc) < 15000]
pbmc.celltype <- sapply(strsplit(colnames(pbmc), "_"), `[[`, 2)
pbmc.libsize <- colSums(pbmc)
rm(pbmc)

pbmc.res <- readRDS("../../output/pbmc_res.rds")
pbmc.genemap <- readRDS("../../output/pbmc_genemap.rds")
plt <- plot.factors(pbmc.res, "pbmc", pbmc.celltype, pbmc.genemap)


montoro <- readRDS("../../data/sp_trachea.rds")
montoro <- expm1(montoro)
montoro <- montoro[, colSums(montoro) < 30000]
montoro.celltype <- sapply(strsplit(colnames(montoro), "_"), `[[`, 3)
montoro.libsize <- colSums(montoro)
rm(montoro)

montoro.res <- readRDS("../../output/montoro_res.rds")
montoro.genemap <- readRDS("../../output/montoro_genemap.rds")
plt <- plot.factors(montoro.res, "montoro", montoro.celltype, montoro.genemap)


# plot.umaps <- function(res,
#                        dataset_name,
#                        cell.types,
#                        libsizes,
#                        max.pt.size = 2) {
#   umap.dat <- res$res %>%
#     group_by(Method) %>%
#     mutate_at(vars(-Method), ~. / max(abs(.))) %>%
#     group_split()
#
#   umap.res <- lapply(umap.dat, function(x) {
#     umap(as.matrix(x %>% select(-Method)), verbose = TRUE)
#   })
#   umap.res <- lapply(umap.res, `[[`, "layout")
#   umap.res <- do.call(rbind, umap.res)
#   colnames(umap.res) = c("umap.x", "umap.y")
#
#   tib <- res$res %>%
#     bind_cols(as_tibble(umap.res)) %>%
#     mutate(Cell.type = rep(cell.types, length.out = nrow(res$res)),
#            Lib.size = rep(log(libsizes), length.out = nrow(res$res)))
#
#   celltype.tib <- tib %>%
#     filter(Method == tib$Method[1]) %>%
#     group_by(Cell.type) %>%
#     summarize(n = n())
#   celltype.tib <- celltype.tib %>%
#     mutate(Pt.size = max.pt.size / sqrt(n / min(n)))
#   tib <- tib %>%
#     left_join(celltype.tib %>% select(Cell.type, Pt.size), by = "Cell.type")
#
#   plt1 <- ggplot(tib, aes(x = umap.x, y = umap.y, col = Cell.type)) +
#     geom_point(size = tib$Pt.size) +
#     scale_color_brewer(palette = "Set3") +
#     theme_void() +
#     facet_wrap(~Method, ncol = 1, scales = "free") +
#     theme(legend.position = "none") +
#     labs(color = "Cell Type")
#
#   plt2 <- ggplot(tib, aes(x = umap.x, y = umap.y, col = Lib.size)) +
#     geom_point() +
#     scale_color_viridis_c() +
#     theme_void() +
#     facet_wrap(~Method, ncol = 1, scales = "free") +
#     theme(legend.position = "none")
#
#   plt <- arrangeGrob(plt1, plt2, ncol = 2, widths = c(4, 4))
#   ggsave(paste0("../../figs/", dataset_name, "_umaps.png"), plot = plt, width = 8, height = 10.5)
#
#   return(plt)
# }

# plt <- plot.umaps(pbmc.res, "pbmc", pbmc.celltype, pbmc.libsize)
# plt <- plot.umaps(montoro.res, "montoro", montoro.celltype, montoro.libsize)
