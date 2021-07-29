library(tidyverse)
library(Matrix)

preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  dat <- t(t(dat) / size.factors)
  gene_cts <- rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]
  return(list(dat = dat, size.factors = size.factors, excluded.genes = gene_cts < min.nzcts))
}

pbmc <- readRDS("../../data/pbmc.rds")
pbmc <- pbmc[, colSums(pbmc) < 15000]
pbmc.pp <- preprocess(pbmc)
pbmc.exprmean <- rowMeans(pbmc)[!pbmc.pp$excluded.genes]

montoro <- readRDS("../../data/sp_trachea.rds")
montoro <- expm1(montoro)
montoro <- montoro[, colSums(montoro) < 30000]
montoro.pp <- preprocess(montoro)
montoro.exprmean <- rowMeans(montoro)[!montoro.pp$excluded.genes]

library(org.Hs.eg.db)

hs <- org.Hs.eg.db
ensembl <- rownames(pbmc.pp$dat)

pbmc.gene.map <- select(
  hs,
  keys = ensembl,
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  group_by(ENSEMBL) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  left_join(tibble(
    ENSEMBL = names(pbmc.exprmean),
    exprmean = pbmc.exprmean
  ))

saveRDS(pbmc.gene.map, "../../output/pbmc_genemap.rds")

montoro.gene.map <- tibble(
  SYMBOL = names(montoro.exprmean),
  exprmean = montoro.exprmean
)

saveRDS(montoro.gene.map, "../../output/montoro_genemap.rds")
