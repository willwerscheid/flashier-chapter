args <- commandArgs(trailingOnly = TRUE)

valid.args <- args %in% c(
  "no-output",
  "no-figs",
  "out-to-file",
  "test"
)
if (!all(valid.args)) {
  stop("Command line argument ",
       min(which(!valid.args)),
       " not recognized.")
}

test <- "test" %in% args


library(tidyverse)
library(ggrepel)
library(grid)
library(gridExtra)
library(ebnm)
library(flashier)
library(Matrix)
library(xtable)


# Output:

if (!("no-output" %in% args)) {
  cat("\nGenerating output...\n\n")

  setwd("./make_output/")

  source("./benchmark_greedy.R")
  source("./benchmark_bf.R")
  source("./benchmark_init.R")
  source("./scRNA_fit_times.R")

  setwd("../")
}

if ("out-to-file" %in% args) {
  out_file <- file("../output/code_output.txt", open = "wt")
  sink(out_file)
  sink(out_file, type = "message")
}

if (!("no-output" %in% args)) {
  cat("\nGenerating output...\n\n")

  setwd("./make_output/")

  source("./misspec_noise.R")
  source("./pseudocount.R")

  setwd("../")
}


# Figures:

if (!("no-figs" %in% args)) {
  cat("\nGenerating figures...\n\n")

  setwd("./make_figs/")

  source("./benchmark_figs.R")
  source("./misspec_noise_fig.R")
  source("./pseudocount_fig.R")
  source("./scRNA_figs.R")
  source("./volcano_plots.R")

  setwd("../")
}


cat("\n\n")
sessionInfo()

