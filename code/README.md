# Code

1. Install the necessary packages:

install.packages(c("tidyverse", "Matrix", "glmpca", "ggrepel", "gridExtra", "xtable", "devtools"))
devtools::install_github("willwerscheid/flashier")
devtools::install_github("stephenslab/ebnm")
devtools::install_github("stephenslab/fastTopics")

2. Run the files in the run_separately directory (from that directory).

3. Run code.R from this directory. Command-line options include:
  no-output: omit output
  no-figs: omit figures
  out-to-file: print output (including progress updates and session Info) to file
  test: changes parameters so that the code runs quicker
