################################################################################
#### Project: Arctic priming
#### Title:   Drake plan
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 July 2020
#### ---------------------------------------------------------------------------


#### PROLOGUE ------------------------------------------------------------------

## Options ----
# remove objects from global environment
rm(list = ls())

# configure default R session options (no factors, bias against scientific #s)
options(stringsAsFactors = F,
scipen = 6)


## Libraries ----
# standard library set
source("packages.r")


## Source scripts ----
# code paths
codePaths <- list.files(
  paste0(getwd(), "/code_clean/"),
  full.names = T
)
# source functions
sapply(
  codePaths, 
  source,
)


#### CLEAN DATA ----------------------------------------------------------------

## Load and clean microbial community data ----
mcc_plan <- drake_plan(
  mcc_dist = "bray",
  raw_mcc = load_mcc(),
  fungi = compile_fungi(raw_mcc = raw_mcc,
                        distance = mcc_dist),
  bacteria = compile_bacteria(raw_mcc = raw_mcc,
                              distance = mcc_dist)
)

## Load and clean other data ----
other_data_plan <- drake_plan(
  mtbs_dist = "cosine",
  other_data = load_other_data(),
  metabolites = clean_metabolites(other_data = other_data,
                                  distance = mtbs_dist),
  priming = clean_priming(other_data = other_data),
  abiotic = clean_physical(other_data = other_data)
)

# collate plan
all_plan <- bind_rows(mcc_plan, 
                      other_data_plan)
make(all_plan)

# visualize workflow
dependencies <- vis_drake_graph(
  all_plan, 
  targets_only = F,
  mode = "all",
  show_output_files = T,
  main = ""
)
dependencies %>%
  visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")
