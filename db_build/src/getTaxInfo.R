
###########################################
## generate tax info table for pathopipe ##
###########################################


## --------------------------------------------------------------------------------
## libraries

library(tidyverse)
library(igraph)
library(foreach)
library(doParallel)


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
threads <- args[1] %>%
    as.integer()


## --------------------------------------------------------------------------------
## process taxDB

cat("__ processing taxonomy DB __\n")

taxInfo <- read_tsv("taxDB", col_names = c("taxId", "taxIdParent", "taxName", "taxRank"), col_types = "cccc")
seqmap <- read_tsv("seqid2taxid.map", col_names = c("accession", "taxId"), col_types = "cc")
dbReport <- read_tsv("database.report.tsv", comment = "#", col_type = "ddddddccc")


## make graph from taxInfo
d1 <- select(taxInfo, taxIdParent, taxId)
gTaxInfo <- graph_from_data_frame(d1, directed = TRUE)

## extract per-genus tax IDs, build summary table and write to file
cat("__ extracting genus taxonomy IDs __\n")

taxIdsGenus <- dbReport %>%
    filter(rank == "genus") %>%
    pull(taxID)

registerDoParallel(threads)
r <- foreach(i = taxIdsGenus) %dopar% {

    r1 <- subcomponent(gTaxInfo, i, "out") %>%
        names() ## all nodes in genus

    r2 <- tibble(taxRank = "genus",
                 taxId = i,
                 taxName = taxInfo$taxName[match(i, taxInfo$taxId)],
                 childIds = r1)
    r2
}
r1 <- bind_rows(r)
write_tsv(r1, file = "taxLists/genus.taxIds.tsv.gz")


## extract per-species tax IDs, build summary table and write to file
cat("__ extracting species taxonomy IDs __\n")

taxIdsSpecies <- dbReport %>%
    filter(rank == "species") %>%
    pull(taxID)

registerDoParallel(threads)
r <- foreach(i = taxIdsSpecies) %dopar% {

    r1 <- subcomponent(gTaxInfo, i, "out") %>%
        names() ## all nodes in genus

    r2 <- tibble(taxRank = "species",
                 taxId = i,
                 taxName = taxInfo$taxName[match(i, taxInfo$taxId)],
                 childIds = r1)
    r2
}
r1 <- bind_rows(r)
write_tsv(r1, file = "taxLists/species.taxIds.tsv.gz")
