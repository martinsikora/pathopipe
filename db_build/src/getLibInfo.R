
##############################################
## generate library info file for pathopipe ##
##############################################


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
## process library file info

dbReport <- read_tsv("database.report.tsv", comment = "#", col_type = "ddddddccc")

r <- dbReport %>%
    filter(rank == "assembly") %>%
    select(taxName, taxID) %>%
    rename(assemblyId = taxName, taxId = taxID) %>%
    mutate(fa = paste("library/", assemblyId, ".masked_sdust.fna", sep = ""),
           mmi = paste("minimap2/", assemblyId, ".mmi", sep = ""))


## get ref sequence lengths
registerDoParallel(threads)
seqComp <- foreach(f = r$fa) %dopar% {

    p <- pipe(paste("seqtk comp", f))
    r1 <- read_tsv(p, col_names = c("contigId", "l", "A", "C", "G", "T", "a2", "a3", "N", "CpG", "tv", "ts", "CpG-ts"),
                   col_types = "cdddddddddddd")
    r2 <- r1 %>%
        summarise(seqLTot = sum(l),
                  seqLMasked = seqLTot - sum(N)) %>%
        mutate(fa = f) %>%
        select(fa, everything())
    r2
}
seqComp <- bind_rows(seqComp)

## parse tax info
taxInfoGenus <- read_tsv("taxLists/genus.taxIds.tsv.gz", col_types = "cccc") %>%
    filter(childIds %in% r$taxId) %>%
    select(childIds, taxId, taxName)
colnames(taxInfoGenus) <- c("taxId", "taxIdGenus", "taxNameGenus")

taxInfoSpecies <- read_tsv("taxLists/species.taxIds.tsv.gz", col_types = "cccc") %>%
    filter(childIds %in% r$taxId) %>%
    select(childIds, taxId, taxName)
colnames(taxInfoSpecies) <- c("taxId", "taxIdSpecies", "taxNameSpecies")

seqInfo <- r %>%
    left_join(taxInfoGenus) %>%
    left_join(taxInfoSpecies) %>%
    left_join(seqComp) %>%
    mutate(taxNameSpecies = gsub(" ", "_", taxNameSpecies)) %>%
    select(assemblyId, taxId, taxIdSpecies, taxNameSpecies, taxIdGenus, taxNameGenus, fa, mmi, seqLTot, seqLMasked)

write_tsv(seqInfo, path = "library.seqInfo.tsv")
