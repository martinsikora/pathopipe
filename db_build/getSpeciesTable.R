
##################################################
## generate lists of taxId for different genera ##
##################################################

## --------------------------------------------------------------------------------
## libraries

library(dplyr)
library(readr)
library(igraph)
library(foreach)
library(doParallel)


## --------------------------------------------------------------------------------
## process taxDB

taxInfo <- read_tsv("taxDB", col_names = c("taxId", "taxIdParent", "taxName", "taxRank"), col_types = "iicc")

## make graph from taxInfo
d1 <- select(taxInfo, taxIdParent, taxId)
gTaxInfo <- graph_from_data_frame(d1, directed = TRUE)
vertices <- V(gTaxInfo) %>% names()

## extract per-genus tax IDs, write them to file and build summary table
d2 <- filter(taxInfo, taxRank == "species")
registerDoParallel(48)
r <- foreach(i = 1:nrow(d2)) %dopar% {
    cat(i, "\r")

    k <- match(d2$taxId[i], vertices) ## index of taxId in edge list
    r1 <- subcomponent(gTaxInfo, k, "out") %>% names() %>% as.integer() %>% sort() ## all nodes in family

    r2 <- tibble(taxRank = "species", taxName = d2$taxName[i], childIds = r1)
    r2
}

r1 <- bind_rows(r)
write_tsv(r1, path = "species.taxIds.tsv.gz")


## extract per-superkingdom tax IDs, write them to file and build summary table
d2 <- filter(taxInfo, taxRank == "superkingdom")
registerDoParallel(48)
r <- foreach(i = 1:nrow(d2)) %dopar% {
    cat(i, "\r")

    k <- match(d2$taxId[i], vertices) ## index of taxId in edge list
    r1 <- subcomponent(gTaxInfo, k, "out") %>% names() %>% as.integer() %>% sort() ## all nodes in family

    r2 <- filter(taxInfo, taxId %in% r1, taxRank == "genus")
    r3 <- tibble(taxRank = "superkingdom", taxName = d2$taxName[i], childIds = r2$taxId, childNames = r2$taxName)
    r3
}


r1 <- bind_rows(r)
write_tsv(r1, path = "superkingdom.taxIds.tsv.gz")


## re-process seqInfo

s <- read_tsv("library.seqInfo.tsv", col_types = "cccccccc")

idxX <- s$assemblyId[duplicated(s$assemblyId)] %>% unique()
s1 <- filter(s, assemblyId %in% idxX)

s2 <- read_tsv("taxLists/genus.taxIds.tsv.gz", col_types = "cccc")
s3 <- left_join(s1, s2, by = c("taxId" = "childIds")) %>% filter(taxIdGenus == taxId.y) %>% mutate(taxIdGenus = taxId.y) %>% select(assemblyId:mmi)

ss <- filter(s, !assemblyId %in% idxX) %>% bind_rows(s3)
write_tsv(ss, path = "library.seqInfo.tsv")
