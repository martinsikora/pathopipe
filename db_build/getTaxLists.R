
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

## ## extract per-genus tax IDs, write them to file and build summary table
## d2 <- filter(taxInfo, taxRank == "genus")
## registerDoParallel(48)
## r <- foreach(i = 1:nrow(d2)) %dopar% {
##     cat(i, "\r")

##     k <- match(d2$taxId[i], vertices) ## index of taxId in edge list
##     r1 <- subcomponent(gTaxInfo, k, "out") %>% names() %>% as.integer() %>% sort() ## all nodes in family

##     r1 <- r1[r1 != d2$taxId[i]]
##     f <- paste("taxLists/", d2$taxId[i], ".genus.taxIds.txt", sep = "")
##     write(r1, file = f, sep = "\n")

##     r2 <- tibble(taxRank = "genus", taxId = d2$taxId[i], taxName = d2$taxName[i], childIds = r1)
##     r2
## }

## r1 <- bind_rows(r)
## write_tsv(r1, path = "taxLists/genus.taxIds.tsv.gz")


## extract per kingdom tax IDs, and build summary table
d2 <- filter(taxInfo, taxRank == "kingdom")
registerDoParallel(48)
r <- foreach(i = 1:nrow(d2)) %dopar% {
    cat(i, "\r")

    k <- match(d2$taxId[i], vertices) ## index of taxId in edge list
    r1 <- subcomponent(gTaxInfo, k, "out") %>% names() %>% as.integer() %>% sort() ## all nodes in family

    r1 <- r1[r1 != d2$taxId[i]]
    f <- paste("taxLists/", d2$taxId[i], ".kingdom.taxIds.txt", sep = "")
    write(r1, file = f, sep = "\n")

    r2 <- tibble(taxRank = "kingdom", taxId = d2$taxId[i], taxName = d2$taxName[i], childIds = r1)
    r2
}

r1 <- bind_rows(r)
write_tsv(r1, path = "taxLists/kingdom.taxIds.tsv.gz")
