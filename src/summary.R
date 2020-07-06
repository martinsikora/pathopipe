
################################
## summary report for species ##
################################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(viridis))
suppressMessages(library(purrr))
suppressMessages(library(Rsamtools))
suppressMessages(library(doParallel))
suppressMessages(library(tidyr))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
dbPath <- args[2]
sampleId <- args[3]
threads <- as.integer(args[4])


## --------------------------------------------------------------------------------
## process tables and summarise

registerDoParallel(threads)

## coverage
f1 <- list.files(path = paste("bam/", sampleId, sep = ""), pattern = "genomecov", full.names = TRUE)
dGc <- foreach(ff = f1) %dopar% {
    r1 <- read_tsv(ff, col_type = "ciddd", col_names = c("contigId", "dp", "count", "l", "p")) %>% filter(contigId != "genome")
    if(nrow(r1) == 0){
        NULL
    } else {
        r2 <- r1 %>% group_by(contigId) %>% summarise(contigL = l[1], coverageAvg = sum(dp * count) / contigL, coverageSd = sqrt(sum((dp - coverageAvg)^2*count)/sum(count)), coverageBp = sum(count[dp>0]), coverageP = 1 - p[1], coverageCv = coverageSd / coverageAvg, coverageEvennessScore = 1 - sum((ceiling(coverageAvg) - dp[dp <= ceiling(coverageAvg)])*count[dp <= ceiling(coverageAvg)] / (ceiling(coverageAvg) * contigL)))
        r3 <- gsub(".*\\/", "", ff) %>% strsplit("\\.")
        r2 <- mutate(r2, genusId = map_chr(r3, 1), assemblyId = paste(map_chr(r3, 2), map_chr(r3, 3), sep = "."))
        select(r2, genusId, assemblyId, contigId, contigL, coverageAvg:coverageEvennessScore)
    }
}
dGc <- bind_rows(dGc)

## edit distance
f1 <- list.files(path = paste("bam/", sampleId, sep = ""), pattern = "filter.bam$", full.names = TRUE)
dEditDist <- foreach(ff = f1) %dopar% {
    p1 <- ScanBamParam(what = "rname", tag=c("NM"), flag = scanBamFlag(isDuplicate = FALSE))
    r1 <- scanBam(ff, param=p1)
    r2 <- tibble(contigId = as.character(r1[[1]]$rname), editDist = r1[[1]]$tag$NM)
    if(nrow(r2) == 0){
        NULL
    } else {
        r3 <- gsub(".*\\/", "", f1)
        r4 <- r2 %>% group_by(contigId) %>% count(editDist) %>% mutate(p = n / sum(n)) %>% select(contigId, editDist, n, p) %>% summarise(nReads = sum(n), editDistMode = editDist[which.max(n)], editDistAvg = sum(n * editDist) / sum(n), editDistAvgDecay = mean(diff(n)) / sum(n)) %>% ungroup()
        r5 <- gsub(".*\\/", "", ff) %>% strsplit("\\.")
        mutate(r4, genusId = map_chr(r5, 1), assemblyId = paste(map_chr(r5, 2), map_chr(r5, 3), sep = ".")) %>% select(genusId, assemblyId, contigId, nReads, editDistMode:editDistAvgDecay)
    }
}
dEditDist <- bind_rows(dEditDist)

## damage
st <- c("C>T", "G>A")
ref <- c("C", "G")

f1 <- list.files(path = paste("mapdamage/", sampleId, sep = ""), pattern = "misincorporation.txt.gz$", recursive = TRUE, full.names = TRUE)
dDamage <- foreach(ff = f1) %dopar% {

    r <- read_tsv(ff, comment = "#", col_types = "cccddddddddddddddddddddddddddd")
    if(nrow(r) == 0){
        NULL
    }
    else {
        r1 <- select(r, Chr:Pos, Total, all_of(c(ref, st))) %>% pivot_longer(all_of(c(st, ref))) %>% mutate(refAllele = substr(name, 1, 1)) %>% group_by(Chr, Pos, End, name, refAllele) %>% summarise(value = sum(value)) %>% group_by(Chr, Pos, End, refAllele) %>% mutate(p = value / value[name == refAllele]) %>% ungroup() %>% filter((name == "C>T" & End == "5p") | (name == "G>A" & End == "3p"))
        r2 <- group_by(r1, Chr) %>% summarise(dam5p = p[Pos == 1 & name == "C>T"], dam3p = p[Pos == 1 & name == "G>A"], dam5pAvgDecay = mean(diff(p[Pos <= 5 & name == "C>T"])), dam3pAvgDecay = mean(diff(p[Pos <= 5 & name == "G>A"]))) %>% ungroup()
        r3 <- strsplit(ff, "\\/") %>% map_chr(3) %>% strsplit("\\.")
        r2 <- mutate(r2, genusId = map_chr(r3, 1), assemblyId = paste(map_chr(r3, 2), map_chr(r3, 3), sep = "."))
        colnames(r2)[1] <- "contigId"
        select(r2, genusId, assemblyId, contigId, dam5p:dam3pAvgDecay)
    }
}
dDamage <- bind_rows(dDamage)

## taxInfo
seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccc")
s1 <- filter(seqInfo, assemblyId %in% dGc$assemblyId) %>% select(assemblyId, taxId, taxIdSpecies, taxNameSpecies)

## final result table
res <- left_join(dGc, dEditDist, by = c("genusId", "assemblyId", "contigId")) %>% left_join(dDamage, by = c("genusId", "assemblyId", "contigId")) %>% mutate(sampleId = sampleId, flag = "") %>% left_join(s1, by = "assemblyId") %>% select(sampleId, genusId, taxIdSpecies, taxNameSpecies, assemblyId:dam3pAvgDecay) %>% filter(!is.na(nReads))
idx <- res$dam5p >= 0.1 & res$dam3p >= 0.1 & res$dam5pAvgDecay < 0
idx1 <- res$editDistAvg <= 1.5 & res$editDistAvgDecay < 0
res$flag[idx & !idx1] <- "damage_0.1"
res$flag[!idx & idx1] <- "editDist_low"
res$flag[idx & idx1] <- "damage_0.1;editDist_low"
write_tsv(res, path = paste("tables/", sampleId, "/", prefix, ".summary.tsv", sep = ""), na = "NaN")
