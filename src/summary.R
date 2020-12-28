
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

## bam stats
f1 <- list.files(path = paste("bam/", sampleId, sep = ""), pattern = "filter.bam$", full.names = TRUE)
dBam <- foreach(ff = f1) %dopar% {
    p1 <- ScanBamParam(what = c("rname", "mapq", "qwidth", "cigar"), tag=c("NM"), flag = scanBamFlag(isDuplicate = FALSE))
    r1 <- scanBam(ff, param=p1)
    if(length(r1[[1]]$rname)== 0){
        NULL
    } else {
        r2 <- as_tibble(r1[[1]][c("rname", "mapq", "qwidth", "cigar")]) %>% mutate(nm = r1[[1]]$tag$NM)

        ## parse cigar string to count number of soft-clipped bases; map function converts all non-S containing cigar fields into NA
        r21 <- gsub("([A-Z])", "\\1.", r2$cigar) %>% strsplit("\\.") %>% map(function(x) as.integer(gsub("S", "", x)) %>% sum(na.rm = TRUE)) %>% as.integer()
        r2 <- mutate(r2, nSoftClip = r21)

        ## generate final summary tibble
        r3 <- gsub(".*\\/", "", f1)

        r41 <- r2 %>% group_by(rname) %>% count(nm) %>% mutate(p = n / sum(n)) %>% select(rname, nm, n, p) %>% summarise(nReads = sum(n), editDistMode = nm[which.max(n)], editDistAvg = sum(n * nm) / sum(n), editDistAvgDecay = mean(diff(n)) / sum(n), editDistDecayEnd = ifelse(is.na(which(diff(n) > 0)[1]), max(nm), which(diff(n) > 0)[1])) %>% ungroup()
        colnames(r41)[1] <- "contigId"
        r42 <- r2 %>% group_by(rname) %>% summarise(readLAvg = mean(qwidth), mqAvg = mean(mapq), nSoftClipAvg = mean(nSoftClip)) %>% ungroup() %>% select(-rname)

        r5 <- gsub(".*\\/", "", ff) %>% strsplit("\\.")
        r6 <- bind_cols(r41, r42) %>% mutate(genusId = map_chr(r5, 1), assemblyId = paste(map_chr(r5, 2), map_chr(r5, 3), sep = ".")) %>% select(genusId, assemblyId, contigId, nReads, editDistMode:nSoftClipAvg)
        r6
    }
}
dBam <- bind_rows(dBam)

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
seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccccc")
s1 <- filter(seqInfo, assemblyId %in% dGc$assemblyId) %>% select(assemblyId, taxId, taxIdSpecies, taxNameSpecies)

## final result table
res <- left_join(dGc, dBam, by = c("genusId", "assemblyId", "contigId")) %>% left_join(dDamage, by = c("genusId", "assemblyId", "contigId")) %>% mutate(sampleId = sampleId, flag = "") %>% left_join(s1, by = "assemblyId") %>% select(sampleId, genusId, taxIdSpecies, taxNameSpecies, assemblyId:dam3pAvgDecay, flag) %>% filter(!is.na(nReads))
idx <- res$dam5p >= 0.1 & res$dam5pAvgDecay < 0
idx1 <- res$editDistAvg <= 1.5 & res$editDistAvgDecay < 0
idx2 <- res$editDistDecayEnd >= 5
idx3 <- res$mqAvg >= 20
res$flag[idx & !idx1 & !idx2] <- "damage_0.1"
res$flag[!idx & idx1 & !idx2] <- "editDistAvg_1.5"
res$flag[!idx & !idx1 & idx2] <- "editDistDecay_5"
res$flag[idx & idx1 & !idx2] <- "damage_0.1;editDistAvg_1.5"
res$flag[idx & !idx1 & idx2] <- "damage_0.1;editDistDecay_5"
res$flag[!idx & idx1 & idx2] <- "editDistAvg_1.5;editDistDecay_5"
res$flag[idx & idx1 & idx2] <- "damage_0.1;editDistAvg_1.5;editDistDecay_5"
res$flag[idx & idx1 & idx2 & idx3] <- "damage_0.1;editDistAvg_1.5;editDistDecay_5;mqAvg_20"

write_tsv(res, path = paste("tables/", sampleId, "/", prefix, ".summary.tsv.gz", sep = ""), na = "NaN")
