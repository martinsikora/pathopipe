
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
## helpers

registerDoParallel(threads)

## helper, extract first field of prefix for matching in file name assemblyId parsing
prefix1 <- strsplit(prefix, "\\.") %>%
    map_chr(1)

seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccccc")


## --------------------------------------------------------------------------------
## genome coverage summary

f1 <- list.files(path = paste("bam/", sampleId, sep = ""),
                 pattern = "genomecov",
                 full.names = TRUE)
dGc <- foreach(ff = f1) %dopar% {
    r1 <- read_tsv(ff,
                   col_type = "ciddd",
                   col_names = c("contigId", "dp", "count", "l", "p"))
    if(nrow(r1) == 0){
        NULL
    } else {

        ## find contigs with coverage
        r11 <- r1 %>%
            filter(dp > 0,
                   contigId != "genome") %>%
            distinct(contigId)

        nContigs <- r1 %>%
            pull(contigId) %>%
            unique() %>%
            length() -1

        ## coverage summary stats
        r2 <- r1 %>%
            filter(contigId %in% r11$contigId) %>%
            group_by(contigId) %>%
            summarise(contigL = l[1],
                      coverageAvg = sum(dp * count) / contigL,
                      coverageSd = sqrt(sum((dp - coverageAvg)^2*count) / sum(count)),
                      coverageBp = sum(count[dp>0]),
                      coverageP = 1 - p[1],
                      coveragePExp = 1 - exp(-coverageAvg),
                      coveragePRatio = coverageP / coveragePExp,
                      coverageCv = coverageSd / coverageAvg,
                      coverageEvennessScore = 1 - sum((ceiling(coverageAvg) - dp[dp <= ceiling(coverageAvg)])*count[dp <= ceiling(coverageAvg)] / (ceiling(coverageAvg) * contigL)))

        r21 <- r1 %>%
            filter(contigId == "genome") %>%
            summarise(contigId = "genome",
                      contigL = l[1],
                      coverageAvg = sum(dp * count) / contigL,
                      coverageSd = sqrt(sum((dp - coverageAvg)^2*count) / sum(count)),
                      coverageBp = sum(count[dp>0]),
                      coverageP = 1 - p[1],
                      coveragePExp = 1 - exp(-coverageAvg),
                      coveragePRatio = coverageP / coveragePExp,
                      coverageCv = coverageSd / coverageAvg,
                      coverageEvennessScore = 1 - sum((ceiling(coverageAvg) - dp[dp <= ceiling(coverageAvg)])*count[dp <= ceiling(coverageAvg)] / (ceiling(coverageAvg) * contigL))) %>%
            mutate(nContigs = nContigs,
                   pContigsCovered = nrow(r11) / nContigs)

        r2 <- bind_rows(r2, r21)

        ## split fields separated by '.' in input file names
        r3 <- ff %>%
            gsub(".*\\/", "", .) %>%
            strsplit("\\.") %>%
            unlist()

        idx <- match(prefix1, r3) ## index of field where prefix starts
        r2 <- r2 %>%
            mutate(genusId = r3[1],
                   assemblyId = paste(r3[2:(idx - 1)], collapse = "."),
                   coveragePRatioCorr = NaN)


        ## add coveragePRatio corrected for non-masked sites
        idx <- which(r2$contigId == "genome")
        s1 <- seqInfo %>%
            filter(assemblyId == r2$assemblyId[1])

        r2$coveragePRatioCorr[idx] <- r2$coverageP[idx] / (r2$coveragePExp[idx] * s1$seqLMasked / s1$seqLTot)
        select(r2, genusId, assemblyId, contigId, contigL, coverageAvg:coveragePRatio, coveragePRatioCorr, coverageCv:pContigsCovered)
    }
}
dGc <- bind_rows(dGc)


## --------------------------------------------------------------------------------
## bam stats

f1 <- list.files(path = paste("bam/", sampleId, sep = ""), pattern = "filter.bam$", full.names = TRUE)
dBam <- foreach(ff = f1) %dopar% {

    p1 <- ScanBamParam(what = c("rname", "mapq", "qwidth", "cigar"),
                       tag=c("NM"),
                       flag = scanBamFlag(isDuplicate = FALSE))
    r1 <- scanBam(ff, param=p1)

    if(length(r1[[1]]$rname)== 0){
        NULL
    } else {
        r2 <- as_tibble(r1[[1]][c("rname", "mapq", "qwidth", "cigar")]) %>%
            mutate(nm = r1[[1]]$tag$NM,
                   rname = as.character(rname))

        ## parse cigar string to count number of soft-clipped bases; map function converts all non-S containing cigar fields into NA
        r21 <- r2 %>%
            pull(cigar) %>%
            gsub("([A-Z])", "\\1.", .) %>%
            strsplit("\\.") %>%
            map(function(x) as.integer(gsub("S", "", x)) %>%
                            sum(na.rm = TRUE)) %>%
            as.integer()
        r2 <- r2 %>%
            mutate(nSoftClip = r21)

        ## generate final summary tibble
        r41 <- r2 %>%
            group_by(rname) %>%
            count(nm) %>%
            mutate(p = n / sum(n)) %>%
            select(rname, nm, n, p) %>%
            summarise(nReads = sum(n),
                      editDistMode = nm[which.max(n)],
                      editDistAvg = sum(n * nm) / sum(n),
                      editDistAvgDecay = mean(diff(n)) / sum(n),
                      editDistDecayEnd = ifelse(is.na(which(diff(n) > 0)[1]), max(nm), which(diff(n) > 0)[1])) %>%
            ungroup()
        colnames(r41)[1] <- "contigId"

        r411 <- r2 %>%
            count(nm) %>%
            mutate(p = n / sum(n)) %>%
            select(nm, n, p) %>%
            summarise(nReads = sum(n),
                      editDistMode = nm[which.max(n)],
                      editDistAvg = sum(n * nm) / sum(n),
                      editDistAvgDecay = mean(diff(n)) / sum(n),
                      editDistDecayEnd = ifelse(is.na(which(diff(n) > 0)[1]), max(nm), which(diff(n) > 0)[1])) %>%
            mutate(contigId = "genome")
        r41 <- bind_rows(r41, r411)

        r42 <- r2 %>%
            group_by(rname) %>%
            summarise(readLAvg = mean(qwidth),
                      mqAvg = mean(mapq),
                      nSoftClipAvg = mean(nSoftClip),
                      ani = 1 - sum(nm)/sum(qwidth)) %>%
            ungroup() %>%
            select(-rname)
        r421 <- r2 %>%
            summarise(readLAvg = mean(qwidth),
                      mqAvg = mean(mapq),
                      nSoftClipAvg = mean(nSoftClip),
                      ani = 1 - sum(nm)/sum(qwidth))
        r42 <- bind_rows(r42, r421)

        ## split fields separated by '.' in input file names
        r5 <- ff %>%
            gsub(".*\\/", "", .) %>%
            strsplit("\\.") %>%
            unlist()

        idx <- match(prefix1, r5) ## index of field where prefix starts
        r6 <- bind_cols(r41, r42) %>%
            mutate(genusId = r5[1],
                   assemblyId = paste(r5[2:(idx - 1)], collapse = ".")) %>%
            select(genusId, assemblyId, contigId, nReads, editDistMode:ani)
        r6
    }
}
dBam <- bind_rows(dBam)

dB1 <- dBam %>%
    filter(contigId == "genome") %>%
    group_by(genusId) %>%
    mutate(aniRank = rank(dplyr::desc(ani))) %>%
    ungroup() %>%
    select(assemblyId, contigId, aniRank)

dB2 <- dBam %>%
    filter(contigId == "genome",
           nReads >= 100) %>%
    group_by(genusId) %>%
    mutate(aniRank100 = rank(dplyr::desc(ani))) %>%
    ungroup() %>%
    select(assemblyId, contigId, aniRank100)

dBam <- dBam %>%
    left_join(dB1) %>%
    left_join(dB2)


## --------------------------------------------------------------------------------
## damage

hdr <- c("contigId", "nReads", "end", "pos", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

f1 <- list.files(path = paste("metadamage/", sampleId, sep = ""), pattern = "bdamage.gz$", full.names = TRUE)
ct <- c("CC", "CA", "CG", "CT")
ga <- c("GG", "GC", "GT", "GA")
dDamage <- foreach(ff = f1) %dopar% {

    bam <- gsub(".bdamage.gz", ".bam", ff) %>%
        gsub("metadamage", "bam", .)

    s1 <- gsub("metadamage/", "", ff) %>%
        strsplit("\\/") %>%
        map_chr(2) %>%
        strsplit("\\.") %>%
        unlist()

    p <- pipe(paste("metadamage print -countout -howmany 25 -bam ", bam, " ", ff))
    r <- read_tsv(p, col_names = hdr, col_types = "cdcidddddddddddddddd", skip = 1)

    if(nrow(r) == 0){
        NULL
    }
    else {
        r1 <- r %>%
            filter(nReads > 0) %>%
            select(contigId, end, pos, all_of(ct), all_of(ga)) %>%
            pivot_longer(c(all_of(ct), all_of(ga))) %>%
            filter((end == "5'" & name %in% ct) | (end == "3'" & name %in% ga))

        r2 <- group_by(r1, contigId, end, pos) %>%
            summarise(dam = value[name %in% c("CT", "GA")] / sum(value)) %>%
            mutate(dam = replace_na(dam, 0)) %>%
            group_by(contigId) %>%
            summarise(dam5p = dam[end == "5'" & pos == 0],
                      dam3p = dam[end == "3'" & pos == 0],
                      dam5pAvgDecay = mean(diff(dam[pos < 5 & end == "5'"])),
                      dam3pAvgDecay = mean(diff(dam[pos < 5 & end == "3'"]))) %>%
            ungroup()

        r21 <- group_by(r1, end, pos) %>%
            summarise(dam = sum(value[name %in% c("CT", "GA")]) / sum(value)) %>%
            mutate(dam = replace_na(dam, 0)) %>%
            ungroup() %>%
            summarise(dam5p = dam[end == "5'" & pos == 0],
                      dam3p = dam[end == "3'" & pos == 0],
                      dam5pAvgDecay = mean(diff(dam[pos < 5 & end == "5'"])),
                      dam3pAvgDecay = mean(diff(dam[pos < 5 & end == "3'"]))) %>%
            mutate(contigId = "genome")

        idx <- match(prefix1, s1) ## index of field where prefix starts
        r2 <- bind_rows(r2, r21) %>%
            mutate(genusId = s1[1],
                   assemblyId = paste(s1[2:(idx - 1)], collapse = "."))
        select(r2, genusId, assemblyId, contigId, dam5p:dam3pAvgDecay)
    }
}
dDamage <- bind_rows(dDamage)


## --------------------------------------------------------------------------------
## krakenUniq results

f1 <- paste("report/", sampleId, ".", prefix, ".krakenuniq.report.tsv.gz", sep = "")

r1 <- read_tsv(f1, c("freq", "nClade", "nTaxon", "nKmer", "dup", "cov", "taxId", "taxRank", "name"), col_types = "ddddddccc", skip = 4)
s2 <- seqInfo %>%
    distinct(taxIdSpecies, taxNameSpecies, taxIdGenus, taxNameGenus)

dKraken <- r1 %>%
    filter(taxRank == "species") %>%
    left_join(s2, by = c("taxId" = "taxIdSpecies")) %>%
    select(taxIdGenus, taxNameGenus, taxId, name, nClade, nKmer, dup) %>%
    group_by(taxIdGenus) %>%
    mutate(kmerRank = rank(-nKmer)) %>%
    ungroup() %>%
    select(taxId, nClade:kmerRank)

colnames(dKraken) <- c("taxIdSpecies", "krakenNClade", "krakenNKmer", "krakenDupRate", "krakenKmerRank")


## --------------------------------------------------------------------------------
## final result table

s1 <- filter(seqInfo, assemblyId %in% dGc$assemblyId) %>%
    select(assemblyId, taxId, taxIdSpecies, taxNameSpecies)

res <- left_join(dGc, dBam, by = c("genusId", "assemblyId", "contigId")) %>%
    left_join(dDamage, by = c("genusId", "assemblyId", "contigId")) %>%
    mutate(sampleId = sampleId, flag = "") %>%
    left_join(s1, by = "assemblyId") %>%
    left_join(dKraken) %>%
    select(sampleId, genusId, taxIdSpecies, taxNameSpecies, assemblyId:dam3pAvgDecay, krakenNClade:krakenKmerRank, flag) %>%
    filter(!is.na(nReads))

idx <- res$dam5p >= 0.1 & res$dam5pAvgDecay < 0
idx1 <- res$aniRank100 < 2
idx2 <- res$coveragePRatioCorr >= 0.8
idx3 <- res$krakenKmerRank < 2 & !is.na(res$krakenKmerRank)

res$flag[idx] <- "damage_0.1"
res$flag[idx1] <- "aniRank100_1"
res$flag[idx2] <- "coveragePRatioCorr_0.8"
res$flag[idx3] <- "krakenKmerRank_1"

res$flag[idx & idx1] <- "damage_0.1;aniRank100_1"
res$flag[idx & idx2] <- "damage_0.1;coveragePRatioCorr_0.8"
res$flag[idx & idx3] <- "damage_0.1;krakenKmerRank_1"
res$flag[idx1 & idx2] <- "aniRank100_1;coveragePRatioCorr_0.8"
res$flag[idx1 & idx3] <- "aniRank100_1;krakenKmerRank_1"
res$flag[idx2 & idx3] <- "coveragePRatio_0.8;krakenKmerRank_1"

res$flag[idx & idx1 & idx2] <- "damage_0.1;aniRank100_1;coveragePRatioCorr_0.8"
res$flag[idx & idx1 & idx3] <- "damage_0.1;aniRank100_1;krakenKmerRank_1"
res$flag[idx & idx2 & idx3] <- "damage_0.1;coveragePRatioCorr_0.8;krakenKmerRank_1"
res$flag[idx1 & idx2 & idx3] <- "aniRank100_1;coveragePRatioCorr_0.8;krakenKmerRank_1"

res$flag[idx & idx1 & idx2 & idx3] <- "damage_0.1;aniRank100_1;coveragePRatioCorr_0.8;krakenKmerRank_1"

write_tsv(res, path = paste("tables/", sampleId, "/", prefix, ".summary.tsv.gz", sep = ""), na = "NaN")
