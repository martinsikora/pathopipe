
##########################
## plot damage patterns ##
##########################


## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(viridis))
suppressMessages(library(tidyr))
suppressMessages(library(purrr))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
dbPath <- args[2]
outPdf <- args[3]
metaDMG <- args[4]
infiles <- args[-1:-4]


## --------------------------------------------------------------------------------
## taxInfo

seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "ccccccccdd")

## helper, extract first field of prefix for matching in file name assemblyId parsing
prefix1 <- prefix %>%
    strsplit("\\.") %>%
    map_chr(1)


## --------------------------------------------------------------------------------
## get damage data

hdr <- c("contigId", "nReads", "end", "pos", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

d <- map_dfr(infiles, ~{
    bam <- gsub(".bdamage.gz", ".bam", .x) %>%
        gsub("metadamage", "bam", .)

    s1 <- gsub("metadamage/", "", .x) %>%
        strsplit("\\/") %>%
        map_chr(2) %>%
        strsplit("\\.") %>%
        unlist()

    idx <- match(prefix1, s1) ## index of field where prefix starts
    assemblyId <- paste(s1[2:(idx - 1)], collapse = ".")

    p <- pipe(paste(metaDMG, " print -howmany 25 ", .x, sep = ""))
    r <- read_tsv(p, col_names = hdr, col_types = "cdcidddddddddddddddd", skip = 1)

    r1 <- r %>%
        pivot_longer(all_of(hdr[-1:-4])) %>%
        mutate(refAllele = substr(name, 1, 1),
               p = value) %>%
        select(end, pos, refAllele, name, p) %>%
        mutate(assemblyId = assemblyId,
               taxNameSpecies = seqInfo$taxNameSpecies[match(assemblyId, seqInfo$assemblyId)])

    r1
})
d <- bind_rows(d)
d$pos[d$end == "3'"] <- -d$pos[d$end == "3'"]


## --------------------------------------------------------------------------------
## plot

d$end <- factor(d$end, levels = c("5'", "3'"))
d <- mutate(d, lab = paste(assemblyId, taxNameSpecies, sep = "\n")) %>%
    filter(!name %in% c("AA", "CC", "GG", "TT") )

st <- c("CT", "GA", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
idxSt <- c("red", "blue", rep("grey", length(st) - 2))
names(idxSt) <- st

h <- 2 + length(unique(d$assemblyId))

pdf(outPdf, width = 10, height = h)
if(nrow(d) == 0){
    dev.off()
} else {
    p <- ggplot(d, aes(x = pos + 1, y = p, colour = name))
    print(p +
      geom_line(size = 0.25) +
      facet_grid(lab ~ end, scales = "free") +
      xlab("Position from read end") +
      ylab("Frequency") +
      coord_cartesian(ylim = c(0, 0.3)) +
      scale_color_manual(values = idxSt) +
      theme_bw() +
      theme(strip.background = element_blank(),
            panel.grid.major = element_line(linetype = "dotted", size = 0.25),
            panel.grid.minor = element_blank(),
            strip.text.y = element_text(angle = 0, hjust = 0),
            legend.position = "none"))
    dev.off()
}
