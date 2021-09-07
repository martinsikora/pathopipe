
############################################
## summary of ancient hits across samples ##
############################################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(viridis))
suppressMessages(library(purrr))
suppressMessages(library(doParallel))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(ggrepel))
suppressMessages(library(scales))
suppressMessages(library(ggh4x))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
targetsPriority <- args[2]
files <- args[-1:-2]


## --------------------------------------------------------------------------------
## process tables and summarise

d <- foreach(f = files) %do% {
    r <- read_tsv(f, col_types = "ccccccddddddddddddddddddddddddddc")
    r
}
d <- bind_rows(d)
write_tsv(d, file = paste("tables/", prefix, ".summary.tsv.gz", sep = ""))

tP <- read_tsv(targetsPriority, col_names = FALSE)

d1 <- d %>%
    filter(genusId %in% tP$X1)
write_tsv(d1, file = paste("tables/", prefix, ".targets_priority.summary.tsv.gz", sep = ""))


## --------------------------------------------------------------------------------
## plot matrix of candidates with aniRank 1 and damage

## find unique species name with ANI rank 1
idxS <- d1 %>%
    filter(aniRank100 < 2,
           dam5p >= 0.1) %>%
    pull(taxNameSpecies) %>%
    unique()

## set up plot data
d3 <- d %>%
    filter(taxNameSpecies %in% idxS) %>%
    mutate(damClass = case_when(
               dam5p >= 0.1 ~ "dam5p_010",
               dam5p >= 0.05 ~ "dam5p_005",
               TRUE ~ "dam5p_low"
           ))

w <- d3$contigId %>%
    unique() %>%
    length() %/% 8 + 2

h <- d3$sampleId %>%
    unique() %>%
    length() %/% 2 + 6

idxC <- c("black", "grey", "white")
names(idxC) <- c("dam5p_010", "dam5p_005", "dam5p_low")

d31 <- d3 %>%
    filter(aniRank100 < 2)

pdf(paste("plots/", prefix, ".targets_priority.matrix.pdf", sep = ""), width = w, height = h)
p <- ggplot(d3, aes(x = contigId,
                    y = sampleId,
                    fill = ani,
                    colour = damClass,
                    size = coverageP))
p +
    geom_point(shape = 22) +
    geom_point(shape = 4, color = "black", stroke = 0.25, data = d31) +
    scale_fill_viridis(name = "ANI", option = "C") +
    scale_colour_manual(name = "Damage", values = idxC) +
    scale_size_continuous("Genome coverage", range = c(0.2, 3)) +
    xlab("Contig") +
    ylab("Sample") +
    facet_nested(. ~ genusId + taxNameSpecies + assemblyId, scales = "free_x", space = "free", nest_line = TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          panel.spacing = unit(0.2, "lines"),
          aspect.ratio = 1)
p +
    geom_point(shape = 22) +
    geom_point(shape = 4, color = "black", stroke = 0.25, data = d31) +
    scale_fill_viridis(name = "ANI", option = "C") +
    scale_colour_manual(name = "Damage", values = idxC) +
    scale_size_continuous("Genome coverage", range = c(0.2, 3), trans = "log10") +
    xlab("Contig") +
    ylab("Sample") +
    facet_nested(. ~ genusId + taxNameSpecies + assemblyId, scales = "free_x", space = "free", nest_line = TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          panel.spacing = unit(0.2, "lines"),
          aspect.ratio = 1)
dev.off()



## --------------------------------------------------------------------------------
## plot individual samples coverage vs damage

xBreaks <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

pdf(paste("plots/", prefix, ".targets_priority.pdf", sep = ""), width = 9, height = 7)
for(x in unique(d3$sampleId)) {
    d4 <- d1 %>%
        filter(sampleId == x, contigId == "genome") %>%
        mutate(isAniTop = aniRank100 < 2 & !is.na(aniRank100))

    d41 <- d4 %>%
        filter(isAniTop)

    d42 <- d41 %>%
        filter(flag == "damage_0.1;aniRank100_1;coveragePRatio_0.9;krakenKmerRank_1")

    p <- ggplot(d4, aes(x = coverageP, y = dam5p))
    print(p +
          geom_hline(yintercept = c(0.05, 0.1), size = 0.25, color = c("grey", "black"), linetype = "dashed") +
          geom_point(aes(fill = ani, size = nReads, alpha = isAniTop), colour = "black", stroke = 0.25, shape = 21) +
          geom_point(aes(size = nReads), colour = "black", stroke = 0.25, shape = 4, data = d42) +
          geom_text_repel(aes(label = taxNameSpecies), size = 1.5, segment.size = 0.25, data = d41) +
          scale_fill_viridis(name = "ANI", option = "C") +
          scale_size_continuous("Number of reads", range = c(1, 5)) +
          scale_alpha_manual("ANI rank 1", values = c(0.1, 1)) +
          scale_x_continuous(breaks = xBreaks, trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
          annotation_logticks(sides = "b") +
          coord_cartesian(xlim = c(1e-4, 1), ylim = c(0, 0.5)) +
          xlab("Fraction of genome covered") +
          ylab("5' C>T damage") +
          theme_bw() +
          theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
                panel.grid.minor = element_blank()
                ) +
          ggtitle(x)
          )
}
dev.off()
