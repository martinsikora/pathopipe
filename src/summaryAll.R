
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
    r <- read_tsv(f, col_types = "ccccccdddddddddddddddddddddddddddddddc", na = c("", "NA", "NaN"))
    r
}
d <- bind_rows(d)

tP <- read_tsv(targetsPriority, col_names = FALSE)
d1 <- d %>%
    filter(genusId %in% tP$X1)


## --------------------------------------------------------------------------------
## plot matrix of candidates with krakenKmerRank < 3 and damage

## find unique species name with krakenKmerRank < 3
idxS <- d1 %>%
    filter(krakenKmerRank < 3,
           dam5p >= 0.1,
           !is.na(aniRank100)) %>%
    pull(taxNameSpecies) %>%
    unique()

## set up plot data
d3 <- d %>%
    filter(taxNameSpecies %in% idxS) %>%
    mutate(damClass = case_when(
               dam5p >= 0.1 ~ "dam5p_010",
               dam5p >= 0.05 ~ "dam5p_005",
               TRUE ~ "dam5p_low"
           ),
           ani1 = case_when(
               ani >= 0.95 ~ ani,
               TRUE ~ NaN
           )
           ) 

d32 <- d3 %>%
    distinct(contigId, contigL) %>%
    arrange(desc(contigL))

d3 <- d3 %>%
    mutate(contigId = factor(contigId, levels = unique(d32$contigId)),
           sampleId = factor(sampleId, levels = unique(d3$sampleId)),
           y = as.integer(sampleId))

h <- d3$sampleId %>%
    unique() %>%
    length() %/% 2 + 6

th <- theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
            axis.text.y = element_text(size = 6),
            panel.spacing = unit(0.2, "lines"),
            aspect.ratio = 1)

idxC <- c("black", "grey", "white")
names(idxC) <- c("dam5p_010", "dam5p_005", "dam5p_low")

yLab <- levels(d3$sampleId)
yBreaks <- 1:length(yLab)

for(g in unique(d3$genusId)){

    d4 <- filter(d3, genusId == g)
    d41 <- d4 %>%
        filter(krakenKmerRank < 3,
               !is.na(aniRank100))

    w <- d4$contigId %>%
        unique() %>%
        length() %/% 8 + 4
    
    
    pdf(paste("plots/targets_priority/", g, ".", prefix, ".targets_priority.matrix.pdf", sep = ""), width = w, height = h)
    p <- ggplot(d4, aes(x = contigId,
                        y = y,
                        fill = ani1,
                        colour = damClass,
                        size = coverageP))
    print(p +
        geom_point(shape = 22) +
        geom_point(shape = 4, color = "black", stroke = 0.25, data = d41) +
        scale_fill_viridis(name = "ANI", option = "C") +
        scale_colour_manual(name = "Damage", values = idxC) +
        scale_size_continuous("Genome coverage", range = c(0.2, 3)) +
        scale_y_continuous(name = "Sample", breaks = yBreaks, labels = yLab, sec.axis = dup_axis()) +
        xlab("Contig") +
        facet_nested(. ~ genusId + taxNameSpecies + assemblyId, scales = "free_x", space = "free", nest_line = TRUE) +
        theme_bw() +
        th)

    print(p +
        geom_point(shape = 22) +
        geom_point(shape = 4, color = "black", stroke = 0.25, data = d41) +
        scale_fill_viridis(name = "ANI", option = "C") +
        scale_colour_manual(name = "Damage", values = idxC) +
        scale_size_continuous("Genome coverage", range = c(0.2, 3), trans = "log10") +
        scale_y_continuous(name = "Sample", breaks = yBreaks, labels = yLab, sec.axis = dup_axis()) +
        xlab("Contig") +
        facet_nested(. ~ genusId + taxNameSpecies + assemblyId, scales = "free_x", space = "free", nest_line = TRUE) +
        theme_bw() +
        th)
    dev.off()
}


## --------------------------------------------------------------------------------
## plot individual samples coverage vs damage

xBreaks <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

pdf(paste("plots/targets_priority/", prefix, ".targets_priority.pdf", sep = ""), width = 9, height = 7)
for(x in unique(d3$sampleId)) {
    d4 <- d1 %>%
        filter(sampleId == x,
               contigId == "genome") %>%
        mutate(isTop = krakenKmerRank < 3 & !is.na(aniRank100) & !is.na(krakenKmerRank) & coveragePRatioCorr >= 0.5,
               ani1 = case_when(
               ani >= 0.95 ~ ani,
               TRUE ~ NaN)
               )

    d41 <- d4 %>%
        filter(isTop)

    d42 <- d41 %>%
        filter(flag == "damage_0.1;aniRank100_1;coveragePRatioCorr_0.8;krakenKmerRank_1")

    p <- ggplot(d4, aes(x = coverageP, y = dam5p))
    print(p +
          geom_hline(yintercept = c(0.05, 0.1), size = 0.25, color = c("grey", "black"), linetype = "dashed") +
          geom_point(aes(fill = ani1, size = nReads, alpha = isTop), colour = "black", stroke = 0.25, shape = 21) +
          geom_point(aes(size = nReads), colour = "black", stroke = 0.25, shape = 4, data = d42) +
          geom_text_repel(aes(label = taxNameSpecies), size = 1.5, segment.size = 0.25, data = d41, max.overlaps = Inf) +
          scale_fill_viridis(name = "ANI", option = "C") +
          scale_size_continuous("Number of reads", range = c(1, 5)) +
          scale_alpha_manual("Kraken kmer rank < 3\ncoverageP ratio > 0.5", values = c(0.2, 1)) +
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


## --------------------------------------------------------------------------------
## write summary tables

write_tsv(d, file = paste("tables/", prefix, ".summary.tsv.gz", sep = ""))
write_tsv(d1, file = paste("tables/", prefix, ".targets_priority.summary.tsv.gz", sep = ""))
