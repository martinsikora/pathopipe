
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


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
targetsPriority <- args[2]
files <- args[-1:-2]


## --------------------------------------------------------------------------------
## process tables and summarise

d <- foreach(f = files) %do% {
    r <- read_tsv(f, col_types = "ccccccdddddddddddddddddddddddddc")
    r
}
d <- bind_rows(d)
write_tsv(d, path = paste("tables/", prefix, ".summary.tsv.gz", sep = ""))

tP <- read_tsv(targetsPriority, col_names = FALSE)

d1 <- d %>%
    filter(genusId %in% tP$X1)
write_tsv(d1, path = paste("tables/", prefix, ".targets_priority.summary.tsv.gz", sep = ""))

## plot aniRank 1 taxa across prioritize samples
idxS <- d1 %>%
    filter(aniRank < 2) %>%
    pull(taxNameSpecies) %>%
    unique()

d2 <- d %>%
    filter(taxNameSpecies %in% idxS,
           coverageP >= 0.001) %>%
    mutate(damClass = case_when(
               dam5p >= 0.1 ~ "dam5p_010",
               dam5p >= 0.05 ~ "dam5p_005",
               TRUE ~ "dam5p_low"
           ))

w <- d2$contigId %>%
    unique() %>%
    length() %/% 8 + 2

h <- d2$sampleId %>%
    unique() %>%
    length() %/% 2 + 4

idxC <- c("black", "grey", "white")
names(idxC) <- c("dam5p_010", "dam5p_005", "dam5p_low")

d21 <- filter(d2, aniRank < 2)

pdf(paste("plots/", prefix, ".targets_priority.matrix.pdf", sep = ""), width = w, height = h)
p <- ggplot(d2, aes(x = contigId,
                    y = sampleId,
                    fill = ani,
                    colour = damClass,
                    size = coverageP))
p +
    geom_point(shape = 22) +
    geom_point(shape = 4, color = "black", stroke = 0.25, data = d21) +
    scale_fill_viridis(name = "ANI", option = "C") +
    scale_colour_manual(name = "Damage", values = idxC) +
    scale_size_continuous("Genome coverage", range = c(0.2, 2.5)) +
    xlab("Contig") +
    ylab("Sample") +
    facet_grid(. ~ taxNameSpecies, scales = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          panel.spacing = unit(0, "lines"),
          aspect.ratio = 1)
p +
    geom_point(shape = 22) +
    geom_point(shape = 4, color = "black", stroke = 0.25, data = d21) +
    scale_fill_viridis(name = "ANI", option = "C") +
    scale_colour_manual(name = "Damage", values = idxC) +
    scale_size_continuous("Genome coverage", range = c(0.2, 3), limits = c(1e-3, 1), trans = "log10") +
    xlab("Contig") +
    ylab("Sample") +
    facet_grid(. ~ taxNameSpecies, scales = "free_x", space = "free") +
    theme_bw() +
    theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          panel.spacing = unit(0, "lines"),
          aspect.ratio = 1)
dev.off()

xBreaks <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)

pdf(paste("plots/", prefix, ".targets_priority.pdf", sep = ""), width = 9, height = 7)
for(x in unique(d2$sampleId)) {
    d3 <- d2 %>%
        filter(sampleId == x, contigId == "genome")

    d31 <- d3 %>%
        filter(nReads >= 100,
               dam5p >= 0.05,
               coverageP >= 0.001)

    p <- ggplot(d3, aes(x = coverageP, y = dam5p))
    print(p +
          geom_hline(yintercept = c(0.05, 0.1), size = 0.25, color = c("grey", "black"), linetype = "dashed") +
          geom_point(aes(fill = ani, size = nReads), colour = "black", stroke = 0.25, shape = 21) +
          geom_text_repel(aes(label = taxNameSpecies), size = 1.5, segment.size = 0.25, data = d31) +
          scale_fill_viridis(name = "ANI", option = "C") +
          scale_size_continuous("Number of reads", range = c(1, 5)) +
          scale_x_continuous(breaks = xBreaks, trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
          annotation_logticks(sides = "b") +
          coord_cartesian(xlim = c(0.0001, 1), ylim = c(0, 0.5)) +
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
