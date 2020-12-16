
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
suppressMessages(library(philentropy))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
files <- args[-1]


## --------------------------------------------------------------------------------
## process tables and summarise

d <- foreach(f = files) %do% {
    r <- read_tsv(f)
    r
}
d <- bind_rows(d)
write_tsv(d, path = paste("tables/", prefix, ".summary.tsv.gz", sep = ""))

flags <- c("editDistAvg_1.5", "editDistDecay_5", "editDistAvg_1.5;editDistDecay_5", "damage_0.1", "damage_0.1;editDistAvg_1.5", "damage_0.1;editDistDecay_5",  "damage_0.1;editDistAvg_1.5;editDistDecay_5")
idxC <- c("grey80", "grey50", "grey30", "gold", "orange1", "orange3", "red")
names(idxC) <- flags

idxS <- filter(d, flag %in% flags[-1:-3]) %>% .$taxNameSpecies %>% unique()
d1 <- filter(d, taxNameSpecies %in% idxS)
w <- d$contigId %>% unique() %>% length() %/% 20 + 2
h <- d$sampleId %>% unique() %>% length() %/% 2 + 2

d1$flag <- factor(d1$flag, levels = flags)
pdf(paste("plots/", prefix, ".hits.pdf", sep = ""), width = w, height = h)
p <- ggplot(d1, aes(x = contigId, y = sampleId, fill = flag, colour = flag, size = nReads))
p + geom_point(shape = 22)  + scale_fill_manual(name = "Flag", values = idxC) + scale_colour_manual(name = "Flag", values = idxC) + scale_size_continuous("Number of reads", range = c(1, 4), trans = "log10") + xlab("Contig") + ylab("Sample") + facet_grid(. ~ taxNameSpecies, scales = "free_x", space = "free") + theme_bw() + theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25), panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 6), panel.spacing = unit(0, "lines"), aspect.ratio = 1)
dev.off()

##
d2 <- filter(d, flag == "damage_0.1;editDistAvg_1.5;editDistDecay_5") %>% mutate(flag1 = 1L) %>% select(sampleId, contigId, flag1) %>% pivot_wider(names_from = contigId, values_from = flag1) %>% column_to_rownames("sampleId") %>% as.matrix()
d2[is.na(d2)] <- 0

d3 <- distance(d2, method = "jaccard")
d4 <- cmdscale(d3) %>% as_tibble() %>% mutate(sampleId = rownames(d2))

pdf(paste("plots/", prefix, ".mds.pdf", sep = ""), width = 8, height = 8)
p <- ggplot(d4, aes(x = V1, y = V2))
p + geom_point(shape = 21, color = "steelblue", fill = "steelblue", alpha = 0.9) + geom_text_repel(aes(label = sampleId), size = 1, segment.size = 0.25) + xlab("MDS 1") + ylab("MDS2") + theme_bw() + theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25), panel.grid.minor = element_blank(), strip.background = element_blank())
dev.off()
