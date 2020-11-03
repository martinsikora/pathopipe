
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


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]
files <- args[-1]


## --------------------------------------------------------------------------------
## process tables and summarise

flags <- c("editDist_low", "damage_0.1", "damage_0.1;editDist_low")
d <- foreach(f = files) %do% {
    r <- read_tsv(f)
    r
}
d <- bind_rows(d)
write_tsv(d, path = paste("tables/", prefix, ".summary.tsv.gz", sep = ""))

d1 <- filter(d, flag %in% flags)
w <- d$contigId %>% unique() %>% length() %/% 5 + 2
h <- d$sampleId %>% unique() %>% length() %/% 4 + 2

idxC <- c("grey", "orange", "red")
names(idxC) <- flags
pdf(paste("plots/", prefix, ".hits.pdf", sep = ""), width = w, height = h)
p <- ggplot(d1, aes(x = contigId, y = sampleId, fill = flag, colour = flag, size = nReads))
p + geom_point(shape = 22)  + scale_fill_manual(name = "Flag", values = idxC) + scale_colour_manual(name = "Flag", values = idxC) + scale_size_continuous("Number of reads", range = c(1, 4), trans = "log10") + xlab("Contig") + ylab("Sample") + facet_grid(. ~ taxNameSpecies, scales = "free_x", space = "free") + theme_bw() + theme(panel.grid.major = element_line(linetype = "dotted", size = 0.25), panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 6), panel.spacing = unit(0, "lines"), aspect.ratio = 1)
dev.off()


