
#########################
## plot edit distances ##
#########################

## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(Rsamtools))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(foreach))
suppressMessages(library(viridis))


## --------------------------------------------------------------------------------
## command line arguments

args <- commandArgs(trailingOnly = TRUE)
dbPath <- args[1]
outPdf <- args[2]
outTsv <- args[3]
bamfiles <- args[-1:-3]


## --------------------------------------------------------------------------------
## taxInfo

seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccccc")


## --------------------------------------------------------------------------------
## get edit distance and write table

d <- foreach(f1 = bamfiles) %do% {
    p1 <- ScanBamParam(tag=c("NM"), flag = scanBamFlag(isDuplicate = FALSE))
    r1 <- scanBam(f1, param=p1)
    r2 <- r1[[1]]$tag$NM %>% as_tibble()

    if(nrow(r2) == 0){
        NULL
    } else {
        r3 <- gsub(".*\\/", "", f1) %>% strsplit("\\.") %>% unlist()
        assemblyId <- paste(r3[2:3], collapse = ".")
        r4 <- count(r2, value) %>% rename(editDist = value) %>% mutate(p = n / sum(n), assemblyId = assemblyId, taxNameSpecies = seqInfo$taxNameSpecies[match(assemblyId, seqInfo$assemblyId)]) %>% select(assemblyId, taxNameSpecies, editDist, n, p)
        r4
    }
}
d <- bind_rows(d)
write_tsv(d, path = outTsv)


## --------------------------------------------------------------------------------
## plot

nBam <- length(bamfiles)

idxC <- viridis(nBam)
idxS <- rep(0:14, length.out = nBam)
th <- theme_bw() + theme(strip.background = element_blank(), panel.grid.major = element_line(linetype = "dotted", size = 0.25), legend.key.size = unit(0.015, "npc"), legend.text = element_text(size = 3), panel.grid.minor = element_blank())

if(nrow(d) == 0){
    pdf(outPdf, width = 10, height = 4)
    dev.off()
} else {
    d <- mutate(d, lab = paste(assemblyId, taxNameSpecies, sep = " | "))
    pdf(outPdf, width = 14, height = 5)
    p <- ggplot(d, aes(x = editDist, y = p, fill = lab, colour = lab, shape = lab))
    print(p + geom_line(size = 0.25, show.legend = FALSE) + geom_point(size = 1) + xlab("Edit distance") + ylab("Fraction of reads") + scale_shape_manual(values = idxS) + scale_colour_manual(values = idxC) + th + guides(shape = guide_legend(ncol = 4, override.aes = list(size = 0.75, stroke = 0.1))))
    print(p + geom_line(size = 0.25, show.legend = FALSE) + geom_point(size = 1) + xlab("Edit distance") + ylab("Fraction of reads") + scale_shape_manual(values = idxS) + scale_colour_manual(values = idxC) + scale_y_log10() + th + guides(shape = guide_legend(ncol = 4, override.aes = list(size = 0.75, stroke = 0.1))))
    dev.off()
}

