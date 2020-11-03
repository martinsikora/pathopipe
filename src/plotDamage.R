
############################
## plot coverage evenness ##
############################


## --------------------------------------------------------------------------------
## libraries

suppressMessages(library(Rsamtools))
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
dbPath <- args[1]
outPdf <- args[2]
outTsv <- args[3]
infiles <- args[-1:-3]


## --------------------------------------------------------------------------------
## taxInfo

seqInfo <- read_tsv(paste(dbPath, "/library.seqInfo.tsv", sep = ""), col_types = "cccccccc")


## --------------------------------------------------------------------------------
## get damage data

st <- c("C>T", "G>A", "A>G", "T>C", "A>C", "A>T", "C>G", "C>A", "T>G", "T>A", "G>C", "G>T", "A>-", "T>-", "C>-", "G>-", "->A", "->T", "->C", "->G", "S")
ref <- c("A", "C", "G", "T")

d <- foreach(f1 = infiles) %do% {
    r <- read_tsv(f1, comment = "#", col_types = "cccddddddddddddddddddddddddddd")
    r1 <- r %>% pivot_longer(all_of(c(st, ref))) %>% mutate(refAllele = substr(name, 1, 1)) %>% filter(refAllele %in% c("A", "C", "G", "T", "S"))
    r21 <- filter(r1, refAllele != "S") %>% group_by(Pos, End, name, refAllele) %>% summarise(value = sum(value)) %>% group_by(Pos, End, refAllele) %>% mutate(p = value / value[name == refAllele]) %>% ungroup()
    r22 <- filter(r1, refAllele == "S") %>% group_by(Pos, End, name, refAllele)  %>% summarise(value = sum(value), Total = sum(Total)) %>% mutate(p = value / Total)
    r2 <- bind_rows(r21, r22)
    r3 <- gsub("mapdamage/", "", f1) %>% strsplit("\\/") %>% map_chr(2) %>% strsplit("\\.") %>% unlist()
    assemblyId <- paste(r3[2:3], collapse = ".")
    mutate(r2, assemblyId = assemblyId, taxNameSpecies = seqInfo$taxNameSpecies[match(assemblyId, seqInfo$assemblyId)])
}
d <- bind_rows(d)
d$Pos[d$End == "3p"] <- -d$Pos[d$End == "3p"]
d <- filter(d, abs(Pos) <= 25)
write_tsv(d, path = outTsv)


## --------------------------------------------------------------------------------
## plot

d$End <- factor(d$End, levels = c("5p", "3p"))
d$name <- factor(d$name, levels = st)
d <- mutate(d, lab = paste(assemblyId, taxNameSpecies, sep = "\n"))

idxSt <- c("red", "blue", rep("grey", length(st) - 3), "orange")
names(idxSt) <- st

h <- 2 + length(infiles)
pdf(outPdf, width = 10, height = h)
p <- ggplot(d, aes(x = Pos, y = p, colour = name))
print(p + geom_line(size = 0.25) + facet_grid(lab ~ End, scales = "free") + xlab("Position from read end") + ylab("Frequency") + coord_cartesian(ylim = c(0, 0.3)) + scale_color_manual(values = idxSt) + theme_bw() + theme(strip.background = element_blank(), panel.grid.major = element_line(linetype = "dotted", size = 0.25), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0), legend.position = "none"))
dev.off()
