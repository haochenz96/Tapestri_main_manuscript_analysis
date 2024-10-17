# devtools::install_github('mskilab/GxG')
library(GxG)
## Bioconductor string libraries
library(Biostrings)
library(BSgenome)
library(Rsamtools)
## allows us to read in 2bit files
library(rtracklayer)

# load genome
library(BSgenome.Hsapiens.1000genomes.hs37d5)
genome <- BSgenome.Hsapiens.1000genomes.hs37d5
######################
# M04 TRIM24-BRAF fusion breakpoints:
# chr7:138189929
# chr7:140481953
# M04 CDKN2A deletion breakpoints:
# chr9:21961008
# chr9:21981193
######################

## define some intervals on the mitochondrial genomem
wins = GRanges(c('7:138189929', '7:140481953 ')) + 100
wins = GRanges(c('9:21961008', '9:21981193')) + 100
seq = getSeq(genome, wins)
gm = homeology(seq, stride = 2, pad = 9, verbose = TRUE, rc  = FALSE)
plot(gm$gt, gm$footprint)

options(repr.plot.width = 8, repr.plot.height = 10)
plot(c(gm$gtrack(clim = c(0, 1), cmap.min = 0)), wins + 100)
