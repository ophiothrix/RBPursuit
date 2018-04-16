# rm(list = ls())
# gc()

### Make up a toy GRanges object ###
require(GenomicRanges)
variants <- GRanges(seqnames = "chr10", IRanges(70931473:70931474, width = 1))
variants$REF <- "C"
variants$ALT <- "T"



#### Extract motif damage for ATtRACT RBP motif database ####

### Convert the ATtRACT PWM matrix to a MEME compatible one ###
## Only needs to be done once
source("./lib/ATtRACT.to.meme.R")
ATtRACT.to.meme(pwm.file = "./data/ATtRACT/pwm.txt", motif.db.file = "./data/ATtRACT/ATtRACT_db.txt", organism = "Homo_sapiens")

### Extract damage scores
source("./lib/motif.damage.annotation.R")
variants <- get.damage.scores.direct(database.path = "./data/ATtRACT/ATtRACT.Homo_sapiens.PWMdb.meme", variants = variants, force.pfm.extract = F)

### Extract motif damage for JASPAR_vertebrates TF motif database ###
variants <- get.damage.scores.direct(database.path = "./data/JASPAR_CORE_2014_vertebrates.meme", variants = variants, force.pfm.extract = F)

### Extract motif damage for a specific motif database ###
variants <- get.damage.scores.direct(database.path = "./data/HOCOMOCOv9.meme", variants = variants, force.pfm.extract = F)


