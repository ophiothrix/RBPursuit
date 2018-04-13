# rm(list = ls())
# gc()

### Make up a toy GRanges object ###
require(GenomicRanges)
variants <- GRanges(seqnames = "chr10", IRanges(70931473:70931474, width = 1))
variants$REF <- "C"
variants$ALT <- "T"

### Convert the ATtRACT PWM matrix to a MEME compatible one ###
## Only needs to be done once
source("./lib/ATtRACT.to.meme.R")
ATtRACT.to.meme(pwm.file = "./data/ATtRACT/pwm.txt", motif.db.file = "./data/ATtRACT/ATtRACT_db.txt", organism = "Homo_sapiens")


### Extract a list of PWMs ###
## Only needs to be done once
source("./lib/motif.damage.annotation.R")
pfm.list <- extract.pfm(database.path = "./data/ATtRACT/ATtRACT.Homo_sapiens.PWMdb.meme", force.run = T)

### Extract motif damage for a specific motif database ###
variants <- get.damage.scores.direct(database.path = "./data/ATtRACT/ATtRACT.Homo_sapiens.PWMdb.meme", variants = variants)

variants
