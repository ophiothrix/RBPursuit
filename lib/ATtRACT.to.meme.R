### Script to convert ATtRACT motf pwm to a MEME compatible format ###

ATtRACT.to.meme <- function(pwm.file, motif.db.file, organism) {
    # pwm.file <- "./data/ATtRACT/pwm.txt"
    # motif.db.file <- "./data/ATtRACT/ATtRACT_db.txt"
    # organism <- "Homo_sapiens"

    ## Load RBP database
    motif.db <- read.table(motif.db.file, sep = "\t", header = T, stringsAsFactors = F)
    head(motif.db)
    
    ## Remove spaces from gene names
    motif.db$Gene_name <- gsub(" ", "", motif.db$Gene_name)
    
    # table(motif.db$Organism)
    ## Subset to human RBPs
    # motif.db <- motif.db$Matrix_id[motif.db$Organism == "Homo_sapiens"]
    # human.motifs <- unique(motif.db$Matrix_id[motif.db$Organism == "Homo_sapiens"])
    # pfm.list <- pfm.list[names(pfm.list) %in% human.motifs]
    # 
    
    
    
    
    
    # Extract PFMs from the motif database -------------------------------------
    motif.file <- readLines(pwm.file)
    # Extract motif IDs and gene symbols
    motifs <- as.data.frame(matrix(unlist(strsplit(grep("^>", motif.file, value = T), "\t")), ncol = 2, byrow = T), stringsAsFactors = F)
    colnames(motifs) <- c("MotifID", "length")
    motifs$MotifID <- gsub(">", "", motifs$MotifID)
    motifs$length <- as.numeric(motifs$length)
    
    # Extract PFM for each motif. Get the line of each motif and the line of the subsequent motif and get the lines between them minus padding. Need some extra care with the last motif. Essentially, add a "pretend" motif at the end of the object). Also need to add an extra line to the object, or it cuts the last PFM short.
    pfms <- c(grep("^>", motif.file), length(motif.file)+1)
    pfm.list <- list()
    
    ## Set the gap between the motifs. If no gap - set to 1
    gap <- 1
    for (i in 1:(length(pfms)-1)){
        pfm.list[[motifs$MotifID[i]]] <- matrix(as.numeric(unlist(strsplit(motif.file[(pfms[i]+1):(pfms[i+1]-gap)], "\t"))), ncol = 4, byrow = T)
    }
    
    # Check that the motif length corresponds to the PFM matrices
    if (!all(unlist(lapply(pfm.list, nrow)) == motifs$length)) {
        stop("Motif length is inconsistent between the declared motif annotation and extracted PFM") 
    }
    
    
    ## Subset to the motif from a desired organism
    human.motifs <- unique(motif.db$Matrix_id[motif.db$Organism == organism])
    pfm.list <- pfm.list[names(pfm.list) %in% human.motifs]
    motif.lengths <- unlist(lapply(pfm.list, nrow))
    
    ## Add MEME file header
    meme.header <- c("MEME version 4",
                     "\nALPHABET= ACGT",
                     "\nstrands: + -",
                     "\nBackground letter frequencies (from file `../HOCOMOCO/bkg.txt'):\nA 0.29182 C 0.20818 G 0.20818 T 0.29182")
    
    ## Construct motif headers
    motif.set <- motif.db[match(names(pfm.list), motif.db$Matrix_id),]
    motif.header <- cbind("\nMOTIF", motif.set$Matrix_id, motif.set$Gene_name, paste0("\n\nletter-probability matrix: alength= 4 w= ", motif.lengths))
    rownames(motif.header) <- motif.set$Matrix_id
    
    ## Combine headers and PWMs
    meme.motifs <- do.call(rbind, lapply(names(pfm.list), function(x) rbind(motif.header[x,], pfm.list[[x]])))
    
    ## Write out MEME compatible motif file
    write.table(rbind(t(meme.header), meme.motifs, NA), paste0("./data/ATtRACT/ATtRACT.", organism, ".PWMdb.meme"), row.names = F, col.names = F, quote = F, na = "")
    
}
