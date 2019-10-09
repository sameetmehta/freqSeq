library(data.table)
library(tidyverse)

# We are going to calculate how many dimers are seen as we inceraase the total
# number of experiments under consideration.  We are going to use Karl's filters
# and then see how many locations are seen.  We want to make sure that across
# all the samples under consideration, each position is counted exactly once.

read.data <- function(fn){
    DT <- fread(fn, stringsAsFactors=F, header=F)
    my.col.names <- c("chr", "damage.pos", "strand", "damage.match", "length",
                      "flags1", "flags2", "CIGAR1", "flags3", "damage.end", "damage.end.match",
                      "flags4", "flags5", "CIGAR2", "flags6", "flags7", "flags8")
    setnames(DT, colnames(DT), my.col.names)
    df <- DT %>% select(1:3, 9, 10, 15:17, 4)
    df
}

get.dat.filtered <- function(fn){
    
    dat <- read.data(fn)
    
    setDF(dat)

    perfect1 <- (dat[,4] == "MD:Z:58")
    perfect1 <- perfect1 & !grepl("58M,0", dat[,7])
    err1 <- grepl("MD:Z:[^0][0-9]*[ACGT][0-9]*$", dat[,4])
    err1 <- err1 & !grepl("58M,[01]", dat[,7])
    dat <- dat[perfect1 | err1,,drop=FALSE]
    
    perfect2 <- (dat[,6] == "MD:Z:66")
    perfect2 <- perfect2 & !grepl("66M,0", dat[,8])
    err2 <- grepl("MD:Z:[^0][0-9]*[ACGT][0-9]*$", dat[,6])
    err2 <- err2 & !grepl("66M,[01]", dat[,8])
    dat <- dat[perfect2 | err2,,drop=FALSE]
    
    unduplicated <- !duplicated(dat[,c(2,5)])
    dat <- dat[unduplicated,,drop=FALSE]
}

make.dimer <- function(vec){
    unlist(lapply(vec, function(x)process.str(x)))
}

process.str <- function(s){
   s <- gsub("\\*", "", s)
   s <- unlist(strsplit(s, ""))[10:11] %>%
        paste0(collapse = "")
}

get.saturation <- function(fl){
    # get the saturation for all the dimers form the given file list.
    for(i in seq_along(fl)){
        dat <- get.dat.filtered(fl[i])
	tot.type <- nrow(dat)
	
    }
}
