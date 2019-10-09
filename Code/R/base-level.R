options(scipen=9999) # we do not want scientific notation
options(stringsAsFactors=F)

library(plyr)
library(data.table) # read and manipulate large files
library(dtplyr)
library(dplyr) # easily summarise data frames.
library(tidyr) # we are going to use the spread function later
library(VennDiagram)
library(RColorBrewer)
library(stringr) # we are going to do something with string manipulation
library(ggplot2)
library(pheatmap)


read.data <- function(fn){
    DT <- fread(fn, stringsAsFactors=F, header=F)
    my.col.names <- c("chr", "damage.pos", "strand", "damage.match", "length",
                      "flags1", "flags2", "CIGAR1", "flags3", "damage.end", "damage.emd.match",
                      "flags4", "flags5", "CIGAR2", "flags6", "flags7", "flags8")
    setnames(DT, colnames(DT), my.col.names)
    return(DT)
}

make.count.table <- function(dt){
    # WE are going to make a count table.  Ultimately we want only two columns
    # in this table, a count column and a id column which we are going to merge
    # on.  I could be open to having to merge on 2 columns instead of 1, lets
    # see how we get about that.

    # furst select columns that we are going to work with.
    use.dt <- dt[chr != "MT"][flags2 > 49][flags7 == "-"]  %>%  # use strict matching, and no extra matches
              select(chr, damage.pos, strand, damage.match)
    use.dt$damage.match <- unlist(lapply(use.dt$damage.match, FUN=function(x){return(paste0(unlist(strsplit(x, ""))[1:21], collapse=""))}))
    my.counts <- plyr::count(use.dt, vars = c("chr", "damage.pos", "strand", "damage.match"))
    my.counts$dinuc <- unlist(lapply(my.counts$damage.match, FUN = function(x){return(paste0(unlist(strsplit(x, ""))[10:12], collapse=""))}))
    my.counts$dinuc <- gsub("\\*", "", my.counts$dinuc) # replace the * with nothing
    setcolorder(my.counts, c("chr", "damage.pos", "strand", "damage.match", "dinuc", "freq"))
    return(my.counts)
}

make.count.table.simple <- function(dt){
    # WE are going to make a count table.  Ultimately we want only two columns
    # in this table, a count column and a id column which we are going to merge
    # on.  I could be open to having to merge on 2 columns instead of 1, lets
    # see how we get about that.

    # furst select columns that we are going to work with.
    use.dt <- dt[chr != "MT"][flags2 > 49][flags7 == "-"]  %>%  # use strict matching, and no extra matches
              select(chr, damage.pos, strand)
    my.counts <- plyr::count(use.dt, vars = c("chr", "damage.pos", "strand"))
    return(my.counts)
}

.make.substring <- function(v){
    # take a vector and make it into a substring.  While we are at it, we will
    # just make it into two parts, the first 20 bases, and the damage
    # dinucleotide itself.

    my.subs <- unlist(lapply(v, FUN=function(x){
                                  return(
				    paste0(unlist(
				      strsplit(x, ""))[1:21], collapse=""))}))
    print(head(my.subs)) # for debugging

    return(my.subs)
}

get.saturation <- function(fl, suff = "PyPy-single-recurrent.txt"){
    nrows.file <- c()
    total.damages <- c()
    unique.damages <- c()
    all.count.table <- data.frame()
    for(i in 1:length(fl)){
	print(my.sat.df) # for debugging and progress reporting
        fn <- paste(fl[i], "/", suff, sep="")
	my.df <- read.data(fn)
	count.table <- make.count.table(my.df)
	colnames(count.table) <- c("chr", "position", "strand", "damage.seq", "dinuc", fl[i])
	nrows.file <- c(nrows.file, nrows(my.df))
	unique.damages <- c(unique.damages, nrow(count.table))
	if(i == 1){
	     all.count.table <- count.table
	} else {
	    all.count.table <- merge(all.count.table, count.table, by = c("chr", "position", "strand", "damage.seq", "dinuc"), all = T)
	}
        total.damages <- c(total.damages, nrow(all.count.table))
    }
    my.sat.df <- data.frame(sample = fl, total.reads = nrows.file, unique.damages = unique.damages, total.damages = total.damages)
    return(my.sat.df)
}

make.merged.table <- function(fl, suff = "PyPy-single-recurrent.txt"){
    # fl is the file list
    my.merged.table <- data.frame()
    for(i in 1:length(fl)){
        print(head(my.merged.table)) # for debugging
        fn <- paste(fl[i], "/", suff, sep="")
	# print(fn) # for debugging
        count.table <- make.count.table(read.data(fn))
	exp.name <- fl[i]
	# print(exp.name) # for debugging
	colnames(count.table) <- c("chr", "position", "strand", "damage.seq", "dinuc", exp.name)
	if(i == 1){
            my.merged.table <- count.table
	} else {
	    my.merged.table <- merge(my.merged.table, count.table, by = c("chr", "position", "strand", "damage.seq", "dinuc"), all=T)
	}
    }
    print(dim(my.merged.table)) # for debugging
    my.merged.table[is.na(my.merged.table)] <- 0
    print(head(my.merged.table)) # for debugging

    row.sum <- rowSums(my.merged.table[, 6:ncol(my.merged.table)])
    my.merged.table$row.sum <- row.sum
      
    my.merged.table <- my.merged.table %>%
                       filter(row.sum > 0)
    
    # damage.context <- unlist(lapply(my.merged.table$damage.seq, FUN=function(x){return(gsub("\\*", "", x))}))
    # damage.dinuc   <- unlist(lapply(my.merged.table$damage.seq, FUN=function(x){return(paste0(unlist(strsplit(x, ""))[10:12], collapse=""))}))
    # damage.dinuc <- gsub("\\*", "", damage.dinuc)
    
    # my.merged.table$damage.context <- damage.context
    # my.merged.table$damage.dinuc   <- damage.dinuc

    # my.merged.table <- my.merged.table[order(as.numeric(my.merged.table$chr), my.merged.table$position), ]
    return(my.merged.table)
}

make.sub.merged <- function(df, samples){
    common.col.indices <- c(1, 2, 3, 4, which(colnames(df) %in% c("damage.context", "damage.dinuc")))
    ch.indices <- which(colnames(df) %in% samples)
    common.sub <- df %>% select(common.col.indices)
    samples.sub <- df %>% select(ch.indices)
    samples.sub$row.sum <- rowSums(samples.sub)

    put.together.df <- cbind(common.sub, samples.sub)
    put.together.df <- put.together.df %>% filter(row.sum > 1)
    return(put.together.df)
}


make.nreads.from.merged <- function(d, my.index, s.name){

}


make.sgr.from.merged.table <- function(d, my.index, s.name){
    # for each sample generate a SGR file, the dinucleotide counts and total
    # number of reads involved.
    # get the columns
    ofl <- list()
    use.index <- 5 + my.index
    sub.df <- d[, c(1:3, 5, use.index), with = FALSE]
    setnames(sub.df, s.name, "freq")
    sub.df <- sub.df %>% filter(freq > 0)
    n.reads <- sum(sub.df$freq)
    print(n.reads) 
    dinuc.counts <- sub.df[, .(counts = sum(freq)), by = .(dinuc)]
    setnames(dinuc.counts, "counts", s.name)
    print(dinuc.counts)
    print(head(sub.df))
    sub.df <- sub.df %>% mutate(window.id = round(position, 2))
    sgr.df <- sub.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
              transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
    if("TG" %in% dinuc.counts$dinuc){
        tg.df <- sub.df %>% filter(dinuc == "TG")
	tg.sgr <- tg.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
              transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
	ofl[["sgr.TG"]] <- tg.sgr
    }
    if("CG" %in% dinuc.counts$dinuc){
        pyA.df <- sub.df %>% filter(grepl("*A", dinuc))
	pyA.sgr <- pyA.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
              transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
	ofl[["sgr.pyA"]] <- pyA.sgr
    }
    print(sgr.df)
    ofl[["nreads"]] <- n.reads
    ofl[["dinuc.counts"]] <- dinuc.counts
    ofl[["sgr.df"]] <- sgr.df

    return(ofl)
}

make.sgr.from.list <- function(count.fn.l, dir.name){  
    print(count.fn.l) # for debugging
    op.l <- list()
    all.dinuc.counts <- data.frame()
    for(j in 1:length(count.fn.l)){
        d <- fread(count.fn.l[j], header=T)
        fl <- colnames(d)[6:(ncol(d)-1)]
	if(length(fl) > 33) fl <- fl[1:33]
	print(fl) # for debugging
        if(unlist(strsplit(count.fn.l[j], "\\."))[2] == "PyPy"){
            dinuc.suff <- "PyPy"
	} else {
            dinuc.suff <- unlist(strsplit(count.fn.l[j], "\\."))[3]
	}
	nn.dinuc.df <- data.frame()
        for(i in 1:length(fl)){
            ofl <- make.sgr.from.merged.table(d, i, fl[i])
	    print(names(ofl)) # for debugging
            suffs <- names(ofl)
	    ofn <- paste(paste(fl[i], dinuc.suff, sep="."), suffs, sep=".")
	    ofn <- file.path(dir.name, ofn)
	    print(ofn) # for debugging
	    
	    for(k in 1:length(ofl)){
                
	    }
	    write.table(ofl[["sgr.df"]], ofn[1], row.names=F, col.names=F, quote=F, sep="\t")
	    if("TG" %in% names(ofl))  write.table(ofl[["sgr.TG"]], ofn[2], row.names=F, col.names=F, quote=F, sep="\t")
	    if("pyA" %in% names(ofl))  write.table(ofl[["sgr.pyA"]], ofn[2], row.names=F, col.names=F, quote=F, sep="\t")
	    
	    if(i == 1){
                nn.dinuc.df <- ofl[["dinuc.counts"]]
	    } else {
                nn.dinuc.df <- merge(nn.dinuc.df, ofl[["dinuc.counts"]], by = "dinuc", all = T)
	    }
        }
	if(j == 1){
            all.dinuc.counts <- nn.dinuc.df
	} else {
	    if(ncol(nn.dinuc.df) > 33){
                all.dinuc.counts <- rbind(all.dinuc.counts, nn.dinuc.df[, 1:33])
	    } else if(ncol(nn.dinuc.df) == 3){
	        all.dinuc.counts <- rbind(all.dinuc.counts, nn.dinuc.df[, 34:36])
	    } else {
                all.dinuc.counts <- rbind(all.dinuc.counts, nn.dinuc.df)
	    }
	}
	all.dinuc.counts[is.na(all.dinuc.counts)] <- 0
	op.l[["all.dinucs"]] <- all.dinuc.counts
	print(all.dinuc.counts)
    }
    return(op.l)
}

#### Based on the the discussion about normalization.  I still have a few
#    doubts but we will figure those out once I have the basic skeleton ready.
#    To the best of my understand, I first need to iterate over all the
#    different combinations per sample, normalize the counts based on the factor
#    I calculate, and then merge all such data-sets together.  The threshold
#    will be rowSum > 0.

make.fn <- function(s.name, suff = "PyPy-single-recurrent.txt"){
   return(file.path(s.name, suff))
}

make.other.counts.vec <- function(dt){
    my.dinuc <- unlist(lapply(dt$damage.match,
                                     FUN = function(x){
                                         return(paste0(unlist(strsplit(x, ""))[10:12], collapse=""))
				     }))
    my.dinuc <- gsub("\\*", "", my.dinuc)
    no.N <- grep("N", my.dinuc, invert = T)
    my.dinuc <- my.dinuc[no.N]
    return(as.data.frame(table(my.dinuc)))
}

make.other.counts <- function(s.name, suff = "PyPu-single-recurrent.txt"){
    dinuc.counts <- data.frame()
    for(i in c("PyPu", "PuPy", "PuPu", "PyPy")){
        suff <- "-single-recurrent.txt"
	suff <- paste(i, suff, sep="")
        my.other.counts <- make.other.counts.vec(read.data(make.fn(s.name, suff = suff)))
	print(my.other.counts) # for debugging
        dinuc.counts <- rbind(dinuc.counts, my.other.counts)
    }
    return(dinuc.counts)
}

make.all.dinucs <- function(fl){
    all.dinucs <- data.frame()
    for(i in 1:length(fl)){
        sname <- fl[i]
	print(sname) # to track progress
	dinuc.count.df <- make.other.counts(sname)
        colnames(dinuc.count.df) <- c("dinuc", sname)
	print(dinuc.count.df) #for debugging
	if(i == 1){
            all.dinucs <- dinuc.count.df
	} else {
            all.dinucs <- merge(all.dinucs, dinuc.count.df, by = "dinuc", all=T)
	    print(nrow(all.dinucs))
	}
    }
    rownames(all.dinucs) <- all.dinucs$dinuc
    all.dinucs[, 1] <- NULL
    all.dinucs <- as.matrix(all.dinucs)
    all.dinucs <- t(all.dinucs)
    return(as.data.frame(all.dinucs))
}

make.normalized.merged.table <- function(fn){
    pypy.counts.table <- get.pypy.counts.table(fn)

}

