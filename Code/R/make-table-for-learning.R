options(scipen=9999) # we do not want scientific notation

library(data.table)
library(dtplyr)
library(dplyr)
library(tidyr)

read.data <- function(fn){
    DT <- fread(fn, stringsAsFactors=F, header=F)
    my.col.names <- c("chr", "damage.pos", "strand", "damage.match", "length",
                      "flags1", "flags2", "CIGAR1", "flags3", "damage.end", "damage.emd.match",
                      "flags4", "flags5", "CIGAR2", "flags6", "flags7", "flags8")
    setnames(DT, colnames(DT), my.col.names)
    return(DT)
}

make.window.column <- function(dt, windows = TRUE, window = 100){
    window.val <- -log10(window)
    dt <- dt[chr != "MT"] %>% setDF() %>% filter(flags2 > 49, flags7 == "-") %>%
                        mutate(window = paste(chr, ":", as.character(round(damage.pos, digits = window.val)), ":", strand, sep = ""), unique.pos = paste(chr, damage.pos, strand, sep=":")) %>%
			select(window, damage.pos, strand, chr)
    setDT(dt)
    return(dt)
}

summarise.dimer <- function(dt){
     tmp.dt <- dt[, .(chr = unique(chr), strand = unique(strand), 
                      damage.min = min(damage.pos), damage.max = max(damage.pos), 
		      count = .N), by = .(window)] %>% setDF() %>%
		      unique()
    setDT(tmp.dt)
    return(tmp.dt)
}

make.sample.table <- function(dname){
    # make a table for one sample, that is merged along the window column. Start
    # with PyPy, and then PuPu, PuPy, and PyPu.

    fn.pref <- c("PyPy", "PuPu", "PuPy", "PyPu")
    sample.dt <- data.frame()
    for(i in 1:length(fn.pref)){
        input.fn <- paste(dname, "/", fn.pref[i], "-single-recurrent.txt", sep="")
	print(input.fn) # for debugging
	dt <- make.window.column(read.data(input.fn))
	dt <- summarise.dimer(dt)
	dt <- dt %>% select(window, count)
	if(i == 1){
	    sample.dt <- dt
	    setDT(sample.dt)
	    setnames(sample.dt, colnames(sample.dt), c("window", paste(dname, fn.pref[i], sep=".")))
	} else {
            sample.dt <- merge(sample.dt, dt, by = "window", all=T)
	    # setDT(sample.dt)
	    setnames(sample.dt, colnames(sample.dt), c("window", paste(dname, fn.pref[1:i], sep=".")))
	}
    }
    sample.dt[is.na(sample.dt)] <- 0 # take care of the introduced NA
    print(dim(sample.dt)) # for debugging
    return(sample.dt)
}

make.all.table <- function(s.names){
    all.tab <- data.frame()
    for(i in 1:length(s.names)){
        sample.table <- make.sample.table(s.names[i])
	if(i == 1){
            all.tab <- sample.table
	} else {
            all.tab <- merge(all.tab, sample.table, by = "window", all = T)
	    print(dim(all.tab)) # for debugging
	}
    }
    all.tab[is.na(all.tab)] <- 0
    print(dim(all.tab)) # for debugging
    return(all.tab)
}
