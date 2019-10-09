options(scipen=9999) # we do not want scientific notation
options(stringsAsFactors=F)

library(tidyverse)
library(plyr)
library(data.table) # read and manipulate large files
library(dtplyr)
library(dplyr) # easily summarise data frames.
library(tidyr) # we are going to use the spread function later
library(stringr) # we are going to do something with string manipulation
library(seqLogo)

read.data <- function(fn){
    DT <- fread(fn, stringsAsFactors=F, header=F)
    my.col.names <- c("chr", "damage.pos", "strand", "damage.match", "length",
                      "flags1", "flags2", "CIGAR1", "flags3", "damage.end", "damage.emd.match",
                      "flags4", "flags5", "CIGAR2", "flags6", "flags7", "flags8")
    setnames(DT, colnames(DT), my.col.names)
    DT
}

get.sub.df <- function(d, my.index, s.name){
    # Create a sub data frame for the selected column.  We only want those
    # entries that are "seen" in that experiment.  
    use.index <- 5 + my.index
    sub.df <- d[, c(1:3, 5, use.index), with = FALSE]
    setnames(sub.df, s.name, "freq")
    sub.df <- sub.df %>% filter(freq > 0)
    return(sub.df)
}

get.saturation.on.windows <- function(fn, window = 100){
    round.no <- -log10(window)
    d <- fread(fn, header = T)
    d <- d %>% mutate(id = paste0(chr, strand, round(position, round.no), sep=""))
    # we are giving it a separator incase we need to re-construct the stuff
    # later.
    d <- d[, c(6:(ncol(d) - 2), ncol(d)), with = FALSE]
    snames <- colnames(d)[1:(ncol(d) - 1)]
    print(snames) # for debugging
    
    unique.windows <- c()
    unique.window.in.sample <- c()
    
    i <- 1
    for(j in 1:length(snames)){
        tmp.df <- d[, c(j, ncol(d)), with = FALSE]
	setnames(tmp.df, snames[j], "freq")
	print(head(tmp.df))
	tmp.df <- tmp.df %>% filter(freq > 0) %>%
	          select(id) %>% unique()
        print(head(tmp.df))
	unique.window.in.sample <- c(unique.window.in.sample, nrow(tmp.df))
        if(i == j){
            sub.df <- d[, c(1, ncol(d)), with = FALSE]
	} else {
            sub.df <- d[, c(i:j, ncol(d)), with = FALSE]
	}
	sub.df <- sub.df %>% mutate(row.sum = rowSums(sub.df[, i:j])) %>%
	          filter(row.sum > 0) %>% select(id) %>% unique()
	print(head(sub.df)) # for debugging
        unique.windows <- c(unique.windows, nrow(sub.df))
    }
    op.df <- data.frame(sample.added = snames,
                        samples = 1:length(snames),
			unique.windows.in.sample = unique.window.in.sample,
			unique.windows = unique.windows)
    return(op.df)
}

get.saturation <- function(fn){
    d <- fread(fn, header=T)
    snames <- colnames(d)[6:(ncol(d) - 1)]
    unique.damages <- c()
    total.damages  <- c()
    unique.damage.in.sample <- c()
    total.damage.in.sample  <- c()
    i <- 1
    for(j in 1:length(snames)){
        print(c(i, j)) # for debugging and tracking progress
	use.index.i <- 5 + i
	use.index.j <- 5 + j
	if(i == j){
	    sub.df <- d[, c(1:3, 5, use.index.i), with = FALSE]
	} else {
            sub.df <- d[, c(1:3, 5, use.index.i:use.index.j), with = FALSE]
	}
	sub.df <- sub.df %>% mutate(row.sum = rowSums(sub.df[, 5:ncol(sub.df)])) %>%
	          filter(row.sum > 0)
	unique.damages <- c(unique.damages, nrow(sub.df))
	total.damages  <- c(total.damages, sum(sub.df$row.sum))
    }
    op.df <- data.frame(sample.added = snames,
                        samples = 1:length(snames), 
                        unique.damages = unique.damages, 
			total.damages  = total.damages,
			unique.damage.perc = (unique.damages / total.damages) * 100,
			recurrent.damage.perc = ((total.damages - unique.damages)/total.damages) * 100)
    return(op.df)
}

make.dinucs.from.merged <- function(d){
    # Requires that each column start with "SP".  But Sanjay is generally naming
    # his samples that way so it should not matter.
    return(d %>% setDF() %>% group_by_("dinuc") %>%
           summarise_at(.cols = vars(starts_with("SP")), .funs = funs("sum")) %>%
           as.data.frame())
}

make.nreads.from.merged <- function(d){
    # return the number of reads in the given file.  This is just sum of
    # frequencies.
    my.nreads <- colSums(d) %>% as.data.frame()
    return(my.nreads)
}

make.sgr.from.merged.table <- function(d, my.index, s.name){
    # for each sample generate a SGR file
    sub.df <- get.sub.df(d, my.index, s.name)
    sub.df <- sub.df %>% mutate(window.id = round(position, digits=-2))
    sgr.df <- sub.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
              transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
    return(sgr.df)
}

make.TG.from.merged.table <- function(d, my.index, s.name){
    # from the PyPu files make SGR file for only the TG events.
    sub.df <- get.sub.df(d, my.index, s.name)
    sub.df <- sub.df %>% mutate(window.id = round(position, digits=-2))
    tg.df <- sub.df %>% filter(grepl("TG", dinuc))
    print(head(tg.df)) # for debugging
    tg.sgr <- tg.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
              transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
    return(tg.sgr)
}

make.pyA.from.merged.table <- function(d, my.index, s.name){
    # from the PyPu files make SGR file for only the pyA events.
    sub.df <- get.sub.df(d, my.index, s.name)
    sub.df <- sub.df %>% mutate(window.id = round(position, digits=-2))
    pyA.df <- sub.df %>% filter(grepl("*A", dinuc))
    print(head(pyA.df)) # for debugging
    pyA.sgr <- pyA.df[, .(freq = base::sum(freq)), by = .(chr, window.id, strand)] %>%
               transmute(chr = paste("chr", chr, sep=""), window.id = window.id, freq = paste(strand, freq, sep=""))
    return(pyA.sgr)
}

make.TGsgr.from.list <- function(count.fl, op.dir){
    for(i in 1:length(count.fl)){
        count.fn <- count.fl[i]
        my.suff <- unlist(strsplit(count.fn, "\\."))[3]
	if(my.suff != "PyPu") next
	d <- fread(count.fn, header=T)
	snames <- colnames(d)[6:(ncol(d) - 1)]
	for(j in 1:length(snames)){
            ofn <- file.path(op.dir, paste(snames[j], "-TG.sgr", sep=""))
	    print(ofn) # for debugging
	    write.table(make.TG.from.merged.table(d, j, snames[j]), ofn, quote=F, row.names=F, col.names=F, sep="\t")
	}
    }
}

make.pyAsgr.from.list <- function(count.fl, op.dir){
    for(i in 1:length(count.fl)){
        count.fn <- count.fl[i]
        my.suff <- unlist(strsplit(count.fn, "\\."))[3]
	if(my.suff != "PyPu") next
	d <- fread(count.fn, header=T)
	snames <- colnames(d)[6:(ncol(d) - 1)]
	for(j in 1:length(snames)){
            ofn <- file.path(op.dir, paste(snames[j], "-pyA.sgr", sep=""))
	    print(ofn) # for debugging
	    write.table(make.pyA.from.merged.table(d, j, snames[j]), ofn, quote=F, row.names=F, col.names=F, sep="\t")
	}
    }
}

make.sgr.from.list <- function(count.fl, op.dir){
    for(i in 1:length(count.fl)){
        count.fn <- count.fl[i]
        my.suff <- unlist(strsplit(count.fn, "\\."))[3]
	d <- fread(count.fn, header=T)
	snames <- colnames(d)[6:(ncol(d) - 1)]
	for(j in 1:length(snames)){
            ofn <- file.path(op.dir, paste(snames[j], "-", my.suff, ".all.sgr", sep=""))
	    print(ofn) # for debugging
	    write.table(make.sgr.from.merged.table(d, j, snames[j]), ofn, quote=F, row.names=F, col.names=F, sep="\t")
	}
    }
}

make.nreads.from.list <- function(count.fl){
    nreads.df <- data.frame()
    suffs <- c()
    for(i in 1:length(count.fl)){
        count.fn <- count.fl[i]
	my.suff  <- unlist(strsplit(count.fn, "\\."))[3]
	suffs <- c(suffs, my.suff) 
        d <- fread(count.fn, header=T)
	nreads <- make.nreads.from.merged(d[, 6:(ncol(d)-1)])
	if(i == 1){
            nreads.df <- nreads
	} else {
            nreads.df <- cbind(nreads.df, nreads)
	}
    }
    return(nreads.df)
}

make.dinucs.from.list <- function(count.fl){
    dinucs.df <- data.frame()
    for(i in 1:length(count.fl)){
        count.fn <- count.fl[i]
	d <- fread(count.fn, header=T)
	snames <- colnames(d)[6:(ncol(d)-1)]
        all.dinucs <- make.dinucs.from.merged(d[, 5:(ncol(d) - 1)])
	if(i == 1){
            dinucs.df <- all.dinucs
	} else {
            dinucs.df <- rbind(dinucs.df, all.dinucs)
	}
    }
    dinucs.df$dinuc <- as.character(dinucs.df$dinuc)
    dinucs.df$dinuc[7] <- "NA"
    dinucs.df[is.na(dinucs.df)] <- 0
    return(dinucs.df)
}

make.pwm <- function(df){
    tmp.df <- df %>% dplyr::select(damage.seq) # select on the sequence column
    p1 <- lapply(tmp.df$damage.seq, function(x){return(unlist(strsplit(x, "")))})
    # make a list of vectors, where each vector is actually all the elements of
    # the string.
    p1 <- lapply(p1, function(x){if(length(x) == 21){return(x)}}) # Make the
                                                                  # list of
								  # vector into
								  # a dataframe.
								  # see
								  # https://stackoverflow.com/questions/4227223/r-list-to-data-frame
    lengths <- unlist(lapply(p1, function(x){return(length(x))})) # get lengths
                                                                  # all the
								  # strings
    my.tab <- as.data.frame(table(lengths)) # make the lengths
    # print(my.tab)  # Into a table. So  that we can choose  the number of rows for the next step

    use.len <- my.tab %>% filter(lengths == 21) %>% dplyr::select(Freq) %>% 
               unlist() %>% as.vector() # select the total number of entries
	                                # that have 
    
    p2 <- data.frame(matrix(unlist(p1), nrow = use.len, byrow=T), stringsAsFactors=FALSE)
    for(i in 1:ncol(p2)){
        tmp.df <- as.data.frame(table(p2[, i]))
	colnames(tmp.df) <- c("nuc", paste("pos", i, sep=""))
        if(i == 1){
	    count.df <- tmp.df
	} else {
	    count.df <- merge(count.df, tmp.df, by = "nuc", all = T)}
        }
    count.df <- count.df %>% filter(nuc %in% c("A", "C", "G", "T"))
    rownames(count.df) <- count.df$nuc
    count.df <- count.df[, -1]  # assumes a 20base sequence, first column is the
                                # Nucleotide
    count.df <- count.df[, -11] # assumes a 20base sequence, the 11th column is
                                # an asterix
    count.df[is.na(count.df)] <- 0 # The merge can generate NAs.
    pwm <- count.df / colSums(count.df)
    return(count.df)
}

make.seq.logo <- function(exp.df){
    if(is.character(exp.df)){
        exp.df <- read.delim(exp.df, header=T, stringsAsFactors=F)
    }
    input.fl <- list.files(pattern = "all.no789.*.txt")
    for(i in 1:length(input.fl)){
	df <- fread(input.fl[i], header = T)
        prefix <- unlist(strsplit(input.fl[i], "\\."))[3]
	my.ids <- unique(exp.df$exp.ind)
        for(j in 1:length(my.ids)){
	    print(j) # for debugging
	    use.exp <- exp.df %>% filter(exp.ind == my.ids[j]) %>%
	               select(experiment.name) %>% unlist() %>% as.vector()
	    exposure <- exp.df %>% filter(exp.ind == my.ids[j]) %>%
	               select(UV) %>% unlist() %>% as.vector() %>% unique()
	    cell.type <- exp.df %>% filter(exp.ind == my.ids[j]) %>%
	                 select(cells) %>% unlist() %>% as.vector() %>% unique()
	    print(c(use.exp, exposure, cell.type))
	    ofn <- paste(cell.type, exposure, prefix, ".pdf", sep="") # for debugging
	    print(ofn) # for debugging
	    
	    sub.df <- df %>% select(damage.seq, which(colnames(df) %in% use.exp))
	    print(head(sub.df))

            sub.df <- sub.df %>% mutate(row.sum = rowSums(sub.df[, 2:ncol(sub.df)])) %>%
	              filter(row.sum > 0)
            print(head(sub.df)) # for debugging
	    op.l <- make.pwm(sub.df)
	    pdf(ofn); seqLogo(op.l[[2]]); dev.off()
	}
    }
}

make.pwm.graduated <- function(exp.df, z.df, df){
    # we are going to first select the columns from the exp.df, then we will
    # slice the rows based on the zscore data frame, and then we will generate a
    # pwm
    op.l <- list()
    counter <- 1
    setnames(z.df, c("Chr", "Start", "End"), c("chr", "start", "end"))
    setnames(z.df, colnames(z.df), gsub(" ", "", colnames(z.df)))

    # df <- fread(fn)
    if(is.character(exp.df)){
        exp.df <- read.delim(exp.df, header=T, stringsAsFactors=F)
    }
    my.ids <- unique(exp.df$exp.ind)
    use.p <- list(c(0, 0.05), c(0.05, 0.25), c(0.25, 0.5), c(0.5, 0.75), c(0.75, 0.95), c(0.95, 1))
    
    for(i in 1:length(my.ids)){
        use.exp <- exp.df %>% filter(exp.ind == my.ids[i]) %>%
	           select(experiment.name) %>% unlist() %>% as.vector()
	exposure <- exp.df %>% filter(exp.ind == my.ids[i]) %>%
                    select(UV) %>% unlist() %>% as.vector() %>% unique()
        cell.type <- exp.df %>% filter(exp.ind == my.ids[i]) %>%
                     select(cells) %>% unlist() %>% as.vector() %>% unique()
	sub.df <- df %>% select(chr, position, damage.seq, which(colnames(df) %in% use.exp))
	sub.df <- sub.df %>% mutate(row.sum = rowSums(sub.df[, 4:ncol(sub.df)])) %>%
		  filter(row.sum > 0) %>% setDF()
	sub.df$chr <- paste("chr", sub.df$chr, sep="")
	sub.df$end <- sub.df$position
	sub.df$start <- sub.df$position
        setDT(sub.df)
	setkey(sub.df, chr, start, end)
	print(head(sub.df)) # for debugging

        for(j in 1:length(use.p)){
	    use.range <- use.p[[j]]
	    use.z <- get.subset.z(z.df, use.range[1], use.range[2])
	    setDT(use.z)
	    setkey(use.z, chr, start, end)
	    common <- foverlaps(sub.df, use.z, nomatch = 0L, type = "any")
	    pwm <- make.pwm(common)
	    colnames(pwm) <- 1:ncol(pwm)
	    op.l[[counter]] <- as.matrix(pwm)
	    counter <- counter + 1
	    print(counter)
	}
    }
    return(op.l)
}

## get.subset.z <- function(z.df, v1, v2){
##     if(v1 == 0){
##         use.z <- z.df %>% filter(Z < quantile(Z, v2))
## 	return(use.z)
##      } else if(v2 == 1){
##          use.z <- z.df %>% filter(Z > quantile(Z, v1))
## 	 return(use.z)
##      } else {
##          use.z <- z.df %>% filter(Z > quantile(Z, v1), Z <= quantile(Z, v2))
## 	 return(use.z)
##      }
## }
## 
get.z.lims <- function(z.df, v1, v2){
    if(v1 == 0){
        return(c(min(z.df$Z), quantile(z.df$Z, v2)))
    } else if (v2 == 1){
        return(c(quantile(z.df$Z, v1), max(z.df$Z)))
    } else {
        return(c(quantile(z.df$Z, v1), quantile(z.df$Z, v2)))
    }
}

get.pypy.lims <- function(z.df, v1, v2){
    print(quantile(z.df$PyPys, c(v1, v2))) # for debugging
    if(v1 == 0){
        return(c(min(z.df$PyPys), quantile(z.df$PyPys, v2)))
    } else if (v2 == 1){
        return(c(quantile(z.df$PyPys, v1), max(z.df$PyPys)))
    } else {
        return(c(quantile(z.df$PyPys, v1), quantile(z.df$PyPys, v2)))
    }
}

## get.subset.cp100b <- function(z.df, v1, v2){
##     if(v1 == 0){
##         use.z <- z.df %>% filter(cp100b < quantile(cp100b, v2))
## 	return(use.z)
##     } else if(v2 == 1){
##         use.z <- z.df %>% filter(cp100b > quantile(cp100b, v1))
##         return(use.z)
##     } else {
##         return(z.df %>% filter(cp100b > quantile(cp100b, v1), cp100b <= quantile(cp100b, v2)))
##     }
## }
## 
## get.subset.pypys <- function(z.df, v1, v2){
##     print(quantile(z.df %>% select(PyPys) %>% unlist() %>% as.vector(), c(v1, v2)))
##     if(v1 == 0){
##         use.z <- z.df %>% filter(PyPys < quantile(PyPys, v2))
## 	print(dim(use.z)) # for debugging
## 	return(use.z)
##      } else if(v2 == 1){
##          use.z <- z.df %>% filter(PyPys > quantile(PyPys, v1))
## 	 return(use.z)
##      } else {
##          use.z <- z.df %>% filter(PyPys > quantile(PyPys, v1), PyPys <= quantile(PyPys, v2))
## 	 return(use.z)
##      }
## }
## 
## get.subset.sepdist <- function(z.df, v1, v2){
##     print(quantile(z.df %>% select(sepration.dist) %>% unlist() %>% as.vector(), c(v1, v2)))
##     if(v1 == 0){
##         use.z <- z.df %>% filter(sepration.dist < quantile(sepration.dist, v2))
## 	return(use.z)
##      } else if(v2 == 1){
##          use.z <- z.df %>% filter(sepration.dist > quantile(sepration.dist, v1))
## 	 return(use.z)
##      } else {
##          use.z <- z.df %>% filter(sepration.dist > quantile(sepration.dist, v1), sepration.dist <= quantile(sepration.dist, v2))
## 	 return(use.z)
##      }
## }




## make.pwm.graduated.2 <- function(read.fn, z.fn, sample.name){
##     # We are going to read in the experiment file, and the file with z scores,
##     # and then we are going to slice the experiment file based on the coordiates
##     # that are selected on the z-score and the pypy window values, and create
##     # pwms.
##     op.l <- list()
##     counter <- c()
##     if(is.character(z.fn)){
##         z.df <- fread(z.fn)
##     } else {
##         z.df <- z.fn
##     }
##     print(z.df) # for debugging
##     # setnames(z.df, c("Chr", "Start", "End", "Strand"), 
##     #                c("chr", "start", "end", "strand"))
##     # setnames(z.df, colnames(z.df), gsub(" ", "", colnames(z.df)))
## 
##     if(is.character(read.fn)){
##         df <- fread(read.fn)
##     } else {
##         df <- read.fn
##     }
##     if("row.sum" %in% colnames(df)){
##         df[, row.sum:=NULL] # remove the row.sum column.
##     }
##     df <- df %>% select(chr, position, strand, damage.seq, sample.name)
##     setnames(df, sample.name, "d.count") # temporarily change the column name so
##                                          # that we can process all the columns
## 					 # or a specific column as rquired.
##     df <- df %>% filter(d.count > 0)
##     df$chr <- paste("chr", df$chr, sep="")
##     df$end <- df$position
##     df$start <- df$position
##     setDT(df)
##     setkey(df, chr, strand, start, end)
##     use.p <- list(c(0, 0.05), c(0.05, 0.25), c(0.25, 0.5), c(0.5, 0.75), c(0.75, 0.95), c(0.95, 1))
##     
##     use.names <- c()
## 
##     for(i in seq_along(use.p)){
##         use.range.1 <- use.p[[i]]
##         use.z.1 <- get.subset.pypys(z.df, use.range.1[1], use.range.1[2])
## 	nums <- range(use.z.1$PyPys)
## 	# print(use.z.1) # for debugging
## 	curr.name.1 <- paste("PyPys (", nums[1], ", ", nums[2], ")", sep="")
##         for(j in seq_along(use.p)){
##             use.range.2 <- use.p[[j]]
## 	    use.z <- get.subset.z(use.z.1, use.range.2[1], use.range.2[2])
##             setDT(use.z)
##             setkey(use.z, chr, strand, start, end)
##             common <- foverlaps(df, use.z, nomatch = 0L, type = "any")
## 	    if(nrow(common) == 0) next # we are only interested in the rows that overlap
## 	    nums2 <- range(use.z$Z)
## 	    curr.name.2 <- paste("Z-score (", nums2[1], ", ", nums2[2], ") N = ", nrow(common), sep="")
## 	    curr.name <- paste(curr.name.1, curr.name.2, sep="\n")
## 	    print(dim(common)) # for debugging
##             pwm <- make.pwm(common)
##             colnames(pwm) <- 1:ncol(pwm)
##             # op.l[[counter]] <- as.matrix(pwm)
##             op.l[[curr.name]] <- as.matrix(pwm)
##             counter <- counter + 1
##             print(counter)
## 	    use.names <- c(use.names, curr.name)
## 	}
## 	print(use.names)
##     }
##     names(op.l) <- use.names
##     op.l
## }

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
    my.counts
}

update.big.table <- function(df.fn, d.l, pattern = "PyPy"){
    # we are going to read in one data-frame.  Then read in one file of a type,
    # and merge on the chr, position, strand, and sequence.
    
    df <- fread(df.fn)
    df$dinuc <- gsub("\\*", "", df$dinuc)
    df[, row.sum := NULL, ] # WE WANT TO REMOVE THE LAST COLUMN
    print(dim(df)) # for debugging
    old.col.names <- colnames(df)
    # it might be better to do several small merges, and one big merge.
    # this will need some testing.
    for(i in 1:length(d.l)){
        fn <- file.path(d.l[i], paste(pattern, "-single-recurrent.txt", sep=""))
	print(fn) # for debugging
	count.table <- make.count.table(read.data(fn))
	setDT(count.table)
	setnames(count.table, c("damage.pos", "freq", "damage.match"), c("position", d.l[i], "damage.seq"))
	setDF(count.table)
	print(head(count.table)) # for debugging
	print(head(df)) # for debugging
	df <- merge(df, count.table, by = c("chr", "position", "strand", "damage.seq", "dinuc"), all = T)
	print(dim(df)) # for updating the current status on screen
    }
    new.col.names <- c(old.col.names, d.l)
    df[is.na(df)] <- 0
    print(dim(df)) # for debugging
    setnames(df, colnames(df), new.col.names)
    print(head(df, 3)) # for debugging
    df
}

make.pwm.graduated.3 <- function(read.fn, z.fn, sample.name){
    # We are going to read in the experiment file, and the file with z scores,
    # and then we are going to slice the experiment file based on the coordiates
    # that are selected on the z-score and the pypy window values, and create
    # pwms.
    op.l <- list()
    counter <- 1
    if(is.character(z.fn)){
        z.df <- fread(z.fn)
    } else {
        z.df <- z.fn
    }
    # print(z.df) # for debugging
    # setnames(z.df, c("Chr", "Start", "End", "Strand"), 
    #                c("chr", "start", "end", "strand"))
    # setnames(z.df, colnames(z.df), gsub(" ", "", colnames(z.df)))

    if(is.character(read.fn)){
        df <- fread(read.fn)
    } else {
        df <- read.fn
    }
    if("row.sum" %in% colnames(df)){
        df[, row.sum:=NULL] # remove the row.sum column.
    }
    df <- df %>% select(chr, position, strand, damage.seq, sample.name)
    setnames(df, sample.name, "d.count") # temporarily change the column name so
                                         # that we can process all the columns
					 # or a specific column as rquired.
    df <- df %>% filter(d.count > 0)
    df$chr <- paste("chr", df$chr, sep="")
    df$end <- df$position
    df$start <- df$position
    setDT(df)
    setkey(df, chr, strand, start, end)
    use.p <- list(c(0, 0.05), c(0.05, 0.25), c(0.25, 0.5), c(0.5, 0.75), c(0.75, 0.95), c(0.95, 1))
    
    use.names <- c()

    for(i in seq_along(use.p)){
        for(j in seq_along(use.p)){
	    lim.z <- get.z.lims(z.df, use.p[[i]][1], use.p[[i]][2])
	    lim.pypy <- get.pypy.lims(z.df, use.p[[j]][1], use.p[[j]][2])
	    print(c(lim.pypy, lim.z))
	    use.z <- z.df %>% filter(Z > lim.z[1], Z <= lim.z[2], PyPys > lim.pypy[1], PyPys <= lim.pypy[2])
	    print(dim(use.z))
	    setDT(use.z)
            setkey(use.z, chr, strand, start, end)
            common <- foverlaps(df, use.z, nomatch = 0L, type = "any")
	    if(nrow(common) == 0){
	        pwm <- matrix(rnorm(80, mean = 1, sd = 0.01), ncol = 20, byrow = T)
		rownames(pwm) <- c("A", "C", "G", "T")
	    } else {# we are only interested in the rows that overlap
                pwm <- make.pwm(common)
	    }
            colnames(pwm) <- 1:ncol(pwm)
	    curr.name.1 <- paste("PyPys (", lim.pypy[1], ", ", lim.pypy[2], ")", sep="")
	    curr.name.2 <- paste("Z-score (", lim.z[1], ", ", lim.z[2], "), N = ", nrow(common), sep="")
	    curr.name <- paste(curr.name.1, curr.name.2, sep="\n")
	    print(curr.name)
	    # print(dim(common)) # for debugging
            # op.l[[counter]] <- as.matrix(pwm)
            op.l[[counter]] <- as.matrix(pwm)
            counter <- counter + 1
            # print(counter)
	    use.names <- c(use.names, curr.name)
	}
	print(use.names)
	print(length(op.l))
    }
    names(op.l) <- use.names
    op.l
}
