options(scipen=9999) # we do not want scientific notation

library(data.table) # read and manipulate large files
library(dplyr) # easily summarise data frames.
library(tidyr) # we are going to use the spread function later
library(VennDiagram)
library(RColorBrewer)
library(RCurl) # for the merge.list function
library(stringr) # we are going to do something with string manipulation
library(ggplot2)

#### coding after Sept092016

# First create a window of 3 kb 
##############

tag.df.with.interval <- function(dt, sizes.df,
                                 length = 3000, offset = 100){
    sizes.df <- sizes.df[2:(nrow(sizes.df) - 1), ] # we do not want anything to do
                                             # Mitochondria for now
    dt <- dt %>% mutate(start = damage.pos, end = damage.pos)
    setkey(dt, chr, start, end)

    new.dt.pos  <- setDT(data.frame())
    new.dt.neg  <- setDT(data.frame())
    
    for(i in 1:nrow(sizes.df)){
        
	my.chr <- gsub("chr", "", sizes.df$Chr[i])
	print(my.chr) ## for debugging
        interval.df <- make.interval.df(i, my.chr, sizes.df, length = length, offset = offset)
        setkey(interval.df, chr, start, end)
        
        dt.neg <- dt[strand == "-"]
	dt.pos <- dt[strand == "+"]

        tmp.pos <- foverlaps(interval.df, dt.pos, type = "any", nomatch = 0L)[, .(mean.dist = mean(diff(damage.pos)), sd.dist = sd(diff(damage.pos)), n.events = .N, chr = chr, start = i.start, end = i.end, damage.pos = paste0(damage.pos, collapse = ",")), by = .(id)][n.events > 1][, .(chr = chr, start = max(start), end = min(end), mean.dist = mean(diff(as.integer(unlist(strsplit(damage.pos, ","))))), n.damages = length(unlist(strsplit(damage.pos, ",")))), by = .(damage.pos)] %>% 
          unique() %>% mutate(length = end - start)
        tmp.neg <- foverlaps(interval.df, dt.neg, type = "any", nomatch = 0L)[, .(mean.dist = mean(diff(damage.pos)), sd.dist = sd(diff(damage.pos)), n.events = .N, chr = chr, start = i.start, end = i.end, damage.pos = paste0(damage.pos, collapse = ",")), by = .(id)][n.events > 1][, .(chr = chr, start = max(start), end = min(end), mean.dist = mean(diff(as.integer(unlist(strsplit(damage.pos, ","))))), n.damages = length(unlist(strsplit(damage.pos, ",")))), by = .(damage.pos)] %>% 
	unique() %>% mutate(length = end - start)
	if(i == 1){
            new.dt.pos <- tmp.pos
            new.dt.neg <- tmp.neg
	} else {
            new.dt.pos <- rbind(new.dt.pos, tmp.pos)
            new.dt.neg <- rbind(new.dt.neg, tmp.neg)
	}
    }
    new.dt.pos$strand <- "+"
    new.dt.neg$strand <- "-"
    new.dt <- rbind(new.dt.pos, new.dt.neg)
    return(new.dt)
}

make.interval.df <- function(i, chr, sizes.df, length = 3000, offset = 100){
    starts <- seq(1, sizes.df$Size[i], by = offset)
    interval.df <- setDT(data.frame(chr = chr, start = starts, 
                                    end = (starts + length - 1), 
				    id = paste(chr, starts, (starts + length - 1), sep = ":"))
		        )
    return(interval.df)
}

box.plot.intercpd <- function(fl){
   all.df <- data.frame()
   for(i in 1:length(fl)){
       fn <- fl[i]
       sample <- unlist(strsplit(fl[i], "-"))[1]
       dt <- fread(fn, header=T)
       dt <- dt %>% select(chr, min.dist, strand, n.damages) %>% mutate(sample = sample)
       if(i == 1){
           all.df <- dt
       } else {
           all.df <- rbind(all.df, dt)
       }
   }
   all.df$chr <- factor(all.df$chr, levels = c(as.character(1:22), "X", "Y"))
   p <- ggplot(data = all.df, aes(x = chr, y = min.dist, fill = sample)) +
        geom_boxplot(outlier.size = 2)
}

hist.plot.intercpd <- function(fl){
    all.df <- .make.all.df(fl)
    p <- ggplot(data = all.df, aes(min.dist, ..density.., color = sample)) +
    geom_freqpoly(binwidth=50)
}

.make.all.df <- function(fl){
   all.df <- data.frame()
   for(i in 1:length(fl)){
       fn <- fl[i]
       sample <- unlist(strsplit(fl[i], "-"))[1]
       dt <- fread(fn, header=T)
       dt <- dt %>% select(chr, mean.dist, min.dist, strand, n.damages, unique.pos) %>% mutate(sample = sample)
       if(i == 1){
           all.df <- dt
       } else {
           all.df <- rbind(all.df, dt)
       }
   }
   return(all.df)
}

.get.min.dist <- function(s){
    return(min(diff(as.integer(unlist(strsplit(s, ","))))))
}

.get.unique <- function(s){
    return(length(unique(unlist(strsplit(s, ",")))))
}

add.unique.damages.min.dist <- function(fn){
    dt <- fread(fn, header=T)
    dt[, min.dist := lapply(damage.pos, FUN=.get.min.dist)]
    dt[, unique.pos := lapply(damage.pos, FUN=.get.unique)]
    dt$min.dist <- as.integer(dt$min.dist)
    dt$unique.pos <- as.integer(dt$unique.pos)
    write.table(dt, fn, row.names=F, col.names=T, sep="\t", quote=F)
}

### the following functions assume that we are working in the Ruddle directory
telescoping <- function(fn.pypy.summary, fn.pupy, fn.pypu, fn.pupu){
    # We want to select those windows that have no other damage events other
    # than PyPy on them.
    summary.dt <- fread(fn.pypy.summary, header=T)
    telescope.dt <- summary.dt
    setkey(telescope.dt, chr, strand, start, end)
    print(dim(telescope.dt))
    for(i in c(fn.pupy, fn.pypu, fn.pupu)){
        tmp.dt <- read.data(i) %>% 
	          select(chr, damage.pos, strand) %>%
		  mutate(start = damage.pos, end = damage.pos) %>%
		  select(chr, strand, start, end) %>% unique() %>%
		  setDT()
	print(dim(tmp.dt)) # for debugging
        setkey(tmp.dt, chr, strand, start, end)
	my.overlaps <- foverlaps(telescope.dt, tmp.dt, type = "any", which = T)
	telescope.dt <- telescope.dt[unique(my.overlaps[is.na(yid)]$xid)]
	print(dim(telescope.dt))
	setkey(telescope.dt, chr, strand, start, end)
    }
    return(telescope.dt)
}

make.fns <- function(s, dir1 = "../analysis-sept132016/", suff2 = "-single-recurrent.txt"){
    pypy <- paste(dir1, s, "-PyPy-hotspots.txt", sep="")
    pypu <- paste(s, "/PyPu", suff2, sep="")
    pupy <- paste(s, "/PuPy", suff2, sep="")
    pupu <- paste(s, "/PuPu", suff2, sep="")
    return(c(pypy, pypu, pupy, pupu))
}

save.pure.pypy <- function(s){
    fns <- make.fns(s)
    ofn <- paste("../analysis-sept282016/", s, "-selected-pure-pypy.txt", sep="")
    telescope.dt <- telescoping(fns[1], fns[2], fns[3], fns[4])
    write.table(telescope.dt, ofn, row.names=F, col.names=T, sep="\t", quote=F)
}

### Functions after oct142016
get.common.windows.overlaps <- function(fl){
    df1 <- fread(fl[1], header=T)
    common.windows <- df1
    for(i in 2:length(fl)){
        s.name <- unlist(strsplit(fl[i], "-"))[1]
        df2 <- fread(fl[i], header=T)
	setkey(df2, chr, strand, start, end)
        common.windows <- foverlaps(common.windows, df2, type = "any", nomatch = 0L)
	new.col.names <- gsub("^i.", paste(s.name, ".", sep=""), colnames(common.windows))
	setnames(common.windows, colnames(common.windows), new.col.names)
	print(dim(common.windows)) # for debugging
    }
    return(common.windows)
}

### Functions for final hotspot detection oct252016
get.merged.files <- function(f){
    my.df <- fread(f, header=T)
    setkey(my.df, chr)

    tmp.2 <- data.frame()
    for(i in unique(my.df$chr)){
        tmp.chr <- my.df[i] %>% select(chr, start, end, strand) %>% arrange(start)
	print(c(i, dim(tmp.chr))) # for debugging
	tmp.2 <- rbind(tmp.2, tmp.chr) # to enable debugging above, the original
                                       # method of doing everything in single
				       # line was most efficent code.
    }
    setDT(tmp.2)

    ofn.pre <- unlist(strsplit(f, "-"))[1]
    write.table(tmp.2[strand == "+"], "tmp.bed", row.names = F, col.names = F, 
                quote= F, sep="\t")
    ofn.pos <- paste(ofn.pre, "-pos-merged.bed", sep="")
    use.mergeBed(ofn.pos)
    write.table(tmp.2[strand == "-"], "tmp.bed", row.names = F, col.names = F, 
                quote= F, sep="\t")
    ofn.neg <- paste(ofn.pre, "-neg-merged.bed", sep="")
    use.mergeBed(ofn.neg)
}

use.mergeBed <- function(fn){
    command <- paste("mergeBed -n -i tmp.bed > ", fn, sep="")
    print(command) # for debugging
    system(command)
}

get.length.distribution <- function(fl){
    my.df <- data.frame()
    my.df.2 <- data.frame()
    my.col.names <- c("chr", "start", "end", "n")
    for(f in fl){
        
	fn.parts <- unlist(strsplit(f, "-"))
	sample.name <- fn.parts[1]
	strand      <- fn.parts[2]

	tmp.df <- fread(f, header=F) # read in the file
	setnames(tmp.df, colnames(tmp.df), my.col.names) # set correct names
	tmp.df[, n := NULL] # delete the last column

	tmp.df <- tmp.df %>% mutate(length = end - start, strand = strand,
	                            sample = sample.name)

	setDT(tmp.df) # we need to make the data frame into data.table
        my.df <- rbind(my.df, 
	               tmp.df[, list(min.window.length = min(length),
	                             mean.window.length = mean(length), 
	                             max.window.length = max(length), 
				     median.window.length = quantile(length, 0.5), 
			             number.of.windows = .N), 
			         # by = list(sample, strand)]
			         by = .(sample)])
	my.df.2 <- rbind(my.df.2, tmp.df)
	                     
    }
    return(list(my.df, my.df.2))
}

### Functions after oct262016
get.sample.names <- function(fl){
    return(unique(unlist(lapply(fl, FUN=function(x){
                                         return(unlist(strsplit(x, "-"))[1])
					}
			       )
		        )
	         )
          )
}

get.df.from.sample <- function(s.name){
    print(s.name)
    my.col.names <- c("chr", "start", "end", "n", "strand")
    tmp.df1 <- fread(paste(s.name, "-pos-merged.bed", sep=""), header=F)
    tmp.df1$strand <- "+"
    tmp.df2 <- fread(paste(s.name, "-neg-merged.bed", sep=""), header=F)
    tmp.df2$strand <- "-"
    my.df <- rbind(tmp.df1, tmp.df2)
    setnames(my.df, colnames(my.df), my.col.names)
    setDT(my.df)
    print(dim(my.df)) # for debugging
    return(my.df)
}

get.overlaps.from.merged <- function(s.names){
    df1 <- get.df.from.sample(s.names[1])
    for(i in 2:length(s.names)){
        s.name <- s.names[i]
        df2 <- get.df.from.sample(s.name)
	setkey(df2, chr, strand, start, end)
	df1 <- foverlaps(df1, df2, type = "any", nomatch = 0L)
	new.col.names <- gsub("^i.", paste(s.name, ".", sep=""), colnames(df1))
	setnames(df1, colnames(df1), new.col.names)
	print(dim(df1)) # for debugging
    }
    return(df1)
}

make.barplot.merged <- function(summary.df, source = "dna", exposure = "exposed", strand = F, ...){
    
    setDT(summary.df)
    setkey(summary.df, source)
    my.names <- summary.df[source]
    setkey(my.names, exposure)
    my.names <- my.names[exposure] %>% select(sample) %>% unlist() %>%
                as.vector()
    print(my.names) # for debugging
    
    merged.df <- get.overlaps.from.merged(my.names)
    print(head(merged.df)) # for debugging

    expected  <- summary.df %>% 
                 filter(source == source, exposure == exposure) %>% 
		 select(total.windows, expected.probability.median, expected.probability.mean) %>% 
		 transmute(n.cells.mean = prod(expected.probability.mean) * sum(total.windows), 
		           n.cells.median = prod(expected.probability.median) * sum(total.windows)) %>% 
	          unique()
    print(expected) # for debugging
    
    ### Start plotting here
    if(strand){
        df.to.plot <- merged.df[, .(n.windows = .N), by = .(chr, strand)] %>%
	              tidyr::spread(strand, n.windows) %>%
		      arrange(as.numeric(chr))
    } else {
        df.to.plot <- merged.df[, .(n.windows = .N), by = .(chr)] %>%
	              arrange(as.numeric(chr))
    }
    
    rownames(df.to.plot) <- df.to.plot$chr
    mat.to.plot <- df.to.plot %>% select(n.windows)
    mat.to.plot <- as.matrix(mat.to.plot)
    print(mat.to.plot)
    

    if(strand){
        barplot(t(mat.to.plot), col = c("blue", "red"), beside = T, ...)
	legend("topright", c("-", "+"), main = "Strand", fill = c("blue", "red"),
	       bty = "n")
    } else {
        barplot(mat.to.plot, names.arg = rownames(mat.to.plot), 
	        xlab = "Chromosomes", ...)
    }

   abline(h = expected[, 1], col = "blue", lty = 2, lwd = 2)
   abline(h = expected[, 2], col = "red", lty = 2, lwd = 2)

    # text(nrow(mat.to.plot), y = expected[, c(1, 2)], col = c("red", "darkred"))
}
