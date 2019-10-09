library(data.table) # read and manipulate large files
library(dplyr) # easily summarise data frames.
library(tidyr) # we are going to use the spread function later
library(VennDiagram)
library(RColorBrewer)
library(RCurl) # for the merge.list function
library(stringr) # we are going to do something with string manipulation
library(ggplot2)
library(pheatmap)

op <- options()
options(scipen = 7)

read.data <- function(fn){
    DT <- fread(fn, stringsAsFactors=F, header=F)
    my.col.names <- c("chr", "damage.pos", "strand", "damage.match", "length",
                      "flags1", "flags2", "CIGAR1", "flags3", "damage.end", "damage.emd.match",
                      "flags4", "flags5", "CIGAR2", "flags6", "flags7", "flags8")
    setnames(DT, colnames(DT), my.col.names)
    return(DT)
}

get.list <- function(pat){
    command = paste("ls", pat)
    print(command) # for debugging
    fl <- system(command, intern=T)
    print(fl) # for debugging
    return(fl)
}

make.vector.for.venn <- function(fname, window, mito=F){
    # Get a filename, read in the data, and convert it to a vector that can be
    # used for a venn diagrama
    dt <- read.data(fname)
    if(mito){
        dt <- dt[chr == "MT"]
    } else {
        dt <- dt[chr != "MT"]
    }
    tmp.dt <- aggregate.on.window(dt, mito = mito, window = window)
    tmp.vec <- make.unique.vector(tmp.dt)
}

make.four.way.venn <- function(pat, window = 100, my.palette = "Dark2", mito=F, ...){
    # make a 4 way venn diagram for each sample
    # This version has all the files read into the memory first and then
    # processed.  WE need to make this more efficient when dealing with bith
    # files
    fl <- get.list(pat)
    
    v1 <- make.vector.for.venn(fl[1], mito=mito, window=window)
    v2 <- make.vector.for.venn(fl[2], mito=mito, window=window)
    v3 <- make.vector.for.venn(fl[3], mito=mito, window=window)
    v4 <- make.vector.for.venn(fl[4], mito=mito, window=window)
    
    total.n <- length(union(v1, union(v2, union(v3, v4))))

    my.cols <- brewer.pal(4, my.palette)
    vd <- venn.diagram(list(v1, v2, v3, v4), 
                       category.names = c("PuPu", "PuPy", "PyPu", "PyPy"),
		       fill = my.cols, col=NA,
		       alpha = c(0.45, 0.45, 0.45, 0.45),
		       filename = NULL, sub=paste("Total N = ", total.n, sep=""),
		       ...)
    return(vd)
    
}

aggregate.on.window <- function(dt, window, strand = FALSE, mito = FALSE){
    window.val <- -log10(window) # calculate the second argument for "round"
    
    if(mito){
        dt <- dt[chr == "MT"]
    } else {
        dt <- dt[chr != "MT"]
    }
    if(strand){
        dt.aggregate <- dt[, .(count = .N), .(Chr = chr, Window =round(damage.pos, digits = window.val), strand = strand)] 
	# aggregate data on window, but keep information on strand on a
	# different row.
    } else {
       dt.aggregate <- dt[, .(count = .N), .(Chr = chr, Window = round(damage.pos, digits = window.val))] #
                                                                        # aggregate the data on the window.
    }
    # rounds it "down", so 459 will become 400, so the window is actually 400 -
    # 500
    return(dt.aggregate)
}

make.unique.vector <- function(dt){
    
    my.vec <- paste(dt$Chr, dt$Window, sep="-")
    return(my.vec)
}

make.ready.for.circos <- function(dt, ofn, window, strand=FALSE){
   # this function will generate a table ready to used in CIRCOS, and write it
   # to the fn.  We are using the "window" only to calculate the bin, we
   # actually do not need the window value itself for anything else
   if(nrow(dt) < 25000){
       if(strand){
           write.table(dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count, strand = strand), ], 
                      ofn, quote=F, row.names=F, col.names=F, sep=" ")
           print(paste("File ", ofn, " written"))
        } else {
           write.table(dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count), ], 
                      ofn, quote=F, row.names=F, col.names=F, sep=" ")
           print(paste("File ", ofn, " written"))
	}
    }
    else { 
        print("Table has too many rows, that cannot be plotted with CURRENT CIRCOS SETTINGS (25000 Max per track)")
    }
}

make.circos.ready.forRMS <- function(fn, ofn, window = 100, mito = FALSE, strand = TRUE, cutoff = 0.9){
    dt <- read.data(fn)
    dt <- aggregate.on.window(dt, window = window, mito = mito, strand = strand)
    dt <- dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count, strand = strand), ]
    dt <- get.correct.number.of.entries(dt)
    write.table(dt, ofn, sep=" ", row.names=F, col.names=F, quote=F)
}

make.circos.ready.forRMS.nocentro <- function(fn, ofn, ofn2, ofn3, window=100, mito=FALSE, strand = TRUE, cutoff = 0.9){
    dt <- read.data(fn)
    centro.df <- read.delim("/home/sm2556/GitHub/CPDAnalysis/Code/resources/centromere.info.txt", as.is=T, header=T, stringsAsFactors=FALSE)
    op.l <- make.no.centro.dt(dt, centro.df)
    dt <- op.l$no.centro
    dt <- aggregate.on.window(dt, window = window, mito = FALSE, strand = strand)
    dt <- dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count, strand = strand), ]
    dt <- get.correct.number.of.entries(dt)
    write.table(dt, ofn, sep=" ", row.names=F, col.names=F, quote=F)
    write.table(op.l$centromeric.df, ofn2, sep="\t", row.names=F, col.names=T, quote=F)
    write.table(op.l$stats.df, ofn3, sep="\t", row.names=F, col.names=T, quote=F)
}

make.no.centro.dt <- function(dt, centro){
    no.centro <- data.frame()
    stats.df  <- data.frame()
    for(i in 1:nrow(centro)){
        my.cen <- centro[i, ] # iterate over each row of the centromere
	my.chr <- gsub("chr", "", my.cen$chr) # our data frame does not contain
	                                      # "chr" in the chromosome
	tmp <- as.data.table(dt) %>% filter(chr == my.chr,
	                     damage.pos < my.cen$bord.start | damage.pos > my.cen$bord.end)
	
	centromeric.count <- as.data.table(dt) %>% filter(chr == my.chr, # note the change in direction of comparison in the next line
	                     damage.pos > my.cen$bord.start,  damage.pos < my.cen$bord.end)
	total.chr         <- nrow(as.data.table(dt) %>% filter(chr == my.chr))
	no.centro <- rbind(no.centro, tmp)
	st.op <- data.frame(chr = my.chr, 
	                    centromeric = nrow(centromeric.count),
			    noncentromeric = nrow(tmp),
			    perc.centromeric = (nrow(centromeric.count)/total.chr)*100)
	stats.df <- rbind(stats.df, st.op)
    }
    return(list(stats.df = stats.df, no.centro = no.centro, centromeric.df = centromeric.count))
}

make.circos.unique.mapping <- function(fn, ofn, window = 100, mito= FALSE, strand=TRUE){
    dt <- read.data(fn)
    dt <- dt %>% filter(flags2 > 49, flags5 > 49)
    dt <- aggregate.on.window(dt, window = window, mito = mito, strand = strand)
    dt <- dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count, strand = strand), ]
    dt <- get.correct.number.of.entries(dt)
    write.table(dt, ofn, sep=" ", row.names=F, col.names=F, quote=F)
}

.make.graph.input.dt <- function(dt, window){
    dt <- dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), mp = Window, count = count, strand = strand), ]
    return(dt)
}

make.graph.input <- function(fn, ofn, window = 100, mito = FALSE, strand = TRUE){
    # This function will generate information in a format that will allow us to
    # plot this in many ways
    dt <- read.data(fn)
    dt.geno <- aggregate.on.window(dt, window = window, mito = mito, strand = strand)
    dt.mito <- aggregate.on.window(dt, window = 10, mito = TRUE, strand = strand)
    dt.geno <- .make.graph.input.dt(dt.geno, 100)   
    dt.mito <- .make.graph.input.dt(dt.mito, 10) 
    write.table(dt.geno, ofn, sep=" ", row.names=F, col.names=T, quote=F)
    write.table(dt.mito, paste(ofn, "mito", sep=""), sep=" ", row.names=F, col.names=T, quote=F)
}

make.graph.input.RMS <- function(pat, suffix = "*-single-recurrent.txt",
                                 window = 100, mito = FALSE, strand = TRUE){

     pat <- paste(pat, suffix, sep="/")
     command <- paste("ls", pat, sep=" ")
     fl <- system(command, intern=T)
     print(fl) # for debugging
     for(i in 1:fl){
         ofn <- paste(unlist(strsplit(fl[i], "."))[1], "forgraph", sep=".")
	 make.graph.input(fl[i], ofn)
     }
}

get.correct.number.of.entries <- function(dt){
    if(nrow(dt) < 25000){
        return(dt)
    }
    for(i in 1:5){
        cutoff <- ((10 ** i) - 1) / (10 ** i)
	print(cutoff)
	print(quantile(dt$count, cutoff))
        new.dt <- dt %>% filter(count >= quantile(count, cutoff))
	print(nrow(new.dt))
	if(nrow(new.dt) < 25000){
            return(new.dt)
        }
    }
}

gather.data.for.sample <- function(pat, window = 10000, mito=FALSE, strand = FALSE){
    # return a data-frame for a pattern and window
    use.pat <- paste("ls ", pat, "*.txt", sep="")
    fl <- get.list(use.pat)
    
    # read in each file
    pupu <- read.data(fl[1]) 
    pupy <- read.data(fl[2])
    pypu <- read.data(fl[3])
    pypy <- read.data(fl[4])
    
    if(mito){
        if(strand){
            pupu.aggregate <- aggregate.on.window(pupu, window, mito=mito, strand=stand)
            pupy.aggregate <- aggregate.on.window(pupy, window, mito=mito, strand=stand)
            pypu.aggregate <- aggregate.on.window(pypu, window, mito=mito, strand=stand)
            pypy.aggregate <- aggregate.on.window(pypy, window, mito=mito, strand=stand)
	} else {
            pupu.aggregate <- aggregate.on.window(pupu, window, mito=mito)
            pupy.aggregate <- aggregate.on.window(pupy, window, mito=mito)
            pypu.aggregate <- aggregate.on.window(pypu, window, mito=mito)
            pypy.aggregate <- aggregate.on.window(pypy, window, mito=mito)

	}
	description = paste(pat, window, "bp, (Mitochondrial genome only)")
    } else {
        if(strand){
            pupu.aggregate <- aggregate.on.window(pupu, window, strand=T)
            pupy.aggregate <- aggregate.on.window(pupy, window, strand=T)
            pypu.aggregate <- aggregate.on.window(pypu, window, strand=T)
            pypy.aggregate <- aggregate.on.window(pypy, window, strand=T)

	} else{
            pupu.aggregate <- aggregate.on.window(pupu, window)
            pupy.aggregate <- aggregate.on.window(pupy, window)
            pypu.aggregate <- aggregate.on.window(pypu, window)
            pypy.aggregate <- aggregate.on.window(pypy, window)

	}
	description = paste(pat, window, "bp, (Genomic DNA only)")
    }
	
    

    op.l <- list(pupu = table(pupu.aggregate$count), 
                 pupy = table(pupy.aggregate$count),
		 pypu = table(pypu.aggregate$count),
		 pypy = table(pypy.aggregate$count))
    
    # now plot these numbers!
    
    quartz()
    layout(matrix(rbind(c(1, 1, 1),
                        c(2, 3, 4)), nrow = 2),
	   heights = c(1, 1), respect=T)
    my.max <- max(c(op.l$pypy, op.l$pupu, op.l$pypu, op.l$pupy))
    my.ylims <- c(1, my.max)
    my.opts <- c(cex.lab = 0.7, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)

    ##TODO!! see if we can actually print the numbers on the top of the bars, so
    ## so that the bars are more readable.  This information will be similar to
    ## 4-way Venn-diagram
    barplot(op.l$pypy, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PyPy", sub = description, cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pupu, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PuPu", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pupy, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PuPy", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pypu, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PyPu", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    return(op.l)
}

make.partner.dist <- function(dt, mito=FALSE){
    if(mito){
        dt.tmp <- dt[chr == "MT"][, .(chr = chr, damage.pos = damage.pos, strand = strand, window.10b = round(damage.pos, -1), window.100b = round(damage.pos, -2), damage.dimer = unlist(lapply(damage.match, FUN=function(x){ gsub("\\*", "", paste0(unlist(strsplit(as.character(x), ""))[10:12], collapse = ""))}))), ]
    } else {
        dt.tmp <- dt[, .(chr = chr, damage.pos = damage.pos, strand = strand, window.100b = round(damage.pos, -2), window.1kb = round(damage.pos, -3), window.10kb = round(damage.pos, -4), damage.dimer = unlist(lapply(damage.match, FUN=function(x){ gsub("\\*", "", paste0(unlist(strsplit(as.character(x), ""))[10:12], collapse = ""))}))), ]
    }
    return(dt.tmp)
}

make.dimer.count.dist <- function(dt, mito=FALSE){
    if(mito){
        dt.tmp <- dt[chr == "MT"][, .(damage.dimer = unlist(lapply(damage.match , FUN=function(x){gsub("\\*", "", paste0(unlist(strsplit(as.character(x), ""))[10:12], collapse = ""))}))), ]
    } else {
        dt.tmp <- dt[chr != "MT"][, .(damage.dimer = unlist(lapply(damage.match , FUN=function(x){gsub("\\*", "", paste0(unlist(strsplit(as.character(x), ""))[10:12], collapse = ""))}))), ]
    }
    return(dt.tmp)
}

make.marginal.dimer.dist <- function(dt, mito=FALSE){
    if(mito){
        dt.tmp <- dt[chr == "MT"][, .(damage.dimer = unlist(lapply(damage.match, .damage.dimer)), marginal = unlist(lapply(damage.match, .context.dimer))), ][, .(freq = .N), by = .(marginal, damage.dimer)]
    } else {
        dt.tmp <- dt[chr != "MT"][, .(damage.dimer = unlist(lapply(damage.match, .damage.dimer)), marginal = unlist(lapply(damage.match, .context.dimer))), ][, .(freq = .N), by = .(marginal, damage.dimer)]
    }
    return(dt.tmp) # return the data-table, with all the necessary columns.
}

## trying to see if external function can be accepted by data.table
.damage.dimer <- function(x){
    return(gsub("\\*", "", paste0(unlist(strsplit(as.character(x), ""))[10:12], collapse="")))
}

.context.dimer <- function(x){
 return(paste0(unlist(strsplit(gsub("\\*", "", paste0(unlist(strsplit(x, ""))[9:13], collapse="")), ""))[c(1, 4)], collapse=", "))
}

get.marginal.dt <- function(fn, mito=FALSE){
    tmp <- make.marginal.dimer.dist(read.data(fn), mito=mito)
    return(tmp)
}

make.marginal.count.df <- function(pat, recurrent=FALSE, mito=FALSE, normalize=FALSE){
    if(recurrent){
        command <- paste(pat, "/*_recurrent.txt", sep="")
    } else {
        command <- paste(pat, "/*_single.txt", sep="")
    }
    command <- paste("ls ", command, sep="")
    fl <- system(command, intern=T)
    if(normalize){
        norm.factor <- get.normalization.factor(pat)
    }
    print(fl) # for debugging
    for(i in 1:length(fl)){
        print(fl[i]) # for debugging
        tmp <- get.marginal.dt(fl[i], mito=mito) %>% 
	        filter(!(grepl("N", damage.dimer))) %>%
		filter(!(grepl("N", marginal)))
	tmp2 <- tmp %>% spread(damage.dimer, freq, fill=0)
	if(i == 1){
            op.df <- tmp2
	} else {
            op.df <- merge(op.df, tmp2, by = "marginal")
	}
    }
    if(normalize){
        op.df <- normalize.df(op.df, norm.factor)
    }
    print(head(op.df)) # for debugging
    return(op.df)
}

normalize.df <- function(df, norm.factor){
    print(norm.factor) # for debugging
    print(dim(df)) # for debugging
    print(head(df)) # for debugging
    tmp <- as.data.frame(df)
    print(head(tmp))
    for(i in 1:ncol(tmp)){
        print(paste(i, class(tmp[i]))) # for debugging
    }
    tmp[2:17] <- lapply(tmp[2:17], FUN=function(x){as.vector(x) / norm.factor})
    print(head(tmp)) # for debugging
    return(tmp)
}

make.dimer.count.vector <- function(pat, suffix = "/*_single.txt", mito=FALSE, normalize=T){
    command <- paste(pat, suffix, sep="") # get sample name, and obtain
                                                  # list of all the files
    print(command) # for debugging
    command <- paste("ls ", command, sep="")
    print(command) # for debugging
    fl <- system(command, intern=T)
    print(fl) # for debugging
    damage.vec <- c()
    dimer.names <- c()
    for(i in 1:length(fl)){
        print(fl[i]) # for debugging
        dt <- read.data(fl[i])
	dt.tmp <- make.dimer.count.dist(dt, mito=mito)
	if(nrow(dt.tmp) == 0){
            dimer.counts <- as.table(c(0, 0, 0, 0))
	    if(i == 1){
                names(dimer.counts) <- c("AA", "AG", "GA", "GG")
	    } else if(i == 2) {
                names(dimer.counts) <- c("AC", "AT", "GC", "GT")
	    } else if(i == 3) {
                names(dimer.counts) <- c("CA", "CG", "TA", "TG")
	    } else {
                names(dimer.counts) <- c("CC", "CT", "TC", "TT")
	    }
	    print(dimer.counts) # for debugging
	} else {
	    dt.tmp <- filter(dt.tmp, !(grepl("N", damage.dimer))) # there were some
	                                                      # cases where
							      # there were "N"
							      # in the damage
							      # dimer. we are
							      # accounting of
							      # that here.
	    dimer.counts <- table(dt.tmp$damage.dimer)
	}
	dimer.names <- c(dimer.names, names(dimer.counts))
        damage.vec <- c(damage.vec, as.vector(dimer.counts))
    }
    names(damage.vec) <- dimer.names
    print(damage.vec) # for debugging
    if(normalize){
        norm.factor <- get.normalization.factor(pat)
	damage.vec <- damage.vec / norm.factor
    }
    print(damage.vec) # for debugging
    if(length(damage.vec) != 16){
        print("Something is not right here")
    }
    return(damage.vec)
}

make.dimer.count.matrix <- function(my.suff = c("04", "05", "06"), normalize=TRUE, mito=FALSE){

    mat <- matrix(ncol = 16) # creates a first row of all NAs
    my.r.names <- c()
    norm.factors <- c()
    for(suff in my.suff){
        my.pat <- paste("SP", suff, "HU", sep="")
	if(normalize){
            norm.fact <- get.normalization.factor(my.pat)
	    norm.factors <- c(norm.factors, norm.fact)
	    print(norm.factors) # for debugging
	dimer.vec <- make.dimer.count.vector(my.pat, mito=mito)
	mat <- rbind(mat, dimer.vec)
	# print(mat) # for debugging # this is working
	my.r.names <- c(my.r.names, my.pat)
	}
    }
    colnames(mat) <- names(dimer.vec)
    mat <- mat[-1, ] # we do not need the first row, its all NAs
    rownames(mat) <- my.r.names
    if(normalize){
        print(mat) # for debugging
        mat <- sweep(mat, 1, norm.factors, "/")
        print(mat) # for debugging
    }
    return(mat)
}

##TODO!! check if there are 16 dimers in the damage events, if not figure out
# correct number of colours for each of the 
check.pupu <- function(dinucs){
    pupu <- c("AA", "AG", "GA", "GG")
    print(paste("PuPu:", sum(pupu %in% dinucs))) # for debugging
    return(sum(pupu %in% dinucs))
}
check.pupy <- function(dinucs){
    pupy <- c("AC", "AT", "GC", "GT")
    print(paste("PuPy:", sum(pupy %in% dinucs))) # for debugging
    return(sum(pupy %in% dinucs))
}
check.pypu <- function(dinucs){
    pypu <- c("CA", "CG", "TA", "TG")
    print(paste("PyPu:", sum(pypu %in% dinucs))) # for debugging
    return(sum(pypu %in% dinucs))
}
check.pypy <- function(dinucs){
    pypy <- c("CC", "CT", "TC", "TT")
    print(paste("PyPy:", sum(pypy %in% dinucs))) # for debugging
    return(sum(pypy %in% dinucs))
}

make.colors <- function(dinucs){
    all.pupu <- check.pupu(dinucs)
}

plot.dimer.matrix.2 <- function(dimer.count.mat, ofn, normalized=TRUE, ...){
    my.cols.1 <- c(brewer.pal(6, "Blues")[3:6],
                   brewer.pal(6, "Greens")[3:6],
                   brewer.pal(6, "Purples")[3:6],
		   brewer.pal(6, "Reds")[3:6]) # reds for pypy because that is
		                               # our primary interest
    
    if(normalized){
        my.lab <- "Frequency (Per million reads)"
    } else {
        my.lab <- "Frequency"
    }
    options("scipen" = 0) # this option is necessary here
    pdf(ofn, width = 10, height = 3.5)
    par(mar=c(3, 5, 1, 1))
    mps <- barplot(t(dimer.count.mat), beside=T, col = my.cols.1,# ylab = my.lab,
            xlab = "", cex.lab=0.9, cex.axis = 0.9, yaxt="n", ...)
    mps <- colMeans(mps)
    print(mps) # for debugging
    axis(2, at = y <- seq(0, 5e+05, by = 100000), las=2, 
    labels = c("0", 
               expression(paste(1, "x", 10^5)), 
	       expression(paste(2, "x", 10^5)),
	       expression(paste(3, "x", 10^5)), 
	       expression(paste(4, "x", 10^5)),
	       expression(paste(5, "x", 10^5))))
    legend("topright", colnames(dimer.count.mat), fill = my.cols.1, ncol=4,
           bty = "n", cex = 1.3, x.intersp=0.2)
    mtext(text = my.lab, side = 2, line = 4, cex = 1.1) # handle the ylab here
                                                      # separately

    dev.off()
}

plot.dimer.matrix <- function(dimer.count.mat, ofn, normalized=T){
    my.cols.1 <- c(brewer.pal(6, "Blues")[3:6],
                   brewer.pal(6, "Greens")[3:6],
                   brewer.pal(6, "Purples")[3:6],
		   brewer.pal(6, "Reds")[3:6]) # reds for pypy because that is
		                               # our primary interest
    
    if(nrow(dimer.count.mat) == 3){
        my.cols.2 <- brewer.pal(6, "Reds")[3:5]
    } else if(nrow(dimer.count.mat) == 2){
        my.cols.2 <- brewer.pal(6, "Reds")[3:4]
    } else {
        my.cols.2 <- brewer.pal(nrow(dimer.count.mat), "Reds")
    }
    if(normalized){
        my.lab <- "Frequency (Per million reads)"
    } else {
        my.lab <- "Frequency"
    }
    if(ncol(dimer.count.mat) < 16){
        dimer.count.mat <- fix.columns.dimer.matrix(dimer.count.mat)
    }
    pdf(ofn)
    split.screen(c(2, 1))
    screen(1)
    barplot(dimer.count.mat, beside=T, col = my.cols.2, ylab = my.lab, 
            xlab = "Dinucleotide", cex.axis=0.8, cex.lab=0.8, cex.names=0.6)
    legend("topleft", rownames(dimer.count.mat), fill = my.cols.2, bty="n",
           cex=0.6, ncol=nrow(dimer.count.mat/3))
    close.screen(1)
    screen(2)
    barplot(t(dimer.count.mat), beside=T, col = my.cols.1, ylab = my.lab,
            xlab = "Experiment", cex.axis=0.6, cex.lab=0.6)
    legend("topright", colnames(dimer.count.mat), fill = my.cols.1, ncol=4, bty="n", cex = 0.5)
    close.screen(2)
    dev.off()
}

fix.columns.dimer.matrix <- function(mat, my.col.names){
    my.col.names <- c("AA", "AG", "GA", "GG", 
                      "AC", "AT", "GC", "GT", 
		      "CA", "CG", "TA", "TG", 
		      "CC", "CT", "TC", "TT")
    ref.mat <- data.frame(dinucs = my.col.names)
    tmp.mat <- as.data.frame(t(mat))
    tmp.mat$dinucs <- rownames(tmp.mat)
    use.mat <- merge(ref.mat, tmp.mat, by="dinucs", all=T)
    print(use.mat) # for debugging
    rownames(use.mat) <- use.mat$dinucs
    use.mat <- use.mat[, -1]
    use.mat <- as.matrix(use.mat)
    print(use.mat) # for debugging
    return(t(use.mat))
}

get.normalization.factor <- function(pat){
##             V1     V2     V3
## 	    1       Linker  Match 896372
## 	    2     Unmapped  56447     NA
## 	    3     Multimap  46662     NA
## 	    4  Wronglength  89469     NA
## 	    5        Other  25459     NA
## 	    6   Deanimated    279     NA
    fn <- paste(pat, "/countlog.txt", sep="")
    df <- read.delim(fn, skip=1, sep="", stringsAsFactors=F, header=F)
    correction.number <- as.numeric(df$V3[1]) - (sum(as.numeric(df$V2[2:6])))    
       
    correction.factor <- correction.number / 1000000 # so many million reads. If
                                                     # total is less than a
						     # million, we will end up
						     # "increasing the number
    ##TODO!! WE NEED TO ADD A FUNCTIONALITY HERE THAT WILL DECIDE ON THE
    # NORMALIZATION factor based on the number itself. This should technically be the
    # same for all n conditions rather than three conditions
    return(correction.factor) 
}

# RMS related function there

make.matrix.from.vec <- function(pat, fn, suff = "-single"){ 
    # fn is the sample-sheet
    # the pat has to be either genomic or mito
    use.samples <- read.delim(fn, as.is=T, header=T)[, 1]
    print(use.samples) # for debugging
    fl <- paste(use.samples, "-dinuc-", pat, suff, ".vec", sep="") # not the
                                                                   # elegantest of solution
    print(fl) # for debugging
    rnames <- c()
    for(i in 1:length(fl)){
        vname <- unlist(strsplit(fl[i], "-"))[1]
        rnames <- c(rnames, vname) 
        tmp.v <- get(load(fl[i]))
        if(i == 1){
            tmp.mat <- matrix(tmp.v, ncol=length(tmp.v))
	    colnames(tmp.mat) <- names(tmp.v)
        } else {
            my.vec <- get(load(fl[i]))
	    tmp.mat <- rbind(tmp.mat, my.vec)
        }
    }
    rownames(tmp.mat) <- rnames
    ofn <- paste(paste(pat, "-dinuc-matrix", suff, sep=""), ".mat", sep="")
    print(ofn) # for debugging
    write.table(tmp.mat, ofn, sep="\t", row.names=T, col.names=T,
			 quote=F)
    return(tmp.mat)
}

read.marginal <- function(pat, recurrent=FALSE, normalized=FALSE){
    # the pat has to be either genomic or mito
    samples <- as.character(read.delim("sample-sheet.txt", header=T, comment.char="#")$sample)
    print(samples) # for debugging
    fl <- paste(samples, "-recurrent-marginal-normalized-", pat, ".txt", sep="")
    ##TODO!! HANDLE ALL THE SUFFIXES.  WE WILL MAKE SURE THAT THE OUTPUT FILES
    # A PATTERNED NAME SO THAT WE CAN WRITE A FUNCTION TO DETERMINE THE PATTERN
    # FOR ALL THE FILES
    for(i in 1:length(samples)){
        vname <- samples[i]
	tmp <- read.delim(fl[i], header=T, stringsAsFactors=F)
	print(tmp) # for debugging
	print(dim(tmp)) # for debugging
        tmp$data <- samples[i]
        if(i == 1){
	    op.df <- tmp
	    print(op.df) # for debugging
	} else {
	    op.df <- rbind(op.df, tmp)
        }
    }
    print(op.df) # for debugging
    return(op.df)
}

plot.marginal.mat <- function(df, my.ylim){
     # we are expecting a sub data frame for a given pair of marinal/context
     # dinucleotides.
    
    if(nrow(df) == 2){
        my.cols <- brewer.pal(3, "Reds")[2:3]
    } else if(nrow(df) == 1){
        my.cols <- brewer.pal(3, "Reds")[3]
    } else {
        my.cols <- brewer.pal(nrow(df), "Reds")
    }
    marginal.txt <- unique(df$marginal) # we expect only one value
    
    plot.mat <- as.matrix(df[, 2:17])
    rownames(plot.mat) <- rownames(df)
    
    opar <- par(no.readonly = T)
    par(mar = c(1, 1, 1, 1), cex.axis=0.3, lheight=0.3, cex.lab = 0.2)
    barplot(plot.mat, beside=T, col = my.cols, cex.names = 0.2, 
            main = marginal.txt, cex.main = 0.65, axisnames=F,
	    ylim = my.ylim) 
    legend("topleft", df$data, fill=my.cols, cex=0.3, bty="n")
}

make.plots <- function(op.df){
    split.screen(c(4, 4))
    marginals <- unique(op.df$marginal)
    marginals <- sort(marginals)
    print(length(marginals)) # for debugging
    my.ylim <- range(as.matrix(op.df[, 2:17]))
    for(i in 1:length(marginals)){
        print(i) # for debugging
        screen(i)
	sub.df <- op.df %>% filter(marginal == marginal[i])
	plot.marginal.mat(sub.df, my.ylim)
	close.screen(i)
    }
}

get.barplot.counts <- function(fn, window=100, mito=FALSE, strand = TRUE){
    print(fn) # for debugging and tracking progress -- working
    dt <- read.data(fn) %>% filter(flags2 > 49, flags7 == "-") # this stringency
                                                               # is hard-coded
							       # here.
    print(head(dt)) # for debugging
    dt.agg <- aggregate.on.window(dt, window=window, mito=mito, strand=strand)
    dt.counts <- as.data.frame(table(dt.agg$count))
    colnames(dt.counts) <- c("RPW", "NumberOfWindows")
    return(dt.counts)
}

get.barplot.counts.save <- function(pat, pat2 = "single", suff="*_single.txt", window = 100, 
                                    mito = FALSE, strand = TRUE){
    command <- paste(pat, suff, sep="/")
    command <- paste("ls", command)
    print(command) # for debugging
    fl <- system(command, intern=T)
    print(fl) # for debugging
    for(i in 1:length(fl)){
        print(fl[i]) # for debugging
	my.fn <- unlist(strsplit(fl[i], "/"))[2]
	my.type <- unlist(strsplit(my.fn, "_"))[1]
	print(my.type) # for debugging
        if(mito){
            ofn <- paste(paste(pat, pat2, window, "mito", my.type, sep="-"), ".counts", sep="")
	} else {
	    ofn <- paste(paste(pat, pat2, window, "genomic", my.type, sep="-"), ".counts", sep="")
	}
	print(ofn) # for debugging
	my.counts <- get.barplot.counts(fl[i], mito=mito, strand = strand, window = window)
	write.table(my.counts, ofn, sep="\t", col.names=T, row.names=F, quote=F)
    }
}

merge.barplot.counts <- function(pat, ofn){
    # Merge all the counts into a single file so that it is easier for Doug to
    # plot.
    fl <- system(paste("ls ", pat, sep=""), intern=T)
    names <- get.names(fl)
    for(i in 1:length(fl)){
        if(i == 1){
            df.1 <- read.delim(fl[i], header=T, stringsAsFactors=F)
	} else {
            df.2 <- read.delim(fl[i], header=T, stringsAsFactors=F)
	    df.1 <- merge(df.1, df.2, by="RPW", all=T)
	}
    }
    colnames(df.1) <- c("RPW", names)
    df.1[is.na(df.1)] <- 0 # lets replace all the NAs with 0s
    write.table(df.1, ofn, row.names=F, col.names=T, sep="\t", quote=F)
}

get.names <- function(v){
    return(unlist(lapply(v, FUN=function(x){return(unlist(strsplit(x, "-"))[1])})))
}

make.matrix <- function(df){
   my.vec <- df[, -1]
   my.tab <- as.table(my.vec)
   names(my.tab) <- as.character(df$RPW)
   return(my.tab)
}

make.window.barplots <- function(pat, ofn, suff = "*_single.txt", 
                                 window = 100, mito = FALSE, strand = TRUE){
    command <- paste(pat, suff, sep="/")
    print(command)
    command <- paste("ls ", command, sep="")
    print(command)
    fl <- system(command, intern=T)
    print(fl)
    op.l <- list()
    for(i in 1:length(fl)){
        counts.tmp <- get.barplot.counts(fl[i], window=window, mito=mito, strand=strand)
	op.l[[i]] <- counts.tmp
    }
    names(op.l) <- c("PuPu", "PuPy", "PyPu", "PyPy")
    if(mito){
        description = paste(pat, window, "bp, (Mitochondrial Genome Only)") 
    } else {
        description = paste(pat, window, "bp, (Genomic DNA only)")
    }
    pdf(ofn)
    layout(matrix(rbind(c(1, 1, 1),
                        c(2, 3, 4)), nrow = 2),
	   heights = c(1, 1), respect = T)
    my.max <- max(c(max(op.l$PuPu$NumberOfWindows),
                    max(op.l$PuPy$NumberOfWindows), 
		    max(op.l$PyPu$NumberOfWindows),
                    max(op.l$PyPy$NumberOfWindows)))
    my.ylims <- c(1, my.max)
    pypy.mat <- make.matrix(op.l$PyPy)
    barplot(pypy.mat, beside=T, 
            ylab = "Number of bins", xlab = "Number of events", main = "PyPy",
	    sub = description, cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8,
	    ylim = my.ylims)
    pupu.mat <- make.matrix(op.l$PuPu)
    barplot(pupu.mat,  beside=T,
            ylab = "Number of bins", xlab = "Number of events", main = "PuPu",
	    cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8,
	    ylim = my.ylims)
    pupy.mat <- make.matrix(op.l$PuPy)
    barplot(pupy.mat, beside=T, 
            ylab = "Number of bins", xlab = "Number of events", main = "PuPy",
	    cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8,
	    ylim = my.ylims)
    pypu.mat <- make.matrix(op.l$PyPu)
    barplot(pypu.mat, beside=T, 
            ylab = "Number of bins", xlab = "Number of events", main = "PyPu",
	    cex.lab = 0.8, cex.axis = 0.6, cex.main = 0.8,
	    ylim = my.ylims)
    dev.off()
}

make.window.column <- function(dt, windows = TRUE, window = 100){
    window.val <- -log10(window)
    dt <- dt[chr != "MT"] %>% filter(flags2 > 49, flags7 == "-") %>%
                        mutate(window = paste(chr, ":", round(damage.pos, digits = window.val), sep = ""),
                        unique.pos = paste(chr, damage.pos, strand, sep=":"))
			# filter(flags2 > 50) #, flags7 == "-")

    return(dt)
}

get.common.entries.windows <- function(dt, common.windows){
    # dt is a data.table
    dt.use <- dt %>% filter(window %in% common.windows) # %>%
    #                aggregate.on.window(mito = F, strand = T, window = 100) %>%
    #	             group_by(Window) %>% summarize(RPW = n())
    return(dt.use)
}

get.common.entries.bases <- function(dt, common.bases){
    dt.tmp <- dt %>% filter(unique.pos %in% common.bases) %>% 
                     transmute(window = window) %>% unique()
    dt.use <- merge(dt.tmp, dt, by = "window") %>% 
              mutate(mutual = ifelse(unique.pos %in% common.bases, "M", "X"))
    setcolorder(dt.use, c(colnames(dt.use)[2:ncol(dt.use)], colnames(dt.use)[1]))
    return(dt.use)
}

get.common.windows.generic <- function(fl, labl){
    common.bases <- c()
    dts.l <- c()
    op.l <- list()
    for(i in 1:length(fl)){
        var.name <- paste("dt", i, sep=".")
	dt <- make.window.column(read.data(fl[i]))
	if(i == 1){
            common.bases <- unique(dt$unique.pos)
	} else {
	    common.bases <- intersect(common.bases, unique(dt$unique.pos))
	    print(length(common.bases)) # for debugging
	}
	assign(var.name, dt)
	dts.l <- c(dts.l, var.name)
    }
    for(i in 1:length(dts.l)){
        var.name <- paste("dt.use", i, sep=".")
	dt.use   <- get.common.entries.bases(get(dts.l[i]), common.bases)
	dt.use$experiment <- labl[i]
	assign(var.name, dt.use)
	op.l[[var.name]] <- dt.use
    }
    return(op.l)
}

get.common.windows <- function(fn1, fn2, lab1, lab2, windows=TRUE){
    # we need information from both the data-sets.  We will have to provide a
    # label1 and label2 so that in the final data-frame we have a column that
    # tells from which data-set a given entry was obtained.  For the small
    # sequencing data-sets this should not be a big problem, for larger
    # data-sets it may be a problem, will need some work.
    dt.1 <- make.window.column(read.data(fn1))
    dt.2 <- make.window.column(read.data(fn2))
    # we want to take only uniquely mapping reads.
    if(windows == TRUE){
        common.windows <- intersect(dt.1$window, dt.2$window)
	print(length(common.windows)) # for debugging
	dt.1.use <- get.common.entries.windows(dt.1)
	dt.2.use <- get.common.entries.windows(dt.2)
    } else {
        common.bases <- intersect(dt.1$unique.pos, dt.2$unique.pos)
	dt.1.use <- get.common.entries.bases(dt.1, common.bases)
	dt.2.use <- get.common.entries.bases(dt.2, common.bases)
    }
    return(merge.commons(dt.1.use, dt.2.use, lab1, lab2))
}

merge.commons <- function(dt.1, dt.2, lab1, lab2){
    # We need to merge the data-sets so that we get all the information at one
    # glance.
    dt.1 <- window_summarize(dt.1, lab1)
}

window_summarize <- function(dt, lab, mito = FALSE){
    # we need to summarize each data table individually!
    # are are going to assume that the window is 100 bp, and that we are
    # ignoring the mitochondria.  We will revisit this later
    if(mito){
        return(NULL)
    } else {
        # create variable names for future use
        var.m <- paste(lab, "M", sep=".")
        var.x <- paste(lab, "X", sep=".")
	var.positions <- paste("positions", lab, sep=".")
        
        dt.tmp <- dt[chr != "MT"][, .(chr = chr, damage.window = round(damage.pos, -2), "M" = sum(which(mutual == "M")), "X" = sum(which(mutual == "X")), positions = paste0(unique(unique.pos), collapse = ",")), by = .(window)] %>%
	    group_by(window) %>% 
	    summarise(chr = unique(chr), damage.window = unique(damage.window), M = unique(M), X = unique(X), positions = unique(positions)) %>%
	    unique() %>% arrange(c(chr, damage.window))
        print(head(dt.tmp))
        setnames(dt.tmp, "M", var.m)
        setnames(dt.tmp, "X", var.x)
        setnames(dt.tmp, "positions", var.positions)
    }
    return(unique(dt.tmp))
} 

make.smooth.sgr <- function(df){
    all.smooth.df <- data.frame()
    for(i in c(1:22, "X", "Y")){
        df.chr <- df[chr == i]
	st <- seq(1, max(df.chr$damage.pos) - 100000, by = 900)
	ends <- st - 1 + 100000
	smooth.df <- data.frame()
	for(j in 1:length(st)){
            use.e <- df.chr %>% filter(damage.pos >= st[j],  damage.pos <= ends[j])
	    use.e <- use.e[, .(count = .N, pos = st[j]), by = .(chr, strand)]
	    smooth.df <- rbind(smooth.df, use.e)
        }
	# print(dim(smooth.df)) # for debugging
	setcolorder(smooth.df, c("chr", "pos", "strand", "count"))
	use.df <- smooth.df %>% transmute(chr = paste("chr", chr, sep=""), 
	                                  pos = pos, 
					  count = paste(strand, count, sep=""))
	all.smooth.df <- rbind(all.smooth.df, use.df)
    }
    return(all.smooth.df)
}

make.smooth.from.RMS <- function(fn){
    dt <- make.window.column(read.data(fn))
    ofn <- paste(fn, ".sgr", sep="")
    smooth.sgr.df <- make.smooth.sgr(dt)
    write.table(smooth.sgr.df, ofn, quote=F, sep="\t", col.names=F, row.names=F)
}

sel.threshold <- function(df, threshold = 1, nrow.limit = 10000){
    if(nrow(df %>% filter(Freq > threshold)) > nrow.limit){
        threshold <- threshold + 1
	sel.threshold(df, threshold = threshold)
	print(threshold)
    } else {
        return(threshold)
    }
}

select.windows <- function(dt, ofn){
    # we want to select only those windows that have more than one damage
    # events.  we are assuming single+recurrent file here.
    freq.tab <- as.data.frame(table(dt$window))
    colnames(freq.tab) <- c("window", "Freq")
    # threshold <- sel.threshold(freq.tab)
    freq.tab <- freq.tab %>% filter(Freq > 1) # only select rows where the
                                              # is greater than 1
    dt.sel <- merge(dt, freq.tab, by = "window") # select all entries from the
                                                 # main table.  We should be
						 # able to do this using dplyr
						 # filter or data.table, we will
						 # worry about it later.
    write.table(transmute(dt.sel, 
                          chr     = paste("chr", chr, sep=""),
			  start   = damage.pos - 50,
			  end     = damage.pos + 50,
			  peak.id = window),
	         ofn, quote=F, sep="\t", row.names=F, col.names=F)
}

call.annotate.peaks <- function(ifn, ofn, stats.fn, genome.go.dir, sample=""){
    # on ruddle call the annotate peaks pipeline.
    prog.name <- "annotatePeaks.pl"
    genome.go <- paste("-genomeOntology ", genome.go.dir)
    annStats  <- paste("-annStats", stats.fn)
    command <- paste(prog.name, ifn, "hg19", annStats, genome.go, ">", ofn)
    print(command)
    system(command)
}

select.genome.features <- function(ifn, basic=F){
    # 
    annot.dt <- fread(ifn, header=T)
    feature.col.names <- c("name", "annotation", "features", "coverage.bp",
                           "avg.feature.size", "overlap.peaks", "overlap.bp",
		  	   "expected.overlap", "log.ratio.enrichment",
			   "log.p.value", "p.value")
    setnames(annot.dt, colnames(annot.dt), feature.col.names)
    if(basic == T){
        annot.dt.sel <- annot.dt
    } else {
        annot.dt.sel <- annot.dt %>% filter(log.p.value < 0)
    }
    return(annot.dt.sel)

}

process.annotation <- function(v, split.char=" "){
    return(unlist(lapply(v, FUN=function(x){
                                   return(unlist(strsplit(x, split.char))[1])})))
}

process.detailed.annotation <- function(v, split.char = "\\|"){
    .process.function.2 <- function(x){
        if(x %in% c("Intergenic", "intron", "exon", "3'", "5'", "TTS", "non-coding", "promoter-TSS")){
            return(x)
        } else {
            return(unlist(strsplit(x, split.char))[3])
        }
    }   
    tmp <- process.annotation(v, split.char = " ")
    tmp <- unlist(lapply(tmp, FUN=.process.function.2))
    return(tmp)
}

.process.function.2 <- function(x, split.char){
    if(x %in% c("Intergenic", "intron", "exon", "3'", "5'", "TTS", "non-coding", "promoter-TSS")){
        return(x)
    } else {
        return(unlist(strsplit(x, split.char))[3])
    }
}

read.annotation.op <- function(ifn){
    # read the annotation
    my.col.names <- c("PeakID", "chr", "start", "end", "strand", "PeakScore",
                      "FocusRatio", "annotation", "detailed.annotation",  
		      "distance.to.tss", "nearest.prom.id", "entrez.id", 
		      "nearest.unigene", "nearest.refseq", "nearest.ensembl", 
		      "gene.name", "gene.alias", "gene.description", 
		      "gene.type")
    dt <- fread(ifn, header=T)
    setnames(dt, colnames(dt), my.col.names)
    dt$annot.brief <- process.annotation(dt$annotation)
    dt$annot2.breif <- process.detailed.annotation(dt$detailed.annotation, split.char="\\|")
    return(dt)
}

gather.data <- function(v){
   df.ret <- data.frame()
   for(i in 1:length(v)){
       fn <- paste(v[i], "annotation", sep=".")
       df <- read.annotation.op(fn)
       freq.df <- as.data.frame(table(df$annot2.breif))
       colnames(freq.df) <- c("Annotation", v[i])
       if(i == 1){
           df.ret <- freq.df
       } else {
           df.ret <- merge(df.ret, freq.df, by = "Annotation", all = T)
       }
   }
   return(df.ret)
}

gather.repeat.data <- function(v, basic=F){
    # df.ret <- data.frame()
    for(i in 1:length(v)){
        fn <- paste(v[i], "-genomeGO/repeats.genomeOntology.txt", sep="")
	df <- select.genome.features(fn) %>% select(name, log.ratio.enrichment)
	colnames(df) <- c("name", paste(gsub("-", ".", v[i]), "log.ratio.enrichment", sep="."))
	print(head(df)) # for debugging
	if(i == 1){
            df.ret <- df
	    print(head(df.ret))
	} else {
	    df.ret <- merge(df.ret, df, by = "name", all=T)  
	    print(head(df.ret))
	}
    }
    return(df.ret)
}

gather.basic.data <- function(v, basic=T){
    # df.ret <- data.frame()
    for(i in 1:length(v)){
        fn <- paste(v[i], "-genomeGO/basic.genomeOntology.txt", sep="")
	df <- select.genome.features(fn, basic = basic) %>% select(name, log.ratio.enrichment)
	colnames(df) <- c("name", paste(gsub("-", ".", v[i]), "log.ratio.enrichment", sep="."))
	print(head(df)) # for debugging
	if(i == 1){
            df.ret <- df
	    print(head(df.ret))
	} else {
	    df.ret <- merge(df.ret, df, by = "name", all=T)  
	    print(head(df.ret))
	}
    }
    return(df.ret)
}

get.all.annotation.RMS <- function(sample){
    sample.fn <- paste(sample, "PyPy-single-recurrent.txt", sep="/")
    sample.bed <- paste(sample, "bed", sep=".")
    stats.fn <- paste(sample, "stats", sep=".")
    annotation.fn <- paste(sample, "annotation", sep=".")
    genome.go.dir <- paste(sample, "genomeGo", sep="-")
    
    dt <- make.window.column(read.data(sample.fn))
    select.windows(dt, sample.bed)
    call.annotate.peaks(sample.bed, annotation.fn, stats.fn, genome.go.dir, sample=sample)
    genome.features <- select.genome.features(paste(genome.go.dir, basic.genomeOntology.txt, sep="/"))
    annot.df <- read.annotation.op(annotation.fn)
    op.l <- list()
    op.l[["annotation"]] <- annot.df
    op.l[["selected.features"]] <- genome.features
    return(op.l)
}


get.windowed.rpw <- function(v, windows, dinuc.v = c("PuPu", "PuPy", "PyPu", "PyPy")){
    for(j in 1:length(windows)){
        for(k in dinuc.v){
	    v.name <- paste(k, "-", windows[j], "-stats.txt", sep="")
	    print(v.name) # for debugging
	    tmp <- data.frame()
            for(i in 1:length(v)){
                fn <- paste(v[i], "/", k, "-single-recurrent.txt", sep="")
                dt <- get.barplot.counts(fn, window=windows[j])
                if(i == 1){
                    setnames(dt, colnames(dt), c("RPW", v[i]))
                    tmp <- dt
                } else {
                    setnames(dt, colnames(dt), c("RPW", v[i]))
                    tmp <- merge(tmp, dt, by = "RPW", all=T)
                }
            }
	    tmp$RPW <- as.numeric(as.character(tmp$RPW))
	    tmp[is.na(tmp)] <- 0
	    print(head(tmp)) # for debugging
	    print(v.name) # for debugging
	    write.table(arrange(tmp, RPW), v.name, sep="\t", row.names=F, col.names=T, quote=F)
        }
    }
}

select.reads.on.rpw <- function(fn, window = 100000, rpw.threshold = 75){
    dt <- read.data(fn) 
    print(dim(dt)) # for debugging
    dt <- dt[chr != "MT"] %>% filter(flags2 > 49,
                                     flags7 == "-") %>%
			      mutate(window = paste(chr, ":", round(damage.pos, digits = -log10(window)), sep = ""))
    print(dim(dt)) # for debugging
    sel.windows <- .select.windows.by.thresh(dt$window, rpw.threshold) %>% 
                   transmute(window = window)
    dt <- dt %>% merge(sel.windows, by = "window")
    print(dim(dt)) # for debugging
    return(dt)
}

.select.windows.by.thresh <- function(v, threshold){
    df <- as.data.frame(table(v))
    colnames(df) <- c("window", "Freq")
    df <- df %>% filter(Freq >= threshold)
    return(df)
}

get.interCPD.dist <- function(fn){
    dt <- read.data(fn) 
    print(dim(dt)) # for debugging
    dt <- dt[chr != "MT"] %>% filter(flags2 > 49, flags7 == "-")
    print(dim(dt)) # for debugging
    my.chrs <- unique(dt$chr)
    new.df <- data.frame()
    for(i in 1:length(my.chrs)){
         tmp <- dt[chr == my.chrs[i]] %>% arrange(damage.pos) %>% transmute(chr, damage.pos)
	 tmp <- tmp %>% transmute(diff.pos = diff(tmp$damage.pos))
	 if(i == 1){
             new.df <- tmp
	 } else {
             new.df <- rbind(new.df, tmp)
	 }
    }
    return(new.df)
}

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
    print(dim(summary.dt))
    print(c(fn.pypy.summary, fn.pupy, fn.pypu, fn.pupu)) # for debugging
    telescope.dt <- summary.dt
    setkey(telescope.dt, chr, strand, start, end)
    print(dim(telescope.dt))
    for(i in c(fn.pupy, fn.pypu, fn.pupu)){
        print(i) # for debugging
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

make.fns.pypu <- function(s, dir1 = "../analysis-sept192016/", suff2 = "-single-recurrent.txt"){
    pypu <- paste(dir1, s, "-PyPu-hotspots.txt", sep="")
    pypy <- paste(s, "/PyPy", suff2, sep="")
    pupy <- paste(s, "/PuPy", suff2, sep="")
    pupu <- paste(s, "/PuPu", suff2, sep="")
    return(c(pypu, pypy, pupy, pupu))
}

make.fns.pupy <- function(s, dir1 = "../analysis-sept192016/", suff2 = "-single-recurrent.txt"){
    pupy <- paste(dir1, s, "-PuPy-hotspots.txt", sep="")
    pypy <- paste(s, "/PyPy", suff2, sep="")
    pypu <- paste(s, "/PyPu", suff2, sep="")
    pupu <- paste(s, "/PuPu", suff2, sep="")
    return(c(pupy, pypu, pypu, pupu))
}

make.fns.pupu <- function(s, dir1 = "../analysis-sept192016/", suff2 = "-single-recurrent.txt"){
    pupu <- paste(dir1, s, "-PuPu-hotspots.txt", sep="")
    pypy <- paste(s, "/PyPy", suff2, sep="")
    pypu <- paste(s, "/PyPu", suff2, sep="")
    pupy <- paste(s, "/PuPy", suff2, sep="")
    return(c(pupu, pypy, pypu, pupy))
}

save.pure <- function(s){
    fns <- make.fns(s)
    ofn <- paste("../analysis-sept282016/", s, "-selected-pure-pypy.txt", sep="")
    telescope.dt <- telescoping(fns[1], fns[2], fns[3], fns[4])
    write.table(telescope.dt, ofn, row.names=F, col.names=T, sep="\t", quote=F)
}

save.pure.others <- function(s){
    fns.2 <- make.fns.pypu(s)
    ofn.2 <- paste("../analysis-nov142016/", s, "-selected-pure-pypu.txt", sep="")
    telescope.dt <- telescoping(fns.2[1], fns.2[2], fns.2[3], fns.2[4])
    write.table(telescope.dt, ofn.2, row.names=F, col.names=T, sep="\t", quote=F)
    
    fns.3 <- make.fns.pupy(s)
    ofn.3 <- paste("../analysis-nov142016/", s, "-selected-pure-pupy.txt", sep="")
    telescope.dt <- telescoping(fns.3[1], fns.3[2], fns.3[3], fns.3[4])
    write.table(telescope.dt, ofn.3, row.names=F, col.names=T, sep="\t", quote=F)

    fns.4 <- make.fns.pupu(s)
    ofn.4 <- paste("../analysis-nov142016/", s, "-selected-pure-pupu.txt", sep="")
    telescope.dt <- telescoping(fns.4[1], fns.4[2], fns.4[3], fns.4[4])
    write.table(telescope.dt, ofn.4, row.names=F, col.names=T, sep="\t", quote=F)
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
	print(dim(commmon.windows)) # for debugging
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

get.length.distribution <- function(fl, dinuc){
    exposure.source.df <- data.frame(sample = c("SP02HU", "SP03HU", "SP04HU",
                                            "SP05HU", "SP06HU", "SP10HU", 
					    "SP11HU", "SP12L16", "SP12L17", 
					    "SP13L16", "SP13L17", "SP14",
					    "SP15", "SP16"),
				 source = c("dna", "dna", "dna", "dna",
				            "dna", "cells", "cells", "cells",
					    "cells", "cells", "cells", "cells",
					    "cells", "cells"),
				 exposure = c("exposed", "exposed", "unexposed",
				              "exposed", "exposed", "unexposed",
					      "exposed", "unexposed",
					      "unexposed", "exposed", "exposed",
					      "unexposed", "exposed", "repair"),
				 uv = c("uvc", "uvc", "no", "uvc", "uvc", "no",
				        "uvc", "no", "no", "uvc", "uvc", "no",
					"uvc", "no")
				)

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
    all.summary <- my.df.2[, .(min.window.length = min(length),
                               perc25.window.length = quantile(length, 0.25),
			       median.window.length = quantile(length, 0.5),
			       perc75.window.length = quantile(length, 0.75),
			       IQR.window.length    = quantile(length, 0.75) - quantile(length, 0.25),
			       mean.window.length   = mean(length),
			       total.windows        = .N),
			       by = .(sample)]
    all.summary <- all.summary %>% mutate(expected.prob.median = (median.window.length / 6000000000) * total.windows,
			                  expected.prob.mean   = (mean.window.length / 6000000000) * total.windows)
    all.summary <- merge(all.summary, exposure.source.df, by = "sample")
    all.summary$sample <- as.character(all.summary$sample)
    all.summary$source <- as.character(all.summary$source)
    all.summary$exposure <- as.character(all.summary$exposure)
    all.summary$dinuc <- dinuc # so that we can aggregate the information for
                               # dinucleotide combinations later if required.
    return(list("summary" = all.summary, "all.data" = my.df.2))
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

.get.names <- function(summary.df, source = "dna", exposure = "exposed"){
    setDT(summary.df)
    setkey(summary.df, "source")
    my.names <- summary.df[source]
    setkey(my.names, "exposure")
    my.names <- my.names[exposure] %>% select(sample) %>% unlist() %>%
                as.vector()
    print(my.names)
    return(my.names)
}

make.barplot.merged <- function(summary.df, source = "dna", exposure = "exposed", strand = F, ...){
    
    my.names <- .get.names(summary.df, source = source, exposure = exposure)
    print(my.names) # for debugging
    
    merged.df <- get.overlaps.from.merged(my.names)
    print(head(merged.df)) # for debugging
    if(nrow(merged.df) == 0){
        print("No common windows detected")
	return()
    } else {

        expected  <- summary.df %>% 
                     filter(source == source, exposure == exposure) %>% 
		     select(total.windows, expected.prob.median, expected.prob.mean) %>% 
		     transmute(n.cells.mean = prod(expected.prob.mean) * sum(total.windows), 
		           n.cells.median = prod(expected.prob.median) * sum(total.windows)) %>% 
	              unique()
        print(expected) # for debugging
    
    ### Start plotting here
        if(strand){
	    df.to.plot <- merged.df %>% select(chr, strand, start, end) %>%
	                  unique()
            df.to.plot <- df.to.plot[, list(n.windows = .N), by = .(chr, strand)] %>%
	                  tidyr::spread(strand, n.windows) %>%
		          arrange(as.numeric(chr))
        } else {
	    df.to.plot <- merged.df %>% select(chr, strand, start, end) %>%
	                  unique()
	    setDT(df.to.plot)
            df.to.plot <- df.to.plot[, .(n.windows = .N), by = .(chr)] %>%
	                  arrange(as.numeric(chr))

        }
    

        if(strand){
	    rownames(df.to.plot) <- df.to.plot$chr
	    df.to.plot <- df.to.plot[, 2:3]
	    mat.to.plot <- as.matrix(df.to.plot)
            b <- barplot(t(mat.to.plot), col = c("blue", "red"), beside = T, 
	                 xlab = "Chromosomes", ...)
	    legend("topright", c("-", "+"), fill = c("blue", "red"),
	           bty = "n")
        } else {
            rownames(df.to.plot) <- df.to.plot$chr
            mat.to.plot <- df.to.plot %>% select(n.windows)
            mat.to.plot <- as.matrix(mat.to.plot)
            b <- barplot(mat.to.plot, names.arg = rownames(mat.to.plot), 
	            xlab = "Chromosomes", col = "grey", beside=T, ...)
        }
       
       print(mat.to.plot)
       
       text(x = 13, y = quantile(mat.to.plot, 0.89), paste("Expected Probability (mean) = ", expected[, 1]), col = "red", pos = 1)
       text(x = 13, y = quantile(mat.to.plot, 0.89), paste("Expected Probability (median) = ", expected[, 2]), col = "blue", pos = 3)
       # abline(h = expected[, 1], col = "blue", type = 2, lwd = 3)
       # abline(h = expected[, 2], col = "red", type = 2, lwd = 3)

    # text(nrow(mat.to.plot), y = expected[, c(1, 2)], col = c("red", "darkred"))
    }
}

## code on Nov 29, 2016
get.hotness <- function(windows.df, s.name){
    # we are only interested in the PyPy single+recurrent data right now.  We
    # are going back to the original outputs and filtering the data, and then
    # determining which window as how many damage events.  We have determined
    # that each window has at least 2 damage events for sure.  We want to now
    # rank the windows according to their hotness.
    dt <- read.data(paste(s.name, "/PyPy-single-recurrent.txt", sep=""))
    dt <- dt %>% transmute(chr = chr, strand = strand, start = damage.pos, end = damage.pos)
    setDT(dt)
    setDT(windows.df)
    setkey(windows.df, chr, strand, start, end)
    tmp.df <- foverlaps(dt, windows.df, type = "any", nomatch = 0L) 
    # the above step essentially labels each damage envent with the window.  We
    # will now count the nuber of times each window is hit.  We exact the final
    # dimensions to have same number of rows as windows.df.

    print(dim(tmp.df)) # for debugging

    tmp.df <- tmp.df %>% 
              mutate(label = paste(chr, ":", start, "-", end, ":", strand, sep=""))
    summary.df <- tmp.df %>% group_by(label) %>% summarise(n.events = n())
    colnames(summary.df) <- c("label", s.name)
    return(summary.df)
}

get.all.hotness <- function(windows.df, s.names){
    main.df <- windows.df
    main.df <- main.df %>%
               mutate(label = paste(chr, ":", start, "-", end, ":", strand, sep = ""))

    for(i in 1:length(s.names)){
        summary.df <- get.hotness(windows.df, s.names[i])
	main.df <- merge(main.df, summary.df, by = "label", all=T)
    }
    main.df <- main.df[, 2:ncol(main.df)]
    print(head(main.df)) # for debugging
    return(main.df)
}

## code in Jan 2017
get.windowed.damages.and.tag <- function(hotwindows.df, pypy.fn, label){
    # We want to get all the damage events that map to a specific window, tag the
    # damage events with the window id, and then whittle down the window to the
    # new minimum and maximum coordinates (strand awar    
    # 1 Read in the file
    damage.df <- read.data(pypy.fn) %>% select(chr, damage.pos, strand) %>%
                 mutate(start = damage.pos, end = damage.pos)

    # 2 Get a subset of damage events in each window
    setkey(hotwindows.df, chr, strand, start, end)
    # 3 Summarise the window to the new start, end 
    tmp.overlaps <- foverlaps(damage.df, hotwindows.df, type = "any", nomatch = 0L)
    tmp.overlaps <- tmp.overlaps %>% select(chr, strand, start, end, damage.pos)
    tmp.overlaps <- tmp.overlaps %>% mutate(experiment = label)
    return(tmp.overlaps)
}

get.whittled.windows <- function(hotwindows.df, s.names){
   # Get damages that map to only the hotwindows from all the experiments.  We
   # are only going to take PyPy-single-recurrent files.  We may change it later
   # to other types of damages.
   all.damages <- data.frame()
   for(i in 1:length(s.names)){
       label   <- s.names[i]
       pypy.fn <- paste(label, "/PyPy-single-recurrent.txt", sep="")
       tmp.overlaps <- get.windowed.damages.and.tag(hotwindows.df, pypy.fn, label)
       all.damages <- rbind(all.damages, tmp.overlaps)
   }
   setDT(all.damages)
   return(all.damages)
}

get.other.damages <- function(hotwindows.df, s.name){
    # Get all the damages from the given experiment that do not map to the given
    # set of coordinates
    pypy.fn <- paste(s.name, "/PyPy-single-recurrent.txt", sep="")
    damage.df <- read.data(pypy.fn) %>% select(chr, damage.pos, strand) %>%
                 mutate(start = damage.pos, end = damage.pos)
    damage.df$id <- 1:nrow(damage.df) # we are just assigning a number to each
                                      # row, we will then choose all the rows
				      # that did not map the given windows.
    setkey(hotwindows.df, chr, strand, start, end)
    tmp.overlaps <- foverlaps(damage.df, hotwindows.df, type = "any", nomatch = 0L)
    non.window.damages <- damage.df[!(id %in% tmp.overlaps$id)][chr != "MT"] %>%  setDF() %>%
                          select(chr, strand, damage.pos) %>%
			  mutate(label = s.name) %>% setDT()
    return(non.window.damages)
}

get.other.damages.all <- function(hotwindows.df, s.names){
    # gather information for all the damage events that lie outside of the
    # windows.

    my.other.damages <- data.frame()
    total.damages <- 0
    for(i in 1:length(s.names)){
        non.window.damages <- get.other.damages(hotwindows.df, s.names[i])
	per.10.kb <- (nrow(non.window.damages) / 6000000000) * 10000
	tmp.df <- data.frame(label = s.names[i], 
	                     n.damages = nrow(non.window.damages), 
			     damage.per.10kb = per.10.kb)
	total.damages <- total.damages + nrow(non.window.damages)
	my.other.damages <- rbind(my.other.damages, tmp.df)
    }
    return(my.other.damages)
}
