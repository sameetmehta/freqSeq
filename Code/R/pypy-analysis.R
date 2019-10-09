options(scipen=9999) # we do not want scientific notation

library(data.table) # read and manipulate large files
library(dplyr) # easily summarise data frames.
library(VennDiagram)
library(RColorBrewer)
library(RCurl) # for the merge.list function
library(stringr) # we are going to do something with string manipulation

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

make.s.aggregate.l <- function(fn, strand=T, mito=F, window=1000){
    # Create a stranded aggregate data.table, that can be used to generated a
    # multipanel Venn diagram
    dt <- aggregate.on.window(read.data(fn), window=window, mito=mito, strand=strand)
    return(dt)
}

make.ready.for.circos <- function(dt, ofn, window, strand=FALSE){
   # this function will generate a table ready to used in CIRCOS, and write it
   # to the fn.  We are using the "window" only to calculate the bin, we
   # actually do not need the window value itself for anything else
   if(nrow(dt) < 25000){
       if(strand){
           write.table(dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count, strand = strand), ], 
                      ofn, quote=F, row.names=F, col.names=F)
           print(paste("File ", ofn, " written"))
        } else {
           write.table(dt[, .(chr = paste("hs", Chr, sep=""), start = ifelse(Window == 0, 0, Window - (window/2)), end = Window + ((window/2) - 1), count = count), ], 
                      ofn, quote=F, row.names=F, col.names=F)
           print(paste("File ", ofn, " written"))
	}
    }
    else { 
        print("Table has too many rows, that cannot be plotted with CURRENT CIRCOS SETTINGS (25000 Max per track)")
    }
}

gather.data <- function(fn, window = 10000, mito = FALSE, strand = FALSE){
    # Same logic as before, we need to refactor this to make it more palatable
    # and runnable on the CLUSTER.  We will deal with one file at a time.
    dt <- read.data(fn)
    dt.aggregate <- aggregate.on.window(dt, window = window, mito = mito, strand=strand)
    return(dt.aggregate)
}

gather.data.for.sample <- function(pat, window = 10000, mito=FALSE, strand = FALSE){
    # return a data-frame for a pattern and window
    fl <- get.list(pat)
    
    pupu.aggregate <- gather.data(fl[1], window = window, mito = mito, strand = strand)
    pupy.aggregate <- gather.data(fl[2], window = window, mito = mito, strand = strand)
    pypu.aggregate <- gather.data(fl[3], window = window, mito = mito, strand = strand)
    pypy.aggregate <- gather.data(fl[4], window = window, mito = mito, strand = strand)

    if(mito){
    	description <- paste(unlist(strsplit(pat, "/"))[1], window, "bp, (Mitochondrial genome only)")
    } else {
	description <- paste(unlist(strsplit(pat, "/"))[1], window, "bp, (Genomic DNA only)")
    }

    op.l <- list(pupu = table(pupu.aggregate$count), 
                 pupy = table(pupy.aggregate$count),
		 pypu = table(pypu.aggregate$count),
		 pypy = table(pypy.aggregate$count),
		 description = description)
    
    # now plot these numbers!
    return(op.l)
}

make.panel.plot <- function(op.l){
    layout(matrix(rbind(c(1, 1, 1),
                        c(2, 3, 4)), nrow = 2),
	   heights = c(1, 1), respect=T)
    my.max <- max(c(op.l$pypy, op.l$pupu, op.l$pypu, op.l$pupy))
    my.ylims <- c(1, my.max)

    ##TODO!! see if we can actually print the numbers on the top of the bars, so
    ## so that the bars are more readable.  This information will be similar to
    ## 4-way Venn-diagram
    barplot(op.l$pypy, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PyPy", sub = op.l$description, cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pupu, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PuPu", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pupy, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PuPy", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
    barplot(op.l$pypu, 
            ylab = "number of bins", xlab = "Reads", 
	    main = "PyPu", cex.lab = 0.8, cex.axis=0.6, cex.main = 0.8, ylim = my.ylims)
}



make.boxplot <- function(l, window, use.chromosome = NA){
   # we want to me a grouped boxplot
   # the "l" contains the list of patterns, so that we can get the files for
   # each pattern as in the venn-diagram function, but we want to stack them
   # into a data-frame that will allow us to make a boxplot.  We are basically
   # going to calculate a fivenum summary

   for(i in 1:length(l)){
       # First gather count data for all 4 compbinations for one pattern
       my.df <- gather.data.for.sample(l[[i]], use.chromosome = use.chromosome, window = window)
   }

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

make.marginal.count.df <- function(pat, mito=FALSE){
    command <- paste(pat, "/*_single.txt", sep="")
    command <- paste("ls ", command, sep="")
    fl <- system(command, intern=T)
    print(fl) # for debugging
    for(i in 1:length(fl)){
        print(fl[i]) # for debugging
        tmp <- get.marginal.dt(fl[i], mito=mito)
	tmp2 <- tidyr::spread(tmp, damage.dimer, freq)
	print(head(tmp2)) # for debugging
	if(i == 1){
            op.df <- tmp2
	} else {
            op.df <- merge(op.df, tmp2, by = "marginal")
	}
    }
    print(head(op.df)) # for debugging
    return(op.df)
}

make.dimer.count.vector <- function(pat, mito=FALSE){
    command <- paste(pat, "/*_single.txt", sep="") # get sample name, and obtain
                                                  # list of all the files
    command <- paste("ls ", command, sep="")
    fl <- system(command, intern=T)
    print(fl) # for debugging
    damage.vec <- c()
    dimer.names <- c()
    for(f in fl){
        dt <- read.data(f)
	dt.tmp <- make.dimer.count.dist(dt, mito=mito)
	dt.tmp <- filter(dt.tmp, !(grepl("N", damage.dimer))) # there were some
	                                                      # cases where
							      # there were "N"
							      # in the damage
							      # dimer. we are
							      # accounting of
							      # that here.
	dimer.counts <- table(dt.tmp$damage.dimer)
	dimer.names <- c(dimer.names, names(dimer.counts))
	damage.vec <- c(damage.vec, as.vector(dimer.counts))
    }
    names(damage.vec) <- dimer.names
    print(damage.vec) # for debugging
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

make.denom <- function(n){
    num.len <- str_length(as.character(n))
    nearest <- num.len - 1
    return(10 ** nearest)
}

get.normalization.factor <- function(pat){
    fn <- paste(pat, "/countlog.txt", sep="")
    df <- read.delim(fn, skip=1, sep="", stringsAsFactors=F, header=F)
    correction.factor <- as.numeric(df$V3[1]) - (sum(as.numeric(df$V2[2:6])))    
    # this is a hack, but right now I dont know a better approach.

    my.denom <- make.denom(correction.factor) # this part will not be used for
                                              # now

    correction.factor <- correction.number / 1000000 # we will right now set it
                                                     # for per million used
						     # reads

    print(correction.factor) # for debugging # this is only to see what kind of numbers we are getting.
    
    ##TODO!! We need to find a way to  to report this information, what the denominator was.
    return(correction.factor) 
}

plot.dimer.matrix <- function(dimer.count.mat, ofn){
    my.cols.1 <- c(brewer.pal(6, "Blues")[3:6], 
                   brewer.pal(6, "Greens")[3:6],
                   brewer.pal(6, "Purples")[3:6],
		   brewer.pal(6, "Reds")[3:6])

    my.cols.2 <- c(brewer.pal(6, "Reds")[3:5])
    pdf(ofn)
    split.screen(c(2, 1))
    screen(1)
    barplot(dimer.count.mat, beside=T, col = my.cols.2, ylab = "Frequency", 
            xlab = "Dimer", cex.axis=0.8, cex.lab=0.8, cex.names=0.6)
    legend("topleft", rownames(dimer.count.mat), fill = my.cols.2, bty="n",
           cex=0.6)
    close.screen(1)
    screen(2)
    barplot(t(dimer.count.mat), beside=T, col = my.cols.1, ylab = "Frequency",
            xlab = "Experiment", cex.axis=0.6, cex.lab=0.6)
    legend("topright", colnames(dimer.count.mat), fill = my.cols.1, ncol=4, bty="n", cex = 0.5)
    close.screen(2)
    dev.off()
}

get.normalization.factor <- function(pat){
## 
##             V1     V2     V3
## 	    1       Linker  Match 896372
## 	    2     Unmapped  56447     NA
## 	    3     Multimap  46662     NA
## 	    4  Wronglength  89469     NA
## 	    5        Other  25459     NA
## 	    6   Deanimated    279     NA
    fn <- paste(pat, "/countlog.txt", sep="")
    df <- read.delim(fn, skip=1, sep="", stringsAsFactors=F, header=F)
    # print(df) # for debugging # now working
    # print(as.numeric(df$V3[1])) # for debugging # now working
    # print(as.numeric(df$V2[2:6])) # for debugging # now working
    correction.factor <- as.numeric(df$V3[1]) - (sum(as.numeric(df$V2[2:6])))    
    # print(correction.number) # for debugging # now working
       
    correction.factor <- correction.number / 1000000 # so many million reads. If
                                                     # total is less than a
						     # million, we will end up
						     # "increasing the number
    ##TODO!! WE NEED TO ADD A FUNCTIONALITY HERE THAT WILL DECIDE ON THE
    # NORMALIZATION factor based on the number itself.  But we want this to be
    # constant for all the three conditions.
    return(correction.factor) 
}


