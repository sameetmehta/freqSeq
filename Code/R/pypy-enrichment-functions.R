source("/home/sameet/GitHub/CPDAnalysis/Code/R/pypy-analysis-functions.R")

read.data.window <- function(fn){
    return(make.window.column(read.data(fn)))
}

select.windows <- function(dt, dinuc, ofn){
    # we want to select only those windows that have more than one damage
    # events.  we are assuming single+recurrent file here.
    freq.tab <- as.data.frame(table(dt$window))
    colnames(freq.tab) <- c("window", "Freq")
    freq.tab <- freq.tab %>% filter(Freq > 1) # only select rows where the
                                              # is greater than 1
    dt.sel <- merge(dt, freq.tab, by = "window") # select all entries from the
                                                 # main table.  We should be
						 # able to do this using dplyr
						 # filter or data.table, we will
						 # worry about it later.
    dt.sel$dinuc <- dinuc
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
