source("~/GitHub/CPDAnalysis/Code/R/pypy-analysis-functions.R")

dl <- system("ls", intern=T)
print(dl) # for debugging


get.comb.num <- function(dl, n, dinuc = "PuPu"){
    # We will generate number of common windows for given combination of files.
    # n is the number of cases to consider together.
    read.files <- c()
    stats.df <- data.frame()
    my.ind <- 1:length(dl)
    print(my.ind) # for debugging
    combn.mat <- combn(my.ind, n)
    for(i in 1:ncol(combn.mat)){
        file.inds <- combn.mat[, i]
	use.files <- paste(dl[file.inds], "/", paste(dinuc, "-single-recurrent.txt", sep=""), sep="")
	print(use.files) # for debugging
        files.used <- c()
        common.windows <- c()
	for(j in 1:length(use.files)){
	    fn <- use.files[j]
	    var.name <- paste(dl[file.inds[j]], dinuc, sep=".")
	    print(fn) # for debugging
	    if(!(var.name %in% read.files)){
                read.files <- c(read.files, var.name)
                windows <- unique(make.window.column(read.data(fn))$window)
		assign(var.name, windows)
	    } else {
                windows <- get(var.name)
	    }
	    if(j == 1){
                common.windows <- c(common.windows, windows)
	    } else {
                common.windows <- intersect(windows, common.windows)
	    }
	    files.used <- c(files.used, var.name)
	}
	n.common.windows <- length(common.windows)
	files.used       <- paste0(files.used, collapse=",")
	print(files.used) # for debugging
	stats.entry <- data.frame(files = files.used, common.windows = n.common.windows)
	stats.df <- rbind(stats.df, stats.entry)
	print(stats.df)
    }
    return(stats.df)
}

## all.df <- data.frame()
## for(i in c(2, 3, 4, 5)){
##     for(dinuc in c("PuPu", "PuPy", "PyPu", "PyPy")){
##         tmp.df <- get.comb.num(dl, i, dinuc=dinuc)
## 	tmp.df$n <- i
## 	tmp.df$dinuc <- dinuc
## 	all.df <- rbind(all.df, tmp.df)
##     }
## }

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
	all.smooth.df <- rbind(all.smooth.df, smooth.df)
    }
    return(all.smooth.df)
}



