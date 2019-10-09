library(data.table)
library(tidyverse)
library(dtplyr)

get.dimer <- function(df){
  dimer.df <- df %>% select(damage.match) %>% unlist() %>%
              as.vector() %>% table() %>% as.data.frame()
  colnames(dimer.df) <- c("dimer", "Freq")
  dimer.df <- dimer.df %>% filter(!grepl("N", dimer))
  sum(dimer.df[, 2])
}

make.fl <- function(suffix = "PyPy", test = 0){
  pattern <- paste(suffix, "-pos.txt", sep = "")
  print(pattern)
  fl <- list.files(pattern = pattern)
  print(fl)

  if(test){
      fl = fl[1:5]
  }

  fl
}

make.saturation.df.all <- function(fl, window = 1, test = T, tmp.file = "intermediate.table.txt"){
  snames <- unlist(lapply(fl, function(x)paste0(unlist(strsplit(x, "-"))[1:2], collapse = "-")))
  for(i in seq_along(fl)){
      sname <- snames[i]
      # print(sname) # for debugging
      
      df <- fread(fl[i], header = T)
      # print(head(df)) # for debugging

      s.dimer.counts <- get.dimer(df)
      # print(dimer.counts) # for debugging
      
      if(window > 1){
          round.v <- -log10(window)
	  df <- df %>% mutate(round.pos = round(damage.pos, round.v))
	  print(head(df)) # for debugginga
	  my.windows <- df %>% select(chr, round.pos, strand) %>% unique() %>% nrow()
          my.unique <- my.windows

      } else {
          my.windows <- nrow(df)
	  df <- unique(df)
	  my.unique <- nrow(df)
      }

      if(i == 1){
          all.df <- df
	  cum.diN.count <- s.dimer.counts
	  all.seen <- nrow(df)
          
	  if(window > 1){
              all.df.op <- df %>% select(chr, round.pos, strand) %>% unique()
	      all.seen <- nrow(all.df.op)
	      # print(head(all.df.op))
	  }

      } else {
          all.df <- merge(all.df, df[, 1:4] %>% unique(), by = c("chr", "damage.pos", "strand", "damage.match"), all = T)[, 1:4]
	  cum.diN.count <- get.dimer(all.df)
	  all.seen <- unique(all.df) %>% nrow()
 	  if(window > 1){
              all.df.op <- rbind(all.df.op, df %>% select(chr, round.pos, strand)) %>% 
	                   unique()
	      all.seen <- nrow(all.df.op)
	 }

      }
      
      my.e <- data.frame(sample = sname,
                         n.windows = my.windows,
		         n.unique  = my.unique,
		         total.seen = all.seen,
		         sample.diN = s.dimer.counts,
		         cum.diN    = cum.diN.count)
       print(my.e)

       if(i == 1){
           all.nums <- my.e
       } else {
          all.nums <- rbind(all.nums, my.e)
       }
      print(all.nums)
      write.table(all.nums, tmp.file, row.names = F, col.names = T, quote = F, sep = "\t")
  }
 
  all.nums
}

make.saturation.df <- function(fl, suffix = "PyPy", window = 1, test = T){
  
  pattern =  paste("-", suffix, "-pos.txt", sep = "")
  
  for(i in seq_along(fl)){
      sname <- gsub(pattern, "", fl[i])
      # print(sname) # for debugging
      
      df <- fread(fl[i], header = T)
      # print(head(df)) # for debugging

      s.dimer.counts <- get.dimer(df)
      # print(dimer.counts) # for debugging
      
      if(window > 1){
          round.v <- -log10(window)
	  df <- df %>% mutate(round.pos = round(damage.pos, round.v))
	  print(head(df)) # for debugginga
	  my.windows <- df %>% select(chr, round.pos, strand) %>% unique() %>% nrow()
          my.unique <- my.windows

      } else {
          my.windows <- nrow(df)
	  df <- unique(df)
	  my.unique <- nrow(df)
      }

      if(i == 1){
          all.df <- df
	  cum.diN.count <- s.dimer.counts
	  all.seen <- nrow(df)
          
	  if(window > 1){
              all.df.op <- df %>% select(chr, round.pos, strand) %>% unique()
	      all.seen <- nrow(all.df.op)
	      # print(head(all.df.op))
	  }

      } else {
          all.df <- merge(all.df, df[, 1:4] %>% unique(), by = c("chr", "damage.pos", "strand", "damage.match"), all = T)[, 1:4]
	  cum.diN.count <- get.dimer(all.df)
	  all.seen <- unique(all.df) %>% nrow()
 	  if(window > 1){
              all.df.op <- rbind(all.df.op, df %>% select(chr, round.pos, strand)) %>% 
	                   unique()
	      all.seen <- nrow(all.df.op)
	 }

      }
      
      my.e <- data.frame(sample = sname,
                         n.windows = my.windows,
		         n.unique  = my.unique,
		         total.seen = all.seen,
		         sample.diN = s.dimer.counts,
		         cum.diN    = cum.diN.count)
       print(my.e)

       if(i == 1){
           all.nums <- my.e
       } else {
          all.nums <- rbind(all.nums, my.e)
       }
      print(all.nums)
  }
 
  all.nums
}
##      sname <- gsub(pattern, "", fl[i])
##      print(sname)
##      # samples <- c(samples, sname)
##      df <- fread(fl[i], header = T)
##      
##      dimer.df <- get.dimer(df)
##      colnames(dimer.df) <- c("dimer", sname)
##      dimer.df <- dimer.df %>% filter(!grepl("N", dimer))
##      dimer.counts <- sum(dimer.df[, 2])
##      # dimer.counts.sample <- c(dimer.counts.sample, dimer.counts)
##      
##      if(window > 1){
##          round.v <- -log10(window)
## 	 print(round.v)
## 	 df <- df %>% mutate(round.pos = round(damage.pos, round.v))
## 	 print(head(df))
## 	 # my.nrows <- c(my.nrows, df %>% select(chr, round.pos) %>% nrow())
## 	 n.windows <- df %>% select(chr, round.pos) %>% nrow()
## 	 # my.uniques <- c(my.uniques, df %>% select(chr, round.pos) %>% unique() %>% nrow())
## 	 n.unique  <- df %>% select(chr, round.pos) %>% unique() %>% nrow()
## 	 print(c(n.windows, n.unique))
##      } else {
##          # my.nrows <- c(my.nrows, nrow(df))
## 	 n.windows <- nrow(df)
## 	 df <- unique(df)
##          n.unique <-  nrow(df)
##      }
##      
##      if(i == 1){
##          all.df <- df
## 	 dimer.counts.cum <- dimer.counts
##          if(window > 1){
## 	     all.df.op <- df %>% select(chr, round.pos) %>% unique()
## 	     print(head(all.df.op))
## 	 }
##      } else {
##          all.df <- merge(all.df, df, by = c("chr", "damage.pos", "strand", "damage.match"), all = T)[, c(1, 2, 3, 4)]
##          cum.dimer <- get.dimer(all.df)
## 	 colnames(cum.dimer) <- c("dimer", "count")
## 	 # dimer.counts.cum <- c(dimer.counts.cum, sum(cum.dimer$count))
## 	 dimer.counts.cum <- sum(cum.dimer$count)
## 	 
## 	 if(window > 1){
##              all.df %>% select(chr, round.pos) %>% rbind(all.df.op) -> all.df.op
## 	     all.df.op %>% unique() -> all.df.op
## 	 }
##      }
##      
##      if(window > 1){
##          # my.totals <- c(my.totals, nrow(all.df.op))
## 	 my.totals <- nrow(all.df.op)
##      } else {
##          my.totals <- nrow(all.df)
##      }
##      
##      print(fl[i])
##      
##      my.e <- data.frame(sample = sname,
##                         n.windows = n.windows,
## 			n.unique = n.unique,
##                         diN.sample = dimer.counts,
##                         total.seen = my.totals,
## 			total.diN = dimer.counts.cum)
##      all.nums <- rbind(all.nums, my.e)
##      print(all.nums)
##                         
##      # print(all.nums)
##   }
##   
## ##   all.nums <- data.frame(sample = samples,
## ##                          n.windows = my.nrows,
## 
## ## 			 unique.windows = my.uniques,
## ## 			 diN.sample = dimer.counts.sample,
## ## 			 total.seen = my.totals,
## ## 			 total.diN  = dimer.counts.cum)
##   all.nums
## }

prefs <- c("PuPu", "PuPy", "PyPu", "PyPy")
wins  <- c(1, 10, 100, 1000, 10000, 100000, 1000000)

## for(i in seq_along(prefs)){
##   for(j in seq_along(wins)){
##     ofn = file.path(paste(prefs[i], "-", wins[j], "-saturation-numbers.txt", sep = ""))
##     print(ofn)
##     fl <- make.fl(prefs[i], test = 1)
##     sat.df <- make.saturation.df(fl, prefs[i], wins[j], test = 1)
## #    write.table(sat.df, ofn, col.names = T, row.names = F, sep = "\t", quote = F)
##   }
## }
