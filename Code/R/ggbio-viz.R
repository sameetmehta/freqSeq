library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggbio)
library(GenomicRanges)

make.karyograph <- function(df, title){
    # make sure that the df has three columns chr, start, and end.
    my.breaks <- fivenum(df$Z)
    data(hg19IdeogramCyto, package = "biovizBase")
    hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
    data(hg19Ideogram, package = "biovizBase")
    # setnames(df, c("Chr", "Start", "End"), c("chr", "start", "end"))
    setnames(df, colnames(df), gsub(" ", "", colnames(df)))
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
    seqlengths(gr) <- seqlengths(hg19Ideogram)[names(seqlengths(gr))]
    p <- ggplot(hg19) + layout_karyogram(cytoband = T)
    # p <- p + layout_karyogram(gr, aes(color = Hotness), ylim = c(0, 10),  geom = "rect") 
    p <- p + layout_karyogram(gr, aes(x = start, y = Z, col = Z), ylim = c(15, 30), alpha = 0.2, geom = "col") 
    p <- p + ggtitle(title) + theme_bw() + ylab("AU") + xlab("Position Along Chromosome") 
    p <- p + theme(strip.text.y = element_text(angle = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.y = element_text(size = 3.5), axis.ticks.y = element_blank())
    # p <- p + scale_color_manual("Z-score",  brewer.pal(5, "Reds"), breaks = my.breaks)
    # + scale_alpha("Hotness", palette = "Reds")
    return(p)
}


make.karyo.dev <- function(df1, df2, df3, title){
    # make sure that the df has three columns chr, start, and end.
    my.breaks <- fivenum(df$Z)
    data(hg19IdeogramCyto, package = "biovizBase")
    hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))
    data(hg19Ideogram, package = "biovizBase")
    setnames(df1, c("Chr", "Start", "End"), c("chr", "start", "end"))
    setnames(df2, c("Chr", "Start", "End"), c("chr", "start", "end"))
    setnames(df3, c("Chr", "Start", "End"), c("chr", "start", "end"))
    setnames(df1, colnames(df1), gsub(" ", "", colnames(df1)))
    setnames(df2, colnames(df2), gsub(" ", "", colnames(df2)))
    setnames(df3, colnames(df3), gsub(" ", "", colnames(df3)))
    gr1 <- makeGRangesFromDataFrame(df1, keep.extra.columns = T)
    gr2 <- makeGRangesFromDataFrame(df2, keep.extra.columns = T)
    gr3 <- makeGRangesFromDataFrame(df3, keep.extra.columns = T)
    seqlengths(gr1) <- seqlengths(hg19Ideogram)[names(seqlengths(gr1))]
    seqlengths(gr2) <- seqlengths(hg19Ideogram)[names(seqlengths(gr2))]
    seqlengths(gr3) <- seqlengths(hg19Ideogram)[names(seqlengths(gr3))]
    p <- ggplot(hg19) + layout_karyogram(cytoband = T)
    p <- p + layout_karyogram(gr1, ylim = c(11, 13),  geom = "rect", col = rgb(1, 0, 0, 0.5)) 
    p <- p + layout_karyogram(gr2, ylim = c(14, 17),  geom = "rect", col = rgb(0, 0, 1, 0.5)) 
    p <- p + layout_karyogram(gr3, ylim = c(18, 20),  geom = "rect", col = rgb(0, 1, 0, 0.5)) 
    # p <- p + layout_karyogram(gr, aes(x = start, y = Z, col = Z), ylim = c(15, 30), alpha = 0.2, geom = "col") 
    p <- p + ggtitle(title) + theme_bw() + ylab("AU") + xlab("Position Along Chromosome") 
    p <- p + theme(strip.text.y = element_text(angle = 0), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.y = element_text(size = 3.5), axis.ticks.y = element_blank())
    # p <- p + scale_color_manual("Z-score",  brewer.pal(5, "Reds"), breaks = my.breaks)
    # + scale_alpha("Hotness", palette = "Reds")
    return(p)
}
