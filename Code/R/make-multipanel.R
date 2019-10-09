library(grid)
library(gridBase)
library(gridExtra)
library(RColorBrewer)

make.stranded.venns <- function(s.agg, exp = "SP04", window = "10 kb", mito = FALSE){# s.agg stands for stranded aggregated data.table
    pushViewport(viewport(x = 0, y = 0, width = unit(1, "npc"), height =
    unit(0.95, "npc"), just = c("left", "bottom"), name = "usablearea"))
    total.venns <- unique(s.agg$count) #total number of Venns to be plotted
    total.venns <- sort(total.venns)
    vname <- paste("vd", 1:length(total.venns), sep=".")
    layout.dim <- ceiling(sqrt(length(total.venns)))
    layout.mat <- matrix(rep(1, layout.dim ** 2), nrow = layout.dim)
    vplay <- grid.layout(layout.dim, layout.dim, 
                         widths = unit(layout.mat/layout.dim, 
			               matrix(rep("npc", layout.dim ** 2), ncol = layout.dim)
				      ),
                         heights = unit(layout.mat/layout.dim,
			               matrix(rep("npc", layout.dim ** 2), ncol = layout.dim)
				      ),
                         respect = matrix(rep(1, layout.dim ** 2), 
			                  ncol = layout.dim)
			)
    pushViewport(viewport(layout=vplay))
    lc <- 1
    lr <- 1
    for(i in 1:length(total.venns)){
        print(vname[i])
        tmp.pos <- make.unique.vector(s.agg[count == total.venns[i]][strand == "+"])
        tmp.neg <- make.unique.vector(s.agg[count == total.venns[i]][strand == "-"])
	tmp.main <- paste(total.venns[i], "RPW")
        tmp.sub <- paste("Total N = ", length(union(tmp.pos, tmp.neg)))
	tmp.vd <- venn.diagram(list(tmp.pos, tmp.neg), 
	                       category.names = c("+", "-"), 
			       fill = c("purple", "orange"), 
			       alpha = c(0.45, 0.45), col=NA,
	                       filename = NULL, main = tmp.main, 
			       sub = tmp.sub, print.mode = "percentage",
			       cex = 0.2, main.cex = 0.4, sub.cex=0.3)
	pushViewport(viewport(layout.pos.col = lc, layout.pos.row = lr))
	grid.rect()
	grid.draw(tmp.vd); popViewport()
	lc <- lc + 1
	if(!(layout.dim %/% lc)&&(layout.dim > 2)){
	    lr <- lr + 1
	    lc <- 1
	}
    }
    popViewport(2)
    if(mito){
        my.title <- paste(exp, " (only Mitochondria) ", "(", window, ")", sep="")
    } else {
        my.title <- paste(exp, " (no Mitrochondria) ", "(", window, ")", sep="")
    }
    pushViewport(viewport(x = 0, y = unit(0.95, "npc"), width = unit(1, "npc"), 
                          height = unit(0.05, "npc"), just = c("left", "bottom"), 
			  name="titlevp"))
    grid.text(my.title, x = unit(0.5, "npc"), y = unit(0.5, "npc"), 
              just = c("center", "center"), gp = gpar(cex = 1.2))
    popViewport()
}
