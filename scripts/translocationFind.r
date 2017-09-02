#This code can be called from another Rscript as e.g. matseg(mat = matrix, lambda = 1.12, method = "BinSeg")
#This will produce a set of segments with their average signal value and corresponding p.value 
library(changepoint)
matseg <- function(mat,lambda,method){ 
	chrAposStart <- c()
	chrBposStart <- c()
	chrBposEnd <- c()
	m.mean <- cpt.mean(mat, method = method)
	for(i in 1:dim(mat)[1]){
		start <- 1
		for(j in 1:length(m.mean[[i]]@cpts)){
			end <- m.mean[[i]]@cpts[j]
			chrAposStart <- c(chrAposStart,i)
			chrBposStart <- c(chrBposStart,start)
			chrBposEnd <- c(chrBposEnd,end)
			start <- m.mean[[i]]@cpts[j]+1	
		}
	}
	mean <- c()
	max <- c()
	percnt_sig <- c()
	pval.mean <- c()
	pval.max <- c()
	for (i in 1:length(chrAposStart)){
		mean <- c(mean,mean(mat[chrAposStart[i],chrBposStart[i]:chrBposEnd[i]]))
		pval.mean <- c(pval.mean,(1-ppois(mean[i],lambda=lambda)))
		max <- c(max,max(mat[chrAposStart[i],chrBposStart[i]:chrBposEnd[i]]))
		pval <- 1-ppois(mat[chrAposStart[i],chrBposStart[i]:chrBposEnd[i]],lambda=lambda)
		percnt_sig <- c(percnt_sig,(length(which(p.adjust(pval,method="BH") <= 0.05))))
		pval.max <- c(pval.max,(1-ppois(max[i],lambda=lambda)))
	}
	pval.mean.adj <- p.adjust(pval.mean,method="BH")
	pval.max.adj <- p.adjust(pval.max,method="BH")
	pval.percnt_sig.adj <- p.adjust(pnorm(-abs(percnt_sig-mean(percnt_sig)/sd(percnt_sig))),method="BH") 	
	info <- data.frame(chrAposStart,chrBposStart,chrBposEnd,mean,pval.mean.adj,pval.percnt_sig.adj,percnt_sig)			
	return (info)
}

args = commandArgs(trailingOnly=TRUE)
u<-read.table(args[1],header=FALSE)
u.mat <- as.matrix(u)*10
u.vec <- as.vector(u.mat)
lambda <- mean(u.vec[u.vec > 0])
chrA <- matseg(u.mat,lambda,"BinSeg")
chrB <- matseg(t(u.mat),lambda,"BinSeg")

dataA <- paste(args[2],"/","dataA.txt",sep="")
dataB <- paste(args[2],"/","dataB.txt",sep="")
write.table(chrA,file=dataA)
write.table(chrB,file=dataB)

filt <- subset(chrA,pval.mean.adj <= 0.05 & pval.percnt_sig.adj <= 0.05)
if (length(filt$chrAposStart) < 2){
	filt <- subset(chrA,pval.mean.adj <= 0.05)
}
chrA.filt <- filt$chrAposStart
chrA.diff <- diff(chrA.filt)
chrA.pval.adj <- p.adjust((1-ppois(chrA.diff,lambda=mean(chrA.diff))),method="BH")
chrA.pval.significant <- c(which(chrA.pval.adj <= 0.05),length(chrA.pval.adj)+1)

chrA.start <- c()
chrA.end <- c()
if (length(chrA.pval.significant) != 1 || (length(chrA.pval.significant) == 2 && chrA.pval.significant[1] != 1)){
	i <- 1
	while (i <= length(chrA.pval.significant)){
	        if (i == 1){
	                start <- 1
	                end <- chrA.pval.significant[i]
			chrA.start <- c(chrA.start,chrA.filt[start])
			chrA.end <- c(chrA.end,chrA.filt[end])
			if (length(chrA.pval.significant) == 1){
	                        start <- chrA.pval.significant[i]+1
	                        end <- length(chrA.pval.adj)+1
				chrA.start <- c(chrA.start,chrA.filt[start])
		                chrA.end <- c(chrA.end,chrA.filt[end])
	                }        
	        }else {
                	start <- chrA.pval.significant[i-1]
                	end <- chrA.pval.significant[i]
                	if (abs(start-end) > 1){
                	        start <- chrA.pval.significant[i-1]+1
                	        end <- chrA.pval.significant[i]
				chrA.start <- c(chrA.start,chrA.filt[start])
	        	        chrA.end <- c(chrA.end,chrA.filt[end])
                	}
        	}	
        	i <- i+1
	}
}else {
	if (length(chrA.pval.significant) == 1){
		chrA.start <- c(chrA.start,chrA.filt[1])                        
		chrA.end <- c(chrA.end,chrA.filt[length(chrA.filt)])
	}else if (length(chrA.pval.significant) == 2 && chrA.pval.significant[1] != 1){
		chrA.start <- c(chrA.start,chrA.filt[2])
                chrA.end <- c(chrA.end,chrA.filt[length(chrA.filt)])
	}
}
filt <- subset(chrB,pval.mean.adj <= 0.05 & pval.percnt_sig.adj <= 0.05)
if (length(filt$chrAposStart) < 2){
        filt <- subset(chrB,pval.mean.adj <= 0.05)
}
chrB.filt <- filt$chrAposStart
chrB.diff <- diff(chrB.filt)
chrB.pval.adj <- p.adjust((1-ppois(chrB.diff,lambda=mean(chrB.diff))),method="BH")
chrB.pval.significant <- c(which(chrB.pval.adj <= 0.05),length(chrB.pval.adj)+1)

chrB.start <- c()
chrB.end <- c()
if (length(chrB.pval.significant) != 1 || (length(chrB.pval.significant) == 2 && chrB.pval.significant[1] != 1)){
	i <- 1
	while (i <= length(chrB.pval.significant)){
	        if (i == 1){
	                start <- 1
	                end <- chrB.pval.significant[i]
			chrB.start <- c(chrB.start,chrB.filt[start])
	                chrB.end <- c(chrB.end,chrB.filt[end])
			if (length(chrB.pval.significant) == 1){
	                        start <- chrB.pval.significant[i]+1
	                        end <- length(chrB.pval.adj)+1
	                        chrB.start <- c(chrB.start,chrB.filt[start])
	                        chrB.end <- c(chrB.end,chrB.filt[end])
	                }
	        }else {
	                start <- chrB.pval.significant[i-1]
	                end <- chrB.pval.significant[i]
	                if (abs(start-end) > 1){
	                        start <- chrB.pval.significant[i-1]+1
	                        end <- chrB.pval.significant[i]
				chrB.start <- c(chrB.start,chrB.filt[start])
	                        chrB.end <- c(chrB.end,chrB.filt[end])
	                }
	        }
	        i <- i+1
	}
}else {
	if (length(chrB.pval.significant) == 1){
		chrB.start <- c(chrB.start,chrB.filt[1])
        	chrB.end <- c(chrB.end,chrB.filt[length(chrB.filt)])
	}else if (length(chrB.pval.significant) == 2 && chrB.pval.significant[1] != 1){
                chrB.start <- c(chrB.start,chrB.filt[2])
                chrB.end <- c(chrB.end,chrB.filt[length(chrB.filt)])
        }
}

if (!is.na(chrA.start) && !is.na(chrB.start)){
	area <- c()
	density <- c()
	chrA.start.pos <- c()
	chrA.end.pos <- c()
	chrB.start.pos <- c()
	chrB.end.pos <- c()
	i <- 1
	while (i <= length(chrA.start)){
		x <- (sqrt((chrA.start[i]-chrA.end[i])^2))
		j <- 1
		while (j <= length(chrB.start)){
			y <- (sqrt((chrB.start[j]-chrB.end[j])^2))
			v <- x*y
			if (v != 0){
				sum <- sum(u.mat[chrA.start[i]:chrA.end[i],chrB.start[j]:chrB.end[j]])
				area <- c(area,v)
				density <- c(density,(sum/v))
				chrA.start.pos <- c(chrA.start.pos,chrA.start[i])
				chrA.end.pos <- c(chrA.end.pos,chrA.end[i])
				chrB.start.pos <- c(chrB.start.pos,chrB.start[j])
	                	chrB.end.pos <- c(chrB.end.pos,chrB.end[j])
			}
			j <- j+1
		}
		i <- i+1
	}
	result <- data.frame(chrA.start.pos,chrA.end.pos,chrB.start.pos,chrB.end.pos,area,density)
	chrA.start <- c()
	chrA.end <- c()
	chrB.start <- c()
	chrB.end <- c()
	chrA.start.log.ratio <- c()
	chrA.end.log.ratio <- c()
	chrB.start.log.ratio <- c()
	chrB.end.log.ratio <- c()
	i <- 1
	while (i <= length(chrA.start.pos)){
		chrA.mid <- as.integer((chrA.start.pos[i]+chrA.end.pos[i])/2)-chrA.start.pos[i]
		chrA.half.start.up <- chrA.start.pos[i]-chrA.mid
		chrA.half.start.dw <- chrA.start.pos[i]+chrA.mid
		chrA.half.end.up <- chrA.end.pos[i]-chrA.mid
	        chrA.half.end.dw <- chrA.end.pos[i]+chrA.mid
		chrB.mid <- as.integer((chrB.start.pos[i]+chrB.end.pos[i])/2)-chrB.start.pos[i]
	        chrB.half.start.up <- chrB.start.pos[i]-chrB.mid
	        chrB.half.start.dw <- chrB.start.pos[i]+chrB.mid
	        chrB.half.end.up <- chrB.end.pos[i]-chrB.mid
	        chrB.half.end.dw <- chrB.end.pos[i]+chrB.mid
	
		if (chrA.half.start.up <= 0){chrA.half.start.up <- 1}
		if (chrA.half.start.dw <= 0){chrA.half.start.dw <- 1}
		if (chrA.half.end.up <= 0){chrA.half.end.up <- 1}
	        if (chrA.half.end.dw <= 0){chrA.half.end.dw <- 1}
		if (chrB.half.start.up <= 0){chrB.half.start.up <- 1}
	        if (chrB.half.start.dw <= 0){chrB.half.start.dw <- 1}
	        if (chrB.half.end.up <= 0){chrB.half.end.up <- 1}
	        if (chrB.half.end.dw <= 0){chrB.half.end.dw <- 1}
	
		if (chrA.half.start.up >= dim(u.mat)[[1]]){chrA.half.start.up <- dim(u.mat)[[1]]}
	        if (chrA.half.start.dw >= dim(u.mat)[[1]]){chrA.half.start.dw <- dim(u.mat)[[1]]}
	        if (chrA.half.end.up >= dim(u.mat)[[1]]){chrA.half.end.up <- dim(u.mat)[[1]]}
	        if (chrA.half.end.dw >= dim(u.mat)[[1]]){chrA.half.end.dw <- dim(u.mat)[[1]]}
	        if (chrB.half.start.up >= dim(u.mat)[[2]]){chrB.half.start.up <- dim(u.mat)[[2]]}
	        if (chrB.half.start.dw >= dim(u.mat)[[2]]){chrB.half.start.dw <- dim(u.mat)[[2]]}
	        if (chrB.half.end.up >= dim(u.mat)[[2]]){chrB.half.end.up <- dim(u.mat)[[2]]}
	        if (chrB.half.end.dw >= dim(u.mat)[[2]]){chrB.half.end.dw <- dim(u.mat)[[2]]}
		
		chrA.start.up.sum <- sum(u.mat[chrA.half.start.up:chrA.start.pos[i],chrB.start.pos[i]:chrB.end.pos[i]])
		chrA.start.dw.sum <- sum(u.mat[chrA.start.pos[i]:chrA.half.start.dw,chrB.start.pos[i]:chrB.end.pos[i]])
		chrA.start.log.ratio <- c(chrA.start.log.ratio,(log2(chrA.start.up.sum/chrA.start.dw.sum)))
		chrA.end.up.sum <- sum(u.mat[chrA.half.end.up:chrA.end.pos[i],chrB.start.pos[i]:chrB.end.pos[i]])
	       	chrA.end.dw.sum <- sum(u.mat[chrA.end.pos[i]:chrA.half.end.dw,chrB.start.pos[i]:chrB.end.pos[i]])
	       	chrA.end.log.ratio <- c(chrA.end.log.ratio,(log2(chrA.end.up.sum/chrA.end.dw.sum)))	
		chrB.start.up.sum <- sum(u.mat[chrA.start.pos[i]:chrA.end.pos[i],chrB.half.start.up:chrB.start.pos[i]])
	       	chrB.start.dw.sum <- sum(u.mat[chrA.start.pos[i]:chrA.end.pos[i],chrB.start.pos[i]:chrB.half.start.dw])
	       	chrB.start.log.ratio <- c(chrB.start.log.ratio,(log2(chrB.start.up.sum/chrB.start.dw.sum)))
	       	chrB.end.up.sum <- sum(u.mat[chrA.start.pos[i]:chrA.end.pos[i],chrB.half.end.up:chrB.end.pos[i]])
	       	chrB.end.dw.sum <- sum(u.mat[chrA.start.pos[i]:chrA.end.pos[i],chrB.end.pos[i]:chrB.half.end.dw])
	       	chrB.end.log.ratio <- c(chrB.end.log.ratio,(log2(chrB.end.up.sum/chrB.end.dw.sum)))
	
		i <- i+1
	}
		
	score <- chrA.start.log.ratio*chrA.end.log.ratio*chrB.start.log.ratio*chrB.end.log.ratio
	result.final <- data.frame(chrA.start.pos,chrA.end.pos,chrB.start.pos,chrB.end.pos,area,density,chrA.start.log.ratio,chrA.end.log.ratio,chrB.start.log.ratio,chrB.end.log.ratio,score)
	data <- paste(args[2],"/","data.txt",sep="")
	write.table(subset(result.final,((abs(score) > 0) | (abs(score) == 0 & (chrA.start.pos <= 10 | chrB.start.pos <= 10 | chrA.end.pos > (dim(u.mat)[[1]]-10) | chrB.end.pos > (dim(u.mat)[[2]]-10))) & area > 5)),file=data)
}else{
	print ("No Translocation Detected")
}
