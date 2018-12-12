##
## This script will calculate the extropy of count distribution inside the translocation
## boxes predicted by HiCtrans (*.Translocation.results) as compared to random boxes.
## If there is a true translocation than count distribution will be more random and should
## have high normalized entropy value as compared to similar sized random boxes.
##

library(optparse)

option_list = list(
        make_option(c("-f", "--fragsFile"),  type="character", help="Fragment file"),
        make_option(c("-t", "--translocationFile"),   type="character", help="Translocation result file"),
	make_option(c("-s", "--sampleNumber"),   type="character", help="Number of random locations from each chromosomes should be compared to generate random entropy upper CI. "),
	make_option(c("-c", "--countValue"),   type="character", help="Process a translocation above this count value."),
	make_option(c("-n", "--noiseRatio"),   type="character", help="Threshold ratio of Box enrtopy by 99% upper confidence interval of random box entropy."),
	make_option(c("-o", "--outfile"),   type="character", help="Output file with translocation box entropy informtation")
)
opt <- parse_args(OptionParser(option_list=option_list))

frags = read.table(as.character(opt$f),header=F)
colnames(frags) = c("chrA","chrA_Mid","chrB","chrB_Mid", "count")

result = read.table(as.character(opt$t),header=T)
###### Abhijit: 12/10/2018 ######
#colnames(result) = c("chrA","BoundaryAS","BreakPointA","BoundaryAE","chrB","BoundaryBS","BreakPointB","BoundaryBE","count")
colnames(result) = c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE","count")
###################

func <- function(as,ae,bs,be,amax,amin,bmax,bmin,count){
	aspan = ae-as
	bspan = be-bs
	arand = as.integer(runif(as.integer(opt$s), min(amin), (max(amax)-aspan)))
	i = 1
	while (i <= nrow(result)){
		span  = result$BoundaryAE[i]-result$BoundaryAS[i]
		span  = span/count
		arand = setdiff(arand,seq((result$BoundaryAS[i]-span):(result$BoundaryAE[i]+span)))
		i = i+1
	}
	brand = as.integer(runif(as.integer(opt$s), min(bmin), (max(bmax)-bspan)))
	i = 1
        while (i <= nrow(result)){
                span  = result$BoundaryBE[i]-result$BoundaryBS[i]
		span  = span/count
                brand = setdiff(brand,seq((result$BoundaryBS[i]-span):(result$BoundaryBE[i]+span)))
                i = i+1
        }
	combination = expand.grid(arand,brand)
	colnames(combination) = c("arand","brand")
	combination = as.data.frame(combination)
	entropy = c()
	i = 1
	while (i <= nrow(combination)){
		astart = combination$arand[i]	
		bstart = combination$brand[i]
		aend = astart + aspan
		bend = bstart + bspan
		frags.filter.tmp = frags[frags$chrA_Mid >= astart & frags$chrA_Mid <= aend & frags$chrB_Mid >= bstart & frags$chrB_Mid <= bend,]$count
		if (length(frags.filter.tmp) < 2){
			entropy = c(entropy,0)
		} else {
			p <- table(frags.filter.tmp)
		        p <- p/sum(p)
        		e <- (sum(-p*log(p)))/log(length(frags.filter.tmp))
			entropy = c(entropy,e)

		}
		i = i+1
	}
	entropy.sd = sd(entropy)
	entropy.len = length(entropy)
	error = qnorm(0.99)*entropy.sd/sqrt(entropy.len)
	uCI = mean(entropy)+error
	if (is.na(uCI)){
		count = count+1
		func(as=as,ae=ae,bs=bs,be=be,amax=amax,amin=amin,bmax=bmax,bmin=bmin,count=count)
	} else {
		return(uCI)
	}
}

if (nrow(result) > 0){
if (length(result[result$count >= as.integer(opt$c),]$count) > 0){
	y = c()
	i = 1
	while (i <= nrow(result)){
		count = 1;
		frags.filter = frags[frags$chrA_Mid >= result$BoundaryAS[i] & frags$chrA_Mid <= result$BoundaryAE[i] & frags$chrB_Mid >= result$BoundaryBS[i] & frags$chrB_Mid <= result$BoundaryBE[i],]
		p = table(frags.filter$count)
		p = p/sum(p)
		box.entropy = (sum(-p*log(p)))/log(length(frags.filter$count))
		random.entropy.99uCI = func(as=result$BoundaryAS[i],ae=result$BoundaryAE[i],bs=result$BoundaryBS[i],be=result$BoundaryBE[i],amax=max(frags$chrA_Mid),amin=min(frags$chrA_Mid),bmax=max(frags$chrB_Mid),bmin=min(frags$chrB_Mid),count=count)
		ratio = box.entropy/random.entropy.99uCI
		x <- cbind(result[i,],box.entropy,random.entropy.99uCI,ratio)
		y <- rbind(y,x)
		i = i+1
	}
	
	y = y[y$count >= as.integer(opt$c),]
	y = y[y$ratio >= as.integer(opt$n),]
	if (nrow(y) > 0){
                colnames(y) = c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE","MaxCount","box.entropy","random.entropy.99uCI","ratio")
		write.table(y,file=as.character(opt$o),row.names=F,quote=F,sep="\t")
	} else {
		cat ("No translocation\n")
	}
} else {
	cat ("No translocation\n")	
}
} else {
        cat ("No translocation\n")
}
