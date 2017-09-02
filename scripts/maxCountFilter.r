args = commandArgs(trailingOnly=TRUE)
frags = read.table(args[1],header=F)
inter = read.table(args[2],header=F)
chromA = args[3]
chromB = args[4]
chromA_total = length(subset(frags,V1==chromA)$V1)
chromB_total = length(subset(frags,V1==chromB)$V1)

counts = inter$V5
inter.pval = c()
n = chromA_total*chromB_total

i = 1
while (i <= length(counts)){
	inter.pval = c(inter.pval,binom.test(counts[i],n,p=(1/n))$p.value)
	i = i+1
}

chromA_name  = inter$V1
chromA_pos   = inter$V2
chromB_name  = inter$V3
chromB_pos   = inter$V4
count_val    = inter$V5
data = data.frame(chromA_name,chromA_pos,chromB_name,chromB_pos,count_val,inter.pval)
data.significant = subset(data,inter.pval <= 0.01 & inter.pval!="NA")
if (length(data.significant$inter.pval) > 0){
	data.significant = data.significant[with(data.significant, order(inter.pval)), ]
}
write.table(data.significant,file=paste(args[5],".tmp.result",sep=""),row.names=F,col.names=F,quote=F)
