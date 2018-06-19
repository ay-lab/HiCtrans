library(data.table)
library(hashmap)
library(optparse)

option_list = list(
        make_option("--bedtools",  type="character", help="Path to bedtools program"),
        make_option("--genomesize",   type="character", help="chromosome wise genome size file. eg. a file like hg19.txt with following information
	chr1    249250621
	chr2    243199373
	chr3    198022430
	chr4    191154276
	chr5    180915260
	chr6    171115067
	chr7    159138663
	chr8    146364022
	chr9    141213431
	chr10   135534747
	chr11   135006516
	chr12   133851895
	chr13   115169878
	chr14   107349540
	chr15   102531392
	chr16   90354753
	chr17   81195210
	chr18   78077248
	chr19   59128983
	chr20   63025520
	chr21   48129895
	chr22   51304566
	chrX    155270560
	chrY    59373566
	chrM    16571
	"),
	make_option("--resolution",  type="integer", help="Hi-C resolution"),
	make_option("--column1",  type="character", help="Column1 chromosome name of juicebox output"),
	make_option("--column2",  type="character", help="Column2 chromosome name of juicebox output"),
	make_option("--juicebox_file",  type="character", help="Juicebox output file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

bedtools    = opt$bedtools
genomeSize  = opt$genomesize
window	    = opt$resolution
chromA_name = opt$column1
chromB_name = opt$column2
juicebox_input = opt$juicebox_file

cmd = paste0(bedtools," makewindows -g ",genomeSize," -w ",window," |awk '{c++; print $0,\"\\t\"c}' > index_",window,".bed")
print (cmd)
system(cmd,wait=T)
window_df = fread(paste0("index_",window,".bed"),header=F)
colnames(window_df) = c("chrom","start","end","index")
head(window_df)
write.table(window_df,file="index.bed",quote=F,col.names=F,row.names=F,sep="\t")

index_chromA = hashmap(keys=window_df[window_df$chrom==as.character(chromA_name),]$start,values=as.character(window_df[window_df$chrom==as.character(chromA_name),]$index))
index_chromB = hashmap(keys=window_df[window_df$chrom==as.character(chromB_name),]$start,values=as.character(window_df[window_df$chrom==as.character(chromB_name),]$index))

juicebox_input = fread(as.character(juicebox_input),header=F)
colnames(juicebox_input) = c("indexA","indexB","count")
head(juicebox_input)

hictrans_matrix = data.frame(
	indexA = index_chromA[[juicebox_input$indexA]],
	indexB = index_chromB[[juicebox_input$indexB]],
	count  = juicebox_input$count
)
head(hictrans_matrix)
window_df = window_df[window_df$chrom==as.character(chromA_name) | window_df$chrom==as.character(chromB_name),]
write.table(window_df,file=paste0(chromA_name,"_",chromB_name,".hictrans.index.bed"),quote=F,sep="\t",row.names=F,col.names=F)
write.table(hictrans_matrix,file=paste0(chromA_name,"_",chromB_name,".hictrans.matrix"),quote=F,sep="\t",row.names=F,col.names=F)
