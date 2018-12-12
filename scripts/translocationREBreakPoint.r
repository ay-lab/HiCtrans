library(data.table)
library(hashmap)
library(caTools)
library(depmixS4)
library(DEoptimR)
library(optparse)
#########################################
#buildmatrix program is taken from 
#HiC-Pro package

buildMatrixPath = "../scripts/buildmatrix"
#########################################

#Process the valid pair file 
mapVPs_on_REs <- function(re.bed,chrom.size,vp.file,prefix,chromA,chromB,id){
 
 system(paste0("grep -w -e ",chromA," -e ",chromB," ",re.bed," > ",prefix,"_",chromA,"_",chromB,"_REfrags.",id,".bed"),wait=T)
 system(paste0("grep -w -e ",chromA," -e ",chromB," ",chrom.size," > ",prefix,"_",chromA,"_",chromB,"_chromSizes.",id,".txt"),wait=T)
 cmd = paste0("./",buildMatrixPath," --binfile ",prefix,"_",chromA,"_",chromB,"_REfrags.",id,".bed --chrsizes ",prefix,"_",chromA,"_",chromB,"_chromSizes.",id,".txt --ifile ",vp.file," --oprefix ",prefix,"_",chromA,"_",chromB,"_",id," --matrix-format complete --detail-progress")
 print (cmd)
 system(cmd,wait=T)
 bed = as.data.frame(fread(as.character(paste0(prefix,"_",chromA,"_",chromB,"_",id,"_abs.bed")),header=F))
 mat = as.data.frame(fread(as.character(paste0(prefix,"_",chromA,"_",chromB,"_",id,".matrix")),header=F))
 colnames(bed) = c("chrom","start","end","index")
 colnames(mat) = c("indexA","indexB","count")
 bed_chrom.hash = hashmap(bed$index,as.character(bed$chrom))
 bed_start.hash = hashmap(bed$index,bed$start)
 bed_end.hash   = hashmap(bed$index,bed$end)
 mat.df = data.frame(
  chromA = bed_chrom.hash[[mat$indexA]],
  startA = bed_start.hash[[mat$indexA]],
  endA   = bed_end.hash[[mat$indexA]],
  chromB = bed_chrom.hash[[mat$indexB]],
  startB = bed_start.hash[[mat$indexB]],
  endB   = bed_end.hash[[mat$indexB]],
  count  = mat$count
 )
 write.table(mat.df,file=paste0(prefix,"_",chromA,"_",chromB,"_",id,".txt"),row.names=F,quote=F,sep="\t")
 return(list(mat.df=mat.df,bed=bed,chromA=chromA,chromB=chromB))
 
}

#Calculate the Cis and trans coverage values
cis_trans_coverage <- function(files,prefix,chrA,chrB,id){

 df = files$mat.df
 df_cis_chromA = df[with(df,chromA==as.character(chrA) & chromB==as.character(chrA)),]
 df_cis_chromB = df[with(df,chromA==as.character(chrB) & chromB==as.character(chrB)),]
 df_trans_chromA = df[with(df,chromA==as.character(chrA) & chromB==as.character(chrB)),]
 df_trans_chromB = df[with(df,chromA==as.character(chrB) & chromB==as.character(chrA)),]
 df_cis_chromA.up = df_cis_chromA[with(df_cis_chromA,startA > endB),]
 df_cis_chromB.up = df_cis_chromB[with(df_cis_chromB,startA > endB),]
 df_cis_chromA.dw = df_cis_chromA[with(df_cis_chromA,endA < startB),]
 df_cis_chromB.dw = df_cis_chromB[with(df_cis_chromB,endA < startB),]

 cis_chromA.up = aggregate(with(df_cis_chromA.up,count ~ startA),FUN = sum)
 cis_chromA.dw = aggregate(with(df_cis_chromA.dw,count ~ startA),FUN = sum)
 cis_chromB.up = aggregate(with(df_cis_chromB.up,count ~ startA),FUN = sum)
 cis_chromB.dw = aggregate(with(df_cis_chromB.dw,count ~ startA),FUN = sum)
 trans_chromA  = aggregate(with(df_trans_chromA,count ~ startA),FUN = sum)
 trans_chromB  = aggregate(with(df_trans_chromB,count ~ startA),FUN = sum)

 cis_chromA.up.hash = hashmap(cis_chromA.up$startA,cis_chromA.up$count)
 cis_chromB.up.hash = hashmap(cis_chromB.up$startA,cis_chromB.up$count)
 cis_chromA.dw.hash = hashmap(cis_chromA.dw$startA,cis_chromA.dw$count)
 cis_chromB.dw.hash = hashmap(cis_chromB.dw$startA,cis_chromB.dw$count)
 trans_chromA.hash  = hashmap(trans_chromA$startA,trans_chromA$count)
 trans_chromB.hash  = hashmap(trans_chromB$startA,trans_chromB$count)

 df = files$bed
 df_chromA = df[with(df,chrom==as.character(chrA)),]
 df_chromB = df[with(df,chrom==as.character(chrB)),]
 df_chromA[,"cis.up"] = cis_chromA.up.hash[[df_chromA$start]]
 df_chromA[,"cis.dw"] = cis_chromA.dw.hash[[df_chromA$start]]
 df_chromB[,"cis.up"] = cis_chromB.up.hash[[df_chromB$start]]
 df_chromB[,"cis.dw"] = cis_chromB.dw.hash[[df_chromB$start]]
 df_chromA[,"trans"]  = trans_chromA.hash[[df_chromA$start]]
 df_chromB[,"trans"]  = trans_chromB.hash[[df_chromB$start]]

 df_chromA$cis.up[is.na(df_chromA$cis.up)] = 0
 df_chromA$cis.dw[is.na(df_chromA$cis.dw)] = 0
 df_chromB$cis.up[is.na(df_chromB$cis.up)] = 0
 df_chromB$cis.dw[is.na(df_chromB$cis.dw)] = 0
 df_chromA$trans[is.na(df_chromA$trans)] = 0
 df_chromB$trans[is.na(df_chromB$trans)] = 0

 df = rbind(df_chromA,df_chromB)
 write.table(df,file=paste0(prefix,"_",chrA,"_",chrB,"_",id,"_Coverage.bed"),col.names=F,row.names=F,quote=F,sep="\t")
}


#Calculate cis and trans coverage and then calculate directionality index. Followed by HMM segmentation (2 states)
di_and_hmm <- function(coverage,chromA,chromA.start,chromA.end,chromB,chromB.start,chromB.end,prefix,id){
 
 bed = as.data.frame(fread(as.character(coverage),header=F))
 colnames(bed) = c("chrom","start","end","index","up","dw","trans")
 chromA_span = chromA.end-chromA.start
 chromB_span = chromB.end-chromB.start
 if (chromA_span < 1e6){
  chromA_span  = 1e6-chromA_span
  chromA_span  = floor(chromA_span/2)
  chromA.start = chromA.start-chromA_span
  chromA.end   = chromA.end+chromA_span 
 }
 if (chromB_span < 1e6){
  chromB_span  = 1e6-chromB_span
  chromB_span  = floor(chromB_span/2)+1
  chromB.start = chromB.start-chromB_span
  chromB.end   = chromB.end+chromB_span
 }

 bed_chromA = bed[with(bed,chrom==chromA),]
 bed_chromB = bed[with(bed,chrom==chromB),]
 bed_chromA_cis_di = list()
 bed_chromA_trans_di = list()
 i = 1
 while (i <= nrow(bed_chromA)){
  up = bed_chromA$up[i]
  dw = bed_chromA$dw[i]
  trans = bed_chromA$trans[i]
  bed_chromA_cis_di[[i]]   = (dw-up)/abs(dw-up)*(((up-mean(c(dw,up)))^2/mean(c(dw,up)))+((dw-mean(c(dw,up)))^2/mean(c(dw,up))))
  bed_chromA_trans_di[[i]] = ((trans-mean(bed_chromA$trans))/abs(trans-mean(bed_chromA$trans)))*((trans-mean(bed_chromA$trans))^2/mean(bed_chromA$trans))
  i = i+1
 }
 
 smooth = 5
 bed_chromA_cis_di = unlist(bed_chromA_cis_di)
 bed_chromA_trans_di = unlist(bed_chromA_trans_di)
 bed_chromA_cis_di[is.na(bed_chromA_cis_di)] = 0
 bed_chromA_trans_di[is.na(bed_chromA_trans_di)] = 0
 bed_chromA_cis_di = runmean(bed_chromA_cis_di,k=smooth,endrule="mean")
 bed_chromA_trans_di = runmean(bed_chromA_trans_di,k=smooth,endrule="mean")
 bed_chromA = cbind(bed_chromA,bed_chromA_cis_di,bed_chromA_trans_di)
 bed_chromA = bed_chromA[with(bed_chromA,start >= chromA.start & end <= chromA.end),]

 bed_chromB_cis_di = list()
 bed_chromB_trans_di = list()
 i = 1
 while (i <= nrow(bed_chromB)){
  up = bed_chromB$up[i]
  dw = bed_chromB$dw[i]
  trans = bed_chromB$trans[i]
  bed_chromB_cis_di[[i]]   = (dw-up)/abs(dw-up)*(((up-mean(c(dw,up)))^2/mean(c(dw,up)))+((dw-mean(c(dw,up)))^2/mean(c(dw,up))))
  bed_chromB_trans_di[[i]] = ((trans-mean(bed_chromB$trans))/abs(trans-mean(bed_chromB$trans)))*((trans-mean(bed_chromB$trans))^2/mean(bed_chromB$trans))
  i = i+1
 }
 
 bed_chromB_cis_di = unlist(bed_chromB_cis_di)
 bed_chromB_trans_di = unlist(bed_chromB_trans_di)
 bed_chromB_cis_di[is.na(bed_chromB_cis_di)] = 0
 bed_chromB_trans_di[is.na(bed_chromB_trans_di)] = 0
 bed_chromB_cis_di = runmean(bed_chromB_cis_di,k=smooth,endrule="mean")
 bed_chromB_trans_di = runmean(bed_chromB_trans_di,k=smooth,endrule="mean")
 bed_chromB = cbind(bed_chromB,bed_chromB_cis_di,bed_chromB_trans_di)
 bed_chromB = bed_chromB[with(bed_chromB,start >= chromB.start & end <= chromB.end),]

 nstates = 2
 print (head(bed_chromA)) 
 chromA.hmm = depmix(response=list(bed_chromA_cis_di~1,bed_chromA_trans_di~1),data=bed_chromA, nstates=nstates,family=list(gaussian(),gaussian()),ntimes=nrow(bed_chromA))
 chromA.hmm.fit = fit(chromA.hmm)
 chromA.hmm.states = posterior(chromA.hmm.fit)$state
 x = summary(chromA.hmm.fit)
 chromA.trans.state = which.max(as.vector(x[1:nstates,3]))
 cis.hash   = hashmap(c(1:nstates),as.vector(x[1:nstates,1]))
 trans.hash = hashmap(c(1:nstates),as.vector(x[1:nstates,3])) 
 chromA_cis_DI = as.vector(x[chromA.trans.state,1])
 bed_chromA = cbind(bed_chromA,chromA.hmm.states,cis.mean=cis.hash[[chromA.hmm.states]],trans.mean=trans.hash[[chromA.hmm.states]])
 chrA.hmmBoundary.start = min(bed_chromA[bed_chromA$chromA.hmm.states==chromA.trans.state,]$start)
 chrA.hmmBoundary.end   = max(bed_chromA[bed_chromA$chromA.hmm.states==chromA.trans.state,]$end)
 cat (chrA.hmmBoundary.start,"\t",chrA.hmmBoundary.end,"\n")
 write.table(bed_chromA,file=paste0(prefix,".",chromA,"_",id,".hmm.states.bed"),col.names=T,row.names=F,sep="\t",quote=F)

 print (head(bed_chromB))
 chromB.hmm = depmix(response=list(bed_chromB_cis_di~1,bed_chromB_trans_di~1),data=bed_chromB, nstates=nstates,family=list(gaussian(),gaussian()),ntimes=nrow(bed_chromB))
 chromB.hmm.fit = fit(chromB.hmm)
 chromB.hmm.states = posterior(chromB.hmm.fit)$state
 x = summary(chromB.hmm.fit)
 chromB.trans.state = which.max(as.vector(x[1:nstates,3]))
 cis.hash   = hashmap(c(1:nstates),as.vector(x[1:nstates,1]))
 trans.hash = hashmap(c(1:nstates),as.vector(x[1:nstates,3]))
 chromB_cis_DI = as.vector(x[chromB.trans.state,1])
 bed_chromB = cbind(bed_chromB,chromB.hmm.states,cis.mean=cis.hash[[chromB.hmm.states]],trans.mean=trans.hash[[chromB.hmm.states]])
 chrB.hmmBoundary.start = min(bed_chromB[bed_chromB$chromB.hmm.states==chromB.trans.state,]$start)
 chrB.hmmBoundary.end   = max(bed_chromB[bed_chromB$chromB.hmm.states==chromB.trans.state,]$end)
 cat (chrB.hmmBoundary.start,"\t",chrB.hmmBoundary.end,"\n")
 write.table(bed_chromB,file=paste0(prefix,".",chromB,"_",id,".hmm.states.bed"),col.names=T,row.names=F,sep="\t",quote=F)
 return(
  list(
   chrA.hmmBoundary.start=chrA.hmmBoundary.start,
   chrA.hmmBoundary.end=chrA.hmmBoundary.end,
   chrB.hmmBoundary.start=chrB.hmmBoundary.start,
   chrB.hmmBoundary.end=chrB.hmmBoundary.end,
   chromA_cis_DI=chromA_cis_DI,
   chromB_cis_DI=chromB_cis_DI
  )
 )

}

#Optimize the HMM segment border (+/-bp) to find the breakpoint using DE optimization
optimization_func <- function(boundary,df,as,bs,ae,be,bk,chromA_cis_DI,chromB_cis_DI){
 
 chromA_span = ae-as
 chromB_span = be-bs
 chromA_pos  = boundary[1]
 chromB_pos  = boundary[2]
 chromA_up   = floor(chromA_pos-(chromA_span/2))
 chromA_dw   = floor(chromA_pos+(chromA_span/2))
 chromB_up   = floor(chromB_pos-(chromB_span/2))
 chromB_dw   = floor(chromB_pos+(chromB_span/2))
 chromA_dw.chromB_up = sum(df[df$startA >= chromA_pos & df$endA <= chromA_dw & df$startB >= chromB_up & df$endB <= chromB_pos,]$count)
 chromA_up.chromB_dw = sum(df[df$startA >= chromA_up & df$endA <= chromA_pos & df$startB >= chromB_pos & df$endB <= chromB_dw,]$count)
 chromAB_up  = sum(df[df$startA >= chromA_up & df$endA <= chromA_pos & df$startB >= chromB_up & df$endB <= chromB_pos,]$count)
 chromAB_dw  = sum(df[df$startA >= chromA_pos & df$endA <= chromA_dw & df$startB >= chromB_pos & df$endB <= chromB_dw,]$count)

 if ((chromA_cis_DI > 0 & chromB_cis_DI < 0) | (chromA_cis_DI < 0 & chromB_cis_DI > 0)){
  if (chromA_cis_DI > 0 & chromB_cis_DI < 0){
   if (bk == 1){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
     } else if (bk == 2){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
     } else if (bk == 3){
       value = (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
     } else if (bk == 4){
       value = (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
     }
  } else if (chromA_cis_DI < 0 & chromB_cis_DI > 0){
     if (bk == 1){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
     } else if (bk == 2){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
     } else if (bk == 3){
       value = (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
     } else if (bk == 4){
       value = (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
     }
  }
 } else if ((chromA_cis_DI > 0 & chromB_cis_DI > 0) | (chromA_cis_DI < 0 & chromB_cis_DI < 0)){
    if (chromA_cis_DI > 0 & chromB_cis_DI > 0){
     if (bk == 1){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
     } else if (bk == 2){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
     } else if (bk == 3){
       value = (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
     } else if (bk == 4){
       value = (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
     }
    } else if (chromA_cis_DI < 0 & chromB_cis_DI < 0){
       if (bk == 1){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
     } else if (bk == 2){
       value = (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
     } else if (bk == 3){
       value = (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
     } else if (bk == 4){
       value = (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
     }
    }
 }
 value
}

bp_optimization <- function(intr.file,boundary,chromA,chromB,ssA,seA,ssB,seB){

 intr.file = as.data.frame(fread(as.character(intr.file),header=T))
 intr.file = intr.file[intr.file$chromA == as.character(chromA) & intr.file$chromB == as.character(chromB),]
 print (head(intr.file))
 print (tail(intr.file))

 as=boundary$chrA.hmmBoundary.start-ssA
 bs=boundary$chrB.hmmBoundary.start-ssB
 ae=boundary$chrA.hmmBoundary.end+seA
 be=boundary$chrB.hmmBoundary.end+seB
 chromA_span = ae-as
 chromB_span = be-bs

 upup = sum(intr.file[with(intr.file,startA >= as & endA <= as+floor(chromA_span/2) & startB >= bs & endB <= bs+floor(chromB_span/2)),]$count)
 dwdw = sum(intr.file[with(intr.file,startA >= as+floor(chromA_span/2) & endA <= ae & startB >= bs+floor(chromB_span/2) & endB <= be),]$count)
 updw = sum(intr.file[with(intr.file,startA >= as & endA <= as+floor(chromA_span/2) & startB >= bs+floor(chromB_span/2) & endB <= be),]$count)
 dwup = sum(intr.file[with(intr.file,startA >= as+floor(chromA_span/2) & endA <= ae & startB >= bs & endB <= bs+floor(chromB_span/2)),]$count)
 mx = which.max(c(upup,dwdw,updw,dwup))
 print (mx)
 chromA_span = ae-as
 chromB_span = be-bs
 if (chromA_span < 1e6){
  chromA_span  = 1e6-chromA_span
  chromA_span  = floor(chromA_span/2)
  as = as-chromA_span
  ae = ae+chromA_span
 }
 if (chromB_span < 1e6){
  chromB_span  = 1e6-chromB_span
  chromB_span  = floor(chromB_span/2)+1
  bs = bs-chromB_span
  be = be+chromB_span
 }

 #Optimization parameters		
 JDEoptim(
  lower=c(as,bs), 
  upper=c(ae,be), 
  fn=optimization_func,
  df=intr.file,
  compare_to="median",
  maxiter=1000,
  NP=100,
  fnscale=1e-5,
  tau_pF=0.1,
  as=as,
  bs=bs,
  ae=ae,
  be=be,
  bk=mx,
  chromA_cis_DI = as.numeric(boundary$chromA_cis_DI),
  chromB_cis_DI = as.numeric(boundary$chromB_cis_DI),
  trace=TRUE
 ) 
}

#Get options here
option_list = list(
 make_option(c("--fragsFile"), type="character", help="Restriction Fragment file\nFormat:
 chr1    0       16007   HIC_chr1_1      0       +
 chr1    16007   24571   HIC_chr1_2      0       +
 chr1    24571   27981   HIC_chr1_3      0       +
 ......"),
 make_option(c("--chromsize"), type="character", help="Chromosome size file\nFormat:
 chr1\t249250621
 chr2\t243199373
 chr3\t198022430
 chr4\t191154276
 ....."),
 make_option(c("--validpair"), type="character", help="Valid pair file of the HiC data\nFormat:
 SRR6213722.1\tchr11\t124331538\t-\tchr11\t124345246\t-
 SRR6213722.2\tchr1\t198436365\t-\tchr1\t199923196\t+
 ....."),      
 make_option(c("--prefix"),    type="character", help="Output files will be saved by this name"),
 make_option(c("--hictrans"),  default=NA,type="character", help="HiCtrans output file i.e. *.Translocation.EntropyFiltered.result file.\n\nIf not given, provide the chromA/startA/endA & chromB/startB/endB regions"),
 make_option(c("--chromA"),    default=NA,type="character", help="Chromosome A to be processed"),
 make_option(c("--startA"),    default=NA,type="integer", help="Chromosome A start position"),
 make_option(c("--endA"),      default=NA,type="integer", help="Chromosome A end position"),
 make_option(c("--chromB"),    default=NA,type="character", help="Chromosome B to be processed"),
 make_option(c("--startB"),    default=NA,type="integer", help="Chromosome B start position"),
 make_option(c("--endB"),      default=NA,type="integer", help="Chromosome B end position\n\nFine tuning parameter:\n"),
 make_option(c("--ssA"),      default=1e5, type="integer", help="Extend -(ve) bp of the 5' HMM segment border of chromosome A for breakpoint identification. Default 100Kb."),
 make_option(c("--seA"),      default=1e5, type="integer", help="Extend +(ve) bp of the 3' HMM segment border of chromosome A for breakpoint identification. Default 100Kb."),
 make_option(c("--ssB"),      default=1e5, type="integer", help="Extend -(ve) bp of the 5' HMM segment border of chromosome B for breakpoint identification. Default 100Kb."),
 make_option(c("--seB"),      default=1e5, type="integer", help="Extend +(ve) bp of the 3' HMM segment border of chromosome B for breakpoint identification. Default 100Kb.")
)

opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$hictrans)){
 if (is.na(opt$chromA) & is.na(opt$startA) & is.na(opt$endA) & is.na(opt$chromB) & is.na(opt$startB) & is.na(opt$endB)){
  stop("Please provide either HiCtrans output file or chromA/startA/endA & chromB/startB/endB regions\n")
 } else {
   df = mapVPs_on_REs(
        re.bed=as.character(opt$fragsFile),
        chrom.size=as.character(opt$chromsize),
        vp.file=as.character(opt$validpair),
        prefix=as.character(opt$prefix),
        chromA=as.character(opt$chromA),
        chromB=as.character(opt$chromB),
	id=1
   )
   cis_trans_coverage(df,as.character(opt$prefix),as.character(opt$chromA),as.character(opt$chromB),id=1)
   boundary = di_and_hmm(
        coverage=paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,"_1_Coverage.bed"),
        chromA=as.character(opt$chromA),
        chromA.start=as.integer(opt$startA),
        chromA.end=as.integer(opt$endA),
        chromB=as.character(opt$chromB),
        chromB.start=as.integer(opt$startB),
        chromB.end=as.integer(opt$endB),
        prefix=opt$prefix,
	id=1
   ) 
   bp_optimized = bp_optimization(
	paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,"_1.txt"),
	boundary,
	as.character(opt$chromA),
	as.character(opt$chromB),
	ssA=as.integer(opt$ssA),
	ssB=as.integer(opt$ssB),
	seA=as.integer(opt$seA),
	seB=as.integer(opt$seB)
	
   )
   chromA_BP_start=head(df$bed[df$bed$chrom==as.character(opt$chromA) & df$bed$start >= floor(bp_optimized$par[1]),]$start,1)
   chromA_BP_end=head(df$bed[df$bed$chrom==as.character(opt$chromA) & df$bed$start >= floor(bp_optimized$par[1]),]$end,1)
   chromB_BP_start=head(df$bed[df$bed$chrom==as.character(opt$chromB) & df$bed$start >= floor(bp_optimized$par[2]),]$start,1)
   chromB_BP_end=head(df$bed[df$bed$chrom==as.character(opt$chromB) & df$bed$start >= floor(bp_optimized$par[2]),]$end,1)
   cat (chromA_BP_start,"\t",chromA_BP_end,"\t",chromB_BP_start,"\t",chromB_BP_end,"\n")	
   param = cbind(
	chrA=as.character(opt$chromA),
	BoundaryAS=as.integer(opt$startA),
 	BoundaryAE=as.integer(opt$endA),
	chrB=as.character(opt$chromB),
	BoundaryBS=as.integer(opt$startB),
        BoundaryBE=as.integer(opt$endB),
	chromA_BP_start=chromA_BP_start,
	chromA_BP_end=chromA_BP_end,
	chromB_BP_start=chromB_BP_start,
	chromB_BP_end=chromB_BP_end
   )
   write.table(param,file=paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,".BreakPoints.txt"),col.names=T,row.names=F,sep="\t",quote=F)
  }
 } else {
   param = read.table(opt$hictrans,header=T)
   chromA_BP_start = list()
   chromA_BP_end=list()
   chromB_BP_start = list()
   chromB_BP_end=list()
   k = 1
   while (k <= nrow(param)){
    opt$chromA = param$chrA[k]
    opt$startA = param$BoundaryAS[k]
    opt$endA   = param$BoundaryAE[k]
    opt$chromB = param$chrB[k]
    opt$startB = param$BoundaryBS[k]
    opt$endB   = param$BoundaryBE[k]

    df  = mapVPs_on_REs(
  	re.bed=as.character(opt$fragsFile),
   	chrom.size=as.character(opt$chromsize),
  	vp.file=as.character(opt$validpair),
  	prefix=as.character(opt$prefix),
   	chromA=as.character(opt$chromA),
 	chromB=as.character(opt$chromB),
 	id=k
    )

    cis_trans_coverage(df,as.character(opt$prefix),as.character(opt$chromA),as.character(opt$chromB),id=k)
    boundary = di_and_hmm(
 	coverage=paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,"_",k,"_Coverage.bed"),
  	chromA=as.character(opt$chromA),
    	chromA.start=as.integer(opt$startA),
 	chromA.end=as.integer(opt$endA),
  	chromB=as.character(opt$chromB),
  	chromB.start=as.integer(opt$startB),
 	chromB.end=as.integer(opt$endB),
       	prefix=opt$prefix,
	id=k
    )
    bp_optimized = bp_optimization(
  	paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,"_",k,".txt"),
   	boundary,
  	as.character(opt$chromA),
  	as.character(opt$chromB),
    	ssA=as.integer(opt$ssA),
	ssB=as.integer(opt$ssB),
  	seA=as.integer(opt$seA),
  	seB=as.integer(opt$seB)
   )
   chromA_BP_start[k]=head(df$bed[df$bed$chrom==as.character(opt$chromA) & df$bed$start >= floor(bp_optimized$par[1]),]$start,1)
   chromA_BP_end[k]=head(df$bed[df$bed$chrom==as.character(opt$chromA) & df$bed$start >= floor(bp_optimized$par[1]),]$end,1)
   chromB_BP_start[k]=head(df$bed[df$bed$chrom==as.character(opt$chromB) & df$bed$start >= floor(bp_optimized$par[2]),]$start,1)
   chromB_BP_end[k]=head(df$bed[df$bed$chrom==as.character(opt$chromB) & df$bed$start >= floor(bp_optimized$par[2]),]$end,1)
   k = k+1
  }
  chromA_BP_start=unlist(chromA_BP_start)
  chromA_BP_end=unlist(chromA_BP_end)
  chromB_BP_start=unlist(chromB_BP_start)
  chromB_BP_end=unlist(chromB_BP_end)
  param = cbind(param,chromA_BP_start,chromA_BP_end,chromB_BP_start,chromB_BP_end)
  write.table(param,file=paste0(opt$prefix,"_",opt$chromA,"_",opt$chromB,".BreakPoints.txt"),col.names=T,row.names=F,sep="\t",quote=F)
}

