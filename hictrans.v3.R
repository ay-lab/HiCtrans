library(data.table)
library(hashmap)
suppressMessages(library(changepoint))
library(hashmap)
library(optparse)
library(Rcpp)
library(caTools)
library(depmixS4)
library(DEoptimR)


############### Get breakpoint calls at restriction level resolution #############

buildMatrixPath <- "buildmatrix"
buildMatrixPath <- normalizePath(buildMatrixPath) 
## Cluster the breakpoints
RE_HClust <- function(box.df,cl.A,cl.B) {

  box.df <- box.df[,c(1:6)]
  colnames(box.df) <- c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE")
  chrA.dist <- matrix(ncol=nrow(box.df),nrow=nrow(box.df))
  chrB.dist <- matrix(ncol=nrow(box.df),nrow=nrow(box.df))

  i <- 1
  while (i <= nrow(box.df)) {
    chrA.AS_i <- box.df$BoundaryAS[i]
    chrA.AE_i <- box.df$BoundaryAE[i]
    chrB.BS_i <- box.df$BoundaryBS[i]
    chrB.BE_i <- box.df$BoundaryBE[i]
    j <- 1
    while (j <= nrow(box.df)) {
      chrA.AS_j <- box.df$BoundaryAS[j]
      chrA.AE_j <- box.df$BoundaryAE[j]
      chrB.BS_j <- box.df$BoundaryBS[j]
      chrB.BE_j <- box.df$BoundaryBE[j]
      chrA.AS_AS.ij <- abs(chrA.AS_i - chrA.AS_j)
      chrA.AE_AS.ij <- abs(chrA.AE_i - chrA.AS_j)
      chrA.AS_AE.ij <- abs(chrA.AS_i - chrA.AE_j)
      chrA.AE_AE.ij <- abs(chrA.AE_i - chrA.AE_j)
      chrB.BS_BS.ij <- abs(chrB.BS_i - chrB.BS_j)
      chrB.BE_BS.ij <- abs(chrB.BE_i - chrB.BS_j)
      chrB.BS_BE.ij <- abs(chrB.BS_i - chrB.BE_j)
      chrB.BE_BE.ij <- abs(chrB.BE_i - chrB.BE_j)
      chrA.dist[i,j] <- min(chrA.AS_AS.ij,chrA.AE_AS.ij,chrA.AS_AE.ij,chrA.AE_AE.ij)
      chrB.dist[i,j] <- min(chrB.BS_BS.ij,chrB.BE_BS.ij,chrB.BS_BE.ij,chrB.BE_BE.ij)
      j <- j + 1
    }
    i <- i + 1
  }

  if (nrow(chrA.dist) > 1 & nrow(chrB.dist) > 1) {
    chrA.hclust <- hclust(as.dist(chrA.dist),method="single")
    chrB.hclust <- hclust(as.dist(chrB.dist),method="single")
    chrA.clust  <- cutree(chrA.hclust,h=cl.A)
    chrB.clust  <- cutree(chrB.hclust,h=cl.B)
    box.df[,"chrA.clus"] <- chrA.clust
    box.df[,"chrB.clus"] <- chrB.clust
  } else {
    chrA.clust <- 1
    chrB.clust <- 1
    box.df[,"chrA.clus"] <- chrA.clust
    box.df[,"chrB.clus"] <- chrB.clust
  }
  clus <- unique(data.frame(chrA.clust,chrB.clust))
  box.cl <- list()
  i <- 1
  while (i <= nrow(clus)) {
    d <- box.df[box.df$chrA.clus==clus$chrA.clust[i] & box.df$chrB.clus==clus$chrB.clust[i],]
    box.cl[[i]] <- data.frame(chrA=d$chrA[1],BoundaryAS=min(d$BoundaryAS),BoundaryAE=max(d$BoundaryAE),
                              chrB=d$chrB[1],BoundaryBS=min(d$BoundaryBS),BoundaryBE=max(d$BoundaryBE))
    i <- i + 1
  }
  box.cl <- as.data.frame(do.call(rbind,box.cl))
  colnames(box.cl) <- c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE")
  return(box.cl)
}

#Process the valid pair file
mapVPs_on_REs <- function(re.bed,chrom.size,vp.file,prefix,chromA,chromB,id){

  system(paste0("grep -w -e ",chromA," -e ",chromB," ",re.bed," > ",prefix,"_",chromA,"_",chromB,"_REfrags.",id,".bed"),wait=T)
  system(paste0("grep -w -e ",chromA," -e ",chromB," ",chrom.size," > ",prefix,"_",chromA,"_",chromB,"_chromSizes.",id,".txt"),wait=T)
  cmd <- paste0(buildMatrixPath," --binfile ",prefix,"_",chromA,"_",chromB,"_REfrags.",id,".bed --chrsizes ",prefix,"_",chromA,"_",chromB,"_chromSizes.",id,".txt --ifile ",vp.file," --oprefix ",prefix,"_",chromA,"_",chromB,"_",id," --matrix-format complete")
  system(cmd,wait=T)
  bed <- as.data.frame(fread(as.character(paste0(prefix,"_",chromA,"_",chromB,"_",id,"_abs.bed")),header=F))
  mat <- as.data.frame(fread(as.character(paste0(prefix,"_",chromA,"_",chromB,"_",id,".matrix")),header=F))
  colnames(bed) <- c("chrom","start","end","index")
  colnames(mat) <- c("indexA","indexB","count")
  bed_chrom.hash <- hashmap(bed$index,as.character(bed$chrom))
  bed_start.hash <- hashmap(bed$index,bed$start)
  bed_end.hash   <- hashmap(bed$index,bed$end)
  mat.df <- data.frame(
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

  df <- files$mat.df
  df_cis_chromA <- df[with(df,chromA==as.character(chrA) & chromB==as.character(chrA)),]
  df_cis_chromB <- df[with(df,chromA==as.character(chrB) & chromB==as.character(chrB)),]
  df_trans_chromA <- df[with(df,chromA==as.character(chrA) & chromB==as.character(chrB)),]
  df_trans_chromB <- df[with(df,chromA==as.character(chrB) & chromB==as.character(chrA)),]
  df_cis_chromA.up <- df_cis_chromA[with(df_cis_chromA,startA > endB),]
  df_cis_chromB.up <- df_cis_chromB[with(df_cis_chromB,startA > endB),]
  df_cis_chromA.dw <- df_cis_chromA[with(df_cis_chromA,endA < startB),]
  df_cis_chromB.dw <- df_cis_chromB[with(df_cis_chromB,endA < startB),]

  cis_chromA.up <- aggregate(with(df_cis_chromA.up,count ~ startA),FUN = sum)
  cis_chromA.dw <- aggregate(with(df_cis_chromA.dw,count ~ startA),FUN = sum)
  cis_chromB.up <- aggregate(with(df_cis_chromB.up,count ~ startA),FUN = sum)
  cis_chromB.dw <- aggregate(with(df_cis_chromB.dw,count ~ startA),FUN = sum)
  trans_chromA  <- aggregate(with(df_trans_chromA,count ~ startA),FUN = sum)
  trans_chromB  <- aggregate(with(df_trans_chromB,count ~ startA),FUN = sum)

  cis_chromA.up.hash <- hashmap(cis_chromA.up$startA,cis_chromA.up$count)
  cis_chromB.up.hash <- hashmap(cis_chromB.up$startA,cis_chromB.up$count)
  cis_chromA.dw.hash <- hashmap(cis_chromA.dw$startA,cis_chromA.dw$count)
  cis_chromB.dw.hash <- hashmap(cis_chromB.dw$startA,cis_chromB.dw$count)
  trans_chromA.hash  <- hashmap(trans_chromA$startA,trans_chromA$count)
  trans_chromB.hash  <- hashmap(trans_chromB$startA,trans_chromB$count)

  df <- files$bed
  df_chromA <- df[with(df,chrom==as.character(chrA)),]
  df_chromB <- df[with(df,chrom==as.character(chrB)),]
  df_chromA[,"cis.up"] <- cis_chromA.up.hash[[df_chromA$start]]
  df_chromA[,"cis.dw"] <- cis_chromA.dw.hash[[df_chromA$start]]
  df_chromB[,"cis.up"] <- cis_chromB.up.hash[[df_chromB$start]]
  df_chromB[,"cis.dw"] <- cis_chromB.dw.hash[[df_chromB$start]]
  df_chromA[,"trans"]  <- trans_chromA.hash[[df_chromA$start]]
  df_chromB[,"trans"]  <- trans_chromB.hash[[df_chromB$start]]

  df_chromA$cis.up[is.na(df_chromA$cis.up)] <- 0
  df_chromA$cis.dw[is.na(df_chromA$cis.dw)] <- 0
  df_chromB$cis.up[is.na(df_chromB$cis.up)] <- 0
  df_chromB$cis.dw[is.na(df_chromB$cis.dw)] <- 0
  df_chromA$trans[is.na(df_chromA$trans)] <- 0
  df_chromB$trans[is.na(df_chromB$trans)] <- 0

  df = rbind(df_chromA,df_chromB)
  write.table(df,file=paste0(prefix,"_",chrA,"_",chrB,"_",id,"_Coverage.bed"),col.names=F,row.names=F,quote=F,sep="\t")
}


#Calculate cis and trans coverage and then calculate directionality index. Followed by HMM segmentation (2 states)
di_and_hmm <- function(coverage,chromA,chromA.start,chromA.end,chromB,chromB.start,chromB.end,prefix,id){

  bed <- as.data.frame(fread(as.character(coverage),header=F))
  colnames(bed) <- c("chrom","start","end","index","up","dw","trans")
  chromA_span <- chromA.end-chromA.start
  chromB_span <- chromB.end-chromB.start
  if (chromA_span < 1e6){
    chromA_span  <- 1e6-chromA_span
    chromA_span  <- floor(chromA_span/2)
    chromA.start <- chromA.start-chromA_span
    chromA.end   <- chromA.end+chromA_span
  }
  if (chromB_span < 1e6){
    chromB_span  <- 1e6-chromB_span
    chromB_span  <- floor(chromB_span/2)+1
    chromB.start <- chromB.start-chromB_span
    chromB.end   <- chromB.end+chromB_span
  }

  bed_chromA <- bed[with(bed,chrom==chromA),]
  bed_chromB <- bed[with(bed,chrom==chromB),]
  bed_chromA_cis_di <- list()
  bed_chromA_trans_di <- list()
  i <- 1
  while (i <= nrow(bed_chromA)){
    up <- bed_chromA$up[i]
    dw <- bed_chromA$dw[i]
    trans <- bed_chromA$trans[i]
    bed_chromA_cis_di[[i]]   <- (dw-up)/abs(dw-up)*(((up-mean(c(dw,up)))^2/mean(c(dw,up)))+((dw-mean(c(dw,up)))^2/mean(c(dw,up))))
    bed_chromA_trans_di[[i]] <- ((trans-mean(bed_chromA$trans))/abs(trans-mean(bed_chromA$trans)))*((trans-mean(bed_chromA$trans))^2/mean(bed_chromA$trans))
    i <- i+1
  }

  smooth <- 5
  bed_chromA_cis_di <- unlist(bed_chromA_cis_di)
  bed_chromA_trans_di <- unlist(bed_chromA_trans_di)
  bed_chromA_cis_di[is.na(bed_chromA_cis_di)] <- 0
  bed_chromA_trans_di[is.na(bed_chromA_trans_di)] <- 0
  bed_chromA_cis_di <- runmean(bed_chromA_cis_di,k=smooth,endrule="mean")
  bed_chromA_trans_di <- runmean(bed_chromA_trans_di,k=smooth,endrule="mean")
  bed_chromA <- cbind(bed_chromA,bed_chromA_cis_di,bed_chromA_trans_di)
  bed_chromA <- bed_chromA[with(bed_chromA,start >= chromA.start & end <= chromA.end),]

  bed_chromB_cis_di <- list()
  bed_chromB_trans_di <- list()
  i <- 1
  while (i <= nrow(bed_chromB)){
    up <- bed_chromB$up[i]
    dw <- bed_chromB$dw[i]
    trans <- bed_chromB$trans[i]
    bed_chromB_cis_di[[i]]   <- (dw-up)/abs(dw-up)*(((up-mean(c(dw,up)))^2/mean(c(dw,up)))+((dw-mean(c(dw,up)))^2/mean(c(dw,up))))
    bed_chromB_trans_di[[i]] <- ((trans-mean(bed_chromB$trans))/abs(trans-mean(bed_chromB$trans)))*((trans-mean(bed_chromB$trans))^2/mean(bed_chromB$trans))
    i <- i+1
  }

  bed_chromB_cis_di <- unlist(bed_chromB_cis_di)
  bed_chromB_trans_di <- unlist(bed_chromB_trans_di)
  bed_chromB_cis_di[is.na(bed_chromB_cis_di)] <- 0
  bed_chromB_trans_di[is.na(bed_chromB_trans_di)] <- 0
  bed_chromB_cis_di <- runmean(bed_chromB_cis_di,k=smooth,endrule="mean")
  bed_chromB_trans_di <- runmean(bed_chromB_trans_di,k=smooth,endrule="mean")
  bed_chromB <- cbind(bed_chromB,bed_chromB_cis_di,bed_chromB_trans_di)
  bed_chromB <- bed_chromB[with(bed_chromB,start >= chromB.start & end <= chromB.end),]

  nstates <- 2
  chromA.hmm <- depmix(response=list(bed_chromA_cis_di~1,bed_chromA_trans_di~1),data=bed_chromA, nstates=nstates,family=list(gaussian(),gaussian()),ntimes=nrow(bed_chromA))
  chromA.hmm.fit <- fit(chromA.hmm)
  chromA.hmm.states <- posterior(chromA.hmm.fit)$state
  x <- summary(chromA.hmm.fit)
  chromA.trans.state <- which.max(as.vector(x[1:nstates,3]))
  cis.hash   <- hashmap(c(1:nstates),as.vector(x[1:nstates,1]))
  trans.hash <- hashmap(c(1:nstates),as.vector(x[1:nstates,3]))
  chromA_cis_DI <- as.vector(x[chromA.trans.state,1])
  bed_chromA <- cbind(bed_chromA,chromA.hmm.states,cis.mean=cis.hash[[chromA.hmm.states]],trans.mean=trans.hash[[chromA.hmm.states]])
  chrA.hmmBoundary.start <- min(bed_chromA[bed_chromA$chromA.hmm.states==chromA.trans.state,]$start)
  chrA.hmmBoundary.end   <- max(bed_chromA[bed_chromA$chromA.hmm.states==chromA.trans.state,]$end)
  write.table(bed_chromA,file=paste0(prefix,".",chromA,"_",id,".hmm.states.bed"),col.names=T,row.names=F,sep="\t",quote=F)

  chromB.hmm <- depmix(response=list(bed_chromB_cis_di~1,bed_chromB_trans_di~1),data=bed_chromB, nstates=nstates,family=list(gaussian(),gaussian()),ntimes=nrow(bed_chromB))
  chromB.hmm.fit <- fit(chromB.hmm)
  chromB.hmm.states <- posterior(chromB.hmm.fit)$state
  x <- summary(chromB.hmm.fit)
  chromB.trans.state <- which.max(as.vector(x[1:nstates,3]))
  cis.hash   <- hashmap(c(1:nstates),as.vector(x[1:nstates,1]))
  trans.hash <- hashmap(c(1:nstates),as.vector(x[1:nstates,3]))
  chromB_cis_DI <- as.vector(x[chromB.trans.state,1])
  bed_chromB <- cbind(bed_chromB,chromB.hmm.states,cis.mean=cis.hash[[chromB.hmm.states]],trans.mean=trans.hash[[chromB.hmm.states]])
  chrB.hmmBoundary.start <- min(bed_chromB[bed_chromB$chromB.hmm.states==chromB.trans.state,]$start)
  chrB.hmmBoundary.end   <- max(bed_chromB[bed_chromB$chromB.hmm.states==chromB.trans.state,]$end)
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

  chromA_span <- ae-as
  chromB_span <- be-bs
  chromA_pos  <- boundary[1]
  chromB_pos  <- boundary[2]
  chromA_up   <- floor(chromA_pos-(chromA_span/2))
  chromA_dw   <- floor(chromA_pos+(chromA_span/2))
  chromB_up   <- floor(chromB_pos-(chromB_span/2))
  chromB_dw   <- floor(chromB_pos+(chromB_span/2))
  chromA_dw.chromB_up <- sum(df[df$startA >= chromA_pos & df$endA <= chromA_dw & df$startB >= chromB_up & df$endB <= chromB_pos,]$count)
  chromA_up.chromB_dw <- sum(df[df$startA >= chromA_up & df$endA <= chromA_pos & df$startB >= chromB_pos & df$endB <= chromB_dw,]$count)
  chromAB_up <- sum(df[df$startA >= chromA_up & df$endA <= chromA_pos & df$startB >= chromB_up & df$endB <= chromB_pos,]$count)
  chromAB_dw <- sum(df[df$startA >= chromA_pos & df$endA <= chromA_dw & df$startB >= chromB_pos & df$endB <= chromB_dw,]$count)

  if ((chromA_cis_DI > 0 & chromB_cis_DI < 0) | (chromA_cis_DI < 0 & chromB_cis_DI > 0)){
    if (chromA_cis_DI > 0 & chromB_cis_DI < 0) {
       if (bk == 1){
          value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
       } else if (bk == 2){
          value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
       } else if (bk == 3){
          value <- (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
       } else if (bk == 4){
          value <- (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
       }
     } else if (chromA_cis_DI < 0 & chromB_cis_DI > 0){
        if (bk == 1){
          value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
        } else if (bk == 2){
          value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
        } else if (bk == 3){
          value <- (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
        } else if (bk == 4){
          value <- (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
        }
      }
   } else if ((chromA_cis_DI > 0 & chromB_cis_DI > 0) | (chromA_cis_DI < 0 & chromB_cis_DI < 0)) {
       if (chromA_cis_DI > 0 & chromB_cis_DI > 0){
         if (bk == 1){
           value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
         } else if (bk == 2){
           value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
         } else if (bk == 3){
           value <- (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
         } else if (bk == 4){
           value <- (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
         }
       } else if (chromA_cis_DI < 0 & chromB_cis_DI < 0){
          if (bk == 1){
            value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_dw
         } else if (bk == 2){
            value <- (chromA_dw.chromB_up+chromA_up.chromB_dw)-chromAB_up
         } else if (bk == 3){
           value <- (chromAB_up+chromAB_dw)-chromA_dw.chromB_up
         } else if (bk == 4){
           value <- (chromAB_up+chromAB_dw)-chromA_up.chromB_dw
         }
       }
   }
   value
}

bp_optimization <- function(intr.file,boundary,chromA,chromB,ssA,seA,ssB,seB){

  intr.file <- as.data.frame(fread(as.character(intr.file),header=T))
  intr.file <- intr.file[intr.file$chromA == as.character(chromA) & intr.file$chromB == as.character(chromB),]

  as <- boundary$chrA.hmmBoundary.start-ssA
  bs <- boundary$chrB.hmmBoundary.start-ssB
  ae <- boundary$chrA.hmmBoundary.end+seA
  be <- boundary$chrB.hmmBoundary.end+seB
  chromA_span <- ae-as
  chromB_span <- be-bs

  upup <- sum(intr.file[with(intr.file,startA >= as & endA <= as+floor(chromA_span/2) & startB >= bs & endB <= bs+floor(chromB_span/2)),]$count)
  dwdw <- sum(intr.file[with(intr.file,startA >= as+floor(chromA_span/2) & endA <= ae & startB >= bs+floor(chromB_span/2) & endB <= be),]$count)
  updw <- sum(intr.file[with(intr.file,startA >= as & endA <= as+floor(chromA_span/2) & startB >= bs+floor(chromB_span/2) & endB <= be),]$count)
  dwup <- sum(intr.file[with(intr.file,startA >= as+floor(chromA_span/2) & endA <= ae & startB >= bs & endB <= bs+floor(chromB_span/2)),]$count)
  mx <- which.max(c(upup,dwdw,updw,dwup))
  chromA_span <- ae-as
  chromB_span <- be-bs
  if (chromA_span < 1e6){
    chromA_span <- 1e6-chromA_span
    chromA_span <- floor(chromA_span/2)
    as <- as-chromA_span
    ae <- ae+chromA_span
  }
  if (chromB_span < 1e6){
    chromB_span <- 1e6-chromB_span
    chromB_span <- floor(chromB_span/2)+1
    bs <- bs-chromB_span
    be <- be+chromB_span
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


############### HiCtrans calling on the binned Hi-C data #########################

## This function is called from HClust to calculate the overlap between the detected segments ##
OverLap <- function(x,y,u,v) {
  if (y < u) {
    return(1)
  } else if (v < x) {
    return(1)
  } else if (y > u & y <= v) {
    return(0)
  } else if (x > u & x <= v) {
    return(0)
  } else if (u > x & v <= y) {
    return(0)
  } else if (v > x & u <= y) {
    return(0)
  } else if (x >= u & y < v) {
    return(0)
  } else if (u >= x & v < y) {
    return(0)
  }
}

## Given distances, this function can cluster a 2D interacting segments ##
HClust <- function(box.df,cl.A,cl.B,colNameA,colNameB,step=0) {

  if (step==0) {
    box.df <- box.df[,c(1:6,8,9,10,12)]
    colnames(box.df) <- c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE","zscore","count","id","resolution")
  }
  chrA.dist <- matrix(ncol=nrow(box.df),nrow=nrow(box.df))
  chrB.dist <- matrix(ncol=nrow(box.df),nrow=nrow(box.df))

  i <- 1
  while (i <= nrow(box.df)) {
    chrA.AS_i <- box.df$BoundaryAS[i]
    chrA.AE_i <- box.df$BoundaryAE[i]
    chrB.BS_i <- box.df$BoundaryBS[i]
    chrB.BE_i <- box.df$BoundaryBE[i]
    j <- 1
    while (j <= nrow(box.df)) {
      chrA.AS_j <- box.df$BoundaryAS[j]
      chrA.AE_j <- box.df$BoundaryAE[j]
      chrB.BS_j <- box.df$BoundaryBS[j]
      chrB.BE_j <- box.df$BoundaryBE[j]
      if (step==0 | step==2) {
        chrA.AS_AS.ij <- abs(chrA.AS_i - chrA.AS_j)
        chrA.AE_AS.ij <- abs(chrA.AE_i - chrA.AS_j)
        chrA.AS_AE.ij <- abs(chrA.AS_i - chrA.AE_j)
        chrA.AE_AE.ij <- abs(chrA.AE_i - chrA.AE_j)
        chrB.BS_BS.ij <- abs(chrB.BS_i - chrB.BS_j)
        chrB.BE_BS.ij <- abs(chrB.BE_i - chrB.BS_j)
        chrB.BS_BE.ij <- abs(chrB.BS_i - chrB.BE_j)
        chrB.BE_BE.ij <- abs(chrB.BE_i - chrB.BE_j)

        chrA.dist[i,j] <- min(chrA.AS_AS.ij,chrA.AE_AS.ij,chrA.AS_AE.ij,chrA.AE_AE.ij)
        chrB.dist[i,j] <- min(chrB.BS_BS.ij,chrB.BE_BS.ij,chrB.BS_BE.ij,chrB.BE_BE.ij)

      } else {

        chrA.dist[i,j] <- OverLap(chrA.AS_i,chrA.AE_i,chrA.AS_j,chrA.AE_j)
        chrB.dist[i,j] <- OverLap(chrB.BS_i,chrB.BE_i,chrB.BS_j,chrB.BE_j)

      }
      j <- j + 1
    }
    i <- i + 1
  }

  if (nrow(chrA.dist) > 1 & nrow(chrB.dist) > 1) {
    chrA.hclust <- hclust(as.dist(chrA.dist),method="single")
    chrB.hclust <- hclust(as.dist(chrB.dist),method="single")
    chrA.clust  <- cutree(chrA.hclust,h=cl.A)
    chrB.clust  <- cutree(chrB.hclust,h=cl.B)
    box.df[,"chrA.clus"] <- chrA.clust
    box.df[,"chrB.clus"] <- chrB.clust
  } else {
    chrA.clust <- 1
    chrB.clust <- 1
    box.df[,"chrA.clus"] <- chrA.clust
    box.df[,"chrB.clus"] <- chrB.clust
  }
  clus <- unique(data.frame(chrA.clust,chrB.clust))
  if (step==1) {
    return(box.df)
  }
  if (step==0) {
    box.cl <- list()
    i <- 1
    while (i <= nrow(clus)) {
      d <- box.df[box.df$chrA.clus==clus$chrA.clust[i] & box.df$chrB.clus==clus$chrB.clust[i],]
      box.cl[[i]] <- data.frame(chrA=d$chrA[1],BoundaryAS=min(d$BoundaryAS),BoundaryAE=max(d$BoundaryAE),
                                chrB=d$chrB[1],BoundaryBS=min(d$BoundaryBS),BoundaryBE=max(d$BoundaryBE),
                                zscore=mean(d$zscore),count=mean(d$count),resolution=d$resolution[1],
    		                id=paste0(d$id,collapse=","),A=clus$chrA.clust[i],B=clus$chrB.clust[i])
      i <- i + 1
    }
    box.cl <- as.data.frame(do.call(rbind,box.cl))
    colnames(box.cl) <- c("chrA","BoundaryAS","BoundaryAE","chrB","BoundaryBS","BoundaryBE","zscore","count","resolution","id",as.character(colNameA),as.character(colNameB))
    return(box.cl)
  }
}

## Finding out the common multi-resolution translocation boxes and breakpoints ##
ZoomIn <- function(d, f, r, l) {

  g <- 1
  clus <- unique(data.frame(d[,c("chrA.clus","chrB.clus")]))
  w <- list()
  i <- 1
  while (i <= nrow(clus)) {
    v <- d[d$chrA.clus==clus$chrA.clus[i] & d$chrB.clus==clus$chrB.clus[i],]
    b <- list()
    j <- 1
    while (j <= length(r)) {
      b[[j]] <- v[v$resolution==r[j],]
      if (nrow(b[[j]]) > 0) {
        b[[j]][,"Level"] <- paste0("L_",j)
      }
      j <- j + 1
    }
    b <- do.call(rbind,b)
    if (nrow(b) >= l) {
      w[[i]] <- b
      w[[i]][,"group"] <- g
      g <- g + 1
    }
    i <- i + 1
  }
  w <- as.data.frame(do.call(rbind, w))
  g <- unique(w$g)
  b <- list()
  h <- list()
  n <- 1
  i <- 1
  while (i <= length(g)) {
    d <- w[w$g==g[i],]
    d <- d[d$Level==d[nrow(d),]$Level,]
    j <- 1
    while (j <= nrow(d)) {
      k <- as.integer(strsplit(as.character(d$id[j]),",")[[1]])
      r <- d$resolution[j]
      b[[n]] <- f[f$resolution==r & f$id %in% k & f$type=="BreakPoint",]
      j <- j + 1
      n <- n + 1
    }
    h[[i]] <- w[w$g==g[i],]
    h[[i]] <- h[[i]][h[[i]]$Level==h[[i]][1,]$Level,]
    i <- i + 1
  }
  b <- as.data.frame(do.call(rbind,b))
  if (nrow(b) > 0) {
    b <- HClust(box.df=b,cl.A=min(b$resolution),cl.B=min(b$resolution),colNameA="resA.cl",colNameB="resB.cl",step=0)
    h <- as.data.frame(do.call(rbind,h))
    b <- b[,c(1:9)]
    h <- h[,c(1:9)]
    b[,"class"] <- "BreakPoint"
    h[,"class"] <- "TranslocationBox"
    df <- as.data.frame(rbind(h,b))
  } else {
    df <- data.frame(chrA=c(),BoundaryAS=c(),BoundaryAE=c(),chrB=c(),BoundaryBS=c(),BoundaryBE=c(),zscore=c(),count=c(),resolution=c(),id=c(),A=c(),B=c(),class=c())
  }
  return(df)
}

## From ijk type data, this function will create a matrix object ##
CreateMatrix <- function(ijk.file, bed.file, chromA, chromB, prefix) {

  ijk <- as.data.frame(fread(ijk.file,h=F))
  bed <- as.data.frame(fread(bed.file,h=F))
  colnames(ijk) <- c("A","B","C")
  colnames(bed) <- c("chr","start","end","index")
  bed <- bed[bed$chr %in% c(chromA,chromB),]
  indexA.conversion <- hashmap(bed[bed$chr==chromA,]$index,c(1:nrow(bed[bed$chr==chromA,])))
  indexB.conversion <- hashmap(bed[bed$chr==chromB,]$index,c(1:nrow(bed[bed$chr==chromB,])))
  ijk <- ijk[ijk$A %in% bed[bed$chr==chromA,]$index,]
  ijk <- ijk[ijk$B %in% bed[bed$chr==chromB,]$index,]
  ijk$A <- indexA.conversion[[ijk$A]]
  ijk$B <- indexB.conversion[[ijk$B]]
  bed[bed$chr==chromA,]$index <- 1:nrow(bed[bed$chr==chromA,])
  bed[bed$chr==chromB,]$index <- 1:nrow(bed[bed$chr==chromB,])
  mat <- matrix(0,nrow=nrow(bed[bed$chr==chromA,]),ncol=nrow(bed[bed$chr==chromB,]))
  print (dim(mat))
  i <- 1
  while (i <= nrow(ijk)) {
    mat[ijk$A[i],ijk$B[i]] <- ijk$C[i]
    i <- i + 1
  }
  write.table(mat,file=paste0(prefix,".",chromA,"_",chromB,".mat.txt"),row.names=F,col.names=F,sep="\t",quote=F)
  return(list(mat=as.matrix(mat),bedA=as.data.frame(bed[bed$chr==chromA,]),bedB=as.data.frame(bed[bed$chr==chromB,])))
}

## This function scans the inter-chromosomal interaction matrix for mean values with a given box size ## 
cppFunction('NumericVector lmat(NumericMatrix mat, int x, int y) {
  int n = mat.nrow();
  int m = mat.ncol();
  int k = 0;
  int s1 = round(n/x);
  int s2 = round(m/y);
  int s3 = (s1 * s2);
  Rcpp::NumericVector M(s3);
  for(int i = 0; i < n; i = i + x) {
    int a = i;
    int b = i + (x-1);
    for (int j = 0; j < m; j = j + y){
      int u = j;
      int v = j + (y-1);
      if (b < mat.nrow()) {
        if (v < mat.ncol()) {
          NumericMatrix F = mat(Range(a,b),Range(u,v));
          M[k] = findMean(F);
          k = k + 1;
        }
      }
    }
  }
  return(M);
}', includes='double findMean(NumericMatrix F) {
  int sum = 0;
  int n = F.nrow();
  int m = F.ncol();
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      sum += F(i,j);
    }
  }
  return (double)sum/(n*m);
}')


## This function shrinks the Hi-C matrix from high to low resolution ##
cppFunction('NumericMatrix ShrinkMatrixCPP(NumericMatrix mat_high, NumericMatrix x_coord, NumericMatrix y_coord) {
  NumericMatrix mat_low(x_coord.nrow(),y_coord.nrow());
  for (int i = 0; i < x_coord.nrow(); i++) {
    int x_0 = x_coord(i,1) - 1;
    int x_1 = x_coord(i,2) - 1;
    for (int j = 0; j < y_coord.nrow(); j++) {
      int y_0 = y_coord(j,1) - 1;
      int y_1 = y_coord(j,2) - 1;
      NumericMatrix sub_mat = mat_high(Range(x_0,x_1),Range(y_0,y_1));
      int sum = 0;
      for (int u = 0; u < sub_mat.nrow(); u++) {
        for (int v = 0; v < sub_mat.ncol(); v++) {
          sum += sub_mat(u,v);
        }
      }
      mat_low(i,j) = sum;
    }
  }
  return(mat_low);
}')


## For a detected translocation box, this function detects the top interacting region ##
EstimateMaxRegion <- function(box, mat, locq) {
 
  b <- list()
  k <- 1
  i <- 1
  while (i <= nrow(box)) {
    s.1 <- box$start1[i]
    e.1 <- box$end1[i]
    s.2 <- box$start2[i]
    e.2 <- box$end2[i]
    m <- mat[s.1:e.1,s.2:e.2]
    rownames(m) <- c(s.1:e.1)
    colnames(m) <- c(s.2:e.2)
    m[m > 0] <- 1
    x.val <- rowSums(m)/sum(rowSums(m))
    y.val <- colSums(m)/sum(colSums(m))
    val.mat <- outer(x.val,y.val)
    val.mat <- mat[s.1:e.1,s.2:e.2] * val.mat
    p <- which(val.mat >= quantile(val.mat, locq), arr.ind = TRUE)
    x.pos <- as.integer(rownames(m)[p[,1]])
    y.pos <- as.integer(colnames(m)[p[,2]])
    j <- 1
    while (j <= nrow(p)) {
      if (x.pos[j] >= 2 & x.pos[j] < nrow(mat) & y.pos[j] >= 2 & y.pos[j] < ncol(mat)) {
        max.count <- max(mat[c(x.pos[j]-1):c(x.pos[j]+1),c(y.pos[j]-1):c(y.pos[j]+1)])
        b[[k]] <- cbind(box[i,],pos.1=x.pos[j],pos.2=y.pos[j],count=max.count)
        k <- k + 1
      }
      j <- j + 1
    }
    i <- i + 1
  }
  b <- as.data.frame(do.call(rbind,b))
  return(b)
}

## This function calls the lmat function described above ##
BackgroundEnrichment <- function(i, mat, f) {
 
  print (f[i,])
  x.s <- f$d1[i]
  y.s <- f$d2[i]
  val <- lmat(mat,x.s,y.s)
  return(data.frame(MEAN=mean(val),SD=sd(val)))

}

## Check for all combination of 1D segments from both the chromosomes, the enrichment of their interaction ##
BoxEnrichment <- function(m, x, y, core, cutoff=0.05) {

  n <- 0
  repeat {
    print (paste0("Estimating box enrichment: Trimmining ",cutoff*1e2,"% of values"))
    d <- list()
    k <- 1
    i <- 1
    while (i <= nrow(x)) {
      a <- x$start[i]
      b <- x$end[i]
      j <- 1
      while (j <= nrow(y)) {
        u <- y$start[j]
        v <- y$end[j]
        d[[k]] <- data.frame(start1=x$start[i],end1=x$end[i],d1=(x$end[i]-x$start[i]),
                            start2=y$start[j],end2=y$end[j],d2=(y$end[j]-y$start[j]),
                            val=mean(m[a:b,u:v],trim=cutoff))
        k <- k + 1
        j <- j + 1
      }
      i <- i + 1
    }
    d <- as.data.frame(do.call(rbind,d))
    d[,"val.center"] <- d$val - mean(d$val)
    d <- d[d$val.center > 0 & d$d1 > 2 & d$d2 > 2,]
    cutoff <- cutoff - 0.04
    if (nrow(d) > 0) {
      break
    }
    if (cutoff < 0) {
      n <- n + 1
      if (n == 1) {
        cutoff <- 0
      } else if (n > 1) {
        print ("No enrichemnt of changepoint boxes detected! Try to run with Lower resolution")
        break
      }
    }
    print ("Lowering trimming threshold")
  }  
  if (nrow(d) > 0) {
    v <- lapply(c(1:nrow(d)), BackgroundEnrichment, m, d)
    v <- as.data.frame(do.call(rbind,v))
    d[,"val.bg.mean"] <- v$MEAN
    d[,"val.bg.sd"]   <- v$SD
    d[,"zscore"]      <- ifelse(is.na(v$SD), NA, (d$val - d$val.bg.mean)/d$val.bg.sd)
    d <- na.omit(d) 
    d <- d[d$zscore >= 1,]
    print (d)
    if (nrow(d) > 0) { 
      return(d)
    } else {
      print (paste0("No enriched segment found"))
      return(data.frame())
    }
  } else {
    print ("No enrichemnt of changepoint boxes detected! Try to run with Lower resolution")
    return(data.frame())
  }
}

## Given a change-point vector, convert the vector to data.frame like object ##
CptsInBedFormat <- function(cpt.obj) {

  start <- list()
  end   <- list()
  val   <- list()
  start[[1]] <- 1
  end[[1]]   <- cpt.obj@cpts[1]
  val[[1]]   <- cpt.obj@param.est$mean[1]
  i <- 2
  while (i <= length(cpt.obj@cpts)) {
    start[[i]] <- cpt.obj@cpts[i-1]
    end[[i]]   <- cpt.obj@cpts[i]
    val[[i]]   <- cpt.obj@param.est$mean[i]
    i <- i + 1
  }
  start <- unlist(start)
  end   <- unlist(end) 
  val   <- unlist(val)
  return(data.frame(start,end,val)) 
}

## This is the main function which call the binary-segmentation algorithm and detects 1D segmentations ##
GetSegments <- function(mat, bed, chromA, chromB, prefix, covq=0.75, locq, mincount, minzscore, resolution, glbq) {

  r.sum <- rowSums(mat)
  c.sum <- rowSums(t(mat))
  r.center <- r.sum - quantile(r.sum, covq)
  c.center <- c.sum - quantile(c.sum, covq)
 
  seg <- ifelse(resolution > 250000, 50, 100) 
  r.cpt <- suppressWarnings(cpt.meanvar(r.center,Q=seg,method="BinSeg"))
  c.cpt <- suppressWarnings(cpt.meanvar(c.center,Q=seg,method="BinSeg"))

  ## Flush changepoint summary ## 
  print ("Row and Column changepoint summary")
  print (r.cpt)
  print (c.cpt)
  cat ("\n")
  ###############################
  
  r.df <- CptsInBedFormat(r.cpt)
  c.df <- CptsInBedFormat(c.cpt)
  r.df <- r.df[r.df$val > 0,]
  c.df <- c.df[c.df$val > 0,]
  
  if (nrow(r.df) > 0 & nrow(c.df) > 0) {
    ## Flush out the changepoints and their mean values ##
    print ("Row segments with positive values")
    print (r.df)
    print ("Column segments with positive values")
    print (c.df)
    ######################################################

    box <- BoxEnrichment(mat,r.df,c.df,core)
    if (nrow(box) > 0) {
      box <- EstimateMaxRegion(box, mat, locq)
      if (nrow(box) > 0) {
        chrA.start.hash <- hashmap(bed[bed$chr==chromA,]$index,bed[bed$chr==chromA,]$start)
        chrB.start.hash <- hashmap(bed[bed$chr==chromB,]$index,bed[bed$chr==chromB,]$start)  
        chrA.end.hash   <- hashmap(bed[bed$chr==chromA,]$index,bed[bed$chr==chromA,]$end)
        chrB.end.hash   <- hashmap(bed[bed$chr==chromB,]$index,bed[bed$chr==chromB,]$end)
        box[,"chr1"] <- chromA
        box[,"x1"]   <- chrA.start.hash[[box$start1]]
        box[,"y1"]   <- chrA.start.hash[[box$end1]] 
        box[,"chr2"] <- chromB
        box[,"x2"]   <- chrB.start.hash[[box$start2]]
        box[,"y2"]   <- chrB.start.hash[[box$end2]]
        box[,"chrA"] <- chromA
        box[,"chrA.Brk.pos1"] <- chrA.start.hash[[box$pos.1]]
        box[,"chrA.Brk.pos2"] <- chrA.end.hash[[box$pos.1]]
        box[,"chrB"] <- chromB
        box[,"chrB.Brk.pos1"] <- chrB.start.hash[[box$pos.2]]
        box[,"chrB.Brk.pos2"] <- chrB.end.hash[[box$pos.2]]
        write.table(box,file=paste0(prefix,".",chromA,"_",chromB,".Details.txt"),row.names=F,sep="\t",quote=F)
   
        wt.threshold <- data.frame(count=as.integer(names(table(as.vector(mat)))),freq=as.vector(table(as.vector(mat)))) 
        wt.threshold <- wt.threshold[-c(1),]
        wt.threshold[,"csum"] <- cumsum(wt.threshold$freq)/sum(wt.threshold$freq)
        wt.threshold <- min(wt.threshold[wt.threshold$csum >= glbq,]$count)
        wt.threshold <- ifelse(wt.threshold < mincount, mincount, wt.threshold)
        print (wt.threshold)
        box[,"color"] <- "0,0,255"
        box  <- box[,c("chr1","x1","y1","chr2","x2","y2","color","zscore","chrA","chrA.Brk.pos1","chrA.Brk.pos2","chrB","chrB.Brk.pos1","chrB.Brk.pos2","count")]
        box  <- box[box$count >= wt.threshold & box$zscore >= minzscore,]
        if (nrow(box) > 0) {
          box[,"id"] <- 1:nrow(box)
          boxA <- box[,c("chr1","x1","y1","chr2","x2","y2","color","zscore","count","id")]
          boxB <- data.frame(chr1=box$chrA,x1=box$chrA.Brk.pos1,y1=box$chrA.Brk.pos2,
  			 chr2=box$chrB,x2=box$chrB.Brk.pos1,y2=box$chrB.Brk.pos2,
    			 color=box$color,zscore=box$zscore,count=box$count,id=box$id) 
          boxB[,"color"] <- "0,0,0"
          boxA[,"type"]  <- "TranslocationBox"
          boxB[,"type"]  <- "BreakPoint"
          box <- rbind(boxA,boxB)
          box[,"resolution"] <- resolution
          write.table(box,file=paste0(prefix,".",chromA,"_",chromB,".Translocations_jcbx.txt"),row.names=F,sep="\t",quote=F)
          return(box)
        } else {
          print ("Minimum count filter not passed. Try to decrease the minimum count to get a breakpoint at this resolution")
          return(data.frame())
        }
      } else {
        print ("Minimum count percentile not passed. Try to decrease the minimum count percentile to get a breakpoint at this resolution")
        return(data.frame())
      }
    } else {
      print ("Minimum count filter not passed. Try to decrease the minimum count to get a breakpoint at this resolution")
      return(data.frame())
    }
  } else {
    if (covq > 0.5) {
      print (paste0("Descresing the quantile value to ",(covq-0.25)," to get positive segments"))
      GetSegments(mat, bed, chromA, chromB, prefix, covq=(covq-0.25), locq, mincount, minzscore, resolution, glbq)
    } else {
      print ("No positive segments from changepoint analysis! Try a lower resolution hic data for this chromosome")
      return(data.frame())
    }
  }
}


CoordinateShrunkage <- function(v,len) {

  d <- list()
  j <- 1
  while(j < length(v)) {
    d[[j]] <- data.frame(p=j,s=v[j],e=(v[j+1]-1))
    j <- j + 1
  }
  if (v[j]==len) {
    d[[j]] <- data.frame(p=j,s=v[j],e=v[j])
  } else {
    d[[j]] <- data.frame(p=j,s=v[j],e=len)
  }
  d <- as.data.frame(do.call(rbind,d))
  return(d)

}

ShrinkMatrix <- function(mat.high,x.coord,y.coord) {

  cat ("Shrinking Matrix to lower resolution\n")
  x.coord <- as.matrix(x.coord)
  y.coord <- as.matrix(y.coord)
  colnames(x.coord) <- NULL
  colnames(y.coord) <- NULL
  mat.low <- ShrinkMatrixCPP(mat.high,x.coord,y.coord)
  return(mat.low)
}

## This function provide a Lower-resolution Hi-C matrix from a Higher-resolution Hi-C matrix ##
HighToLowResolutionCoversion <- function(mat, high.resolution, low.resolution, bed.x, bed.y, bed, prefix) {

  mat_list     <- list()
  bed_list     <- list()
  x.coord_list <- list()
  y.coord_list <- list()
  mat_list[[1]] <- mat
  bed_list[[1]] <- bed
  i <- 1
  while (i <= length(low.resolution)) {
    x.start <- hashmap(c(bed.x$index),c(bed.x$start))
    x.end   <- hashmap(c(bed.x$index),c(bed.x$end))
    y.start <- hashmap(c(bed.y$index),c(bed.y$start))
    y.end   <- hashmap(c(bed.y$index),c(bed.y$end))
    scale.factor <- low.resolution[i]/high.resolution
    x.dim <- nrow(mat)
    y.dim <- ncol(mat)
    x.divison <- seq(1,x.dim,by=scale.factor)
    y.divison <- seq(1,y.dim,by=scale.factor)
    x.coord <- CoordinateShrunkage(x.divison,x.dim)
    y.coord <- CoordinateShrunkage(y.divison,y.dim)
    bedA <- data.frame(chr=rep(as.character(bed.x$chr[1]),nrow(x.coord)),start=x.start[[x.coord$s]],end=x.end[[x.coord$e]],index=c(1:nrow(x.coord)))
    bedB <- data.frame(chr=rep(as.character(bed.y$chr[1]),nrow(y.coord)),start=y.start[[y.coord$s]],end=y.end[[y.coord$e]],index=c(1:nrow(y.coord)))   
    x.coord_list[[i]] <- cbind(bedA,x.coord)
    y.coord_list[[i]] <- cbind(bedB,y.coord)
    mat_list[[i+1]]   <- ShrinkMatrix(mat,x.coord,y.coord)
    bed_list[[i+1]]   <- data.frame(rbind(bedA,bedB))

    write.table(mat_list[[i+1]],file=paste0(prefix,"_",as.integer(low.resolution[i]),".",bedA$chr[1],"_",bedB$chr[1],".matrix"),row.names=F,col.names=F,sep="\t",quote=F)
    write.table(bed_list[[i+1]],file=paste0(prefix,"_",as.integer(low.resolution[i]),".",bedA$chr[1],"_",bedB$chr[1],"_abs.bed"),row.names=F,col.names=F,sep="\t",quote=F)

    i <- i + 1
  }
  d <- list(mat=mat_list,bed=bed_list,x.coord=x.coord_list,y.coord=y.coord_list)
  return(d)
}

option_list = list(
  make_option(c("--mat"), type="character", help="An upper triangular Hi-C sparse matrix\n\t\tIt should have the following columns\n
 \t\t<indexA> <indexB> <count>\n
 \t\t1 1 300
 \t\t1 2 30
 \t\t1 3 10
 \t\t2 2 200
 \t\t2 3 20
 \t\t3 3 200
 \t\t....\n"),

  make_option(c("--bed"), type="character", help="Bed file with index information\n\t\tIt should have the following columns\n\n
 \t\t<chr> <start> <end> <index>\n
 \t\tchr1 1 40000 1
 \t\tchr1 40000 80000 2
 \t\tchr1 80000 120000 3
 \t\t....\n"),

  make_option(c("--chrA"), type="character", help="Chromosome A name. It will represent the rows in the inter-chromosomal matrix. It should be the <indexA> chromosome.\n"),

  make_option(c("--chrB"), type="character", help="Chromosome B name. It will represent the columns in the inter-chromosomal matrix. It should be the <indexB> chromosome.\n"),

  make_option(c("--prefix"), type="character", help="Prefix of the output file <prefix.chrA_chrB>. All the output files and folders will be generated with this prefix.\n"),

  make_option(c("--covq"), type="numeric", default=0.1, help="Quantile value to be subtracted from one dimensional trans-coverage profile [trans.coverage - quantile(trans.coverage, covq)] [default 0.10].\n
  \t\tBins with very low coverage values are removed with this filter. 
  \t\tIncreasing <covq> value will keep only the most stringent bins in the two chromosome.\n"),

  make_option(c("--minzscore"), type="numeric", default=1, help="Minimum Zscore of a possible translocation box to be retained [default is 1].\n
  \t\tHiCtrans will find enriched boxes within the inter-chromosomal matrix as potential translocation box.
  \t\tThe enrichment is calculated as Z-score against a background with all possible similar sized boxes in the inter-chromosomal matrix.
  \t\tIncreasing <minzscore> value will keep the most enriched trans interacting boxes.\n"),

  make_option(c("--locq"), type="numeric", default=0.1, help="Top percentile to be reported as possible breakpoints within a translocation box [default top 0.1%]\n
  \t\tFor each enriched translocated boxes, HiCtrans will report top <locq>% interacting pairs (Weighted by the frequency of total interaction).
  \t\tDecreasing <locq> will reduce the number of reported breakpoints within an enriched trans interacting box. A <locq> value of 0 will report only the top interacting pair\n"),

  make_option(c("--mincount"), type="numeric", default=10, help="Minimum count of a possible breakpoint to be retained when compared to all possible chrA-chrB interaction [default cutoff is 10].\n
  \t\tThis is an absolute minimum count cutoff to filter out breakpoints detected at any resolution.
  \t\tIncreasing <mincount> value will keep only the most stringent interacting pair.\n"),

  make_option(c("--glbq"), type="numeric", default=0.1, help="Percentile value for minimum count cutoff at each resolution [default cutoff is at top 0.1% of the count distribution].\n
  \t\tThis is a relative count cutoff based on the inter-chromosomal count distribution determined for each resolution independently.   
  \t\tIncreasing <glbq> value keep only the most stringent interacting pair.\n"),

  make_option(c("--resolutions"), type="character", default="2,4,5,10", help="Comma separated list of integers to be multiplied with the starting Hi-C resolution to get the lower resolutions [default 2,4,5,10].\n
  \t\tHiCtrans will search for enriched trans-interacting boxes and breakpoint finding at different resolutions.
  \t\tProvide only integer values in a comma separated list.\n"),

  make_option(c("--multires"), type="numeric", default=2, help="Number of Hi-C resolutions at which the breakpoint should be supported with [default is at least 2 different resolutions].\n
  \t\tHiCtrans wiil find enriched trans-interacting boxes and subsequent breakpoints at different resolutions. The ultimate goal is to find a true translocation supported by multiple resolutions.
  \t\tIncreasing <multires> value will keep only the enriched boxes and breakpoints supported by at least <multires> number of different resolution.\n"),

  make_option(c("--maxres"), type="numeric", default=3, help="Maximum resolution upto which the breakpoint is kept after multi-resolution filtering [default is 3 X user provided Hi-C resolution].\n
  \t\tHiCtrans wiil find enriched trans-interacting boxes and subsequent breakpoints at different resolutions. The ultimate goal is to find a true translocation supported by the highest resolutions.
  \t\tIncreasing <maxres> value will keep only the enriched boxes and breakpoints supported by upto <maxres> X starting Hi-C resolutions.\n"),

  make_option(c("--relevel"), type="character", default="No", help="Should the breakpoints be refined upto restriction-level resolution [default is 'No'; If 'Yes', the following parameters are MUST]\n"),

  make_option(c("--fragsFile"), type="character", help="Restriction Fragment file [MUST].\n
 \t\tchr1    0       16007   HIC_chr1_1      0       +
 \t\tchr1    16007   24571   HIC_chr1_2      0       +
 \t\tchr1    24571   27981   HIC_chr1_3      0       +
 \t\t......\n"),

  make_option(c("--chromsize"), type="character", help="Chromosome size file [MUST].\n
 \t\tchr1\t249250621
 \t\tchr2\t243199373
 \t\tchr3\t198022430
 \t\tchr4\t191154276
 \t\t.....\n"),

  make_option(c("--validpair"), type="character", help="Valid pair file of the HiC data [MUST].\n
 \t\tSRR6213722.1\tchr11\t124331538\t-\tchr11\t124345246\t-
 \t\tSRR6213722.2\tchr1\t198436365\t-\tchr1\t199923196\t+
 \t\t.....\n"),
 
  make_option(c("--clusdist"),  default=1e6,type="integer", help="Distance threshold in basepairs to cluster the nearby breakpoints obtained from multi-resolution filtered (MultiResolution_Filtered.Translocation.txt) or individual Translocations_jcbx.txt files [Default 1Mb]"),
  make_option(c("--ssA"), default=1e5, type="integer", help="Extend -(ve) bp of the 5' HMM segment border of chromosome A for breakpoint identification. Default 100Kb."),
  make_option(c("--seA"), default=1e5, type="integer", help="Extend +(ve) bp of the 3' HMM segment border of chromosome A for breakpoint identification. Default 100Kb."),
  make_option(c("--ssB"), default=1e5, type="integer", help="Extend -(ve) bp of the 5' HMM segment border of chromosome B for breakpoint identification. Default 100Kb."),
  make_option(c("--seB"), default=1e5, type="integer", help="Extend +(ve) bp of the 3' HMM segment border of chromosome B for breakpoint identification. Default 100Kb.")

)
opt = parse_args(OptionParser(option_list=option_list))

prefix <- as.character(opt$prefix)

bed <- as.character(opt$bed)
bed <- read.table(bed,h=F)
colnames(bed) <- c("chr","start","end","index")

high.resolution <- bed$end[1] - bed$start[1]
low.resolution  <- as.integer(strsplit(as.character(opt$resolutions),",")[[1]]) * high.resolution
high.resolution <- as.integer(high.resolution)
low.resolution  <- as.integer(low.resolution)
multires <- as.integer(opt$multires)
maxres <- as.integer(opt$maxres) * high.resolution
chromA <- as.character(opt$chrA)
chromB <- as.character(opt$chrB)

if (!file.exists(paste0(prefix,"_hictrans"))) {
  dir.create(paste0(prefix,"_hictrans"))
}
if (!file.exists(paste0(prefix,"_hictrans/",prefix,"_hictrans_",chromA,"_",chromB,"_",as.integer(high.resolution)))) {
  dir.create(paste0(prefix,"_hictrans/",prefix,"_hictrans_",chromA,"_",chromB,"_",as.integer(high.resolution)))
  mat.dir.name <- dirname(opt$mat)
  bed.dir.name <- dirname(opt$bed)

  mat_pat <- normalizePath(opt$mat)
  bed_pat <- normalizePath(opt$bed)
 
  setwd(paste0(prefix,"_hictrans/",prefix,"_hictrans_",chromA,"_",chromB,"_",as.integer(high.resolution)))
  cmd <- paste0("ln -s ",mat_pat)
  system(cmd, wait=T)
  cmd <- paste0("ln -s ",bed_pat)
  system(cmd, wait=T)
  
  if (opt$relevel == "Yes") {
    fragsFile_path <- normalizePath(opt$fragsFile)
    chromsize_path <- normalizePath(opt$chromsize)
    validpair_path <- normalizePath(opt$validpair)
    cmd <- paste0("ln -s ",fragsFile_path)
    system(cmd, wait=T)
    cmd <- paste0("ln -s ",chromsize_path)
    system(cmd, wait=T) 
    cmd <- paste0("ln -s ",validpair_path)
    system(cmd, wait=T)
  }
 
} else {

  mat.dir.name <- dirname(opt$mat)
  bed.dir.name <- dirname(opt$bed)

  mat_pat <- normalizePath(opt$mat)
  bed_pat <- normalizePath(opt$bed)

  setwd(paste0(prefix,"_hictrans/",prefix,"_hictrans_",chromA,"_",chromB,"_",as.integer(high.resolution)))
  cmd <- paste0("ln -s ",mat_pat)
  system(cmd, wait=T)
  cmd <- paste0("ln -s ",bed_pat)
  system(cmd, wait=T)

  if (opt$relevel == "Yes") {
    fragsFile_path <- normalizePath(opt$fragsFile)
    chromsize_path <- normalizePath(opt$chromsize)
    validpair_path <- normalizePath(opt$validpair)
    cmd <- paste0("ln -s ",fragsFile_path)
    system(cmd, wait=T)
    cmd <- paste0("ln -s ",chromsize_path)
    system(cmd, wait=T)
    cmd <- paste0("ln -s ",validpair_path)
    system(cmd, wait=T)
  }
}

bed <- bed[bed$chr %in% c(chromA,chromB),]
bed[bed$chr==chromA,]$index <- 1:nrow(bed[bed$chr==chromA,])
bed[bed$chr==chromB,]$index <- 1:nrow(bed[bed$chr==chromB,])

mat <- CreateMatrix(as.character(opt$mat),
	       as.character(opt$bed),
	       as.character(opt$chrA),
	       as.character(opt$chrB), 
	       as.character(opt$prefix)) 

mat <- HighToLowResolutionCoversion(mat$mat,high.resolution,low.resolution, mat$bedA, mat$bedB, bed, prefix)

all.resolution <- sort(c(high.resolution,low.resolution))
translocation.boxes <- list()

sink(paste0(prefix,".",chromA,"_",chromB,".log.txt"))
print(paste0("Started ",Sys.time()))
for (i in 1:length(all.resolution)) {
    
  covq <- opt$covq
  locq <- opt$locq
  glbq <- opt$glbq
  locq <- 1 - (locq/100)
  glbq <- 1 - (glbq/100)  
  mincount <- as.integer(opt$mincount)
  minzscore<- as.integer(opt$minzscore)
  prefix_resolution <- paste0(prefix,"_",as.integer(all.resolution[i]))
  print (paste0("Running for ",all.resolution[i]," resolution"))
  translocation.boxes[[i]] <- GetSegments(mat$mat[[i]], mat$bed[[i]], chromA, chromB, prefix_resolution, covq, locq, mincount, minzscore, as.integer(all.resolution[i]), glbq)

}
translocation.boxes <- as.data.frame(do.call(rbind, translocation.boxes))
write.table(translocation.boxes, file=paste0(prefix,"_hictrans.",chromA,"_",chromB,".preCluster.txt"),row.names=F,sep="\t",quote=F)
print ("Clustering the Translocation boxes at resolution level")

print ("Filtering breakpoint from multi-resolution data")
if (nrow(translocation.boxes) > 0) {
  resolution <- sort(unique(as.integer(translocation.boxes$resolution)),decreasing=T)
  if (length(resolution) >= opt$multires) {
    data <- list()
    for(i in 1:length(resolution)) {
      box.tmp  <- translocation.boxes[translocation.boxes$resolution==resolution[i] & translocation.boxes$type=="TranslocationBox",]
      if (nrow(box.tmp) > 0) {
        data[[i]] <- HClust(box.df=box.tmp,cl.A=0,cl.B=0,colNameA="resA.cl",colNameB="resB.cl",step=0)
      }
    }
    data <- as.data.frame(do.call(rbind,data))
    data <- HClust(box.df=data,cl.A=0,cl.B=0,colNameA="boxA.cl",colNameB="boxB.cl",step=1)
    data <- ZoomIn(data, translocation.boxes, resolution, multires)
    data <- data[data$resolution <= maxres,]
    if (nrow(data) > 0) {
      write.table(data, file=paste0(prefix,"_hictrans.",chromA,"_",chromB,".MultiResolution_Filtered.Translocation.txt"),row.names=F,sep="\t",quote=F)
    }
    ## Call restriction resolution breakpoint filtering here ## 
    if (opt$relevel == "Yes"){
      param <- data
      param <- param[param$class=="BreakPoint",]
      param <- RE_HClust(param,cl.A=opt$clusdist,cl.B=opt$clusdist)
      chromA_BP_start <- list()
      chromA_BP_end <- list()
      chromB_BP_start <- list()
      chromB_BP_end <- list()
      k <- 1
      while (k <= nrow(param)){
        opt$chrA <- param$chrA[k]
        opt$startA <- param$BoundaryAS[k]
        opt$endA   <- param$BoundaryAE[k]
        opt$chrB <- param$chrB[k]
        opt$startB <- param$BoundaryBS[k]
        opt$endB   <- param$BoundaryBE[k]

        df <- mapVPs_on_REs(
              re.bed=as.character(opt$fragsFile),
              chrom.size=as.character(opt$chromsize),
              vp.file=as.character(opt$validpair),
              prefix=as.character(opt$prefix),
              chromA=as.character(opt$chrA),
              chromB=as.character(opt$chrB),
              id=k
         )

        cis_trans_coverage(df,as.character(opt$prefix),as.character(opt$chrA),as.character(opt$chrB),id=k)
        boundary <- di_and_hmm(
              coverage=paste0(opt$prefix,"_",opt$chrA,"_",opt$chrB,"_",k,"_Coverage.bed"),
              chromA=as.character(opt$chrA),
              chromA.start=as.integer(opt$startA),
              chromA.end=as.integer(opt$endA),
              chromB=as.character(opt$chrB),
              chromB.start=as.integer(opt$startB),
              chromB.end=as.integer(opt$endB),
              prefix=opt$prefix,
              id=k
         )

         bp_optimized <- bp_optimization(
              paste0(opt$prefix,"_",opt$chrA,"_",opt$chrB,"_",k,".txt"),
              boundary,
              as.character(opt$chrA),
              as.character(opt$chrB),
              ssA=as.integer(opt$ssA),
              ssB=as.integer(opt$ssB),
              seA=as.integer(opt$seA),
              seB=as.integer(opt$seB)
         )
         chromA_BP_start[k] <- head(df$bed[df$bed$chrom==as.character(opt$chrA) & df$bed$start >= floor(bp_optimized$par[1]),]$start,1)
         chromA_BP_end[k]   <- head(df$bed[df$bed$chrom==as.character(opt$chrA) & df$bed$start >= floor(bp_optimized$par[1]),]$end,1)
         chromB_BP_start[k] <- head(df$bed[df$bed$chrom==as.character(opt$chrB) & df$bed$start >= floor(bp_optimized$par[2]),]$start,1)
         chromB_BP_end[k]   <- head(df$bed[df$bed$chrom==as.character(opt$chrB) & df$bed$start >= floor(bp_optimized$par[2]),]$end,1)
         k = k+1
      }
      chromA_BP_start <- unlist(chromA_BP_start)
      chromA_BP_end   <- unlist(chromA_BP_end)
      chromB_BP_start <- unlist(chromB_BP_start)
      chromB_BP_end   <- unlist(chromB_BP_end)
      param <- cbind(param,chromA_BP_start,chromA_BP_end,chromB_BP_start,chromB_BP_end)
      write.table(param,file=paste0(opt$prefix,"_",opt$chrA,"_",opt$chrB,".RE_BreakPoints.txt"),col.names=T,row.names=F,sep="\t",quote=F)
    }
  } else {
    print ("Breakpoints have not enough resolution support!")
  }
} else {
  print ("No enriched translocation box found")
}

## Orginizing results ##
system("mkdir Lower_Resolution_HiC_Data", wait=T)
files_matrix <- list.files(".", pattern=paste0(chromA,"_",chromB,".matrix"))
files_bed    <- list.files(".", pattern=paste0(chromA,"_",chromB,"_abs.bed"))
for(i in 1:length(files_matrix)) {
  system(paste0("mv ",files_matrix[i]," Lower_Resolution_HiC_Data/"),wait=T)
  system(paste0("mv ",files_bed[i]," Lower_Resolution_HiC_Data/"),wait=T)
}

system("mkdir Translocations", wait=T)
system("mkdir Translocations/juicebox_files", wait=T)
system("mkdir Translocations/Details", wait=T)
f <- list.files(".",pattern="preCluster.txt")
if (file.exists(f[1])) {
  system("mv *.preCluster.txt Translocations/")
  system("mv *.Details.txt Translocations/Details/")
}

f <- list.files(".",pattern="Translocations_jcbx.txt")
if (length(f) > 0) {
  system("mv *.Translocations_jcbx.txt Translocations/juicebox_files/")
}

f <- list.files(".",pattern="MultiResolution_Filtered.Translocation.txt")
if (file.exists(f[1])) {
  system("mkdir MultiResolution_supported_Translocations")
  system("mv *.MultiResolution_Filtered.Translocation.txt  MultiResolution_supported_Translocations/")
}

f <- list.files(".",pattern="_REfrags")
if (length(f) > 0) {
  system("mkdir RE_Level_Translocation")
  system(paste0("mv ",prefix,"_",chromA,"_",chromB,"* RE_Level_Translocation/"))
  system(paste0("mv *.hmm.states.bed RE_Level_Translocation/"))
 
}

print(paste0("Ended ",Sys.time()))

sink()
