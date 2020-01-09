findfilter <- function(datatmp,maxmean,maxsd){
  suitmp <- kmeans(t(log2(as.matrix(datatmp)+1)),2)
  mittmp <- data.frame(suitmp$cluster)
  mittmp$allcell <- rownames(mittmp)
  pos <- which(mittmp$suitmp.cluster == 1)
  numusedtmp1 <- length(pos)
  numusedtmp2 <- length(mittmp$suitmp.cluster) - length(pos)
  if(numusedtmp1 == 1 | numusedtmp2 == 1){
#    choosenum = 1
    if(numusedtmp1 == 1){
      site_u <- pos
      cellclu1 <- mittmp$allcell[pos]
      pos <- which(colnames(datatmp) %in% cellclu1)
      mean1 <- datatmp[,pos]
      s1 <- sqrt(mean1)
      ss1 <- sum(s1)
      ssm1 <- sum(s1)
      
      cellclu2 <- mittmp[-site_u,]
      pos <- which(colnames(datatmp) %in% cellclu2$allcell)
      datatmp2 <- datatmp[,pos]
      site2 <- NULL
      posmedia2 <- NULL
      for(n in 1:length(rownames(datatmp2))){
        tmp <- log2(as.numeric(as.matrix(datatmp2[n,]))+1)
        sd1 <- sd(tmp)
        meidatmp <- median(tmp)
        site2 <- c(site2,sd1)
        posmedia2 <- c(posmedia2,meidatmp)
      }
      mean2 <- rowMeans(log2(as.matrix(datatmp2)+1))
      s2 <- sqrt(mean2)
      ss2 <- sum(s2)
      possd2 <- length(which(site2 >= maxsd))
      posmean2 <- length(which(mean2 >= maxmean))
      sm2 <- sqrt(posmedia2)
      ssm2 <- sum(sm2)
      choosenum <- 1
      if(ssm1 > ssm2){choosenum <- 2}else{
        if(ssm1 == ssm2){
          if(ss1 > ss2){
            choosenum <- 2
          }
        }
      }
    }else{
      site_u <- pos
      cellclu2 <- mittmp$allcell[-pos]
      pos <- which(colnames(datatmp) %in% cellclu2)
      mean2 <- datatmp[,pos]
      s2 <- sqrt(mean2)
      ss2 <- sum(s2)
      ssm2 <- sum(s2)
      
      cellclu1 <- mittmp[site_u,]
      pos <- which(colnames(datatmp) %in% cellclu1$allcell)
      datatmp1 <- datatmp[,pos]
      site1 <- NULL
      posmedia1 <- NULL
      for(m in 1:length(rownames(datatmp1))){
        tmp <- log2(as.numeric(as.matrix(datatmp1[m,]))+1)
        sd1 <- sd(tmp)
        meidatmp <- median(tmp)
        site1 <- c(site1,sd1)
        posmedia1 <- c(posmedia1,meidatmp)
      }
      mean1 <- rowMeans(log2(as.matrix(datatmp1)+1))
      s1 <- sqrt(mean1)
      ss1 <- sum(s1)
      possd1 <- length(which(site1 >= maxsd))
      posmean1 <- length(which(mean1 >= maxmean))
      sm1 <- sqrt(posmedia1)
      ssm1 <- sum(sm1)
      choosenum <- 2
      if(ssm2 > ssm1){choosenum <- 1}else{
        if(ssm1 == ssm2){
          if(ss2 > ss1){
            choosenum <- 1
          }
        }
      }
    }
    statout <- data.frame(sdstat = c(0,0),meanstat= c(0,0),cellnum = c(numusedtmp1,numusedtmp2),sq = c(0,0))
  }else{
    
    cellclu1 <- mittmp[pos,]
    cellclu2 <- mittmp[-pos,]
    pos <- which(colnames(datatmp) %in% cellclu1$allcell)
    datatmp1 <- datatmp[,pos]
    pos <- which(colnames(datatmp) %in% cellclu2$allcell)
    datatmp2 <- datatmp[,pos]
    
    site1 <- NULL
    posmedia1 <- NULL
    for(m in 1:length(rownames(datatmp1))){
      tmp <- log2(as.numeric(as.matrix(datatmp1[m,]))+1)
      sd1 <- sd(tmp)
      meidatmp <- median(tmp)
      site1 <- c(site1,sd1)
      posmedia1 <- c(posmedia1,meidatmp)
    }
    mean1 <- rowMeans(log2(as.matrix(datatmp1)+1))
    s1 <- sqrt(mean1)
    ss1 <- sum(s1)
    possd1 <- length(which(site1 >= maxsd))
    posmean1 <- length(which(mean1 >= maxmean))
    sm1 <- sqrt(posmedia1)
    ssm1 <- sum(sm1)
    
    site2 <- NULL
    posmedia2 <- NULL
    for(n in 1:length(rownames(datatmp2))){
      tmp <- log2(as.numeric(as.matrix(datatmp2[n,]))+1)
      sd1 <- sd(tmp)
      meidatmp <- median(tmp)
      site2 <- c(site2,sd1)
      posmedia2 <- c(posmedia2,meidatmp)
    }
    mean2 <- rowMeans(log2(as.matrix(datatmp2)+1))
    s2 <- sqrt(mean2)
    ss2 <- sum(s2)
    possd2 <- length(which(site2 >= maxsd))
    posmean2 <- length(which(mean2 >= maxmean))
    sm2 <- sqrt(posmedia2)
    ssm2 <- sum(sm2)
    statout <- data.frame(sdstat = c(possd1,possd2),meanstat= c(posmean1,posmean2),cellnum = c(length(colnames(datatmp1)),length(colnames(datatmp2))),sq = c(ss1,ss2))
    rownames(statout) <- c(1,2)
    choosenum <- 1
    if(ssm1 > ssm2){choosenum <- 2}else{
      if(ssm1 == ssm2){
        if(ss1 > ss2){
          choosenum <- 2
        }
      }
    }
    
  }
  tmpout <- list(clu = mittmp,statclu = statout,num_u = choosenum)
  return(tmpout)
}

filter <- function(matall,cell.min,maxsd,maxmean){
  suialltmp <- findfilter(matall,maxsd,maxmean)
  cluall <- suialltmp$clu
  clu_usedtmp <- suialltmp$clu
  statall <- suialltmp$statclu
  num_0 <- suialltmp$num_u
  num_in <- suialltmp$num_u
  sdpan <- statall$sdstat[num_0]
  meanpan <-  statall$meanstat[num_0]
  cellpan <- statall$cellnum[num_0]
  num <-  3
  sui_used_tmp <- NULL
  while ( (sdpan > 0 | meanpan > 0) & cellpan > cell.min){
    pos <- which(clu_usedtmp$suitmp.cluster == num_0)
    cellused <- rownames(clu_usedtmp)[pos]
    pos <- which(colnames(matall) %in% cellused)
    matall_used <- matall[,pos];
    
    sui_used_tmp <- findfilter(matall_used,maxsd,maxmean)
    
    num_0 <- sui_used_tmp$num_u
    clu_usedtmp <- sui_used_tmp$clu
    stat_all <- sui_used_tmp$statclu
    sdpan <- stat_all$sdstat[num_0]
    meanpan <-  stat_all$meanstat[num_0]
    cellpan <- stat_all$cellnum[num_0]
      pos <- which(clu_usedtmp$suitmp.cluster == num_0)
      cellnoused <- rownames(clu_usedtmp)[-pos]
      pos <- which(rownames(cluall) %in% cellnoused)
      cluall$suitmp.cluster[pos] <- num
      num <- num + 1
  }
  res <- list(cluinfo = cluall,cl = num_in)
  return (res)
}
