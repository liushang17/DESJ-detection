# Authors:              Shang Liu
# Contact information:  liushang@genomics.cn
# Date:                 2019-12-16
# Version:		1.0

######################## Overview ########################

# Description: DESJ-detection, a novel computation framework, is designed to address the main barriersin differential splicing
#              between groups at single cell levels, such as technical noise and low coverage depth. It applies iterative K-means
#              and a new normalization method to relieve the above diffculty. 
#

##########################################################
library(optparse)
library(edgeR)
library(limma)
library(parallel)
library(pheatmap)
source('./filter.cell.r')

##### get options

option_list = list(
  make_option(c("--junction_matrix"), type="character", default=NULL,metavar="character",
              help="The junction_cell count matrix"),
  make_option(c("--col_ann"), type="character", default=NULL,metavar="character",
	      help="The text file, containing the clustering info, possess two columns, including cell name and clusters name"),
  make_option(c("--min_pct"), type="numeric", default=0.2,metavar="numeric",
              help="The minimum scale in either of the two groups of cells, is able to detect the junctions of any genes"),
  make_option(c("--cellinfo"), type="character", default=NULL,metavar="character",
              help="The text file, containing the uniquely mapped reads number, possess two colums, including cell name and
		reads number"),
  make_option(c("--row_ann"), type="character", default=NULL,metavar="character",
              help="The text file, containing the junctions info, including three columns, including junction id, gene id 
		and gene symbol"),
  make_option(c("--clus"), type="character", default=NULL,metavar="character",
              help="The name of the two clusters, will be detected for differential junctions. for example: clusterA,clusterB"),
  make_option(c("--outdir"), type="character", default=NULL,metavar="character",
              help="The pathway of the outdir"),
  make_option(c("--mincell"), type="numeric", default=10,metavar="numeric",
              help="The minimum cell number of the cluster determined by iterative k-means would be filtered"),
  make_option(c("--maxsd"), type="numeric", default=0.2,metavar="numeric",
              help="The maximum standard deviation of the expression of the junctions of the cells in the 
		filtered cluster determined by iterative k-means"),
  make_option(c("--maxmean"), type="numeric", default=0.2,metavar="numeric",
              help="The maximum mean of the expression of the junctions of the cells in the
		filtered cluster determined by iterative k-means"),
  make_option(c("--cpu"), type="numeric", default=6,metavar="numeric",
              help="The cpu number")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$junction_matrix)){
	print_help(opt_parser)
	stop("The junction count matrix must be supplied", call.=FALSE)	
}
if (is.null(opt$col_ann)){
	print_help(opt_parser)
	stop("The cluster information must be supplied", call.=FALSE)
}
if (is.null(opt$cellinfo)){
	print_help(opt_parser)
	stop("The uniquely mapped reads number information must be supplied", call.=FALSE)
}
if (is.null(opt$row_ann)){
	print_help(opt_parser)
	stop("The junction information must be supplied", call.=FALSE)
}
if (is.null(opt$clus)){
	print_help(opt_parser)
	stop("The name of two clusters must be supplied", call.=FALSE)
}
if (is.null(opt$outdir)){
	print_help(opt_parser)
	stop("The path of outdir must be supplied", call.=FALSE)
}

## loading data
mat <- read.table(opt$junction_matrix,sep = '\t',header= T,row.names=1)
colann <- read.table(opt$col_ann,sep = '\t',header = T)
min.pct <- as.numeric(opt$min_pct)
cellstat <- read.table(opt$cellinfo,sep = '\t',header = T)
ann <- read.table(opt$row_ann,sep = '\t')
clus <- opt$clus
outdir <- opt$outdir
mincell <- as.numeric(opt$mincell)
maxsd <- as.numeric(opt$maxsd)
maxmean <- as.numeric(opt$maxmean)

## annotate the data
mat <- data.frame(t(mat))
colnames(colann) <- c('allcell','majorCluster')
clus <- as.vector(strsplit(clus,',')[[1]])
pos <- which(colann$majorCluster %in% clus)
colann <- colann[pos,]
pos <- which(colann$majorCluster %in% clus[1])
clu1.mincell <- length(pos) * min.pct
pos <- which(colann$majorCluster %in% clus[2])
clu2.mincell <- length(pos) * min.pct

clu.mincell <- clu1.mincell
if(clu1.mincell > clu2.mincell){clu.mincell <- clu2.mincell}

system(paste("mkdir -p",outdir,sep = " "))
## find the gene
ann1 <- ann
frq <- data.frame(table(as.character(ann1$V3)))
pos <- which(frq$Freq >= 2)
frq1 <- frq[pos,]
geneall <- unique(as.character(frq1$Var1))
pos <- which(geneall %in% "")
if(length(pos) == 1){geneall <- geneall[-1]}

## filter the low expression cell and calculate the logFC and adjust.p.value
geneall <- geneall
resall <- NULL
cl <- makeCluster(opt$cpu,type="FORK")
resall = parLapply(cl,geneall,function(x){
  pos <- which(ann$V3 %in% x)
  anntmp <- ann[pos,]
  pos <- which(rownames(mat) %in% anntmp$V1)
  mattmp <- mat[pos,]

  tempOutput <- NULL
  mattmp2 <- data.frame(t(mattmp))
  mattmp2 <- mattmp2[1:length(rownames(mattmp2)),]
  colnames(mattmp2) <- rownames(mattmp)
  mattmp2$allcell <- rownames(mattmp2)
  mattmp3 <- merge(mattmp2,cellstat,by="allcell")
  num <- length(colnames(mattmp2))
  for(i1 in 2:num){
    mattmp3[,i1] <- as.numeric(as.character(mattmp3[,i1])) / as.numeric(as.character(mattmp3$uniq)) * 1000000
  }

  mattmp4 <- mattmp3[,1:num]
  mattmp5 <- mattmp4[,2:num]
  mattmp6 <- data.frame(t(mattmp5))
  colnames(mattmp6) <- mattmp4$allcell
  rownames(mattmp6) <- colnames(mattmp5)

  mittmp <- mattmp
  mattmp <- mattmp6
  pos <- which(colnames(mattmp) %in% colann$allcell)
  mattmp <- mattmp[,pos]
  
  mattmp_test = rowSums(mattmp)
  pos <- which(mattmp_test > 0)
  
  if(length(pos) > 0){
  
	  sui <- filter(mattmp,mincell,maxsd,maxmean)
	  numfilter <- sui$cl
	  cluallinfo <- sui$cluinfo
	  pos <- which(cluallinfo$suitmp.cluster == numfilter)
	  clufilter <- cluallinfo[-pos,]
	  pos <- which(colnames(mattmp) %in% rownames(clufilter))
	  mattmp <- mattmp[,pos]

	  if(length(pos) > clu.mincell){
		  mattmp <- mattmp[,colSums(mattmp)!=0]
		  pos <- which(colnames(mittmp) %in% colnames(mattmp))
		  mittmp <- mittmp[,pos]

		  numsum <- colSums(mittmp)
		  
		  mittmp2 <- data.frame(t(mittmp))
		  mittmp2 <- mittmp2[1:length(rownames(mittmp2)),]
		  colnames(mittmp2) <- rownames(mittmp)
		  mittmp2$allcell <- rownames(mittmp2)
		  mittmp2$allstat <- numsum
		  mittmp3 <- merge(mittmp2,cellstat,by="allcell")
		  num <- length(colnames(mittmp2))
		  for(i1 in 2:(num-1)){
			mittmp3[,i1] <- as.numeric(as.character(mittmp3[,i1])) / as.numeric(as.character(mittmp3$allstat)) * 100
		  }
		  
		  mittmp4 <- mittmp3[,1:(num-1)]
		  mittmp5 <- mittmp4[,2:(num-1)]
		  mittmp6 <- data.frame(t(mittmp5))
		  colnames(mittmp6) <- mittmp4$allcell
		  rownames(mittmp6) <- colnames(mittmp5)
		  
		  
		  mittmp <- mittmp6
		  mattmp <- mittmp
		  
		  juncnum <- length(rownames(mattmp))
		  clu_used <- clus[1]
		  pos <- which(colann$majorCluster %in% clu_used)
		  pos <- which(colnames(mattmp) %in% colann$allcell[pos])
		  numused <- length(pos)
		  mat_used <- mattmp[,pos]

		  pos <- which(colann$majorCluster %in% clus[2])
		  pos <- which(colnames(mattmp) %in% colann$allcell[pos])
		  numfilter <- length(pos)
		  mat_filter = mattmp[,pos]
		  if(numused > clu.mincell & numfilter > clu.mincell){
			mat_used=mat_used[,colSums(mat_used)!=0]
			mat_filter = mat_filter[,colSums(mat_filter)!=0]

			matexnex=cbind(mat_filter,mat_used)
			matexnex=matexnex[rowSums(matexnex)!=0,]
			keep <- rowSums(matexnex>0) >= length(colnames(mat_used)) * 0.01
			mat_used<-mat_used[keep,]
			mat_filter <- mat_filter[keep,]
		
			numused2 <- length(which(colSums(mat_used)!=0))
			mat_used=mat_used[,colSums(mat_used)!=0]
			numfilter2 <- length(which(colSums(mat_filter)!=0))
			mat_filter = mat_filter[,colSums(mat_filter)!=0]
			if(numused2 > clu.mincell & numfilter2 > clu.mincell){
			  matexnex=cbind(mat_filter,mat_used)
			
			  alll=length(colnames(matexnex))
			  numrow = length(rownames(matexnex))
			  if(length(rownames(matexnex)) < 5){
				 for(nrowi in (numrow+1):5){
				 matexnex[nrowi,] = 100
			    }
			  matexnex[6,] = 10000
			  }
			  dge <- DGEList(counts=matexnex)
			  group1=1
			  for (i3 in 1:length(colnames(mat_filter))){group1[i3]="EX"}
			  for (j3 in (length(colnames(mat_filter))+1):(length(colnames(matexnex)))){group1[j3]="NEX"}
			  group<-factor(group1)
			  design <- model.matrix(~group)
			  logCPM <- cpm(dge, log=TRUE, prior.count=1)
			  fit <- lmFit(logCPM, design)
			  fit <- eBayes(fit, trend=TRUE)
			  tempOutput = topTable(fit, coef=ncol(design))

			  tempOutput = data.frame(tempOutput)
			  if(numrow < 5){
			    numtmp <- (numrow+1):6
			    pos <- which(rownames(tempOutput) %in% numtmp)
			    if(length(pos) >= 1){tempOutput <- tempOutput[-pos,]}
			  }
			
			  pos <- which(rownames(matexnex) %in% rownames(tempOutput))
			  matexnex1 <- matexnex[pos,]
			  kares <- matrix(nrow = length(rownames(matexnex1)),ncol = 4)
			  junc <- rownames(tempOutput)
			  for(i2 in 1:length(junc)){
			    chr_used <- junc[i2]
			    pos <- which(rownames(mat_filter) %in% chr_used)
			    tmp_filter <- as.numeric(as.character(as.matrix(mat_filter[pos,])))
			    tmp_filter1 <- length(which(tmp_filter > 0))
			    tmp_filter0 <- length(which(tmp_filter == 0))
			  
			    pos <- which(rownames(mat_used) %in% chr_used)
			    tmp_used <- as.numeric(as.character(as.matrix(mat_used[pos,])))
			  
			    mean_used <- mean(log2(tmp_used+1))
			    mean_filter <- mean(log2(tmp_filter +1))
			    tmp_used1 <- length(which(tmp_used > 0))
			    tmp_used0 <- length(which(tmp_used == 0))
			    test_all <- matrix(c(tmp_filter1,tmp_filter0,tmp_used1,tmp_used0),nrow = 2,ncol = 2,byrow = T)
			    frqall <- tmp_filter1 + tmp_filter0 + tmp_used1 + tmp_used0
			    t_theory <- (tmp_filter1 + tmp_filter0) * (tmp_used1 + tmp_used0) / frqall
			    expe <- chisq.test(test_all)$expected
			    if(expe[1,1] > 5 & expe[1,2] > 5 & expe[2,1] > 5 & expe[2,2] >5){
				ts <- chisq.test(test_all)
			    }else{
				ts <- fisher.test(test_all)
			    }
			  
			    kares[i2,1] <- chr_used
			    kares[i2,2] <- ts$p.value
			    kares[i2,3] <- mean_used
			    kares[i2,4] <- mean_filter
			  }
			  tempOutput$ka.p.value <- kares[,2]
			  tempOutput$genename <- rep(x,length(rownames(tempOutput)))
			  tempOutput$used <- kares[,3]
			  tempOutput$filter <- kares[,4]
			}
		  }
	  }
    }
  return(tempOutput)
})

resall <- do.call(rbind,resall)
stopCluster(cl)

filename <- paste(outdir,'first.res.xls',sep = '/')
write.table(resall,file=filename,sep = '\t',quote=F)

pos <- which(abs(resall$logFC) > 1 & as.numeric(resall$adj.P.Val) < 0.01)
resall1 <- resall[pos,]
genefilter <- as.character(unique(resall1$genename))

outdir = paste(outdir,'res.pdf',sep = '/')
system(paste("mkdir -p",outdir,sep = " "))
write.table(resall1,file=paste(outdir,'res.xls',sep = '/'),sep = '\t',quote = F)
for(i in 1:length(genefilter)){
  x <- genefilter[i]
  pos <- which(ann$V3 %in% x)
  anntmp <- ann[pos,]
  pos <- which(rownames(mat) %in% anntmp$V1)
  mattmp <- mat[pos,]

  mattmp2 <- data.frame(t(mattmp))
  mattmp2 <- mattmp2[1:length(rownames(mattmp2)),]
  colnames(mattmp2) <- rownames(mattmp)
  mattmp2$allcell <- rownames(mattmp2)
  mattmp3 <- merge(mattmp2,cellstat,by="allcell")
  num <- length(colnames(mattmp2))
  for(i1 in 2:num){
    mattmp3[,i1] <- as.numeric(as.character(mattmp3[,i1])) / as.numeric(as.character(mattmp3$uniq)) * 1000000
  }

  mattmp4 <- mattmp3[,1:num]
  mattmp5 <- mattmp4[,2:num]
  mattmp6 <- data.frame(t(mattmp5))
  colnames(mattmp6) <- mattmp4$allcell
  rownames(mattmp6) <- colnames(mattmp5)
  
  mittmp <- mattmp
  mattmp <- mattmp6
  pos <- which(colnames(mattmp) %in% colann$allcell)
  mattmp <- mattmp[,pos]
  
  sui <- filter(mattmp,mincell,maxsd,maxmean)
  numfilter <- sui$cl
  cluallinfo <- sui$cluinfo
  pos <- which(cluallinfo$suitmp.cluster == numfilter)
  clufilter <- cluallinfo[-pos,]
  pos <- which(colnames(mattmp) %in% rownames(clufilter))
  mattmp <- mattmp[,pos]
  
  pos <- which(colnames(mittmp) %in% colnames(mattmp))
  mittmp <- mittmp[,pos]

  numsum <- colSums(mittmp)
  
  mittmp2 <- data.frame(t(mittmp))
  mittmp2 <- mittmp2[1:length(rownames(mittmp2)),]
  colnames(mittmp2) <- rownames(mittmp)
  mittmp2$allcell <- rownames(mittmp2)
  mittmp2$allstat <- numsum
  mittmp3 <- merge(mittmp2,cellstat,by="allcell")
  num <- length(colnames(mittmp2))
  for(i1 in 2:(num-1)){
    mittmp3[,i1] <- as.numeric(as.character(mittmp3[,i1])) / as.numeric(as.character(mittmp3$allstat)) * 100
  }
  
  mittmp4 <- mittmp3[,1:(num-1)]
  mittmp5 <- mittmp4[,2:(num-1)]
  mittmp6 <- data.frame(t(mittmp5))
  colnames(mittmp6) <- mittmp4$allcell
  rownames(mittmp6) <- colnames(mittmp5)
  
  
  mittmp <- mittmp6
  mattmp <- mittmp
  
  clu_used <- clus[1]
  pos <- which(colann$majorCluster %in% clu_used)
  colann_used <- colann[pos,]
  pos <- which(colnames(mattmp) %in% colann_used$allcell)
  mat_used <- mattmp[,pos]
  
  pos <- which(colann$majorCluster %in% clus[1])
  colann_filter <- colann[-pos,]
  pos <- which(colnames(mattmp) %in% colann_filter$allcell)
  mat_filter = mattmp[,pos]
  hgroup <- data.frame(clu=c(rep('used',length(colnames(mat_used))),rep('filter',length(colnames(mat_filter)))))
  mat_h <- cbind(mat_used,mat_filter)
  rownames(hgroup) <- colnames(mat_h)
  pdfname <- paste(outdir,x,sep = '/')
  pdfname <- paste(pdfname,'.pdf',sep ='')
  pdf(pdfname)
  pheatmap(log2(as.matrix(mat_h)+1),annotation_col = hgroup,cluster_cols = F,show_colnames = F,cluster_rows = F)
  dev.off()
}
