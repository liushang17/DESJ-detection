# Authors:              Shang Liu
# Contact information:  liushang@genomics.cn
# Date:                 2019-12-16
# Version:		1.0

######################## Overview ########################

# Description: Plot.cell.filter.r, a novel computation framework, is designed to address the main barriersin differential splicing
#              between groups at single cell levels, such as technical noise and low coverage depth. It applies iterative K-means
#              and a new normalization method to relieve the above diffculty. 
#

##########################################################
library(optparse)
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
              help="The cpu number"),
  make_option(c("--genename"), type="character", default=NULL,metavar="character",
              help="The name of gene is used to show how to filter cells")
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

## find the gene
ann1 <- ann
frq <- data.frame(table(as.character(ann1$V3)))
pos <- which(frq$Freq >= 2)
frq1 <- frq[pos,]
geneall <- unique(as.character(frq1$Var1))
pos <- which(geneall %in% "")
if(length(pos) == 1){geneall <- geneall[-1]}

##plot the gene
  filterpdf <- paste(outdir,'filter.pdf',sep = '/')
  system(paste("mkdir -p ",filterpdf))

  x <- opt$genename

  pos <- which(ann$V3 %in% x)
  anntmp <- ann[pos,]
  pos <- which(rownames(mat) %in% anntmp$V1)
  mattmp <- mat[pos,]

  numsum <- colSums(mattmp)
  pos <- which(numsum == 0)
  numsum[pos] <- 1

  mattmp2 <- data.frame(t(mattmp))
  mattmp2 <- mattmp2[1:length(rownames(mattmp2)),]
  colnames(mattmp2) <- rownames(mattmp)
  mattmp2$allcell <- rownames(mattmp2)
  mattmp2$allstat <- numsum
  mattmp3 <- merge(mattmp2,cellstat,by="allcell")
  num <- length(colnames(mattmp2))
  for(i1 in 2:num){
    mattmp3[,i1] <- as.numeric(as.character(mattmp3[,i1])) / as.numeric(as.character(mattmp3$uniq)) * 1000000
  }

  mattmp4 <- mattmp3[,1:(num-1)]
  mattmp5 <- mattmp4[,2:(num-1)]
  mattmp6 <- data.frame(t(mattmp5))
  colnames(mattmp6) <- mattmp4$allcell
  rownames(mattmp6) <- colnames(mattmp5)

  mattmp <- mattmp6
  pos <- which(colnames(mattmp) %in% colann$allcell)
  mattmp <- mattmp[,pos]

  mattmp_test = rowSums(mattmp)
  pos <- which(mattmp_test > 0)
  if(length(pos) >= 1){
   tes <- filter(mattmp,mincell,maxsd,maxmean)
   sui <- tes
   numfilter <- sui$cl
   cluallinfo <- sui$cluinfo
   cluallinfo$clu = 'filter'
   pos <- which(cluallinfo$suitmp.cluster != numfilter)
   cluallinfo$clu[pos] <- cluallinfo$suitmp.cluster[pos]
   mit <- cluallinfo
   mit1 <- mit[order(mit$suitmp.cluster),]
   mit2 <- data.frame(clu = mit1$clu)
   rownames(mit2) <- rownames(mit1)
   mattmp1 <- mattmp[,rownames(mit2)]
   filename <- paste(filterpdf,x,sep = '/')
   filename <- paste(filename,'.pdf',sep = '')
   pdf(filename)
   pheatmap(log2(as.matrix(mattmp1)+1),annotation_col=mit2,cluster_cols = F,show_colnames = F,cluster_rows = F)
   dev.off()
  }
