# Authors:              Shang Liu
# Contact information:  liushang@genomics.cn
# Date:                 2019-12-16
# Version:              1.0

######################## Overview ########################

# Description: merge.sj.r, a script used to detect junction satisfy an expression threshold of read count > min.read and an cell th#		  reshold of cell count > min.cell 
#

##########################################################

library(optparse)
##### get options

option_list = list(
  make_option(c("--sj_dir"), type="character", default=NULL,metavar="character",
              help="The path of splicing junction outputed by STAR"),
  make_option(c("--outdir"), type="character", default=NULL,metavar="character",
              help="The pathway of outdir"),
  make_option(c("--min_read"), type="numeric", default=4,metavar="numeric",
              help="The read count threshold"),
  make_option(c("--min_cell"), type="numeric", default=10,metavar="numeric",
              help="The cell count threshold"),
  make_option(c("--cpu"), type="numeric", default=6,metavar="numeric",
              help="The cpu number")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(ggplot2)
library(parallel)
if (!is.null(opt$sj_dir)){
	sjdir <- opt$sj_dir
}else{
	print_help(opt_parser)
	stop("The path of splicing junction outputed by STAR should be supplied", call.=FALSE)
}

if (!is.null(opt$outdir)){
	outdir <- opt$outdir
}else{
	print_help(opt_parser)
	stop("The path of output dir by STAR should be supplied", call.=FALSE)
}

cl <- makeCluster(opt$cpu,type="FORK")

junctions <- NULL
cell<- list.files(sjdir)

## merge junction information and the reads number
junctions = parLapply(cl,cell,function(x){
	filename <- paste(sjdir,x,sep = "/")
	filename <- paste(filename,x,sep = "/")
	filename <- paste(filename,'SJ.out.tab',sep = "")
	mat <- read.table(filename,sep = '\t')
	mat$cell <- x
	id <- paste(mat$V1,mat$V2,sep = "_")
	id <- paste(id,mat$V3,sep = "_")
	id <- paste(id,mat$V4,sep = "_")
	tmp <- data.frame(junction_id = id,readnum = mat$V7)
	return(tmp)
})
junctions <- do.call(rbind,junctions)

filename1 <- paste(outdir,'Alljunction.list.reads.number.xls',sep = '/')
write.table(junctions,file=filename1,sep = '\t',quote=F,row.names=F,col.names=F)

## filter junctions of which the reads number is less than 4
pos <- which(junctions$readnum >= 11)
readstmp <- junctions
readstmp$readnum[pos] <- '>= 11'

freqnum <- table(as.character(readstmp$readnum))
freqnum <- data.frame(freqnum)
freqnum <- freqnum[order(freqnum$Var1),]

filename2 <- paste(outdir,'Alljunction.reads.number.distribution.pdf',sep = '/')
pdf(filename2)
ggplot(freqnum, aes(x=factor(Var1,levels=c("0","1","2","3","4","5","6","7","8","9","10",">= 11")), y=as.numeric(Freq))) +
  geom_bar(fill="#4DAF4A",alpha = .9, stat="identity",width=0.8) +
  geom_text(aes(x=Var1,y=Freq+20,label=Freq))+
  guides(fill=FALSE)+
  theme(legend.key = element_blank(),legend.title = element_blank())+
  xlab("Read Count")+ylab("Frequency") +ggtitle("Read count distribution")
dev.off()

pos <- which(junctions$readnum >= opt$min_read)
junctionsfilter <- junctions[pos,]

## filter junctions of which the cell numer is less than 10
freqcell <- data.frame(table(as.character(junctionsfilter$junction_id)))
freqcelltmp <- freqcell
pos <- which(as.numeric(as.character(freqcell$Freq)) > 10)
freqcell$Freq[pos] <- '>= 11'
freqcellnum <- data.frame(table(as.character(freqcell$Freq)))

filename3 <- paste(outdir,'Junctions.cell.number.distribution.pdf',sep = '/')
pdf(filename2)
ggplot(freqcellnum, aes(x=factor(Var1,levels=c("0","1","2","3","4","5","6","7","8","9","10",">= 11")), y=as.numeric(Freq))) +
  geom_bar(fill="#4DAF4A",alpha = .9, stat="identity",width=0.8) +
  geom_text(aes(x=Var1,y=Freq+20,label=Freq))+
  guides(fill=FALSE)+
  theme(legend.key = element_blank(),legend.title = element_blank()
  )+  xlab("Cell number")+ylab("Frequency") +ggtitle("Cell number distribution")
dev.off()

pos <- which(as.numeric(as.character(freqcelltmp$Freq)) > opt$min_cell)
freqcelltmpfilter <- freqcelltmp[pos,]
filename4 <- paste(outdir,'Alljunction.filter.list.xls',sep = '/')
head(freqcelltmpfilter)
freqcelltmpfilter1 <- data.frame(junction = freqcelltmpfilter$Var1)
write.table(freqcelltmpfilter1,file=filename4,sep = "\t",quote =F,row.names=F,col.names=F)

##DONE
