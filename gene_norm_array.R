# GENE Normalization of microarray data
# Sunny

gene_norm_m<-function(FileName, Name)
{
# load the gff file (saved as txt by excel)
out=read.table('~/Documents/LAB/data_analysis/SK1_annotation/SK1_annotation_modified.txt')

# load the microarray data run by FE2 (FE_SK1)
#Name='Red1'
#FileName='red1'
data=read.table(FileName, header=T)
plat<-read.table("SK1rosetta_Keeney.txt", header=T, sep="\t")
chr<-function(x)
{
  if (x<=9)
  {
    return (which(plat[,'chr']==paste('chr0',x,sep='')))
  }
  if (x>9)
  {
    return (which(plat[,'chr']==paste('chr',x,sep='')))
  }
}

chrorder<-function(x)
{
  return (order(plat[chr(x),'start']))
}

red1=list()
for (i in 1:16)
{
  red1[[i]]=cbind(plat[chr(i),][chrorder(i),'start'],data[,4][chr(i)][chrorder(i)])
  colnames(red1[[i]])=c('start','array')
}


# extract the GENEs + 650bp on either side
newmatrix=matrix(NA,ncol=12,nrow=1)
for (i in 1:(nrow(out)-1))
{
  if (strsplit(as.character(out[i,10]),'')[[1]][1]=='Y')
  {
    print (i)
    newline=out[i,]
    newmatrix=rbind(newmatrix,newline)
  }
}
newmatrix=newmatrix[2:nrow(newmatrix),]
mid=newmatrix


midmatrix=as.data.frame(matrix(NA,ncol=4,nrow=1))
for (i in 1:nrow(mid))
{
  newline=c(mid[i,1],as.numeric(mid[i,4])-650,as.numeric(mid[i,5])+650,mid[i,7])
  midmatrix=rbind(midmatrix,newline)
}
midmatrix=midmatrix[2:nrow(midmatrix),]


# ORF regions
midmatrix=cbind(midmatrix, 2000/(as.numeric(midmatrix[,3])-as.numeric(midmatrix[,2])))

# STRAND==+
midmatrix_w=midmatrix[which(midmatrix[,4]=='+'),]

# STRAND==-
midmatrix_c=midmatrix[which(midmatrix[,4]=='-'),]


# normalize genes on the watson strand
alldata_mid_w=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  ## strand
  chr=midmatrix_w[which(midmatrix_w[,1]==c),]
  if (nrow(chr)!=0)
  {wigmatrix=matrix(NA,ncol=2,nrow=1)
   for (i in 1:nrow(chr))
   {
     start=as.numeric(chr[i,2])
     end=as.numeric(chr[i,3])
     if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
     {
       tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,5]
       tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
       tmp=cbind(tmp1,tmp2)
       wigmatrix=rbind(wigmatrix,tmp)
     }
   }
   wigmatrix=wigmatrix[2:nrow(wigmatrix),]
   
   # mean
   wig_sort=wigmatrix[order(wigmatrix[,1]),]
   wig_sort[,1]=round(wig_sort[,1])
   mean_matrix_w=as.data.frame(matrix(NA,ncol=2,nrow=2001))
   for (i in 0:2000)
   {
     mean_matrix_w[i+1,1]=i
     mean_matrix_w[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
   }
   #plot(mean_matrix, pch=16,cex=0.6,col='blue')
   
   # smooth
   bp=5
   sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
   for (i in seq(1,2000, bp))
   {
     sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
     sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
   }
   alldata_mid_w[[k]]=sm_matrix
   #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
  }
}



# normalize genes on the crick strand
alldata_mid_c=list()
for (k in 1:16)
{
  print (k)
  if (k <= 9) {
    c<-paste("chr0",k,sep="")
  } else {
    c<-paste("chr",k,sep="")
  }
  ## strand
  chr=midmatrix_c[which(midmatrix_c[,1]==c),]
  if (nrow(chr)!=0)
  {wigmatrix=matrix(NA,ncol=2,nrow=1)
   for (i in 1:nrow(chr))
   {
     start=as.numeric(chr[i,2])
     end=as.numeric(chr[i,3])
     if ((length(which(red1[[k]][,1]>=start & red1[[k]][,1]<=end)))!=0)
     {
       tmp1=(red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),1]-start)*chr[i,5]
       tmp2=red1[[k]][which(red1[[k]][,1]>=start & red1[[k]][,1]<=end),2]
       tmp=cbind(tmp1,tmp2)
       wigmatrix=rbind(wigmatrix,tmp)
     }
   }
   wigmatrix=wigmatrix[2:nrow(wigmatrix),]
   
   # mean
   wig_sort=wigmatrix[order(wigmatrix[,1]),]
   wig_sort[,1]=round(wig_sort[,1])
   mean_matrix_c=as.data.frame(matrix(NA,ncol=2,nrow=2001))
   for (i in 0:2000)
   {
     mean_matrix_c[i+1,1]=i
     mean_matrix_c[i+1,2]=mean(wig_sort[which(wig_sort[,1]==i),2])
   }
   #plot(mean_matrix, pch=16,cex=0.6,col='blue')
   
   # smooth
   bp=5
   sm_matrix=as.data.frame(matrix(NA,ncol=2,nrow=2000/bp))
   for (i in seq(1,2000, bp))
   {
     sm_matrix[(i+bp-1)/bp,1]=i+(bp/2)
     sm_matrix[(i+bp-1)/bp,2]=mean(wig_sort[which(wig_sort[,1]>=i & wig_sort[,1]<i+bp),2])
   }
   alldata_mid_c[[k]]=sm_matrix
   #plot(sm_matrix, pch=16,col='blue',xlab=paste('Chromosome',k,' Position (kb)', sep=""), ylab='Red1 Signal')
  }
}


# Combine both strands
alldata_mid=list()
for (k in 1:16)
{
  start=alldata_mid_w[[k]][,1]
  w=alldata_mid_w[[k]][,2]
  c=alldata_mid_c[[k]][nrow(alldata_mid_c[[k]]):1,2]
  temp=as.data.frame(matrix(NA,ncol=2,nrow=1))
  for (i in 1:400)
  {
    newline=c(start[i],(as.numeric(w[i])+as.numeric(c[i]))/2)
    temp=rbind(temp,newline)
  }
  alldata_mid[[k]]=temp[2:nrow(temp),]
}

#
mean_matrix_all=cbind(mean_matrix_w, mean_matrix_c[nrow(mean_matrix_c):1,2])
mean_matrix_all=cbind(mean_matrix_all, rowMeans(mean_matrix_all[,2:3]))
mean_matrix=mean_matrix_all[,c(1,4)]

# plot each chromosome
color_s=c('red','yellow','orange','blue','purple','springgreen','brown','burlywood','pink','grey','skyblue','turquoise','dark green','violetred','magenta','gold')
color=color_s
pdf(file=paste('~/Desktop/', FileName, '-', '_gene_norm_plot.pdf', sep=""), width=10, height=8)
plot(alldata_mid[[1]][,1],alldata_mid[[1]][,2], pch=16,ylim=c(-1,5),col=color[1],
     xlab='Normalized 5\'-3\' ORF (bp)', ylab=Name)
for (i in 2:16)
{
  points(alldata_mid[[i]][,1],alldata_mid[[i]][,2],col=color[i], pch=16)
}
legend(x=0, y=4, col=color, pch=16, legend=1:16, horiz=T, cex=0.85)
dev.off()


# plot average signals
pdf(file=paste('~/Desktop/', FileName, '-', '_gene_norm_ave_plot.pdf', sep=""), width=10, height=8)
plot(mean_matrix, pch=16, cex=0.6, col='blue', xlab='Normalized 5\'-3\' ORF (bp)', ylab=Name,
     ylim=c(-1,5))
dev.off()
}