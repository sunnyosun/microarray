## Feature Extraction of micriarray data using Keeney's sequences##
# by Sunny #

# Install "marray" package the first time you run this function, run the following two lines:
## source("http://bioconductor.org/biocLite.R")
## biocLite("marray")

# Before running the function, make sure to put the "SK1rosetta_Keeney.txt" file into your working 
# directory, the "SK1rosetta_Keeney.txt" file can be found in WD_Drive/R_files/FE folder.

# The four variables in the FE function are FileName, cy0 and strain:
## 1. "FileName" is the name of your .txt file outputed by Agilent Feature Extractor. 
## You can use Tab to find them in your working directory.
## 2. "cy0" is the cyanine dye you use for the control (t=0). e.g. if you use cy5, just put the
##     number 5.
## 3. "strain" is the name of your strain, which will appear as the title of your plot.
## 4. "range" is the y axis limit of your plots. If not given, it runs as the maximum value of your 
##     signal.

# e.g. If you want to look at file "US93503718_251523910299_S01_CGH_107_Sep09_2_Hochwagen_1_4.txt",
#      and the dye for the control is cy5, strain name is "WT", and the y axis limit is 10.
#   Your FE function will work as: 
#   step1: source the function.
#   step2: type: 
#               FE("US93503718_251523910299_S01_CGH_107_Sep09_2_Hochwagen_1_1.txt", 5, "WT", 10)
#          into the console
#   step3: press return.

# The output of this function will be:
# 1. A 16 ".pdf" files of your microarray signal plot along chromosomes will be automatically saved
#    in your working directory.
# 2. The plots will also appear on new windows after you runing the function.

# The FE_PAIR function will output pairwise plots between two strains you want to compare. You
# will need to give information of two strains.

# If you want all the signal data of your array. Please run FE2. The three variables in FE2 are 
# the same as FE. A ".txt" file will be saved in your working directory with 4 columns: chr, start, 
# end, Log2Ratio.

# e.g. If you want to look at file "US93503718_251523910299_S01_CGH_107_Sep09_2_Hochwagen_1_4.txt",
#      and the dye for the control is cy5, strain name is "WT".
#      your FE2 function will work as: 
#   step1: source the function.
#   step2: type: 
#               FE2("US93503718_251523910299_S01_CGH_107_Sep09_2_Hochwagen_1_1.txt", 5, "WT")
#          into the console
#   step3: press return.


# function FE.
FE<-function(FileName, cy0, strain, range=NULL)
{
  plat<-read.table("SK1rosetta_Keeney.txt", header=T, sep="\t")
  # read the platform information into your workspace
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
  # chrorder is a function to select information for specific chromosomes, and ordered by
  # the start positions along chromosome. Use as chrorder(x), x is the chromosome number.
  library(marray)
  maData<-read.Agilent(fnames=FileName, name.Rf="rMeanSignal", name.Gf="gMeanSignal",
                       name.Rb="rBGMeanSignal", name.Gb="gBGMeanSignal", sep="\t")
  # marray is the package to extract the signal data.
  maNorm<-maNorm(maData,norm="printTipLoess")
  # normalize the data by "printTipLoess" method.
  if (cy0==5)
  {
    lr=-maNorm@maM
  }
  if (cy0==3)
  {
    lr=maNorm@maM
  }
  # Different read-out depends on what control dye you use.
  if (is.null(range)==FALSE)
  {
    for (i in 1:16)
    {
      pdf(file=paste(strain,'_chr', i, '.pdf', sep=""))
      plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr[chr(i)][chrorder(i)]),
           xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
           ylab='Signal', frame.plot=F,pch=16, col='blue', main=strain, ylim=c(0,range))
      dev.off()
    }
  }
  else
  {
  for (i in 1:16)
  {
    pdf(file=paste(strain,'_chr', i, '.pdf', sep=""))
    plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr[chr(i)][chrorder(i)]),
         xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
         ylab='Signal', frame.plot=F,pch=16, col='blue', main=strain)
    dev.off()
  }
  }
  # save the plots into a pdf file into your working directory.
}


# function FE_PAIR.
FE_PAIR<-function(FileName_1, cy0_1, strain_1, FileName_2, cy0_2, strain_2, range=NULL)
{
  plat<-read.table("SK1rosetta_Keeney.txt", header=T, sep="\t")
  # read the platform information into your workspace
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
  # chrorder is a function to select information for specific chromosomes, and ordered by
  # the start positions along chromosome. Use as chrorder(x), x is the chromosome number.
  library(marray)
  maData1<-read.Agilent(fnames=FileName_1, name.Rf="rMeanSignal", name.Gf="gMeanSignal",
                        name.Rb="rBGMeanSignal", name.Gb="gBGMeanSignal", sep="\t")
  maData2<-read.Agilent(fnames=FileName_2, name.Rf="rMeanSignal", name.Gf="gMeanSignal",
                        name.Rb="rBGMeanSignal", name.Gb="gBGMeanSignal", sep="\t")
  # marray is the package to extract the signal data.
  maNorm1<-maNorm(maData1,norm="printTipLoess")
  maNorm2<-maNorm(maData2,norm="printTipLoess")
  # normalize the data by "printTipLoess" method.
  if (cy0_1==5)
  {
    lr1<--maNorm1@maM
  }
  if (cy0_1==3)
  {
    lr1<-maNorm1@maM
  }
  if (cy0_2==5)
  {
    lr2<--maNorm2@maM
  }
  if (cy0_2==3)
  {
    lr2<-maNorm2@maM
  }
  # Different read-out depends on what control dye you use.
  if (is.null(range)==FALSE)
  {
    for (i in 1:16)
    {
      pdf(file=paste(strain_1, '-', strain_2, '_chr', i, '.pdf', sep=""))
      par(mfrow=c(2,1))
      plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr1[chr(i)][chrorder(i)]),
           xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
           ylab='Signal', frame.plot=F,pch=16, col='blue', main=strain_1, ylim=c(0,range))
      plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr2[chr(i)][chrorder(i)]),
           xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
           ylab='Signal', frame.plot=F,pch=16, col='red', main=strain_2, ylim=c(0,range))
      dev.off()
    }
  }
  else
  {
  for (i in 1:16)
  {
    pdf(file=paste(strain_1, '-', strain_2, '_chr', i, '.pdf', sep=""))
    par(mfrow=c(2,1))
    plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr1[chr(i)][chrorder(i)]),
         xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
         ylab='Signal', frame.plot=F,pch=16, col='blue', main=strain_1)
    plot(plat[chr(i),][chrorder(i),'start']/1000,2^(lr2[chr(i)][chrorder(i)]),
         xlab=paste('Chromosome',i,' Position (kb)', sep=""), 
         ylab='Signal', frame.plot=F,pch=16, col='red', main=strain_2)
    dev.off()
  }
  }
  # save the plots into a pdf file into your working directory.
}




# function FE2.
FE2<-function(FileName, cy0, strain)
{
  library(marray) 
  plat<-read.table("SK1rosetta_Keeney.txt", header=T, sep="\t")
  maData<-read.Agilent(fnames=FileName, name.Rf="rMeanSignal", name.Gf="gMeanSignal",
                       name.Rb="rBGMeanSignal", name.Gb="gBGMeanSignal", sep="\t")
  maNorm<-maNorm(maData,norm="printTipLoess")
  if (cy0==5)
  {
    lr<--maNorm@maM
  }
  if (cy0==3)
  {
    lr<-maNorm@maM
  }
  data<-as.data.frame(matrix(0,nrow=nrow(plat),ncol=4))
  data[,1]<-plat[,'chr']
  data[,2]<-plat[,'start']
  data[,3]<-plat[,'stop']
  data[,4]<-lr[,1]
  colnames(data)<-c('chr','start', 'end', "Log2Ratio")
  write.table(data,strain,col.names=T, row.names=F, quote=F, sep="\t")
  # save the signal data into a txt file in your working directory.
}


