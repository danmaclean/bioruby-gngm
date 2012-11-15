# functions for dealing with generic chastity belt formation and manipulation
# Revision: 0.5.0


library("RColorBrewer")
#library("hopach")   # only called when k=NA

filterMaqSNPs <- function(snps,minDepth=6,maxDepth=250,neighQ=40,bestQ=20,conQ=20)
  {
    snps <- snps[snps$V6 >= minDepth & snps$V6 <= maxDepth & snps$V9 >= neighQ & snps$V8 >= bestQ & snps$V5 >= conQ,]
    return(snps)
  }

filterBwaSNPs <- function(snps,minDepth=6,maxDepth=250,neighQ=40,bestQ=20,conQ=20)
  {
    snps <- snps[snps$V6 >= minDepth & snps$V6 <= maxDepth & snps$V9 >= neighQ & snps$V8 >= bestQ & snps$V7 >= conQ,]
    return(snps)
  }

# produce a plot of chastity belts and return info on partitions from a set of ChD annotated SNPs
ChBand <- function(snps,k,mutant,adj,win,slide,chr,plot=TRUE)   # assume snps from a single chromosome here
  {
    snps <- data.frame(snps[,2],snps[,dim(snps)[2]])  # positions are 2nd column and chastityNR last column
    names(snps) <- c("pos","Cnr")
    a <- densWindow(snps,win,slide,adj)  # snps, window, slide
    b <- getY(a,pch="y")
    rn <- rownames(b)
    b <- apply(b,2,as.numeric)
    rownames(b) <- rn
    c <- cor(t(b))
    if (is.na(k))
      {
        library(hopach)
        k <- silcheck(c,diss=TRUE)[1]
      }
    d <- kmeans(c,k,nstart=1000)
    e <- factor(as.numeric(d$cluster))
    clrs <- e
    levels(clrs) <- brewer.pal(k,"Paired")
    clrs <- as.character(clrs)
    if (plot)
      {
        for (x in 1:dim(b)[1])
          {
         if (x == 1) {plot(a[[1]],b[x,],t="l",col=clrs[x],xlab="Position (bp)",ylab="Density",main=paste(mutant,"       All phases: k = ",k),xaxp=c(1,chrLengths[chr],3),xlim=c(1,chrLengths[chr]),cex.main=2,cex.lab=2,cex.axis=2)  }
            else        {   lines(a[[1]],b[x,],t="l",col=clrs[x])     }
          }
      }
    return(list(a,b,c,e,clrs))
  }

# return list of density curves over sliding window of chastity values - helper to ChBands
densWindow <- function(snps,winSize,slideSize,adj)
  {
    result <- list()
    for (x in seq(0.2,1,by=slideSize))
      {
        maxpos <- max(snps$pos)
        hits <- snps$pos[snps$Cnr >= x & snps$Cnr < x + winSize]
        if (length(hits) > 0)
          {
            tmp <- density(hits, n=240, from=0, to=maxpos,kernel="gaussian",adjust=adj)
            result <- c(result,tmp,Cnr=x)
          }
      }
    return(result)
  }

# get a list of yvalues for the density curves retrieved by interval_split - helper to ChBands
getY <- function(mat,pch="y")
  {
  u <- which(names(mat) == pch)
  v <- which(names(mat) == "Cnr")
  if (length(u) != length(v))
    {
      stop(paste("error: x values and Cnr values unequal in getY():",as.character(u), " != ",as.character(v)))
    }
  cvals <- list()
  for (x in v)
    {
      cvals <- c(cvals,mat[[x]])
    }
  cvals <- as.character(cvals)
  result <- list()
  for (x in u)
    {
      result <- rbind(result,mat[[x]])
    }
  rownames(result) <- cvals
  return(result)
  }

# produce a plot of the two chastity belts that encompass the mono and biphasic signals
ChSplits <- function(ChList,mono,bi,chromosome)
{
  a <- ChList[[1]]
  b <- ChList[[2]]
  c <- ChList[[3]]
  e <- ChList[[4]]
  clrs <- ChList[[5]]
  # next plot splitting clusters to compare
  botNum <- e[which(as.numeric(rownames(b)) == bi)]   # the biphasic
  bottom <- b[which(e == botNum),]    # cluster under the top curve
  bc <- clrs[which(e == botNum)]
  br <- c(min(as.numeric(rownames(bottom))),max(as.numeric(rownames(bottom))))
  topNum <- e[which(as.numeric(rownames(b)) == mono)]  # the monophasic
  top <- b[which(e == topNum),]
  tc <- clrs[which(e == topNum)]
  tr <- c(min(as.numeric(rownames(top))),max(as.numeric(rownames(top))))
  newClrs <- c(bc,tc)
  # get ratio values
  ratios <- rbind(bottom,top)
#  lab <- paste("Homozygous: ",as.character(tr[1])," < ChD < ",as.character(tr[2]),"   ","Heterozygous: ",as.character(br[1])," < ChD < ", as.character(br[2]))
  lab <- "Homozygous and Heterozygous Chastity Belts"
  for (x in 1:dim(ratios)[1])
    {
      if (x == 1)   {   plot(a[[1]],ratios[x,],t="l",col=newClrs[x],xlab="Position (bp)",ylab="Density",main=lab,cex.main=2,cex.lab=2,cex.axis=2,xaxp=c(1,chrLengths[chromosome],3),xlim=c(1,chrLengths[chromosome]),ylim=c(0,max(max(ratios[[1]]),max(ratios[[2]])) * 2))      }
      else          {   lines(a[[1]],ratios[x,],t="l",col=newClrs[x])  }
    }
}

# plot the ratio of the two clusters encompassing the mono and biphasic belts
ChRatio <- function(ChList,title,mono,bi,chr)
  {
    a <- ChList[[1]]
    b <- ChList[[2]]
    c <- ChList[[3]]
    e <- ChList[[4]]
    clrs <- ChList[[5]]
    # next plot splitting clusters to compare
    botNum <- e[which(as.numeric(rownames(b)) == bi)]   # the biphasic
    bottom <- b[which(e == botNum),]    # cluster under the top curve
    bc <- clrs[which(e == botNum)]
    br <- c(min(as.numeric(rownames(bottom))),max(as.numeric(rownames(bottom))))
    topNum <- e[which(as.numeric(rownames(b)) == mono)]  # the monophasic
    top <- b[which(e == topNum),]
    tc <- clrs[which(e == topNum)]
    tr <- c(min(as.numeric(rownames(top))),max(as.numeric(rownames(top))))
    newClrs <- c(bc,tc)
    # plot ratio values
    ratios <- rbind(bottom,top)
    tp <- apply(top,2,function(x) mean(as.numeric(x),trim=0.05))
    bt <- apply(bottom,2,function(x) mean(as.numeric(x),trim=0.05))
    xvals <- a[[1]]
    lab1 <- paste("Homozygous: ",as.character(tr[1])," < ChD < ",as.character(tr[2]))
    lab2 <- paste("Heterozygous: ",as.character(br[1])," < ChD < ", as.character(br[2]))
    plot(xvals,tp/bt,t="l",main=title,xlab="Position (bp)",ylab="Ratio",cex.main=2,cex.axis=2,cex.lab=2,xaxp=c(1,chrLengths[chr],3),xlim=c(1,chrLengths[chr]))
#    plot(xvals,tp/bt,t="l",main=title,xlab="Position (bp)",ylab="Ratio",cex.main=1,cex.axis=0.5,cex.lab=0.5,xaxp=c(1,chrLengths[chromosome],3),xlim=c(1,chrLengths[chromosome]))
    mtext(lab1, cex=1.25, line=-2, at=250000,adj=0)
    mtext(lab2, cex=1.25, line=-4, at=250000,adj=0)
#    mtext(lab1, cex=0.5, line=-2, at=250000,adj=0)
#    mtext(lab2, cex=0.5, line=-4, at=250000,adj=0)
    return(list(tp,bt,xvals))
  }



## Functions for returning and plotting the candidate mutation region
# add a shading box to suspected mutant area
addBox <- function(a)
  {
    a1 <- a[[1]] / a[[2]]
    x <- which(a1 > (range(a1)[2] - range(a1)[1]) / 4)
    y <- a1[x]
    x2 <- a[[3]][x]
    polygon(c(x2[1],x2,x2[length(x)]),c(-1,y,-1),density=25,col="grey20")
  }
# return the region one sd above the mean of the ratio amplitude
getRegion <- function(a)
  {
    a1 <- a[[1]] / a[[2]]
    x <- which(a1 > (range(a1)[2] - range(a1)[1]) / 4)
    x2 <- a[[3]][x]
    return(list(min=min(x2),max=max(x2)))
  }



# plot the ratio of the two clusters encompassing the mono and biphasic belts
ChRatioNoPlot <- function(ChList,mono,bi)
  {
    a <- ChList[[1]]
    b <- ChList[[2]]
    c <- ChList[[3]]
    e <- ChList[[4]]
    clrs <- ChList[[5]]
    botNum <- e[which(as.numeric(rownames(b)) == bi)]   # the biphasic
    bottom <- b[which(e == botNum),]    # cluster under the top curve
    bc <- clrs[which(e == botNum)]
    br <- c(min(as.numeric(rownames(bottom))),max(as.numeric(rownames(bottom))))
    topNum <- e[which(as.numeric(rownames(b)) == mono)]  # the monophasic
    top <- b[which(e == topNum),]
    tc <- clrs[which(e == topNum)]
    tr <- c(min(as.numeric(rownames(top))),max(as.numeric(rownames(top))))
    newClrs <- c(bc,tc)
    ratios <- rbind(bottom,top)
    tp <- apply(top,2,function(x) mean(as.numeric(x),trim=0.05))
    bt <- apply(bottom,2,function(x) mean(as.numeric(x),trim=0.05))
    xvals <- a[[1]]
    rat <- tp/bt
    return(list(rat,xvals))
  }


# plot a  descriptive plot detailing the ratio of actual SNP histograms and frequencies
hister <- function(snps,chr,homo,hetLow,hetHi)
  {
    snpF <- snps[snps$V1 == chr,]
    snpHomo <- snpF[snpF$V5 > homo,]
    snpHet <- snpF[snpF$V5 > hetLow && snpF$V5 < hetHi,]
    par(mfrow=c(3,1))
    homo <- hist(snpHomo$V2,breaks=100)
    het <- hist(snpHet$V2,breaks=100)
    plot(het$counts/homo$counts,t="l",col="red")
  }
