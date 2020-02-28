## ----setup, include=FALSE-----------------------------------------------------
library(selectSNPs)
knitr::opts_chunk$set(echo = TRUE)

## ----MAF, echo=FALSE----------------------------------------------------------
x<-seq(0,1,length.out=100)
y<-(-1)*(x*log2(x)+(1-x)*log2(1-x))

plot(y~x,type="p",xlab="Allele A frequency",ylab="E score (t1=1)")
abline(v=0.5)

y100<-y^100
plot(y100~x,type="p",xlab="Allele A frequency",ylab="E score (t1=100)")
abline(v=0.5)

## ----example_data_objects, echo=TRUE------------------------------------------
snp1<-new("Locus",
          Name="SNP1",
          Chromosome="1",
          Position=600000,
          Maf=0.20,
          Type="C",
          Status=as.integer(0))
str(snp1)

chr1<-new("Chrom",
          Chromosome="1",
          Name=paste("SNP",1:10,sep=""),
          Position=seq(600000,11000000,length.out=10),
          Maf=runif(10,0,0.5),
          Type=rep("C",10),
          Status=as.integer(rep(0,10)))
str(chr1)

## ----bov80K, echo=TRUE--------------------------------------------------------
data(bov80K)
summary(bov80K)
plotMap(bov80K)
u80K<-scoreU(bov80K)
e80K<-scoreE(bov80K)

## ----bov80K_continued, echo=TRUE----------------------------------------------
tmpdat<-as.data.frame(bov80K)
head(tmpdat)
tmpdat<-tmpdat[(tmpdat$Chromosome%in%c(1,10,15)),]
map1a<-as.Map(tmpdat)
str(map1a)

## ----map1, echo=TRUE----------------------------------------------------------
map1<-selectLocalOptimalSNPs(bov80K,n=2000,w1=0.5,w2=0.5)
summary(map1)
plotMap(map1)
u1<-scoreU(map1)
e1<-scoreE(map1)

## ----map2, echo=TRUE----------------------------------------------------------
map2<-selectLocalOptimalSNPs(bov80K,n=2000,w1=0,w2=1)
plotMap(map2)
u2<-scoreU(map2)
e2<-scoreE(map2)

## ----map3, echo=TRUE----------------------------------------------------------
map3<-selectLocalOptimalSNPs(bov80K,n=2000,w1=1,w2=0)
plotMap(map3)
u3<-scoreU(map3)
e3<-scoreE(map3)

