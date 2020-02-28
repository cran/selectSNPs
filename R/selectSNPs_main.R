
#' @title Select locally-optimized SNPs on one or more chromosomes
#'
#' @description This function is used to select locally-optimal SNPs using an unified, local function 
#' uder varied scenarios. The idea is to divide each chromosome into equal segments, and within each 
#' segment, an local optimal SNP is selected as the one having the largest weighted score between the 
#' uniformmness (U score) and a function (average Shannon entropy, or E score) of minor allele frequency. 
#' If the weight for the E score is zero, it selects (approximately) uniformly-distributed SNPs. If the
#' weight for the U score is zero, it selects SNPs with the largest minor allele frequencies within 
#' pre-defined local chromosome regions or bins.  
#'
#' @param object An input object, which can be a vector of map positions, a Chrom object, a Map object,
#' or a data frame containing map information for all the candidate SNPs.
#' @param n NUmber of SNP positions to be filled.
#' @param f A vector of minor allele frequencies for all the candidate SNPs.
#' @param t1 A penalty parameter (0~100) for the E score, which by default is 1.
#' @param t2 A shrinkage parameter (0-100) for the U score, which by default is 1.
#' @param w1 A weight for the impact of MAF on selecting markers, which by default is 0.5.
#' @param w2 A weight for the impact of map position on selecting markers, which by default is 0.5.
#' @param ... Extra input parameters, as needed.
#' @return A \code{Map} object with locally-optimal SNPs when the input is a \code{Map} object, or a
#' \code{Chrom} object with locally-optimal SNPs when the input is a \code{Chrom} object, or a vector
#' of map positions of locally-optimal SNPs when the input are a vector of map positions of SNPs.
#' 
#' @examples
#' \donttest{
#' data(bov80K)
#' map1<-selectLocalOptimalSNPs(bov80K,bin=1,n=1000,w1=0,w2=1)
#' plot(map1)
#' scoreU(map1)
#' scoreE(map1)
#' 
#' map2<-selectLocalOptimalSNPs(bov80K,bin=2,n=1000,w1=0,w2=1)
#' plot(map2)
#' scoreU(map2)
#' scoreE(map2)
#' 
#' map3<-selectLocalOptimalSNPs(bov80K,bin=4,n=1000,w1=0.5,w2=0.5)
#' plot(map3)
#' scoreU(map3)
#' scoreE(map3)
#' } 
#'
#' @export
#' @docType methods
#' @rdname selectLocalOptimalSNPs-methods
#' @author Nick X-L Wu
#' @references {
#'     1. Wu X-L, Li H, Ferretti R, Simpson B, Walker J, Parham J, Mastro L, Qiu J, Schultz T, Tait
#' RG. Jr., and Bauck S. (2020) A unified local objective function for optimally selecting SNPs on arrays
#'  for agricultural genomics applications. Anim. Genet. (in press)
#'
#'     2. Wu XL, Xu J, Feng G, Wiggans GR, Taylor JF, He J, Qian C, Qiu J, Simpson B, Walker J,
#'  Bauck S. Optimal Design of Low-Density SNP Arrays for Genomic Prediction: Algorithm and Applications.
#'  PLoS One. 2016, 11(9):0161719. doi: 10.1371/journal.pone.0161719. eCollection 2016.
#'  }
#'
setGeneric("selectLocalOptimalSNPs",def=function(object,...){
  standardGeneric("selectLocalOptimalSNPs")
})

#' @rdname selectLocalOptimalSNPs-methods
#' @aliases selectLocalOptimalSNPs,numeric-method
setMethod("selectLocalOptimalSNPs",signature(object="numeric"),
          function(object,f,n,t1=1,t2=1,w1=0.5,w2=0.5){
            if(missing(f)){
              f<-rep(0,length(object))
            }
            wkdat<-data.frame(x=object,f=f)
            wkdat<-wkdat[order(wkdat$x),]

            nx<-length(wkdat$x)
            if(nx>3){
              z<-combn(wkdat$x[2:(nx-1)],n)

              nz<-ncol(z)
              u<-numeric()
              wt<-numeric()
              for (j in 1:nz){
                zj<-c(wkdat$x[1],z[,j],wkdat$x[nx])
                fj=wkdat$f[(wkdat$x%in%c(z[,j]))]

                wtj<-localScore(f=fj,x=zj,t1=t1,t2=t2,w1=w1,w2=w2)
                wt<-c(wt,wtj)
              }
              out<-z[,which.max(wt)]
            }else{
              out<-wkdat$x[2]
            }
            return(out)
          })

#' @rdname selectLocalOptimalSNPs-methods
#' @aliases selectLocalOptimalSNPs,Chrom-method
setMethod("selectLocalOptimalSNPs",signature(object="Chrom"),
          function(object,n,t1=1,t2=1,w1=0.5,w2=0.5){
            x<-object@Position
            f<-object@Maf
            wkdat<-data.frame(x=object@Position,
                              f=object@Maf)
            wkdat<-wkdat[(order(wkdat$x)),]

            low<-min(wkdat$x)
            upp<-max(wkdat$x)
            selLst<-c(low,upp)

            y<-seq(from=low,to=upp,length.out=(n-1))
            nj<-rep(1,times=(length(y)-1))

            go<-TRUE
            while(go){
              # number of SNPs in each bin
              nBin<-length(nj)
              z<-numeric()
              nsize<-integer()
              for (j in 1:nBin){
                datj<-subset(wkdat,
                             subset=((wkdat$x>=y[j]) & (wkdat$x<=y[(j+1)])))
                nsize<-c(nsize,nrow(datj))
                if((j==1) || (j==nBin)){
                  nsize[j]<-nsize[j]-1
                }
              }
            
              if(any(nsize<nj)){
                idx<-(1:nBin)[(nsize<nj)]
                idx[idx==1]<-idx[idx==1]+1
              
                if(length(idx)>0){
                  y<-y[-c(idx)]
                  dy<-diff(y)
                  nj<-floor((n-2)*dy/sum(dy)+0.5)
                }
              }
              
              nsize<-nsize[(nsize!=0)]
              go=any(nsize<nj)
            }

            # selecct locally optimal SNPs
            xsel<-numeric()
            nBin=length(nj)
            for (j in 1:nBin){
              datj<-subset(wkdat,
                           subset=((wkdat$x>=y[j]) & (wkdat$x<=y[(j+1)])))
              if((nrow(datj)>0) && (nj[j]>0)){
                tks<-data.frame(x=c(y[j],y[(j+1)]),
                                f=0)
                datj<-rbind(datj,tks)
                datj<-datj[order(datj$x),]

                selj<-selectLocalOptimalSNPs(object=datj$x,
                                             f=datj$f,
                                             n=nj[j],
                                             t1=t1,
                                             t2=t2,
                                             w1=w1,
                                             w2=w2)
                #cat(j,nj[j],length(selj),"......\n")
                xsel<-c(xsel,selj)
              }
            }

            Position<-sort(unique(c(min(x),xsel,max(x))))
            Name<-object@Name[(object@Position%in%Position)]
            Chromosome<-object@Chromosome
            Maf<-object@Maf[(object@Position%in%Position)]
            Type<-object@Type[(object@Position%in%Position)]
            Status<-object@Status[(object@Position%in%Position)]

            out<-new("Chrom",
                     Name = as.character(Name),
                     Chromosome = as.character(Chromosome),
                     Position = as.numeric(Position),
                     Maf = as.numeric(Maf),
                     Type = as.character(Type),
                     Status = as.integer(as.character(Status)))

            return(out)
          })

#' @rdname selectLocalOptimalSNPs-methods
#' @aliases selectLocalOptimalSNPs,Map-method
setMethod("selectLocalOptimalSNPs",signature(object="Map"),
          function(object,n,t1=1,t2=1,w1=0.5,w2=0.5){
            map<-cleanMap(object)
            nj<-nSNPByChrom(map,n=n)
            nChr<-length(map)
            mapLst<-vector("list",nChr)

            cat("Selecting SNPs by chromosomes ......\n")
            for (j in 1:nChr){
              chrj<-(map[[j]])
              outj<-selectLocalOptimalSNPs(chrj,
                                           n=nj[j],
                                           t1=t1,
                                           t2=t2,
                                           w1=w1,
                                           w2=w2)
              mapLst[[j]]<-outj
              cat(j,chrj@Chromosome,length(chrj@Position),nj[j],length(outj@Position),"......\n")
            }
            out<-new("Map",.Data=mapLst)
            return(out)
          })


