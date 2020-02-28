#' @importFrom graphics abline axis barplot image par plot segments title
#' @importFrom grDevices topo.colors
#' @importFrom methods new
#' @importFrom stats sd na.omit
#' @importFrom utils combn read.table
#' @importFrom Rcpp sourceCpp
#' @useDynLib selectSNPs
#' 
#' @title Convert a \code{Map} Object to a Data Frame
#' @description This inherited S3 function is devised to convert a Map object to a data frame. The Map
#' object has six slots: Marker name (Name), Chromosome, Positin, minor alllel frequency (Maf), Tyep,
#' and Status.
#' @param x An input \code{Map} object.
#' @return A data frame having the map informtion.
#' 
#' @rdname as.data.frame
#' @author Nick X-L Wu
#' @examples
#' data(bov80K)
#' class(bov80K)
#' tmpdat<-as.data.frame(bov80K)
#' class(tmpdat)
#'
#' @export
as.data.frame<-function(x){
  UseMethod("as.data.frame",x)
}

#' @export
as.data.frame.Map<-function(x){
  nChr<-length(x)
  for (j in 1:nChr){
    Xj<-(x[[j]])
    outj<-data.frame(Name=Xj@Name,
                     Chromosome=Xj@Chromosome,
                     Position=Xj@Position,
                     Maf=Xj@Maf,
                     Type=Xj@Type,
                     Status=Xj@Status)
    if(j==1){
      out<-outj
    }else{
      out<-rbind(out,outj)
    }
  }
  out$Name<-as.character(out$Name)
  out$Chromosome<-as.character(out$Chromosome)
  out$Position<-as.numeric(out$Position)
  out$Maf<-as.numeric(out$Maf)
  out$Type<-as.character(out$Type)
  out$Status<-as.integer(as.character(out$Status))
  return(out)
}

#' @title Convert a Data Frame to a \code{Map} Object.
#'
#' @description This S4 function convert or creates a \code{Map} object. The input can be either a data frame or
#' a list of \code{Chrom} objects. If the data for minor allele frequency (Maf) are missing in the
#' input, their values are to be NAs. The default values for "Type" and "Status" are "C" and 0,
#' repectivelyl, if they are missing in the input.
#'
#' @param object An input object, which can be either a data frame or a list of Chrom objects.
#' @return An \code{Map} object. 
#' 
#' @export
#' @docType methods
#' @rdname as.Map-methods
#' @author Nick X-L Wu
#' @examples
#' tmpdat<-as.data.frame(bov80K)
#' class(tmpdat)
#' tmpmap<-as.Map(tmpdat)
#' class(tmpmap)
#'
setGeneric("as.Map",def=function(object){
  standardGeneric("as.Map")
})

#' @rdname as.Map-methods
#' @aliases as.Map,list-method
setMethod("as.Map",signature(object="list"),
          function(object){
            if(class(object[[1]])!="Chrom"){
              stop("Only accept a list of Chrom objects as the input!\n")
            }
            out<-new("Map",.Data = object)
            return(out)
          })

#' @rdname as.Map-methods
#' @aliases as.Map,data.frame-method
setMethod("as.Map",signature(object="data.frame"),
          function(object){
            object<-object[order(object$Chromosome,object$Position),]
            chrs<-unique(object$Chromosome)
            nChr<-length(chrs)

            map<-vector('list',nChr)
            for (j in 1:nChr){
              this.chr<-chrs[j]
              Xj<-object[(object$Chromosome==this.chr),]

              colNam<-colnames(Xj)
              if(!"MAF"%in%toupper(colNam)){
                Xj$Maf<-NA
              }
              if(!"Type"%in%toupper(colNam)){
                Xj$Type<-"C"
              }
              if(!"Status"%in%toupper(colNam)){
                Xj$Status<-0
              }

              chrj<-new("Chrom",
                        Chromosome=as.character(this.chr),
                        Name=as.character(Xj$Name),
                        Position=as.numeric(Xj$Position),
                        Maf=as.numeric(Xj$Maf),
                        Type=as.character(Xj$Type),
                        Status=as.integer(as.character(Xj$Status)))
              map[[j]]<-chrj
            }
            out<-new("Map",.Data = map)
            return(out)
          })



#' @title Summary statiscs for a Map object.
#' @description This inherited S3 function is devised specifically for compute summary statistics for
#' a Map object. The output includes Chromosome names or indexes (chrom), numuber of SNPs on each
#' chromosome (nLoci), minimum map position on each chromosome (Min), maximum map position on each
#' chromosome (Max), chromosome length (Length), average spacing or bin width (mu.bw), standard
#' deviation of average spacing (sd.bw), minimum spacing (min.bw), and maximum spacing or gap
#' length (max.bw).
#' @param x An input map object.
#' @return Summary statistics of the \code{Map} object.
#' 
#' @rdname summary
#' @author Nick X-L Wu
#' @examples
#' data("bov80K")
#' summary(bov80K)
#'
#' @export
summary <- function (x) {
  UseMethod("summary", x)
}
#' @export
summary.Map<-function(x){
  map<-cleanMap(x)
  map<-as.data.frame(map)
  chrs<-sort(unique(map$Chromosome))
  out<-numeric()
  for (j in chrs){
    # summary of this map
    chrj <- map[(map$Chromosome==j),]
    jnum <- nrow(chrj)
    jmin<-min(chrj$Position,na.rm=TRUE)
    jmax<-max(chrj$Position,na.rm=TRUE)
    range<-jmax-jmin
    # summary of map bins
    ncol<-nrow(chrj)
    bin<-numeric(ncol)
    for (k in 2:ncol){
      bin[(k-1)]<-chrj$Position[k]-chrj$Position[(k-1)]
    }
    # output as one row
    out<-rbind(out,c(j,jnum,jmin,jmax,range,
                     mean(bin),sd(bin),min(bin),max(bin)))

  }

  class(out)<-"numeric"
  colnames(out)<-c("chrom","nLoci","Min","Max","Length",
                   "mu.bw","sd.bw","min.bw","max.bw")
  out<-as.data.frame(out)
  out<-out[order(out$chrom),]
  rownames(out)<-1:nrow(out)
  return(out)
}

#' @title Plotting a Map object.
#' @description This S4 generic function is used to plot a Map object. 
#' The Map object needs to have at leaset three slots: Marker names (Name), 
#' Chromosome, and Position.
#' @param object An input Map object.
#' @param plotMaf A logical value. If TRUE, MAF density will be plotted; otherwise, SNPs are plotted 
#' based on their map positions. 
#' plotted based on their chromosomal po
#' @param bin Length of the bin in which SNPs are counted. The default length is 1M bp.
#' @param maf.base  The baseline MAF, which by default is 0.3. A MAF density above maf.base
#'  is trucated to be 1.
#' @param ... Some extra input parameter values; optional.
#' @return None.
#' 
#' @export
#' @docType methods
#' @rdname plotMap-methods
#' @author Nick X-L Wu
#'
setGeneric("plotMap",def=function(object,...){
  standardGeneric("plotMap")
})

#' @rdname plotMap-methods
#' @aliases plotMap,Map-method
setMethod("plotMap",signature(object="Map"),
          function(object,plotMaf=FALSE,bin=1000000,maf.base=0.3){
            map<-as.data.frame(object)
            map<-cleanMap(map)

            chrs<-sort(unique(map$Chromosome))
            nChr<-length(chrs)  
  
            if(plotMaf){
              pos.max<-ceiling(max(map$Position)/bin)
    
              den<-matrix(0,nrow=nChr,ncol=pos.max)
              for (j in 1:nChr){
                chrj<-chrs[j]
                posj<-map$Position[(map$Chromosome==chrj)]
                mafj<-map$Maf[(map$Chromosome==chrj)]
      
                for (k in 1:pos.max){
                  upp<-k*bin
                  low<-upp-bin
        
                  mafj.k<-mafj[(posj<upp)]
                  posj.k<-posj[(posj<upp)]
        
                  mafj.k<-mafj.k[(posj.k>low)]
                  den[j,k]<-mean(mafj.k)
                }
              }
    
              den[den==0]<-NA
    
              if(missing(maf.base)){
                maf.base<-0.5
              }
              den<-den/maf.base
              den[den>1]<-1
    
              x<-1:nChr
              y <- 1:pos.max
              den2<-1-den
              image(x, y, den2, col = topo.colors(n=100), axes = T)
            }else{
              # pos to list
              pos<-vector("list")

              idx<-1
              for (j in chrs){
                tmpdat<-subset(map,map$Chromosome==j)
                posj<-tmpdat$Position
                names(posj)<-tmpdat$SNP
                pos[[idx]]<-posj
                idx<-idx+1
              }

              xlab <- "Chromosome"
              ylab <- "Position (bp)"

              n.chr <- length(pos)
              chrpos <- 1:n.chr
              thelim <- range(chrpos) + c(-0.5, 0.5)
              maxlen <- max(unlist(lapply(pos, max)))

              xlim <- thelim
              ylim <- c(0,maxlen)

              old.xpd <- par("xpd")
              old.las <- par("las")
              par(xpd = TRUE, las = 1)
              on.exit(par(xpd = old.xpd, las = old.las))

              plot(0, 0, type = "n", ylim = ylim, xlim = xlim,
              xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
              a <- par("usr")
              for (i in 1:n.chr) {
                segments(chrpos[i], min(pos[[i]]), chrpos[i], max(pos[[i]]),col="red")
                segments(chrpos[i] - 0.25, pos[[i]], chrpos[i] + 0.25, pos[[i]],col="red")
              }

              if(n.chr>15){
                even <- seq(2, length(chrpos), by = 2)
                for (i in even) {
                  axis(side = 1, at = chrpos[i], labels = "")
                  axis(side = 1, at = chrpos[i], labels = names(pos)[i], line = +0.4, tick = FALSE)
                }
              }else{
                chrs <- seq(1, length(chrpos), by = 1)
                for (i in chrs) {
                  axis(side = 1, at = chrpos[i], labels = "")
                  axis(side = 1, at = chrpos[i], labels = names(pos)[i], line = -0.4, tick = FALSE)
                }
              }
              invisible()
            }
          })



#' @title Compute the Uniformness(U) score
#' @description This function computes the (U) scores, which measures the uniformness of all the markers
#' on one or more chromsoomes.
#'
#' @param object An input object, which can be a vector of map position, a Chrom object, or a Map object.
#' @param chrom Names of the chromosomes for which the U score is to be computed; Optional.
#' @param ... Extra input parameters, as needed.
#' @return The computed U score of the input object. 
#'
#' @export scoreU
#' @docType methods
#' @rdname  scoreU-methods
#' @author Nick X-L Wu
#' @examples
#' data(bov80K)
#' scoreU(bov80K)
#' 
#' @references {
#' 1. Wu XL, Li H, Ferretti R, Simpson B, Walker J, Parham J, Mastro L, Qiu J, Schultz T, Tait RG Jr, 
#' Bauck S, (2010). A unified local objective function for optimally selecting SNPs on arrays for 
#' agricultural genomics applications. Anim Genet. 2020 Jan 31. doi: 10.1111/age.12916. 
#'
#' 2. Wu XL, Xu J, Feng G, Wiggans GR, Taylor JF, He J, Qian C, Qiu J, Simpson B, Walker J, Bauck S. 
#' (2016). Optimal Design of Low-Density SNP Arrays for Genomic Prediction: Algorithm and Applications.
#' PLoS One. 11(9):0161719. doi: 10.1371/journal.pone.0161719. 
#' }

setGeneric("scoreU",def=function(object,...) {
  standardGeneric("scoreU")})

#' @rdname scoreU-methods
#' @aliases scoreU,Chrom-method
setMethod("scoreU",signature(object="Chrom"),
          function(object){
            u<-round(UScore(object@Position),4)
            return(u)
          })

#' @rdname scoreU-methods
#' @aliases scoreU,Map-method
setMethod("scoreU",signature(object="Map"),
          function(object,chrom){
            
            map<-cleanMap(object)
            if(missing(chrom)){
              chrom<-namesOfChrom(map)
            }

            u<-numeric()
            ch<-character()
            nChr<-length(map)
            for (j in 1:nChr){
              chrj<-(map[[j]])
              if((chrj@Chromosome)%in%chrom){
                u<-c(u,UScore(chrj@Position))
                ch<-c(ch,(chrj@Chromosome))
              }
            }

            names(u)<-ch
            mu<-round(mean(u),3)
            std<-round(sd(u),3)
            barplot(u,
                    xlab="Chromosome",
                    ylim=c(0,1),
                    main=paste("U-score: Mean (SD) =",mu," (",std,")",sep=""),
                    col="purple",
                    names.arg=ch)
            abline(h=1,col="red")
            return(u)
          })


#' @title  Compute the E score
#'
#' @description This function computes the average Shannon entropy (E score) as a measure of information
#'  for all the markers on one or more chromosome.
#'
#' @param object An input object, which can be a vector of map positions, a Chrom object, or a Map object.
#' @param chrom Names of the chromosomes for which the E score is to be computed; Optional.
#' @param ... Extra input parameters, as needed.
#' @return The computed E score of the input object. 
#'
#' @export scoreE
#' @docType methods
#' @rdname scoreE-methods
#' @author Nick X-L Wu
#' @examples
#' data(bov80K)
#' scoreE(bov80K)
#' 
#' @references{
#' 1. Wu XL, Li H, Ferretti R, Simpson B, Walker J, Parham J, Mastro L, Qiu J, Schultz T, Tait RG Jr, Bauck S, 
#' (2010). A unified local objective function for optimally selecting SNPs on arrays for agricultural genomics 
#' applications. Anim Genet. 2020 Jan 31. doi: 10.1111/age.12916.
#'
#' 2. Wu XL, Xu J, Feng G, Wiggans GR, Taylor JF, He J, Qian C, Qiu J, Simpson B, Walker J, Bauck S. (2016) 
#' Optimal Design of Low-Density SNP Arrays for Genomic Prediction: Algorithm and Applications. PLoS One. 
#' 11(9):0161719. doi: 10.1371/journal.pone.0161719.
#' }
#' 
setGeneric("scoreE",def=function(object,...) {
  standardGeneric("scoreE")})

#' @rdname scoreE-methods
#' @aliases scoreE,Chrom-methods
setMethod("scoreE",signature(object="Chrom"),
          function(object){
            f<-object@Maf
            mu<-round(mean(EScore(f),na.rm=T),4)
            return(mu)
          })

#' @rdname scoreE-methods
#' @aliases scoreE,Map-method
setMethod("scoreE",signature(object="Map"),
          function(object,chrom){
            
            map<-cleanMap(object)
            if(missing(chrom)){
              chrom<-namesOfChrom(map)
            }

            out<-numeric()
            ch<-character()
            nChr<-length(map)
            for (j in 1:nChr){
              chrj<-(map[[j]])
              if((chrj@Chromosome)%in%chrom){
                ej<-round(mean(EScore(chrj@Maf),na.rm = T),4)
                out<-c(out,ej)
                ch<-c(ch,(chrj@Chromosome))
              }
            }

            names(out)<-ch
            mu<-round(mean(out),3)
            std<-round(sd(out),3)
            barplot(out,
                    xlab="Chromosome",
                    ylim=c(0,1),
                    main=paste("E-score: Mean (SD) =",mu," (",std,")",sep=""),
                    col="purple",
                    names.arg=ch)
            abline(h=1,col="red")

            return(out)
          })


#' @title Compute the weighted local score
#' @description This function computes the weighted Local score, which is the weighted average between
#' the U score and the E socre, with the weights given from the inputs.
#'
#' @param object A vector of minor allele frequency values.
#' @param chrom Names of chromosomes, for which local scores are to be computed; optioonal.
#' @param t1 A shrinkage parameter on the E score, which by default is 1. 
#' @param t2 A shrinkage parameter on the U score, which by defaukt us 1.
#' @param w1 A weight for MAF, wchich by default is 0.5.
#' @param w2 A weight for map position, which by default is 0.5.
#' @param ... Extra input parameters, as needed.
#' @return The computed local score of the input object.
#'
#' @export
#' @docType methods
#' @rdname scoreLocal-methods
#' @author Nick X-L Wu
#' @examples
#' data("bov80K")
#' scoreLocal(bov80K)
#' 
#' @references Wu XL, Li H, Ferretti R, Simpson B, Walker J, Parham J, Mastro L, Qiu J, Schultz T, 
#' Tait RG Jr, Bauck S, (2010). A unified local objective function for optimally selecting SNPs on 
#' arrays for agricultural genomics applications. Anim Genet. 2020 Jan 31. doi: 10.1111/age.12916.
#'
setGeneric("scoreLocal",def=function(object,...) {
  standardGeneric("scoreLocal")})

#' @rdname scoreLocal-methods
#' @aliases scoreLocal,Map-method
setMethod("scoreLocal",signature(object="Chrom"),
          function(object,t1=1,t2=1,w1=0.5,w2=0.5){
            f<-object@Maf
            x<-object@Position

            local<-localScore(f=f,x=x,w1=w1,w2=w2,t1=t1,t2=t2)
            return(round(local,4))
          })

#' @rdname scoreLocal-methods
#' @aliases scoreLocal,Map-method
setMethod("scoreLocal",signature(object="Map"),
          function(object,chrom,t1=1,t2=1,w1=0.5,w2=0.5){
            
            map<-cleanMap(object)
            if(missing(chrom)){
              chrom<-namesOfChrom(map)
            }

            local<-numeric()
            ch<-character()
            nChr<-length(map)
            for (j in 1:nChr){
              chrj<-(map[[j]])
              if((chrj@Chromosome)%in%chrom){
                local<-c(local,scoreLocal(chrj,t1=t1,t2=t2,w1=w1,w2=w2))
                ch<-c(ch,(chrj@Chromosome))
              }
            }

            names(local)<-ch
            mu<-round(mean(local),3)
            std<-round(sd(local),3)
            barplot(local,
                    xlab="Chromosome",
                    ylim=c(0,1),
                    main=paste("Local score: Mean (SD) =",mu," (",std,")",sep=""),
                    col="purple",
                    names.arg=ch)
            abline(h=1,col="red")
          })

#' @title Extract chromosome indexes or names
#'
#' @description This S4 function is used to extract chromosome indexes or names for \code{Chrom} object,
#' a \code{Map} object, or a data frame.
#'
#' @param object An input object, which can be a \code{Chrom} object, a \code{Map} object, or a data frame.
#' @return Chromosome names of the input object.
#'
#' @export
#' @docType methods
#' @rdname namesOfChrom-methods
#' @author Nick X-L Wu
#' @examples
#' data("bov80K")
#' namesOfChrom(bov80K)
#'
setGeneric("namesOfChrom",def=function(object){
  standardGeneric("namesOfChrom")
})

#' @rdname namesOfChrom-methods
#' @aliases namesOfChrom,Chrom
setMethod("namesOfChrom",signature(object="Chrom"),
          function(object){
            return(object@Chromosome)
          })

#' @rdname namesOfChrom-methods
#' @aliases namesOfChrom,Map
setMethod("namesOfChrom",signature(object="Map"),
          function(object){
            nChr<-length(object)
            chrs<-character()
            for (j in 1:nChr){
              chrj<-(object[[j]])
              chrs<-c(chrs,chrj@Chromosome)
            }
            return(chrs)
          })


#' @title Extract the number of loci on each chromosome
#'
#' @description This S4 function is used to extract the number of loci on each chromsome.
#'
#' @param object An input object, which can be a \code{Chrom} object, a \code{Map} object, or a data frame.
#' @return Number of SNP loci by chrommosomes.
#'
#' @export
#' @docType methods
#' @rdname numberOfLoci-methods
#' @author Nick X-L Wu
#' @examples
#' data("bov80K")
#' numberOfLoci(bov80K)
#'
setGeneric("numberOfLoci",def=function(object){
  standardGeneric("numberOfLoci")
})

#' @rdname numberOfLoci-methods
#' @aliases numberOfLoci,Chrom
setMethod("numberOfLoci",signature(object="Chrom"),
          function(object){
            numLoci<-length(object@Position)
            names(numLoci)<-object@Chromosome
            return(numLoci)
          })

#' @rdname numberOfLoci-methods
#' @aliases numberOfLoci,Map
setMethod("numberOfLoci",signature(object="Map"),
          function(object){
            nChr<-length(object)
            chrs<-character()
            numLoci<-integer()
            for (j in 1:nChr){
              chrj<-(object[[j]])
              chrs<-c(chrs,chrj@Chromosome)
              numLoci<-c(numLoci,length(chrj@Position))
            }
            names(numLoci)<-chrs
            return(numLoci)
          })

#' @title Create a new \code{Map} object.
#'
#' @description This S4 function is used to create a new \code{Map} object. The input information is taken
#' from an input DNA report in matrix format, plus a map or manifest file.
#'
#' @param map An input map or manifest file. A typical manifest file has the following columns: Index,
#' 	Name,	Chromosome,	Position,	GenTrain Score, SNP, ILMN Strand, Customer Strand, and NormID.
#' @param freq An input object containing genotypes, which can be either the name of the input DNA report
#'  file name object, or a data frame containing genotypes. A DNA report file has 9 header lines (which
#'  will be skipped when reading). The column names are on the tenth lines and the genotypes start from the
#'  eleventh lines, separated by tabs.
#' @param type The types of a SNP, which by default is "C".
#' @param status The status of a SNP, which by default is all 0 (i.e., non-obligatory SNPs).
#' @param ... Extra parameters as needed.
#' @return An \code{Map} object.
#'
#' @export
#' @docType methods
#' @rdname newMap-methods
#' @author Nick X-L Wu
#'
setGeneric("newMap",def=function(map,freq,...){
  standardGeneric("newMap")
})

#' @rdname newMap-methods
#' @aliases newMap,data.frame
setMethod("newMap",signature(map="data.frame",freq="data.frame"),
          function(map,freq){
            map<-map[,c("Name","Chromosome","Position")]
            
            if(!"Type"%in%colnames(freq)){
              freq$Type="C"
            }
            if(!"Status"%in%colnames(freq)){
              freq$Status=0
            }
            out<-merge(map,freq,by="Name")
            out<-as.Map(out)
          return(out)})


#' @rdname newMap-methods
#' @aliases newMap,character
setMethod("newMap",signature(map="character",freq="character"),
          function(map,freq,type="C",status=0){
            mapdat<-read.table(file=map,
                               header=T,
                               sep="\t",
                               stringsAsFactors = F)
            mapdat<-mapdat[,c("Name","Chromosome","Position")]
            
            # compute MAF  
            wkdat<-read.table(file=freq,
                             sep="\t",
                             header=T,
                             skip=9,
                             stringsAsFactors=F)

            snpNam <- wkdat[, 1]
            wkdat <- wkdat[, -1]
            wkdat <- as.matrix(wkdat)
            
            system.time({
              cat("compute genotype frequencies...\n\n")
              wt <- apply(X = wkdat, MARGIN = 1, FUN = function(x) {
                ts = summary(as.factor(na.omit(x)))
                return(ts/sum(ts))
              })
            })
            
            system.time({
              n <- length(wt)
              maf<-numeric(n)
              for (j in 1:n) {
                if(j%%10000==0){
                  cat(j,"......\n")
                }
                wtj<-wt[[j]]
                wtj<-wtj[!(names(wtj)%in%c("--","NA"))]
                wtj<-wtj/sum(wtj)
                
                freqAA=0
                if(any(names(wtj)%in%"AA")){
                  freqAA<-wtj[(names(wtj)%in%"AA")]
                }
                freqAB=0
                if(any(names(wtj)%in%"AB")){
                  freqAB<-wtj[(names(wtj)=="AB")]
                }
                
                maf[j]<-freqAA + 0.5 * freqAB
              }
            })
            maf[maf>0.5]<-1-maf[(maf>0.5)]
            
            mafdat <- data.frame(Name=snpNam,
                                 Maf = maf)
            mafdat$Type=type
            mafdat$Status=status

            # generate a new Map object
            out<-newMap(map=mapdat,freq=mafdat)
            return(out)})

