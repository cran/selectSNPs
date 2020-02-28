setGeneric("cleanMap",def=function(object) {
  standardGeneric("cleanMap")})

setMethod("cleanMap",signature(object="data.frame"),
          function(object){
            map<-object
            # remove NA data
            map<-map[(!is.na(map$Name)),]
            map<-map[(!is.na(map$Chromosome)),]
            map<-map[(!is.na(map$Position)),]
            map<-map[(!map$Chromosome%in%c("0","Y","MT")),]

            # remove SNPs without position
            map$Position<-as.numeric(as.character(map$Position))
            map<-map[(map$Position>0),]

            # remove duplicated SNPs
            map<-map[(!duplicated(map$Name)),]

            # re-number chromosomes
            chrs<-unique(map$Chromosome)
            chrs.int<-as.integer(as.character(chrs[(chrs%in%c(0:1000))]))
            chrs.str<-chrs[(!chrs%in%chrs.int)]

            nchr.max<-max(chrs.int)
            nToNum<-length(chrs.str)

            for (j in 1:nToNum){
              nj<-nchr.max+j
              this.chr<-chrs.str[j]
              map$Chromosome[(map$Chromosome==this.chr)]<-nj
            }

            map$Chromosome<-as.integer(as.character(map$Chromosome))
            map<-map[order(map$Chromosome,map$Position),]

            return(map)
          })

setMethod("cleanMap",signature(object="Map"),
          function(object){
            daf<-as.data.frame(object)
            daf<-cleanMap(daf)
            daf<-as.Map(daf)
            return(daf)
          })

setGeneric("nSNPByChrom",def=function(object,...){
  standardGeneric("nSNPByChrom")
})

setMethod("nSNPByChrom",signature(object="numeric"),
          function(object,n){
            z<-object
            p<-z/sum(z)

            bin<-length(z)

            cnt<-n
            t<-0
            while(t<n){
              nj<-floor(cnt*p)
              t<-sum(nj)
              #cat(n,t,"......\n")
              cnt<-cnt+0.1
            }

            idx2<-integer()
            if(t>n){
              df<-t-n
              idx2<-(1:bin)[nj>1]
              idxSel<-sample(idx2,size=df,replace=F)
              nj[idxSel]<-nj[idxSel]-1
            }

            idx0<-integer()
            idx1<-integer()
            if(t<n){
              df<-n-t
              idx0<-(1:bin)[nj==0]
              l0<-length(idx0)

              idx1-(1:bin)[nj==1]
              l1<-length(idx1)

              if(l0<df){
                d0<-l1-l0
                idx0<-c(idx0+idx1[1:d0])
              }
              idxSel<-sample(idx0,size=df,replace=F)
              nj[idxSel]<-nj[idxSel]+1
            }
            return(nj)
          })


setMethod("nSNPByChrom",signature(object="Chrom"),
          function(object,n){
            x<-object@Position

            nj<-nSNPByChrom(object=x,n=n)
            return(nj)
          })

setMethod("nSNPByChrom",signature(object="Map"),
          function(object,n){
            map<-cleanMap(object)
            summ<-summary(map)

            nj<-nSNPByChrom(object=summ$Length,n=n)
            return(nj)
          })

setMethod("nSNPByChrom",signature(object="data.frame"),
          function(object,n){
            map<-cleanMap(object)
            summ<-summary(map)

            nj<-nSNPByChrom(object=summ$Length,n=n)
            return(nj)
          })

