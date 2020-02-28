#' @title The Locus class
#'
#' @description This is a S4 \code{Locus} classs for a marker (e.g., SNP) locus.
#' @details This \code{Locus} class abstracts the information for a marker or SNP locus. It consists of six slots,
#' as follows.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{Name}:}{A character string as the name of this marker locus.}
#'    \item{\code{Chromosome}:}{A character string as the  name of the chromosome for this marker.}
#'    \item{\code{Position}:}{A numeric number  as the chromosomal position of this marker.}
#'    \item{\code{Maf}:}{A numeric number for the minor allele frquency at this locus.}
#'    \item{\code{Type}:}{A character specifying the type of this locus, which can be "A", "B", or "C".
#'    The default value is "C" when the input is missing.}
#'    \item{\code{Status}:}{A integer specifying the status of tihs locus, wihch can be 1 for obligatory
#'    markers or 0 otherwise. The default value is 0 when missing the input.}
#'  }
#'
#' @name Locus
#' @rdname Locus
#' @aliases Locus-class
#' @exportClass Locus
#' @author Nick X-L Wu
#'
setClass("Locus",
         representation(Name = "character",
                        Chromosome = "character",
                        Position = "numeric",
                        Maf="numeric",
                        Type="character",
                        Status = "integer"))


#' @title The Chromosome class
#'
#' @description This is a S4 \code{Chrom} classs for a single chromosome
#' @details This \code{Chrom} class abstracts the information for all the markers on a single chromosome. It has
#' the following slots.
#'
#'@section Slots:
#'  \describe{
#'    \item{\code{Chromosome}:}{A character string as the chromosome name or index.}
#'    \item{\code{Name}:}{A character vector containing the names for all the markers on this chromosome.}
#'    \item{\code{Position}:}{A numeric vector containing the positions of all the markers on this chromosome.}
#'    \item{\code{Maf}:}{A numeric vector containing the minor allele frequencies of all the markers on
#'     this chromosome.}
#'    \item{\code{Type}:}{A character vector containing the types for all the markers on this chromosome.}
#'    \item{\code{Status}:}{An integer vector containing the statuses for all the markers on this chromosome.}
#'  }
#'
#' @name Chrom
#' @rdname Chrom
#' @aliases Chrom-class
#' @exportClass Chrom
#' @author Nick X-L Wu
#'
setClass("Chrom",
         representation(Chromosome = "character",
                        Name = "character",
                        Position = "numeric",
                        Maf="numeric",
                        Type="character",
                        Status = "integer"))


#' @title The Map class
#'
#' @description This is a S4 \code{Map} class representing a genome.
#' @details This \code{Map} class has information of all the markers on a map, which represents either entirely
#' or partially a genome. It consists of multiple chromosomes, though a \code{Map} object can have only one
#' chromosome.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{.Data}:}{A list of \code{Chrom} objects.}
#'  }
#'
#' @name Map
#' @rdname Map
#' @aliases Map-class
#' @exportClass Map
#' @author Nick X-L Wu
#'
setClass("Map",
         contains = c("list"))


#' A \code{Map} object with 80K bovine SNPs
#'
#' This example \code{Map} object contains map information for 76,694 SNPs on 30 chromosomes, with chromosome X
#' represnted by 30. Each \code{Chrom} object has six slots: Locus names (Name), Chromosome, Position, minor
#' allele frquency (Maf), Type and Status.
#'
#' @docType data
#' @keywords datasets
#' @name bov80K
#' @usage data(bov80K)
#' @format A Map object with a list of 30 Chrom objects.
NULL
