% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectSNPs_class.R
\docType{class}
\name{Locus}
\alias{Locus}
\alias{Locus-class}
\title{The Locus class}
\description{
This is a S4 \code{Locus} classs for a marker (e.g., SNP) locus.
}
\details{
This \code{Locus} class abstracts the information for a marker or SNP locus. It consists of six slots,
as follows.
}
\section{Slots}{

 \describe{
   \item{\code{Name}:}{A character string as the name of this marker locus.}
   \item{\code{Chromosome}:}{A character string as the  name of the chromosome for this marker.}
   \item{\code{Position}:}{A numeric number  as the chromosomal position of this marker.}
   \item{\code{Maf}:}{A numeric number for the minor allele frquency at this locus.}
   \item{\code{Type}:}{A character specifying the type of this locus, which can be "A", "B", or "C".
   The default value is "C" when the input is missing.}
   \item{\code{Status}:}{A integer specifying the status of tihs locus, wihch can be 1 for obligatory
   markers or 0 otherwise. The default value is 0 when missing the input.}
 }
}

\author{
Nick X-L Wu
}
