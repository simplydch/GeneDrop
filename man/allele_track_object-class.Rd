% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allele_tracking_functions.R
\docType{class}
\name{allele_track_object-class}
\alias{allele_track_object-class}
\alias{allele_track_object}
\alias{get_allele_dam}
\alias{get_allele_sire}
\title{An S4 class to store to store allele tracking output}
\usage{
get_allele_dam(allele_track_object)

get_allele_sire(allele_track_object)
}
\description{
An S4 class to store to store allele tracking output
}
\section{Functions}{
\itemize{
\item \code{get_allele_dam}: A method to access the dam allele section of the allele_track_object

\item \code{get_allele_sire}: A method to access the sire allele section of the allele_track_object
}}

\section{Slots}{

\describe{
\item{\code{id}}{A numeric or character value which corresponds to ID in pedigree of gene_drop_object}

\item{\code{locus}}{A numeric  value which corresponds to locus position in of gene_drop_object}

\item{\code{allele_dam}}{A matrix to contain allele tracking information for the dam}

\item{\code{allele_sire}}{A matrix to contain allele tracking information for the sire}
}}

\examples{

}
