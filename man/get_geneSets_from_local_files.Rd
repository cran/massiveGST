\name{get_geneSets_from_local_files}

\alias{get_geneSets_from_local_files}

\title{Load the gene-sets collection from local gmt files}

\description{Load the gene-sets collection from local gmt files}

\usage{get_geneSets_from_local_files(ffiles)}

\arguments{
  \item{ffiles}{a character string or a list of a character pointing to local files}
}

\value{A vector list of gene-sets}

\author{Stefano M. Pagnotta}

\seealso{\link{write_geneSets_to_gmt}
}

\examples{library(massiveGST)

# getting one collection
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "h.all.v2024.1.Hs.symbols.gmt")
length(geneSets <- get_geneSets_from_local_files(fname))
head(geneSets)

# getting two or more collections
geneSets <- get_geneSets_from_local_files(c(fname, fname))
length(geneSets)

}
