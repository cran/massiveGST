\name{save_as_tsv}
\alias{save_as_tsv}

\title{Save the results in tab-separeted value file}

\description{Save the data frame coming from the massiveGST function as tab-separeted value.
}

\usage{
save_as_tsv(x, file_name = "massiveGST.tsv", sep = "\t", ...)
}
\arguments{
  \item{x}{a data frame of "mGST" class coming from massiveGST function.}
  \item{file_name}{a character value ("massiveGST.tsv" as default)}
  \item{sep}{a character value}
  \item{...}{Arguments to be passed to methods}
}


\value{No return value.}

\author{Stefano M. Pagnotta}

\seealso{\link{massiveGST}}

\examples{

library(massiveGST)

# get the gene profile
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "pre_ranked_list.txt")
geneProfile <- get_geneProfile(fname)

# get the gene-sets
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "h.all.v2024.1.Hs.symbols.gmt")
geneSets <- get_geneSets_from_local_files(fname)

# run the function
ans <- massiveGST(geneProfile, geneSets, alternative = "two.sided")

# save the results
fname <- file.path(tempdir(), "massiveGST_results.tsv")
save_as_tsv(ans, file_name = fname)

}