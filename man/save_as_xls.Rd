\name{save_as_xls}
\alias{save_as_xls}

\title{Save the results in xls file format}
\description{Save the data frame coming from the massiveGST function as Excel 2003 (XLS) or Excel 2007 (XLSX) files}
\usage{
save_as_xls(x, file_name = "massiveGST.xls", ...)
}

\arguments{
  \item{x}{a data frame of "mGST" class coming from massiveGST function.}
  \item{file_name}{a character value ("massiveGST.xls" as default)}
  \item{...}{Arguments to be passed to methods}
}

\value{No return value.}

\author{Stefano M. Pagnotta}

\seealso{\link[WriteXLS:WriteXLS]{WriteXLS::WriteXLS()}, \link{massiveGST}
}

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
fname <- file.path(tempdir(), "massiveGST_results.xls")
save_as_xls(ans, file_name = fname)

}