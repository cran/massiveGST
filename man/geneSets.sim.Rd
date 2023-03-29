\name{geneSets.sim}
\alias{geneSets.sim}
\title{Compute the similarities between a collection of gene sets.}
\description{Compute the similarities between a collection of gene sets using a convex function of the Jaccard and overlap indeces.}
\usage{geneSets.sim(gs, eps = 0.25)}
\arguments{
  \item{gs}{a character vector of gene-sets.}
  \item{eps}{a real value between 0.0 and 1.0 controlling the contribution of the Jaccard and overlap similaties to their convex combination; eps = 0.25 (default), see details.}
}
\details{The similarity between the gene-set is computed a convex combination of the Jaccard and overlap similarities. See the reference for further details.}
\value{returns an object of class "dist", where the values are the similaties between gene sets.
}
\references{Cerulo, Pagnotta (2022) \doi{10.3390/e24050739}}

\author{Stefano M. Pagnotta}

\seealso{\link{dist}}
\examples{
library(massiveGST)

# get the gene-sets
geneSets <- get_geneSets_from_msigdbr(category = "H", what = "gene_symbol")[1:5]

# compute the similarities
geneSets.sim(geneSets)
ssim <- geneSets.sim(geneSets)
ssim <- as.matrix(ssim)
diag(ssim) <- 1
ssim
}
