\name{massiveORT}
\alias{massiveORT}
\title{A wrapper to fisher.test to get over representation analysis of gene sets.}
\description{The function massiveORT essentially is a wrapper to the function fisher.test in charge to 1) arrange the input to feed fisher.test in sequence for each gene set, 2) arrange the output in a data frame compatible with the other function of the package, and 3) compute the universe of genes for the analysis. }

\usage{massiveORT(gene_list, gene_sets, universe = NULL, 
          alternative = c("greater", "less", "two.sided"))}

\arguments{
  \item{gene_list}{a list of gene names, or gene ids that have to match the corresponidng in the gene-set.}
  \item{gene_sets}{a character vector of gene-sets.}
  \item{universe}{a list of gene, or gene ids, that defines the universe for the analysis (see details); NULL by default.}
  \item{alternative}{a character string specifying the alternative hypothesis of the fisher.test; must be one of "two.sided", "greater" (default) or "less".}
}
\details{This function allows to define externally or compute the universe of reference of the analysis. By default (universe = NULL), the universe is computed starting from the gene names contributing at least once in each gene set. }
\value{A data frame with columns
\item{universe_size}{size of the universe of genes.}
\item{geneList_size}{size of intersection between the gene list and the universe.}
\item{geneSet_size}{size of intersection between the gene set and the universe.}
\item{geneList_in_GenesSet_size}{size of the intersection between the geneList and the geneSet.}
\item{odds_ratio}{odd ratio coming from the fisher.test}
\item{log2_odds_ratio}{log2 transformation of odds_ratio.}
\item{p.value}{p-values associated with the gene-set coming from the fisher.test}
\item{BH.value}{Benijamini and Hockberg adjustment of the p.values}
\item{B.value}{Bonferroni adjustment of the p.values}
\item{relevance}{marginal ordering of the table.}
}
\references{Cerulo, Pagnotta (2022) \doi{10.3390/e24050739}}

\author{Stefano M. Pagnotta}

\seealso{\link{fisher.test}, \link{cut_by_significance}}

\examples{
library(massiveGST)

# get the gene profile
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "pre_ranked_list.txt")
geneProfile <- get_geneProfile(fname)
geneList <- names(head(geneProfile, 1000))

# get the gene-sets
geneSets <- get_geneSets_from_msigdbr(category = "C5", subcategory = "CC", what = "gene_symbol")
geneSets <- geneSets[1:250]

# run the function
ans <- massiveORT(geneList, geneSets)
cut_by_significance(ans)

plot(cut_by_significance(ans), geneSets,as.network = TRUE)
}