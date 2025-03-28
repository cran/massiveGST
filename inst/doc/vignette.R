## ----eval =FALSE--------------------------------------------------------------
# install.packages("massiveGST")

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("stefanoMP/massiveGST")

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(massiveGST, quietly = TRUE))

## -----------------------------------------------------------------------------
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "pre_ranked_list.txt")
geneProfile <- get_geneProfile(fname)
class(geneProfile)
head(geneProfile)
tail(geneProfile)

## -----------------------------------------------------------------------------
fname <- system.file("extdata", package="massiveGST")
fname <- file.path(fname, "h.all.v2024.1.Hs.symbols.gmt")
geneSets <- get_geneSets_from_local_files(fname)
class(geneSets)
head(names(geneSets))

## -----------------------------------------------------------------------------
fname <- file.path(tempdir(), "hallmarks.gmt")
write_geneSets_to_gmt(geneSets, fileName = fname)

## -----------------------------------------------------------------------------
fname <- file.path(tempdir(), "hallmarks.gmt")
tmp <- get_geneSets_from_local_files(fname)

class(geneSets)
head(names(geneSets))

## ----eval = TRUE--------------------------------------------------------------
system.time({ans <- massiveGST(geneProfile, geneSets, alternative = "two.sided")})
class(ans)
ans[1:6,]

## ----eval=FALSE---------------------------------------------------------------
# fname <- file.path(tempdir(), "massiveGST_results.tsv")
# save_as_tsv(ans, file_name = fname)
# 
# fname <- file.path(tempdir(), "massiveGST_results.xls")
# save_as_xls(ans, file_name = fname)

## -----------------------------------------------------------------------------
summary(ans)[1:10,]

## -----------------------------------------------------------------------------
summary(ans, order_by = "NES")[1:10,]
summary(ans, order_by = "p.value")[1:10,]
summary(ans, order_by = "bonferroni")[1:10,]

## -----------------------------------------------------------------------------
(tmp <- summary(ans, order_by = "p.value", cols_to_remove = c("BH.value", "B.value"))[1:10,])

## -----------------------------------------------------------------------------
summary(ans, as.formattable = TRUE)
summary(ans, order_by = "p.value", cols_to_remove = c("BH.value", "B.value"), as.formattable = TRUE)

## -----------------------------------------------------------------------------
tmp

## ----eval = TRUE--------------------------------------------------------------
summary(cut_by_significance(ans), as.formattable = TRUE)

## -----------------------------------------------------------------------------
summary(cut_by_significance(ans, level_of_significance = 0.01, where = "bonferroni"),
        cols_to_remove = c("BH.value", "NES", "size"), 
        order_by = "logit2NES",
        as.formattable = TRUE)

## -----------------------------------------------------------------------------
tmp <- cut_by_significance(ans)
summary(cut_by_logit2NES(tmp), as.formattable = TRUE, order_by = "NES")
summary(cut_by_NES(tmp), as.formattable = TRUE, order_by = "NES")

## ----bar-plot, eval = TRUE----------------------------------------------------
plot(ans)

## -----------------------------------------------------------------------------
plot(cut_by_significance(ans), top = 30)

## ----network------------------------------------------------------------------
plot(cut_by_significance(ans), gene_sets = geneSets, as.network = TRUE)

## ----massive OR-test, echo=TRUE-----------------------------------------------
#geneSets <- get_geneSets_from_msigdbr(category = "C5", subcategory = "CC", what = "gene_symbol")

#fname <- system.file("extdata", package="massiveGST")
#fname <- file.path(fname, "pre_ranked_list.txt")
#geneProfile <- get_geneProfile(fname)

## ----echo=TRUE----------------------------------------------------------------
geneList <- names(head(geneProfile, 500))

## ----echo=TRUE----------------------------------------------------------------
ans <- massiveORT(geneList, geneSets)

## ----echo=TRUE----------------------------------------------------------------
summary(cut_by_significance(ans), as.formattable = TRUE)

## ----echo=TRUE----------------------------------------------------------------
plot(cut_by_significance(ans))

## -----------------------------------------------------------------------------
plot(ans, gene_sets = geneSets, as.network = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

