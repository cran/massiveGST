## ---- eval =FALSE-------------------------------------------------------------
#  install.packages("massiveGST")

## ---- eval=FALSE--------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#  BiocManager::install("stefanoMP/massiveGST")

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
geneSets <- get_geneSets_from_msigdbr(category = "H", what = "gene_symbol")
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

## ---- eval = TRUE-------------------------------------------------------------
system.time({ans <- massiveGST(geneProfile, geneSets, alternative = "two.sided")})
class(ans)
ans[1:6,]

## ---- eval=FALSE--------------------------------------------------------------
#  fname <- file.path(tempdir(), "massiveGST_results.tsv")
#  save_as_tsv(ans, file_name = fname)
#  
#  fname <- file.path(tempdir(), "massiveGST_results.xls")
#  save_as_xls(ans, file_name = fname)

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

## ---- eval = TRUE-------------------------------------------------------------
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

## ---- eval = TRUE-------------------------------------------------------------
system.time({C5BP_gs <- get_geneSets_from_msigdbr(category = "C5", 
                                                  subcategory = "BP", 
                                                  what = "gene_symbol")})

system.time({C5MF_gs <- get_geneSets_from_msigdbr(category = "C5", 
                                                  subcategory = "MF", 
                                                  what = "gene_symbol")})


system.time({C5CC_gs <- get_geneSets_from_msigdbr(category = "C5", 
                                                  subcategory = "CC", 
                                                  what = "gene_symbol")})

system.time({H_gs <- get_geneSets_from_msigdbr(category = "H", 
                                               what = "gene_symbol")})

# merging gene-sets collections
geneSets <- c(C5MF_gs, C5CC_gs, C5BP_gs, H_gs)
length(geneSets)

# running the analysis
system.time({ans <- massiveGST(geneProfile, geneSets, 
                               alternative = "greater")})

# removing non significant results
ans <- cut_by_significance(ans, 
                           level_of_significance = 0.05, 
                           where = "bonferroni")

# Tabular results
summary(ans, as.formattable = TRUE, cols_to_remove = "BH.value")

# Saving the table
fname <- file.path(tempdir(), "massiveGST_results.tsv")
save_as_xls(ans, file_name = fname)

# Inspecting the network of gene-sets
plot(ans, gene_sets = geneSets, as.network = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

