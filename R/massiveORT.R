massiveORT <-
function(gene_list, gene_sets, universe = NULL,
                       alternative = c("greater", "less", "two.sided")) {
  
  if(is.null(universe)) {
    universe <- unlist(gene_sets)
    names(universe) <- NULL
    universe <- unique(universe)
  }
  
  message(alternative <- match.arg(alternative))
  
  collection <- unlist(sapply(gene_sets, function(x) attr(x, "collection")))
  if(is.null(collection)) collection <- rep(NA, length(gene_sets))
  link <- unlist(sapply(gene_sets, function(x) attr(x, "link")))
  size <- unlist(sapply(gene_sets, length))
 
  gene_sets <- lapply(gene_sets, function(x) intersect(x, universe))
  gene_list <- intersect(gene_list, universe)
  actualSize <- unlist(sapply(gene_sets, length))
  
  
  results <- data.frame(
    row.names = names(gene_sets),
    collection, 
    universe_size = length(universe),
    geneList_size = length(gene_list),
    geneSet_size = actualSize,
    geneList_in_GeneSet = NA,
    odds_ratio = NA, 
    log2_odds_ratio = NA, 
    p.value = NA
    )
  
  theGeneSet_template <- rep(FALSE, length(universe))
  names(theGeneSet_template) <- universe
  
  drawn <- theGeneSet_template
  drawn[intersect(gene_list, universe)] <- TRUE
  
  k <- 1
  for(k in 1:length(gene_sets)) {
    theGeneSet <- theGeneSet_template
    theGeneSet[gene_sets[[k]]] <- TRUE
    ans <- fisher.test(tt <- table(drawn, theGeneSet), alternative =  alternative)
    results$geneList_in_GeneSet[k] <- tt[2, 2]
    results$odds_ratio[k] <- ans$estimate
    results$p.value[k] <-  ans$p.value
    
  }
  results$log2_odds_ratio <- log2(results$odds_ratio)
  results$BH.value <- p.adjust(results$p.value, method = "BH")
  results$B.value <- p.adjust(results$p.value, method = "bonferroni")
  
  
  positive <- which(results$log2_odds_ratio >= 0)
  order_p <- rank(results[positive, "log2_odds_ratio"]) + rank(results[positive, "geneSet_size"]) + rank(1-results[positive, "p.value"])
  names(order_p) <- rownames(results)[positive]
  
  negative <- which(results$log2_odds_ratio < 0)
  order_n <- rank(-results[negative, "log2_odds_ratio"]) + rank(results[negative, "geneSet_size"]) + rank(1 - results[negative, "p.value"])
  names(order_n) <- rownames(results)[negative]
  
  results$relevance <- relevance <- c(order_p, -order_n)[rownames(results)]
  
  if(alternative == "less") result$relevance <- -result$relevance
  
  results$link <- as.character(link[rownames(results)])
  class(results$link) <- "hyperlink"
  
  results <- results[order(results$relevance, decreasing = TRUE),]
  
  class(results) <- c("mORT", "data.frame")
  attr(results, "alternative") <- alternative
  
  invisible(results)
}
