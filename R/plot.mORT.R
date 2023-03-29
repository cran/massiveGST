plot.mORT <-
function(x, gene_sets = NULL, order_by = "log2_odds_ratio", top = 30, 
                      eps = 0.25, 
                      as.network = FALSE, 
                      similarity_threshold = 1/3, 
                      manipulation = FALSE, 
                      autoResize = TRUE, ...) {
  
  x$NES <- pnorm(x$log2_odds_ratio)
  x$logit2NES <- x$log2_odds_ratio
  x$actualSize <- x$geneSet_size
  alternative <- attr(x, "alternative")
  if(alternative == "two.sided") x$abs_logit2NES <- abs(x$logit2NES)
  if(order_by == "log2_odds_ratio") order_by = "logit2NES"
  if(order_by == "odds_ratio") order_by = "NES"
  class(x)[1] <- "mGST"
  plot(x, gene_sets = gene_sets, order_by = order_by, top = top, 
                      eps = eps, 
                      as.network = as.network, 
                      similarity_threshold = similarity_threshold, 
                      manipulation = manipulation, 
                      autoResize = autoResize, ...)
}
