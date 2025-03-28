summary.mORT <-
function (object, cols_to_remove = "link", order_by = c("relevance", 
    "odds_ratio", "log2_odds_ratio", "p.value", "BH.value", "bonferroni"), 
    top = NULL, as.formattable = FALSE, ...) 
{
    alternative <- attr(object, "alternative")
    order_by <- match.arg(order_by)
    if (order_by == "relevance") 
        object$order <- order(object$relevance, decreasing = TRUE)
    if (order_by == "odds_ratio") 
        object$order <- order(object$NES, decreasing = TRUE)
    if (order_by == "log2_odds_ratio") 
        object$order <- order(object$NES, decreasing = TRUE)
    if (order_by == "p.value") 
        object$order <- order(object$p.value)
    if (order_by == "BH.value") 
        object$order <- order(object$BH.value)
    if (order_by == "bonferroni") 
        object$order <- order(object$B.value)
    tmp <- object[object$order, ]
    tmp$log2_odds_ratio <- round(tmp$log2_odds_ratio, 4)
    tmp$odds_ratio <- round(tmp$odds_ratio, 4)
    tmp$p.value <- format(tmp$p.value, digits = 2)
    tmp$BH.value <- format(tmp$BH.value, digits = 2)
    tmp$B.value <- format(tmp$B.value, digits = 2)
    cols_to_remove <- c(cols_to_remove, "order")
    if (order_by != "relevance") 
        cols_to_remove <- c(cols_to_remove, "relevance")
    if (!is.null(cols_to_remove)) 
        tmp <- tmp[, -which(colnames(tmp) %in% cols_to_remove)]
    if (!is.null(top) & is.numeric(top)) {
        top <- as.integer(top)
        above <- min(top, sum(tmp$log2_odds_ratio >= 0))
        below <- min(top, sum(tmp$log2_odds_ratio < 0))
        nnrow <- nrow(tmp)
        if (alternative == "two.sided") {
            top <- trunc(top/2)
            tmp <- tmp[c(1:above, (nnrow - below + 1):nnrow), 
                ]
        }
        else if (alternative == "greater") 
            tmp <- tmp[1:above, ]
        else tmp <- tmp[(nnrow - below + 1):nnrow, ]
    }
    if (as.formattable) {
        colnames(tmp)[2] <- "universe\nsize"
        colnames(tmp)[3] <- "geneList\nsize"
        colnames(tmp)[4] <- "geneSet\nsize"
        colnames(tmp)[5] <- "geneList in\ngeneSet"
        colnames(tmp)[7] <- "log2\nodds_ratio"
        formattable(tmp, list(log2_odds_ratio = formatter("span", 
            style = ~style(color = ifelse(log2_odds_ratio >= 
                0, "red", "green")))))
    }
    else {
        class(tmp) <- class(tmp)[-1]
        return(tmp)
    }
}
