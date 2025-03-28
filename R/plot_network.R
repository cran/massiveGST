plot_network <-
function (ttable, gs, eps = 0.25, similarity_threshold = 1/3, 
    manipulation = FALSE, autoResize = TRUE, use.clustering = FALSE, 
    edges.color = "#9E9E9E", ...) 
{
    if (use.clustering & !("clustering" %in% colnames(ttable))) {
        message("Error: no clustering found.")
        return(NULL)
    }
    similarity_matrix <- as.matrix(geneSets.sim(gs[rownames(ttable)], 
        eps = eps))
    colnames(similarity_matrix) <- rownames(similarity_matrix) <- rownames(ttable)
    diag(similarity_matrix) <- 1
    similarity_matrix_tmp <- as.dist(similarity_matrix)
    similarity_matrix <- similarity_matrix > similarity_threshold
    min_node_size <- 10
    ttable$size.tmp <- sapply(ttable$actualSize, function(x) max(min_node_size, 
        x))
    ttable$size.tmp <- floor(sapply(ttable$size.tmp, function(x) min_node_size + 
        50 * (2 * (pnorm(x, sd = 250) - 0.5))))
    ttable$color <- 1 + round(99 * pnorm(ttable$logit2NES, sd = 0.29999999999999999))
    ttable$color <- colorRampPalette(c("green2", "red"))(100)[ttable$color]
    if (use.clustering) 
        ttable$color <- ttable$clustering
    nodes <- data.frame(id = 1:nrow(ttable), label = rownames(ttable), 
        size = ttable$size.tmp, color = ttable$color)
    K <- nrow(ttable)
    edges <- data.frame(from = rep(2:K, (K - 1):1), to = rep(1:(K - 
        1), (K - 1):1))
    for (k in 2:K) edges[which(edges[, "to"] == k - 1), "from"] <- k:K
    wwhich <- which(unclass(as.dist(similarity_matrix)) > 0)
    edges <- edges[wwhich, ]
    edges$color <- edges.color
    similarity_matrix_tmp <- similarity_matrix_tmp[wwhich]
    similarity_matrix_tmp <- (similarity_matrix_tmp - min(similarity_matrix_tmp))/(max(similarity_matrix_tmp) - 
        min(similarity_matrix_tmp))
    similarity_matrix_tmp <- cut(similarity_matrix_tmp, include.lowest = TRUE, 
        breaks = 0:10/10)
    edges$width <- as.numeric(similarity_matrix_tmp)
    ans <- visNetwork(nodes, edges, ...) %>% visOptions(manipulation = manipulation, 
        autoResize = autoResize) %>% visIgraphLayout()
    return(ans)
}
