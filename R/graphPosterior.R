graphPosterior <- function(tpost) {
  n.post <- dim(tableObjs(tpost))[3]
  n.gene <- dim(perturbationObj(tpost))[1]
  n.sample <- dim(perturbationObj(tpost))[2]
  parent.post.table <- matrix(0, n.gene, n.gene + 1)
  for (i in seq_len(n.gene)) {
    for (j in seq_len(n.post)) {
      gr <- graphObjs(tpost)[, , j]
      pp <- gr[1:maxDegree(inputParams(tpost)), i]
      if (all(pp == 0)) {
        parent.post.table[i, 1] <- parent.post.table[i, 1] + 1
      } else {
        for (k in seq_len(length(pp))) {
          if (pp[k] != 0) parent.post.table[i, pp[k] + 1] <- parent.post.table[i, pp[k] + 1] + 1
        }
      }
    }
  }
  parent.post.table <- parent.post.table / n.post
  if (is.null(rownames(perturbationObj(tpost)))) {
    rownames(parent.post.table) <- seq_len(nrow(parent.post.table))
    colnames(parent.post.table) <- 0:(ncol(parent.post.table) -
      1)
  }
  else {
    rownames(parent.post.table) <- rownames(perturbationObj(tpost))
    colnames(parent.post.table) <- c("-", rownames(perturbationObj(tpost)))
  }
  return(parent.post.table)
}
