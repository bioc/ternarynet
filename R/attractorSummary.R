attractorSummary <- function(tpost, post.prob.limit = 0.01, wildtype = TRUE) {
  perturbationObj(tpost) <- perturbationObj(tpost) + 2
  steadyStateObj(tpost) <- round(steadyStateObj(tpost)) + 2

  nOutcomes <- 3
  n.genes <- nrow(perturbationObj(tpost))
  n.sample <- ncol(perturbationObj(tpost))
  RowNames <- seq_len(n.genes)
  ColNames <- seq_len(n.sample)
  forced.gene.list <- NULL
  for (i in seq_len(n.sample)) {
    fg <- seq_len(n.genes)[perturbationObj(tpost)[, i] != 2]
    fs <- perturbationObj(tpost)[fg, i]
    forced.gene.list <- append(
      forced.gene.list,
      list(list(forced.genes = fg, forced.states = fs))
    )
  }
  initial.state.obj <- list(forced.gene.list = forced.gene.list, initial.states = steadyStateObj(tpost))
  attr <- attractor.posterior.2(tpost, initial.state.obj, wt = wildtype)
  post.attr.summary.list <- NULL
  for (i in seq_len(n.sample)) {
    post.attr.summary.list <- append(post.attr.summary.list, list(list(obs = steadyStateObj(tpost)[, i] - nOutcomes + 1, attr = condense.attractor(attr[[i]], post.prob.limit))))
  }
  post.attr.summary.table <- NULL
  for (i in seq_len(n.sample)) {
    m <- rbind(post.attr.summary.list[[i]]$obs, post.attr.summary.list[[i]]$attr$uni.attr)
    m <- cbind(rep(i, dim(m)[1]), m, c(0, post.attr.summary.list[[i]]$attr$p.attr))
    post.attr.summary.table <- rbind(post.attr.summary.table, m)
    colnames(post.attr.summary.table) <- c("index", RowNames, "PostProb")
  }
  post.attr.summary.table
}
