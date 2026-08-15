# bootstrapPruning2.R
#
# Weakest-link (cost-complexity style) pruning for the dipolar survival trees.
#
# Selection rule: build the nested sequence of pruned trees T1 > T2 > ... > Tm
# (T1 = full input tree, Tm = single node) together with the cost g_h at which
# each successive tree is reached. The cost sequence is non-decreasing along the
# sequence (weakest link removed first). Given a complexity threshold `alpha`,
# we prune away every link whose cost is <= alpha and keep the smallest tree
# whose cheapest surviving link still exceeds alpha:
#
#   - alpha < min(cost)  -> no link is worth pruning, return the full tree
#   - alpha >= max(cost) -> prune everything, return the smallest tree
#   - otherwise          -> first tree in the sequence whose cost exceeds alpha
#
# `time` and `censor` are passed as arguments (defaults "stop"/"status") so the
# leaf KM recomputation does not depend on those names being bound in the
# calling scope.

bootstrapPruning <- function(intree, dipolarmodel, alpha,
                             time = "stop", censor = "status") {
  
  # Pruning criterion for a subtree rooted at `node`:
  # average internal-node log-rank statistic, g_h = GT_h / S_h.
  pruningstat <- function(node) {
    GTh <- sum(node$Get("lrstat", filterFun = isNotLeaf))
    Sh  <- node$totalCount - node$leafCount
    GTh / Sh
  }
  
  alldata <- dipolarmodel$traindata[strtoi(rownames(intree$data)), ]
  
  # Recompute the leaf-level Kaplan-Meier fit (stored as node$KMest).
  KMcompute <- function(node) {
    subsetX    <- rownames(node$data)
    datasubset <- alldata[subsetX, ]
    Y <- Surv(datasubset[[time]], datasubset[[censor]] == 1)
    survfit(Y ~ 1)
  }
  
  # Attach leaf KM fits to whatever tree we ultimately return.
  finalize <- function(tree) {
    tree$Do(function(node) node$KMest <- KMcompute(node), filterFun = isLeaf)
    tree
  }
  
  # Prune both children (the "l" and "r" subtrees) of a named node.
  prunenodesfun <- function(opttree, prunenode) {
    Prune(opttree, pruneFun = function(x) x$name != paste0(prunenode, "l"))
    Prune(opttree, pruneFun = function(x) x$name != paste0(prunenode, "r"))
  }
  
  # Build the nested sequence of pruned trees and the cost at which each is reached.
  lrpruningfun <- function(treetoprune) {
    pstats  <- c()
    opttree <- list()
    opttree[[1]] <- Clone(treetoprune)
    opttree[[1]]$Do(function(node) node$prnstat <- pruningstat(node),
                    filterFun = isNotLeaf)
    pstats[1] <- min(opttree[[1]]$Get("prnstat", filterFun = isNotLeaf))
    i <- 0
    while (opttree[[i + 1]]$totalCount != 1) {
      i <- i + 1
      prunenode <- names(which.min(opttree[[i]]$Get("prnstat", filterFun = isNotLeaf)))
      prunenodesfun(opttree[[i]], prunenode)
      opttree[[i + 1]] <- Clone(opttree[[i]])
      opttree[[i + 1]]$Do(function(node) node$prnstat <- pruningstat(node),
                          filterFun = isNotLeaf)
      if (opttree[[i]]$totalCount != 1) {
        pstats[i + 1] <- min(opttree[[i + 1]]$Get("prnstat", filterFun = isNotLeaf))
      }
    }
    opttree <- opttree[-(i + 1)]   # drop the trailing single-node stump
    list(pstats, opttree)
  }
  
  # --- single-node input: nothing to prune ---
  if (intree$totalCount == 1) {
    return(list("No pruning, input tree is single node",
                "No pruning, input tree is single node",
                1L, finalize(intree)))
  }
  
  listprune     <- lrpruningfun(intree)
  listprunestat <- listprune[[1]]
  listprunetree <- c(intree, listprune[[2]])   # index 1 = full tree; larger index = smaller tree
  
  # --- alpha below the cheapest link: keep the full tree ---
  if (alpha < min(listprunestat)) {
    return(list(listprunestat, listprunetree, 1L, finalize(listprunetree[[1]])))
  }
  
  # --- selection: smallest tree whose cheapest surviving link exceeds alpha ---
  prune_away <- listprunestat <= alpha
  if (all(prune_away)) {
    index <- length(listprunetree)        # alpha large: smallest tree
  } else {
    index <- max(which(prune_away)) + 1L  # first tree whose cost exceeds alpha
  }
  index <- max(1L, min(index, length(listprunetree)))
  
  finaltree <- finalize(listprunetree[[index]])
  list(listprunestat, listprunetree, index, finaltree)
}