## Niche optima of a species pool, evolved by Brownian motion along a random
## phylogeny and rescaled to the gradient [0, 100]. Taken from VirtualCom
## (Münkemüller & Gallien 2015), reduced to the case used here.
draw_niche_optima <- function(n.species) {
  ## a non-zero seed stops geiger from reseeding R from the clock
  tree <- geiger::sim.bdtree(0.1, 0, stop = "taxa", n = n.species,
                             seed = sample.int(.Machine$integer.max, 1L))
  tree$edge[tree$edge <= n.species] <- seq_len(n.species)
  tree$tip.label <- seq_len(n.species)
  trait <- ape::rTraitCont(tree, model = "BM", sigma = 0.01, root.value = 0)
  optima <- (trait - min(trait)) * 100 / abs(max(trait) - min(trait))
  names(optima) <- seq_len(n.species)
  ## The phylogenetic signal is not used, but VirtualCom computes it and its
  ## randomisation draws random numbers: kept so that the worlds are those of
  ## the published grids.
  invisible(picante::phylosignal(optima, tree, reps = 1))
  optima
}
