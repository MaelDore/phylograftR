
### Test script used to benchmark the time needed to compute combined range and state constrained grafting

## Result: we witness an inflation of combined allowed states, but computation time remains reasonable.

## Benchmark of ranges + traits (graft_missing_taxa with one range trait): synthetic Markov maps on random coalescent trees; timings in seconds (stage and full runs).
## Run after installing the package (R_LIBS pointing to it). Second part (larger configs) = edit `cfgs` at the end.

suppressMessages(library(phylograftR)); library(ape)
ns <- asNamespace("phylograftR")
set.seed(1)

# --- Markov simulator of discrete maps on a tree (root to tips) ---
sim_map <- function(phy, states, rate, neighbour = NULL, root_state = NULL) {
  ne <- nrow(phy$edge); n <- Ntip(phy); maps <- vector("list", ne); end_state <- character(n + phy$Nnode)
  end_state[n + 1] <- if (is.null(root_state)) sample(states, 1) else root_state
  ord <- seq_len(ne) # edges are in cladewise order => parents before children
  for (e in ord) {
    p <- phy$edge[e, 1]; ch <- phy$edge[e, 2]; len <- phy$edge.length[e]
    cur <- end_state[p]; t <- 0; nm <- cur; du <- numeric(0)
    repeat { w <- rexp(1, rate); if (t + w >= len) break
      du <- c(du, w); t <- t + w
      cands <- if (is.null(neighbour)) setdiff(states, cur) else neighbour(cur)
      cur <- sample(cands, 1); nm <- c(nm, cur) }
    du <- c(du, len - t); names(du) <- nm; maps[[e]] <- du; end_state[ch] <- cur }
  s <- phy; s$maps <- maps; class(s) <- c("simmap", "phylo"); s
}
# range states: all subsets of up to maxr areas among k areas
make_ranges <- function(k, maxr) { a <- LETTERS[1:k]; unlist(lapply(1:maxr, function(m) combn(a, m, paste, collapse = ""))) }
nb_range <- function(R) { ar <- strsplit(R, "")[[1]]; a <- LETTERS[1:7]
  out <- c(); for (x in setdiff(a, ar)) out <- c(out, paste(sort(c(ar, x)), collapse = "")); if (length(ar) > 1) for (x in ar) out <- c(out, paste(setdiff(ar, x), collapse = ""))
  unique(out) }

run_case <- function(n_tips, k_areas, maxr, n_missing, trait_states, n_normal, rate_r = 0.2, rate_t = 0.1, seed = 1) {
  set.seed(seed)
  phy <- rcoal(n_tips); phy$edge.length <- phy$edge.length / max(node.depth.edgelength(phy)) * 50
  rs <- make_ranges(k_areas, maxr)
  nb <- function(R) { ar <- strsplit(R, "")[[1]]; a <- LETTERS[1:k_areas]; out <- c()
    for (x in setdiff(a, ar)) if (length(ar) < maxr) out <- c(out, paste(sort(c(ar, x)), collapse = ""))
    if (length(ar) > 1) for (x in ar) out <- c(out, paste(setdiff(ar, x), collapse = "")); out }
  biogeo <- sim_map(phy, rs, rate_r, nb, root_state = "A")
  traits <- list(biogeo = biogeo)
  for (j in seq_len(n_normal)) traits[[paste0("trait", j)]] <- sim_map(phy, paste0("s", seq_len(trait_states)), rate_t)
  obs_r <- length(unique(unlist(lapply(biogeo$maps, names))))
  # missing taxa: target clade = 2 tips of a random internal node
  nodes <- sample((n_tips + 1):(n_tips + phy$Nnode), n_missing, replace = TRUE)
  tl <- lapply(nodes, function(nd) { d <- phytools::getDescendants(phy, nd); d <- d[d <= n_tips]; phy$tip.label[sample(d, min(2, length(d)))] })
  tl <- lapply(tl, function(x) if (length(x) < 2) c(x, x) else x)
  tl <- lapply(tl, unique)
  mtl <- list(Taxon_name = paste0("X", seq_len(n_missing)), Target_clade = stats::setNames(tl, paste0("X", seq_len(n_missing))))
  # allowed: a single-area range present in the map (broad), no restriction on traits
  allowed_r <- sort(unique(unlist(strsplit(unique(unlist(lapply(biogeo$maps, names))), ""))))
  res <- list(n_tips = n_tips, range_states = length(rs), observed_range_states = obs_r)
  tm <- function(expr) { t0 <- proc.time()[3]; v <- try(suppressWarnings(expr), silent = TRUE); if (inherits(v, "try-error")) cat("   ERROR:", substr(conditionMessage(attr(v, "condition")), 1, 200), "\n"); c(time = unname(proc.time()[3] - t0), failed = inherits(v, "try-error")) }
  # stage timings
  t_merge <- system.time(merged <- ns$merge_trait_simmaps(phy = {p <- phy; p}, trait_maps = structure(traits, class = "phylograftR_trait_maps")))[3]
  res$composite_states <- length(unique(unlist(lapply(merged$simmap$maps, names))))
  res$pieces <- sum(lengths(merged$simmap$maps))
  res$t_merge <- unname(t_merge)
  lb <- ns$label_backbone_nodes(phy)$phy
  tr <- ns$start_state_tracking(lb, ns$merge_trait_simmaps(lb, structure(traits, class = "phylograftR_trait_maps"))$simmap)
  res$t_resolve <- unname(system.time(r <- suppressWarnings(ns$resolve_trait_allowed_ranges(allowed_states = c(list(allowed_r), rep(list(NULL), n_normal)), traits = merged$traits, composite_states = tr$states, taxon_names = mtl$Taxon_name, range_trait = 1)))[3])
  res$inherited <- length(r$inherited_states)
  # full runs
  full <- function(sm, ...) tm(graft_missing_taxa(phy, mtl, simmap = sm, method = "uniform", verbose = FALSE, seed = 1, ...))
  res$t_full_ranges_traits <- full(traits, treat_states_as_ranges = c(TRUE, rep(FALSE, n_normal)), allowed_states = c(list(allowed_r), rep(list(NULL), n_normal)))["time"]
  res$fail_ranges_traits <- full(traits, treat_states_as_ranges = c(TRUE, rep(FALSE, n_normal)), allowed_states = c(list(allowed_r), rep(list(NULL), n_normal)))["failed"]
  res$t_full_ranges_only <- full(biogeo, treat_states_as_ranges = TRUE, allowed_states = allowed_r)["time"]
  res$t_full_traits_noranges <- full(traits, allowed_states = NULL)["time"]
  as.data.frame(res)
}
cfgs <- list(
  list(n_tips = 100, k_areas = 4, maxr = 2, n_missing = 50, trait_states = 3, n_normal = 1),
  list(n_tips = 100, k_areas = 5, maxr = 3, n_missing = 50, trait_states = 3, n_normal = 1),
  list(n_tips = 100, k_areas = 7, maxr = 3, n_missing = 50, trait_states = 3, n_normal = 1),
  list(n_tips = 100, k_areas = 7, maxr = 3, n_missing = 50, trait_states = 3, n_normal = 2),
  list(n_tips = 300, k_areas = 7, maxr = 3, n_missing = 100, trait_states = 3, n_normal = 2)
)
out <- do.call(rbind, lapply(cfgs, function(cf) { cat("running", unlist(cf), "\n"); do.call(run_case, cf) }))
cfg_df <- do.call(rbind, lapply(cfgs, as.data.frame)); print(cbind(cfg_df[c("n_tips","k_areas","maxr","n_missing","trait_states","n_normal")], out[-1]), digits = 3, row.names = FALSE)
