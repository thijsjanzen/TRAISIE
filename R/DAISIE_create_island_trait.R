#' @keywords internal
DAISIE_create_island_trait <- function(
    stt_table,
    total_time,
    island_spec,
    mainland,
    trait_pars,
    num_observed_states,
    num_hidden_states
) {
  # ---- helpers --------------------------------------------------------------
  cnames <- c("Species",
              "Mainland Ancestor",
              "Colonisation time (BP)",
              "Species type",
              "branch_code",
              "branching time (BP)",
              "Anagenetic_origin",
              "trait_state",
              "connection")

  coerce_island_spec <- function(x) {
    if (is.null(x) || length(x) == 0) {
      out <- as.data.frame(matrix(nrow = 0, ncol = length(cnames)))
      names(out) <- cnames
      return(out)
    }
    if (is.vector(x) && !is.list(x)) x <- matrix(x, nrow = 1)
    if (is.matrix(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (ncol(x) < length(cnames)) x[(ncol(x) + 1):length(cnames)] <- NA
    names(x) <- cnames
    x
  }

  mainland_total <- if (is.list(mainland)) sum(unlist(mainland)) else as.integer(mainland)

  root_from_mainland <- function(anc_id) {
    if (is.list(mainland)) {
      counts <- unlist(mainland)
      cum <- cumsum(counts)
      k <- length(counts)
      idx <- findInterval(as.numeric(anc_id), c(0, cum))
      one <- rep(0L, k)
      if (idx >= 1 && idx <= k) one[idx] <- 1L
      return(one)
    } else {
      return(1L) # single mainland pool
    }
  }

  # ---- coerce & sort --------------------------------------------------------
  isdf <- coerce_island_spec(island_spec)

  # sort by branch_code if present/numeric
  if ("branch_code" %in% names(isdf)) {
    o <- order(suppressWarnings(as.numeric(isdf$branch_code)), na.last = TRUE)
    isdf <- isdf[o, , drop = FALSE]
  }

  # ---- empty island: trivial return -----------------------------------------
  if (nrow(isdf) == 0) {
    return(list(
      list(island_age = total_time, not_present = mainland_total),
      list(branching_times = total_time, stac = 0, missing_species = 0, stt_table = stt_table)
    ))
  }

  # ---- times to BP (counting back from present) -----------------------------
  isdf[["branching time (BP)"]]    <- total_time - suppressWarnings(as.numeric(isdf[["branching time (BP)"]]))
  isdf[["Colonisation time (BP)"]] <- total_time - suppressWarnings(as.numeric(isdf[["Colonisation time (BP)"]]))

  # ---- single-colonist fast path --------------------------------------------
  if (mainland_total == 1L) {
    return(DAISIE:::DAISIE_ONEcolonist(total_time, isdf, stt_table))
  }

  # ---- build clades per ancestor & colonisation time ------------------------
  clades <- list()
  idx <- 0L

  mas <- sort(unique(as.numeric(isdf[["Mainland Ancestor"]])))
  for (ma in mas) {
    sub_ma <- isdf[as.numeric(isdf[["Mainland Ancestor"]]) == ma, , drop = FALSE]
    # distinct colonisation times for this ancestor
    col_times <- sort(unique(suppressWarnings(as.numeric(sub_ma[["Colonisation time (BP)"]]))))
    for (ct in col_times) {
      sub_ct <- sub_ma[suppressWarnings(as.numeric(sub_ma[["Colonisation time (BP)"]])) == ct, , drop = FALSE]
      if (!is.data.frame(sub_ct) || nrow(sub_ct) == 0) next

      idx <- idx + 1L
      res <- DAISIE:::DAISIE_ONEcolonist(
        total_time,
        island_spec = sub_ct,
        stt_table = NULL
      )

      # Attach extras safely
      res$root_state <- root_from_mainland(ma)
      res$sampling_fraction <- rep(1, num_observed_states)
      res$traits <- suppressWarnings(as.numeric(sub_ct$trait_state))

      if (nrow(sub_ct) > 1 && exists("build_phylo_tree_from_island_spec", mode = "function")) {
        res$phylogeny <- tryCatch(build_phylo_tree_from_island_spec(sub_ct),
                                  error = function(e) NA)
      } else {
        res$phylogeny <- NA
      }

      res$stt_table <- NULL
      clades[[idx]] <- res
    }
  }

  clades <- Filter(Negate(is.null), clades)

  island_header <- list(
    island_age  = total_time,
    not_present = max(0L, mainland_total - length(unique(mas)))
  )

  island <- append(list(island_header), clades)
  return(island)
}
