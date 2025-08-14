#' Reconstruct a Phylogenetic Tree from DAISIE Trait-Dependent Simulation Output
#'
#' This function reconstructs a phylogenetic tree (as an \code{ape::phylo} object) from a matrix
#' of species connection paths and branching times produced by a trait-dependent DAISIE simulation.
#' It processes string-encoded ancestral paths, aligns them with branching times, builds parent-child
#' edges, and constructs a valid rooted ultrametric tree.
#'
#' @param island_spec_matrix A character vector where each element encodes the ancestry of a species,
#'        typically taken from the \code{connection} column of the island specification matrix. Each
#'        string is a space-separated path of branching times (e.g., \code{"NA 1.7 1.3"}).
#' @param branch_codes A character vector giving the labels for each tip of the tree, typically
#'        corresponding to the \code{branch_code} of each species (e.g., \code{c("AA", "B", "AB")}).
#' @param traits_vec A character or numeric vector indicating the trait state of each species,
#'        typically from column 8 of the island spec matrix. Not used directly, but may be useful
#'        for downstream annotation.
#' @param branching_times A numeric or character vector of all valid branching times in the island.
#'        Used to validate and clean each ancestral path.
#'
#' @return A \code{phylo} object (from the \pkg{ape} package), containing:
#' \itemize{
#'   \item \code{edge}: A matrix of parent-child node pairs.
#'   \item \code{edge.length}: Numeric vector of branch lengths.
#'   \item \code{tip.label}: Vector of labels for each tip (species).
#'   \item \code{Nnode}: Number of internal nodes.
#' }
#'
#' @details
#' - Each species' ancestry is encoded as a path of branching times. These paths are first cleaned
#'   by removing non-matching times.
#' - Nodes are assigned unique integer labels and parent-child relationships are extracted.
#' - Edge lengths are computed as differences between successive times.
#' - Duplicate edges and non-matching times are removed.
#' - Node numbers are remapped to meet the \code{ape} standard (tips first, then internal nodes).
#'
#' The function assumes the input paths are sorted from root to tip. Tip labels must match
#' the order of \code{island_spec_matrix}.
#'
#'

build_phylo_tree_from_island_spec <- function(island_spec) {


  # Ensure required package is loaded
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The 'ape' package is required but not installed. Please install it with install.packages('ape').")
  }
  island_spec <- island_spec[order(island_spec[, 5]), ]
  # Convert colonisation time column to numeric
  colonisation_times <- as.numeric( island_spec [, "Colonisation time (BP)"])

  # Check if there are at least 2 unique colonisation times
  if (length(unique(colonisation_times)) > 1) {
    # Find the smallest colonisation time
    min_col_time <- min(colonisation_times, na.rm = TRUE)

    # Keep only the rows where colonisation time is NOT the smallest
    island_spec  <-  island_spec [colonisation_times != min_col_time, , drop = FALSE]
  }



  # Inputs
  history_vec = island_spec[, 9]  # e.g., character vector like c("NA 1.7 1.3", ...)
  bt_vec = island_spec[, 6]
  col_time = island_spec[, 3] # e.g., numeric or character vector
  branch_codes = island_spec[,5]
  traits_vec = island_spec[, 8]


  # Convert bt_vec to numeric for comparison
  bt_vec_num <- suppressWarnings(as.numeric(bt_vec))
  tolerance <- 1e-8

  # Initialize cleaned vector
  cleaned_history <- character(length(history_vec))

  # Loop through each element of history_vec
  for (i in seq_along(history_vec)) {

    # Split the string into elements
    vals <- unlist(strsplit(history_vec[i], " "))

    # Initialize list to keep valid values
    valid_vals <- character(0)

    for (val in vals) {
      # Keep "NA" as-is
      if (val == "NA") {
        valid_vals <- c(valid_vals, "NA")
      } else {
        num_val <- suppressWarnings(as.numeric(val))
        # Keep if close to any value in bt_vec_num
        if (!is.na(num_val)) {
          if (any(abs(num_val - bt_vec_num) < tolerance, na.rm = TRUE)) {
            valid_vals <- c(valid_vals, val)
          }
        }
      }
    }

    # Join the valid values back into a string
    cleaned_history[i] <- paste(valid_vals, collapse = " ")
  }
  # Step 2: Create value → label map
  all_vals <- unique(unlist(strsplit(cleaned_history, " ")))
  all_vals <- all_vals[all_vals != "NA" & all_vals != ""]
  label_map <- setNames(seq_along(all_vals), all_vals)

  # Step 3: Recode each row with unique terminal node
  used_numbers <- c()
  next_label <- max(label_map) + 1

  # Step 3: Recode each row with unique terminal node
  translated_result <- vector("character", length(cleaned_history))
  used_numbers <- c()
  next_label <- max(label_map) + 1

  for (i in seq_along(cleaned_history)) {
    row_str <- cleaned_history[i]
    vals <- unlist(strsplit(row_str, " "))
    vals <- vals[vals != "NA" & vals != ""]

    labels <- as.character(label_map[vals])
    current_numbers <- as.numeric(labels)

    while (next_label %in% current_numbers || next_label %in% used_numbers) {
      next_label <- next_label + 1
    }

    final_numbers <- c(current_numbers, next_label)
    used_numbers <- c(used_numbers, next_label)
    next_label <- next_label + 1

    translated_result[i] <- paste(final_numbers, collapse = " ")
  }

  # Final output
  translated_history <- translated_result

  # Step 4: Create node vectors
  vecs <- lapply(translated_history, function(str_row) {
    as.integer(unlist(strsplit(str_row, " ")))
  })

  # Step 5: Extract edges
  edges <- do.call(rbind, lapply(vecs, function(path) {
    if (length(path) < 2) return(NULL)
    cbind(parent = path[-length(path)], child = path[-1])
  }))

  # Step 6: Identify and remove duplicate edges
  duplicated_rows <- duplicated(edges)
  duplicate_indices <- which(duplicated_rows)
  edges <- edges[!duplicated_rows, ]

  # Step 7: Compute edge lengths
  time_vecs_with_zero <- lapply(cleaned_history, function(row_str) {
    vals <- unlist(strsplit(as.character(row_str), " "))
    vals <- vals[vals != "NA" & vals != ""]
    numeric_vals <- as.numeric(vals)
    c(numeric_vals, 0)
  })

  edge_lengths <- unlist(lapply(time_vecs_with_zero, function(times) {
    if (length(times) < 2) return(numeric(0))
    head(times, -1) - tail(times, -1)
  }))

  if (length(duplicate_indices) > 0) {
    edge_lengths <- edge_lengths[-duplicate_indices]
  }

  # Step 8: Remap edges to comply with ape format
  all_nodes <- unique(as.vector(edges))
  parent_nodes <- unique(edges[, 1])
  child_nodes  <- unique(edges[, 2])
  tip_nodes <- sort(setdiff(child_nodes, parent_nodes))
  ntip <- length(tip_nodes)
  internal_nodes <- sort(setdiff(all_nodes, tip_nodes))

  node_order <- c(tip_nodes, internal_nodes)
  node_map <- setNames(seq_along(node_order), node_order)

  edge_remapped <- matrix(
    c(node_map[as.character(edges[, 1])], node_map[as.character(edges[, 2])]),
    ncol = 2
  )
  # branch_codes <- as.vector(sim[[2]][[22]]$island_spec[[1]][, 5])
  # traits_vec <- as.vector(sim[[2]][[22]]$island_spec[[1]][, 8])
  # Step 9: Build phylo object
  tree <- list()
  tree$edge <- matrix(as.integer(edge_remapped), ncol = 2)
  tree$edge.length <- as.numeric(edge_lengths)
  tree$Nnode <- length(internal_nodes)
  tree$tip.label <- branch_codes
  class(tree) <- "phylo"

  # Optional validation
  # ape::checkValidPhylo(tree)

  return(tree)
}



