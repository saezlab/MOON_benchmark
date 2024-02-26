#' Decompress Moon Result
#'
#' This function decompresses the results obtained from moon analysis by incorporating 
#' node signatures and handling duplicated parents. It merges these details with the 
#' provided meta network data and returns a comprehensive data frame.
#'
#' @param moon_res A data frame containing the results of a moon analysis.
#' @param meta_network_compressed_list A list containing compressed meta network details, 
#'        including node signatures and duplicated parents.
#' @param meta_network A data frame representing the original meta network.
#'
#' @return A data frame which merges the moon analysis results with the meta network data,
#'         including additional details about node signatures and handling of duplicated parents.
#'
#' @examples
#' # Example usage (requires appropriate data structures for moon_res, 
#' # meta_network_compressed_list, and meta_network)
#' # decompressed_result <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network)
#'
#' @export
decompress_moon_result <- function(moon_res, meta_network_compressed_list, meta_network) {
  # Extract node_signatures and duplicated_parents from the list
  node_signatures <- meta_network_compressed_list$node_signatures
  duplicated_parents <- meta_network_compressed_list$duplicated_signatures
  
  # Create a dataframe for duplicated parents
  duplicated_parents_df <- data.frame(duplicated_parents)
  duplicated_parents_df$source_original <- row.names(duplicated_parents_df)
  names(duplicated_parents_df)[1] <- "source"
  
  # Create a dataframe for addons
  addons <- data.frame(names(node_signatures)[-which(names(node_signatures) %in% duplicated_parents_df$source_original)]) 
  names(addons)[1] <- "source"
  addons$source_original <- addons$source
  
  # Get final leaves
  final_leaves <- meta_network[!(meta_network$target %in% meta_network$source),"target"]
  final_leaves <- as.data.frame(cbind(final_leaves,final_leaves))
  names(final_leaves) <- names(addons)
  
  # Combine addons and final leaves
  addons <- as.data.frame(rbind(addons,final_leaves))
  
  # Create mapping table by combining duplicated parents and addons
  mapping_table <- as.data.frame(rbind(duplicated_parents_df,addons))
  
  mapping_table <- unique(mapping_table)
  # Merge the moon_res data frame with the mapping table
  moon_res <- merge(moon_res, mapping_table, by = "source")
  
  # Return the merged data frame
  return(moon_res)
}

prepare_for_roc = function(df, filter_tn = FALSE) {
  res = df %>%
    dplyr::mutate(response = case_when(.data$source == .data$target ~ 1,
                                       .data$source != .data$target ~ 0),
                  predictor = .data$score)
  res$response = factor(res$response, levels = c(1, 0))
  
  if (filter_tn == TRUE) {
    z = intersect(res$source, res$target)
    res = res %>%
      filter(.data$source %in% z, .data$target %in% z)
  }
  res %>%
    dplyr::select(.data$source, .data$id, .data$response, .data$predictor)
}

#' Get Moon Scoring Network
#'
#' This function analyzes a given meta network based on moon scores and an upstream node. 
#' It filters and processes the network by controlling and observing neighbours 
#' according to specified parameters. The function returns a list containing a filtered 
#' network and updated moon scores.
#'
#' @param upstream_node The node from which the network analysis starts.
#' @param meta_network The complete network data.
#' @param moon_scores Scores associated with each node in the network.
#' @param keep_upstream_node_peers Logical; whether to keep peers of the upstream node. Default is FALSE.
#' 
#' @return A list with two elements: 
#'   - `SIF`: A data frame representing the filtered meta network.
#'   - `ATT`: A data frame representing the updated moon scores.
#' 
#' @examples
#' # Example usage (requires appropriate data structures for meta_network and moon_scores)
#' # result <- get_moon_scoring_network(upstream_node, meta_network, moon_scores)
#' 
#' @export
get_moon_scoring_network <- function(upstream_node,
                                     meta_network,
                                     moon_scores,
                                     keep_upstream_node_peers = F){
  
  # Determine the number of steps from the upstream node based on moon score level.
  n_steps <- moon_scores[moon_scores$source == upstream_node,"level"]
  
  # If level peers of the upstream node are not to be kept, filter out these nodes.
  if(!keep_upstream_node_peers)
  {
    moon_scores <- moon_scores[!(moon_scores$level == n_steps & moon_scores$source != upstream_node),]
  }
  
  # Filter the meta network to keep only controllable neighbours of the upstream node.
  meta_network_filtered <- cosmosR:::keep_controllable_neighbours(network = meta_network, n_steps = n_steps,input_nodes = upstream_node)
  
  # Identify downstream inputs from the moon scores.
  downstream_inputs <- moon_scores[which(moon_scores$level == 0 & moon_scores$source %in% meta_network_filtered$target),"source"]
  
  # Further filter the network to keep only observable neighbours.
  meta_network_filtered <- cosmosR:::keep_observable_neighbours(network = meta_network_filtered, n_steps = n_steps,observed_nodes = downstream_inputs)
  
  # Update moon scores to include only those present in the filtered network.
  moon_scores <- moon_scores[moon_scores$source %in% meta_network_filtered$source |
                               moon_scores$source %in% meta_network_filtered$target,]
  
  # Filter the network to include only those connections present in moon scores.
  meta_network_filtered <- meta_network_filtered[meta_network_filtered$source %in% moon_scores$source &
                                                   meta_network_filtered$target %in% moon_scores$source,]
  
  # Return a list containing the filtered network (SIF) and the updated moon scores (ATT).
  return(list("SIF" = meta_network_filtered, "ATT" = moon_scores))
}