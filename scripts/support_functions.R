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

get_moon_scoring_network <- function(upstream_node,
                                 meta_network,
                                 moon_scores,
                                 keep_upstream_node_peers = F){
  
  
  n_steps <- moon_scores[moon_scores$source == upstream_node,"level"]
  
  if(!keep_upstream_node_peers)
  {
    moon_scores <- moon_scores[!(moon_scores$level == n_steps & moon_scores$source != upstream_node),]
  }
  
  meta_network_filtered <- cosmosR:::keep_controllable_neighbours(network = meta_network, n_steps = n_steps,input_nodes = upstream_node)
  downstream_inputs <- moon_scores[which(moon_scores$level == 0 & moon_scores$source %in% meta_network_filtered$target),"source"]
  meta_network_filtered <- cosmosR:::keep_observable_neighbours(network = meta_network_filtered, n_steps = n_steps,observed_nodes = downstream_inputs)
  
  moon_scores <- moon_scores[moon_scores$source %in% meta_network_filtered$source |
                               moon_scores$source %in% meta_network_filtered$target,]
  
  
  meta_network_filtered <- meta_network_filtered[meta_network_filtered$source %in% moon_scores$source &
                                                   meta_network_filtered$target %in% moon_scores$source,]
  
  return(list("SIF" = meta_network_filtered, "ATT" = moon_scores))
}
