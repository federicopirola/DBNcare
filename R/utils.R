#' rename_nodes_unroll
#' @description
#' Renames a list of nodes, making nodes names compliant to the dbn unrolled structure
#' 
#' @param time_slice an integer representing a time slice
#' @param nodes_names a list of nodes names 
#'
#' @export
#'
rename_nodes_unroll <- function(time_slice, nodes_names) {
  renamed_nodes <- c()
  for (name_n in nodes_names) {
    #does the node end with t-1
    end_t_1 <- grepl("_t-1$", name_n)
    end_with_t <- grepl("_t$", name_n)
    if (end_t_1) {
      previous_time_slice <- time_slice - 1
      root_node_name <- get_generic_node_name_rex(name_n)
      new_node_name <-
        paste(root_node_name,
              "_",
              as.character(previous_time_slice),
              sep = "")
    }
    if (end_with_t) {
      root_node_name <- get_generic_node_name_rex(name_n)
      new_node_name <-
        paste(root_node_name, "_", as.character(time_slice), sep = "")
    }
    renamed_nodes <- c(renamed_nodes, new_node_name)
  }
  return(renamed_nodes)
}

#' get_node_edges
#' @description
#' This function gets all the edges of a node
#' 
#' @param dbn_transition a dbn transition
#' @param node name of the node in the dbn transition
#' @param time_slice an integer representing a time slice
#'
#' @return all the edges of a specific node in a time slice
#' @export
#'
get_node_edges <- function(dbn_transition, node, time_slice) {
  if (class(dbn_transition) != "bn.fit") {
    stop("dbn must be a bn.fit object!")
  }
  if (is.character(node) == FALSE) {
    stop("node must be a character!")
  }
  if (is.numeric(time_slice)==FALSE){
    stop("time slice must be an integer!")
  }
  node_edges <- c()
  root_node_name_n_t <- get_generic_node_name_rex(node)
  new_node_name_n_t <-
    paste(root_node_name_n_t, "_", as.character(time_slice), sep = "")
  node_parents <- get_parents(dbn_transition, node)
  
  if (length(node_parents) > 0) {
    renamed_parents <- rename_nodes_unroll(time_slice, node_parents)
    
    for (n_renamed_parent in renamed_parents) {
      #using "from" "to" standard
      single_edge <- c(n_renamed_parent, new_node_name_n_t)
      node_edges <- c(node_edges, single_edge)
    }
  }
  return(node_edges)
}

#' get_unrolled_dbn
#' 
#' @description 
#' this function generates an unrolled dynamic bayesian network
#' @param dbn_fitted a dbn.fit class object
#' @param slices number of time-slice of the unrolled network
#'
#' @return a bn.fit object 
#' @export
#'
#' @examples get_unrolled_dbn(my_fitted_dbn, 4)
get_unrolled_dbn <- function(dbn_fitted, slices) {
  if (is.character(slices)) {
    stop("slices must be an integer!")
  }
  if (slices < 1){
    stop("slices must be greater than 0!")
  }
  if (class(dbn_fitted) != "dbn.fit") {
    stop("dbn_fitted must be a dbn.fit object")
  }
  nodes_bn <- c()
  edges_bn <- c()
  cpt_bn <- list()
  my_dbn_transition <-
    from_fitted_DBN_to_fitted_G_transition(dbn_fitted)
  nodes_in_transition <- get_nodes_t(my_dbn_transition)
  my_dbn_g_0 <- from_fitted_DBN_to_fitted_G_0(dbn_fitted)
  nodes_in_0 <- bnlearn::node.ordering(my_dbn_g_0)
  #iterate on nodes and get,for each one, edges and probability tables
  for (n_0 in nodes_in_0) {
    #add the node in bn_nodes
    nodes_bn <- c(nodes_bn, n_0)
    cpt_bn[[n_0]] <- my_dbn_g_0[[n_0]][["prob"]]
    parents_n_0 <- my_dbn_g_0[[n_0]][["parents"]]
    if (length(parents_n_0) > 0) {
      for (p_n_0 in parents_n_0) {
        #using "from" "to" standard
        edge_n_0 <- c(p_n_0, n_0)
        edges_bn <- c(edges_bn, edge_n_0)
      }
    }
  }
  for (time_slice in seq(slices)) {
    for (n_t in nodes_in_transition) {
      cpt_n_t <- my_dbn_transition[[n_t]][["prob"]]
      names_cpt_n_t <- names(dimnames(cpt_n_t))
      #rename each single node if node has t rename t with time_slice
      #if name has t-1 subtract 1 to time_slice.
      renamed_nodes <- rename_nodes_unroll(time_slice, names_cpt_n_t)
      root_node_name_n_t <- get_generic_node_name_rex(n_t)
      new_node_name_n_t <-
        paste(root_node_name_n_t, "_", as.character(time_slice), sep = "")
      names(dimnames(cpt_n_t)) <- renamed_nodes
      cpt_bn[[new_node_name_n_t]] <- cpt_n_t
      #adding node
      nodes_bn <- c(nodes_bn, new_node_name_n_t)
      #adding edges
      edges_time_slice <-
        get_node_edges(my_dbn_transition, n_t, time_slice)
      if (length(edges_time_slice) > 0) {
        edges_bn <- c(edges_bn, edges_time_slice)
        
      }
    }
  }
  dag_unrolled <- bnlearn::empty.graph(nodes = nodes_bn)
  arc_unrolled.set <- matrix(
    edges_bn,
    byrow = TRUE,
    ncol = 2,
    dimnames = list(NULL, c("from", "to"))
  )
  bnlearn::arcs(dag_unrolled) <- arc_unrolled.set
  dbn_unrolled <- bnlearn::custom.fit(dag_unrolled, cpt_bn)
  return(dbn_unrolled)
}
