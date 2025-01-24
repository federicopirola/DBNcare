#' Removes the suffix _t-1 from the node's name
#'
#' @param node_name a string representing the name of the node
#'
#' @return the modified string
#' @export
#'
#' @examples remove_suffix_t_minus_1("A_t-1")
remove_suffix_t_minus_1 <- function(node_name) {
  if(is.character(node_name)==FALSE){
    stop("node_name must be a character")
  }
  gsub("-[0-9]+$", "", node_name)
}

#' Removes the suffix _t from the node's name
#'
#' @param node_name a string representing the name of the node
#'
#' @return the modified string
#' @export
#'
#' @examples remove_suffix_t("A_t")
remove_suffix_t <- function(node_name) {
  if(is.character(node_name)==FALSE){
    stop("node_name must be a character")
  }
  modified_node_name <- sub("_t$", "", node_name)
  return(modified_node_name)
}

#' Defines a level for every node in the network at time t
#'
#' @description
#' Calculates the levels in which the nodes will be distributed when plotting
#' the structure. This level is defined by their parent nodes: a node with no
#' parents will always be in the level 0. Subsequently, the level of a node
#' will be one more of the maximum level of its parents.
#' @param net an object of class DBN representing a dynamic bayesian network
#' @param order a topological order of the nodes at t, with the orphan nodes
#' in the first place. See \code{\link{node.ordering}}
#' @param lvl current level being processed
#' @param acc accumulator of the nodes already processed
#' @return a matrix with the names of the nodes in the first row and their
#' levels in the second
#' @export
node_levels <- function(net, order, lvl = 1, acc = NULL){
  if(class(net)!="DBN"){
    stop("net must be of class DBN")
  }
  ret <- acc
  if(length(order) > 0){
    clean_order <- remove_suffix_t(order[1])
  
    parent_n_order <-  net[['nodes']][[clean_order]][['t']][['parents']]
    parent_n_order <- select_nodes_names_with_t(parent_n_order)
    
    pa <- which(acc[1,] %in% parent_n_order)
  
    if(length(pa) > 0 && max(as.numeric(cbind(c("_","0"),acc[,pa])[2,])) == lvl)
      lvl <- lvl + 1
    ret <- node_levels(net, order[-1], lvl, cbind(acc, c(order[1], lvl)))
  }
  return(ret)
}

#' Defines a level for every node in the network g_0
#'
#' @description
#' Calculates the levels in which the nodes will be distributed when plotting
#' the structure. This level is defined by their parent nodes: a node with no
#' parents will always be in the level 0. Subsequently, the level of a node
#' will be one more of the maximum level of his parents.
#' @param net the object of class bn G_0 
#' @param order a topological order of the nodes, with the orphan nodes
#' in the first place. See \code{\link{node.ordering}}
#' @param lvl current level being processed
#' @param acc accumulator of the nodes already processed
#' @return a matrix with the names of the nodes in the first row and their
#' level in the second
#' @export
node_levels_g_0 <- function(net, order, lvl = 1, acc = NULL){
  if(class(net)!="bn"){
    stop("net must be of class bn")
  }
  ret <- acc
  if(length(order) > 0){
    pa <- which(acc[1,] %in% bnlearn::parents(net, order[1]))
    
    if(length(pa) > 0 && max(as.numeric(cbind(c("_","0"),acc[,pa])[2,])) == lvl)
      lvl <- lvl + 1
    ret <- node_levels_g_0(net, order[-1], lvl, cbind(acc, c(order[1], lvl)))
  }
  return(ret)
}


#' Plots the network G_0
#'
#' @param dbn a list representing a dbn
#'
#' @export
#'
plot_g0 <- function(dbn){
  library(dplyr)
  if(class(dbn)!="DBN"){
    stop("dbn must be of class DBN")
  }
  g_0 <- from_DBN_to_G_0(dbn)
  nodes_uniq <- bnlearn::node.ordering(g_0)
  nodes <- data.frame(id = nodes_uniq,
                      label = nodes_uniq,
                      #font.size = 24,
                      level = node_levels_g_0(g_0, nodes_uniq)[2,],
                      shadow = FALSE,
                      physics = FALSE)
  
  edges <- data.frame(from = bnlearn::arcs(g_0)[,1], #refactor using $arcs
                      to = bnlearn::arcs(g_0)[,2], #refactor using $arcs
                      arrows = "to",
                      #smooth = TRUE, # visNetwork's bug
                      shadow = FALSE
  )
  
  ret <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visHierarchicalLayout(levelSeparation = 100) %>%
    visNetwork::visOptions(nodesIdSelection = T)
  eval(ret)
}


#' Returns a vector with the number of consecutive nodes in each level
#'
#' @description
#' This method processes the vector of node levels to get the position of
#' each node inside the level. E.g. c(1,1,1,2,2,3,4,4,5,5) turns into 
#' c(1,2,3,1,2,1,1,2,1,2)
#' @param nodes a vector with the level of each node
#' @param res the accumulative results of the sub successions
#' @param prev the level of the previous node processed
#' @param acc the accumulator of the index in the current sub successions
#' @return the vector of sub successions in each level
#' @export
#' @keywords internal
acc_successions <- function(nodes, res = NULL, prev = 0, acc = 0){
  if(length(nodes) == 0)
    return(res)
  else if(prev == nodes[1])
    return(acc_successions(nodes[-1], c(res, acc+1), prev, acc+1))
  else
    return(acc_successions(nodes[-1], c(res, 1), nodes[1], 1))
}


#' Filters out all the nodes names that do not end up with t
#'
#' @param nodes_list a list containing nodes names
#'
#' @return a filtered list with nodes names that end up with t
#' @export
#'
select_nodes_names_with_t <- function(nodes_list) {
  if(is.character(nodes_list)==FALSE){
    stop("nodes_list must be a character!")
  }
  filtered <- nodes_list[grep("_t$", nodes_list)]
  return(filtered)
}

#' Filters out all the nodes names that do not end up with t-n
#'
#' @param nodes_list a list containing nodes names
#' @param n time point
#'
#' @export
#'
select_nodes_names_with_t_minus_n <- function(nodes_list, n) {
  if(is.character(nodes_list)==FALSE){
    stop("nodes_list must be a character!")
  }
  if(is.numeric(n)==FALSE){
    stop("n must be numeric!")
  }
  # Regular expression pattern to match strings ending with t-n 
  pattern <- paste0("\\w+_t-", n, "$")  
  filtered_strings <- grep(pattern, nodes_list, value = TRUE)
  return(filtered_strings)
}

#' Gets the topological ordered nodes list
#'
#' @description 
#' This method gets the structure of a DBN, isolates the nodes of the 
#' time slice t and then gives a topological ordering of them.
#' @param structure the structure of the network.
#' @return the ordered nodes at time t
#' 
#' @export
dynamic_ordering <- function(dbn){
  if(class(dbn)!="DBN"){
    stop("dbn must be of class DBN")
  }
  g_transition <- from_DBN_to_G_transition(dbn)
  nodes_t <- bnlearn::node.ordering(g_transition)
  nodes_t <- select_nodes_names_with_t(nodes_t)
  return(nodes_t)
}

#' Plots the network G_transition
#'
#' @param dbn a list representing a dbn 
#'
#' @export
#'
plot_g_transition <- function(dbn){
  library(dplyr)
  if(class(dbn)!="DBN"){
    stop("dbn must be of class DBN")
  }
  g_transition <- from_DBN_to_G_transition(dbn)
  ordered_node_at_t <- dynamic_ordering(dbn)
  nodes_at_t_levels <- node_levels(dbn,ordered_node_at_t)
  x_axis_positions <- acc_successions(as.numeric(nodes_at_t_levels[2,]))*100
  nodes_list <-  bnlearn::node.ordering(g_transition)
  nodes_list <- select_nodes_names_with_t(nodes_list)
  edges <- data.frame(g_transition$arcs)
  y_axis_positions <-  as.integer(nodes_at_t_levels[2,])*100
  nodes_df <- list()
  for(i in 1:length(nodes_list)){
    nodes_df[['id']] <- c(nodes_df[['id']],nodes_list[i])
    nodes_df[['label']] <- c(nodes_df[['label']],nodes_list[i])
    nodes_df[['x']] <- c(nodes_df[['x']],x_axis_positions[i])
    nodes_df[['y']] <-c(nodes_df[['y']],y_axis_positions[i])
  }
  nodes_t_minus_n <- bnlearn::node.ordering(g_transition)
  nodes_t_minus_n <- select_nodes_names_with_t_minus_n(nodes_t_minus_n, 1)
  nodes_df_tmp <- data.frame(nodes_df)
  for(i in 1:length(nodes_t_minus_n)){
    string_with_no_time <- remove_suffix_t_minus_1(nodes_t_minus_n[i])
    row_index <- nodes_df_tmp[, "label"] == string_with_no_time
    value_of_x <- nodes_df_tmp[row_index, "x"]-450 #set default offset taking into account  the max number of node for the same level
    value_of_y <- nodes_df_tmp[row_index, "y"]
    nodes_df[['id']] <- c(nodes_df[['id']],nodes_t_minus_n[i])
    nodes_df[['label']] <- c(nodes_df[['label']],nodes_t_minus_n[i])
    nodes_df[['x']] <- c(nodes_df[['x']], value_of_x)
    nodes_df[['y']] <- c(nodes_df[['y']], value_of_y)
  }
  nodes_df <- data.frame(nodes_df)
  nodes_df$physics <- FALSE
  
  visNetwork::visNetwork(nodes_df, edges, width = "100%", height = 1000) %>% 
    visNetwork::visEdges(arrows = "to")
  
}
