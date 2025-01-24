#' get_node_name_vector
#' @description
#' the function transforms the node name to the vector format A_t-1 -> c("A", "t-1")
#'
#' @param node_name a string representing the name of a node
#'
#' @return a vector representing the node name.
#' @export
#'
#' @examples get_node_name_vector("A_0")
get_node_name_vector <- function(node_name) {
  if (is.character(node_name) == FALSE) {
    stop("node_name must be a character!")
  }
  #getting last char of the node
  last_char_node_name <-
    substr(node_name, nchar(node_name), nchar(node_name))
  
  if (last_char_node_name == "1") {
    root_node_name <- get_generic_node_name_rex(node_name)
    return(c(root_node_name, "t-1"))
  }
  if (last_char_node_name == "t") {
    root_node_name <- get_generic_node_name_rex(node_name)
    return(c(root_node_name, "t"))
  }
  if (last_char_node_name == "0") {
    root_node_name <- get_generic_node_name_0(node_name)
    return(c(root_node_name, "t_0"))
  }
}

#' generate_dbn_nodes_names
#'
#' @description
#' the function returns the nodes names for g_0 and g_transition.
#'
#' @param nodes_names a vector containing the names of nodes.
#'
#' @return a list containing the nodes names for each time slice.
#' @export
#'
#' @examples generate_dbn_nodes_names(c("A", "B", "C"))
generate_dbn_nodes_names <- function(nodes_names) {
  if (!is.vector(nodes_names) ||
      !is.character(nodes_names) || length(nodes_names) <= 1) {
    stop("nodes_names must be a vector containing names of nodes!")
  }
  g_0_nodes_names <- c()
  g_t_1_nodes_names <- c()
  g_t_nodes_names <- c()
  
  for (n in nodes_names) {
    g_0_nodes_names <- c(g_0_nodes_names, paste(n, "_0", sep = ""))
    g_t_1_nodes_names <-
      c(g_t_1_nodes_names, paste(n, "_t-1", sep = ""))
    g_t_nodes_names <- c(g_t_nodes_names, paste(n, "_t", sep = ""))
    
  }
  list_of_nodes_names <-
    list(g_0 = g_0_nodes_names, g_t_1 = g_t_1_nodes_names, g_t = g_t_nodes_names)
  return(list_of_nodes_names)
}

#' generate_dbn_random_structure
#'
#' @description
#' Generates a random structure of a DBN. This function uses random.graph
#' from the bnlearn R package. Generates a DBN structure consistent with the
#' ordering of the nodes provided by the user; in each arc, the node at the
#' tail precedes the node at the head in the node set.
#' Arcs are samples independently with a probability of inclusion specified by
#' the prob_edge argument. For more details please check:
#' https://www.bnlearn.com/examples/dag/.
#'
#' @param nodes_names a vector containing the names of nodes.
#' @param is_same a boolean, if TRUE g_0 and g_transition will have the same
#' edges except for the temporal ones
#' @param prob_edge_intraslice a double parameter representing the probability of an edge to be present in the intraslice.
#' @param prob_edge_interslice a double parameter representing the probability of an edge to be present in the interslice.
#'
#' @return an object of class 'DBN'
#' @export
#'
#' @examples generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5, 0.4)
generate_dbn_random_structure <-
  function(nodes_names, is_same, prob_edge_intraslice, prob_edge_interslice) {
    if (!is.vector(nodes_names) ||
        !is.character(nodes_names) || length(nodes_names) <= 1) {
      stop("nodes_names must be a vector containing names of nodes!")
    }
    if (is.logical(is_same) == FALSE) {
      stop("is_same must be logical!")
    }
    if (!(0 < prob_edge_intraslice && prob_edge_intraslice <= 1) || !(0 < prob_edge_interslice && prob_edge_interslice <= 1) ) {
      stop("prob_edge must be greater than 0 and less than or equal to 1.")
    }
    list_of_nodes_names <- generate_dbn_nodes_names(nodes_names)
    g_0_nodes_names <- list_of_nodes_names[['g_0']]
    g_t_1_nodes_names <- list_of_nodes_names[['g_t_1']]
    g_t_nodes_names <- list_of_nodes_names[['g_t']]
    intra_slice_edges <-
      as.vector(t(bnlearn::arcs(
        bnlearn::random.graph(g_t_nodes_names, prob = prob_edge_intraslice)
      )))
    inter_slice_edges <- c()
    g_0_edges <- c()
    for (n_t_1 in g_t_1_nodes_names) {
      for (n_t in g_t_nodes_names) {
        add_edge <-
          sample(
            c(0, 1),
            size = 1,
            prob = c(1 - prob_edge_interslice, prob_edge_interslice),
            replace = FALSE
          )
        #if the arc between X_t and X_t-1 must be added put and if statement and draw an arc
        if (add_edge == 1) {
          inter_slice_edges <- c(inter_slice_edges, c(n_t_1, n_t))
        }
      }
    }
    if (is_same) {
      for (n_t in intra_slice_edges) {
        g_0_edge <- paste0(substr(n_t, 1, nchar(n_t) - 1), "0")
        g_0_edges <- c(g_0_edges, g_0_edge)
      }
    } else{
      g_0_edges <-
        as.vector(t(bnlearn::arcs(
          bnlearn::random.graph(g_0_nodes_names, prob = prob_edge_intraslice)
        )))
    }
    all_edges <- c(g_0_edges, intra_slice_edges, inter_slice_edges)
    
    generated_dbn <-
      empty_DBN(dynamic_nodes = nodes_names, markov_order = 1)
    
    for (i in seq(1, length(all_edges), by = 2)) {
      node_from <- all_edges[i]
      node_to <- all_edges[i + 1]
      generated_dbn <- add_arc_DBN(DBN = generated_dbn,
                                   from = get_node_name_vector(node_from),
                                   to = get_node_name_vector(node_to))
    }
    return(generated_dbn)
  }

#' generate_dbn_nodes_distributions
#'
#' @description
#' This function generates a set of multinomial distributions for each node in a network,
#' using the Dirichlet distribution to ensure uniform sampling of probability values.
#' The generated distributions can either have a fixed cardinality for all nodes or a
#' variable cardinality for each node, depending on the specified parameters.
#'
#'
#' @param generated_dbn an object of class DBN
#' @param fixed_cardinality A logical parameter (`TRUE` or `FALSE`). If `TRUE`,
#' all nodes will have the same cardinality as specified by
#' `max_variables_cardinality`. If `FALSE`, the cardinality for each node is randomly
#' determined, with a minimum of 2 and a maximum defined by `max_variables_cardinality`.
#' @param max_variables_cardinality An integer specifying the maximum cardinality to be used for the nodes.
#' If `fixed_cardinality` is `TRUE`, this parameter sets the exact cardinality for all nodes. If `fixed_cardinality`
#' is `FALSE`, this parameter defines the upper limit of the range from which the cardinality is randomly selected.
#'
#' @return a DBN.fit object
#' @export
#'
#' @examples generate_dbn_nodes_distributions(generated_dbn,c("A","B","C"), TRUE, 2)
generate_dbn_nodes_distributions <-
  function(generated_dbn,
           fixed_cardinality,
           max_variables_cardinality) {
    if (class(generated_dbn) != "DBN") {
      stop("generated_dbn must be a DBN object!")
    }
    if (!is.logical(fixed_cardinality)) {
      stop("fixed_cardinality must be TRUE or FALSE.")
    }
    if (max_variables_cardinality < 2) {
      stop("max_variables_cardinality must be greater than or equal to 2.")
    }
    nodes_names <- names(generated_dbn$nodes)
    list_of_nodes_names <- generate_dbn_nodes_names(nodes_names)
    g_0_nodes_names <- list_of_nodes_names[['g_0']]
    g_t_1_nodes_names <- list_of_nodes_names[['g_t_1']]
    g_t_nodes_names <- list_of_nodes_names[['g_t']]
    
    #generating levels for each variable
    list_of_levels <- list()
    #fixed_cardinality == TRUE all the nodes must have a fixed cardinality
    #the cardinality is equal to cardinality parameter
    if (fixed_cardinality == TRUE) {
      for (n in nodes_names) {
        list_of_levels[[n]] <- 0:(max_variables_cardinality - 1)
      }
    } else{
      #fixed_cardinality == FALSE for each node we will have a random cardinality
      # 2<= cardinality <=max_variables_cardinality
      for (n in nodes_names) {
        if (max_variables_cardinality == 2) {
          cardinality_n <- 1
        }
        else{
          cardinality_n <- sample(2:max_variables_cardinality, 1) - 1
        }
        list_of_levels[[n]] <- 0:cardinality_n
      }
    }
    #generate distribution for t_0
    list_cpt <- list()
    for (n_0 in g_0_nodes_names) {
      root_node_name <- get_generic_node_name_0(n_0)
      dim_names_n_0 <- list()
      cardinality_cpt_n_0 <- c()
      dim_names_n_0[[n_0]] <- list_of_levels[[root_node_name]]
      cardinality_cpt_n_0 <-
        c(cardinality_cpt_n_0, length(list_of_levels[[root_node_name]]))
      parents_n_0 <-
        generated_dbn[["nodes"]][[root_node_name]][["t_0"]][["parents"]]
      
      #adding each node's and parent's levels and creating the list of cardinality
      n_distributions <- 1
      
      for (p_n_0 in parents_n_0) {
        root_node_name_p <- get_generic_node_name_0(p_n_0)
        dim_names_n_0[[p_n_0]] <- list_of_levels[[root_node_name_p]]
        cardinality_cpt_n_0 <-
          c(cardinality_cpt_n_0, length(list_of_levels[[root_node_name_p]]))
        n_distributions <-
          n_distributions * length(list_of_levels[[root_node_name_p]])
      }
      k <- length(list_of_levels[[root_node_name]])
      alpha <-
        rep(1, k)  # Parameters for Dirichlet distribution, can be adjusted
      vector_random_distributions <- c()
      for (i in seq(1, n_distributions)) {
        random_probabilities <- MCMCpack::rdirichlet(1, alpha)
        distribution_n_0 <- as.vector(random_probabilities)
        vector_random_distributions <-
          c(vector_random_distributions, distribution_n_0)
      }
      list_cpt[[n_0]] <-
        array(vector_random_distributions,
              cardinality_cpt_n_0,
              dim_names_n_0)
    }
    #Generate distribution for t
    for (n_t in g_t_nodes_names) {
      root_node_name <- get_generic_node_name_rex(n_t)
      dim_names_n_t <- list()
      cardinality_cpt_n_t <- c()
      dim_names_n_t[[n_t]] <- list_of_levels[[root_node_name]]
      cardinality_cpt_n_t <-
        c(cardinality_cpt_n_t, length(list_of_levels[[root_node_name]]))
      parents_n_t <-
        generated_dbn[["nodes"]][[root_node_name]][["t"]][["parents"]]
      
      #adding each node's and parent's levels and creating the list of cardinality
      n_distributions <- 1
      
      for (p_n_t in parents_n_t) {
        root_node_name_p <- get_generic_node_name_rex(p_n_t)
        dim_names_n_t[[p_n_t]] <- list_of_levels[[root_node_name_p]]
        cardinality_cpt_n_t <-
          c(cardinality_cpt_n_t, length(list_of_levels[[root_node_name_p]]))
        n_distributions <-
          n_distributions * length(list_of_levels[[root_node_name_p]])
      }
      k <- length(list_of_levels[[root_node_name]])
      alpha <-
        rep(1, k)  # Parameters for Dirichlet distribution, can be adjusted
      vector_random_distributions <- c()
      for (i in seq(1, n_distributions)) {
        random_probabilities <- MCMCpack::rdirichlet(1, alpha)
        distribution_n_t <- as.vector(random_probabilities)
        vector_random_distributions <-
          c(vector_random_distributions, distribution_n_t)
      }
      list_cpt[[n_t]] <-
        array(vector_random_distributions,
              cardinality_cpt_n_t,
              dim_names_n_t)
    }
    CPTs <- list_cpt
    fitted_dbn <- DBN_parameters(DBN = generated_dbn, CPTs = list_cpt)
    
    return(fitted_dbn)
  }

#' random_dbn
#'
#' @description
#' Generates a random structure of a Dynamic Bayesian Network using the order
#' of the variable nodes_names
#'
#'
#'
#' @param nodes_names a list of strings representing the names of nodes in the
#' network
#' @param is_same TRUE if you want to keep the same edges defined in G_0
#' in G_transition. FALSE if you DO NOT want G_transition to have all the relationships
#' defined in G_0.
#' @param prob_edge_intraslice a double parameter representing the probability of an edge to be present in the intraslice.
#' @param prob_edge_interslice a double parameter representing the probability of an edge to be present in the interslice.
#'
#' @return a DBN object.
#' @export
#'
#' @examples random_dbn(c("A","B","C"), TRUE, 0.6, TRUE, 2)
random_dbn <- function(nodes_names, is_same, prob_edge_intraslice, prob_edge_interslice) {
  #generating the random structure of the dbn
  generated_dbn <-
    generate_dbn_random_structure(nodes_names, is_same, prob_edge_intraslice, prob_edge_interslice)
  return(generated_dbn)
}
#' fit_random_dbn
#'
#' @param generated_dbn a DBN class object
#' @param fixed_cardinality A logical parameter (`TRUE` or `FALSE`). If `TRUE`,
#' all nodes will have the same cardinality as specified by
#' `max_variables_cardinality`. If `FALSE`, the cardinality for each node is randomly
#' determined, with a minimum of 2 and a maximum defined by `max_variables_cardinality`.
#' @param max_variables_cardinality An integer specifying the maximum cardinality to be used for the nodes.
#' If `fixed_cardinality` is `TRUE`, this parameter sets the exact cardinality for all nodes. If `fixed_cardinality`
#' is `FALSE`, this parameter defines the upper limit of the range from which the cardinality is randomly selected.
#'
#' @return a DBN.fit object
#' @export
#'
#' @examples fit_random_dbn(generated_dbn, TRUE, 2)
fit_random_dbn <-
  function(generated_dbn,
           fixed_cardinality = FALSE,
           max_variables_cardinality) {
    #try to get nodes_names from the generated_dbn otherwise check that nodes_names is present in the dbn
    fitted_random_dbn <-
      generate_dbn_nodes_distributions(generated_dbn,
                                       fixed_cardinality,
                                       max_variables_cardinality)
    return(fitted_random_dbn)
  }