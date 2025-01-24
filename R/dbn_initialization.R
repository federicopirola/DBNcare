library(dplyr)
# CREATING AN EMPTY DBN GIVEN THE SETS OF STATIC NODES AND DYNAMIC NODES AND THE MARKOVIAN ORDER

#' Function to create a Dynamic Bayesian Network
#'
#' @param static_nodes character vector listing the time-independent nodes
#' @param dynamic_nodes character vector listing the time-dependent nodes
#' @param markov_order integer value (> 1): markovian order of the process
#'
#' @return object of class 'DBN'
#' @export
#'
#' @examples
#' DBN_example <- empty_DBN(static_nodes = c("K"), dynamic_nodes = c("A", "S", "E", "O", "R", "T", "B"), markov_order = 1)
#' class(DBN_example) # 'DBN'
empty_DBN <-
  function(static_nodes = c(),
           dynamic_nodes,
           markov_order) {
    if (is.null(static_nodes)) {
      if (length(intersect(static_nodes, dynamic_nodes)) != 0) {
        stop('A node could not be both Static and Dynamic!!!')
      }
    }
    if (markov_order < 1) {
      stop('Markov order of the process must be 1 or higher!!!')
    }
    black_list_fc <-
      function(static_nodes,
               dynamic_nodes,
               markov_order) {
        if (length(c(dynamic_nodes, static_nodes)) > 0) {
          perm <-
            gtools::permutations(length(c(dynamic_nodes, static_nodes)),
                                 2,
                                 c(dynamic_nodes, static_nodes),
                                 repeats.allowed = TRUE)
          blacklist <-
            matrix(
              ncol = 2,
              nrow = 0,
              byrow = TRUE,
              dimnames = list(character(0), c("from", "to"))
            )
          for (i in 1:nrow(perm)) {
            from <- perm[i, 1]
            to <- perm[i, 2]
            if (from %in% dynamic_nodes & to %in% dynamic_nodes) {
              blacklist <-
                rbind(blacklist, c(gsub(" ", "", paste(
                  from, '_0'
                )), gsub(" ", "", paste(to, '_t'))))
              blacklist <-
                rbind(blacklist, c(gsub(" ", "", paste(
                  from, '_t'
                )), gsub(" ", "", paste(to, '_0'))))
              for (i in 1:markov_order) {
                blacklist <-
                  rbind(blacklist, c(gsub(" ", "", paste(
                    from, '_t-', i
                  )), gsub(" ", "", paste(
                    to, '_0'
                  ))))
                blacklist <-
                  rbind(blacklist, c(gsub(" ", "", paste(
                    from, '_0'
                  )), gsub(" ", "", paste(
                    to, '_t-', i
                  ))))
                blacklist <-
                  rbind(blacklist, c(gsub(" ", "", paste(
                    from, '_t'
                  )), gsub(" ", "", paste(
                    to, '_t-', i
                  ))))
                for (j in 1:markov_order) {
                  blacklist <-
                    rbind(blacklist, c(gsub(
                      " ", "", paste(from, '_t-', i)
                    ), gsub(" ", "", paste(
                      to, '_t-', j
                    ))))
                }
              }
            }
            else if (from %in% dynamic_nodes &
                     to %in% static_nodes) {
              blacklist <-
                rbind(blacklist, c(gsub(" ", "", paste(
                  from, '_0'
                )), gsub(" ", "", paste(to, '_0'))))
              blacklist <-
                rbind(blacklist, c(gsub(" ", "", paste(
                  from, '_t'
                )), gsub(" ", "", paste(to, '_0'))))
              for (i in 1:markov_order) {
                blacklist <-
                  rbind(blacklist, c(gsub(" ", "", paste(
                    from, '_t-', i
                  )), gsub(" ", "", paste(
                    to, '_0'
                  ))))
              }
            }
            else if (from %in% static_nodes &
                     to %in% dynamic_nodes) {
              blacklist <-
                rbind(blacklist, c(gsub(" ", "", paste(
                  from, '_0'
                )), gsub(" ", "", paste(to, '_t'))))
              for (i in 1:markov_order) {
                blacklist <-
                  rbind(blacklist, c(gsub(" ", "", paste(
                    from, '_0'
                  )), gsub(" ", "", paste(
                    to, '_t-', i
                  ))))
              }
            }
          }
        }
        blacklist
      }
    DBN = list(
      learning = list(
        whitelist = NULL,
        blacklist = black_list_fc(static_nodes, dynamic_nodes, markov_order),
        test = 'none',
        ntests = 0,
        algo = 'empty',
        args = list()
      ),
      markov_order = markov_order,
      arcs = matrix(
        ncol = 2,
        nrow = 0,
        byrow = TRUE,
        dimnames = list(character(0), c("from", "to"))
      ),
      nodes = list()
    )
    details = list(
      mb = character(0),
      nbr = character(0),
      parents = character(0),
      children = character(0)
    )
    details_2 = list(type = 'Dynamic', 't_0' = details)
    for (i in 0:markov_order) {
      details_2[[ifelse(i == 0, 't', gsub(" ", "", paste("t-", i)))]] = details
    }
    for (i in dynamic_nodes) {
      DBN$nodes[[i]] <- details_2
    }
    for (i in static_nodes) {
      DBN$nodes[[i]][['type']] = 'Static'
      DBN$nodes[[i]][['t_0']] <- details
    }
    class(DBN) <- "DBN"
    DBN
  }

# ADDING NODE TO DBN

#' Function for node addiction in Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#' @param node object of class 'character': name of the node to add
#' @param type 'Dynamic' (default) or 'Static' node
#'
#' @return object of class 'DBN'
#' @export
#'
#' @examples
#' DBN_example <- add_node_DBN(DBN=DBN_example, node="F", type='Dynamic')
#' DBN_example <- add_node_DBN(DBN=DBN_example, node="K", type='Static')
add_node_DBN <- function(DBN, node, type = 'Dynamic') {
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  if (!is.character(node))
    stop("ERROR: node is not a character")
  if (node %in% names(DBN$nodes)) {
    stop("ERROR: node name already exists")
  }
  if (type == 'Dynamic') {
    details_2 = list(
      type = 'Dynamic',
      't_0' = list(
        mb = character(0),
        nbr = character(0),
        parents = character(0),
        children = character(0)
      )
    )
    for (i in 0:DBN$markov_order) {
      details_2[[ifelse(i == 0, 't', gsub(" ", "", paste("t-", i)))]] = list(
        mb = character(0),
        nbr = character(0),
        parents = character(0),
        children = character(0)
      )
    }
    DBN$nodes[[node]] <- details_2
    
  }
  else if (type == 'Static') {
    DBN$nodes[[node]][['type']] = 'Static'
    DBN$nodes[[node]][['t_0']] <-
      list(
        mb = character(0),
        nbr = character(0),
        parents = character(0),
        children = character(0)
      )
  }
  else{
    stop("ERROR: type must be 'Dynamic' or 'Static")
  }
  DBN
}


# ADDING ARC FUNCTION

#' Function for arc addiction in Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'DBN'
#' @export
#'
#' @examples
#' DBN_example <- add_arc_DBN(DBN=DBN_example,from=c('A','t_0'),to=c('R','t_0'))
#' DBN_example <- add_arc_DBN(DBN=DBN_example,from=c('S','t'),to=c('O','t'))
#' DBN_example <- add_arc_DBN(DBN=DBN_example,from=c('S','t-1'),to=c('S','t'))
add_arc_DBN <- function(DBN, from, to, cycle_OK = TRUE) {
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  if (!is.character(from))
    stop("ERROR: from is not a character")
  if (!is.character(to))
    stop("ERROR: to is not a character")
  if (!(length(from) == 2 & length(to) == 2)) {
    stop("ERROR: from and to must be vectors of type (variable, time)")
  }
  
  if (!(from[1] %in% names(DBN$nodes) &
        to[1] %in% names(DBN$nodes))) {
    stop("ERROR: Defined node not in DBN!!!")
  }
  
  if ((from[1] == to[1] & from[2] == to[2])) {
    stop("ERROR: The arc defines a loop!!!")
  }
  
  if (to[2] == 't_0') {
    if (from[2] != 't_0') {
      stop("ERROR: Prior arcs must be between nodes at 't_0'")
    }
    else{
      from_ <- gsub(" ", "", paste(from[1], "_", 0))
      to_ <- gsub(" ", "", paste(to[1], "_", 0))
    }
  }
  else if (to[2] == 't') {
    if (!from[2] == 't' &
        (nchar(from[2]) <= 2 |
         !(substring(from[2], 1, 2) == 't-' &
           !is.na(as.numeric(
             substring(from[2], 3)
           ))))) {
      stop("ERROR: Parent nodes in DBN must be at time 't_0', 't' or 't-k' (with k <= Markov Order)")
    }
    else if (!from[2] == 't' &
             as.numeric(substring(from[2], 3)) > DBN$markov_order) {
      stop("ERROR: Transition arcs could not be of order higher than DBN's Markov Order")
    }
    else{
      from_ <- gsub(" ", "", paste(from[1], "_", from[2]))
      to_ <- gsub(" ", "", paste(to[1], "_", to[2]))
    }
  }
  else{
    stop("ERROR: Children nodes in DBN must be at time 't' or 't_0'")
  }
  
  if ((DBN[['nodes']][[from[1]]][['type']] == 'Dynamic' &
       DBN[['nodes']][[to[1]]][['type']] == 'Static')) {
    stop("ERROR: A Dynamic Node could not be parent of a Static Node!!!")
  }
  
  if ((from[2] != 't_0' &
       DBN[['nodes']][[from[1]]][['type']] == 'Static')) {
    stop(
      "ERROR: Static Node have not order higher than 0, thus they are not included in transition arcs!!!"
    )
  }
  
  if ((to[2] != 't_0' &
       DBN[['nodes']][[to[1]]][['type']] == 'Static')) {
    stop(
      "ERROR: Static Node have not order higher than 0, thus they are not included in transition arcs!!!"
    )
  }
  
  if (cycle_OK == FALSE &
      !identical(DBN[['nodes']][[from[1]]][[from[2]]][['parents']], character(0)) &
      !identical(DBN[['nodes']][[to[1]]][[to[2]]][['children']], character(0))) {
    chi <- DBN[['nodes']][[to[1]]][[to[2]]][['children']]
    s = TRUE
    while (s) {
      if (!identical(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['parents']], chi), character(0)) &
          !is.null(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['parents']], chi))) {
        stop("ERROR: The arc create a cycle!!!")
      }
      else{
        chi1 <- c()
        for (j in chi) {
          j_star <- strsplit(j, '_')[[1]]
          j_1 <-
            paste(j_star[1:(length(j_star) - 1)], collapse = '_')
          temp_j <-
            ifelse(j_star[length(j_star)] == '0', 't_0', j_star[length(j_star)])
          chi1 <-
            c(chi1, DBN[['nodes']][[j_1]][[temp_j]][['children']])
        }
        chi <- chi1
        if (identical(chi, character(0)) | is.null(chi)) {
          s = FALSE
        }
      }
    }
  }
  
  if (!to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']]) {
    DBN[['nodes']][[from[1]]][[from[2]]][['children']] <-
      c(DBN[['nodes']][[from[1]]][[from[2]]][['children']], to_)
    DBN[['nodes']][[to[1]]][[to[2]]][['parents']] <-
      c(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], from_)
    if (!to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['nbr']]) {
      DBN[['nodes']][[from[1]]][[from[2]]][['nbr']] <-
        c(DBN[['nodes']][[from[1]]][[from[2]]][['nbr']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['nbr']] <-
        c(DBN[['nodes']][[to[1]]][[to[2]]][['nbr']], from_)
    }
    if (!to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['mb']]) {
      DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
        c(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
        c(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_)
    }
    if (!is.na(DBN$arcs[DBN$arcs[, 'to'] == to_,]['from'])) {
      for (i in DBN$arcs[DBN$arcs[, 'to'] == to_,]['from']) {
        i_star <- strsplit(i, '_')[[1]]
        i_1 <- paste(i_star[1:(length(i_star) - 1)], collapse = '_')
        temp_i <-
          ifelse(i_star[length(i_star)] == '0', 't_0', i_star[length(i_star)])
        if (!i %in% DBN[['nodes']][[from[1]]][[temp_i]][['mb']]) {
          DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
            c(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], i)
          DBN[['nodes']][[i_1]][[temp_i]][['mb']] <-
            c(DBN[['nodes']][[i_1]][[temp_i]][['mb']], from_)
        }
      }
    }
    DBN[['arcs']] <- rbind(DBN$arcs, c(from_, to_))
  }
  
  DBN
}

# REMOVING ARC FUNCTION

#' Function for arc deletion in Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'DBN'
#' @export
#'
#' @examples
#' DBN_example <- delete_arc_DBN(DBN=DBN_example,from=c('A','t_0'),to=c('R','t_0'))
#' DBN_example <- delete_arc_DBN(DBN=DBN_example,from=c('S','t'),to=c('O','t'))
#' DBN_example <- delete_arc_DBN(DBN=DBN_example,from=c('S','t-1'),to=c('S','t'))
delete_arc_DBN <- function(DBN, from, to) {
  if (!class(DBN) == 'DBN')
    stop("Error: DBN argument is not of class 'DBN'")
  if (!is.character(from))
    stop("Error: from is not a character")
  if (!is.character(to))
    stop("Error: to is not a character")
  delete_element <- function(x, element) {
    if (class(x)[1] == 'character') {
      x[!x == element]
    }
    else{
      if (dim(x)[1] != 2) {
        x[-prodlim::row.match(element, x),]
      }
      else{
        matrix(
          x[-prodlim::row.match(element, x),],
          nrow = 1,
          ncol = 2,
          dimnames = list(NULL, c('from', 'to'))
        )
      }
    }
  }
  from_ <-
    gsub(" ", "", paste(from[1], "_", ifelse(from[2] == 't_0', '0', from[2])))
  to_ <-
    gsub(" ", "", paste(to[1], "_", ifelse(to[2] == 't_0', '0', to[2])))
  if (to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']]) {
    DBN[['arcs']] <- delete_element(DBN$arcs, c(from_, to_))
    DBN[['nodes']][[from[1]]][[from[2]]][['children']] <-
      delete_element(DBN[['nodes']][[from[1]]][[from[2]]][['children']], to_)
    DBN[['nodes']][[to[1]]][[to[2]]][['parents']] <-
      delete_element(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], from_)
    if (!to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['parents']]) {
      DBN[['nodes']][[from[1]]][[from[2]]][['nbr']] <-
        delete_element(DBN[['nodes']][[from[1]]][[from[2]]][['nbr']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['nbr']] <-
        delete_element(DBN[['nodes']][[to[1]]][[to[2]]][['nbr']], from_)
      DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
        delete_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
        delete_element(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_)
    }
    else if (identical(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['children']], DBN[['nodes']][[to[1]]][[to[2]]][['children']]), character(0))) {
      DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
        delete_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], to_)
      DBN[['nodes']][[to[1]]][[to[2]]][['mb']] <-
        delete_element(DBN[['nodes']][[to[1]]][[to[2]]][['mb']], from_)
    }
    if (!identical(DBN[['nodes']][[to[1]]][[to[2]]][['parents']], character(0))) {
      for (i in DBN[['nodes']][[to[1]]][[to[2]]][['parents']]) {
        i_star <- strsplit(i, '_')[[1]]
        i_1 <- paste(i_star[1:(length(i_star) - 1)], collapse = '_')
        temp_i <-
          ifelse(i_star[length(i_star)] == '0', 't_0', i_star[length(i_star)])
        
        
        if (identical(intersect(DBN[['nodes']][[from[1]]][[from[2]]][['children']], DBN[['nodes']][[i_1]][[temp_i]][['children']]), character(0))) {
          DBN[['nodes']][[from[1]]][[from[2]]][['mb']] <-
            delete_element(DBN[['nodes']][[from[1]]][[from[2]]][['mb']], i)
          DBN[['nodes']][[i_1]][[temp_i]][['mb']] <-
            delete_element(DBN[['nodes']][[i_1]][[temp_i]][['mb']], from_)
        }
      }
    }
  }
  else{
    stop("ERROR: Arc do not exists!!!")
  }
  DBN
}

# REVERSING ARC FUNCTION

#' Function for arc reversal in Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#' @param from two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param to two-element character vector c(var_name, time), with time is \eqn{t_0}, \eqn{t} or \eqn{t-k} (with \eqn{0 <= k <= markovian order})
#' @param cycle_OK boolean value, TRUE if cycle are admitted (default)
#'
#' @return object of class 'DBN'
#' @export
#'
#' @examples
#' DBN_example <- add_arc_DBN(DBN=DBN_example,from=c('A','t-1'),to=c('R','t'))
reverse_arc_DBN <- function(DBN, from, to, cycle_OK = TRUE) {
  if (!class(DBN) == 'DBN')
    stop("Error: DBN argument is not of class 'DBN'")
  if (!is.character(from))
    stop("Error: from is not a character")
  if (!is.character(to))
    stop("Error: to is not a character")
  from_ <-
    gsub(" ", "", paste(from[1], "_", ifelse(from[2] == 't_0', '0', from[2])))
  to_ <-
    gsub(" ", "", paste(to[1], "_", ifelse(to[2] == 't_0', '0', to[2])))
  if (!(to_ %in% DBN[['nodes']][[from[1]]][[from[2]]][['children']])) {
    stop("ERROR: An arc that does not exist could not be reversed!!!")
  }
  if (any(apply(DBN$learning$blacklist, 1, function(x)
    (x[['from']] == to_ & x[['to']] == from_)))) {
    stop("ERROR: Temporal arcs are irreversible in DBNs!!!")
  }
  DBN <-
    DynamicBayesianNetwork::delete_arc_DBN(DBN = DBN, from = from, to = to)
  DBN <-
    DynamicBayesianNetwork::add_arc_DBN(
      DBN = DBN,
      from = to,
      to = from,
      cycle_OK = TRUE
    )
  DBN
}

# DBN SUMMARY FUNCTION

#' Function returning a summary of input Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#'
#' @details
#' The function returns from terminal:
#' - the model string corresponding to the DAG at G_0
#' - the model string corresponding to the DAG at G_transition
#' - the Markovian Order of the process
#' - the number of nodes, divided between Dynamic nodes and Static nodes
#' - the number of arcs, divided between arcs of G_0, inner-slice arcs and intra-slice arcs of G_transition
#' - the average Markov Blanket size of the nodes
#' - the average neighboorhood size of the nodes
#'
#' @export
#'
#' @examples
#' summary(DBN_example)
summary_DBN <- function(DBN, test = FALSE) {
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  mo = DBN$markov_order
  nodes <- length(DBN$nodes)
  dynamic_nodes = 0
  mbs <- c()
  neighbs <- c()
  for (i in names(DBN$nodes)) {
    dynamic_nodes = dynamic_nodes + ifelse(DBN$nodes[[i]][['type']] == 'Dynamic', 1, 0)
    mbs <-
      c(mbs, length(DBN$nodes[[i]][['t_0']][['mb']]), length(DBN$nodes[[i]][['t']][['mb']]))
    neighbs <-
      c(neighbs, length(DBN$nodes[[i]][['t_0']][['nbr']]), length(DBN$nodes[[i]][['t']][['nbr']]))
  }
  static_nodes <- nodes - dynamic_nodes
  arcs = nrow(DBN$arcs)
  intra_slice_arcs = 0
  inner_slice_arcs = length(unlist(sapply(DBN$arcs[, 'from'], function(x) {
    quanteda::char_select(x, '*_t')
  })))
  prior_arcs = length(unlist(sapply(DBN$arcs[, 'to'], function(x) {
    quanteda::char_select(x, '*_0')
  })))
  for (j in 1:mo) {
    intra_slice_arcs <-
      intra_slice_arcs + length(unlist(sapply(DBN$arcs[, 'from'], function(x) {
        quanteda::char_select(x, gsub(" ", "", paste('*_t-', j)))
      })))
  }
  if (test == FALSE) {
    cat('Dynamic Bayesian Network \n\n\n')
    if (bnlearn::directed(from_DBN_to_G_0(DBN))) {
      cat('Prior Network Model\n ',
          bnlearn::modelstring(from_DBN_to_G_0(DBN)),
          '\n\n')
    }
    else{
      cat('Prior Network Model\n  [partially directed graph]\n\n')
    }
    if (bnlearn::directed(from_DBN_to_G_transition(DBN))) {
      cat(
        'Trantision Network Model\n ',
        bnlearn::modelstring(from_DBN_to_G_transition(DBN)),
        '\n\n'
      )
    }
    else{
      cat('Trantision Network Model\n  [partially directed graph]\n\n')
    }
    cat("Markovian Order:           \t\t", mo, '\n\n')
    cat('Nodes:                     \t\t', nodes, '\n')
    cat('\tDynamic Nodes:           \t', dynamic_nodes, '\n')
    cat('\tStatic Nodes:            \t', static_nodes, '\n\n')
    cat('Arcs:                      \t\t', nrow(DBN$arcs), '\n')
    cat('\tPrior Arcs:              \t', prior_arcs, '\n')
    cat('\tInner-slice Arcs:        \t', inner_slice_arcs, '\n')
    cat('\tIntra-slice Arcs:        \t', intra_slice_arcs, '\n\n')
    cat('Avg Markov Blanket size:   \t\t', round(mean(mbs), digits = 4), '\n')
    cat('Avg Neighboorhood size:    \t\t', round(mean(neighbs), digits = 4), '\n')
  }
  result_list <-
    list(
      markov_order = mo,
      n_nodes = nodes,
      n_dynamic_nodes = dynamic_nodes,
      n_static_nodes = static_nodes,
      n_arcs = nrow(DBN$arcs),
      n_prior_arcs = prior_arcs,
      n_inner_slice_arcs = inner_slice_arcs,
      n_intra_slice_arcs = intra_slice_arcs,
      avg_mb = round(mean(mbs), digits = 4),
      avg_nbr = round(mean(neighbs), digits = 4)
    )
  result_list
}

# DBN TO G_TRANSITION (OF CLASS BN)

#' Transformation of a DBN in a G_transition network
#'
#' @param DBN object of class 'DBN'
#'
#' @return object of class 'bn'
#' @export
#'
#' @examples
#' G_transition <- from_DBN_to_G_transition(DBN_example)
from_DBN_to_G_transition <- function(DBN) {
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  TN = list(
    learning = list(
      whitelist = NULL,
      blacklist = NULL,
      test = 'none',
      ntests = 0,
      algo = 'empty',
      args = list()
    ),
    arcs = matrix(
      ncol = 2,
      nrow = 0,
      byrow = TRUE,
      dimnames = list(character(0), c("from", "to"))
    ),
    nodes = list()
  )
  for (i in 0:DBN$markov_order) {
    for (j in names(DBN$nodes)) {
      if (DBN[['nodes']][[j]][['type']] == 'Dynamic') {
        TN[['nodes']][[gsub(" ", "", paste(j, '_', ifelse(i == 0, 't', gsub(
          " ", "", paste('t-', i)
        ))))]] = DBN[['nodes']][[j]][[ifelse(i == 0, 't', gsub(" ", "", paste('t-', i)))]]
        for (k in DBN[['nodes']][[j]][[ifelse(i == 0, 't', gsub(" ", "", paste('t-', i)))]][['parents']]) {
          TN[['arcs']] <-
            rbind(TN$arcs, c(k, gsub(" ", "", paste(
              j, '_', ifelse(i == 0, 't', gsub(" ", "", paste('t-', i)))
            ))))
        }
      }
    }
  }
  class(TN) <- "bn"
  TN
}

# DBN TO G_0 (OF CLASS BN)

#' Transformation of a DBN in a G_0 network
#'
#' @param DBN object of class 'DBN'
#'
#' @return object of class 'bn'
#' @export
#'
#' @examples
#' G_0 <- from_DBN_to_G_0(DBN_example)
from_DBN_to_G_0 <- function(DBN) {
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  PN = list(
    learning = list(
      whitelist = NULL,
      blacklist = NULL,
      test = 'none',
      ntests = 0,
      algo = 'empty',
      args = list()
    ),
    arcs = matrix(
      ncol = 2,
      nrow = 0,
      byrow = TRUE,
      dimnames = list(character(0), c("from", "to"))
    ),
    nodes = list()
  )
  for (j in names(DBN$nodes)) {
    PN[['nodes']][[gsub(" ", "", paste(j, '_', 0))]] = DBN[['nodes']][[j]][['t_0']]
    for (k in DBN[['nodes']][[j]][['t_0']][['parents']]) {
      PN[['arcs']] <-
        rbind(PN$arcs, c(k, gsub(" ", "", paste(j, '_', 0))))
    }
  }
  class(PN) <- "bn"
  PN
}

# DEFINE CPTs FROM TERMINAL

#' Function for Conditionally Probability Table definition by terminal
#'
#' @param var_name target variable to be conditioned (optional - default: ' ')
#' @param DBN object of class 'DBN' (optional - default: an empty DBN with no nodes)
#' @param nodes_level named list with the listed levels for each node
#'
#' @return multi-dimensional numerical vector
#' @export
#'
#' @examples
#' CPT <- define_CPT() #no information about the variables
#' B_0.prob <- define_CPT(var_name='B_0', DBN=DBN_example) #definition of the CPT for target variable B_0 given the DBN
define_CPT <-
  function(var_name = ' ',
           DBN = empty_DBN(dynamic_nodes = c(), markov_order = 1),
           nodes_level = list()) {
    if (!class(DBN) == 'DBN')
      stop("ERROR: DBN argument is not of class 'DBN'")
    if (!is.character(var_name))
      stop("ERROR: var_name is not a character")
    if (var_name == ' ' & length(DBN$nodes) == 0) {
      var_name <- readline(prompt = "target variable name: ")
      n_vars_conditioning_set <-
        readline(prompt = "number of variables conditioning set: ")
      conditioning_set <- c()
      if (n_vars_conditioning_set < 0) {
        stop("ERROR: number of variables in the conditioning set cannot be negative")
      }
      else if (n_vars_conditioning_set > 0) {
        for (j in 1:max(1, n_vars_conditioning_set)) {
          conditioning <-
            readline(prompt = paste("conditional variable", j, ': '))
          if (conditioning %in% c(var_name, conditioning_set)) {
            stop("ERROR: variable name already used")
          }
          conditioning_set <- c(conditioning_set, conditioning)
        }
      }
    }
    else{
      var_star <-
        var_star <- strsplit(var_name, '_')[[1]]
      var_1 <-
        paste(var_star[1:(length(var_star) - 1)], collapse = '_')
      temp_var <-
        ifelse(var_star[length(var_star)] == '0', 't_0', var_star[length(var_star)])
      if (class(DBN) == 'DBN' &
          var_1 %in% names(DBN$nodes) &
          temp_var %in% names(DBN$nodes[[var_1]])) {
        cat('target variable name:', var_name)
        n_vars_conditioning_set <-
          length(DBN$nodes[[var_1]][[temp_var]][['parents']])
        cat('\nconditioning set:', DBN$nodes[[var_1]][[temp_var]][['parents']])
        conditioning_set <- c()
        for (v in DBN$nodes[[var_1]][[temp_var]][['parents']]) {
          conditioning_set <- c(conditioning_set, v)
        }
      }
      else{
        stop("ERROR: var_name must be a DBN node")
      }
    }
    if (var_name %in% names(nodes_level)) {
      var_n_levels <- length(nodes_level[[var_name]])
      variable_set <- list()
      variable_set[[var_name]] <- nodes_level[[var_name]]
      cat('\n', var_name, "levels:", nodes_level[[var_name]], '\n')
    }
    else{
      var_n_levels <-
        as.numeric(readline(prompt = paste(var_name, "number of levels: ")))
      if (var_n_levels <= 1) {
        stop("ERROR: target variable must have at least 2 levels")
      }
      var_levels <- c()
      for (i in 1:var_n_levels) {
        var_levels <-
          c(var_levels, readline(prompt = paste("level", i, 'of', var_name, ': ')))
      }
      variable_set <- list()
      variable_set[[var_name]] <- var_levels
    }
    conditioning_set_levels <- c(var_n_levels)
    if (n_vars_conditioning_set > 0) {
      for (j in 1:n_vars_conditioning_set) {
        conditioning <- conditioning_set[j]
        if (conditioning %in% names(nodes_level)) {
          conditioning_n_levels <- length(nodes_level[[conditioning]])
          conditioning_set_levels <-
            c(conditioning_set_levels, conditioning_n_levels)
          variable_set[[conditioning]] <-
            nodes_level[[conditioning]]
          cat(conditioning, "levels:", nodes_level[[conditioning]], '\n')
        }
        else{
          conditioning_n_levels <-
            as.numeric(readline(prompt = paste(conditioning, "number of levels: ")))
          if (conditioning_n_levels <= 1) {
            stop("ERROR: variables in the conditioning set must have at least 2 levels")
          }
          conditioning_set_levels <-
            c(conditioning_set_levels, conditioning_n_levels)
          conditioning_levels <- c()
          for (k in 1:conditioning_n_levels) {
            conditioning_levels <-
              c(conditioning_levels, readline(prompt = paste(
                "level", k, 'of', conditioning, ': '
              )))
          }
          variable_set[[conditioning]] <- conditioning_levels
        }
      }
    }
    probabilities <- c()
    for (row in 1:nrow(expand.grid(variable_set))) {
      if (numbers::rem(row, var_n_levels) == 1) {
        prob_sum = 0
      }
      for (var in 1:length(variable_set)) {
        if (var == 1) {
          string <-
            ifelse(
              n_vars_conditioning_set > 0,
              paste("P(", var_name, "=", expand.grid(variable_set)[row, 1], '|'),
              paste("P(", var_name, "=", expand.grid(variable_set)[row, 1], ') =')
            )
        }
        else if (var == length(variable_set)) {
          string <-
            paste(
              string,
              conditioning_set[var - 1],
              '=',
              expand.grid(variable_set)[row, var],
              ifelse(((prob_sum == 1) | (numbers::rem(row, var_n_levels) == 0)
              ), ')', ') =')
            )
        }
        else{
          string <-
            paste(string,
                  conditioning_set[var - 1],
                  '=',
                  expand.grid(variable_set)[row, var],
                  ',')
        }
      }
      if (prob_sum == 1) {
        cat(paste(string, 'automatically set to 0\n'))
        probabilities = c(probabilities, 0)
        next
      }
      else if (numbers::rem(row, var_n_levels) == 0) {
        cat(paste(
          string,
          'automatically set to',
          as.character(1 - prob_sum),
          '\n'
        ))
        probabilities = c(probabilities, 1 - prob_sum)
        next
      }
      prob = as.numeric(readline(prompt = string))
      if (prob < 0 | prob > 1) {
        stop("ERROR: Probabilities must vary in [0,1]")
      }
      prob_sum = prob_sum + prob
      if (prob_sum > 1) {
        stop("ERROR: The sum of probabilities given the conditioning set exceed 1")
      }
      probabilities <- c(probabilities, prob)
    }
    cat('\n')
    assign(
      gsub(' ', '', paste(var_name, '.prob')),
      array(probabilities, dim = conditioning_set_levels, dimnames = variable_set)
    )
    get(gsub(' ', '', paste(var_name, '.prob')))
  }

#' Function for Conditionally Probability Tables definition by terminal given the Dynamic Bayesian Network
#'
#' @param DBN object of class 'DBN'
#'
#' @return list of multi-dimensional vector (CPTs for each DBN node)
#' @export
#'
#' @examples
#' define_CPTs(DBN_example)
define_CPTs <-
  function(DBN = empty_DBN(dynamic_nodes = c(), markov_order = 1)) {
    static_nodes <- bnlearn::node.ordering(from_DBN_to_G_0(DBN))
    dynamic_nodes <-
      quanteda::char_select(bnlearn::node.ordering(from_DBN_to_G_transition(DBN)),
                            "*t",
                            valuetype = "glob")
    CPTs <- list()
    defined_levels <- list()
    for (i in static_nodes) {
      CPT <-
        define_CPT(var_name = i,
                   DBN = DBN,
                   nodes_level = defined_levels)
      CPTs[[i]] <- CPT
      def_lev <-
        c(dimnames(CPT),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(
                x
              ) - 2)), '_t')) # generalize for var_name length > 1
            })
          ))),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(
                x
              ) - 2)), '_t-1')) # generalize for var_name length > 1
            })
          ))))
      defined_levels <-
        c(defined_levels, def_lev[setdiff(names(def_lev), names(defined_levels))])
    }
    for (j in dynamic_nodes) {
      CPT <-
        define_CPT(var_name = j,
                   DBN = DBN,
                   nodes_level = defined_levels)
      CPTs[[j]] <- CPT
    }
    CPTs
  }