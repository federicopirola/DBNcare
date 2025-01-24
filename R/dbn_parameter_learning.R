# DEFINE FUNCTION FOR dbn.fit OBJECT CREATION -> NODES (with childrens, parents and CPTs)

#' Function for parameter learning in Dynamic Bayesian Networks
#'
#' @param DBN object of class 'DBN'
#' @param CPTs list of multi-dimensional vector (CPTs for each DBN node)
#' @param data data.frame object
#' @param replace.unidentifiable If TRUE conditional probabilities for unobserved parents combinations (unidentifiable parameters) are replaced by uniform conditional probabilities, if FALSE (default) they are set as NA
#'
#' @return object of class 'DBN.fit'
#' @export
#'
#' @examples
#' learned_dbn <- DBN_parameters(DBN = DBN_example, data = sampling_set)
DBN_parameters <- function(DBN,
                           CPTs = list(),
                           data = data.frame(),
                           replace.unidentifiable = FALSE) {
  suppressPackageStartupMessages(library(dplyr))
  if (!class(DBN) == 'DBN')
    stop("ERROR: DBN argument is not of class 'DBN'")
  if (!class(CPTs) == 'list')
    stop("ERROR: CPTs argument is not of class 'list'")
  if (!class(data) == 'data.frame')
    stop("ERROR: data argument is not of class 'data.frame'")
  if (any(is.na(DBN)))
    stop("ERROR: missing data detected")
  static_nodes <- bnlearn::node.ordering(from_DBN_to_G_0(DBN))
  dynamic_nodes <-
    quanteda::char_select(bnlearn::node.ordering(from_DBN_to_G_transition(DBN)),
                          "*t",
                          valuetype = "glob")
  nodes <- c(static_nodes, dynamic_nodes)
  vars <-
    colnames(data)[!colnames(data) %in% c('Time', 'Sample_id')]
  lvs <- list()
  dbn_fitted = list()
  for (i in vars) {
    var <- gsub(" ", "", paste(i, '_0'))
    var_1 <- gsub(" ", "", paste(i, '_t'))
    var_2 <- gsub(" ", "", paste(i, '_t-1'))
    lvs[[var]] <- sort(as.array(levels(factor(data[[i]]))))
    lvs[[var_1]] <- sort(as.array(levels(factor(data[[i]]))))
    lvs[[var_2]] <- sort(as.array(levels(factor(data[[i]]))))
  }
  if (length(CPTs) == 0 & nrow(data) > 0) {
    if (!setequal(names(DBN$nodes), setdiff(colnames(data), c('Sample_id', 'Time'))))
      stop("ERROR: nodes in DBN and variables in dataframe do not match")
    #check if each combination have cardinality one
    if (!(all(dplyr::count(
      data %>% dplyr::group_by(Sample_id, Time)
    )$n == 1)))
      stop("ERROR: One or mode combinations of IDs and time slices is repeated")
    #check presence of all combinations
    if (any(is.na(data %>% tidyr::complete(Sample_id, Time))))
      stop("ERROR: One or more sample/individual present an incomplete temporal sequences")
    df_0 <- data[data$Time == 0,]
    names(df_0) <-
      lapply(names(data), function(x)
        ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_0'))))
    df_transition <- data
    names(df_transition) <-
      lapply(names(data), function(x)
        ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_t-1'))))
    for (k in dynamic_nodes) {
      df_transition <-
        df_transition %>% dplyr::group_by(Sample_id) %>% dplyr::mutate(!!k := lead(get(gsub(
          " ", "", paste(k, '-1')
        )), n = 1, default = NA))
    }
    df_transition <-
      df_transition %>% dplyr::ungroup() %>% stats::na.omit()
    for (f in static_nodes) {
      f_star <- strsplit(f, '_')[[1]]
      f_1 <- paste(f_star[1:(length(f_star) - 1)], collapse = '_')
      temp_f <-
        ifelse(f_star[length(f_star)] == '0', 't_0', f_star[length(f_star)])
      parents <- DBN[['nodes']][[f_1]][[temp_f]][['parents']]
      children <- DBN[['nodes']][[f_1]][[temp_f]][['children']]
      if (length(parents) > 0) {
        pr <-
          df_0[, c(rev(parents), f)] %>% dplyr::group_by_all() %>% dplyr::count() %>% dplyr::ungroup() %>% tidyr::complete(!!! rlang::syms(c(rev(parents), f))) %>% replace(is.na(.), 0) %>% group_by(across(rev(parents))) %>% reframe(!!sym(f), prob = n/sum(n)) %>% dplyr::ungroup() %>% dplyr::arrange_all(.vars = c(rev(parents), f))
        #prob_vec <-
        #  unlist(lapply(c(1:length(pr$n)), function(x)
        #    ifelse(
        #      numbers::rem(x, 2) == 1,
        #      pr$n[x] / (pr$n[x] + pr$n[x + 1]),
        #      pr$n[x] / (pr$n[x] + pr$n[x - 1])
        #    )))
        if (any(is.na(pr$prob))){
          if (replace.unidentifiable){
            pr$prob <- replace(pr$prob,is.na(pr$prob),1/length(lvs[[f]]))
          }
          else{
            warning(
              "WARNING: Probabilities of the conditioning set equal to 0: Relative frequency is NULL")  
          }
        }
        dbn_fitted[[f]] <-
          list(
            node = f,
            parents = parents,
            children = children,
            prob = array(
              pr$prob,
              dim = unname(unlist(lapply(
                lvs[c(f, parents)], length
              ))),
              dimnames = lvs[c(f, parents)]
            )
          )
      }
      else{
        prob_vec <-
          unlist(lapply(lvs[[f]], function(x)
            nrow(df_0[df_0[, f] == x,]) / nrow(df_0)))
        dbn_fitted[[f]] <-
          list(
            node = f,
            parents = parents,
            children = children,
            prob = array(
              prob_vec,
              dim = unname(unlist(lapply(
                lvs[c(f, parents)], length
              ))),
              dimnames = lvs[c(f, parents)]
            )
          )
      }
    }
    for (g in dynamic_nodes) {
      g_star <- strsplit(g, '_')[[1]]
      g_1 <- paste(g_star[1:(length(g_star) - 1)], collapse = '_')
      temp_g <-
        ifelse(g_star[length(g_star)] == '0', 't_0', g_star[length(g_star)])
      parents <- DBN[['nodes']][[g_1]][[temp_g]][['parents']]
      children <- DBN[['nodes']][[g_1]][[temp_g]][['children']]
      if (length(parents) > 0) {
        pr <-
          df_transition[, c(rev(parents), g)] %>% dplyr::group_by_all() %>% dplyr::count() %>% dplyr::ungroup() %>% tidyr::complete(!!! rlang::syms(c(rev(parents), g))) %>% replace(is.na(.), 0) %>% group_by(across(rev(parents))) %>% reframe(!!sym(g),prob = n/sum(n)) %>% dplyr::ungroup() %>%  dplyr::arrange_all(.vars = c(rev(parents), g))
        #prob_vec <-
        #  unlist(lapply(c(1:length(pr$n)), function(x)
        #    ifelse(
        #      numbers::rem(x, 2) == 1,
        #      pr$n[x] / (pr$n[x] + pr$n[x + 1]),
        #      pr$n[x] / (pr$n[x] + pr$n[x - 1])
        #    )))
        if (any(is.na(pr$prob))){
          if (replace.unidentifiable){
            pr$prob <- replace(pr$prob,is.na(pr$prob),1/length(lvs[[g]]))
          }
          else{
            warning(
              "WARNING: Probabilities of the conditioning set equal to 0: Relative frequency is NULL")  
          }
        }
        dbn_fitted[[g]] <-
          list(
            node = g,
            parents = parents,
            children = children,
            prob = array(
              pr$prob,
              dim = unname(unlist(lapply(
                lvs[c(g, parents)], length
              ))),
              dimnames = lvs[c(g, parents)]
            )
          )
      }
      else{
        prob_vec <-
          unlist(lapply(lvs[[g]], function(x)
            nrow(df_transition[df_transition[, g] == x,]) / nrow(df_transition)))
        dbn_fitted[[g]] <-
          list(
            node = g,
            parents = parents,
            children = children,
            prob = array(
              prob_vec,
              dim = unname(unlist(lapply(
                lvs[c(g, parents)], length
              ))),
              dimnames = lvs[c(g, parents)]
            )
          )
        
      }
    }
    class(dbn_fitted) <- "dbn.fit"
    dbn_fitted
  }
  else if (length(CPTs) > 0 & nrow(data) == 0) {
    if (!setequal(names(CPTs), nodes)) {
      stop("ERROR: nodes in DBN and variables in CPTs do not match")
    }
    defined_levels <- list()
    nodes_cpt <- names(DBN$nodes)
    nodes_info <- list()
    if (!(any(lapply(CPTs, class)  %in% c('matrix', 'array')))) {
      stop("ERROR: CPT must be of class 'matrix' or 'array'")
    }
    for (i in nodes) {
      CPT <- CPTs[[i]]
      if (!(all(lapply(CPT, class) == 'numeric'))) {
        stop("ERROR: Probabilities must be numeric")
      }
      if (length(dim(CPT)) > 1) {
        l <- length(dim(CPT))
        idx_target <- which(names(dimnames(CPT)) == i)
        if (!(all(apply(
          CPT, setdiff(1:l, idx_target), sum
        ) == as.character(1)))) {
          stop("ERROR: Probabilities for each conditioning set must sum to 1")
        }
      }
      else{
        if (sum(CPT) != as.character(1)) {
          stop("ERROR: Probabilities for each conditioning set must sum to 1")
        }
      }
      i_star <- strsplit(i, '_')[[1]]
      i_1 <- paste(i_star[1:(length(i_star) - 1)], collapse = '_')
      temp_i <-
        ifelse(i_star[length(i_star)] == '0', 't_0', i_star[length(i_star)])
      parents <- DBN[['nodes']][[i_1]][[temp_i]][['parents']]
      children <- DBN[['nodes']][[i_1]][[temp_i]][['children']]
      if (!setequal(setdiff(names(dimnames(CPT)), i), parents)) {
        stop("ERROR: CPTs do not match parents set")
      }
      def_lev <-
        c(dimnames(CPT),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(x)-2)), '_t'))
            })
          ))),
          setNames(dimnames(CPT), array(unlist(
            sapply(names(dimnames(CPT)) , function(x) {
              gsub(" ", "", paste(substring(x, 1, (nchar(x)-2)), '_t-1'))
            })
          ))))
      int_nodes <- intersect(names(def_lev), names(defined_levels))
      for (c in int_nodes) {
        if (!(setequal(def_lev[[c]], defined_levels[[c]]))) {
          stop("ERROR: Inconsistency in node's levels")
        }
      }
      defined_levels <-
        c(defined_levels, def_lev[setdiff(names(def_lev), names(defined_levels))])
      nodes_info[[i]] <-
        list(
          node = i,
          parents = parents,#setdiff(names(dimnames(CPT)), i),
          children = children,
          prob = aperm(CPT, c(i, parents))
        )
    }
    class(nodes_info) <- "dbn.fit"
    nodes_info
  }
  else{
    stop(
      "ERROR: Only one between data or CPTs must be defined to learn the parameters of the DBN!!"
    )
  }
}

# DEFINE FUNCTION FROM dbn.fit TO G_0 bn.fit

#' Function for G_0 parameters set extraction
#'
#' @param DBN_fitted object of class 'dbn.fit'
#'
#' @return object of class 'bn.fit'
#' @export
#'
#' @examples
#' fitted_0 <- from_fitted_DBN_to_fitted_G_0(fitted_DBN)
from_fitted_DBN_to_fitted_G_0 <- function(DBN_fitted) {
  if (!class(DBN_fitted) == 'dbn.fit')
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  BN_0_fitted <-
    DBN_fitted[quanteda::char_select(names(DBN_fitted), "*0", valuetype = "glob")]
  class(BN_0_fitted) <- "bn.fit"
  BN_0_fitted
}

# DEFINE FUNCTION FROM dbn.fit TO G_transition bn.fit

#' Function for G_transition parameters set extraction
#'
#' @param DBN_fitted object of class 'dbn.fit'
#'
#' @return object of class 'bn.fit'
#' @export
#'
#' @examples
#' fitted_transition <- from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
from_fitted_DBN_to_fitted_G_transition <- function(DBN_fitted) {
  if (!class(DBN_fitted) == 'dbn.fit')
    stop("ERROR: DBN_fitted argument is not of class 'dbn.fit'")
  BN_transition_fitted <-
    DBN_fitted[quanteda::char_select(names(DBN_fitted), "*t", valuetype = "glob")]
  class(BN_transition_fitted) <- "bn.fit"
  BN_transition_fitted
}
