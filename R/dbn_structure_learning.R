#' A function returning the blacklist given a DBN object
#'
#' @param DBN object of class 'DBN'
#'
#' @return data.frame of the arcs in the blacklist
#' @export
#'
blacklist_from_DBN <- function(DBN) {
  bl <- as.data.frame(DBN$learning$blacklist)
  bl
}

#' A function returning the whitelist given a DBN object
#'
#' @param DBN object of class 'DBN'
#'
#' @return data.frame of the arcs in the whitelist
#' @export
#'
whitelist_from_DBN <- function(DBN) {
  wl <- as.data.frame(DBN$learning$whitelist)
  wl
}

#' It returns a blacklist matrix represented by the pair of temporal nodes which are forbidden  
#'
#' @param static_nodes character vector of time independent variables
#' @param dynamic_nodes character vector of time dependent variables
#' @param markov_order integer
#'
#' @return matrix of node pairs which edge is forbidden
#' @export
#'
#' @examples
#' blacklist.DBN(dynamic_nodes = ('A','B'))
blacklist.DBN <-
  function(static_nodes = c(),
           dynamic_nodes,
           markov_order = 1) {
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

#' Given a data.frame it returns a blacklist matrix represented by the pair of temporal nodes which are forbidden
#'
#' @param data a data.frame with column 'Time', 'Sample_id' and the other columns with the values of the time series of 'Sample_id' at time 'Time' 
#'
#' @return matrix of node pairs which edge is forbidden
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' forbidden_edges <- blacklist.temp.DBN(data)
blacklist.temp.DBN <- function(data){
  bl <- blacklist.DBN(dynamic_nodes=setdiff(names(data), c("Time", "Sample_id")))
  bl[which(arr.ind = T, grepl("^.+(t|t-1)$",bl[,'from']) & grepl("^.+(t-1)$",bl[,'to'])),]
}

#' Preprocessing function for Structure Learning algorithm. It scans for the consistency of structure learning parameters in Dynamic Bayesian Networks
#'
#' @param data object of type data.frame to be given as an input to structure learning algorithm
#' @param blacklist matrix with 'from' and 'to' columns of forbidden edges in the Dynamic Bayesian Network (default is NULL)
#' @param whitelist matrix with 'from' and 'to' columns of fixed edges in the Dynamic Bayesian Network (default is NULL)
#' @param method a string defining the structure learning method class, possible values are 'constraint' (default), 'score' and 'hybrid'
#' @param test a string defining the test for constraint based or hybrid algorithm (default is `Mutual Information` mi)
#' @param score a string defining the score for score based or hybrid algorithm (default is `Bayesian Information Criterion` bic )
#' @param max.sx maximum number of parent per node (default is no limit)
#'
#' @return character vector made of ordered data to fit G_0, data to fit G_transition, blacklist to fit G_0, blacklist to fit G_transition, whitelist to fit G_0, whitelist to fit G_transition
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' StructureLearning.Preprocess(data)
StructureLearning.Preprocess <- function(data, blacklist = NULL, whitelist = NULL, method = 'constraint', test = 'mi', score = 'bic', max.sx = NULL){
  is.naturalnumber <-
    function(x, tol = .Machine$double.eps^0.5)  x > tol & abs(x - round(x)) < tol
  
  if (!is.null(max.sx)){
    if(!is.naturalnumber(max.sx)){
      stop('ERROR: max.sx must be a positive integer')
    }
  }
  if (method == 'constraint'){
    if (!is.null(test)){
      if (!test %in% c("mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh",  "x2-adf", "mc-x2", "smc-x2", "sp-x2")){
        stop('ERROR: Test must be one the following: "mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh",  "x2-adf", "mc-x2", "smc-x2", "sp-x2"')
      }
    }
  }
  else if (method == 'score'){
    if (!is.null(score)){
      if (!score %in% c("loglik", "aic", "bic", "ebic", "bde",  "bds", "mbde", "bdla", "k2", "fnml", "qnml", "nal", "pnal")){
        stop('ERROR: Score must be one the following: "loglik", "aic", "bic", "ebic", "pred-loglik", "bde",  "bds", "mbde", "bdla", "k2", "fnml", "qnml", "nal", "pnal"')
      }
    }
  }
  else if (method == 'hybrid'){
    if (!is.null(test)){
      if (!test %in% c("mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh",  "x2-adf", "mc-x2", "smc-x2", "sp-x2")){
        stop('ERROR: Test must be one the following: "mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh",  "x2-adf", "mc-x2", "smc-x2", "sp-x2"')
      }
    }
    if (!is.null(score)){
      if (!score %in% c("loglik", "aic", "bic", "ebic", "pred-loglik", "bde",  "bds", "mbde", "bdla", "k2", "fnml", "qnml", "nal", "pnal")){
        stop('ERROR: Score must be one the following: "loglik", "aic", "bic", "ebic", "pred-loglik", "bde",  "bds", "mbde", "bdla", "k2", "fnml", "qnml", "nal", "pnal"')
      }
    }
  }
  else{
    stop("ERROR: method must be 'score', 'constraint' or 'hybrid'")
  }
  if (!all(c('Time','Sample_id') %in% names(data))){
    stop("ERROR: Sample_id and Time are not in data columns")
  }
  if (!is.character(data$Sample_id)){
    stop("ERROR: Sample_id must be character")
  }
  if (!is.numeric(data$Time) & !is.character(data$Time)){
    stop("ERROR: Time must be numeric or character")
  }
  if (any(is.na(data)))
    stop("ERROR: Missing data detected")
  if (!(all(dplyr::count(
    data %>% dplyr::group_by(Sample_id, Time)
  )$n == 1)))
    stop("ERROR: One or mode combinations of IDs and time slices is repeated")
  #check presence of all combinations
  if (any(is.na(data %>% tidyr::complete(Sample_id, Time))))
    stop("ERROR: Missing data detected: One or more samples present an incomplete temporal sequences")
  if (!"matrix" %in% class(blacklist) & !is.null(blacklist)){
    stop("ERROR: blacklist must be a matrix")
  }
  else if(!is.null(blacklist)){
    if (!all(colnames(blacklist) == c("from","to"))){
      stop("ERROR: blacklist must be a matrix with columns 'from' and 'to'")
    }
    if (any(duplicated(as.data.frame(blacklist)))){
      stop("ERROR: blacklist could not contain duplicates")
    }
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t|t-1|0)$',sep = ''),blacklist))){
      stop("ERROR: blacklist entries must be in the form `var`_0 `var`_t or `var`_t-1, where var is a column of data")
    }
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_0$',sep = ''),blacklist[,'to']) == grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_0$',sep = ''),blacklist[,'from']))){
      stop("ERROR: blacklist ('from', 'to') must be pairs in the form (`var1`_0, `var2`_0), (`var1`_t-1, `var2`_t) or (`var1`_t, `var2`_t), where var1 and var2 are columns of data")
    } 
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t)$',sep = ''),blacklist[,'to']) == grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t|t-1)$',sep = ''),blacklist[,'from']))){
      stop("ERROR: blacklist ('from', 'to') must be pairs in the form (`var1`_0, `var2`_0), (`var1`_t-1, `var2`_t) or (`var1`_t, `var2`_t), where var1 and var2 are columns of data")
    }
  }
  if (!"matrix" %in% class(whitelist) & !is.null(whitelist)){
    stop("ERROR: whitelist must be a matrix")
  }
  else if(!is.null(whitelist)){
    if (!all(colnames(whitelist) == c("from","to"))){
      stop("ERROR: whitelist must be a matrix with columns 'from' and 'to'")
    }
    if (any(duplicated(as.data.frame(whitelist)))){
      stop("ERROR: whitelist could not contain duplicates")
    }
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t|t-1|0)$',sep = ''),whitelist))){
      stop("ERROR: whitelist entries must be in the form `var`_0 `var`_t or `var`_t-1, where var is a column of data")
    }
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_0$',sep = ''),whitelist[,'to']) == grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_0$',sep = ''),whitelist[,'from']))){
      stop("ERROR: whitelist ('from', 'to') must be pairs in the form (`var1`_0, `var2`_0), (`var1`_t-1, `var2`_t) or (`var1`_t, `var2`_t), where var1 and var2 are columns of data")
    } 
    if (!all(grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t)$',sep = ''),whitelist[,'to']) == grepl(paste('^(',paste(setdiff(names(data), c("Time", "Sample_id")), collapse = "|"),')_(t|t-1)$',sep = ''),whitelist[,'from']))){
      stop("ERROR: whitelist ('from', 'to') must be pairs in the form (`var1`_0, `var2`_0), (`var1`_t-1, `var2`_t) or (`var1`_t, `var2`_t), where var1 and var2 are columns of data")
    }
  }
  if (any(duplicated(as.data.frame(rbind(blacklist,whitelist))))){
    stop("ERROR: an element could not be both in the whitelist and in the blacklist")
  }
  
  blacklist <- unique(rbind(blacklist, blacklist.temp.DBN(data)))
  
  data[setdiff(names(data),c('Sample_id','Time'))] <- lapply(data[setdiff(names(data),c('Sample_id','Time'))],as.factor)
  
  data_0 <- data[data$Time == 0,]
  excluded_vars <- c("Time", "Sample_id")
  data_0 <- select(data_0, -one_of(excluded_vars))
  if (!all(sapply(data_0, nlevels) >=2)){
    stop("ERROR: variables at time 0 must have at least 2 levels")
  }
  names(data_0) <-
    lapply(names(data_0), function(x)
      ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_0'))))
  
  data_transition <- data
  names(data_transition) <-
    lapply(names(data), function(x)
      ifelse(x %in% c('Sample_id', 'Time'), x, gsub(" ", "", paste(x, '_t-1'))))
  for (k in unlist(lapply(names(data)[!(names(data) %in% excluded_vars)], function(x) gsub(" ", "", paste(x, '_t'))))) {
    data_transition <-
      data_transition %>% dplyr::group_by(Sample_id) %>% dplyr::mutate(!!k := lead(get(gsub(
        " ", "", paste(k, '-1')
      )), n = 1, default = NA))
  }
  data_transition <-
    data_transition %>% dplyr::ungroup() %>% stats::na.omit(na.action=na.pass)
  data_transition <- select(data_transition, -one_of(excluded_vars))
  data_transition <- as.data.frame(data_transition)
  if (!all(sapply(data_transition, nlevels) >=2)){
    stop("ERROR: variables at time t must have at least 2 levels")
  }
  blacklist_0 <-  blacklist[which(arr.ind = T, grepl("^.+(0)$",blacklist[,'from'])),] 
  blacklist_transition <- blacklist[which(arr.ind = T, grepl("^.+(t|t-1)$",blacklist[,'from'])),] 
  whitelist_0 <- whitelist[which(arr.ind = T, grepl("^.+(0)$",whitelist[,'from'])),] 
  whitelist_transition <- whitelist[which(arr.ind = T, grepl("^.+(t|t-1)$",whitelist[,'from'])),] 
  list(data_0,data_transition,blacklist_0, blacklist_transition, whitelist_0, whitelist_transition)
}

#' Generation of a Dynamic Bayesian Network given two compatible networks \eqn{G_0} and \eqn{G_{transition}}
#'
#' @param PN \eqn{G_0}, an object of class 'bn'
#' @param TN \eqn{G_{transition}}, an object of class 'bn'
#' @param nodes character vector of node names
#'
#' @return An object of class 'DBN' with the characteristics of \eqn{G_0} and \eqn{G_{transition}}
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' c(data0, dataTR, ., blacklistTR, ., .) %<-% StructureLearning.Preprocess(data = data)
#' bn_sampled_0 <- pc.stable(data0)
#' bn_sampled_TR <- pc.stable(dataTR, blacklist = blacklistTR)
#' StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
#' 
StructureLearning.Data <- function(PN, TN, nodes) {
  if (!class(PN) == 'bn')
    stop("ERROR: PN argument is not of class 'bn'")
  if (!class(TN) == 'bn')
    stop("ERROR: TN argument is not of class 'bn'")
  DBN = list(
    learning = list(
      whitelist = rbind(PN[['learning']][['whitelist']], TN[['learning']][['whitelist']]),
      blacklist = unique(rbind(blacklist.DBN(dynamic_nodes = nodes), PN[['learning']][['blacklist']], TN[['learning']][['blacklist']])),
      test = list(G_0 = PN[['learning']][['test']], G_transition = TN[['learning']][['test']]),
      ntests = list(G_0 = PN[['learning']][['ntests']], G_transition = TN[['learning']][['ntests']]),
      algo = list(G_0 = PN[['learning']][['algo']], G_transition = TN[['learning']][['algo']]),
      args = list(G_0 = PN[['learning']][['args']], G_transition = TN[['learning']][['args']])
    ),
    markov_order = 1,
    arcs = matrix(
      ncol = 2,
      nrow = 0,
      byrow = TRUE,
      dimnames = list(character(0), c("from", "to"))
    ),
    nodes = list()
  )
  vars <- unlist(lapply(names(PN$nodes), function(x) strsplit(x, '_')[[1]][[1]]))
  for (j in vars) {
    DBN[['nodes']][[j]][['t_0']] = PN[['nodes']][[gsub(" ", "", paste(j, '_0'))]]
    DBN[['nodes']][[j]][['type']] = 'Dynamic'
    DBN[['nodes']][[j]][['t']] = TN[['nodes']][[gsub(" ", "", paste(j, '_t'))]]
    DBN[['nodes']][[j]][['t-1']] = TN[['nodes']][[gsub(" ", "", paste(j, '_t-1'))]]
  }
  DBN[['arcs']] <- rbind(PN[['arcs']],TN[['arcs']]) 
  class(DBN) <- "DBN"
  DBN
}


#' Learn the equivalence class of DBN via PC-stable algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via PC stable algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' pc.stable.DBN(data)
#' 
pc.stable.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- pc.stable(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- pc.stable(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Grow-Shrink algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Grow-Shrink algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' gs.DBN(data)
#' 
gs.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- gs(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- gs(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Incremental Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Incremental Association algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' iamb.DBN(data)
#' 
iamb.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Fast Incremental Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Fast Incremental Association algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' fast.iamb.DBN(data)
#' 
fast.iamb.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- fast.iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- fast.iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Interleaved Association algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Interleaved Association algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' inter.iamb.DBN(data)
#' 
inter.iamb.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- inter.iamb(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- inter.iamb(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Incremental Association algorithm with FDR
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Incremental Association algorithm with FDR
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' iamb.fdr.DBN(data)
#' 
iamb.fdr.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- iamb.fdr(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- iamb.fdr(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Max-Min Parents and Children algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Max-Min Parents and Children algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' mmpc.DBN(data)
#' 
mmpc.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- mmpc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- mmpc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Semi-Interleaved HITON-PC algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Semi-Interleaved HITON-PC algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' si.hiton.pc.DBN(data)
#' 
si.hiton.pc.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- si.hiton.pc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- si.hiton.pc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Hybrid Parents and Children algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters for the chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Hybrid Parents and Children algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' hpc.DBN(data)
#' 
hpc.DBN <- function(data, test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'constraint', test = test, max.sx = max.sx)
  bn_sampled_0 <- hpc(data0,  test = test, max.sx = min(max.sx,length(names(data0))-1), blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- hpc(dataTR,  test = test, max.sx = min(max.sx,length(names(dataTR))-1), blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Hill Climbing algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function
#'
#' @return An object of class 'DBN', which is a DBN learned via Hill Climbing algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' hc.DBN(data)
#' 
hc.DBN <- function(data, score = 'bic', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'score', score = score, max.sx = max.sx)
  bn_sampled_0 <- hc(data0, score = score, maxp = min(max.sx,length(names(data0))-1) , blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- hc(dataTR, score = score, maxp = min(max.sx,length(names(dataTR))-1) , blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Tabu Search algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function
#'
#' @return An object of class 'DBN', which is a DBN learned via Tabu Search algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' tabu.DBN(data)
#' 
tabu.DBN <- function(data, score = 'bic', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'score', score = score, max.sx = max.sx)
  bn_sampled_0 <- tabu(data0, score = score, maxp = min(max.sx,length(names(data0))-1) , blacklist = blacklist0, whitelist = whitelist0,...)
  bn_sampled_TR <- tabu(dataTR, score = score, maxp = min(max.sx,length(names(dataTR))-1) , blacklist = blacklistTR, whitelist = whitelistTR,...)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Hybrid HPC algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via Hybrid HPC algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' h2pc.DBN(data)
#'
h2pc.DBN <- function(data, score = 'bic', test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'hybrid', test = test, score = score, max.sx = max.sx)
  bn_sampled_0 <- h2pc(data0, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- h2pc(dataTR, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via Max-Min Hill Climbing algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via  Max-Min Hill Climbing algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' mmhc.DBN(data)
#'
mmhc.DBN <- function(data, score = 'bic', test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'hybrid', test = test, score = score, max.sx = max.sx)
  bn_sampled_0 <- mmhc(data0, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- mmhc(dataTR, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Learn the equivalence class of DBN via 2-phase Restricted Maximization algorithm
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via 2-phase Restricted Maximization algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' rsmax2.DBN(data)
#'
rsmax2.DBN <- function(data, restrict = 'pc.stable', maximize = 'hc', score = 'bic', test = 'mi', max.sx = NULL, blacklist = NULL, whitelist = NULL,...){
  c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% StructureLearning.Preprocess(data = data, blacklist = blacklist, whitelist = whitelist, method = 'hybrid', test = test, score = score, max.sx = max.sx)
  bn_sampled_0 <- rsmax2(data0, restrict = restrict, maximize = maximize, maximize.args = list(score = score, maxp = min(max.sx,length(names(data0))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(data0))-1), ...), blacklist = blacklist0, whitelist = whitelist0)
  bn_sampled_TR <- rsmax2(dataTR, restrict = restrict, maximize = maximize, maximize.args = list(score = score, maxp = min(max.sx,length(names(dataTR))-1)), restrict.args = list(test = test, max.sx = min(max.sx,length(names(dataTR))-1), ...), blacklist = blacklistTR, whitelist = whitelistTR)
  StructureLearning.Data(PN = bn_sampled_0, TN =  bn_sampled_TR, nodes = setdiff(names(data), c('Sample_id','Time')))
}

#' Structure Learning function for Dynamic Bayesian Networks
#'
#' @param data a data.frame to be given as an input the to structure learning algorithm
#' @param algorithm a character string, the label of the structure learning algorithm to be used (default is `hc`)
#' @param algorithm.res a character string, the label of the structure learning algorithm to be used for the restriction phase of the hybrid algorithms (default is `pc.stable`)
#' @param algorithm.max a character string, the label of the structure learning algorithm to be used for the maximization phase of the hybrid algorithms (default is `hc`)
#' @param blacklist a matrix with columns 'from' and 'to' defining the set of arcs that are forbidden
#' @param whitelist a matrix with columns 'from' and 'to' defining the set of arcs that are fixed
#' @param test a character string, the label of the conditional independence test to be used in the algorithm test (default is `Mutual Information` mi)
#' @param score a character string, the label of the scoring function to be used in the algorithm test (default is `Bayesian Information Criterion` bic)
#' @param max.sx an integer, the maximum number of parents for each node (default is no limit) 
#' @param ... additional parameters of the scoring function and/or chosen test 
#'
#' @return An object of class 'DBN', which is a DBN learned via 2-phase Restricted Maximization algorithm
#' @export
#'
#' @examples
#' DBN <- generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5)
#' DBN_fitted <- generate_dbn_nodes_distributions(DBN,c("A","B","C"), TRUE, 2)
#' data <- dbn_sampling(DBN_fitted, 2000, 5)
#' StructureLearning.Network(data, algorithm='h2pc',score='aic') 
StructureLearning.Network <- function(data, algorithm = 'hc', algorithm.res = 'pc.stable', algorithm.max = 'hc', blacklist = NULL, whitelist = NULL, test = NULL, score = NULL, max.sx = NULL, ...){
  if (algorithm %in% c('hc','tabu')){
    get(paste(algorithm,'.DBN',sep=''))(data = data, score = score, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')){
    get(paste(algorithm,'.DBN',sep=''))(data = data, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm %in% c('h2pc','mmhc') & algorithm.res %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr') & algorithm.max %in% c('hc','tabu')){
    get(paste(algorithm,'.DBN',sep=''))(data = data, score = score, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else if (algorithm == 'rsmax2' & algorithm.res %in% c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr') & algorithm.max %in% c('hc','tabu')){
    get(paste(algorithm,'.DBN',sep=''))(data = data, restrict = algorithm.res, maximize = algorithm.max, score = score, test = test, blacklist = blacklist, whitelist = whitelist, max.sx = max.sx, ...)
  }
  else{
    stop("ERROR: Algorithm must be one of the following: 'hc','tabu','pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr', 'h2pc','mmhc', 'rsmax2'")
  }
}
