#' Return the variables indexes by t in G_transition
#'
#' @param G_transition a G_transition graph.
#' @returns nodes indexed by t only e.g. A_t, B_t.
#' @export
#' @examples
#' get_nodes_t(G_transition)
get_nodes_t <- function(G_transition) {
  if (class(G_transition) != "bn.fit") {
    stop("G_transition must be a bn.fit object!")
  }
  nodes_dbn <- bnlearn::node.ordering(G_transition)
  nodes_t <- c()
  for (n in nodes_dbn) {
    ends_with_t <- substr(n, nchar(n), nchar(n)) == "t"
    if (ends_with_t == TRUE) {
      nodes_t <- c(nodes_t, n)
    }
  }
  return(nodes_t)
}
#' Get the parents of a node
#'
#' @param G A bn.fit object
#' @param n name of node n
#' @returns parents of the node n.
#' @export
#' @examples
#' get_parents(G, "A_t")
get_parents <- function(G_transition, n) {
  if (class(G_transition) != "bn.fit") {
    stop("G must be a bn.fit object!")
  }
  #given a node, returns its parents
  parents_n <- G_transition[[n]][["parents"]]
  return(parents_n)
}
#' Get the time index of a node
#'
#' @param n name of node n
#' @returns time index of a node n.
#' @export
#' @examples
#' get_index_regular_expression("A_t-1")
get_index_regular_expression <- function(n) {
  if (is.character(n) == FALSE) {
    stop("node name must be a string!")
  }
  time_index <- substr(n, nchar(n), nchar(n))
  if (time_index == "t") {
    return(0)
  } else{
    # Define a regular expression pattern
    pattern <- ".*_t-(\\d+)$"
    # Extract the last n_characters from each string
    time_index <- sub(pattern, "\\1", n)
    
    #if time index isn't found will return n
    if (time_index == n) {
      string_error <- paste("Time index can't be found for node ", n)
      stop(string_error)
    } else{
      # Convert the result to numeric
      time_index <- as.numeric(time_index)
      return(time_index)
    }
  }
}
#' Get the generic name of a node belonging to the transition graph
#'
#' @param n name of node n
#' @returns name of node with no time index.
#' @export
#' @examples
#' get_generic_node_name_rex("A_t-1")
get_generic_node_name_rex_old<- function(n) {
  if (is.character(n) == FALSE) {
    stop("node name must be a string!")
  }
  # Define a regular expression pattern
  pattern <- "^(.*)_t(?:-(\\d+))?$"
  # Extract the substring preceding _t_number
  generic_node_name <- sub(pattern, "\\1", n)
  if (generic_node_name == n) {
    string_error <- paste("Generic node name can't be found for ", n)
    stop(string_error)
  } else{
    return(generic_node_name)
  }
}
#' Get the generic name of a node belonging to the transition graph
#'
#' @param n name of node n
#' @returns name of node with no time index.
#' @export
#' @examples
#' get_generic_node_name_rex("A_t-1")
get_generic_node_name_rex<- function(n) {
  if (is.character(n) == FALSE) {
    stop("node name must be a string!")
  }
  # Define a regular expression pattern
  pattern <- "^(.*)_t(?:-(\\d+))?$"
  # Extract the substring preceding _t_number
  generic_node_name <- stringi::stri_replace_first_regex(n, pattern, "$1")
  if (generic_node_name == n) {
    string_error <- paste("Generic node name can't be found for ", n)
    stop(string_error)
  } else{
    return(generic_node_name)
  }
}
#' Get the generic name of a node belonging to the graph G_0
#'
#' @param n name of node n
#' @returns name of node with no time index.
#' @export
#' @examples
#' get_generic_node_name_0("A_0")
get_generic_node_name_0 <- function(n) {
  if (is.character(n) == FALSE) {
    stop("Node name must be a string!")
  }
  last_two_chars <- substring(n, nchar(n) - 1, nchar(n))
  if (last_two_chars == n) {
    string_error <- paste("Generic node name can't be found for ", n)
    stop(string_error)
  }
  if (!last_two_chars == "_0") {
    stop("Inconsistency found in node name. Name must end with _0")
  }
  generic_node_name <- substr(n, 1, nchar(n) - 2)
  return(generic_node_name)
}
#' Get the time point to retrieve the values of variables in the time series
#'
#' @param t time point
#' @param index_n time index of node n
#' @returns time-point in the series used to retrieve the value of the node n
#' @export
#' @examples
#' get_time_point(3,1)
get_time_point <- function(t, index_n) {
  if (is.character(t) | is.character(index_n)) {
    stop("t and index_n must be numeric")
  }
  # compute t required to get the value of the variable
  # one unit is added cause the series starts from zero
  time_index <- as.integer(t) - as.integer(index_n) + 1
  return(time_index)
}
#' Filter the cpt of node n given the value of its parents
#'
#' @param df_cpt_n a dataframe representing the cpt
#' @param parents_values a list with the parents of n and their values
#' @returns the filtered cpt
#' @export
#' @examples
#' filter_cpt(df_cp_n, parents_values)
filter_cpt_old <- function(df_cpt_n, parents_values) {
  if (is.data.frame(df_cpt_n) == FALSE) {
    stop("df_cpt_n must be a dataframe!")
  }
  if (is.list(parents_values) == FALSE) {
    stop("parents_values must be a list!")
  }
  filtered_df <- df_cpt_n
  for (key in names(parents_values)) {
    filtered_df <-
      filtered_df[filtered_df[key] == parents_values[key],]
  }
  return(filtered_df)
}

#' Filter the cpt of node n given the value of its parents
#'
#' @param df_cpt_n a dataframe representing the cpt
#' @param parents_values a list with the parents of n and their values
#' @returns the filtered cpt
#' @export
#' @examples
#' filter_cpt(df_cp_n, parents_values)
filter_cpt<- function(df_cpt_n, parents_values) {
  # Check if df_cpt_n is a dataframe
  if (!is.data.frame(df_cpt_n)) {
    stop("df_cpt_n must be a dataframe!")
  }
  
  #Check if parents_values is a list
  if (!is.list(parents_values)) {
    stop("parents_values must be a list!")
  }
  
  # Create a logical index vector for all conditions at once
  filter_index <- rep(TRUE, nrow(df_cpt_n))
  
  for (key in names(parents_values)) {
    filter_index <- filter_index & (df_cpt_n[[key]] == parents_values[[key]])
  }
  
  # Apply the combined logical index to filter the dataframe in one go
  filtered_df <- df_cpt_n[filter_index, ]
  
  return(filtered_df)
}


#' Generate a sampling dataset
#'
#' @param fitted_DBN an object of class 'DBN'
#' @param N_samples number of samples
#' @param Time time series length
#' @returns the generated dataframe
#' @export
#' @examples
#' dbn_sampling(DBN_example, N_samples, Time)
dbn_sampling_old <- function(fitted_DBN, N_samples, Time) {
  if (is.character(Time)) {
    stop("Time must be an integer!")
  }
  if (Time < 1){
    stop("Time must be greater than 0!")
  }
  if (is.character(N_samples)) {
    stop("N_samples must be an integer!")
  }
  if (N_samples < 1){
    stop("N_samples must be greater than 0!")
  }
  if (class(fitted_DBN) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }
  bn_0 <- from_fitted_DBN_to_fitted_G_0(fitted_DBN)
  bn_transition <-
    from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
  # numbers in input are of type double
  N_samples <- as.integer(N_samples)
  Time <- as.integer(Time)
  timeseries_dict <- list(Time = c())
  for (obs in seq(N_samples)) {
    #when t=0
    timeseries_dict[["Time"]] <- c(timeseries_dict[["Time"]], 0)
    timeseries_dict[["Sample_id"]] <-
      c(timeseries_dict[["Sample_id"]],
        paste("sample", obs, sep = ""))
    for (n in  bnlearn::node.ordering(bn_0)) {
      parents_n <- get_parents(bn_0, n)
      #has parents
      if (length(parents_n) > 0) {
        #this isn't a root node we need to check for the values of its parents
        #create a dictionary for parents' values
        parents_values <- list()
        #iterate over parents_n
        for (p in parents_n) {
          #remove _0 index
          p_with_no_index <- get_generic_node_name_0(p)
          #get the value of the parent from the timeseries
          value_p <-
            timeseries_dict[[p_with_no_index]][(obs - 1) * (Time + 1) + 1]
          parents_values[[p]] <- value_p
        }
        #select and filtering cpt
        df_cp_n <- data.frame(as.table(bn_0[[n]][["prob"]]))
        filtered_cpt_n <- filter_cpt_old(df_cp_n, parents_values)
        #sampling n using frequencies
        sampling_n <-
          as.character(sample(filtered_cpt_n[[n]], size = 1,
                              prob = filtered_cpt_n[["Freq"]]))
        #adding generated sample to the dataset
        n_with_no_index <- get_generic_node_name_0(n)
        timeseries_dict[[n_with_no_index]] <-
          c(timeseries_dict[[n_with_no_index]],
            sampling_n)
      } else{
        #the current node is a root node
        # remove _0 index
        n_with_no_index <- get_generic_node_name_0(n)
        # get the cpt of the node n
        df_prob_n <- data.frame(as.table(bn_0[[n]][["prob"]]))
        sampled_value <-
          sample(df_prob_n[[n]], size = 1, prob = df_prob_n[["Freq"]])
        timeseries_dict[[n_with_no_index]] <-
          c(timeseries_dict[[n_with_no_index]],
            as.character(sampled_value))
        
      }
    }
    #when t>0
    #iterate over time
    for (t in 1:Time) {
      #get the nodes for the current order
      nodes <- get_nodes_t(bn_transition)
      for (n in nodes) {
        #retrieve the parents for the node-th
        parents_n = get_parents(bn_transition, n)
        #two cases could occur parents_n == 0 or parents_n > 0
        if (length(parents_n) > 0) {
          #create a dictionary for parents' values
          parents_values <- list()
          #iterate over parents_n
          for (p in parents_n) {
            #get the index of the parent to retrieve its value
            index_p <- get_index_regular_expression(p)
            p_with_no_index <- get_generic_node_name_rex_old(p)
            #get the time point for the parent in order to retrieve the
            #value from the time series
            time_point <- get_time_point(t, index_p)
            value_p <-
              timeseries_dict[[p_with_no_index]][(obs - 1) * (Time + 1) + time_point]
            parents_values[[p]] <- value_p
          }
          #select and filtering cpt
          df_cp_n <-
            data.frame(as.table(bn_transition[[n]][["prob"]]))
          #data.frame replace - with dot
          #replacing dot with dash in dataframe columns
          colnames(df_cp_n) <- gsub("\\.", "-", colnames(df_cp_n))
          filtered_cpt_n <- filter_cpt_old(df_cp_n, parents_values)
          #sampling n using frequencies
          sampling_n <-
            as.character(sample(filtered_cpt_n[[n]], size = 1,
                                prob = filtered_cpt_n[["Freq"]]))
          #adding the sample to the dataset
          n_with_no_index <- get_generic_node_name_rex_old(n)
          timeseries_dict[[n_with_no_index]] <-
            c(timeseries_dict[[n_with_no_index]],
              sampling_n)
        } else {
          #get the cpt of the node
          df_cp_n <-
            data.frame(as.table(bn_transition[[n]][["prob"]]))
          #sampling n using frequencies
          sampling_n <- as.character(sample(df_cp_n[[n]], size = 1,
                                            prob = df_cp_n[["Freq"]]))
          #adding the sample to the dataset
          n_with_no_index <- get_generic_node_name_rex_old(n)
          timeseries_dict[[n_with_no_index]] <-
            c(timeseries_dict[[n_with_no_index]],
              sampling_n)
        }
      }
    }
    timeseries_dict[["Time"]] <-
      c(timeseries_dict[["Time"]], seq(Time))
    timeseries_dict[["Sample_id"]] <-
      c(timeseries_dict[["Sample_id"]],
        rep(paste("sample", obs, sep = ""), Time))
  }
  df_timeseries <- data.frame(timeseries_dict)
  return(df_timeseries)
}





#' Generate a sampling dataset
#'
#' @param fitted_DBN an object of class 'DBN'
#' @param N_samples number of samples
#' @param Time time series length
#' @returns the generated dataframe
#' @export
#' @examples
#' dbn_sampling(DBN_example, N_samples, Time)
dbn_sampling <- function(fitted_DBN, N_samples, Time) {
  if (is.character(Time)) {
    stop("Time must be an integer!")
  }
  if (Time < 1){
    stop("Time must be greater than 0!")
  }
  if (is.character(N_samples)) {
    stop("N_samples must be an integer!")
  }
  if (N_samples < 1){
    stop("N_samples must be greater than 0!")
  }
  if (class(fitted_DBN) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }
  bn_0 <- from_fitted_DBN_to_fitted_G_0(fitted_DBN)
  bn_transition <-
    from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
  # numbers in input are of type double
  N_samples <- as.integer(N_samples)
  Time <- as.integer(Time)
  timeseries_dict <- list(Time = c())
  for (obs in seq(N_samples)) {
    #when t=0
    timeseries_dict[["Time"]] <- c(timeseries_dict[["Time"]], 0)
    timeseries_dict[["Sample_id"]] <-
      c(timeseries_dict[["Sample_id"]],
        paste("sample", obs, sep = ""))
    for (n in  bnlearn::node.ordering(bn_0)) {
      parents_n <- get_parents(bn_0, n)
      #has parents
      if (length(parents_n) > 0) {
        #this isn't a root node we need to check for the values of its parents
        #create a dictionary for parents' values
        parents_values <- list()
        #iterate over parents_n
        for (p in parents_n) {
          #remove _0 index
          p_with_no_index <- get_generic_node_name_0(p)
          #get the value of the parent from the timeseries
          value_p <-
            timeseries_dict[[p_with_no_index]][(obs - 1) * (Time + 1) + 1]
          parents_values[[p]] <- value_p
        }
        #select and filtering cpt
        df_cp_n <- data.frame(as.table(bn_0[[n]][["prob"]]))
        filtered_cpt_n <- filter_cpt(df_cp_n, parents_values)
        #sampling n using frequencies
        sampling_n <-
          as.character(sample(filtered_cpt_n[[n]], size = 1,
                              prob = filtered_cpt_n[["Freq"]]))
        #adding generated sample to the dataset
        n_with_no_index <- get_generic_node_name_0(n)
        timeseries_dict[[n_with_no_index]] <-
          c(timeseries_dict[[n_with_no_index]],
            sampling_n)
      } else{
        #the current node is a root node
        # remove _0 index
        n_with_no_index <- get_generic_node_name_0(n)
        # get the cpt of the node n
        df_prob_n <- data.frame(as.table(bn_0[[n]][["prob"]]))
        sampled_value <-
          sample(df_prob_n[[n]], size = 1, prob = df_prob_n[["Freq"]])
        timeseries_dict[[n_with_no_index]] <-
          append(timeseries_dict[[n_with_no_index]],
            as.character(sampled_value))
        
      }
    }
    #when t>0
    #iterate over time
    for (t in 1:Time) {
      #get the nodes for the current order
      nodes <- get_nodes_t(bn_transition)
      for (n in nodes) {
        #retrieve the parents for the node-th
        parents_n = get_parents(bn_transition, n)
        #two cases could occur parents_n == 0 or parents_n > 0
        if (length(parents_n) > 0) {
          #create a dictionary for parents' values
          parents_values <- list()
          #iterate over parents_n
          for (p in parents_n) {
            #get the index of the parent to retrieve its value
            index_p <- get_index_regular_expression(p)
            p_with_no_index <- get_generic_node_name_rex(p)
            #get the time point for the parent in order to retrieve the
            #value from the time series
            time_point <- get_time_point(t, index_p)
            value_p <-
              timeseries_dict[[p_with_no_index]][(obs - 1) * (Time + 1) + time_point]
            parents_values[[p]] <- value_p
          }
          #select and filtering cpt
          df_cp_n <-
            data.frame(as.table(bn_transition[[n]][["prob"]]))
          #data.frame replace - with dot
          #replacing dot with dash in dataframe columns
          colnames(df_cp_n) <- gsub("\\.", "-", colnames(df_cp_n))
          filtered_cpt_n <- filter_cpt(df_cp_n, parents_values)
          #sampling n using frequencies
          sampling_n <-
            as.character(sample(filtered_cpt_n[[n]], size = 1,
                                prob = filtered_cpt_n[["Freq"]]))
          #adding the sample to the dataset
          n_with_no_index <- get_generic_node_name_rex(n)
          timeseries_dict[[n_with_no_index]] <-
            append(timeseries_dict[[n_with_no_index]],
              sampling_n)
        } else {
          #get the cpt of the node
          df_cp_n <-
            data.frame(as.table(bn_transition[[n]][["prob"]]))
          #sampling n using frequencies
          sampling_n <- as.character(sample(df_cp_n[[n]], size = 1,
                                            prob = df_cp_n[["Freq"]]))
          #adding the sample to the dataset
          n_with_no_index <- get_generic_node_name_rex(n)
          timeseries_dict[[n_with_no_index]] <-
            c(timeseries_dict[[n_with_no_index]],
              sampling_n)
        }
      }
    }
    timeseries_dict[["Time"]] <-
      c(timeseries_dict[["Time"]], seq(Time))
    timeseries_dict[["Sample_id"]] <-
      c(timeseries_dict[["Sample_id"]],
        rep(paste("sample", obs, sep = ""), Time))
  }
  df_timeseries <- data.frame(timeseries_dict)
  return(df_timeseries)
}



#' Generate a sampling dataset
#'
#' @param fitted_DBN an object of class 'DBN'
#' @param observations a list containing the observations for each variable
#' @param timepoints number of timepoints to be forecasted
#' @returns a dataframe containing the forecasting
#' @export
#' @examples
#' dbn_forecasting(fitted_DBN, observations, timepoints)
dbn_forecasting <- function(fitted_DBN, observations, timepoints) {
  if (is.character(timepoints)) {
    stop("timepoints must be an integer!")
  }
  if (timepoints < 1){
    stop("timepoints must be greater than 0!")
  }
  #check that observation is a list
  if (class(fitted_DBN) != "dbn.fit") {
    stop("fitted_DBN must be a dbn.fit object")
  }
  bn_transition <-
    from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
  #when t>0
  #iterate over time
  for (t in 1:timepoints) {
    observations[["Time"]] <-
      c(observations[["Time"]], observations[["Time"]][1]+t)
    #get the nodes for the current order
    nodes <- get_nodes_t(bn_transition)
    for (n in nodes) {
      #retrieve the parents for the node-th
      parents_n = get_parents(bn_transition, n)
      #two cases could occur parents_n == 0 or parents_n > 0
      if (length(parents_n) > 0) {
        #create a dictionary for parents' values
        parents_values <- list()
        #iterate over parents_n
        for (p in parents_n) {
          #get the index of the parent to retrieve its value
          index_p <- get_index_regular_expression(p)
          p_with_no_index <- get_generic_node_name_rex(p)
          #get the time point for the parent in order to retrieve the
          #value from the time series
          time_point <- get_time_point(t, index_p)
          value_p <-
            observations[[p_with_no_index]][(1 - 1) * (timepoints + 1) + time_point]
          parents_values[[p]] <- value_p
        }
        #select and filtering cpt
        df_cp_n <-
          data.frame(as.table(bn_transition[[n]][["prob"]]))
        #data.frame replace - with dot
        #replacing dot with dash in dataframe columns
        colnames(df_cp_n) <- gsub("\\.", "-", colnames(df_cp_n))
        filtered_cpt_n <- filter_cpt(df_cp_n, parents_values)
        #sampling n using frequencies
        sampling_n <-
          as.character(sample(filtered_cpt_n[[n]], size = 1,
                              prob = filtered_cpt_n[["Freq"]]))
        #adding the sample to the dataset
        n_with_no_index <- get_generic_node_name_rex(n)
        observations[[n_with_no_index]] <-
          append(observations[[n_with_no_index]],
                 sampling_n)
      } else {
        #get the cpt of the node
        df_cp_n <-
          data.frame(as.table(bn_transition[[n]][["prob"]]))
        #sampling n using frequencies
        sampling_n <- as.character(sample(df_cp_n[[n]], size = 1,
                                          prob = df_cp_n[["Freq"]]))
        #adding the sample to the dataset
        n_with_no_index <- get_generic_node_name_rex(n)
        observations[[n_with_no_index]] <-
          c(observations[[n_with_no_index]],
            sampling_n)
      }
    }
  }
df_timeseries <- data.frame(observations)
return(df_timeseries)
}






