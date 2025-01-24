
library(dplyr)
test_that("get_node_name_vector returns the correct output", {
  expect_equal(get_node_name_vector("A_t-1"), c("A", "t-1"))
  expect_equal(get_node_name_vector("B_0"), c("B", "t_0"))
  expect_equal(get_node_name_vector("C_t"), c("C", "t"))
})
test_that("get_node_name_vector raises error when wrong input is given", {
  expect_error(get_node_name_vector(1))
})
test_that("generate_dbn_nodes_names returns the correct output", {
  expect_equal(generate_dbn_nodes_names(c("A", "B", "C")),
               list(
                 g_0 = c("A_0", "B_0", "C_0"),
                 g_t_1 = c("A_t-1", "B_t-1", "C_t-1"),
                 g_t = c("A_t", "B_t", "C_t")
               ))
})
test_that("generate_dbn_nodes_names raises error when wrong input is given",
          {
            #input not a vector
            expect_error(generate_dbn_nodes_names("A"))
            expect_error(generate_dbn_nodes_names(c(1, 2)))
          })
test_that("generate_dbn_random_structure returns the correct output", {
  check_correspondence <- function(dbn_edges) {
    # Extract pairs ending with _0 and _t
    pairs_0 <-
      unique(dbn_edges[grepl("_0$", dbn_edges[, 1]) &
                         grepl("_0$", dbn_edges[, 2]),])
    pairs_t <-
      unique(dbn_edges[grepl("_t$", dbn_edges[, 1]) &
                         grepl("_t$", dbn_edges[, 2]),])
    
    # Function to create a set of pairs without suffix
    remove_suffix <- function(pairs) {
      unique(paste0(
        sub("_0$|_t$", "", pairs[, 1]),
        "_",
        sub("_0$|_t$", "", pairs[, 2])
      ))
    }
    
    # Create sets of pairs without suffix
    set_0 <- remove_suffix(pairs_0)
    set_t <- remove_suffix(pairs_t)
    
    # Check if all elements in set_0 are in set_t
    all(set_0 %in% set_t)
  }
  #g_0 does contain the same edges of g_transition when is_same = TRUE
  dbn_random_generated <-
    generate_dbn_random_structure(c("A", "B", "C", "D"), TRUE, 0.5, 0.5)
  dbn_edges <- dbn_random_generated$arcs
  check_result <- check_correspondence(dbn_edges)
  expect_equal(check_result, TRUE)
  #g_0 does not contain the same edges of g_transition when is_same = FALSE
  dbn_random_generated <-
    generate_dbn_random_structure(c("A", "B", "C", "D"), FALSE, 0.5, 0.5)
  dbn_edges <- dbn_random_generated$arcs
  check_result <- check_correspondence(dbn_edges)
  expect_equal(check_result, FALSE)
  
  #the function returns a DBN object
  expect_equal(class(generate_dbn_random_structure(c("A", "B"), TRUE, 0.5, 0.5)), "DBN")
  
  #no error when the input parameters are correct
  expect_no_error(generate_dbn_random_structure(c("A", "B", "C", "D"), TRUE, 0.6, 0.6))
  expect_no_error(generate_dbn_random_structure(c("A", "B", "C", "D"), FALSE, 0.6, 0.6))
})
test_that("generate_dbn_random_structure raises error when wrong input is given",
          {
            #input not a vector
            expect_error(generate_dbn_random_structure("A", TRUE, 0.5, 0.5))
            #input is not a char vector
            expect_error(generate_dbn_random_structure(c(1, 2), TRUE, 0.5, 0.5))
            #is_same is not a boolean
            expect_error(generate_dbn_random_structure(c("A", "B"), "Not a boolean", 0.5, 0.5))
            #edge_prob is higher than 1
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 10, 10))
            #edge_prob is zero
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 0, 0))
            #edge_prob is higher than 1
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 0.5, 10))
            #edge_prob is zero
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 0.5, 0))
            #edge_prob is higher than 1
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 10, 0.5))
            #edge_prob is zero
            expect_error(generate_dbn_random_structure(c("A", "B"), TRUE, 0, 0.5))
            
          })

test_that("generate_dbn_nodes_distributions raises error when wrong input is given",
          {
            #generated_dbn is not a DBN object
            expect_error(generate_dbn_nodes_distributions("A", TRUE, 2))
            generated_random_dbn <-
              generate_dbn_random_structure(c("A", "B"), TRUE, 0.5, 0.5)
            #fixed_cardinality not a boolean
            expect_error(generate_dbn_nodes_distributions(generated_random_dbn, "Not boolean", 2))
            #max_variables_cardinality <2
            expect_error(generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 1))
            
          })

test_that("generate_dbn_nodes_distributions correctly istantiates nodes' distribution",
          {
            #no errors is produced when the input is correct
            generated_random_dbn <-
              generate_dbn_random_structure(c("A", "B", "C"), TRUE, 0.6, 0.6)
            expect_no_error(generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 2))
            #variables levels are aligned with max_cardinality when fixed_cardinality = TRUE
            generated_random_dbn <-
              generate_dbn_random_structure(c("A", "B", "C"), TRUE, 0.6, 0.6)
            fitted_random_dbn <-
              generate_dbn_nodes_distributions(generated_random_dbn, TRUE, 2)
            all_variables <- names(fitted_random_dbn)
            correct_instantiation <- TRUE
            for (var in all_variables) {
              if (length(dimnames(fitted_random_dbn$B_t$prob)[[1]]) != 2) {
                correct_instantiation <- FALSE
              }
            }
            expect_equal(correct_instantiation, TRUE)
            #variables levels are aligned with max_cardinality when fixed_cardinality = FALSE
            fitted_random_dbn <-
              generate_dbn_nodes_distributions(generated_random_dbn, FALSE, 4)
            all_variables <- names(fitted_random_dbn)
            correct_instantiation <- TRUE
            for (var in all_variables) {
              if (!(length(dimnames(fitted_random_dbn[[var]][["prob"]])[[1]]) >= 2 &
                    length(dimnames(fitted_random_dbn[[var]][["prob"]])[[1]]) <= 4)) {
                correct_instantiation <- FALSE
              }
            }
            expect_equal(correct_instantiation, TRUE)
          })
