
# DBN definition
DBN_example <- empty_DBN(dynamic_nodes = c("A"), markov_order = 1)

# Node Addition
DBN_example <-
  add_node_DBN(DBN = DBN_example, node = "B", type = 'Dynamic')

# Arcs Addition
DBN_example <-
  add_arc_DBN(DBN = DBN_example,
              from = c('A', 't'),
              to = c('B', 't'))
DBN_example <-
  add_arc_DBN(DBN = DBN_example,
              from = c('A', 't_0'),
              to = c('B', 't_0'))
DBN_example <-
  add_arc_DBN(DBN = DBN_example,
              from = c('A', 't-1'),
              to = c('A', 't'))



# CPTs external definition
A_lv <- c('yes', 'no')
B_lv <- c('high', 'low')
A_0.prob <- array(c(0.2, 0.8),
                 dim = length(A_lv),
                 dimnames = list(A_0 = A_lv))
B_0.prob <- array(
  c(0.25, 0.75, 0.6, 0.4),
  dim = c(length(B_lv), length(A_lv)),
  dimnames = list(B_0 = B_lv, A_0 = A_lv)
)

dims <- list(A_t = A_lv)
dims[['A_t-1']] <- A_lv
A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                 dim = c(length(A_lv), length(A_lv)),
                 dimnames = dims)
B_t.prob <- array(
  c(0.2, 0.8, 0.6, 0.4),
  dim = c(length(B_lv), length(A_lv)),
  dimnames = list(B_t = B_lv, A_t = A_lv)
)
CPTs_toy <- list(
  A_0 = A_0.prob,
  B_0 = B_0.prob,
  A_t = A_t.prob,
  B_t = B_t.prob
)
# fitting the dbn network dbn.fit object
fitted_DBN <- DBN_parameters(DBN = DBN_example, CPTs = CPTs_toy)

# G_0 bn.fit object
bn_0 <-  from_fitted_DBN_to_fitted_G_0(fitted_DBN)
# G_transition bn.fit object
bn_transition <- from_fitted_DBN_to_fitted_G_transition(fitted_DBN)

nodes_time_t <- get_nodes_t(bn_transition)

# testing get_nodes_t
test_that("get_nodes_t returns a list of strings", {
  expect_true(all(sapply(nodes_time_t, is.character)))
})
test_that("get_nodes_t can extract specific nodes form G_transition",{
  nodes_list_t <- c("A_t","B_t")
  expect_equal(nodes_time_t,nodes_list_t)
})
test_that("get_nodes_t raise an error if G_transition isn't a bn.fit object",{
  transition_network <- c("A", "B")
  expect_error(get_nodes_t(transition_network))
})

# testing get_parents
test_that("get_parents returns the right parents set for a node", {
  expect_equal(get_parents(bn_0, "B_0"), c("A_0"))
  expect_equal(get_parents(bn_transition, "B_t"), c("A_t"))
  expect_equal(get_parents(bn_0, "A_0"), character(0))
})
test_that("get_parents raises an error if G is not a bn.fit object", {
  transition_network <- c("A", "B")
  expect_error(get_parents(transition_network, "A_t"))
})
test_that("get_parents returns NULL if the node does not exist", {
  expect_equal(get_parents(bn_transition, "H_t"), NULL)
})

# testing get_index_regular_expression
test_that("get_index_regular_expression returns the correct time index", {
  expect_equal(get_index_regular_expression("A_t-30"), 30)
  expect_equal(get_index_regular_expression("A_t-1"), 1)
  expect_equal(get_index_regular_expression("A_t"), 0)
  expect_equal(get_index_regular_expression("A_t30_t-1"), 1)
  expect_error(get_index_regular_expression("A"))
  expect_error(get_index_regular_expression(10))
})

# testing get_generic_node_name_rex
test_that("get_generic_node_name_rex returns the correct generic node name", {
  expect_equal(get_generic_node_name_rex("A_t-30"), "A")
  expect_equal(get_generic_node_name_rex("A_t30_t-20"), "A_t30")
  expect_equal(get_generic_node_name_rex("B_t-1"), "B")
  expect_equal(get_generic_node_name_rex("C_t"), "C")
  expect_error(get_generic_node_name_rex(10))
  expect_error(get_generic_node_name_rex("A_b"))
  expect_error(get_generic_node_name_rex("B"))
})

# testing get_generic_node_name_0
test_that("get_generic_node_name_0 returns the correct generic node name",{
  expect_equal(get_generic_node_name_0("H_0"), "H")
  expect_equal(get_generic_node_name_0("H_0N2_0"), "H_0N2")
  expect_error(get_generic_node_name_0("A_N"))
  expect_error(get_generic_node_name_0("A"))
})

# testing get_time_point
test_that("get_time_point returns the correct value", {
  expect_equal(get_time_point(10, 1), 10)
  expect_equal(get_time_point(10, 0), 11)
  expect_equal(get_time_point(10, 2), 9)
  expect_error(get_time_point("10", 2))
  expect_error(get_time_point(10, "2"))
  expect_error(get_time_point("10", "2"))
})

# testing filter_cpt
test_that("filter_cpt filters a cpt correctly", {
  cpt_to_filter <-
    data.frame(as.table(bn_transition[["B_t"]][["prob"]]))
  parents_values_list <- list("A_t" = c("yes"))
  filtered_cpt <- filter_cpt(cpt_to_filter, parents_values_list)
  expected_cpt <- cpt_to_filter[cpt_to_filter$A_t == "yes", ]
  expect_equal(filtered_cpt, expected_cpt)
  expect_error(filter_cpt(cpt_to_filter, c("A_t")))
  expect_error(filter_cpt(c("A_t"), parents_values_list))
})

# testing dbn_sampling
test_that("dbn_sampling raises error in case of wrong inputs", {
  expect_error(dbn_sampling(bn_0, bn_transition, "5", 4))
  expect_error(dbn_sampling(bn_0, bn_transition, 5, "4"))
  expect_error(dbn_sampling(bn_0, bn_transition, "5", "4"))
  expect_error(dbn_sampling(c("dbn"), 10, 10))
  expect_error(dbn_sampling(fitted_DBN, 10, 0))
  expect_error(dbn_sampling(fitted_DBN, 0, 2))
})

test_that(
  "dbn_sampling produces a dataset that accurately reflects the process from which it originates",
  {
    tolerance <- 0.02
    
    sampled_dataset <-
      dbn_sampling(fitted_DBN = fitted_DBN, 1e+4, 5)
    
    fitted_dbn_with_param_learning <-
      DBN_parameters(DBN = DBN_example, data = sampled_dataset)
    compare_with_tolerance <-
      function(value1, value2, tolerance)
        abs(value1 - value2) <= tolerance
    
    #testing for A_0
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_0$prob[["yes"]],
                             fitted_dbn_with_param_learning$A_0$prob[["yes"]],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_0$prob[["no"]],
                             fitted_dbn_with_param_learning$A_0$prob[["no"]],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing for B_0
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_0$prob["high", "yes"],
                             fitted_dbn_with_param_learning$B_0$prob["high", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_0$prob["high", "no"],
                             fitted_dbn_with_param_learning$B_0$prob["high", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_0$prob["low", "yes"],
                             fitted_dbn_with_param_learning$B_0$prob["low", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_0$prob["low", "no"],
                             fitted_dbn_with_param_learning$B_0$prob["low", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing for A_t
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_t$prob["yes", "yes"],
                             fitted_dbn_with_param_learning$A_t$prob["yes", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_t$prob["yes", "no"],
                             fitted_dbn_with_param_learning$A_t$prob["yes", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_t$prob["no", "yes"],
                             fitted_dbn_with_param_learning$A_t$prob["no", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$A_t$prob["no", "no"],
                             fitted_dbn_with_param_learning$A_t$prob["no", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing for B_t
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_t$prob["high", "yes"],
                             fitted_dbn_with_param_learning$B_t$prob["high", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_t$prob["high", "no"],
                             fitted_dbn_with_param_learning$B_t$prob["high", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_t$prob["low", "yes"],
                             fitted_dbn_with_param_learning$B_t$prob["low", "yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_DBN$B_t$prob["low", "no"],
                             fitted_dbn_with_param_learning$B_t$prob["low", "no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
  }
)
test_that(
  "dbn_sampling produces a dataset that accurately reflects the process from which it originates part.2",{
    # Creating a DBN
    my_dbn <-
      empty_DBN(dynamic_nodes = c("A", "B", "C"),
                markov_order = 1)
    
    # Creating G_0
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('A', 't_0'),
                  to = c('B', 't_0'))
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('A', 't_0'),
                  to = c('C', 't_0'))
    
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('B', 't_0'),
                  to = c('C', 't_0'))
    
    # Creating G_transition
    my_dbn <- add_arc_DBN(DBN = my_dbn,
                          from = c('A', 't'),
                          to = c('B', 't'))
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('A', 't'),
                  to = c('C', 't'))
    
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('B', 't'),
                  to = c('C', 't'))
    
    my_dbn <- add_arc_DBN(DBN = my_dbn,
                          from = c('A', 't-1'),
                          to = c('A', 't'))
    
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('B', 't-1'),
                  to = c('B', 't'))
    
    my_dbn <-
      add_arc_DBN(DBN = my_dbn,
                  from = c('C', 't-1'),
                  to = c('C', 't'))
    
    #Defining CPT G_0
    A_lv <- c('yes', 'no')
    B_lv <- c('high', 'low')
    C_lv <- c('medium', 'high')
    A_0.prob <- array(c(0.2, 0.8),
                      dim = length(A_lv),
                      dimnames = list(A_0 = A_lv))
    B_0.prob <- array(
      c(0.25, 0.75, 0.6, 0.4),
      dim = c(length(B_lv), length(A_lv)),
      dimnames = list(B_0 = B_lv, A_0 = A_lv)
    )
    C_0.prob <- array(
      c(0.1, 0.9, 0.5, 0.5, 0.3, 0.7, 0.6, 0.4),
      dim = c(length(C_lv), length(B_lv), length(A_lv)),
      dimnames = list(C_0 = C_lv, B_0 = B_lv, A_0 = A_lv)
    )
    
    # Defining CPT G_transition
    dims_A_t <- list(A_t = A_lv)
    dims_A_t[['A_t-1']] = A_lv
    # defining dims for B
    dims_B_t <- list(B_t = B_lv, A_t = A_lv)
    dims_B_t[['B_t-1']] = B_lv
    
    # defining dims for C
    dims_C_t <- list(C_t = C_lv, A_t = A_lv, B_t = B_lv)
    dims_C_t[['C_t-1']] = C_lv
    
    #defining A_t CPT
    A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                     dim = c(length(A_lv), length(A_lv)),
                     dimnames = dims_A_t)
    # defining B_t CPT
    B_t.prob = array(
      c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.1, 0.9),
      dim = c(length(B_lv), length(A_lv), length(B_lv)),
      dimnames = dims_B_t
    )
    #defining C_t CPT
    C_t.prob = array(
      c(
        0.1,
        0.9,
        0.2,
        0.8,
        0.12,
        0.88,
        0.9,
        0.1,
        0.3,
        0.7,
        0.4,
        0.6,
        0.25,
        0.75,
        0.45,
        0.55
      ),
      dim = c(length(C_lv), length(A_lv), length(B_lv), length(C_lv)),
      dimnames = dims_C_t
    )
    
    CPTs_mydbn = list(
      A_0 = A_0.prob,
      B_0 = B_0.prob,
      C_0 = C_0.prob,
      A_t = A_t.prob,
      B_t = B_t.prob,
      C_t = C_t.prob
    )
    #defining the two dbns
    fitted_my_dbn <- DBN_parameters(my_dbn, CPTs_mydbn)
    sampled_dataset <-
      dbn_sampling(fitted_DBN = fitted_my_dbn, 15e+3, 5)
    fitted_dbn_with_param_learning <-
      DBN_parameters(DBN = my_dbn, data = sampled_dataset)
    
    tolerance <- 0.04
    compare_with_tolerance <-
      function(value1, value2, tolerance)
        abs(value1 - value2) <= tolerance
    
    #testing A_0
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_0$prob[["yes"]],
                             fitted_dbn_with_param_learning$A_0$prob[["yes"]],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_0$prob[["no"]],
                             fitted_dbn_with_param_learning$A_0$prob[["no"]],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing B_0
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_0$prob["high","yes"],
                             fitted_dbn_with_param_learning$B_0$prob["high","yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_0$prob["high","no"],
                             fitted_dbn_with_param_learning$B_0$prob["high","no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_0$prob["low","yes"],
                             fitted_dbn_with_param_learning$B_0$prob["low","yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_0$prob["low","no"],
                             fitted_dbn_with_param_learning$B_0$prob["low","no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    
    #testing C_0 when A_0 = yes
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["medium","yes", "high"],
                             fitted_dbn_with_param_learning$C_0$prob["medium","yes","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["medium","yes","low"],
                             fitted_dbn_with_param_learning$C_0$prob["medium","yes","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["high","yes","high"],
                             fitted_dbn_with_param_learning$C_0$prob["high","yes","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["high","yes","low"],
                             fitted_dbn_with_param_learning$C_0$prob["high","yes","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing C_0 when A_0 = no
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["medium","no","high"],
                             fitted_dbn_with_param_learning$C_0$prob["medium","no","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["medium","no","low"],
                             fitted_dbn_with_param_learning$C_0$prob["medium","no","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["high","no","high"],
                             fitted_dbn_with_param_learning$C_0$prob["high","no","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_0$prob["high","no","low"],
                             fitted_dbn_with_param_learning$C_0$prob["high","no","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #testing G_transition
    
    #testing A_t
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_t$prob["no","no"],
                             fitted_dbn_with_param_learning$A_t$prob["no","no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_t$prob["no","yes"],
                             fitted_dbn_with_param_learning$A_t$prob["no","yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_t$prob["yes","no"],
                             fitted_dbn_with_param_learning$A_t$prob["yes","no"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$A_t$prob["yes","yes"],
                             fitted_dbn_with_param_learning$A_t$prob["yes","yes"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    # testing B_t
    # when B_t-1 = high
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["high","no","high"],
                             fitted_dbn_with_param_learning$B_t$prob["high","no","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["high","yes","high"],
                             fitted_dbn_with_param_learning$B_t$prob["high","yes","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["low","no","high"],
                             fitted_dbn_with_param_learning$B_t$prob["low","no","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["low","yes","high"],
                             fitted_dbn_with_param_learning$B_t$prob["low","yes","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    #when B_t-1 =low
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["high","no","low"],
                             fitted_dbn_with_param_learning$B_t$prob["high","no","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["high","yes","low"],
                             fitted_dbn_with_param_learning$B_t$prob["high","yes","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["low","no","low"],
                             fitted_dbn_with_param_learning$B_t$prob["low","no","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$B_t$prob["low","yes","low"],
                             fitted_dbn_with_param_learning$B_t$prob["low","yes","low"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    # testing C_t
    # when B_t = high and C_t-1=medium
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","yes","high","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","yes","high","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","no","high","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","no","high","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","yes","high","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["high","yes","high","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","no","high","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["high","no","high","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    # when B_t = low and C_t-1=medium
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","yes","low","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","yes","low","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","no","low","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","no","low","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","yes","low","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["high","yes","low","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","no","low","medium"],
                             fitted_dbn_with_param_learning$C_t$prob["high","no","low","medium"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    # when B_t = high, C_t-1 = high
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","yes","high","high"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","yes","high","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","no","high","high"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","no","high","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","yes","high","high"],
                             fitted_dbn_with_param_learning$C_t$prob["high","yes","high","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","no","high","high"],
                             fitted_dbn_with_param_learning$C_t$prob["high","no","high","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    # when B_t = low, C_t-1 = high
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","yes","low","high"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","yes","low","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["medium","no","low","high"],
                             fitted_dbn_with_param_learning$C_t$prob["medium","no","low","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","yes","low","high"],
                             fitted_dbn_with_param_learning$C_t$prob["high","yes","low","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
    
    is_equal <-
      compare_with_tolerance(fitted_my_dbn$C_t$prob["high","no","low","high"],
                             fitted_dbn_with_param_learning$C_t$prob["high","no","low","high"],
                             tolerance)
    expect_equal(is_equal, TRUE)
  })
