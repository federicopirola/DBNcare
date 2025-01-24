
my_dbn <-
  empty_DBN(dynamic_nodes = c("A", "B", "C", "D", "E","F"),
            markov_order = 1)

# Creating G_0


my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('A', 't_0'),
              to = c('D', 't_0'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('A', 't_0'),
              to = c('C', 't_0'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('B', 't_0'),
              to = c('C', 't_0'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('B', 't_0'),
              to = c('D', 't_0'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('D', 't_0'),
              to = c('E', 't_0'))
my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('E', 't_0'),
              to = c('F', 't_0'))


# Creating G_transition

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('A', 't'),
                      to = c('D', 't'))
my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('A', 't'),
              to = c('C', 't'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('B', 't'),
              to = c('C', 't'))

my_dbn <-
  add_arc_DBN(DBN = my_dbn,
              from = c('B', 't'),
              to = c('D', 't'))

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('D', 't'),
                      to = c('E', 't'))

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('E', 't'),
                      to = c('F', 't'))



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

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('D', 't-1'),
                      to = c('D', 't'))

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('E', 't-1'),
                      to = c('E', 't'))

my_dbn <- add_arc_DBN(DBN = my_dbn,
                      from = c('F', 't-1'),
                      to = c('F', 't'))

# test remove_suffix_t_minus_1
test_that("remove_suffix_t_minus_1 returns the correct output",{
  expect_equal(remove_suffix_t_minus_1("A_t-1"), "A_t")
  expect_equal(remove_suffix_t_minus_1("A_t_t-1"), "A_t_t")
  expect_equal(remove_suffix_t_minus_1("B_t-1"),"B_t")
  expect_equal(remove_suffix_t_minus_1(c("B_t-1","A_t-1")),c("B_t","A_t"))
  expect_error(remove_suffix_t_minus_1(1))
  
})

# test remove_suffix_t
test_that("remove_suffix_t returns the correct output",{
  expect_equal(remove_suffix_t(c("A_t","B_t")), c("A","B"))
  expect_equal(remove_suffix_t("A_t_t"), "A_t")
  expect_equal(remove_suffix_t("B_t"), "B")
  expect_error(remove_suffix_t(1))
})

# test select_nodes_names_with_t
test_that("select_nodes_names_with_t returns the right output",{
  expect_equal(select_nodes_names_with_t(c("A_t", "A_t-1")), c("A_t"))  
  expect_equal(select_nodes_names_with_t(c("A_t_t", "A_t-1")), c("A_t_t"))  
  expect_error(select_nodes_names_with_t(c(1,2,3,4)))
  
})

# test select_nodes_names_with_t_minus_n
test_that("select_nodes_names_with_t_minus_n returns the right output",{
  expect_equal(select_nodes_names_with_t_minus_n(c("A_t","A_t-1"),1), c("A_t-1"))
  expect_equal(select_nodes_names_with_t_minus_n(c("A_t","B_t-1"),1), c("B_t-1"))
  expect_equal(select_nodes_names_with_t_minus_n(c("A_t","B_t-2"),2), c("B_t-2"))
  expect_equal(select_nodes_names_with_t_minus_n(c("A_t","B_t_t-2"),2), c("B_t_t-2"))
  expect_error(select_nodes_names_with_t_minus_n(c("A_t","B_t-1"),"1"))
})

# test dynamic_ordering
test_that("dynamic_ordering returns the right output",{
  result_ordering <- c("A_t", "B_t","C_t","D_t","E_t","F_t")
  expect_equal(dynamic_ordering(my_dbn), result_ordering)
  expect_error(dynamic_ordering(c("a")))
})

# test node_levels
test_that("node_levels returns the correct levels for each node",{
  ordered_nodes_t <- dynamic_ordering(my_dbn)
  nodes_levels_t <- node_levels(my_dbn,ordered_nodes_t)
  expected_result <- matrix(c("A_t", "B_t", "C_t", "D_t", "E_t", "F_t",
                              "1", "1", "2", "2", "3", "4"),
                            nrow = 2,      
                            ncol = 6,      
                            byrow = TRUE) 
  expect_equal(nodes_levels_t, expected_result)
  expect_error(node_levels(c("not the right input"), ordered_nodes_t))
})

# test nodes_levels_g_0 
test_that("node_levels_g_0 returns the correct levels for each node",{
  g_0 <- from_DBN_to_G_0(my_dbn)
  nodes_uniq <- bnlearn::node.ordering(g_0)
  nodes_levels_g_0 <- node_levels_g_0(g_0, nodes_uniq)[2,]
  expected_result <- c("1", "1", "2", "2", "3", "4")
  expect_equal(nodes_levels_g_0, expected_result)
  
  expect_error(node_levels_g_0(c("wrong input"), nodes_uniq))
})

# test plot_g_0
test_that("plot_g0 generates error in case of wrong input",{
  expect_error(plot_g0(list(a=c("wrong input"))))
})

#plot_g_transition
test_that("plot_g_transition generates error in case of wrong input",{
  expect_error(plot_g_transition(list(a=c("wrong input"))))
})

#test acc_succession
test_that("acc_successions returns the correct x_axis coordinates",{
  ordered_nodes_t <- dynamic_ordering(my_dbn)
  nodes_levels_t <- node_levels(my_dbn,ordered_nodes_t)
  expected_result <- c(1, 2, 1, 2, 1, 1)
  x_positions <- acc_successions(as.numeric(nodes_levels_t[2,]))
  expect_equal(x_positions, expected_result)
  
})

