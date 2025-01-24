test_that("empty_DBN returns the expected S3 class and values", {
  d_nodes <- c('A', 'B', 'C')
  s_nodes <- c('D', 'E')
  m_os <- 5
  for (m_o in 1:m_os) {
    expect_error(empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = m_o
    ),
    NA)
    dbn <-
      empty_DBN(
        dynamic_nodes = d_nodes,
        static_nodes = s_nodes,
        markov_order = m_o
      )
    expect_equal(class(dbn), 'DBN')
    expect_equal(names(dbn$nodes), c(d_nodes, s_nodes))
    temps <- c('type', 't_0', 't')
    expect_equal(dbn$markov_order, m_o)
    for (i in 1:m_o) {
      temps <- c(temps, gsub(" ", "", paste("t-", i)))
    }
    for (node in d_nodes) {
      expect_equal(names(dbn$nodes[[node]]), temps)
    }
    for (node in s_nodes) {
      expect_equal(names(dbn$nodes[[node]]), c('type', 't_0'))
    }
  }
  expect_error(
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 0
    ),
    'Markov order of the process must be 1 or higher!!!'
  )
})

test_that("add_node_DBN returns the expected output", {
  d_nodes <- c('A', 'B', 'C')
  s_nodes <- c('D', 'E')
  added_dynamic <- 'F'
  added_static <- 'G'
  dbn <-
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 1
    )
  expect_error(add_node_DBN(DBN = dbn, node = added_dynamic), NA)
  dbn <- add_node_DBN(DBN = dbn, node = added_dynamic)
  expect_equal(names(dbn$nodes), c(d_nodes, s_nodes, added_dynamic))
  expect_equal(dbn$nodes[[added_dynamic]]$type, 'Dynamic')
  expect_equal(names(dbn$nodes[[added_dynamic]]), c('type', 't_0', 't', 't-1'))
  expect_error(add_node_DBN(DBN = dbn, node = added_static, type = 'Static'), NA)
  dbn <-
    add_node_DBN(DBN = dbn, node = added_static, type = 'Static')
  expect_equal(names(dbn$nodes),
               c(d_nodes, s_nodes, added_dynamic, added_static))
  expect_equal(dbn$nodes[[added_static]]$type, 'Static')
  expect_equal(names(dbn$nodes[[added_static]]), c('type', 't_0'))
  expect_error(add_node_DBN(DBN = dbn, node = 1), 'ERROR: node is not a character')
  expect_error(add_node_DBN(DBN = c(1, 2), node = 1),
               "ERROR: DBN argument is not of class 'DBN'")
  expect_error(add_node_DBN(DBN = dbn, node = 'A'),
               'ERROR: node name already exists')
})

test_that("add_arc_DBN returns the expected output", {
  d_nodes <- c('A', 'B', 'C')
  s_nodes <- c('D', 'E')
  dbn <-
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 1
    )
  expect_error(add_arc_DBN(
    DBN = dbn,
    from = c('A', 't'),
    to = c('B', 't')
  ), NA)
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('A', 't'),
                     to = c('B', 't'))
  expect_equal(class(dbn), 'DBN')
  expect_equal(names(dbn$nodes), c(d_nodes, s_nodes))
  expect_equal(dbn$nodes$A$t$children, 'B_t')
  expect_equal(dbn$nodes$A$t$nbr, 'B_t')
  expect_equal(dbn$nodes$A$t$mb, 'B_t')
  expect_equal(dbn$nodes$B$t$parents, 'A_t')
  expect_equal(dbn$nodes$B$t$nbr, 'A_t')
  expect_equal(dbn$nodes$B$t$mb, 'A_t')
  expect_equal(unname(dbn$arcs[1,]), c('A_t', 'B_t'))
  expect_error(add_arc_DBN(
    DBN = dbn,
    from = c('K', 't_0'),
    to = c('B', 't')
  ),
  "ERROR: Defined node not in DBN!!!")
  expect_error(add_arc_DBN(
    DBN = dbn,
    from = c('A', 't_0'),
    to = c('A', 't_0')
  ),
  "ERROR: The arc defines a loop!!!")
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('C', 't'),
      to = c('A', 't_0')
    ),
    "ERROR: Prior arcs must be between nodes at 't_0'"
  )
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('A', 't'),
      to = c('B', 't-1')
    ),
    "ERROR: Children nodes in DBN must be at time 't' or 't_0'"
  )
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('D', 't'),
      to = c('B', 't')
    ),
    "ERROR: Static Node have not order higher than 0, thus they are not included in transition arcs"
  )
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('A', 't_0'),
      to = c('D', 't_0')
    ),
    "ERROR: A Dynamic Node could not be parent of a Static Node!!!"
  )
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('A', 't'),
      to = c('C', 't+1')
    ),
    "ERROR: Children nodes in DBN must be at time 't' or 't_0'"
  )
  expect_error(add_arc_DBN(
    DBN = dbn,
    from = c('A', 't+1'),
    to = c('C', 't')
  ))#,"ERROR: Parent nodes in DBN must be at time 't_0', 't' or 't-k' (with k <= Markov Order)")
  expect_error(
    add_arc_DBN(
      DBN = dbn,
      from = c('A', 't-4'),
      to = c('C', 't')
    ),
    "ERROR: Transition arcs could not be of order higher than DBN's Markov Order"
  )
})

test_that("remove_arc_DBN returns the expected output", {
  d_nodes <- c('A', 'B', 'C')
  s_nodes <- c('D', 'E')
  dbn <-
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 1
    )
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('A', 't'),
                     to = c('B', 't'))
  expect_error(delete_arc_DBN(
    DBN = dbn,
    from = c('A', 't'),
    to = c('B', 't')
  ), NA)
  dbn <-
    delete_arc_DBN(DBN = dbn,
                   from = c('A', 't'),
                   to = c('B', 't'))
  expect_equal(class(dbn), 'DBN')
  expect_equal(names(dbn$nodes), c(d_nodes, s_nodes))
  expect_equal(dbn$nodes$A$t$children, character(0))
  expect_equal(dbn$nodes$A$t$nbr, character(0))
  expect_equal(dbn$nodes$A$t$mb, character(0))
  expect_equal(dbn$nodes$B$t$parents, character(0))
  expect_equal(dbn$nodes$B$t$nbr, character(0))
  expect_equal(dbn$nodes$B$t$mb, character(0))
  expect_equal(length(unname(dbn$arcs)), 0)
  expect_error(delete_arc_DBN(
    DBN = dbn,
    from = c('A', 't_0'),
    to = c('B', 't_0')
  ),
  "ERROR: Arc do not exists!!!")
})

test_that("reverse_arc_DBN returns the expected output", {
  d_nodes <- c('A', 'B', 'C', 'E')
  s_nodes <- c('D')
  dbn <-
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 1
    )
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('A', 't'),
                     to = c('B', 't'))
  expect_error(reverse_arc_DBN(
    DBN = dbn,
    from = c('A', 't'),
    to = c('B', 't')
  ), NA)
  dbn <-
    reverse_arc_DBN(DBN = dbn,
                    from = c('A', 't'),
                    to = c('B', 't'))
  expect_equal(class(dbn), 'DBN')
  expect_equal(names(dbn$nodes), c(d_nodes, s_nodes))
  expect_equal(dbn$nodes$B$t$children, 'A_t')
  expect_equal(dbn$nodes$B$t$nbr, 'A_t')
  expect_equal(dbn$nodes$B$t$mb, 'A_t')
  expect_equal(dbn$nodes$A$t$parents, 'B_t')
  expect_equal(dbn$nodes$A$t$nbr, 'B_t')
  expect_equal(dbn$nodes$A$t$mb, 'B_t')
  expect_equal(unname(dbn$arcs[1,]), c('B_t', 'A_t'))
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('E', 't-1'),
                     to = c('C', 't'))
  expect_error(
    reverse_arc_DBN(
      DBN = dbn,
      from = c('E', 't-1'),
      to = c('C', 't')
    ),
    "ERROR: Temporal arcs are irreversible in DBNs!!!"
  )
})

test_that("G_0 is generated correctly", {
  d_nodes <- c("A", "B")
  nodes_0 <-
    dbn <- empty_DBN(dynamic_nodes = d_nodes, markov_order = 1)
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('B', 't'),
                to = c('A', 't'))
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('A', 't_0'),
                to = c('B', 't_0'))
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('A', 't-1'),
                to = c('A', 't'))
  dbn <-
    reverse_arc_DBN(DBN = dbn,
                    from = c('B', 't'),
                    to = c('A', 't'))
  expect_error(from_DBN_to_G_0(dbn), NA)
  G_0 <- from_DBN_to_G_0(dbn)
  expect_equal(class(G_0), 'bn')
  expect_equal(sort(names(G_0$nodes)), sort(unlist(lapply(d_nodes, function(x) {
    gsub(' ', '', paste(x, '_0'))
  }))))
  expect_equal(G_0$nodes[['A_0']]$parents, dbn$nodes$A[['t_0']]$parents)
  expect_equal(G_0$nodes[['A_0']]$children, dbn$nodes$A[['t_0']]$children)
  expect_equal(G_0$nodes[['A_0']]$mb, dbn$nodes$A[['t_0']]$mb)
  expect_equal(G_0$nodes[['A_0']]$nbr, dbn$nodes$A[['t_0']]$nbr)
  expect_equal(G_0$nodes[['B_0']]$parents, dbn$nodes$B[['t_0']]$parents)
  expect_equal(G_0$nodes[['B_0']]$children, dbn$nodes$B[['t_0']]$children)
  expect_equal(G_0$nodes[['B_0']]$mb, dbn$nodes$B[['t_0']]$mb)
  expect_equal(G_0$nodes[['B_0']]$nbr, dbn$nodes$B[['t_0']]$nbr)
})

test_that("G_transition is generated correctly", {
  d_nodes <- c("A", "B")
  dbn <- empty_DBN(dynamic_nodes = d_nodes, markov_order = 1)
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('B', 't'),
                to = c('A', 't'))
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('A', 't_0'),
                to = c('B', 't_0'))
  dbn <-
    add_arc_DBN(DBN = dbn,
                from = c('A', 't-1'),
                to = c('A', 't'))
  dbn <-
    reverse_arc_DBN(DBN = dbn,
                    from = c('B', 't'),
                    to = c('A', 't'))
  expect_error(from_DBN_to_G_transition(dbn), NA)
  G_transition <- from_DBN_to_G_transition(dbn)
  expect_equal(class(G_transition), 'bn')
  expect_equal(sort(names(G_transition$nodes)), sort(unlist(lapply(d_nodes, function(x) {
    c(gsub(' ', '', paste(x, '_t')), gsub(' ', '', paste(x, '_t-1')))
  }))))
  expect_equal(G_transition$nodes$A_t$parents, dbn$nodes$A$t$parents)
  expect_equal(G_transition$nodes$A_t$children, dbn$nodes$A$t$children)
  expect_equal(G_transition$nodes$A_t$mb, dbn$nodes$A$t$mb)
  expect_equal(G_transition$nodes$A_t$nbr, dbn$nodes$A$t$nbr)
  expect_equal(G_transition$nodes[['A_t-1']]$parents, dbn$nodes$A[['t-1']]$parents)
  expect_equal(G_transition$nodes[['A_t-1']]$children, dbn$nodes$A[['t-1']]$children)
  expect_equal(G_transition$nodes[['A_t-1']]$mb, dbn$nodes$A[['t-1']]$mb)
  expect_equal(G_transition$nodes[['A_t-1']]$nbr, dbn$nodes$A[['t-1']]$nbr)
  expect_equal(G_transition$nodes$B_t$parents, dbn$nodes$B$t$parents)
  expect_equal(G_transition$nodes$B_t$children, dbn$nodes$B$t$children)
  expect_equal(G_transition$nodes$B_t$mb, dbn$nodes$B$t$mb)
  expect_equal(G_transition$nodes$B_t$nbr, dbn$nodes$B$t$nbr)
  expect_equal(G_transition$nodes[['B_t-1']]$parents, dbn$nodes$B[['t-1']]$parents)
  expect_equal(G_transition$nodes[['B_t-1']]$children, dbn$nodes$B[['t-1']]$children)
  expect_equal(G_transition$nodes[['B_t-1']]$mb, dbn$nodes$B[['t-1']]$mb)
  expect_equal(G_transition$nodes[['B_t-1']]$nbr, dbn$nodes$B[['t-1']]$nbr)
})

test_that("summary is generated correctly", {
  d_nodes <- c('A', 'B', 'C')
  s_nodes <- c('D', 'E')
  added_dynamic <- 'F'
  added_static <- 'G'
  dbn <-
    empty_DBN(
      dynamic_nodes = d_nodes,
      static_nodes = s_nodes,
      markov_order = 2
    )
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('A', 't'),
                     to = c('B', 't'))
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('A', 't-1'),
                     to = c('B', 't'))
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('D', 't_0'),
                     to = c('C', 't_0'))
  dbn <- add_arc_DBN(DBN = dbn,
                     from = c('B', 't'),
                     to = c('C', 't'))
  expect_error(summary_DBN(dbn, test = TRUE), NA)
  summ <- summary_DBN(dbn, test = TRUE)
  expect_equal(class(summ), 'list')
  expect_equal(summ$n_nodes, length(c(d_nodes, s_nodes)))
  expect_equal(summ$n_dynamic_nodes, length(d_nodes))
  expect_equal(summ$n_static_nodes, length(s_nodes))
  expect_equal(summ$n_arcs, 4)
  expect_equal(summ$n_inner_slice_arcs, 2)
  expect_equal(summ$n_intra_slice_arcs, 1)
  expect_equal(summ$n_prior_arcs, 1)
})