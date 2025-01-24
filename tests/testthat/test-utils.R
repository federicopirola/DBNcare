test_that("get_unrolled_dbn returns the correct output", {
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
  #fitting the dbn
  fitted_my_dbn <- DBN_parameters(my_dbn, CPTs_mydbn)
  dbn_unrolled <- get_unrolled_dbn(fitted_my_dbn, 4)
  all_nodes <- bnlearn::nodes(dbn_unrolled)
  all_arcs <- bnlearn::arcs(dbn_unrolled)
  expected_nodes <-
    c(
      "A_0",
      "B_0",
      "C_0",
      "A_1",
      "B_1",
      "C_1",
      "A_2",
      "B_2",
      "C_2",
      "A_3",
      "B_3",
      "C_3",
      "A_4",
      "B_4",
      "C_4"
    )
  expected_edges <-
    c(
      "A_0",
      "B_0",
      "A_0",
      "C_0",
      "A_0",
      "A_1",
      "B_0",
      "C_0",
      "B_0",
      "B_1",
      "C_0",
      "C_1",
      "A_1",
      "B_1",
      "A_1",
      "C_1",
      "A_1",
      "A_2",
      "B_1",
      "C_1",
      "B_1",
      "B_2",
      "C_1",
      "C_2",
      "A_2",
      "B_2",
      "A_2",
      "C_2",
      "A_2",
      "A_3",
      "B_2",
      "C_2",
      "B_2",
      "B_3",
      "C_2",
      "C_3",
      "A_3",
      "B_3",
      "A_3",
      "C_3",
      "A_3",
      "A_4",
      "B_3",
      "C_3",
      "B_3",
      "B_4",
      "C_3",
      "C_4",
      "A_4",
      "B_4",
      "A_4",
      "C_4",
      "B_4",
      "C_4"
    )
  expected_edges <- matrix(
    expected_edges,
    byrow = TRUE,
    ncol = 2,
    dimnames = list(NULL, c("from", "to"))
  )
  expected_cpts_dbn_unrolled <- list()
  
  
  expect_equal(all_nodes, expected_nodes)
  expect_equal(all_arcs, expected_edges)
  
  #sanity check - CPTs 
  #unfortunately given to changing of order of the cpts in the fitted_my_dbn
  #the check can't be automated with a for loop.
  #checking cpts for timeslice 0
  expect_equal(as.vector(dbn_unrolled$A_0$prob), as.vector(fitted_my_dbn$A_0$prob))
  expect_equal(as.vector(dbn_unrolled$B_0$prob), as.vector(fitted_my_dbn$B_0$prob))
  expect_equal(as.vector(dbn_unrolled$C_0$prob), as.vector(fitted_my_dbn$C_0$prob))
  
  #cheking cpts for A, timeslice 1,2,3,4
  expect_equal(as.vector(dbn_unrolled$A_1$prob), as.vector(fitted_my_dbn$A_t$prob))
  expect_equal(as.vector(dbn_unrolled$A_2$prob), as.vector(fitted_my_dbn$A_t$prob))
  expect_equal(as.vector(dbn_unrolled$A_3$prob), as.vector(fitted_my_dbn$A_t$prob))
  expect_equal(as.vector(dbn_unrolled$A_4$prob), as.vector(fitted_my_dbn$A_t$prob))
  
  #checking cpts for B timeslice 1
  expect_equal(dbn_unrolled$B_1$prob["high","high","yes"],
              fitted_my_dbn$B_t$prob["high","yes", "high"])
  expect_equal(dbn_unrolled$B_1$prob["low","high","yes"],
              fitted_my_dbn$B_t$prob["low","yes","high"])
  expect_equal(dbn_unrolled$B_1$prob["high","low","yes"],
              fitted_my_dbn$B_t$prob["high","yes","low"])
  expect_equal(dbn_unrolled$B_1$prob["low","low","yes"],
              fitted_my_dbn$B_t$prob["low","yes","low"])
  
  #
  expect_equal(dbn_unrolled$B_1$prob["high","high","no"],
              fitted_my_dbn$B_t$prob["high","no", "high"])
  expect_equal(dbn_unrolled$B_1$prob["low","high","no"],
              fitted_my_dbn$B_t$prob["low","no","high"])
  expect_equal(dbn_unrolled$B_1$prob["high","low","no"],
              fitted_my_dbn$B_t$prob["high","no","low"])
  expect_equal(dbn_unrolled$B_1$prob["low","low","no"],
              fitted_my_dbn$B_t$prob["low","no","low"])
  
  
  
  #checking cpts for B timeslice 2
  expect_equal(dbn_unrolled$B_2$prob["high","high","yes"],
               fitted_my_dbn$B_t$prob["high","yes", "high"])
  expect_equal(dbn_unrolled$B_2$prob["low","high","yes"],
               fitted_my_dbn$B_t$prob["low","yes","high"])
  expect_equal(dbn_unrolled$B_2$prob["high","low","yes"],
               fitted_my_dbn$B_t$prob["high","yes","low"])
  expect_equal(dbn_unrolled$B_2$prob["low","low","yes"],
               fitted_my_dbn$B_t$prob["low","yes","low"])
  
  #
  expect_equal(dbn_unrolled$B_2$prob["high","high","no"],
               fitted_my_dbn$B_t$prob["high","no", "high"])
  expect_equal(dbn_unrolled$B_2$prob["low","high","no"],
               fitted_my_dbn$B_t$prob["low","no","high"])
  expect_equal(dbn_unrolled$B_2$prob["high","low","no"],
               fitted_my_dbn$B_t$prob["high","no","low"])
  expect_equal(dbn_unrolled$B_2$prob["low","low","no"],
               fitted_my_dbn$B_t$prob["low","no","low"])
  
  #checking C time slice 1

  #dbn_cpt_slice<-[bn_cpt_slice[1],bn_cpt_slice[3],bn_cpt_slice[4],bn_cpt_slice[2]]
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'medium', 'yes','high'],
              fitted_my_dbn$C_t$prob['medium','yes','high','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'medium', 'yes','high'],
              fitted_my_dbn$C_t$prob['high','yes','high','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'high', 'yes','high'],
              fitted_my_dbn$C_t$prob['medium','yes','high','high'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'high', 'yes','high'],
              fitted_my_dbn$C_t$prob['high','yes','high','high'])
  
  #
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'medium', 'no','high'],
               fitted_my_dbn$C_t$prob['medium','no','high','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'medium', 'no','high'],
               fitted_my_dbn$C_t$prob['high','no','high','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'high', 'no','high'],
               fitted_my_dbn$C_t$prob['medium','no','high','high'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'high', 'no','high'],
               fitted_my_dbn$C_t$prob['high','no','high','high'])
  
  #
  expect_equal(dbn_unrolled$C_1$prob['medium', 'medium', 'yes','low'],
               fitted_my_dbn$C_t$prob['medium','yes','low','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'medium', 'yes','low'],
               fitted_my_dbn$C_t$prob['high','yes','low','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'high', 'yes','low'],
               fitted_my_dbn$C_t$prob['medium','yes','low','high'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'high', 'yes','low'],
               fitted_my_dbn$C_t$prob['high','yes','low','high'])
  
  #
  expect_equal(dbn_unrolled$C_1$prob['medium', 'medium', 'no','low'],
               fitted_my_dbn$C_t$prob['medium','no','low','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'medium', 'no','low'],
               fitted_my_dbn$C_t$prob['high','no','low','medium'])
  
  expect_equal(dbn_unrolled$C_1$prob['medium', 'high', 'no','low'],
               fitted_my_dbn$C_t$prob['medium','no','low','high'])
  
  expect_equal(dbn_unrolled$C_1$prob['high', 'high', 'no','low'],
               fitted_my_dbn$C_t$prob['high','no','low','high'])
  
  
  # checking C_2
  expect_equal(dbn_unrolled$C_2$prob['medium', 'medium', 'yes','high'],
               fitted_my_dbn$C_t$prob['medium','yes','high','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'medium', 'yes','high'],
               fitted_my_dbn$C_t$prob['high','yes','high','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['medium', 'high', 'yes','high'],
               fitted_my_dbn$C_t$prob['medium','yes','high','high'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'high', 'yes','high'],
               fitted_my_dbn$C_t$prob['high','yes','high','high'])
  
  #
  
  expect_equal(dbn_unrolled$C_2$prob['medium', 'medium', 'no','high'],
               fitted_my_dbn$C_t$prob['medium','no','high','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'medium', 'no','high'],
               fitted_my_dbn$C_t$prob['high','no','high','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['medium', 'high', 'no','high'],
               fitted_my_dbn$C_t$prob['medium','no','high','high'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'high', 'no','high'],
               fitted_my_dbn$C_t$prob['high','no','high','high'])
  
  #
  expect_equal(dbn_unrolled$C_2$prob['medium', 'medium', 'yes','low'],
               fitted_my_dbn$C_t$prob['medium','yes','low','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'medium', 'yes','low'],
               fitted_my_dbn$C_t$prob['high','yes','low','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['medium', 'high', 'yes','low'],
               fitted_my_dbn$C_t$prob['medium','yes','low','high'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'high', 'yes','low'],
               fitted_my_dbn$C_t$prob['high','yes','low','high'])
  
  #
  expect_equal(dbn_unrolled$C_2$prob['medium', 'medium', 'no','low'],
               fitted_my_dbn$C_t$prob['medium','no','low','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'medium', 'no','low'],
               fitted_my_dbn$C_t$prob['high','no','low','medium'])
  
  expect_equal(dbn_unrolled$C_2$prob['medium', 'high', 'no','low'],
               fitted_my_dbn$C_t$prob['medium','no','low','high'])
  
  expect_equal(dbn_unrolled$C_2$prob['high', 'high', 'no','low'],
               fitted_my_dbn$C_t$prob['high','no','low','high'])
  
  # testing different parameters for G_0
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
  A_0.prob <- array(c(0.3, 0.7),
                    dim = length(A_lv),
                    dimnames = list(A_0 = A_lv))
  B_0.prob <- array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  C_0.prob <- array(
    c(0.2, 0.8, 0.5, 0.5, 0.5, 0.5, 0.6, 0.4),
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
  #fitting the dbn
  fitted_my_dbn <- DBN_parameters(my_dbn, CPTs_mydbn)
  dbn_unrolled <- get_unrolled_dbn(fitted_my_dbn, 4)
  
  expect_equal(as.vector(dbn_unrolled$A_0$prob), as.vector(fitted_my_dbn$A_0$prob))
  expect_equal(as.vector(dbn_unrolled$B_0$prob), as.vector(fitted_my_dbn$B_0$prob))
  expect_equal(as.vector(dbn_unrolled$C_0$prob), as.vector(fitted_my_dbn$C_0$prob))
})
test_that("get_unrolled_dbn raises error in case of wrong inputs", {
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
  A_0.prob <- array(c(0.3, 0.7),
                    dim = length(A_lv),
                    dimnames = list(A_0 = A_lv))
  B_0.prob <- array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  C_0.prob <- array(
    c(0.2, 0.8, 0.5, 0.5, 0.5, 0.5, 0.6, 0.4),
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
  #fitting the dbn
  fitted_my_dbn <- DBN_parameters(my_dbn, CPTs_mydbn)
  expect_error(get_unrolled_dbn(c(),4))
  expect_error(get_unrolled_dbn(fitted_my_dbn, "3"))
})
test_that(" get_node_edges error in case of wrong inputs", {
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
  A_0.prob <- array(c(0.3, 0.7),
                    dim = length(A_lv),
                    dimnames = list(A_0 = A_lv))
  B_0.prob <- array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  C_0.prob <- array(
    c(0.2, 0.8, 0.5, 0.5, 0.5, 0.5, 0.6, 0.4),
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
  #fitting the dbn
  fitted_my_dbn <- DBN_parameters(my_dbn, CPTs_mydbn)
  g_transition<- from_fitted_DBN_to_fitted_G_transition(fitted_my_dbn)
  expect_error(get_node_edges(c(),"A_t", 3))
  expect_error(get_node_edges(g_transition,3,3))
  expect_error(get_node_edges(g_transition,"A_t","3"))
})
