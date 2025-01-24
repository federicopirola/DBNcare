test_that("CPTs are generated correctly", {
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
  G_0 <- from_DBN_to_G_0(dbn)
  
  G_transition <- from_DBN_to_G_transition(dbn)
  
  A_lv = c('yes', 'no')
  B_lv = c('high', 'low')
  A_0.prob = array(c(0.2, 0.8),
                   dim = length(A_lv),
                   dimnames = list(A_0 = A_lv))
  B_0.prob = array(
    c(0.25, 0.75, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  dims <- list(A_t = A_lv)
  dims[['A_t-1']] = A_lv
  A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                   dim = c(length(A_lv), length(A_lv)),
                   dimnames = dims)
  B_t.prob = array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_t = B_lv, A_t = A_lv)
  )
  
  CPTs_toy = list(
    A_0 = A_0.prob,
    B_0 = B_0.prob,
    A_t = A_t.prob,
    B_t = B_t.prob
  )
  
  fitted_DBN <-
    DynamicBayesianNetwork::DBN_parameters(DBN = dbn, CPTs = CPTs_toy)
  expect_equal(class(fitted_DBN), 'dbn.fit')
  expect_equal(fitted_DBN$A_0$node, 'A_0')
  expect_equal(fitted_DBN$B_0$node, 'B_0')
  expect_equal(fitted_DBN$A_t$node, 'A_t')
  expect_equal(fitted_DBN$B_t$node, 'B_t')
  expect_equal(fitted_DBN$A_0$parents, dbn$nodes$A$t_0$parents)
  expect_equal(fitted_DBN$B_0$parents, dbn$nodes$B$t_0$parents)
  expect_equal(fitted_DBN$A_t$parents, dbn$nodes$A$t$parents)
  expect_equal(fitted_DBN$B_t$parents, dbn$nodes$B$t$parents)
  expect_equal(fitted_DBN$A_0$children, dbn$nodes$A$t_0$children)
  expect_equal(fitted_DBN$B_0$children, dbn$nodes$B$t_0$children)
  expect_equal(fitted_DBN$A_t$children, dbn$nodes$A$t$children)
  expect_equal(fitted_DBN$B_t$children, dbn$nodes$B$t$children)
  expect_equal(fitted_DBN$A_0$prob, CPTs_toy$A_0)
  expect_equal(fitted_DBN$B_0$prob, CPTs_toy$B_0)
  expect_equal(fitted_DBN$A_t$prob, CPTs_toy$A_t)
  expect_equal(fitted_DBN$B_t$prob, CPTs_toy$B_t)
  expect_equal(sum(fitted_DBN$A_0$prob), 1)
  expect_equal(all(apply(
    fitted_DBN$B_0$prob, setdiff(1:length(dimnames(CPTs_toy$B_0)), which(names(
      dimnames(CPTs_toy$B_0)
    ) == 'B_0')), sum
  ) == 1), TRUE)
  expect_equal(all(apply(
    fitted_DBN$A_t$prob, setdiff(1:length(dimnames(CPTs_toy$A_t)), which(names(
      dimnames(CPTs_toy$A_t)
    ) == 'A_t')), sum
  ) == 1), TRUE)
  expect_equal(all(apply(
    fitted_DBN$B_t$prob, setdiff(1:length(dimnames(CPTs_toy$B_t)), which(names(
      dimnames(CPTs_toy$B_t)
    ) == 'B_t')), sum
  ) == 1), TRUE)
  
})

test_that("Learned G_0 has been generated correctly", {
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
  G_0 <- from_DBN_to_G_0(dbn)
  
  G_transition <- from_DBN_to_G_transition(dbn)
  
  A_lv = c('yes', 'no')
  B_lv = c('high', 'low')
  A_0.prob = array(c(0.2, 0.8),
                   dim = length(A_lv),
                   dimnames = list(A_0 = A_lv))
  B_0.prob = array(
    c(0.25, 0.75, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  dims <- list(A_t = A_lv)
  dims[['A_t-1']] = A_lv
  A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                   dim = c(length(A_lv), length(A_lv)),
                   dimnames = dims)
  B_t.prob = array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_t = B_lv, A_t = A_lv)
  )
  
  CPTs_toy = list(
    A_0 = A_0.prob,
    B_0 = B_0.prob,
    A_t = A_t.prob,
    B_t = B_t.prob
  )
  
  fitted_DBN <-
    DynamicBayesianNetwork::DBN_parameters(DBN = dbn, CPTs = CPTs_toy)
  fitted_G_0 <-
    DynamicBayesianNetwork::from_fitted_DBN_to_fitted_G_0(fitted_DBN)
  expect_equal(class(fitted_G_0), 'bn.fit')
  expect_equal(names(fitted_G_0), c('A_0', 'B_0'))
})

test_that("Learned Transition Network has been generated correctly", {
  d_nodes <- c("A", "B")
  dbn <-
    DynamicBayesianNetwork::empty_DBN(dynamic_nodes = d_nodes, markov_order = 1)
  dbn <-
    DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                        from = c('B', 't'),
                                        to = c('A', 't'))
  dbn <-
    DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                        from = c('A', 't_0'),
                                        to = c('B', 't_0'))
  dbn <-
    DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                        from = c('A', 't-1'),
                                        to = c('A', 't'))
  dbn <-
    DynamicBayesianNetwork::reverse_arc_DBN(DBN = dbn,
                                            from = c('B', 't'),
                                            to = c('A', 't'))
  G_0 <- DynamicBayesianNetwork::from_DBN_to_G_0(dbn)
  
  G_transition <-
    DynamicBayesianNetwork::from_DBN_to_G_transition(dbn)
  
  A_lv = c('yes', 'no')
  B_lv = c('high', 'low')
  A_0.prob = array(c(0.2, 0.8),
                   dim = length(A_lv),
                   dimnames = list(A_0 = A_lv))
  B_0.prob = array(
    c(0.25, 0.75, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_0 = B_lv, A_0 = A_lv)
  )
  dims <- list(A_t = A_lv)
  dims[['A_t-1']] = A_lv
  A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                   dim = c(length(A_lv), length(A_lv)),
                   dimnames = dims)
  B_t.prob = array(
    c(0.2, 0.8, 0.6, 0.4),
    dim = c(length(B_lv), length(A_lv)),
    dimnames = list(B_t = B_lv, A_t = A_lv)
  )
  
  CPTs_toy = list(
    A_0 = A_0.prob,
    B_0 = B_0.prob,
    A_t = A_t.prob,
    B_t = B_t.prob
  )
  
  fitted_DBN <-
    DynamicBayesianNetwork::DBN_parameters(DBN = dbn, CPTs = CPTs_toy)
  fitted_G_transition <-
    DynamicBayesianNetwork::from_fitted_DBN_to_fitted_G_transition(fitted_DBN)
  expect_equal(class(fitted_G_transition), 'bn.fit')
  expect_equal(names(fitted_G_transition), c('A_t', 'B_t'))
})

test_that("CPTs learned from randomly generate data are almost equal to original CPTs",
          {
            d_nodes <- c("A", "B")
            dbn <-
              DynamicBayesianNetwork::empty_DBN(dynamic_nodes = d_nodes, markov_order = 1)
            dbn <-
              DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                                  from = c('A', 't'),
                                                  to = c('B', 't'))
            dbn <-
              DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                                  from = c('A', 't_0'),
                                                  to = c('B', 't_0'))
            dbn <-
              DynamicBayesianNetwork::add_arc_DBN(DBN = dbn,
                                                  from = c('A', 't-1'),
                                                  to = c('A', 't'))
            G_0 <- DynamicBayesianNetwork::from_DBN_to_G_0(dbn)
            
            G_transition <-
              DynamicBayesianNetwork::from_DBN_to_G_transition(dbn)
            
            #B_lv = c('yes','no')
            A_lv = c('yes', 'no')
            B_lv = c('high', 'low')
            A_0.prob = array(c(0.2, 0.8),
                             dim = length(A_lv),
                             dimnames = list(A_0 = A_lv))
            B_0.prob = array(
              c(0.25, 0.75, 0.6, 0.4),
              dim = c(length(B_lv), length(A_lv)),
              dimnames = list(B_0 = B_lv, A_0 = A_lv)
            )
            dims <- list(A_t = A_lv)
            dims[['A_t-1']] = A_lv
            A_t.prob = array(c(0.8, 0.2, 0.05, 0.95),
                             dim = c(length(A_lv), length(A_lv)),
                             dimnames = dims)
            B_t.prob = array(
              c(0.2, 0.8, 0.6, 0.4),
              dim = c(length(B_lv), length(A_lv)),
              dimnames = list(B_t = B_lv, A_t = A_lv)
            )
            
            CPTs_toy = list(
              A_0 = A_0.prob,
              B_0 = B_0.prob,
              A_t = A_t.prob,
              B_t = B_t.prob
            )
            
            fitted_DBN <-
              DynamicBayesianNetwork::DBN_parameters(DBN = dbn, CPTs = CPTs_toy)
            N_samples <- 5000
            Time <- 4
            sampling_set <-
              DynamicBayesianNetwork::dbn_sampling(fitted_DBN, N_samples, Time)
            learned_dbn <-
              DynamicBayesianNetwork::DBN_parameters(DBN = dbn, data = sampling_set)
            expect_equal(sum(learned_dbn$A_0$prob), 1)
            expect_equal(all(apply(
              learned_dbn$B_0$prob, setdiff(1:length(dimnames(CPTs_toy$B_0)), which(names(
                dimnames(CPTs_toy$B_0)
              ) == 'B_0')), sum
            ) == 1), TRUE)
            expect_equal(all(apply(
              learned_dbn$A_t$prob, setdiff(1:length(dimnames(CPTs_toy$A_t)), which(names(
                dimnames(CPTs_toy$A_t)
              ) == 'A_t')), sum
            ) == 1), TRUE)
            expect_equal(all(apply(
              learned_dbn$B_t$prob, setdiff(1:length(dimnames(CPTs_toy$B_t)), which(names(
                dimnames(CPTs_toy$B_t)
              ) == 'B_t')), sum
            ) == 1), TRUE)
            expect_equal(array(learned_dbn$A_0$prob),
                         array(fitted_DBN$A_0$prob[order(A_lv)]),
                         tolerance = 0.05)
            expect_equal(array(learned_dbn$A_t$prob),
                         array(fitted_DBN$A_t$prob[order(A_lv), order(A_lv)]),
                         tolerance = 0.05)
            expect_equal(array(learned_dbn$B_0$prob),
                         array(fitted_DBN$B_0$prob[order(B_lv), order(A_lv)]),
                         tolerance = 0.05)
            expect_equal(array(learned_dbn$B_t$prob),
                         array(fitted_DBN$B_t$prob[order(B_lv), order(A_lv)]),
                         tolerance = 0.05)
            d_nodes_2 <- c('H_2_0', 'Infected', 'Tumor_size')
            dbn_2 <-
              DynamicBayesianNetwork::empty_DBN(dynamic_nodes = d_nodes_2, markov_order = 1)
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('H_2_0', 't'),
                to = c('Infected', 't')
              )
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('Infected', 't-1'),
                to = c('Infected', 't')
              )
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('Tumor_size', 't-1'),
                to = c('Tumor_size', 't')
              )
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('Infected', 't-1'),
                to = c('Tumor_size', 't')
              )
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('H_2_0', 't_0'),
                to = c('Infected', 't_0')
              )
            dbn_2 <-
              DynamicBayesianNetwork::add_arc_DBN(
                DBN = dbn_2,
                from = c('Infected', 't_0'),
                to = c('Tumor_size', 't_0')
              )
            Infected_lv <- c('yes', 'no')
            H_2_0_lv <- c('high', 'low')
            Tumor_size_lv <- c('big', 'small')
            H_2_0_0.prob = array(c(0.2, 0.8),
                                 dim = length(H_2_0_lv),
                                 dimnames = list(H_2_0_0 = H_2_0_lv))
            
            Infected_0.prob = array(
              c(0.22, 0.78, 0.6, 0.4),
              dim = c(length(Infected_lv), length(H_2_0_lv)),
              dimnames = list(Infected_0 = Infected_lv, H_2_0_0 = H_2_0_lv)
            )
            
            Tumor_size_0.prob = array(
              c(0.9, 0.1, 0.05, 0.95),
              dim = c(length(Tumor_size_lv), length(Infected_lv)),
              dimnames = list(Tumor_size_0 = Tumor_size_lv, Infected_0 = Infected_lv)
            )
            H_2_0_t.prob = array(c(0.2, 0.8),
                                 dim = length(H_2_0_lv),
                                 dimnames = list(H_2_0_t = H_2_0_lv))
            
            Infected_t.prob = array(
              c(0.80, 0.20, 0.98, 0.02, 0.03, 0.97, 0.10, 0.90),
              dim = c(length(Infected_lv), length(H_2_0_lv), length(Infected_lv)),
              dimnames = list(
                Infected_t = Infected_lv,
                H_2_0_t = H_2_0_lv,
                `Infected_t-1` = Infected_lv
              )
            )
            
            Tumor_size_t.prob = array(
              c(0.83, 0.17, 0.45, 0.55, 0.80, 0.20, 0.31, 0.69),
              dim = c(
                length(Tumor_size_lv),
                length(Tumor_size_lv),
                length(Infected_lv)
              ),
              dimnames = list(
                Tumor_size_t = Tumor_size_lv,
                `Tumor_size_t-1` = Tumor_size_lv,
                `Infected_t-1` = Infected_lv
              )
            )
            cpts_2 = list(
              H_2_0_0 = H_2_0_0.prob,
              Infected_0 = Infected_0.prob,
              Tumor_size_0 = Tumor_size_0.prob,
              H_2_0_t = H_2_0_t.prob,
              Infected_t = Infected_t.prob,
              Tumor_size_t = Tumor_size_t.prob
            )
            dbn_from_cpts <-
              DynamicBayesianNetwork::DBN_parameters(DBN = dbn_2, CPTs = cpts_2)
            
            sampling_set_2 <-
              DynamicBayesianNetwork::dbn_sampling(dbn_from_cpts, N_samples, Time)
            
            dbn_from_random_data <-
              DynamicBayesianNetwork::DBN_parameters(DBN = dbn_2, data = sampling_set_2)
            
            expect_equal(
              array(dbn_from_random_data$H_2_0_0$prob),
              array(dbn_from_cpts$H_2_0_0$prob[order(H_2_0_lv)]),
              tolerance = 0.05
            )
            expect_equal(
              array(dbn_from_random_data$Infected_0$prob),
              array(dbn_from_cpts$Infected_0$prob[order(Infected_lv), order(H_2_0_lv)]),
              tolerance = 0.05
            )
            expect_equal(
              array(dbn_from_random_data$Tumor_size_0$prob),
              array(dbn_from_cpts$Tumor_size_0$prob[order(Tumor_size_lv), order(Infected_lv)]),
              tolerance = 0.05
            )
            expect_equal(
              array(dbn_from_random_data$H_2_0_t$prob),
              array(dbn_from_cpts$H_2_0_t$prob[order(H_2_0_lv)]),
              tolerance = 0.05
            )
            expect_equal(
              array(dbn_from_random_data$Infected_t$prob),
              array(dbn_from_cpts$Infected_t$prob[order(Infected_lv), order(H_2_0_lv), order(Infected_lv)]),
              tolerance = 0.05
            )
            expect_equal(
              array(dbn_from_random_data$Tumor_size_t$prob),
              array(dbn_from_cpts$Tumor_size_t$prob[order(Tumor_size_lv), order(Tumor_size_lv), order(Infected_lv)]),
              tolerance = 0.05
            )
          })