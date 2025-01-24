#first step: a random dbn with 5 nodes
library(DynamicBayesianNetwork)
library(dplyr)

dbn_at_random<-random_dbn(LETTERS[1:10],TRUE,0.6)
#second step: fit all nodes distributions generating a multinomial distribution with a cardinality = 2.
dbn_at_random_fitted<-fit_random_dbn(dbn_at_random, TRUE, 2)

dbn_sampling(dbn_at_random_fitted, 1000, 6)

system.time({
  dbn_sampling(dbn_at_random_fitted, 1000, 6)
})
system.time({
  dbn_sampling_old(dbn_at_random_fitted, 1000, 6)
})

set.seed(123)


library(microbenchmark)


microbenchmark(
  original = get_generic_node_name_rex_old("A_t-1"),
  optimized = get_generic_node_name_rex("A_t-1")
)



