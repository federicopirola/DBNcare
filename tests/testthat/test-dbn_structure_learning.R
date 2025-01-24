DBN <- DynamicBayesianNetwork::generate_dbn_random_structure(c("A","B","C"), TRUE, 0.5, 0.5)
DBN_fitted <- DynamicBayesianNetwork::generate_dbn_nodes_distributions(DBN, TRUE, 2)
data <- DynamicBayesianNetwork::dbn_sampling(DBN_fitted, 2000, 5)
c(data0, dataTR, blacklist0, blacklistTR, whitelist0, whitelistTR) %<-% DynamicBayesianNetwork::StructureLearning.Preprocess(data)
test_that("check DBN Preprocessing Outcome",{
  StructureLearning.Preprocess(data)
})
DBN_learned <- DynamicBayesianNetwork::StructureLearning.Network(data)
DynamicBayesianNetwork::plot_g_transition(DBN_learned)
DynamicBayesianNetwork::plot_g_transition(DBN)
#test_that(DBN_learned$arcs
for (i in c("loglik", "aic", "bic", "ebic", "bde",  "bds", "mbde", "bdla", "k2", "fnml", "qnml", "nal", "pnal")){
  print(paste('score',i))
  DynamicBayesianNetwork::StructureLearning.Network(data, score=i)
}

for (i in c("mi", "mi-adf", "mc-mi", "smc-mi", "sp-mi", "mi-sh",  "x2-adf", "mc-x2", "smc-x2", "sp-x2")){
  print(paste('test',i))
  DynamicBayesianNetwork::StructureLearning.Network(data, algorithm = 'pc.stable', score=i)
}

for (i in c('hc','tabu')){
  print(paste('algorithm',i))
  DynamicBayesianNetwork::StructureLearning.Network(data, algorithm = i)
}

for (i in c('pc.stable','gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')){
  print(paste('algorithm',i))
  DynamicBayesianNetwork::StructureLearning.Network(data, algorithm = i)
}

for (i in c('h2pc','mmhc','rsmax2')){
  print(paste('algorithm',i))
  DBN_learned <- DynamicBayesianNetwork::StructureLearning.Network(data, algorithm = i)
}
