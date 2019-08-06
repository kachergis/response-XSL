source("fitting_methods.R")
require(DEoptim) # require(pso)

final_hum = .3449 # 18AFC accuracy for humans at final test

#load("cogsci2014_hypothesis_expanded.RData")
#load("cogsci2014_hypothesis.RData")
#fin18afc = aggregate(Correct ~ Word, mean, data=subset(test, Cond=="hypothesis")) 
#final18afc_acc = fin18afc$Correct

load("training_responses_by_final_acc.RData") # gg - avg resp trajectories for learned vs. unlearned items
ord = read.table("orig_4x4.txt", header=F)

#graph_response_by_correct(gg, file="human_training_responses_by_final_acc")

fit_and_graph_model <- function(modelname, K=NA, lower=c(.00001,.01,.5,0), upper=c(10,30,1,1)) {
  source(paste(modelname,".R",sep=''))
  #bestk1 = DEoptim(fn=fit_sampling, ord=ord, human=gg, K=1, name="model", lower=c(.00001,.01,.5,0), upper=c(10,30,1,1), DEoptim.control(reltol=.001, steptol=50, itermax=200, trace=10))
  best = DEoptim(fn=loglik, ord=ord, human=gg, K=K, name="model", lower=lower, upper=upper, 
                 DEoptim.control(reltol=.001, steptol=60, itermax=200, trace=10))
  print(best)
  par = best$optim$bestmem # X=learn rate, B=unc/PK bias, C=memory decay, ent_thresh=when to say "Don't Know"
  print(par)
  mod = simulate_subjects(par, ord, name="model", K=K, Nsubj=10000)
  save(best, mod, file=paste(modelname,"_k",K,"_best_fit.RData",sep=''))
  graph_response_by_correct(mod$split, file=paste(modelname,"_k",K,"_training_responses_by_final_acc",sep=''))
}

#fit_and_graph_model("model_detailed_sampling_str", K=1)

fit_and_graph_model("model_detailed_sampling_str", K=3)

#fit_and_graph_model("model_detailed_sampling", K=2)
# loglik = 32.23 par=c(0.17420, 0.01820, 0.99183, 0.13067) 

#fit_and_graph_model("model_detailed_sampling", K=3)
# loglik = 32.93 par=c(0.13703, 0.01363, 0.97935, 0.12406) 
# BIC = log(54)*4 + 2*32.93 = 81.82

#fit_and_graph_model("model_detailed_sampling", K=1)
# loglik = 31.60 par = (0.1172661, 0.1800034, 0.9359698, 0.1521185) # final acc: .37
# BIC = log(54)*4 + 2*31.60 = 79.16

#fit_and_graph_model("model_detailed_sampling", K=4)
# loglik = 33.17 par=(0.1317579, 21.1838329, 0.9980332, 0.1146410) # final acc: .33
# BIC = log(54)*4 + 2*33.17 = 82.30

#fit_and_graph_model("trueswell2012_model_detailed_uncertain", K=NA, lower=c(0,0,0), upper=c(1,1,1))
# loglik = 37.52 tru_par = c(0.14823280, 0.68238302, 0.05281296) 18AFCacc=.41 final item-level acc cor: .42
# BIC = log(54)*3 + 2*37.52 = 87.01

#fit_and_graph_model("hypoth_model_detailed", K=NA, lower=c(0,0,0), upper=c(1,1,1))
# loglik = 37.16 hyppar = c(0.1748480, 0.3726927, 0.3738978) # 18AFCacc = .36  item-level final acc cor: -.08
# BIC = ln(n)*k - 2*ln(L), k = # free params, n = # data points
# L = maximized value of likelihood function
# BIC = log(54)*3 + 2*37.16 = 86.29

#mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)

#load("model_detailed_sampling_k3_best_fit.RData")

#load("hypoth_model_detailed_kNA_best_fit.RData")

#load("trueswell2012_model_detailed_uncertain_kNA_best_fit.RData")

# fit_sampling(park1, ord, gg, 2) # .42
# loglik(par, ord, gg, 1) # 92.24

# BIC = -2*logLik + k*ln(n), k=#pars, n=number of observations in the data set
