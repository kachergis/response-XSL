source("fitting_methods.R")
require(DEoptim) # require(pso)

final_hum = .3449 # 18AFC accuracy for humans at final test

load("training_responses_by_final_acc.RData") # gg - avg resp trajectories for learned vs. unlearned items
gg$Training = factor(gg$Training, levels=c("Uncertain","Incorrect","Correct"))
gg = gg[order(gg$Training,gg$Test,gg$Occurrence),] # to match the order of the model data.frame
ord = read.table("orig_4x4.txt", header=F)


#graph_response_by_correct(gg, file="human_training_responses_by_final_acc")

fit_and_graph_model <- function(modelname, K=NA, lower=c(.00001,.01,.5,0), upper=c(10,30,1,1)) {
  source(paste(modelname,".R",sep=''))
  #bestk1 = DEoptim(fn=fit_sampling, ord=ord, human=gg, K=1, name="model", lower=c(.00001,.01,.5,0), upper=c(10,30,1,1), DEoptim.control(reltol=.001, steptol=50, itermax=200, trace=10))
  best = DEoptim(fn=loglik, ord=ord, human=gg, K=K, name="model", lower=lower, upper=upper, 
                 DEoptim.control(reltol=.001, steptol=50, itermax=200, trace=10))
  print(best)
  par = best$optim$bestmem # X=learn rate, B=unc/PK bias, C=memory decay, ent_thresh=when to say "Don't Know"
  print(par)
  mod = simulate_subjects(par, ord, name="model", K=K, Nsubj=10000)
  save(best, mod, file=paste(modelname,"_k",K,"_best_fit.RData",sep=''))
  graph_response_by_correct(mod$split, file=paste(modelname,"_k",K,"_training_responses_by_final_acc",sep=''))
}


fit_and_graph_model("model_detailed_sampling", K=1)
# G^2 = 19.96 par = c(0.125119, 0.123456, 0.927544, 0.149003) # final acc: .37 SSE=.42
loglik(par, ord, gg, K=1, name=name) # 102.65
oldloglik(par, ord, gg, K=1, name=name) # 85.09
# G^2 = 1.06 par = c(.174647, 4.711791, .849637, .202020) 18AFCacc=.71 (not done fitting yet)
loglik(par, ord, gg, K=1, name=name) # 158.5
oldloglik(par, ord, gg, K=1, name=name) # 136.4
# BIC = 19.96 / -62 + 4*log(37) = 14.12


#fit_and_graph_model("model_detailed_sampling", K=2)
# G^2 = 36.347072 par=c(0.160752, 0.012067, 0.847887, 0.156473) final acc: .69

#fit_and_graph_model("model_detailed_sampling", K=3)
# G^2 = 34.334731 par=c(0.032850, 8.321289, 0.555843, 0.196116) final acc: .84

fit_and_graph_model("model_detailed_sampling", K=4)
# G^2 = 38.02 par = c(0.03344791, 6.62286674, 0.52243956, 0.18556479) # final acc: .88  SSE=1.134193
loglik(par, ord, gg, K=4, name=name) # 151.98
oldloglik(par, ord, gg, K=4, name=name) # 184.46
# G^2 = 11.00 par = c(.033499, 7.215464, .520277, .190457) 18AFCacc=.86 (not done fitting yet)

fit_and_graph_model("trueswell2012_model_detailed_uncertain", K=NA, lower=c(0,0,0), upper=c(1,1,1))
# G^2 = 605.04 tru_par = c(0.4148646, 0.8686193, 0.1232299)  18AFCacc=.73  SSE=.59 with G^2 loss function
# G^2 = 560.05 tru_par = c(0.5028206, 0.9717068, 0.1430175)  18AFCacc=.81
# BIC = 560.05 / -62 + 3*log(37) = 1.80

fit_and_graph_model("hypoth_model_detailed", K=NA, lower=c(0,0,0), upper=c(1,1,1))
# G^2 = 338.32 hyppar = c(0.095662, 0.005485, 0.257254) # 18afc acc=.01
# G^2 = 329.51 hyppar = c(0.07598824, 0.49927647, 0.45697262) 18AFCacc=
# BIC = 329.51 / -62 + 3*log(37) = 5.52

#mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)

#load("model_detailed_sampling_k3_best_fit.RData")

#load("hypoth_model_detailed_kNA_best_fit.RData")

#load("trueswell2012_model_detailed_uncertain_kNA_best_fit.RData")

# fit_sampling(park1, ord, gg, 2) # .42
# loglik(par, ord, gg, 1) # 92.24

# BIC = -2*logLik + k*ln(n), k=#pars, n=number of observations in the data set
