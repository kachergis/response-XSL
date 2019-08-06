

BIC <- function(SSE, n, k) {
  # Schwarz, 1978 (aka SIC or SBC)
  # n = sample size, k = # of predictors
  # proportional to:
  return(n*log(SSE/n) + log(n)*k)
}


graph_response_by_correct <- function(ld, file='') {
  require(ggplot2)
  ld$Occurrence = factor(ld$Occurrence)
  aesth = aes(color=Training, shape=Training, linetype=Training)
  palette = c("blue1","red3","green3") # uncertain, incorrect, correct
  d <- ggplot(data=ld, aes(x=Occurrence, y=Proportion, group=Training)) + geom_point(aesth) + 
    geom_line(aesth) + facet_grid(. ~ Test) + scale_color_manual(values=palette) + 
    xlab("Training Occurrence") + ylim(0,1.0) + ylab("Response Proportion") + theme_bw()
  if(file!='') {
    ggsave(filename=paste(file,".pdf",sep=''), plot=d, height=4, width=6)
  } else {
    print(d)
  }
}


fit_trues2 <- function(par, ord, human, final_hum, Nsubj=1000) {
  mresps = matrix(0, nrow=Nsubj, ncol=length(human$resps)+1)
  for(s in 1:Nsubj) {
    mp = model(par, ord)
    mresps[s,] = c(colMeans(mp$resps), mp$final)
  }
  mod_means = colMeans(mresps)
  mod_sd = apply(mresps, 2, sd)
  sse = sum((mod_means-c(human$resps,human$final))^2)
  print(c("Par:",par,"SSE:",sse))
  print(mod_means)
  print(mod_sd)
  return(sse)
}

fit_trues <- function(par, ord, human, final_hum, Nsubj=1000) {
  mresps = rep(0, length(human$resps))
  mfin = 0
  for(s in 1:Nsubj) {
    mp = model(par, ord)
    mresps = mresps + colMeans(mp$resps)
    mfin = mfin + mp$final
  }
  mresps = mresps/Nsubj
  mfin = mfin/Nsubj
  sse = sum((mresps-human$resps)^2)
  sse = sse + (mfin-human$final)^2
  print(c("Par:",par,"SSE:",sse))
  print(c(mresps,mfin))
  return(sse)
}


fit_sampling <- function(par, ord, human, K, name="") { # , final_hum -- could try to fit final proportion correct...
  mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)
  sse = sum((mdat$split$Proportion-human$Proportion)^2)
  sse = sse + (mdat$all$Correct[7]-final_hum)^2 
  #print(sse)
  return(sse)
}

fit_sampling_bias_final <- function(par, ord, human, K) { 
  mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)
  sse = sum((mdat$split$Proportion-human$Proportion)^2)
  sse = sse + 18*(mdat$all$Correct[7]-final_hum)^2 # focus mostly on getting 18AFC accuracy
  #print(sse)
  return(sse)
}

loglik <- function(par, ord, human, K, name="") {
  # The discrepancy between predicted and empirical proportions was assessed using the log-likelihood measure, 
  # G2=2N sum_i fi log(fi/mi) ,where N is the number of subjects, fi is the human choice proportion,  
  # mi is the model's predicted proportion, and the index i runs over all 132 data points.
  mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)
  mdat$split$Proportion[which(mdat$split$Proportion==0)] = 1e-9 # for non-inf vals
  finacc = final_hum * log(final_hum / mdat$all$Correct[7])
  Gsq = 2*62 * (sum(human$Proportion * log(human$Proportion / mdat$split$Proportion)) + finacc)
  return(Gsq)
}

oldloglik <- function(par, ord, human, K, name="") {
  # The discrepancy between predicted and empirical proportions was assessed using the log-likelihood measure, 
  # G2=2N sum_i fi log(fi/mi) ,where N is the number of subjects, fi is the human choice proportion,  
  # mi is the model's predicted proportion, and the index i runs over all 132 data points.
  mdat = simulate_subjects(par, ord, K, name, Nsubj=2000)
  mdat$split$Proportion[which(mdat$split$Proportion==0)] = 1e-9 # for non-inf vals
  Gsq = 2*62 * sum(human$Proportion * log(human$Proportion / mdat$split$Proportion))
  return(Gsq)
}

simulate_subjects <- function(par, ord, K, name, Nsubj=1000) {
  dat = data.frame()
  for(s in 1:Nsubj) {
    subjdat = data.frame(model(par, ord, K, name)$resps)
    dat = rbind(dat, subjdat)
  }
  datsum = data.frame(list("Uncertain"=NA,"Incorrect"=NA,"Correct"=NA))
  corsum = data.frame(list("Uncertain"=NA,"Incorrect"=NA,"Correct"=NA))
  incorsum = data.frame(list("Uncertain"=NA,"Incorrect"=NA,"Correct"=NA))
  cor = dat[which(dat[,"Final"]==1),]
  incor = dat[which(dat[,"Final"]==0),] # was -1 (but that's Don't Know)
  for(c in names(dat)) {
    dat[,c] = factor(dat[,c], levels=c("-1","0","1"))
    cor[,c] = factor(cor[,c], levels=c("-1","0","1"))
    incor[,c] = factor(incor[,c], levels=c("-1","0","1"))
    prop = table(dat[,c]) / length(dat[,c])
    prop_cor = table(cor[,c]) / length(cor[,c])
    prop_incor = table(incor[,c]) / length(incor[,c])
    #names(prop)
    datsum = rbind(datsum, prop)
    if(c!="Final") {
      corsum = rbind(corsum, prop_cor)
      incorsum = rbind(incorsum, prop_incor)
    }
  }
  datsum = na.omit(datsum)
  corsum = na.omit(corsum)
  incorsum = na.omit(incorsum)
  datsum$Occurrence = c(paste("Occ.",1:6), "Final")
  corsum$Occurrence = 1:6
  incorsum$Occurrence = 1:6
  incorsum$Test = "Finally Incorrect"
  corsum$Test = "Finally Correct"
  split = rbind(corsum, incorsum)
  library(reshape2)
  ld = melt(split, id=c("Test","Occurrence"))
  ld$Training = ld$variable
  ld$variable = NULL
  ld$Proportion = ld$value
  ld$value = NULL
  return(list(all=datsum, split=ld))
}


get_response_acc_by_final_correctness <- function(par, ord, Nsubj=10000) {
  mresps_cor = matrix(0, nrow=Nsubj, ncol=7)
  mresps_incor = matrix(0, nrow=Nsubj, ncol=7)
  for(s in 1:Nsubj) {
    mp = model(par, ord)
    cor = mp$resps[which(mp$resps[,7]==1),]
    incor = mp$resps[which(mp$resps[,7]==0),]
    if(is.matrix(cor) & nrow(as.matrix(cor))>1) {
      mresps_cor[s,] = colMeans(cor)
    } else if(is.vector(cor)) {
      mresps_cor[s,] = cor
    } 
    if(is.matrix(incor) & nrow(incor)>1) {
      mresps_incor[s,] = colMeans(incor)
    } else if(is.vector(incor)) {
      mresps_incor[s,] = incor
    } 
  }
  mod_means_cor = colMeans(mresps_cor[which(rowSums(mresps_cor)!=0),1:6])
  mod_means_incor = colMeans(mresps_incor[which(rowSums(mresps_incor)!=0),1:6])
  mod_sd_cor = apply(mresps_cor, 2, sd)[1:6]
  mod_sd_incor = apply(mresps_incor, 2, sd)[1:6]
  d = data.frame(Response="Correct",Occurrence=1:6,Proportion=mod_means_cor, SD=mod_sd_cor)
  d = rbind(d, data.frame(Response="Incorrect",Occurrence=1:6,Proportion=mod_means_incor, SD=mod_sd_incor))
  return(d)
}