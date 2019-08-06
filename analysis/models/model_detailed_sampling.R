# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  gkacherg@indiana.edu  June 10, 2011

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}

update_known <- function(m, tr) {
	startval = .01
		
	for(i in tr) {
		for(c in 1:dim(m)[2]) {
			if(sum(m[,c]>0) & m[i,c]==0) {
				m[i,c] = startval
				m[c,i] = startval
			}
		}
		for(j in tr) {
			if(m[i,j]==0) m[i,j] = startval
			if(m[j,i]==0) m[j,i] = startval
			}
		}
		return(m)
	}


model <- function(params, ord=c(), K=1, name="model", save_traj=FALSE, print_matrix=FALSE) {
  if(K=="max") K = dim(ord)[2] 
  #K = of assocs to update per word
  X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	ent_thresh <- params[4] # should decrease as entropy of an item's associates decreases..let's just multiply by ent
	
  if(length(ord)==0) ord = read.table("orig_4x4.txt", header=F)
	
	voc_sz = max(ord, na.rm=TRUE) # vocabulary size
	mean_ent = c()
	m <- matrix(0, voc_sz, voc_sz) # association matrix
	trial_sz = dim(ord)[2]
  
	# item x occurrence matrix, to be filled in during training 
	resps = matrix(0, voc_sz, 7) # final column is 18AFC correct
	freq = rep(0,voc_sz) # number of occurrences per pair, so far (to index the resps matrix)
  
	perf = c() # mean perf at end of each training ord
	# training
	for(t in 1:dim(ord)[1]) { 
		#print(format(m, digits=3))
		tr = as.integer(ord[t,])
		tr = tr[!is.na(tr)]
		freq[tr] = freq[tr] + 1
		m = update_known(m, tr) # what's been seen so far?
		ent = rep(0, voc_sz)
		for(w in tr) { ent[w] = shannon.entropy(m[w,]) }
		ent = exp(B*ent) 
		
		temp_wts = matrix(0, voc_sz, voc_sz)
		temp_wts[tr,tr] = m[tr,tr] # use these weights to calculate entropy
		if(name=="model") {
		  temp_wts = temp_wts * (ent %*% t(ent)) 
		} else if(name=="unbiased") {
		  temp_wts[tr,tr] = 1
		}
		# if strength-biased, just use tmp_wts
		
		chosen_assocs = matrix(0, voc_sz, voc_sz)
		#print(tr)
		for(w in tr) { 
		  x <- sample(1:voc_sz, K, replace=FALSE, prob=temp_wts[w,]) 
		  chosen_assocs[w,x] = m[w,x] # PK for chosen
		}
		
		if(name=="model") {
		  denom = sum(chosen_assocs * (ent %*% t(ent)))
		  chosen_assocs = (X * chosen_assocs * (ent %*% t(ent))) / denom
		} else if(name=="strength-biased" | name=="unbiased") {
		  denom = sum(chosen_assocs)
		  chosen_assocs = (X * chosen_assocs) / denom
		} 
		  
		m = m*C # decay everything
		m = m + chosen_assocs
		
		if(print_matrix) print(m)
		
		# store knowledge for Nth appearance of each pair on this trial
		mm = m 
		for(s in tr) {
		  probCor = mm[s,s] / sum(mm[s,tr])
		  
		  if(runif(1) <= ent_thresh*shannon.entropy(mm[s,tr])) {
		    resps[s,freq[s]] = -1 # Don't Know
		  } else { # sample from associations
		    response = sample(tr, 1, prob=mm[s,tr])
		    if(response==s) {
		      resps[s,freq[s]] = 1 # Correct
		    } else {
		      resps[s,freq[s]] = 0 # Incorrect
		    }
		  }
		}
    
	  }
	  perf = c(perf, mean(diag(m) / rowSums(m)))
  resp_prob = diag(m) / rowSums(m)
  fin_cor = runif(length(resp_prob))<resp_prob
  fin_cor[which(fin_cor==F)] = 0
  fin_cor[which(fin_cor==T)] = 1
  colnames(resps) = c("Occ.1","Occ.2","Occ.3","Occ.4","Occ.5","Occ.6","Final")
  resps[,7] = fin_cor
	return(list(resps=resps, final=perf))
	}

