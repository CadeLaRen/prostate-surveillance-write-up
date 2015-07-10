
cat("model {

###PRIORS FOR LATENT CLASS MODEL

#flat prior on probability of latent class membership for all subjects
p_eta ~ dbeta(1,1)


###PRIORS FOR MIXED MODEL
#model correlated random effects distribution

for (index in 1:d.Z) {
	xi[index]~dunif(0,100) #scale parameter, same across classes
	for(k in 1:K){
		mu_raw[index,k]~dnorm(0, 0.01) 
		mu[index,k]<-xi[index] *mu_raw[index,k]}  }
		
for(k in 1:K){ #save this iteration of mu with clearer labels. This is not used in model, just to get posteriors more easily.
	mu_int[k] <- mu[1,k]  
	mu_slope[k] <- mu[2,k]} 

#same covariance matrix (Sigma_B_raw) across latent classes
Tau_B_raw ~ dwish(I_d.Z[,], (d.Z+1))  #this is unscaled covariance matrix
Sigma_B_raw[1:d.Z, 1:d.Z] <- inverse(Tau_B_raw[1:d.Z, 1:d.Z])	
for (index in 1:d.Z){
		sigma[index]<-xi[index]*sqrt(Sigma_B_raw[index,index]) } #take into account scaling when saving elements of the covariance matrix 

sigma_int <- sigma[1] # again, this is to track elements of the cov matrix with easier labels
sigma_slope <- sigma[2] 

rho_int_slope <- Sigma_B_raw[1,2]/sqrt(Sigma_B_raw[1,1] * Sigma_B_raw[2,2])
cov_int_slope<- rho_int_slope * sigma_int * sigma_slope * xi[1] * xi[2] #again, take into account scaling (xi) when saving the covariance  


##residual variance, independent of correlated random effects, same across classes
sigma_res ~ dunif(0, 1) 
tau_res <- pow(sigma_res,-2)

##fixed effects
for(index in 1:d.X){
	beta[index] ~ dnorm(0,0.01)}


###PRIORS FOR OUTCOME MODEL
#last element in each gamma is coefficient for class membership eta=1
for(index in 1:(d.W.BX+1)){gamma.BX[index] ~ dnorm(0,0.01)}
for(index in 1:(d.W.RC+1)){gamma.RC[index] ~ dnorm(0,0.01)}
for(index in 1:(d.W.RRP+2)){gamma.RRP[index] ~ dnorm(0,0.01)} #+2 because includes an interaction


###LIKELIHOOD

##latent variable for true cancer state
for(i in 1:n_eta_known){
	eta.data[i] ~ dbern(p_eta)
	eta[i] <- eta.data[i] + 1} #this is for those with path reports from RRP, eta known 

for(i in (n_eta_known+1):n){
	eta.hat[(i-n_eta_known)] ~ dbern(p_eta)
	eta[i] <- eta.hat[(i-n_eta_known)] + 1}  #for those without RRP

##linear mixed effects model for PSA 
#generate random intercept and slope for individual given latent class
for (i in 1:n) {
	B_raw[i,1:d.Z] ~ dmnorm(mu_raw[1:d.Z,eta[i]], Tau_B_raw[1:d.Z, 1:d.Z])
	for(index in 1:d.Z){b.vec[i,index] <- xi[index]*B_raw[i,index]} }

#fit LME
for(j in 1:n_obs_psa){ 
	mu_obs_psa[j] <- inprod(b.vec[subj_psa[j],1:d.Z], Z[j,1:d.Z])  + (beta[1]*X[j,1:d.X]) 
	Y[j] ~ dnorm(mu_obs_psa[j], tau_res) }


##all biopsy data

#logistic regression for biopsy (prob of observation)
for(j in 1:n_bx){
	logit(p_bx[j]) <-inprod(gamma.BX[1:d.W.BX], W.BX[j,1:d.W.BX]) + gamma.BX[(d.W.BX+1)]*equals(eta[subj_bx[j]],2)  
	BX[j] ~ dbern(p_bx[j]) }

#logistic regression for reclassification 	
for(j in 1:n_rc){
	logit(p_rc[j]) <-inprod(gamma.RC[1:d.W.RC], W.RC[j,1:d.W.RC]) + gamma.RC[(d.W.RC+1)]*equals(eta[subj_rc[j]],2) 
	RC[j] ~ dbern(p_rc[j]) }

#logistic regression for surgery 	
for(j in 1:n_rrp){
	logit(p_rrp[j]) <-inprod(gamma.RRP[1:d.W.RRP], W.RRP[j,1:d.W.RRP]) + gamma.RRP[(d.W.RRP+1)]*equals(eta[subj_rrp[j]],2)  + gamma.RRP[(d.W.RRP+2)]*equals(eta[subj_rrp[j]],2)*W.RRP[j,d.W.RRP]  
	RRP[j] ~ dbern(p_rrp[j]) }

	
 }", fill=TRUE, file="model-for-jags-inf-obs.txt")