
cat("model {

###PRIORS FOR LATENT CLASS MODEL

#flat prior on probability of latent class membership for all subjects
p_eta ~ dbeta(1,1)


###PRIORS FOR MIXED MODEL
#model correlated random effects distribution

for (index in 1:d.Z) {for(k in 1:K){
		mu[index,k]~dnorm(0, 0.01) }  }

for(k in 1:K){
	mu_int[k] <- mu[1,k] 
	mu_slope[k] <- mu[2,k]
	mu_spline[k] <- mu[3,k] }

#same covariance matrix (Sigma_B) across latent classes
Tau_B ~ dwish(I_d.Z[,], (d.Z+1)) 
Sigma_B[1:d.Z, 1:d.Z] <- inverse(Tau_B[1:d.Z, 1:d.Z])
for (index in 1:d.Z) {
	sigma[index] <- sqrt(Sigma_B[index,index])}

sigma_int <- sigma[1] 
sigma_slope <- sigma[2] 
sigma_spline <- sigma[3] 

#the rest is just saving the correlation and elements of cov matrix
for (row in 1:d.Z) { for (col in 1:d.Z) {
		rho[row,col] <- Sigma_B[row, col]/sqrt(Sigma_B[row, row] * Sigma_B[col, col])} } 

rho_int_slope <- rho[1,2]
rho_int_spline <- rho[1,3]
rho_slope_spline <- rho[2,3]

cov_int_slope<- rho_int_slope*sigma_int*sigma_slope
cov_int_spline<- rho_int_spline*sigma_int*sigma_spline
cov_slope_spline<- rho_slope_spline*sigma_slope*sigma_spline


##residual variance, independent of correlated random effects, same across classes
sigma_res ~ dunif(0, 3) 
tau_res <- pow(sigma_res,-2)

##fixed effects
#stratified by latent class, so 2 coefficients for each predictor in X
for(index in 1:(2*d.X)){
	beta[index] ~ dnorm(0,0.01)}


###PRIORS FOR OUTCOME MODEL
for(index in 1:(d.W.RC+1)){gamma.RC[index] ~ dnorm(0,0.01)}
#last element in gamma is coefficient for class membership eta=1

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
	b.vec[i,1:d.Z] ~ dmnorm(mu[1:d.Z,eta[i]], Tau_B[1:d.Z,1:d.Z])}

#fit LME
for(j in 1:n_obs_psa){ 
	mu_obs_psa[j] <- inprod(b.vec[subj_psa[j],1:d.Z], Z[j,1:d.Z])  + (beta[1]*X[j,1:d.X]) + (beta[2]*X[j,1:d.X])*equals(eta[subj_psa[j]],2) 
	Y[j] ~ dnorm(mu_obs_psa[j], tau_res) }


##logistic regression for reclassification 	
for(j in 1:n_obs_bx){
	logit(p_rc[j]) <-inprod(gamma.RC[1:d.W.RC], W.RC[j,1:d.W.RC]) + gamma.RC[(d.W.RC+1)]*equals(eta[subj_bx[j]],2) 
	R[j] ~ dbern(p_rc[j]) }
	
 }", fill=TRUE, file="prediction-model.txt")