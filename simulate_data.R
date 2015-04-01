library(MASS)

#Hard coded
p_eta = .4 
d.Z =  3
I_d.Z=diag(d.Z)
K = 2
d.X = 1 # ??
n=896
d.W.RC=3 #W.RC=as.matrix(cbind(rep(1,n_obs_bx),  data.use$age.std, data.use$bx.time))  ??

n_eta_known= 162
n_obs_psa= 10 # ?? (per subj??)

# n_obs_bx = ??
# subj_bx = ??




###PRIORS FOR MIXED MODEL
#model correlated random effects distribution
mu<-matrix(NA,d.Z,K)
mu_int<-
mu_slope<-
mu_spline<-rep(NA,K)
for (index in 1:d.Z) {for(k in 1:K){
		mu[index,k]<-rnorm(1,0, 0.01) }  }
for(k in 1:K){
	mu_int[k] <- mu[1,k] 
	mu_slope[k] <- mu[2,k]
	mu_spline[k] <- mu[3,k] }

#same covariance matrix (Sigma_B) across latent classes
Tau_B <- (d.Z+1)*I_d.Z #mean of rWishart(1, (d.Z+1), I_d.Z)  ??
Sigma_B  <- solve(Tau_B[1:d.Z, 1:d.Z])
sigma <- sqrt(diag(Sigma_B))
sigma_int <- sigma[1] 
sigma_slope <- sigma[2] 
sigma_spline <- sigma[3] 
#the rest is just saving the correlation and elements of cov matrix
rho <- matrix(NA,d.Z,d.Z)
for (row in 1:d.Z) { for (col in 1:d.Z) {
		rho[row,col] <- Sigma_B[row, col]/sqrt(Sigma_B[row, row] * Sigma_B[col, col])} } 
rho_int_slope <- rho[1,2]
rho_int_spline <- rho[1,3]
rho_slope_spline <- rho[2,3]
cov_int_slope<- rho_int_slope*sigma_int*sigma_slope
cov_int_spline<- rho_int_spline*sigma_int*sigma_spline
cov_slope_spline<- rho_slope_spline*sigma_slope*sigma_spline
##residual variance, independent of correlated random effects, same across classes
sigma_res <- 1.5 #mean of runif(0, 3) 
tau_res <- sigma_res^-2
##fixed effects
#stratified by latent class, so 2 coefficients for each predictor in X
beta<- rnorm(2*d.X,0,0.01)
###PRIORS FOR OUTCOME MODEL
gamma.RC <- rnorm((d.W.RC+1),0,0.01)
#last element in gamma is coefficient for class membership eta=1





###LIKELIHOOD

###### latent variable for true cancer state
#this is for those with path reports from RRP, eta known 
eta.data <- c(rbinom(n_eta_known,1,p_eta) ,  
	rep(NA,n-n_eta_known))
eta.hat <- rbinom(n-n_eta_known,1,p_eta)#for those without RRP
eta <- c(eta.data[1:n_eta_known],eta.hat) + 1  #add 1 to use as index later on?
######


##linear mixed effects model for PSA 
#generate random intercept and slope for individual given latent class
b.vec<-matrix(NA,n,d.Z)
for (i in 1:n){
	b.vec[i,] <- mvrnorm(1,mu[1:d.Z,eta[i]], Tau_B)
}


###################################
###################################
# not sure how to proceed from here??

#fit LME
for(j in 1:n_obs_psa){ 
	#subj_psa is some kind of ordering??
	mu_obs_psa[j] <- inprod(b.vec[subj_psa[j],1:d.Z], Z[j,1:d.Z])  + (beta[1]*X[j,1:d.X]) + (beta[2]*X[j,1:d.X])*equals(eta[subj_psa[j]],2) 
	Y[j] ~ dnorm(mu_obs_psa[j], tau_res) }
##logistic regression for reclassification 	
for(j in 1:n_obs_bx){
	logit(p_rc[j]) <-inprod(gamma.RC[1:d.W.RC], W.RC[j,1:d.W.RC]) + gamma.RC[(d.W.RC+1)]*equals(eta[subj_bx[j]],2) 
	R[j] ~ dbern(p_rc[j]) }
	
 }



