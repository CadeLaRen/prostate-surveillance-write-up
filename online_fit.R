library(dplyr)

invLogit <- function(x)
return(exp(x)/(1+exp(x)))

# oo<-readRDS('posterior_full_25k.rds')$sims.list
# oo<-oo[!names(oo) %in% c('b.vec','eta.hat')] #don't use previous subj's posteriors.. at least not at this point in the project


##############
# Observed data for subject_star
# star <- 1 #star is defined in the previous script

psa.data.star<-filter(psa.data.full, subj == star)

Y_star <- psa.data.star$log.psa
X_star <- psa.data.star$std.vol
Z_star <- cbind(1, psa.data.star$std.age, psa.data.star$age.basis)

bx.data_star <- filter(bx.data.full, subj==star) 

R_star <- bx.data_star$rc
W.RC_star <-cbind(1,bx.data_star$std.age, bx.data_star$bx.time)
##############





K<-2  # two latent classes
P<-length(oo$p_eta) # number of particles

mu<-array(NA,dim=c(P,3,K)) #indeces: #[particle, ind/slope/spline, k in {1,2}=class]
mu[,1,]<-oo$mu_int
mu[,2,]<-oo$mu_slope
mu[,3,]<-oo$mu_spline

#Get random draws for eta
eta<-rbinom(P,1,prob=oo$p_eta) 
#Get P random draws for b.vec
b.vec.star <- matrix(NA,P,3)
cov_for_bvec <- list()
for(p in 1:P){
	cov_for_bvec[[p]] <- diag(c(oo$sigma_int[p]^2,oo$sigma_slope[p]^2,oo$sigma_spline[p]^2))
	cov_for_bvec[[p]][1,2]<-
	cov_for_bvec[[p]][2,1]<-oo$cov_int_slope[p]
	cov_for_bvec[[p]][1,3]<-
	cov_for_bvec[[p]][3,1]<-oo$cov_int_spline[p]
	cov_for_bvec[[p]][2,3]<-
	cov_for_bvec[[p]][3,2]<-oo$cov_slope_spline[p]
	b.vec.star[p,]<-mvrnorm(1,mu=mu[p,,eta[p]+1], Sigma=cov_for_bvec[[p]])
}

#Y_star is length = n_psa_star
#X_star is dim(n_psa_star,d.X) =  dim(n_psa_star,1)
#Z_star is dim(n_psa_star,d.Z) =  dim(n_psa_star,3)
#R_star is a (n_bx_star) length vector of reclassifications
#W.RC_star<-matrix(runif(d.W.RC_star*n_bx_star),n_bx_star,d.W.RC)

# Likelihood of a particle, given subj_star's data (Y, Z, X, and W),
  # and a proposed set of random effects (b-vec & eta) and 
  # hyper-params (beta, gamma.RC)
W <- rep(NA,P)
for(p in 1:P){
	 #likelihood of PSA data
	mu_obs_psa <- Z_star %*% b.vec.star[p,]  + 
		( oo$beta[p,1] + oo$beta[p,2]*(eta[p]==1) ) * X_star
	L_Y <- prod(dnorm(Y_star,mean=mu_obs_psa, sd=oo$sigma_res[p]))

	#likelihood of reclassifications
	if(!is.null(R_star)){
		logit_p_rc<-cbind(W.RC_star,eta[p]) %*% oo$gamma.RC[p,1:(d.W.RC+1)] 
		L_R <- prod(dbinom(x=R_star,size=1,prob=c(invLogit(logit_p_rc))))
	}else{ 
		L_R <- 1
	}
	W[p] <- L_Y * L_R
}
W <- W/sum(W)
plot(density(eta,weights=W))
# Our eta here is the eta.hat from the main JAGS model.

sum(eta*W)
hist(colMeans(oo$eta.hat),breaks=30,main='Distribution of mean posterior eta_hat,\nacross people')
abline(v=mean(oo$eta.hat),lwd=3)
abline(v=mean(oo$eta.hat[,star+1]),lwd=2, col='black')
abline(v=sum(eta*W),lwd=3,col='blue',lty=2)





