#setwd("/Users/aaronfisher/Dropbox/Future Projects/inHealth Prostate Screening/repo")

library(dplyr)
library(MASS)

invLogit <- function(x)
return(exp(x)/(1+exp(x)))






########
#get full data
psa.data.full<-read.csv("simulation-data/psa-data-sim.csv")
pt.ordered.full<-read.csv("simulation-data/pt-ordered-sim.csv") 
bx.data.full <- read.csv("simulation-data/bx-data-sim.csv")
eta.data.full<-read.csv("simulation-data/eta-data-sim.csv") #

#Output from leave-one-Out JAGS model  (output-out -> oo)
	# For now just use the output from the full model as an approximate replacement
oo <- readRDS('posterior_full_100k.rds')$sims.list 



################
#Posterior estimates from full JAGS model
of<-readRDS('posterior_full_100k.rds')$sims.list
known_etas<-eta.data.full[!is.na(eta.data.full[,2]),2]
eta_true_or_jags<-c(known_etas,colMeans(of$eta.hat))

# Expand these posterior draws by repeating them nreps times.
	# Reduce monte-carlo error.
# This gives an unweighted particle distribution for latent variables for next subject.
# We'll then reweight these particles individually, for each subject, based on that subject's likelihood.

K<-2  # two latent classes
P<-length(oo$p_eta) # number of particles from JAGS output
nreps<-20 # extend all arrays by repeating them nreps times
P*nreps # total number of particles we'll be weighting over

#Assign by cycle over nreps times (two vectors of different lengths)
mu<-array(NA,dim=c(nreps*P,3,K))
	#indeces: #[particle, int/slope/spline, k in {1,2}=class]
mu[,1,1]<-oo$mu_int[,1]
mu[,1,2]<-oo$mu_int[,2]
mu[,2,1]<-oo$mu_slope[,1]
mu[,2,2]<-oo$mu_slope[,2]
mu[,3,1]<-oo$mu_spline[,1]
mu[,3,2]<-oo$mu_spline[,2]
all(mu[1,,]==mu[1+P,,]) #if TRUE then cycling worked!

#Get random draws for eta
eta<-rbinom(P*nreps,1,prob=rep(c(oo$p_eta),times=nreps))
	# Our eta here is analogous to eta.hat from the main JAGS model.
#Get P random draws for b.vec
b.vec.star <- matrix(NA,P*nreps,3)
cov_for_bvec <- list() #Also consider *not* storing this, since it's not used after b.vec.star is generated.
(expand_time<-system.time({
pb_sim <- txtProgressBar(min = 0, max = P*nreps, char = "=", style=3)
for(r in 1:nreps){ # nreps = number of times we expand the particles from oo.
for(oo_ind in 1:P){ #index for oo
	p <- (r-1)*P+oo_ind #index along the extended particles we're creating.
	setTxtProgressBar(pb_sim,p)

	cov_for_bvec[[p]] <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2,oo$sigma_spline[oo_ind]^2))
	cov_for_bvec[[p]][1,2]<-
	cov_for_bvec[[p]][2,1]<-oo$cov_int_slope[oo_ind]
	cov_for_bvec[[p]][1,3]<-
	cov_for_bvec[[p]][3,1]<-oo$cov_int_spline[oo_ind]
	cov_for_bvec[[p]][2,3]<-
	cov_for_bvec[[p]][3,2]<-oo$cov_slope_spline[oo_ind]
	b.vec.star[p,]<-mvrnorm(1,mu=mu[p,,eta[p]+1], Sigma=cov_for_bvec[[p]])
		#note, eta[i] is *not* nessecarily the same as eta[i+P] due to random draws.
}}})) #~ 2 min for nreps=50, P=25000
################



################
#Get weights for each subject, save weighted etas.

n_subj_to_est <- dim(oo$eta.hat)[2]
etas_weighted <- rep(NA,1000) #will contain posterior mean estimates (eta_hats) from online fits. 
# One for each subject. Some entries will be NAs at the end of this script, 
# to maintain indexing system of: subj_star's eta = etas_weighted[star]

#To test for star = some random number, use
	#star <- round(runif(1,1000-n_subj_to_est+1,1000)) 


pb <- txtProgressBar(min = 0, max = n_subj_to_est, char = "=", style=3)
for(star in (1000-n_subj_to_est+1):1000){
	setTxtProgressBar(pb, star-(1000-n_subj_to_est+1))

	########
	#data for subj_star
	psa.data.star<-filter(psa.data.full, subj == star)

	Y_star <- psa.data.star$log.psa
	X_star <- psa.data.star$std.vol
	Z_star <- cbind(1, psa.data.star$std.age, psa.data.star$age.basis)

	bx.data_star <- filter(bx.data.full, subj==star) 

	R_star <- bx.data_star$rc
	W.RC_star <-cbind(1,bx.data_star$std.age, bx.data_star$bx.time)

	d.W.RC<-dim(W.RC_star)[2]
	d.Z<-dim(Z_star)[2]
	#########


	#################
	# Get weigths based on likelihood of a particle, 
	# given subj_star's data (Y, Z, X, and W),
	# and a proposed set of random effects (b-vec & eta) and 
	# hyper-params (beta, gamma.RC)
	W <- rep(NA,P*nreps)
	fit_time<-system.time({  # ~ 35k / second
	for(r in 1:nreps){ 
	for(oo_ind in 1:P){ #index for oo
		p <- (r-1)*P+oo_ind #index for our particle set

		 #likelihood of PSA data
		mu_obs_psa <- Z_star %*% b.vec.star[p,]  + 
			( oo$beta[oo_ind,1] + oo$beta[oo_ind,2]*(eta[p]==1) ) * X_star
		L_Y <- prod(dnorm(Y_star,mean=mu_obs_psa, sd=oo$sigma_res[oo_ind]))

		#likelihood of reclassifications
		if(!is.null(R_star)){
			logit_p_rc<-cbind(W.RC_star,eta[p]) %*% oo$gamma.RC[oo_ind,1:(d.W.RC+1)] 
			L_R <- prod(dbinom(x=R_star,size=1,prob=c(invLogit(logit_p_rc))))
		}else{ 
			L_R <- 1
		}

		W[p] <- L_Y * L_R #weights
	}}})
	#################

	W <- W/sum(W)
	etas_weighted[star] <- crossprod(eta,W)
	# plot(density(eta,weights=W))
}

abs_error <- etas_weighted-eta_true_or_jags
per_error <- abs_error/eta_true_or_jags

hist(star_abs_error,breaks=20)
range(star_abs_error, na.rm=TRUE)

hist(star_per_error)
range(star_per_error, na.rm=TRUE)

plot(eta_true_or_jags,etas_weighted)
abline(0,1)


#Choose an arbitrary subject to plot
star <- round(runif(1,1000-n_subj_to_est+1,1000))

hist(colMeans(of$eta.hat),breaks=30,main=paste0('Distribution of mean posterior eta_hat,\nacross people (star=',star,')'),xlab='Population Distribution')
abline(v=mean(of$eta.hat),lwd=3)
abline(v=eta_true_or_jags[star],lwd=2, col='darkgreen')
abline(v=etas_weighted[star],lwd=3,col='blue',lty=2)
legend('topright',c('overall pop mean','subj* online','subj* full model'),lty=c(1,2,1),lwd=c(3,2,2),col=c('black','blue','darkgreen'))

#dev.copy2pdf(file=paste0('plots/',Sys.Date(),'subj_',star,'_plot.pdf'))

length(Y_star)
length(X_star)

fit_time
nreps*P

