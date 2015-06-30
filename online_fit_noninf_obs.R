#setwd("/Users/aaronfisher/Dropbox/Future Projects/inHealth Prostate Screening/repo")




######### Workflow: 
# Load data - get the posterior from leaving one subject out. Also get the posterior from fitting on the entire dataset, as something to compare against.
# Generate particles or candidate draws for the posterior for each subject. (`gen_particles`)
# These particles are weighted or accepted/rejected to get an estimated posterior for each subject. Subjects share the same candidate draws, but have different weights (for importance sampling (IS)) or have different values accepted (for acceptance/rejection sampling (RS)).
# Compare results to assess performance. 
######### 


########################
# Packages and basic internal functions
library(dplyr)
library(MASS)
library(ggplot2)

invLogit <- function(x)
	return(exp(x)/(1+exp(x)))
########################








###############
# Load Data


psa.data.full<-read.csv("simulation-data/psa-data-sim.csv")
pt.ordered.full<-read.csv("simulation-data/pt-data-sim.csv")
bx.data.full <- read.csv("simulation-data/bx-data-sim.csv")

### Output from leave-one-Out JAGS model  (output-out -> oo)
	# For now just use the output from the full model as an approximate replacement
oo <- readRDS('posterior_full_100k.rds')$sims.list
#oo <- readRDS('2015-06-05_posterior_full_100k_seed_5.rds')$sims.list #on cluster
P<-length(oo$p_eta)
nreps <- 50
(P*nreps)

### Collect posterior estimates from full JAGS model
of<-readRDS('posterior_full_100k.rds')$sims.list
#of <- readRDS('2015-06-05_posterior_full_100k_seed_5.rds')$sims.list
known_etas<-dplyr::filter(pt.ordered.full,obs.eta==1) %>%
	dplyr::select(eta.true)
eta_true_or_jags<-c(known_etas$eta.true,colMeans(of$eta.hat))

###############


###############
# Generate candidate particles to weight, or accept/reject. 

gen_particles<-function(oo,nreps,talk=TRUE){
#set to getOption('verbose') if packaging
#oo = output from leave-one-out JAGS object
#nreps = number of times to expand each posterior draw

	# Expand the posterior draws in oo by repeating them nreps times.
		# Reduces the monte-carlo error.

	K<-2  # two latent classes
	P<-length(oo$p_eta) # number of particles from JAGS output
	
	P*nreps # total number of particles we'll be weighting over

	# Rather than looping operations nreps times,
	# we'll create a new posterior of size nreps * P,
	# to more easily vectorize the operations for increased speed.

	#Assign by cycle over nreps times (two vectors of different lengths)
	mu<-array(NA,dim=c(nreps*P,3,K))
		#indeces: #[particle, int/slope/spline, k in {1,2}=class]
	mu[,1,1]<-oo$mu_int[,1]
	mu[,1,2]<-oo$mu_int[,2]
	mu[,2,1]<-oo$mu_slope[,1]
	mu[,2,2]<-oo$mu_slope[,2]
	mu[,3,1]<-oo$mu_spline[,1]
	mu[,3,2]<-oo$mu_spline[,2]

	#expand beta
	beta_exp <- matrix(NA,P*nreps,2)
	beta_exp[,1]<-oo$beta[,1]
	beta_exp[,2]<-oo$beta[,2]

	#expand sigma_res
	sigma_res_exp<-rep(NA,P*nreps)
	sigma_res_exp[]<-oo$sigma_res

	#expand gamma
	gamma.RC.exp <- matrix(NA,P*nreps,dim(oo$gamma.RC)[2])
	gamma.RC.exp[,1]<-oo$gamma.RC[,1]
	gamma.RC.exp[,2]<-oo$gamma.RC[,2]
	gamma.RC.exp[,3]<-oo$gamma.RC[,3]
	gamma.RC.exp[,4]<-oo$gamma.RC[,4]


	if(nreps>1){
		#if TRUE then cycling worked!
		all(sigma_res_exp[1:10]==sigma_res_exp[1:10+P])
		all(beta_exp[1:10,]==beta_exp[1:10+P,]) 
		all(mu[1:10,,]==mu[1:10+P,,]) 
	}

	#Get random draws for eta (re-used for each subject)
	eta<-rbinom(P*nreps,1,prob=rep(c(oo$p_eta),times=nreps))
		# Our eta here is analogous to eta.hat from the main JAGS model.
	#Get P random draws for b.vec
	b.vec.star <- matrix(NA,P*nreps,3)
	cov_for_bvec <- list() #Also consider *not* storing this, since it's not used after b.vec.star is generated.
	(expand_time<-system.time({
	if(talk) pb_sim <- txtProgressBar(min = 0, max = P*nreps, char = "=", style=3)
	for(r in 1:nreps){ # nreps = number of times we expand the particles from oo.
	for(oo_ind in 1:P){ #index for oo
		p <- (r-1)*P+oo_ind #index along the extended particles we're creating.

		cov_for_bvec[[p]] <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2,oo$sigma_spline[oo_ind]^2))
		cov_for_bvec[[p]][1,2]<-
		cov_for_bvec[[p]][2,1]<-oo$cov_int_slope[oo_ind]
		cov_for_bvec[[p]][1,3]<-
		cov_for_bvec[[p]][3,1]<-oo$cov_int_spline[oo_ind]
		cov_for_bvec[[p]][2,3]<-
		cov_for_bvec[[p]][3,2]<-oo$cov_slope_spline[oo_ind]
		b.vec.star[p,]<-mvrnorm(1,mu=mu[p,,eta[p]+1], Sigma=cov_for_bvec[[p]])
			#note, eta[i] is *not* nessecarily the same as eta[i+P] due to random draws.

		if(talk) setTxtProgressBar(pb_sim,p)
	}}})) #~ 2 min for nreps=50, P=25000


	return(list(
		eta = eta,
		sigma_res_exp = sigma_res_exp,
		beta_exp = beta_exp,
		mu = mu,
		gamma.RC.exp = gamma.RC.exp,
		cov_for_bvec = cov_for_bvec,
		b.vec.star = b.vec.star
	))
}




#' Get likelihood of each particle, 
#' given subj_star's data (Y, Z, X, and W),
#' and a proposed set of random effects (b-vec & eta) and 
#' hyper-params (beta, gamma.RC)
#'
#' @param ps particle set (list). Output from get_particles.
#' @param psa.data.star psa data for subject of interest
#' @param bx.data_star biopsy data for subject of interest
get_likelihood<-function(ps,psa.data.star,bx.data_star){

	PP<-length(ps[[1]]) #!! change later to just P and P_pre

	######### Setup
	#Abbreviate particle variables
	eta <- ps$eta
	sigma_res_exp <- ps$sigma_res_exp
	beta_exp <- ps$beta_exp
	mu <- ps$mu
	gamma.RC.exp <- ps$gamma.RC.exp
	cov_for_bvec <- ps$cov_for_bvec
	b.vec.star <- ps$b.vec.star

	#Setup subject data
	Y_star <- psa.data.star$log.psa
	X_star <- psa.data.star$std.vol
	Z_star <- cbind(1, psa.data.star$age.std)

	R_star <- bx.data_star$rc
	W.RC_star <-cbind(1,bx.data_star$age.std,bx.data_star$time, bx.data_star$time.ns, bx.data_star$sec.time.std)

	d.W.RC<-dim(W.RC_star)[2]
	d.Z<-dim(Z_star)[2]
	#########


	#########
	#likelihood of PSA data
	#To vectorize likelihood fits, we use expanded vectors
		# Expanded vectors are have suffix `_exp` or `.exp`
	# 1) repeat data vectors PP times
	# 2) get likelihood of each visit
	# 3) add a group variable p_ind that groups visits by the particle

	Z_star_X_bvec<-tcrossprod(b.vec.star,Z_star)
	beta_eta<-beta_exp[,1] + beta_exp[,2]*(eta==1)
	beta_exp_eta_Xstar<-tcrossprod(beta_eta,X_star)
	mu_obs_psa_exp <-c(t(Z_star_X_bvec +beta_exp_eta_Xstar))
	Y_star_exp<-rep(NA,length(Y_star)*PP)
	Y_star_exp[]<-Y_star #cycle Y_star up
	p_ind <- rep(1:(PP),each=length(Y_star)) #particle index to group Y_star likelihoods
	sigma_res_exp2<-rep(sigma_res_exp,each=length(Y_star))
	L_Y_j <- dnorm(Y_star_exp,mean=mu_obs_psa_exp, sd=sigma_res_exp2) #the likelihood for each visit, grouped by particle.
	L_Y_frame<-data.frame('L_Y_all'=L_Y_j,'p_ind'=as.factor(p_ind))%>%
		group_by(p_ind) %>%
		summarize(prod=prod(L_Y_all)) 
	L_Y <- L_Y_frame$prod


	if(!is.null(R_star)){
		logit_p_rc_exp<-gamma.RC.exp[,1:3] %*% t(W.RC_star) + gamma.RC.exp[,4]*eta 

		R_star_exp<-rep(NA,length(R_star)*PP)
		R_star_exp[]<-R_star
		prob_R_exp<-c(t(invLogit(logit_p_rc_exp)))
		L_R_j <- matrix(dbinom(x=R_star_exp,size=1,prob=prob_R_exp),PP,length(R_star),byrow=TRUE)
		L_R_frame <- data.frame(L_R_all=c(t(L_R_j)),ind=rep(1:(PP),each=length(R_star))) %>%
			group_by(ind) %>%
			summarize(prod=prod(L_R_all))
		L_R <- L_R_frame$prod
	}else{
		L_R <- 1
	}

	return(L_R*L_Y)


	# #####
	# Alternate (older) version, looping nreps times instead of vectorizing.
	# fit_time<-system.time({  # ~ 35k / second
	# L_Y2_log<-
	# L_R2_log<-
	# W2 <- rep(NA,PP)
	# for(r in 1:nreps){ 
	# for(oo_ind in 1:P){ #index for oo
	# 	p <- (r-1)*P+oo_ind #index for our particle set

	# 	 #likelihood of PSA data
	# 	mu_obs_psa <- Z_star %*% b.vec.star[p,]  + 
	# 		( oo$beta[oo_ind,1] + oo$beta[oo_ind,2]*(eta[p]==1) ) * X_star
	# 	L_Y2 <- prod(dnorm(Y_star,mean=mu_obs_psa, sd=oo$sigma_res[oo_ind]))
	# 	L_Y2_log[p]<-L_Y2

	# 	#likelihood of reclassifications
	# 	if(!is.null(R_star)){
	# 		logit_p_rc<-cbind(W.RC_star,eta[p]) %*% oo$gamma.RC[oo_ind,1:(d.W.RC+1)] 
	# 		L_R2 <- prod(dbinom(x=R_star,size=1,prob=c(invLogit(logit_p_rc))))
	# 		L_R2_log[p]<-L_R2
	# 	}else{ 
	# 		L_R2 <- 1
	# 	}

	# 	W2[p] <- L_Y2 * L_R2 #weights
	# }}})
	# #################

	# W2 <- W2/sum(W2)
	# if(any(W!=W2)) browser()
	# all(L_R==L_R2_log)
	# all(L_Y==L_Y2_log)

	# etas_IS[star] <- crossprod(eta,W2)
	######
	######
	######
}





##############################
# Generate particles:

seed<-101
set.seed(seed)

#`ps` = particle set
ps <- gen_particles(oo,nreps=nreps,talk=TRUE)
runifs <- runif(P*nreps) #needed for rejection sampling later





##############################
# Estimate posteriors for each subject.

# For each subject, calcualte and store:
	# for Importance sampling
		# * the weighted eta posterior mean
		# * the effective sample size for the number of posterior draws
	# for Rejection sampling
		# * average eta over accepted draws
		# * number of accepted draws from the posterior



###### Vectors to store results
N <- max(psa.data.full$subj) #number of subjects

# store:
# * effective sample for importance weighting
# * number of proposals accepted by rejection sampling 
# * mean among proposals accepted by RS
# * posterior mean estimates for IS
effective_ss <- 
num_accepted <- 
etas_RS <-  
etas_IS <- rep(NA,N)

# Some entries will be NAs at the end of this script, 
# to maintain indexing system of:
	# subj_star's eta = etas_IS[star]
######



###### Function to get posterior means for each subject

#' A function with references to data objects in parent environment. 
#' useful to write it this way so we can re-update subject's estimates.
#' We have this depend on the particle set (ps) so that we can
#' easily redo it for specific subjects with a different particle set.
#' @param star the index of the subject to fit.
posterior_star<-function(star,ps,rej_const=NULL){
	
	#data for subj_star

	likelihood <- get_likelihood(
		ps=ps,
		psa.data.star=filter(psa.data.full, subj == star),
		bx.data_star=filter(bx.data.full, subj==star)
		)
	
	## Importance Weighting ## 

	W <- likelihood/sum(likelihood)
	etas_IS_star <<- crossprod(W,ps$eta)
	effective_ss_star <<- 1/crossprod(W)
	if(FALSE) plot(density(ps$eta))
	if(FALSE) plot(density(ps$eta,weights=W))


	## Rejection Sampling ##

	#we don't know integrating constant. Does that matter if we just set this to max(W)??
	if(is.null(rej_const)) rej_const <- max(likelihood)
	accept_ind <- (likelihood/rej_const) >= runifs
	num_accepted_star <<- sum(accept_ind)
	etas_RS_star <<- mean(ps$eta[accept_ind])

	return(list(
		etas_IS_star = etas_IS_star,
		effective_ss_star = effective_ss_star,
		num_accepted_star = num_accepted_star,
		etas_RS_star = etas_RS_star,
		W=W
	))
}
######


###### Run loop over subjects


#To test for star = some random number, use
	#star <- round(runif(1,N-n_subj_to_est+1,N)) 
n_subj_to_est <- dim(oo$eta.hat)[2]
subj2fit <- (N-n_subj_to_est+1):N #subjects to estimate eta

pb <- txtProgressBar(min = 0, max = n_subj_to_est, char = "=", style=3)
for(star in subj2fit){
	
	# rej_const<-0.0007 #!!?? need to get better version for. If ratio is the likelihood, then rej_const is the value at the MLE. We should know this from what we generated from? Equal to the latent value or something?
	# if(star>1) rej_const <- max(rej_const,likelihood) #this tends to get too small.
	rej_const=NULL
	fit_time<-system.time({
		post_star <- posterior_star(star,ps,rej_const=rej_const)
	})

	etas_IS[star] <- post_star$etas_IS_star
	effective_ss[star] <- post_star$effective_ss_star
	num_accepted[star] <- post_star$num_accepted_star
	etas_RS[star] <- post_star$etas_RS_star

	setTxtProgressBar(pb, star-(N-n_subj_to_est+1))
}


#save.image(file=paste0(Sys.Date(),'_seed_',seed,'_online_fit_results.RData'))
#load('2015-06-08_seed_101_online_fit_results_with_refitv2.RData')




######################################
# Performance

# IS appears to do better than RS

 #NAs for etas_IS are removed from eta_true_or_jags in plot()
plot(etas_IS,etas_RS,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)
plot(eta_true_or_jags,etas_RS,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)
plot(eta_true_or_jags,etas_IS,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)

skip <- is.na(etas_IS)
cor(eta_true_or_jags[!skip],etas_IS[!skip])
cor(eta_true_or_jags[!skip],etas_RS[!skip])

errors_IS <- (eta_true_or_jags - etas_IS)^2
errors_RS <- (eta_true_or_jags - etas_RS)^2
(RMSE_weighted<- sqrt(mean(errors_IS,na.rm=TRUE)))
(RMSE_accepted<- sqrt(mean(errors_RS,na.rm=TRUE)))
sd(eta_true_or_jags[!skip])
quantile(errors_IS[!skip],c(.95,.99,.995,.999,.9999,1))
quantile(errors_RS[!skip],c(.95,.99,.995,.999,.9999,1))


####################
# Why does IS and RS do so badly for some subjects?
# It appears to be because, for those subjects, few proposals being accepted. (Or small effective sample size for IS)

# Which people do we have bad estimates (for IS or RS)
bad_IS <- abs(eta_true_or_jags-etas_IS)>.05
bad_RS <- abs(eta_true_or_jags-etas_RS)>.05
which(bad_IS)
which(bad_RS)
# Previously you saw, for the importance weighted people with high error, for which(bad_IS), you got
#seed = 101 gives *265* 499 *514* *532* *619* ***972***
#seed = 0 gives 254 *514* *532* *619* 893 **972**
# (for unset seed): *265* 460 463 577 637 **972**
# (unset seed): **972**, 303, *514*, *619*, 809



######## Rejection sampling & Num accepted

quantile(num_accepted,c(0,.005,.01,.02,.05,.1,.5),na.rm=TRUE)

quantile(num_accepted,na.rm=TRUE)
mean(num_accepted<100,na.rm=TRUE)

plot(eta_true_or_jags,etas_RS,cex=.5,xlim=0:1,ylim=0:1,col=(num_accepted<100)+1)
quantile(abs(eta_true_or_jags-etas_RS)[num_accepted>100],
	c(.90,.95,.99,.999),
	na.rm=TRUE) #yes! this appears to identify they well!

error_df<-data.frame(error=sqrt(errors_RS),num_accepted=num_accepted)
ggplot(error_df,aes(x=error,y=num_accepted)) + geom_point(alpha=.4)+ coord_trans(y = "log10")


######### Importance weighting

hist(effective_ss,breaks=100)
hist(effective_ss,breaks=200,xlim=c(1,10000))
quantile(effective_ss,c(0,.005,.01,.02,.05,.1,.5),na.rm=TRUE)

eff_ss_cutoff <- 500

plot(eta_true_or_jags,etas_IS,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1)
abline(0,1)

quantile(sqrt(errors_IS)[effective_ss>eff_ss_cutoff],
	c(.90,.95,.99,.999),
	na.rm=TRUE)
quantile(sqrt(errors_IS),
	c(.90,.95,.99,.999),
	na.rm=TRUE) #yup, makes a huge difference!

error_df<-data.frame(error=sqrt(errors_IS),effective_sample_size=effective_ss)
ggplot(error_df,aes(x=error,y=effective_sample_size)) + geom_point(alpha=.4)+ coord_trans(y = "log10")


##############################
# REDO IS & RS based on effective sample size in IS

#How many do we need to refit?
sum(effective_ss < eff_ss_cutoff, na.rm=TRUE)

subj2refit_v2<- effective_ss < eff_ss_cutoff

etas_IS_v2 <-  etas_IS
effective_ss_v2 <-  effective_ss
num_accepted_v2 <-  num_accepted
etas_RS_v2 <-  etas_RS

nreps_v2 <- 200

#important! re-gen the particles!
seed2 <- 99
set.seed(seed2)
ps2 <- gen_particles(oo,nreps=nreps_v2,talk=TRUE)

for(star in which(subj2refit_v2)){
	
	# rej_const<-0.0007 #!!?? need to get better version for. If ratio is the likelihood, then rej_const is the value at the MLE. We should know this from what we generated from? Equal to the latent value or something?
	# if(star>1) rej_const <- max(rej_const,likelihood) #this tends to get too small.
	rej_const=NULL
	fit_time_v2<-system.time({
		post_star <- posterior_star(star,ps2,rej_const=rej_const)
	})

	etas_IS_v2[star] <- post_star$etas_IS_star
	effective_ss_v2[star] <- post_star$effective_ss_star
	num_accepted_v2[star] <- post_star$num_accepted_star
	etas_RS_v2[star] <- post_star$etas_RS_star

	setTxtProgressBar(pb, star-(N-n_subj_to_est+1))
}

# save.image(file=paste0(Sys.Date(),'_seed_',seed,'_online_fit_results_with_refitv2.RData'))
# load('2015-06-08_seed_101_online_fit_results_with_refitv2.RData')

sum( effective_ss_v2 < eff_ss_cutoff ,na.rm=TRUE) #how many do we still have to refit?

hist(effective_ss_v2,breaks=100)
hist(effective_ss_v2,breaks=200,xlim=c(1,10000))

quantile(effective_ss,c(0,.005,.01,.02,.05,.1,.5),na.rm=TRUE)
quantile(effective_ss_v2,c(0,.005,.01,.02,.05,.1,.5),na.rm=TRUE)

of2<-readRDS('posterior_full_100k.rds')$sims.list
eta_true_or_jags2<-c(known_etas,colMeans(of2$eta.hat))


png(file=paste0('plots/',Sys.Date(),'_compare_fits.png'),height=1200,width=500,pointsize=20)
par(mfrow=c(3,1),mar=c(4,4,3,1))
plot(eta_true_or_jags,etas_IS,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='First attempt',xlab='JAGS',ylab='Importance Sampling')
legend('topleft',c('Flagged','Normal'),col=c('red','black'),pch=1)
abline(0,1,lty=3)
plot(eta_true_or_jags,etas_IS_v2,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='Refitting',xlab='JAGS',ylab='Importance Sampling 2')
abline(0,1,lty=3)
plot(eta_true_or_jags,eta_true_or_jags2,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='References Coherence',xlab='JAGS 1',ylab='JAGS 2')
abline(0,1,lty=3)

dev.off()



errors_IS_v2 <- (eta_true_or_jags - etas_IS_v2)^2
errors_RS_v2 <- (eta_true_or_jags - etas_RS_v2)^2
sqrt(mean(errors_IS_v2,na.rm=TRUE))
sqrt(mean(errors_RS_v2,na.rm=TRUE))

quantile(sqrt(errors_IS_v2)[effective_ss_v2>eff_ss_cutoff],
	c(.90,.95,.99,.999,1),
	na.rm=TRUE)
max(sqrt(errors_IS_v2),na.rm=TRUE)
max(sqrt(errors_IS),na.rm=TRUE)

cor(eta_true_or_jags[!skip],etas_IS_v2[!skip])


#Choose an arbitrary subject to plot
star <- round(runif(1,N-n_subj_to_est+1,N))

hist(colMeans(of$eta.hat),breaks=30,main=paste0('Distribution of mean posterior eta_hat,\nacross people (star=',star,')'),xlab='Population Distribution')
abline(v=mean(of$eta.hat),lwd=3)
abline(v=eta_true_or_jags[star],lwd=2, col='darkgreen')
abline(v=etas_IS[star],lwd=3,col='blue',lty=2)
legend('topright',c('overall pop mean','subj* online','subj* full model'),lty=c(1,2,1),lwd=c(3,2,2),col=c('black','blue','darkgreen'))

#dev.copy2pdf(file=paste0('plots/',Sys.Date(),'subj_',star,'_plot.pdf'))




length(Y_star)
length(X_star)


nreps*P





####################################
####################################
####################################
####################################
#Workspace

if(FALSE){

#Vectorized version

 #likelihood of PSA data
system.time({
Z_star_X_bvec<-tcrossprod(b.vec.star,Z_star)
beta_eta<-beta_exp[,1] + beta_exp[,2]*(eta==1)
beta_exp_eta_Xstar<-tcrossprod(beta_eta,X_star)
mu_obs_psa_exp <-c(t(Z_star_X_bvec +beta_exp_eta_Xstar))
Y_star_exp<-rep(NA,length(Y_star)*nreps*P)
Y_star_exp[]<-Y_star #cycle Y_star up
p_ind <- rep(1:(nreps*P),each=length(Y_star))
sigma_res_exp2<-rep(sigma_res_exp,each=length(Y_star))
L_Y_j <- dnorm(Y_star_exp,mean=mu_obs_psa_exp, sd=sigma_res_exp2) #the likelihood for each visit, grouped by particle.
L_Y_frame<-data.frame('L_Y_all'=L_Y_j,'p_ind'=as.factor(p_ind))%>%
	group_by(p_ind) %>%
	summarize(prod=prod(L_Y_all)) 
L_Y2<-L_Y_frame$prod

})


system.time({
gamma.RC.exp <- matrix(NA,P*nreps,dim(oo$gamma.RC)[2])
gamma.RC.exp[,1]<-oo$gamma.RC[,1]
gamma.RC.exp[,2]<-oo$gamma.RC[,2]
gamma.RC.exp[,3]<-oo$gamma.RC[,3]
gamma.RC.exp[,4]<-oo$gamma.RC[,4]


logit_p_rc<-gamma.RC.exp[,1:3] %*% t(W.RC_star) + gamma.RC.exp[,4]*eta 

L_R_j <- matrix(dbinom(x=R_star,size=1,prob=c(invLogit(logit_p_rc))),nreps*P,length(R_star))
L_R_frame <- data.frame(L_R_all=c(t(L_R_j)),ind=rep(1:(nreps*P),each=length(R_star))) %>%
	group_by(ind) %>%
	summarize(prod=prod(L_R_all))
L_R2<-L_R_frame$prod

})


#LOOP VERSION

L_Y_log<-
L_R_log<-
W <- rep(NA,P*nreps)
fit_time<-system.time({  # ~ 35k / second
for(r in 1:nreps){ 
for(oo_ind in 1:P){ #index for oo
	p <- (r-1)*P+oo_ind #index for our particle set

	 #likelihood of PSA data
	mu_obs_psa <- Z_star %*% b.vec.star[p,]  + 
		( oo$beta[oo_ind,1] + oo$beta[oo_ind,2]*(eta[p]==1) ) * X_star
	L_Y <- prod(dnorm(Y_star,mean=mu_obs_psa, sd=oo$sigma_res[oo_ind]))
	L_Y_log[p]<-L_Y

	#likelihood of reclassifications
	if(!is.null(R_star)){
		logit_p_rc<-cbind(W.RC_star,eta[p]) %*% oo$gamma.RC[oo_ind,1:(d.W.RC+1)] 
		L_R <- prod(dbinom(x=R_star,size=1,prob=c(invLogit(logit_p_rc))))
		L_R_log[p]<-L_R
	}else{ 
		L_R <- 1
	}

	W[p] <- L_Y * L_R #weights
}}})
#################

W <- W/sum(W)
etas_IS[star] <- crossprod(eta,W)






}