

######################### 
##     Workflow:

# Load data - get the posterior from leaving one subject out. Also get the posterior from fitting on the entire dataset, as something to compare against.
# Generate particles or candidate draws for the posterior for each subject. (`gen_particles`)
# These particles are weighted or accepted/rejected to get an estimated posterior for each subject. 
# All new subjects share the same candidate draws, but have different weights (for importance sampling (IS)) or have different values accepted (for acceptance/rejection sampling (RS)).
# Compare results to assess performance.
#########################


########################
# Packages and basic internal functions
library(dplyr)
library(MASS)
library(ggplot2)

invLogit <- function(x)
	return(exp(x)/(1+exp(x)))
########################


psa_data_full<-read.csv("simulation-data/psa-data-sim.csv")
pt_data_full<-read.csv("simulation-data/pt-data-sim.csv") #true eta is now contained here, in patient data (pt data)
bx_data_full <- read.csv("simulation-data/bx-data-sim.csv")

#Splines for surg and bx
couldve_had_biopsy_all <- !is.na(bx_data_full$bx_here)
did_have_biopsy_all <- couldve_had_biopsy_all & bx_data_full$bx_here




### Output from leave-one-Out JAGS model  (output-out is abbreviated as oo)
### Then collect posterior estimates from output from full JAGS model (of)
#In practice we should separately fit an MCMC leaving each patient out. However, for simplicity of reproducing this computational example, we instead use the posterior after running MCMC on all subjects. 
# In this example, we only use the population parameter posteriors from `oo`. We assume that the posterior for the population parameters will not be substancially different if one subject is left out.
oo <- readRDS('posterior_SEED-0_star-0_crop-FALSE.rds')$sims.list
of <- readRDS('posterior_SEED-0_star-0_crop-FALSE.rds')$sims.list ##this seems to be a repeat of the above object?

#nreps parameter lets you get more out of a limited posterior draw of the population parameters.
#For each population parameter draw, we draw several (nreps) draws for the random effect.
nreps <- 40 




n_post<-length(oo$p_eta)
(P<-n_post*nreps)



missing_etas <- which(is.na(pt_data_full$obs_eta))#=1 if obs and aggressive, 0 if obs and not, or NA if not observed
eta_true_or_jags<-pt_data_full$obs_eta
eta_true_or_jags[missing_etas]<-of$eta_hat_means #indexed by subject

###############


###############
# Generate candidate particles to weight, or accept/reject.

gen_particles<-function(oo,nreps,talk=TRUE){
#set to getOption('verbose') if packaging
#oo = Output from leave-one-Out JAGS object
#nreps = number of times to expand each posterior draw

	# Expand the posterior draws in oo by repeating them nreps times.
		# Reduces the monte-carlo error.

	K<-2  # two latent classes
	n_post<-length(oo$p_eta) # number of particles from JAGS output
	
	P<-n_post*nreps # total number of particles we'll be weighting over

	# Rather than looping operations nreps times,
	# we'll create a new posterior of size nreps * n_post,
	# to more easily vectorize the operations for increased speed_

	#Assign by cycle over nreps times (two vectors of different lengths)
	mu<-array(NA,dim=c(P,2,K))
		#indeces: #[particle, int/slope, k in {1,2}=class]
		#Random effects are intercepts and slopes for age
	mu[,1,1]<-oo$mu_int[,1]
	mu[,1,2]<-oo$mu_int[,2]
	mu[,2,1]<-oo$mu_slope[,1]
	mu[,2,2]<-oo$mu_slope[,2]

	#expand beta
	beta_exp <- matrix(NA,P,1) #beta is no longer beta_k?? No longer 2 columns, just one??
	beta_exp[,1]<-oo$beta[,1] #beta is no longer beta_k??

	#expand sigma_res
	sigma_res_exp<-rep(NA,P)
	sigma_res_exp[]<-oo$sigma_res

	#expand gamma
	gamma_RC_exp <- matrix(NA,P,dim(oo$gamma_RC)[2])
	for(d in 1:dim(oo$gamma_RC)[2])
		gamma_RC_exp[,d]<-oo$gamma_RC[,d]
	
	

	#Get random candidate draws for eta (re-used for each new subject)
	eta<-rbinom(P,1,prob=rep(c(oo$p_eta),times=nreps))
		# Our eta here is analogous to eta_hat from the main JAGS model.
	#Get n_post random draws for b_vec
	b_vec_star <- matrix(NA,P,2)

	(expand_time<-system.time({
	if(talk) pb_sim <- txtProgressBar(min = 0, max = P, char = "=", style=3)
	for(r in 1:nreps){ # nreps = number of times we expand the particles from oo.
	for(oo_ind in 1:n_post){ #index for oo
		p <- (r-1)*n_post+oo_ind #index along the extended particles we're creating.

		cov_for_bvec_p <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2))
		cov_for_bvec_p[1,2]<-
		cov_for_bvec_p[2,1]<-oo$cov_int_slope[oo_ind]
		b_vec_star[p,]<-mvrnorm(1,mu=mu[p,,eta[p]+1], Sigma=cov_for_bvec_p)
			#note, eta[i] is *not* nessecarily the same as eta[i+n_post] due to random draws.
		# we don't store cov_for_bvec_p because it's not used after b_vec_star is generated_
		if(talk) setTxtProgressBar(pb_sim,p)
	}}})) #~ 2 min for nreps=50, n_post=25000

	return(list( #all items in this list are length nreps * n_post.
		#You previously had "exp" suffixes, but not you're ommiting those.
		eta = eta,
		sigma_res = sigma_res_exp,
		beta = beta_exp,
		mu = mu,
		gamma_RC = gamma_RC_exp,
		b_vec_star = b_vec_star
	))
}



#' Get likelihood of each particle,
#' given subj_star's data (Y, Z, X, and V),
#' and a proposed set of random effects (b-vec & eta) and
#' hyper-params (beta, gamma_RC)
#'
#' @param ps particle set (list). Output from get_particles.
#' @param psa_data_star psa data for subject of interest
#' @param bx_data_star biopsy data for subject of interest
get_likelihood<-function(ps, psa_data_star, bx_data_star, verbose=getOption('verbose')){

	P<-length(ps$sigma_res)

	#Setup subject data
	Y_star <- psa_data_star$log_psa
	X_star <- psa_data_star$vol_std #What if this is also a risk factor??
	Z_star <- cbind(1, psa_data_star$age_std) #right!???

	RC_star <- dplyr::filter(bx_data_star, bx_here==1)$rc
	prev_pos_biopsy <- dplyr::filter(bx_data_star, !is.na(bx_here))$prev_G7 #do I have this right, or do I need to exclude periods where bx_here = NA?? Does bx_here=NA mean they've exited the study (reclassified)?


	#Where should we select covariate variables from:
	couldve_had_biopsy_star <- !(is.na(bx_data_star$bx_here))
	did_have_biopsy_star <- couldve_had_biopsy_star & bx_data_star$bx_here

	V_RC_star <- dplyr::mutate(bx_data_star, intercept=1) %>%
		dplyr::select(intercept, age_std, time, time_ns, sec_time_std) %>%
		dplyr::filter(did_have_biopsy_star)

	d_V_RC<-dim(V_RC_star)[2] #Note, these don't include eta yet
	d_Z<-dim(Z_star)[2]

	
	#########
	#likelihood of PSA data
	#To vectorize likelihood fits, we use expanded vectors
		# Expanded vectors are have suffix `_exp` or `_exp`
	# 1) expand data and parameters so that they're grouped by particle, and each particle is associated with a full copy of the dataset.
		# let I = length of data vector, of # of visits
		# a) repeat I-length data vectors P times
		# b) repeate P-length parameter vectors with each = I
	# 2) get likelihood of each visit
	# 3) add a group variable p_ind that groups visits by the particle
	# 4) get the log-likelihood of each particle group

##this is confusing because you are inconsistent with the _exp !!

	if(length(Y_star)==0){
		LL_Y<-0
	}else{
		Z_star_X_bvec<-tcrossprod(ps$b_vec_star,Z_star)
		beta_X_star<-tcrossprod(ps$beta,X_star)
		mu_obs_psa_exp <-c(t(Z_star_X_bvec +beta_X_star)) #take matrix of means, and turn it into a vector.
		Y_star_exp<-rep(Y_star,times=P) #expand Y_star for each particle
		p_ind <- rep(1:P,each=length(Y_star)) #particle index to group Y_star likelihoods
		sigma_res_exp<-rep(ps$sigma_res,each=length(Y_star))
		LL_Y_j <- log(dnorm(Y_star_exp,mean=mu_obs_psa_exp, sd=sigma_res_exp)) #the likelihood for each visit, grouped by particle.
		LL_Y<-( #logLik of all visits
			data_frame('LL_Y_all'=LL_Y_j,'p_ind'=as.factor(p_ind))%>%
			group_by(p_ind) %>%
			summarize(sum=sum(LL_Y_all))
			)$sum
	}

	


	#' For RC, get the likelihood for all visits from a subject.
	#' This returns a vector of length P, with the joint log-likelihood of all visits, under each particle
	#' It works under the model that for any particle p,
	#' we have a bernoulli outcome with:
	#'  logit(mean)= V %*% gamma[p,1:dim(V)] + eta * gamma[p,dim(V)+1]
	get_joint_LL_measurements<-function(V,outcomes,gamma,eta){ 
		
		if(length(outcomes)==0) return(0)

		nVisits <- length(outcomes)
		if(nVisits != dim(V)[1]) error('Outcome length does not match covariate length')

		P <- dim(gamma)[1] #shouldn't need to be redefined?
		d_V <- dim(V)[2]
		p_ind <- rep(1:P,each=nVisits) #particle index to group

		eta_gamma <- gamma[,d_V+1] * eta
		eta_gamma_exp <- rep(eta_gamma, each = nVisits)
		V_gamma_exp <-  c(tcrossprod(as.matrix(V), gamma[,1:d_V] ))
		logit_p_exp <- eta_gamma_exp + V_gamma_exp

		p_exp<-c(t(invLogit(logit_p_exp)))  ##why c(t()) ? its already a vector
		outcomes_exp<-rep(outcomes,times=P)		

		LL_j <- log(dbinom(x=outcomes_exp,size=1,prob=p_exp))
		LL <- ( #logLik of all visits
			data_frame(LL_all=c(t(LL_j)),ind=p_ind) %>% #confused by c(t()) again
			group_by(ind) %>%
			summarize(sum=sum(LL_all))
			)$sum

		return(LL)
	}


	LL_RC <- get_joint_LL_measurements(
		V=V_RC_star,
		outcomes=RC_star,
		gamma=ps$gamma_RC,
		eta=ps$eta)

	W <- exp(LL_Y + LL_RC) #I may have messed up on my earlier revisions (changing W cov mats to V) since you named the likelihood for each particle W


	return(W)
}






##############
# Estimate posteriors for each subject.

# For each subject, calcualte and store:
	# for importance sampling
		# * the weighted eta posterior mean
		# * the effective sample size for the number of posterior draws
	# for Rejection sampling
		# * average eta over accepted draws
		# * number of accepted draws from the posterior



###### Function to get posterior means for each subject

#' A function with references to data objects in parent environment.
#' useful to write it this way so we can re-update subject's estimates.
#' We have this depend on the particle set (ps) so that we can
#' easily redo it for specific subjects with a different particle set.
#' @param star the index of the subject to fit.
posterior_star<-function(star,ps,runifs,rej_const=NULL){
	
	#data for subj_star

	likelihood <- get_likelihood(
		ps=ps, #!! GLOBAL REFS !!
		psa_data_star=filter(psa_data_full, subj == star),
		bx_data_star=filter(bx_data_full, subj==star)
		)
	
	## Importance Weighting ##
	W <- likelihood/sum(likelihood)
	etas_IS_star <- crossprod(W,ps$eta) 
	effective_ss_star <- 1/crossprod(W)
	

	## Rejection Sampling ##

	#we don't know integrating constant. Does that matter if we just set this to max(likelihood)??
	if(is.null(rej_const)) rej_const <- max(likelihood)
	accept_ind <- (likelihood/rej_const) >= runifs
	num_accepted_star <- sum(accept_ind)
	etas_RS_star <- mean(ps$eta[accept_ind])

	return(list(
		etas_IS_star = etas_IS_star,
		effective_ss_star = effective_ss_star,
		num_accepted_star = num_accepted_star,
		etas_RS_star = etas_RS_star,
		W=W
	))
}
######
















##############################
####### Run Functions
##############################




##### Generate particles:

seed<-101
set.seed(seed)

#`ps` = particle set
ps <- gen_particles(oo,nreps=nreps,talk=TRUE)
runifs <- runif(P) #needed for rejection sampling later





######################## 
# Run loop over subjects


###### Vectors to store results
N <- max(psa_data_full$subj) #number of subjects

# store:
# * effective sample for importance weighting
# * number of proposals accepted by rejection sampling
# * mean among proposals accepted by RS
# * posterior mean estimates for IS
fit_time <-
effective_ss <-
num_accepted <-
etas_RS <-
etas_IS <- rep(NA,N)


# Some entries will be NAs at the end of this script,
# to maintain indexing system of:
	# subj_star's eta = etas_IS[star]
######



#To test for star = a random number, use
	#star <- round(runif(1,N-n_subj_to_est+1,N))

pb <- txtProgressBar(min = min(missing_etas), max = max(missing_etas), char = "=", style=3)
for(star in missing_etas){

	rej_const=NULL
	fit_time[star]<-system.time({
		post_star <- posterior_star(star,ps,rej_const=rej_const,runifs=runifs)
	})['elapsed']

	etas_IS[star] <- post_star$etas_IS_star
	effective_ss[star] <- post_star$effective_ss_star
	num_accepted[star] <- post_star$num_accepted_star
	etas_RS[star] <- post_star$etas_RS_star

	setTxtProgressBar(pb, star)
}

save.image(file='checkpoint_IS_fit.RData')
# load('checkpoint_IS_fit.RData')








##############################
####### Performance & Agreement
##############################


#Show the distribution (over patients) for effective sample sizes for IS.
if(interactive()) 
	ggplot(as_data_frame(effective_ss))+ geom_histogram(aes(x=effective_ss))+scale_x_log10() + labs(title='Effective Sample Size')


#IS and RS give very similar results
plotTitle <- 'Agreement between Estimated Probabilities\nof Having Aggressive Cancer'

if(interactive()){
	plot(etas_IS,etas_RS,cex=.5,xlim=0:1,ylim=0:1,main=plotTitle,xlab='Importance Sampling', ylab='Rejection Sampling')
	abline(0,1)
}

#However, IS appears to have a slightly closer resemblance to MCMC draws
squared_errors_IS <- (etas_IS - eta_true_or_jags)^2
squared_errors_RS <- (etas_RS - eta_true_or_jags)^2
mean(squared_errors_IS,na.rm=TRUE) #IS mean squared error
mean(squared_errors_RS,na.rm=TRUE) #RS mean squared error








######### Save plots of results

eff_ss_error_data<-data_frame(etas_IS,eta_true_or_jags,squared_errors_IS,squared_errors_RS,effective_ss)
tail(eff_ss_error_data)

png('effective_SS_and_error.png',width = 6, height = 6, type='cairo')
ggplot(eff_ss_error_data) + 
	geom_point(aes(y=sqrt(squared_errors_IS),x=effective_ss),alpha=.5)+ 
	labs(title='Effective Sample Size v. Absolute Error\n(plotted on log scale)',x='Effective Sample Size for IS',y='Absolute Error')+
	theme(text=element_text(size=18)) +
	scale_x_continuous(breaks=10^(1:6)/2) +
	scale_y_continuous(breaks=c(.03,.01,.005,.001)) +
	coord_trans(y = "log10",x= "log10")
dev.off()

png('agreement_MCMC_IS.png',width = 5.5, height = 5.5, type='cairo')
plot(x=etas_IS,y=eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1,ylab='Estimates from MCMC',xlab='Estimates from Importance Sampling',main='Agreement between Estimated Probabilities\nof Having Aggressive Cancer')
abline(0,1)
dev.off()
