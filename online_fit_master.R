


######################### 
##     Workflow:

# Load data - get the posterior from leaving one subject out. Also get the posterior from fitting on the entire dataset, as something to compare against.
# Generate particles or candidate draws for the posterior for each subject. (`gen_particles`)
# These particles are weighted or accepted/rejected to get an estimated posterior for each subject. 
# All new subjects share the same candidate draws, but have different weights (for importance sampling (IS)) or have different values accepted (for acceptance/rejection sampling (RS)).
# Compare results to assess performance.
#########################



# !! go over order of params (e.g. the 2nd dimension of the W's to make sure you coded it in the right order).





########################
# Packages and basic internal functions
library(dplyr)
library(MASS)
library(ggplot2)
#library(splines)

invLogit <- function(x)
	return(exp(x)/(1+exp(x)))
########################


#IOP (Informative observation process)
IOP_BX <- TRUE #Informative observation process for biopsy
IOP_SURG <- TRUE #Informative observation process for surgery
leave_one_out <- FALSE

IOPs<-paste0(
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG')
batch_path<-paste0(
	'batches/',
	IOPs,
	'/leave_one_out_',leave_one_out,
	'/batch-1/')
posterior_path<-paste0(batch_path,
	'concatenated_posterior.rds')
# posterior_path2<-paste0('batches/',
# 	IOPs,
# 	'/2/concatenated_posterior.rds')

###############
# Load Data


psa_data_full<-read.csv("simulation-data/psa-data-sim.csv")
pt_data_full<-read.csv("simulation-data/pt-data-sim.csv") #true eta is now contained here, in patient data (pt data)
bx_data_full <- read.csv("simulation-data/bx-data-sim.csv")

#Splines for surg and bx
couldve_had_biopsy_all <- !is.na(bx_data_full$bx_here)
did_have_biopsy_all <- couldve_had_biopsy_all & bx_data_full$bx_here

if(IOP_SURG|IOP_BX){
	sum(couldve_had_biopsy_all)
	sum(did_have_biopsy_all)
}





### Output from leave-one-Out JAGS model  (output-out -> oo)
	# For now just use the output from the full model as an approximate replacement
### Then collect posterior estimates from full JAGS model
# of <- readRDS('batches/inf-obs/1/2015-07-24_JAGS_full_sample_summarized_12.rds')
oo <- readRDS(posterior_path)
of <- readRDS(posterior_path)
of2 <- readRDS(posterior_path2)

nreps <- 1



n_post<-length(oo$p_eta)
(P<-n_post*nreps)



missing_etas <- which(is.na(pt_data_full$obs_eta)) #=1 if observed and aggressive, 0 if observed and not aggressive, or NA if not observed
eta_true_or_jags2<-
eta_true_or_jags<-pt_data_full$obs_eta
eta_true_or_jags[missing_etas]<-colMeans(of$eta_hat_means) #indexed by subject
eta_true_or_jags2[missing_etas]<-colMeans(of2$eta_hat_means) #indexed by subject
#sum(missing_etas)==length(colMeans(of$eta_hat))

range(abs(eta_true_or_jags2-eta_true_or_jags)) #largest error 
sqrt(mean((eta_true_or_jags2- eta_true_or_jags)^2, na.rm=TRUE)) #rMSD
cor(colMeans(of2$eta_hat_means),
	colMeans( of$eta_hat_means))
###############


###############
# Generate candidate particles to weight, or accept/reject.

gen_particles<-function(oo,nreps,verbose=TRUE){
#set to getOption('verbose') if packaging
#oo = Output from leave-one-Out JAGS object
#nreps = number of times to expand each posterior draw

	# Expand the posterior draws in oo by repeating them nreps times.
	# This reduces the monte-carlo error.
	# Loop first over n_post, then increment nreps by 1, then loop over n_post again etc.
	# In this way there is as much diversity as possible in the *head* of each posterior draw. 
	# This frontloaded diversity will be useful when we adaptively draw a number of particles to get a good effective sample size.


	K<-2  # two latent classes
	n_post<-length(oo$p_eta) # number of particles from JAGS output
	
	P<-n_post*nreps # total number of particles we'll be weighting over

	# Rather than looping operations nreps times,
	# we'll create a new posterior of size nreps * n_post,
	# to more easily vectorize the operations for increased speed_

	#Assign by cycle over nreps times (two vectors of different lengths)
	p_eta_exp<-rep(c(oo$p_eta),times=nreps)

	mu<-array(NA,dim=c(P,2,K))
		#indeces: #[particle, int/slope, k in {1,2}=class]
		#Random effects are intercepts and slopes for age
	mu[,1,1]<-oo$mu_int[,1]
	mu[,1,2]<-oo$mu_int[,2]
	mu[,2,1]<-oo$mu_slope[,1]
	mu[,2,2]<-oo$mu_slope[,2]

	#expand beta
	beta_exp <- matrix(NA,P,1) #beta is no longer beta_k?? No longer 2 columns, just one?? ##yes, just one beta, assumed constant over classes. (separate beta was not identifiable)
	beta_exp[,1]<-oo$beta[,1] #beta is no longer beta_k??

	#expand sigma_res
	sigma_res_exp<-rep(NA,P)
	sigma_res_exp[]<-oo$sigma_res

	#expand gamma
	gamma_RC_exp <- matrix(NA,P,dim(oo$gamma_RC)[2])
	for(d in 1:dim(oo$gamma_RC)[2])
		gamma_RC_exp[,d]<-oo$gamma_RC[,d]
	
	nu_BX_exp <-
	omega_SURG_exp <- NULL

	if(IOP_BX){
		nu_BX_exp <- matrix(NA,P,dim(oo$nu_BX)[2])
		for(d in 1:dim(oo$nu_BX)[2])
			nu_BX_exp[,d]<-oo$nu_BX[,d]
	}
	
	if(IOP_SURG){
		omega_SURG_exp <- matrix(NA,P,dim(oo$omega_SURG)[2])
		for(d in 1:dim(oo$omega_SURG)[2])
			omega_SURG_exp[,d]<-oo$omega_SURG[,d]
	}
	
	
	##############
	# !!WORKCHECK!!
	# if TRUE then recycling for nreps worked!
	if(FALSE){
		all(sigma_res_exp[1:10]==sigma_res_exp[1:10+n_post])
		all(beta_exp[1:10,]==beta_exp[1:10+n_post,])
		all(mu[1:10,,]==mu[1:10+n_post,,])
		all(gamma_RC_exp[1:10,]==gamma_RC_exp[1:10+n_post,])
		all(nu_BX_exp[1:10,]==nu_BX_exp[1:10+n_post,])
		all(omega_SURG_exp[1:10,]==omega_SURG_exp[1:10+n_post,])
	}
	##############

	#Get random candidate draws for eta (re-used for each new subject)
	eta<-rbinom(P,1,prob=p_eta_exp)
		# Our eta here is analogous to eta_hat from the main JAGS model.
	#Get n_post random draws for b_vec
	b_vec_star <- matrix(NA,P,2)

	(expand_time<-system.time({
	if(verbose) pb_sim <- txtProgressBar(min = 0, max = P, char = "=", style=3)
	for(r in 1:nreps){ # nreps = number of times we expand the particles from oo.
	for(oo_ind in 1:n_post){ #index for oo
		p <- (r-1)*n_post+oo_ind #index along the extended particles we're creating.

		cov_for_bvec_p <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2))
		cov_for_bvec_p[1,2]<-
		cov_for_bvec_p[2,1]<-oo$cov_int_slope[oo_ind]
		b_vec_star[p,]<-mvrnorm(1,mu=mu[p,,eta[p]+1], Sigma=cov_for_bvec_p)
			#note, eta[i] is *not* nessecarily the same as eta[i+n_post] due to random draws.
		# we don't store cov_for_bvec_p because it's not used after b_vec_star is generated_
		if(verbose) setTxtProgressBar(pb_sim,p)
	}}})) #~ 2 min for nreps=50, n_post=25000


	return(list( #all items in this list are length nreps * n_post.
		#You previously had "exp" suffixes, but not you're ommiting those.
		eta = eta,
		p_eta = p_eta_exp,
		sigma_res = sigma_res_exp,
		beta = beta_exp,
		mu = mu,
		gamma_RC = gamma_RC_exp,
		nu_BX = nu_BX_exp,
		omega_SURG = omega_SURG_exp,
		b_vec_star = b_vec_star
	))
}



#' Get likelihood of each particle,
#' given subj_star's data (Y, Z, X, and W),
#' and a proposed set of random effects (b-vec & eta) and
#' hyper-params (beta, gamma_RC)
#'
#' @param ps particle set (list). Output from get_particles.
#' @param psa_data_star psa data for subject of interest
#' @param bx_data_star biopsy data for subject of interest
#' @param ns_BX_star splines for subject of interest
#' @param ns_BX_star splines for subject of interest
get_likelihood<-function(ps, psa_data_star, bx_data_star, ns_BX_star, ns_SURG_star, verbose=getOption('verbose')){

	P<-length(ps$sigma_res)

	#Setup subject data
	Y_star <- psa_data_star$log_psa
	X_star <- psa_data_star$vol_std #What if this is also a risk factor?? # we do not expect prostate volume to be a risk factor for prostate cancer. even increases in volume due to tumor growth are negligible in comparison to measurement error in prostate colume
	Z_star <- cbind(1, psa_data_star$age_std) #right?? ##yes, this is correct

	RC_star <- dplyr::filter(bx_data_star, bx_here==1)$rc
	BX_star <- dplyr::filter(bx_data_star, !is.na(bx_here))$bx_here
	SURG_star <- bx_data_star$surg
	prev_pos_biopsy <- dplyr::filter(bx_data_star, !is.na(bx_here))$prev_G7 #do I have this right, or do I need to exclude periods where bx_here = NA?? Does bx_here=NA mean they've exited the study (reclassified)? 
	##bx_here=NA does mean that the pt has reclassified. prev_G7=1 only in intervals with bx_here=1 & rc=1 (that is, the interval where reclassification occurred) or later intervals where bx_here=NA


	#Only take the covariates for which we have actual biopsies, otherwise R will do some bad recycling later on.
	couldve_had_biopsy_star <- !(is.na(bx_data_star$bx_here))
	did_have_biopsy_star <- couldve_had_biopsy_star & bx_data_star$bx_here

	V_RC_star <- dplyr::mutate(bx_data_star, intercept=1) %>%
		dplyr::select(intercept, 
			contains("rc_time_ns"), 
			contains("rc_date_ns"), 
			rc_age_std ) %>%
		dplyr::filter(did_have_biopsy_star)

	d_V_RC<-dim(V_RC_star)[2] #Note, these don't include eta yet
	d_Z<-dim(Z_star)[2]

	if(IOP_BX){
		U_BX_star <- bx_data_star %>% 
			dplyr::filter(couldve_had_biopsy_star) %>%#covariates just for BX
			dplyr::mutate(intercept=1) %>%
			dplyr::select(intercept, 
					contains("bx_time_ns"),
					contains("bx_date_ns"),
					contains("bx_age_ns"),
					contains("bx_num_prev_bx_ns") )
	}

	if(IOP_SURG){
		W_SURG_star<- 	bx_data_star %>%
		dplyr::mutate(intercept=1) %>%
		dplyr::select(intercept,
					contains("surg_time_ns"),
					contains("surg_date_ns"),
					contains("surg_age_ns"),
					surg_num_prev_bx_ns_std,
					prev_G7)
	}




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
			data.frame('LL_Y_all'=LL_Y_j,'p_ind'=as.factor(p_ind))%>%
			group_by(p_ind) %>%
			summarize(sum=sum(LL_Y_all))
			)$sum
	}

	

	# The remaining likelihoods are all binomial, with a
	# common form, so we can calculate them more easily with a
	# common function


	#' For either SURG, BX, or RC, get the likelihood for all visits from a subject.
	#' This returns a vector of length P, with the joint log-likelihood of all visits, under each particle
	#' It works under the assumption that for any particle p,
	#' we have a bernoulli outcome with:
	#'  logit(mean)= W %*% coeffs[p,1:dim(W)] + eta * coeffs[p,dim(W)+1] + interact_SURG * eta * W[dim(W)] * coeffs[p,dim(W)+2]
	#' Where `interact_SURG` is an indicator that adds an interaction term for the SURG regression
	get_joint_LL_measurements<-function(W, outcomes, coeffs, eta, interact_SURG=FALSE){

		if(length(outcomes)==0) return(0)

		nVisits <- length(outcomes)
		if(nVisits != dim(W)[1]) error('Outcome length does not match covariate length')

		P <- dim(coeffs)[1]
		d_W <- dim(W)[2]
		p_ind <- rep(1:P,each=nVisits) #particle index to us in dplyr grouping

		#below, exp = expanded for using in dplyr
		W_coeffs_exp <-  c(tcrossprod(as.matrix(W), coeffs[,1:d_W] ))

		eta_coeffs_dWp1 <- coeffs[,d_W+1] * eta #dWp1 indicates d_W+1
		eta_coeffs_dWp1_exp <- rep(eta_coeffs_dWp1, each = nVisits) 

		eta_coeffs_dWp2_G7_exp <- rep(0,P*nVisits)
		if(interact_SURG){
			eta_coeffs_dWp2 <- coeffs[,d_W+2] * eta
			eta_coeffs_dWp2_G7_mat <- tcrossprod(W[,d_W],eta_coeffs_dWp2) #each column corresponds to 1 particle.
			eta_coeffs_dWp2_G7_exp <- c(eta_coeffs_dWp2_G7_mat) #expand by stacking columns together
		}

		logit_p_exp <-  W_coeffs_exp + eta_coeffs_dWp1_exp + eta_coeffs_dWp2_G7_exp

		p_exp<-c(t(invLogit(logit_p_exp)))
		outcomes_exp<-rep(outcomes,times=P)		

		LL_j <- log(dbinom(x=outcomes_exp,size=1,prob=p_exp))
		LL <- ( #logLik of all visits
			data.frame(LL_all=c(t(LL_j)),ind=p_ind) %>%
			group_by(ind) %>%
			summarize(sum=sum(LL_all))
			)$sum

		return(LL)
	}


	LL_RC <- get_joint_LL_measurements(
		W=V_RC_star,
		outcomes=RC_star,
		coeffs=ps$gamma_RC,
		eta=ps$eta,
		interact_SURG=FALSE)


	LL_BX <-LL_SURG <- 0 #zero if not included in model
	if(IOP_BX){
		LL_BX <- get_joint_LL_measurements(
		W=U_BX_star,
		outcomes=BX_star,
		coeffs=ps$nu_BX,
		eta=ps$eta,
		interact_SURG=FALSE)
		
	}
	if(IOP_SURG){
		LL_SURG <- get_joint_LL_measurements(
			W=W_SURG_star,
			outcomes=SURG_star,
			coeffs=ps$omega_SURG,
			eta=ps$eta,
			interact_SURG=TRUE)
	}

	W <- exp(LL_Y + LL_BX + LL_RC + LL_SURG)


	return(W)

}






##############################
# Estimate posteriors for each subject.

# For each subject, calcualte and store:
	# for Importance sampling
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
#' @param data_star a list of dataframes for subject star
posterior_star<-function( data_star, ps, runifs, rej_const=NULL,	e_ss_threshold=800,	n_draws_init = min(length(ps[[1]]), e_ss_threshold*2)){

	P <- length(ps[[1]])
	
	if(e_ss_threshold > P) stop("Effective sample size exceeds the number of particles")

	#Get seq starting at 0, ending at P, and increasing by factors of 2, from n_draws_init
	l2P <- log(P,base=2)
	l2n1 <- log(n_draws_init,base=2)
	log_breaks <- unique(c(-Inf,seq(from=l2n1,to=l2P,by=1),l2P))
	breaks <- round(2^log_breaks)


	#Function to select a subset (inds) of the particle set
	#Particle index is always stored in the first dimension
	select_first_dim <- function(x, inds){

		ldx<-length(dim(x))

		if(ldx<=1) return(x[inds]) #for vectors, ldx = 0.
		if(ldx==2) return(x[inds,])
		if(ldx==3) return(x[inds,,])
		if(ldx==4) return(x[inds,,,]) #levels beyond here aren't necessary. No array is that large for us.
		if(ldx==5) return(x[inds,,,,])
	}

	#Calculate likelihood of particles.
	#If effective sample size not met, double size of particle set and recalculate.
	likelihood <- c()
	for(i in 2:(length(breaks))){

		inds_i <- (breaks[i-1]+1):breaks[i]

		psi <- lapply(ps,function(x){
			select_first_dim(x,inds=inds_i)
		})

		likelihood_i <- get_likelihood(
			ps=psi, #!! GLOBAL REFS !!
			psa_data_star=data_star$PSA,
			bx_data_star=data_star$BX
			)
		
		likelihood <- c(likelihood, likelihood_i)

		# Below, save copies without altering likelihood, so likelihood can be appended if needed in next stage.
		W <- likelihood / sum(likelihood) 
		effective_ss <- 1/crossprod(W)
		
		last_ind_used <- breaks[i]

		if(effective_ss >= e_ss_threshold) break

	}

	if(effective_ss < e_ss_threshold) warning('Even with full particle set, the maximum effective sample size required was not met.')

	ps_cumulative <- lapply(ps,function(x){
			select_first_dim(x,inds=1:last_ind_used)
		})

	## Importance Weighting ##
	etas_IS_star <- crossprod(W,ps_cumulative$eta)

	## Rejection Sampling ##
	#we don't know integrating constant. Does that matter if we just set this to max(likelihood)??
	if(is.null(rej_const)) rej_const <- max(likelihood)
	accept_ind <- (likelihood/rej_const) >= runifs[1:last_ind_used]
	num_accepted_star <<- sum(accept_ind)
	etas_RS_star <<- mean(ps_cumulative$eta[accept_ind])

	return(list(
		#Importance Sampling:
		W=W,
		etas_IS_star = etas_IS_star,
		effective_ss_star = effective_ss,
		particle_draws = last_ind_used,
		#Rejection Sampling:
		etas_RS_star = etas_RS_star,
		num_accepted_star = num_accepted_star
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
ps <- gen_particles(oo,nreps=nreps,verbose=TRUE)
runifs <- runif(P) #needed for rejection sampling later





###### Run loop over subjects


###### Vectors to store results
N <- max(psa_data_full$subj) #number of subjects

# store:
# * effective sample for importance weighting
# * number of proposals accepted by rejection sampling
# * mean among proposals accepted by RS
# * posterior mean estimates for IS
fit_time <-
effective_ss <-
particle_draws <-
num_accepted <-
etas_RS <-
etas_IS <- rep(NA,N)


# Some entries will be NAs at the end of this script,
# to maintain indexing system of:
	# subj_star's eta = etas_IS[star]
######



#To test for star = some random number, use
	#star <- sample(missing_etas)[1]

pb <- txtProgressBar(min = min(missing_etas), max = max(missing_etas), char = "=", style=3)
for(star in missing_etas){


	# rej_const<-0.0007 #!!?? need to get better version for. If ratio is the likelihood, then rej_const is the value at the MLE. We should know this from what we generated from? Equal to the latent value or something?
	# if(star>1) rej_const <- max(rej_const,likelihood) #this method of saving values across subjects tends to get too small.
	data_star <- list(
		PSA=filter(psa_data_full, subj == star),
		BX=filter(bx_data_full, subj==star),
		ns_BX=filter(ns_BX, subj==star),
		ns_SURG=filter(ns_SURG, subj==star)
	)

	rej_const=NULL
	fit_time[star]<-system.time({
		post_star <- posterior_star(
			data_star,ps,
			runifs=runifs,
			rej_const=rej_const,
			e_ss_threshold = 1000,
			n_draws_init = 9000)
	})['elapsed']

	etas_IS[star] <- post_star$etas_IS_star
	effective_ss[star] <- post_star$effective_ss_star
	num_accepted[star] <- post_star$num_accepted_star
	etas_RS[star] <- post_star$etas_RS_star
	particle_draws[star] <- post_star$particle_draws

	setTxtProgressBar(pb, star)

	# save.image(file=paste0(Sys.Date(),'_seed_',seed,'_online_fit_results_incremental.RData'))
}

save('seed', 'nreps','posterior_path',
	'fit_time',
	'effective_ss',
	'particle_draws',
	'num_accepted',
	'etas_RS',
	'etas_IS',
	file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_P-',P,'_online_fit_results_EffSS_done.RData'))

# load(...)

###########
#INF-OBS
# load('batches/IOP_BX-IOP_SURG/1/2015-09-14_seed_101_IOP_BX-TRUE_IOP_SURG-TRUE_P-625000_online_fit_results_done.RData')

# load('batches/IOP_BX-IOP_SURG/1/2015-09-15_seed_102_IOP_BX-TRUE_IOP_SURG-TRUE_P-625000_online_fit_results_done.RData')

# load('batches/IOP_BX-IOP_SURG/1/2015-09-14_seed_101_IOP_BX-TRUE_IOP_SURG-TRUE_P-62500_online_fit_results_done.RData')
###########


###########
#NON-INF

# load('batches/NIOP_BX-NIOP_SURG/1/2015-08-23_seed_101_IOP_BX-FALSE_IOP_SURG-FALSE_P-625000_online_fit_results_done.RData') #On AJF's computer, with nreps=10

# load('batches/NIOP_BX-NIOP_SURG/1/2015-08-28_seed_101_IOP_BX-FALSE_IOP_SURG-FALSE_P-62500_online_fit_results_done.RData') #version with nreps =1.
###########



######################################
# Performance


range(fit_time, na.rm=TRUE)
range(effective_ss, na.rm=TRUE)
hist(particle_draws, na.rm=TRUE)
hist(effective_ss, na.rm=TRUE)


library(ggplot2)

ggplot(as.data.frame(effective_ss))+ geom_histogram(aes(x=effective_ss))+scale_x_log10()
mean(effective_ss,na.rm=TRUE)
# IS appears to do better than RS

 #NAs for etas_IS are removed from eta_true_or_jags in plot()
plot(etas_IS,etas_RS,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)
plot(etas_RS,eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)

# png(file=paste0('plots/',Sys.Date(),'_',IOPs,'_agreement_MCMC_IS.png'),pointsize=18,width=520,height=520)
plot(x=etas_IS,y=eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1,ylab='Estimates from MCMC',xlab='Estimates from Importance Sampling',main='Agreement between Estimated Probabilities\nof Having Aggressive Cancer')
abline(0,1,lty=2)
dev.off()



#However, IS appears to have a slightly closer resemblance to MCMC draws
squared_errors_IS <- (etas_IS - eta_true_or_jags)^2
squared_errors_RS <- (etas_RS - eta_true_or_jags)^2
sqrt(mean(squared_errors_IS,na.rm=TRUE)) #IS root mean squared error
sqrt(mean(squared_errors_RS,na.rm=TRUE)) #RS root mean squared error
cor.na.rm<-function(x,y){
	nas<-is.na(x)
	cor(x[!nas],y[!nas])
}
cor.na.rm(etas_IS,eta_true_or_jags)

quantile(sqrt(squared_errors_IS),probs=seq(.9,1,by=.01),na.rm=TRUE)



eff_ss_error_data <- data.frame(etas_IS,eta_true_or_jags,squared_errors_IS,squared_errors_RS,effective_ss,nreps=nreps,P=P)
tail(eff_ss_error_data)

# saveRDS(eff_ss_error_data,file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_P-',P,'_effective_ss_plotframe.RData'))

# png(paste0('plots/',Sys.Date(),'_effective_SS_and_error.png'),pointsize=17,width = 520, height = 520,)
ggplot(eff_ss_error_data) + 
	geom_point(aes(y=sqrt(squared_errors_IS),x=effective_ss),alpha=.5) + 
	labs(title='Effective Sample Size v. Absolute Error\n(plotted on log scale)',x='Effective Sample Size for IS',y='Absolute Difference') +
	theme(text=element_text(size=18)) +
	scale_x_continuous(breaks=10^(1:6)/2) +
	scale_y_continuous(breaks=c(.1,.01,.001)) +
	coord_trans(y = "log10",x= "log10")
dev.off()





#FACET output from TWO runs with different values of nreps

effSS_error_nreps1<-readRDS('batches/IOP_BX-IOP_SURG/1/2015-09-14_seed_101_P-62500_effective_ss_plotframe.RData')
effSS_error_nreps10<-readRDS('batches/IOP_BX-IOP_SURG/1/2015-09-14_seed_101_P-625000_effective_ss_plotframe.RData')
effSS_error_compare<-rbind(effSS_error_nreps1,effSS_error_nreps10)


lf<-function(value,variable){
	out<-rep(NA,length(variable))
	out[variable==10] <- 'x10 draws'
	out[variable==1] <- 'x1 draw'
	out
}


#Two different ways to lay this plot out
	# overlaid with color
	# faceted
# png(paste0('plots/',Sys.Date(),'_effective_ss_facet.png'),pointsize=17,width = 520, height = 720)
# png(file=paste0('plots/',Sys.Date(),'_',IOPs,'_effective_ss_overlaid_png'),width = 620, height = 520)
ggplot(effSS_error_compare) + 
	geom_point(aes(y=sqrt(squared_errors_IS),x=effective_ss,col=as.factor(P))) + 
	labs(title='Effective Sample Size v. Absolute Error\n(on log scale)',x='Effective Sample Size for IS',y='Absolute Difference', col='# of\nParticles') +
	theme(text=element_text(size=18)) +
	scale_x_continuous(breaks=10^(1:6)/2) +
	scale_y_continuous(breaks=c(.1,.01,.001)) +
	# facet_grid(nreps~.,labeller=lf)+
	coord_trans(x= "log10",y="log10")
dev.off()







countSubj <- group_by(bx_data_full,subj)%>%
	summarise(
		nCouldBX=sum(!is.na(bx_here)),
		nDidBX=sum(bx_here==1,na.rm=TRUE),
		nCouldSURG=n(),
		everRC=sum(rc)
		#No point in checking everSURG, because if they did, we would've be fitting them here.
		)
countSubj$nPSA <- (group_by(psa_data_full,subj) %>%
	summarise(nPSA=n()))$nPSA
countSubj$RCfactor<-factor(countSubj$everRC,labels=c('Never RC','Eventually RC'))



plotData<-data.frame(etas_IS,eta_true_or_jags,effective_ss,fit_time,countSubj)
tail(plotData)

# png(paste0('plots/',Sys.Date(),'_',IOPs,'_agreement_eff_ss.png'),pointsize=17,width = 530, height = 480,)
ggplot(plotData) + geom_point(aes(x=etas_IS,y=eta_true_or_jags,color=effective_ss<500),pch=1) + labs(x='IS',y='MCMC', title='Agreement between MCMC and\n Importance Sampling (IS)') + geom_abline(intercept=0, slope=1, size = .5, lty=2)
dev.off()


qplot(fit_time,nPSA+nCouldBX,data=plotData, main='Fit time v. amount of data')

ggplot(plotData) + geom_point(aes(x=etas_IS,y=eta_true_or_jags,color=RCfactor,size=effective_ss),pch=1) + labs(color='Ever Reclassify',x='IS',y='MCMC',size='Effective Sample Size', title='Agreement between MCMC and\n Importance Sampling (IS)') + geom_abline(intercept=0, slope=1, size = .5, lty=2) + scale_size(trans="log10")

# png(paste0('plots/',Sys.Date(),'_',IOPs,'_agreement.png'),pointsize=17,width = 530, height = 480,)
ggplot(plotData) + geom_point(aes(x=etas_IS,y=eta_true_or_jags,color=RCfactor),pch=1) + labs(color='Whether\npatients\nreclassify (RC)',x='Importance Sampling',y='Markov chain Monte Carlo (via JAGS)',title='Agreement of Model Fitting Methods') + geom_abline(intercept=0, slope=1, size = .5, lty=2) + theme(text = element_text(size=15))
dev.off()
# deva('agreement_with_size_color')

ggplot(plotData) + geom_boxplot(aes(x=as.factor(nDidBX),y=etas_IS-eta_true_or_jags)) + facet_grid(~RCfactor) + labs(x='Number of biopsies',y='IS - MCMC')
# deva('number_of_BX_distribution')

boxplot((etas_IS-eta_true_or_jags)~nDidBX,data=plotData[plotData$everRC==0,],xlab='Reclassify Ever',ylab='IS - MCMC', main='Error in predictions, by RC')


error <-  etas_IS - eta_true_or_jags

plot(abs(error),log(effective_ss))

plot(plotData$nCouldBX,error)
plot(plotData$nDidBX,error)
plot(plotData$nCouldSURG,error)

boxplot(error~everRC,data=plotData,xlab='Reclassify Ever',ylab='IS - MCMC', main='Error in predictions, by RC')

boxplot((etas_IS-eta_true_or_jags)~nDidBX,data=plotData[plotData$everRC==0,],xlab='Reclassify Ever',ylab='IS - MCMC', main='Error in predictions, by RC')


#OK!! It appears to be related to the RC being positive! Maybe there is a weird issue here that's somewhat due to the simulation?
#Or why else might this over-estimate for positive people??? Hm...


check<- eta_true_or_jags - etas_IS< -.2
which(check)
#check subj numbers:
# 215 229 236 246 281 296 309 366 369 416 429 431
#  455 507 545 644 646 670 672 689 695 703 709 771
#  782 833 861 865 890 919 939 950
filter(bx_data_full,subj%in%head(which(check),3))
filter(bx_data_full,subj%in%head(which(!check),3))



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
# It appears to be because, for those subjects, few proposals being accepted_ (Or small effective sample size for IS)

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
ggplot(error_df,aes(x=error,y=effective_sample_size)) + geom_point(alpha=.4)+ coord_trans(y = "log10") + geom_abline(intercept=10000,slope=0)
deva()


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
ps2 <- gen_particles(oo,nreps=nreps_v2,verbose=TRUE)
runifs2 <- runif(length(oo$p_eta)*nreps_v2)

for(star in which(subj2refit_v2)){
	
	# rej_const<-0.0007 #!!?? need to get better version for. If ratio is the likelihood, then rej_const is the value at the MLE. We should know this from what we generated from? Equal to the latent value or something?
	# if(star>1) rej_const <- max(rej_const,likelihood) #this tends to get too small.
	rej_const=NULL
	fit_time_v2<-system.time({
		post_star <- posterior_star(star,ps2,rej_const=rej_const,runifs=runifs2)
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
eta_true_or_jags2<-c(known_etas,colMeans(of2$eta_hat))


png(file=paste0('plots/',Sys.Date(),'_compare_fits.png'),height=1200,width=500,pointsize=20)
par(mfrow=c(3,1),mar=c(4,4,3,1))
plot(eta_true_or_jags,etas_IS,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='First attempt',xlab='MCMC',ylab='Importance Sampling')
legend('topleft',c('Flagged','Normal'),col=c('red','black'),pch=1)
abline(0,1,lty=3)
plot(eta_true_or_jags,etas_IS_v2,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='Refitting',xlab='MCMC',ylab='Importance Sampling 2')
abline(0,1,lty=3)
plot(eta_true_or_jags,eta_true_or_jags2,cex=.5,xlim=0:1,ylim=0:1,col=(effective_ss<eff_ss_cutoff)+1,main='References Coherence',xlab='MCMC 1',ylab='MCMC 2')
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

hist(colMeans(of$eta_hat),breaks=30,main=paste0('Distribution of mean posterior eta_hat,\nacross people (star=',star,')'),xlab='Population Distribution')
abline(v=mean(of$eta_hat),lwd=3)
abline(v=eta_true_or_jags[star],lwd=2, col='darkgreen')
abline(v=etas_IS[star],lwd=3,col='blue',lty=2)
legend('topright',c('overall pop mean','subj* online','subj* full model'),lty=c(1,2,1),lwd=c(3,2,2),col=c('black','blue','darkgreen'))

#dev.copy2pdf(file=paste0('plots/',Sys.Date(),'subj_',star,'_plot.pdf'))




length(Y_star)
length(X_star)


P





####################################
####################################
####################################
####################################
#Workspace

if(FALSE){

#Vectorized version

 #likelihood of PSA data
system.time({
Z_star_X_bvec<-tcrossprod(b_vec_star,Z_star)
beta_eta<-beta_exp[,1] + beta_exp[,2]*(eta==1)
beta_exp_eta_Xstar<-tcrossprod(beta_eta,X_star)
mu_obs_psa_exp <-c(t(Z_star_X_bvec +beta_exp_eta_Xstar))
Y_star_exp<-rep(NA,length(Y_star)*P)
Y_star_exp[]<-Y_star #cycle Y_star up
p_ind <- rep(1:(P),each=length(Y_star))
sigma_res_exp2<-rep(sigma_res_exp,each=length(Y_star))
L_Y_j <- dnorm(Y_star_exp,mean=mu_obs_psa_exp, sd=sigma_res_exp2) #the likelihood for each visit, grouped by particle.
L_Y_frame<-data.frame('L_Y_all'=L_Y_j,'p_ind'=as.factor(p_ind))%>%
	group_by(p_ind) %>%
	summarize(prod=prod(L_Y_all))
L_Y2<-L_Y_frame$prod

})


system.time({
gamma_RC_exp <- matrix(NA,n_post*nreps,dim(oo$gamma_RC)[2])
gamma_RC_exp[,1]<-oo$gamma_RC[,1]
gamma_RC_exp[,2]<-oo$gamma_RC[,2]
gamma_RC_exp[,3]<-oo$gamma_RC[,3]
gamma_RC_exp[,4]<-oo$gamma_RC[,4]


logit_p_rc<-gamma_RC_exp[,1:3] %*% t(V_RC_star) + gamma_RC_exp[,4]*eta

L_RC_j <- matrix(dbinom(x=RC_star,size=1,prob=c(invLogit(logit_p_rc))),P,length(RC_star))
L_RC_frame <- data.frame(L_RC_all=c(t(L_RC_j)),ind=rep(1:(P),each=length(RC_star))) %>%
	group_by(ind) %>%
	summarize(prod=prod(L_RC_all))
L_R2<-L_RC_frame$prod

})


#LOOP VERSION

L_Y_save<-
L_RC_save<-
W <- rep(NA,n_post*nreps)
fit_time<-system.time({  # ~ 35k / second
for(r in 1:nreps){
for(oo_ind in 1:n_post){ #index for oo
	p <- (r-1)*n_post+oo_ind #index for our particle set

	 #likelihood of PSA data
	mu_obs_psa <- Z_star %*% b_vec_star[p,]  +
		( oo$beta[oo_ind,1] + oo$beta[oo_ind,2]*(eta[p]==1) ) * X_star
	L_Y <- prod(dnorm(Y_star,mean=mu_obs_psa, sd=oo$sigma_res[oo_ind]))
	L_Y_save[p]<-L_Y

	#likelihood of reclassifications
	if(!is.null(RC_star)){
		logit_p_rc<-cbind(V_RC_star,eta[p]) %*% oo$gamma_RC[oo_ind,1:(d_V_RC+1)]
		L_R <- prod(dbinom(x=RC_star,size=1,prob=c(invLogit(logit_p_rc))))
		L_RC_save[p]<-L_R
	}else{
		L_R <- 1
	}

	W[p] <- L_Y * L_R #weights
}}})
#################

W <- W/sum(W)
etas_IS[star] <- crossprod(eta,W)






}