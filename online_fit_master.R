


######################### 
##     Workflow:

# Load data - get the posterior from leaving one subject out. Also get the posterior from fitting on the entire dataset, as something to compare against.
# Generate particles or candidate draws for the posterior for each subject. (`gen_particles`)
# These particles are weighted or accepted/rejected to get an estimated posterior for each subject. 
# All new subjects share the same candidate draws, but have different weights (for importance sampling (IS)) or have different values accepted (for acceptance/rejection sampling (RS)).
# Compare results to assess performance.
#########################



# !! go over order of params (e.g. the 2nd dimension of the W's to make sure you coded it in the right order).


# Ranges from 1-21 (1 for small, 2-21 for big data frame)
# qsub -N fitIS -t 1-21 -V -l mf=10G,h_vmem=10G -cwd -b y 'R CMD BATCH --no-save online_fit_master.R fitIS.Rout'
# tail fitIS.Rout -n 10
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))

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
# of2 <- readRDS(posterior_path2)

nreps <- 1

n_post<-length(oo$p_eta)
(P<-n_post*nreps)



missing_etas <- which(is.na(pt_data_full$obs_eta)) #=1 if observed and aggressive, 0 if observed and not aggressive, or NA if not observed
# eta_true_or_jags2<-
eta_true_or_jags<-pt_data_full$obs_eta
eta_true_or_jags[missing_etas]<-colMeans(of$eta_hat_means) #indexed by subject
# eta_true_or_jags2[missing_etas]<-colMeans(of2$eta_hat_means) #indexed by subject
#sum(missing_etas)==length(colMeans(of$eta_hat))

# range(abs(eta_true_or_jags2 - eta_true_or_jags)) #largest error
# sqrt(mean((eta_true_or_jags2 - eta_true_or_jags)^2, na.rm=TRUE)) #rMSD
# cor(colMeans(of2$eta_hat_means),
# 	colMeans( of$eta_hat_means))
###############


###############
# Generate candidate particles to weight, or accept/reject.

gen_particles_space_efficient<-function(oo,nreps,verbose=TRUE){
#oo = Output from leave-one-Out JAGS object
#nreps = number of times to expand each posterior draw

	# Here, same idea as gen_particles, but store it differently.
	# Don't repeat redundant information
	# Only expand the subj-specific params
	# Store these as *arrays* with:
		#first index denoting 1:n_post
		#second index denoting 1:nreps
	# Also reshape the mu param to matrix form, but do not expand

	K<-2  # two latent classes
	n_post<-length(oo$p_eta) # number of particles from JAGS output
	
	P<-n_post*nreps # total number of particles we'll be weighting over

	mu<-array(NA,dim=c(n_post,2,K)) #No recycling here.
	mu[,1,1]<-oo$mu_int[,1]
	mu[,1,2]<-oo$mu_int[,2]
	mu[,2,1]<-oo$mu_slope[,1]
	mu[,2,2]<-oo$mu_slope[,2]
	
	out<-oo[names(oo)!='eta_hat_means'] #mostly just keep existing posterior information
	out$mu_int <- out$mu_slope <- NULL
	out$mu <- mu

	#Get random candidate draws for eta (reuse this for each new subject)
	#dim = n_post by nreps
	eta <- array( 
		rbinom(n_post*nreps,1,prob=rep(oo$p_eta,times=nreps)),dim=c(n_post,nreps)
		)
	
	
	(expand_time<-system.time({
	b_vec_star <- array(NA,dim=c(n_post,nreps,2))
	if(verbose) pb_sim <- txtProgressBar(min = 0, max = P, char = "=", style=3)
	for(r in 1:nreps){ # nreps = number of times we expand the particles from oo.
	for(oo_ind in 1:n_post){ #index for oo
		#

		cov_for_this_bvec <- diag(c(oo$sigma_int[oo_ind]^2,oo$sigma_slope[oo_ind]^2))
		cov_for_this_bvec[1,2]<-
		cov_for_this_bvec[2,1]<-oo$cov_int_slope[oo_ind]

		b_vec_star[oo_ind,r,]<-mvrnorm(1,mu=mu[oo_ind,,eta[oo_ind,r]+1], Sigma=cov_for_this_bvec)
			#note, eta[i] is *not* nessecarily the same as eta[i+n_post] due to random draws.
		# we don't store cov_for_bvec_p because it's not used after b_vec_star is generated_
		if(verbose) setTxtProgressBar(pb_sim,(r-1)*n_post+oo_ind)
	}}})) #~ 2 min for nreps=50, n_post=25000

	out$b_vec_star <- b_vec_star
	out$eta <- eta
	return(out)
}




#' Get likelihood of each particle,
#' given subj_star's data (Y, Z, X, and W),
#' and a proposed set of random effects (b-vec & eta) and
#' hyper-params (beta, gamma_RC)
#'
#' @param ps particle set (list). Output from get_particles.
#' @param psa_data_star psa data for subject of interest
#' @param bx_data_star biopsy data for subject of interest
get_likelihood_space_efficient<-function(ps, psa_data_star, bx_data_star, verbose=getOption('verbose')){

	n_post<-dim(ps$eta)[1]
	nreps<-dim(ps$eta)[2]

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
		#effect of eta has already been added in gen_particles function. Eta determines the betas that are selected.

		#Tensors are defined in this section are general of dimension:
			# 1:n_visits x 1:n_post x 1:nreps
			# tensors in `ps` have first two dimensions: 1:n_post x 1:nreps x ...

		n_visits <- length(Y_star)

		#Need final result in terms of 
		mu_obs_psa <- aperm(
			tensor(#Z_star_X_bvec - result is 1:n_post x 1:nreps x 1:n_visits
				ps$b_vec_star, Z_star,
				alongA=3, alongB=2
			) + 
			tensor( #beta_X_star
				array(rep(ps$beta, times=nreps),dim=c(n_post,nreps,1)),# beta_tensor 
				as.matrix(X_star), 
				alongA=3, alongB=2
			),c(3,1,2)) #need to get n_visits on first dim, so Y_star recycles correctly.

		sigma_obs_psa <- aperm(
			array(ps$sigma_res, dim=c(n_post,nreps,n_visits)),
			c(3,1,2)
			)

		#Is there a way we can avoid storing these two intermediate matrices?
		LL_Y_j <- log(dnorm(Y_star,mean=mu_obs_psa, sd=sigma_obs_psa))
		rm(mu_obs_psa, sigma_obs_psa)

		LL_Y<-apply(LL_Y_j,MARGIN=c(2,3),FUN=sum)
		rm(LL_Y_j)
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

		n_visits <- length(outcomes)
		if(n_visits != dim(W)[1]) error('Outcome length does not match covariate length')

		d_W <- dim(W)[2]
		n_post <- dim(eta)[1]
		nreps <- dim(eta)[2]

		#For each calculate the sum of 3 terms:
			# coeffs_W; coef_eta; and (optionally) coef_eta_interact_SURG
			# Each will have dimension 1:n_post x 1:nreps x 1:n_visits

		#below, exp = expanded for using in dplyr
		coeffs_W <-  tensor(
			array(coeffs[,1:d_W],dim=c(n_post,1,d_W)),
			as.matrix(W),
			alongA=3,
			alongB=2
			)
			#results in 1:n_post x 1 x 1:n_visits
				# same for all 1:nreps because no person-specific effects here.

		#
		eta_coeffs_dWp1 <- coeffs[,d_W+1] * array(eta,dim=c(dim(eta),1))
			#dWp1 indicates d_W+1 
			#results in 1:n_post x 1:nreps x 1
				#same for all 1:n_visits because eta is constant over time.
			
		eta_coeffs_dWp2_G7 <- 0
		if(interact_SURG){
			eta_coeffs_dWp2 <- coeffs[,d_W+2] * array(eta,dim=c(dim(eta),1))
				#n_post x nreps x 1 (same for all visits)
			eta_coeffs_dWp2_G7 <- tensor(
				eta_coeffs_dWp2,
				as.matrix(W[,d_W]),
				alongA=3,
				alongB=2
				) #n_post x nreps x n_visits
		}

		#PLACEHOLDER!!?? Try to find better way to add these arrays that uses recycling instead of creating new objects
			#read tensor code, and think about outer.
		logit_p_exp <- aperm(  #Need to put n_visits on the first dimension for use in dbinom. This will get summed out anyway.
			aperm(tensor(coeffs_W,matrix(1,nrow=nreps),2,2), c(1,3,2)) + 
			tensor(eta_coeffs_dWp1,matrix(1,nrow=n_visits),3,2) +
			eta_coeffs_dWp2_G7,
			c(3,1,2))

		p_exp<- invLogit(logit_p_exp)


		LL_j <- log(dbinom(x=outcomes,size=1,prob=invLogit(logit_p_exp)))

		LL <- apply(LL_j,MARGIN=c(2,3),sum)

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
posterior_star<-function( data_star, ps, runifs, rej_const=NULL,	e_ss_threshold=800,	n_draws_init = min(length(ps$eta), e_ss_threshold*2),get_zhenkes_approach=TRUE){

	P <- length(ps$eta)
	n_post <- dim(ps$eta)[1]
	nreps <- dim(ps$eta)[2]
	
	if(e_ss_threshold > P) stop("Effective sample size exceeds the number of particles")
	if(n_draws_init > P) stop('n_draws_init must be less than or equal to the number of particles')

	#Get seq starting at 0, ending at P, and increasing by factors of 2, from n_draws_init
	#making sure we include a midpoint at n_post ensures that we're always slicing away squares from array(1:P,n_post,nreps). This follows because all elements in seq2 can always be written as n_post * 2^z for some integer z. So elements in seq2 are all multiples of n_post.
		#proof: 2^(l2n_post+z) = 2^(l2n_post)*2^z = n_post * 2^z
		#Or if we use l2n1_big, we instead get some multiple of n_post
	

	l2P <- log(P,base=2)
	l2n1 <- log(n_draws_init,base=2)
	l2n1_big <- log(ceiling(n_draws_init/n_post)*n_post,base=2)
	l2n_post <- log(n_post,base=2)

	seq1<-c()
	seq2<-seq(from=l2n1_big,to=l2P,by=1) #if(n_draws_init > n_post)
	if(l2n1 <= l2n_post){
		seq1 <- c(seq(from=l2n1,to=l2n_post,by=1),l2n_post)
		seq2 <- seq(from=l2n_post,to=l2P,by=1)
	}

	log_breaks <- unique( c(-Inf, seq1, seq2, l2P) )
	breaks <- round(2^log_breaks)

	# Double checks -- these errors shouldn't ever be triggered, but are here as a redundancy.
	if(any(breaks > n_post & breaks %% n_post!=0)) stop('error occured')

	#Function to select a subset (inds) of the particle set
	#n_post index is always stored in the first dimension
	select_nth_dim <- function(x, n=1, inds, verbose=FALSE, dim_intact=TRUE){

	    d <- dim(x)
	    ld <- length(d)
	    if(is.null(d)) ld <- 1 #if we have a vector

	    if(n>ld) stop('n must be less than dimension length of x')

	    if (ld == 1) return(x[inds])

	    prefix <-
	    suffix <- ''

	    if(dim_intact){
	    	dOut <- d
	    	dOut[n] <- length(inds)

	    	prefix <- 'array('
	    	suffix <- ',dim=dOut)'
	    }

	    text_output <- 
	    	paste0(prefix,"x[",
	    		paste0(rep(",",n-1),collapse=''),
	    		'inds',
	    		paste0(rep(",",ld-n),collapse=''), 
	    		"]", suffix
	    	)
	    if(verbose) message(paste('returning',text_output))
	    eval(parse(text = text_output))
	}

	# (x<-array(1:100,dim=c(4,5,2)))
	# select_nth_dim(x,n=1,inds=1:2)
	# select_nth_dim(x,n=2,inds=3:4)
	# select_nth_dim(1:10,n=1,inds=3:4)

	subject_specific <- c('eta','b_vec_star')

	#Calculate likelihood of particles.
	#If effective sample size not met, double size of particle set and recalculate.
	time_IS<-system.time({
	likelihood <- matrix(NA,n_post,nreps)
	
	#First count up through firt column
	for(i in 2:(length(breaks))){

		these_inds_vector<-(breaks[i-1]+1):breaks[i]
		inds_1 <- unique(these_inds_vector %% n_post)
		inds_1[which(inds_1==0)] <- n_post
			#if we passed n_post, then store it as n_post rather than 0
		inds_2 <- unique(ceiling(these_inds_vector/n_post))


		if(max(inds_2)>1 & length(inds_1)!=n_post) stop('dimension mismatch has occured') #this should never happen, but putting in a redundancy just in case.


		psi <- list()
		for(j in 1:length(ps)){
			#subset by n_post
			part_j <- select_nth_dim(ps[[j]],n=1,inds=inds_1,dim_intact=TRUE)
			psi[[j]]<-part_j
			if(names(ps)[j] %in% subject_specific){
				#Further subset by nreps
				psi[[j]] <- select_nth_dim(part_j,n=2,inds=inds_2,dim_intact=TRUE)
			}
			names(psi)[j]<-names(ps)[j]
		}

		system.time({
		likelihood_i <- get_likelihood_space_efficient(
			ps=psi, #!! NEED TO UPDATE
			psa_data_star=data_star$PSA,
			bx_data_star=data_star$BX
			)})
		

		likelihood[inds_1,inds_2] <- likelihood_i

		# Below, save copies without altering likelihood, so likelihood can be appended if needed in next stage.
		W <- c(likelihood) / sum(likelihood,na.rm=TRUE)
		W <- W[!is.na(W)]
		effective_ss <- 1/crossprod(W)
		
		last_ind_used <- breaks[i]
		if(effective_ss >= e_ss_threshold) break

	}

	if(effective_ss < e_ss_threshold) warning('Even with full particle set, the maximum effective sample size required was not met.')

	inds_cum <- which(!is.na(likelihood),arr.ind=TRUE)
	ps_cumulative_eta <- select_nth_dim(
			ps[[j]],
			n=1,
			inds=unique(inds_cum[,1]),
			dim_intact=TRUE) %>%
			select_nth_dim(.,n=2,
				inds=unique(inds_cum[,2]),
				dim_intact=TRUE)
	
	## Importance Weighting ##
	etas_IS_star <- crossprod(W,c(ps_cumulative_eta))

	})


	## Rejection Sampling ##
	#we don't know integrating constant. Does that matter if we just set this to max(likelihood)??
	if(is.null(rej_const)) rej_const <- max(likelihood,na.rm=TRUE)
	accept_ind <- which(
		(c(likelihood[!is.na(likelihood)])/rej_const) >= runifs[1:last_ind_used]
		)
	num_accepted_star <<- length(accept_ind)
	etas_RS_star <<- mean(ps_cumulative_eta[accept_ind])

	## Zhenke's Approach ##
	n_missing_zhenke <- 
	etas_zhenke_star <- NA
	if(get_zhenkes_approach){
		psi_eta_1 <-
		psi_eta_0 <- psi
		psi_eta_1$eta[] <- 1
		psi_eta_0$eta[] <- rep(0,length(psi$eta))

		likelihood_eta_1 <- get_likelihood_space_efficient( #denote as l1
				ps=psi_eta_1, #!! GLOBAL REFS !!
				psa_data_star=data_star$PSA,
				bx_data_star=data_star$BX
				) 
		likelihood_eta_0 <- get_likelihood_space_efficient( #denote as l0
				ps=psi_eta_0, #!! GLOBAL REFS !!
				psa_data_star=data_star$PSA,
				bx_data_star=data_star$BX
				)

		#Want to calculate 
			# l1 * p_eta / (l0 *(1-p_eta) + l1 *p_eta)
		# and then take the mean
		# Can also work with logs, but the calculations work out to be the same
		# at least with the examples you tried in R.
		W_zhenke <- likelihood_eta_1 * psi$p_eta / 
			(likelihood_eta_1 * psi$p_eta + likelihood_eta_0 * (1-psi$p_eta))

		n_missing_zhenke <- sum(is.na(W_zhenke)) #for some particles, we end up with 0 / 0, because likelihood_eta_1 & likelihood_eta_2 both end up
		etas_zhenke_star <- mean(W_zhenke,na.rm=TRUE)
	}

	return(list(
		#Importance Sampling:
		W=W,
		etas_IS_star = etas_IS_star,
		time_IS= time_IS,
		effective_ss_star = effective_ss,
		particle_draws = last_ind_used,
		#Rejection Sampling:
		etas_RS_star = etas_RS_star,
		num_accepted_star = num_accepted_star,
		#Zhenke Sampling
		etas_zhenke_star= etas_zhenke_star,
		ns_missing_zhenke = n_missing_zhenke
	))
}








###### Run loop over subjects
posteriors_all_subjects <- function(missing_etas,
	ps, runifs, e_ss_threshold, n_draws_init, get_zhenkes_approach=FALSE,
	psa_data_full=psa_data_full,
	bx_data_full=bx_data_full,
	verbose=TRUE
	){


N <- length(missing_etas)

namesOut <-c(
	'subj',
	'time_fit_total', #time spent on subject star
	'time_IS', #time spent total doing IS (possibly dynamically)
	'effective_ss', #effective sample for importance weighting
	'particle_draws', #number of proposals used
	'num_accepted', #number of proposals accepted by rejection sampling
	'etas_RS', #mean among proposals accepted by RS
	'etas_IS', #posterior mean estimates for IS
	'etas_zhenke',
	'ns_missing_zhenke')#how many times does this method get wonky
out<-data.frame(matrix(NA,N,length(namesOut)))
names(out)<-namesOut

time_fit_total<-rep(NA,N)


if(verbose) pb <- txtProgressBar(min = 1, max = length(missing_etas), char = "=", style=3)
for(i in 1:length(missing_etas)){

	#To test for star = some random number, use
	#star <- sample(missing_etas)[1]


	# lowest E_SS is at star = 205
	#! note - worst patient is at star= 957
		# Changing the seed for gen_particles doesn't seem to change anything
		# Neither does increasing the # of particles to a huge number.

	# rej_const<-0.0007 #!!?? need to get better version for. If ratio is the likelihood, then rej_const is the value at the MLE. We should know this from what we generated from? Equal to the latent value or something?
	# if(star>1) rej_const <- max(rej_const,likelihood) #this method of saving values across subjects tends to get too small.
	data_star <- list(
		PSA=filter(psa_data_full, subj == missing_etas[i]),
		BX=filter(bx_data_full, subj == missing_etas[i])
	)

	rej_const=NULL
	time_fit_total<-system.time({
		post_star <- posterior_star(
			data_star,ps,
			runifs=runifs,
			rej_const=rej_const,
			e_ss_threshold = e_ss_threshold,
			n_draws_init = n_draws_init,
			get_zhenkes_approach=get_zhenkes_approach)
	})['elapsed']
	#eta_true_or_jags[missing_etas[i]]

	names(post_star)<-sapply(names(post_star),function(x){
		strsplit(x,'_star')[[1]][1]
	})
	
	vec_star<-c('subj'=missing_etas[i],
		unlist(post_star[names(post_star) %in% namesOut])
		)
	out[i,]<- vec_star[namesOut]
	out[i,'time_IS']<- post_star$time_IS['elapsed']
	out[i,'time_fit_total']<- time_fit_total
	

	if(verbose) setTxtProgressBar(pb, i)

	# save.image(file=paste0(Sys.Date(),'_seed_',seed,'_online_fit_results_incremental.RData'))
}


return(out)

}




##############################
####### Run Functions
##############################



##### Generate particles:

seed<-101
set.seed(seed)

#`ps` = particle set
# ps1 <- gen_particles(oo,nreps=nreps,verbose=TRUE)
ps2 <- gen_particles_space_efficient(oo,nreps=nreps,verbose=TRUE) #only expanding the stuff that needs expanding.
# ps1$b_vec_star <- cbind(c(ps2$b_vec_star[,,1]),c(ps2$b_vec_star[,,2]))
runifs <- runif(P) #needed for rejection sampling later


###### Vectors to store results
N <- max(psa_data_full$subj) #number of subjects
e_ss_threshold <- 1000
n_draws_init <- 50000

if(taskID==1){
	dynamic <- posteriors_all_subjects(
		missing_etas=missing_etas[2:5],
		ps=ps2,
		runifs=runifs, 
		e_ss_threshold=e_ss_threshold,
		n_draws_init=n_draws_init,
		get_zhenkes_approach=TRUE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)

	small <- posteriors_all_subjects(
		missing_etas=missing_etas,
		ps=ps2,
		runifs=runifs,
		e_ss_threshold=0,
		n_draws_init=50000,
		get_zhenkes_approach=TRUE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)


	save('seed', 'nreps','posterior_path',
		'small',
		'dynamic',
		file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_P-',P,'_online_fit_results.RData'))
}

if(taskID>1){

	cc <- cut(1:length(missing_etas),breaks=20)
	taskInds <- which(as.numeric(cc)==taskID-1)

	big <- posteriors_all_subjects(
		missing_etas=missing_etas[taskInds],
		ps=ps2,
		runifs=runifs, 
		e_ss_threshold=0,
		n_draws_init=P,
		get_zhenkes_approach=FALSE,
		psa_data_full=psa_data_full,
		bx_data_full=bx_data_full,
		verbose=TRUE
		)

	saveRDS(big,
		file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_IOP_BX-',IOP_BX,'_IOP_SURG-',IOP_SURG,'_P-',P,'_online_fit_big_taskID-',taskID,'.rds'))
}

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

if(FALSE){


range(small$time_fit_total, na.rm=TRUE)
quantile(time_fit_total,na.rm=TRUE)
range(effective_ss, na.rm=TRUE)
hist(particle_draws, na.rm=TRUE)
hist(effective_ss, na.rm=TRUE)
quantile(effective_ss,na.rm=TRUE)
quantile(effective_ss,prob=c(0,.001,.01,.05,.1,.5),na.rm=TRUE)

sum(effective_ss<e_ss_threshold,na.rm=TRUE)
mean(effective_ss<e_ss_threshold,na.rm=TRUE)

library(ggplot2)

ggplot(as.data.frame(effective_ss))+ geom_histogram(aes(x=effective_ss))+scale_x_log10()
mean(effective_ss,na.rm=TRUE)
# IS appears to do better than RS

plot(etas_IS,etas_zhenke,cex=.5,xlim=0:1,ylim=0:1)


 #NAs for etas_IS are removed from eta_true_or_jags in plot()
plot(etas_IS,etas_RS,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)
plot(etas_RS,eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1)
abline(0,1)

# png(file=paste0('plots/',Sys.Date(),'_',IOPs,'_agreement_MCMC_IS.png'),pointsize=18,width=520,height=520)
plot(x=etas_IS,y=eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1,ylab='Estimates from MCMC',xlab='Estimates from Importance Sampling',main='Agreement between Estimated Probabilities\nof Having Aggressive Cancer')
abline(0,1,lty=2)
dev.off()


plot(x=etas_zhenke,y=eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1,ylab='Estimates from MCMC',xlab='Estimates from Importance Sampling',main='Agreement between Estimated Probabilities\nof Having Aggressive Cancer')
abline(0,1,lty=2)


plot(x=etas_IS_big,y=eta_true_or_jags,cex=.5,xlim=0:1,ylim=0:1,ylab='Estimates from MCMC',xlab='Estimates from Importance Sampling',main='Agreement between Estimated Probabilities\nof Having Aggressive Cancer')


#However, IS appears to have a slightly closer resemblance to MCMC draws
squared_errors_IS <- (dynamic$etas_IS - eta_true_or_jags)^2
squared_errors_IS_big <- (etas_IS_big - eta_true_or_jags)^2
squared_errors_RS <- (etas_RS - eta_true_or_jags)^2
sqrt(mean(squared_errors_IS,na.rm=TRUE)) #IS root mean squared error
sqrt(mean(squared_errors_RS,na.rm=TRUE)) #RS root mean squared error
cor.na.rm<-function(x,y){
	nas<-is.na(x)
	cor(x[!nas],y[!nas])
}
cor.na.rm(etas_IS,eta_true_or_jags)

cbind(quantile(sqrt(squared_errors_IS),probs=seq(.9,1,by=.01),na.rm=TRUE),
quantile(sqrt(squared_errors_IS_big),probs=seq(.9,1,by=.01),na.rm=TRUE))



eff_ss_error_data <- data.frame(etas_IS,eta_true_or_jags,squared_errors_IS,squared_errors_IS_big,effective_ss_big,squared_errors_RS,effective_ss,particle_draws,time_fit_total,nreps=nreps,P=P)
tail(eff_ss_error_data)

# saveRDS(eff_ss_error_data,file=paste0(batch_path,Sys.Date(),'_seed_',seed,'_P-',P,'_effective_ss_plotframe.RData'))

# png(paste0('plots/',Sys.Date(),'_effective_SS_and_error.png'),pointsize=17,width = 520, height = 520,)
ggplot(eff_ss_error_data) + 
	geom_point(aes(y=sqrt(squared_errors_IS_big),x=effective_ss_big, color=particle_draws)) + 
	labs(title='Effective Sample Size v. Absolute Error',x='Effective Sample Size for IS',y='Absolute Difference (log spacing)', color='Number of\nparticle\ndraws\n(log spacing)\n') +
	theme(text=element_text(size=18)) +
	# scale_x_log10(breaks=c(10^(1:6))) +
	scale_y_log10(breaks=c(.1,.05,.01,.005,.001))+ 
	scale_colour_gradient(trans='log',limits=c(n_draws_init,P),breaks=c(n_draws_init*2^c(0,2,4,6),P),low='#7fcdbb', high='#081d58')+	
	geom_vline(xintercept=e_ss_threshold)

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



plotData<-data.frame(etas_IS,eta_true_or_jags,effective_ss,time_fit_total,countSubj)
tail(plotData)

# png(paste0('plots/',Sys.Date(),'_',IOPs,'_agreement_eff_ss.png'),pointsize=17,width = 530, height = 480,)
ggplot(plotData) + geom_point(aes(x=etas_IS,y=eta_true_or_jags,color=effective_ss<500),pch=1) + labs(x='IS',y='MCMC', title='Agreement between MCMC and\n Importance Sampling (IS)') + geom_abline(intercept=0, slope=1, size = .5, lty=2)
dev.off()



}

####################################
####################################
####################################
####################################
#Workspace

if(FALSE){

# Non vectorized version of above code

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
time_fit_total<-system.time({  # ~ 35k / second
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