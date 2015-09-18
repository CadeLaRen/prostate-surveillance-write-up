
# User-specific, non-generalizable code
# rm(list=ls())
# setwd("/Users/ryc/Dropbox/inhealth/prediction-model-final/sim-data")
# setwd("/Users/ryc/GitHub/prostate_surveillance")
# setwd("/Users/aaronfisher/Dropbox/Future Projects/inHealth Prostate Screening/repo")
# setwd("/home/bst/student/afisher/inHealth_prostate")



# !! When you do the importance sampling comparison for *new data* for person with existing data, don't just do last PSA, go back to last biopsy, and treat all data from that day on (including the biopsy) as "newly acquired" data_
# !! Use age field (not age_std) to separate newly acquired from previously acquired data_
# Make a plot of "new patient" and plot of "new data on existing patient".
# Run on actual data eventually.


#load necessary packages
library("lme4")
library("bayesm")
library("R2jags")
library("dplyr")





#######################################
#######################################
#Load & filter data

pt_data<-read_csv("simulation-data/pt-data-sim.csv")
psa_data_all<-read_csv("simulation-data/psa-data-sim.csv")
data_use_all<-read_csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring
#data frames with '_all' suffix will be cropped in some cases to for "leave one out" JAGS fits.

#Crop out data for "leave-one-out" model fits.
data_use_star<- dplyr::filter(data_use_all, subj==star)
bx_inds <- which(data_use_star$bx_here==1)
last_bx_ind <- bx_inds[length(bx_inds)-1] #Take all data since the 2nd to last biopsy as "new data_"
last_age <- data_use_star$age[last_bx_ind]


if(!crop | length(last_age)==0 ) last_age <- 0 #there could be some PSA measurements before you're enrolled in this biopsy study.
psa_data <- dplyr::filter(psa_data_all, !(subj==star & age>=last_age))
data_use <- dplyr::filter(data_use_all, !(subj==star & age>=last_age))




(n<-dim(pt_data)[1]) #there are 1000 patients
#This matrix will always be fully intact, but psa and bx data may not be.




#######################################
#######################################
#Generate covariates and datavectors to be sent to JAGS





#get observed latent class
eta_data<-pt_data$obs_eta
table(eta_data) #107 in each
(n_eta_known<-sum(!is.na(eta_data))) #214


#########
# PSA outcome model

(n_obs_psa<-dim(psa_data)[1])
Y<-psa_data$log_psa
summary(Y)
subj_psa<-psa_data$subj

#covariates with random effects
Z_data<-as.matrix(cbind(rep(1,n_obs_psa), psa_data$age_std))
(d_Z<-dim(Z_data)[2])
round(apply(Z_data,2,summary),2)

#covariates with only fixed effects
X_data<-as.matrix(cbind(psa_data$vol_std))
(d_X<-dim(X_data)[2])
summary(X_data)
#########


#########
# Reclassification (RC) observation model, logistic regression


bx_data<-
n_bx<-
BX<-
subj_bx<-
U_BX_data<-
d_U_BX<- NULL

if(IOP_BX){
	bx_data<-data_use[!is.na(data_use$bx_here),] #remove patients who have already had RC observed but haven't had surgery or been censored
	(n_bx<-dim(bx_data)[1])
	BX<-as.numeric(bx_data$bx_here) #indicator of bx
	subj_bx<-bx_data$subj


	U_BX_data<-as.matrix(cbind(rep(1,n_bx), bx_data$age_std, bx_data$age_ns, bx_data$time, bx_data$time_ns, bx_data$sec_time_std, bx_data$sec_time_ns, bx_data$num_prev_bx ))
	(d_U_BX<-dim(U_BX_data)[2]) #12
	round(apply(U_BX_data,2,summary),2)
}


#########
# Reclassification (RC) outcome model, logistic regression

rc_data<-data_use[data_use$bx_here==1 & !is.na(data_use$bx_here),] #only use records where a biopsy occurred
(n_rc<-dim(rc_data)[1])
RC<-as.numeric(rc_data$rc)
subj_rc<-rc_data$subj

#covariates influencing risk of reclassification
V_RC_data<-as.matrix(cbind(rep(1,n_rc),  rc_data$age_std, rc_data$time, rc_data$time_ns, rc_data$sec_time_std ))
(d_V_RC<-dim(V_RC_data)[2])
round(apply(V_RC_data,2,summary) ,2)


#########
#Observation model for surgery (SURG), logistic regression

#this uses all records, because patients always at risk of choosing surgery

SURG<-
n_surg<-
subj_surg<-
W_SURG_data<- NULL

if(IOP_SURG){
	SURG<-as.numeric(data_use$surg)
	(n_surg<-dim(data_use)[1])
	subj_surg<-data_use$subj

	W_SURG_data<-as.matrix(cbind(rep(1,n_surg), data_use$age_std, data_use$age_ns, data_use$time, data_use$time_ns, data_use$sec_time_std, data_use$sec_time_ns, data_use$num_prev_bx_surg, data_use$prev_G7))
	(d_W_SURG<-dim(W_SURG_data)[2]) #12	 
	round(apply(W_SURG_data,2,summary) ,2)
}

d_W_SURG<- surg needs to be defined!!??

##get starting values, other functions necessary for call to JAGS





#######################################
#######################################
# Initialize parameters for JAGS
# 	(and other JAGS helper functions)

#lmer fit for initializing variance parameter in JAGS
mod_lmer<-lmer(log_psa~ vol_std + (1+ age_std |id), data=psa_data)
(var_vec <- apply(coef(mod_lmer)$id, 2, var)[1:d_Z])
(var_vec <- c(var_vec[2], var_vec[1]))

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2
jags_data<-list(
	K=K, n=n, eta_data=eta_data, n_eta_known=n_eta_known, 
	Y=Y,n_obs_psa=n_obs_psa, subj_psa=subj_psa, Z=Z_data, X=X_data, d_Z=d_Z, d_X=d_X, I_d_Z=diag(d_Z),
	RC=RC, n_rc=n_rc, subj_rc=subj_rc, V_RC=V_RC_data, d_V_RC=d_V_RC
	)
 
if(IOP_BX) jags_data <- c(jags_data, list(BX=BX, n_bx=n_bx, subj_bx=subj_bx, U_BX=U_BX_data, d_U_BX=d_U_BX))

if(IOP_SURG) jags_data <- c(jags_data, list(SURG=SURG, n_surg=n_surg, subj_surg=subj_surg, W_SURG=W_SURG_data, d_W_SURG=d_W_SURG))


#initialize model
#this is to set initial values of parameters
#note that not all "parameters" need to be initialized_ specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent variables that are not observed (here, eta)

inits <- function(){
		
	p_eta<-rbeta(1,1,1)

	eta_hat<-pt_data$rc[is.na(eta_data)]

	xi<-c(min(rlnorm(1),100), min(rlnorm(1),100))
	mu_raw<-as.matrix(cbind(rnorm(d_Z),rnorm(d_Z)))
	Tau_B_raw<-rwishart((d_Z+1),diag(d_Z)*var_vec)$W
	sigma_res<-min(rlnorm(1),1)

	beta<-rnorm(d_X)

	gamma_RC<-rnorm((d_V_RC+1), mean=0, sd=0.1) #ditto


	out <- list(p_eta=p_eta, eta_hat=eta_hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, gamma_RC=gamma_RC)

	
	if(IOP_BX){
		nu_BX<-rnorm((d_U_BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
		out$nu_BX <- nu_BX
	}
	if(IOP_SURG){
		omega_SURG<-c(rnorm((d_W_SURG+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1
		out$omega_SURG <- omega_SURG
	}

	out
}



# parameters to track
params <- c("p_eta", "eta_hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b_vec", "beta", "p_rc","gamma_RC",)

if(IOP_BX) params <- c(params, "nu_BX", "p_bx",)
if(IOP_SURG) params <- c(params, "omega_SURG", "p_surg")

#you may not need to monitor p_bx, p_rc, and p_surg. taking them out of the list should improve computing time a bit



# MCMC settings
#ni, nb, nt, and nc are now set in separate files.


model_file_IOP <- paste0(
	'model-for-jags-',
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG',
	'.txt'
	)

do_one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data,
		inits=inits,
		parameters.to.save=params,
		model.file=model_file_IOP,
		n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

	return(outj$BUGSoutput)
}










