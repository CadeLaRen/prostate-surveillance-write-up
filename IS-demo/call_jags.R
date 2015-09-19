



SEED <- 0

crop <- FALSE # If desired
star <- 0 #!!


# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1
#ni <- 1000; nb <- 20; nt <- 5; nc <- 1
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1



#load necessary packages
library("lme4")
library("bayesm")
library("R2jags")
library("dplyr")





#######################################
#######################################
#Load & filter data

pt_data<-read.csv("simulation-data/pt-data-sim.csv")
psa_data_all<-read.csv("simulation-data/psa-data-sim.csv")
data_use_all<-read.csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring
#data frames with '_all' suffix will be cropped down in some cases (depending on `crop` and `star` options) for "leave-one-subject-out" or "leave-one-visit-out" MCMC models.

#Crop out data for "leave-one-out" model fits.
#If star==0, no data is removed_
#In this script, we use star=0 as an approximate illustration of the method that is more easily reproduced across systems.
data_use_star<- filter(data_use_all, subj==star)
bx_inds <- which(data_use_star$bx_here==1)
last_bx_ind <- bx_inds[length(bx_inds)-1] #Take all data since the 2nd to last biopsy as "new data_"
last_age <- data_use_star$age[last_bx_ind]


if(!crop | length(last_age)==0 ) last_age <- 0 #there can PSA measurements before you're enrolled in this biopsy study.
psa_data <- filter(psa_data_all, !(subj==star & age>=last_age))
data_use <- filter(data_use_all, !(subj==star & age>=last_age))




(n<-dim(pt_data)[1]) #there are 1000 patients
#This matrix will always be fully intact, but  matrices psa and bx data might not be, if star!=0.




#######################################
#######################################
#Generate covariates and data vectors to be sent to JAGS




#get observed latent class
#NA when no surgery, true cancer state not observed
eta_data<-pt_data$obs_eta
table(eta_data) #107 in each
(n_eta_known<-sum(!is.na(eta_data))) #214


#########
# PSA model

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
# Reclassification (RC) logistic regression outcome model

rc_data<-data_use[data_use$bx_here==1 & !is.na(data_use$bx_here),] #only use records where a biopsy occurred
(n_rc<-dim(rc_data)[1])
RC<-as.numeric(rc_data$rc)
subj_rc<-rc_data$subj

#covariates influencing risk of reclassification
V_RC_data<-as.matrix(cbind(rep(1,n_rc),  rc_data$age_std, rc_data$time, rc_data$time_ns, rc_data$sec_time_std ))
(d_V_RC<-dim(V_RC_data)[2])
round(apply(V_RC_data,2,summary) ,2)
#########




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

jags_data<-list(K=K, n=n, eta_data=eta_data, n_eta_known=n_eta_known,
	Y=Y, subj_psa=subj_psa, Z=Z_data, X=X_data, d_Z=d_Z, d_X=d_X, I_d_Z=diag(d_Z), n_obs_psa=n_obs_psa,
	RC=RC, n_rc=n_rc, subj_rc=subj_rc, V_RC=V_RC_data, d_V_RC=d_V_RC)


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

	gamma_RC<-rnorm((d_V_RC+1), mean=0, sd=0.1) #last coefficient is effect of eta=1

	
	out<-list(p_eta=p_eta, eta_hat=eta_hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta,
		gamma_RC=gamma_RC)
	
	out

}



# parameters to track
params <- c("p_eta", "eta_hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b_vec", "beta", "gamma_RC", "p_rc")

#Monitoring p_rc is optional. Taking it out of the list should improve computing time a bit




#######################################
#######################################
# Call JAGS

do_one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data,
		inits=inits,
		parameters.to.save=params,
		model.file='model_for_jags.txt',
		n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

	return(outj$BUGSoutput)
}



outJAGS <- do_one(seed=SEED) #Output from JAGS



#######################################
#######################################
# Create a smaller filesize version of the output (outJAGS) to save for future use.



#Vector of parameters to ignore and not save
params2ignore<-c('b_vec','rho_int_slope','eta_hat','p_rc')


outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore] #New sims.list argument to add to outJAGS2
outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.array')] #Remove redundant information
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta_hat,2,mean) #rather than saving all estimates of latent variables, just store the posterior mean for each subject.
if(crop){ #if `crop`, save subj-specific posterior for subj star. Otherwise do not.
	outJAGS2_sl$b_vec_star <- outJAGS$sims.list$b_vec[,star,]
	outJAGS2_sl$eta_hat_star <- NULL
	if(star > 214)
		#subject 214 is the last subject for which we know their true latent state.
		outJAGS2_sl$eta_hat_star <- outJAGS$sims.list$eta_hat[,star - 214]
}

outJAGS2$sims.list <- outJAGS2_sl


str(outJAGS2$sims.list)

saveRDS(outJAGS2,file=paste0('posterior_SEED-',SEED,'_star-',star,'_crop-',crop,'.rds'))










