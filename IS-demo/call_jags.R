



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
library("splines")
library("dplyr")





#######################################
#######################################
#Load & filter data

pt.data<-read.csv("simulation-data/pt-data-sim.csv")
psa.data.all<-read.csv("simulation-data/psa-data-sim.csv")
data.use.all<-read.csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring
#data frames with '.all' suffix will be cropped down in some cases (depending on `crop` and `star` options) for "leave-one-subject-out" or "leave-one-visit-out" MCMC models.

#Crop out data for "leave-one-out" model fits.
#If star==0, no data is removed.
#In this script, we use star=0 as an approximate illustration of the method that is more easily reproduced across systems.
data.use.star<- filter(data.use.all, subj==star)
bx.inds <- which(data.use.star$bx.here==1)
last.bx.ind <- bx.inds[length(bx.inds)-1] #Take all data since the 2nd to last biopsy as "new data."
last.age <- data.use.star$age[last.bx.ind]


if(!crop | length(last.age)==0 ) last.age <- 0 #there can PSA measurements before you're enrolled in this biopsy study.
psa.data <- filter(psa.data.all, !(subj==star & age>=last.age))
data.use <- filter(data.use.all, !(subj==star & age>=last.age))




(n<-dim(pt.data)[1]) #there are 1000 patients
#This matrix will always be fully intact, but  matrices psa and bx data might not be, if star!=0.




#######################################
#######################################
#Generate covariates and data vectors to be sent to JAGS




#get observed latent class
eta.data<-pt.data$obs.eta
table(eta.data) #107 in each
(n_eta_known<-sum(!is.na(eta.data))) #214


#########
# PSA model

(n_obs_psa<-dim(psa.data)[1])
Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj

#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std))
(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2)

#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$std.vol))
(d.X<-dim(X.data)[2])
summary(X.data)
#########


#########
# Reclassification (RC) logistic regression outcome model

rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),] #only use records where a biopsy occurred
(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc)
subj_rc<-rc.data$subj

#covariates influencing risk of reclassification
V.RC.data<-as.matrix(cbind(rep(1,n_rc),  rc.data$age.std, rc.data$time, rc.data$time.ns, rc.data$sec.time.std ))
(d.V.RC<-dim(V.RC.data)[2])
round(apply(V.RC.data,2,summary) ,2)
#########




#######################################
#######################################
# Initialize parameters for JAGS
# 	(and other JAGS helper functions)

#lmer fit for initializing variance parameter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ age.std |id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z])
(var_vec <- c(var_vec[2], var_vec[1]))

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2

jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known,
	Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), n_obs_psa=n_obs_psa,
	RC=RC, n_rc=n_rc, subj_rc=subj_rc, V.RC=V.RC.data, d.V.RC=d.V.RC)


#initialize model
#this is to set initial values of parameters
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent variables that are not observed (here, eta)

inits <- function(){
		
	p_eta<-rbeta(1,1,1)

	eta.hat<-pt.data$rc[is.na(eta.data)]

	xi<-c(min(rlnorm(1),100), min(rlnorm(1),100))
	mu_raw<-as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
	Tau_B_raw<-rwishart((d.Z+1),diag(d.Z)*var_vec)$W
	sigma_res<-min(rlnorm(1),1)

	beta<-rnorm(d.X)

	gamma.RC<-rnorm((d.V.RC+1), mean=0, sd=0.1) #last coefficient is effect of eta=1

	
	out<-list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta,
		gamma.RC=gamma.RC)
	
	out

}



# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b.vec", "beta", "gamma.RC", "p_rc")

#Monitoring p_rc is optional. Taking it out of the list should improve computing time a bit




#######################################
#######################################
# Call JAGS

do.one<-function(seed){
	set.seed(seed)	
	outj<-jags(jags_data,
		inits=inits,
		parameters.to.save=params,
		model.file='model_for_jags.txt',
		n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

	return(outj$BUGSoutput)
}



outJAGS <- do.one(seed=SEED) #Output from JAGS



#######################################
#######################################
# Create a smaller filesize version of the output (outJAGS) to save for future use.



#Vector of parameters to ignore and not save
params2ignore<-c('b.vec','rho_int_slope','eta.hat','p_rc')


outJAGS2_sl <- outJAGS$sims.list[!names(outJAGS$sims.list) %in% params2ignore] #New sims.list argument to add to outJAGS2
outJAGS2 <- outJAGS[!names(outJAGS) %in% c('sims.matrix','sims.array')] #Remove redundant information
outJAGS2_sl$eta_hat_means <- apply(outJAGS$sims.list$eta.hat,2,mean) #rather than saving all estimates of latent variables, just store the posterior mean for each subject.
if(crop){ #if `crop`, save subj-specific posterior for subj star. Otherwise do not.
	outJAGS2_sl$b.vec.star <- outJAGS$sims.list$b.vec[,star,]
	outJAGS2_sl$eta.hat.star <- NULL
	if(star > 214)
		#subject 214 is the last subject for which we know their true latent state.
		outJAGS2_sl$eta.hat.star <- outJAGS$sims.list$eta.hat[,star - 214]
}

outJAGS2$sims.list <- outJAGS2_sl


str(outJAGS2$sims.list)

saveRDS(outJAGS2,file=paste0('posterior_SEED-',SEED,'_star-',star,'_crop-',crop,'.rds'))










