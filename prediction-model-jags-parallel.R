rm(list=ls())
# setwd("/Users/ryc/Dropbox/inhealth/prediction-model")


#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
if(is.na(SEED)) SEED<-4

#load necessary packages
#library("splines")
library("lme4")
library("bayesm")
library("R2jags")



#get data
psa.data<-read.csv("simulation-data/psa-data-sim.csv") #this has data on psa observations for all the individuals in the analysis. there is one record per test. data includes unique pt id, date of test, total PSA


pt.ordered<-read.csv("simulation-data/pt-ordered-sim.csv") #this is a dataset that has one record per person. variables include unique pt id, diagnosis date, and if any reclassification is observed. data set is ordered following eta.data (below)
(n<-dim(pt.ordered)[1]) #1000


bx.data<-read.csv("simulation-data/bx-data-sim.csv") #biopsy data, one record per person per post-dx biopsy, includes pt id, time of biopsy, results


eta.data<-read.csv("simulation-data/eta-data-sim.csv") #this dataset has all the observed gleason scores (from surgery) and NA for those without surgery. data is ordered based on value (0,1,NA) and the order corresponds to the "subj" variable in all the other data sets
eta.data<-eta.data[,2]
(n_eta_known<-sum(!is.na(eta.data)))



#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

#latent class model
#no regression model here

#PSA model
(n_obs_psa<-dim(psa.data)[1])
	#this is the number of PSA observations we have, so it is >>number subjects
	#I will loop through all 1:n_obs_psa observations in the JAGS code

Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj #unique id that corresponds to the vector eta.data and will be used to index the random effects


#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$std.age, psa.data$age.basis)) 

#this is the design matrix for the random effects
#age.std is standardized age, i.e., centered at mean of all ages for PSA observations and divided by the sd of those ages
#psa.age.basis is the 3rd basis fcn, h_3(), in the write-up. it corresponds to the "curvature"

(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2) #this is just here to check that I defined things correctly, have pulled in the right data


#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$std.vol)) 
(d.X<-dim(X.data)[2])
summary(X.data)


#outcome model (logistic regression for reclassification)
(n_obs_bx<-dim(bx.data)[1])
R<-as.numeric(bx.data$rc)
#table(R) #290

subj_bx<-bx.data$subj

#covariates influencing risk of reclassification
W.RC.data<-as.matrix(cbind(rep(1,n_obs_bx),  bx.data$std.age, bx.data$bx.time))  
(d.W.RC<-dim(W.RC.data)[2])
round(apply(W.RC.data,2,summary) ,2)


##get starting values, other functions necessary for call to JAGS


#lmer fit for initializing parameters
#do this to get the starting value for a variance paramter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ std.age + age.basis|id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z])


#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2
jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), n_obs_bx=n_obs_bx, R=R, subj_bx=subj_bx, W.RC=W.RC.data, d.W.RC=d.W.RC) 


#initialize model
#this is to set initial values of parameters
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent variables that are not observed (here, eta)

inits <- function() {
	
p_eta<-rbeta(1,1,1)

eta.hat<-pt.ordered$rc[is.na(eta.data)]

mu<-as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
Tau_B<-rwishart((d.Z+1),diag(d.Z)*var_vec)$W
sigma_res<-min(rlnorm(1),3)

beta<-rnorm(d.X*2)

gamma.RC<-rnorm((d.W.RC+1),mean=0,sd=0.25)

list(p_eta=p_eta, eta.hat=eta.hat, mu=mu, Tau_B=Tau_B, sigma_res=sigma_res, beta=beta, gamma.RC=gamma.RC) } 


# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "mu_spline", "sigma_spline", "sigma_res", "rho_int_slope", "rho_int_spline", "rho_slope_spline", "cov_int_slope", "cov_int_spline", "cov_slope_spline", "b.vec", "beta", "gamma.RC") 

# MCMC settings
#ni <- 250; nb <- 50; nt <- 5; nc <- 1 
#ni <- 25000; nb <- 5000; nt <- 20; nc <- 1 #mixing usually good by here
ni <- 100000; nb <- 50000; nt <- 20; nc <- 1


source("prediction-model.R") 


#seed<-SEED


do.one<-function(seed,return_R_obj=FALSE, save_output=TRUE){
set.seed(seed)	
outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="prediction-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput

if(save_output){
	for(j in 1:length(out$sims.list)){
		write.csv(out$sims.list[[j]], paste("jags-prediction-",names(out$sims.list)[j],"-",seed,".csv",sep=""))}}

if(return_R_obj) return(out)
}

do.one(seed=SEED)



# WORKSPACE
if(FALSE){
	out<-do.one(seed=SEED,return_R_obj=TRUE, save_output=FALSE)
	#saveRDS(out,file='posterior_full.rds')
	str(out$sims.list)
	summary(out$sims.list$mu_spline)
}



