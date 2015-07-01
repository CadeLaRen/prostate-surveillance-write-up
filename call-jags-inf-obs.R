#rm(list=ls())
#setwd("/Users/ryc/Dropbox/inhealth/prediction-model-final/sim-data")
#setwd("/Users/ryc/GitHub/prostate_surveillance")
# setwd("/Users/aaronfisher/Dropbox/Future Projects/inHealth Prostate Screening/repo")
# setwd("/home/bst/student/afisher/inHealth_prostate")

#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
if(is.na(SEED)) SEED <- 0


#load necessary packages
library("lme4")
library("bayesm")
library("R2jags")
library("splines")

#get data
psa.data<-read.csv("simulation-data/psa-data-sim.csv")
pt.data<-read.csv("simulation-data/pt-data-sim.csv")

data.use<-read.csv("simulation-data/bx-data-sim.csv")
#this contains one record per annual interval for each patient until surgery or censoring



#function to get natrual spline basis
get.ns.basis<-function(obs.data,knots){
#	knots<-quantile(obs.data,p=c(0.25,0.5,0.75))
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}


#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

(n<-dim(pt.data)[1]) #1000

#get observed latent class
eta.data<-pt.data$obs.eta
table(eta.data) #107 in each
(n_eta_known<-sum(!is.na(eta.data))) #214



#PSA model
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


###bx data

#observation model
bx.data<-data.use[!is.na(data.use$bx.here),] #remove patients who have already had RC observed but haven't had surgery or been censored
(n_bx<-dim(bx.data)[1])
BX<-as.numeric(bx.data$bx.here) #indicator of bx
subj_bx<-bx.data$subj


W.BX.data<-as.matrix(cbind(rep(1,n_bx), bx.data$age.std, bx.data$age.ns, ns(bx.data$time,4), bx.data$num.prev.bx, ns(bx.data$sec.time.std,4)  ))
(d.W.BX<-dim(W.BX.data)[2]) #12
round(apply(W.BX.data,2,summary),2)



#outcome model (logistic regression for reclassification)
rc.data<-data.use[data.use$bx.here==1 & !is.na(data.use$bx.here),] #only use records where a biopsy occurred
(n_rc<-dim(rc.data)[1])
RC<-as.numeric(rc.data$rc)
subj_rc<-rc.data$subj

#covariates influencing risk of reclassification
W.RC.data<-as.matrix(cbind(rep(1,n_rc),  rc.data$age.std, rc.data$time, rc.data$time.ns, rc.data$sec.time.std ))  
(d.W.RC<-dim(W.RC.data)[2])
round(apply(W.RC.data,2,summary) ,2)



#logistic regression for RRP
#this uses all records, because patients always at risk of choosing surgery 
RRP<-as.numeric(data.use$rrp)
(n_rrp<-dim(data.use)[1])
subj_rrp<-data.use$subj


W.RRP.data<-as.matrix(cbind(rep(1,n_rrp), data.use$age.std, data.use$age.ns, ns(data.use$time,4), ns(data.use$sec.time.std,3) , data.use$num.prev.bx.rrp, data.use$prev.G7)) #
(d.W.RRP<-dim(W.RRP.data)[2])
round(apply(W.RRP.data,2,summary) ,2)
 

##get starting values, other functions necessary for call to JAGS

#lmer fit for initializing parameters
#do this to get the starting value for a variance paramter in JAGS
mod.lmer<-lmer(log.psa~ std.vol + (1+ age.std |id), data=psa.data)
(var_vec <- apply(coef(mod.lmer)$id, 2, var)[1:d.Z])
(var_vec <- c(var_vec[2], var_vec[1]))

#bundle data for call to JAGS
#this is observed data and constant variables that have already been assigned values (e.g. number of class K=2, number of subjects n, etc.)
K<-2
jags_data<-list(K=K, n=n, eta.data=eta.data, n_eta_known=n_eta_known, n_obs_psa=n_obs_psa, Y=Y, subj_psa=subj_psa, Z=Z.data, X=X.data, d.Z=d.Z, d.X=d.X, I_d.Z=diag(d.Z), BX=BX, n_bx=n_bx, subj_bx=subj_bx, W.BX=W.BX.data, d.W.BX=d.W.BX, RC=RC, n_rc=n_rc, subj_rc=subj_rc, W.RC=W.RC.data, d.W.RC=d.W.RC, RRP=RRP, n_rrp=n_rrp, subj_rrp=subj_rrp, W.RRP=W.RRP.data, d.W.RRP=d.W.RRP) 


#initialize model
#this is to set initial values of parameters
#note that not all "parameters" need to be initialized. specifically, don't initialize random effects, but do need to set initial values for mean and covariance of random effects
# also need to set initial values for latent variables that are not observed (here, eta)

inits <- function() {
	
p_eta<-rbeta(1,1,1)

eta.hat<-pt.data$rc[is.na(eta.data)]

xi<-rlnorm(d.Z)
mu_raw<-as.matrix(cbind(rnorm(d.Z),rnorm(d.Z)))
Tau_B_raw<-rwishart((d.Z+1),diag(d.Z)*var_vec)$W
sigma_res<-min(rlnorm(1),1)

beta<-rnorm(d.X)

gamma.BX<-rnorm((d.W.BX+1), mean=0, sd=0.1) #last coefficient is effect of eta=1
gamma.RC<-rnorm((d.W.RC+1), mean=0, sd=0.1) #ditto
gamma.RRP<-c(rnorm((d.W.RRP+2), mean=0, sd=0.01))  #here, include interaction with last prediction and eta=1

list(p_eta=p_eta, eta.hat=eta.hat, xi=xi, mu_raw=mu_raw, Tau_B_raw=Tau_B_raw, sigma_res=sigma_res, beta=beta, gamma.BX=gamma.BX, gamma.RC=gamma.RC, gamma.RRP=gamma.RRP) } 



# parameters to track
params <- c("p_eta", "eta.hat", "mu_int", "mu_slope", "sigma_int", "sigma_slope", "sigma_res", "rho_int_slope", "cov_int_slope", "b.vec", "beta", "gamma.BX", "gamma.RC", "gamma.RRP", "p_bx", "p_rc", "p_rrp")  #you may not need to monitor p_bx, p_rc, and p_rrp. taking them out of the list should improve computing time a bit

# MCMC settings
#ni <- 100; nb <- 20; nt <- 5; nc <- 1 
ni <- 50000; nb <- 25000; nt <- 20; nc <- 1 
#note, I am needing fewer sampling iterations here because I've solved mixing problems

source("model-for-jags-inf-obs.R")

do.one<-function(seed, return_R_obj=FALSE, save_output=TRUE){
set.seed(seed)	
outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="model-for-jags-inf-obs.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput

if(save_output){
	for(j in 1:length(out$sims.list)){
	write.csv(out$sims.list[[j]], paste("jags-prediction-inf-obs-", names(out$sims.list)[j],"-",seed,".csv",sep=""))}
}

if(return_R_obj) return(out)
}

do.one(seed=SEED)



# out<-do.one(seed=SEED,return_R_obj=TRUE, save_output=FALSE)
# len.sim<-length(out$sims.list$p_eta)
# saveRDS(out,file=paste0(Sys.Date(),'_posterior_full_nsim-',len.sim,'.rds'))
# str(out$sims.list)
# summary(out$sims.list$mu_slope)










