#rm(list=ls())
#setwd("/Users/ryc/Dropbox/inhealth/prediction-model")

#import environment variable, used for running multiple chains in parallel
(SEED<-as.numeric(Sys.getenv("SGE_TASK_ID")))
#SEED<-4

#load necessary packages
#library("splines")
library("lme4")
library("bayesm")
library("R2jags")


#function to get natrual spline basis
get.ns.basis<-function(obs.data,knots){
	od.k1<- obs.data-knots[1]
	od.k1[od.k1<0]<-0
	od.k2<- obs.data-knots[2]
	od.k2[od.k2<0]<-0
	od.k3<- obs.data-knots[3]
	od.k3[od.k3<0]<-0
	return(as.vector((od.k1^3 - od.k3^3)/(knots[3]-knots[1]) - (od.k2^3 - od.k3^3)/(knots[3]-knots[2])))}

#get data
	#What are these different files??
psa.data<-read.csv("psa-data-for-prediction.csv")
pt.data<-read.csv("pt-data-for-prediction.csv")
data.use<-read.csv("data-to-use-prediction.csv")
tx.data<-read.csv("tx-data-for-prediction.csv")


#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

(n<-dim(pt.data)[1]) #896
length(unique(psa.data$id)) 
	#what are these values??
length(unique(data.use$id)) 
	#what are these values??


#get observed latent class
eta.data<-vector(length=n)
for(i in 1:n){
	if(sum(pt.data$id[i]%in%tx.data$id)==1){
		eta.data[i]<-as.numeric(tx.data$RRP_Gleason[tx.data$id==pt.data$id[i]]>=7)}
	else{eta.data[i]<-NA}}
(n_eta_known<-sum(!is.na(eta.data)))

#79+83 #162 #number eta observed

#order data based on whether eta was observed, necessary for JAGS
ordered<-order(eta.data)
eta.data<-eta.data[ordered]
pt.ordered<-pt.data[ordered,]

#need a list of subjects (1) ordered by eta observed and (2) consecutive numbering
	# I'm having trouble following this part, what do psa.data$id and data.use$id represent??
ids<-unique(pt.ordered$id)
psa.data$subj<-rep(0,dim(psa.data)[1])
for(j in 1:dim(psa.data)[1]){psa.data$subj[j]<-c(1:n)[ids==psa.data$id[j]]}
data.use$subj<-rep(0,dim(data.use)[1])
for(j in 1:dim(data.use)[1]){data.use$subj[j]<-c(1:n)[ids==data.use$id[j]]}
ids<-NULL


#latent class model
#no regression model here

#PSA model
(n_obs_psa<-dim(psa.data)[1])
	#this is > the number subjects??
Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj 
	#what is this??

#get natural spline basis for age
(knots.psa<-quantile(psa.data$age.std,p=c(0.25,0.5,0.75)))
psa.age.basis<- get.ns.basis( obs.data = psa.data$age.std, knots = knots.psa)

#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std, psa.age.basis)) 
(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2)

#covariates with only fixed effects
X.data<-as.matrix(cbind(psa.data$vol.std)) 
(d.X<-dim(X.data)[2])
summary(X.data)



#outcome model (logistic regression for reclassification)
data.use<-data.use[data.use$time>0,]
(n_obs_bx<-dim(data.use)[1])
R<-as.numeric(data.use$rc)
subj_bx<-data.use$subj

#natural spline basis for age
data.use$age.std<-(data.use$age-mean(data.use$age))/sd(data.use$age)

#covariates influencing risk of reclassification
W.RC.data<-as.matrix(cbind(rep(1,n_obs_bx),  data.use$age.std, data.use$bx.time))  
(d.W.RC<-dim(W.RC.data)[2])
round(apply(W.RC.data,2,summary) ,2)



##get starting values, other functions necessary for call to JAGS


#lmer fit for initializing parameters
#do this to get the starting value for a variance paramter in JAGS
mod.lmer<-lmer(log.psa~ vol.std + (1+ age.std + psa.age.basis|id), data=psa.data)
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


do.one<-function(seed){
set.seed(seed)	
outj<-jags(jags_data, inits=inits, parameters.to.save=params, model.file="prediction-model.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

out<-outj$BUGSoutput

for(j in 1:length(out$sims.list)){
	write.csv(out$sims.list[[j]], paste("jags-prediction-",names(out$sims.list)[j],"-",seed,".csv",sep=""))}
}

do.one(seed=SEED)
