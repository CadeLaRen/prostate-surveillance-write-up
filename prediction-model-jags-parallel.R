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
psa.data<-read.csv("psa-data-for-prediction.csv") #this has data on psa observations for all the individuals in the analysis. there is one record per test. data includes unique pt id, date of test, total PSA

pt.data<-read.csv("pt-data-for-prediction.csv") #this is a dataset that has one record per person. variables include unique pt id, diagnosis date, and if any reclassification is observed

data.use<-read.csv("data-to-use-prediction.csv") #this is a dataset that I have set-up for the reclassification logistic regression model. it has one record per 6 month interval with a biopsy per person. includes pt id, time of biopsy, biopsy results, and bx history up until that biopsy

tx.data<-read.csv("tx-data-for-prediction.csv") #this dataset has a record for each patient who had their prostate removed, includes "true" gleason score



#Before call to JAGS, get the data into simple matrices and vectors to send to JAGS

(n<-dim(pt.data)[1]) #896
length(unique(psa.data$id)) 
	#what are these values??
	#here, I am just checking that I have 896 unique patients in my psa data
length(unique(data.use$id)) 
	#what are these values??
	#here, I am just checking that I have 896 unique patients in my biopsy data
	#I just had these checks in here bc I had done so much data tidying, I wanted to make sure I pulled the right files


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
	
	# $id is the unique id field for both of these data sets
	#here I am basically re-assigning unique ids to meet the two objectives mentioned in the above comment. 
	#consecutive numbering- I need a list of unique ids that doesn't have any skips (i.e., subj ids 1-896, not 1-2000 with 896 people somewhere among those 2000) because I index the random effects by this subject id. 
	#ordered by true cancer state- it is not so important that those with indolent cancer have subject ids lower than those with aggressive, but what is important is that we know true cancer state for subjects 1-162 and we do not know cancer state for subjects 163-896. it's important for the JAGS model, because I define eta (the latent state) differently for these two groups (for the first, eta is known data, for the second, it is a latent variable I will sample)
	#$id will refer to the unique ids that I originally got in the data from Dr. Carter's lab, so this id does have skips (e.g., there is no id=8) and this id is not ordered with respect to true cancer state. I need to save this id to refer back to all the other files I have that use this unique identifier 
	#"subj" will be this consecutive, ordered list of unique identifiers
	
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
	#this is the number of PSA observations we have, so it is >>number subjects
	#I will loop through all 1:n_obs_psa observations in the JAGS code
Y<-psa.data$log.psa
summary(Y)
subj_psa<-psa.data$subj 
	#what is this??
	#this variable has length=n_obs_psa, i.e., it is defined for every PSA observation
	#I have matched the newly created unique subject variable with the psa data (which previously just had $id), so I will use it in the JAGS code to refer to a patient's random effects and latent state

#get natural spline basis for age
(knots.psa<-quantile(psa.data$age.std,p=c(0.25,0.5,0.75)))
psa.age.basis<- get.ns.basis( obs.data = psa.data$age.std, knots = knots.psa)

#covariates with random effects
Z.data<-as.matrix(cbind(rep(1,n_obs_psa), psa.data$age.std, psa.age.basis)) 

#this is the design matrix for the random effects
#age.std is standardized age, i.e., centered at mean of all ages for PSA observations and divided by the sd of those ages
#psa.age.basis is the 3rd basis fcn, h_3(), in the write-up. it corresponds to the "curvature"

(d.Z<-dim(Z.data)[2])
round(apply(Z.data,2,summary),2) #this is just here to check that I defined things correctly, have pulled in the right data

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
