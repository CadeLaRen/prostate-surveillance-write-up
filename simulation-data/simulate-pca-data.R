### This code simulates data in order to demonstrate the joint modeling approach outlined in Coley et al (2015) as well as the importance sampling algorithm outlined in Fisher et al (2015).
### Data is generated using posterior estimates from fitting the joint model to data from the Johns Hopkins Active Surveillance cohort.

rm(list=ls())

### LOAD PACKAGES
library(MASS)


### SET SEED
set.seed(1)

### DEFINE FUNCTIONS
expit<-function(x){return(exp(x)/(1+exp(x)))}

#function to get natural spline basis with 3 knots_ (Just an alternate definition. See Ch 11 of Wakefield (2013))
get_ns_basis<-function(obs_data,knots){
	od_k1<- obs_data-knots[1]
	od_k1[od_k1<0]<-0
	od_k2<- obs_data-knots[2]
	od_k2[od_k2<0]<-0
	od_k3<- obs_data-knots[3]
	od_k3[od_k3<0]<-0
	return(as.vector((od_k1^3 - od_k3^3)/(knots[3]-knots[1]) - (od_k2^3 - od_k3^3)/(knots[3]-knots[2])))}


### DEFINE PARAMETER VALUES (similar to posterior estimates)
p_eta <- 0.22

mu_int <- c(1.36, 1.61)
mu_slope <- c(0.26,0.51)
mu_mat <- as.matrix(rbind(mu_int, mu_slope))

Sigma <- matrix(c(0.55^2, 0.04, 0.04, 0.4^2), nrow=2, ncol=2)
sigma_res <- 0.3

beta <- c(0.31)

nu_bx<-c(0, 0.5, -0.15, 1, -0.1, 0.7, -0.2, 0.2 , -0.5)
#int, age, age_ns, time, time_ns, sec_time, sec_time,ns, # previous biopsies, eta

gam_rc<-c(-2, 0.5, -0.2, -0.1, 0.25, 2)
#int, age, time, time_ns, sec_time, eta

omega_surg<-c(-4, -0.4, -1, 1, -0.1,  0.8, -0.3, -0.2, 1.5, 0.6, 2.5)
#int, age, age_ns, time, time_ns, sec_time, sec_time_ns, # previous biopsies, previous reclassification (Gleason >=7), eta, interaction with previous G7 and eta

#from real data, for design matrices for biopsy data
mean_age_bx<-69.4
sd_age_bx<-6.5
knots_age_bx<- c(-0.5, 0.1, 0.7)

knots_time_bx<- c(1.3, 3.2, 5.8)

mean_sec_time_bx<-4.5
sd_sec_time_bx<-4.1
knots_sec_time_bx<-c(-0.5,0.3,0.9)



### SIMULATE DATA
n <- 1000
id <- c(1:n)

ages_dx <- rnorm(n, mean=65.5, sd=5.5) #from real data
sec_time_dx <- rnorm(n, mean=1.6, sd=4.3)  #secular time, in relation to 2005


pt_data<-as.data.frame(cbind(id,ages_dx, sec_time_dx))
names(pt_data) <- c("id","age_dx","sec_time_dx")

#latent class
pt_data$eta_true <- eta_true <- rbinom(n,1,p_eta)
table(pt_data$eta_true)




### all biopsy data

times<-seq(1,10,1)

ids<-c(rep(1,10))
for(i in 2:n){ids<-c(ids,rep(i,10))}

bx_sim<-as.data.frame(cbind(ids, rep(times,n)))
names(bx_sim)<-c("id","time")

(N<-dim(bx_sim)[1])

bx_sim$eta<-rep(0,N)
for(i in 1:n){
	bx_sim$eta[bx_sim$id==i]<-pt_data$eta_true[pt_data$id==i]}

bx_sim$age<-bx_sim$sec_time<-rep(0,N)
for(i in 1:n){
	bx_sim$age[bx_sim$id==i]<-pt_data$age_dx[i] + bx_sim$time[bx_sim$id==1] + 0.5
	bx_sim$sec_time[bx_sim$id==i]<-pt_data$sec_time_dx[i] + bx_sim$time[bx_sim$id==1] + 0.5}

bx_sim$age_std<-(bx_sim$age-mean_age_bx)/sd_age_bx
bx_sim$age_ns<-get_ns_basis(bx_sim$age_std, knots_age_bx)

bx_sim$time_ns<-get_ns_basis(bx_sim$time,knots_time_bx)

bx_sim$sec_time_std<-(bx_sim$sec_time-mean_sec_time_bx)/sd_sec_time_bx
bx_sim$sec_time_ns<-get_ns_basis(bx_sim$sec_time, knots_sec_time_bx)

bx_sim$rm<-rep(0,N)


##biopsies
bx_sim$bx_here<-rep(0,N)
bx_sim$num_prev_bx<-rep(1,N)

bx_sub<-bx_sim[bx_sim$time==1,]
(n_bx<-dim(bx_sub)[1])
U_BX<-as.matrix(cbind( rep(1,n_bx), bx_sub$age_std, bx_sub$age_ns, bx_sub$time, bx_sim$time_ns,  bx_sub$sec_time_std, bx_sub$sec_time_ns, bx_sub$num_prev_bx, bx_sub$eta  ))
summary(as.vector(expit(U_BX%*%nu_bx)))

bx_sim$bx_here[bx_sim$time==1]<-rbinom(n,1,as.vector(expit(U_BX%*%nu_bx)))
#table(bx_sim$bx_here[bx_sim$time==1])

for(j in 2:10){
	for(i in 1:n){
		bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time==j]<-sum(bx_sim$bx_here[bx_sim$id==i & bx_sim$time<j]) + 1}

	bx_sub<-bx_sim[bx_sim$time==j,]
	(n_bx<-dim(bx_sub)[1])
	U_BX<-as.matrix(cbind( rep(1,n_bx), bx_sub$age_std, bx_sub$age_ns, bx_sub$time, bx_sim$time_ns,  bx_sub$sec_time_std, bx_sub$sec_time_ns, bx_sub$num_prev_bx, bx_sub$eta  ))
	
	bx_sim$bx_here[bx_sim$time==j]<-rbinom(n,1,as.vector(expit(U_BX%*%nu_bx)))}
table(bx_sim$bx_here)	

#reclassifications
bx_sim$rc<-bx_sim$prev_G7<-rep(0,N)
rc_sub<-bx_sim[bx_sim$bx_here==1,]
(n_rc<-dim(rc_sub)[1])
V_RC<-as.matrix(cbind(rep(1,n_rc), rc_sub$age_std, rc_sub$time, rc_sub$time_ns, rc_sub$sec_time_std, rc_sub$eta))

bx_sim$rc[bx_sim$bx_here==1]<-rbinom(n_rc,1,as.vector(expit(V_RC%*%gam_rc)))

for(i in 1:n){
	if(sum(bx_sim$rc[bx_sim$id==i]==1)>0){
		rc_time<-min(bx_sim$time[bx_sim$rc==1 & bx_sim$id==i])
		bx_sim$rc[bx_sim$id==i & bx_sim$time>rc_time]<-0
		bx_sim$bx_here[bx_sim$id==i & bx_sim$time>rc_time]<-0
		bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time>rc_time]<-(bx_sim$num_prev_bx[bx_sim$id==i & bx_sim$time==rc_time] + 1)
		bx_sim$prev_G7[bx_sim$id==i & bx_sim$time>=rc_time]<-1
		bx_sim$rm[bx_sim$id==i & bx_sim$time>(rc_time+2)]<-1}}


# surgery
bx_sim$surg<-rep(0,N) 
bx_sim$num_prev_bx_surg <- bx_sim$num_prev_bx + bx_sim$bx_here

W_SURG<-as.matrix(cbind(rep(1,N), bx_sim$age_std, bx_sim$age_ns, bx_sim$time, bx_sim$time_ns, bx_sim$sec_time_std, bx_sim$sec_time_ns,bx_sim$num_prev_bx_surg, bx_sim$prev_G7, bx_sim$eta, (bx_sim$prev_G7*bx_sim$eta) ))

bx_sim$surg<-rbinom(N,1,as.vector(expit(W_SURG%*%omega_surg)))

#messes up design matrices to delete columns earlier
bx_sim<-bx_sim[bx_sim$rm==0,]
(N<-dim(bx_sim)[1])


pt_data$rc<-pt_data$surg<-rep(0,n)

for(i in 1:n){
	if(sum(bx_sim$surg[bx_sim$id==i])>0){
		surg_time<-min(bx_sim$time[bx_sim$id==i & bx_sim$surg==1])
		bx_sim$rm[bx_sim$id==i & bx_sim$time>surg_time]<-1	
		pt_data$surg[pt_data$id==i]<-1}	}
table(pt_data$surg)


bx_sim<-bx_sim[bx_sim$rm==0,]
(N<-dim(bx_sim)[1])

for(i in 1:n){
	pt_data$rc[i]<-sum(bx_sim$rc[bx_sim$id==pt_data$id[i]])}
table(pt_data$rc) 

pt_data$obs_eta<-rep(NA,n)
pt_data$obs_eta[pt_data$surg==1]<-pt_data$eta_true[pt_data$surg==1]
table(pt_data$obs_eta)

for(i in 1:n){
	if(max(bx_sim$rc[bx_sim$id==i])==1){
		rc_time<-bx_sim$time[bx_sim$rc==1 & bx_sim$id==i]
		bx_sim$bx_here[bx_sim$id==i & bx_sim$time>rc_time]<-NA	} }

table(bx_sim$bx_here)



##psa data

psa_time<-seq(-1, max(bx_sim$time[bx_sim$id==1]),0.5)
psa_id<-rep(1, length(psa_time))

for(i in 2:n){
	psa_add<-seq(-1, max(bx_sim$time[bx_sim$id==i]), 0.5)
	psa_time<-c(psa_time,psa_add)
	psa_id<-c(psa_id, rep(i, length(psa_add)))}	
	
psa_data<-as.data.frame(cbind(psa_id, psa_time))
names(psa_data)<-c("id","psa_time")
(n_obs_psa<-dim(psa_data)[1])

psa_data$psa_time<-psa_data$psa_time + runif(n_obs_psa, min=-0.25, max=0.25)
psa_data$age<-vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	psa_data$age[j] <- psa_data$psa_time[j] + pt_data$age[pt_data$id==psa_data$id[j]]}
	

mean(psa_data$age) #69.53362
sd(psa_data$age) #6.624688
psa_data$age_std<-(psa_data$age-mean(psa_data$age))/sd(psa_data$age)

pt_data$vol_std<-rnorm(n,0,1)
psa_data$vol_std<-vector(length=n_obs_psa)
for(i in 1:n){
	psa_data$vol_std[psa_data$id==i] <- pt_data$vol_std[i]}

b.vec <- matrix(nrow=n, ncol=2)
for(i in 1:n){
	b.vec[i,] <- mvrnorm(n=1, mu=mu_mat[,(pt_data$eta_true[pt_data$id==i]+1)], Sigma=Sigma)}



psa_data$log_psa <- vector(length=n_obs_psa)
for(j in 1:n_obs_psa){
	lin_pred <- NULL
	lin_pred <- sum(b.vec[pt_data$id==psa_data$id[j],] * c(1, psa_data$age_std[j])) + beta[1]*psa_data$vol_std[j]
	psa_data$log_psa[j] <- rnorm(1, mean=lin_pred, sd=sigma_res)}
summary(psa_data$log_psa)




#get ordered subject variable

pt_data<-pt_data[order(pt_data$obs_eta),]
pt_data$subj<-c(1:n)
psa_data$subj<-rep(0,n_obs_psa)
for(i in 1:n){psa_data$subj[psa_data$id==pt_data$id[i]]<-pt_data$subj[i]}
bx_sim$subj<-rep(0,N)
for(i in 1:n){bx_sim$subj[bx_sim$id==pt_data$id[i]]<-pt_data$subj[i]}


write.csv(psa_data,"psa-data-sim.csv")
write.csv(pt_data,"pt-data-sim.csv")
write.csv(bx_sim,"bx-data-sim.csv")







