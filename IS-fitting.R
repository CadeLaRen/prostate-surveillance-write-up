


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
# qsub -N fitIS -t 1-21 -V -l mf=10G,h_vmem=10G -cwd -b y 'R CMD BATCH --no-save IS-fitting.R fitIS.Rout'
# rm fitIS.e*
# rm fitIS.o*
# tail fitIS.Rout -n 10
taskID <- as.numeric(Sys.getenv("SGE_TASK_ID"))

########################
# Packages and basic internal functions
library(dplyr)
library(MASS)
library(ggplot2)
library(tensor)
#library(splines)

invLogit <- function(x)
	return(exp(x)/(1+exp(x)))
########################


#IOP (Informative observation process)
IOP_BX <- TRUE #Informative observation process for biopsy
IOP_SURG <- TRUE #Informative observation process for surgery
leave_one_out <- FALSE
crop <- FALSE

IOPs<-paste0(
	c('N')[!IOP_BX],'IOP_BX-',
	c('N')[!IOP_SURG],'IOP_SURG')
batch_path<-paste0(
	'batches/',
	IOPs,
	'/leave_one_out_',leave_one_out,
	'/crop_',crop,
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



##############################
####### Run Functions
##############################

source('IS-functions-setup.R')

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
n_draws_init <- n_post/10

if(taskID==1){
	dynamic <- posteriors_all_subjects(
		missing_etas=missing_etas,
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