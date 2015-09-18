# Codebook

### Abbreviations and Acronyms
* PSA: prostate-specific antigen
* PCa: prostate cancer
* BX: biopsy
* RC: reclassification
* SURG: surgical removal of prostate (radical retropubic prostatectomy)
* PT: patient
* DX: diagnosis

### PT data
1 record per patient
* id: pt-specific unique identifier, assigned for data generation
* age_dx: age at PCa diagnosis
* sec_time_dx: “secular time”, calendar date of diagnosis. 0 is January 1, 2005.
* eta_true: latent state, either 0 or 1, corresponding to indolent (Gleason <=6) or aggressive (Gleason >=7) cancer (respectively); generated for all patients for data simulation but only observed on a subset with surgery (rap)
* obs_eta: eta_true is surg=1, NA otherwise
* surg: does this subject ever have surgery?
* rc: does this subject ever reclassify, i.e. have a biopsy with a Gleason score 7 or higher
* vol_std: prostate volume, mean and sd standardized
* subj: pt-specific unique identifier, assigned after data simulation to correspond to data sorting. data sorted based on obs_eta as follows: eta observed and =0, eta observed and =1, eta unobserved. (sorting by obs_eta makes posterior sampling in JAGS much easier)

### PSA data
1 record per person per PSA observation
* psa: outcome measured in continuous time. Measurements may occur before initial PCa diagnosis. 
* log_psa: log-transformed PSA used for regression
* age: age at PSA measurement
* age_std: mean and std dev-standardized age at psa (standardized within this dataset); used in X
* vol_std: prostate volume, standardized within psa dataset; used in Z


### BX Data
1 record per annual interval where pt eligible for biopsy or (after rc) surgery.
Patients eligible for BX until RC (as per active surveillance protocol); it is only possible to reclassify once.
Patients eligible for surgery prior to RC or up to 2 years after RC
Patients censored after surgery, 2 years post-RC, or 10 years
No record for initial diagnostic bx 


* eta: true state (for all patients, used to generate data)
* bx_here: biopsy occurs in this interval; NA after a pt reclassifies
* rc: grade reclassification occurs at this annual biopsy
* surg: surgery performed in this annual interval

* time, time_ns: years since dx for bx, natural spline representation of years since dx
* age, age_std, age_ns: age at bx, standardized, and ns representation
* sec_time, sec_time_std: calendar time, standardized (see above)
* num_prev_bx: number of prior bx at beginning of interval (used to predict P(BX) in this interval); everyone has one prior bx (diagnostic) at time=1 (
* num_prev_bx_surg: number of prior bx + bx in this interval (used to predict P(SURG) in this interval)
* prev_G7: bx grade of Gleason 7 or higher (i.e. reclassification) in this interval or previous intervals (used to predict P(SURG) in this interval)

* rm: indicator for removing a record from dataset (used in data generation); not used in analysis


