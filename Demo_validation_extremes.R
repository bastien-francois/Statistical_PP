rm(list=ls())
library(scoringRules)
############################################################
#### twCRPS with indicator weight function
############################################################
#### Adaptation of scoringRules::crps_sample() function using chaining function
#### from Allen2022:"Evaluating forecasts for high-impact events using transformed kernel scores"

#Calculate scores (twCRPS) given observations and draws from predictive distributions.
# y: vector of realized values.
# dat: vector or matrix (depending on y; see details) of simulation draws from forecast distribution.
# t: threshold value

twcrps_sample<-function(y, dat, t){
  chaining_indic<-function(x,t){ 
    return(pmax(x,t))
  }
  res=crps_sample(chaining_indic(y,t), chaining_indic(dat,t))
  return(res)
}

### Other option can be using gaussian weight function (Allen2022)
# chaining_gauss<-function(x,t){
#   return((x-t)*pnorm(x, t, 1)+dnorm(x,t,1))
# }

############################################################
### Example of applications
############################################################
### Retrieving the results of Table 4 of Lerch2017. See paper below
# https://projecteuclid.org/journals/statistical-science/volume-32/issue-1/Forecasters-Dilemma-Extreme-Events-and-Forecast-Evaluation/10.1214/16-STS588.full
set.seed(42)
sigma2<-2/3
sd_<-sqrt(1-sigma2)
mu_<-rnorm(10000, 0, sd_)

### Generate some obs. and predictive distrib.
y<-rnorm(10000, mu_, sqrt(sigma2))
n_sample=1000
perfect=unconditional=extremist=matrix(NaN, ncol=n_sample, nrow=length(y))#rnorm(10000, mu_, sqrt(sigma2))
for(i in 1:length(y)){
  perfect[i,]<-rnorm(n_sample,mu_[i], sqrt(sigma2))
  unconditional[i,]=rnorm(n_sample, 0, 1)
  extremist[i,]=rnorm(n_sample,mu_[i]+5/2, sqrt(sigma2))
}

#### Define twCRPS threshold
t_=1.64

### Results of Table 4
mean(twcrps_sample(y, perfect, t_)) #0.018
mean(twcrps_sample(y, unconditional, t_)) #0.019
mean(twcrps_sample(y, extremist, t_))  #0.573

#### Some specific cases of twcrps ###
### Choice of t is crucial. t can be station-based (e.g., quantiles of climato. ref)
### When t is very large (outside support of predictive distrib.), then twCRPS=0 
### whatever the value of y
mean(twcrps_sample(y, perfect, 100))

### if t is above support of predictive distrib. 
### and slightly above an extreme obs., then twcrps=0
twcrps_sample(5, rnorm(10000), 5.1)





############################################################
#### extremeIndex (Taillardat et al., 2022)
############################################################
#Paper: https://www.sciencedirect.com/science/article/pii/S0169207022001017
#### Verification of forecasts for extreme events
rm(list=ls())
library(extremeIndex)
library(evmix)
library(pracma)
library(eva)

### Modification of the indexfore() function 
### because of numerical instabilities. 
### Similar to the code used in the paper available in the supplementary
indexforeb<-function (score_fore, clim)
{
  stopifnot(class(clim) == "indexclim")
  testdata = function(yy) {
    exc = score_fore[which(clim$obs > yy)]
    if (clim$estim_xi == TRUE) {
      pa = evir::gpd(clim$obs, threshold = yy, method = "pwm")$par.ests[c(1,
                                                                          2)]
    }
    else { #### changes here
      pa = c(clim$xi, 1+clim$xi*yy)
    }#### end changes
    
    ##### we condition the crps to obs > yy, it happens that crps < yy and so issues in pgpd function : all the crps < yy are set to yy :
    exc[exc<yy]=yy
    U=sapply(exc, function(q, xi, mu, beta) {
      if (xi > 0.01) {
        (1 - (1 + (xi * (q - mu))/beta)^(-1/xi))
      }
      else {
        1 - exp(-(q - mu)/beta)
      }
    }, xi = pa[1], mu = yy, beta = pa[2])
    
    U = sort(U)
    n = length(U)
    k = seq_len(n)
    omega2 = 1/(12 * n) + sum((U - (2 * k - 1)/(2 * n))^2)
    return(omega2)
  }
  cvm = sapply(clim$quantiles, function(q) testdata(q))
  result = list(quantiles = clim$quantiles, index = 1 - clim$index/cvm,
                obs = clim$obs, clim = clim, score = score_fore, estim_xi = clim$estim_xi, omega=cvm)
  class(result) = "indexfore"
  return(result)
}

set.seed(42)
#### Reproducing results of Taillardat2022
gamma_=1/4
n_sample=100000
delta=rgamma(n_sample,1/gamma_,1/gamma_)
y=rexp(n_sample, delta)

### Closed form of CRPS for extremist forecaster
CRPS_extremist<-function(y_, v_, delta_){
  res=y_+2*v_/delta_*exp((-1)*(delta_*y_)/v_)-(3*v_)/(2*delta_)
  return(res)
}

### Incomplete Gamma function
incogam=function(u){
  sapply(u,function(x) incgam(x,0))
}

### Closed form of CRPS for lambda-informed forecaster
CRPS_l_informed<-function(y, l, d){
  res=y+2*l/d*(exp(-d*y)-1)+8*(1-l)/3*(64/(y+4)**3 -1) +l*l/(2*d) + 4/7*(1-l)**2 + 8/3*l*(1-l)*(1-2*d+8*d*d-32*exp(4*d)*d*d*d*incogam(4*d)) 
  return(res)
}


### Values in Table 2 in Taillardat2022
### Ideal
ref=mean(crps_exp(y, delta))
ref
ref/ref*100 #100

#### Extremist 1.1
mean(CRPS_extremist(y, 1.1, delta))/ref*100 #100.509
#### Lambda-informed 0.75
mean(CRPS_l_informed(y, 0.75, delta))/ref*100 #100.9259
#### Lambda-informed 0.5
mean(CRPS_l_informed(y, 0.5, delta))/ref*100 #103.6508
#### Extremist 1.4
mean(CRPS_extremist(y, 1.4, delta))/ref*100 #106.7813
#### Lambda-informed 0.25
mean(CRPS_l_informed(y, 0.25, delta))/ref*100 #108.1746
#### Climato
mean(crps_gpd(y, gamma_))/ref*100 #114.4973
#### Extremist 1.8 
mean(CRPS_extremist(y, 1.8, delta))/ref*100 #123.0484




##########################################
### Application of the extremeIndex ###
##########################################
### Warning: some assumptions are needed to use extremeIndex (see below)

###0. Generate some obs. y
set.seed(42)
gamma_=1/4
n_sample=100000
delta=rgamma(n_sample,shape=1/gamma_,scale=gamma_)
y=rexp(n_sample, delta)

### 1. Goodness of fit test for obs. y with Cramer-von Mises Test
#### H0: y follows GPD(0,sigma,shape)
tmp_gof=gpdCvm(y)
if(tmp_gof$p.value<=0.05){
  print("Rejection, obs. is not GPD at 5% -> STOP: index cannot be used")
}else{
  print("Obs. is GPD at 5% conf. level, -> index can be used")
}

### 2. Estimation of shape parameter xi of y. 
###!Warning: xi should be >0 so that the package works
### also, sigma is assumed to be equal to 1 in the package: To potentially modify
### Reminder: Threshold stability: if Y follows GPD(0, sigma, xi), then Y-u|Y>u follows GPD(0, sigma+xi*u)
tmp_xi=0.25#otherwise, do evir::gpd(y, threshold = 0, method = "ml")$par.ests[1]
tmp_scale=1#otherwise, do evir::gpd(y, threshold = 0, method = "ml")$par.ests[2]
if(tmp_xi<0){
  print("Warning! xi<0, package cannot be used.")
}else{
  print("xi>0, package can be used.")
}

#### 3. Defining the thresholds u of interest (as in Taillardat's paper):
probs_thresh=c(0.5,0.8,0.9,0.95,0.975,0.99,0.995,0.999,0.9995,0.9999,0.99995)
thresh=evd::qgpd(probs_thresh,loc=0,scale=tmp_scale,shape=tmp_xi)
### In practice, thresholds can be determined using parameter stability plots
#with choosethres
#see Papastathopoulos et al., 2013 "Extended generalised Pareto models for tail estimation"
#Choice of u such that conf. interval of k_u contains 1.
choosethres(y,quantile(y, probs=c(0.75,0.8,0.85,0.9)))

estim_xi=rep(tmp_xi, length(thresh)) 

### 4. Compute CRPS of climatological ref (with closed form here)
crps_climato_ref=crps_gpd(y, shape=gamma_,location=0,scale=1)

### 5. Compute CRPS of the different forecasters (see Taillardat, 2022)
print("Warning, calibration of forecasters needs to be assessed before computing extremeIndex")
### Ideal 
crps_ideal=crps_exp(y, delta)
### Extremist (not calibrated)
crps_extremist1.1=CRPS_extremist(y, 1.1, delta)
crps_extremist1.4=CRPS_extremist(y, 1.4, delta)
crps_extremist1.8=CRPS_extremist(y, 1.8, delta)
### lambda_informed  (calibrated)
crps_informed0.75=CRPS_l_informed(y, 0.75, delta)
crps_informed0.5=CRPS_l_informed(y, 0.5, delta)
crps_informed0.25=CRPS_l_informed(y, 0.25, delta)

### 6. Compute the extremeIndex for the forecasters
T_ideal=rep(NaN, length(thresh))
T_climato=rep(NaN, length(thresh))
T_extremist1.1=rep(NaN, length(thresh))
T_extremist1.4=rep(NaN, length(thresh))
T_extremist1.8=rep(NaN, length(thresh))
T_informed0.75=rep(NaN, length(thresh))
T_informed0.5=rep(NaN, length(thresh))
T_informed0.25=rep(NaN, length(thresh))

for(i in 1:length(thresh)){
  #### Compute indexclim object (incl. CvM values) between 
  ####CRPS values when obs. above threshold
  #### and the fitted GPD of y above threshold 
  tmp_indexclim=indexclim(y, thresh=thresh[i], score_clim=crps_climato_ref, 
                          xi=estim_xi[i], estim_xi=FALSE)
  #### Compute the index T between CRPS of forecasters (1-ratio of CvM) 
  T_ideal[i]=indexforeb(crps_ideal,tmp_indexclim)$index
  T_climato[i]=indexforeb(crps_climato_ref,tmp_indexclim)$index
  T_extremist1.1[i]=indexforeb(crps_extremist1.1,tmp_indexclim)$index
  T_extremist1.4[i]=indexforeb(crps_extremist1.4,tmp_indexclim)$index
  T_extremist1.8[i]=indexforeb(crps_extremist1.8,tmp_indexclim)$index
  T_informed0.75[i]=indexforeb(crps_informed0.75,tmp_indexclim)$index
  T_informed0.5[i]=indexforeb(crps_informed0.5,tmp_indexclim)$index
  T_informed0.25[i]=indexforeb(crps_informed0.25,tmp_indexclim)$index
}
#### Under calibration, the higher the T value, the greater the skills for extremes


#### Plot Figure 2 in Taillardat2022
plot(thresh, T_ideal, ylim=c(0,1), col="blue", 
     type="l", ylab="Index", xlab="Threshold", lwd=2, main="Under calibration, the higher, the better")
lines(thresh,T_informed0.75, col="plum", lwd=2)
lines(thresh,T_informed0.5, col="darkviolet", lwd=2)
lines(thresh,T_informed0.25, col="indianred3", lwd=2)
lines(thresh,T_climato, col="red", lwd=2)
lines(thresh,T_extremist1.1, col="cyan", lwd=2)
lines(thresh,T_extremist1.4, col="darkgreen", lwd=2)
lines(thresh,T_extremist1.8, col="green", lwd=2)

legend("topleft", c("ideal","0.75-informed","0.5-informed",
                    "0.25-informed","climatological","1.1-extremist",
                    "1.4-extremist","1.8-extremist"),
       col=c("blue","plum","darkviolet","indianred3","red",
             "cyan","darkgreen","green"), lty=rep(1,8), lwd=rep(2,8))
### Ideal is the best among the calibrated forecasters (as expected)
### Second best is 0.75-informed (as expected)
### Extremists are not calibrated -> should not be considered.
