# Code to apply the DistribForest method (Schlosser et al., 2019).
# Author: Bastien Francois, KNMI, 2023

rm(list=ls())
path_code  <- "/usr/people/francois/Documents/KNMI/SPENS/Code/"

setwd(path_code)
#Load useful functions to perform post-processing
source("helper_functions_main.R")
source("helper_functions_DistribForest.R")

#Set a random seed
set.seed(42)

#Train validation and test dataset are pre-processed to produce one RData file
#by lead time, run hour and season.
#Local NWP data used: ECMWF-IFS ensemble data (3 years) over NL.
#November 2018 to Sept. 2022.
#Resolution: 18km x 18km
#Number of members: 51
#Lead time: from 6 to 240.
#Season: winter (Oct. to Mar.) and summer (Apr. to Sept.)
#Init. time: 00, 12

#Local reference data used: stations (for wind gusts)
#and upscaled gauge-adjusted radar dataset (for precipitation)

#Data are divided into training/validation dataset (3 years):
#2018-2019, 2019-2020, 2020-2021
#using 4, 8, ..., 28 days for validation
#Test set is 2021-2022 year (e.g.,for winter).

#Format of data: data.frame of size (n_timestep x features)
#containing in column observations, covariates and information
#on the data such as station id, data, validTimes etc.

rh_ <- "00"
season_ <- "winter"
lt_=6

print("############# DistribForest ########################")
#Load training and test dataset for a specific run hour, a lead time and a specific season.
setwd("/nobackup/users/francois/SPENS/TrainValTest")
load(file = paste0("df_kf_windgust_TrainValTest_RH", rh_, "_LT", lt_, "_", season_, ".RData"))
  
data_train_and_val=get(paste0("df_train_and_val"))
data_test=get(paste0("df_test"))

#Remove dataframes
rm(df_train_and_val)
rm(df_train)
rm(df_val)
rm(df_test)

dim(data_train_and_val)
dim(data_test)

#Define the names of predictors to use
names_pred=names(data_train_and_val)
#Clear the names that are not useful or should not be used
#(such as stations id, observations, ensemble members and quantiles of climatology)
names_pred=names_pred[-which(names_pred %in% c('sta_id','gp_id', 'date', 'time', 'hour', 'FX', 'obs', 
                                                'validTimes','p0', 'p03', paste0("ens", 1:51), 
                                                names_pred[grep("climato_q", names_pred)]))]
#Add back some names to the list (e.g., median of the climatology)
names_pred=c(names_pred, "climato_q050")
#Delete some covariates
names_pred=names_pred[!names_pred %in% grep(paste0("_q10", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_q25", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_q50", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_q75", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_q90", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_min", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_max", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("_kurt", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("gp", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("skew", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("wh_mean", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("IQR", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("nao", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("pev", collapse = "|"), names_pred, value = T)]
names_pred=names_pred[!names_pred %in% grep(paste0("sf", collapse = "|"), names_pred, value = T)]
names_pred=names_pred
#Add for gbex
names_pred=c(names_pred, "x10fg6_max", "climato_q100")
#Add for + consistent with Taillardat
names_pred=c(names_pred, "x10fg6_q50", "x10fg6_IQR", "x10fg6_q10", "x10fg6_q90")
length(names_pred)
names_pred

##################################################
#### Application of DistribForest (Schlosser et al., 2019)
##################################################
#Set the hyper-parameters to use
gen.trun(par=c(0),family="LO", name="0", type="left")
family_=LO0()

param_control_disttree=disttree_control(type.tree="ctree", 
                                decorrelate="none", # type of decorrelation for the empirical estimating functions (or scores)
                                teststat = "quad", # type of the test to be applied for variable selection. 
                                splitstat = "quad", # type of the test to be applied for splitpoint selection. 
                                testtype = "Univ", # specifying how to compute the distribution of the test statistic.
                                method = "L-BFGS-B",#"Nelder-Mead",#"L-BFGS-B",
                                intersplit = FALSE, # =FALSE: splits in numeric variables are simple
                                mincriterion = 0, # the value of the test statistic or 1 - p-value that must be exceeded in order to implement a split. 
                                minsplit = 20, #minimum number of observations in a node to split.
                                minbucket = 20,
                                maxdepth=999) #minimum number of observations in terminal node.

param_DistribForest=list(names_pred=names_pred,
                    family=family_,
                    nb_trees=100,
                    mtry=trunc(length(names_pred)/3),
                    ctree_ctrl=param_control_disttree)

data_train_and_val=data_train_and_val[1:100,]
data_test=data_test[1:100,]
    
#Train DistForest over the training+validation dataset
DistForest_trained=DistribForest_train_forest(data_train_and_val, param_DistribForest)

#Predict for test set
PPtest_DistribForest=DistribForest_predict_forest(DistForest_trained$forest_train, data_test, param_DistribForest, probs_to_draw=1:51/52, add_info=TRUE)
      


plot(density(as.numeric(PPtest_DistribForest$df_PP[1,paste0("ens", 1:51)])))
abline(v=as.numeric(PPtest_DistribForest$df_PP[1,"obs"]))


 
