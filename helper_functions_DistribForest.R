# Code containing functions for postprocessing ensemble forecasts with DistribForest (Schlosser et al., 2019)
# Author: Bastien Francois, KNMI, 2023
# Original code from R package disttree

library(partykit)
library(disttree)
library(gamlss.dist)
library(crch)
library(gamlss)
library(gamboostLSS)
library(scoringRules)
library(gamlss.tr)
library(gamlss.cens)


# Function to train DistribForest 
DistribForest_train_forest<-function(df_train, list_hp){
    ###-----------------------------------------------------------------------------
    ###Input
    #df_train.........Train data (data.frame) including predictors and obs. information.  (n_train x (n_preds + 1) data.frame)
    #list_hp..........List of hyperparameters for the DistribForest model.
    #.................The list should at least contain the names of predictors to use (vector names_pred), 
    #.................the distribution family to fit (family_),
    #................. the number of trees in the forest (nb_trees)
    #.................the number of variables to possibly split at in each node (mtry),
    #.................the minimal node size to consider further split (minsplit).
    #.................the minimal terminal node size allowed (min_node_size)
    #See package disttree for further details.
    ###-----------------------------------------------------------------------------
    ###Output
    #res..............List containing:
    #df_train.........Train data used as inputs.
    #forest_train.....The fitted forest (distforest object)
    #list_param.......List of hyperparameters for the distforest model.
    ###------------------------------------                                                           
  
    
    # Define formula for regression
    {
        formula_char='obs ~ ens_mean'
    for(pred_ in list_hp$names_pred){
        if(pred_ !="obs" & pred_ !="ens_mean" ){
            formula_char=paste0(formula_char,  paste0(" + ", pred_))
        }
    }
    # tree and forest formula
    eval(parse(text=paste0('df_formula <- as.formula(formula_char)')))
    }

    #Define Parametric Distrib used.
    if(list_hp$family$family[1]=="LO0"){
        gen.trun(par=c(0),family="LO", name="0", type="left")
        tmp_family= LO0()
    }
    if(list_hp$family$family[1]=="GA0"){
        gen.trun(par=c(0),family="GA", name="0", type="left")
        tmp_family= GA0()
    }
    if(list_hp$family$family[1]=="NO0"){
        gen.trun(par=c(0),family="NO", name="0", type="left")
        tmp_family= NO0()
    }
    if(list_hp$family$family[1]=="BCT0"){
        gen.trun(par=c(0),family="BCT", name="0", type="left")
        tmp_family= BCT0()
    }

    print('Training...')
    print(list_hp$family$family[1])
    #Fit Distributional Forest
    fitted_DistForest <- distforest(df_formula, 
                 data = df_train, 
                 family = tmp_family, 
                 type.tree = "ctree", 
                 ntree = list_hp$nb_trees, 
                 mtry = list_hp$mtry, 
                 control = list_hp$ctree_ctrl, 
                 trace=TRUE)

  return(list(df_train=df_train, 
              forest_train=fitted_DistForest, 
              list_param=list_hp))
}



# Function to predict with forest_train from DistribForest_train_forest() 
DistribForest_predict_forest<-function(forest_train_object, df_test, list_hp, probs_to_draw=1:51/52, add_info=TRUE){
    ###-----------------------------------------------------------------------------
    ###Input
    #forest_train_object.........Trained forest object obtained from DistribForest_train_forest().
    #df_test.........Test data used as inputs.
    #list_hp..........List of hyperparameters for the DistribForest model.
    #.................The list should at least contain the names of predictors to use (vector names_pred), 
    #.................the distribution family to consider (family_),
    #................. the number of trees in the forest (nb_trees)
    #.................the number of variables to possibly split at in each node (mtry),
    #.................the minimal node size to consider further split (minsplit).
    #.................the minimal terminal node size allowed (min_node_size)
    #See package disttree for further details.
    #probs_to_draw............Probabilities for which quantiles have to be estimated.
    #add_info.................Boolean. Should information from the test dataset be
    #.........................added to the final output? =FALSE if 
    #.........................different data than those from KNMI are used (general case).
    ###-----------------------------------------------------------------------------
    ###Output
    #res..............List containing:
    #df_PP.................Data.frame containing DistribForest post-processed data in ens1,...,ens51.
    #df_test...............Test data used as input.
    #list_hp...............List of hyperparameters for the DistribForest
    ###-----------------------------------------------------------------------------
    print('Predict...')
    print(list_hp$family$family[1])
    nb_members=51
    # Compute outputs for test set
    tmp_res=matrix(NaN, ncol=nb_members, nrow=length(df_test[,1]))

    ### Prediction on test set
    #new
    #Select predictors
    predictor_df_test=(df_test[,list_hp$names_pred])
    #endnew
    #pred.hat <- predict(forest_train_object, df_test, type = "parameter")
    pred.hat <- predict(forest_train_object, predictor_df_test, type = "parameter")
    
    if(list_hp$family$family[1] %in% c("LO0", "GA0", "NO0")){
        fitted_mu <- pred.hat$mu
        fitted_sigma <- pred.hat$sigma
        for(i in 1:nrow(predictor_df_test)){
            if(list_hp$family$family[1]=="LO0"){tmp_res[i,]<-qLO0(probs_to_draw, mu=fitted_mu[i], sigma=fitted_sigma[i])}
            if(list_hp$family$family[1]=="GA0"){tmp_res[i,]<-qGA0(probs_to_draw, mu=fitted_mu[i], sigma=fitted_sigma[i])}
            if(list_hp$family$family[1]=="NO0"){tmp_res[i,]<-qNO0(probs_to_draw, mu=fitted_mu[i], sigma=fitted_sigma[i])}
        }
    }
    if(list_hp$family$family[1] %in% c("BCT0")){
        fitted_mu <- pred.hat$mu
        fitted_sigma <- pred.hat$sigma
        fitted_nu <- pred.hat$nu
        fitted_tau <- pred.hat$tau
        for(i in 1:nrow(predictor_df_test)){
            if(list_hp$family$family[1]=="BCT0"){
                tmp_res[i,]<-qBCT0(probs_to_draw, mu=fitted_mu[i], sigma=fitted_sigma[i],nu=fitted_nu[i], tau=fitted_tau[i])
            }
        }
    }
    tmp_res=as.data.frame(tmp_res)
    colnames(tmp_res)=paste0("ens", 1:length(probs_to_draw))
    if(add_info==TRUE){
        final_output<-add_info_to_dataframe(tmp_res, df_test)
    }else{
        final_output <- tmp_res
    }

    return(list(df_PP=final_output,
              df_test=df_test, 
              list_param=list_hp))
}




