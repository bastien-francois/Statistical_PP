# Code containing general purpose functions
# Author: Bastien Francois, KNMI, 2023

library(dplyr)
#Function to add information such as stations id, date...
#to a data.frame containing post-processed data. 
#The use of this function is specifically design for the local data used at KNMI.
#It should not be used in other applications if the format of data are different.
add_info_to_dataframe<-function(df_pp_without_info, df_test_with_info){
  ###-----------------------------------------------------------------------------
  ###Input
  #df_pp_without_info.........Data frame containing the post-processed dataset.
  #df_test_with_info..........Data frame of the dataset that has been post-processed. It should include information
  #...........................It should contain information such as stations id, date... to add to the post-processed data.
  #Warning! df_pp_without_info and df_test_with_info should be ordered similarly (in rows) 
  #so that the added information matches the correct predictive distribution.
  ###-----------------------------------------------------------------------------
  ###Output
  #res..............Data.frame containing post-processed data and the desired information from the test dataset.
  ###-----------------------------------------------------------------------------
  
  tmp_df=data.frame("obs"=as.numeric(df_test_with_info$obs),
                    "gp_id"=df_test_with_info$gp_id,
                    "date"=df_test_with_info$date,
                    "sta_id"=df_test_with_info$sta_id,
                    "time"=df_test_with_info$time,
                    "validTimes"=df_test_with_info$validTimes,
                    "climato_q050"=df_test_with_info$climato_q050,
                    "climato_q060"=df_test_with_info$climato_q060,
                    "climato_q070"=df_test_with_info$climato_q070,
                    "climato_q080"=df_test_with_info$climato_q080,
                    "climato_q090"=df_test_with_info$climato_q090,
                    "climato_q092"=df_test_with_info$climato_q092,
                    "climato_q094"=df_test_with_info$climato_q094,
                    "climato_q096"=df_test_with_info$climato_q096,
                    "climato_q098"=df_test_with_info$climato_q098,
                    "climato_q100"=df_test_with_info$climato_q100)
  
  df_pp_without_info$validTimes=df_test_with_info$validTimes
  df_pp_without_info$sta_id=df_test_with_info$sta_id
  df_pp_without_info$gp_id=df_test_with_info$gp_id
  df_pp_without_info$obs=df_test_with_info$obs
  res=merge(tmp_df, df_pp_without_info, by=c("sta_id","gp_id","validTimes", "obs"))
  res <- res %>% arrange(sta_id, gp_id, validTimes)
  return(res)
}
