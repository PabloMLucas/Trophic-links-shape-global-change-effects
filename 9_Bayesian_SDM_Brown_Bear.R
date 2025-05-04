

############################################################################################################
#Readme:
############################################################################################################
#R code to fit bayesian models to explain brown bear distribution a funtion of biotic and/or abiotic variables
#Authors: Pablo M. Lucas
#Last update: 

#DESCRIPTION (number in parenthesis correspond with the first level of the schema for the script):
  #This script fits, validate, visualize, compare and predict the distribution (in raster maps and in KLM for  
  #visualization in Google Earth) of differents SDMs for brown bear using Bayesian statistics:
  #fits a SDM with historical distribution and historical climate for the brown bear at range scale 
  # at a gross resolution (10x10km) for all Eurasia(4); fits a SDM for brown bear with only current data at 
  #high resolution (1x1 km) with climate, land use and biotic variables (Energy variables from food species)(5);
  #fits a Bayesian hierarchical model using results for climate variables from the SDM historical as priors 
  #for climate variables and current data, land use and biotic variables use non-informative priors (6);
  #fits differents combinations of models hierarchical/non hierarechical and with different combinations of variables land use,
  #climate, biotic (7-12). Abiotic variables refers to climate and land use variables. 
  #The selected model (based in the best fit with WAIC) was the 
  #"Bayesian hierarchical model with current data using results of Bayesian model with historical data as priors" (6).

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run script:
  #/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_model_all_5km.Rdata
  #/7_Construction_of_Historical_database/presences_absences_eurasia_hist.RData
  #/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1901/bioclimatic_var_1901_cea.RData
  #/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1901/DF_bio1901_range.RData
  #/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1981/bioclimatic_var_1981_cea.RData
  #/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_model_all_5km.Rdata
  #/8_Univarable_SDM_brown_bear_comparing_biotic_variables/DF_AIC_univar_models_habitat_land_use.RData # We check the best land use variables
  #/8_Univarable_SDM_brown_bear_comparing_biotic_variables/DF_AIC_univar_models_habitat_energy2.RData # We check the best biotic variables
  #/0_Construction_of_the_Spatial_Database/data_matrix.RData #We load the data of coordinates and ID pixels
  #/4_SDM_Food_species/vec_extrap_raster.RData #We load the vector of extrapolation
  #/4_SDM_Food_species/vec_subpop_raster.RData #We load the delimitation of the subpopulations which are used as random factors in the GLMMs       
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP26.RData
  #/6_Calculation_of_Biotic_variables/merged_DF_Variables_Energy_RCP26.RData # #We load the scenario for current biotic variables 
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP1_RCP26.RData We load the vector of water
  #/0_Construction_of_the_Spatial_Database/vec_gaus_pr_sel_natural_SSP1_RCP26.Rdata # We load the vector of Frag           
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP60.RData
  #/6_Calculation_of_Biotic_variables/merged_DF_Variables_Energy_RCP60.RData
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP3_RCP70.RData
  #/0_Construction_of_the_Spatial_Database/vec_gaus_pr_sel_natural_SSP3_RCP70.Rdata
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP85.RData
  #/6_Calculation_of_Biotic_variables/merged_DF_Variables_Energy_RCP85.RData
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP5_RCP85.RData
  #/12_Construction_Database_Brown_Bear_occurrences/DF_Scen_current_all_noNA.RData      
  #/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_validation_all_5km.Rdata
  #shapefile('data_0.shp') # Brown bear data of current range download from the IUCN
  #raster("PROTECTED_AREAS_TER_COAST_POINT_WORLD_REC.rst") # Rasterized terrestrial protected areas
#Data output to other scripts:
  #No outputs to other scripts

##############################################################################################                
#Schema
##############################################################################################                
#9_Bayesian_SDM_Brown_Bear
  #9.1  Software installation and load datasets
    #9.1.1 Install required packages 
    #9.1.2 We load the data
  #9.2  Variable selection and data preparation
        #9.2.1 We drop variables highly correlated based in VIF 
        #9.2.2 We select the best four variables based on univariable GLM including the quadratic effect
        #9.2.3 We create a subset with the selected variables and the presences/absences
        #9.2.4 We are going to divide the variables to have a more homogeneous range of data among them
        #9.2.5 We are going to calculate the cuadrateic effect from each variable as in Stan we need to include it as a variable
  #9.3 Sampling historical data 
    #9.3.1 Stratified sampling
  #9.4 Bayesian model at range scale with historical data 
    #9.4.1 We define the priors location (mean) scale (standard deviation)
    #9.4.2 We write the model  
    #9.4.3 Summaries of priors and estimates
    #9.4.4 Model evaluation
    #9.4.5 We write a bayesian null model 
    #9.4.6 Visualization/Checking the chain #9.8.3.6 CHECKING THE CHAIN
    #9.4.7 Visualization of variable response
    #9.4.8 Spatial Prediction  
      #9.4.8.1 Spatial Prediction  for historic climate
      #9.4.8.2 Spatial Prediction  for current conditions
      #9.4.8.3 Creation of the rasters
        #9.4.8.3.1 Historical prediction
        #9.4.8.3.2 Current prediction
  #9.5 Bayesian model with only current data (no information from historical model for the priors)
    #9.5.1 WE define the priors location (mean) scale (standard deviation)
    #9.5.2 We write the model  
    #9.5.3 Summaries of priors and estimates
    #9.5.4 Model evaluation
    #9.5.5 We write a bayesian null model 
    #9.5.6 Visualization/Checking the chain #9.8.3.6 CHECKING THE CHAIN
    #9.5.7 Visualization of variable response
  #9.6 Bayesian hierarchical model with current data using results of Bayesian model with historical data as priors
    #9.6.1 We define the priors location (mean) scale (standard deviation)
    #9.6.2 We write the model  
    #9.6.3 Summaries of priors and estimates
    #9.6.4 Model evaluation
    #9.6.5 We write a bayesian null model 
    #9.6.6 Visualization/Checking the chain #9.8.3.6 CHECKING THE CHAIN
    #9.6.7 Visualization of variable response  
    #9.6.8 Prediction Bayesian hierarchical
      #9.6.8.1 Preparation of normal future scenarios (current was already prepared)
      #9.6.8.2 Scalation of normal scenarios (current and future)
      #9.6.8.3 Combination of scenarios to simulate scenario with only change in  biotic and scenario with only change in abiotic
      #9.6.8.4 We calculate the prediction for each scenario
      #9.6.8.5 Creation of the rasters with the prediction of habitat and uncertainity
    #9.6.9 Validation of the Bayesian hierarchical
      #9.6.9.1 Validation of the model selecting a threshold including the 90% of presences
      #9.6.9.2 Validation of the model selecting a threshold which maxize the TSS
        #9.6.9.2.1 A unique threshold for all Europe
        #9.6.9.2.2 A different threshold by subpopulation
        #9.6.9.2.3 We binarize the predictions into 0 and 1. We are going to select a threshold that maximize the TSS (A unique threshold for all Europe obtained in #9.6.9.2.1)
    #9.6.10 We asses the predictions for each scenario
        #9.6.10.1 We calcualte the suitable habitat fragments and a variable to calculate the distance in idrisi
      #9.6.11.2 We prepare the data of current range in our study area and of protected areas
        #9.6.11.2.1 current range in our study area 
        #9.6.11.2.2 Protected areas
      #9.6.11.3 We calculate the suitable habitat measures for the habitat predictions
        #9.6.11.3.1 We are going to asses the change in the variables from the current scenario to future  
        #9.6.11.3.2 We are going to asses the change in the variables from the current scenario to future for each subpopulation  
      #9.6.11.4 We calculate the suitable habitat measures for the currently occupied area
        #9.6.11.4.1 We are going to asses the change in the variables from the current scenario to future  
        #9.6.11.4.2 We are going to asses the change in the variables from the current scenario to future for each subpopulation  
      #9.6.11.5 We calculate rasters of habitat change for each scenario
  #9.7 Bayesian hierarchical model ONLY CLIMATE with current data using results of Bayesian model with historical data as priors
  #9.8 Bayesian hierarchical model ONLY LAND USE with current data using results of Bayesian model with historical data as priors
  #9.9 Bayesian hierarchical model ONLY ABIOTIC with current data using results of Bayesian model with historical data as priors
  #9.10 Bayesian hierarchical model ONLY BIOTIC with current data using results of Bayesian model with historical data as priors
  #9.11 Bayesian NOT hierarchical model ONLY ABIOTIC with current data 
  #9.12 Bayesian model ONLY CLIMATE with ONLY current data
  #9.13 Checking the error propagation in the model from 9.6 Bayesian hierarchical model with current data using results of Bayesian model with historical data as priors


################################################################################################    
#9.1  Software installation and data preparation
################################################################################################   

## Run in R 4.02
  rm(list=ls()) 
  folder<-"D:/Project_name/Spatial_Database/R_analysis/"
  setwd(paste0(folder,"Modified/"))
  folder_data<-paste0(folder,"Historical/")
  
  install.packages(c("coda","mvtnorm","devtools","loo"))
library(devtools)
devtools::install_github("rmcelreath/rethinking")
  
  
  
  #9.1.1 Install required packages 
    if(!requireNamespace("spatialEco", quietly=TRUE))
      install.packages("spatialEco", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("rgdal", quietly=TRUE))
      install.packages("rgdal", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("sp", quietly=TRUE))
      install.packages("sp", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("usdm", quietly=TRUE))
      install.packages("usdm", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("MuMIn", quietly=TRUE))
      install.packages("MuMIn", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("splitstackshape", quietly=TRUE))
      install.packages("splitstackshape", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("rstan", quietly=TRUE))
      install.packages("rstan", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("raster", quietly=TRUE))
      install.packages("raster", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("loo", quietly=TRUE))
      install.packages("loo", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("rstanarm", quietly=TRUE))
      install.packages("rstanarm", quiet=TRUE, dependencies=TRUE)
     if(!requireNamespace("xlsx", quietly=TRUE))
      install.packages("xlsx", quiet=TRUE, dependencies=TRUE)
    if(!requireNamespace("bayesplot", quietly=TRUE))
      install.packages("bayesplot", quiet=TRUE, dependencies=TRUE)
    library(spatialEco)
    library(rgdal) 
    library(sp)
    library(usdm)
    library(MuMIn)
    library(splitstackshape)
    library(rethinking)
    library(rstan)
    library(raster)
    library(stats)    
    library(loo)    
    library(rstanarm)  
    library(xlsx)    
    library(rstan)
    library(bayesplot)
    library(usdm)
    library(MuMIn)
    #devtools::install_github("stan-dev/loo")
    #install.packages("splitstackshape")
    library(splitstackshape)
    #install.packages("stats")
    library(stats)   
  

  #9.1.2 We load the data
    #We load the data of presences/absences and variables for historical distribution at 50x50 km scale:
    load(paste0(folder_data,"presences_absences_eurasia_hist.RData"))
    data.frame(names(presences_absences_eurasia_hist))
    presences_absences_eurasia_hist$vec_his<-presences_absences_eurasia_hist$Range
    data.frame(names(presences_absences_eurasia_hist))

################################################################################################    
#9.2  Variable selection of climate variables based in historical range data
################################################################################################   
      #9.2.1 We drop variables highly correlated based in VIF 
        vif_result_clim<-vifstep(presences_absences_eurasia_hist[,8:26], th=10) 
        save(vif_result_clim, file="vif_result_clim.RData")
        names_select_vif_clim<-vif_result_clim @ results $Variables   
        save(names_select_vif_clim, file="names_select_vif_clim.RData")
        data_for_mod_clim <- presences_absences_eurasia_hist[, names(presences_absences_eurasia_hist) %in% names_select_vif_clim, drop = F]
        data_for_mod_clim$vec_his<- presences_absences_eurasia_hist$vec_his
        n_var_average_clim<-length(names_select_vif_clim)

      #9.2.2 We select the best four variables based on univariable GLM including the quadratic effect
        DF_AIC_univar_models_clim<-c()
        for (v in 1:n_var_average_clim){
          var_in_loop_clim<-names_select_vif_clim[v]
          write_model_prefix_clim<-c("(vec_his)~")
          fixed_term_clim<-var_in_loop_clim
          write_model_var_in_loop_clim_1<-paste(var_in_loop_clim,sep="")# Lineal
          write_model_var_in_loop_clim_2<-paste("I(",var_in_loop_clim,"^2",")",sep="")# Quadratic
          write_model_var_in_loop_clim_3<-paste(var_in_loop_clim," + ","I(",var_in_loop_clim,"^2",")",sep="")# Lineal and quadratic
      
          write_model_1_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_1)# Lineal
          write_model_2_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_2)# Quadratic
          write_model_3_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_3)# Lineal and quadratic
          
          from_mod_1_clim<-as.formula(write_model_1_clim)
          from_mod_2_clim<-as.formula(write_model_2_clim)
          from_mod_3_clim<-as.formula(write_model_3_clim)
          
          mod_1_clim<- glm(formula = from_mod_1_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail") 
          mod_2_clim<- glm(formula = from_mod_2_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail") 
          mod_3_clim<- glm(formula = from_mod_3_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail") 
          #mod_null<- glm(formula = (sp_1_europa)~1, data = data_for_mod, family = "binomial", na.action = "na.fail") 
          AICc_mod_1_clim<-AICc(mod_1_clim)
          AICc_mod_2_clim<-AICc(mod_2_clim)
          AICc_mod_3_clim<-AICc(mod_3_clim)
          AICc_mod_clim_123<-cbind(AICc_mod_1_clim,AICc_mod_2_clim,AICc_mod_3_clim)
          AICc_mod_clim_123_min<-min(AICc_mod_clim_123[1,])
          DF_AIC_univar_models_clim<-rbind(DF_AIC_univar_models_clim,AICc_mod_clim_123_min)
        }
        rownames(DF_AIC_univar_models_clim)<-names_select_vif_clim
        sorted_variables_clim<-sort(DF_AIC_univar_models_clim[,1])
        #We select the 4 better
        number_variables_selected_for_models<-4
        variables_less_AICc_clim<-sorted_variables_clim[1:(number_variables_selected_for_models)]
        variables_less_AICc_clim<-data.frame(variables_less_AICc_clim)
        variables_less_AICc_clim$variables_names<-names(variables_less_AICc_clim)
        variables_names_less_AICc_clim<-variables_less_AICc_clim
   
        
        colnames(variables_names_less_AICc_clim)<-c("Model/Variable", "AICc")
        save(variables_names_less_AICc_clim, file="variables_names_less_AICc_clim.RData")
        DF_AIC_univar_models_clim<-as.data.frame(DF_AIC_univar_models_clim)
        DF_AIC_univar_models_clim$Model_Variable<-rownames(DF_AIC_univar_models_clim)
        DF_AIC_univar_models_clim<-DF_AIC_univar_models_clim[,c(2,1)]
        rownames(DF_AIC_univar_models_clim)<-c()
        colnames(DF_AIC_univar_models_clim)<-c("Model/Variable","AICc")
        save(DF_AIC_univar_models_clim, file="DF_AIC_univar_models_clim.RData") # Supplementary Table 61 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            write.csv(DF_AIC_univar_models_clim, file='DF_AIC_univar_models_clim.csv')# # Supplementary Table 62 # Supplementary Table 63 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        print(paste0("OUTPUT: Supplementary Table 61#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 61. Comparasion of univarible models for historical bioclimatic variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        
        
      #9.2.3 We create a subset with the selected variables and the presences/absences
        sub_presences_absences_eurasia_hist<-presences_absences_eurasia_hist[,c("bio03","bio04","bio08","bio09","vec_his")]#,"LC_1","LC_2","LC_1"
        names(sub_presences_absences_eurasia_hist)[names(sub_presences_absences_eurasia_hist) == "vec_his"] <- "presence"
        names(sub_presences_absences_eurasia_hist)[names(sub_presences_absences_eurasia_hist) == "bio03"] <- "Clim_3" 
        names(sub_presences_absences_eurasia_hist)[names(sub_presences_absences_eurasia_hist) == "bio04"] <- "Clim_4"
        names(sub_presences_absences_eurasia_hist)[names(sub_presences_absences_eurasia_hist) == "bio08"] <- "Clim_8"
        names(sub_presences_absences_eurasia_hist)[names(sub_presences_absences_eurasia_hist) == "bio09"] <- "Clim_9"

      #9.2.4 We are going to divide the variables to have a more homogeneous range of data among them
        summary(sub_presences_absences_eurasia_hist)
        sub_presences_absences_eurasia_hist$Clim_3<-sub_presences_absences_eurasia_hist$Clim_3/100
        sub_presences_absences_eurasia_hist$Clim_4<-sub_presences_absences_eurasia_hist$Clim_4/1000
        sub_presences_absences_eurasia_hist$Clim_8<-sub_presences_absences_eurasia_hist$Clim_8/10
        sub_presences_absences_eurasia_hist$Clim_9<-sub_presences_absences_eurasia_hist$Clim_9/10
        summary(sub_presences_absences_eurasia_hist)
        
      #9.2.5 We are going to calculate the cuadratic effect from each variable as in Stan we need to include it as a variable
        sub_presences_absences_eurasia_hist$Clim_3_c<-sub_presences_absences_eurasia_hist$Clim_3^2
        sub_presences_absences_eurasia_hist$Clim_4_c<-sub_presences_absences_eurasia_hist$Clim_4^2
        sub_presences_absences_eurasia_hist$Clim_8_c<-sub_presences_absences_eurasia_hist$Clim_8^2
        sub_presences_absences_eurasia_hist$Clim_9_c<-sub_presences_absences_eurasia_hist$Clim_9^2
        summary(sub_presences_absences_eurasia_hist)
   
   

################################################################################################    
#9.3 Sampling historical data 
################################################################################################   
  #9.3.1 Stratified sampling
    #We select the environmental linear variables
    env<- sub_presences_absences_eurasia_hist[,c(1:4)]
    head(env)
    strata<-c()
    #Here we transform the original data of each variable in three categories. We are going to apply a trick to have a unique value of
    #category which is "10^var" this makes that the first category be in the range cof decimals, the second varible in range of hundredds,
    #the third in the range of thousands... this is an adaptation from the book "Habitat Suitabililty and Distributions Models with applications in R" 
    #from Antoine Guisan, W. Thuiller and N. E. Zimmermann, 7.4 Sampling Design and Data Collection.
    for (var in 1:ncol(env)){
      #We obatain a sequence of four numbers ranging from minimum to the maximunm of the variable
      s<-seq(from=min(env[,var]),to=max(env[,var]),length.out =4)
      df_var1<-as.data.frame(env[,var])
      colnames(df_var1)<-c("Value_original")
      str(df_var1)
      #we combine values to from a vector
      m <- c(s[1], s[2], 1*10^var,  s[2], s[3], 2*10^var,  s[3], s[4], 3*10^var)
      #We obtain a matrix with the range of the variable and with a new integer calue for each class
      rclmat <- matrix(m, ncol=3, byrow=TRUE)
      df_rclmat<-as.data.frame(rclmat)
  
     df3 <- transform(df_var1, Rec_value=df_rclmat$V3[findInterval(Value_original, df_rclmat$V1)])     
     strata<-cbind(strata,df3$Rec_value)     
    }   
    #We have a data frame with three categories for each of the four variables
    strata_df<-as.data.frame(strata)
    save(strata_df, file="strata_df.RData")
    #With this we create a variables with a combination of all the categories
    str<-rowSums(strata_df)
    X11()
    hist(str, main="Histogram values")
    #We see the number of groups and the pixels in each group  
    table(str)
    n_groups<-unique(str)
    str(n_groups)
    sub_presences_absences_eurasia_hist_str<-cbind(sub_presences_absences_eurasia_hist,str)
    head(sub_presences_absences_eurasia_hist_str)
    #We divide betwen presences and absences
          sub2_presences_hist<-subset(sub_presences_absences_eurasia_hist_str,presence==1)
          summary(sub2_presences_hist)
          str(sub2_presences_hist)
  
          sub2_absences_hist<-subset(sub_presences_absences_eurasia_hist_str,presence==0)
          summary(sub2_absences_hist)
          str(sub2_absences_hist)
    #We do a stratified sampling of the presences and absences
    set.seed(1)

    
    out_sub2_presences_hist <- stratified(sub2_presences_hist, c("str"), 50)
    str(out_sub2_presences_hist)
    out_sub2_absences_hist <- stratified(sub2_absences_hist, c("str"), 50)
    str(out_sub2_absences_hist)
    #WE merge presences and basences again
    dat_50km<-rbind(out_sub2_presences_hist,out_sub2_absences_hist)
    dat_50km<-as.data.frame(dat_50km)
    str(dat_50km)
    save(dat_50km,file="dat_50km.Rdata")
   

################################################################################################    
#9.4 Bayesian model at range scale with historical data 
################################################################################################    
  ## Run in R 4.02
    load("dat_50km.Rdata")
  #We say to R that detect the number of cores of the computers
  options(mc.cores = parallel::detectCores())
  
   #install.packages("rstanarm") 
   library(rstanarm)
    
  #9.4.1 WE define the priors location (mean) scale (standard deviation)
  my_prior <- normal(location = c(0,0,0,0,0,0,0,0), scale = c(10,10,10,10,10,10,10,10))
  my_prior_intercept<- normal(location = c(0), scale = c(10))
  
  #9.4.2 We write the model  
  m10.4_rstanarm<-stan_glm(presence~ Clim_3+ Clim_4+Clim_8+Clim_9 +Clim_3_c+Clim_4_c+Clim_8_c+Clim_9_c,
    data = dat_50km, family = "binomial",
    prior_intercept =my_prior_intercept,
    prior = my_prior,
    chains = 4,
    iter = 8000)
  save (m10.4_rstanarm, file='m10.4_rstanarm.RData')
  #load(file="m10.4_rstanarm.RData")
   
 #9.4.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(m10.4_rstanarm)
    
    #A summary of the model  
    summary_m10.4_rstanarm<- as.data.frame(summary(m10.4_rstanarm))
    save (summary_m10.4_rstanarm, file='summary_m10.4_rstanarm.RData')# Supplementary Table 62 # Supplementary Table 63
    write.csv(summary_m10.4_rstanarm, file='summary_m10.4_rstanarm.csv')# # Supplementary Table 62 # Supplementary Table 63 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 62#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 62. Results for the species distribution model at range scale, SDMRange based on bayesian GLM predicting presence (using historical distribution in Eurasia) as a function of bioclimatic variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_m10.4_rstanarm [1:9,]

        print(paste0("OUTPUT: Supplementary Table 63#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 63. Results for the species distribution model at range scale, SDMRange based on bayesian GLM predicting presence (using historical distribution in Eurasia) as a function of bioclimatic variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_m10.4_rstanarm [10,]
        
    #We compare the bayesian model with a frequentist model
    m10.4_glm<-glm(presence~ Clim_3+ Clim_4+Clim_8+Clim_9 +Clim_3_c+Clim_4_c+Clim_8_c+Clim_9_c, data = dat_50km, family = "binomial")
    summary(m10.4_glm)

  #9.4.4 Model evaluation
  #Widely Applicable Information Criterion (WAIC). 
  #References: 
    #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
    #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
    #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
    #and WAIC for evaluating fitted Bayesian models.
    print(paste0("OUTPUT: Supplementary Table 64#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 64. Estimates and standard error (SE) for the expected log pointwise predictive density (elpd_waic), the effective number of parameters (p_waic) and the information criterion waic (which is just -2 * elpd_waic, i.e., converted to deviance scale) for the species distribution model at range scale, SDMRange.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    waic(m10.4_rstanarm) # # Supplementary Table 64 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #From https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
    
    library(loo)
    library(bayesplot)
        library(xlsx)

    loo1 <- loo(m10.4_rstanarm, save_psis = TRUE)#
    X11()  
      plot(loo1)
    pareto_k_table(loo1)
    pareto_k_values(loo1)
    yrep <- posterior_predict(m10.4_rstanarm)
    X11()
    ppc_loo_pit_overlay(
      y = dat_50km$presence,
      yrep = yrep,
      lw = weights(loo1$psis_object)
    )
    X11()
    mcmc_areas(as.matrix(m10.4_rstanarm),prob_outer = .99)

    #9.4.4.1 mcmc_pairs
      X11()
      mcmc_pairs(m10.4_rstanarm)
         
    #9.4.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
      posterior_df <- as.data.frame(m10.4_rstanarm)
      cor_poterior_m10.4_rstanarm<-cor(posterior_df)
      write.xlsx(cor_poterior_m10.4_rstanarm, file='cor_poterior_m10.4_rstanarm.xlsx')
      print(paste0("OUTPUT: Supplementary Table 66#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 66. Correlation of the posterior samples among the predictors used in the SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(cor_poterior_m10.4_rstanarm)
      

  #9.4.5 We write a bayesian null model 
    #WE define the priors location (mean) scale (standard deviation)
    #my_prior <- normal(location = c(0,0,0,0,0,0,0,0), scale = c(10,10,10,10,10,10,10,10))#, autoscale = FALSE
    my_prior_intercept<- normal(location = c(0), scale = c(10))#, autoscale = FALSE

    m10.4_rstanarm_null<-stan_glm(presence~ 1,
      data = dat_50km, family = "binomial",
      prior_intercept =my_prior_intercept,
      #prior = my_prior,
      chains = 4,
      iter = 8000)
    save (m10.4_rstanarm_null, file='m10.4_rstanarm_null.RData')
      #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
        #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
        #and WAIC for evaluating fitted Bayesian models.
    
    print(paste0("OUTPUT: Supplementary Table 65#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 65. Estimates and standard error (SE) for the expected log pointwise predictive density (elpd_waic), the effective number of parameters (p_waic) and the information criterion waic (which is just -2 * elpd_waic, i.e., converted to deviance scale) for a null species distribution model at range scale, SDMRange_Null~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    waic(m10.4_rstanarm_null) # # Supplementary Table 65 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    save(m10.4_rstanarm_null,file="m10.4_rstanarm_null.Rdata")

  #9.4.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(m10.4_rstanarm)
    dim(posterior)
    str(posterior)
    
    print(paste0("OUTPUT: Supplementary Figure 14#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 14.Chains for the Species Distribution Model at range scale, the SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 14 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 1, strip.position = "left"))
    X11()
    plot(m10.4_rstanarm)  

  #9.4.7 Visualization of variable response
    #We are gouing to calculate to see how is the response to changes for each variable.
    #We are going to change values of a variable meanwhile we mantein constant with a mean value the other variables
    #Prepare new counterfactual data R code 5.9
    load("m10.4_rstanarm.Rdata")
    
    sub_variables<-dat_50km[,c(1:4,6:9)]
    head(dat_50km)
    head(sub_variables)
    nvariables<-length(sub_variables)
    variables_names<-colnames(sub_variables)
    
    Clim_3_median<-median(dat_50km$Clim_3)    
    Clim_4_median<-median(dat_50km$Clim_4) 
    Clim_8_median<-median(dat_50km$Clim_8)    
    Clim_9_median<-median(dat_50km$Clim_9)    

    print(paste0("OUTPUT: Supplementary Figure 7a~~~~Documents\Bear\Paper\Figures\Dibujos~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 7a Response plot of the three Bayesian models explaining the distribution of the brown bear. a Bayesian model at range scale using the historical distribution of brown bear and historical climate variables as predictors.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    #X11() # PDF 10 x 8 inches
    par(mfrow=c(3,4))
    par(mar=c(2,2,2,0.5))
    
    all_lines<-c()
    all_shades<-c()
    all_var_in_loop_seq<-c()
    
      for (var in 1:4){
        #We create a data frame with median values. For see the response of a variable we will fixed the values of all other 
        #variables using these median values
        pred.data<- data.frame(
          Clim_3=Clim_3_median,
          Clim_4=Clim_4_median,
          Clim_8=Clim_8_median,
          Clim_9=Clim_9_median,
          Clim_3_c=Clim_3_median^2,  
          Clim_4_c=Clim_4_median^2,
          Clim_8_c=Clim_8_median^2,  
          Clim_9_c=Clim_9_median^2
        )
        #We create a dataframe with 1000 rows repeating the data
        pred.data<-pred.data[rep(seq_len(nrow(pred.data)), each = 1000), ]
      
        #We see the varriables in the loop. The linear and the quadratic names  
        var_in_loop<-sub_variables[,c(var)]  
        var_in_loop_c<-sub_variables[,c(var+4)] 
        var_name_in_loop<-variables_names[c(var)]  
        var_name_in_loop_c<-variables_names[c(var+4)] 
     
        #We create a sequence of values for the variables ranging from the minimum to tje maximum of all the values oberved and with a length of 1000
        var_in_loop_seq<-seq(from=min(dat_50km[,var_name_in_loop]), to=max(dat_50km[,var_name_in_loop]), length.out = 1000)
        var_c_in_loop_seq<-var_in_loop_seq^2
        #We sustitute the constant mean values dataframe with 1000 rows repeating the data
        # with the sequences ranging from min-max of the variables in the loop, created in the two above lines
        pred.data[,c(var)]<-var_in_loop_seq
        pred.data[,c(var+4)]<-var_c_in_loop_seq
        #We extract the posterior
        stan.mu <- posterior_epred(m10.4_rstanarm, newdata = pred.data)
        ## R code 4.56
        # summarize the distribution of mu
        mu.mean <- apply( stan.mu , 2 , mean )
        mu.HPDI <- apply( stan.mu , 2 , HPDI , prob=0.95 )
        ## R code 4.57
        # plot raw data
        # fading out points to make line and interval more visible
        plot( presence ~ dat_50km[,c(var)],  data=dat_50km , col=col.alpha(rangi2,0.5), cex=0,pch = 20, main=(var_name_in_loop))
        # plot the MAP line, aka the mean mu for each weight
        lines( var_in_loop_seq , mu.mean )
        # plot a shaded region for 89% HPDI
        shade( mu.HPDI , var_in_loop_seq )
        #Add points of presences/absences:
       #  points(x = pirates$height[pirates$sex == "male"],
       # y = pirates$weight[pirates$sex == "male"],
       # pch = 16,
       # col = transparent("coral2", trans.val = .8))
        # points(x=dat_50km[c(var_name_in_loop)][[1]],
        #        y=dat_50km[c("presence")][[1]],
        #         pch = 16,
        #        
        #         col = alpha("red",0.01))
      }
    

    
  #9.4.8 Spatial Prediction  
    #9.4.8.1 Spatial Prediction  for historic climate
      #We load the data for bioclimatic variables in the world for 1981
          load(file="F:/G/Project_name/Historical_range/CRU_TS_v_402/Climate_1901/bioclimatic_var_1901_cea.RData") #borarlinea
      load(file="bioclimatic_var_1901_cea.RData") #borarlinea
      load(paste0(my_dir,"/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1901/bioclimatic_var_1901_cea.RData")) #borarlinea
      bioclimatic_var_1901_cea<-as.data.frame(bioclimatic_var_1901_cea)
      names(bioclimatic_var_1901_cea)[names(bioclimatic_var_1901_cea) == "bio3"] <- "Clim_3" 
      names(bioclimatic_var_1901_cea)[names(bioclimatic_var_1901_cea) == "bio4"] <- "Clim_4"
      names(bioclimatic_var_1901_cea)[names(bioclimatic_var_1901_cea) == "bio8"] <- "Clim_8"
      names(bioclimatic_var_1901_cea)[names(bioclimatic_var_1901_cea) == "bio9"] <- "Clim_9"
      
      sub_pred_bioclimatic_var_1901_cea<-bioclimatic_var_1901_cea[,c("Clim_3","Clim_4","Clim_8","Clim_9")]
      summary(sub_pred_bioclimatic_var_1901_cea)
      sub_pred_bioclimatic_var_1901_cea$Clim_3<-sub_pred_bioclimatic_var_1901_cea$Clim_3/100
      sub_pred_bioclimatic_var_1901_cea$Clim_4<-sub_pred_bioclimatic_var_1901_cea$Clim_4/1000
      sub_pred_bioclimatic_var_1901_cea$Clim_8<-sub_pred_bioclimatic_var_1901_cea$Clim_8/10
      sub_pred_bioclimatic_var_1901_cea$Clim_9<-sub_pred_bioclimatic_var_1901_cea$Clim_9/10
      
      sub_pred_bioclimatic_var_1901_cea$Clim_3_c<-sub_pred_bioclimatic_var_1901_cea$Clim_3^2
      sub_pred_bioclimatic_var_1901_cea$Clim_4_c<-sub_pred_bioclimatic_var_1901_cea$Clim_4^2
      sub_pred_bioclimatic_var_1901_cea$Clim_8_c<-sub_pred_bioclimatic_var_1901_cea$Clim_8^2
      sub_pred_bioclimatic_var_1901_cea$Clim_9_c<-sub_pred_bioclimatic_var_1901_cea$Clim_9^2
      #We add and ID to merge later the predictions
      sub_pred_bioclimatic_var_1901_cea$ID<-seq(1:201600)
      #We omit the NA values to be able to obtain posterior
      sub_pred_bioclimatic_var_1901_cea_no_NA<-na.omit(sub_pred_bioclimatic_var_1901_cea)
      
      #We extract the posterior
      link.m10.4_1901 <- posterior_epred(m10.4_rstanarm, newdata = sub_pred_bioclimatic_var_1901_cea_no_NA)
      str(link.m10.4_1901)  
      save(link.m10.4_1901,file="link.m10.4_1901.Rdata")
      # summarize, we calculate the mean of all the columns of the matrix link.m10.4
      sum_mean_1901<-apply( link.m10.4_1901 , 2 , mean)
      #We compute the percentile intervals, we select prob=0.95
      #sum_PI<-apply( link.m10.4 , 2 , PI)
      #This is a modification of above original code for data with Na values
      sum_PI_1901 <- apply( link.m10.4_1901 , 2 , HPDI , prob=0.95 )

      df_sum_mean_1901<-as.data.frame(sum_mean_1901)
      df_sum_PI_1901<-as.data.frame(sum_PI_1901)
      df_sum_PI_t_1901<-t(df_sum_PI_1901)  
      df_sum_PI_t_1901<-as.data.frame(df_sum_PI_t_1901)
  
      ##Prediction (I wrote Europe but is predicted for all the world)
      df_prediciton_europe_1901<-cbind(df_sum_mean_1901,df_sum_PI_t_1901)
      df_prediciton_europe_1901$Uncertainity<-df_prediciton_europe_1901[,c(3)]-df_prediciton_europe_1901[,c(2)]
      df_prediciton_europe_1901<-cbind(df_prediciton_europe_1901,sub_pred_bioclimatic_var_1901_cea_no_NA$ID)
      colnames(df_prediciton_europe_1901)<-c("sum_mean","|0.95","0.95|","Uncertainity","ID")
      merge_Edf_prediciton_europe_1901<-merge(sub_pred_bioclimatic_var_1901_cea,df_prediciton_europe_1901, by.x="ID", by.y = "ID",all.x = T)
      save(merge_Edf_prediciton_europe_1901,file="merge_Edf_prediciton_europe_1901.Rdata")
      
      #9.4.8.1.1 We are going to calculate the cutoff for retain the 95% of the historical range for Eurasia (Which is the area that we used to fit the model)
      load("F:/G/Project_name/Historical_range/CRU_TS_v_402/Climate_1901/DF_bio1901_range.RData")#borarlinea
      load(paste0(my_dir,"/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1901/DF_bio1901_range.RData")) #borarlinea
        v_sum_mean<-merge_Edf_prediciton_europe_1901$sum_mean
        DF_bio1901_range2<-cbind(DF_bio1901_range,v_sum_mean)
           #Analysis with Eurasian brown bear
                  #We subset areas based in latitude and longitude
                  DF_bio1901_range2_eurasia<-subset(DF_bio1901_range2, x>-7000000 & x<20000000 & y>0 & y< 6300000)
                  head(DF_bio1901_range2_eurasia)
                  DF_1901_observed_predicted<-DF_bio1901_range2_eurasia[,c("vec_his","v_sum_mean")]
                  colnames(DF_1901_observed_predicted)<-c("Observed","Predicted")
          
        sub_DF_predicted_observed_of_training<-DF_1901_observed_predicted[DF_1901_observed_predicted$Observed==1,]
        summary(sub_DF_predicted_observed_of_training)
        order_data<-sub_DF_predicted_observed_of_training[order(sub_DF_predicted_observed_of_training$Predicted),]
        str(order_data)
        quantile(order_data$Predicted, probs = c(0, 0.05, 0.9, 1)) # quartile 5% = 0.29388912 # 0.38685954


    #9.4.8.2 Spatial Prediction  for current conditions
      #We load the data for bioclimatic variables in the world for 1981
      load(file="F:/G/Project_name/Historical_range/CRU_TS_v_402/Climate_1981/bioclimatic_var_1981_cea.RData")
      load(paste0(my_dir,"/7_Construction_of_Historical_database/CRU_TS_v_402/Climate_1981/bioclimatic_var_1981_cea.RData")) 
      bioclimatic_var_1981_cea<-as.data.frame(bioclimatic_var_1981_cea)
      names(bioclimatic_var_1981_cea)[names(bioclimatic_var_1981_cea) == "bio3"] <- "Clim_3" 
      names(bioclimatic_var_1981_cea)[names(bioclimatic_var_1981_cea) == "bio4"] <- "Clim_4"
      names(bioclimatic_var_1981_cea)[names(bioclimatic_var_1981_cea) == "bio8"] <- "Clim_8"
      names(bioclimatic_var_1981_cea)[names(bioclimatic_var_1981_cea) == "bio9"] <- "Clim_9"
      
      sub_pred_bioclimatic_var_1981_cea<-bioclimatic_var_1981_cea[,c("Clim_3","Clim_4","Clim_8","Clim_9")]
      summary(sub_pred_bioclimatic_var_1981_cea)
      sub_pred_bioclimatic_var_1981_cea$Clim_3<-sub_pred_bioclimatic_var_1981_cea$Clim_3/100
      sub_pred_bioclimatic_var_1981_cea$Clim_4<-sub_pred_bioclimatic_var_1981_cea$Clim_4/1000
      sub_pred_bioclimatic_var_1981_cea$Clim_8<-sub_pred_bioclimatic_var_1981_cea$Clim_8/10
      sub_pred_bioclimatic_var_1981_cea$Clim_9<-sub_pred_bioclimatic_var_1981_cea$Clim_9/10
      
      sub_pred_bioclimatic_var_1981_cea$Clim_3_c<-sub_pred_bioclimatic_var_1981_cea$Clim_3^2
      sub_pred_bioclimatic_var_1981_cea$Clim_4_c<-sub_pred_bioclimatic_var_1981_cea$Clim_4^2
      sub_pred_bioclimatic_var_1981_cea$Clim_8_c<-sub_pred_bioclimatic_var_1981_cea$Clim_8^2
      sub_pred_bioclimatic_var_1981_cea$Clim_9_c<-sub_pred_bioclimatic_var_1981_cea$Clim_9^2
      #We add and ID to merge later the predictions
      sub_pred_bioclimatic_var_1981_cea$ID<-seq(1:201600)
      #We omit the NA values to be able to obtain posterior
      sub_pred_bioclimatic_var_1981_cea_no_NA<-na.omit(sub_pred_bioclimatic_var_1981_cea)
      #We extract the posterior
      link.m10.4_1981 <- posterior_epred(m10.4_rstanarm, newdata = sub_pred_bioclimatic_var_1981_cea_no_NA)
      str(link.m10.4_1981)  
      save(link.m10.4_1981,file="link.m10.4_1981.Rdata")
      #summarize, we calculate the mean of all the columns of the matrix link.m10.4
      sum_mean_1981<-apply( link.m10.4_1981 , 2 , mean)
      #We compute the percentile intervals, we select prob=0.95
      sum_PI_1981 <- apply( link.m10.4_1981 , 2 , HPDI , prob=0.95 )

      df_sum_mean_1981<-as.data.frame(sum_mean_1981)
      df_sum_PI_1981<-as.data.frame(sum_PI_1981)
      df_sum_PI_t_1981<-t(df_sum_PI_1981)  
      df_sum_PI_t_1981<-as.data.frame(df_sum_PI_t_1981)

      ##Prediction (I wrote Europe but is predicted for all the world)
      df_prediciton_europe_1981<-cbind(df_sum_mean_1981,df_sum_PI_t_1981)
      df_prediciton_europe_1981$Uncertainity<-df_prediciton_europe_1981[,c(3)]-df_prediciton_europe_1981[,c(2)]
      df_prediciton_europe_1981<-cbind(df_prediciton_europe_1981,sub_pred_bioclimatic_var_1981_cea_no_NA$ID)
      colnames(df_prediciton_europe_1981)<-c("sum_mean","|0.95","0.95|","Uncertainity","ID")
      merge_Edf_prediciton_europe_1981<-merge(sub_pred_bioclimatic_var_1981_cea,df_prediciton_europe_1981, by.x="ID", by.y = "ID",all.x = T)
      save(merge_Edf_prediciton_europe_1981,file="merge_Edf_prediciton_europe_1981.Rdata")

      
    #9.4.8.3 Creation of the rasters
      #We create rasters in img format with the mean values of the predicion and the range of the confidence intervals  
      #We change the R version, close and open R in earlier version (3.5.3) because a strange error:
      print("We change the R version, close and open R in earlier version (3.5.3)")
      ## Run in R 3.5.3
      rm(list=ls()) 
      my_dir <-"writehereyourpath" # path of the folder where your data are stored
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)

      #9.4.8.3.1 Historical prediction
        #Install required packages 
        if(!requireNamespace("raster", quietly=TRUE))
          install.packages("raster", quiet=TRUE, dependencies=TRUE)
        if(!requireNamespace("rgdal", quietly=TRUE))
          install.packages("rgdal", quiet=TRUE, dependencies=TRUE)
        if(!requireNamespace("sp", quietly=TRUE))
          install.packages("sp", quiet=TRUE, dependencies=TRUE)
        if(!requireNamespace("GDAL", quietly=TRUE))
          install.packages("GDAL", quiet=TRUE, dependencies=TRUE)
        if(!requireNamespace("PROJ", quietly=TRUE))
          install.packages("PROJ", quiet=TRUE, dependencies=TRUE)
        library(raster)
        library(rgdal) 
        library (sp)
        library (PROJ)
        load(file="merge_Edf_prediciton_europe_1901.Rdata")
        #Raster of reference: 
        cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" 
        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data
        raster_cea <- raster(xmn=-20000000, xmx=20000000, ymn=-6300000, ymx=6300000,  ncols=800, nrows=252)
        projection(raster_cea) <- cea_proj
        raster_habitat_bear_cea_1901<-raster_cea
        values(raster_habitat_bear_cea_1901)<-round(merge_Edf_prediciton_europe_1901$sum_mean,digits=2)
        save(raster_habitat_bear_cea_1901,file="raster_habitat_bear_cea_1901.Rdata")
        writeRaster(raster_habitat_bear_cea_1901, "raster_habitat_bear_cea_1901.img", overwrite=TRUE) #Supplementary Figure 15a
        print(paste0("OUTPUT: Supplementary Figure 15a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 15a. Map at World scale with mean predicted probabilities brown bear distribution using model at range scale, SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()
        plot(raster_habitat_bear_cea_1901) # Supplementary Figure 15a #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        raster_uncert_bear_cea_1901<-raster_cea
        values(raster_uncert_bear_cea_1901)<-round(merge_Edf_prediciton_europe_1901$Uncertainity,digits=3)
        save(raster_uncert_bear_cea_1901,file="raster_uncert_bear_cea_1901.Rdata")
        writeRaster(raster_uncert_bear_cea_1901, "raster_uncert_bear_cea_1901.img", overwrite=TRUE) 
        print(paste0("OUTPUT: Supplementary Figure 15b~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 15b. Map at World scale with mean predicted probabilities brown bear distribution using model at range scale, SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()
        plot(raster_uncert_bear_cea_1901) # Supplementary Figure 15a #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      #9.4.8.3.2 Current prediction
        load(file="merge_Edf_prediciton_europe_1981.Rdata")
        #Raster of reference: 
        cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" 
        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data
        raster_cea <- raster(xmn=-20000000, xmx=20000000, ymn=-6300000, ymx=6300000,  ncols=800, nrows=252)
        projection(raster_cea) <- cea_proj
        raster_habitat_bear_cea_1981<-raster_cea
        values(raster_habitat_bear_cea_1981)<-round(merge_Edf_prediciton_europe_1981$sum_mean,digits=2)
        save(raster_habitat_bear_cea_1981,file="raster_habitat_bear_cea_1981.Rdata")
        writeRaster(raster_habitat_bear_cea_1981, "raster_habitat_bear_cea_1981.img", overwrite=TRUE) #Supplementary Figure 15c
        print(paste0("OUTPUT: Supplementary Figure 15c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 15c. Map at World scale with mean predicted probabilities brown bear distribution using model at range scale, SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()
        plot(raster_habitat_bear_cea_1981) # Supplementary Figure 15a #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        raster_uncert_bear_cea_1981<-raster_cea
        values(raster_uncert_bear_cea_1981)<-round(merge_Edf_prediciton_europe_1981$Uncertainity,digits=3)
        save(raster_uncert_bear_cea_1981,file="raster_uncert_bear_cea_1981.Rdata")
        writeRaster(raster_uncert_bear_cea_1981, "raster_uncert_bear_cea_1981.img", overwrite=TRUE) #Supplementary Figure 15d
        print(paste0("OUTPUT: Supplementary Figure 15d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 15d. Map at World scale with mean predicted probabilities brown bear distribution using model at range scale, SDMRange~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()
        plot(raster_uncert_bear_cea_1981) # Supplementary Figure 15d #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
######################################################################################################################  
#9.5 Bayesian model with only current data (no information from historical model for the priors)
######################################################################################################################  

    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)
 #We load the presences
    load("RANDOM_presences_absences_for_model_all_5km.Rdata") 
  
## Run in R 4.02
rm(list=ls()) 
    my_dir <-"writehereyourpath"
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)
    #We load the presences
    load(paste0(folder_working,"/8_Univarable_SDM_brown_bear_comparing_biotic_variables/RANDOM_presences_absences_for_model_all_5km.Rdata"))
    presences_absences_5km<-RANDOM_presences_absences_for_model_all_5km
    

    #We check the best land use variables
    load(paste0(folder_working,"/8_Univarable_SDM_brown_bear_comparing_biotic_variables/DF_AIC_univar_models_habitat_land_use.RData"))
    DF_AIC_univar_models_habitat_land_use
    #We check the best biotic variables
    load(paste0(folder_working,"/8_Univarable_SDM_brown_bear_comparing_biotic_variables/DF_AIC_univar_models_habitat_energy2.RData"))
    DF_AIC_univar_models_habitat_energy2
    
    #We are going to scale the variables that we use with uninformative priors (the 4 land use variables and the 4 biotic variables):
      var_model_uninformative_priors_only_current<-presences_absences_5km[,c(23,25,26,31,40,41,43,44)]
      str(var_model_uninformative_priors_only_current)
      save(var_model_uninformative_priors_only_current, file="var_model_uninformative_priors_only_current.RData")
      #See https://stats.stackexchange.com/questions/89172/how-to-scale-new-observations-for-making-predictions-when-the-model-was-fitted-w
      scaled.var_model_uninformative_priors_only_current <- scale(var_model_uninformative_priors_only_current, scale=TRUE)
      save(scaled.var_model_uninformative_priors_only_current, file="scaled.var_model_uninformative_priors_only_current.RData")
      summary(scaled.var_model_uninformative_priors_only_current)
      data.frame(names(presences_absences_5km))
      
    #This dataframe include the scaled and the non scaled values
    #presences_absences_5km_scaled_only_current<-as.data.frame(cbind(presences_absences_5km[,c(1:3,6,7,11,13,118:121)],scaled.var_model_uninformative_priors_only_current,presences_absences_5km[,c("factor_subpop","Europe_presences_0_1")]))
    presences_absences_5km_scaled_only_current<-as.data.frame(cbind(presences_absences_5km[,c(1:3,6,7,11,12,69:72)],scaled.var_model_uninformative_priors_only_current,presences_absences_5km[,c("factor_subpop","Europe_presences_0_1")]))
    str(presences_absences_5km_scaled_only_current)
    save(presences_absences_5km_scaled_only_current, file="presences_absences_5km_scaled_only_current.RData")###########################################################################################################################################################################################################################
    #We say to R that detect the number of cores of the computers
    options(mc.cores = parallel::detectCores())
  
  #9.5.1 WE define the priors location (mean) scale (standard deviation)

  #9.5.2 We write the model  
  #Time 7334.11 seconds (Total)    
    mod_stan_glmer.2<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag + ene_Wild_rep_plant + ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = normal(), prior_intercept = normal())
    save(mod_stan_glmer.2,file="mod_stan_glmer.2.Rdata")
    #load(file="mod_stan_glmer.2.RData")

  #9.5.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_glmer.2)
    #A summary of the model  
    summary_mod_stan_glmer.2<- as.data.frame(summary(mod_stan_glmer.2))
    save (summary_mod_stan_glmer.2, file='summary_mod_stan_glmer.2.RData')
    write.csv(summary_mod_stan_glmer.2, file='summary_mod_stan_glmer.2.csv') # Supplementary Table 45 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 45#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 45.Results for the simple Bayesian model (no hierarchical) using abiotic and biotic factors to explain brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_glmer.2 [1:31,]

        print(paste0("OUTPUT: Supplementary Table 43 mean_PPD Model SHMABC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_glmer.2 [33,]
        
        
  #9.5.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      
      print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHMABC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_glmer.2)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_glmer.2, save_psis = TRUE)
      X11()
      plot(loo1)
      pareto_k_table(loo1)
      pareto_k_values(loo1)

  #9.5.4.1 mcmc_pairs
    #X11()
    #mcmc_pairs(mod_stan_glmer.2)
       
  #9.5.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_glmer.2)
    str(posterior_df)
    cor_poterior_mod_stan_glmer.2<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_glmer.2, file='cor_poterior_mod_stan_glmer.2.xlsx')
    print(paste0("OUTPUT: Supplementary Table 50#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 50. Correlation of the posterior samples among the predictors used in the simple Bayesian model (no hierarchical) using abiotic and biotic factors to explain brown bear distribution~~~"))
    print(cor_poterior_mod_stan_glmer.2)
    
  #9.5.5 We write a bayesian null model 
    #WE define the priors location (mean) scale (standard deviation)
    mod_stan_glmer.2_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = normal(), prior_intercept = normal())
    save (mod_stan_glmer.2_null, file='mod_stan_glmer.2_null.RData')
    load(file="mod_stan_glmer.2_null.RData")
    
    summary_mod_stan_glmer.2_null<- as.data.frame(summary(mod_stan_glmer.2_null))
    save (summary_mod_stan_glmer.2_null, file='summary_mod_stan_glmer.2_null.RData')
    write.csv(summary_mod_stan_glmer.2_null, file='summary_mod_stan_glmer.2_null.csv')

    print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHM_Null #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    waic(mod_stan_glmer.2_null)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #9.5.6 Visualization/Checking the chain
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_glmer.2)
    dim(posterior)
    str(posterior)
    print(paste0("OUTPUT: Supplementary Figure 9 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 9.Chains for the Bayesian model of brown bear habitat with abiotic and biotic factors using current data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 9 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("Clim_3" "Clim_4" "Clim_8"), 
    facet_args = list(ncol = 6, strip.position = "left"))
    X11()
    plot(mod_stan_glmer.2)  

  #9.5.7 Visualization of variable response
    #We calculate the median values for the variables of the model    
    #We are gouing to calculate to see how is the response to changes for each variable.
    #We are going to change values of a variable meanwhile we mantein constant with a mean value the other variables
    #Prepare new counterfactual data R code 5.9
    load("mod_stan_glmer.2.Rdata")
    load("dat_50km.Rdata") 
    load("presences_absences_5km_scaled_only_current.RData") 

    sub_variables<-dat_50km[,c(1:4,6:9)]
    head(sub_variables)
    nvariables<-length(sub_variables)
    variables_names<-colnames(sub_variables)
    
    Clim_3_median_h<-median(dat_50km$Clim_3)    
    Clim_4_median_h<-median(dat_50km$Clim_4) 
    Clim_8_median_h<-median(dat_50km$Clim_8)    
    Clim_9_median_h<-median(dat_50km$Clim_9)    

    Clim_3_median<-median(presences_absences_5km_scaled_only_current$Clim_3)    
    Clim_4_median<-median(presences_absences_5km_scaled_only_current$Clim_4) 
    Clim_8_median<-median(presences_absences_5km_scaled_only_current$Clim_8)    
    Clim_9_median<-median(presences_absences_5km_scaled_only_current$Clim_9)    

    LC_1_median<-median(presences_absences_5km_scaled_only_current$LC_1)
    LC_3_median<-median(presences_absences_5km_scaled_only_current$LC_3)
    LC_4_median<-median(presences_absences_5km_scaled_only_current$LC_4)
    Frag_median<-median(presences_absences_5km_scaled_only_current$Frag)

    ene_Wild_rep_plant_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_rep_plant)
    ene_Wild_invertebrates_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_invertebrates)
    ene_Wild_unk_plant_oth_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_unk_plant_oth)
    ene_Wild_vertebrates_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_vertebrates)

    print(paste0("OUTPUT: Supplementary Figure 7b~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 7b. Response plot of the three Bayesian models explaining the distribution of the brown bear.b) Simple Bayesian model (no hierarchical) using abiotic and biotic factors to explain brown bear distribution.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    #X11() # PDF 8 x 7 inches
    par(mfrow=c(3,4))
    par(mar=c(2,2,2,0.5))
    
    all_lines<-c()
    all_shades<-c()
    for (var in 1:12){
    #for (var in 1:1){
            #This is a dataframe with constant values for all variables/predictors using median values from the historical data (Data with all the range of climate variables across the distribution)
            pred.data_h<- data.frame(
            Clim_3=Clim_3_median_h,
            Clim_4=Clim_4_median_h,
            Clim_8=Clim_8_median_h,
            Clim_9=Clim_9_median_h,
            Clim_3_c=Clim_3_median_h^2,  
            Clim_4_c=Clim_4_median_h^2,
            Clim_8_c=Clim_8_median_h^2,  
            Clim_9_c=Clim_9_median_h^2,
            LC_1=LC_1_median,
            LC_3=LC_3_median,
            LC_4=LC_4_median,
            Frag=Frag_median,
            ene_Wild_rep_plant=ene_Wild_rep_plant_median,
            ene_Wild_unk_plant_oth=ene_Wild_unk_plant_oth_median,
            ene_Wild_invertebrates =ene_Wild_invertebrates_median,
            ene_Wild_vertebrates =ene_Wild_vertebrates_median,
            factor_subpop=as.integer(1)
          )
  
      #This is a dataframe with constant values for all variables/predictors  using median values from the current data (Data with a small part of the range of climate variables across the distribution)
      pred.data.current<- data.frame(
            Clim_3=Clim_3_median,
            Clim_4=Clim_4_median,
            Clim_8=Clim_8_median,
            Clim_9=Clim_9_median,
            Clim_3_c=Clim_3_median^2,  
            Clim_4_c=Clim_4_median^2,
            Clim_8_c=Clim_8_median^2,  
            Clim_9_c=Clim_9_median^2,
            LC_1=LC_1_median,
            LC_3=LC_3_median,
            LC_4=LC_4_median,
            Frag=Frag_median,
            ene_Wild_rep_plant=ene_Wild_rep_plant_median,
            ene_Wild_unk_plant_oth=ene_Wild_unk_plant_oth_median,
            ene_Wild_invertebrates =ene_Wild_invertebrates_median,
            ene_Wild_vertebrates =ene_Wild_vertebrates_median,
            factor_subpop=as.integer(1)
          )
      variables_names<-colnames(pred.data.current)
      variables_names_h<-colnames(pred.data_h)
        #We create a dataframe with 1000 rows repeating the data
        pred.data.current<-pred.data.current[rep(seq_len(nrow(pred.data.current)), each = 1000), ]
        #We create a dataframe with 1000 rows repeating the data
        pred.data_h<-pred.data_h[rep(seq_len(nrow(pred.data_h)), each = 1000), ]
        
     #Now for each variable (var) we are going to substitute the values of pred.data.current with the median
      # by a sequence of values  
        if(var<5){ 
          var<-var
          print(paste("Variable in loop ", var))
          var_in_loop<-pred.data_h[,c(var)]  
          var_in_loop_c<-pred.data_h[,c(var+4)] 
          var_name_in_loop<-variables_names_h[c(var)]  
          var_name_in_loop_c<-variables_names_h[c(var+4)] 
          print(paste("Linear Variable in loop    ", var_name_in_loop))
          print(paste("Quadratic Variable in loop ", var_name_in_loop_c))
          #If we plot using the range of values of the historical range:
            #var_in_loop_seq<-seq(from=min(dat_50km[,var_name_in_loop]), to=max(dat_50km[,var_name_in_loop]), length.out = 1000)
          #If we plot using the range of values of current range:
            var_in_loop_seq<-seq(from=min(dat_50km[,var_name_in_loop]), to=max(dat_50km[,var_name_in_loop]), length.out = 1000)
          #We calculate the sequence of values for the quadratic effect:  
            var_c_in_loop_seq<-var_in_loop_seq^2
          #We sustitute the sequence of values for the variables in the data frame with the values to predict
          pred.data_h[,c(var)]<-var_in_loop_seq
          pred.data_h[,c(var+4)]<-var_c_in_loop_seq
          plot( presence ~ dat_50km[,c(var)],  data=dat_50km , col=col.alpha(rangi2,0.5), cex=0,pch = 20, main=(var_name_in_loop))
        }  
        if(var>4){ 
          var<-var+4
          print(paste("Variable in loop ", var))
          var_in_loop<-pred.data.current[,c(var)] 
          
          var_name_in_loop<-variables_names[c(var)]  
          print(paste("Variable >4 in loop    ", var_name_in_loop))
  
          var_in_loop_seq<-seq(from=min(presences_absences_5km_scaled_only_current[,var_name_in_loop]), to=max(presences_absences_5km_scaled_only_current[,var_name_in_loop]), length.out = 1000)
          pred.data.current[,c(var)]<-var_in_loop_seq
          plot( Europe_presences_0_1 ~ presences_absences_5km_scaled_only_current[,var_name_in_loop],  data=presences_absences_5km_scaled_only_current , col="blue", cex=0,pch = 20, main=(var_name_in_loop))
          pred.data_h<-pred.data.current
          }  
        #We extract the posterior
        stan.mu <- posterior_epred(mod_stan_glmer.2, newdata = pred.data_h)
        #we obtain a matrix with the number of columns as the number of values with which we fed it with 
        #this number of values "length.out = 1000"
    # summarize the distribution of mu
    mu.mean <- apply( stan.mu , 2 , mean )
    mu.HPDI <- apply( stan.mu , 2 , HPDI , prob=0.95 )
    # plot the MAP line, aka the mean mu for each weight
    rect(min(presences_absences_5km_scaled_only_current[,var_name_in_loop]),0,max(presences_absences_5km_scaled_only_current[,var_name_in_loop]),1,col = "lightblue1")
    # plot a shaded region for 89% HPDI
    shade( mu.HPDI , var_in_loop_seq ,col = "grey80")
    lines( var_in_loop_seq , mu.mean )
    
    # points(x=presences_absences_5km_scaled_only_current[c(var_name_in_loop)][[1]],
    #            y=presences_absences_5km_scaled_only_current[c("Europe_presences_0_1")][[1]],
    #             pch = 19,
    # 
    #             col = alpha("red",0.002))

    }#Save as 8x 7 inches and reduce 83.4%

######################################################################################################################  
#9.6 Bayesian hierarchical model with current data using results of Bayesian model with historical data as priors
###################################################################################################################### 
 
  ## Run in R 4.02
  rm(list=ls()) 
  my_dir <-"writehereyourpath" 
  folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
  setwd(folder_working)

   #We load the presences
   load("presences_absences_5km_scaled_only_current.RData") 
    data.frame(names(presences_absences_5km_scaled_only_current))
    summary(presences_absences_5km_scaled_only_current)
    
    
    
   #We load the previous fitted model at range scale ussing historical data (we need to calculate the informative priors)
   load("summary_m10.4_rstanarm.RData")
   #We say to R that detect the number of cores of the computers
   options(mc.cores = parallel::detectCores())
  
  #9.6.1 We define the priors location (mean) scale (standard deviation)
    my_prior_integrated <- normal(location = c(summary_m10.4_rstanarm[2,1],summary_m10.4_rstanarm[3,1],summary_m10.4_rstanarm[4,1],summary_m10.4_rstanarm[5,1],summary_m10.4_rstanarm[6,1],summary_m10.4_rstanarm[7,1],summary_m10.4_rstanarm[8,1],summary_m10.4_rstanarm[9,1],0,0,0,0,0,0,0,0), scale = c(summary_m10.4_rstanarm[2,3],summary_m10.4_rstanarm[3,3],summary_m10.4_rstanarm[4,3],summary_m10.4_rstanarm[5,3],summary_m10.4_rstanarm[6,3],summary_m10.4_rstanarm[7,3],summary_m10.4_rstanarm[8,3],summary_m10.4_rstanarm[9,3],2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5))#, autoscale = FALSE
    my_prior_intercept_integrated<- normal(location = c(summary_m10.4_rstanarm[1,1]), scale = c(summary_m10.4_rstanarm[1,3]))#, autoscale = FALSE
    #my_prior_intercept_random_current<- normal(location = c(0), scale = c(10))#, autoscale = FALSE
      #See https://mc-stan.org/rstanarm/reference/priors.html#examples  
      #Covariance
      #~ decov(reg. = 1, conc. = 1, shape = 1, scale = 1)  
  
  #9.6.2 We write the model  
    #Time 37045.3 seconds (Total)  
    #mod_glmer<-glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_10 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_10_c  + LC_1 + LC_3 + LC_4 + LC_7 + ene_Wild_rep_plant+ene_Wild_invertebrates+ene_Wild_unk_plant_oth+ene_Wild_vertebrates +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial")
    #summary(mod_glmer)
    mod_stan_integrated2<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag + ene_Wild_rep_plant + ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = my_prior_integrated, prior_intercept = my_prior_intercept_integrated)
    #N=  39852 obs..#  25845.1 seconds (Total) with Frag
    save(mod_stan_integrated2,file="mod_stan_integrated2.Rdata")

    
  #9.6.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_integrated2)
    #A summary of the model  
    summary_mod_stan_integrated2<- as.data.frame(summary(mod_stan_integrated2))
    save (summary_mod_stan_integrated2, file='summary_mod_stan_integrated2.RData')
    write.csv(summary_mod_stan_integrated2, file='summary_mod_stan_integrated2.csv')# Supplementary Table 44 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 44#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 44. Results for the Hiarerchical Bayesian model using abiotic and biotic factors to explain brown bear distribution Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated2 [1:31,]

        data.frame(rownames(summary_mod_stan_integrated2))
        
        print(paste0("OUTPUT: Supplementary Table 43 mean_PPD Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated2 [33,]

  #9.6.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_integrated2)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated2, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)

  #9.6.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_integrated2)
       
  #9.6.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_integrated2)
    str(posterior_df)
    cor_poterior_mod_stan_integrated2<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_integrated2, file='cor_poterior_mod_stan_integrated2.xlsx')
    print(paste0("OUTPUT: Supplementary Table 49#  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 49. Correlation of the posterior samples among the predictors used in the Bayesian hierarchical model using abiotic and biotic factors to explain brown bear distribution  Model SHM_ABI~~~~~"))
    print(cor_poterior_mod_stan_integrated2)

  #9.6.5 We write a bayesian null model 
  print("This was calculated on section 9.5.5")
  
  #9.6.6 Visualization/Checking the chain
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_integrated2)
    dim(posterior)
    str(posterior)
    print(paste0("OUTPUT: Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 8.Chains for the Bayesian model of brown bear habitat with abiotic and biotic factors combining data.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    X11()
    plot(mod_stan_integrated2)  

  #9.6.7 Visualization of variable response    
    #We are gouing to calculate to see how is the response to changes for each variable.
    #We are going to change values of a variable meanwhile we mantein constant with a mean value the other variables
    #Prepare new counterfactual data R code 5.9
    load(file="dat_50km.Rdata")
    sub_variables<-dat_50km[,c(1:4,6:9)]
    head(sub_variables)
    nvariables<-length(sub_variables)

    Clim_3_median_h<-median(dat_50km$Clim_3)    
    Clim_4_median_h<-median(dat_50km$Clim_4) 
    Clim_8_median_h<-median(dat_50km$Clim_8)    
    Clim_9_median_h<-median(dat_50km$Clim_9)    

    
    Clim_3_median<-median(presences_absences_5km_scaled_only_current$Clim_3)    
    Clim_4_median<-median(presences_absences_5km_scaled_only_current$Clim_4) 
    Clim_8_median<-median(presences_absences_5km_scaled_only_current$Clim_8)    
    Clim_9_median<-median(presences_absences_5km_scaled_only_current$Clim_9)    

    LC_1_median<-median(presences_absences_5km_scaled_only_current$LC_1)
    LC_3_median<-median(presences_absences_5km_scaled_only_current$LC_3)
    LC_4_median<-median(presences_absences_5km_scaled_only_current$LC_4)
    Frag_median<-median(presences_absences_5km_scaled_only_current$Frag)

    ene_Wild_rep_plant_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_rep_plant)
    ene_Wild_invertebrates_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_invertebrates)
    ene_Wild_unk_plant_oth_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_unk_plant_oth)
    ene_Wild_vertebrates_median<-median(presences_absences_5km_scaled_only_current$ene_Wild_vertebrates)

    print(paste0("OUTPUT: Figure 5 MAIN TEXT & Supplementary Figure 7C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Fig. 5. Partial response plots of brown bear distribution to both abiotic and biotic variables. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    #X11() # PDF 8 x 7 inches and 83.4%
    par(mfrow=c(3,4))
    par(mar=c(2,2,2,0.5))
    
    all_lines<-c()
    all_shades<-c()
    
    for (var in 1:12){
      
            pred.data_h<- data.frame(
            Clim_3=Clim_3_median_h,
            Clim_4=Clim_4_median_h,
            Clim_8=Clim_8_median_h,
            Clim_9=Clim_9_median_h,
            Clim_3_c=Clim_3_median_h^2,  
            Clim_4_c=Clim_4_median_h^2,
            Clim_8_c=Clim_8_median_h^2,  
            Clim_9_c=Clim_9_median_h^2,
            LC_1=LC_1_median,
            LC_3=LC_3_median,
            LC_4=LC_4_median,
            Frag=Frag_median,
            ene_Wild_rep_plant=ene_Wild_rep_plant_median,
            ene_Wild_unk_plant_oth=ene_Wild_unk_plant_oth_median,
            ene_Wild_invertebrates =ene_Wild_invertebrates_median,
            ene_Wild_vertebrates =ene_Wild_vertebrates_median,
            factor_subpop=as.integer(1)
          )
    
  
      #We first use constant data for all variables/predictors  
      pred.data.current<- data.frame(
            Clim_3=Clim_3_median,
            Clim_4=Clim_4_median,
            Clim_8=Clim_8_median,
            Clim_9=Clim_9_median,
            Clim_3_c=Clim_3_median^2,  
            Clim_4_c=Clim_4_median^2,
            Clim_8_c=Clim_8_median^2,  
            Clim_9_c=Clim_9_median^2,
            LC_1=LC_1_median,
            LC_3=LC_3_median,
            LC_4=LC_4_median,
            Frag=Frag_median,
            ene_Wild_rep_plant=ene_Wild_rep_plant_median,
            ene_Wild_unk_plant_oth=ene_Wild_unk_plant_oth_median,
            ene_Wild_invertebrates =ene_Wild_invertebrates_median,
            ene_Wild_vertebrates =ene_Wild_vertebrates_median,
            factor_subpop=as.integer(1)
          )
      variables_names<-colnames(pred.data.current)
      variables_names_h<-colnames(pred.data_h)
        #We create a dataframe with 1000 rows repeating the data
        pred.data.current<-pred.data.current[rep(seq_len(nrow(pred.data.current)), each = 1000), ]
        #We create a dataframe with 1000 rows repeating the data
        pred.data_h<-pred.data_h[rep(seq_len(nrow(pred.data_h)), each = 1000), ]
        
     #Now for each variable (var) we are going to substitute the values of pred.data.current with the median
      # by a sequence of values  
        if(var<5){ 
          var<-var
          print(paste("Variable in loop ", var))
          var_in_loop<-pred.data_h[,c(var)]  
          var_in_loop_c<-pred.data_h[,c(var+4)] 
      
          var_name_in_loop<-variables_names_h[c(var)]  
          var_name_in_loop_c<-variables_names_h[c(var+4)] 
          print(paste("Linear Variable in loop    ", var_name_in_loop))
          print(paste("Quadratic Variable in loop ", var_name_in_loop_c))
          
          #If we plot using the range of values of the historical range:
            #var_in_loop_seq<-seq(from=min(dat_50km[,var_name_in_loop]), to=max(dat_50km[,var_name_in_loop]), length.out = 1000)
          #If we plot using the range of values of current range:
            var_in_loop_seq<-seq(from=min(dat_50km[,var_name_in_loop]), to=max(dat_50km[,var_name_in_loop]), length.out = 1000)
          #We calculate the sequence of values for the quadratic effect:  
            var_c_in_loop_seq<-var_in_loop_seq^2
          #We sustitute the sequence of values for the variables in the data frame with the values to predict
          pred.data_h[,c(var)]<-var_in_loop_seq
          pred.data_h[,c(var+4)]<-var_c_in_loop_seq
          plot( presence ~ dat_50km[,c(var)],  data=dat_50km , col=col.alpha(rangi2,0.5), cex=0,pch = 20, main=(var_name_in_loop))
        }  
        if(var>4){ 
          var<-var+4
          print(paste("Variable in loop ", var))
          var_in_loop<-pred.data.current[,c(var)] 
          
          var_name_in_loop<-variables_names[c(var)]  
          print(paste("Variable >4 in loop    ", var_name_in_loop))
  
          var_in_loop_seq<-seq(from=min(presences_absences_5km_scaled_only_current[,var_name_in_loop]), to=max(presences_absences_5km_scaled_only_current[,var_name_in_loop]), length.out = 1000)
          pred.data.current[,c(var)]<-var_in_loop_seq
          plot( Europe_presences_0_1 ~ presences_absences_5km_scaled_only_current[,var_name_in_loop],  data=presences_absences_5km_scaled_only_current , col="blue", cex=0,pch = 20, main=(var_name_in_loop))
          pred.data_h<-pred.data.current
          }  
        #We extract the posterior
        stan.mu <- posterior_epred(mod_stan_integrated2, newdata = pred.data_h)
    # summarize the distribution of mu
    mu.mean <- apply( stan.mu , 2 , mean )
    mu.HPDI <- apply( stan.mu , 2 , HPDI , prob=0.95 )
    ## R code 4.57
    # plot raw data
    # fading out points to make line and interval more visible
    # plot the MAP line, aka the mean mu for each weight
    rect(min(presences_absences_5km_scaled_only_current[,var_name_in_loop]),0,max(presences_absences_5km_scaled_only_current[,var_name_in_loop]),1,col = "lightblue1")
    shade( mu.HPDI , var_in_loop_seq ,col = "grey80")
    lines( var_in_loop_seq , mu.mean )
    
    points(x=presences_absences_5km_scaled_only_current[c(var_name_in_loop)][[1]],
               y=presences_absences_5km_scaled_only_current[c("Europe_presences_0_1")][[1]],
                pch = 19,

                col = alpha("red",0.002))
    
    }
    #Save as PDF 8 x 7 inches and 83.4%
    
    
  #9.6.8 Prediction Bayesian hierarchical
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)
    
    #We load the fitted integrated model
    load(file="mod_stan_integrated2.RData")
    dfmod<-(mod_stan_integrated2$glmod$fr)
    variables_model<-colnames(dfmod)[2:17]
 
    #See https://stats.stackexchange.com/questions/89172/how-to-scale-new-observations-for-making-predictions-when-the-model-was-fitted-w
    #We load the scaled variables (Climate variables are not scaled) and their attributes 
    load("scaled.var_model_uninformative_priors_only_current.RData")
    
    #9.6.8.1 Preparation of normal future scenarios (current was already prepared)
      #Common files to all scenarios:    
        #We load the data of coordinates and ID pixels
        #We load the data of coordinates and ID pixels
        load(paste0(my_dir,"/0_Construction_of_the_Spatial_Database/data_matrix.RData"))
        #We load the vector of extrapolation
        load(paste0(my_dir,"/4_SDM_Food_species/vec_extrap_raster.RData"))
        #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
        load(paste0(my_dir,"/4_SDM_Food_species/vec_subpop_raster.RData"))
      
        for (sc in 1:3){          
          #For RCP26
          if (sc==1){
            #We load the scenario current abiotic
            load("F:/G/Project_name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP26.RData")
            #We load the scenario current biotic/energy variables NO  energy standarizeD, calculated summer 2020
            load("F:/G/Project_name/Results_Biomod/R_analysis/Energy_variables/merged_DF_Variables_Energy_RCP26.RData")
            #We load the vector of cober
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_cober_future_2050_SSP1_RCP26.RData")
            #We load the vector of water
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_water_future_2050_SSP1_RCP26.RData")
            #We load the vector of Frag
            load("F:/G/Project_name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_gaus_pr_sel_natural_SSP1_RCP26.Rdata")
            #We change the names  
            #The name of the scenario
            scenario_name<-"RCP26"
            #For the data frame with the abiotic variables
            Scenario<-Scenario_RCP26
            #For the data frame with the biotic variables
            Energy<-merged_DF_Variables_Energy_RCP26
            #Variable cober
            cober<-vec_cober_future_2050_SSP1_RCP26
            #Variable water
            water<-vec_water_future_2050_SSP1_RCP26
            #Variable Frag
            Frag<-vec_gaus_pr_sel_natural_SSP1_RCP26
          }   
          #For RCP60
          if (sc==2){
            #We load the scenario current abiotic
            load("F:/G/Project_name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP60.RData")
            #We load the scenario current biotic/energy variables NO  energy standarizeD, calculated summer 2020
            load("F:/G/Project_name/Results_Biomod/R_analysis/Energy_variables/merged_DF_Variables_Energy_RCP60.RData")
            #We load the vector of cober
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_cober_future_2050_SSP3_RCP70.RData")
            #We load the vector of water
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_water_future_2050_SSP3_RCP70.RData")
            #We load the vector of Frag
            load("F:/G/Project_name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_gaus_pr_sel_natural_SSP3_RCP70.Rdata")
            #We change the names  
            #The name of the scenario
            scenario_name<-"RCP60"
            #For the data frame with the abiotic variables
            Scenario<-Scenario_RCP60
            #For the data frame with the biotic variables
            Energy<-merged_DF_Variables_Energy_RCP60
            #Variable cober
            cober<-vec_cober_future_2050_SSP3_RCP70
            #Variable water
            water<-vec_water_future_2050_SSP3_RCP70
           #Variable Frag
            Frag<-vec_gaus_pr_sel_natural_SSP3_RCP70
          }   
          #For RCP85
          if (sc==3){
            #We load the scenario current abiotic
            load("F:/G/Project_name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP85.RData")
            #We load the scenario current biotic/energy variables NO  energy standarizeD, calculated summer 2020
            load("F:/G/Project_name/Results_Biomod/R_analysis/Energy_variables/merged_DF_Variables_Energy_RCP85.RData")
            #We load the vector of cober
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_cober_future_2050_SSP5_RCP85.RData")
            #We load the vector of water
            load("C:/Users/Author/Documents/Bear/Analisis_t/vec_water_future_2050_SSP5_RCP85.RData")
            #We change the names  
            #The name of the scenario
            scenario_name<-"RCP85"
            #For the data frame with the abiotic variables
            Scenario<-Scenario_RCP85
            #For the data frame with the biotic variables
            Energy<-merged_DF_Variables_Energy_RCP85
            #Variable cober
            cober<-vec_cober_future_2050_SSP5_RCP85
            #Variable water
            water<-vec_water_future_2050_SSP5_RCP85
          }   
          #This code is common to the three scenarios    
          Energy <- Energy[ ,-c(1) ]
          #We merge all the data frames
              DF_Scen_all<-as.data.frame(cbind(data_matrix,
              Scenario,
              Energy,
              cober,
              water, 
              Frag,  
              vec_extrap_raster,
              vec_subpop_raster))
          #We select the areas inside the subpopulations
          DF_Scen_all_subpopul<-subset(DF_Scen_all,vec_subpop_raster >= 1 & vec_subpop_raster < 15)
          #We select the areas with data of land use   
          DF_Scen_all_subpopul2<-subset(DF_Scen_all_subpopul,LC_1 >= 0 & LC_1 < 2)
          DF_Scen_all_subpopul2$factor_subpop<-as.factor(DF_Scen_all_subpopul2$vec_subpop_raster)
          #We select the areas without water 
          DF_Scen_all_noNA<-subset(DF_Scen_all_subpopul2,water == 0 )

          #First we apply the conversions to make equal the data of CRU TS used for historical and the data of CHELSA used for current
          #We need to multiply CHELSA data bioclimatic variables from Bio1 to Bio11 (both included) by 10:
          #DF_Scen_all_noNA$Clim_3<-DF_Scen_all_noNA$Clim_3/10 #For some unknown reason Clim_was not in the same scale #We need to multiply each RCP$Clim_3 by 10 to be equal to original Scenario_current
          DF_Scen_all_noNA$Clim_4<-DF_Scen_all_noNA$Clim_4/10
          DF_Scen_all_noNA$Clim_8<-DF_Scen_all_noNA$Clim_8/10
          DF_Scen_all_noNA$Clim_9<-DF_Scen_all_noNA$Clim_9/10
          #Now we need to apply the same transformation to scalar moreless the data around 1
          DF_Scen_all_noNA$Clim_3<-DF_Scen_all_noNA$Clim_3/100
          DF_Scen_all_noNA$Clim_4<-DF_Scen_all_noNA$Clim_4/1000
          DF_Scen_all_noNA$Clim_8<-DF_Scen_all_noNA$Clim_8/10
          DF_Scen_all_noNA$Clim_9<-DF_Scen_all_noNA$Clim_9/10
      
          DF_Scen_all_noNA$Clim_3_c<-DF_Scen_all_noNA$Clim_3^2
          DF_Scen_all_noNA$Clim_4_c<-DF_Scen_all_noNA$Clim_4^2
          DF_Scen_all_noNA$Clim_8_c<-DF_Scen_all_noNA$Clim_8^2
          DF_Scen_all_noNA$Clim_9_c<-DF_Scen_all_noNA$Clim_9^2
          #We assign a name with the scenario to the data frame and we save it
          assign(paste0("DF_Scen_",scenario_name,"_all_noNA"),DF_Scen_all_noNA)
          save(list=paste0("DF_Scen_",scenario_name,"_all_noNA"), file=paste0("DF_Scen_",scenario_name,"_all_noNA.Rdata"))   
    }  

    #9.6.8.2 Scalation of normal scenarios (current and future)
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)
      #See https://stats.stackexchange.com/questions/89172/how-to-scale-new-observations-for-making-predictions-when-the-model-was-fitted-w
      #We load the scaled variables and their attributes 
      load("scaled.var_model_uninformative_priors_only_current.RData")

      #FOR EACH SCENARIO
      for (sc in 1:4){          
        #For current
        if (sc==1){
          #We load the scenario (This is the scenario with all Europe, not only the presence/pdeudoabsences)
          load("DF_Scen_current_all_noNA.RData")
          scenario_in_loop<-DF_Scen_current_all_noNA
          scenario_name<-"current"
        }
        #For RCP26
        if (sc==2){
          #We load the scenario 
          load("DF_Scen_RCP26_all_noNA.RData")
          scenario_in_loop<-DF_Scen_RCP26_all_noNA
          scenario_name<-"RCP26"
        }
        #For RCP60
        if (sc==3){
          #We load the scenario 
          load("DF_Scen_RCP60_all_noNA.RData")
          scenario_in_loop<-DF_Scen_RCP60_all_noNA
          scenario_name<-"RCP60"
        }
        #For RCP85
        if (sc==4){
          #We load the scenario 
          load("DF_Scen_RCP85_all_noNA.RData")
          scenario_in_loop<-DF_Scen_RCP85_all_noNA
          scenario_name<-"RCP85"
        }
        #This code is common to the four scenarios    
        #Create a subset of land use and biotic variables
        var_uninformative<-scenario_in_loop[,c("LC_1","LC_3","LC_4","Frag",
            "ene_Wild_rep_plant","ene_Wild_invertebrates","ene_Wild_unk_plant_oth","ene_Wild_vertebrates")] 
          #we are going to use the scaled attributes of the variables of the fitted model to scale the scenario  
          scaled.new <- scale(var_uninformative, attr(scaled.var_model_uninformative_priors_only_current, "scaled:center"), attr(scaled.var_model_uninformative_priors_only_current, "scaled:scale"))
          scaled.new_DF<-as.data.frame(scaled.new)
          #We create a dataframe with all variables the rescaled variables of land use and biotic variables and the climate variables
          var_clim_0<-scenario_in_loop[,c("ID_pixel","Clim_3","Clim_4","Clim_8","Clim_9",
            "Clim_3_c","Clim_4_c","Clim_8_c","Clim_9_c")]  
          factor_subpop<-scenario_in_loop$factor_subpop
          scaled_variables_for_prediction<-as.data.frame(cbind(var_clim_0,scaled.new_DF,factor_subpop))
          #We create an integer variable to divide all the
          scaled_variables_for_prediction$inte_row<-floor((c(1:nrow(scaled_variables_for_prediction))/100000))
          #We assign a name with the scenario to the data frame and we save it
          assign(paste0("scaled_variables_for_prediction_",scenario_name),scaled_variables_for_prediction)
          save(list=paste0("scaled_variables_for_prediction_", scenario_name), file=paste0("scaled_variables_for_prediction_", scenario_name, ".Rdata"))   
      }
    
      
    #9.6.8.3 Combination of scenarios to simulate scenario with only change in  biotic and scenario with only change in abiotic
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)

      load("scaled_variables_for_prediction_current.RData")
      current_abiotic<-scaled_variables_for_prediction_current[,c(1:13,18)]
      current_biotic<-scaled_variables_for_prediction_current[,c(1,14:18)]

      #For the RCP26
      load("scaled_variables_for_prediction_RCP26.RData")
      RCP26_abiotic<-scaled_variables_for_prediction_RCP26[,c(1:13)]
      RCP26_biotic<-scaled_variables_for_prediction_RCP26[,c(1,14:17)]

      #For the RCP60
      load("scaled_variables_for_prediction_RCP60.RData")
      RCP60_abiotic<-scaled_variables_for_prediction_RCP60[,c(1:13)]
      RCP60_biotic<-scaled_variables_for_prediction_RCP60[,c(1,14:17)]

      #For the RCP85
      load("scaled_variables_for_prediction_RCP85.RData")
      RCP85_abiotic<-scaled_variables_for_prediction_RCP85[,c(1:13)]
      RCP85_biotic<-scaled_variables_for_prediction_RCP85[,c(1,14:17)]

      #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
      RCP26_change_onlybiotic<-merge(current_abiotic,RCP26_biotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP26_change_onlybiotic, file='RCP26_change_onlybiotic.RData')
      
      RCP60_change_onlybiotic<-merge(current_abiotic,RCP60_biotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP60_change_onlybiotic, file='RCP60_change_onlybiotic.RData')
      
      RCP85_change_onlybiotic<-merge(current_abiotic,RCP85_biotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP85_change_onlybiotic, file='RCP85_change_onlybiotic.RData')
      
      #Scenario with only change in abiotic (We are goint to take the biotic in the current scenario)
      RCP26_change_onlyabiotic<-merge(current_biotic,RCP26_abiotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP26_change_onlyabiotic, file='RCP26_change_onlyabiotic.RData')
      
      RCP60_change_onlyabiotic<-merge(current_biotic,RCP60_abiotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP60_change_onlyabiotic, file='RCP60_change_onlyabiotic.RData')
      
      RCP85_change_onlyabiotic<-merge(current_biotic,RCP85_abiotic, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
      save(RCP85_change_onlyabiotic, file='RCP85_change_onlyabiotic.RData')

    #9.6.8.4 We calculate the prediction for each scenario
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)

      #We load the data of coordinates and ID pixels
      load(paste0(my_dir,"/0_Construction_of_the_Spatial_Database/data_matrix.RData"))
      #We load the fitted integrated model
      load(file="mod_stan_integrated2.RData")
      dfmod<-(mod_stan_integrated2$glmod$fr)
      variables_model<-colnames(dfmod)[2:17]

      #FOR EACH SCENARIO NORMAL AND SIMULATED
      for (sc in 1:10){   
        print(paste0("#######################Scenario in loop ",sc,"#######################################"))
        #For current
        if (sc==1){
        #For current
          #We load the scenario 
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          scenario_name<-"current"
        }
        #For RCP26
        if (sc==2){
          #We load the scenario 
          load("scaled_variables_for_prediction_RCP26.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP26
          scenario_name<-"RCP26"
        }
        #For RCP60
        if (sc==3){
          #We load the scenario 
          load("scaled_variables_for_prediction_RCP60.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP60
          scenario_name<-"RCP60"
        }
        #For RCP85
        if (sc==4){
          #We load the scenario 
          load("scaled_variables_for_prediction_RCP85.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP85
          scenario_name<-"RCP85"
        }
        ############SIMULATED SCENARIOS################
        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==5){
          #We load the scenario 
          load("RCP26_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlybiotic
          scenario_name<-"RCP26_change_onlybiotic"
        }
        #For RCP60
        if (sc==6){
          #We load the scenario 
          load("RCP60_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlybiotic
          scenario_name<-"RCP60_change_onlybiotic"
        }
        #For RCP85
        if (sc==7){
          #We load the scenario 
          load("RCP85_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlybiotic
          scenario_name<-"RCP85_change_onlybiotic"
        }

        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==8){
          #We load the scenario 
          load("RCP26_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlyabiotic
          scenario_name<-"RCP26_change_onlyabiotic"
        }
        #For RCP60
        if (sc==9){
          #We load the scenario 
          load("RCP60_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlyabiotic
          scenario_name<-"RCP60_change_onlyabiotic"
        }
        #For RCP85
        if (sc==10){
          #We load the scenario 
          load("RCP85_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlyabiotic
          scenario_name<-"RCP85_change_onlyabiotic"
        }
 
        if (sc>4){
         scaled_variables_for_prediction$inte_row<-floor((c(1:nrow(scaled_variables_for_prediction))/100000))
        }

        max_integer<-max(scaled_variables_for_prediction$inte_row)
        df_pred_integr_all<-c()
        for (i in 0:max_integer){
          print(paste0("Loop i = ",i," #####"))
          sub_in_loop<-scaled_variables_for_prediction[scaled_variables_for_prediction$inte_row==i,]  
          #We extract the posterior
          link.mod_stan_integrated2 <- posterior_epred(mod_stan_integrated2, newdata = sub_in_loop)
          # summarize, we calculate the mean of all the columns of the matrix
          sum_mean_integrated<-apply( link.mod_stan_integrated2 , 2 , mean)
          #We compute the percentile intervals, we select prob=0.95
          sum_PI_integrated <- apply( link.mod_stan_integrated2 , 2 , HPDI , prob=0.95 )
    
          df_sum_mean_integrated<-as.data.frame(sum_mean_integrated)
          df_sum_PI_integrated<-as.data.frame(sum_PI_integrated)
          df_sum_PI_t_integrated<-t(df_sum_PI_integrated)  
          df_sum_PI_t_integrated<-as.data.frame(df_sum_PI_t_integrated)
    
          df_pred_integr<-cbind(df_sum_mean_integrated,df_sum_PI_t_integrated)
          df_pred_integr$Uncertainity<-df_pred_integr[,c(3)]-df_pred_integr[,c(2)]
          df_pred_integr<-cbind(df_pred_integr,sub_in_loop$ID_pixel)
          colnames(df_pred_integr)<-c("sum_mean","|0.95","0.95|","Uncertainity","ID")
          df_pred_integr_all<-rbind(df_pred_integr_all,df_pred_integr)
        }
        #This is the prediction incuding different factor for each subpopulation
        merge_prediciton_integrated<-merge(data_matrix,df_pred_integr_all, by.x="ID_pixel", by.y = "ID",all.x = T)
        assign(paste0("merge_prediciton_integrated_",scenario_name),merge_prediciton_integrated)
        save(list=paste0("merge_prediciton_integrated_", scenario_name), file=paste0("merge_prediciton_integrated_", scenario_name, ".Rdata"))   

        #9.6.8.5 Creation of the rasters with the prediction of habitat and uncertainity
        #We create rasters in img format with the mean values of the predicion and the range of the confidence intervals  
        #We change the R version, close and open R in earlier version (3.5.3) because a strange error:
          
        #We save into a raster file to see in Arcgis
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
        #We are going to use Europe Albers Equal Area Conic  
        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
        #We defined the projection of our raster:
        projection(new_raster) <- newproj
    
        #For the habitat
          raster_habitat_integrated<-new_raster
          values(raster_habitat_integrated)<-round(merge_prediciton_integrated$sum_mean,digits=4)
          assign(paste0("raster_habitat_integrated_",scenario_name),raster_habitat_integrated)
          save(list=paste0("raster_habitat_integrated_", scenario_name), file=paste0("raster_habitat_integrated_", scenario_name, ".Rdata"))   
          print(paste0("OUTPUT: Figure 4 & Electronic Material 1-9, Raster files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          writeRaster(raster_habitat_integrated, paste0("raster_habitat_integrated_", scenario_name, ".img"), overwrite=TRUE) ##################################################################################### Supplementary Electronic Material 1-9

        #For the uncertainity
          raster_uncertainity_integrated<-new_raster
          values(raster_uncertainity_integrated)<-round(merge_prediciton_integrated$Uncertainity,digits=5)
          assign(paste0("raster_uncertainity_integrated_",scenario_name),raster_uncertainity_integrated)
          save(list=paste0("raster_uncertainity_integrated_", scenario_name), file=paste0("raster_uncertainity_integrated_", scenario_name, ".Rdata"))   
          print(paste0("OUTPUT: Electronic Material 10-18, Raster files with the uncertainity of the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          writeRaster(raster_uncertainity_integrated, paste0("raster_uncertainity_integrated_", scenario_name, ".img"), overwrite=TRUE) ############################################################################ Supplementary Electronic Material 10-18

        #We save it in a kml file to see in google earth      
          #We crate a new rster i lat long: 
          new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
          crs(new_raster_2) <- CRS('+init=EPSG:4326')
          raster_habitat_integrated_lat_long <- projectRaster(raster_habitat_integrated, new_raster_2, method='bilinear')
          print(paste0("OUTPUT: Electronic Material 19-27, Klm files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          KML(raster_habitat_integrated_lat_long, file=paste0("raster_habitat_integrated_kml_", scenario_name, ".kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)######################### Supplementary Electronic Material 19-27
      }

  #9.6.9 Validation of the Bayesian hierarchical
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    load(paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_validation_all_5km.Rdata"))
    load("merge_prediciton_integrated_current.RData")
    load("RANDOM_presences_absences_for_validation_all_5km.RData")
    presences_absences_5km<-RANDOM_presences_absences_for_validation_all_5km
    merge_prediciton_integrated_observed<-merge(presences_absences_5km,merge_prediciton_integrated_current, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
    merge_prediciton_integrated_observed2<-merge_prediciton_integrated_observed[,c("Europe_presences_0_1","factor_subpop","sum_mean","|0.95","0.95|","Uncertainity")]         
    colnames(merge_prediciton_integrated_observed2)<-c("Observed","factor_subpop","Predicted","|0.95","0.95|","Uncertainity")
    save (merge_prediciton_integrated_observed2, file='merge_prediciton_integrated_observed2.RData')
    head(merge_prediciton_integrated_observed2)
    
    #9.6.9.1 Validation of the model selecting a threshold including the 90% of presences
      #We binarize the predictions into 0 and 1. We are going to select a threshold to include the 90% of the presences
      sub_DF_predicted_observed_of_training<-merge_prediciton_integrated_observed2[merge_prediciton_integrated_observed2$Observed==1,]
      order_data<-sub_DF_predicted_observed_of_training[order(sub_DF_predicted_observed_of_training$Predicted),]
      quantiles<-data.frame(quantile(order_data$Predicted, probs = c(0, 0.1, 0.9, 1))) 
      colnames(quantiles)<-c("thresholds") 
      quantiles
      Observed_predicted_50<-merge_prediciton_integrated_observed2
      Observed_predicted_50$code_PREDIC = ifelse(Observed_predicted_50$Predicted < quantiles["10%",],0,1)
  
      tableGLMM<- table(Observed_predicted_50$Observed,Observed_predicted_50$code_PREDIC)
      dftableGLMM<-as.data.frame(tableGLMM)
      dftableGLMM$Measure<-c("True negative","False negative","False positive","True positive")
      save (tableGLMM, file='tableGLMM.RData')
      save (dftableGLMM, file='dftableGLMM.RData')
  
      #SPECIFITY is the proportion of true-negatives
      Specifity<-dftableGLMM[1,3]/((dftableGLMM[1,3])+(dftableGLMM[3,3]))
      #SENSITIVITY is the proportion of true-positives
      Sensivity<-dftableGLMM[4,3]/((dftableGLMM[4,3])+(dftableGLMM[2,3]))
      #Accuracy
      Accuracy<-(dftableGLMM[4,3]+dftableGLMM[1,3])/((dftableGLMM[1,3])+(dftableGLMM[2,3])+(dftableGLMM[3,3])+(dftableGLMM[4,3]))
  
      #By subpopulations
      validation_subp_all<-c()
      for (subp in 1:14){
        print(paste0("SUBPOPULATION IN LOOP ",subp))
        sub_pop_Observed_predicted_50<-subset(Observed_predicted_50, factor_subpop==subp)
      
        tableGLMM_subp<- table(sub_pop_Observed_predicted_50$Observed,sub_pop_Observed_predicted_50$code_PREDIC)
        dftableGLMM_subp<-as.data.frame(tableGLMM_subp)
        dftableGLMM_subp$Measure<-c("True negative","False negative","False positive","True positive")
        
        #SPECIFITY is the proportion of true-negatives
        Specifity_subp<-dftableGLMM_subp[1,3]/((dftableGLMM_subp[1,3])+(dftableGLMM_subp[3,3]))
        Specifity_subp
        #SENSITIVITY is the proportion of true-positives
        Sensivity_subp<-dftableGLMM_subp[4,3]/((dftableGLMM_subp[4,3])+(dftableGLMM_subp[2,3]))
        Sensivity_subp
        #Accuracy
        Accuracy_subp<-(dftableGLMM_subp[4,3]+dftableGLMM_subp[1,3])/((dftableGLMM_subp[1,3])+(dftableGLMM_subp[2,3])+(dftableGLMM_subp[3,3])+(dftableGLMM_subp[4,3]))
        Accuracy_subp
        validation_subp<-c(Specifity_subp,Sensivity_subp,Accuracy_subp)
        validation_subp_all<-cbind(validation_subp_all,validation_subp)  
      }
      rownames(validation_subp_all)<-c("Specifity_subp","Sensivity_subp","Accuracy_subp")
      colnames(validation_subp_all)<-c(1:14)
      validation<-as.data.frame(c(Specifity,Sensivity,Accuracy))
      rownames(validation)<-c("Specifity","Sensivity","Accuracy")
      colnames(validation)<-c("All study area")
      print(paste0("OUTPUT: Supplementary Table 54. Results for the validation of the best Bayesian model (BM) explaining brown bear distribution, the Bayesian hierarchical model using abiotic and biotic factors #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 54. We show the values to correctly classify the pseudo-absences of brown bear (true negative rate; TNR), to correctly classify the presences of brown bear (true positive rate; TPR) and classification accuracy (Acc.) at European scale and by subpopulation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      validation
      validation_subp_all
      save (validation, file='validation.RData')
      save (validation_subp_all, file='validation_subp_all.RData')
      write.xlsx(validation, file='validation.xlsx')
      write.xlsx(validation_subp_all, file='validation_subp_all.xlsx')
   
    #9.6.9.2 Validation of the model selecting a threshold which maxize the TSS
      library(ecospat)
      
      #9.6.9.2.1 A unique threshold for all Europe
      tss_threshold_max<-ecospat.max.tss(merge_prediciton_integrated_observed2$Predicted, merge_prediciton_integrated_observed2$Observed)       
      write.xlsx(tss_threshold_max, file='tss_threshold_max.xlsx')
      
      #9.6.9.2.2 A different threshold by subpopulation
      tss_threshold_max_subp_all<-c()
      for (subp in 1:14){
        print(paste0("SUBPOPULATION IN LOOP ",subp))
        sub_pop<-subset(merge_prediciton_integrated_observed2, factor_subpop==subp)
        sub_pop_tss_threshold_max<-ecospat.max.tss(sub_pop$Predicted, sub_pop$Observed)       
        tss_in_loop<-c(sub_pop_tss_threshold_max$max.TSS,sub_pop_tss_threshold_max$max.threshold)
        tss_threshold_max_subp_all<-rbind(tss_threshold_max_subp_all,tss_in_loop)  
      }
      colnames(tss_threshold_max_subp_all)<-c("max.TSS","max.threshold")
      tss_threshold_max_subp_all_df<-as.data.frame(tss_threshold_max_subp_all)
      summary(tss_threshold_max_subp_all_df)
      mean(tss_threshold_max_subp_all[,1])
      table(merge_prediciton_integrated_observed2$factor_subpop)
      write.xlsx(tss_threshold_max, file='tss_threshold_max.xlsx')

      #9.6.9.2.3 We binarize the predictions into 0 and 1. We are going to select a threshold that maximize the TSS (A unique threshold for all Europe obtained in #9.6.9.2.1)
      sub_DF_predicted_observed_of_training2<-merge_prediciton_integrated_observed2[merge_prediciton_integrated_observed2$Observed==1,]
      order_data2<-sub_DF_predicted_observed_of_training2[order(sub_DF_predicted_observed_of_training2$Predicted),]
      quantiles<-data.frame(quantile(order_data2$Predicted, probs = c(0, 0.1, 0.9, 1))) 
      colnames(quantiles)<-c("thresholds") 
      quantiles
      Observed_predicted_maxTSS<-merge_prediciton_integrated_observed2
      Observed_predicted_maxTSS$code_PREDIC = ifelse(Observed_predicted_maxTSS$Predicted < tss_threshold_max$max.threshold,0,1)
  
      tableGLMM_maxTSS<- table(Observed_predicted_maxTSS$Observed,Observed_predicted_maxTSS$code_PREDIC)
      dftableGLMM_maxTSS<-as.data.frame(tableGLMM_maxTSS)
      dftableGLMM_maxTSS$Measure<-c("True negative","False negative","False positive","True positive")
      save (tableGLMM_maxTSS, file='tableGLMM_maxTSS.RData')
      save (dftableGLMM_maxTSS, file='dftableGLMM_maxTSS.RData')
  
      #SPECIFITY is the proportion of true-negatives
      Specifity<-dftableGLMM_maxTSS[1,3]/((dftableGLMM_maxTSS[1,3])+(dftableGLMM_maxTSS[3,3]))
      #SENSITIVITY is the proportion of true-positives
      Sensivity<-dftableGLMM_maxTSS[4,3]/((dftableGLMM_maxTSS[4,3])+(dftableGLMM_maxTSS[2,3]))
      #Accuracy
      Accuracy<-(dftableGLMM_maxTSS[4,3]+dftableGLMM_maxTSS[1,3])/((dftableGLMM_maxTSS[1,3])+(dftableGLMM_maxTSS[2,3])+(dftableGLMM_maxTSS[3,3])+(dftableGLMM_maxTSS[4,3]))
  
      #By subpopulations
      validation_subp_all_maxTSS<-c()
      for (subp in 1:14){
        print(paste0("SUBPOPULATION IN LOOP ",subp))
        sub_pop_Observed_predicted_50<-subset(Observed_predicted_50, factor_subpop==subp)
      
        tableGLMM_subp<- table(sub_pop_Observed_predicted_50$Observed,sub_pop_Observed_predicted_50$code_PREDIC)
        dftableGLMM_maxTSS_subp<-as.data.frame(tableGLMM_subp)
        dftableGLMM_maxTSS_subp$Measure<-c("True negative","False negative","False positive","True positive")
        
        #SPECIFITY is the proportion of true-negatives
        Specifity_subp<-dftableGLMM_maxTSS_subp[1,3]/((dftableGLMM_maxTSS_subp[1,3])+(dftableGLMM_maxTSS_subp[3,3]))
        Specifity_subp
        #SENSITIVITY is the proportion of true-positives
        Sensivity_subp<-dftableGLMM_maxTSS_subp[4,3]/((dftableGLMM_maxTSS_subp[4,3])+(dftableGLMM_maxTSS_subp[2,3]))
        Sensivity_subp
        #Accuracy
        Accuracy_subp<-(dftableGLMM_maxTSS_subp[4,3]+dftableGLMM_maxTSS_subp[1,3])/((dftableGLMM_maxTSS_subp[1,3])+(dftableGLMM_maxTSS_subp[2,3])+(dftableGLMM_maxTSS_subp[3,3])+(dftableGLMM_maxTSS_subp[4,3]))
        Accuracy_subp
        validation_subp<-c(Specifity_subp,Sensivity_subp,Accuracy_subp)
        validation_subp_all_maxTSS<-cbind(validation_subp_all_maxTSS,validation_subp)  
      }
      rownames(validation_subp_all_maxTSS)<-c("Specifity_subp","Sensivity_subp","Accuracy_subp")
      colnames(validation_subp_all_maxTSS)<-c(1:14)
      validation_maxTSS<-as.data.frame(c(Specifity,Sensivity,Accuracy))
      rownames(validation_maxTSS)<-c("Specifity","Sensivity","Accuracy")
      colnames(validation_maxTSS)<-c("All study area")
      validation_maxTSS
      validation_subp_all_maxTSS
      save (validation_maxTSS, file='validation_maxTSS.RData')
      save (validation_subp_all_maxTSS, file='validation_subp_all_maxTSS.RData')
      write.xlsx(validation_maxTSS, file='validation_maxTSS.xlsx')
      write.xlsx(validation_subp_all_maxTSS, file='validation_subp_all_maxTSS.xlsx')

    #9.6.10 We asses the predictions for each scenario
      #9.6.10.1 We calculate the suitable habitat fragments and a variable to calculate the distance in idrisi
        rm(list=ls()) 
        my_dir <-"writehereyourpath" 
        folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
        setwd(folder_working)

        code_all<-c()
        code_group_all<-c()
        #FOR EACH SCENARIO NORMAL AND SIMULATED
        for (sc in 1:10){   
          print(paste0("#######################Scenario in loop ",sc,"#######################################"))
          #For current
          if (sc==1){
          #For current
            #We load the scenario 
            load(file="raster_habitat_integrated_current.RData")
            #load(file="raster_habitat_integrated_current.RData")
            raster_habitat_integrated<-raster_habitat_integrated_current
            scenario_name<-"current"
          }
          #For RCP26
          if (sc==2){
            #We load the scenario 
            load("raster_habitat_integrated_RCP26.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP26
            scenario_name<-"RCP26"
          }
          #For RCP60
          if (sc==3){
            #We load the scenario 
            load("raster_habitat_integrated_RCP60.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP60
            scenario_name<-"RCP60"
          }
          #For RCP85
          if (sc==4){
            #We load the scenario 
            load("raster_habitat_integrated_RCP85.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP85
            scenario_name<-"RCP85"
          }
          ############SIMULATED SCENARIOS################
          #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
          #For RCP26
          if (sc==5){
            #We load the scenario 
            load("raster_habitat_integrated_RCP26_change_onlybiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlybiotic
            scenario_name<-"RCP26_change_onlybiotic"
          }
          #For RCP60
          if (sc==6){
            #We load the scenario 
            load("raster_habitat_integrated_RCP60_change_onlybiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlybiotic
            scenario_name<-"RCP60_change_onlybiotic"
          }
          #For RCP85
          if (sc==7){
            #We load the scenario 
            load("raster_habitat_integrated_RCP85_change_onlybiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlybiotic
            scenario_name<-"RCP85_change_onlybiotic"
          }
  
          #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
          #For RCP26
          if (sc==8){
            #We load the scenario 
            load("raster_habitat_integrated_RCP26_change_onlyabiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlyabiotic
            scenario_name<-"RCP26_change_onlyabiotic"
          }
          #For RCP60
          if (sc==9){
            #We load the scenario 
            load("raster_habitat_integrated_RCP60_change_onlyabiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlyabiotic
            scenario_name<-"RCP60_change_onlyabiotic"
          }
          #For RCP85
          if (sc==10){
            #We load the scenario 
            load("raster_habitat_integrated_RCP85_change_onlyabiotic.RData")
            raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlyabiotic
            scenario_name<-"RCP85_change_onlyabiotic"
          }
          
         #We create a map of the habitat in binary
         raster_habitat_integrated_binary<-raster_habitat_integrated
         #We reclass human to integer values to calculate water areas
         raster_habitat_integrated_binary[raster_habitat_integrated_binary >= 0.5 & raster_habitat_integrated_binary <= 1 ] <- 1
         raster_habitat_integrated_binary[raster_habitat_integrated_binary < 0.5 ] <- 0
         writeRaster(raster_habitat_integrated_binary, paste0("raster_habitat_integrated_binary_", scenario_name, ".rst"), datatype='INT4S', overwrite=TRUE)
         #We calculate the patches     
         #patches0<- get_patches(raster_habitat_integrated_binary)
         #We select the patches of suitable habitat (1)     
         #patches<-patches0$`1`
         #assign(paste0("raster_patches_",scenario_name),patches)
         #save(list=paste0("raster_patches_", scenario_name), file=paste0("raster_patches_", scenario_name, ".Rdata"))  
         #raster_NO_habitat_binary<-raster_habitat_integrated
         #raster_NO_habitat_binary[raster_NO_habitat_binary >= 0.5 & raster_NO_habitat_binary <= 1 ] <- 0
         #raster_NO_habitat_binary[raster_NO_habitat_binary > 0 & raster_NO_habitat_binary < 0.5 ] <- 1
         #writeRaster(raster_NO_habitat_binary, paste0("raster_NO_habitat_binary_", scenario_name, ".rst"), datatype='INT4S', overwrite=TRUE)
         #Go to idrisi and calculate distance to raster_NO_habitat_binary
         #code<-paste0("distance x ",paste0(folder_working,"/raster_NO_habitat_binary_", scenario_name, ".rst"),"*",paste0(folder_working,"/raster_NO_habitat_binary_", scenario_name, "_dis.rst"))
         #code_all<-rbind(code_all,code)
        }
       #macro_idrisi_distance<-as.data.frame(code_all)
       #write.csv(file="macro_idrisi_distance.csv", x=macro_idrisi_distance)

    #9.6.11.2 We prepare the data of current range in our study area and of protected areas
      #9.6.11.2.1 current range in our study area 
        rm(list=ls()) 
        my_dir <-"writehereyourpath" 
        folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
        setwd(folder_working)
        
          mapa_bear <- shapefile('C:/Users/Author/Documents/Bear/Current_distribution/data_0.shp')
          sub_pres_cur = subset(mapa_bear, PRESENCE=="1")
          #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
              new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
              #We are going to use Europe Albers Equal Area Conic  
              newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
              #We defined the projection of our raster:
              projection(new_raster) <- newproj
          #We project the shapefile into CEA:
           sub_pres_cur_aea <- spTransform(sub_pres_cur, newproj) 
          #Rasterize the data
          sub_pres_cur_raster<-rasterize(sub_pres_cur_aea, new_raster, field=1)# field="presence"
          #We extract the values of the raster to a vector:
          vec_range_current_brown_bear_1km<-values(sub_pres_cur_raster)
          save(vec_range_current_brown_bear_1km,file="vec_range_current_brown_bear_1km.RData")

      #9.6.9.2.2 Protected areas
        rm(list=ls()) 
        my_dir <-"writehereyourpath" 
        folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
        setwd(folder_working)
 
          #1.1.2 The protected areas
          protected_areas<-raster("PROTECTED_AREAS_TER_COAST_POINT_WORLD_REC.rst")
          projection(protected_areas) <- newproj
        
          #Raster of reference: 
          cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" 
          #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data
          projection(protected_areas) <- cea_proj

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
          #We are going to cut/crop the map for reduce to Europe
          #We stablish an extension for crop the map:
          e <- extent(-2000000,7500000,3000000,6000000)#Extension que utlizamos de normal 
          #We use the extension to crop the original map    
          m1_crop <- crop(protected_areas,e) 
          m1_crop[m1_crop >= 0.5 & m1_crop <= 2 ] <- 1
          m1_crop[m1_crop > 0 & m1_crop < 0.5 ] <- 0
          #We project our data to the raster created with the new projection
          protected_areas_aea <- projectRaster(m1_crop, new_raster, method='bilinear')
          protected_areas_aea[protected_areas_aea >= 0.5 & protected_areas_aea <= 2 ] <- 1
          protected_areas_aea[protected_areas_aea > -2 & protected_areas_aea < 0.5 ] <- 0
          #We extract the values of the raster to a vector:
          vec_protected_areas<-values(protected_areas_aea)
          save(vec_protected_areas,file="vec_protected_areas.RData")
          X11()
          plot(protected_areas_aea)

    #9.6.11.3 We calculate the suitable habitat measures for the habitat predictions
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)
      
      load(paste0(my_dir,"/0_Construction_of_the_Spatial_Database/data_matrix.RData"))     
      load("vec_protected_areas.RData")      
      load("vec_range_current_brown_bear_1km.RData")      
      #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
      load(paste0(my_dir,"/4_SDM_Food_species/vec_subpop_raster.RData"))

      CV=function(x) (sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))

      measure_in_europe_all<-c()
      measure_in_loop_all<-c()
      #FOR EACH SCENARIO NORMAL AND SIMULATED
      for (sc in 1:10){   
        print(paste0("#######################Scenario in loop ",sc,"#######################################"))
        #For current
        if (sc==1){
        #For current
          #We load the scenario 
          scenario_name<-"current"
          load(file="raster_habitat_integrated_current.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_current
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_current.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_current
          #raster_dis<-raster("raster_NO_habitat_binary_current_dis.rst")
        }
        #For RCP26
        if (sc==2){
          #We load the scenario 
          scenario_name<-"RCP26"
          load("raster_habitat_integrated_RCP26.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_RCP26
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP26
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_dis.rst")
        }
        #For RCP60
        if (sc==3){
          #We load the scenario 
          scenario_name<-"RCP60"
          load("raster_habitat_integrated_RCP60.RData")
          load("scaled_variables_for_prediction_RCP60.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP60
          raster_habitat_integrated<-raster_habitat_integrated_RCP60
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP60
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_dis.rst")
        }
        #For RCP85
        if (sc==4){
          #We load the scenario 
          scenario_name<-"RCP85"
          load("raster_habitat_integrated_RCP85.RData")
          load("scaled_variables_for_prediction_RCP85.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP85
          raster_habitat_integrated<-raster_habitat_integrated_RCP85
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP85#26058
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_dis.rst")
        }
        ############SIMULATED SCENARIOS################
        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==5){
          #We load the scenario 
          scenario_name<-"RCP26_change_onlybiotic"
          load("raster_habitat_integrated_RCP26_change_onlybiotic.RData")
          load("RCP26_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlybiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP26_change_onlybiotic
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlybiotic_dis.rst")
        }
        #For RCP60
        if (sc==6){
          #We load the scenario 
          scenario_name<-"RCP60_change_onlybiotic"
          load("raster_habitat_integrated_RCP60_change_onlybiotic.RData")
          load("RCP60_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlybiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP60_change_onlybiotic
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlybiotic_dis.rst")
        }
        #For RCP85
        if (sc==7){
          #We load the scenario 
          scenario_name<-"RCP85_change_onlybiotic"
          load("raster_habitat_integrated_RCP85_change_onlybiotic.RData")
          load("RCP85_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlybiotic.rst")
          #raster_patches<-raster("raster_patches_RCP85_change_onlybiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP85_change_onlybiotic
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlybiotic_dis.rst")
        }

        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==8){
          #We load the scenario 
          scenario_name<-"RCP26_change_onlyabiotic"
          load("raster_habitat_integrated_RCP26_change_onlyabiotic.RData")
          load("RCP26_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP26_change_onlyabiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP26_change_onlyabiotic
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlyabiotic_dis.rst")
        }
        #For RCP60
        if (sc==9){
          #We load the scenario 
          scenario_name<-"RCP60_change_onlyabiotic"
          load("raster_habitat_integrated_RCP60_change_onlyabiotic.RData")
          load("RCP60_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP60_change_onlyabiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP60_change_onlyabiotic
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlyabiotic_dis.rst")
        }
        #For RCP85
        if (sc==10){
          #We load the scenario 
          scenario_name<-"RCP85_change_onlyabiotic"
          load("raster_habitat_integrated_RCP85_change_onlyabiotic.RData")
          load("RCP85_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP85_change_onlyabiotic.rst")
          #load(paste0("raster_patches_", scenario_name, ".Rdata"))#raster_patches_R_file<-
          #raster_patches<-raster_patches_RCP85_change_onlyabiotic#42492
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlyabiotic_dis.rst")
        }

        if (sc>4){
          scaled_variables_for_prediction$inte_row<-floor((c(1:nrow(scaled_variables_for_prediction))/100000))
        }
        vector_habitat<-values(raster_habitat_integrated)
        vector_habitat_binary<-values(raster_habitat_integrated_binary)
        #table(vector_habitat_binary)
        #vector_patches<-values(raster_patches)
        #vector_dis<-values(raster_dis)
    
        #DF_descriptors_variables<-as.data.frame(cbind(data_matrix,vector_habitat,vector_habitat_binary,
        #  vector_patches,vector_dis,vec_range_current_brown_bear_1km,vec_protected_areas))
       DF_descriptors_variables<-as.data.frame(cbind(data_matrix,vector_habitat,vector_habitat_binary,
          vec_range_current_brown_bear_1km,vec_protected_areas))
         
        merge_descriptors_variables0<-merge(DF_descriptors_variables,scaled_variables_for_prediction, by.x="ID_pixel", by.y = "ID_pixel",all.y = T)
        #merge_descriptors_variables<-merge_descriptors_variables0[,c("ID_pixel","x","y","vector_habitat","vector_habitat_binary","vector_patches","vector_dis","vec_range_current_brown_bear_1km","vec_protected_areas","factor_subpop")]
        merge_descriptors_variables<-merge_descriptors_variables0[,c("ID_pixel","x","y","vector_habitat","vector_habitat_binary","vec_range_current_brown_bear_1km","vec_protected_areas","factor_subpop")]
        merge_descriptors_variables_habitable<-subset(merge_descriptors_variables,vector_habitat_binary==1)

        #We calculate the area of the pathes
        #area_in_pixels_vector_patches0<- as.data.frame(table(merge_descriptors_variables_habitable$vector_patches))
        #colnames(area_in_pixels_vector_patches0)<-c("Patch_ID","Area")
        #area_in_pixels_vector_patches0$Patch_ID_num<-as.numeric(as.character(area_in_pixels_vector_patches0$Patch_ID))
        #colnames(area_in_pixels_vector_patches0)<-c("Patch_ID","Area","Patch_ID_num")
        #area_in_pixels_vector_patches<-subset(area_in_pixels_vector_patches0,Patch_ID_num>=0)
        #Measures of habitat to include
          #Area of habitat
          area_habitat<- sum((vector_habitat_binary==1))
          #Number of fragments
          #n_frag<-length(area_in_pixels_vector_patches$Patch_ID)
          #CV area of patches
          #cv_area_patches<-CV(area_in_pixels_vector_patches$Area)
          #Mean area by fragment
          #mean_Area_frag<-area_habitat/n_frag
          #Median area by fragment
          #median_Area_frag_vec<-median(area_in_pixels_vector_patches$Area)
          #Max area of fragment
          #max_Area_frag_vec<-max(area_in_pixels_vector_patches$Area)
          #Heterogeneity. Proportion of the total range area represented by the largest fragment.
          #It ranges from close to 0 (similiar fragment size) to close to 1 (very
          #different fragment size).
          #heterogeneity<-max_Area_frag_vec/area_habitat 
          #Distance to the border measures
          #Mean distance to the border
          #mean_distance<-mean(merge_descriptors_variables_habitable$vector_dis)
          #Median distance to the border
          #median_distance<-median(merge_descriptors_variables_habitable$vector_dis)
          #Core habitat
          #core_habitat_dis_big2km<-subset(merge_descriptors_variables_habitable,vector_dis>2000)
          #area_core_habitat_dis_big2km<- sum(core_habitat_dis_big2km$vector_habitat_binary)
        #Mesures for fragmetns bigger than a fragment for have a stable population size
          #subset_area_bigger_than3000km2<-subset(area_in_pixels_vector_patches,Area>3000)
            #Area of habitat
            #area_habitat_bigger_than3000km2<- sum(subset_area_bigger_than3000km2$Area)
            #Number of fragments
            #n_frag_vec_bigger_than3000km2<-length(subset_area_bigger_than3000km2$Patch_ID)
        #Measures of distribution and protected areas  
          #Habitat in protected areas   
          subset_habitat_protected<-subset(merge_descriptors_variables_habitable,vec_protected_areas>0.5)
            #Area of habitat protected
            area_habitat_protected_km2<- sum(subset_habitat_protected$vector_habitat_binary)
            percentage_in_protected_areas<-area_habitat_protected_km2/area_habitat*100
          #Occupied habitat
          subset_habitat_occupied<-subset(merge_descriptors_variables_habitable,vec_range_current_brown_bear_1km>0.5)
            #Area of habitat Occupied
            area_habitat_occupied_km2<- sum(subset_habitat_occupied$vector_habitat_binary)
            percentage_occupied<-area_habitat_occupied_km2/area_habitat*100
        #We bind the data at European level of the sc loop    
        measure_in_europe<-cbind(
          sc,
          area_habitat,
          #n_frag,
          #cv_area_patches,
          #mean_Area_frag,
          #median_Area_frag_vec,
          #max_Area_frag_vec,
          #heterogeneity,
          #mean_distance,
          #median_distance,
          #area_core_habitat_dis_big2km,
          #area_habitat_bigger_than3000km2,
          #n_frag_vec_bigger_than3000km2,
          area_habitat_protected_km2,
          percentage_in_protected_areas,
          area_habitat_occupied_km2,
          percentage_occupied)
    
        #We bind the data at European level of the sc loop  to the data of all sc  
        measure_in_europe_all<-as.data.frame(rbind(measure_in_europe_all,measure_in_europe))
    #subpop=1
        for (subpop in 1:14){
          merge_descriptors_variables_habitable_in_loop<-subset(merge_descriptors_variables_habitable,factor_subpop==subpop) 
          if(nrow(merge_descriptors_variables_habitable_in_loop)>=1){
          #We calculate the area of the pathes
          #area_in_pixels_vector_patches0_in_loop<- as.data.frame(table(merge_descriptors_variables_habitable_in_loop$vector_patches))
          #colnames(area_in_pixels_vector_patches0_in_loop)<-c("Patch_ID","Area")
          #area_in_pixels_vector_patches0_in_loop$Patch_ID_num<-as.numeric(as.character(area_in_pixels_vector_patches0_in_loop$Patch_ID))
          #colnames(area_in_pixels_vector_patches0_in_loop)<-c("Patch_ID","Area","Patch_ID_num")
          #area_in_pixels_vector_patches_in_loop<-subset(area_in_pixels_vector_patches0_in_loop,Patch_ID_num>=0)
          #Measures of habitat to include
            #Area of habitat
            area_habitat_in_loop<- sum(merge_descriptors_variables_habitable_in_loop$vector_habitat_binary==1)
            #Number of fragments
            #n_frag_in_loop<-length(area_in_pixels_vector_patches_in_loop$Patch_ID)
            #CV area of patches
            #cv_area_patches_in_loop<-CV(area_in_pixels_vector_patches_in_loop$Area)
            #Mean area by fragment
            #mean_Area_frag_in_loop<-area_habitat_in_loop/n_frag_in_loop
            #Median area by fragment
            #median_Area_frag_vec_in_loop<-median(area_in_pixels_vector_patches_in_loop$Area)
            #Max area of fragment
            #max_Area_frag_vec_in_loop<-max(area_in_pixels_vector_patches_in_loop$Area)
            #Heterogeneity. Proportion of the total range area represented by the largest fragment.
            #It ranges from close to 0 (similiar fragment size) to close to 1 (very
            #different fragment size).
            #heterogeneity_in_loop<-max_Area_frag_vec_in_loop/area_habitat_in_loop 
            #Distance to the border measures
            #Mean distance to the border
            #mean_distance_in_loop<-mean(merge_descriptors_variables_habitable_in_loop$vector_dis)
            #Median distance to the border
            #median_distance_in_loop<-median(merge_descriptors_variables_habitable_in_loop$vector_dis)
            #Core habitat
            #core_habitat_dis_big2km_in_loop<-subset(merge_descriptors_variables_habitable_in_loop,vector_dis>2000)
            #area_core_habitat_dis_big2km_in_loop<- sum(core_habitat_dis_big2km_in_loop$vector_habitat_binary)
          #Mesures for fragmetns bigger than a fragment for have a stable population size
            #subset_area_bigger_than3000km2_in_loop<-subset(area_in_pixels_vector_patches_in_loop,Area>3000)
              #Area of habitat
              #area_habitat_bigger_than3000km2_in_loop<- sum(subset_area_bigger_than3000km2_in_loop$Area)
              #Number of fragments
              #n_frag_vec_bigger_than3000km2_in_loop<-length(subset_area_bigger_than3000km2_in_loop$Patch_ID)
          #Measures of distribution and protected areas  
            #Habitat in protected areas   
            subset_habitat_protected_in_loop<-subset(merge_descriptors_variables_habitable_in_loop,vec_protected_areas>0.5)
              #Area of habitat protected
              area_habitat_protected_km2_in_loop<- sum(subset_habitat_protected_in_loop$vector_habitat_binary)
              percentage_in_protected_areas_in_loop<-area_habitat_protected_km2_in_loop/area_habitat_in_loop*100
            #Occupied habitat
            subset_habitat_occupied_in_loop<-subset(merge_descriptors_variables_habitable_in_loop,vec_range_current_brown_bear_1km>0.5)
              #Area of habitat Occupied
              area_habitat_occupied_km2_in_loop<- sum(subset_habitat_occupied_in_loop$vector_habitat_binary)
              percentage_occupied_in_loop<-area_habitat_occupied_km2_in_loop/area_habitat_in_loop*100
          #We bind the data at subpopulation level of the sc loop    
            measure_in_loop<-cbind(
              sc,
              subpop,  
              area_habitat_in_loop,
              #n_frag_in_loop,
              #cv_area_patches_in_loop,
              #mean_Area_frag_in_loop,
              #median_Area_frag_vec_in_loop,
              #max_Area_frag_vec_in_loop,
              #heterogeneity_in_loop,
              #mean_distance_in_loop,
              #median_distance_in_loop,
              #area_core_habitat_dis_big2km_in_loop,
              #area_habitat_bigger_than3000km2_in_loop,
              #n_frag_vec_bigger_than3000km2_in_loop,
              area_habitat_protected_km2_in_loop,
              percentage_in_protected_areas_in_loop,
              area_habitat_occupied_km2_in_loop,
              percentage_occupied_in_loop)
          }else{
           measure_in_loop<-c(sc,subpop,0,NA,NA,NA)  
          }
            #We bind the data at subpopulation level of the sc loop with all data at subpopulation from all sc   
            measure_in_loop_all<-as.data.frame(rbind(measure_in_loop_all,measure_in_loop))
          }
      }
      save (measure_in_europe_all, file='measure_in_europe_all.RData')
      write.xlsx(measure_in_europe_all, file='measure_in_europe_all.xlsx')
      print(paste0("OUTPUT: Supplementary Table 55 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 55.Current potential habitat area, % of occupation, %protected areas in Europe~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      measure_in_europe_all
      
      save (measure_in_loop_all, file='measure_in_loop_all.RData')
      write.xlsx(measure_in_loop_all, file='measure_in_loop_all.xlsx')
      print(paste0("OUTPUT: Supplementary Table 56 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 56.Current potential habitat area, % of occupation, %protected areas by subpopulation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      measure_in_loop_all

      #9.6.11.3.1 We are going to asses the change in the variables from the current scenario to future  
        rm(list=ls()) 
        my_dir <-"writehereyourpath" 
        folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
        setwd(folder_working)

        load("measure_in_europe_all.RData")      
        measure_in_europe_current<-measure_in_europe_all[1,]
    
        percentage_change_sc_all_europe<-c()
        real_sc_europe_all<-c()
        for (sc in 1:10){
          measure_in_europe_sc<-measure_in_europe_all[sc,]  
          percentage_change_sc_europe<-measure_in_europe_sc/measure_in_europe_current*100
          percentage_change_sc_europe$sc<-percentage_change_sc_europe$sc/100
          percentage_change_sc_all_europe<-as.data.frame(rbind(percentage_change_sc_all_europe, percentage_change_sc_europe[,2:6]))
          real_sc_europe<-as.data.frame(cbind("sce"=sc, "subpopulation"=0,"subpopulation_scenario"=paste0("Subp_0_Sc_",sc)))
          real_sc_europe_all<-as.data.frame(rbind(real_sc_europe_all, real_sc_europe))
        }
         
        df_percentage_change_sc_europe<-as.data.frame(cbind(real_sc_europe_all,percentage_change_sc_all_europe))
        save (df_percentage_change_sc_europe, file='df_percentage_change_sc_europe.RData')
        write.xlsx(df_percentage_change_sc_europe, file='df_percentage_change_sc_europe.xlsx')
        print(paste0("OUTPUT: Supplementary Table 57 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 57.Change of current potential habitat area, % of occupation, %protected areas in Europe~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        df_percentage_change_sc_europe
 
      #9.6.11.3.2 We are going to asses the change in the variables from the current scenario to future for each subpopulation  
        measure_in_loop_all_current<-measure_in_loop_all[1:14,]
        percentage_change_sc_all<-c()
        real_sc_subpop_all<-c()
        for (subpopulation in 1:14){
           measure_in_loop_sce_subpop<-subset(measure_in_loop_all,subpop==subpopulation)
           measure_in_loop_current_subpop<-subset(measure_in_loop_all_current,subpop==subpopulation)
           for (sce in 1:10){
           measure_in_loop_sce<-subset(measure_in_loop_sce_subpop,sc==sce)
           percentage_change_sc<-measure_in_loop_sce/measure_in_loop_current_subpop*100
           percentage_change_sc_all<-as.data.frame(rbind(percentage_change_sc_all, percentage_change_sc[,3:7]))
           real_sc_subpop<-as.data.frame(cbind(sce, subpopulation,"subpopulation_scenario"=paste0("Subp_",subpopulation,"_Sc_",sce)))
           real_sc_subpop_all<-as.data.frame(rbind(real_sc_subpop_all, real_sc_subpop))
         }
       }
       df_percentage_change_sc_subpopulations<-as.data.frame(cbind(real_sc_subpop_all,percentage_change_sc_all))
       head(df_percentage_change_sc_subpopulations,14)
       save (df_percentage_change_sc_subpopulations, file='df_percentage_change_sc_subpopulations.RData')
        
       #write.xlsx(df_percentage_change_sc_subpopulations, file='df_percentage_change_sc_subpopulations.xlsx')
       str(df_percentage_change_sc_europe)
       str(df_percentage_change_sc_subpopulations) 
       variables_nam<-colnames(df_percentage_change_sc_europe) 
       colnames(df_percentage_change_sc_subpopulations)<-variables_nam 
            write.xlsx(df_percentage_change_sc_subpopulations, file='df_percentage_change_sc_subpopulations.xlsx')

       df_percentage_change_sc_europe_subpopulations<-as.data.frame(rbind(df_percentage_change_sc_europe,df_percentage_change_sc_subpopulations))
       save (df_percentage_change_sc_europe_subpopulations, file='df_percentage_change_sc_europe_subpopulations.RData')
       load("df_percentage_change_sc_europe_subpopulations.RData")     
       save (measure_in_loop_all, file='measure_in_loop_all.RData')
       print(paste0("OUTPUT: Supplementary Table 58 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
       print(paste0("Supplementary Table 58. Change from current potential habitat area, % of occupation, %protected areas by subpopulation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
       df_percentage_change_sc_europe_subpopulations
        write.xlsx(df_percentage_change_sc_europe_subpopulations, file='df_percentage_change_sc_europe_subpopulations.xlsx')
       
    #9.6.11.4 We calculate the suitable habitat measures for the currently occupied area
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)
      
      load("data_matrix.RData")      
      load("vec_protected_areas.RData")      
      load("vec_range_current_brown_bear_1km.RData")      
      #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
      load("vec_subpop_raster.RData")
      CV=function(x) (sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))

      measure_in_europe_all<-c()
      measure_in_loop_all<-c()
      #FOR EACH SCENARIO NORMAL AND SIMULATED
      for (sc in 1:10){   
        print(paste0("#######################Scenario in loop ",sc,"#######################################"))
        #For current
        if (sc==1){
        #For current
          #We load the scenario 
          load(file="raster_habitat_integrated_current.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_current
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_current.rst")
          raster_patches<-raster("raster_patches_current.rst")
          raster_dis<-raster("raster_NO_habitat_binary_current_dis.rst")
          scenario_name<-"current"
        }
        #For RCP26
        if (sc==2){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_RCP26
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26.rst")
          raster_patches<-raster("raster_patches_RCP26.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP26_dis.rst")
          scenario_name<-"RCP26"
        }
        #For RCP60
        if (sc==3){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60.RData")
          load("scaled_variables_for_prediction_RCP60.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP60
          raster_habitat_integrated<-raster_habitat_integrated_RCP60
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60.rst")
          raster_patches<-raster("raster_patches_RCP60.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP60_dis.rst")
          scenario_name<-"RCP60"
        }
        #For RCP85
        if (sc==4){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85.RData")
          load("scaled_variables_for_prediction_RCP85.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP85
          raster_habitat_integrated<-raster_habitat_integrated_RCP85
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85.rst")
          raster_patches<-raster("raster_patches_RCP85.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP85_dis.rst")
          scenario_name<-"RCP85"
        }
        ############SIMULATED SCENARIOS################
        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==5){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26_change_onlybiotic.RData")
          load("RCP26_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlybiotic.rst")
          raster_patches<-raster("raster_patches_RCP26_change_onlybiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlybiotic_dis.rst")
          scenario_name<-"RCP26_change_onlybiotic"
        }
        #For RCP60
        if (sc==6){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60_change_onlybiotic.RData")
          load("RCP60_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlybiotic.rst")
          raster_patches<-raster("raster_patches_RCP60_change_onlybiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlybiotic_dis.rst")
          scenario_name<-"RCP60_change_onlybiotic"
        }
        #For RCP85
        if (sc==7){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85_change_onlybiotic.RData")
          load("RCP85_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlybiotic.rst")
          raster_patches<-raster("raster_patches_RCP85_change_onlybiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlybiotic_dis.rst")
          scenario_name<-"RCP85_change_onlybiotic"
        }

        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==8){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26_change_onlyabiotic.RData")
          load("RCP26_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlyabiotic.rst")
          raster_patches<-raster("raster_patches_RCP26_change_onlyabiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP26_change_onlyabiotic"
        }
        #For RCP60
        if (sc==9){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60_change_onlyabiotic.RData")
          load("RCP60_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlyabiotic.rst")
          raster_patches<-raster("raster_patches_RCP60_change_onlyabiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP60_change_onlyabiotic"
        }
        #For RCP85
        if (sc==10){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85_change_onlyabiotic.RData")
          load("RCP85_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlyabiotic.rst")
          raster_patches<-raster("raster_patches_RCP85_change_onlyabiotic.rst")
          raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP85_change_onlyabiotic"
        }

        if (sc>4){
          scaled_variables_for_prediction$inte_row<-floor((c(1:nrow(scaled_variables_for_prediction))/100000))
        }
        vector_habitat<-values(raster_habitat_integrated)
        vector_habitat_binary<-values(raster_habitat_integrated_binary)
        vector_patches<-values(raster_patches)
        vector_dis<-values(raster_dis)
    
        DF_descriptors_variables<-as.data.frame(cbind(data_matrix,vector_habitat,vector_habitat_binary,
          vector_patches,vector_dis,vec_range_current_brown_bear_1km,vec_protected_areas))
        merge_descriptors_variables0<-merge(DF_descriptors_variables,scaled_variables_for_prediction, by.x="ID_pixel", by.y = "ID_pixel",all.y = T)
        merge_descriptors_variables<-merge_descriptors_variables0[,c("ID_pixel","x","y","vector_habitat","vector_habitat_binary","vector_patches","vector_dis","vec_range_current_brown_bear_1km","vec_protected_areas","factor_subpop")]
        
        #We select the area of current habitat
        merge_descriptors_variables_habitable<-subset(merge_descriptors_variables,vector_habitat_binary==1)
        #We select the area of current habitat which is occupied
        merge_descriptors_variables_habitable_occupied<-subset(merge_descriptors_variables_habitable,vec_range_current_brown_bear_1km==1)

        #We calculate the area of the pathes
        area_in_pixels_vector_patches0<- as.data.frame(table(merge_descriptors_variables_habitable_occupied$vector_patches))
        colnames(area_in_pixels_vector_patches0)<-c("Patch_ID","Area")
        area_in_pixels_vector_patches0$Patch_ID_num<-as.numeric(as.character(area_in_pixels_vector_patches0$Patch_ID))
        colnames(area_in_pixels_vector_patches0)<-c("Patch_ID","Area","Patch_ID_num")
        area_in_pixels_vector_patches<-subset(area_in_pixels_vector_patches0,Patch_ID_num>=0)
        #Measures of habitat to include
          #Area of habitat
          area_habitat<- sum(area_in_pixels_vector_patches$Area)
          #Number of fragments
          n_frag<-length(area_in_pixels_vector_patches$Patch_ID)
          #CV area of patches
          cv_area_patches<-CV(area_in_pixels_vector_patches$Area)
          #Mean area by fragment
          mean_Area_frag<-area_habitat/n_frag
          #Median area by fragment
          median_Area_frag_vec<-median(area_in_pixels_vector_patches$Area)
          #Max area of fragment
          max_Area_frag_vec<-max(area_in_pixels_vector_patches$Area)
          #Heterogeneity. Proportion of the total range area represented by the largest fragment.
          #It ranges from close to 0 (similiar fragment size) to close to 1 (very
          #different fragment size).
          heterogeneity<-max_Area_frag_vec/area_habitat 
          #Distance to the border measures
          #Mean distance to the border
          mean_distance<-mean(merge_descriptors_variables_habitable_occupied$vector_dis)
          #Median distance to the border
          median_distance<-median(merge_descriptors_variables_habitable_occupied$vector_dis)
          #Core habitat
          core_habitat_dis_big2km<-subset(merge_descriptors_variables_habitable_occupied,vector_dis>2000)
          area_core_habitat_dis_big2km<- sum(core_habitat_dis_big2km$vector_habitat_binary)
        #Mesures for fragmetns bigger than a fragment for have a stable population size
          subset_area_bigger_than3000km2<-subset(area_in_pixels_vector_patches,Area>3000)
            #Area of habitat
            area_habitat_bigger_than3000km2<- sum(subset_area_bigger_than3000km2$Area)
            #Number of fragments
            n_frag_vec_bigger_than3000km2<-length(subset_area_bigger_than3000km2$Patch_ID)
        #Measures of distribution and protected areas  
          #Habitat in protected areas   
          subset_habitat_protected<-subset(merge_descriptors_variables_habitable_occupied,vec_protected_areas>0.5)
            #Area of habitat protected
            area_habitat_protected_km2<- sum(subset_habitat_protected$vector_habitat_binary)
            percentage_in_protected_areas<-area_habitat_protected_km2/area_habitat*100
          #Occupied habitat
          subset_habitat_occupied<-subset(merge_descriptors_variables_habitable_occupied,vec_range_current_brown_bear_1km>0.5)
            #Area of habitat Occupied
            area_habitat_occupied_km2<- sum(subset_habitat_occupied$vector_habitat_binary)
            percentage_occupied<-area_habitat_occupied_km2/area_habitat*100
        #We bind the data at European level of the sc loop    
        measure_in_europe<-cbind(
          sc,
          area_habitat,
          n_frag,
          cv_area_patches,
          mean_Area_frag,
          median_Area_frag_vec,
          max_Area_frag_vec,
          heterogeneity,
          mean_distance,
          median_distance,
          area_core_habitat_dis_big2km,
          area_habitat_bigger_than3000km2,
          n_frag_vec_bigger_than3000km2,
          area_habitat_protected_km2,
          percentage_in_protected_areas,
          area_habitat_occupied_km2,
          percentage_occupied)
    
        #We bind the data at European level of the sc loop  to the data of all sc  
        measure_in_europe_all<-as.data.frame(rbind(measure_in_europe_all,measure_in_europe))
    
        for (subpop in 1:14){
          merge_descriptors_variables_habitable_occupied_in_loop<-subset(merge_descriptors_variables_habitable_occupied,factor_subpop==subpop) 
          if(nrow(merge_descriptors_variables_habitable_occupied_in_loop)>=1){
          #We calculate the area of the pathes
          area_in_pixels_vector_patches0_in_loop<- as.data.frame(table(merge_descriptors_variables_habitable_occupied_in_loop$vector_patches))
          colnames(area_in_pixels_vector_patches0_in_loop)<-c("Patch_ID","Area")
          area_in_pixels_vector_patches0_in_loop$Patch_ID_num<-as.numeric(as.character(area_in_pixels_vector_patches0_in_loop$Patch_ID))
          colnames(area_in_pixels_vector_patches0_in_loop)<-c("Patch_ID","Area","Patch_ID_num")
          area_in_pixels_vector_patches_in_loop<-subset(area_in_pixels_vector_patches0_in_loop,Patch_ID_num>=0)
          #Measures of habitat to include
            #Area of habitat
            area_habitat_in_loop<- sum(area_in_pixels_vector_patches_in_loop$Area)
            #Number of fragments
            n_frag_in_loop<-length(area_in_pixels_vector_patches_in_loop$Patch_ID)
            #CV area of patches
            cv_area_patches_in_loop<-CV(area_in_pixels_vector_patches_in_loop$Area)
            #Mean area by fragment
            mean_Area_frag_in_loop<-area_habitat_in_loop/n_frag_in_loop
            #Median area by fragment
            median_Area_frag_vec_in_loop<-median(area_in_pixels_vector_patches_in_loop$Area)
            #Max area of fragment
            max_Area_frag_vec_in_loop<-max(area_in_pixels_vector_patches_in_loop$Area)
            #Heterogeneity. Proportion of the total range area represented by the largest fragment.
            #It ranges from close to 0 (similiar fragment size) to close to 1 (very
            #different fragment size).
            heterogeneity_in_loop<-max_Area_frag_vec_in_loop/area_habitat_in_loop 
            #Distance to the border measures
            #Mean distance to the border
            mean_distance_in_loop<-mean(merge_descriptors_variables_habitable_occupied_in_loop$vector_dis)
            #Median distance to the border
            median_distance_in_loop<-median(merge_descriptors_variables_habitable_occupied_in_loop$vector_dis)
            #Core habitat
            core_habitat_dis_big2km_in_loop<-subset(merge_descriptors_variables_habitable_occupied_in_loop,vector_dis>2000)
            area_core_habitat_dis_big2km_in_loop<- sum(core_habitat_dis_big2km_in_loop$vector_habitat_binary)
          #Mesures for fragmetns bigger than a fragment for have a stable population size
            subset_area_bigger_than3000km2_in_loop<-subset(area_in_pixels_vector_patches_in_loop,Area>3000)
              #Area of habitat
              area_habitat_bigger_than3000km2_in_loop<- sum(subset_area_bigger_than3000km2_in_loop$Area)
              #Number of fragments
              n_frag_vec_bigger_than3000km2_in_loop<-length(subset_area_bigger_than3000km2_in_loop$Patch_ID)
          #Measures of distribution and protected areas  
            #Habitat in protected areas   
            subset_habitat_protected_in_loop<-subset(merge_descriptors_variables_habitable_occupied_in_loop,vec_protected_areas>0.5)
              #Area of habitat protected
              area_habitat_protected_km2_in_loop<- sum(subset_habitat_protected_in_loop$vector_habitat_binary)
              percentage_in_protected_areas_in_loop<-area_habitat_protected_km2_in_loop/area_habitat_in_loop*100
            #Occupied habitat
            subset_habitat_occupied_in_loop<-subset(merge_descriptors_variables_habitable_occupied_in_loop,vec_range_current_brown_bear_1km>0.5)
              #Area of habitat Occupied
              area_habitat_occupied_km2_in_loop<- sum(subset_habitat_occupied_in_loop$vector_habitat_binary)
              percentage_occupied_in_loop<-area_habitat_occupied_km2_in_loop/area_habitat_in_loop*100
          #We bind the data at subpopulation level of the sc loop    
            measure_in_loop<-cbind(
              sc,
              subpop,  
              area_habitat_in_loop,
              n_frag_in_loop,
              cv_area_patches_in_loop,
              mean_Area_frag_in_loop,
              median_Area_frag_vec_in_loop,
              max_Area_frag_vec_in_loop,
              heterogeneity_in_loop,
              mean_distance_in_loop,
              median_distance_in_loop,
              area_core_habitat_dis_big2km_in_loop,
              area_habitat_bigger_than3000km2_in_loop,
              n_frag_vec_bigger_than3000km2_in_loop,
              area_habitat_protected_km2_in_loop,
              percentage_in_protected_areas_in_loop,
              area_habitat_occupied_km2_in_loop,
              percentage_occupied_in_loop)
          }else{
           measure_in_loop<-c(sc,subpop,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)  
          }
            #We bind the data at subpopulation level of the sc loop with all data at subpopulation from all sc   
            measure_in_loop_all<-as.data.frame(rbind(measure_in_loop_all,measure_in_loop))
          }
      }
    measure_in_europe_all_habitat_occupied<-measure_in_europe_all
    measure_in_loop_all_habitat_occupied<-measure_in_loop_all
    save (measure_in_europe_all_habitat_occupied, file='measure_in_europe_all_habitat_occupied.RData')
    save (measure_in_loop_all_habitat_occupied, file='measure_in_loop_all_habitat_occupied.RData')

     #9.6.11.4.2 We are going to asses the change in the variables from the current scenario to future for each subpopulation  
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)

      load("measure_in_europe_all_habitat_occupied.RData")  
      
       measure_in_europe_current<-measure_in_europe_all_habitat_occupied[1,]
       percentage_change_sc_all_europe<-c()
       real_sc_europe_all<-c()
       for (sc in 1:10){
       measure_in_europe_sc<-measure_in_europe_all_habitat_occupied[sc,]  
       percentage_change_sc_europe<-measure_in_europe_sc/measure_in_europe_current*100
       percentage_change_sc_europe$sc<-percentage_change_sc_europe$sc/100
       percentage_change_sc_all_europe<-as.data.frame(rbind(percentage_change_sc_all_europe, percentage_change_sc_europe[,2:17]))
       real_sc_europe<-as.data.frame(cbind("sce"=sc, "subpopulation"=0,"subpopulation_scenario"=paste0("Subp_0_Sc_",sc)))
       real_sc_europe_all<-as.data.frame(rbind(real_sc_europe_all, real_sc_europe))
       }
       df_percentage_change_sc_europe_habitat_occupied<-as.data.frame(cbind(real_sc_europe_all,percentage_change_sc_all_europe))
       save (df_percentage_change_sc_europe_habitat_occupied, file='df_percentage_change_sc_europe_habitat_occupied.RData')
   
      #9.6.11.4.1 We are going to asses the change in the variables from the current scenario to future  
        load("measure_in_loop_all_habitat_occupied.RData")      
        measure_in_loop_all_current<-measure_in_loop_all_habitat_occupied[1:14,]
        percentage_change_sc_all<-c()
        real_sc_subpop_all<-c()

       for (subpopulation in 1:14){
         measure_in_loop_sce_subpop<-subset(measure_in_loop_all_habitat_occupied,subpop==subpopulation)
         measure_in_loop_current_subpop<-subset(measure_in_loop_all_current,subpop==subpopulation)
         for (sce in 1:10){
         measure_in_loop_sce<-subset(measure_in_loop_sce_subpop,sc==sce)
         percentage_change_sc<-measure_in_loop_sce/measure_in_loop_current_subpop*100
         percentage_change_sc_all<-as.data.frame(rbind(percentage_change_sc_all, percentage_change_sc[,3:18]))
         real_sc_subpop<-as.data.frame(cbind(sce, subpopulation,"subpopulation_scenario"=paste0("Subp_",subpopulation,"_Sc_",sce)))
         real_sc_subpop_all<-as.data.frame(rbind(real_sc_subpop_all, real_sc_subpop))
        }
        }
       df_percentage_change_sc_subpopulations_habitat_occupied<-as.data.frame(cbind(real_sc_subpop_all,percentage_change_sc_all))
       save (df_percentage_change_sc_subpopulations_habitat_occupied, file='df_percentage_change_sc_subpopulations_habitat_occupied.RData')
       variables_nam<-colnames(df_percentage_change_sc_europe_habitat_occupied) 
       colnames(df_percentage_change_sc_subpopulations_habitat_occupied)<-variables_nam 
       df_percentage_change_sc_europe_subpopulations_habitat_occupied<-as.data.frame(rbind(df_percentage_change_sc_europe_habitat_occupied,df_percentage_change_sc_subpopulations_habitat_occupied))
       save (df_percentage_change_sc_europe_subpopulations_habitat_occupied, file='df_percentage_change_sc_europe_subpopulations_habitat_occupied.RData')
       df_percentage_change_sc_europe_subpopulations_current_RCP65_habitat_occupied<-subset(df_percentage_change_sc_europe_subpopulations_habitat_occupied, sce==1|sce==3|sce==6|sce==9)

    #9.6.11.5 We calculate rasters of habitat change for each scenario
      rm(list=ls()) 
      my_dir <-"writehereyourpath" 
      folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
      setwd(folder_working)
 
      #load("data_matrix.RData")      
      #load("vec_protected_areas.RData")      
      #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
      #load("vec_subpop_raster.RData")
      #CV=function(x) (sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
      raster_habitat_integrated_binary_current<-raster("raster_habitat_integrated_binary_current.rst")
      #measure_in_europe_all<-c()
      #measure_in_loop_all<-c()
      #FOR EACH SCENARIO NORMAL AND SIMULATED
      for (sc in 2:10){   
        print(paste0("#######################Scenario in loop ",sc,"#######################################"))
        #For current
        if (sc==1){
        #For current
          #We load the scenario 
          load(file="raster_habitat_integrated_current.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_current
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_current.rst")
          #raster_patches<-raster("raster_patches_current.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_current_dis.rst")
          scenario_name<-"current"
        }
        #For RCP26
        if (sc==2){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26.RData")
          load("scaled_variables_for_prediction_current.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_current
          raster_habitat_integrated<-raster_habitat_integrated_RCP26
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26.rst")
          #raster_patches<-raster("raster_patches_RCP26.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_dis.rst")
          scenario_name<-"RCP26"
        }
        #For RCP60
        if (sc==3){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60.RData")
          load("scaled_variables_for_prediction_RCP60.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP60
          raster_habitat_integrated<-raster_habitat_integrated_RCP60
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60.rst")
          #raster_patches<-raster("raster_patches_RCP60.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_dis.rst")
          scenario_name<-"RCP60"
        }
        #For RCP85
        if (sc==4){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85.RData")
          load("scaled_variables_for_prediction_RCP85.RData")
          scaled_variables_for_prediction<-scaled_variables_for_prediction_RCP85
          raster_habitat_integrated<-raster_habitat_integrated_RCP85
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85.rst")
          #raster_patches<-raster("raster_patches_RCP85.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_dis.rst")
          scenario_name<-"RCP85"
        }
        ############SIMULATED SCENARIOS################
        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==5){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26_change_onlybiotic.RData")
          load("RCP26_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlybiotic.rst")
          #raster_patches<-raster("raster_patches_RCP26_change_onlybiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlybiotic_dis.rst")
          scenario_name<-"RCP26_change_onlybiotic"
        }
        #For RCP60
        if (sc==6){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60_change_onlybiotic.RData")
          load("RCP60_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlybiotic.rst")
          #raster_patches<-raster("raster_patches_RCP60_change_onlybiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlybiotic_dis.rst")
          scenario_name<-"RCP60_change_onlybiotic"
        }
        #For RCP85
        if (sc==7){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85_change_onlybiotic.RData")
          load("RCP85_change_onlybiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlybiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlybiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlybiotic.rst")
          #raster_patches<-raster("raster_patches_RCP85_change_onlybiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlybiotic_dis.rst")
          scenario_name<-"RCP85_change_onlybiotic"
        }
        #Scenario with only change in biotic (We are goint to take the abiotic in the current scenario)
        #For RCP26
        if (sc==8){
          #We load the scenario 
          load("raster_habitat_integrated_RCP26_change_onlyabiotic.RData")
          load("RCP26_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP26_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP26_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP26_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP26_change_onlyabiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP26_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP26_change_onlyabiotic"
        }
        #For RCP60
        if (sc==9){
          #We load the scenario 
          load("raster_habitat_integrated_RCP60_change_onlyabiotic.RData")
          load("RCP60_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP60_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP60_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP60_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP60_change_onlyabiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP60_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP60_change_onlyabiotic"
        }
        #For RCP85
        if (sc==10){
          #We load the scenario 
          load("raster_habitat_integrated_RCP85_change_onlyabiotic.RData")
          load("RCP85_change_onlyabiotic.RData")
          scaled_variables_for_prediction<-RCP85_change_onlyabiotic
          raster_habitat_integrated<-raster_habitat_integrated_RCP85_change_onlyabiotic
          raster_habitat_integrated_binary<-raster("raster_habitat_integrated_binary_RCP85_change_onlyabiotic.rst")
          #raster_patches<-raster("raster_patches_RCP85_change_onlyabiotic.rst")
          #raster_dis<-raster("raster_NO_habitat_binary_RCP85_change_onlyabiotic_dis.rst")
          scenario_name<-"RCP85_change_onlyabiotic"
        }
        if (sc>4){
          scaled_variables_for_prediction$inte_row<-floor((c(1:nrow(scaled_variables_for_prediction))/100000))
        }
        habitat_change<-raster_habitat_integrated_binary_current-raster_habitat_integrated_binary
        NAvalue(habitat_change) <- -32768
        habitat_change[habitat_change == -32768 ] <- -1000
        writeRaster(habitat_change, paste0("raster_habitat_change_binary_", scenario_name, ".rst"), datatype='INT4S', overwrite=TRUE)
      }

######################################################################################################################  
#9.7 Bayesian hierarchical model ONLY CLIMATE with current data using results of Bayesian model with historical data as priors
###################################################################################################################### 
  #We change the R version, to 4.02
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the model ata range scale with historical data
    load("summary_m10.4_rstanarm.RData")
    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())
  
  #9.7.1 We define the priors location (mean) scale (standard deviation)
    my_prior_integrated <- normal(location = c(summary_m10.4_rstanarm[2,1],summary_m10.4_rstanarm[3,1],summary_m10.4_rstanarm[4,1],summary_m10.4_rstanarm[5,1],summary_m10.4_rstanarm[6,1],summary_m10.4_rstanarm[7,1],summary_m10.4_rstanarm[8,1],summary_m10.4_rstanarm[9,1]), scale = c(summary_m10.4_rstanarm[2,3],summary_m10.4_rstanarm[3,3],summary_m10.4_rstanarm[4,3],summary_m10.4_rstanarm[5,3],summary_m10.4_rstanarm[6,3],summary_m10.4_rstanarm[7,3],summary_m10.4_rstanarm[8,3],summary_m10.4_rstanarm[9,3]))#, autoscale = FALSE
    my_prior_intercept_integrated<- normal(location = c(summary_m10.4_rstanarm[1,1]), scale = c(summary_m10.4_rstanarm[1,3]))#, autoscale = FALSE
    #my_prior_intercept_random_current<- normal(location = c(0), scale = c(10))#, autoscale = FALSE
      #See https://mc-stan.org/rstanarm/reference/priors.html#examples  
      #Covariance
      #~ decov(reg. = 1, conc. = 1, shape = 1, scale = 1)  
  
  #9.7.2 We write the model  
    #Time 37045.3 seconds (Total)  
    mod_stan_integrated_only_climate<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_10 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_10_c +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = my_prior_integrated, prior_intercept = my_prior_intercept_integrated)
    save(mod_stan_integrated_only_climate,file="mod_stan_integrated_only_climate.Rdata")

  #9.7.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_integrated_only_climate)
    #A summary of the model  
    summary_mod_stan_integrated_only_climate<- as.data.frame(summary(mod_stan_integrated_only_climate))
    save (summary_mod_stan_integrated_only_climate, file='summary_mod_stan_integrated_only_climate.RData')
    write.csv(summary_mod_stan_integrated_only_climate, file='summary_mod_stan_integrated_only_climate.csv') # Supplementary Table FINALLY NOT INCLUDED IN METHODS & SUPPLEMENTARY #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: FINALLY NOT INCLUDED IN METHODS & SUPPLEMENTARY#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("FINALLY NOT INCLUDED IN METHODS & SUPPLEMENTARY.Results for the simple Bayesian model (no hierarchical) using abiotic and biotic factors to explain brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_only_climate [1:24,]

        print(paste0("OUTPUT: FINALLY NOT INCLUDED IN METHODS & SUPPLEMENTARY mean_PPD Model SHMABC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table FINALLY NOT INCLUDED IN METHODS & SUPPLEMENTARY. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_only_climate [25,]

  #9.7.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      waic(mod_stan_integrated_only_climate)
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated_only_climate, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)
      #pareto_k_values(loo1)

  #9.7.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_integrated_only_climate)
       
  #9.7.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_integrated_only_climate)
    str(posterior_df)
    cor_poterior_mod_stan_integrated_only_climate<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_integrated_only_climate, file='cor_poterior_mod_stan_integrated_only_climate.xlsx')

  #9.7.5 We write a bayesian null model GO TO 5.5 SAME MODEL
    #WE define the priors location (mean) scale (standard deviation)
    #mod_stan_integrated_only_climate_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop_integer), data = presences_absences_20km, family = "binomial",  prior = normal(), prior_intercept = normal())
    #save (mod_stan_integrated_only_climate_null, file='mod_stan_integrated_only_climate_null.RData')
    #waic(mod_stan_integrated_only_climate_null)
    #save(mod_stan_integrated_only_climate_null,file="mod_stan_integrated_only_climate_null.Rdata")

  #9.7.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_integrated_only_climate)
    dim(posterior)
    str(posterior)
    X11()
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    plot(mod_stan_integrated_only_climate)  

######################################################################################################################  
#9.8 Bayesian hierarchical model ONLY LAND USE with current data using results of Bayesian model with historical data as priors
###################################################################################################################### 
  #We change the R version, to 4.02
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the model ata range scale with historical data
    load("summary_m10.4_rstanarm.RData")
    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())

  #9.8.1 We define the priors location (mean) scale (standard deviation)

  #9.8.2 We write the model  
    #Time 37045.3 seconds (Total)  
    mod_stan_integrated_only_land_use<-stan_glmer(Europe_presences_0_1~ LC_1 + LC_3 + LC_4 + Frag +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = normal(), prior_intercept = normal())
    save(mod_stan_integrated_only_land_use,file="mod_stan_integrated_only_land_use.Rdata")

  #9.8.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_integrated_only_land_use)
    #A summary of the model  
    summary_mod_stan_integrated_only_land_use<- as.data.frame(summary(mod_stan_integrated_only_land_use))
    save (summary_mod_stan_integrated_only_land_use, file='summary_mod_stan_integrated_only_land_use.RData')
    write.csv(summary_mod_stan_integrated_only_land_use, file='summary_mod_stan_integrated_only_land_use.csv')

  #9.8.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
    waic(mod_stan_integrated_only_land_use)
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated_only_land_use, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)
      #pareto_k_values(loo1)

  #9.8.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_integrated_only_land_use)
    X11()
    pairs(mod_stan_integrated_only_land_use, pars = c("LC_1", "LC_3", "LC_4", "Frag"))

  #9.8.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_integrated_only_land_use)
    str(posterior_df)
    cor_poterior_mod_stan_integrated_only_land_use<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_integrated_only_land_use, file='cor_poterior_mod_stan_integrated_only_land_use.xlsx')

  #9.8.5 We write a bayesian null model GO TO 5.5 SAME MODEL
    #WE define the priors location (mean) scale (standard deviation)
    #mod_stan_integrated_only_land_use_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop_integer), data = presences_absences_20km, family = "binomial",  prior = normal(), prior_intercept = normal())
    #save (mod_stan_integrated_only_land_use_null, file='mod_stan_integrated_only_land_use_null.RData')
    #waic(mod_stan_integrated_only_land_use_null)
    #save(mod_stan_integrated_only_land_use_null,file="mod_stan_integrated_only_land_use_null.Rdata")

  #9.8.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_integrated_only_land_use)
    dim(posterior)
    str(posterior)
    X11()
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    plot(mod_stan_integrated_only_land_use)  

######################################################################################################################  
#9.9 Bayesian hierarchical model ONLY ABIOTIC with current data using results of Bayesian model with historical data as priors
###################################################################################################################### 
  #We change the R version, to 4.02
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the model ata range scale with historical data
    load("summary_m10.4_rstanarm.RData")
    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())

  #9.9.1 We define the priors location (mean) scale (standard deviation)
    my_prior_integrated <- normal(location = c(summary_m10.4_rstanarm[2,1],summary_m10.4_rstanarm[3,1],summary_m10.4_rstanarm[4,1],summary_m10.4_rstanarm[5,1],summary_m10.4_rstanarm[6,1],summary_m10.4_rstanarm[7,1],summary_m10.4_rstanarm[8,1],summary_m10.4_rstanarm[9,1],0,0,0,0), scale = c(summary_m10.4_rstanarm[2,3],summary_m10.4_rstanarm[3,3],summary_m10.4_rstanarm[4,3],summary_m10.4_rstanarm[5,3],summary_m10.4_rstanarm[6,3],summary_m10.4_rstanarm[7,3],summary_m10.4_rstanarm[8,3],summary_m10.4_rstanarm[9,3],2.5,2.5,2.5,2.5))#, autoscale = FALSE
    my_prior_intercept_integrated<- normal(location = c(summary_m10.4_rstanarm[1,1]), scale = c(summary_m10.4_rstanarm[1,3]))#, autoscale = FALSE
      #See https://mc-stan.org/rstanarm/reference/priors.html#examples  

  #9.9.2 We write the model  
    #Time 37045.3 seconds (Total)  
    #mod_stan_integrated_abiotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag +(1 | factor_subpop),
    #data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = my_prior_integrated, prior_intercept = my_prior_intercept_integrated)
    
    mod_stan_integrated_abiotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag + 
                                                        (1 | factor_subpop),
                                                        data = presences_absences_5km_scaled_only_current,
                                                        family = "binomial",
                                                        prior = my_prior_integrated,
                                                        prior_intercept = my_prior_intercept_integrated,
                                                        iter = 6000,
                                                        adapt_delta=0.999)
    save(mod_stan_integrated_abiotic,file="mod_stan_integrated_abiotic.Rdata")


    
  #9.9.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_integrated_abiotic)
    #A summary of the model  
    summary_mod_stan_integrated_abiotic<- as.data.frame(summary(mod_stan_integrated_abiotic))
    save (summary_mod_stan_integrated_abiotic, file='summary_mod_stan_integrated_abiotic.RData')
    write.csv(summary_mod_stan_integrated_abiotic, file='summary_mod_stan_integrated_abiotic.csv')# Supplementary Table 46 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 46#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 46. Results for the Bayesian hierarchical model using abiotic factors to explain brown bear distribution SHM_AI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_abiotic [1:13,]

        print(paste0("OUTPUT: Supplementary Table 43 mean_PPD Model SHM_AI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_abiotic [29,]

        

  #9.9.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      
      print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHMABC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_integrated_abiotic)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated_abiotic, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)
      #pareto_k_values(loo1)

  #9.9.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_integrated_abiotic)
       
  #9.9.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_integrated_abiotic)
    str(posterior_df)
    cor_poterior_mod_stan_integrated_abiotic<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_integrated_abiotic, file='cor_poterior_mod_stan_integrated_abiotic.xlsx')
    print(paste0("OUTPUT: Supplementary Table 51#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 51. Correlation of the posterior samples among the predictors used in the Bayesian hierarchical model using abiotic factors to explain brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(cor_poterior_mod_stan_integrated_abiotic)

  #9.9.5 We write a bayesian null model GO TO 5.5 SAME MODEL
 
  #9.9.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_integrated_abiotic)
    dim(posterior)
    str(posterior)
    print(paste0("OUTPUT: Supplementary Figure 10 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 10. Chains for the Bayesian model of brown bear habitat with abiotic factors combining data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 10 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    plot(mod_stan_integrated_abiotic)  

######################################################################################################################  
#9.10 Bayesian hierarchical model ONLY BIOTIC with current data using results of Bayesian model with historical data as priors
###################################################################################################################### 
  #We change the R version, to 4.02
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the model ata range scale with historical data
    load("summary_m10.4_rstanarm.RData")
    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())
  
  #9.10.1 We define the priors location (mean) scale (standard deviation)

  #9.10.2 We write the model  
    #Time 37045.3 seconds (Total)  
    #mod_stan_integrated_biotic<-stan_glmer(Europe_presences_0_1~ ene_Wild_rep_plant + ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +
                                                        #(1 | factor_subpop),
                                                        #data = presences_absences_5km_scaled_only_current,
                                                        #family = "binomial",  
                                                        #prior = normal(),
                                                        #prior_intercept = normal())
    mod_stan_integrated_biotic<-stan_glmer(Europe_presences_0_1~ ene_Wild_rep_plant + ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +
                                                        (1 | factor_subpop),
                                                        data = presences_absences_5km_scaled_only_current,
                                                        family = "binomial",
                                                        prior = normal(),
                                                        prior_intercept = normal(),
                                                        iter = 6000,
                                                        adapt_delta=0.999)
    save(mod_stan_integrated_biotic,file="mod_stan_integrated_biotic.Rdata")
load("mod_stan_integrated_biotic.Rdata")
mean(presences_absences_5km_scaled_only_current$Europe_presences_0_1)
  #9.10.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_integrated_biotic)
    #A summary of the model  
    summary_mod_stan_integrated_biotic<- as.data.frame(summary(mod_stan_integrated_biotic))
    save (summary_mod_stan_integrated_biotic, file='summary_mod_stan_integrated_biotic.RData')
    write.csv(summary_mod_stan_integrated_biotic, file='summary_mod_stan_integrated_biotic.csv') # Supplementary Table 48 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 48#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 48.Results for the simple Bayesian model (no hierarchical) using biotic factors to explain brown bear distribution. Model SHM_BC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_biotic [1:5,]

        print(paste0("OUTPUT: Supplementary Table 43 mean_PPD Model SHM_BC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_integrated_biotic [21,]

  #9.10.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHM_BC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution. Model SHM_BC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_integrated_biotic)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated_biotic, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)

  #9.10.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_integrated_biotic)
       
  #9.10.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_integrated_biotic)
    str(posterior_df)
    cor_poterior_mod_stan_integrated_biotic<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_integrated_biotic, file='cor_poterior_mod_stan_integrated_biotic.xlsx')
    print(paste0("OUTPUT: Supplementary Table 53#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 53. Correlation of the posterior samples among the predictors used in the simple Bayesian model (no hierarchical) using biotic factors to explain brown bear distribution. Model SHM_BC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(cor_poterior_mod_stan_integrated_biotic)

  #9.10.5 We write a bayesian null model GO TO 5.5 SAME MODEL
    #WE define the priors location (mean) scale (standard deviation)
    #mod_stan_integrated_biotic_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop_integer), data = presences_absences_20km, family = "binomial",  prior = normal(), prior_intercept = normal())
    #save (mod_stan_integrated_biotic_null, file='mod_stan_integrated_biotic_null.RData')
    #waic(mod_stan_integrated_biotic_null)
    #save(mod_stan_integrated_biotic_null,file="mod_stan_integrated_biotic_null.Rdata")

  #9.10.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_integrated_biotic)
    dim(posterior)
    str(posterior)
    print(paste0("OUTPUT: Supplementary Figure 12 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 12.Chains for the simple Bayesian model (no hierarchical) using biotic factors to explain brown bear distribution. Model SHM_BC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 12 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    plot(mod_stan_integrated_biotic)  

######################################################################################################################  
#9.11 Bayesian NOT hierarchical model ONLY ABIOTIC with current data 
###################################################################################################################### 
  #We change the R version, to 4.02
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())

  #9.11.1 We define the priors location (mean) scale (standard deviation)

  #9.11.2 We write the model  
    #Time 37045.3 seconds (Total)  
    # mod_stan_not_integrated_abiotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag +
    #                                               (1 | factor_subpop),
    #                                             data = presences_absences_5km_scaled_only_current,
    #                                             family = "binomial",
    #                                             prior = normal(),
    #                                             prior_intercept = normal())
    mod_stan_not_integrated_abiotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag +
                                                (1 | factor_subpop),
                                                data = presences_absences_5km_scaled_only_current,
                                                family = "binomial",
                                                prior = normal(),
                                                prior_intercept = normal(),
                                                iter = 6000,
                                                adapt_delta=0.999)  
    save(mod_stan_not_integrated_abiotic,file="mod_stan_not_integrated_abiotic.Rdata")
    load("mod_stan_not_integrated_abiotic.Rdata")
  #9.11.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_not_integrated_abiotic)
    #A summary of the model  
    summary_mod_stan_not_integrated_abiotic<- as.data.frame(summary(mod_stan_not_integrated_abiotic))
    save (summary_mod_stan_not_integrated_abiotic, file='summary_mod_stan_not_integrated_abiotic.RData')
    write.csv(summary_mod_stan_not_integrated_abiotic, file='summary_mod_stan_not_integrated_abiotic.csv')
        print(paste0("OUTPUT: Supplementary Table 47#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 47.Results for the simple Bayesian model using abiotic factors to explain brown bear distribution Model SHM_AC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_not_integrated_abiotic [1:13,]

        print(paste0("OUTPUT: Supplementary Table 43 mean_PPD Model SHM_AC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution Model SHM_AC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        summary_mod_stan_not_integrated_abiotic [29,]

  #9.11.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      print(paste0("OUTPUT: Supplementary Table 43 WAIC, elpd_WAIC and p_WAIC for Model SHM_AC #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 43. Performance for bayesian models (BMs) explaining brown bear distribution Model SHM_AC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_not_integrated_abiotic)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_not_integrated_abiotic, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)

  #9.11.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_not_integrated_abiotic)
       
  #9.11.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_not_integrated_abiotic)
    str(posterior_df)
    cor_poterior_mod_stan_not_integrated_abiotic<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_not_integrated_abiotic, file='cor_poterior_mod_stan_not_integrated_abiotic.xlsx')
    print(paste0("OUTPUT: Supplementary Table 52#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 52. Correlation of the posterior samples among the predictors used in the simple Bayesian model using abiotic factors to explain brown bear distribution Model SHM_AC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(cor_poterior_mod_stan_not_integrated_abiotic)

  #9.11.5 We write a bayesian null model GO TO 5.5 SAME MODEL
    #WE define the priors location (mean) scale (standard deviation)
    #mod_stan_not_integrated_abiotic_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop_integer), data = presences_absences_20km, family = "binomial",  prior = normal(), prior_intercept = normal())
    #save (mod_stan_not_integrated_abiotic_null, file='mod_stan_not_integrated_abiotic_null.RData')
    #waic(mod_stan_not_integrated_abiotic_null)
    #save(mod_stan_not_integrated_abiotic_null,file="mod_stan_not_integrated_abiotic_null.Rdata")

  #9.11.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_not_integrated_abiotic)
    dim(posterior)
    str(posterior)
    print(paste0("OUTPUT: Supplementary Figure 11 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 11. Chains for the simple Bayesian model using abiotic factors to explain brown bear distribution Model SHM_AC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    X11()# # Supplementary Figure 11 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    
    plot(mod_stan_not_integrated_abiotic)  

######################################################################################################################  
#9.12 Bayesian model ONLY CLIMATE with ONLY current data
###################################################################################################################### 
  #We set the working directory
    rm(list=ls()) 
    my_dir <-"writehereyourpath" 
    folder_working<-paste0(my_dir,"/9_Bayesian_SDM_Brown_Bear")
    setwd(folder_working)

    #We load the model ata range scale with historical data
    load("summary_m10.4_rstanarm.RData")
    #We load the current presences
    load("presences_absences_5km_scaled_only_current.RData") 
    #We say to R that detect the number of cores of the computer
    options(mc.cores = parallel::detectCores())
  
  #9.12.1 We define the priors location (mean) scale (standard deviation)

  #9.12.2 We write the model  
    mod_stan_current_only_climate<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_10 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_10_c +(1 | factor_subpop), data = presences_absences_5km_scaled_only_current, family = "binomial",  prior = normal(), prior_intercept = normal())
    #25845.1 seconds (Total) 
    save(mod_stan_current_only_climate,file="mod_stan_current_only_climate.Rdata")

  #9.12.3 Summaries of priors and estimates
    #A summary of the priors used in the model  
    prior_summary(mod_stan_current_only_climate)
    #A summary of the model  
    summary_mod_stan_current_only_climate<- as.data.frame(summary(mod_stan_current_only_climate))
    save (summary_mod_stan_current_only_climate, file='summary_mod_stan_current_only_climate.RData')
    write.csv(summary_mod_stan_current_only_climate, file='summary_mod_stan_current_only_climate.csv')

  #9.12.4 Model evaluation
    #Widely Applicable Information Criterion (WAIC). 
    #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      waic(mod_stan_current_only_climate)
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_current_only_climate, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)
      #pareto_k_values(loo1)

  #9.12.4.1 mcmc_pairs
    X11()
    mcmc_pairs(mod_stan_current_only_climate)
       
  #9.12.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
    posterior_df <- as.data.frame(mod_stan_current_only_climate)
    str(posterior_df)
    cor_poterior_mod_stan_current_only_climate<-cor(posterior_df)
    write.xlsx(cor_poterior_mod_stan_current_only_climate, file='cor_poterior_mod_stan_current_only_climate.xlsx')

  #9.12.5 We write a bayesian null model GO TO 5.5 SAME MODEL
    #WE define the priors location (mean) scale (standard deviation)
    #mod_stan_current_only_climate_null<-stan_glmer(Europe_presences_0_1~ (1 | factor_subpop_integer), data = presences_absences_20km, family = "binomial",  prior = normal(), prior_intercept = normal())
    #save (mod_stan_current_only_climate_null, file='mod_stan_current_only_climate_null.RData')
    #waic(mod_stan_current_only_climate_null)
    #save(mod_stan_current_only_climate_null,file="mod_stan_current_only_climate_null.Rdata")

  #9.12.6 Visualization/Checking the chain #8.3.6 CHECKING THE CHAIN
    #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
    #order, joined by a line.
    ## R code 8.12  
    posterior <- as.array(mod_stan_current_only_climate)
    dim(posterior)
    str(posterior)
    X11()
    color_scheme_set("mix-blue-red")
    mcmc_trace(posterior, #pars = c("wt", "sigma"), 
           facet_args = list(ncol = 6, strip.position = "left"))
    plot(mod_stan_current_only_climate)  

    
    
#9.13 Checking the error propagation in the model from 9.6 Bayesian hierarchical model with current data using results of Bayesian model with historical data as priors

#Schema:
  #9.13.1 Point estimate of mean
  #9.13.2 Replicas of models using the modified data with upper and lower bounds of biotic data
    #9.13.2.1 Replica of model using lower bound of biotic variables
      #9.13.2.1.1 We define the priors location (mean) scale (standard deviation)
      #9.13.2.1.2 We write the model  
      #9.13.2.1.3 Summaries of priors and estimates    
      #9.13.2.1.4 Model evaluation    
        #9.13.2.1.4.1 mcmc_pairs
        #9.13.2.1.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters  
      #9.13.2.1.5 We write a bayesian null model 
      #9.13.2.1.6 Visualization/Checking the chain
    #9.13.2.2 Replica of model using upper bound of biotic variables    
      #9.13.2.2.1 We define the priors location (mean) scale (standard deviation)
      #9.13.2.2.2 We write the model  
      #9.13.2.2.3 Summaries of priors and estimates    
      #9.13.2.2.4 Model evaluation    
        #9.13.2.2.4.1 mcmc_pairs
        #9.13.2.2.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters 
      #9.13.2.2.5 We write a bayesian null model 
      #9.13.2.2.6 Visualization/Checking the chain

  #9.13.3 Compassion of replicas of models with original model
    #9.13.3.1 Correlation among models
    


  #9.13.1 Point estimate of mean
  #We load the presences and the environmental variables:
  load("presences_absences_5km_scaled_only_current.RData") 
  data.frame(names(presences_absences_5km_scaled_only_current))
  summary(presences_absences_5km_scaled_only_current)    
    
  df_all_high<-c()   
  df_all_low<-c()   
  for (i in 1:4){
    column<-i+15  
    data<-presences_absences_5km_scaled_only_current[,column]  
    total <- length(data)  
    favourable <- mean(data, na.rm = TRUE)
    s <- sd(data)
    # calculate margin of error
    margin <- qt(0.975,df=total-1)*s/sqrt(total)
    # calculate lower and upper bounds of confidence interval
    low <- favourable - margin
    print(low)
    high <- favourable + margin
    print(high)
    #Apply the lower and upper bounds over the data
    data_high<- data+high    
    data_low<- data+low    
    #We bind the data creating a dataframe with high values and a dataframe with low values
    df_all_high<-cbind(df_all_high,data_high)
    df_all_low<-cbind(df_all_low,data_low)
  }  
  df_all_high2<-data.frame(df_all_high) 
  colnames(df_all_high2)<-names(presences_absences_5km_scaled_only_current)[16:19]
  df_all_low2<-data.frame(df_all_low) 
  colnames(df_all_low2)<-names(presences_absences_5km_scaled_only_current)[16:19]
  
  #Data frame with modified variables lower bound:
  presences_absences_5km_scaled_only_current_low_biotic<-data.frame(presences_absences_5km_scaled_only_current[,-c(16:19)],df_all_low2)
  save (presences_absences_5km_scaled_only_current_low_biotic, file='presences_absences_5km_scaled_only_current_low_biotic.RData')
  
  #Data frame with modified variables high bound:
  presences_absences_5km_scaled_only_current_high_biotic<-data.frame(presences_absences_5km_scaled_only_current[,-c(16:19)],df_all_high2)
  save (presences_absences_5km_scaled_only_current_high_biotic, file='presences_absences_5km_scaled_only_current_high_biotic.RData')



  #9.13.2 Replicas of models using the modified data with upper and lower bounds of biotic data
    # We are going to run two models with the same structure as 9.6, our hierarchical model using abiotic and biotic variables:

    #9.13.2.1 Replica of model using lower bound of biotic variables
      #We load the previous fitted model at range scale using historical data (we need to calculate the informative priors)
      load("summary_m10.4_rstanarm.RData")
      #We say to R that detect the number of cores of the computers
      options(mc.cores = parallel::detectCores())

      #9.13.2.1.1 We define the priors location (mean) scale (standard deviation)
      #We define the priors location (mean) scale (standard deviation)
      my_prior_integrated <- normal(location = c(summary_m10.4_rstanarm[2,1],summary_m10.4_rstanarm[3,1],summary_m10.4_rstanarm[4,1],summary_m10.4_rstanarm[5,1],summary_m10.4_rstanarm[6,1],summary_m10.4_rstanarm[7,1],summary_m10.4_rstanarm[8,1],summary_m10.4_rstanarm[9,1],0,0,0,0,0,0,0,0), scale = c(summary_m10.4_rstanarm[2,3],summary_m10.4_rstanarm[3,3],summary_m10.4_rstanarm[4,3],summary_m10.4_rstanarm[5,3],summary_m10.4_rstanarm[6,3],summary_m10.4_rstanarm[7,3],summary_m10.4_rstanarm[8,3],summary_m10.4_rstanarm[9,3],2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5))#, autoscale = FALSE
      my_prior_intercept_integrated<- normal(location = c(summary_m10.4_rstanarm[1,1]), scale = c(summary_m10.4_rstanarm[1,3]))#, autoscale = FALSE
      
      #9.13.2.1.2 We write the model  
      mod_stan_integrated2_low_bound_biotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag + ene_Wild_rep_plant +
                                                        ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +(1 | factor_subpop),
                                                        data = presences_absences_5km_scaled_only_current_low_biotic,
                                                        family = "binomial",
                                                        prior = my_prior_integrated,
                                                        prior_intercept = my_prior_intercept_integrated,
                                                        iter = 6000,
                                                        adapt_delta=0.999)
      save(mod_stan_integrated2_low_bound_biotic,file="mod_stan_integrated2_low_bound_biotic.Rdata")
      #  load(mod_stan_integrated2_low_bound_biotic.Rdata)
      
      #9.13.2.1.3 Summaries of priors and estimates    
      summary_mod_stan_integrated2_low_bound_biotic<- as.data.frame(summary(mod_stan_integrated2_low_bound_biotic))
      save (summary_mod_stan_integrated2_low_bound_biotic, file='summary_mod_stan_integrated2_low_bound_biotic.RData')
      write.csv(summary_mod_stan_integrated2_low_bound_biotic, file='summary_mod_stan_integrated2_low_bound_biotic.csv')# Supplementary Table 44 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          print(paste0("OUTPUT: Supplementary Table 70#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          print(paste0("Supplementary Table 70. Results for the Hiarerchical Bayesian model using abiotic and biotic factors to explain brown bear distribution Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          summary_mod_stan_integrated2_low_bound_biotic [1:17,]
  
          data.frame(rownames(summary_mod_stan_integrated2_low_bound_biotic))
          
          print(paste0("OUTPUT: Supplementary Table 69 mean_PPD Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          print(paste0("Supplementary Table 69. Performance for bayesian models (BMs) explaining brown bear distribution  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          summary_mod_stan_integrated2_low_bound_biotic [33,]

      #9.13.2.1.4 Model evaluation    
      #Widely Applicable Information Criterion (WAIC). 
      #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      print(paste0("OUTPUT: Supplementary Table 69 WAIC, elpd_WAIC and p_WAIC for Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 69. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_integrated2_low_bound_biotic)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated2_low_bound_biotic, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)

        #9.13.2.1.4.1 mcmc_pairs
        X11()
        mcmc_pairs(mod_stan_integrated2_low_bound_biotic)
       
        #9.6.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
        posterior_df <- as.data.frame(mod_stan_integrated2_low_bound_biotic)
        str(posterior_df)
        cor_poterior_mod_stan_integrated2_low_bound_biotic<-cor(posterior_df)
        write.xlsx(cor_poterior_mod_stan_integrated2_low_bound_biotic, file='cor_poterior_mod_stan_integrated2_low_bound_biotic.xlsx')
        print(paste0("OUTPUT: Supplementary Table 72#  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 72. Correlation of the posterior samples among the predictors used in the Bayesian hierarchical model using abiotic and biotic factors to explain brown bear distribution  Model SHM_ABI~~~~~"))
        print(cor_poterior_mod_stan_integrated2_low_bound_biotic)
    
      #9.13.2.1.5 We write a bayesian null model 
      print("This was calculated on section 9.5.5")
      
      #9.13.2.1.6 Visualization/Checking the chain
        #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
        #order, joined by a line.
        ## R code 8.12  
        posterior <- as.array(mod_stan_integrated2_low_bound_biotic)
        print(paste0("OUTPUT: Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 8.Chains for the Bayesian model of brown bear habitat with abiotic and biotic factors combining data.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()# # Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        color_scheme_set("mix-blue-red")
        mcmc_trace(posterior, #pars = c("wt", "sigma"), 
               facet_args = list(ncol = 6, strip.position = "left"))
     
       
      #9.13.2.1.7 We calculate the prediction for scenario current
        #We load the data of coordinates and ID pixels
        str(data_matrix)
        load(paste0(my_dir,"/0_Construction_of_the_Spatial_Database/data_matrix.RData"))
        #We load the fitted integrated model
        dfmod<-(mod_stan_integrated2_low_bound_biotic$glmod$fr)
        variables_model<-colnames(dfmod)[2:17]
  
        #FOR EACH SCENARIO NORMAL AND SIMULATED
        #for (sc in 1:10){  
        #sc=1
          print(paste0("#######################Scenario in loop ",sc,"#######################################"))
          #For current
          #if (sc==1){
          #For current
            #We load the scenario 
            load("scaled_variables_for_prediction_current.RData")
            scaled_variables_for_prediction<-scaled_variables_for_prediction_current
            scenario_name<-"current"
          #}
  
          max_integer<-max(scaled_variables_for_prediction$inte_row)
          df_pred_integr_all<-c()
          for (i in 0:max_integer){
            print(paste0("Loop i = ",i," #####"))
            sub_in_loop<-scaled_variables_for_prediction[scaled_variables_for_prediction$inte_row==i,]  
            #We extract the posterior
            link.mod_stan_integrated2_low_bound_biotic <- posterior_epred(mod_stan_integrated2_low_bound_biotic, newdata = sub_in_loop)
            # summarize, we calculate the mean of all the columns of the matrix
            sum_mean_integrated<-apply( link.mod_stan_integrated2_low_bound_biotic , 2 , mean)
            #We compute the percentile intervals, we select prob=0.95
            sum_PI_integrated <- apply( link.mod_stan_integrated2_low_bound_biotic , 2 , HPDI , prob=0.95 )
      
            df_sum_mean_integrated<-as.data.frame(sum_mean_integrated)
            df_sum_PI_integrated<-as.data.frame(sum_PI_integrated)
            df_sum_PI_t_integrated<-t(df_sum_PI_integrated)  
            df_sum_PI_t_integrated<-as.data.frame(df_sum_PI_t_integrated)
      
            df_pred_integr<-cbind(df_sum_mean_integrated,df_sum_PI_t_integrated)
            df_pred_integr$Uncertainity<-df_pred_integr[,c(3)]-df_pred_integr[,c(2)]
            df_pred_integr<-cbind(df_pred_integr,sub_in_loop$ID_pixel)
            colnames(df_pred_integr)<-c("sum_mean","|0.95","0.95|","Uncertainity","ID")
            df_pred_integr_all<-rbind(df_pred_integr_all,df_pred_integr)
          }
          #This is the prediction incuding different factor for each subpopulation
          merge_prediciton_integrated<-merge(data_matrix,df_pred_integr_all, by.x="ID_pixel", by.y = "ID",all.x = T)
          assign(paste0("merge_prediciton_integrated_biotic_lower_bound_",scenario_name),merge_prediciton_integrated)
          save(list=paste0("merge_prediciton_integrated_biotic_lower_bound_", scenario_name), file=paste0("merge_prediciton_integrated_biotic_lower_bound_", scenario_name, ".Rdata"))   
  
          #9.6.8.5 Creation of the rasters with the prediction of habitat and uncertainity
          #We create rasters in img format with the mean values of the predicion and the range of the confidence intervals  
          #We change the R version, close and open R in earlier version (3.5.3) because a strange error:
            
          #We save into a raster file to see in Arcgis
          new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
          #We are going to use Europe Albers Equal Area Conic  
          newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          #We defined the projection of our raster:
          projection(new_raster) <- newproj
      
          #For the habitat
            raster_habitat_integrated<-new_raster
            values(raster_habitat_integrated)<-round(merge_prediciton_integrated$sum_mean,digits=4)
            assign(paste0("raster_habitat_integrated_biotic_lower_bound_",scenario_name),raster_habitat_integrated)
            save(list=paste0("raster_habitat_integrated_biotic_lower_bound_", scenario_name), file=paste0("raster_habitat_integrated_biotic_lower_bound_", scenario_name, ".Rdata"))   
            print(paste0("OUTPUT: Figure 4 & Electronic Material 1-9, Raster files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            writeRaster(raster_habitat_integrated, paste0("raster_habitat_integrated_biotic_lower_bound_", scenario_name, ".img"), overwrite=TRUE) ##################################################################################### Supplementary Electronic Material 1-9
  
          #For the uncertainity
            raster_uncertainity_integrated<-new_raster
            values(raster_uncertainity_integrated)<-round(merge_prediciton_integrated$Uncertainity,digits=5)
            assign(paste0("raster_uncertainity_integrated_biotic_lower_bound_",scenario_name),raster_uncertainity_integrated)
            save(list=paste0("raster_uncertainity_integrated_biotic_lower_bound_", scenario_name), file=paste0("raster_uncertainity_integrated_biotic_lower_bound_", scenario_name, ".Rdata"))   
            print(paste0("OUTPUT: Electronic Material 10-18, Raster files with the uncertainity of the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            writeRaster(raster_uncertainity_integrated, paste0("raster_uncertainity_integrated_biotic_lower_bound_", scenario_name, ".img"), overwrite=TRUE) ############################################################################ Supplementary Electronic Material 10-18
  
          #We save it in a kml file to see in google earth      
            #We crate a new rster i lat long: 
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            raster_habitat_integrated_lat_long <- projectRaster(raster_habitat_integrated, new_raster_2, method='bilinear')
            print(paste0("OUTPUT: Electronic Material 19-27, Klm files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            KML(raster_habitat_integrated_lat_long, file=paste0("raster_habitat_integrated_kml_biotic_lower_bound_", scenario_name, ".kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)######################### Supplementary Electronic Material 19-27
        #}
      ########################################################################################
      ########################################################################################
      ########################################################################################
            
      
  #9.13.2.2 Replica of model using upper bound of biotic variables
      #We load the previous fitted model at range scale using historical data (we need to calculate the informative priors)
      load("summary_m10.4_rstanarm.RData")
      #We say to R that detect the number of cores of the computers
      options(mc.cores = parallel::detectCores())

      #9.13.2.2.1 We define the priors location (mean) scale (standard deviation)
      #We define the priors location (mean) scale (standard deviation)
      my_prior_integrated <- normal(location = c(summary_m10.4_rstanarm[2,1],summary_m10.4_rstanarm[3,1],summary_m10.4_rstanarm[4,1],summary_m10.4_rstanarm[5,1],summary_m10.4_rstanarm[6,1],summary_m10.4_rstanarm[7,1],summary_m10.4_rstanarm[8,1],summary_m10.4_rstanarm[9,1],0,0,0,0,0,0,0,0), scale = c(summary_m10.4_rstanarm[2,3],summary_m10.4_rstanarm[3,3],summary_m10.4_rstanarm[4,3],summary_m10.4_rstanarm[5,3],summary_m10.4_rstanarm[6,3],summary_m10.4_rstanarm[7,3],summary_m10.4_rstanarm[8,3],summary_m10.4_rstanarm[9,3],2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5))#, autoscale = FALSE
      my_prior_intercept_integrated<- normal(location = c(summary_m10.4_rstanarm[1,1]), scale = c(summary_m10.4_rstanarm[1,3]))#, autoscale = FALSE
      
      #9.13.2.2 Replica of model using upper bound of biotic variables    
      mod_stan_integrated2_high_bound_biotic<-stan_glmer(Europe_presences_0_1~ Clim_3+ Clim_4 + Clim_8+ Clim_9 + Clim_3_c + Clim_4_c + Clim_8_c + Clim_9_c + LC_1 + LC_3 + LC_4 + Frag + ene_Wild_rep_plant +
                                                        ene_Wild_invertebrates + ene_Wild_unk_plant_oth + ene_Wild_vertebrates +(1 | factor_subpop),
                                                        data = presences_absences_5km_scaled_only_current_high_biotic,
                                                        family = "binomial",
                                                        prior = my_prior_integrated,
                                                        prior_intercept = my_prior_intercept_integrated,
                                                        iter = 6000,
                                                        adapt_delta=0.999)
      save(mod_stan_integrated2_high_bound_biotic,file="mod_stan_integrated2_high_bound_biotic.Rdata")
      #  load(mod_stan_integrated2_high_bound_biotic.Rdata)
      
      #9.13.2.2.3 Summaries of priors and estimates    
      summary_mod_stan_integrated2_high_bound_biotic<- as.data.frame(summary(mod_stan_integrated2_high_bound_biotic))
      save (summary_mod_stan_integrated2_high_bound_biotic, file='summary_mod_stan_integrated2_high_bound_biotic.RData')
      write.csv(summary_mod_stan_integrated2_high_bound_biotic, file='summary_mod_stan_integrated2_high_bound_biotic.csv')# Supplementary Table 44 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          print(paste0("OUTPUT: Supplementary Table 71#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          print(paste0("Supplementary Table 71. Results for the Hiarerchical Bayesian model using abiotic and biotic factors to explain brown bear distribution Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          summary_mod_stan_integrated2_high_bound_biotic [1:17,]
  
          data.frame(rownames(summary_mod_stan_integrated2_high_bound_biotic))
          
          print(paste0("OUTPUT: Supplementary Table 69 mean_PPD Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          print(paste0("Supplementary Table 69. Performance for bayesian models (BMs) explaining brown bear distribution  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
          summary_mod_stan_integrated2_high_bound_biotic [33,]

      #9.13.2.2.4 Model evaluation    
      #Widely Applicable Information Criterion (WAIC). 
      #References: 
      #Watanabe, S. 2010. Asymptotic equivalence of Bayes cross validation and Widely Applicable 
      #Information Criterion in singular learning theory. Journal of Machine Learning Research 11:3571-3594.
      #Vehtari, A., A. Gelman, and J. Gabry. 2015. Efficient implementation of leave-one-out cross-validation 
      #and WAIC for evaluating fitted Bayesian models.
      print(paste0("OUTPUT: Supplementary Table 69 WAIC, elpd_WAIC and p_WAIC for Model SHM_ABI #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 69. Performance for bayesian models (BMs) explaining brown bear distribution~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      waic(mod_stan_integrated2_high_bound_biotic)# Supplementary Table 43 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #FRom https://cran.r-project.org/web/packages/loo/vignettes/loo2-example.html
      loo1 <- loo(mod_stan_integrated2_high_bound_biotic, save_psis = TRUE)#
      print(loo1)
      X11()
      plot(loo1)
      pareto_k_table(loo1)

        #9.13.2.2.4.1 mcmc_pairs
        X11()
        mcmc_pairs(mod_stan_integrated2_high_bound_biotic)
       
        #9.6.4.2 run cor() on the posterior samples to get an estimate of correlations among parameters.  
        posterior_df <- as.data.frame(mod_stan_integrated2_high_bound_biotic)
        str(posterior_df)
        cor_poterior_mod_stan_integrated2_high_bound_biotic<-cor(posterior_df)
        write.xlsx(cor_poterior_mod_stan_integrated2_high_bound_biotic, file='cor_poterior_mod_stan_integrated2_high_bound_biotic.xlsx')
        print(paste0("OUTPUT: Supplementary Table 73#  Model SHM_ABI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 73. Correlation of the posterior samples among the predictors used in the Bayesian hierarchical model using abiotic and biotic factors to explain brown bear distribution  Model SHM_ABI~~~~~"))
        print(cor_poterior_mod_stan_integrated2_high_bound_biotic)
    
      #9.13.2.2.5 We write a bayesian null model 
      print("This was calculated on section 9.5.5")
      
      #9.13.2.2.6 Visualization/Checking the chain
        #The most useful tool for diagnosing malfunction is a trace plot. A trace plot merely plots the samples in sequencial 
        #order, joined by a line.
        ## R code 8.12  
        posterior <- as.array(mod_stan_integrated2_high_bound_biotic)
        print(paste0("OUTPUT: Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Figure 8.Chains for the Bayesian model of brown bear habitat with abiotic and biotic factors combining data.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        X11()# # Supplementary Figure 8 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        color_scheme_set("mix-blue-red")
        mcmc_trace(posterior, #pars = c("wt", "sigma"), 
               facet_args = list(ncol = 6, strip.position = "left"))
           
      
  ################################################################################    
        #9.13.2.1.7 We calculate the prediction for scenario current
        #We load the data of coordinates and ID pixels
        str(data_matrix)
        load(paste0(my_dir,"/0_Construction_of_the_Spatial_Database/data_matrix.RData"))
        #We load the fitted integrated model
        dfmod<-(mod_stan_integrated2_high_bound_biotic$glmod$fr)
        variables_model<-colnames(dfmod)[2:17]
  
        #FOR EACH SCENARIO NORMAL AND SIMULATED
        #for (sc in 1:10){  
        #sc=1
          print(paste0("#######################Scenario in loop ",sc,"#######################################"))
          #For current
          #if (sc==1){
          #For current
            #We load the scenario 
            load("scaled_variables_for_prediction_current.RData")
            scaled_variables_for_prediction<-scaled_variables_for_prediction_current
            scenario_name<-"current"
          #}
  
          max_integer<-max(scaled_variables_for_prediction$inte_row)
          df_pred_integr_all<-c()
          for (i in 0:max_integer){
            print(paste0("Loop i = ",i," #####"))
            sub_in_loop<-scaled_variables_for_prediction[scaled_variables_for_prediction$inte_row==i,]  
            #We extract the posterior
            link.mod_stan_integrated2_high_bound_biotic <- posterior_epred(mod_stan_integrated2_high_bound_biotic, newdata = sub_in_loop)
            # summarize, we calculate the mean of all the columns of the matrix
            sum_mean_integrated<-apply( link.mod_stan_integrated2_high_bound_biotic , 2 , mean)
            #We compute the percentile intervals, we select prob=0.95
            sum_PI_integrated <- apply( link.mod_stan_integrated2_high_bound_biotic , 2 , HPDI , prob=0.95 )
      
            df_sum_mean_integrated<-as.data.frame(sum_mean_integrated)
            df_sum_PI_integrated<-as.data.frame(sum_PI_integrated)
            df_sum_PI_t_integrated<-t(df_sum_PI_integrated)  
            df_sum_PI_t_integrated<-as.data.frame(df_sum_PI_t_integrated)
      
            df_pred_integr<-cbind(df_sum_mean_integrated,df_sum_PI_t_integrated)
            df_pred_integr$Uncertainity<-df_pred_integr[,c(3)]-df_pred_integr[,c(2)]
            df_pred_integr<-cbind(df_pred_integr,sub_in_loop$ID_pixel)
            colnames(df_pred_integr)<-c("sum_mean","|0.95","0.95|","Uncertainity","ID")
            df_pred_integr_all<-rbind(df_pred_integr_all,df_pred_integr)
          }
          #This is the prediction incuding different factor for each subpopulation
          merge_prediciton_integrated<-merge(data_matrix,df_pred_integr_all, by.x="ID_pixel", by.y = "ID",all.x = T)
          assign(paste0("merge_prediciton_integrated_biotic_higher_bound_",scenario_name),merge_prediciton_integrated)
          save(list=paste0("merge_prediciton_integrated_biotic_higher_bound_", scenario_name), file=paste0("merge_prediciton_integrated_biotic_higher_bound_", scenario_name, ".Rdata"))   
  
          #9.6.8.5 Creation of the rasters with the prediction of habitat and uncertainity
          #We create rasters in img format with the mean values of the predicion and the range of the confidence intervals  
          #We change the R version, close and open R in earlier version (3.5.3) because a strange error:
            
          #We save into a raster file to see in Arcgis
          new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
          #We are going to use Europe Albers Equal Area Conic  
          newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          #We defined the projection of our raster:
          projection(new_raster) <- newproj
      
          #For the habitat
            raster_habitat_integrated<-new_raster
            values(raster_habitat_integrated)<-round(merge_prediciton_integrated$sum_mean,digits=4)
            assign(paste0("raster_habitat_integrated_biotic_higher_bound_",scenario_name),raster_habitat_integrated)
            save(list=paste0("raster_habitat_integrated_biotic_higher_bound_", scenario_name), file=paste0("raster_habitat_integrated_biotic_higher_bound_", scenario_name, ".Rdata"))   
            print(paste0("OUTPUT: Figure 4 & Electronic Material 1-9, Raster files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            writeRaster(raster_habitat_integrated, paste0("raster_habitat_integrated_biotic_higher_bound_", scenario_name, ".img"), overwrite=TRUE) ##################################################################################### Supplementary Electronic Material 1-9
  
          #For the uncertainity
            raster_uncertainity_integrated<-new_raster
            values(raster_uncertainity_integrated)<-round(merge_prediciton_integrated$Uncertainity,digits=5)
            assign(paste0("raster_uncertainity_integrated_biotic_higher_bound_",scenario_name),raster_uncertainity_integrated)
            save(list=paste0("raster_uncertainity_integrated_biotic_higher_bound_", scenario_name), file=paste0("raster_uncertainity_integrated_biotic_higher_bound_", scenario_name, ".Rdata"))   
            print(paste0("OUTPUT: Electronic Material 10-18, Raster files with the uncertainity of the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            writeRaster(raster_uncertainity_integrated, paste0("raster_uncertainity_integrated_biotic_higher_bound_", scenario_name, ".img"), overwrite=TRUE) ############################################################################ Supplementary Electronic Material 10-18
  
          #We save it in a kml file to see in google earth      
            #We crate a new rster i lat long: 
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            raster_habitat_integrated_lat_long <- projectRaster(raster_habitat_integrated, new_raster_2, method='bilinear')
            print(paste0("OUTPUT: Electronic Material 19-27, Klm files with the habitat prediction for current and future scenarios #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
            KML(raster_habitat_integrated_lat_long, file=paste0("raster_habitat_integrated_kml_biotic_higher_bound_", scenario_name, ".kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)######################### Supplementary Electronic Material 19-27
        #}
      ########################################################################################
      ########################################################################################
      ########################################################################################
           
       
      
    load(file=paste0("merge_prediciton_integrated_biotic_lower_bound_", scenario_name, ".Rdata"))  
    load(file=paste0("merge_prediciton_integrated_biotic_higher_bound_", scenario_name, ".Rdata"))  
    load(file=paste0("merge_prediciton_integrated_", scenario_name, ".Rdata"))  
      
      
     merge_prediciton_integrated_biotic_lower_bound_current 
     merge_prediciton_integrated_biotic_higher_bound_current 
     merge_prediciton_integrated_current 
      str(merge_prediciton_integrated_current)
      str(merge_prediciton_integrated_biotic_lower_bound_current)
      summary(merge_prediciton_integrated_biotic_lower_bound_current)
      summary(merge_prediciton_integrated_current)
     #Correlation with lower bound prediction
     cor(merge_prediciton_integrated_current$sum_mean, merge_prediciton_integrated_biotic_lower_bound_current$sum_mean, )
      
     #Correlation with lower bound prediction
      cor(merge_prediciton_integrated_current$sum_mean, merge_prediciton_integrated_biotic_lower_bound_current$sum_mean, use = "pairwise.complete.obs")
      # correlation = 0.9999991
      
     #Correlation with higher bound prediction
      cor(merge_prediciton_integrated_current$sum_mean, merge_prediciton_integrated_biotic_higher_bound_current$sum_mean, use = "pairwise.complete.obs")
      # correlation = 0.9999994
      
      
      
      
      
      
      
      
      
      
      
