

############################################################################################################
#Readme:
############################################################################################################
#R code for calculate the SDMs for the species in the brown bear diet
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
 #/0_Construction_of_the_Spatial_Database/Scenario_current.RData
 #/0_Construction_of_the_Spatial_Database/Scenario_RCP26.RData
 #/0_Construction_of_the_Spatial_Database/Scenario_RCP60.RData
 #/0_Construction_of_the_Spatial_Database/Scenario_RCP85.RData
 #/3_Download_occurrences_from_GBIF/df_datos_1_276_description_GBIF.RData
 #/0_Construction_of_the_Spatial_Database/data_matrix.RData
 # For each species:
  #/Buffer_10000/Buffer_10000_sps_",i,".rst
  #/Buffer_2000/Buffer_2000_sps_",i,".rst
  #/Rasters_in_vector/Vector_sps_", list_sps_in_loop$i, ".Rdata
#Data output:
  #/eval_DF_all3.xlsx
  #/3_Download_occurrences_from_GBIF/df_datos_1_276_description_GBIF.RData
  #Ensemble Forecasting of SDMs for Food species for Current and future Scenarios

##############################################################################################                
#Schema
##############################################################################################   
#4_SDM_Food_species
  #4.1 Prepare databases with environmental predictors
    #4.1.1 For scenario current
    #4.1.2 For scenario RCP26
    #4.1.3 For scenario RCP60
    #4.1.4 For scenario RCP85
  #4.2 Variable selection for habitat models
    #4.2.1 We load the data and stablish pseudoabsence parameters
    #4.2.2 Loop for each species in the diet    
      #4.2.2.1 Selection of variables based in univariable GLM and correlation among variables
        #4.2.1.1.1 Selection of important climate variables 
        #4.2.1.1.2 Selection of important land use variables 
  #4.3 We create a vector with the six selected best variables (three from clim and three from land use)     
  #4.4 Biomod analysis
        #4.4.1 We formating the data for biomod
        #4.4.2 Modeling options 
        #4.4.3 Building models 
        #4.4.4 Ensemble Modeling
        #4.4.5 Ensemble Modeling Evaluation and variable importance  
          #4.5.5.1 Evaluation
          #4.5.5.2 Variable importance
  #4.5 Ensemble Forecasting  
    #4.5.1 Ensemble Forecasting of Current Scenario 
      #4.5.1.1 We merge the data of the matrix to the Current conditions          
      #4.5.1.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      #4.5.1.3 Application of BIOMOD_EnsembleForecasting to Current Scenario    
      #4.5.1.4 Code for create a ID for the cells of raster
      #4.5.1.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      #4.5.1.6 We merge for each subpopulations the predictions for all species
      #4.5.1.7 We divide the dataframe with all species in smaller dataframes for do the merge
      #4.5.1.8 We do the merge for each small data frame
      #4.5.1.9 Create the raster files with the predicitons
      #4.5.1.10 Code for copy the results to folders
      	#4.5.1.10.1 Code for copy the results to folders for WILD species
      #4.5.1.11 Code for copy the results to unique tables with all species
      #4.5.1.11.1 Code for copy the results to unique tables with ONLY WILD species
    #4.5.2 Ensemble Forecasting of RCP_25 Scenario 
      #4.5.2.1 We merge the data of the matrix to the RCP_25 conditions          
      #4.5.2.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      #4.5.2.3 Application of BIOMOD_EnsembleForecasting to RCP_25 Scenario    
      #4.5.2.4 Code for create a ID for the cells of raster
      #4.5.2.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      #4.5.2.6 We merge for each subpopulations the predictions for all species
      #4.5.2.7 We divide the dataframe with all species in smaller dataframes for do the merge
      #4.5.2.8 We do the merge for each small data frame
      #4.5.2.9 Create the raster files with the predicitons
      #4.5.2.10 Code for copy the results to folders
      	#4.5.2.10.1 Code for copy the results to folders for WILD species
    #4.5.3 Ensemble Forecasting of RCP_60 Scenario 
      #4.5.3.1 We merge the data of the matrix to the RCP_60 conditions          
      #4.5.3.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      #4.5.3.3 Application of BIOMOD_EnsembleForecasting to RCP_60 Scenario    
      #4.5.3.4 Code for create a ID for the cells of raster
      #4.5.3.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      #4.5.3.6 We merge for each subpopulations the predictions for all species
      #4.5.3.7 We divide the dataframe with all species in smaller dataframes for do the merge
      #4.5.3.8 We do the merge for each small data frame
      #4.5.3.9 Create the raster files with the predicitons
      #4.5.3.10 Code for copy the results to folders
      	#4.5.3.10.1 Code for copy the results to folders for WILD species
      #4.5.3.11 Code for copy the results to unique tables with all species
      #4.5.3.11.1 Code for copy the results to unique tables with ONLY WILD species
    #4.5.4 Ensemble Forecasting of RCP_85 Scenario 
      #4.5.4.1 We merge the data of the matrix to the RCP_85 conditions          
      #4.5.4.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      #4.5.4.3 Application of BIOMOD_EnsembleForecasting to RCP_85 Scenario    
      #4.5.4.4 Code for create a ID for the cells of raster
      #4.5.4.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      #4.5.4.6 We merge for each subpopulations the predictions for all species
      #4.5.4.7 We divide the dataframe with all species in smaller dataframes for do the merge
      #4.5.4.8 We do the merge for each small data frame
      #4.5.4.9 Create the raster files with the predicitons
      #4.5.4.10 Code for copy the results to folders
      	#4.5.4.10.1 Code for copy the results to folders for WILD species
      #4.5.4.11 Code for copy the results to unique tables with all species
      #4.5.4.11.1 Code for copy the results to unique tables with ONLY WILD species
  #4.6 Calculation of change in habitat from current to future RCPs scenarios
    #4.6.1 Code for calculate and save the delta for each species
    #4.6.2 Code for load the information of the delta for each species

##############################################################################################                
#4.1 Prepare databases with environmental predictors
##############################################################################################                
  rm(list=ls()) 
  setwd("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2")
  load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/Scenario_current.RData")# Copy here to try to be faster
  load("H:/G/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP26.RData")
  load("H:/G/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP60.RData")
  load("H:/G/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP85.RData")
  str(Scenario_current)
  str(Scenario_RCP26)
  
  #4.1.1 For scenario current
    sub_Scenario_current<-Scenario_current[,c(1,7,12,15,21:27)]      
    str(sub_Scenario_current)
    save(sub_Scenario_current, file='sub_Scenario_current.RData')

  #4.1.2 For scenario RCP26
    sub_Scenario_RCP26<-Scenario_RCP26[,c(1,7,12,15,21:27)]      
    str(sub_Scenario_RCP26)
    save(sub_Scenario_RCP26, file='sub_Scenario_RCP26.RData')
  
  #4.1.3 For scenario RCP60
    sub_Scenario_RCP60<-Scenario_RCP60[,c(1,7,12,15,21:27)]      
    str(sub_Scenario_RCP60)
    save(sub_Scenario_RCP60, file='sub_Scenario_RCP60.RData')

  #4.1.4 For scenario RCP85
    sub_Scenario_RCP85<-Scenario_RCP85[,c(1,7,12,15,21:27)]      
    str(sub_Scenario_RCP85)
    save(sub_Scenario_RCP85, file='sub_Scenario_RCP85.RData')
    
      
##############################################################################################       
#4.2 Variable selection for habitat models
##############################################################################################      

  #4.2.1 We load the data and stablish pseudoabsence parameters
    rm(list=ls()) 
    library(biomod2) 
    library(slam)
    library(MASS)
    library(usdm)
    library(rgdal) 
    library(biomod2) 
    library(caret)
    library(MuMIn)
    library(raster)
    path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
    path_save = paste0(path_data,"SDM_MOD/")
    dir.create(path_save)
    setwd(path_save)
    dir.create(path_save)
    setwd(path_save)
    #Parameters of modelling:
    radio_min_pseudoabsences<-2000
    radio_max_pseudoabsences<-10000
    number_variables_selected_for_models<-6
    parameters_save = paste0(path_save,"Parameters/")
    dir.create(parameters_save)
       save(radio_min_pseudoabsences, file=paste0(parameters_save,"radio_min_pseudoabsences.RData"))
       save(radio_max_pseudoabsences, file=paste0(parameters_save,"radio_max_pseudoabsences.RData"))
       save(number_variables_selected_for_models, file=paste0(parameters_save,"number_variables_selected_for_models.RData"))
        #Data of list of presences of species in raster maps; 
        load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        #Data of coordinates and Ids of pixels of the raster
        load("D:/Project_Name/Results_Biomod/Data_for_analysis/data_matrix.RData")# Copy here to try to be faster
        #Scenario current (Current Climate and Land use variables)
        load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_current.RData")# Copy here to try to be faster
 
  #4.2.2 Loop for each species in the diet    
      print("############################################################################################################################")
      print("############################################################################################################################")
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print(paste0("             ######## STARTING LOOP FOR SELECTION OF VARIABLES ",Sys.time(),"########          "))
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print("############################################################################################################################")
      print("############################################################################################################################")
      DF_n_presencesabsesences_all<-c()
        for (i in 1:276){#The number of species to model Originally 276
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          #Condition of a minimum number of presences below n=150 we consider that we have not enought data
          list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
          if ((list_sps_in_loop$N_pixels_presence)>=50){
            path_save_sps = paste0(path_save,"sps_",list_sps_in_loop$i,"/")
            dir.create(path_save_sps)
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Starting Variable selection for species ",list_sps_in_loop$Species,"  ###   "))
            #We select the species that will be run in the loop. This is a column with 1 for presences and Na for the other areas of the map (including sea and pseudbasences)            
              file_load_sps=load(file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters_in_vector/Vector_sps_", list_sps_in_loop$i, ".Rdata",sep=""))       
              sp_1_europa = get(file_load_sps)
              buf_max<-raster(paste0("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Buffer_10000/Buffer_10000_sps_",i,".rst"))
              buf_min<-raster(paste0("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Buffer_2000/Buffer_2000_sps_",i,".rst"))

    #4.2.2.1 Selection of variables based in univariable GLM and correlation among variables
                    #We extract the values of the raster to a vector:
                        vec_buf_min<-values(buf_min)
                        vec_buf_max<-values(buf_max)
                        #We change the NA values of the buffer by 0
                        vec_buf_min[is.na(vec_buf_min)] <- 0
                        vec_buf_max[is.na(vec_buf_max)] <- 0
                #We merge data of Species presence, data about the coordinates, data of climate and data of land cover
                sp_1_europa[is.na(sp_1_europa)]<-0
                DF_LC_current_data_presence<-as.data.frame(cbind(sp_1_europa,data_matrix,sub_Scenario_current))#,vec_density_with_NAs_scaled))
                #We add the data of buffers minimum and maximum:
                data_pre<-cbind(DF_LC_current_data_presence,vec_buf_min,vec_buf_max)
                #We exclude areas outside the buffer
                data_pre2<-subset(data_pre, vec_buf_max==1 )
                #We exclude the values that are in the sea (NA):
                data <- na.exclude(data_pre2)
                save(data, file=paste0(path_save_sps,"data.RData"))
                #We select the areas of presences
                data_presences = subset(data, sp_1_europa==1)
                n_presences<-length(data_presences$sp_1_europa)
                #We select the areas of absences
                data_in_buf_absences_circle<- subset(data, vec_buf_max==1 & vec_buf_min==0)
                n_absesences<-length(data_in_buf_absences_circle$sp_1_europa)
                DF_n_presencesabsesences<-as.data.frame(cbind(n_presences,n_absesences))
                save(DF_n_presencesabsesences, file=paste0(path_save_sps,"DF_n_presencesabsesences.RData"))
                    #We calculate the pseudoabsences:
                        data_absences<-data_in_buf_absences_circle[sample(nrow(data_in_buf_absences_circle), n_presences), ]#20*n_presences
                        #As these are pseudoabsences we return the value of 0 to this data for sp_1_europa
                        data_absences$sp_1_europa<-0
                        #We merge the presences and pseudoabsences calculated:
                            presences_absences<-rbind(data_presences,data_absences)
                            path_save_Univariable_models = paste0(path_save_sps,"Univariable_models/")
                            dir.create(path_save_Univariable_models)
                            save(presences_absences, file=paste0(path_save_Univariable_models,"presences_absences.RData"))
                            
      #4.2.1.1.1 Selection of important climate variables 
        vif_result_clim<-vifstep(presences_absences[,5:8], th=10) #We have change from 31 to 32 to add variable frag
        save(vif_result_clim, file=paste0(path_save_Univariable_models,"vif_result_clim.RData"))
        names_select_vif_clim<-vif_result_clim @ results $Variables   
        save(names_select_vif_clim, file=paste0(path_save_Univariable_models,"names_select_vif_clim.RData"))
        data_for_mod_clim <- presences_absences[, names(presences_absences) %in% names_select_vif_clim, drop = F]
        data_for_mod_clim$sp_1_europa<- presences_absences$sp_1_europa
        n_var_average_clim<-length(names_select_vif_clim)
                      
        DF_AIC_univar_models_clim<-c()
        for (v in 1:n_var_average_clim){
          var_in_loop_clim<-names_select_vif_clim[v]
          write_model_prefix_clim<-c("(sp_1_europa)~")
          fixed_term_clim<-var_in_loop_clim
          write_model_var_in_loop_clim_1<-paste(var_in_loop_clim,sep="")# Efecto lineal
          write_model_var_in_loop_clim_2<-paste("I(",var_in_loop_clim,"^2",")",sep="")# Efecto cuadratico
          write_model_var_in_loop_clim_3<-paste(var_in_loop_clim," + ","I(",var_in_loop_clim,"^2",")",sep="")# Efecto lineal y cuadratico,
      
          write_model_1_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_1)# Efecto lineal
          write_model_2_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_2)# Efecto cuadratico
          write_model_3_clim<-paste(write_model_prefix_clim,write_model_var_in_loop_clim_3)# Efecto lineal y cuadratico
          
          from_mod_1_clim<-as.formula(write_model_1_clim)
          from_mod_2_clim<-as.formula(write_model_2_clim)
          from_mod_3_clim<-as.formula(write_model_3_clim)
          
          mod_1_clim<- glm(formula = from_mod_1_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          mod_2_clim<- glm(formula = from_mod_2_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          mod_3_clim<- glm(formula = from_mod_3_clim, data = data_for_mod_clim, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          
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
        #We select the 3 better
        variables_less_AICc_clim<-sorted_variables_clim[1:(number_variables_selected_for_models/2)]
        variables_names_less_AICc_clim<-names(variables_less_AICc_clim)
        save(variables_names_less_AICc_clim, file=paste0(path_save_Univariable_models,"variables_names_less_AICc_clim.RData"))
        save(DF_AIC_univar_models_clim, file=paste0(path_save_Univariable_models,"DF_AIC_univar_models_clim.RData"))
      
      #4.2.1.1.2 Selection of important land use variables 
        vif_result_use<-vifstep(presences_absences[,9:15], th=10) #We have change from 31 to 32 to add variable frag
        save(vif_result_use, file=paste0(path_save_Univariable_models,"vif_result_use.RData"))
        names_select_vif_use<-vif_result_use @ results $Variables   
        save(names_select_vif_use, file=paste0(path_save_Univariable_models,"names_select_vif_use.RData"))
        data_for_mod_use <- presences_absences[, names(presences_absences) %in% names_select_vif_use, drop = F]
        data_for_mod_use$sp_1_europa<- presences_absences$sp_1_europa
        n_var_average_use<-length(names_select_vif_use)
                      
        DF_AIC_univar_models_use<-c()
        for (r in 1:n_var_average_use){
          var_in_loop_use<-names_select_vif_use[r]
          write_model_prefix_use<-c("(sp_1_europa)~")
          fixed_term_use<-var_in_loop_use
          write_model_var_in_loop_use_1<-paste(var_in_loop_use,sep="")# Efecto lineal
          write_model_var_in_loop_use_2<-paste("I(",var_in_loop_use,"^2",")",sep="")# Efecto cuadratico
          write_model_var_in_loop_use_3<-paste(var_in_loop_use," + ","I(",var_in_loop_use,"^2",")",sep="")# Efecto lineal y cuadratico,
      
          write_model_1_use<-paste(write_model_prefix_use,write_model_var_in_loop_use_1)# Efecto lineal
          write_model_2_use<-paste(write_model_prefix_use,write_model_var_in_loop_use_2)# Efecto cuadratico
          write_model_3_use<-paste(write_model_prefix_use,write_model_var_in_loop_use_3)# Efecto lineal y cuadratico
          
          from_mod_1_use<-as.formula(write_model_1_use)
          from_mod_2_use<-as.formula(write_model_2_use)
          from_mod_3_use<-as.formula(write_model_3_use)
          
          mod_1_use<- glm(formula = from_mod_1_use, data = data_for_mod_use, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          mod_2_use<- glm(formula = from_mod_2_use, data = data_for_mod_use, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          mod_3_use<- glm(formula = from_mod_3_use, data = data_for_mod_use, family = "binomial", na.action = "na.fail")  #,weight=presences_absences[,"weight"] 
          
          #mod_null<- glm(formula = (sp_1_europa)~1, data = data_for_mod, family = "binomial", na.action = "na.fail") 
          AICc_mod_1_use<-AICc(mod_1_use)
          AICc_mod_2_use<-AICc(mod_2_use)
          AICc_mod_3_use<-AICc(mod_3_use)
          AICc_mod_use_123<-cbind(AICc_mod_1_use,AICc_mod_2_use,AICc_mod_3_use)
          AICc_mod_use_123_min<-min(AICc_mod_use_123[1,])
          DF_AIC_univar_models_use<-rbind(DF_AIC_univar_models_use,AICc_mod_use_123_min)
        }
        rownames(DF_AIC_univar_models_use)<-names_select_vif_use
        sorted_variables_use<-sort(DF_AIC_univar_models_use[,1])
        #We select the 3 better
        variables_less_AICc_use<-sorted_variables_use[1:(number_variables_selected_for_models/2)]
        variables_names_less_AICc_use<-names(variables_less_AICc_use)
        save(variables_names_less_AICc_use, file=paste0(path_save_Univariable_models,"variables_names_less_AICc_use.RData"))
        save(DF_AIC_univar_models_use, file=paste0(path_save_Univariable_models,"DF_AIC_univar_models_use.RData"))
       
  #4.3 We create a vector with the six selected best variables (three from clim and three from land use)     
     variables_names_less_AICc<-c(variables_names_less_AICc_clim,variables_names_less_AICc_use)
     save(variables_names_less_AICc, file=paste0(path_save_Univariable_models,"variables_names_less_AICc.RData"))
 
            print(paste0("N presences ",DF_n_presencesabsesences$n_presences,"## N absences ",DF_n_presencesabsesences$n_absesences,"########          "))
            print(paste0(" ######## Selected variables for model ########          "))
            print(paste0(variables_names_less_AICc))
            print(paste0("Relation N absences/N presences ",DF_n_presencesabsesences$n_absesences/DF_n_presencesabsesences$n_presences,"  ##  Should be bigger than 2 "))
            print(paste0("             ######## Variable selection finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   FINISHED    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          rm(file_load_sps)
          rm(variables_names_less_AICc_clim)
          rm(variables_names_less_AICc_use)
          rm(variables_names_less_AICc)
          rm(DF_AIC_univar_models_clim)
          rm(DF_AIC_univar_models_use)
          rm(names_select_vif_clim)
          rm(names_select_vif_use)
          rm(sp_1_europa)
          rm(data)
          rm(presences_absences)
          rm(DF_LC_current_data_presence)
          }#End of condition of a minimum number of presences    
        }#End of the loop for each species in the diet         


###################################################################################################################              
#4.4 Biomod analysis
###################################################################################################################              
  rm(list=ls()) 
  library(biomod2) 
  library(slam)
  library(MASS)
  library(usdm)
  library(rgdal) 
  library(biomod2) 
  library(caret)
  library(MuMIn)
  library(raster)
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
  path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
  path_save = paste0(path_data,"SDM_MOD/")
  parameters_save = paste0(path_save,"Parameters/")
  load(paste0(parameters_save,"radio_min_pseudoabsences.RData"))
  load(paste0(parameters_save,"radio_max_pseudoabsences.RData"))
  load(paste0(parameters_save,"number_variables_selected_for_models.RData"))
  setwd(path_save)
  print("############################################################################################################################")
  print("############################################################################################################################")
  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
  print(paste0("             ######## STARTING LOOP FOR BIOMOD ANALYSIS ",Sys.time(),"########          "))
  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
  print("############################################################################################################################")
  print("############################################################################################################################")
  for (i in 1:276){#The number of species to model Originally
    list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
    list_sps_in_loop$n_presences_sps<-(list_sps_in_loop$N_pixels_presence)
    list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
    if ((list_sps_in_loop$N_pixels_presence)>=50){
      path_save_sps = paste0(path_save,"sps_",list_sps_in_loop$i,"/")
      print("             ##############################################################          ")
      print("             ########           Starting a new species             ########          ")
      print("             ##############################################################          ")
      print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
      path_save_Univariable_models = paste0(path_save_sps,"Univariable_models/")
      load(file=paste0(path_save_Univariable_models,"variables_names_less_AICc.RData"))
      load(file=paste0(path_save_sps,"data.RData"))
      #We select a subset of columns reducing to 5 explanatory variables (5 Climatic)
      data_for_mod <- data[, names(data) %in% variables_names_less_AICc, drop = F]
      data_for_mod$sp_1_europa<-data$sp_1_europa
      data_for_mod$x<-data$x
      data_for_mod$y<-data$y
      data_for_mod$vec_buf_min<-data$vec_buf_min
      #We select the areas of presences
      data_for_mod_presences <- subset(data_for_mod, sp_1_europa==1)
      n_presences<-length(data_for_mod_presences$sp_1_europa)
      #We select the areas of absences
      data_in_buf_absences_circle <- subset(data_for_mod, vec_buf_min==0)
      n_absences_available<-length(data_in_buf_absences_circle$sp_1_europa)
        #Random selection:
        data_absences<-data_in_buf_absences_circle[sample(nrow(data_in_buf_absences_circle), n_presences*5), ]
      #We merge the presences and pseudoabsences calculated:
      presences_absences<-rbind(data_for_mod_presences,data_absences)
      #We convert the "0" of the column of presences of the species to NA
      sp_1<-presences_absences[,c("sp_1_europa")]
      sp_1[sp_1=="0"] <- NA
      sum_presences<-sum(sp_1,na.rm =T)                
      path_save_sps_Biomod = paste0(path_save_sps,"Biomod/")
      dir.create(path_save_sps_Biomod)
      setwd(path_save_sps_Biomod)
      
      #4.4.1 We formating the data for biomod
        format_data<-BIOMOD_FormatingData(resp.var=sp_1,expl.var= presences_absences[,c(1:number_variables_selected_for_models)],resp.xy=presences_absences[,c("x","y")],resp.name = "resp_name",#Change data_for_mod for presences_absences
          PA.nb.rep=2,PA.nb.absences=sum_presences, PA.strategy="random")         
        save(format_data, file=paste0(path_save_sps_Biomod,"format_data.RData"))
        save(presences_absences, file=paste0(path_save_sps_Biomod,"presences_absences.RData"))
      #4.4.2 Modeling options 
        #We customize the set of parameters and options of the models that are going to be build   
        ptm_Model_opt <- proc.time()  
        Model_opt<-BIOMOD_ModelingOptions(
          GLM=list(type="quadratic",interaction.level=0),
          GBM=list(n.trees=3000),
          RF = list(ntree = 750),
        )
        save(Model_opt, file=paste0(path_save_sps_Biomod,"Model_opt.RData"))

      #4.4.3 Building models 
        Sps_models<-BIOMOD_Modeling(
          data=format_data,
          models=c("GLM","GBM","RF"),
          models.options = Model_opt,
          NbRunEval = 2,
          DataSplit = 70,
          VarImport = 4,#VarImport = 0,
          models.eval.meth = c("TSS"),
          do.full.models=F,
          modeling.id = "ex2")
        save(Sps_models, file=paste0(path_save_sps_Biomod,"Sps_models.RData"))

      #4.4.4 Ensemble Modeling
        Sps_essemble_models<-BIOMOD_EnsembleModeling(modeling.output=Sps_models,
                                                     em.by = "all", #'PA_dataset+repet',
                                                     eval.metric = 'TSS',
                                                     eval.metric.quality.threshold = 0.2,
                                                     models.eval.meth = c('TSS'),
                                                     prob.mean = FALSE,
                                                     prob.cv = FALSE,
                                                     committee.averaging = TRUE,
                                                     prob.mean.weight = FALSE,
                                                     VarImport = 4)#VarImport = 0)
        save(Sps_essemble_models, file=paste0(path_save_sps_Biomod,"Sps_essemble_models.RData"))
      
      #4.5.5 Ensemble Modeling Evaluation and variable importance  
        #4.5.5.1 Evaluation
          eval<-get_evaluations(Sps_essemble_models,as.data.frame=T)
          save(eval, file=paste0(path_save_sps_Biomod,"eval.RData"))
          write.csv(eval, file=(paste0(path_save_sps_Biomod,'/eval.csv')))

        #4.5.5.2 Variable importance
          var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
          var_imp_df_mean<-apply(var_imp_df,1,mean)
          save (var_imp_df, file=(paste0(path_save_sps_Biomod,'/var_imp_df.RData')))
          write.csv(var_imp_df, file=(paste0(path_save_sps_Biomod,'/var_imp_df.csv')))
          save (var_imp_df_mean, file=(paste0(path_save_sps_Biomod,'/var_imp_df_mean.RData')))
          write.csv(var_imp_df_mean, file=(paste0(path_save_sps_Biomod,'/var_imp_df_mean.csv')))
      rm(variables_names_less_AICc)
      rm(data)
      print(paste0("             ######## Biomod analysis finished at ",Sys.time(),"########          "))
      print("############################################################################################################################")
      print("############################################################################################################################")
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   FINISHED    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
    }#End of condition of a minimum number of presences    
  }#End of the loop for each species in the diet   


##########################################################################################################
##########################################################################################################
##########################################################################################################
#4.5 Ensemble Forecasting  
##########################################################################################################
##########################################################################################################
##########################################################################################################

  ##########################################################################################################
  #4.5.1 Ensemble Forecasting of Current Scenario 
  ##########################################################################################################
  
    ##########################################################################################################
    #4.5.1.1 We merge the data of the matrix to the Current conditions 
    ##########################################################################################################
      rm(list=ls()) 
      path_data = "D:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("D:/Test_Biomod")
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_current.RData")# Copy here to try to be faster
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/data_matrix.RData")# Copy here to try to be faster
      DF_Scen_current_Matrix<-as.data.frame(cbind(data_matrix,sub_Scenario_current))
      save(DF_Scen_current_Matrix, file="D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix.RData")
      DF_Scen_current_Matrix_noNA<-na.omit(DF_Scen_current_Matrix)
      save(DF_Scen_current_Matrix_noNA, file="D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA.RData")

    ##########################################################################################################
    #4.5.1.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
    ##########################################################################################################
      rm(list=ls()) 
      path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("H:/D/Test_Biomod")
      library(rgdal) 
      library(raster) 
      mapa_extrapolacion <- shapefile('C:/Users/First_Author_Name/Documents/Project_Name/species_41688/area_extrap_200.shp')
      #Raster of reference: 
      new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
      #We are going to use Europe Albers Equal Area Conic  
      newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
      #We defined the projection of our raster:
      projection(new_raster) <- newproj
      #Rasterize the data
      extrap_raster<-rasterize(mapa_extrapolacion, new_raster, field=1)
      e2<- extent(-2800000,3100000,-800000,4800000) 
      vec_extrap_raster<-extract(extrap_raster,e2)
      save(vec_extrap_raster,file="vec_extrap_raster.RData")
      #Map of subpopulations
      mapa_subpopulations <- shapefile('H:/G/Project_Name/Maps/Extrapolation limits/Subpopulations3.shp')
      mapa_subpopulations_df<-as.data.frame(mapa_subpopulations)
      #Raster of reference: 
      new_raster2 <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
      #We defined the projection of our raster:
      projection(new_raster2) <- newproj
      #Rasterize the data
      subpop_raster<-rasterize(mapa_subpopulations, new_raster2, field="code_num")
      X11()
      plot(subpop_raster, col=rainbow(14), legend =T)
      vec_subpop_raster<-extract(subpop_raster,e2)
      table(vec_subpop_raster)
      save(vec_subpop_raster,file="vec_subpop_raster.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix.RData")# Copy here to try to be faster
      DF_Scen_current_Matrix_2<-as.data.frame(cbind(DF_Scen_current_Matrix,vec_extrap_raster,vec_subpop_raster))
      str(DF_Scen_current_Matrix_2)
      DF_Scen_current_Matrix_noNA_2<-na.omit(DF_Scen_current_Matrix_2)
      save(DF_Scen_current_Matrix_noNA_2, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")

    ##########################################################################################################
    #4.5.1.3 Application of BIOMOD_EnsembleForecasting to Current Scenario   
    ##########################################################################################################
      rm(list=ls()) 
      library(biomod2) 
      library(slam)
      library(MASS)
      library(usdm)
      library(rgdal) 
      library(caret)
      library(MuMIn)
      library(raster)
      path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("G:/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      merge_list_diet_energy<-merge(df_datos_1_276_description_GBIF,subpopulations_diet, by.x="Species", by.y = "species",all.x = T)
      merge_list_diet_energy2<-merge(subpopulations_diet,df_datos_1_276_description_GBIF, by.x="species", by.y = "Species",all.x = T)
      summary(merge_list_diet_energy2)
      summary(merge_list_diet_energy)
      sub_Na_merge_list_diet_energy<-merge_list_diet_energy[is.na(merge_list_diet_energy$central),]
      sub_Na_merge_list_diet_energy2<-merge_list_diet_energy2[is.na(merge_list_diet_energy2$N_pixels_presence),]
      #Descripcion de los NAs:
                #Para sub_Na_merge_list_diet_energy tenemos tres especies con NA que son Fagus orientalis, Lonicera caucasica y Pahseolus vulgaris que se les 
                #ha considerado rF_original = 0 en el estudio por lo que la REDEC sale NA, aunque la especie est? descrita coo parte de la dieta
                
                #Para sub_Na_merge_list_diet_energy2 tenemos un caso que es el Sus scrofa domesticus que est? en el estudio de Naves y no lo tenemos invluido como 
                #especie/subes pecie o taxon, con este sumarian 277 taxones
                
                #save(df_datos_1_276_description_GBIF, file="G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      table_range_environmental_data_current_scenario_all<-c()
      X11()
        for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=50 we consider that we have not enought data
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i,"/")
          path_save_Univariable_models = paste0(path_save_sps,"/Univariable_models/")
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod/")
          setwd(path_save_sps_Biomod)
          #This is the dataset with the predictors used in the model "Predictors_current" these data are unique for each species
          #We select from the dataframe with the extension to project the variables included in Predictors_current"
          load(file=paste0(path_save_sps_Biomod,"Sps_essemble_models.RData"))  
          variables_names<-Sps_essemble_models@ expl.var.names
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
          pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
          vec_populations_with_species<-colnames(pop_diet)
          vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
          vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
          vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
          vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
          vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
          vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
          vec_populations_with_species[vec_populations_with_species=="central"] <- 7
          vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
          vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
          vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
          vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
          vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
          vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
          vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
          vec_populations_with_species<-as.numeric(vec_populations_with_species)           
          data_pre_for_current_projection_0<-DF_Scen_current_Matrix_noNA_2[DF_Scen_current_Matrix_noNA_2$vec_subpop_raster %in% vec_populations_with_species,]
          data_pre_for_current_projection<-data_pre_for_current_projection_0[,variables_names]
          predicted_area_absolute<-NROW(data_pre_for_current_projection)
          var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
          var_imp_df_mean<-data.frame(apply(var_imp_df,1,mean))
          colnames(var_imp_df_mean)<-c("Var_importance")
          #Code for see the range of the data of variables used to fit the model and compare with the range of data of each scenario:
          #We load the data used to fit the model of the species
          load(file=paste0(path_save_sps_Biomod,"format_data.RData"))
          df_data_fit<-format_data @ data.env.var
          data_var_all<-c()
          pdf(file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Plots_environ_current/Range_env_data_sps_",i,".pdf",sep=""), width = 0.394 * 64, height = 0.394 * 8)
          par(mfrow = c(1, 6), cex = 0.7)
            for (col in 1:ncol(df_data_fit)) {
            d_all_europe <- density(data_pre_for_current_projection[,col])
            plot(d_all_europe, main=colnames(df_data_fit[col]))
            d <- density(df_data_fit[,col])
            lines(d, main=colnames(df_data_fit[col]))
            polygon(d, col="red", border="red")
            abline(v = max(df_data_fit[,col]), lty=2, col="blue")
            abline(v = min(df_data_fit[,col]), lty=2, col="blue")
            lines(d_all_europe, main=colnames(df_data_fit[col]))
            range_var_max<-max(df_data_fit[,col])
            range_var_min<-min(df_data_fit[,col])
            range_var<-max(df_data_fit[,col])-min(df_data_fit[,col])
            range_var_plus5per<-range_var_min + range_var*1.2
            range_var_less5per<-range_var_min-range_var*0.2
            data_in_range<-with(data_pre_for_current_projection, data_pre_for_current_projection[,col] <= range_var_plus5per & 
               data_pre_for_current_projection[,col] >= range_var_less5per)
            tab_var_europe<-table(data_in_range)
            df_tab_var_europe<-as.data.frame(tab_var_europe)
            per_n_cells_in_range<-df_tab_var_europe[which(df_tab_var_europe$data_in_range=="TRUE"),"Freq"]/6252207*100
            data_var<-data.frame(i,list_sps_in_loop$Species,colnames(df_data_fit[col]),range_var_min,range_var_max,range_var,range_var_less5per,range_var_plus5per,per_n_cells_in_range)
            data_var_all<-rbind(data_var_all,data_var)
            }
          colnames(data_var_all)<-c("ID_species","Species","Variable","range_var_min","range_var_max","range_var",
          "range_var_less5per","range_var_plus5per","per_n_cells_in_range")
          dev.off()  
          #Using the range of data of table_range_environmental_data we are going to convert the data ouside the range to NA
            for (j in 1:6) {
            data_pre_for_current_projection[,j][data_pre_for_current_projection[,j]>data_var_all[j,"range_var_plus5per"]]<-NA
            data_pre_for_current_projection[,j][data_pre_for_current_projection[,j]<data_var_all[j,"range_var_less5per"]]<-NA
            }
          #BIOMOD_EnsembleForecasting Current copy from WILFRED corrections
          Sps_esemble_curr_proj<-BIOMOD_EnsembleForecasting(EM.output=Sps_essemble_models,
          new.env=data_pre_for_current_projection,
          proj.name="Current",
          selected.models = 'all', 
          binary.meth = "TSS", 
          do.stack=FALSE,
          build.clamping.mask = F)
          save(Sps_esemble_curr_proj, file="Sps_esemble_curr_proj.RData")
          not_predicted<-sum(is.na(Sps_esemble_curr_proj@ proj @ val))
          predicted_area<-(predicted_area_absolute-not_predicted)/predicted_area_absolute*100
          save(predicted_area, file=paste0(path_save_sps_Biomod,"predicted_area_sps",i,".RData"))
          table_range_environmental_data_current_scenario<-data.frame(data_var_all,var_imp_df_mean,"Predict_area"=predicted_area)
          save(table_range_environmental_data_current_scenario, file=paste0(path_save_sps_Biomod,"table_range_environmental_data_current_scenario_sps_",i,".RData"))
          table_range_environmental_data_current_scenario_all<-rbind(table_range_environmental_data_current_scenario_all,table_range_environmental_data_current_scenario)   
          print(paste0("Precentage of predicted area for species ",list_sps_in_loop$Species," = ",predicted_area))
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        }
      setwd(path_data)
      save(table_range_environmental_data_current_scenario_all, file=paste0(path_save,"table_range_environmental_data_current_scenario_all.RData"))
      summary(table_range_environmental_data_current_scenario_all)
      
    ##########################################################################################################
    #4.5.1.4 Code for create a ID for the cells of raster
    ##########################################################################################################
      rm(list=ls()) 
      path_data = "D:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix.RData")# Copy here to try to be faster
      sub_DF_Scen_current_Matrix<-DF_Scen_current_Matrix[,c(2:3)]
      #We create a small dataframe with the ID of the cells of the raster
      save(sub_DF_Scen_current_Matrix, file="D:/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_current_Matrix.RData")

    #4.5.1.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      rm(list=ls())
      library(biomod2) 
      library(slam)
      library(MASS)
      library(usdm)
      library(rgdal) 
      library(caret)
      library(MuMIn)
      library(raster)
      library(rgdal) 
      library(raster) 
      path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("G:/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_current<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
          pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
          vec_populations_with_species<-colnames(pop_diet)
          vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
          vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
          vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
          vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
          vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
          vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
          vec_populations_with_species[vec_populations_with_species=="central"] <- 7
          vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
          vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
          vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
          vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
          vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
          vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
          vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
          vec_populations_with_species<-as.numeric(vec_populations_with_species)           
          DF_subpopulation$subpopulations_species_detected_not_detected<-c(NA)
          DF_subpopulation$subpopulations_species_detected_not_detected[DF_subpopulation$vec_subpop_raster %in% vec_populations_with_species] <- 1               
          #We select the subpopulations where the species is detected: 
          DF_supop_diet_present<- subset(DF_subpopulation, subpopulations_species_detected_not_detected==1)
          DF_supop_diet_NO_present<- subset(DF_subpopulation, is.na(subpopulations_species_detected_not_detected))
          DF_supop_diet_NO_present$habitat_continuous<-c(NA)
          supopulations_present<- unique(DF_supop_diet_present$vec_subpop_raster)
          supopulations_NO_present<- unique(DF_supop_diet_NO_present$vec_subpop_raster)
          #For ALL the subpopulations where it has been dtectd the species:
          #We load the prediction of habitat for current scenario:
          #We load the essemble forecast in continuous:
          load(file="Sps_esemble_curr_proj.RData")
          load(file="Sps_essemble_models.RData")
          #str(Sps_essemble_models)
          Model_Evaluation_DF<-(Sps_essemble_models@ em.models$resp.name_EMcaByTSS_mergedAlgo_mergedRun_mergedData@ model_evaluation)
          save(Model_Evaluation_DF, file="Model_Evaluation_DF.RData")
          habitat_continuous<-(Sps_esemble_curr_proj@ proj@ val[,1])
          DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          dir.create(path_habitat)
          dir.create(path_scenario)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:NROW(supopulations_present)){#The number of species to model
            subpopulation_in_loop<-supopulations_present[j]  
            print("             ########           Starting a new Subpopulation             ########          ")
            print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
            sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
            assign(paste0("sub_habitat_subpop_",subpopulation_in_loop),sub_habitat_subpop$habitat_continuous)
            save(list=paste0("sub_habitat_subpop_", subpopulation_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_in_loop, ".Rdata"))   
          }  
          #For each the subpopulations where NO has been detected the species:
          for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
            subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
            print("             ########           Starting a new Subpopulation             ########          ")
            print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
            sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
            sub_habitat_subpop<-sub_habitat_subpop_no_present
            assign(paste0("sub_habitat_subpop_",subpopulation_no_present_in_loop),sub_habitat_subpop$habitat_continuous)
            save(list=paste0("sub_habitat_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
          }  
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 
                    
      
    #4.5.1.6 We merge for each subpopulations the predictions for all species
      rm(list=ls())
      path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("G:/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_current<-c()
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sub_habitat_subpop_1_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
      sub_habitat_subpop_2_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
      sub_habitat_subpop_3_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
      sub_habitat_subpop_4_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
      sub_habitat_subpop_5_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
      sub_habitat_subpop_6_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
      sub_habitat_subpop_7_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
      sub_habitat_subpop_8_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
      sub_habitat_subpop_9_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
      sub_habitat_subpop_10_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
      sub_habitat_subpop_11_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
      sub_habitat_subpop_12_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
      sub_habitat_subpop_13_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
      sub_habitat_subpop_14_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
      sub_habitat_list_sps_loop<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=50 we consider that we have not enought data
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        if ((list_sps_in_loop$N_pixels_presence)>=50){
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          sub_habitat_list_sps_loop<-cbind(sub_habitat_list_sps_loop,i)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:14){#The number of species to model
            load(file=paste0(path_scenario,"/sub_habitat_subpop_", j, ".Rdata"))
          }  
          sub_habitat_subpop_1_ALL<-cbind(sub_habitat_subpop_1_ALL, sub_habitat_subpop_1)
          sub_habitat_subpop_2_ALL<-cbind(sub_habitat_subpop_2_ALL, sub_habitat_subpop_2)
          sub_habitat_subpop_3_ALL<-cbind(sub_habitat_subpop_3_ALL, sub_habitat_subpop_3)
          sub_habitat_subpop_4_ALL<-cbind(sub_habitat_subpop_4_ALL, sub_habitat_subpop_4)
          sub_habitat_subpop_5_ALL<-cbind(sub_habitat_subpop_5_ALL, sub_habitat_subpop_5)
          sub_habitat_subpop_6_ALL<-cbind(sub_habitat_subpop_6_ALL, sub_habitat_subpop_6)
          sub_habitat_subpop_7_ALL<-cbind(sub_habitat_subpop_7_ALL, sub_habitat_subpop_7)
          sub_habitat_subpop_8_ALL<-cbind(sub_habitat_subpop_8_ALL, sub_habitat_subpop_8)
          sub_habitat_subpop_9_ALL<-cbind(sub_habitat_subpop_9_ALL, sub_habitat_subpop_9)
          sub_habitat_subpop_10_ALL<-cbind(sub_habitat_subpop_10_ALL, sub_habitat_subpop_10)
          sub_habitat_subpop_11_ALL<-cbind(sub_habitat_subpop_11_ALL, sub_habitat_subpop_11)
          sub_habitat_subpop_12_ALL<-cbind(sub_habitat_subpop_12_ALL, sub_habitat_subpop_12)
          sub_habitat_subpop_13_ALL<-cbind(sub_habitat_subpop_13_ALL, sub_habitat_subpop_13)
          sub_habitat_subpop_14_ALL<-cbind(sub_habitat_subpop_14_ALL, sub_habitat_subpop_14)
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 
      sub_habitat_list_sps_loop<-c("ID_pixel",sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_1_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_2_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_3_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_4_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_5_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_6_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_7_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_8_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_9_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_10_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_11_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_12_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_13_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_14_ALL)<-c(sub_habitat_list_sps_loop)
      sub_habitat_Europe_current<-rbind(sub_habitat_subpop_1_ALL,
                sub_habitat_subpop_2_ALL,
                sub_habitat_subpop_3_ALL,
                sub_habitat_subpop_4_ALL,
                sub_habitat_subpop_5_ALL,
                sub_habitat_subpop_6_ALL,
                sub_habitat_subpop_7_ALL,
                sub_habitat_subpop_8_ALL,
                sub_habitat_subpop_9_ALL,
                sub_habitat_subpop_10_ALL,
                sub_habitat_subpop_11_ALL,
                sub_habitat_subpop_12_ALL,
                sub_habitat_subpop_13_ALL,
                sub_habitat_subpop_14_ALL)
      setwd(path_data)
      save(sub_habitat_Europe_current, file=paste0(path_save,"sub_habitat_Europe_current.RData"))

    #4.5.1.7 We divide the dataframe with all species in smaller dataframes for do the merge
      rm(list=ls())
      path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load(file=paste0(path_save,"sub_habitat_Europe_current.RData"))# Copy here to try to be faster
      str(sub_habitat_Europe_current)
      load("D:/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_current_Matrix.RData")# Copy here to try to be faster
      DF_ID_pixel<-as.data.frame(sub_DF_Scen_current_Matrix$ID_pixel)
      colnames(DF_ID_pixel)<-c("ID_pixel")
      save(DF_ID_pixel, file=paste0(path_save,"DF_ID_pixel.RData"))
      sub_habitat_Europe_current_1<- sub_habitat_Europe_current[,1:41]
      sub_habitat_Europe_current_2<- sub_habitat_Europe_current[,c(1,42:82)]
      sub_habitat_Europe_current_3<- sub_habitat_Europe_current[,c(1,83:123)]
      sub_habitat_Europe_current_4<- sub_habitat_Europe_current[,c(1,124:164)]
      sub_habitat_Europe_current_5<- sub_habitat_Europe_current[,c(1,165:205)]
      sub_habitat_Europe_current_6<- sub_habitat_Europe_current[,c(1,206:237)]
      save(sub_habitat_Europe_current_1, file=paste0(path_save,"sub_habitat_Europe_current_1.RData"))
      save(sub_habitat_Europe_current_2, file=paste0(path_save,"sub_habitat_Europe_current_2.RData"))
      save(sub_habitat_Europe_current_3, file=paste0(path_save,"sub_habitat_Europe_current_3.RData"))
      save(sub_habitat_Europe_current_4, file=paste0(path_save,"sub_habitat_Europe_current_4.RData"))
      save(sub_habitat_Europe_current_5, file=paste0(path_save,"sub_habitat_Europe_current_5.RData"))
      save(sub_habitat_Europe_current_6, file=paste0(path_save,"sub_habitat_Europe_current_6.RData"))

    #4.5.1.8 We do the merge for each small data frame
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_1.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_1<-merge(DF_ID_pixel, sub_habitat_Europe_current_1, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_1)
          head(merged_habitat_Europe_Current_1)
          save(merged_habitat_Europe_Current_1, file=paste0(path_save,"merged_habitat_Europe_Current_1.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_2.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_2<-merge(DF_ID_pixel, sub_habitat_Europe_current_2, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_2)
          head(merged_habitat_Europe_Current_2)
          save(merged_habitat_Europe_Current_2, file=paste0(path_save,"merged_habitat_Europe_Current_2.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_3.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_3<-merge(DF_ID_pixel, sub_habitat_Europe_current_3, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_3)
          head(merged_habitat_Europe_Current_3)
          save(merged_habitat_Europe_Current_3, file=paste0(path_save,"merged_habitat_Europe_Current_3.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_4.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_4<-merge(DF_ID_pixel, sub_habitat_Europe_current_4, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_4)
          head(merged_habitat_Europe_Current_4)
          save(merged_habitat_Europe_Current_4, file=paste0(path_save,"merged_habitat_Europe_Current_4.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_5.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_5<-merge(DF_ID_pixel, sub_habitat_Europe_current_5, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_5)
          head(merged_habitat_Europe_Current_5)
          save(merged_habitat_Europe_Current_5, file=paste0(path_save,"merged_habitat_Europe_Current_5.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_current_6.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_Current_6<-merge(DF_ID_pixel, sub_habitat_Europe_current_6, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_habitat_Europe_Current_6)
          head(merged_habitat_Europe_Current_6)
          save(merged_habitat_Europe_Current_6, file=paste0(path_save,"merged_habitat_Europe_Current_6.RData"))
        
    ##########################################################################################
    #4.5.1.9 Create the raster files with the predicitons
    ##########################################################################################
        rm(list=ls())
          library(dplyr)
          library(biomod2) 
          library(slam)
          library(MASS)
          library(usdm)
          library(rgdal) 
          library(caret)
          library(MuMIn)
          library(raster)
          library(rgdal) 
          library(raster) 
          path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          #We load the satellite image for background maps created in part 3  
          load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
          #As these files are super heavy we will load one, then run the code with the species of the dataframe,
          #then close R and load other dataframe: 
            #load(file=paste0(path_save,"merged_habitat_Europe_Current_1.RData"))# 
            #load(file=paste0(path_save,"merged_habitat_Europe_Current_2.RData"))# 
            #load(file=paste0(path_save,"merged_habitat_Europe_Current_3.RData"))# From 92 to 138
            #load(file=paste0(path_save,"merged_habitat_Europe_Current_4.RData"))# From 139 to 188
            #load(file=paste0(path_save,"merged_habitat_Europe_Current_5.RData"))# From 189 to 240
            load(file=paste0(path_save,"merged_habitat_Europe_Current_6.RData"))#  From 242 to 276
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("G:/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sps_in_loop<-c()
      X11()
      for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          if ((i)<46){
          merged_habitat_Europe_Current<-merged_habitat_Europe_Current_1
          }
          if (i>=46 && i<92){ 
          merged_habitat_Europe_Current<-merged_habitat_Europe_Current_2
          }
         if (i>=92 && i<139){ 
          merged_habitat_Europe_Current<-merged_habitat_Europe_Current_3
          }
         if (i>=139 && i<189){  
          merged_habitat_Europe_Current<-merged_habitat_Europe_Current_4
          }
         if (i>=189 && i<241){
         merged_habitat_Europe_Current<-merged_habitat_Europe_Current_5
          }
         if (i>=241 && i<277){  
          merged_habitat_Europe_Current<-merged_habitat_Europe_Current_6
          }
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
            #We save into a raster file to see in Arcgis
                new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                #We defined the projection of our raster:
                projection(new_raster) <- newproj
                values(new_raster)<-merged_habitat_Europe_Current[,as.character(i)]
                writeRaster(new_raster, paste0(path_scenario,"/habitat.img"), overwrite=TRUE)
           #We plot a map with the google satellite image and we save it in jpeg         
            try({
            plotRGB(pr2_crop) 
            })
            plot(new_raster, add = T,col=(gray.colors(12)), colNA=NA)  
            try({#col = gray.colors(12)
            dev.copy(pdf,file=paste0(path_scenario,"/Habitat_Current_Species_",i,".pdf",sep=""))
            })
            dev.off()
          #We save it in a kml file to see in google earth      
            #RWe crate a new rster i lat long: 
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            habitat_lat_long <- projectRaster(new_raster, new_raster_2, method='bilinear')
            KML(habitat_lat_long, file=paste0(path_scenario,"/habitat.kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)
                  print(paste0("             ######## Maps Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
                  }        
            } 
      
    ##########################################################################################
    #4.5.1.10 Code for copy the results to folders
    ##########################################################################################
      rm(list=ls())
      install.packages("pdftools")
      library(pdftools)
      destination_folder_eval = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Eval"
      destination_folder_var_imp_mean = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Var_imp_mean"
      destination_folder_var_imp = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Var_imp"
      destination_folder_kml = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_kml"
      destination_folder_jpeg = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_jpeg"
      destination_folder_pdf = "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_pdf"
      path_data = "G:/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("G:/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.kml"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.png"), paste0(destination_folder_kml,"/")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".img")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_current_sps",i,".img.aux.xml")) 
          file.rename(paste0(destination_folder_kml,"/habitat.kml"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".kml")) 
          file.rename(paste0(destination_folder_kml,"/habitat.png"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".png")) 
          file.copy(paste0(path_save_sps_Biomod,"/eval.csv"), paste0(destination_folder_eval,"/")) 
          file.rename(paste0(destination_folder_eval,"/eval.csv"), paste0(destination_folder_eval,"/eval_sps_",i,".csv")) 
          file.copy(paste0(path_save_sps_Biomod,"/var_imp_df.csv"), paste0(destination_folder_var_imp,"/")) 
          file.rename(paste0(destination_folder_var_imp,"/var_imp_df.csv"), paste0(destination_folder_var_imp,"/var_imp_sps_",i,".csv")) 
          file.copy(paste0(path_save_sps_Biomod,"/var_imp_df_mean.csv"), paste0(destination_folder_var_imp_mean,"/")) 
          file.rename(paste0(destination_folder_var_imp_mean,"/var_imp_df_mean.csv"), paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")) 
        setwd(destination_folder_pdf)
        pdf_convert(paste0(path_scenario,"/Habitat_Current_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_Current_Species_",i,".jpeg"),
        dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        

          }
       }
       
      ###############################################################################################################
      #4.5.1.10.1 Code for copy the results to folders for WILD species
      ###############################################################################################################
        rm(list=ls())
        library(pdftools)
        destination_folder_eval = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Eval"
        destination_folder_var_imp_mean = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Var_imp_mean"
        destination_folder_var_imp = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Var_imp"
        destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_kml"
        destination_folder_jpeg = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_jpeg"
        destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_pdf"
        path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD")
        setwd(path_data)
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        # We load the information of wild/human species  
        load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
        #We change the value of Human/Wild for the reindeer
        merged_sps_list_diet_category3[214,]
        merged_sps_list_diet_category3[214,4] <- 0
        merged_sps_list_diet_category3[214,]
        #We change the value of Human/Wild for Ficus carica
        merged_sps_list_diet_category3[276,4] <- 1
        merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)
         for (i in 1:276){#The number of species to model 276
         list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
         #Condition of a minimum number of presences below n=150 we consider that we have not enought data
         #Condition of a minimum number of presences below n=150 we consider that we have not enought data
         list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
         #Here we are going to load the data in what supopulations the species is eaten and 
         #we will use the subpoplations to predict only in these areas:
         sps_in_loop<-factor(list_sps_in_loop$Species)
         sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
         #We have some species with energy as zero
         try({
         row_sum<-rowSums(sub_subpopulations_diet)
         })
         #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
         #the species when there is no energy of the species in any subpopulation. This happen with some species
         #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
         list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
         list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
         #Conditional for skype species with errors and the brown bear in the models:
         if ((i==109|i==164|i==257)){  
         list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
          }
            if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
            list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
            path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
            path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
            #setwd(path_save_sps_Biomod)
            path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
            path_scenario = paste0(path_habitat,"/Current")
             file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
             file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
             file.copy(paste0(path_scenario,"/habitat.kml"), paste0(destination_folder_kml,"/")) 
             file.copy(paste0(path_scenario,"/habitat.png"), paste0(destination_folder_kml,"/")) 
               file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".img")) 
               file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_current_sps",i,".img.aux.xml")) 
               file.rename(paste0(destination_folder_kml,"/habitat.kml"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".kml")) 
               file.rename(paste0(destination_folder_kml,"/habitat.png"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".png")) 
            #file.copy(paste0(path_scenario,"/Habitat_Current_Species_",i,".pdf"), paste0(destination_folder_pdf,"/habitat_current_sps_",i,".pdf")) 
            file.copy(paste0(path_save_sps_Biomod,"/eval.csv"), paste0(destination_folder_eval,"/")) 
               file.rename(paste0(destination_folder_eval,"/eval.csv"), paste0(destination_folder_eval,"/eval_sps_",i,".csv")) 
             file.copy(paste0(path_save_sps_Biomod,"/var_imp_df.csv"), paste0(destination_folder_var_imp,"/")) 
               file.rename(paste0(destination_folder_var_imp,"/var_imp_df.csv"), paste0(destination_folder_var_imp,"/var_imp_sps_",i,".csv")) 
             file.copy(paste0(path_save_sps_Biomod,"/var_imp_df_mean.csv"), paste0(destination_folder_var_imp_mean,"/")) 
               file.rename(paste0(destination_folder_var_imp_mean,"/var_imp_df_mean.csv"), paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")) 
          setwd(destination_folder_pdf)
          pdf_convert(paste0(path_scenario,"/Habitat_Current_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_Current_Species_",i,".jpeg"),
          dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
            }
         }
         
  
    ###############################################################################################################
    #4.5.1.11 Code for copy the results to unique tables with all species
    ###############################################################################################################
      rm(list=ls())
      destination_folder_eval = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Eval"
      destination_folder_var_imp_mean = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Var_imp_mean"
      destination_folder_var_imp = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Var_imp"
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_kml"
      destination_folder_jpeg = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_jpeg"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Current_prediction_pdf"
      destination_folder_Results_summary_tables = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Results_summary_tables"
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_current.RData")# Copy here to try to be faster
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      merged_var_imp_mean<-as.data.frame(colnames(sub_Scenario_current))
      colnames(merged_var_imp_mean)<-c("Variable")      
      eval_DF_all<-c()     
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")
          var_imp_mean_DF <- read.csv(paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")) # header = TRUE, na.strings = "<NA>"
          colnames(var_imp_mean_DF)<-c("Variable",i)
          merged_var_imp_mean<-merge(merged_var_imp_mean, var_imp_mean_DF, by.x="Variable", by.y = "Variable",all.x = T)
          eval_DF <- read.csv(paste0(destination_folder_eval,"/eval_sps_",i,".csv")) # header = TRUE, na.strings = "<NA>"
          eval_DF<-cbind("i"=i,eval_DF)
          eval_DF_all<-rbind(eval_DF_all,eval_DF[1,c(1,3:9)])
          }
       }
      merged_var_imp_mean2<-merged_var_imp_mean[,c(2:237)]   
      rownames(merged_var_imp_mean2)<-merged_var_imp_mean[,c(1)] 
      Means_merged_var_imp_mean2<-as.data.frame(rowMeans(merged_var_imp_mean2, na.rm = TRUE))
      colnames(Means_merged_var_imp_mean2)<-c("Mean_importance")
      Means_merged_var_imp_mean3<-as.matrix(Means_merged_var_imp_mean2)
      sorted_variables<-sort(Means_merged_var_imp_mean3[,1],decreasing = T)
      DF_sorted_variables<-as.data.frame(sorted_variables)
      DF_summary_eval_DF_all<-summary(eval_DF_all)    
      summary_eval<-DF_summary_eval_DF_all[,c(4,6:8)]
      #Variable importance
      merged_var_imp_mean3<-t(merged_var_imp_mean2)
      DF_i<-as.data.frame(as.numeric(rownames(merged_var_imp_mean3)))
      colnames(DF_i)<-c("i")
      merged_var_imp_mean4<-cbind(merged_var_imp_mean3,DF_i)
      merged_var_imp_mean4_sps<-merge(merged_var_imp_mean4, df_datos_1_276_description_GBIF, by.x="i", by.y = "i",all.x = T)
      merged_var_imp_mean5_sps<-merged_var_imp_mean4_sps[,c(1,13,2:12)]
      write.xlsx(merged_var_imp_mean5_sps, file=paste0(destination_folder_Results_summary_tables,"/merged_var_imp_mean5_sps.xlsx"))# Supplemenry Table 35. Variable importance for all species
      print(paste0("OUTPUT: Supplemenry Table 35 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplemenry Table 35. Supplementary Table 35. Variable importance for the fitted species distribution models of food species. NA means that this variable was not 
        included in the model. We show the identificator of the species (i). ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

      merged_var_imp_mean5_sps[is.na(merged_var_imp_mean5_sps)]<-0
      Means_var_imp<-as.data.frame(colMeans(merged_var_imp_mean5_sps[,c(3:13)], na.rm = TRUE))
      colnames(Means_var_imp)<-c("Variable importance")
      sum(Means_var_imp[1:4,1])
      sum(Means_var_imp[5:11,1])
      write.xlsx(Means_var_imp, file=paste0(destination_folder_Results_summary_tables,'/Means_var_imp.xlsx'))# Supplemenry Table 36. Summary of variable importance for all species

      #Evaluation
      eval_DF_all2<-merge(eval_DF_all, df_datos_1_276_description_GBIF, by.x="i", by.y = "i",all.x = T)
      eval_DF_all3<-eval_DF_all2[,c(1,9,3,4,6:8)]
      write.xlsx(eval_DF_all3, file=paste0(destination_folder_Results_summary_tables,'/eval_DF_all3.xlsx'))# Supplemenry Table 33.  Evaluation of model for all  species
      print(paste0("OUTPUT: Supplemenry Table 33.   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplemenry Table 33. Supplemenry Table 33. List of the 236 species in the diet (including wild and human origin species) of the brown bear for which it has been fitted 
      a species distribution model (notice that we finally only considered  wild species for calculation of biotic interactions). We show the identificator of species (i), the scientific name 
        of the species (Species), the TSS value, and the sensivity and specifity of the modes.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

      summary_eval<-summary(eval_DF_all3[,c(4:7)])
      write.xlsx(summary_eval, file=paste0(destination_folder_Results_summary_tables,'/summary_eval.xlsx'))# Supplemenry Table 34. Summary of evaluation of model for all species
      print(paste0("OUTPUT: Supplemenry Table 34. #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplemenry Table 34.Summary of evaluation of models for all species  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))


    ############################################################################################        
    #4.5.1.11.1 Code for copy the results to unique tables with ONLY WILD species
    ############################################################################################        
      rm(list=ls())
      library(xlsx)
      destination_folder_eval = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Eval"
      destination_folder_var_imp_mean = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Var_imp_mean"
      destination_folder_var_imp = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Var_imp"
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_kml"
      destination_folder_jpeg = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_jpeg"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_pdf"
      destination_folder_Results_summary_tables = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Results_summary_tables"
      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_current.RData")# Copy here to try to be faster
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)
      merged_var_imp_mean<-as.data.frame(colnames(sub_Scenario_current))
      colnames(merged_var_imp_mean)<-c("Variable")      
      eval_DF_all<-c()     
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/Current")
          paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")
          var_imp_mean_DF <- read.csv(paste0(destination_folder_var_imp_mean,"/var_imp_mean_sps_",i,".csv")) # header = TRUE, na.strings = "<NA>"
          colnames(var_imp_mean_DF)<-c("Variable",i)
          merged_var_imp_mean<-merge(merged_var_imp_mean, var_imp_mean_DF, by.x="Variable", by.y = "Variable",all.x = T)
          eval_DF <- read.csv(paste0(destination_folder_eval,"/eval_sps_",i,".csv")) # header = TRUE, na.strings = "<NA>"
          eval_DF<-cbind("i"=i,eval_DF)
          eval_DF_all<-rbind(eval_DF_all,eval_DF[1,c(1,3:9)])
          }
       }
      merged_var_imp_mean2<-merged_var_imp_mean[,c(2:205)]   
      rownames(merged_var_imp_mean2)<-merged_var_imp_mean[,c(1)] 
      Means_merged_var_imp_mean2<-as.data.frame(rowMeans(merged_var_imp_mean2, na.rm = TRUE))
      colnames(Means_merged_var_imp_mean2)<-c("Mean_importance")
      Means_merged_var_imp_mean3<-as.matrix(Means_merged_var_imp_mean2)
      sorted_variables<-sort(Means_merged_var_imp_mean3[,1],decreasing = T)
      DF_sorted_variables<-as.data.frame(sorted_variables)
      DF_summary_eval_DF_all<-summary(eval_DF_all)    
      summary_eval<-DF_summary_eval_DF_all[,c(4,6:8)]
      
      #Variable importance for only wild species
      merged_var_imp_mean3<-t(merged_var_imp_mean2)
      DF_i<-as.data.frame(as.numeric(rownames(merged_var_imp_mean3)))
      colnames(DF_i)<-c("i")
      merged_var_imp_mean4<-cbind(merged_var_imp_mean3,DF_i)
      merged_var_imp_mean4_sps<-merge(merged_var_imp_mean4, df_datos_1_276_description_GBIF, by.x="i", by.y = "i",all.x = T)
      merged_var_imp_mean5_sps<-merged_var_imp_mean4_sps[,c(1,13,2:12)]
      merged_var_imp_mean5_sps_wild<-merged_var_imp_mean5_sps
      write.xlsx(merged_var_imp_mean5_sps_wild, file=paste0(destination_folder_Results_summary_tables,"/merged_var_imp_mean5_sps_wild.xlsx"))# 
      
        #Summary of variable importance for only wild species
        merged_var_imp_mean5_sps[is.na(merged_var_imp_mean5_sps)]<-0
        Means_var_imp<-as.data.frame(colMeans(merged_var_imp_mean5_sps[,c(3:13)], na.rm = TRUE))
        colnames(Means_var_imp)<-c("Variable importance")
        sum(Means_var_imp[1:4,1])
        sum(Means_var_imp[5:11,1])
        Means_var_imp_wild<-Means_var_imp
        write.xlsx(Means_var_imp_wild, file=paste0(destination_folder_Results_summary_tables,'/Means_var_imp_wild.xlsx'))#
        
      #Evaluation of models for only wild species
      eval_DF_all2<-merge(eval_DF_all, df_datos_1_276_description_GBIF, by.x="i", by.y = "i",all.x = T)
      eval_DF_all3<-eval_DF_all2[,c(1,9,3,4,6:8)]
      eval_DF_all3_wild<-eval_DF_all3
      write.xlsx(eval_DF_all3_wild, file=paste0(destination_folder_Results_summary_tables,'/eval_DF_all3_wild.xlsx'))# 
      
        #Summary of evaluation of models for only wild species
        summary_eval_wild<-summary(eval_DF_all3[,c(4:7)])
        write.xlsx(summary_eval_wild, file=paste0(destination_folder_Results_summary_tables,'/summary_eval_wild.xlsx'))# 
        
      
  ##########################################################################################################
  ##########################################################################################################
  #4.5.2 Ensemble Forecasting of RCP_25 Scenario 
  ##########################################################################################################
  ##########################################################################################################
        
        
    ##########################################################################################################
    #4.5.2.1 We merge the data of the matrix to the RCP_25 conditions 
    ##########################################################################################################

      #4.5.2.1 We merge the data of the matrix to the RCP_25 conditions          
        rm(list=ls()) 
        path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD/")
        setwd("H:/D/Test_Biomod")
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_RCP26.RData")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/data_matrix.RData")
        DF_Scen_RCP26_Matrix<-as.data.frame(cbind(data_matrix,sub_Scenario_RCP26))
        save(DF_Scen_RCP26_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix.RData")
        DF_Scen_RCP26_Matrix_noNA<-na.omit(DF_Scen_RCP26_Matrix)
        save(DF_Scen_RCP26_Matrix_noNA, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix_noNA.RData")

      #4.5.2.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
        rm(list=ls()) 
        path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD/")
        setwd("H:/D/Test_Biomod")
        #Map extrapolacion
        load("H:/D/Test_Biomod/vec_extrap_raster.RData")
        #Map of subpopulations
        load("H:/D/Test_Biomod/vec_subpop_raster.RData")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix.RData")# Copy here to try to be faster
        DF_Scen_RCP26_Matrix_2<-as.data.frame(cbind(DF_Scen_RCP26_Matrix,vec_extrap_raster,vec_subpop_raster))
        DF_Scen_RCP26_Matrix_noNA_2<-na.omit(DF_Scen_RCP26_Matrix_2)
        save(DF_Scen_RCP26_Matrix_noNA_2, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix_noNA_2.RData")

      #4.5.2.3 Application of BIOMOD_EnsembleForecasting to RCP26 Scenario    
        rm(list=ls()) 
        library(biomod2) 
        library(slam)
        library(MASS)
        library(usdm)
        library(rgdal) 
        library(caret)
        library(MuMIn)
        library(raster)
        path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD/")
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix_noNA_2.RData")# Copy here to try to be faster
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        merge_list_diet_energy<-merge(df_datos_1_276_description_GBIF,subpopulations_diet, by.x="Species", by.y = "species",all.x = T)
        merge_list_diet_energy2<-merge(subpopulations_diet,df_datos_1_276_description_GBIF, by.x="species", by.y = "Species",all.x = T)
        sub_Na_merge_list_diet_energy<-merge_list_diet_energy[is.na(merge_list_diet_energy$central),]
        sub_Na_merge_list_diet_energy2<-merge_list_diet_energy2[is.na(merge_list_diet_energy2$N_pixels_presence),]
        #Descripcion de los NAs:
                  #Para sub_Na_merge_list_diet_energy tenemos tres especies con NA que son Fagus orientalis, Lonicera caucasica y Pahseolus vulgaris que se les 
                  #ha considerado rF_original = 0 en el estudio por lo que la REDEC sale NA, aunque la especie est? descrita coo parte de la dieta
                  #Para sub_Na_merge_list_diet_energy2 tenemos un caso que es el Sus scrofa domesticus que est? en el estudio de Naves y no lo tenemos invluido como 
                  #especie/subes pecie o taxon, con este sumarian 277 taxones
        table_range_environmental_data_RCP26_scenario_all<-c()
        X11()
          for (i in 1:276){#The number of species to model
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          #Condition of a minimum number of presences below n=150 we consider that we have not enought data
          #Condition of a minimum number of presences below n=150 we consider that we have not enought data
          #Condition for i=109 # With this species we have no model ! No models kept due to threshold filtering... Ensemble Modeling was skipped!
          list_sps_in_loop$N_pixels_presence[(list_sps_in_loop$i)==109]<-0
          list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sps_in_loop<-factor(list_sps_in_loop$Species)
          sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
          #We have some species with energy as zero
          try({
            row_sum<-rowSums(sub_subpopulations_diet)
          })
          #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
          #the species when there is no energy of the species in any subpopulation. This happen with some species
          #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
          list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
          list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
            if ((list_sps_in_loop$N_pixels_presence)>=50){  
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
            path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i,"/")
            path_save_Univariable_models = paste0(path_save_sps,"/Univariable_models/")
            path_save_sps_Biomod = paste0(path_save_sps,"/Biomod/")
            setwd(path_save_sps_Biomod)
            #This is the dataset with the predictors used in the model "Predictors_RCP26" these data are unique for each species
            #We select from the dataframe with the extension to project the variables included in Predictors_RCP26"
            #load(file=paste(path_save_Univariable_models,"variables_names_less_AICc.RData"))
            load(file=paste0(path_save_sps_Biomod,"Sps_essemble_models.RData"))  
            variables_names<-Sps_essemble_models@ expl.var.names
            #Here we are going to load the data in what supopulations the species is eaten and 
            #we will use the subpoplations to predict only in these areas:
            sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
            pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
            vec_populations_with_species<-colnames(pop_diet)
            vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
            vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
            vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
            vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
            vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
            vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
            vec_populations_with_species[vec_populations_with_species=="central"] <- 7
            vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
            vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
            vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
            vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
            vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
            vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
            vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
            vec_populations_with_species<-as.numeric(vec_populations_with_species)           
            data_pre_for_RCP26_projection_0<-DF_Scen_RCP26_Matrix_noNA_2[DF_Scen_RCP26_Matrix_noNA_2$vec_subpop_raster %in% vec_populations_with_species,]
            data_pre_for_RCP26_projection<-data_pre_for_RCP26_projection_0[,variables_names]
            predicted_area_absolute_RCP26<-NROW(data_pre_for_RCP26_projection)
            var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
            var_imp_df_mean<-data.frame(apply(var_imp_df,1,mean))
            colnames(var_imp_df_mean)<-c("Var_importance")
            #Code for see the range of the data of variables used to fit the model and compare with the range of data of each scenario:
            #We load the data used to fit the model of the species
            load(file=paste0(path_save_sps_Biomod,"format_data.RData"))
            df_data_fit<-format_data @ data.env.var
            data_var_all<-c()
            pdf(file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Plots_environ_RCP26/Range_env_data_sps_",i,".pdf",sep=""), width = 0.394 * 64, height = 0.394 * 8)
            par(mfrow = c(1, 6), cex = 0.7)
              for (col in 1:ncol(df_data_fit)) {
              d_all_europe <- density(data_pre_for_RCP26_projection[,col])
              plot(d_all_europe, main=colnames(df_data_fit[col]))
              d <- density(df_data_fit[,col])
              lines(d, main=colnames(df_data_fit[col]))
              polygon(d, col="red", border="red")
              abline(v = max(df_data_fit[,col]), lty=2, col="blue")
              abline(v = min(df_data_fit[,col]), lty=2, col="blue")
              lines(d_all_europe, main=colnames(df_data_fit[col]))
              range_var_max<-max(df_data_fit[,col])
              range_var_min<-min(df_data_fit[,col])
              range_var<-max(df_data_fit[,col])-min(df_data_fit[,col])
              range_var_plus5per<-range_var_min + range_var*1.2
              range_var_less5per<-range_var_min-range_var*0.2
              data_in_range<-with(data_pre_for_RCP26_projection, data_pre_for_RCP26_projection[,col] <= range_var_plus5per & 
                data_pre_for_RCP26_projection[,col] >= range_var_less5per)
              tab_var_europe<-table(data_in_range)
              df_tab_var_europe<-as.data.frame(tab_var_europe)
              per_n_cells_in_range<-df_tab_var_europe[which(df_tab_var_europe$data_in_range=="TRUE"),"Freq"]/6252207*100
              data_var<-data.frame(i,list_sps_in_loop$Species,colnames(df_data_fit[col]),range_var_min,range_var_max,range_var,range_var_less5per,range_var_plus5per,per_n_cells_in_range)
              data_var_all<-rbind(data_var_all,data_var)
              }
            colnames(data_var_all)<-c("ID_species","Species","Variable","range_var_min","range_var_max","range_var",
            "range_var_less5per","range_var_plus5per","per_n_cells_in_range")
            dev.off()  
            #Using the range of data of table_range_environmental_data we are going to convert the data ouside the range to NA
              for (j in 1:6) {
              data_pre_for_RCP26_projection[,j][data_pre_for_RCP26_projection[,j]>data_var_all[j,"range_var_plus5per"]]<-NA
              data_pre_for_RCP26_projection[,j][data_pre_for_RCP26_projection[,j]<data_var_all[j,"range_var_less5per"]]<-NA
              }
            #BIOMOD_EnsembleForecasting RCP26 copy from WILFRED corrections
            Sps_esemble_RCP26_proj<-BIOMOD_EnsembleForecasting(EM.output=Sps_essemble_models,
            new.env=data_pre_for_RCP26_projection,
            proj.name="RCP26",
            selected.models = 'all', 
            binary.meth = "TSS", 
            do.stack=FALSE,
            build.clamping.mask = F)
            save(Sps_esemble_RCP26_proj, file="Sps_esemble_RCP26_proj.RData")
            not_predicted_RCP26<-sum(is.na(Sps_esemble_RCP26_proj@ proj @ val))
            predicted_area_RCP26<-(predicted_area_absolute_RCP26-not_predicted_RCP26)/predicted_area_absolute_RCP26*100
            save(predicted_area_RCP26, file=paste0(path_save_sps_Biomod,"predicted_area_RCP26_sps",i,".RData"))
            table_range_environmental_data_RCP26_scenario<-data.frame(data_var_all,var_imp_df_mean,"Predict_area"=predicted_area_RCP26)
            save(table_range_environmental_data_RCP26_scenario, file=paste0(path_save_sps_Biomod,"table_range_environmental_data_RCP26_scenario_sps_",i,".RData"))
            table_range_environmental_data_RCP26_scenario_all<-rbind(table_range_environmental_data_RCP26_scenario_all,table_range_environmental_data_RCP26_scenario)   
            print(paste0("Precentage of predicted area for species ",list_sps_in_loop$Species," = ",predicted_area_RCP26))
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
            }        
          }
        setwd(path_data)
        save(table_range_environmental_data_RCP26_scenario_all, file=paste0(path_save,"table_range_environmental_data_RCP26_scenario_all.RData"))
        View(table_range_environmental_data_RCP26_scenario_all)
        
      #4.5.2.4 Code for create a ID for the cells of raster
        rm(list=ls()) 
        path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD/")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix.RData")# Copy here to try to be faster
        sub_DF_Scen_RCP26_Matrix<-DF_Scen_RCP26_Matrix[,c(2:3)]
        #We create a small dataframe with the ID of the cells of the raster
        save(sub_DF_Scen_RCP26_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP26_Matrix.RData")

     #4.5.2.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
      rm(list=ls())
      library(biomod2) 
      library(slam)
      library(MASS)
      library(usdm)
      library(rgdal) 
      library(caret)
      library(MuMIn)
      library(raster)
      library(rgdal) 
      library(raster) 
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      DF_subpopulation<-as.data.frame(DF_Scen_RCP26_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_RCP26<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
          pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
          vec_populations_with_species<-colnames(pop_diet)
          vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
          vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
          vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
          vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
          vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
          vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
          vec_populations_with_species[vec_populations_with_species=="central"] <- 7
          vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
          vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
          vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
          vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
          vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
          vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
          vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
          vec_populations_with_species<-as.numeric(vec_populations_with_species)           
          DF_subpopulation$subpopulations_species_detected_not_detected<-c(NA)
          DF_subpopulation$subpopulations_species_detected_not_detected[DF_subpopulation$vec_subpop_raster %in% vec_populations_with_species] <- 1               
          #We select the subpopulations where the species is detected: 
          DF_supop_diet_present<- subset(DF_subpopulation, subpopulations_species_detected_not_detected==1)
          DF_supop_diet_NO_present<- subset(DF_subpopulation, is.na(subpopulations_species_detected_not_detected))
          DF_supop_diet_NO_present$habitat_continuous<-c(NA)
          supopulations_present<- unique(DF_supop_diet_present$vec_subpop_raster)
          supopulations_NO_present<- unique(DF_supop_diet_NO_present$vec_subpop_raster)
          #For ALL the subpopulations where it has been dtectd the species:
          #We load the prediction of habitat for RCP26 scenario:
          #We load the essemble forecast in continuous:
          load(file="Sps_esemble_RCP26_proj.RData")
          load(file="Sps_essemble_models.RData")
          habitat_continuous<-(Sps_esemble_RCP26_proj@ proj@ val[,1])
          #Now we are going to replace "habitat_continuous" by "habitat_binary"
          DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP26")
          #dir.create(path_habitat)
          dir.create(path_scenario)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:NROW(supopulations_present)){#The number of species to model
            subpopulation_in_loop<-supopulations_present[j]  
            print("             ########           Starting a new Subpopulation             ########          ")
            print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
            sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
            assign(paste0("sub_habitat_subpop_",subpopulation_in_loop),sub_habitat_subpop$habitat_continuous)
            save(list=paste0("sub_habitat_subpop_", subpopulation_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_in_loop, ".Rdata"))   
          }  
          #For each the subpopulations where NO has been detected the species:
          for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
            subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
            print("             ########           Starting a new Subpopulation             ########          ")
            print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
            sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
            sub_habitat_subpop<-sub_habitat_subpop_no_present
            assign(paste0("sub_habitat_subpop_",subpopulation_no_present_in_loop),sub_habitat_subpop$habitat_continuous)
            save(list=paste0("sub_habitat_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
          }  
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 

    #4.5.2.6 We merge for each subpopulations the predictions for all species
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP26_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      DF_subpopulation<-as.data.frame(DF_Scen_RCP26_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_RCP26<-c()
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sub_habitat_subpop_1_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
      sub_habitat_subpop_2_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
      sub_habitat_subpop_3_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
      sub_habitat_subpop_4_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
      sub_habitat_subpop_5_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
      sub_habitat_subpop_6_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
      sub_habitat_subpop_7_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
      sub_habitat_subpop_8_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
      sub_habitat_subpop_9_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
      sub_habitat_subpop_10_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
      sub_habitat_subpop_11_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
      sub_habitat_subpop_12_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
      sub_habitat_subpop_13_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
      sub_habitat_subpop_14_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
      sub_habitat_list_sps_loop<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        # list_sps_in_loop$N_pixels_presence[(i==109)]<-0#with this code we skipe this species with cero energy
        #list_sps_in_loop$N_pixels_presence[(i=164)]<-0#with this code we skipe this species with cero energy
        if ((list_sps_in_loop$N_pixels_presence)>=50){
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP26")
          sub_habitat_list_sps_loop<-cbind(sub_habitat_list_sps_loop,i)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:14){#The number of species to model
            load(file=paste0(path_scenario,"/sub_habitat_subpop_", j, ".Rdata"))
          }  
          sub_habitat_subpop_1_ALL<-cbind(sub_habitat_subpop_1_ALL, sub_habitat_subpop_1)
          sub_habitat_subpop_2_ALL<-cbind(sub_habitat_subpop_2_ALL, sub_habitat_subpop_2)
          sub_habitat_subpop_3_ALL<-cbind(sub_habitat_subpop_3_ALL, sub_habitat_subpop_3)
          sub_habitat_subpop_4_ALL<-cbind(sub_habitat_subpop_4_ALL, sub_habitat_subpop_4)
          sub_habitat_subpop_5_ALL<-cbind(sub_habitat_subpop_5_ALL, sub_habitat_subpop_5)
          sub_habitat_subpop_6_ALL<-cbind(sub_habitat_subpop_6_ALL, sub_habitat_subpop_6)
          sub_habitat_subpop_7_ALL<-cbind(sub_habitat_subpop_7_ALL, sub_habitat_subpop_7)
          sub_habitat_subpop_8_ALL<-cbind(sub_habitat_subpop_8_ALL, sub_habitat_subpop_8)
          sub_habitat_subpop_9_ALL<-cbind(sub_habitat_subpop_9_ALL, sub_habitat_subpop_9)
          sub_habitat_subpop_10_ALL<-cbind(sub_habitat_subpop_10_ALL, sub_habitat_subpop_10)
          sub_habitat_subpop_11_ALL<-cbind(sub_habitat_subpop_11_ALL, sub_habitat_subpop_11)
          sub_habitat_subpop_12_ALL<-cbind(sub_habitat_subpop_12_ALL, sub_habitat_subpop_12)
          sub_habitat_subpop_13_ALL<-cbind(sub_habitat_subpop_13_ALL, sub_habitat_subpop_13)
          sub_habitat_subpop_14_ALL<-cbind(sub_habitat_subpop_14_ALL, sub_habitat_subpop_14)
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 
      sub_habitat_list_sps_loop<-c("ID_pixel",sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_1_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_2_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_3_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_4_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_5_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_6_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_7_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_8_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_9_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_10_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_11_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_12_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_13_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_14_ALL)<-c(sub_habitat_list_sps_loop)
      sub_habitat_Europe_RCP26<-rbind(sub_habitat_subpop_1_ALL,
                sub_habitat_subpop_2_ALL,
                sub_habitat_subpop_3_ALL,
                sub_habitat_subpop_4_ALL,
                sub_habitat_subpop_5_ALL,
                sub_habitat_subpop_6_ALL,
                sub_habitat_subpop_7_ALL,
                sub_habitat_subpop_8_ALL,
                sub_habitat_subpop_9_ALL,
                sub_habitat_subpop_10_ALL,
                sub_habitat_subpop_11_ALL,
                sub_habitat_subpop_12_ALL,
                sub_habitat_subpop_13_ALL,
                sub_habitat_subpop_14_ALL)
      setwd(path_data)
      save(sub_habitat_Europe_RCP26, file=paste0(path_save,"sub_habitat_Europe_RCP26.RData"))

    #4.5.2.7 We divide the dataframe with all species in smaller dataframes for do the merge
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load(file=paste0(path_save,"sub_habitat_Europe_RCP26.RData"))# Copy here to try to be faster
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP26_Matrix.RData")# Copy here to try to be faster
      DF_ID_pixel<-as.data.frame(sub_DF_Scen_RCP26_Matrix$ID_pixel)
      colnames(DF_ID_pixel)<-c("ID_pixel")
      save(DF_ID_pixel, file=paste0(path_save,"DF_ID_pixel.RData"))
      sub_habitat_Europe_RCP26_1<- sub_habitat_Europe_RCP26[,1:41]
      sub_habitat_Europe_RCP26_2<- sub_habitat_Europe_RCP26[,c(1,42:82)]
      sub_habitat_Europe_RCP26_3<- sub_habitat_Europe_RCP26[,c(1,83:123)]
      sub_habitat_Europe_RCP26_4<- sub_habitat_Europe_RCP26[,c(1,124:164)]
      sub_habitat_Europe_RCP26_5<- sub_habitat_Europe_RCP26[,c(1,165:205)]
      sub_habitat_Europe_RCP26_6<- sub_habitat_Europe_RCP26[,c(1,206:237)]
      save(sub_habitat_Europe_RCP26_1, file=paste0(path_save,"sub_habitat_Europe_RCP26_1.RData"))
      save(sub_habitat_Europe_RCP26_2, file=paste0(path_save,"sub_habitat_Europe_RCP26_2.RData"))
      save(sub_habitat_Europe_RCP26_3, file=paste0(path_save,"sub_habitat_Europe_RCP26_3.RData"))
      save(sub_habitat_Europe_RCP26_4, file=paste0(path_save,"sub_habitat_Europe_RCP26_4.RData"))
      save(sub_habitat_Europe_RCP26_5, file=paste0(path_save,"sub_habitat_Europe_RCP26_5.RData"))
      save(sub_habitat_Europe_RCP26_6, file=paste0(path_save,"sub_habitat_Europe_RCP26_6.RData"))

    #4.5.2.8 We do the merge for each small data frame
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_1.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_1<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_1, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_1, file=paste0(path_save,"merged_habitat_Europe_RCP26_1.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_2.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_2<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_2, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_2, file=paste0(path_save,"merged_habitat_Europe_RCP26_2.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_3.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_3<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_3, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_3, file=paste0(path_save,"merged_habitat_Europe_RCP26_3.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_4.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_4<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_4, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_4, file=paste0(path_save,"merged_habitat_Europe_RCP26_4.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_5.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_5<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_5, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_5, file=paste0(path_save,"merged_habitat_Europe_RCP26_5.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP26_6.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP26_6<-merge(DF_ID_pixel, sub_habitat_Europe_RCP26_6, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP26_6, file=paste0(path_save,"merged_habitat_Europe_RCP26_6.RData"))

    #4.5.2.9 Create the raster files with the predicitons
          rm(list=ls())
          library(raster) 
          library(sp) 
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          #We load the satellite image for background maps created in part 3  
          load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
          #As these files are super heavy we will load one, then run the code with the species of the dataframe,
          #then close R and load other dataframe: 
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP26_1.RData"))# From 1 to 45 
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP26_2.RData"))# From 46 to 91
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP26_3.RData"))# From 92 to 138
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP26_4.RData"))# From 139 to 188
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP26_5.RData"))# From 189 to 240
            load(file=paste0(path_save,"merged_habitat_Europe_RCP26_6.RData"))#  From 242 to 276
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sps_in_loop<-c()
      X11()
      for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP26")
          if ((i)<46){
          merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_1
          }
          if (i>=46 && i<92){ 
          merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_2
          }
         if (i>=92 && i<139){ 
          merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_3
          }
         if (i>=139 && i<189){  
          merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_4
          }
         if (i>=189 && i<241){
         merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_5
          }
         if (i>=241 && i<277){  
          merged_habitat_Europe_RCP26<-merged_habitat_Europe_RCP26_6
          }
            #We save into a raster file to see in Arcgis
                new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                #We defined the projection of our raster:
                projection(new_raster) <- newproj
                values(new_raster)<-merged_habitat_Europe_RCP26[,as.character(i)]
                writeRaster(new_raster, paste0(path_scenario,"/habitat.img"), overwrite=TRUE)
            #We plot a map with the google satellite image and we save it in jpeg         
            try({
            plotRGB(pr2_crop) 
            })
            plot(new_raster, add = T,col=(gray.colors(12)), colNA=NA)  
            try({#col = gray.colors(12)
            dev.copy(pdf,file=paste0(path_scenario,"/Habitat_RCP26_Species_",i,".pdf",sep=""))
            })
            dev.off()
            #We save it in a kml file to see in google earth      
            #We crate a new rster i lat lonH:/G 
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            habitat_lat_long <- projectRaster(new_raster, new_raster_2, method='bilinear')
            KML(habitat_lat_long, file=paste0(path_scenario,"/habitat.kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)
                  print(paste0("             ######## Maps Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
                  }        
            } 

    #4.5.2.10 Code for copy the results to folders
      rm(list=ls())
      install.packages("pdftools")
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP26_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP26_prediction_pdf"
      dir.create(destination_folder_kml)
      dir.create(destination_folder_pdf)
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)

       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          #setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP26")

           file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
           file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP26_sps_",i,".img")) 
           file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP26_sps",i,".img.aux.xml")) 
           file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP26_sps_",i,".kmz")) 
        setwd(destination_folder_pdf)
        pdf_convert(paste0(path_scenario,"/Habitat_RCP26_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP26_Species_",i,".jpeg"),
        dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }
      
    #4.5.2.10.1 Code for copy the results to folders for WILD species
      rm(list=ls())
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP26_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP26_prediction_pdf"
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors and the brown bear in the models:
       if ((i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP26")
           file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
             file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP26_sps_",i,".img")) 
             file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP26_sps",i,".img.aux.xml")) 
             file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP26_sps_",i,".kmz")) 
        setwd(destination_folder_pdf)
        pdf_convert(paste0(path_scenario,"/Habitat_RCP26_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP26_Species_",i,".jpeg"),
        dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }


  ###################################################################################################################    
  ###################################################################################################################    
  ##4.5.3 Ensemble Forecasting of RCP60 Scenario 
  ###################################################################################################################    
  ###################################################################################################################    
  
    #4.5.3.1 We merge the data of the matrix to the RCP_25 conditions          
      rm(list=ls()) 
      path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("H:/D/Test_Biomod")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_RCP60.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/data_matrix.RData")
      DF_Scen_RCP60_Matrix<-as.data.frame(cbind(data_matrix,sub_Scenario_RCP60))
      save(DF_Scen_RCP60_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix.RData")
      DF_Scen_RCP60_Matrix_noNA<-na.omit(DF_Scen_RCP60_Matrix)
      save(DF_Scen_RCP60_Matrix_noNA, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA.RData")

    #4.5.3.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      rm(list=ls()) 
      path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("H:/D/Test_Biomod")
      #Map extrapolacion
      load("H:/D/Test_Biomod/vec_extrap_raster.RData")
      #Map of subpopulations
      load("H:/D/Test_Biomod/vec_subpop_raster.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix.RData")# Copy here to try to be faster
      DF_Scen_RCP60_Matrix_2<-as.data.frame(cbind(DF_Scen_RCP60_Matrix,vec_extrap_raster,vec_subpop_raster))
      DF_Scen_RCP60_Matrix_noNA_2<-na.omit(DF_Scen_RCP60_Matrix_2)
      save(DF_Scen_RCP60_Matrix_noNA_2, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA_2.RData")

    #4.5.3.3 Application of BIOMOD_EnsembleForecasting to RCP60 Scenario    
      rm(list=ls()) 
      library(biomod2) 
      library(slam)
      library(MASS)
      library(usdm)
      library(rgdal) 
      library(caret)
      library(MuMIn)
      library(raster)
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      #load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/df_datos_1_276_description_GBIF.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      merge_list_diet_energy<-merge(df_datos_1_276_description_GBIF,subpopulations_diet, by.x="Species", by.y = "species",all.x = T)
      merge_list_diet_energy2<-merge(subpopulations_diet,df_datos_1_276_description_GBIF, by.x="species", by.y = "Species",all.x = T)
      sub_Na_merge_list_diet_energy<-merge_list_diet_energy[is.na(merge_list_diet_energy$central),]
      sub_Na_merge_list_diet_energy2<-merge_list_diet_energy2[is.na(merge_list_diet_energy2$N_pixels_presence),]
      #Descripcion de los NAs:
                #Para sub_Na_merge_list_diet_energy tenemos tres especies con NA que son Fagus orientalis, Lonicera caucasica y Pahseolus vulgaris que se les 
                #ha considerado rF_original = 0 en el estudio por lo que la REDEC sale NA, aunque la especie est? descrita coo parte de la dieta
                #Para sub_Na_merge_list_diet_energy2 tenemos un caso que es el Sus scrofa domesticus que est? en el estudio de Naves y no lo tenemos invluido como 
                #especie/subes pecie o taxon, con este sumarian 277 taxones
      table_range_environmental_data_RCP60_scenario_all<-c()
      X11()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        #Condition of a minimum number of presences below n=150 we consider that we have not enought data
        #Condition for i=109 # With this species we have no model ! No models kept due to threshold filtering... Ensemble Modeling was skipped!
        list_sps_in_loop$N_pixels_presence[(list_sps_in_loop$i)==109]<-0
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i,"/")
          path_save_Univariable_models = paste0(path_save_sps,"/Univariable_models/")
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod/")
          setwd(path_save_sps_Biomod)
          #This is the dataset with the predictors used in the model "Predictors_RCP60" these data are unique for each species
          #We select from the dataframe with the extension to project the variables included in Predictors_RCP60"
          load(file=paste0(path_save_sps_Biomod,"Sps_essemble_models.RData"))  
          variables_names<-Sps_essemble_models@ expl.var.names
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
          pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
          vec_populations_with_species<-colnames(pop_diet)
          vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
          vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
          vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
          vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
          vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
          vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
          vec_populations_with_species[vec_populations_with_species=="central"] <- 7
          vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
          vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
          vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
          vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
          vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
          vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
          vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
          vec_populations_with_species<-as.numeric(vec_populations_with_species)           
          data_pre_for_RCP60_projection_0<-DF_Scen_RCP60_Matrix_noNA_2[DF_Scen_RCP60_Matrix_noNA_2$vec_subpop_raster %in% vec_populations_with_species,]
          data_pre_for_RCP60_projection<-data_pre_for_RCP60_projection_0[,variables_names]
          predicted_area_absolute_RCP60<-NROW(data_pre_for_RCP60_projection)
          var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
          var_imp_df_mean<-data.frame(apply(var_imp_df,1,mean))
          colnames(var_imp_df_mean)<-c("Var_importance")
          #Code for see the range of the data of variables used to fit the model and compare with the range of data of each scenario:
          #We load the data used to fit the model of the species
          load(file=paste0(path_save_sps_Biomod,"format_data.RData"))
          df_data_fit<-format_data @ data.env.var
          data_var_all<-c()
          pdf(file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Plots_environ_RCP60/Range_env_data_sps_",i,".pdf",sep=""), width = 0.394 * 64, height = 0.394 * 8)
          par(mfrow = c(1, 6), cex = 0.7)
            for (col in 1:ncol(df_data_fit)) {
            d_all_europe <- density(data_pre_for_RCP60_projection[,col])
            plot(d_all_europe, main=colnames(df_data_fit[col]))
            d <- density(df_data_fit[,col])
            lines(d, main=colnames(df_data_fit[col]))
            polygon(d, col="red", border="red")
            abline(v = max(df_data_fit[,col]), lty=2, col="blue")
            abline(v = min(df_data_fit[,col]), lty=2, col="blue")
            lines(d_all_europe, main=colnames(df_data_fit[col]))
            range_var_max<-max(df_data_fit[,col])
            range_var_min<-min(df_data_fit[,col])
            range_var<-max(df_data_fit[,col])-min(df_data_fit[,col])
            range_var_plus5per<-range_var_min + range_var*1.2
            range_var_less5per<-range_var_min-range_var*0.2
            data_in_range<-with(data_pre_for_RCP60_projection, data_pre_for_RCP60_projection[,col] <= range_var_plus5per & 
               data_pre_for_RCP60_projection[,col] >= range_var_less5per)
            tab_var_europe<-table(data_in_range)
            df_tab_var_europe<-as.data.frame(tab_var_europe)
            per_n_cells_in_range<-df_tab_var_europe[which(df_tab_var_europe$data_in_range=="TRUE"),"Freq"]/6252207*100
            data_var<-data.frame(i,list_sps_in_loop$Species,colnames(df_data_fit[col]),range_var_min,range_var_max,range_var,range_var_less5per,range_var_plus5per,per_n_cells_in_range)
            data_var_all<-rbind(data_var_all,data_var)
            }
          colnames(data_var_all)<-c("ID_species","Species","Variable","range_var_min","range_var_max","range_var",
          "range_var_less5per","range_var_plus5per","per_n_cells_in_range")
          dev.off()  
          #Using the range of data of table_range_environmental_data we are going to convert the data ouside the range to NA
            for (j in 1:6) {
            data_pre_for_RCP60_projection[,j][data_pre_for_RCP60_projection[,j]>data_var_all[j,"range_var_plus5per"]]<-NA
            data_pre_for_RCP60_projection[,j][data_pre_for_RCP60_projection[,j]<data_var_all[j,"range_var_less5per"]]<-NA
            }
          #BIOMOD_EnsembleForecasting RCP60 copy from WILFRED corrections
          Sps_esemble_RCP60_proj<-BIOMOD_EnsembleForecasting(EM.output=Sps_essemble_models,
          new.env=data_pre_for_RCP60_projection,
          proj.name="RCP60",
          selected.models = 'all', 
          binary.meth = "TSS", 
          do.stack=FALSE,
          build.clamping.mask = F)
          #proc.time() - ptm3 
          save(Sps_esemble_RCP60_proj, file="Sps_esemble_RCP60_proj.RData")
          not_predicted_RCP60<-sum(is.na(Sps_esemble_RCP60_proj@ proj @ val))
          predicted_area_RCP60<-(predicted_area_absolute_RCP60-not_predicted_RCP60)/predicted_area_absolute_RCP60*100
          save(predicted_area_RCP60, file=paste0(path_save_sps_Biomod,"predicted_area_RCP60_sps",i,".RData"))
          table_range_environmental_data_RCP60_scenario<-data.frame(data_var_all,var_imp_df_mean,"Predict_area"=predicted_area_RCP60)
          save(table_range_environmental_data_RCP60_scenario, file=paste0(path_save_sps_Biomod,"table_range_environmental_data_RCP60_scenario_sps_",i,".RData"))
          table_range_environmental_data_RCP60_scenario_all<-rbind(table_range_environmental_data_RCP60_scenario_all,table_range_environmental_data_RCP60_scenario)   
          print(paste0("Precentage of predicted area for species ",list_sps_in_loop$Species," = ",predicted_area_RCP60))
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        }
      setwd(path_data)
      save(table_range_environmental_data_RCP60_scenario_all, file=paste0(path_save,"table_range_environmental_data_RCP60_scenario_all.RData"))

   ##4.5.3.4 Code for create a ID for the cells of raster
     rm(list=ls()) 
     path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
     path_save = paste0(path_data,"SDM_MOD/")
     load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix.RData")# Copy here to try to be faster
     sub_DF_Scen_RCP60_Matrix<-DF_Scen_RCP60_Matrix[,c(2:3)]
     #We create a small dataframe with the ID of the cells of the raster
     save(sub_DF_Scen_RCP60_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP60_Matrix.RData")

   #4.5.3.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
    rm(list=ls())
    library(biomod2) 
    library(slam)
    library(MASS)
    library(usdm)
    library(rgdal) 
    library(caret)
    library(MuMIn)
    library(raster)
    library(rgdal) 
    library(raster) 
    path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
    path_save = paste0(path_data,"SDM_MOD")
    load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA_2.RData")# Copy here to try to be faster
    load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
    levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
    df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
    df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
    DF_subpopulation<-as.data.frame(DF_Scen_RCP60_Matrix_noNA_2$vec_subpop_raster)
    colnames(DF_subpopulation)<-c("vec_subpop_raster")
    sps_in_loop<-c()
    DF_habitat_continuous_RCP60<-c()
    #Species with errors
    for (i in 1:276){#The number of species to model
      list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
      #Condition of a minimum number of presences below n=150 we consider that we have not enought data
      #Condition of a minimum number of presences below n=150 we consider that we have not enought data
      list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
      #Here we are going to load the data in what supopulations the species is eaten and 
      #we will use the subpoplations to predict only in these areas:
      sps_in_loop<-factor(list_sps_in_loop$Species)
      sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
      #We have some species with energy as zero
      try({
        row_sum<-rowSums(sub_subpopulations_diet)
      })
      #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
      #the species when there is no energy of the species in any subpopulation. This happen with some species
      #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
      list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
      list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
      #Conditional for skype species with errors in the models:
      if ((i==109|i==164)){  
        list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
      }
      if ((list_sps_in_loop$N_pixels_presence)>=50){  
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        print("             ##############################################################          ")
        print("             ########           Starting a new species             ########          ")
        print("             ##############################################################          ")
        print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
        path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
        path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
        setwd(path_save_sps_Biomod)
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
        pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
        vec_populations_with_species<-colnames(pop_diet)
        vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
        vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
        vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
        vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
        vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
        vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
        vec_populations_with_species[vec_populations_with_species=="central"] <- 7
        vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
        vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
        vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
        vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
        vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
        vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
        vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
        vec_populations_with_species<-as.numeric(vec_populations_with_species)           
        DF_subpopulation$subpopulations_species_detected_not_detected<-c(NA)
        DF_subpopulation$subpopulations_species_detected_not_detected[DF_subpopulation$vec_subpop_raster %in% vec_populations_with_species] <- 1               
        #We select the subpopulations where the species is detected: 
        DF_supop_diet_present<- subset(DF_subpopulation, subpopulations_species_detected_not_detected==1)
        DF_supop_diet_NO_present<- subset(DF_subpopulation, is.na(subpopulations_species_detected_not_detected))
        DF_supop_diet_NO_present$habitat_continuous<-c(NA)
        supopulations_present<- unique(DF_supop_diet_present$vec_subpop_raster)
        supopulations_NO_present<- unique(DF_supop_diet_NO_present$vec_subpop_raster)
        #For ALL the subpopulations where it has been dtectd the species:
        #We load the prediction of habitat for RCP60 scenario:
        #We load the essemble forecast in continuous:
        load(file="Sps_esemble_RCP60_proj.RData")
        load(file="Sps_essemble_models.RData")
        habitat_continuous<-(Sps_esemble_RCP60_proj@ proj@ val[,1])
        #Now we are going to replace "habitat_continuous" by "habitat_binary"
        DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
        path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
        path_scenario = paste0(path_habitat,"/RCP60")
        #dir.create(path_habitat)
        dir.create(path_scenario)
        #For each the subpopulations where it has been detected the species:
        for (j in 1:NROW(supopulations_present)){#The number of species to model
          subpopulation_in_loop<-supopulations_present[j]  
          print("             ########           Starting a new Subpopulation             ########          ")
          print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
          sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
          assign(paste0("sub_habitat_subpop_",subpopulation_in_loop),sub_habitat_subpop$habitat_continuous)
          save(list=paste0("sub_habitat_subpop_", subpopulation_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_in_loop, ".Rdata"))   
        }  
        #For each the subpopulations where NO has been detected the species:
        for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
          subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
          print("             ########           Starting a new Subpopulation             ########          ")
          print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
          sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
          sub_habitat_subpop<-sub_habitat_subpop_no_present
          assign(paste0("sub_habitat_subpop_",subpopulation_no_present_in_loop),sub_habitat_subpop$habitat_continuous)
          save(list=paste0("sub_habitat_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
        }  
        print(paste0("             ######## Finished at ",Sys.time(),"########          "))
        print("############################################################################################################################")
        print("############################################################################################################################")
        print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      }        
    } 

    #4.5.3.6 We merge for each subpopulations the predictions for all species
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      DF_subpopulation<-as.data.frame(DF_Scen_RCP60_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_RCP60<-c()
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sub_habitat_subpop_1_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
      sub_habitat_subpop_2_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
      sub_habitat_subpop_3_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
      sub_habitat_subpop_4_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
      sub_habitat_subpop_5_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
      sub_habitat_subpop_6_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
      sub_habitat_subpop_7_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
      sub_habitat_subpop_8_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
      sub_habitat_subpop_9_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
      sub_habitat_subpop_10_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
      sub_habitat_subpop_11_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
      sub_habitat_subpop_12_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
      sub_habitat_subpop_13_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
      sub_habitat_subpop_14_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
      sub_habitat_list_sps_loop<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        if ((list_sps_in_loop$N_pixels_presence)>=50){
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP60")
          sub_habitat_list_sps_loop<-cbind(sub_habitat_list_sps_loop,i)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:14){#The number of species to model
            load(file=paste0(path_scenario,"/sub_habitat_subpop_", j, ".Rdata"))
          }  
          sub_habitat_subpop_1_ALL<-cbind(sub_habitat_subpop_1_ALL, sub_habitat_subpop_1)
          sub_habitat_subpop_2_ALL<-cbind(sub_habitat_subpop_2_ALL, sub_habitat_subpop_2)
          sub_habitat_subpop_3_ALL<-cbind(sub_habitat_subpop_3_ALL, sub_habitat_subpop_3)
          sub_habitat_subpop_4_ALL<-cbind(sub_habitat_subpop_4_ALL, sub_habitat_subpop_4)
          sub_habitat_subpop_5_ALL<-cbind(sub_habitat_subpop_5_ALL, sub_habitat_subpop_5)
          sub_habitat_subpop_6_ALL<-cbind(sub_habitat_subpop_6_ALL, sub_habitat_subpop_6)
          sub_habitat_subpop_7_ALL<-cbind(sub_habitat_subpop_7_ALL, sub_habitat_subpop_7)
          sub_habitat_subpop_8_ALL<-cbind(sub_habitat_subpop_8_ALL, sub_habitat_subpop_8)
          sub_habitat_subpop_9_ALL<-cbind(sub_habitat_subpop_9_ALL, sub_habitat_subpop_9)
          sub_habitat_subpop_10_ALL<-cbind(sub_habitat_subpop_10_ALL, sub_habitat_subpop_10)
          sub_habitat_subpop_11_ALL<-cbind(sub_habitat_subpop_11_ALL, sub_habitat_subpop_11)
          sub_habitat_subpop_12_ALL<-cbind(sub_habitat_subpop_12_ALL, sub_habitat_subpop_12)
          sub_habitat_subpop_13_ALL<-cbind(sub_habitat_subpop_13_ALL, sub_habitat_subpop_13)
          sub_habitat_subpop_14_ALL<-cbind(sub_habitat_subpop_14_ALL, sub_habitat_subpop_14)
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 
      sub_habitat_list_sps_loop<-c("ID_pixel",sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_1_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_2_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_3_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_4_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_5_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_6_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_7_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_8_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_9_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_10_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_11_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_12_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_13_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_14_ALL)<-c(sub_habitat_list_sps_loop)
      sub_habitat_Europe_RCP60<-rbind(sub_habitat_subpop_1_ALL,
                sub_habitat_subpop_2_ALL,
                sub_habitat_subpop_3_ALL,
                sub_habitat_subpop_4_ALL,
                sub_habitat_subpop_5_ALL,
                sub_habitat_subpop_6_ALL,
                sub_habitat_subpop_7_ALL,
                sub_habitat_subpop_8_ALL,
                sub_habitat_subpop_9_ALL,
                sub_habitat_subpop_10_ALL,
                sub_habitat_subpop_11_ALL,
                sub_habitat_subpop_12_ALL,
                sub_habitat_subpop_13_ALL,
                sub_habitat_subpop_14_ALL)
      setwd(path_data)
      save(sub_habitat_Europe_RCP60, file=paste0(path_save,"sub_habitat_Europe_RCP60.RData"))

    #4.5.3.7 We divide the dataframe with all species in smaller dataframes for do the merge
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load(file=paste0(path_save,"sub_habitat_Europe_RCP60.RData"))# Copy here to try to be faster
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP60_Matrix.RData")# Copy here to try to be faster
      DF_ID_pixel<-as.data.frame(sub_DF_Scen_RCP60_Matrix$ID_pixel)
      colnames(DF_ID_pixel)<-c("ID_pixel")
      save(DF_ID_pixel, file=paste0(path_save,"DF_ID_pixel.RData"))
      sub_habitat_Europe_RCP60_1<- sub_habitat_Europe_RCP60[,1:41]
      sub_habitat_Europe_RCP60_2<- sub_habitat_Europe_RCP60[,c(1,42:82)]
      sub_habitat_Europe_RCP60_3<- sub_habitat_Europe_RCP60[,c(1,83:123)]
      sub_habitat_Europe_RCP60_4<- sub_habitat_Europe_RCP60[,c(1,124:164)]
      sub_habitat_Europe_RCP60_5<- sub_habitat_Europe_RCP60[,c(1,165:205)]
      sub_habitat_Europe_RCP60_6<- sub_habitat_Europe_RCP60[,c(1,206:237)]
      save(sub_habitat_Europe_RCP60_1, file=paste0(path_save,"sub_habitat_Europe_RCP60_1.RData"))
      save(sub_habitat_Europe_RCP60_2, file=paste0(path_save,"sub_habitat_Europe_RCP60_2.RData"))
      save(sub_habitat_Europe_RCP60_3, file=paste0(path_save,"sub_habitat_Europe_RCP60_3.RData"))
      save(sub_habitat_Europe_RCP60_4, file=paste0(path_save,"sub_habitat_Europe_RCP60_4.RData"))
      save(sub_habitat_Europe_RCP60_5, file=paste0(path_save,"sub_habitat_Europe_RCP60_5.RData"))
      save(sub_habitat_Europe_RCP60_6, file=paste0(path_save,"sub_habitat_Europe_RCP60_6.RData"))

    #4.5.3.8 We do the merge for each small data frame
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_1.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_1<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_1, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_1, file=paste0(path_save,"merged_habitat_Europe_RCP60_1.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_2.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_2<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_2, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_2, file=paste0(path_save,"merged_habitat_Europe_RCP60_2.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_3.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_3<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_3, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_3, file=paste0(path_save,"merged_habitat_Europe_RCP60_3.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_4.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_4<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_4, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_4, file=paste0(path_save,"merged_habitat_Europe_RCP60_4.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_5.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_5<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_5, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_5, file=paste0(path_save,"merged_habitat_Europe_RCP60_5.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP60_6.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP60_6<-merge(DF_ID_pixel, sub_habitat_Europe_RCP60_6, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP60_6, file=paste0(path_save,"merged_habitat_Europe_RCP60_6.RData"))
          
    #4.5.3.9 Create the raster files with the predicitons
        rm(list=ls())
          library(raster) 
          library(sp) 
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          #We load the satellite image for background maps created in part 3  
          load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
          #As these files are super heavy we will load one, then run the code with the species of the dataframe,
          #then close R and load other dataframe: 
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP60_1.RData"))# From 1 to 45 
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP60_2.RData"))# From 46 to 91
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP60_3.RData"))# From 92 to 138
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP60_4.RData"))# From 139 to 188
            load(file=paste0(path_save,"merged_habitat_Europe_RCP60_5.RData"))# From 189 to 240
            load(file=paste0(path_save,"merged_habitat_Europe_RCP60_6.RData"))#  From 242 to 276
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sps_in_loop<-c()
      X11()
      for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP60")
          if ((i)<46){
          merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_1
          }
          if (i>=46 && i<92){ 
          merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_2
          }
         if (i>=92 && i<139){ 
          merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_3
          }
         if (i>=139 && i<189){  
          merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_4
          }
         if (i>=189 && i<241){
         merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_5
          }
         if (i>=241 && i<277){  
          merged_habitat_Europe_RCP60<-merged_habitat_Europe_RCP60_6
          }
            #We save into a raster file to see in Arcgis
                new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                #We defined the projection of our raster:
                projection(new_raster) <- newproj
                values(new_raster)<-merged_habitat_Europe_RCP60[,as.character(i)]
                #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
                writeRaster(new_raster, paste0(path_scenario,"/habitat.img"), overwrite=TRUE)
           #We plot a map with the google satellite image and we save it in jpeg         
            try({
            plotRGB(pr2_crop) 
            })
            plot(new_raster, add = T,col=(gray.colors(12)), colNA=NA)  
            try({#col = gray.colors(12)
            dev.copy(pdf,file=paste0(path_scenario,"/Habitat_RCP60_Species_",i,".pdf",sep=""))
            })
            dev.off()
          #We save it in a kml file to see in google earth      
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            habitat_lat_long <- projectRaster(new_raster, new_raster_2, method='bilinear')
            KML(habitat_lat_long, file=paste0(path_scenario,"/habitat.kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)
                  print(paste0("             ######## Maps Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
                  }        
            } 

    #4.5.3.10 Code for copy the results to folders
      rm(list=ls())
      install.packages("pdftools")
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP60_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP60_prediction_pdf"
      dir.create(destination_folder_kml)
      dir.create(destination_folder_pdf)
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP60")
          file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP60_sps_",i,".img")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP60_sps",i,".img.aux.xml")) 
          file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP60_sps_",i,".kmz")) 
          setwd(destination_folder_pdf)
          pdf_convert(paste0(path_scenario,"/Habitat_RCP60_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP60_Species_",i,".jpeg"),
          dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }

    #4.5.3.10.1 Code for copy the results to folders for WILD species
      rm(list=ls())
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP60_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP60_prediction_pdf"
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)
      for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors and the brown bear in the models:
       if ((i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP60")
          file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP60_sps_",i,".img")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP60_sps",i,".img.aux.xml")) 
          file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP60_sps_",i,".kmz")) 
          setwd(destination_folder_pdf)
          pdf_convert(paste0(path_scenario,"/Habitat_RCP60_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP60_Species_",i,".jpeg"),
          dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }

  #########################################################################################################
  #4.5.4 Ensemble Forecasting of RCP85 Scenario 
  #########################################################################################################
      
    #4.5.4.1 We merge the data of the matrix to the RCP85 conditions          
      rm(list=ls()) 
      path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("H:/D/Test_Biomod")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_RCP85.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/data_matrix.RData")
      DF_Scen_RCP85_Matrix<-as.data.frame(cbind(data_matrix,sub_Scenario_RCP85))
      save(DF_Scen_RCP85_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix.RData")
      DF_Scen_RCP85_Matrix_noNA<-na.omit(DF_Scen_RCP85_Matrix)
      save(DF_Scen_RCP85_Matrix_noNA, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix_noNA.RData")
    #4.5.4.2 Code for create a reduced area for extrapolate models and for include the subpopulations asigned to the diet         
      rm(list=ls()) 
      path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      setwd("H:/D/Test_Biomod")
      #Map extrapolacion
      load("H:/D/Test_Biomod/vec_extrap_raster.RData")
      #Map of subpopulations
      load("H:/D/Test_Biomod/vec_subpop_raster.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix.RData")# Copy here to try to be faster
      DF_Scen_RCP85_Matrix_2<-as.data.frame(cbind(DF_Scen_RCP85_Matrix,vec_extrap_raster,vec_subpop_raster))
      DF_Scen_RCP85_Matrix_noNA_2<-na.omit(DF_Scen_RCP85_Matrix_2)
      save(DF_Scen_RCP85_Matrix_noNA_2, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix_noNA_2.RData")

    #4.5.4.3 Application of BIOMOD_EnsembleForecasting to RCP85 Scenario    
      rm(list=ls()) 
      library(biomod2) 
      library(slam)
      library(MASS)
      library(usdm)
      library(rgdal) 
      library(caret)
      library(MuMIn)
      library(raster)
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD/")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      merge_list_diet_energy<-merge(df_datos_1_276_description_GBIF,subpopulations_diet, by.x="Species", by.y = "species",all.x = T)
      merge_list_diet_energy2<-merge(subpopulations_diet,df_datos_1_276_description_GBIF, by.x="species", by.y = "Species",all.x = T)
      sub_Na_merge_list_diet_energy<-merge_list_diet_energy[is.na(merge_list_diet_energy$central),]
      sub_Na_merge_list_diet_energy2<-merge_list_diet_energy2[is.na(merge_list_diet_energy2$N_pixels_presence),]
      #Descripcion de los NAs:
                #Para sub_Na_merge_list_diet_energy tenemos tres especies con NA que son Fagus orientalis, Lonicera caucasica y Pahseolus vulgaris que se les 
                #ha considerado rF_original = 0 en el estudio por lo que la REDEC sale NA, aunque la especie est? descrita coo parte de la dieta
                #Para sub_Na_merge_list_diet_energy2 tenemos un caso que es el Sus scrofa domesticus que est? en el estudio de Naves y no lo tenemos invluido como 
                #especie/subes pecie o taxon, con este sumarian 277 taxones
      table_range_environmental_data_RCP85_scenario_all<-c()
      X11()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        #Condition for i=109 # With this species we have no model ! No models kept due to threshold filtering... Ensemble Modeling was skipped!
        list_sps_in_loop$N_pixels_presence[(list_sps_in_loop$i)==109]<-0
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i,"/")
          path_save_Univariable_models = paste0(path_save_sps,"/Univariable_models/")
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod/")
          setwd(path_save_sps_Biomod)
          #This is the dataset with the predictors used in the model "Predictors_RCP85" these data are unique for each species
          #We select from the dataframe with the extension to project the variables included in Predictors_RCP85"
          #load(file=paste(path_save_Univariable_models,"variables_names_less_AICc.RData"))
          load(file=paste0(path_save_sps_Biomod,"Sps_essemble_models.RData"))  
          variables_names<-Sps_essemble_models@ expl.var.names
          #Here we are going to load the data in what supopulations the species is eaten and 
          #we will use the subpoplations to predict only in these areas:
          sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
          pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
          vec_populations_with_species<-colnames(pop_diet)
          vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
          vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
          vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
          vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
          vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
          vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
          vec_populations_with_species[vec_populations_with_species=="central"] <- 7
          vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
          vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
          vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
          vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
          vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
          vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
          vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
          vec_populations_with_species<-as.numeric(vec_populations_with_species)           
          data_pre_for_RCP85_projection_0<-DF_Scen_RCP85_Matrix_noNA_2[DF_Scen_RCP85_Matrix_noNA_2$vec_subpop_raster %in% vec_populations_with_species,]
          data_pre_for_RCP85_projection<-data_pre_for_RCP85_projection_0[,variables_names]
          predicted_area_absolute_RCP85<-NROW(data_pre_for_RCP85_projection)
          var_imp_df<- get_variables_importance(Sps_essemble_models,as.data.frame=T) 
          var_imp_df_mean<-data.frame(apply(var_imp_df,1,mean))
          colnames(var_imp_df_mean)<-c("Var_importance")
          #Code for see the range of the data of variables used to fit the model and compare with the range of data of each scenario:
          #We load the data used to fit the model of the species
          load(file=paste0(path_save_sps_Biomod,"format_data.RData"))
          df_data_fit<-format_data @ data.env.var
          data_var_all<-c()
          pdf(file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Plots_environ_RCP85/Range_env_data_sps_",i,".pdf",sep=""), width = 0.394 * 64, height = 0.394 * 8)
          par(mfrow = c(1, 6), cex = 0.7)
            for (col in 1:ncol(df_data_fit)) {
            d_all_europe <- density(data_pre_for_RCP85_projection[,col])
            plot(d_all_europe, main=colnames(df_data_fit[col]))
            d <- density(df_data_fit[,col])
            lines(d, main=colnames(df_data_fit[col]))
            polygon(d, col="red", border="red")
            abline(v = max(df_data_fit[,col]), lty=2, col="blue")
            abline(v = min(df_data_fit[,col]), lty=2, col="blue")
            lines(d_all_europe, main=colnames(df_data_fit[col]))
            range_var_max<-max(df_data_fit[,col])
            range_var_min<-min(df_data_fit[,col])
            range_var<-max(df_data_fit[,col])-min(df_data_fit[,col])
            range_var_plus5per<-range_var_min + range_var*1.2
            range_var_less5per<-range_var_min-range_var*0.2
            data_in_range<-with(data_pre_for_RCP85_projection, data_pre_for_RCP85_projection[,col] <= range_var_plus5per & 
            data_pre_for_RCP85_projection[,col] >= range_var_less5per)
            #This is a condition to apply when all pixels are outside the range
            if (any(data_in_range)){  
              tab_var_europe<-table(data_in_range)
              df_tab_var_europe<-as.data.frame(tab_var_europe)
              per_n_cells_in_range<-df_tab_var_europe[which(df_tab_var_europe$data_in_range=="TRUE"),"Freq"]/6252207*100
            }  else   {
              per_n_cells_in_range<-0
            }  
            data_var<-data.frame(i,list_sps_in_loop$Species,colnames(df_data_fit[col]),range_var_min,range_var_max,range_var,range_var_less5per,range_var_plus5per,per_n_cells_in_range)
            data_var_all<-rbind(data_var_all,data_var)
            }
          colnames(data_var_all)<-c("ID_species","Species","Variable","range_var_min","range_var_max","range_var",
          "range_var_less5per","range_var_plus5per","per_n_cells_in_range")
          dev.off()  
          #Using the range of data of table_range_environmental_data we are going to convert the data ouside the range to NA
            for (j in 1:6) {
            data_pre_for_RCP85_projection[,j][data_pre_for_RCP85_projection[,j]>data_var_all[j,"range_var_plus5per"]]<-NA
            data_pre_for_RCP85_projection[,j][data_pre_for_RCP85_projection[,j]<data_var_all[j,"range_var_less5per"]]<-NA
            }
          #BIOMOD_EnsembleForecasting RCP85 copy from WILFRED corrections
          Sps_esemble_RCP85_proj<-BIOMOD_EnsembleForecasting(EM.output=Sps_essemble_models,
          new.env=data_pre_for_RCP85_projection,
          proj.name="RCP85",
          selected.models = 'all', 
          binary.meth = "TSS", 
          do.stack=FALSE,
          build.clamping.mask = F)
          #proc.time() - ptm3 
          save(Sps_esemble_RCP85_proj, file="Sps_esemble_RCP85_proj.RData")
          not_predicted_RCP85<-sum(is.na(Sps_esemble_RCP85_proj@ proj @ val))
          predicted_area_RCP85<-(predicted_area_absolute_RCP85-not_predicted_RCP85)/predicted_area_absolute_RCP85*100
          save(predicted_area_RCP85, file=paste0(path_save_sps_Biomod,"predicted_area_RCP85_sps",i,".RData"))
          table_range_environmental_data_RCP85_scenario<-data.frame(data_var_all,var_imp_df_mean,"Predict_area"=predicted_area_RCP85)
          save(table_range_environmental_data_RCP85_scenario, file=paste0(path_save_sps_Biomod,"table_range_environmental_data_RCP85_scenario_sps_",i,".RData"))
          table_range_environmental_data_RCP85_scenario_all<-rbind(table_range_environmental_data_RCP85_scenario_all,table_range_environmental_data_RCP85_scenario)   
          print(paste0("Precentage of predicted area for species ",list_sps_in_loop$Species," = ",predicted_area_RCP85))
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        }
      setwd(path_data)
      save(table_range_environmental_data_RCP85_scenario_all, file=paste0(path_save,"table_range_environmental_data_RCP85_scenario_all.RData"))

   #4.5.4.4 Code for create a ID for the cells of raster
     rm(list=ls()) 
     path_data = "H:/D/Project_Name/Results_Biomod/R_analysis/"
     path_save = paste0(path_data,"SDM_MOD/")
     load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix.RData")# Copy here to try to be faster
     sub_DF_Scen_RCP85_Matrix<-DF_Scen_RCP85_Matrix[,c(2:3)]
     #We create a small dataframe with the ID of the cells of the raster
     save(sub_DF_Scen_RCP85_Matrix, file="H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP85_Matrix.RData")

   #4.5.4.5 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
    rm(list=ls())
    library(biomod2) 
    library(slam)
    library(MASS)
    library(usdm)
    library(rgdal) 
    library(caret)
    library(MuMIn)
    library(raster)
    library(rgdal) 
    library(raster) 
    path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
    path_save = paste0(path_data,"SDM_MOD")
    load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix_noNA_2.RData")# Copy here to try to be faster
    load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
    levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
    df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
    df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
    DF_subpopulation<-as.data.frame(DF_Scen_RCP85_Matrix_noNA_2$vec_subpop_raster)
    colnames(DF_subpopulation)<-c("vec_subpop_raster")
    sps_in_loop<-c()
    DF_habitat_continuous_RCP85<-c()
    for (i in 1:276){#The number of species to model
      list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
      list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
      #Here we are going to load the data in what supopulations the species is eaten and 
      #we will use the subpoplations to predict only in these areas:
      sps_in_loop<-factor(list_sps_in_loop$Species)
      sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
      #We have some species with energy as zero
      try({
        row_sum<-rowSums(sub_subpopulations_diet)
      })
      #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
      #the species when there is no energy of the species in any subpopulation. This happen with some species
      #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
      list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
      list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
      #Conditional for skype species with errors in the models:
      if ((i==109|i==164)){  
        list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
      }
      if ((list_sps_in_loop$N_pixels_presence)>=50){  
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        print("             ##############################################################          ")
        print("             ########           Starting a new species             ########          ")
        print("             ##############################################################          ")
        print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
        path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
        path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
        setwd(path_save_sps_Biomod)
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
        pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
        vec_populations_with_species<-colnames(pop_diet)
        vec_populations_with_species[vec_populations_with_species=="alpine"] <- 1
        vec_populations_with_species[vec_populations_with_species=="baltic"] <- 2
        vec_populations_with_species[vec_populations_with_species=="cantabrian"] <- 3
        vec_populations_with_species[vec_populations_with_species=="carpath_east"] <- 4
        vec_populations_with_species[vec_populations_with_species=="carpath_west"] <- 5
        vec_populations_with_species[vec_populations_with_species=="caucasian"] <- 6
        vec_populations_with_species[vec_populations_with_species=="central"] <- 7
        vec_populations_with_species[vec_populations_with_species=="dinaric_east"] <- 8
        vec_populations_with_species[vec_populations_with_species=="dinaric_south"] <- 9
        vec_populations_with_species[vec_populations_with_species=="dinaric_west"] <- 10
        vec_populations_with_species[vec_populations_with_species=="karelian"] <- 11
        vec_populations_with_species[vec_populations_with_species=="pirynees"] <- 12
        vec_populations_with_species[vec_populations_with_species=="scandinavian"] <- 13
        vec_populations_with_species[vec_populations_with_species=="turkey"] <- 14
        vec_populations_with_species<-as.numeric(vec_populations_with_species)           
        DF_subpopulation$subpopulations_species_detected_not_detected<-c(NA)
        DF_subpopulation$subpopulations_species_detected_not_detected[DF_subpopulation$vec_subpop_raster %in% vec_populations_with_species] <- 1               
        #We select the subpopulations where the species is detected: 
        DF_supop_diet_present<- subset(DF_subpopulation, subpopulations_species_detected_not_detected==1)
        DF_supop_diet_NO_present<- subset(DF_subpopulation, is.na(subpopulations_species_detected_not_detected))
        DF_supop_diet_NO_present$habitat_continuous<-c(NA)
        supopulations_present<- unique(DF_supop_diet_present$vec_subpop_raster)
        supopulations_NO_present<- unique(DF_supop_diet_NO_present$vec_subpop_raster)
        #For ALL the subpopulations where it has been dtectd the species:
        #We load the prediction of habitat for RCP85 scenario:
        #We load the essemble forecast in continuous:
        load(file="Sps_esemble_RCP85_proj.RData")
        load(file="Sps_essemble_models.RData")
        habitat_continuous<-(Sps_esemble_RCP85_proj@ proj@ val[,1])
        DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
        path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
        path_scenario = paste0(path_habitat,"/RCP85")
        dir.create(path_scenario)
        #For each the subpopulations where it has been detected the species:
        for (j in 1:NROW(supopulations_present)){#The number of species to model
          subpopulation_in_loop<-supopulations_present[j]  
          print("             ########           Starting a new Subpopulation             ########          ")
          print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
          sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
          assign(paste0("sub_habitat_subpop_",subpopulation_in_loop),sub_habitat_subpop$habitat_continuous)
          save(list=paste0("sub_habitat_subpop_", subpopulation_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_in_loop, ".Rdata"))   
        }  
        #For each the subpopulations where NO has been detected the species:
        for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
          subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
          print("             ########           Starting a new Subpopulation             ########          ")
          print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
          sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
          sub_habitat_subpop<-sub_habitat_subpop_no_present
          assign(paste0("sub_habitat_subpop_",subpopulation_no_present_in_loop),sub_habitat_subpop$habitat_continuous)
          save(list=paste0("sub_habitat_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario,"/sub_habitat_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
        }  
        print(paste0("             ######## Finished at ",Sys.time(),"########          "))
        print("############################################################################################################################")
        print("############################################################################################################################")
        print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      }        
    } 
                  
    #4.5.4.6 We merge for each subpopulations the predictions for all species
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP85_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      DF_subpopulation<-as.data.frame(DF_Scen_RCP85_Matrix_noNA_2$vec_subpop_raster)
      colnames(DF_subpopulation)<-c("vec_subpop_raster")
      sps_in_loop<-c()
      DF_habitat_continuous_RCP85<-c()
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sub_habitat_subpop_1_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
      sub_habitat_subpop_2_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
      sub_habitat_subpop_3_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
      sub_habitat_subpop_4_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
      sub_habitat_subpop_5_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
      sub_habitat_subpop_6_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
      sub_habitat_subpop_7_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
      sub_habitat_subpop_8_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
      sub_habitat_subpop_9_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
      sub_habitat_subpop_10_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
      sub_habitat_subpop_11_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
      sub_habitat_subpop_12_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
      sub_habitat_subpop_13_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
      sub_habitat_subpop_14_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
      sub_habitat_list_sps_loop<-c()
      for (i in 1:276){#The number of species to model
        list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
        list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
        #Here we are going to load the data in what supopulations the species is eaten and 
        #we will use the subpoplations to predict only in these areas:
        sps_in_loop<-factor(list_sps_in_loop$Species)
        sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
        #We have some species with energy as zero
        try({
          row_sum<-rowSums(sub_subpopulations_diet)
        })
        #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
        #the species when there is no energy of the species in any subpopulation. This happen with some species
        #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
        list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
        list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
        #Conditional for skype species with errors in the models:
        if ((i==109|i==164|i==19)){  
          list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
        if ((list_sps_in_loop$N_pixels_presence)>=50){
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP85")
          sub_habitat_list_sps_loop<-cbind(sub_habitat_list_sps_loop,i)
          #For each the subpopulations where it has been detected the species:
          for (j in 1:14){#The number of species to model
            load(file=paste0(path_scenario,"/sub_habitat_subpop_", j, ".Rdata"))
          }  
          sub_habitat_subpop_1_ALL<-cbind(sub_habitat_subpop_1_ALL, sub_habitat_subpop_1)
          sub_habitat_subpop_2_ALL<-cbind(sub_habitat_subpop_2_ALL, sub_habitat_subpop_2)
          sub_habitat_subpop_3_ALL<-cbind(sub_habitat_subpop_3_ALL, sub_habitat_subpop_3)
          sub_habitat_subpop_4_ALL<-cbind(sub_habitat_subpop_4_ALL, sub_habitat_subpop_4)
          sub_habitat_subpop_5_ALL<-cbind(sub_habitat_subpop_5_ALL, sub_habitat_subpop_5)
          sub_habitat_subpop_6_ALL<-cbind(sub_habitat_subpop_6_ALL, sub_habitat_subpop_6)
          sub_habitat_subpop_7_ALL<-cbind(sub_habitat_subpop_7_ALL, sub_habitat_subpop_7)
          sub_habitat_subpop_8_ALL<-cbind(sub_habitat_subpop_8_ALL, sub_habitat_subpop_8)
          sub_habitat_subpop_9_ALL<-cbind(sub_habitat_subpop_9_ALL, sub_habitat_subpop_9)
          sub_habitat_subpop_10_ALL<-cbind(sub_habitat_subpop_10_ALL, sub_habitat_subpop_10)
          sub_habitat_subpop_11_ALL<-cbind(sub_habitat_subpop_11_ALL, sub_habitat_subpop_11)
          sub_habitat_subpop_12_ALL<-cbind(sub_habitat_subpop_12_ALL, sub_habitat_subpop_12)
          sub_habitat_subpop_13_ALL<-cbind(sub_habitat_subpop_13_ALL, sub_habitat_subpop_13)
          sub_habitat_subpop_14_ALL<-cbind(sub_habitat_subpop_14_ALL, sub_habitat_subpop_14)
          print(paste0("             ######## Finished at ",Sys.time(),"########          "))
          print("############################################################################################################################")
          print("############################################################################################################################")
          print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        }        
      } 
      sub_habitat_list_sps_loop<-c("ID_pixel",sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_1_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_2_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_3_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_4_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_5_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_6_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_7_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_8_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_9_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_10_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_11_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_12_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_13_ALL)<-c(sub_habitat_list_sps_loop)
      colnames(sub_habitat_subpop_14_ALL)<-c(sub_habitat_list_sps_loop)
      sub_habitat_Europe_RCP85<-rbind(sub_habitat_subpop_1_ALL,
                sub_habitat_subpop_2_ALL,
                sub_habitat_subpop_3_ALL,
                sub_habitat_subpop_4_ALL,
                sub_habitat_subpop_5_ALL,
                sub_habitat_subpop_6_ALL,
                sub_habitat_subpop_7_ALL,
                sub_habitat_subpop_8_ALL,
                sub_habitat_subpop_9_ALL,
                sub_habitat_subpop_10_ALL,
                sub_habitat_subpop_11_ALL,
                sub_habitat_subpop_12_ALL,
                sub_habitat_subpop_13_ALL,
                sub_habitat_subpop_14_ALL)
      setwd(path_data)
      save(sub_habitat_Europe_RCP85, file=paste0(path_save,"sub_habitat_Europe_RCP85.RData"))

     #4.5.4.7 We divide the dataframe with all species in smaller dataframes for do the merge
      rm(list=ls())
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load(file=paste0(path_save,"sub_habitat_Europe_RCP85.RData"))# Copy here to try to be faster
      load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP85_Matrix.RData")# Copy here to try to be faster
      DF_ID_pixel<-as.data.frame(sub_DF_Scen_RCP85_Matrix$ID_pixel)
      colnames(DF_ID_pixel)<-c("ID_pixel")
      save(DF_ID_pixel, file=paste0(path_save,"DF_ID_pixel.RData"))
      sub_habitat_Europe_RCP85_1<- sub_habitat_Europe_RCP85[,1:41]
      sub_habitat_Europe_RCP85_2<- sub_habitat_Europe_RCP85[,c(1,42:82)]
      sub_habitat_Europe_RCP85_3<- sub_habitat_Europe_RCP85[,c(1,83:123)]
      sub_habitat_Europe_RCP85_4<- sub_habitat_Europe_RCP85[,c(1,124:164)]
      sub_habitat_Europe_RCP85_5<- sub_habitat_Europe_RCP85[,c(1,165:205)]
      sub_habitat_Europe_RCP85_6<- sub_habitat_Europe_RCP85[,c(1,206:236)]
      save(sub_habitat_Europe_RCP85_1, file=paste0(path_save,"sub_habitat_Europe_RCP85_1.RData"))
      save(sub_habitat_Europe_RCP85_2, file=paste0(path_save,"sub_habitat_Europe_RCP85_2.RData"))
      save(sub_habitat_Europe_RCP85_3, file=paste0(path_save,"sub_habitat_Europe_RCP85_3.RData"))
      save(sub_habitat_Europe_RCP85_4, file=paste0(path_save,"sub_habitat_Europe_RCP85_4.RData"))
      save(sub_habitat_Europe_RCP85_5, file=paste0(path_save,"sub_habitat_Europe_RCP85_5.RData"))
      save(sub_habitat_Europe_RCP85_6, file=paste0(path_save,"sub_habitat_Europe_RCP85_6.RData"))

      #4.5.4.8 We do the merge for each small data frame
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_1.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_1<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_1, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_1, file=paste0(path_save,"merged_habitat_Europe_RCP85_1.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_2.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_2<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_2, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_2, file=paste0(path_save,"merged_habitat_Europe_RCP85_2.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_3.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_3<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_3, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_3, file=paste0(path_save,"merged_habitat_Europe_RCP85_3.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_4.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_4<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_4, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_4, file=paste0(path_save,"merged_habitat_Europe_RCP85_4.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_5.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_5<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_5, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_5, file=paste0(path_save,"merged_habitat_Europe_RCP85_5.RData"))
          #WE do the merge
          rm(list=ls())
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          load(file=paste0(path_save,"sub_habitat_Europe_RCP85_6.RData"))# Copy here to try to be faster
          load(file=paste0(path_save,"DF_ID_pixel.RData"))
          merged_habitat_Europe_RCP85_6<-merge(DF_ID_pixel, sub_habitat_Europe_RCP85_6, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_habitat_Europe_RCP85_6, file=paste0(path_save,"merged_habitat_Europe_RCP85_6.RData"))

    #4.5.4.9 Create the raster files with the predicitons
          rm(list=ls())
          library(raster) 
          library(sp) 
          path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
          path_save = paste0(path_data,"SDM_MOD")
          setwd(path_data)
          #We load the satellite image for background maps created in part 3  
          load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
          #As these files are super heavy we will load one, then run the code with the species of the dataframe,
          #then close R and load other dataframe: 
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP85_1.RData"))# From 1 to 45 /46
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP85_2.RData"))# From 47 to 92
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP85_3.RData"))# From 92 to 138
            #load(file=paste0(path_save,"merged_habitat_Europe_RCP85_4.RData"))# From 139 to 188
            load(file=paste0(path_save,"merged_habitat_Europe_RCP85_5.RData"))# From 189 to 240
            load(file=paste0(path_save,"merged_habitat_Europe_RCP85_6.RData"))#  From 242 to 276
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      sps_in_loop<-c()
      X11()
      for (i in 242:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          setwd(path_save_sps_Biomod)
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP85")
          if ((i)<47){
          merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_1
          }
          if (i>=47 && i<93){ 
          merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_2
          }
         if (i>=93 && i<140){ 
          merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_3
          }
         if (i>=140 && i<190){  
          merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_4
          }
         if (i>=190 && i<243){
         merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_5
          }
         if (i>=243 && i<277){  
          merged_habitat_Europe_RCP85<-merged_habitat_Europe_RCP85_6
          }
            #We save into a raster file to see in Arcgis
                new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                #We defined the projection of our raster:
                projection(new_raster) <- newproj
                values(new_raster)<-merged_habitat_Europe_RCP85[,as.character(i)]
                #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
                writeRaster(new_raster, paste0(path_scenario,"/habitat.img"), overwrite=TRUE)
           #We plot a map with the google satellite image and we save it in jpeg         
            try({
            plotRGB(pr2_crop) 
            })
            plot(new_raster, add = T,col=(gray.colors(12)), colNA=NA)  
            try({#col = gray.colors(12)
            dev.copy(pdf,file=paste0(path_scenario,"/Habitat_RCP85_Species_",i,".pdf",sep=""))
            })
            dev.off()
          #We save it in a kml file to see in google earth      
            #RWe crate a new rster i lat lonH:/G 
            new_raster_2 <- raster(xmn=-20, xmx=105, ymn=15, ymx=75, ncols=5900, nrows=5600)
            crs(new_raster_2) <- CRS('+init=EPSG:4326')
            habitat_lat_long <- projectRaster(new_raster, new_raster_2, method='bilinear')
            KML(habitat_lat_long, file=paste0(path_scenario,"/habitat.kml"),col=(gray.colors(12)), colNA=NA, maxpixels=33050000, overwrite=TRUE)
                  print(paste0("             ######## Maps Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
                  }        
            } 

    #4.5.4.10 Code for copy the results to folders
      rm(list=ls())
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP85_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/RCP85_prediction_pdf"
      dir.create(destination_folder_kml)
      dir.create(destination_folder_pdf)
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors in the models:
       if ((i==109|i==164)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP85")
          file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
          file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP85_sps_",i,".img")) 
          file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP85_sps",i,".img.aux.xml")) 
          file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP85_sps_",i,".kmz")) 
          setwd(destination_folder_pdf)
          pdf_convert(paste0(path_scenario,"/Habitat_RCP85_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP85_Species_",i,".jpeg"),
          dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }

    #4.5.4.10.1 Code for copy the results to folders for WILD species
      rm(list=ls())
      library(pdftools)
      destination_folder_kml = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP85_prediction_kml"
      destination_folder_pdf = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP85_prediction_pdf"
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
      setwd(path_data)
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors and the brown bear in the models:
       if ((i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }
          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          path_scenario = paste0(path_habitat,"/RCP85")
           file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
           file.copy(paste0(path_scenario,"/habitat.kmz"), paste0(destination_folder_kml,"/")) 
             file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_RCP85_sps_",i,".img")) 
             file.rename(paste0(destination_folder_kml,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/habitat_RCP85_sps",i,".img.aux.xml")) 
             file.rename(paste0(destination_folder_kml,"/habitat.kmz"), paste0(destination_folder_kml,"/habitat_RCP85_sps_",i,".kmz")) 
        setwd(destination_folder_pdf)
        pdf_convert(paste0(path_scenario,"/Habitat_RCP85_Species_",i,".pdf"), format = "jpeg", pages = NULL, filenames = paste0("Habitat_RCP85_Species_",i,".jpeg"),
        dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE)        
          }
       }

     
        
####################################################################################################################################################################################
#4.6 Calculation of change in habitat from current to future RCPs scenarios
####################################################################################################################################################################################
  # We are going to calculate for each species the delta of change in habitat suitability comparing with the current scenario

    #4.6.1 Code for calculate and save the delta for each species
      rm(list=ls())
         library(dplyr)
          library(biomod2) 
          library(slam)
          library(MASS)
          library(usdm)
          library(rgdal) 
          library(caret)
          library(MuMIn)
          library(raster)
          library(rgdal) 
          library(raster) 

      library(pdftools)
      destination_folder_kml_current = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_kml"
      destination_folder_kml_RCP26 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP26_prediction_kml"
      destination_folder_kml_RCP60 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP60_prediction_kml"
      destination_folder_kml_RCP85 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP85_prediction_kml"
      
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
          
      setwd(path_data)
      
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
      load("H:/D/Test_Biomod/vec_subpop_raster.RData")      
      
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)

      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)

       for (i in 20:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors and the brown bear in the models:
       if ((i==19|i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }

          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
          #path_save_sps = paste0(path_save,"/sps_",list_sps_in_loop$i)
          #path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
          #setwd(path_save_sps_Biomod)
          #path_habitat = paste0(path_save_sps_Biomod,"/Habitat_predicition")
          #path_scenario = paste0(path_habitat,"/Current")

          # file.copy(paste0(path_scenario,"/habitat.img"), paste0(destination_folder_kml,"/")) 
          # file.copy(paste0(path_scenario,"/habitat.img.aux.xml"), paste0(destination_folder_kml,"/")) 
          # file.copy(paste0(path_scenario,"/habitat.kml"), paste0(destination_folder_kml,"/")) 
          # file.copy(paste0(path_scenario,"/habitat.png"), paste0(destination_folder_kml,"/")) 
          # file.rename(paste0(destination_folder_kml,"/habitat.img"), paste0(destination_folder_kml,"/habitat_current_sps_",i,".img"))
             
   
         habitat_current<-raster(paste0(destination_folder_kml_current,"/habitat_current_sps_",i,".img"))    
         habitat_RCP26<-raster(paste0(destination_folder_kml_RCP26,"/habitat_RCP26_sps_",i,".img"))    
         habitat_RCP60<-raster(paste0(destination_folder_kml_RCP60,"/habitat_RCP60_sps_",i,".img"))    
         habitat_RCP85<-raster(paste0(destination_folder_kml_RCP85,"/habitat_RCP85_sps_",i,".img"))    
         
        #Delta_habitat_RCP26_current<-habitat_RCP26-habitat_current
        #writeRaster(Delta_habitat_RCP26_current, "Delta_habitat_RCP26_current.img", overwrite=TRUE)
      
        #Delta_habitat_RCP60_current<-habitat_RCP60-habitat_current
        #writeRaster(Delta_habitat_RCP60_current, "Delta_habitat_RCP60_current.img", overwrite=TRUE)
      
        #Delta_habitat_RCP85_current<-habitat_RCP85-habitat_current
        #writeRaster(Delta_habitat_RCP85_current, "Delta_habitat_RCP85_current.img", overwrite=TRUE)

        #We extract the values of the raster to a vector:
        vec_habitat_current<-values(habitat_current)
        vec_habitat_RCP26<-values(habitat_RCP26)
        vec_habitat_RCP60<-values(habitat_RCP60)
        vec_habitat_RCP85<-values(habitat_RCP85)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_habitat_current,vec_habitat_RCP26,vec_habitat_RCP60,vec_habitat_RCP85,vec_subpop_raster))
        
       #We calculate the sums at subpopulation scale 
          sum_habitat_current_by_subpopulation<-tapply(data_sps_in_loop$vec_habitat_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_habitat_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_habitat_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_habitat_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_habitat_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_habitat_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_habitat_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_habitat_RCP26_current_by_subpopulation<- (sum_habitat_RCP26_by_subpopulation - sum_habitat_current_by_subpopulation)/sum_habitat_current_by_subpopulation*100
          Delta_habitat_RCP60_current_by_subpopulation<- (sum_habitat_RCP60_by_subpopulation - sum_habitat_current_by_subpopulation)/sum_habitat_current_by_subpopulation*100
          Delta_habitat_RCP85_current_by_subpopulation<- (sum_habitat_RCP85_by_subpopulation - sum_habitat_current_by_subpopulation)/sum_habitat_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_habitat_current_Europe<-sum(vec_habitat_current,na.rm=T)
          sum_habitat_RCP26_Europe<-sum(vec_habitat_RCP26,na.rm=T)
          sum_habitat_RCP60_Europe<-sum(vec_habitat_RCP60,na.rm=T)
          sum_habitat_RCP85_Europe<-sum(vec_habitat_RCP85,na.rm=T)
  
          Delta_habitat_RCP26_current_Europe <- (sum_habitat_RCP26_Europe - sum_habitat_current_Europe)/sum_habitat_current_Europe*100
          Delta_habitat_RCP60_current_Europe <- (sum_habitat_RCP60_Europe - sum_habitat_current_Europe)/sum_habitat_current_Europe*100
          Delta_habitat_RCP85_current_Europe <- (sum_habitat_RCP85_Europe - sum_habitat_current_Europe)/sum_habitat_current_Europe*100
          Delta_habitat_RCPs_current_Europe<-c(Delta_habitat_RCP26_current_Europe,Delta_habitat_RCP60_current_Europe,Delta_habitat_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_habitat_RCPs_current<-as.data.frame(rbind(Delta_habitat_RCP26_current_by_subpopulation,Delta_habitat_RCP60_current_by_subpopulation,Delta_habitat_RCP85_current_by_subpopulation))
          Delta_habitat_RCPs_current2<-as.data.frame(cbind(i,Delta_habitat_RCPs_current,Delta_habitat_RCPs_current_Europe))        
    
        #We save the data of delta for each species
          save(Delta_habitat_RCPs_current2, file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Delta_habitat/Delta_habitat_RCPs_current2_", i,".RData", sep=""))     
          }
       }
  
    #4.6.2 Code for load the information of the delta for each species
      rm(list=ls())
         library(dplyr)
          library(biomod2) 
          library(slam)
          library(MASS)
          library(usdm)
          library(rgdal) 
          library(caret)
          library(MuMIn)
          library(raster)
          library(rgdal) 
          library(raster) 

      #install.packages("pdftools")
      library(pdftools)
      destination_folder_kml_current = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Current_prediction_kml"
      destination_folder_kml_RCP26 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP26_prediction_kml"
      destination_folder_kml_RCP60 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP60_prediction_kml"
      destination_folder_kml_RCP85 = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/RCP85_prediction_kml"
      
      path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
      path_save = paste0(path_data,"SDM_MOD")
          
      setwd(path_data)
      
      load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      #We load the subpopulations used for diet and that we will use ofr subset data and run GLMM
      load("H:/D/Test_Biomod/vec_subpop_raster.RData")      
      
      levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
      df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
      df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)

      # We load the information of wild/human species  
      load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
      #We change the value of Human/Wild for the reindeer
      merged_sps_list_diet_category3[214,]
      merged_sps_list_diet_category3[214,4] <- 0
      merged_sps_list_diet_category3[214,]
      #We change the value of Human/Wild for Ficus carica
      merged_sps_list_diet_category3[276,4] <- 1
      
      merged_df_datos_1_276_description_GBIF_categories_diet<-merge(df_datos_1_276_description_GBIF, merged_sps_list_diet_category3[-1], by.x="i", by.y = "i",all.x = T)

      Delta_habitat_RCPs_current2_all<-c()
       for (i in 1:276){#The number of species to model 276
       list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       #Condition of a minimum number of presences below n=150 we consider that we have not enought data
       list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
       #Here we are going to load the data in what supopulations the species is eaten and 
       #we will use the subpoplations to predict only in these areas:
       sps_in_loop<-factor(list_sps_in_loop$Species)
       sub_subpopulations_diet <- subset(subpopulations_diet, species==as.character(sps_in_loop) ,select=central:caucasian)
       #We have some species with energy as zero
       try({
       row_sum<-rowSums(sub_subpopulations_diet)
       })
       #We "modify the value of "N_pixels_presence" in the species in loop for do not calculate
       #the species when there is no energy of the species in any subpopulation. This happen with some species
       #like Martes foina which is reported in frequency but with Volume 0 in the data of CIUC study
       list_sps_in_loop$N_pixels_presence[(nrow(sub_subpopulations_diet)<1)]<-0#with this code we skipe this species with cero energy
       list_sps_in_loop$N_pixels_presence[(row_sum==0)]<-0#with this code we skipe this species with cero energy
       #Conditional for skype species with errors and the brown bear in the models:
       if ((i==19|i==109|i==164|i==257)){  
       list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
        }

          if ((list_sps_in_loop$N_pixels_presence)>=50 & (list_sps_in_loop$Human_origin)==0){  
          list_sps_in_loop<-merged_df_datos_1_276_description_GBIF_categories_diet[c(i),]
          print("             ##############################################################          ")
          print("             ########           Starting a new species             ########          ")
          print("             ##############################################################          ")
          print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
      
          load(file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Wild/Delta_habitat/Delta_habitat_RCPs_current2_", i,".RData", sep=""))  
          Delta_habitat_RCPs_current2_all<-rbind(Delta_habitat_RCPs_current2_all,Delta_habitat_RCPs_current2)
          }
       }
    rownames(Delta_habitat_RCPs_current2_all)<-c()      
    Delta_habitat_RCPs_current2_all_2<-merge(Delta_habitat_RCPs_current2_all, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)
    Delta_vector<-rep(c("RCP26","RCP60","RCP85"),times=203)      
    Delta_habitat_RCPs_current2_all_2$Delta_RCP<-Delta_vector
    
    Delta_habitat_RCPs_current2_all_3<-(Delta_habitat_RCPs_current2_all_2[,c(1,17,20,21,16,2:15)])
    mapa_subpopulations <- shapefile('H:/G/Project_Name/Maps/Extrapolation limits/Subpopulations3.shp')
    mapa_subpopulations_df<-as.data.frame(mapa_subpopulations)
    mapa_subpopulations_df_names<-mapa_subpopulations_df$code
    mapa_subpopulations_df_names2<-mapa_subpopulations_df_names[c(1:5,13,6:12,14)]
    colnames(Delta_habitat_RCPs_current2_all_3)<-c("i","Species","Diet_category2","Delta_RCP","Delta_habitat_RCPs_current_Europe",mapa_subpopulations_df_names2)
    write.csv(Delta_habitat_RCPs_current2_all_3,file="Delta_habitat_RCPs_current2_all_3.csv", row.names = F)# Supplementary Table 37
    print(paste0("OUTPUT: Supplemenry Table 37. Variable importance for all Wild species #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("#Supplementary Table 37. Habitat suitability change for wild food species. For each species we show a numeric identificator (i), 
    their scientific name (Species), the most frequents diet categories in all subpopulations (Diet_category2), the change in habitat suitability 
    in all subpopulations where the species is consumed by brown bear (Delta_habitat_RCPs_current_Europe), and the change in habitat suitability for
    each of the subpopulations where the species is consumed by brown bear (columns 6-19). For each species there are three rows, one for each of the
    three future shared socioeconomic parthways (SSPs) scenarios for land use and climate, namely SSP1-2.6, SSP3-6.0 and SSP5-8.5.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

   #Mean delta in habitat for each RCPs for all species in Europe  
   Delta_habitat_mean_all<-tapply(Delta_habitat_RCPs_current2_all_3$Delta_habitat_RCPs_current_Europe,Delta_habitat_RCPs_current2_all_3$Delta_RCP,mean, na.rm=TRUE)
   write.csv(Delta_habitat_mean_all,file="Delta_habitat_mean_all.csv", row.names = F)# Supplementary Table 37
    print(paste0("OUTPUT: Supplemenry Table 37 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("#Supplementary Table 37. Habitat suitability change for wild food species. For each species we show a numeric identificator (i), 
    their scientific name (Species), the most frequents diet categories in all subpopulations (Diet_category2), the change in habitat suitability 
    in all subpopulations where the species is consumed by brown bear (Delta_habitat_RCPs_current_Europe), and the change in habitat suitability for
    each of the subpopulations where the species is consumed by brown bear (columns 6-19). For each species there are three rows, one for each of the
    three future shared socioeconomic parthways (SSPs) scenarios for land use and climate, namely SSP1-2.6, SSP3-6.0 and SSP5-8.5.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
   path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
   setwd(path_data)
   Delta_habitat_RCPs_current2_all_3<-read.csv("Delta_habitat_RCPs_current2_all_3.csv")     
   #Mean delta in habitat  for each RCPs and diet category for all species in Europe 
   Delta_habitat_RCPs_current2_all_3$Group_Delta_RCP_Diet_category2<-paste0(Delta_habitat_RCPs_current2_all_3$Delta_RCP,Delta_habitat_RCPs_current2_all_3$Diet_category2)
   Delta_habitat_mean_all_diet_category<-tapply(Delta_habitat_RCPs_current2_all_3$Delta_habitat_RCPs_current_Europe,Delta_habitat_RCPs_current2_all_3$Group_Delta_RCP_Diet_category2,mean, na.rm=TRUE)
   write.csv(Delta_habitat_mean_all_diet_category,file="Delta_habitat_mean_all_diet_category.csv", row.names = T)#Supplementary Table 38 # 
    print(paste0("OUTPUT: Supplemenry Table 38 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 38. Mean habitat suitability change (in percentage) for each food category and scenario.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

   #Mean dekta in habitat for each RCPs for all species by subpoplation 
   Delta_habitat_mean_all_by_subpopulation<-aggregate(Delta_habitat_RCPs_current2_all_3[,6:19],by=list(Delta_habitat_RCPs_current2_all_3$Delta_RCP), mean, na.rm=TRUE)
   write.csv(Delta_habitat_mean_all_by_subpopulation,file="Delta_habitat_mean_all_by_subpopulation.csv", row.names = T)#Supplementary Table 39# 
    print(paste0("OUTPUT: Supplemenry Table 39 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("#Supplementary Table 39. Supplementary Table 39. Mean habitat suitability change (in percentage) for all species by subpopulation and scenario.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

