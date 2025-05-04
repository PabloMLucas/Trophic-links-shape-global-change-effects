

############################################################################################################
#Readme:
############################################################################################################
#R code for calculate biotic variables for each species
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
  #/3_Download_occurrences_from_GBIF/df_datos_1_276_description_GBIF.RData
  #/4_SDM_Food_species/DF_Scen_current_Matrix_noNA_2.RData 
  #/2_Representative_Food-web_and_Food_species_list/subpopulations_diet.RData
  #/eval_DF_all3.xlsx
  #/3_Download_occurrences_from_GBIF/df_datos_1_276_description_GBIF.RData
  #/4_SDM_Food_species/Ensemble Forecasting of SDMs for Food species for Current and future Scenarios
#Data output:
    #/SDM_MODsub_energy_Europe_current_option_13_1_1.RData
    #/SDM_MODsub_energy_Europe_RCP26_option_13_1_1.RData")
    #/SDM_MODsub_energy_Europe_RCP60_option_13_1_1.RData")
    #/SDM_MODsub_energy_Europe_RCP85_option_13_1_1.RData")
    #/SDM_MODsub_energy_Europe_current_option_13_1_1_homogeneous.RData"

##############################################################################################                
#Schema
############################################################################################## 
#5_Calculation_of_Biotic_variables_for_each_sps
  #5.1 Calculation of energy spatial variables for each species
    #5.1.1 Calculation of energy variables from all species for Current Scenario 
      #5.1.1.1 Taking into account all species (wild and domestic)    
        #5.1.1.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
        #5.1.1.1.2 We merge for each subpopulations the values of spatial energy for all species
    #5.1.2 Calculation of energy variables from all species for RCP26 Scenario 
      #5.1.2.1 Taking into account all species (wild and domestic)    
        #5.1.2.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
        #5.1.2.1.2 We merge for each subpopulations the values of spatial energy for all species
    #5.1.3 Calculation of energy variables from all species for RCP60 Scenario 
      #5.1.3.1 Taking into account all species (wild and domestic)    
        #5.1.3.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
        #5.1.3.1.2 We merge for each subpopulations the values of spatial energy for all species
    #5.1.4 Calculation of energy variables from all species for RCP85 Scenario 
      #5.1.4.1 Taking into account all species (wild and domestic)    
        #5.1.4.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
        #5.1.4.1.2 We merge for each subpopulations the values of spatial energy for all species
  #5.2 Calculation of biotic interactions without account the variation in rEDEC among species
    #5.2.1 Calculation of energy variables from all species for Current Scenario 
      #5.2.1.1 Taking into account all species (wild and domestic)    
        #5.2.1.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
        #5.2.1.1.2 We merge for each subpopulations the values of spatial energy for all species
      

####################################################################################################################################################################################
####################################################################################################################################################################################
#5.1 Calculation of energy spatial variables for each species
####################################################################################################################################################################################
####################################################################################################################################################################################

  ####################################################################################################################################################################################
  #5.1.1 Calculation of energy variables from all species for Current Scenario 
  ####################################################################################################################################################################################

    #5.1.1.1 Taking into account all species (wild and domestic)    
      #5.1.1.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
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
      setwd(path_data)

      load("F:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
      load("F:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
      load("F:/G/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
      load("F:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_Scenario_current.RData")# Copy here to try to be faster

        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
        colnames(DF_subpopulation)<-c("vec_subpop_raster")
        #We create a Data frame numeric for species energy matrix
        subpopulations_diet_sub_pop_numeric<-subpopulations_diet
        colnames(subpopulations_diet_sub_pop_numeric)<-c("species","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
        sps_in_loop<-c()
        DF_habitat_continuous_current<-c()
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
          if ((i==40|i==44|i==95)){  #40,44, 95            
            list_sps_in_loop$N_pixels_presence<-0#with this code we skipe this species with cero energy
          }
          if ((list_sps_in_loop$N_pixels_presence)>=50){  
            list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Starting Species ",list_sps_in_loop$Species,"  ###  i = ",list_sps_in_loop$i))
            path_save_sps = paste0(path_save_new,"/sps_",list_sps_in_loop$i)
            path_save_sps_Biomod = paste0(path_save_sps,"/Biomod")
            setwd(path_save_sps_Biomod)
            #Here we are going to load the data in what supopulations the species is eaten and 
            #we will use the subpoplations to predict only in these areas:
            sub_subpopulations_diet[sub_subpopulations_diet == 0] <- NA
            pop_diet<-sub_subpopulations_diet[colSums(!is.na(sub_subpopulations_diet)) > 0]
            vec_populations_with_species<-colnames(pop_diet)
            #we are going to change the names of the populations to numbers
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
            habitat_continuous<-(Sps_esemble_curr_proj@ proj@ val[,1])
            DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/Current")
            energy_of_species_by_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop))
            energy_of_species_by_subpop<-energy_of_species_by_subpop[,2:15]
            dir.create(path_energy_sps)
            dir.create(path_scenario_energy_sps)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:NROW(supopulations_present)){#The number of species to model
              subpopulation_in_loop<-supopulations_present[j]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
              sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
              energy_of_species_by_subpop<-energy_of_species_by_subpop[,c(4,11,2,5,13,14,1,6,10,7,12,9,3,8)]
              colnames(energy_of_species_by_subpop)<-c(1:14)
              sub_energy_subpop_diet <- subset(energy_of_species_by_subpop, select=subpopulation_in_loop)
              vec_sub_energy_subpop_diet<-as.vector(sub_energy_subpop_diet[1,1])
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous*vec_sub_energy_subpop_diet*100
              assign(paste0("sub_energy_subpop_",subpopulation_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_in_loop, ".Rdata"))   
            }  
            #For each the subpopulations where NO has been detected the species:
            for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
              subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
              sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
              sub_habitat_subpop<-sub_habitat_subpop_no_present
              #energy_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop) ,select=subpopulation_no_present_in_loop)
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous
              #sub_energy_subpop[is.na(sub_energy_subpop)] <- 0
              assign(paste0("sub_energy_subpop_",subpopulation_no_present_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
            }  
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 

      #5.1.1.1.2 We merge for each subpopulations the values of spatial energy for all species
        rm(list=ls())
        path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD")
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
        colnames(DF_subpopulation)<-c("vec_subpop_raster")
        sps_in_loop<-c()
        DF_habitat_continuous_current<-c()
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        sub_energy_subpop_1_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
        sub_energy_subpop_2_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
        sub_energy_subpop_3_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
        sub_energy_subpop_4_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
        sub_energy_subpop_5_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
        sub_energy_subpop_6_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
        sub_energy_subpop_7_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
        sub_energy_subpop_8_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
        sub_energy_subpop_9_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
        sub_energy_subpop_10_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
        sub_energy_subpop_11_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
        sub_energy_subpop_12_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
        sub_energy_subpop_13_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
        sub_energy_subpop_14_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
        sub_energy_list_sps_loop<-c()
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
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/Current")
            sub_energy_list_sps_loop<-cbind(sub_energy_list_sps_loop,i)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:14){#The number of species to model
              load(file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", j, ".Rdata"))
            }  
            sub_energy_subpop_1_ALL<-cbind(sub_energy_subpop_1_ALL, sub_energy_subpop_1)
            sub_energy_subpop_2_ALL<-cbind(sub_energy_subpop_2_ALL, sub_energy_subpop_2)
            sub_energy_subpop_3_ALL<-cbind(sub_energy_subpop_3_ALL, sub_energy_subpop_3)
            sub_energy_subpop_4_ALL<-cbind(sub_energy_subpop_4_ALL, sub_energy_subpop_4)
            sub_energy_subpop_5_ALL<-cbind(sub_energy_subpop_5_ALL, sub_energy_subpop_5)
            sub_energy_subpop_6_ALL<-cbind(sub_energy_subpop_6_ALL, sub_energy_subpop_6)
            sub_energy_subpop_7_ALL<-cbind(sub_energy_subpop_7_ALL, sub_energy_subpop_7)
            sub_energy_subpop_8_ALL<-cbind(sub_energy_subpop_8_ALL, sub_energy_subpop_8)
            sub_energy_subpop_9_ALL<-cbind(sub_energy_subpop_9_ALL, sub_energy_subpop_9)
            sub_energy_subpop_10_ALL<-cbind(sub_energy_subpop_10_ALL, sub_energy_subpop_10)
            sub_energy_subpop_11_ALL<-cbind(sub_energy_subpop_11_ALL, sub_energy_subpop_11)
            sub_energy_subpop_12_ALL<-cbind(sub_energy_subpop_12_ALL, sub_energy_subpop_12)
            sub_energy_subpop_13_ALL<-cbind(sub_energy_subpop_13_ALL, sub_energy_subpop_13)
            sub_energy_subpop_14_ALL<-cbind(sub_energy_subpop_14_ALL, sub_energy_subpop_14)
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 
        sub_energy_list_sps_loop<-c("ID_pixel",sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_1_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_2_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_3_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_4_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_5_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_6_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_7_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_8_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_9_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_10_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_11_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_12_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_13_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_14_ALL)<-c(sub_energy_list_sps_loop)
        sub_energy_Europe_current<-rbind(sub_energy_subpop_1_ALL,
                  sub_energy_subpop_2_ALL,
                  sub_energy_subpop_3_ALL,
                  sub_energy_subpop_4_ALL,
                  sub_energy_subpop_5_ALL,
                  sub_energy_subpop_6_ALL,
                  sub_energy_subpop_7_ALL,
                  sub_energy_subpop_8_ALL,
                  sub_energy_subpop_9_ALL,
                  sub_energy_subpop_10_ALL,
                  sub_energy_subpop_11_ALL,
                  sub_energy_subpop_12_ALL,
                  sub_energy_subpop_13_ALL,
                  sub_energy_subpop_14_ALL)
        setwd(path_data)
        sub_energy_Europe_current_option_13_1_1<-sub_energy_Europe_current
        save(sub_energy_Europe_current_option_13_1_1, file=paste0(path_save,"sub_energy_Europe_current_option_13_1_1.RData"))

  ####################################################################################################################################################################################
  #5.1.2 Calculation of energy variables from all species for RCP26 Scenario         
  ####################################################################################################################################################################################
  
    #5.1.2.1 Taking into account all species (wild and domestic)    
      #5.1.2.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
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
        #We create a Data frame numeric for species energy matrix
        subpopulations_diet_sub_pop_numeric<-subpopulations_diet
        colnames(subpopulations_diet_sub_pop_numeric)<-c("species","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
        sps_in_loop<-c()
        DF_habitat_continuous_RCP26<-c()
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
            #We load the prediction of habitat for RCP26 scenario:
            #We load the essemble forecast in continuous:
            load(file="Sps_esemble_RCP26_proj.RData")
            habitat_continuous<-(Sps_esemble_RCP26_proj@ proj@ val[,1])
            DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP26")
            energy_of_species_by_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop))
            energy_of_species_by_subpop<-energy_of_species_by_subpop[,2:15]
            dir.create(path_scenario_energy_sps)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:NROW(supopulations_present)){#The number of species to model
              subpopulation_in_loop<-supopulations_present[j]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
              sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
              energy_of_species_by_subpop<-energy_of_species_by_subpop[,c(4,11,2,5,13,14,1,6,10,7,12,9,3,8)]
              colnames(energy_of_species_by_subpop)<-c(1:14)
              sub_energy_subpop_diet <- subset(energy_of_species_by_subpop, select=subpopulation_in_loop)
              vec_sub_energy_subpop_diet<-as.vector(sub_energy_subpop_diet[1,1])
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous*vec_sub_energy_subpop_diet*100
              assign(paste0("sub_energy_subpop_",subpopulation_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_in_loop, ".Rdata"))   
            }  
            #For each the subpopulations where NO has been detected the species:
            for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
              subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
              sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
              sub_habitat_subpop<-sub_habitat_subpop_no_present
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous
              assign(paste0("sub_energy_subpop_",subpopulation_no_present_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
            }  
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 

      #5.1.2.1.2 We merge for each subpopulations the values of spatial energy for all species
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
        sub_energy_subpop_1_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
        sub_energy_subpop_2_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
        sub_energy_subpop_3_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
        sub_energy_subpop_4_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
        sub_energy_subpop_5_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
        sub_energy_subpop_6_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
        sub_energy_subpop_7_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
        sub_energy_subpop_8_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
        sub_energy_subpop_9_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
        sub_energy_subpop_10_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
        sub_energy_subpop_11_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
        sub_energy_subpop_12_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
        sub_energy_subpop_13_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
        sub_energy_subpop_14_ALL<-subset(DF_Scen_RCP26_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
        sub_energy_list_sps_loop<-c()
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
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            #path_scenario = paste0(path_habitat,"/RCP26")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP26")
            sub_energy_list_sps_loop<-cbind(sub_energy_list_sps_loop,i)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:14){#The number of species to model
              load(file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", j, ".Rdata"))
            }  
            sub_energy_subpop_1_ALL<-cbind(sub_energy_subpop_1_ALL, sub_energy_subpop_1)
            sub_energy_subpop_2_ALL<-cbind(sub_energy_subpop_2_ALL, sub_energy_subpop_2)
            sub_energy_subpop_3_ALL<-cbind(sub_energy_subpop_3_ALL, sub_energy_subpop_3)
            sub_energy_subpop_4_ALL<-cbind(sub_energy_subpop_4_ALL, sub_energy_subpop_4)
            sub_energy_subpop_5_ALL<-cbind(sub_energy_subpop_5_ALL, sub_energy_subpop_5)
            sub_energy_subpop_6_ALL<-cbind(sub_energy_subpop_6_ALL, sub_energy_subpop_6)
            sub_energy_subpop_7_ALL<-cbind(sub_energy_subpop_7_ALL, sub_energy_subpop_7)
            sub_energy_subpop_8_ALL<-cbind(sub_energy_subpop_8_ALL, sub_energy_subpop_8)
            sub_energy_subpop_9_ALL<-cbind(sub_energy_subpop_9_ALL, sub_energy_subpop_9)
            sub_energy_subpop_10_ALL<-cbind(sub_energy_subpop_10_ALL, sub_energy_subpop_10)
            sub_energy_subpop_11_ALL<-cbind(sub_energy_subpop_11_ALL, sub_energy_subpop_11)
            sub_energy_subpop_12_ALL<-cbind(sub_energy_subpop_12_ALL, sub_energy_subpop_12)
            sub_energy_subpop_13_ALL<-cbind(sub_energy_subpop_13_ALL, sub_energy_subpop_13)
            sub_energy_subpop_14_ALL<-cbind(sub_energy_subpop_14_ALL, sub_energy_subpop_14)
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 
        sub_energy_list_sps_loop<-c("ID_pixel",sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_1_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_2_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_3_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_4_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_5_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_6_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_7_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_8_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_9_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_10_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_11_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_12_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_13_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_14_ALL)<-c(sub_energy_list_sps_loop)
        sub_energy_Europe_RCP26<-rbind(sub_energy_subpop_1_ALL,
                  sub_energy_subpop_2_ALL,
                  sub_energy_subpop_3_ALL,
                  sub_energy_subpop_4_ALL,
                  sub_energy_subpop_5_ALL,
                  sub_energy_subpop_6_ALL,
                  sub_energy_subpop_7_ALL,
                  sub_energy_subpop_8_ALL,
                  sub_energy_subpop_9_ALL,
                  sub_energy_subpop_10_ALL,
                  sub_energy_subpop_11_ALL,
                  sub_energy_subpop_12_ALL,
                  sub_energy_subpop_13_ALL,
                  sub_energy_subpop_14_ALL)
        setwd(path_data)
        sub_energy_Europe_RCP26_option_13_1_1<-sub_energy_Europe_RCP26
        save(sub_energy_Europe_RCP26_option_13_1_1, file=paste0(path_save,"sub_energy_Europe_RCP26_option_13_1_1.RData"))

  ####################################################################################################################################################################################
  #5.1.3 Calculation of energy variables from all species for RCP60 Scenario         
  ####################################################################################################################################################################################
        
    #5.1.3.1 Taking into account all species (wild and domestic)    
      #5.1.3.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
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
        #path_relation = "H:/G/Project_Name/Relation_Energy_Habitat"
        #load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/sub_DF_Scen_RCP60_Matrix.RData")# Copy here to try to be faster
        #ID_pixel_Scen_RCP60_Matrix_noNA_2<-DF_Scen_RCP60_Matrix_noNA_2$ID_pixel
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        head(df_datos_1_276_description_GBIF)
        
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_RCP60_Matrix_noNA_2.RData")# Copy here to try to be faster
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        DF_subpopulation<-as.data.frame(DF_Scen_RCP60_Matrix_noNA_2$vec_subpop_raster)
        colnames(DF_subpopulation)<-c("vec_subpop_raster")
        #We create a Data frame numeric for species energy matrix
        subpopulations_diet_sub_pop_numeric<-subpopulations_diet
        colnames(subpopulations_diet_sub_pop_numeric)<-c("species","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
        sps_in_loop<-c()
        DF_habitat_continuous_RCP60<-c()
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
            #We load the prediction of habitat for RCP60 scenario:
            #We load the essemble forecast in continuous:
            load(file="Sps_esemble_RCP60_proj.RData")
            habitat_continuous<-(Sps_esemble_RCP60_proj@ proj@ val[,1])
            DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP60")
            energy_of_species_by_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop))
            energy_of_species_by_subpop<-energy_of_species_by_subpop[,2:15]
            dir.create(path_scenario_energy_sps)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:NROW(supopulations_present)){#The number of species to model
              subpopulation_in_loop<-supopulations_present[j]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
              sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
              energy_of_species_by_subpop<-energy_of_species_by_subpop[,c(4,11,2,5,13,14,1,6,10,7,12,9,3,8)]
              colnames(energy_of_species_by_subpop)<-c(1:14)
              sub_energy_subpop_diet <- subset(energy_of_species_by_subpop, select=subpopulation_in_loop)
              vec_sub_energy_subpop_diet<-as.vector(sub_energy_subpop_diet[1,1])
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous*vec_sub_energy_subpop_diet*100
              assign(paste0("sub_energy_subpop_",subpopulation_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_in_loop, ".Rdata"))   
            }  
            #For each the subpopulations where NO has been detected the species:
            for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
              subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
              sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
              sub_habitat_subpop<-sub_habitat_subpop_no_present
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous
              assign(paste0("sub_energy_subpop_",subpopulation_no_present_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
            }  
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 

      #5.1.3.1.2 We merge for each subpopulations the values of spatial energy for all species
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
        sub_energy_subpop_1_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
        sub_energy_subpop_2_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
        sub_energy_subpop_3_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
        sub_energy_subpop_4_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
        sub_energy_subpop_5_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
        sub_energy_subpop_6_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
        sub_energy_subpop_7_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
        sub_energy_subpop_8_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
        sub_energy_subpop_9_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
        sub_energy_subpop_10_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
        sub_energy_subpop_11_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
        sub_energy_subpop_12_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
        sub_energy_subpop_13_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
        sub_energy_subpop_14_ALL<-subset(DF_Scen_RCP60_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
        sub_energy_list_sps_loop<-c()
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
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP60")
            sub_energy_list_sps_loop<-cbind(sub_energy_list_sps_loop,i)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:14){#The number of species to model
              load(file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", j, ".Rdata"))
            }  
            sub_energy_subpop_1_ALL<-cbind(sub_energy_subpop_1_ALL, sub_energy_subpop_1)
            sub_energy_subpop_2_ALL<-cbind(sub_energy_subpop_2_ALL, sub_energy_subpop_2)
            sub_energy_subpop_3_ALL<-cbind(sub_energy_subpop_3_ALL, sub_energy_subpop_3)
            sub_energy_subpop_4_ALL<-cbind(sub_energy_subpop_4_ALL, sub_energy_subpop_4)
            sub_energy_subpop_5_ALL<-cbind(sub_energy_subpop_5_ALL, sub_energy_subpop_5)
            sub_energy_subpop_6_ALL<-cbind(sub_energy_subpop_6_ALL, sub_energy_subpop_6)
            sub_energy_subpop_7_ALL<-cbind(sub_energy_subpop_7_ALL, sub_energy_subpop_7)
            sub_energy_subpop_8_ALL<-cbind(sub_energy_subpop_8_ALL, sub_energy_subpop_8)
            sub_energy_subpop_9_ALL<-cbind(sub_energy_subpop_9_ALL, sub_energy_subpop_9)
            sub_energy_subpop_10_ALL<-cbind(sub_energy_subpop_10_ALL, sub_energy_subpop_10)
            sub_energy_subpop_11_ALL<-cbind(sub_energy_subpop_11_ALL, sub_energy_subpop_11)
            sub_energy_subpop_12_ALL<-cbind(sub_energy_subpop_12_ALL, sub_energy_subpop_12)
            sub_energy_subpop_13_ALL<-cbind(sub_energy_subpop_13_ALL, sub_energy_subpop_13)
            sub_energy_subpop_14_ALL<-cbind(sub_energy_subpop_14_ALL, sub_energy_subpop_14)
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 
        sub_energy_list_sps_loop<-c("ID_pixel",sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_1_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_2_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_3_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_4_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_5_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_6_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_7_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_8_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_9_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_10_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_11_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_12_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_13_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_14_ALL)<-c(sub_energy_list_sps_loop)
        sub_energy_Europe_RCP60<-rbind(sub_energy_subpop_1_ALL,
                  sub_energy_subpop_2_ALL,
                  sub_energy_subpop_3_ALL,
                  sub_energy_subpop_4_ALL,
                  sub_energy_subpop_5_ALL,
                  sub_energy_subpop_6_ALL,
                  sub_energy_subpop_7_ALL,
                  sub_energy_subpop_8_ALL,
                  sub_energy_subpop_9_ALL,
                  sub_energy_subpop_10_ALL,
                  sub_energy_subpop_11_ALL,
                  sub_energy_subpop_12_ALL,
                  sub_energy_subpop_13_ALL,
                  sub_energy_subpop_14_ALL)
        setwd(path_data)
        sub_energy_Europe_RCP60_option_13_1_1<-sub_energy_Europe_RCP60
        save(sub_energy_Europe_RCP60_option_13_1_1, file=paste0(path_save,"sub_energy_Europe_RCP60_option_13_1_1.RData"))

  ####################################################################################################################################################################################
  #5.1.4 Calculation of energy variables from all species for RCP85 Scenario         
  ####################################################################################################################################################################################
    #5.1.4.1 Taking into account all species (wild and domestic)    
      #5.1.4.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
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
        #We create a Data frame numeric for species energy matrix
        subpopulations_diet_sub_pop_numeric<-subpopulations_diet
        colnames(subpopulations_diet_sub_pop_numeric)<-c("species","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
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
            habitat_continuous<-(Sps_esemble_RCP85_proj@ proj@ val[,1])
            DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP85")
            energy_of_species_by_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop))
            energy_of_species_by_subpop<-energy_of_species_by_subpop[,2:15]
            dir.create(path_scenario_energy_sps)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:NROW(supopulations_present)){#The number of species to model
              subpopulation_in_loop<-supopulations_present[j]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
              sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
              energy_of_species_by_subpop<-energy_of_species_by_subpop[,c(4,11,2,5,13,14,1,6,10,7,12,9,3,8)]
              colnames(energy_of_species_by_subpop)<-c(1:14)
              sub_energy_subpop_diet <- subset(energy_of_species_by_subpop, select=subpopulation_in_loop)
              vec_sub_energy_subpop_diet<-as.vector(sub_energy_subpop_diet[1,1])
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous*vec_sub_energy_subpop_diet*100
              assign(paste0("sub_energy_subpop_",subpopulation_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_in_loop, ".Rdata"))   
            }  
            #For each the subpopulations where NO has been detected the species:
            for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
              subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
              sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
              sub_habitat_subpop<-sub_habitat_subpop_no_present
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous
              assign(paste0("sub_energy_subpop_",subpopulation_no_present_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
            }  
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 

      #5.1.4.1.2 We merge for each subpopulations the values of spatial energy for all species
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
        sub_energy_subpop_1_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
        sub_energy_subpop_2_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
        sub_energy_subpop_3_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
        sub_energy_subpop_4_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
        sub_energy_subpop_5_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
        sub_energy_subpop_6_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
        sub_energy_subpop_7_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
        sub_energy_subpop_8_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
        sub_energy_subpop_9_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
        sub_energy_subpop_10_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
        sub_energy_subpop_11_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
        sub_energy_subpop_12_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
        sub_energy_subpop_13_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
        sub_energy_subpop_14_ALL<-subset(DF_Scen_RCP85_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
        sub_energy_list_sps_loop<-c()
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
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/RCP85")
            sub_energy_list_sps_loop<-cbind(sub_energy_list_sps_loop,i)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:14){#The number of species to model
              load(file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", j, ".Rdata"))
            }  
            sub_energy_subpop_1_ALL<-cbind(sub_energy_subpop_1_ALL, sub_energy_subpop_1)
            sub_energy_subpop_2_ALL<-cbind(sub_energy_subpop_2_ALL, sub_energy_subpop_2)
            sub_energy_subpop_3_ALL<-cbind(sub_energy_subpop_3_ALL, sub_energy_subpop_3)
            sub_energy_subpop_4_ALL<-cbind(sub_energy_subpop_4_ALL, sub_energy_subpop_4)
            sub_energy_subpop_5_ALL<-cbind(sub_energy_subpop_5_ALL, sub_energy_subpop_5)
            sub_energy_subpop_6_ALL<-cbind(sub_energy_subpop_6_ALL, sub_energy_subpop_6)
            sub_energy_subpop_7_ALL<-cbind(sub_energy_subpop_7_ALL, sub_energy_subpop_7)
            sub_energy_subpop_8_ALL<-cbind(sub_energy_subpop_8_ALL, sub_energy_subpop_8)
            sub_energy_subpop_9_ALL<-cbind(sub_energy_subpop_9_ALL, sub_energy_subpop_9)
            sub_energy_subpop_10_ALL<-cbind(sub_energy_subpop_10_ALL, sub_energy_subpop_10)
            sub_energy_subpop_11_ALL<-cbind(sub_energy_subpop_11_ALL, sub_energy_subpop_11)
            sub_energy_subpop_12_ALL<-cbind(sub_energy_subpop_12_ALL, sub_energy_subpop_12)
            sub_energy_subpop_13_ALL<-cbind(sub_energy_subpop_13_ALL, sub_energy_subpop_13)
            sub_energy_subpop_14_ALL<-cbind(sub_energy_subpop_14_ALL, sub_energy_subpop_14)
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 
        sub_energy_list_sps_loop<-c("ID_pixel",sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_1_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_2_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_3_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_4_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_5_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_6_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_7_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_8_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_9_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_10_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_11_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_12_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_13_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_14_ALL)<-c(sub_energy_list_sps_loop)
        sub_energy_Europe_RCP85<-rbind(sub_energy_subpop_1_ALL,
                  sub_energy_subpop_2_ALL,
                  sub_energy_subpop_3_ALL,
                  sub_energy_subpop_4_ALL,
                  sub_energy_subpop_5_ALL,
                  sub_energy_subpop_6_ALL,
                  sub_energy_subpop_7_ALL,
                  sub_energy_subpop_8_ALL,
                  sub_energy_subpop_9_ALL,
                  sub_energy_subpop_10_ALL,
                  sub_energy_subpop_11_ALL,
                  sub_energy_subpop_12_ALL,
                  sub_energy_subpop_13_ALL,
                  sub_energy_subpop_14_ALL)
        setwd(path_data)
        sub_energy_Europe_RCP85_option_13_1_1<-sub_energy_Europe_RCP85
        save(sub_energy_Europe_RCP85_option_13_1_1, file=paste0(path_save,"sub_energy_Europe_RCP85_option_13_1_1.RData"))

####################################################################################################################################################################################
####################################################################################################################################################################################
#5.2 Calculation of biotic interactions without account the variation in rEDEC among species
####################################################################################################################################################################################
####################################################################################################################################################################################
  # This idea is the  concept in how biotic interactions are usually calculated, they assum no variation in the importance of the species. For develop this scenario/method/variables we will calculate
  # we  will assume the same interaction for all species (no account for variaion in rEDEC). Thus in summary, we are going to sum
  # all habitat predictions within each subpopulation
        
  ####################################################################################################################################################################################
  #5.2.1 Calculation of energy variables from all species for Current Scenario 
  ####################################################################################################################################################################################
    #5.2.1.1 Taking into account all species (wild and domestic)    
      #5.2.1.1.1 We are going to create for each species a DF for each subpopulations with the predictions or with NA values
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
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
        colnames(DF_subpopulation)<-c("vec_subpop_raster")
        #We create a Data frame numeric for species energy matrix
        subpopulations_diet_sub_pop_numeric<-subpopulations_diet
        colnames(subpopulations_diet_sub_pop_numeric)<-c("species","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
        sps_in_loop<-c()
        DF_habitat_continuous_current<-c()
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
            #we are going to change the names of the populations to numbers
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
            habitat_continuous<-(Sps_esemble_curr_proj@ proj@ val[,1])
            DF_habitat_species_detected<-cbind(DF_supop_diet_present,habitat_continuous)
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_homogeneous_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/Current")
            energy_of_species_by_subpop <- subset(subpopulations_diet_sub_pop_numeric, species==as.character(sps_in_loop))
            energy_of_species_by_subpop<-energy_of_species_by_subpop[,2:15]
            dir.create(path_energy_sps)
            dir.create(path_scenario_energy_sps)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:NROW(supopulations_present)){#The number of species to model
              subpopulation_in_loop<-supopulations_present[j]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop j = ",j," ## Starting Subpopulation ",subpopulation_in_loop,"  ###  "))
              sub_habitat_subpop<- subset(DF_habitat_species_detected, vec_subpop_raster==subpopulation_in_loop)
              energy_of_species_by_subpop<-energy_of_species_by_subpop[,c(4,11,2,5,13,14,1,6,10,7,12,9,3,8)]
              colnames(energy_of_species_by_subpop)<-c(1:14)
              sub_energy_subpop_diet <- subset(energy_of_species_by_subpop, select=subpopulation_in_loop)
              vec_sub_energy_subpop_diet<-as.vector(sub_energy_subpop_diet[1,1])
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous#*vec_sub_energy_subpop_diet*100
              assign(paste0("sub_energy_subpop_",subpopulation_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_in_loop, ".Rdata"))   
            }  
            #For each the subpopulations where NO has been detected the species:
            for (k in 1:NROW(supopulations_NO_present)){#The number of species to model
              subpopulation_no_present_in_loop<-supopulations_NO_present[k]  
              print("             ########           Starting a new Subpopulation             ########          ")
              print(paste0("Loop k = ",k," ## Starting Subpopulation ",subpopulation_no_present_in_loop,"  ###  "))
              sub_habitat_subpop_no_present<- subset(DF_supop_diet_NO_present, vec_subpop_raster==subpopulation_no_present_in_loop)
              sub_habitat_subpop<-sub_habitat_subpop_no_present
              sub_energy_subpop<- sub_habitat_subpop$habitat_continuous
              assign(paste0("sub_energy_subpop_",subpopulation_no_present_in_loop),sub_energy_subpop)
              save(list=paste0("sub_energy_subpop_", subpopulation_no_present_in_loop), file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", subpopulation_no_present_in_loop, ".Rdata"))   
            }  
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 

      #5.2.1.1.2 We merge for each subpopulations the values of spatial energy for all species
        rm(list=ls())
        path_data = "H:/G/Project_Name/Results_Biomod/R_analysis/"
        path_save = paste0(path_data,"SDM_MOD")
        load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/DF_Scen_current_Matrix_noNA_2.RData")# Copy here to try to be faster
        load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
        DF_subpopulation<-as.data.frame(DF_Scen_current_Matrix_noNA_2$vec_subpop_raster)
        colnames(DF_subpopulation)<-c("vec_subpop_raster")
        sps_in_loop<-c()
        DF_habitat_continuous_current<-c()
        levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
        df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
        df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
        sub_energy_subpop_1_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==1, select=c(ID_pixel))
        sub_energy_subpop_2_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==2, select=c(ID_pixel))
        sub_energy_subpop_3_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==3, select=c(ID_pixel))
        sub_energy_subpop_4_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==4, select=c(ID_pixel))
        sub_energy_subpop_5_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==5, select=c(ID_pixel))
        sub_energy_subpop_6_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==6, select=c(ID_pixel))
        sub_energy_subpop_7_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==7, select=c(ID_pixel))
        sub_energy_subpop_8_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==8, select=c(ID_pixel))
        sub_energy_subpop_9_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==9, select=c(ID_pixel))
        sub_energy_subpop_10_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==10, select=c(ID_pixel))
        sub_energy_subpop_11_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==11, select=c(ID_pixel))
        sub_energy_subpop_12_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==12, select=c(ID_pixel))
        sub_energy_subpop_13_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==13, select=c(ID_pixel))
        sub_energy_subpop_14_ALL<-subset(DF_Scen_current_Matrix_noNA_2,vec_subpop_raster==14, select=c(ID_pixel))
        sub_energy_list_sps_loop<-c()
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
            path_energy_sps = paste0(path_save_sps_Biomod,"/Energy_homogeneous_predicition_sps")
            path_scenario_energy_sps = paste0(path_energy_sps,"/Current")
            sub_energy_list_sps_loop<-cbind(sub_energy_list_sps_loop,i)
            #For each the subpopulations where it has been detected the species:
            for (j in 1:14){#The number of species to model
              load(file=paste0(path_scenario_energy_sps,"/sub_energy_subpop_", j, ".Rdata"))
            }  
            sub_energy_subpop_1_ALL<-cbind(sub_energy_subpop_1_ALL, sub_energy_subpop_1)
            sub_energy_subpop_2_ALL<-cbind(sub_energy_subpop_2_ALL, sub_energy_subpop_2)
            sub_energy_subpop_3_ALL<-cbind(sub_energy_subpop_3_ALL, sub_energy_subpop_3)
            sub_energy_subpop_4_ALL<-cbind(sub_energy_subpop_4_ALL, sub_energy_subpop_4)
            sub_energy_subpop_5_ALL<-cbind(sub_energy_subpop_5_ALL, sub_energy_subpop_5)
            sub_energy_subpop_6_ALL<-cbind(sub_energy_subpop_6_ALL, sub_energy_subpop_6)
            sub_energy_subpop_7_ALL<-cbind(sub_energy_subpop_7_ALL, sub_energy_subpop_7)
            sub_energy_subpop_8_ALL<-cbind(sub_energy_subpop_8_ALL, sub_energy_subpop_8)
            sub_energy_subpop_9_ALL<-cbind(sub_energy_subpop_9_ALL, sub_energy_subpop_9)
            sub_energy_subpop_10_ALL<-cbind(sub_energy_subpop_10_ALL, sub_energy_subpop_10)
            sub_energy_subpop_11_ALL<-cbind(sub_energy_subpop_11_ALL, sub_energy_subpop_11)
            sub_energy_subpop_12_ALL<-cbind(sub_energy_subpop_12_ALL, sub_energy_subpop_12)
            sub_energy_subpop_13_ALL<-cbind(sub_energy_subpop_13_ALL, sub_energy_subpop_13)
            sub_energy_subpop_14_ALL<-cbind(sub_energy_subpop_14_ALL, sub_energy_subpop_14)
            print(paste0("             ######## Finished at ",Sys.time(),"########          "))
            print("############################################################################################################################")
            print("############################################################################################################################")
            print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
          }        
        } 
        sub_energy_list_sps_loop<-c("ID_pixel",sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_1_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_2_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_3_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_4_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_5_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_6_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_7_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_8_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_9_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_10_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_11_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_12_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_13_ALL)<-c(sub_energy_list_sps_loop)
        colnames(sub_energy_subpop_14_ALL)<-c(sub_energy_list_sps_loop)
        sub_energy_Europe_current<-rbind(sub_energy_subpop_1_ALL,
                  sub_energy_subpop_2_ALL,
                  sub_energy_subpop_3_ALL,
                  sub_energy_subpop_4_ALL,
                  sub_energy_subpop_5_ALL,
                  sub_energy_subpop_6_ALL,
                  sub_energy_subpop_7_ALL,
                  sub_energy_subpop_8_ALL,
                  sub_energy_subpop_9_ALL,
                  sub_energy_subpop_10_ALL,
                  sub_energy_subpop_11_ALL,
                  sub_energy_subpop_12_ALL,
                  sub_energy_subpop_13_ALL,
                  sub_energy_subpop_14_ALL)
        setwd(path_data)
        sub_energy_Europe_current_option_13_1_1_homogeneous<-sub_energy_Europe_current
        save(sub_energy_Europe_current_option_13_1_1_homogeneous, file=paste0(path_save,"sub_energy_Europe_current_option_13_1_1_homogeneous.RData"))
