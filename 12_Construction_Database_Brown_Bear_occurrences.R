
############################################################################################################
#Readme:
############################################################################################################
#R code to construct the database of brown bear occurrences
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run script:
  #/4_SDM_Food_species/eval_DF_all3.xlsx
  #/10_Descriptive_analysis_Brown_Bear_Food-web/merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF.xlsx
#Data output to other scripts:
  #No outputs to other scripts

##############################################################################################                
#Schema
##############################################################################################                


#12.1 Obtain a unique map with presences/absences for all Europe
#12.2 Calculation of buffers to extract pseudoabsences
#12.3 We merge all datasets (presences/pseudoabsences, buffers, matrix, abiotic variables, biotic variables, subpopulation, others)
#12.4 Do subsampling of presences/pseudoabsences


############################################################################################################
#12.1 Obtain a unique database of presences/absences for all Europe
############################################################################################################
  rm(list=ls()) 
  my_dir <-"writehereyourpath" # path of the folder where your data are stored

  folder_working<-paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences")
  setwd(folder_working)
  library(rgdal) 
  library(raster) 
  library(readxl)    
  
  #Load data
  load("DF_Europe_presences_1.RData")            
  load("DF_Europe_presences_2.RData")    
  load("DF_Europe_presences_3.RData")    
  load("DF_n_presences_1_15.RData")
  load("DF_n_presences_16_35.RData")
  load("DF_n_presences_36.RData")
  load("vec_subpop_raster.RData")

  #Calculation of number of pixels (1x1 km) with brown bear presence for each research group
  n_pixels_with_presence_by_research_group<-rbind(DF_n_presences_1_15,DF_n_presences_16_35,DF_n_presences_36)# 
  write.xlsx(n_pixels_with_presence_by_research_group, file="n_pixels_with_presence_by_research_group.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
  print(paste0("OUTPUT: Supplementary Table 1 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  print(paste0("Supplementary Table 1 (Column Number of pixels by group)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  n_pixels_with_presence_by_research_group

  #Calculation of number of pixels (1x1 km) for each dataset (a pixel from a dataset may overlap with another pixel from another dataset)
  DF_Europe_presences<-cbind(DF_Europe_presences_1,DF_Europe_presences_2,DF_Europe_presences_3)  
  n_pixels_with_presence_by_dataset<-colSums(DF_Europe_presences,na.rm=TRUE)# Supplementary Table 1 (Number of  presences by dataset)
  print(paste0("OUTPUT: Supplementary Table 1 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  print(paste0("Supplementary Table 1 (Column Number of  presences by dataset)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  n_pixels_with_presence_by_dataset
  write.xlsx(n_pixels_with_presence_by_dataset, file="n_pixels_with_presence_by_dataset.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   

  Europe_presences_sum<-rowSums(DF_Europe_presences,na.rm=TRUE)
  Europe_presences_sum2<-Europe_presences_sum
  Europe_presences_sum2[Europe_presences_sum2 >= 1] <- 1
  print(paste0("OUTPUT: Number of presences (1x1km) in the database #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  sum(na.omit(Europe_presences_sum2))# 
  
  Europe_presences_0_1<-Europe_presences_sum2
  save(Europe_presences_0_1,file="Europe_presences_0_1.RData")
  
  #We check the data by subpopulation
  df_presences_subpopulation<-data.frame(Europe_presences_0_1,vec_subpop_raster)
  str(df_presences_subpopulation)
  table_presences_by_subpopulation_initial<-table(df_presences_subpopulation$Europe_presences_0_1,df_presences_subpopulation$vec_subpop_raster)
  write.xlsx(table_presences_by_subpopulation_initial, file="table_presences_by_subpopulation_initial.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
  print(paste0("OUTPUT: Supplementary Table 3. Column n presences in terrestrial #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  print(paste0("Supplementary Table 3. Number of pixels of 1x1 km with presence of brown bear by subpopulation (n presences all systems),Column n presences all systems~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  table_presences_by_subpopulation_initial # # Supplementary Figure 3 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #We are going to use Europe Albers Equal Area Conic  
  newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
  #Raster of reference: 
  Europe_presences_0_1_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
  #We defined the projection of our raster:
  projection(Europe_presences_0_1_raster) <- newproj
  values(Europe_presences_0_1_raster)<-Europe_presences_0_1
  save(Europe_presences_0_1_raster,file="Europe_presences_0_1_raster.RData")
  #We exported the raster
  writeRaster(Europe_presences_0_1_raster, 'Europe_presences_0_1_raster.img',overwrite=TRUE) 

############################################################################################################
#12.2 Calculation of buffers to extract pseudoabsences
############################################################################################################
  rm(list=ls()) 
  my_dir <-"writehereyourpath" # path of the folder where your data are stored
  folder_working<-paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences")
  setwd(folder_working)
  library(rgdal) 
  library(raster) 
  library(readxl)  
  
  #We load the raster of presences
  load("Europe_presences_0_1_raster.RData")
  Europe_presences_0_1_raster2[Europe_presences_0_1_raster2==0]<-NA
  e2<- extent(-2800000,3100000,-800000,4800000) 
  #Buffer of 20km
  buf20km_presences<-buffer(Europe_presences_0_1_raster2, width=20000, dissolve=TRUE) 
  save(buf20km_presences, file="buf20km_presences.Rdata")  
  vec_buf20km_presences<-extract(buf20km_presences,e2)
  save(vec_buf20km_presences, file="vec_buf20km_presences.Rdata")  
  #Buffer of 15km
  buf15km_presences<-buffer(Europe_presences_0_1_raster2, width=15000, dissolve=TRUE) 
  save(buf15km_presences, file="buf15km_presences.Rdata")  
  vec_buf15km_presences<-extract(buf15km_presences,e2)
  save(vec_buf15km_presences, file="vec_buf15km_presences.Rdata")  
  #Buffer of 10km
  buf10km_presences<-buffer(Europe_presences_0_1_raster2, width=10000, dissolve=TRUE) 
  save(buf10km_presences, file="buf10km_presences.Rdata")  
  vec_buf10km_presences<-extract(buf10km_presences,e2)
  save(vec_buf10km_presences, file="vec_buf10km_presences.Rdata")  
  #Buffer of 5km
  buf5km_presences<-buffer(Europe_presences_0_1_raster2, width=5000, dissolve=TRUE) 
  save(buf5km_presences, file="buf5km_presences.Rdata")  
  vec_buf5km_presences<-extract(buf5km_presences,e2)
  save(vec_buf5km_presences, file="vec_buf5km_presences.Rdata")  
  DF_buffers<-as.data.frame(cbind(vec_buf5km_presences,vec_buf10km_presences,vec_buf15km_presences,vec_buf20km_presences))
  save(DF_buffers, file="DF_buffers.Rdata")  

############################################################################################################
#12.3 We merge all datasets (presences/pseudoabsences, buffers, matrix, abiotic variables, biotic variables, subpopulation, others)
############################################################################################################
  rm(list=ls()) 
  my_dir <-"writehereyourpath" # path of the folder where your data are stored
  folder_working<-paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences")
  setwd(folder_working)
 
  library(rgdal) 
  library(raster) 
  library(readxl)    

  #We load the presences
  load("Europe_presences_0_1.RData")
  #We load the buffers
  load("DF_buffers.Rdata")
  #We load the data of coordinates and ID pixels
  load("data_matrix.RData")
  #We load the scenario current abiotic
  load("Scenario_current.RData")
  #We load the scenario current biotic/energy variables NO  energy standarizeD, calculated summer 2020
  load("merged_DF_Variables_Energy_current.RData")
  load("merged_DF_Variables_Energy_current_homogeneous.RData")
  #We load the vector of extrapolation
  load("vec_extrap_raster.RData")
  #We load the subpopulations
  load("vec_subpop_raster.RData")
  load("vec_LC_2_dis.Rdata")
  load("vec_cober_current.RData")
  load("vec_water_current.RData")
  #We change some variable names
  cober<-vec_cober_current
  water<-vec_water_current
  #We elimnate not necessary columns:
  merged_DF_Variables_Energy_current <- merged_DF_Variables_Energy_current[ ,-c(1) ]
  merged_DF_Variables_Energy_current_homogeneous <- merged_DF_Variables_Energy_current_homogeneous[ ,-c(1) ]
  #We create a dataframe with all data
  DF_Scen_current_all<-as.data.frame(cbind(data_matrix,
    Scenario_current,
    merged_DF_Variables_Energy_current,
    merged_DF_Variables_Energy_current_homogeneous,
    vec_extrap_raster,
    vec_subpop_raster,
    Europe_presences_0_1,
    DF_buffers,
    vec_LC_2_dis,
    cober,
    water))
  #We change the NA values of the buffer by 0
  DF_Scen_current_all$vec_buf5km_presences[is.na(DF_Scen_current_all$vec_buf5km_presences)] <- 0
  DF_Scen_current_all$vec_buf10km_presences[is.na(DF_Scen_current_all$vec_buf10km_presences)] <- 0
  DF_Scen_current_all$vec_buf15km_presences[is.na(DF_Scen_current_all$vec_buf15km_presences)] <- 0
  DF_Scen_current_all$vec_buf20km_presences[is.na(DF_Scen_current_all$vec_buf20km_presences)] <- 0
  DF_Scen_current_all$ene_Human_invertebrates[is.na(DF_Scen_current_all$ene_Human_invertebrates)] <- 0
  save(DF_Scen_current_all, file="DF_Scen_current_all.RData")

  #We select the areas inside the subpopulations  
  DF_Scen_current_all_subpopul<-subset(DF_Scen_current_all,vec_subpop_raster >= 1 & vec_subpop_raster < 15)
  #We select the areas with data of land use   
  DF_Scen_current_all_subpopul1<-subset(DF_Scen_current_all_subpopul,LC_1 >= 0 & LC_1 < 2)
  DF_Scen_current_all_subpopul1$factor_subpop<-as.factor(DF_Scen_current_all_subpopul1$vec_subpop_raster)
  #We select the areas without water 
  DF_Scen_current_all_subpopul2<-subset(DF_Scen_current_all_subpopul1,water == 0 )
  #We change the NA values of the buffer by 0
  DF_Scen_current_all_subpopul2$ene_Human_invertebrates[is.na(DF_Scen_current_all_subpopul2$ene_Human_invertebrates)] <- 0
  DF_Scen_current_all_subpopul2$ene_Human_invertebrates_h[is.na(DF_Scen_current_all_subpopul2$ene_Human_invertebrates_h)] <- 0
  DF_Scen_current_all_noNA<- DF_Scen_current_all_subpopul2
  #First we apply the conversions to make equal the data of CRU TS used for historical and the data of CHELSA used for current
    #We need to divide CHELSA data current bioclimatic variables by 10:
    DF_Scen_current_all_noNA$Clim_3<-DF_Scen_current_all_noNA$Clim_3/10
    DF_Scen_current_all_noNA$Clim_4<-DF_Scen_current_all_noNA$Clim_4/10
    DF_Scen_current_all_noNA$Clim_8<-DF_Scen_current_all_noNA$Clim_8/10
    DF_Scen_current_all_noNA$Clim_9<-DF_Scen_current_all_noNA$Clim_9/10
    #Now we need to apply the same transformation used in historical data to scalar moreless the data around 1
    DF_Scen_current_all_noNA$Clim_3<-DF_Scen_current_all_noNA$Clim_3/100
    DF_Scen_current_all_noNA$Clim_4<-DF_Scen_current_all_noNA$Clim_4/1000
    DF_Scen_current_all_noNA$Clim_8<-DF_Scen_current_all_noNA$Clim_8/10
    DF_Scen_current_all_noNA$Clim_9<-DF_Scen_current_all_noNA$Clim_9/10
 
    #Calculation of quadratic variables
    DF_Scen_current_all_noNA$Clim_3_c<-DF_Scen_current_all_noNA$Clim_3^2
    DF_Scen_current_all_noNA$Clim_4_c<-DF_Scen_current_all_noNA$Clim_4^2
    DF_Scen_current_all_noNA$Clim_8_c<-DF_Scen_current_all_noNA$Clim_8^2
    DF_Scen_current_all_noNA$Clim_9_c<-DF_Scen_current_all_noNA$Clim_9^2
  save(DF_Scen_current_all_noNA, file="DF_Scen_current_all_noNA.RData")

############################################################################################################
#12.4 Do subsampling of presences/pseudoabsences
############################################################################################################
  rm(list=ls()) 
  my_dir <-"writehereyourpath" # path of the folder where your data are stored
  folder_working<-paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences")
  setwd(folder_working)
  library(rgdal) 
  library(raster) 
  library(readxl) 
  library(xlsx)
  #library(splitstackshape)
  
  #We load all the merged data
  load("DF_Scen_current_all_noNA.RData")
  #We see the presence by wach subpopulation
  table_presences_by_subpopulation<-table(DF_Scen_current_all_noNA$Europe_presences_0_1,DF_Scen_current_all_noNA$vec_subpop_raster)
  table(DF_Scen_current_all_noNA$Europe_presences_0_1,DF_Scen_current_all_noNA$factor_subpop)
  DF_Scen_current_all_noNA$factor_subpop_integer<-as.integer(DF_Scen_current_all_noNA$factor_subpop)
  write.xlsx(table_presences_by_subpopulation, file="table_presences_by_subpopulation.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
  print(paste0("OUTPUT: Supplementary Table 3. Column n presences in terrestrial #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  print(paste0("Supplementary Table 3. Number of pixels of 1x1 km with presence of brown bear by subpopulation in terrestrial systems (n presences in terrestrial),Column n presences in terrestrial~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  table_presences_by_subpopulation # # Supplementary Figure 3 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
  #We subset the data within the area of each buffer
  DF_Scen_current_5km<-subset(DF_Scen_current_all_noNA,vec_buf5km_presences==1)          
  DF_Scen_current_10km<-subset(DF_Scen_current_all_noNA,vec_buf10km_presences==1)          
  DF_Scen_current_15km<-subset(DF_Scen_current_all_noNA,vec_buf15km_presences==1)          
  DF_Scen_current_20km<-subset(DF_Scen_current_all_noNA,vec_buf20km_presences==1)          

 list_of_df <- list(DF_Scen_current_5km,DF_Scen_current_10km,DF_Scen_current_15km, DF_Scen_current_20km)
 vec_scales<-c("5km","10km","15km","20km")
 
 #make this example reproducible
  set.seed(1)
 
  for (scale in 1:4){ # We apply 4 different scales one for each buffer
    df_scale_in_loop <- as.data.frame(list_of_df[[scale]])
    print(paste0("SCALE ",vec_scales[scale]," #########################################################"))
    RANDOM_presences_absences_for_model_all<-c()
    RANDOM_presences_absences_for_validation_all<-c()
    
    for (subp in 1:14){# We subset the data within each subpopulation, if for a subpopulation we have less than 2000 presences we do not subset
      DF_Scen_current_all_noNA_subp<-subset(df_scale_in_loop,vec_subpop_raster==subp)
      DF_Scen_current_all_noNA_subp_absences<-subset(DF_Scen_current_all_noNA_subp,Europe_presences_0_1==0) 
      DF_Scen_current_all_noNA_subp_presences<-subset(DF_Scen_current_all_noNA_subp,Europe_presences_0_1==1)
      n_presences_subp<-nrow(DF_Scen_current_all_noNA_subp_presences)
      n_absences_subp<-nrow(DF_Scen_current_all_noNA_subp_absences)
      min_absences_presences<-min(n_presences_subp,n_absences_subp)
      if(min_absences_presences>2000){
        min_absences_presences=2000 
      }
      ##Over the cells with absence we select a sample equals to the numer of presences
      absences_balanced_subp<-DF_Scen_current_all_noNA_subp_absences[sample(nrow(DF_Scen_current_all_noNA_subp_absences), min_absences_presences), ]
      presences_balanced_subp<-DF_Scen_current_all_noNA_subp_presences[sample(nrow(DF_Scen_current_all_noNA_subp_presences), min_absences_presences), ]
      ##We select 80% of data for analysis and 20% for validation:
      split_data_absences_DF_FOR_GLMM<-split(absences_balanced_subp, sample(1:nrow(absences_balanced_subp) > round(nrow(absences_balanced_subp) * .20)))
      RANDOM_absences_for_validation<-split_data_absences_DF_FOR_GLMM$'FALSE'
      RANDOM_absences_for_model<-split_data_absences_DF_FOR_GLMM$'TRUE'
      #Presences
      split_data_presences_DF_FOR_GLMM<-split(presences_balanced_subp, sample(1:nrow(presences_balanced_subp) > round(nrow(presences_balanced_subp) * .20)))
      RANDOM_presences_for_validation<-split_data_presences_DF_FOR_GLMM$'FALSE'
      RANDOM_presences_for_model<-split_data_presences_DF_FOR_GLMM$'TRUE'
      #Absences
      RANDOM_presences_absences_for_model<-rbind(RANDOM_presences_for_model,RANDOM_absences_for_model)
      RANDOM_presences_absences_for_validation<-rbind(RANDOM_presences_for_validation,RANDOM_absences_for_validation)
      #A counter of the number of presences for model
      n_RANDOM_presences_for_model_subp<-nrow(RANDOM_presences_for_model)
      #A counter of the number of absences for model
      n_RANDOM_absences_for_model_subp<-nrow(RANDOM_absences_for_model)
        print(paste0("Subpopulation number ",subp,": Number of presences for model = ",n_RANDOM_presences_for_model_subp,"# Number of absences for model = ",n_RANDOM_absences_for_model_subp))   
      RANDOM_presences_absences_for_model_all<-rbind(RANDOM_presences_absences_for_model_all,RANDOM_presences_absences_for_model)
      RANDOM_presences_absences_for_validation_all<-rbind(RANDOM_presences_absences_for_validation_all,RANDOM_presences_absences_for_validation)
    }  
    assign( paste("RANDOM_presences_absences_for_model_all_",vec_scales[scale], sep=""),  RANDOM_presences_absences_for_model_all)
    save(list=paste("RANDOM_presences_absences_for_model_all_",vec_scales[scale], sep=""), file=paste("RANDOM_presences_absences_for_model_all_",vec_scales[scale],".RData", sep=""))     
    assign( paste("RANDOM_presences_absences_for_validation_all_",vec_scales[scale], sep=""),  RANDOM_presences_absences_for_validation_all)
    save(list=paste("RANDOM_presences_absences_for_validation_all_",vec_scales[scale], sep=""), file=paste("RANDOM_presences_absences_for_validation_all_",vec_scales[scale],".RData", sep=""))     
  }      

