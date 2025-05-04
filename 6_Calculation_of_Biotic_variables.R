

############################################################################################################
#Readme:
############################################################################################################
#R code for calculate biotic variables. We calculate two kind of biotic variables:
  # Biotic variables quantitatives (6.2 in the schema)
  # Biotic variables binary (6.4 in the schema)
#Here we just stack/sum the pixels from individual (for each species) calculations from "5_Calculation_of_Biotic_variables_for_each_sps"
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
 #/3_Download_occurrences_from_GBIF/df_datos_1_276_description_GBIF.RData
  #/2_Representative_Food-web_and_Food_species_list/Species_list_clean_binary.xlsx # The list of species
  #/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names.xlsx
  #For biotic quantitative / accounting the rEDEC:
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/SDM_MODsub_energy_Europe_current_option_13_1_1.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/SDM_MODsub_energy_Europe_RCP26_option_13_1_1.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/SDM_MODsub_energy_Europe_RCP60_option_13_1_1.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/SDM_MODsub_energy_Europe_RCP85_option_13_1_1.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #  load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
  #For biotic binary / homogeneous variable:
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/SDM_MODsub_energy_Europe_current_option_13_1_1_homogeneous.RData")
    #/5_Calculation_of_Biotic_variables_for_each_sps/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
#Data output:
  #For biotic quantitative / accounting the rEDEC:
    # merged_DF_Variables_Energy_current.RData
    # Rasters of energy variables current
    # merged_DF_Variables_Energy_RCP26.RData
    # Rasters of energy variables RCP26
    # merged_DF_Variables_Energy_RCP60.RData
    # Rasters of energy variables RCP60
  #Raster of change in the different RCPS in the biotic quantitative variables:
  #For biotic binary / homogeneous variables:
    # DF_Variables_Energy_current_homogeneous.RData
    # Rasters of energy binary variables current

##############################################################################################                
#Schema
############################################################################################## 
#6_Calculation_of_Biotic_variables
  #6.1 We prepare some tables summarizing/grouping the data
  #6.2 We create subsets of energy by category # QUANTITATIVE option
    #6.2.1 For Current scenario
      #6.2.1.1 We are going to do subset of species based in human/wild origin and in categories of diet
          #6.2.1.1.1 Human and wild species (all species excluding the brown bear)
          #6.2.1.1.2 Human species
          #6.2.1.1.3 Wild species
    #6.2.2 For RCP26 scenario
      #6.2.2.1 We are going to do subset of species based in human/wild origin and in categories of diet
          #6.2.2.1.1 Human and wild species (all species excluding the brown bear)
          #6.2.2.1.2 Human species
          #6.2.2.1.3 Wild species
    #6.2.3 For RCP60 scenario
      #6.2.3.1 We are going to do subset of species based in human/wild origin and in categories of diet
          #6.2.3.1.1 Human and wild species (all species excluding the brown bear)
          #6.2.3.1.2 Human species
          #6.2.3.1.3 Wild species
    #6.2.4 For RCP85 scenario
      #6.2.4.1 We are going to do subset of species based in human/wild origin and in categories of diet
          #6.2.4.1.1 Human and wild species (all species excluding the brown bear)
          #6.2.4.1.2 Human species
          #6.2.4.1.3 Wild species
  #6.3 #CHANGE OF BIOTIC QUANTITATIVE VARIABLES/MAPS # CHANGE FOR THE QUANTITATIVE option from current to future scenarios
    #6.3.1 Biotic delta of change by percentage from current and absolute values summary by subpopulation    
      #6.3.1.1 For ras_current_ene_Wild.img   
      #6.3.1.2 For ras_current_ene_Wild_unk_plant_oth.img   
      #6.3.1.3 For ras_current_ene_Wild_rep_plant.img   
      #6.3.1.4 For ras_current_ene_Wild_veg_plant.img   
      #6.3.1.5 For ras_current_ene_Wild_vertebrates.img   
      #6.3.1.6 For ras_current_ene_Wild_invertebrates_invertebrates.img   
      #6.3.1.7 We bind all groups 
  #6.4 We create subsets of energy by category HOMOGENEOUS OPTION # BINARY option
       #6.4.1.1 We are going to do subset of species based in human/wild origin and in categories of diet
          #6.4.1.1.1 Human and wild species (all species excluding the brown bear)
          #6.4.1.1.2 Human species
          #6.4.1.1.3 Wild species

#######################################################################################################################
#######################################################################################################################
#6.1 We prepare some tables summarizing/grouping the data
#######################################################################################################################
#######################################################################################################################
  rm(list=ls())
  setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
  library(readxl)
  load("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
  levels(df_datos_1_276_description_GBIF$Species) <- c(levels(df_datos_1_276_description_GBIF$Species),"Ficus carica")
  df_datos_1_276_description_GBIF[86,"Species"]<-"Ficus carica"
  df_datos_1_276_description_GBIF$Species<-factor(df_datos_1_276_description_GBIF$Species)
  #WE add the list of species which contain information abour the human origin 
  sheet<-read_excel("H:/G/Project_Name/Database_diet/R_analysis/Descriptive_results/Species_list_clean_binary.xlsx", sheet = 1)
  data_taxo<-as.data.frame(sheet)
  #Data of the reassignated energy
  sheet2<-read_excel("H:/G/Project_Name/Database_diet/Reasignation_of_diet/database_original_with_energy_imputed2_ALL_corrected_names.xlsx", sheet = 1)
  table_all_imputed_energy<-as.data.frame(sheet2)
  table_all_imputed_energy_2<-subset(table_all_imputed_energy,table_all_imputed_energy$Species_corrected_names!="NA")
  table_diet2<-table(table_all_imputed_energy_2$Species_corrected_names,table_all_imputed_energy_2$diet2) 
  matrix_diet2<-as.data.frame.matrix(table_diet2)
  max_values<-colnames(matrix_diet2)[apply(matrix_diet2,1,which.max)]
  df_list_sps_categories_diet2_mode<-as.data.frame(cbind(rownames(matrix_diet2),max_values))
  colnames(df_list_sps_categories_diet2_mode)<-c("species","Diet_category")
  merged_sps_list_diet_category<-merge(df_datos_1_276_description_GBIF, df_list_sps_categories_diet2_mode, by.x="Species", by.y = "species",all.x = T)
  data_taxo2<-data_taxo[,c("Species","Human_origin")]
  merged_sps_list_diet_category2<-merge(merged_sps_list_diet_category, data_taxo2, by.x="Species", by.y = "Species",all.x = T)
  save(merged_sps_list_diet_category2, file="merged_sps_list_diet_category2.RData")
  #We create two similar datafrabes 
    #merged_sps_list_diet_category3 contain less columns which we are not going to need
    merged_sps_list_diet_category3<-merged_sps_list_diet_category2[,c(1,2,15,16)]
    #merged_sps_list_diet_category3b contain all columns which we are not going to need
    merged_sps_list_diet_category3b<-merged_sps_list_diet_category2
  #merged_sps_list_diet_category3 
  merged_sps_list_diet_category3$Diet_category2<-merged_sps_list_diet_category3$Diet_category
  levels(merged_sps_list_diet_category3$Diet_category2) <- c(levels(merged_sps_list_diet_category3$Diet_category2),"unknown_plant_material_and_others")
  merged_sps_list_diet_category3$Diet_category2[merged_sps_list_diet_category3$Diet_category2 == 'fungi_lichens_bryophytes_algae'] <- 'unknown_plant_material_and_others'
  merged_sps_list_diet_category3$Diet_category2[merged_sps_list_diet_category3$Diet_category2 == 'unknown_plant_material'] <- 'unknown_plant_material_and_others'
  merged_sps_list_diet_category3$Diet_category2<-factor(merged_sps_list_diet_category3$Diet_category2)
  save(merged_sps_list_diet_category3, file="merged_sps_list_diet_category3.RData")
  ##merged_sps_list_diet_category3b
  merged_sps_list_diet_category3b$Diet_category2<-merged_sps_list_diet_category3b$Diet_category
  levels(merged_sps_list_diet_category3b$Diet_category2) <- c(levels(merged_sps_list_diet_category3b$Diet_category2),"unknown_plant_material_and_others")
  merged_sps_list_diet_category3b$Diet_category2[merged_sps_list_diet_category3b$Diet_category2 == 'fungi_lichens_bryophytes_algae'] <- 'unknown_plant_material_and_others'
  merged_sps_list_diet_category3b$Diet_category2[merged_sps_list_diet_category3b$Diet_category2 == 'unknown_plant_material'] <- 'unknown_plant_material_and_others'
  merged_sps_list_diet_category3b$Diet_category2<-factor(merged_sps_list_diet_category3b$Diet_category2)
  save(merged_sps_list_diet_category3b, file="merged_sps_list_diet_category3b.RData")
  
  table(merged_sps_list_diet_category3$Diet_category)
  (67+57+86)/276*100 #Percentage plants
  (67)/276*100 #Reprod plant material
  (86)/276*100 #Reprod plant material
  (57)/276*100 #Reprod plant material
      
  (28+36)/276*100 #Percentage plants
  (36)/276*100 #Reprod plant material
  (28)/276*100 #Reprod plant material

  table(merged_sps_list_diet_category3$Diet_category2,merged_sps_list_diet_category3$Human_origin)

#######################################################################################################################
#######################################################################################################################
#6.2 We create subsets of energy by category
#######################################################################################################################
#######################################################################################################################
  
  ##########################################################################################################            
  #6.2.1 For Current scenario
  ##########################################################################################################            
    rm(list=ls())
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    library(readxl)
    #We load the data frame with the energy by pixel of each species
    load("H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODsub_energy_Europe_current_option_13_1_1.RData")
    load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #We change the value of Human/Wild for the reindeer
    merged_sps_list_diet_category3[214,4] <- 0

    colnames_sub_energy_Europe_current_option_13_1_1<-colnames(sub_energy_Europe_current_option_13_1_1[,2:237])
    colnames_sub_energy_Europe_current_option_13_1_1<-as.data.frame(as.numeric(colnames_sub_energy_Europe_current_option_13_1_1))
    colnames(colnames_sub_energy_Europe_current_option_13_1_1)<-c("i")
    merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet<-merge(colnames_sub_energy_Europe_current_option_13_1_1, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)
    str(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet)
    
    #6.2.1.1 We are going to do subset of species based in human/wild origin and in categories of diet
        #6.2.1.1.1 Human and wild species (all species excluding the brown bear)
          merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_All<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet[!merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_All$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_All$Human_origin)
          vector_sub_merged_sps_list_diet_category3_All<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_All$i
          str(vector_sub_merged_sps_list_diet_category3_All)
          vector_sub_merged_sps_list_diet_category3_All_character<-as.character(vector_sub_merged_sps_list_diet_category3_All)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_All<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_All_character]
          sum_sub_energy_Europe_current_option_13_1_1_All<-rowSums(sub_energy_Europe_current_option_13_1_1_All[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_All)
          head(sum_sub_energy_Europe_current_option_13_1_1_All)
          energy_current_All<-sum_sub_energy_Europe_current_option_13_1_1_All
          save(energy_current_All, file="energy_current_All.RData")
  
        #6.2.1.1.2 Human species
          merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet, Human_origin==1)
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human$i
          str(vector_sub_merged_sps_list_diet_category3_Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_origin)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_Human_origin<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_origin_character]
          sum_sub_energy_Europe_current_option_13_1_1_Human_origin<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_origin[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_Human_origin)
          head(sum_sub_energy_Europe_current_option_13_1_1_Human_origin)
          energy_current_Human<-sum_sub_energy_Europe_current_option_13_1_1_Human_origin
          save(energy_current_Human, file="energy_current_Human.RData")

          #Human_invertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_invertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Human_invertebrates<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character]
            #sum_sub_energy_Europe_current_option_13_1_1_Human_invertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_invertebrates[,], na.rm = T)
            sum_sub_energy_Europe_current_option_13_1_1_Human_invertebrates<-sub_energy_Europe_current_option_13_1_1_Human_invertebrates
            str(sub_energy_Europe_current_option_13_1_1_Human_invertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_Human_invertebrates)
            energy_current_Human_invertebrates<-sum_sub_energy_Europe_current_option_13_1_1_Human_invertebrates
            save(energy_current_Human_invertebrates, file="energy_current_Human_invertebrates.RData")
          
          #Human_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_reproductive_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material)
            energy_current_Human_rep_plant<-sum_sub_energy_Europe_current_option_13_1_1_Human_reproductive_plant_material
            save(energy_current_Human_rep_plant, file="energy_current_Human_rep_plant.RData")

          #Human_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vegetative_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material)
            energy_current_Human_veg_plant<-sum_sub_energy_Europe_current_option_13_1_1_Human_vegetative_plant_material
            save(energy_current_Human_veg_plant, file="energy_current_Human_veg_plant.RData")

          #Human_vertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_vertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Human_vertebrates<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_Human_vertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_vertebrates[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Human_vertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_Human_vertebrates)
            energy_current_Human_vertebrates<-sum_sub_energy_Europe_current_option_13_1_1_Human_vertebrates
            save(energy_current_Human_vertebrates, file="energy_current_Human_vertebrates.RData")

          #Human_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others)
            energy_current_Human_unk_plant_oth<-sum_sub_energy_Europe_current_option_13_1_1_Human_unknown_plant_material_and_others
            save(energy_current_Human_unk_plant_oth, file="energy_current_Human_unk_plant_oth.RData")

        #6.2.1.1.3 Wild species
          merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet, Human_origin==0)
          merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild[!merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Wild<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild$i
          str(vector_sub_merged_sps_list_diet_category3_Wild)
          vector_sub_merged_sps_list_diet_category3_Wild_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_Wild<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_character]
          sum_sub_energy_Europe_current_option_13_1_1_Wild<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_Wild)
          head(sum_sub_energy_Europe_current_option_13_1_1_Wild)
          energy_current_Wild<-sum_sub_energy_Europe_current_option_13_1_1_Wild
          save(energy_current_Wild, file="energy_current_Wild.RData")

          #Wild_invertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_invertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Wild_invertebrates<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_Wild_invertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild_invertebrates[,], na.rm = T)
            #sum_sub_energy_Europe_current_option_13_1_1_Wild_invertebrates<-sub_energy_Europe_current_option_13_1_1_Wild_invertebrates
            str(sub_energy_Europe_current_option_13_1_1_Wild_invertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_Wild_invertebrates)
            energy_current_Wild_invertebrates<-sum_sub_energy_Europe_current_option_13_1_1_Wild_invertebrates
            save(energy_current_Wild_invertebrates, file="energy_current_Wild_invertebrates.RData")
          
          #Wild_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_reproductive_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material)
            energy_current_Wild_rep_plant<-sum_sub_energy_Europe_current_option_13_1_1_Wild_reproductive_plant_material
            save(energy_current_Wild_rep_plant, file="energy_current_Wild_rep_plant.RData")

          #Wild_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vegetative_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material)
            energy_current_Wild_veg_plant<-sum_sub_energy_Europe_current_option_13_1_1_Wild_vegetative_plant_material
            save(energy_current_Wild_veg_plant, file="energy_current_Wild_veg_plant.RData")

          #Wild_vertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_vertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Wild_vertebrates<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_Wild_vertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild_vertebrates[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Wild_vertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_Wild_vertebrates)
            energy_current_Wild_vertebrates<-sum_sub_energy_Europe_current_option_13_1_1_Wild_vertebrates
            save(energy_current_Wild_vertebrates, file="energy_current_Wild_vertebrates.RData")

          #Wild_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_current_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others<- sub_energy_Europe_current_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others)
            energy_current_Wild_unk_plant_oth<-sum_sub_energy_Europe_current_option_13_1_1_Wild_unknown_plant_material_and_others
            save(energy_current_Wild_unk_plant_oth, file="energy_current_Wild_unk_plant_oth.RData")

      #We sum the energy for all species in each pixel:
        DF_Variables_Energy_current<- as.data.frame(cbind(sub_energy_Europe_current_option_13_1_1[,1],
          energy_current_All,
          energy_current_Human,
          energy_current_Human_invertebrates,
          energy_current_Human_rep_plant,
          energy_current_Human_veg_plant,
          energy_current_Human_vertebrates,
          energy_current_Human_unk_plant_oth,
          energy_current_Wild,
          energy_current_Wild_invertebrates,
          energy_current_Wild_rep_plant,
          energy_current_Wild_veg_plant,
          energy_current_Wild_vertebrates,
          energy_current_Wild_unk_plant_oth))
        head(DF_Variables_Energy_current)
        colnames(DF_Variables_Energy_current)<-c("ID_pixel",
          "ene_All",
          "ene_Human",
          "ene_Human_invertebrates",
          "ene_Human_rep_plant",
          "ene_Human_veg_plant",
          "ene_Human_vertebrates",
          "ene_Human_unk_plant_oth",
          "ene_Wild",
          "ene_Wild_invertebrates",
          "ene_Wild_rep_plant",
          "ene_Wild_veg_plant",
          "ene_Wild_vertebrates",
          "ene_Wild_unk_plant_oth")
        head(DF_Variables_Energy_current)
        save(DF_Variables_Energy_current, file="DF_Variables_Energy_current.RData")

        #WE do the merge
          rm(list=ls())
          setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
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
          load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
          load(file="DF_Variables_Energy_current.RData")
          merged_DF_Variables_Energy_current<-merge(DF_ID_pixel, DF_Variables_Energy_current, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_DF_Variables_Energy_current, file="merged_DF_Variables_Energy_current.RData")

          #We save into a raster file to see in Arcgis
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          
            #For ene_All
              ras_current_ene_All <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_All) <- newproj
              values(ras_current_ene_All)<-merged_DF_Variables_Energy_current$ene_All
              writeRaster(ras_current_ene_All, "ras_current_ene_All.img", overwrite=TRUE)
              
            #For ene_Human
              ras_current_ene_Human <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human) <- newproj
              values(ras_current_ene_Human)<-merged_DF_Variables_Energy_current$ene_Human
              writeRaster(ras_current_ene_Human, "ras_current_ene_Human.img", overwrite=TRUE)

            #For ene_Human_invertebrates
              ras_current_ene_Human_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_invertebrates) <- newproj
              values(ras_current_ene_Human_invertebrates)<-merged_DF_Variables_Energy_current$ene_Human_invertebrates
              writeRaster(ras_current_ene_Human_invertebrates, "ras_current_ene_Human_invertebrates.img", overwrite=TRUE)
              
             #For ene_Human_rep_plant
              ras_current_ene_Human_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_rep_plant) <- newproj
              values(ras_current_ene_Human_rep_plant)<-merged_DF_Variables_Energy_current$ene_Human_rep_plant
              writeRaster(ras_current_ene_Human_rep_plant, "ras_current_ene_Human_rep_plant.img", overwrite=TRUE)
              
            #For ene_Human_veg_plant
              ras_current_ene_Human_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_veg_plant) <- newproj
              values(ras_current_ene_Human_veg_plant)<-merged_DF_Variables_Energy_current$ene_Human_veg_plant
              writeRaster(ras_current_ene_Human_veg_plant, "ras_current_ene_Human_veg_plant.img", overwrite=TRUE)
             
              
             #For ene_Human_vertebrates
              ras_current_ene_Human_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_vertebrates) <- newproj
              values(ras_current_ene_Human_vertebrates)<-merged_DF_Variables_Energy_current$ene_Human_vertebrates
              writeRaster(ras_current_ene_Human_vertebrates, "ras_current_ene_Human_vertebrates.img", overwrite=TRUE)
              
             #For ene_Human_unk_plant_oth
              ras_current_ene_Human_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_unk_plant_oth) <- newproj
              values(ras_current_ene_Human_unk_plant_oth)<-merged_DF_Variables_Energy_current$ene_Human_unk_plant_oth
              writeRaster(ras_current_ene_Human_unk_plant_oth, "ras_current_ene_Human_unk_plant_oth.img", overwrite=TRUE)
              
            #For ene_Wild
              ras_current_ene_Wild <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild) <- newproj
              values(ras_current_ene_Wild)<-merged_DF_Variables_Energy_current$ene_Wild
              writeRaster(ras_current_ene_Wild, "ras_current_ene_Wild.img", overwrite=TRUE) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Supplementary Figure 6b #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Supplementary Figure 6b. Maps with change in biotic variables for future scenarios. Bio_All_species for current ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

            #For ene_Wild_invertebrates
              ras_current_ene_Wild_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_invertebrates) <- newproj
              values(ras_current_ene_Wild_invertebrates)<-merged_DF_Variables_Energy_current$ene_Wild_invertebrates
              writeRaster(ras_current_ene_Wild_invertebrates, "ras_current_ene_Wild_invertebrates.img", overwrite=TRUE)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Figure 1d and Supplementary Figure 6f#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Figure 1d and Supplementary Figure 6f. Maps with change in biotic variables for future scenarios. Bio_Invertebrates for current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

             #For ene_Wild_rep_plant
              ras_current_ene_Wild_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_rep_plant) <- newproj
              values(ras_current_ene_Wild_rep_plant)<-merged_DF_Variables_Energy_current$ene_Wild_rep_plant
              writeRaster(ras_current_ene_Wild_rep_plant, "ras_current_ene_Wild_rep_plant.img", overwrite=TRUE)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Figure 1d and Supplementary Figure 6c#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Figure 1d and Supplementary Figure 6c. Maps with change in biotic variables for future scenarios. Bio_Reproductive_plant for current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

            #For ene_Wild_veg_plant
              ras_current_ene_Wild_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_veg_plant) <- newproj
              values(ras_current_ene_Wild_veg_plant)<-merged_DF_Variables_Energy_current$ene_Wild_veg_plant
              writeRaster(ras_current_ene_Wild_veg_plant, "ras_current_ene_Wild_veg_plant.img", overwrite=TRUE)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Figure 1d and Supplementary Figure 6d#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Figure 1d and Supplementary Figure 6d. Maps with change in biotic variables for future scenarios. Bio_Vegetative_plant for current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

             #For ene_Wild_vertebrates
              ras_current_ene_Wild_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_vertebrates) <- newproj
              values(ras_current_ene_Wild_vertebrates)<-merged_DF_Variables_Energy_current$ene_Wild_vertebrates
              writeRaster(ras_current_ene_Wild_vertebrates, "ras_current_ene_Wild_vertebrates.img", overwrite=TRUE)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Figure 1d and Supplementary Figure 6g#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Figure 1d and Supplementary Figure 6g. Maps with change in biotic variables for future scenarios. Bio_Vertebrates for current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

             #For ene_Wild_unk_plant_oth
              ras_current_ene_Wild_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_unk_plant_oth) <- newproj
              values(ras_current_ene_Wild_unk_plant_oth)<-merged_DF_Variables_Energy_current$ene_Wild_unk_plant_oth
              writeRaster(ras_current_ene_Wild_unk_plant_oth, "ras_current_ene_Wild_unk_plant_oth.img", overwrite=TRUE)#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              print(paste0("OUTPUT: Figure 1d and Supplementary Figure 6e#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
              print(paste0("Figure 1d and Supplementary Figure 6e. Maps with change in biotic variables for future scenarios. Bio_Unknown_plant for current~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

   
  ##########################################################################################################            
  #6.2.2 For RCP26 scenario
  ##########################################################################################################            
    rm(list=ls())
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    library(readxl)
    #We load the data frame with the energy by pixel of each species
    load("H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODsub_energy_Europe_RCP26_option_13_1_1.RData")
    load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #We change the value of Human/Wild for the reindeer
    View(merged_sps_list_diet_category3)
    merged_sps_list_diet_category3[214,]
    merged_sps_list_diet_category3[214,4] <- 0
    merged_sps_list_diet_category3[214,]
    colnames_sub_energy_Europe_RCP26_option_13_1_1<-colnames(sub_energy_Europe_RCP26_option_13_1_1[,2:237])
    colnames_sub_energy_Europe_RCP26_option_13_1_1<-as.data.frame(as.numeric(colnames_sub_energy_Europe_RCP26_option_13_1_1))
    colnames(colnames_sub_energy_Europe_RCP26_option_13_1_1)<-c("i")
    merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet<-merge(colnames_sub_energy_Europe_RCP26_option_13_1_1, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)

    #6.2.2.1 We are going to do subset of species based in human/wild origin and in categories of diet  
        #6.2.2.1.1 Human and wild species (all species excluding the brown bear)
          merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_All<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet[!merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_All$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_All$Human_origin)
          vector_sub_merged_sps_list_diet_category3_All<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_All$i
          vector_sub_merged_sps_list_diet_category3_All_character<-as.character(vector_sub_merged_sps_list_diet_category3_All)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP26_option_13_1_1_All<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_All_character]
          sum_sub_energy_Europe_RCP26_option_13_1_1_All<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_All[,], na.rm = T)
          energy_RCP26_All<-sum_sub_energy_Europe_RCP26_option_13_1_1_All
          save(energy_RCP26_All, file="energy_RCP26_All.RData")
  
        #6.2.2.1.2 Human species
          merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet, Human_origin==1)
          table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human$i
          vector_sub_merged_sps_list_diet_category3_Human_origin_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_origin)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP26_option_13_1_1_Human_origin<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_origin_character]
          sum_sub_energy_Europe_RCP26_option_13_1_1_Human_origin<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_origin[,], na.rm = T)
          energy_RCP26_Human<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_origin
          save(energy_RCP26_Human, file="energy_RCP26_Human.RData")

          #Human_invertebrates species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_invertebrates$i
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character]
            #sum_sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates[,], na.rm = T)
            sum_sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates<-sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates
            energy_RCP26_Human_invertebrates<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_invertebrates
            save(energy_RCP26_Human_invertebrates, file="energy_RCP26_Human_invertebrates.RData")
          
          #Human_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_reproductive_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Human_reproductive_plant_material<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Human_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_reproductive_plant_material[,], na.rm = T)
            energy_RCP26_Human_rep_plant<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_reproductive_plant_material
            save(energy_RCP26_Human_rep_plant, file="energy_RCP26_Human_rep_plant.RData")

          #Human_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vegetative_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Human_vegetative_plant_material<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Human_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_vegetative_plant_material[,], na.rm = T)
            energy_RCP26_Human_veg_plant<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_vegetative_plant_material
            save(energy_RCP26_Human_veg_plant, file="energy_RCP26_Human_veg_plant.RData")

          #Human_vertebrates species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_vertebrates$i
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Human_vertebrates<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Human_vertebrates<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_vertebrates[,], na.rm = T)
            energy_RCP26_Human_vertebrates<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_vertebrates
            save(energy_RCP26_Human_vertebrates, file="energy_RCP26_Human_vertebrates.RData")

          #Human_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$i
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Human_unknown_plant_material_and_others<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Human_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Human_unknown_plant_material_and_others[,], na.rm = T)
            energy_RCP26_Human_unk_plant_oth<-sum_sub_energy_Europe_RCP26_option_13_1_1_Human_unknown_plant_material_and_others
            save(energy_RCP26_Human_unk_plant_oth, file="energy_RCP26_Human_unk_plant_oth.RData")

        #6.2.2.1.3 Wild species
          merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet, Human_origin==0)
          merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild[!merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Wild<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild$i
          vector_sub_merged_sps_list_diet_category3_Wild_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP26_option_13_1_1_Wild<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_character]
          sum_sub_energy_Europe_RCP26_option_13_1_1_Wild<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild[,], na.rm = T)
          energy_RCP26_Wild<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild
          save(energy_RCP26_Wild, file="energy_RCP26_Wild.RData")

          #Wild_invertebrates species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_invertebrates$i
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Wild_invertebrates<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_invertebrates<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild_invertebrates[,], na.rm = T)
            energy_RCP26_Wild_invertebrates<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_invertebrates
            save(energy_RCP26_Wild_invertebrates, file="energy_RCP26_Wild_invertebrates.RData")
          
          #Wild_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_reproductive_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Wild_reproductive_plant_material<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild_reproductive_plant_material[,], na.rm = T)
            energy_RCP26_Wild_rep_plant<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_reproductive_plant_material
            save(energy_RCP26_Wild_rep_plant, file="energy_RCP26_Wild_rep_plant.RData")

          #Wild_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vegetative_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Wild_vegetative_plant_material<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild_vegetative_plant_material[,], na.rm = T)
            energy_RCP26_Wild_veg_plant<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_vegetative_plant_material
            save(energy_RCP26_Wild_veg_plant, file="energy_RCP26_Wild_veg_plant.RData")

          #Wild_vertebrates species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_vertebrates$i
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Wild_vertebrates<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_vertebrates<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild_vertebrates[,], na.rm = T)
            energy_RCP26_Wild_vertebrates<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_vertebrates
            save(energy_RCP26_Wild_vertebrates, file="energy_RCP26_Wild_vertebrates.RData")

          #Wild_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP26_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$i
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP26_option_13_1_1_Wild_unknown_plant_material_and_others<- sub_energy_Europe_RCP26_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP26_option_13_1_1_Wild_unknown_plant_material_and_others[,], na.rm = T)
            energy_RCP26_Wild_unk_plant_oth<-sum_sub_energy_Europe_RCP26_option_13_1_1_Wild_unknown_plant_material_and_others
            save(energy_RCP26_Wild_unk_plant_oth, file="energy_RCP26_Wild_unk_plant_oth.RData")

      #We sum the energy for all species in each pixel:
        DF_Variables_Energy_RCP26<- as.data.frame(cbind(sub_energy_Europe_RCP26_option_13_1_1[,1],
          energy_RCP26_All,
          energy_RCP26_Human,
          energy_RCP26_Human_invertebrates,
          energy_RCP26_Human_rep_plant,
          energy_RCP26_Human_veg_plant,
          energy_RCP26_Human_vertebrates,
          energy_RCP26_Human_unk_plant_oth,
          energy_RCP26_Wild,
          energy_RCP26_Wild_invertebrates,
          energy_RCP26_Wild_rep_plant,
          energy_RCP26_Wild_veg_plant,
          energy_RCP26_Wild_vertebrates,
          energy_RCP26_Wild_unk_plant_oth))
        colnames(DF_Variables_Energy_RCP26)<-c("ID_pixel",
          "ene_All",
          "ene_Human",
          "ene_Human_invertebrates",
          "ene_Human_rep_plant",
          "ene_Human_veg_plant",
          "ene_Human_vertebrates",
          "ene_Human_unk_plant_oth",
          "ene_Wild",
          "ene_Wild_invertebrates",
          "ene_Wild_rep_plant",
          "ene_Wild_veg_plant",
          "ene_Wild_vertebrates",
          "ene_Wild_unk_plant_oth")
        save(DF_Variables_Energy_RCP26, file="DF_Variables_Energy_RCP26.RData")

        #WE do the merge
          load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
          merged_DF_Variables_Energy_RCP26<-merge(DF_ID_pixel, DF_Variables_Energy_RCP26, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          save(merged_DF_Variables_Energy_RCP26, file="merged_DF_Variables_Energy_RCP26.RData")
          
        #We close and open R in 3.5.3 for write rasters
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
          setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
          load(file="merged_DF_Variables_Energy_RCP26.RData")

          #We save into a raster file to see in Arcgis
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          
            #For ene_All
              ras_RCP26_ene_All <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_All) <- newproj
              values(ras_RCP26_ene_All)<-merged_DF_Variables_Energy_RCP26$ene_All
              writeRaster(ras_RCP26_ene_All, "ras_RCP26_ene_All.img", overwrite=TRUE)
  
            #For ene_Human
              ras_RCP26_ene_Human <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human) <- newproj
              values(ras_RCP26_ene_Human)<-merged_DF_Variables_Energy_RCP26$ene_Human
              writeRaster(ras_RCP26_ene_Human, "ras_RCP26_ene_Human.img", overwrite=TRUE)

            #For ene_Human_invertebrates
              ras_RCP26_ene_Human_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human_invertebrates) <- newproj
              values(ras_RCP26_ene_Human_invertebrates)<-merged_DF_Variables_Energy_RCP26$ene_Human_invertebrates
              writeRaster(ras_RCP26_ene_Human_invertebrates, "ras_RCP26_ene_Human_invertebrates.img", overwrite=TRUE)
              
             #For ene_Human_rep_plant
              ras_RCP26_ene_Human_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human_rep_plant) <- newproj
              values(ras_RCP26_ene_Human_rep_plant)<-merged_DF_Variables_Energy_RCP26$ene_Human_rep_plant
              writeRaster(ras_RCP26_ene_Human_rep_plant, "ras_RCP26_ene_Human_rep_plant.img", overwrite=TRUE)
              
            #For ene_Human_veg_plant
              ras_RCP26_ene_Human_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human_veg_plant) <- newproj
              values(ras_RCP26_ene_Human_veg_plant)<-merged_DF_Variables_Energy_RCP26$ene_Human_veg_plant
              writeRaster(ras_RCP26_ene_Human_veg_plant, "ras_RCP26_ene_Human_veg_plant.img", overwrite=TRUE)
             
             #For ene_Human_vertebrates
              ras_RCP26_ene_Human_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human_vertebrates) <- newproj
              values(ras_RCP26_ene_Human_vertebrates)<-merged_DF_Variables_Energy_RCP26$ene_Human_vertebrates
              writeRaster(ras_RCP26_ene_Human_vertebrates, "ras_RCP26_ene_Human_vertebrates.img", overwrite=TRUE)
              
             #For ene_Human_unk_plant_oth
              ras_RCP26_ene_Human_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Human_unk_plant_oth) <- newproj
              values(ras_RCP26_ene_Human_unk_plant_oth)<-merged_DF_Variables_Energy_RCP26$ene_Human_unk_plant_oth
              writeRaster(ras_RCP26_ene_Human_unk_plant_oth, "ras_RCP26_ene_Human_unk_plant_oth.img", overwrite=TRUE)
              
            #For ene_Wild
              ras_RCP26_ene_Wild <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild) <- newproj
              values(ras_RCP26_ene_Wild)<-merged_DF_Variables_Energy_RCP26$ene_Wild
              writeRaster(ras_RCP26_ene_Wild, "ras_RCP26_ene_Wild.img", overwrite=TRUE)

            #For ene_Wild_invertebrates
              ras_RCP26_ene_Wild_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild_invertebrates) <- newproj
              values(ras_RCP26_ene_Wild_invertebrates)<-merged_DF_Variables_Energy_RCP26$ene_Wild_invertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP26_ene_Wild_invertebrates, "ras_RCP26_ene_Wild_invertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_rep_plant
              ras_RCP26_ene_Wild_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild_rep_plant) <- newproj
              values(ras_RCP26_ene_Wild_rep_plant)<-merged_DF_Variables_Energy_RCP26$ene_Wild_rep_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP26_ene_Wild_rep_plant, "ras_RCP26_ene_Wild_rep_plant.img", overwrite=TRUE)
              
            #For ene_Wild_veg_plant
              ras_RCP26_ene_Wild_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild_veg_plant) <- newproj
              values(ras_RCP26_ene_Wild_veg_plant)<-merged_DF_Variables_Energy_RCP26$ene_Wild_veg_plant
              writeRaster(ras_RCP26_ene_Wild_veg_plant, "ras_RCP26_ene_Wild_veg_plant.img", overwrite=TRUE)

             #For ene_Wild_vertebrates
              ras_RCP26_ene_Wild_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild_vertebrates) <- newproj
              values(ras_RCP26_ene_Wild_vertebrates)<-merged_DF_Variables_Energy_RCP26$ene_Wild_vertebrates
              writeRaster(ras_RCP26_ene_Wild_vertebrates, "ras_RCP26_ene_Wild_vertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_unk_plant_oth
              ras_RCP26_ene_Wild_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP26_ene_Wild_unk_plant_oth) <- newproj
              values(ras_RCP26_ene_Wild_unk_plant_oth)<-merged_DF_Variables_Energy_RCP26$ene_Wild_unk_plant_oth
              writeRaster(ras_RCP26_ene_Wild_unk_plant_oth, "ras_RCP26_ene_Wild_unk_plant_oth.img", overwrite=TRUE)

  ##########################################################################################################            
  #6.2.3 For RCP60 scenario 
  ##########################################################################################################            
    rm(list=ls())
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    library(readxl)
    #We load the data frame with the energy by pixel of each species
    load("H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODsub_energy_Europe_RCP60_option_13_1_1.RData")
    load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #We change the value of Human/Wild for the reindeer
    merged_sps_list_diet_category3[214,]
    merged_sps_list_diet_category3[214,4] <- 0
    merged_sps_list_diet_category3[214,]
    colnames_sub_energy_Europe_RCP60_option_13_1_1<-colnames(sub_energy_Europe_RCP60_option_13_1_1[,2:237])
    colnames_sub_energy_Europe_RCP60_option_13_1_1<-as.data.frame(as.numeric(colnames_sub_energy_Europe_RCP60_option_13_1_1))
    colnames(colnames_sub_energy_Europe_RCP60_option_13_1_1)<-c("i")
    merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet<-merge(colnames_sub_energy_Europe_RCP60_option_13_1_1, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)

      #6.2.3.1 We are going to do subset of species based in human/wild origin and in categories of diet
        #6.2.3.1.1 Human and wild species (all species excluding the brown bear)
          merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_All<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet[!merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_All$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_All$Human_origin)
          vector_sub_merged_sps_list_diet_category3_All<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_All$i
          vector_sub_merged_sps_list_diet_category3_All_character<-as.character(vector_sub_merged_sps_list_diet_category3_All)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP60_option_13_1_1_All<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_All_character]
          sum_sub_energy_Europe_RCP60_option_13_1_1_All<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_All[,], na.rm = T)
          energy_RCP60_All<-sum_sub_energy_Europe_RCP60_option_13_1_1_All
          save(energy_RCP60_All, file="energy_RCP60_All.RData")
  
        #6.2.3.1.2 Human species
          merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet, Human_origin==1)
          table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human$i
          vector_sub_merged_sps_list_diet_category3_Human_origin_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_origin)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP60_option_13_1_1_Human_origin<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_origin_character]
          sum_sub_energy_Europe_RCP60_option_13_1_1_Human_origin<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Human_origin[,], na.rm = T)
          energy_RCP60_Human<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_origin
          save(energy_RCP60_Human, file="energy_RCP60_Human.RData")

          #Human_invertebrates species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Human_invertebrates<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Human_invertebrates<-sub_energy_Europe_RCP60_option_13_1_1_Human_invertebrates
            energy_RCP60_Human_invertebrates<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_invertebrates
            save(energy_RCP60_Human_invertebrates, file="energy_RCP60_Human_invertebrates.RData")
          
          #Human_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_reproductive_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Human_reproductive_plant_material<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Human_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Human_reproductive_plant_material[,], na.rm = T)
            energy_RCP60_Human_rep_plant<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_reproductive_plant_material
            save(energy_RCP60_Human_rep_plant, file="energy_RCP60_Human_rep_plant.RData")

          #Human_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vegetative_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Human_vegetative_plant_material<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Human_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Human_vegetative_plant_material[,], na.rm = T)
            energy_RCP60_Human_veg_plant<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_vegetative_plant_material
            save(energy_RCP60_Human_veg_plant, file="energy_RCP60_Human_veg_plant.RData")

          #Human_vertebrates species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_vertebrates$i
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Human_vertebrates<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Human_vertebrates<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Human_vertebrates[,], na.rm = T)
            energy_RCP60_Human_vertebrates<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_vertebrates
            save(energy_RCP60_Human_vertebrates, file="energy_RCP60_Human_vertebrates.RData")

          #Human_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$i
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Human_unknown_plant_material_and_others<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Human_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Human_unknown_plant_material_and_others[,], na.rm = T)
            energy_RCP60_Human_unk_plant_oth<-sum_sub_energy_Europe_RCP60_option_13_1_1_Human_unknown_plant_material_and_others
            save(energy_RCP60_Human_unk_plant_oth, file="energy_RCP60_Human_unk_plant_oth.RData")

        ##6.2.3.1.3 Wild species
          merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet, Human_origin==0)
          merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild[!merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Wild<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild$i
          vector_sub_merged_sps_list_diet_category3_Wild_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP60_option_13_1_1_Wild<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_character]
          sum_sub_energy_Europe_RCP60_option_13_1_1_Wild<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild[,], na.rm = T)
          energy_RCP60_Wild<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild
          save(energy_RCP60_Wild, file="energy_RCP60_Wild.RData")

          #Wild_invertebrates species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_invertebrates$i
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Wild_invertebrates<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_invertebrates<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild_invertebrates[,], na.rm = T)
            energy_RCP60_Wild_invertebrates<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_invertebrates
            save(energy_RCP60_Wild_invertebrates, file="energy_RCP60_Wild_invertebrates.RData")
          
          #Wild_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_reproductive_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Wild_reproductive_plant_material<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild_reproductive_plant_material[,], na.rm = T)
            energy_RCP60_Wild_rep_plant<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_reproductive_plant_material
            save(energy_RCP60_Wild_rep_plant, file="energy_RCP60_Wild_rep_plant.RData")

          #Wild_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vegetative_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Wild_vegetative_plant_material<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild_vegetative_plant_material[,], na.rm = T)
            energy_RCP60_Wild_veg_plant<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_vegetative_plant_material
            save(energy_RCP60_Wild_veg_plant, file="energy_RCP60_Wild_veg_plant.RData")
          #Wild_vertebrates species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_vertebrates$i
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Wild_vertebrates<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_vertebrates<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild_vertebrates[,], na.rm = T)
            energy_RCP60_Wild_vertebrates<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_vertebrates
            save(energy_RCP60_Wild_vertebrates, file="energy_RCP60_Wild_vertebrates.RData")

          #Wild_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP60_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$i
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP60_option_13_1_1_Wild_unknown_plant_material_and_others<- sub_energy_Europe_RCP60_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP60_option_13_1_1_Wild_unknown_plant_material_and_others[,], na.rm = T)
            energy_RCP60_Wild_unk_plant_oth<-sum_sub_energy_Europe_RCP60_option_13_1_1_Wild_unknown_plant_material_and_others
            save(energy_RCP60_Wild_unk_plant_oth, file="energy_RCP60_Wild_unk_plant_oth.RData")

      #We sum the energy for all species in each pixel:
        DF_Variables_Energy_RCP60<- as.data.frame(cbind(sub_energy_Europe_RCP60_option_13_1_1[,1],
          energy_RCP60_All,
          energy_RCP60_Human,
          energy_RCP60_Human_invertebrates,
          energy_RCP60_Human_rep_plant,
          energy_RCP60_Human_veg_plant,
          energy_RCP60_Human_vertebrates,
          energy_RCP60_Human_unk_plant_oth,
          energy_RCP60_Wild,
          energy_RCP60_Wild_invertebrates,
          energy_RCP60_Wild_rep_plant,
          energy_RCP60_Wild_veg_plant,
          energy_RCP60_Wild_vertebrates,
          energy_RCP60_Wild_unk_plant_oth))
        colnames(DF_Variables_Energy_RCP60)<-c("ID_pixel",
          "ene_All",
          "ene_Human",
          "ene_Human_invertebrates",
          "ene_Human_rep_plant",
          "ene_Human_veg_plant",
          "ene_Human_vertebrates",
          "ene_Human_unk_plant_oth",
          "ene_Wild",
          "ene_Wild_invertebrates",
          "ene_Wild_rep_plant",
          "ene_Wild_veg_plant",
          "ene_Wild_vertebrates",
          "ene_Wild_unk_plant_oth")
        save(DF_Variables_Energy_RCP60, file="DF_Variables_Energy_RCP60.RData")

        #WE do the merge
          load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
          merged_DF_Variables_Energy_RCP60<-merge(DF_ID_pixel, DF_Variables_Energy_RCP60, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_DF_Variables_Energy_RCP60)
          head(merged_DF_Variables_Energy_RCP60)
          save(merged_DF_Variables_Energy_RCP60, file="merged_DF_Variables_Energy_RCP60.RData")

        #We close and open R in 3.5.3 for write rasters
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
          setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
          load(file="merged_DF_Variables_Energy_RCP60.RData")

          #We save into a raster file to see in Arcgis
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          
            #For ene_All
              ras_RCP60_ene_All <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_All) <- newproj
              values(ras_RCP60_ene_All)<-merged_DF_Variables_Energy_RCP60$ene_All
              writeRaster(ras_RCP60_ene_All, "ras_RCP60_ene_All.img", overwrite=TRUE)
  
            #For ene_Human
              ras_RCP60_ene_Human <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human) <- newproj
              values(ras_RCP60_ene_Human)<-merged_DF_Variables_Energy_RCP60$ene_Human
              writeRaster(ras_RCP60_ene_Human, "ras_RCP60_ene_Human.img", overwrite=TRUE)

            #For ene_Human_invertebrates
              ras_RCP60_ene_Human_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human_invertebrates) <- newproj
              values(ras_RCP60_ene_Human_invertebrates)<-merged_DF_Variables_Energy_RCP60$ene_Human_invertebrates
              writeRaster(ras_RCP60_ene_Human_invertebrates, "ras_RCP60_ene_Human_invertebrates.img", overwrite=TRUE)
              
             #For ene_Human_rep_plant
              ras_RCP60_ene_Human_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human_rep_plant) <- newproj
              values(ras_RCP60_ene_Human_rep_plant)<-merged_DF_Variables_Energy_RCP60$ene_Human_rep_plant
              writeRaster(ras_RCP60_ene_Human_rep_plant, "ras_RCP60_ene_Human_rep_plant.img", overwrite=TRUE)
              
            #For ene_Human_veg_plant
              ras_RCP60_ene_Human_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human_veg_plant) <- newproj
              values(ras_RCP60_ene_Human_veg_plant)<-merged_DF_Variables_Energy_RCP60$ene_Human_veg_plant
              writeRaster(ras_RCP60_ene_Human_veg_plant, "ras_RCP60_ene_Human_veg_plant.img", overwrite=TRUE)

             #For ene_Human_vertebrates
              ras_RCP60_ene_Human_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human_vertebrates) <- newproj
              values(ras_RCP60_ene_Human_vertebrates)<-merged_DF_Variables_Energy_RCP60$ene_Human_vertebrates
              writeRaster(ras_RCP60_ene_Human_vertebrates, "ras_RCP60_ene_Human_vertebrates.img", overwrite=TRUE)
              
             #For ene_Human_unk_plant_oth
              ras_RCP60_ene_Human_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Human_unk_plant_oth) <- newproj
              values(ras_RCP60_ene_Human_unk_plant_oth)<-merged_DF_Variables_Energy_RCP60$ene_Human_unk_plant_oth
              writeRaster(ras_RCP60_ene_Human_unk_plant_oth, "ras_RCP60_ene_Human_unk_plant_oth.img", overwrite=TRUE)
              
            #For ene_Wild
              ras_RCP60_ene_Wild <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild) <- newproj
              values(ras_RCP60_ene_Wild)<-merged_DF_Variables_Energy_RCP60$ene_Wild
              writeRaster(ras_RCP60_ene_Wild, "ras_RCP60_ene_Wild.img", overwrite=TRUE)

            #For ene_Wild_invertebrates
              ras_RCP60_ene_Wild_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild_invertebrates) <- newproj
              values(ras_RCP60_ene_Wild_invertebrates)<-merged_DF_Variables_Energy_RCP60$ene_Wild_invertebrates
              writeRaster(ras_RCP60_ene_Wild_invertebrates, "ras_RCP60_ene_Wild_invertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_rep_plant
              ras_RCP60_ene_Wild_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild_rep_plant) <- newproj
              values(ras_RCP60_ene_Wild_rep_plant)<-merged_DF_Variables_Energy_RCP60$ene_Wild_rep_plant
              writeRaster(ras_RCP60_ene_Wild_rep_plant, "ras_RCP60_ene_Wild_rep_plant.img", overwrite=TRUE)
              
            #For ene_Wild_veg_plant
              ras_RCP60_ene_Wild_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild_veg_plant) <- newproj
              values(ras_RCP60_ene_Wild_veg_plant)<-merged_DF_Variables_Energy_RCP60$ene_Wild_veg_plant
              writeRaster(ras_RCP60_ene_Wild_veg_plant, "ras_RCP60_ene_Wild_veg_plant.img", overwrite=TRUE)

             #For ene_Wild_vertebrates
              ras_RCP60_ene_Wild_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild_vertebrates) <- newproj
              values(ras_RCP60_ene_Wild_vertebrates)<-merged_DF_Variables_Energy_RCP60$ene_Wild_vertebrates
              writeRaster(ras_RCP60_ene_Wild_vertebrates, "ras_RCP60_ene_Wild_vertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_unk_plant_oth
              ras_RCP60_ene_Wild_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP60_ene_Wild_unk_plant_oth) <- newproj
              values(ras_RCP60_ene_Wild_unk_plant_oth)<-merged_DF_Variables_Energy_RCP60$ene_Wild_unk_plant_oth
              writeRaster(ras_RCP60_ene_Wild_unk_plant_oth, "ras_RCP60_ene_Wild_unk_plant_oth.img", overwrite=TRUE)
               
                  
  ##########################################################################################################            
  #6.2.4 For RCP85 scenario
  ##########################################################################################################            
    rm(list=ls())
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    library(readxl)
    #We load the data frame with the energy by pixel of each species
    load("H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODsub_energy_Europe_RCP85_option_13_1_1.RData")
    load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #We change the value of Human/Wild for the reindeer
    merged_sps_list_diet_category3[214,]
    merged_sps_list_diet_category3[214,4] <- 0
    merged_sps_list_diet_category3[214,]
    colnames_sub_energy_Europe_RCP85_option_13_1_1<-colnames(sub_energy_Europe_RCP85_option_13_1_1[,2:236])
    colnames_sub_energy_Europe_RCP85_option_13_1_1<-as.data.frame(as.numeric(colnames_sub_energy_Europe_RCP85_option_13_1_1))
    colnames(colnames_sub_energy_Europe_RCP85_option_13_1_1)<-c("i")
    merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet<-merge(colnames_sub_energy_Europe_RCP85_option_13_1_1, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)

      #6.2.4.1 We are going to do subset of species based in human/wild origin and in categories of diet
        #6.2.4.1.1 Human and wild species (all species excluding the brown bear)
          merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_All<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet[!merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_All$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_All$Human_origin)
          vector_sub_merged_sps_list_diet_category3_All<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_All$i
          vector_sub_merged_sps_list_diet_category3_All_character<-as.character(vector_sub_merged_sps_list_diet_category3_All)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP85_option_13_1_1_All<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_All_character]
          sum_sub_energy_Europe_RCP85_option_13_1_1_All<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_All[,], na.rm = T)
          energy_RCP85_All<-sum_sub_energy_Europe_RCP85_option_13_1_1_All
          save(energy_RCP85_All, file="energy_RCP85_All.RData")
  
        #6.2.4.1.2 Human species
          merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet, Human_origin==1)
          table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human$i
          vector_sub_merged_sps_list_diet_category3_Human_origin_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_origin)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP85_option_13_1_1_Human_origin<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_origin_character]
          sum_sub_energy_Europe_RCP85_option_13_1_1_Human_origin<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Human_origin[,], na.rm = T)
          energy_RCP85_Human<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_origin
          save(energy_RCP85_Human, file="energy_RCP85_Human.RData")

          #Human_invertebrates species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_invertebrates$i
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Human_invertebrates<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Human_invertebrates<-sub_energy_Europe_RCP85_option_13_1_1_Human_invertebrates
            energy_RCP85_Human_invertebrates<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_invertebrates
            save(energy_RCP85_Human_invertebrates, file="energy_RCP85_Human_invertebrates.RData")
          
          #Human_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_reproductive_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Human_reproductive_plant_material<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Human_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Human_reproductive_plant_material[,], na.rm = T)
            energy_RCP85_Human_rep_plant<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_reproductive_plant_material
            save(energy_RCP85_Human_rep_plant, file="energy_RCP85_Human_rep_plant.RData")

          #Human_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vegetative_plant_material$i
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Human_vegetative_plant_material<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Human_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Human_vegetative_plant_material[,], na.rm = T)
            energy_RCP85_Human_veg_plant<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_vegetative_plant_material
            save(energy_RCP85_Human_veg_plant, file="energy_RCP85_Human_veg_plant.RData")

          #Human_vertebrates species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_vertebrates$i
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Human_vertebrates<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Human_vertebrates<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Human_vertebrates[,], na.rm = T)
            energy_RCP85_Human_vertebrates<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_vertebrates
            save(energy_RCP85_Human_vertebrates, file="energy_RCP85_Human_vertebrates.RData")

          #Human_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Human_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others)
            energy_RCP85_Human_unk_plant_oth<-sum_sub_energy_Europe_RCP85_option_13_1_1_Human_unknown_plant_material_and_others
            save(energy_RCP85_Human_unk_plant_oth, file="energy_RCP85_Human_unk_plant_oth.RData")

        #6.2.4.1.3 Wild species
          merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet, Human_origin==0)
          merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild[!merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Wild<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild$i
          str(vector_sub_merged_sps_list_diet_category3_Wild)
          vector_sub_merged_sps_list_diet_category3_Wild_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_RCP85_option_13_1_1_Wild<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_character]
          sum_sub_energy_Europe_RCP85_option_13_1_1_Wild<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild[,], na.rm = T)
          str(sub_energy_Europe_RCP85_option_13_1_1_Wild)
          head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild)
          energy_RCP85_Wild<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild
          save(energy_RCP85_Wild, file="energy_RCP85_Wild.RData")

          #Wild_invertebrates species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_invertebrates<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates[,], na.rm = T)
            #sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates<-sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates
            str(sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates)
            energy_RCP85_Wild_invertebrates<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_invertebrates
            save(energy_RCP85_Wild_invertebrates, file="energy_RCP85_Wild_invertebrates.RData")
          
          #Wild_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_reproductive_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material)
            energy_RCP85_Wild_rep_plant<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_reproductive_plant_material
            save(energy_RCP85_Wild_rep_plant, file="energy_RCP85_Wild_rep_plant.RData")

          #Wild_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vegetative_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material)
            energy_RCP85_Wild_veg_plant<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vegetative_plant_material
            save(energy_RCP85_Wild_veg_plant, file="energy_RCP85_Wild_veg_plant.RData")

          #Wild_vertebrates species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vertebrates<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_vertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates[,], na.rm = T)
            str(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates)
            energy_RCP85_Wild_vertebrates<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_vertebrates
            save(energy_RCP85_Wild_vertebrates, file="energy_RCP85_Wild_vertebrates.RData")

          #Wild_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_RCP85_option_13_1_1_categories_diet_Wild_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others<- sub_energy_Europe_RCP85_option_13_1_1[,vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others)
            energy_RCP85_Wild_unk_plant_oth<-sum_sub_energy_Europe_RCP85_option_13_1_1_Wild_unknown_plant_material_and_others
            save(energy_RCP85_Wild_unk_plant_oth, file="energy_RCP85_Wild_unk_plant_oth.RData")


      #We sum the energy for all species in each pixel:
        DF_Variables_Energy_RCP85<- as.data.frame(cbind(sub_energy_Europe_RCP85_option_13_1_1[,1],
          energy_RCP85_All,
          energy_RCP85_Human,
          energy_RCP85_Human_invertebrates,
          energy_RCP85_Human_rep_plant,
          energy_RCP85_Human_veg_plant,
          energy_RCP85_Human_vertebrates,
          energy_RCP85_Human_unk_plant_oth,
          energy_RCP85_Wild,
          energy_RCP85_Wild_invertebrates,
          energy_RCP85_Wild_rep_plant,
          energy_RCP85_Wild_veg_plant,
          energy_RCP85_Wild_vertebrates,
          energy_RCP85_Wild_unk_plant_oth))
        head(DF_Variables_Energy_RCP85)
        colnames(DF_Variables_Energy_RCP85)<-c("ID_pixel",
          "ene_All",
          "ene_Human",
          "ene_Human_invertebrates",
          "ene_Human_rep_plant",
          "ene_Human_veg_plant",
          "ene_Human_vertebrates",
          "ene_Human_unk_plant_oth",
          "ene_Wild",
          "ene_Wild_invertebrates",
          "ene_Wild_rep_plant",
          "ene_Wild_veg_plant",
          "ene_Wild_vertebrates",
          "ene_Wild_unk_plant_oth")
        head(DF_Variables_Energy_RCP85)
        save(DF_Variables_Energy_RCP85, file="DF_Variables_Energy_RCP85.RData")

        #WE do the merge
          load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
          merged_DF_Variables_Energy_RCP85<-merge(DF_ID_pixel, DF_Variables_Energy_RCP85, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_DF_Variables_Energy_RCP85)
          head(merged_DF_Variables_Energy_RCP85)
          save(merged_DF_Variables_Energy_RCP85, file="merged_DF_Variables_Energy_RCP85.RData")

        #We close and open R in 3.5.3 for write rasters
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
          setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
          load(file="merged_DF_Variables_Energy_RCP85.RData")

          #We save into a raster file to see in Arcgis
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          
            #For ene_All
              ras_RCP85_ene_All <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_All) <- newproj
              values(ras_RCP85_ene_All)<-merged_DF_Variables_Energy_RCP85$ene_All
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_All, "ras_RCP85_ene_All.img", overwrite=TRUE)
  
            #For ene_Human
              ras_RCP85_ene_Human <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human) <- newproj
              values(ras_RCP85_ene_Human)<-merged_DF_Variables_Energy_RCP85$ene_Human
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human, "ras_RCP85_ene_Human.img", overwrite=TRUE)

            #For ene_Human_invertebrates
              ras_RCP85_ene_Human_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human_invertebrates) <- newproj
              values(ras_RCP85_ene_Human_invertebrates)<-merged_DF_Variables_Energy_RCP85$ene_Human_invertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human_invertebrates, "ras_RCP85_ene_Human_invertebrates.img", overwrite=TRUE)
              
             #For ene_Human_rep_plant
              ras_RCP85_ene_Human_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human_rep_plant) <- newproj
              values(ras_RCP85_ene_Human_rep_plant)<-merged_DF_Variables_Energy_RCP85$ene_Human_rep_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human_rep_plant, "ras_RCP85_ene_Human_rep_plant.img", overwrite=TRUE)
              
            #For ene_Human_veg_plant
              ras_RCP85_ene_Human_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human_veg_plant) <- newproj
              values(ras_RCP85_ene_Human_veg_plant)<-merged_DF_Variables_Energy_RCP85$ene_Human_veg_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human_veg_plant, "ras_RCP85_ene_Human_veg_plant.img", overwrite=TRUE)

             #For ene_Human_vertebrates
              ras_RCP85_ene_Human_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human_vertebrates) <- newproj
              values(ras_RCP85_ene_Human_vertebrates)<-merged_DF_Variables_Energy_RCP85$ene_Human_vertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human_vertebrates, "ras_RCP85_ene_Human_vertebrates.img", overwrite=TRUE)
              
             #For ene_Human_unk_plant_oth
              ras_RCP85_ene_Human_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Human_unk_plant_oth) <- newproj
              values(ras_RCP85_ene_Human_unk_plant_oth)<-merged_DF_Variables_Energy_RCP85$ene_Human_unk_plant_oth
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Human_unk_plant_oth, "ras_RCP85_ene_Human_unk_plant_oth.img", overwrite=TRUE)
              
            #For ene_Wild
              ras_RCP85_ene_Wild <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild) <- newproj
              values(ras_RCP85_ene_Wild)<-merged_DF_Variables_Energy_RCP85$ene_Wild
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild, "ras_RCP85_ene_Wild.img", overwrite=TRUE)

            #For ene_Wild_invertebrates
              ras_RCP85_ene_Wild_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild_invertebrates) <- newproj
              values(ras_RCP85_ene_Wild_invertebrates)<-merged_DF_Variables_Energy_RCP85$ene_Wild_invertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild_invertebrates, "ras_RCP85_ene_Wild_invertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_rep_plant
              ras_RCP85_ene_Wild_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild_rep_plant) <- newproj
              values(ras_RCP85_ene_Wild_rep_plant)<-merged_DF_Variables_Energy_RCP85$ene_Wild_rep_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild_rep_plant, "ras_RCP85_ene_Wild_rep_plant.img", overwrite=TRUE)
              
            #For ene_Wild_veg_plant
              ras_RCP85_ene_Wild_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild_veg_plant) <- newproj
              values(ras_RCP85_ene_Wild_veg_plant)<-merged_DF_Variables_Energy_RCP85$ene_Wild_veg_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild_veg_plant, "ras_RCP85_ene_Wild_veg_plant.img", overwrite=TRUE)

             #For ene_Wild_vertebrates
              ras_RCP85_ene_Wild_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild_vertebrates) <- newproj
              values(ras_RCP85_ene_Wild_vertebrates)<-merged_DF_Variables_Energy_RCP85$ene_Wild_vertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild_vertebrates, "ras_RCP85_ene_Wild_vertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_unk_plant_oth
              ras_RCP85_ene_Wild_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_RCP85_ene_Wild_unk_plant_oth) <- newproj
              values(ras_RCP85_ene_Wild_unk_plant_oth)<-merged_DF_Variables_Energy_RCP85$ene_Wild_unk_plant_oth
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_RCP85_ene_Wild_unk_plant_oth, "ras_RCP85_ene_Wild_unk_plant_oth.img", overwrite=TRUE)

                            
#######################################################################################################################
#######################################################################################################################
#6.3 #CHANGE OF BIOTIC VARIABLES/MAPS
#######################################################################################################################
#######################################################################################################################
    
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
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    #load(file="merged_DF_Variables_Energy_RCP85.RData")
              

    #For ras_current_ene_Wild.img   
    ras_current_ene_Wild<-raster("ras_current_ene_Wild.img")          
    ras_RCP26_ene_Wild<-raster("ras_RCP26_ene_Wild.img")          
    ras_RCP60_ene_Wild<-raster("ras_RCP60_ene_Wild.img")          
    ras_RCP85_ene_Wild<-raster("ras_RCP85_ene_Wild.img")          

    Delta_ene_Wild_RCP26_current<-ras_RCP26_ene_Wild-ras_current_ene_Wild
    writeRaster(Delta_ene_Wild_RCP26_current, "Delta_ene_Wild_RCP26_current.img", overwrite=TRUE)

    Delta_ene_Wild_RCP60_current<-ras_RCP60_ene_Wild-ras_current_ene_Wild
    writeRaster(Delta_ene_Wild_RCP60_current, "Delta_ene_Wild_RCP60_current.img", overwrite=TRUE)

    Delta_ene_Wild_RCP85_current<-ras_RCP85_ene_Wild-ras_current_ene_Wild
    writeRaster(Delta_ene_Wild_RCP85_current, "Delta_ene_Wild_RCP85_current.img", overwrite=TRUE)
      
      
    #For ras_current_ene_Wild_unk_plant_oth.img   
    ras_current_ene_Wild_unk_plant_oth<-raster("ras_current_ene_Wild_unk_plant_oth.img")          
    ras_RCP26_ene_Wild_unk_plant_oth<-raster("ras_RCP26_ene_Wild_unk_plant_oth.img")          
    ras_RCP60_ene_Wild_unk_plant_oth<-raster("ras_RCP60_ene_Wild_unk_plant_oth.img")          
    ras_RCP85_ene_Wild_unk_plant_oth<-raster("ras_RCP85_ene_Wild_unk_plant_oth.img")          

    Delta_ene_Wild_unk_plant_oth_RCP26_current<-ras_RCP26_ene_Wild_unk_plant_oth-ras_current_ene_Wild_unk_plant_oth
    writeRaster(Delta_ene_Wild_unk_plant_oth_RCP26_current, "Delta_ene_Wild_unk_plant_oth_RCP26_current.img", overwrite=TRUE)

    Delta_ene_Wild_unk_plant_oth_RCP60_current<-ras_RCP60_ene_Wild_unk_plant_oth-ras_current_ene_Wild_unk_plant_oth
    writeRaster(Delta_ene_Wild_unk_plant_oth_RCP60_current, "Delta_ene_Wild_unk_plant_oth_RCP60_current.img", overwrite=TRUE)

    Delta_ene_Wild_unk_plant_oth_RCP85_current<-ras_RCP85_ene_Wild_unk_plant_oth-ras_current_ene_Wild_unk_plant_oth
    writeRaster(Delta_ene_Wild_unk_plant_oth_RCP85_current, "Delta_ene_Wild_unk_plant_oth_RCP85_current.img", overwrite=TRUE)
      
    #For ras_current_ene_Wild_rep_plant.img   
    ras_current_ene_Wild_rep_plant<-raster("ras_current_ene_Wild_rep_plant.img")          
    ras_RCP26_ene_Wild_rep_plant<-raster("ras_RCP26_ene_Wild_rep_plant.img")          
    ras_RCP60_ene_Wild_rep_plant<-raster("ras_RCP60_ene_Wild_rep_plant.img")          
    ras_RCP85_ene_Wild_rep_plant<-raster("ras_RCP85_ene_Wild_rep_plant.img")          

    Delta_ene_Wild_rep_plant_RCP26_current<-ras_RCP26_ene_Wild_rep_plant-ras_current_ene_Wild_rep_plant
    writeRaster(Delta_ene_Wild_rep_plant_RCP26_current, "Delta_ene_Wild_rep_plant_RCP26_current.img", overwrite=TRUE)

    Delta_ene_Wild_rep_plant_RCP60_current<-ras_RCP60_ene_Wild_rep_plant-ras_current_ene_Wild_rep_plant
    writeRaster(Delta_ene_Wild_rep_plant_RCP60_current, "Delta_ene_Wild_rep_plant_RCP60_current.img", overwrite=TRUE)

    Delta_ene_Wild_rep_plant_RCP85_current<-ras_RCP85_ene_Wild_rep_plant-ras_current_ene_Wild_rep_plant
    writeRaster(Delta_ene_Wild_rep_plant_RCP85_current, "Delta_ene_Wild_rep_plant_RCP85_current.img", overwrite=TRUE)


    #For ras_current_ene_Wild_veg_plant.img   
    ras_current_ene_Wild_veg_plant<-raster("ras_current_ene_Wild_veg_plant.img")          
    ras_RCP26_ene_Wild_veg_plant<-raster("ras_RCP26_ene_Wild_veg_plant.img")          
    ras_RCP60_ene_Wild_veg_plant<-raster("ras_RCP60_ene_Wild_veg_plant.img")          
    ras_RCP85_ene_Wild_veg_plant<-raster("ras_RCP85_ene_Wild_veg_plant.img")          

    Delta_ene_Wild_veg_plant_RCP26_current<-ras_RCP26_ene_Wild_veg_plant-ras_current_ene_Wild_veg_plant
    writeRaster(Delta_ene_Wild_veg_plant_RCP26_current, "Delta_ene_Wild_veg_plant_RCP26_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_veg_plant_RCP60_current<-ras_RCP60_ene_Wild_veg_plant-ras_current_ene_Wild_veg_plant
    writeRaster(Delta_ene_Wild_veg_plant_RCP60_current, "Delta_ene_Wild_veg_plant_RCP60_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_veg_plant_RCP85_current<-ras_RCP85_ene_Wild_veg_plant-ras_current_ene_Wild_veg_plant
    writeRaster(Delta_ene_Wild_veg_plant_RCP85_current, "Delta_ene_Wild_veg_plant_RCP85_current.img", overwrite=TRUE)
    
    #For ras_current_ene_Wild_vertebrates.img   
    ras_current_ene_Wild_vertebrates<-raster("ras_current_ene_Wild_vertebrates.img")          
    ras_RCP26_ene_Wild_vertebrates<-raster("ras_RCP26_ene_Wild_vertebrates.img")          
    ras_RCP60_ene_Wild_vertebrates<-raster("ras_RCP60_ene_Wild_vertebrates.img")          
    ras_RCP85_ene_Wild_vertebrates<-raster("ras_RCP85_ene_Wild_vertebrates.img")          

    Delta_ene_Wild_vertebrates_RCP26_current<-ras_RCP26_ene_Wild_vertebrates-ras_current_ene_Wild_vertebrates
    writeRaster(Delta_ene_Wild_vertebrates_RCP26_current, "Delta_ene_Wild_vertebrates_RCP26_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_vertebrates_RCP60_current<-ras_RCP60_ene_Wild_vertebrates-ras_current_ene_Wild_vertebrates
    writeRaster(Delta_ene_Wild_vertebrates_RCP60_current, "Delta_ene_Wild_vertebrates_RCP60_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_vertebrates_RCP85_current<-ras_RCP85_ene_Wild_vertebrates-ras_current_ene_Wild_vertebrates
    writeRaster(Delta_ene_Wild_vertebrates_RCP85_current, "Delta_ene_Wild_vertebrates_RCP85_current.img", overwrite=TRUE)
      
    #For ras_current_ene_Wild_invertebrates.img   
    ras_current_ene_Wild_invertebrates<-raster("ras_current_ene_Wild_invertebrates.img")          
    ras_RCP26_ene_Wild_invertebrates<-raster("ras_RCP26_ene_Wild_invertebrates.img")          
    ras_RCP60_ene_Wild_invertebrates<-raster("ras_RCP60_ene_Wild_invertebrates.img")          
    ras_RCP85_ene_Wild_invertebrates<-raster("ras_RCP85_ene_Wild_invertebrates.img")          

    Delta_ene_Wild_invertebrates_RCP26_current<-ras_RCP26_ene_Wild_invertebrates-ras_current_ene_Wild_invertebrates
    writeRaster(Delta_ene_Wild_invertebrates_RCP26_current, "Delta_ene_Wild_invertebrates_RCP26_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_invertebrates_RCP60_current<-ras_RCP60_ene_Wild_invertebrates-ras_current_ene_Wild_invertebrates
    writeRaster(Delta_ene_Wild_invertebrates_RCP60_current, "Delta_ene_Wild_invertebrates_RCP60_current.img", overwrite=TRUE)
    
    Delta_ene_Wild_invertebrates_RCP85_current<-ras_RCP85_ene_Wild_invertebrates-ras_current_ene_Wild_invertebrates
    writeRaster(Delta_ene_Wild_invertebrates_RCP85_current, "Delta_ene_Wild_invertebrates_RCP85_current.img", overwrite=TRUE)
  
     
    #For ras_current_ene_Human_vertebrates.img   
    ras_current_ene_Human_vertebrates<-raster("ras_current_ene_Human_vertebrates.img")          
    ras_RCP26_ene_Human_vertebrates<-raster("ras_RCP26_ene_Human_vertebrates.img")          
    ras_RCP60_ene_Human_vertebrates<-raster("ras_RCP60_ene_Human_vertebrates.img")          
    ras_RCP85_ene_Human_vertebrates<-raster("ras_RCP85_ene_Human_vertebrates.img")          

    Delta_ene_Human_vertebrates_RCP26_current<-ras_RCP26_ene_Human_vertebrates-ras_current_ene_Human_vertebrates
    writeRaster(Delta_ene_Human_vertebrates_RCP26_current, "Delta_ene_Human_vertebrates_RCP26_current.img", overwrite=TRUE)
    
    Delta_ene_Human_vertebrates_RCP60_current<-ras_RCP60_ene_Human_vertebrates-ras_current_ene_Human_vertebrates
    writeRaster(Delta_ene_Human_vertebrates_RCP60_current, "Delta_ene_Human_vertebrates_RCP60_current.img", overwrite=TRUE)
    
    Delta_ene_Human_vertebrates_RCP85_current<-ras_RCP85_ene_Human_vertebrates-ras_current_ene_Human_vertebrates
    writeRaster(Delta_ene_Human_vertebrates_RCP85_current, "Delta_ene_Human_vertebrates_RCP85_current.img", overwrite=TRUE)
 
    
  #6.3.1 Biotic delta of change by percentage from current and absolute values summary by subpopulation    
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
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    load("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy/subpopulations_diet.RData")
    load("H:/D/Test_Biomod/vec_subpop_raster.RData")      

    mapa_subpopulations <- shapefile('H:/G/Project_Name/Maps/Extrapolation limits/Subpopulations3.shp')
    mapa_subpopulations_df<-as.data.frame(mapa_subpopulations)
    mapa_subpopulations_df_names<-mapa_subpopulations_df$code
    mapa_subpopulations_df_names2<-mapa_subpopulations_df_names[c(1:5,13,6:12,14)]
    Delta_vector<-c("RCP26","RCP60","RCP85")     
    
    #6.3.1.1 For ras_current_ene_Wild.img   
        ras_current_ene_Wild<-raster("ras_current_ene_Wild.img")          
        ras_RCP26_ene_Wild<-raster("ras_RCP26_ene_Wild.img")          
        ras_RCP60_ene_Wild<-raster("ras_RCP60_ene_Wild.img")          
        ras_RCP85_ene_Wild<-raster("ras_RCP85_ene_Wild.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_current<-values(ras_current_ene_Wild)
        vec_ene_Wild_RCP26<-values(ras_RCP26_ene_Wild)
        vec_ene_Wild_RCP60<-values(ras_RCP60_ene_Wild)
        vec_ene_Wild_RCP85<-values(ras_RCP85_ene_Wild)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_current,vec_ene_Wild_RCP26,vec_ene_Wild_RCP60,vec_ene_Wild_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_RCP26_current_by_subpopulation<- (sum_ene_Wild_RCP26_by_subpopulation - sum_ene_Wild_current_by_subpopulation)/sum_ene_Wild_current_by_subpopulation*100
          Delta_ene_Wild_RCP60_current_by_subpopulation<- (sum_ene_Wild_RCP60_by_subpopulation - sum_ene_Wild_current_by_subpopulation)/sum_ene_Wild_current_by_subpopulation*100
          Delta_ene_Wild_RCP85_current_by_subpopulation<- (sum_ene_Wild_RCP85_by_subpopulation - sum_ene_Wild_current_by_subpopulation)/sum_ene_Wild_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_current_Europe<-sum(vec_ene_Wild_current,na.rm=T)
          sum_ene_Wild_RCP26_Europe<-sum(vec_ene_Wild_RCP26,na.rm=T)
          sum_ene_Wild_RCP60_Europe<-sum(vec_ene_Wild_RCP60,na.rm=T)
          sum_ene_Wild_RCP85_Europe<-sum(vec_ene_Wild_RCP85,na.rm=T)
  
          Delta_ene_Wild_RCP26_current_Europe <- (sum_ene_Wild_RCP26_Europe - sum_ene_Wild_current_Europe)/sum_ene_Wild_current_Europe*100
          Delta_ene_Wild_RCP60_current_Europe <- (sum_ene_Wild_RCP60_Europe - sum_ene_Wild_current_Europe)/sum_ene_Wild_current_Europe*100
          Delta_ene_Wild_RCP85_current_Europe <- (sum_ene_Wild_RCP85_Europe - sum_ene_Wild_current_Europe)/sum_ene_Wild_current_Europe*100
          Delta_ene_Wild_RCPs_current_Europe<-c(Delta_ene_Wild_RCP26_current_Europe,Delta_ene_Wild_RCP60_current_Europe,Delta_ene_Wild_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_RCP26_current_by_subpopulation,Delta_ene_Wild_RCP60_current_by_subpopulation,Delta_ene_Wild_RCP85_current_by_subpopulation))
          Delta_ene_Wild_RCPs_current2<-as.data.frame(cbind("All_species",Delta_ene_Wild_RCPs_current,Delta_ene_Wild_RCPs_current_Europe))        
          Delta_ene_Wild_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_RCPs_current2_3<-(Delta_ene_Wild_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_RCPs_current2_3, file="Delta_ene_Wild_RCPs_current2_3.RData")  
          
        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_current_Europe<-mean(vec_ene_Wild_current, na.rm=TRUE)
          
          
    
    #6.3.1.2 For ras_current_ene_Wild_unk_plant_oth.img   
        ras_current_ene_Wild_unk_plant_oth<-raster("ras_current_ene_Wild_unk_plant_oth.img")          
        ras_RCP26_ene_Wild_unk_plant_oth<-raster("ras_RCP26_ene_Wild_unk_plant_oth.img")          
        ras_RCP60_ene_Wild_unk_plant_oth<-raster("ras_RCP60_ene_Wild_unk_plant_oth.img")          
        ras_RCP85_ene_Wild_unk_plant_oth<-raster("ras_RCP85_ene_Wild_unk_plant_oth.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_unk_plant_oth_current<-values(ras_current_ene_Wild_unk_plant_oth)
        vec_ene_Wild_unk_plant_oth_RCP26<-values(ras_RCP26_ene_Wild_unk_plant_oth)
        vec_ene_Wild_unk_plant_oth_RCP60<-values(ras_RCP60_ene_Wild_unk_plant_oth)
        vec_ene_Wild_unk_plant_oth_RCP85<-values(ras_RCP85_ene_Wild_unk_plant_oth)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_unk_plant_oth_current,vec_ene_Wild_unk_plant_oth_RCP26,vec_ene_Wild_unk_plant_oth_RCP60,vec_ene_Wild_unk_plant_oth_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_unk_plant_oth_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_unk_plant_oth_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_unk_plant_oth_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_unk_plant_oth_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_unk_plant_oth_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_unk_plant_oth_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_unk_plant_oth_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_unk_plant_oth_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_unk_plant_oth_RCP26_current_by_subpopulation<- (sum_ene_Wild_unk_plant_oth_RCP26_by_subpopulation - sum_ene_Wild_unk_plant_oth_current_by_subpopulation)/sum_ene_Wild_unk_plant_oth_current_by_subpopulation*100
          Delta_ene_Wild_unk_plant_oth_RCP60_current_by_subpopulation<- (sum_ene_Wild_unk_plant_oth_RCP60_by_subpopulation - sum_ene_Wild_unk_plant_oth_current_by_subpopulation)/sum_ene_Wild_unk_plant_oth_current_by_subpopulation*100
          Delta_ene_Wild_unk_plant_oth_RCP85_current_by_subpopulation<- (sum_ene_Wild_unk_plant_oth_RCP85_by_subpopulation - sum_ene_Wild_unk_plant_oth_current_by_subpopulation)/sum_ene_Wild_unk_plant_oth_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_unk_plant_oth_current_Europe<-sum(vec_ene_Wild_unk_plant_oth_current,na.rm=T)
          sum_ene_Wild_unk_plant_oth_RCP26_Europe<-sum(vec_ene_Wild_unk_plant_oth_RCP26,na.rm=T)
          sum_ene_Wild_unk_plant_oth_RCP60_Europe<-sum(vec_ene_Wild_unk_plant_oth_RCP60,na.rm=T)
          sum_ene_Wild_unk_plant_oth_RCP85_Europe<-sum(vec_ene_Wild_unk_plant_oth_RCP85,na.rm=T)
  
          Delta_ene_Wild_unk_plant_oth_RCP26_current_Europe <- (sum_ene_Wild_unk_plant_oth_RCP26_Europe - sum_ene_Wild_unk_plant_oth_current_Europe)/sum_ene_Wild_unk_plant_oth_current_Europe*100
          Delta_ene_Wild_unk_plant_oth_RCP60_current_Europe <- (sum_ene_Wild_unk_plant_oth_RCP60_Europe - sum_ene_Wild_unk_plant_oth_current_Europe)/sum_ene_Wild_unk_plant_oth_current_Europe*100
          Delta_ene_Wild_unk_plant_oth_RCP85_current_Europe <- (sum_ene_Wild_unk_plant_oth_RCP85_Europe - sum_ene_Wild_unk_plant_oth_current_Europe)/sum_ene_Wild_unk_plant_oth_current_Europe*100
          Delta_ene_Wild_unk_plant_oth_RCPs_current_Europe<-c(Delta_ene_Wild_unk_plant_oth_RCP26_current_Europe,Delta_ene_Wild_unk_plant_oth_RCP60_current_Europe,Delta_ene_Wild_unk_plant_oth_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_unk_plant_oth_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_unk_plant_oth_RCP26_current_by_subpopulation,Delta_ene_Wild_unk_plant_oth_RCP60_current_by_subpopulation,Delta_ene_Wild_unk_plant_oth_RCP85_current_by_subpopulation))
          Delta_ene_Wild_unk_plant_oth_RCPs_current2<-as.data.frame(cbind("Bio_unknown_plant",Delta_ene_Wild_unk_plant_oth_RCPs_current,Delta_ene_Wild_unk_plant_oth_RCPs_current_Europe))        
          Delta_ene_Wild_unk_plant_oth_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_unk_plant_oth_RCPs_current2_3<-(Delta_ene_Wild_unk_plant_oth_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_unk_plant_oth_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_unk_plant_oth_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_unk_plant_oth_RCPs_current2_3, file="Delta_ene_Wild_unk_plant_oth_RCPs_current2_3.RData")     
    
        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_unk_plant_oth_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_unk_plant_oth_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_unk_plant_oth_current_Europe<-mean(vec_ene_Wild_unk_plant_oth_current, na.rm=TRUE)
          
    #6.3.1.3 For ras_current_ene_Wild_rep_plant.img   
        ras_current_ene_Wild_rep_plant<-raster("ras_current_ene_Wild_rep_plant.img")          
        ras_RCP26_ene_Wild_rep_plant<-raster("ras_RCP26_ene_Wild_rep_plant.img")          
        ras_RCP60_ene_Wild_rep_plant<-raster("ras_RCP60_ene_Wild_rep_plant.img")          
        ras_RCP85_ene_Wild_rep_plant<-raster("ras_RCP85_ene_Wild_rep_plant.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_rep_plant_current<-values(ras_current_ene_Wild_rep_plant)
        vec_ene_Wild_rep_plant_RCP26<-values(ras_RCP26_ene_Wild_rep_plant)
        vec_ene_Wild_rep_plant_RCP60<-values(ras_RCP60_ene_Wild_rep_plant)
        vec_ene_Wild_rep_plant_RCP85<-values(ras_RCP85_ene_Wild_rep_plant)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_rep_plant_current,vec_ene_Wild_rep_plant_RCP26,vec_ene_Wild_rep_plant_RCP60,vec_ene_Wild_rep_plant_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_rep_plant_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_rep_plant_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_rep_plant_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_rep_plant_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_rep_plant_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_rep_plant_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_rep_plant_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_rep_plant_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_rep_plant_RCP26_current_by_subpopulation<- (sum_ene_Wild_rep_plant_RCP26_by_subpopulation - sum_ene_Wild_rep_plant_current_by_subpopulation)/sum_ene_Wild_rep_plant_current_by_subpopulation*100
          Delta_ene_Wild_rep_plant_RCP60_current_by_subpopulation<- (sum_ene_Wild_rep_plant_RCP60_by_subpopulation - sum_ene_Wild_rep_plant_current_by_subpopulation)/sum_ene_Wild_rep_plant_current_by_subpopulation*100
          Delta_ene_Wild_rep_plant_RCP85_current_by_subpopulation<- (sum_ene_Wild_rep_plant_RCP85_by_subpopulation - sum_ene_Wild_rep_plant_current_by_subpopulation)/sum_ene_Wild_rep_plant_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_rep_plant_current_Europe<-sum(vec_ene_Wild_rep_plant_current,na.rm=T)
          sum_ene_Wild_rep_plant_RCP26_Europe<-sum(vec_ene_Wild_rep_plant_RCP26,na.rm=T)
          sum_ene_Wild_rep_plant_RCP60_Europe<-sum(vec_ene_Wild_rep_plant_RCP60,na.rm=T)
          sum_ene_Wild_rep_plant_RCP85_Europe<-sum(vec_ene_Wild_rep_plant_RCP85,na.rm=T)
  
          Delta_ene_Wild_rep_plant_RCP26_current_Europe <- (sum_ene_Wild_rep_plant_RCP26_Europe - sum_ene_Wild_rep_plant_current_Europe)/sum_ene_Wild_rep_plant_current_Europe*100
          Delta_ene_Wild_rep_plant_RCP60_current_Europe <- (sum_ene_Wild_rep_plant_RCP60_Europe - sum_ene_Wild_rep_plant_current_Europe)/sum_ene_Wild_rep_plant_current_Europe*100
          Delta_ene_Wild_rep_plant_RCP85_current_Europe <- (sum_ene_Wild_rep_plant_RCP85_Europe - sum_ene_Wild_rep_plant_current_Europe)/sum_ene_Wild_rep_plant_current_Europe*100
          Delta_ene_Wild_rep_plant_RCPs_current_Europe<-c(Delta_ene_Wild_rep_plant_RCP26_current_Europe,Delta_ene_Wild_rep_plant_RCP60_current_Europe,Delta_ene_Wild_rep_plant_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_rep_plant_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_rep_plant_RCP26_current_by_subpopulation,Delta_ene_Wild_rep_plant_RCP60_current_by_subpopulation,Delta_ene_Wild_rep_plant_RCP85_current_by_subpopulation))
          Delta_ene_Wild_rep_plant_RCPs_current2<-as.data.frame(cbind("Bio_reprod_plant",Delta_ene_Wild_rep_plant_RCPs_current,Delta_ene_Wild_rep_plant_RCPs_current_Europe))        
          Delta_ene_Wild_rep_plant_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_rep_plant_RCPs_current2_3<-(Delta_ene_Wild_rep_plant_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_rep_plant_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_rep_plant_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_rep_plant_RCPs_current2_3, file="Delta_ene_Wild_rep_plant_RCPs_current2_3.RData")     

        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_rep_plant_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_rep_plant_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_rep_plant_current_Europe<-mean(vec_ene_Wild_rep_plant_current, na.rm=TRUE)
              
    #6.3.1.4 For ras_current_ene_Wild_veg_plant.img   
        ras_current_ene_Wild_veg_plant<-raster("ras_current_ene_Wild_veg_plant.img")          
        ras_RCP26_ene_Wild_veg_plant<-raster("ras_RCP26_ene_Wild_veg_plant.img")          
        ras_RCP60_ene_Wild_veg_plant<-raster("ras_RCP60_ene_Wild_veg_plant.img")          
        ras_RCP85_ene_Wild_veg_plant<-raster("ras_RCP85_ene_Wild_veg_plant.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_veg_plant_current<-values(ras_current_ene_Wild_veg_plant)
        vec_ene_Wild_veg_plant_RCP26<-values(ras_RCP26_ene_Wild_veg_plant)
        vec_ene_Wild_veg_plant_RCP60<-values(ras_RCP60_ene_Wild_veg_plant)
        vec_ene_Wild_veg_plant_RCP85<-values(ras_RCP85_ene_Wild_veg_plant)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_veg_plant_current,vec_ene_Wild_veg_plant_RCP26,vec_ene_Wild_veg_plant_RCP60,vec_ene_Wild_veg_plant_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_veg_plant_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_veg_plant_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_veg_plant_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_veg_plant_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_veg_plant_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_veg_plant_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_veg_plant_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_veg_plant_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_veg_plant_RCP26_current_by_subpopulation<- (sum_ene_Wild_veg_plant_RCP26_by_subpopulation - sum_ene_Wild_veg_plant_current_by_subpopulation)/sum_ene_Wild_veg_plant_current_by_subpopulation*100
          Delta_ene_Wild_veg_plant_RCP60_current_by_subpopulation<- (sum_ene_Wild_veg_plant_RCP60_by_subpopulation - sum_ene_Wild_veg_plant_current_by_subpopulation)/sum_ene_Wild_veg_plant_current_by_subpopulation*100
          Delta_ene_Wild_veg_plant_RCP85_current_by_subpopulation<- (sum_ene_Wild_veg_plant_RCP85_by_subpopulation - sum_ene_Wild_veg_plant_current_by_subpopulation)/sum_ene_Wild_veg_plant_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_veg_plant_current_Europe<-sum(vec_ene_Wild_veg_plant_current,na.rm=T)
          sum_ene_Wild_veg_plant_RCP26_Europe<-sum(vec_ene_Wild_veg_plant_RCP26,na.rm=T)
          sum_ene_Wild_veg_plant_RCP60_Europe<-sum(vec_ene_Wild_veg_plant_RCP60,na.rm=T)
          sum_ene_Wild_veg_plant_RCP85_Europe<-sum(vec_ene_Wild_veg_plant_RCP85,na.rm=T)
  
          Delta_ene_Wild_veg_plant_RCP26_current_Europe <- (sum_ene_Wild_veg_plant_RCP26_Europe - sum_ene_Wild_veg_plant_current_Europe)/sum_ene_Wild_veg_plant_current_Europe*100
          Delta_ene_Wild_veg_plant_RCP60_current_Europe <- (sum_ene_Wild_veg_plant_RCP60_Europe - sum_ene_Wild_veg_plant_current_Europe)/sum_ene_Wild_veg_plant_current_Europe*100
          Delta_ene_Wild_veg_plant_RCP85_current_Europe <- (sum_ene_Wild_veg_plant_RCP85_Europe - sum_ene_Wild_veg_plant_current_Europe)/sum_ene_Wild_veg_plant_current_Europe*100
          Delta_ene_Wild_veg_plant_RCPs_current_Europe<-c(Delta_ene_Wild_veg_plant_RCP26_current_Europe,Delta_ene_Wild_veg_plant_RCP60_current_Europe,Delta_ene_Wild_veg_plant_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_veg_plant_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_veg_plant_RCP26_current_by_subpopulation,Delta_ene_Wild_veg_plant_RCP60_current_by_subpopulation,Delta_ene_Wild_veg_plant_RCP85_current_by_subpopulation))
          Delta_ene_Wild_veg_plant_RCPs_current2<-as.data.frame(cbind("Bio_veget_plant",Delta_ene_Wild_veg_plant_RCPs_current,Delta_ene_Wild_veg_plant_RCPs_current_Europe))        
          Delta_ene_Wild_veg_plant_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_veg_plant_RCPs_current2_3<-(Delta_ene_Wild_veg_plant_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_veg_plant_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_veg_plant_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_veg_plant_RCPs_current2_3, file="Delta_ene_Wild_veg_plant_RCPs_current2_3.RData")     
    
        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_veg_plant_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_veg_plant_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_veg_plant_current_Europe<-mean(vec_ene_Wild_veg_plant_current, na.rm=TRUE)
              
    #6.3.1.5 For ras_current_ene_Wild_vertebrates.img   
        ras_current_ene_Wild_vertebrates<-raster("ras_current_ene_Wild_vertebrates.img")          
        ras_RCP26_ene_Wild_vertebrates<-raster("ras_RCP26_ene_Wild_vertebrates.img")          
        ras_RCP60_ene_Wild_vertebrates<-raster("ras_RCP60_ene_Wild_vertebrates.img")          
        ras_RCP85_ene_Wild_vertebrates<-raster("ras_RCP85_ene_Wild_vertebrates.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_vertebrates_current<-values(ras_current_ene_Wild_vertebrates)
        vec_ene_Wild_vertebrates_RCP26<-values(ras_RCP26_ene_Wild_vertebrates)
        vec_ene_Wild_vertebrates_RCP60<-values(ras_RCP60_ene_Wild_vertebrates)
        vec_ene_Wild_vertebrates_RCP85<-values(ras_RCP85_ene_Wild_vertebrates)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_vertebrates_current,vec_ene_Wild_vertebrates_RCP26,vec_ene_Wild_vertebrates_RCP60,vec_ene_Wild_vertebrates_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_vertebrates_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_vertebrates_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_vertebrates_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_vertebrates_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_vertebrates_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_vertebrates_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_vertebrates_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_vertebrates_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_vertebrates_RCP26_current_by_subpopulation<- (sum_ene_Wild_vertebrates_RCP26_by_subpopulation - sum_ene_Wild_vertebrates_current_by_subpopulation)/sum_ene_Wild_vertebrates_current_by_subpopulation*100
          Delta_ene_Wild_vertebrates_RCP60_current_by_subpopulation<- (sum_ene_Wild_vertebrates_RCP60_by_subpopulation - sum_ene_Wild_vertebrates_current_by_subpopulation)/sum_ene_Wild_vertebrates_current_by_subpopulation*100
          Delta_ene_Wild_vertebrates_RCP85_current_by_subpopulation<- (sum_ene_Wild_vertebrates_RCP85_by_subpopulation - sum_ene_Wild_vertebrates_current_by_subpopulation)/sum_ene_Wild_vertebrates_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_vertebrates_current_Europe<-sum(vec_ene_Wild_vertebrates_current,na.rm=T)
          sum_ene_Wild_vertebrates_RCP26_Europe<-sum(vec_ene_Wild_vertebrates_RCP26,na.rm=T)
          sum_ene_Wild_vertebrates_RCP60_Europe<-sum(vec_ene_Wild_vertebrates_RCP60,na.rm=T)
          sum_ene_Wild_vertebrates_RCP85_Europe<-sum(vec_ene_Wild_vertebrates_RCP85,na.rm=T)
  
          Delta_ene_Wild_vertebrates_RCP26_current_Europe <- (sum_ene_Wild_vertebrates_RCP26_Europe - sum_ene_Wild_vertebrates_current_Europe)/sum_ene_Wild_vertebrates_current_Europe*100
          Delta_ene_Wild_vertebrates_RCP60_current_Europe <- (sum_ene_Wild_vertebrates_RCP60_Europe - sum_ene_Wild_vertebrates_current_Europe)/sum_ene_Wild_vertebrates_current_Europe*100
          Delta_ene_Wild_vertebrates_RCP85_current_Europe <- (sum_ene_Wild_vertebrates_RCP85_Europe - sum_ene_Wild_vertebrates_current_Europe)/sum_ene_Wild_vertebrates_current_Europe*100
          Delta_ene_Wild_vertebrates_RCPs_current_Europe<-c(Delta_ene_Wild_vertebrates_RCP26_current_Europe,Delta_ene_Wild_vertebrates_RCP60_current_Europe,Delta_ene_Wild_vertebrates_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_vertebrates_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_vertebrates_RCP26_current_by_subpopulation,Delta_ene_Wild_vertebrates_RCP60_current_by_subpopulation,Delta_ene_Wild_vertebrates_RCP85_current_by_subpopulation))
          Delta_ene_Wild_vertebrates_RCPs_current2<-as.data.frame(cbind("Bio_vertebrates",Delta_ene_Wild_vertebrates_RCPs_current,Delta_ene_Wild_vertebrates_RCPs_current_Europe))        
          Delta_ene_Wild_vertebrates_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_vertebrates_RCPs_current2_3<-(Delta_ene_Wild_vertebrates_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_vertebrates_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_vertebrates_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_vertebrates_RCPs_current2_3, file="Delta_ene_Wild_vertebrates_RCPs_current2_3.RData")     

        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_vertebrates_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_vertebrates_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_vertebrates_current_Europe<-mean(vec_ene_Wild_vertebrates_current, na.rm=TRUE)
              
    #6.3.1.6 For ras_current_ene_Wild_invertebrates_invertebrates.img   
        ras_current_ene_Wild_invertebrates<-raster("ras_current_ene_Wild_invertebrates.img")          
        ras_RCP26_ene_Wild_invertebrates<-raster("ras_RCP26_ene_Wild_invertebrates.img")          
        ras_RCP60_ene_Wild_invertebrates<-raster("ras_RCP60_ene_Wild_invertebrates.img")          
        ras_RCP85_ene_Wild_invertebrates<-raster("ras_RCP85_ene_Wild_invertebrates.img")          
   
       #We extract the values of the raster to a vector:
        vec_ene_Wild_invertebrates_current<-values(ras_current_ene_Wild_invertebrates)
        vec_ene_Wild_invertebrates_RCP26<-values(ras_RCP26_ene_Wild_invertebrates)
        vec_ene_Wild_invertebrates_RCP60<-values(ras_RCP60_ene_Wild_invertebrates)
        vec_ene_Wild_invertebrates_RCP85<-values(ras_RCP85_ene_Wild_invertebrates)
        
        data_sps_in_loop<-as.data.frame(cbind(vec_ene_Wild_invertebrates_current,vec_ene_Wild_invertebrates_RCP26,vec_ene_Wild_invertebrates_RCP60,vec_ene_Wild_invertebrates_RCP85,vec_subpop_raster))

       #We calculate the sums at subpopulation scale 
          sum_ene_Wild_invertebrates_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_invertebrates_current,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_invertebrates_RCP26_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_invertebrates_RCP26,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_invertebrates_RCP60_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_invertebrates_RCP60,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          sum_ene_Wild_invertebrates_RCP85_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_invertebrates_RCP85,data_sps_in_loop$vec_subpop_raster,sum, na.rm=TRUE)
          
          Delta_ene_Wild_invertebrates_RCP26_current_by_subpopulation<- (sum_ene_Wild_invertebrates_RCP26_by_subpopulation - sum_ene_Wild_invertebrates_current_by_subpopulation)/sum_ene_Wild_invertebrates_current_by_subpopulation*100
          Delta_ene_Wild_invertebrates_RCP60_current_by_subpopulation<- (sum_ene_Wild_invertebrates_RCP60_by_subpopulation - sum_ene_Wild_invertebrates_current_by_subpopulation)/sum_ene_Wild_invertebrates_current_by_subpopulation*100
          Delta_ene_Wild_invertebrates_RCP85_current_by_subpopulation<- (sum_ene_Wild_invertebrates_RCP85_by_subpopulation - sum_ene_Wild_invertebrates_current_by_subpopulation)/sum_ene_Wild_invertebrates_current_by_subpopulation*100
  
       #We calculate the sums at Europe scale 
          sum_ene_Wild_invertebrates_current_Europe<-sum(vec_ene_Wild_invertebrates_current,na.rm=T)
          sum_ene_Wild_invertebrates_RCP26_Europe<-sum(vec_ene_Wild_invertebrates_RCP26,na.rm=T)
          sum_ene_Wild_invertebrates_RCP60_Europe<-sum(vec_ene_Wild_invertebrates_RCP60,na.rm=T)
          sum_ene_Wild_invertebrates_RCP85_Europe<-sum(vec_ene_Wild_invertebrates_RCP85,na.rm=T)
  
          Delta_ene_Wild_invertebrates_RCP26_current_Europe <- (sum_ene_Wild_invertebrates_RCP26_Europe - sum_ene_Wild_invertebrates_current_Europe)/sum_ene_Wild_invertebrates_current_Europe*100
          Delta_ene_Wild_invertebrates_RCP60_current_Europe <- (sum_ene_Wild_invertebrates_RCP60_Europe - sum_ene_Wild_invertebrates_current_Europe)/sum_ene_Wild_invertebrates_current_Europe*100
          Delta_ene_Wild_invertebrates_RCP85_current_Europe <- (sum_ene_Wild_invertebrates_RCP85_Europe - sum_ene_Wild_invertebrates_current_Europe)/sum_ene_Wild_invertebrates_current_Europe*100
          Delta_ene_Wild_invertebrates_RCPs_current_Europe<-c(Delta_ene_Wild_invertebrates_RCP26_current_Europe,Delta_ene_Wild_invertebrates_RCP60_current_Europe,Delta_ene_Wild_invertebrates_RCP85_current_Europe)  
        
        #We bind the data at subpopulation and at European scale
          Delta_ene_Wild_invertebrates_RCPs_current<-as.data.frame(rbind(Delta_ene_Wild_invertebrates_RCP26_current_by_subpopulation,Delta_ene_Wild_invertebrates_RCP60_current_by_subpopulation,Delta_ene_Wild_invertebrates_RCP85_current_by_subpopulation))
          Delta_ene_Wild_invertebrates_RCPs_current2<-as.data.frame(cbind("Bio_invertebrates",Delta_ene_Wild_invertebrates_RCPs_current,Delta_ene_Wild_invertebrates_RCPs_current_Europe))        
          Delta_ene_Wild_invertebrates_RCPs_current2$Delta_RCP<-Delta_vector
          Delta_ene_Wild_invertebrates_RCPs_current2_3<-(Delta_ene_Wild_invertebrates_RCPs_current2[,c(1,17,16,2:15)])
          colnames(Delta_ene_Wild_invertebrates_RCPs_current2_3)<-c("Biotic_variable","Delta_RCP","Delta_Biotic_variable_RCPs_current_Europe",mapa_subpopulations_df_names2)
          rownames(Delta_ene_Wild_invertebrates_RCPs_current2_3)<-c()   
          
        #We save the data of delta for each species
          save(Delta_ene_Wild_invertebrates_RCPs_current2_3, file="Delta_ene_Wild_invertebrates_RCPs_current2_3.RData")     
   
        #We calculate the sums at subpopulation scale 
          mean_ene_Wild_invertebrates_current_by_subpopulation<-tapply(data_sps_in_loop$vec_ene_Wild_invertebrates_current,data_sps_in_loop$vec_subpop_raster,mean, na.rm=TRUE)
          mean_ene_Wild_invertebrates_current_Europe<-mean(vec_ene_Wild_invertebrates_current, na.rm=TRUE)
              
    #6.3.1.7 We bind all groups 
          
      #The change of biotic variables    
        Delta_ene_RCPs_current<-rbind(Delta_ene_Wild_RCPs_current2_3,
          Delta_ene_Wild_rep_plant_RCPs_current2_3,
          Delta_ene_Wild_veg_plant_RCPs_current2_3,
          Delta_ene_Wild_unk_plant_oth_RCPs_current2_3,
          Delta_ene_Wild_invertebrates_RCPs_current2_3,
          Delta_ene_Wild_vertebrates_RCPs_current2_3)

        write.csv(Delta_ene_RCPs_current,file="Delta_ene_RCPs_current.csv", row.names = F)
    
        
      #The current values
        mean_ene_current_by_subpopulation<-as.data.frame(rbind(
          mean_ene_Wild_current_by_subpopulation,
          mean_ene_Wild_rep_plant_current_by_subpopulation,
          mean_ene_Wild_veg_plant_current_by_subpopulation,
          mean_ene_Wild_unk_plant_oth_current_by_subpopulation,
          mean_ene_Wild_invertebrates_current_by_subpopulation,
          mean_ene_Wild_vertebrates_current_by_subpopulation
        ))
        
        mean_ene_current_by_subpopulation$Biotic_variable<-c("All_species","Bio_reprod_plant","Bio_veget_plant","Bio_unknown_plant","Bio_invertebrates","Bio_vertebrates")
        mean_ene_current_by_subpopulation2<-(mean_ene_current_by_subpopulation[,c(15,1:14)])
        rownames(mean_ene_current_by_subpopulation2)<-c()
        
        mean_ene_current_by_subpopulation2
        
        
        
        
        mean_ene_current_by_Europe<-rbind(
          mean_ene_Wild_current_Europe,
          mean_ene_Wild_rep_plant_current_Europe,
          mean_ene_Wild_veg_plant_current_Europe,
          mean_ene_Wild_unk_plant_oth_current_Europe,
          mean_ene_Wild_invertebrates_current_Europe,
          mean_ene_Wild_vertebrates_current_Europe
        )
        
        mean_ene_current_by_subpopulation2_Europe<-cbind(mean_ene_current_by_Europe,mean_ene_current_by_subpopulation2)
        mean_ene_current_by_subpopulation2_Europe2<-(mean_ene_current_by_subpopulation2_Europe[,c(2,1,3:16)])
        rownames(mean_ene_current_by_subpopulation2_Europe2)<-c()
        colnames(mean_ene_current_by_subpopulation2_Europe2)<-c("Biotic_variable","Mean_ene_current_for_Europe",mapa_subpopulations_df_names2)
        
        write.csv(mean_ene_current_by_subpopulation2_Europe2,file="mean_ene_current_by_subpopulation2_Europe2.csv", row.names = F)
        
            
#######################################################################################################################
#######################################################################################################################
#6.4 We create subsets of energy by category HOMOGENEOUS OPTION binary option
#######################################################################################################################
#######################################################################################################################
    
  #6.4.1 For Current scenario
    rm(list=ls())
    setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
    library(readxl)
    #We load the data frame with the energy by pixel of each species
    load("H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODsub_energy_Europe_current_option_13_1_1_homogeneous.RData")
    load("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables/merged_sps_list_diet_category3.RData")
    #View(merged_sps_list_diet_category3)
    #We change the value of Human/Wild for the reindeer
    merged_sps_list_diet_category3[214,4] <- 0

    colnames_sub_energy_Europe_current_option_13_1_1_homogeneous<-colnames(sub_energy_Europe_current_option_13_1_1_homogeneous[,2:237])
    colnames_sub_energy_Europe_current_option_13_1_1_homogeneous<-as.data.frame(as.numeric(colnames_sub_energy_Europe_current_option_13_1_1_homogeneous))
    colnames(colnames_sub_energy_Europe_current_option_13_1_1_homogeneous)<-c("i")
    merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet<-merge(colnames_sub_energy_Europe_current_option_13_1_1_homogeneous, merged_sps_list_diet_category3, by.x="i", by.y = "i",all.x = T)
    str(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet)
    
      #6.4.1.1 We are going to do subset of species based in human/wild origin and in categories of diet
        #6.4.1.1.1 Human and wild species (all species excluding the brown bear)
          merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_All<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet[!merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_All$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_All$Human_origin)
          vector_sub_merged_sps_list_diet_category3_All<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_All$i
          str(vector_sub_merged_sps_list_diet_category3_All)
          vector_sub_merged_sps_list_diet_category3_All_character<-as.character(vector_sub_merged_sps_list_diet_category3_All)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_homogeneous_All<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_All_character]
          sum_sub_energy_Europe_current_option_13_1_1_homogeneous_All<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_All[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_homogeneous_All)
          head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_All)
          energy_current_All_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_All
          save(energy_current_All_homogeneous, file="energy_current_All_homogeneous.RData")
  
        #6.4.1.1.2 Human species
          merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet, Human_origin==1)
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human$i
          str(vector_sub_merged_sps_list_diet_category3_Human_origin)
          vector_sub_merged_sps_list_diet_category3_Human_origin_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_origin)
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_origin_character]
          sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin)
          head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin)
          energy_current_Human_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_origin
          save(energy_current_Human_homogeneous, file="energy_current_Human_homogeneous.RData")

          #Human_invertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_invertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_invertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_invertebrates_character]
            #sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates[,], na.rm = T)
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates<-sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates
            str(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates)
            energy_current_Human_invertebrates_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_invertebrates
            save(energy_current_Human_invertebrates_homogeneous, file="energy_current_Human_invertebrates_homogeneous.RData")
          
          #Human_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_reproductive_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_reproductive_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material)
            energy_current_Human_rep_plant_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_reproductive_plant_material
            save(energy_current_Human_rep_plant_homogeneous, file="energy_current_Human_rep_plant_homogeneous.RData")

          #Human_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vegetative_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_vegetative_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material)
            energy_current_Human_veg_plant_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vegetative_plant_material
            save(energy_current_Human_veg_plant_homogeneous, file="energy_current_Human_veg_plant_homogeneous.RData")

          #Human_vertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_vertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_vertebrates)
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_vertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates)
            energy_current_Human_vertebrates_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_vertebrates
            save(energy_current_Human_vertebrates_homogeneous, file="energy_current_Human_vertebrates_homogeneous.RData")

          #Human_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Human_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Human_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others)
            energy_current_Human_unk_plant_oth_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Human_unknown_plant_material_and_others
            save(energy_current_Human_unk_plant_oth_homogeneous, file="energy_current_Human_unk_plant_oth_homogeneous.RData")

        #6.4.1.1.3 Wild species
          merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet, Human_origin==0)
          merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild[!merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild$Species == "Ursus arctos", ]
          table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild$Human_origin)
          vector_sub_merged_sps_list_diet_category3_Wild<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild$i
          str(vector_sub_merged_sps_list_diet_category3_Wild)
          vector_sub_merged_sps_list_diet_category3_Wild_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild)
      
          #We subset the data frame of all pixels in rows and species in the columns
          sub_energy_Europe_current_option_13_1_1_homogeneous_Wild<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_character]
          sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild[,], na.rm = T)
          str(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild)
          head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild)
          energy_current_Wild_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild
          save(energy_current_Wild_homogeneous, file="energy_current_Wild_homogeneous.RData")

          #Wild_invertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_invertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild, Diet_category2=="invertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_invertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_invertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_invertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_invertebrates)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_invertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates[,], na.rm = T)
            #sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates<-sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates
            str(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates)
            energy_current_Wild_invertebrates_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_invertebrates
            save(energy_current_Wild_invertebrates_homogeneous, file="energy_current_Wild_invertebrates_homogeneous.RData")
          
          #Wild_reproductive_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_reproductive_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild, Diet_category2=="reproductive_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_reproductive_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_reproductive_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_reproductive_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_reproductive_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material)
            energy_current_Wild_rep_plant_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_reproductive_plant_material
            save(energy_current_Wild_rep_plant_homogeneous, file="energy_current_Wild_rep_plant_homogeneous.RData")

          #Wild_vegetative_plant_material species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vegetative_plant_material<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild, Diet_category2=="vegetative_plant_material")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vegetative_plant_material$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vegetative_plant_material$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vegetative_plant_material$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
            vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_vegetative_plant_material_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material)
            energy_current_Wild_veg_plant_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vegetative_plant_material
            save(energy_current_Wild_veg_plant_homogeneous, file="energy_current_Wild_veg_plant_homogeneous.RData")

          #Wild_vertebrates species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vertebrates<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild, Diet_category2=="vertebrates")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vertebrates$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vertebrates$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_vertebrates$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
            vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_vertebrates)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_vertebrates_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates)
            energy_current_Wild_vertebrates_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_vertebrates
            save(energy_current_Wild_vertebrates_homogeneous, file="energy_current_Wild_vertebrates_homogeneous.RData")

          #Wild_unknown_plant_material_and_others species
            merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_unknown_plant_material_and_others<-subset(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild, Diet_category2=="unknown_plant_material_and_others")
            table(merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_unknown_plant_material_and_others$Diet_category2,merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_unknown_plant_material_and_others$Human_origin)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others<-merged_colnames_sub_energy_Europe_current_option_13_1_1_homogeneous_categories_diet_Wild_unknown_plant_material_and_others$i
            str(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
            vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character<-as.character(vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others)
        
            #We subset the data frame of all pixels in rows and species in the columns
            sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others<- sub_energy_Europe_current_option_13_1_1_homogeneous[,vector_sub_merged_sps_list_diet_category3_Wild_unknown_plant_material_and_others_character]
            sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others<-rowSums(sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others[,], na.rm = T)
            str(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others)
            head(sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others)
            energy_current_Wild_unk_plant_oth_homogeneous<-sum_sub_energy_Europe_current_option_13_1_1_homogeneous_Wild_unknown_plant_material_and_others
            save(energy_current_Wild_unk_plant_oth_homogeneous, file="energy_current_Wild_unk_plant_oth_homogeneous.RData")


      #We sum the energy for all species in each pixel:
        DF_Variables_Energy_current_homogeneous<- as.data.frame(cbind(sub_energy_Europe_current_option_13_1_1_homogeneous[,1],
          energy_current_All_homogeneous,
          energy_current_Human_homogeneous,
          energy_current_Human_invertebrates_homogeneous,
          energy_current_Human_rep_plant_homogeneous,
          energy_current_Human_veg_plant_homogeneous,
          energy_current_Human_vertebrates_homogeneous,
          energy_current_Human_unk_plant_oth_homogeneous,
          energy_current_Wild_homogeneous,
          energy_current_Wild_invertebrates_homogeneous,
          energy_current_Wild_rep_plant_homogeneous,
          energy_current_Wild_veg_plant_homogeneous,
          energy_current_Wild_vertebrates_homogeneous,
          energy_current_Wild_unk_plant_oth_homogeneous))
        head(DF_Variables_Energy_current_homogeneous)
        colnames(DF_Variables_Energy_current_homogeneous)<-c("ID_pixel",
          "ene_All_h",
          "ene_Human_h",
          "ene_Human_invertebrates_h",
          "ene_Human_rep_plant_h",
          "ene_Human_veg_plant_h",
          "ene_Human_vertebrates_h",
          "ene_Human_unk_plant_oth_h",
          "ene_Wild_h",
          "ene_Wild_invertebrates_h",
          "ene_Wild_rep_plant_h",
          "ene_Wild_veg_plant_h",
          "ene_Wild_vertebrates_h",
          "ene_Wild_unk_plant_oth_h")
        head(DF_Variables_Energy_current_homogeneous)
        save(DF_Variables_Energy_current_homogeneous, file="DF_Variables_Energy_current_homogeneous.RData")

        #WE do the merge
          rm(list=ls())
          setwd("H:/G/Project_Name/Results_Biomod/R_analysis/Energy_variables")
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
           #We load the satellite image for background maps created in part 3  
          #load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
          #load(file=paste0(path_save,"DF_Variables_Energy_Europe_current_option_13_1_1.RData"))
          load(file="H:/G/Project_Name/Results_Biomod/R_analysis/SDM_MODDF_ID_pixel.RData")
          load(file="DF_Variables_Energy_current_homogeneous.RData")
          merged_DF_Variables_Energy_current_homogeneous<-merge(DF_ID_pixel, DF_Variables_Energy_current_homogeneous, by.x="ID_pixel", by.y = "ID_pixel",all.x = T)
          str(merged_DF_Variables_Energy_current_homogeneous)
          head(merged_DF_Variables_Energy_current_homogeneous)
          save(merged_DF_Variables_Energy_current_homogeneous, file="merged_DF_Variables_Energy_current_homogeneous.RData")

          #We save into a raster file to see in Arcgis
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          
            #For ene_All
              ras_current_ene_All <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_All) <- newproj
              values(ras_current_ene_All)<-merged_DF_Variables_Energy_current$ene_All
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_All, "ras_current_ene_All.img", overwrite=TRUE)
  
            #For ene_Human
              ras_current_ene_Human <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human) <- newproj
              values(ras_current_ene_Human)<-merged_DF_Variables_Energy_current$ene_Human
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human, "ras_current_ene_Human.img", overwrite=TRUE)

            #For ene_Human_invertebrates
              ras_current_ene_Human_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_invertebrates) <- newproj
              values(ras_current_ene_Human_invertebrates)<-merged_DF_Variables_Energy_current$ene_Human_invertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human_invertebrates, "ras_current_ene_Human_invertebrates.img", overwrite=TRUE)
              
             #For ene_Human_rep_plant
              ras_current_ene_Human_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_rep_plant) <- newproj
              values(ras_current_ene_Human_rep_plant)<-merged_DF_Variables_Energy_current$ene_Human_rep_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human_rep_plant, "ras_current_ene_Human_rep_plant.img", overwrite=TRUE)
              
            #For ene_Human_veg_plant
              ras_current_ene_Human_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_veg_plant) <- newproj
              values(ras_current_ene_Human_veg_plant)<-merged_DF_Variables_Energy_current$ene_Human_veg_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human_veg_plant, "ras_current_ene_Human_veg_plant.img", overwrite=TRUE)
             
              
             #For ene_Human_vertebrates
              ras_current_ene_Human_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_vertebrates) <- newproj
              values(ras_current_ene_Human_vertebrates)<-merged_DF_Variables_Energy_current$ene_Human_vertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human_vertebrates, "ras_current_ene_Human_vertebrates.img", overwrite=TRUE)
              
             #For ene_Human_unk_plant_oth
              ras_current_ene_Human_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Human_unk_plant_oth) <- newproj
              values(ras_current_ene_Human_unk_plant_oth)<-merged_DF_Variables_Energy_current$ene_Human_unk_plant_oth
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Human_unk_plant_oth, "ras_current_ene_Human_unk_plant_oth.img", overwrite=TRUE)
              
            #For ene_Wild
              ras_current_ene_Wild <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild) <- newproj
              values(ras_current_ene_Wild)<-merged_DF_Variables_Energy_current$ene_Wild
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild, "ras_current_ene_Wild.img", overwrite=TRUE)

            #For ene_Wild_invertebrates
              ras_current_ene_Wild_invertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_invertebrates) <- newproj
              values(ras_current_ene_Wild_invertebrates)<-merged_DF_Variables_Energy_current$ene_Wild_invertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild_invertebrates, "ras_current_ene_Wild_invertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_rep_plant
              ras_current_ene_Wild_rep_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_rep_plant) <- newproj
              values(ras_current_ene_Wild_rep_plant)<-merged_DF_Variables_Energy_current$ene_Wild_rep_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild_rep_plant, "ras_current_ene_Wild_rep_plant.img", overwrite=TRUE)
              
            #For ene_Wild_veg_plant
              ras_current_ene_Wild_veg_plant <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_veg_plant) <- newproj
              values(ras_current_ene_Wild_veg_plant)<-merged_DF_Variables_Energy_current$ene_Wild_veg_plant
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild_veg_plant, "ras_current_ene_Wild_veg_plant.img", overwrite=TRUE)
             
              
             #For ene_Wild_vertebrates
              ras_current_ene_Wild_vertebrates <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_vertebrates) <- newproj
              values(ras_current_ene_Wild_vertebrates)<-merged_DF_Variables_Energy_current$ene_Wild_vertebrates
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild_vertebrates, "ras_current_ene_Wild_vertebrates.img", overwrite=TRUE)
              
             #For ene_Wild_unk_plant_oth
              ras_current_ene_Wild_unk_plant_oth <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600) 
              #We defined the projection of our raster:
              projection(ras_current_ene_Wild_unk_plant_oth) <- newproj
              values(ras_current_ene_Wild_unk_plant_oth)<-merged_DF_Variables_Energy_current$ene_Wild_unk_plant_oth
              #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
              writeRaster(ras_current_ene_Wild_unk_plant_oth, "ras_current_ene_Wild_unk_plant_oth.img", overwrite=TRUE)
