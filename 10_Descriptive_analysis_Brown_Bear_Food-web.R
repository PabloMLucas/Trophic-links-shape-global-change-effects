

############################################################################################################
#Readme:
############################################################################################################
#R code for do a descriptive analysis of brown bear food web
#Authors: Pablo M. Lucas
#Last update: 13/03/2025

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run script:
  #/6_Calculation_of_Biotic_variables/merged_sps_list_diet_category3b.RData
  #/2_Representative_Food-web_and_Food_species_list/subpopulations_diet2.csv
  #/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names2.csv
#Data output to other scripts:
  #No outputs to other scripts

##############################################################################################                
#Schema
##############################################################################################   
#10_Descriptive_analysis_Brown_Bear_Food-web
  #10.0 Load data
  #10.1 Data of available energy at species level in each subpopulation
  #10.2 Data of energy for species for all subpopulations
  #10.3 We are going to see the spatial variation of the energy by groups of diet
  #10.4 We are going to see the spatial variation of the energy by groups of diet
  #10.5 Plots of diet by subpopulations
    #10.5.1  Energy by diet group
      #10.5.1.1 Analysis of subpopulations_diet.csv the representative diet for all subpopulations.
      #10.5.1.2 Plot
    #10.5.2  Energy human/wild
    #10.5.3  Plot energy identified at species level
    #10.5.4 The three plots of Figure 2                                

################################################################################################################################################################
################################################################################################################################################################

  #10.0 Load data
    rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/10_Descriptive_analysis_Brown_Bear_Food-web")
    setwd(folder_working)
    
    library(readxl)   
    library(sp)  
    library(dplyr)
    library(bipartite) 
    library(xlsx)  
    #We load the data which contain the classification of species in the diet, the GBIF data and human classification
      # READ: We differentiated in this list of species between subspecies of 
      #Canis lupus familiaris and Canis lupus lupus, but we didn't include the separation between the subspecies of 
      #Sus scrofa and Sus scrofa domesticus due to the difficulties of differentiate the GBIF occurrences between domestic
      #and wild subspecies. In addition we include in this list the study Ambarli et al. 2015 which haven't data 
      #of imputed energy at species level and was then excluded from the other calculations but it includes three 
      #new species (Fagus orientalis, Lonicera caucasica and Phaseolus vulgaris). 
      #Overall we considered 276 species/subspecies items.
    load(paste0(my_dir,"/6_Calculation_of_Biotic_variables/merged_sps_list_diet_category3b.RData"))
    #We load the matrix of species and subpoplations
    subpopulations_diet<-read.csv(paste0(my_dir,"/2_Representative_Food-web_and_Food_species_list/subpopulations_diet2.csv"))
    subpopulations_diet<-as.data.frame(subpopulations_diet)
    #Data of the reassignated energy with the information of taxonomy 
    sheet2 <- read.csv(paste0(my_dir,"/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names2.csv")) # 
    table_all_imputed_energy<-as.data.frame(sheet2)
    res <- subset(table_all_imputed_energy, subset = !duplicated(table_all_imputed_energy[c("Species_corrected_names")]),)
    res2<-res[,c(14:19,50)]
    #Merge the list of species with GBIF data withh taxonomy
    merged_sps_list_diet_category3b_taxonomy<-merge(merged_sps_list_diet_category3b,res2, by.x="Species", by.y = "Species_corrected_names",all.x = T)
    # This include all spcies, we have one species with NA in energy because was present but in a very low percentage in a study and we have the brown bear
    merged_sps_list_diet_category3b_taxonomy2<-merged_sps_list_diet_category3b_taxonomy
    write.xlsx(merged_sps_list_diet_category3b_taxonomy2, file="merged_sps_list_diet_category3b_taxonomy2.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
    save(merged_sps_list_diet_category3b_taxonomy2, file="merged_sps_list_diet_category3b_taxonomy2.RData")

  #10.1 Data of available energy at species level in each subpopulation
    subpopulations_energy_avalable_diet_sps<-as.data.frame(colSums(subpopulations_diet[,2:15]))
    colnames(subpopulations_energy_avalable_diet_sps)<-c("Energy_sum")
    summary(subpopulations_energy_avalable_diet_sps)
    subpopulations_energy_avalable_diet_sps2<-subpopulations_energy_avalable_diet_sps 
    write.csv(file="subpopulations_energy_avalable_diet_sps2.csv", x=subpopulations_energy_avalable_diet_sps2, row.names = T)# Supplemenry Table 10. 
    print(paste0("OUTPUT: Supplemenry Table 10.  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 10. Sum of rates of Estimated Dietary Energy Content (rEDEC) described at the species level for each subpopulation (rEDECSubp).~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(subpopulations_energy_avalable_diet_sps2)

  #10.2 Data of energy for species for all subpopulations
    subpopulations_diet_sps<-as.data.frame(rowSums(subpopulations_diet[,2:15]))
    colnames(subpopulations_diet_sps)<-c("Energy_sum")
    species<-as.data.frame((subpopulations_diet[,1]))
    colnames(species)<-c("Species")
    Energy_sum_DF<-as.data.frame(cbind(species,subpopulations_diet_sps))
    format(Energy_sum_DF, scientific=F)
    subset(Energy_sum_DF, Energy_sum_DF$Energy_sum<0.000000000000000000000000001)# 6 Species with zero energy 
    merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF<-merge(merged_sps_list_diet_category3b_taxonomy,Energy_sum_DF, by.x="Species", by.y = "Species",all.x = T)
    write.xlsx(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF, file="merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
    merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2<-merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF
    #Supplementary Table 8:
    write.xlsx(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2, file="merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
    save(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2, file="merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2.RData")# Supplemenry Table 8. 
    print(paste0("OUTPUT: Supplemenry Table 8.  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplemenry Table 8. List of brown bear food species. For each species we show its scientific name (Species), a numeric identificator (i), 
    the GBIF key (GBIF_key), the total number of occurrences in GBIF, the number of GBIF occurrences filtered using Coordinatecleaner package (Zizka et al 2019),
    the number of remaining occurrences (Filter_occur), the number of pixels of 1x1 km with presence of the species which will be used to model the species 
    distribution (N_pixels_presence), the most frequents diet categories in all subpopulations (Diet_category and Diet_category2), a description of the human
    origin (Human_origin, 1 indicating human and 0 wild), the taxonomy of the species, and the sum of the rEDEC of the species among all subpopulations. We differentiated
    in this list of species between subspecies of Canis lupus familiaris and Canis lupus lupus, but we didn't differentiated between subspecies of Sus scrofa
    (e.g. Sus scrofa domesticus) due to the difficulties to differentiate the GBIF occurrences between domestic and wild subspecies. In addition we include in this list
    the study Ambarli et al. 2015 which haven't data of imputed energy at species level and was then excluded from the other calculations but it includes three new species
    (Fagus orientalis, Lonicera caucasica and Phaseolus vulgaris). Overall we considered 276 species/subspecies items.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2)

    #We have three species with no energy 
    new_DF <- merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF[is.na(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Energy_sum),]  
    #We ommit species with zero frequency and zero energy 
    merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2<-na.omit(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF)
    Energy_by_cat_DF<-tapply(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum, merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Diet_category2, FUN=sum)
    #Percentage of energy by category
      Energy_by_cat_DF[1]/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum)*100
      Energy_by_cat_DF[2]/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum)*100
      Energy_by_cat_DF[3]/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum)*100
      Energy_by_cat_DF[4]/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum)*100
      Energy_by_cat_DF[5]/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2$Energy_sum)*100
    #Frequency of species by category     
    print(paste0("OUTPUT: Fig. 1.  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Fig. 1. Diagram showing our model system to assess the importance of biotic interactions in understanding the consequences
      of global change for biodiversity. a, Construction of a database with detailed explicit knowledge of biotic interactions
      (in our model system brown bear food species in Europe) based on a literature review which accounts for the spatial 
      variability of interactions. .~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    table(merged_sps_list_diet_category3b$Diet_category)
    table(merged_sps_list_diet_category3b$Diet_category2,merged_sps_list_diet_category3b$Human_origin)
    subset(merged_sps_list_diet_category3b, Human_origin==1)
    (67+57+86)/276*100 #Percentage plants
      (67)/276*100 #Reprod plant material
      (86)/276*100 #vegetative_plant_material 
      (57)/276*100 #unknown_plant_material
    #########################################
    (28+36)/276*100 #Percentage animals
      (36)/276*100 # Vertebrates
      (28)/276*100 # Invertebrates

  #10.3 We are going to see the spatial variation of the energy by groups of diet
    merged_subpop_diet_taxo_groups<-merge(subpopulations_diet, merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF, by.x="species", by.y = "Species",all.x = T)
    Energy_by_cat_subpop_1_DF<-tapply(merged_subpop_diet_taxo_groups[,(2)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum) / subpopulations_energy_avalable_diet_sps[1,1] *100
    Energy_by_cat_subpop_2_DF<-tapply(merged_subpop_diet_taxo_groups[,(3)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum) / subpopulations_energy_avalable_diet_sps[2,1] *100
    Energy_by_cat_subpop_3_DF<-tapply(merged_subpop_diet_taxo_groups[,(4)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[3,1] *100
    Energy_by_cat_subpop_4_DF<-tapply(merged_subpop_diet_taxo_groups[,(5)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[4,1] *100
    Energy_by_cat_subpop_5_DF<-tapply(merged_subpop_diet_taxo_groups[,(6)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[5,1] *100
    Energy_by_cat_subpop_6_DF<-tapply(merged_subpop_diet_taxo_groups[,(7)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[6,1] *100
    Energy_by_cat_subpop_7_DF<-tapply(merged_subpop_diet_taxo_groups[,(8)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[7,1] *100
    Energy_by_cat_subpop_8_DF<-tapply(merged_subpop_diet_taxo_groups[,(9)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[8,1] *100
    Energy_by_cat_subpop_9_DF<-tapply(merged_subpop_diet_taxo_groups[,(10)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[9,1] *100
    Energy_by_cat_subpop_10_DF<-tapply(merged_subpop_diet_taxo_groups[,(11)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[10,1] *100
    Energy_by_cat_subpop_11_DF<-tapply(merged_subpop_diet_taxo_groups[,(12)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[11,1] *100
    Energy_by_cat_subpop_12_DF<-tapply(merged_subpop_diet_taxo_groups[,(13)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[12,1] *100
    Energy_by_cat_subpop_13_DF<-tapply(merged_subpop_diet_taxo_groups[,(14)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[13,1] *100
    Energy_by_cat_subpop_14_DF<-tapply(merged_subpop_diet_taxo_groups[,(15)], merged_subpop_diet_taxo_groups$Diet_category2, FUN=sum)/ subpopulations_energy_avalable_diet_sps[14,1] *100
    Energy_by_cat_subpop_ALL_DF_relative <- rbind(Energy_by_cat_subpop_1_DF,Energy_by_cat_subpop_2_DF,Energy_by_cat_subpop_3_DF,Energy_by_cat_subpop_4_DF,
    Energy_by_cat_subpop_5_DF,Energy_by_cat_subpop_6_DF,Energy_by_cat_subpop_7_DF,Energy_by_cat_subpop_8_DF,Energy_by_cat_subpop_9_DF,
    Energy_by_cat_subpop_10_DF,Energy_by_cat_subpop_11_DF,Energy_by_cat_subpop_12_DF,Energy_by_cat_subpop_13_DF,Energy_by_cat_subpop_14_DF)  
    rownames(Energy_by_cat_subpop_ALL_DF_relative)<-colnames(merged_subpop_diet_taxo_groups[,c(2:15)])  
    
    # Supplementary Table 11:
    write.xlsx(Energy_by_cat_subpop_ALL_DF_relative, file="Energy_by_cat_subpop_ALL_DF_relative.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)  # Supplemenry Table 11. 
    print(paste0("OUTPUT: Supplemenry Table 11.  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 11. Sum of rates of Estimated Dietary Energy Content, rEDEC, described at the species level for each subpopulation by diet group.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2) 
    
  #10.4 We are going to see the spatial variation of the energy by groups of diet
    head(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF)  
    Energy_by_cat_subpop_1_DF<-tapply(merged_subpop_diet_taxo_groups[,(2)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum) / subpopulations_energy_avalable_diet_sps[1,1] *100
    Energy_by_cat_subpop_2_DF<-tapply(merged_subpop_diet_taxo_groups[,(3)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum) / subpopulations_energy_avalable_diet_sps[2,1] *100
    Energy_by_cat_subpop_3_DF<-tapply(merged_subpop_diet_taxo_groups[,(4)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[3,1] *100
    Energy_by_cat_subpop_4_DF<-tapply(merged_subpop_diet_taxo_groups[,(5)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[4,1] *100
    Energy_by_cat_subpop_5_DF<-tapply(merged_subpop_diet_taxo_groups[,(6)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[5,1] *100
    Energy_by_cat_subpop_6_DF<-tapply(merged_subpop_diet_taxo_groups[,(7)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[6,1] *100
    Energy_by_cat_subpop_7_DF<-tapply(merged_subpop_diet_taxo_groups[,(8)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[7,1] *100
    Energy_by_cat_subpop_8_DF<-tapply(merged_subpop_diet_taxo_groups[,(9)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[8,1] *100
    Energy_by_cat_subpop_9_DF<-tapply(merged_subpop_diet_taxo_groups[,(10)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[9,1] *100
    Energy_by_cat_subpop_10_DF<-tapply(merged_subpop_diet_taxo_groups[,(11)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[10,1] *100
    Energy_by_cat_subpop_11_DF<-tapply(merged_subpop_diet_taxo_groups[,(12)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[11,1] *100
    Energy_by_cat_subpop_12_DF<-tapply(merged_subpop_diet_taxo_groups[,(13)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[12,1] *100
    Energy_by_cat_subpop_13_DF<-tapply(merged_subpop_diet_taxo_groups[,(14)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[13,1] *100
    Energy_by_cat_subpop_14_DF<-tapply(merged_subpop_diet_taxo_groups[,(15)], merged_subpop_diet_taxo_groups$Human_origin, FUN=sum)/ subpopulations_energy_avalable_diet_sps[14,1] *100
  
    Energy_by_cat_subpop_ALL_DF_relative_human <- rbind(Energy_by_cat_subpop_1_DF,Energy_by_cat_subpop_2_DF,Energy_by_cat_subpop_3_DF,Energy_by_cat_subpop_4_DF,
    Energy_by_cat_subpop_5_DF,Energy_by_cat_subpop_6_DF,Energy_by_cat_subpop_7_DF,Energy_by_cat_subpop_8_DF,Energy_by_cat_subpop_9_DF,
    Energy_by_cat_subpop_10_DF,Energy_by_cat_subpop_11_DF,Energy_by_cat_subpop_12_DF,Energy_by_cat_subpop_13_DF,Energy_by_cat_subpop_14_DF)  
    colnames(Energy_by_cat_subpop_ALL_DF_relative_human)<-c("Wild","Human")  
    rownames(Energy_by_cat_subpop_ALL_DF_relative_human)<-colnames(merged_subpop_diet_taxo_groups[,c(2:15)])  
    # Supplementary Table 12:
    write.xlsx(Energy_by_cat_subpop_ALL_DF_relative_human, file="Energy_by_cat_subpop_ALL_DF_relative_human.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
    print(paste0("OUTPUT: Supplemenry Table 12.  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 12. Rates (from the total described at the species level for each subpopulation) of Estimated Dietary Energy Content, rEDEC, by origin (Wild/Human) for each subpopulation.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF2) 

  #10.5 Plots of diet by subpopulations
    rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/10_Descriptive_analysis_Brown_Bear_Food-web")
    setwd(folder_working)

    library(sp)
    library(readxl) 
    library(dplyr)
    library(bipartite)
    library(xlsx)

    #10.5.1  Energy by diet group
      Energy_by_cat_subpop_ALL_DF_relative<-read.xlsx("Energy_by_cat_subpop_ALL_DF_relative.xlsx", sheetName = "Sheet1",row.names = TRUE)
      Energy_by_cat_subpop_ALL_DF_relative_traposed<-as.data.frame(t(Energy_by_cat_subpop_ALL_DF_relative))
      colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed<-as.data.frame(colnames(Energy_by_cat_subpop_ALL_DF_relative_traposed))               
      colnames(colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed)<-c("Subpopulation") 
      meta_data_2<-read.xlsx(paste0(folder_diet_data,"/meta_data_2b.xlsx"), sheetName = "Sheet1")
      meta_data_3<-subset(meta_data_2, Included=="yes")    
      meta_data_3_unique<-meta_data_3[!duplicated(meta_data_3$Subpopulation),]
      merged_subpopulations<-merge(colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed, meta_data_3_unique, by.x="Subpopulation", by.y = "Subpopulation",all.x = T,all.y = F)
      merged_subpopulations2<-merged_subpopulations[,c(1,3)]
      merged_subpopulations3<-merged_subpopulations2[with(merged_subpopulations2, order(latitude_dd)), ]
      rownames(merged_subpopulations3)<-merged_subpopulations3$Subpopulation
      merged_subpopulations3_t<-t(merged_subpopulations3)
      Energy_by_cat_subpop_ALL_DF_relative_traposed2<-Energy_by_cat_subpop_ALL_DF_relative_traposed[colnames(merged_subpopulations3_t)]  
      Energy_by_cat_subpop_ALL_DF_relative_traposed3<-Energy_by_cat_subpop_ALL_DF_relative_traposed2[c("reproductive_plant_material","unknown_plant_material_and_others","vegetative_plant_material","invertebrates","vertebrates"),]  
  
      #10.5.1.1 Analysis of subpopulations_diet.csv the representative diet for all subpopulations.
      # Statistical analysis about differences in % among subpopulations
        Ener_subpop_DF<-Energy_by_cat_subpop_ALL_DF_relative
        Ener_subpop_DF$Subpopulation<-rownames(Ener_subpop_DF)
        merged_subpopulations_analysis<-merge(Ener_subpop_DF, merged_subpopulations3, by.x="Subpopulation", by.y = "Subpopulation",all = T)
        merged_subpopulations_analysis <- merged_subpopulations_analysis[order(merged_subpopulations_analysis$latitude_dd),]  
  
      #10.5.1.2 Plot
        #Figure 2d
          X11()
          #dev.off()
          # PDF 8 x 3.96 inches
          plotweb(Energy_by_cat_subpop_ALL_DF_relative_traposed3, method="normal",col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90, high.lablength=0,low.lablength=0)    
          #Figure 2d with labels
          X11()
          plotweb(Energy_by_cat_subpop_ALL_DF_relative_traposed3, method="normal",col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90)    
     
    #10.5.2  Energy human/wild
      #Plot energy human/wild
      Energy_by_cat_subpop_ALL_DF_relative_human<-read.xlsx("Energy_by_cat_subpop_ALL_DF_relative_human.xlsx", sheetName = "Sheet1",row.names = TRUE)
      Energy_by_cat_subpop_ALL_DF_relative_human
      Energy_by_cat_subpop_ALL_DF_relative_human_traposed<-as.data.frame(t(Energy_by_cat_subpop_ALL_DF_relative_human))
      colnames_Energy_by_cat_subpop_ALL_DF_relative_human_traposed<-as.data.frame(colnames(Energy_by_cat_subpop_ALL_DF_relative_human_traposed))               
      colnames(colnames_Energy_by_cat_subpop_ALL_DF_relative_human_traposed)<-c("Subpopulation") 
      Energy_by_cat_subpop_ALL_DF_relative_human_traposed2<-Energy_by_cat_subpop_ALL_DF_relative_human_traposed[colnames(merged_subpopulations3_t)]  
       
      #10.5.2.1 Analysis of subpopulations_diet.csv the representative diet for all subpopulations
        # Statistical analysis about differences in % among subpopulations
        Ener_subpop_DF<-Energy_by_cat_subpop_ALL_DF_relative_human
        Ener_subpop_DF$Subpopulation<-rownames(Ener_subpop_DF)
        merged_subpopulations_analysis<-merge(Ener_subpop_DF, merged_subpopulations3, by.x="Subpopulation", by.y = "Subpopulation",all = T)
        merged_subpopulations_analysis <- merged_subpopulations_analysis[order(merged_subpopulations_analysis$latitude_dd),]  
      
        #Figure 2c
          X11()
            # PDF 8 x 3.96 inches
          plotweb(Energy_by_cat_subpop_ALL_DF_relative_human_traposed2, method="normal",empty=T, col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90, high.lablength=0,low.lablength=0)    
          #rowMeans(Energy_by_cat_subpop_ALL_DF_relative_human_traposed2)
          #Figure 2c (with labels)
            X11()
            plotweb(Energy_by_cat_subpop_ALL_DF_relative_human_traposed2, method="normal",empty=T, col.interaction="lightgreen",
            col.high= "grey90",col.low= "grey90",text.rot=90)    
      
    #10.5.3  Plot energy identified at species level
      Energy_energy_avalable_subpop_ALL_DF_relative<-read.csv("subpopulations_energy_avalable_diet_sps2.csv")
      colnames(Energy_energy_avalable_subpop_ALL_DF_relative)<-c("Subpopulation","Energy_sum") 
      rownames(Energy_energy_avalable_subpop_ALL_DF_relative)<-Energy_energy_avalable_subpop_ALL_DF_relative$Subpopulation 
          
      Energy_energy_avalable_subpop_ALL_DF_relative$Energy_sum_no_available<-round(1-Energy_energy_avalable_subpop_ALL_DF_relative$Energy_sum,2)
      Energy_energy_avalable_subpop_ALL_DF_relative$Energy_sum<-round(Energy_energy_avalable_subpop_ALL_DF_relative$Energy_sum,2)
          
      Energy_energy_avalable_subpop_ALL_DF_relative2<-Energy_energy_avalable_subpop_ALL_DF_relative[,c(2,3)]
      #rownames(colnames_Energy_energy_avalable_subpop_ALL_DF_relative_traposed)<-c("Subpopulation") 
      Energy_energy_avalable_subpop_ALL_DF_relative_traposed<-as.data.frame(t(Energy_energy_avalable_subpop_ALL_DF_relative2))
      colnames_Energy_energy_avalable_subpop_ALL_DF_relative_traposed<-as.data.frame(colnames(Energy_energy_avalable_subpop_ALL_DF_relative_traposed))               
      colnames(colnames_Energy_energy_avalable_subpop_ALL_DF_relative_traposed)<-c("Subpopulation") 
         
      Energy_energy_avalable_subpop_ALL_DF_relative_traposed2<-Energy_energy_avalable_subpop_ALL_DF_relative_traposed[colnames(merged_subpopulations3_t)]  
      rownames(Energy_energy_avalable_subpop_ALL_DF_relative_traposed2)<-c("Known","Not_Known")
      #Figure 2b
        X11()
        plotweb(Energy_energy_avalable_subpop_ALL_DF_relative_traposed2, method="normal",empty=T, col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90, high.lablength=0,low.lablength=0)    
        #Figure 2b (with labels)
          X11()
          plotweb(Energy_energy_avalable_subpop_ALL_DF_relative_traposed2, method="normal",empty=T, col.interaction="limegreen",
              col.high= "grey90",col.low= "grey90",text.rot=90)  
          
    #10.5.4 The three plots of Figure 2 
      print(paste0("OUTPUT: Fig. 2 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Fig. 2. Brown bear food web in Europe. b, Relative estimated dietary energy content (rEDEC, a proxy for the relative importance of each species)
        identified or not identified at the species level for each of the 14 brown bear subpopulations in Europe. c, Proportion of the rEDEC for wild species and those
        of human origin. d, Proportion of the rEDEC identified at the species level for each food category (U for unknown plant material, V for vegetative plant material,
        I for invertebrates), for each of the 14 brown bear subpopulations in Europe.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      X11()
      par(mfrow=c(3,1)) 
        #Figure 2b (with labels)
        plotweb(Energy_energy_avalable_subpop_ALL_DF_relative_traposed2, method="normal",empty=T, col.interaction="#A5E8C6",
        col.high= "grey90",col.low= "grey90",text.rot=90)    
        #Figure 2c (with labels)
        plotweb(Energy_by_cat_subpop_ALL_DF_relative_human_traposed2, method="normal",empty=T, col.interaction="#A5E8C6",
        col.high= "grey90",col.low= "grey90",text.rot=90)    
        #Figure 2d (with labels)
        plotweb(Energy_by_cat_subpop_ALL_DF_relative_traposed3, method="normal",col.interaction="#A5E8C6", #col.interaction="limegreen",
        col.high= "grey90",col.low= "grey90",text.rot=90)    
        
        
      #10.5.1  Energy by diet group
      Energy_by_cat_subpop_ALL_DF_relative<-read.xlsx("Energy_by_cat_subpop_ALL_DF_relative.xlsx", sheetName = "Sheet1",row.names = TRUE)
      Energy_by_cat_subpop_ALL_DF_relative_traposed<-as.data.frame(t(Energy_by_cat_subpop_ALL_DF_relative))
      colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed<-as.data.frame(colnames(Energy_by_cat_subpop_ALL_DF_relative_traposed))               
      colnames(colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed)<-c("Subpopulation") 
      meta_data_2<-read.xlsx(paste0("C:/Users/Pablo/Documents/Bear/Paper/Submission/R_code/diet/data","/meta_data_2b.xlsx"), sheetName = "Sheet1")
      meta_data_3<-subset(meta_data_2, Included=="yes")    
      meta_data_3_unique<-meta_data_3[!duplicated(meta_data_3$Subpopulation),]
      merged_subpopulations<-merge(colnames_Energy_by_cat_subpop_ALL_DF_relative_traposed, meta_data_3_unique, by.x="Subpopulation", by.y = "Subpopulation",all.x = T,all.y = F)
      merged_subpopulations2<-merged_subpopulations[,c(1,3)]
      merged_subpopulations3<-merged_subpopulations2[with(merged_subpopulations2, order(latitude_dd)), ]
      rownames(merged_subpopulations3)<-merged_subpopulations3$Subpopulation
      merged_subpopulations3_t<-t(merged_subpopulations3)
      Energy_by_cat_subpop_ALL_DF_relative_traposed2<-Energy_by_cat_subpop_ALL_DF_relative_traposed[colnames(merged_subpopulations3_t)]  
      Energy_by_cat_subpop_ALL_DF_relative_traposed3<-Energy_by_cat_subpop_ALL_DF_relative_traposed2[c("reproductive_plant_material","unknown_plant_material_and_others","vegetative_plant_material","invertebrates","vertebrates"),]  
  
      #10.5.1.1 Analysis of subpopulations_diet.csv the representative diet for all subpopulations.
      # Statistical analysis about differences in % among subpopulations
        Ener_subpop_DF<-Energy_by_cat_subpop_ALL_DF_relative
        Ener_subpop_DF$Subpopulation<-rownames(Ener_subpop_DF)
        merged_subpopulations_analysis<-merge(Ener_subpop_DF, merged_subpopulations3, by.x="Subpopulation", by.y = "Subpopulation",all = T)
        merged_subpopulations_analysis <- merged_subpopulations_analysis[order(merged_subpopulations_analysis$latitude_dd),]  
  
      #10.5.1.2 Plot
        #Figure 2d
          X11()
          #dev.off()
          # PDF 8 x 3.96 inches
          plotweb(Energy_by_cat_subpop_ALL_DF_relative_traposed3, method="normal",col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90, high.lablength=0,low.lablength=0)    
          #Figure 2d with labels
          X11()
          plotweb(Energy_by_cat_subpop_ALL_DF_relative_traposed3, method="normal",col.interaction="#A5E8C6", 
          col.high= "grey90",col.low= "grey90",text.rot=90)    
    
        
    # Load required libraries
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    

# Create the stacked bar plot for food categories
    # Load required libraries
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(viridis)  # Colorblind-friendly palettes
    
    # Convert your matrix/dataframe to a long format
    Energy_long <- Energy_by_cat_subpop_ALL_DF_relative_traposed3 %>%
      tibble::rownames_to_column(var = "Category") %>%  # Convert row names to a column
      pivot_longer(cols = -Category, names_to = "Region", values_to = "Percentage")
    
    # Ensure that "Region" maintains its original order
    Energy_long$Region <- factor(Energy_long$Region, levels = colnames(Energy_by_cat_subpop_ALL_DF_relative_traposed3))
    # Create the stacked bar plot with a colorblind-friendly palette
    ggplot(Energy_long, aes(x = Region, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      labs(x = "Region", y = "Percentage", fill = "Diet Category") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), 
      axis.ticks.x = element_line(),  # Keep the ticks
      panel.grid = element_blank(),  # Remove grid lines
      ) +  # Rotate x-axis labels for readability
      scale_fill_viridis_d(option = "plasma")  # Use "plasma" (or "viridis", "cividis", etc.)
           
    #Without legend        
    
    # Create the stacked bar plot with a colorblind-friendly palette and NO legend
    stack_cat<-ggplot(Energy_long, aes(x = Region, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      labs(x = NULL, y = NULL) +  # Remove x-axis label
      theme_minimal() +
      theme(axis.text.x = element_blank(),  # Remove x-axis labels
            #axis.ticks.x = element_line(),  # Keep the ticks for x
            #axis.ticks.length = unit(0.2, "cm"),  # Make x-axis ticks closer
            axis.text.y = element_blank(), 
            #axis.ticks.y = element_line(),  # Add ticks to y-axis
            #panel.grid = element_blank(),  # Remove grid lines
            panel.grid.major = element_line(color = "grey", size = 0.5),  # Main grid lines
            legend.position = "none")+ # Remove legend
            #axis.line = element_line(color = "black")) +  # Make axis lines visible
            scale_fill_viridis_d(option = "plasma")#+  # Colorblind-friendly palette
  #scale_y_continuous(limits = c(0, NA)) +  # Set the y-axis to start at 0
  #coord_cartesian(ylim = c(0, NA))  # Ensure y-axis starts at 0 and does not cut off bars
      #Save as PDF 3.43 x 1.2 inches 27/02/2025 in F:/G/Proyect_name/Database_diet/R_analysis/Imputation of dietary energy
    stack_cat
        # Set up the PDF device
    pdf("Stacked_plot_diet_categories.pdf", width = 3.33, height = 1.5)  # Create PDF file with dimensions 3.43 x 1.3 inches
    stack_cat
    # Close the PDF device
    dev.off()    
    
# Create the stacked bar plot for human/wild categories
  Energy_by_cat_subpop_ALL_DF_relative_human<-read.xlsx("Energy_by_cat_subpop_ALL_DF_relative_human.xlsx", sheetName = "Sheet1",row.names = TRUE)
      Energy_by_cat_subpop_ALL_DF_relative_human
      Energy_by_cat_subpop_ALL_DF_relative_human_traposed<-as.data.frame(t(Energy_by_cat_subpop_ALL_DF_relative_human))
      colnames_Energy_by_cat_subpop_ALL_DF_relative_human_traposed<-as.data.frame(colnames(Energy_by_cat_subpop_ALL_DF_relative_human_traposed))               
      colnames(colnames_Energy_by_cat_subpop_ALL_DF_relative_human_traposed)<-c("Subpopulation") 
      Energy_by_cat_subpop_ALL_DF_relative_human_traposed2<-Energy_by_cat_subpop_ALL_DF_relative_human_traposed[colnames(merged_subpopulations3_t)]  
colSums(Energy_by_cat_subpop_ALL_DF_relative_human_traposed2)
      
Energy_by_cat_subpop_ALL_DF_relative_human_traposed2


# Scale each column so that it sums to 100
Energy_scaled <- Energy_by_cat_subpop_ALL_DF_relative_human_traposed2

# Apply scaling: for each column, divide by the column sum and multiply by 100
Energy_scaled <- apply(Energy_scaled, 2, function(x) x / sum(x) * 100)

# Check if the columns now sum to 100
colSums(Energy_scaled)  # This should return 100 for each column
str(Energy_scaled)
Energy_scaled<-as.data.frame(Energy_scaled)


    # Convert your matrix/dataframe to a long format
    Energy_long <- Energy_scaled %>%
      tibble::rownames_to_column(var = "Category") %>%  # Convert row names to a column
      pivot_longer(cols = -Category, names_to = "Region", values_to = "Percentage")
    
# Set the order of the 'Region' factor explicitly to match the column names order in the data
Energy_long$Region <- factor(Energy_long$Region, levels = colnames(Energy_scaled))
    
    # Create the stacked bar plot with a colorblind-friendly palette and NO legend
    stack_cat<-ggplot(Energy_long, aes(x = Region, y = Percentage, fill = Category)) +
      geom_bar(stat = "identity") +
      labs(x = NULL, y = NULL) +  # Remove x-axis label
      theme_minimal() +
      theme(axis.text.x = element_blank(),  # Remove x-axis labels
            #axis.ticks.x = element_line(),  # Keep the ticks for x
            #axis.ticks.length = unit(0.2, "cm"),  # Make x-axis ticks closer
            axis.text.y = element_blank(), 
            #axis.ticks.y = element_line(),  # Add ticks to y-axis
            #panel.grid = element_blank(),  # Remove grid lines
            panel.grid.major = element_line(color = "grey", size = 0.5),  # Main grid lines
            legend.position = "none")+ # Remove legend
            #axis.line = element_line(color = "black")) +  # Make axis lines visible
            #scale_fill_viridis_d(option = "Set2")#+  # Colorblind-friendly palette
     scale_fill_manual(values = c("#0072B2", "#D55E00")) #+  # Blue & Orange (safe)
    
  #scale_y_continuous(limits = c(0, NA)) +  # Set the y-axis to start at 0
  #coord_cartesian(ylim = c(0, NA))  # Ensure y-axis starts at 0 and does not cut off bars
      #Save as PDF 3.43 x 1.2 inches 27/02/2025 in F:/G/Proyect_name/Database_diet/R_analysis/Imputation of dietary energy
stack_cat
    # Set up the PDF device
pdf("Stacked_plot_diet_human_wild.pdf", width = 3.33, height = 1.5)  # Create PDF file with dimensions 3.43 x 1.3 inches
stack_cat
# Close the PDF device
dev.off()    


        