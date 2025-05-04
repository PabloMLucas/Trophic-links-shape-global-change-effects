

############################################################################################################
#Readme:
############################################################################################################
#R code for construct the representative diet/food-web of brown bear in each subpopulation
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
  #/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names2.csv
  #/1_Construction_of_the_Trophic_Database/meta_data_2b.xlsx
  #Species_list_clean_binary.xlsx # A modified version of the list of species
  #/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2.RData
#Data output:
  #subpopulations_diet2.RData # Supplementary Table 9
  #subpopulations_diet2.csv # Supplementary Table 9
  #Species_list_clean_binary.xlsx # The list of species

##############################################################################################                
#Schema
##############################################################################################                
#2_Representative_Food-web_and_Food_species_list
  #2.1 Construction of the representative diet/food-web of brown bear in each subpopulation
  #2.2 Construction of the list of species

############################################################################################################################
#2.1 Construction of the representative diet/food-web of brown bear in each subpopulation
############################################################################################################################
  
rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/2_Representative_Food-web_and_Food_species_list")
    setwd(folder_working)

  # Required libraries ------------------------------------------------------------------
    #Install required packages 
    if(!requireNamespace("readxl", quietly=TRUE))
      install.packages("readxl", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("dplyr", quietly=TRUE))
      install.packages("dplyr", quiet=TRUE, dependencies=TRUE) 
    if(!requireNamespace("bipartite", quietly=TRUE))
      install.packages("bipartite", quiet=TRUE, dependencies=TRUE)      
    library(readxl)
    library(dplyr)
    library(bipartite)

  #Load the database modified with the reassignation of energy (Supplementary Table 7)
    database_diet <- read.csv(paste0(my_dir,"/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names2.csv")) # header = TRUE, na.strings = "<NA>"
        database_diet$rEDEC_final<-as.numeric(as.character(database_diet$rEDEC_final))
        database_diet_selected_studies<-database_diet
        database_diet_selected_studies$code<-factor(database_diet_selected_studies$code)
        str(levels(database_diet_selected_studies$code))#The list of all itmes have 31 studies 
        str(database_diet)#1362 diet items
        database_diet_selected_studies$kingdom[database_diet_selected_studies$kingdom=="fungi"]<-"Fungi"
        database_diet_selected_studies$kingdom[database_diet_selected_studies$kingdom=="plantae"]<-"Plantae"
        database_diet_selected_studies$kingdom<-factor(database_diet_selected_studies$kingdom)#7 phylums
        database_diet_selected_studies$Species_corrected_names<-factor(database_diet_selected_studies$Species_corrected_names)#7 phylums
      #WE load the metadata of all diet studies in our database (Supplementary Table 4)
        sheet<-read_excel(paste0(my_dir,"/1_Construction_of_the_Trophic_Database/meta_data_2b.xlsx"), sheet = 1)
        metaData<-as.data.frame(sheet)
        metaData<-metaData[,-1]
        metaData_included<-subset(metaData,Included=="yes")
    length(metaData_included$code)
    length(unique(database_diet_selected_studies$code))
    sort(metaData_included$code)
    sort(unique(database_diet_selected_studies$code))
    ori_imp_met<-merge(database_diet_selected_studies, metaData, by.x="code", by.y = "code",all.x = T)
    ori_imp_met_selected<-subset(ori_imp_met,Included=="yes")
    subpopulations<-unique(ori_imp_met_selected$Subpopulation)#14 subpopulations
    subpopulations_character<-as.character(subpopulations)
    species_list_data_studies<- as.data.frame(unique(factor(ori_imp_met_selected$Species_corrected_names)))
    species_list_data_studies<-na.omit(species_list_data_studies)
    colnames(species_list_data_studies)<-c("species")
    species_list_data_studies <- as.data.frame(species_list_data_studies[order(species_list_data_studies$species),]) 
    colnames(species_list_data_studies)<-c("species")
 
    for (i in 1:14){
      subpopulation_in_loop<-subpopulations[i] 
      data_studies_in_loop<-subset(ori_imp_met_selected,Subpopulation==subpopulation_in_loop)
      studies_in_loop<-unique(data_studies_in_loop$code)
      n_studies_in_loop<-length(studies_in_loop)
      sum_n_data_subpopulation<-0
      vec_rEDEC_proportional_0 <- c()
      species_list_data_studies_in_loop<- as.data.frame(unique(factor(data_studies_in_loop$Species_corrected_names)))
      colnames(species_list_data_studies_in_loop)<-c("species")
      species_list_data_studies_in_loop_no_NA<-na.omit(species_list_data_studies_in_loop)
       for (j in 1:n_studies_in_loop){
         study_code_in_loop<-studies_in_loop[j]
         study_in_loop<-subset(data_studies_in_loop,code==study_code_in_loop)
         list_sps_study<-as.data.frame(factor(study_in_loop$Species_corrected_names))
         colnames(list_sps_study)<-c("species")
         rEDEC_proportional_0<-as.data.frame(study_in_loop$rEDEC_final)*(study_in_loop$n_data)
         colnames(rEDEC_proportional_0)<-paste0("rEDEC_proportional_",study_code_in_loop)
         dfff<-cbind(list_sps_study,rEDEC_proportional_0)
         dfff<-na.omit(dfff)
         dfff$species <-factor(dfff$species)
         apply_sum<-as.data.frame(tapply(dfff[,2], dfff$species, sum))
         apply_sum<-data.frame(rownames(apply_sum),apply_sum)
         colnames(apply_sum)<-c("species",paste0("rEDEC_proportional_",study_code_in_loop))
         species_list_data_studies_in_loop_no_NA<-merge(species_list_data_studies_in_loop_no_NA,apply_sum,by.x="species", by.y = "species",all.x = T)
         #This will sum the n of samples from all the studies in that subpopulation
         sum_n_data_subpopulation<-study_in_loop[1,"n_data"]+sum_n_data_subpopulation
       }
      
       for (j in 1:n_studies_in_loop){
       species_list_data_studies_in_loop_no_NA[,c(1+j)]<-species_list_data_studies_in_loop_no_NA[,c(1+j)]/sum_n_data_subpopulation
       }
      
       species_list_data_studies_in_loop_no_NA[is.na(species_list_data_studies_in_loop_no_NA)] <- 0
       energy_in_Pop<-0
       for (j in 1:n_studies_in_loop){ 
       energy_in_Pop<-energy_in_Pop+species_list_data_studies_in_loop_no_NA[,c(j+1)]
       }
       
     species_list_data_studies_in_loop_no_NA<-data.frame(species_list_data_studies_in_loop_no_NA,energy_in_Pop)
     sub_species_list_wiht_energy_in_Pop<-subset(species_list_data_studies_in_loop_no_NA,energy_in_Pop>0)
     DF_loop<-sub_species_list_wiht_energy_in_Pop[,c("species","energy_in_Pop")]
     colnames(DF_loop)<-c("species",subpopulation_in_loop)
     species_list_data_studies<-merge(species_list_data_studies,DF_loop,by.x="species", by.y = "species",all.x = T)
    }    
      colnames(species_list_data_studies)<-c("species",subpopulations_character)
      head(species_list_data_studies)
      na_values<- rowSums(is.na(species_list_data_studies))
      subpopulations_diet<-species_list_data_studies
      subpopulations_diet[is.na(subpopulations_diet)] <- 0

      colSums(subpopulations_diet[,c(2:15)])
      print(paste0("OUTPUT: Supplementary Table 9  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 9. Matrix showing the summarized data of rEDEC for each species by subpopulation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      # Supplementary Table 9
      save(subpopulations_diet, file=paste0(folder_working,"/subpopulations_diet2.RData"))
      write.csv(file=paste0(folder_working,"/subpopulations_diet2.csv"), x=subpopulations_diet, row.names = F)
      subpopulations_diet
      
############################################################################################################################
#2.2 Construction of the list of species
############################################################################################################################

rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/2_Representative_Food-web_and_Food_species_list")
    setwd(folder_working)
    folder_data_input<-paste0(folder_working,"/Data_input")
    folder_data_ouput<-paste0(folder_working,"/Data_output")
  
  # Required libraries ------------------------------------------------------------------
    #Install required packages 
    if(!requireNamespace("readxl", quietly=TRUE))
      install.packages("readxl", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("dplyr", quietly=TRUE))
      install.packages("dplyr", quiet=TRUE, dependencies=TRUE) 
    if(!requireNamespace("bipartite", quietly=TRUE))
      install.packages("bipartite", quiet=TRUE, dependencies=TRUE)      
    library(readxl)
    library(dplyr)
    library(bipartite)

  #Load the database modified with the reassignation of energy (Supplementary Table 7)
      database_diet <- read.csv(paste0(my_dir,"/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2_ALL_corrected_names2.csv")) 
        database_diet$rEDEC_final<-as.numeric(as.character(database_diet$rEDEC_final))
        database_diet_selected_studies<-database_diet
        database_diet_selected_studies$code<-factor(database_diet_selected_studies$code)
        str(levels(database_diet_selected_studies$code))#The list of all itmes have 31 studies 
        str(database_diet)#1362 diet items
        database_diet_selected_studies$kingdom[database_diet_selected_studies$kingdom=="fungi"]<-"Fungi"
        database_diet_selected_studies$kingdom[database_diet_selected_studies$kingdom=="plantae"]<-"Plantae"
        database_diet_selected_studies$kingdom<-factor(database_diet_selected_studies$kingdom)#7 phylums
        database_diet_selected_studies$Species_corrected_names<-factor(database_diet_selected_studies$Species_corrected_names)#7 phylums
      #WE load the metadata of all diet studies in our database 
        sheet<-read_excel(paste0(my_dir,"/1_Construction_of_the_Trophic_Database/meta_data_2b.xlsx"), sheet = 1)
        metaData<-as.data.frame(sheet)
        metaData<-metaData[,-1]
        metaData_included<-subset(metaData,Included=="yes")
    length(metaData_included$code)
    length(unique(database_diet_selected_studies$code))
    sort(metaData_included$code)
    sort(unique(database_diet_selected_studies$code))
    ori_imp_met<-merge(database_diet_selected_studies, metaData, by.x="code", by.y = "code",all.x = T)
    ori_imp_met_selected<-subset(ori_imp_met,Included=="yes")
    subpopulations<-unique(ori_imp_met_selected$Subpopulation)#14 subpopulations
    subpopulations_character<-as.character(subpopulations)
    species_list_data_studies<- as.data.frame(unique(factor(ori_imp_met_selected$Species_corrected_names)))
    species_list_data_studies<-na.omit(species_list_data_studies)
    colnames(species_list_data_studies)<-c("species")
  
    #We are going to add the species detected in Ambarli et al:
    load(file=paste0(my_dir,"/1_Construction_of_the_Trophic_Database/database_original_with_energy_imputed2.RData"))
    database_original_with_energy_imputed2_AMBA<-subset(database_original_with_energy_imputed2,code=="AMBA")  
    species_list_data_studies_AMBA<- as.data.frame(unique(factor(database_original_with_energy_imputed2_AMBA$species)))
    species_list_data_studies_AMBA<-na.omit(species_list_data_studies_AMBA)
    colnames(species_list_data_studies_AMBA)<-c("species")
    species_list_data_studies_AMBA <- as.data.frame(species_list_data_studies_AMBA[order(species_list_data_studies_AMBA$species),]) 
    colnames(species_list_data_studies_AMBA)<-c("species")
    species_list_data_studies2<-rbind(species_list_data_studies,species_list_data_studies_AMBA)
    species_list_data_studies3<- as.data.frame(unique(factor(species_list_data_studies2$species)))
    colnames(species_list_data_studies3)<-c("species")
    species_list_data_studies3$species<-factor(as.character(species_list_data_studies3$species))
    species_list_data_studies4 <- as.data.frame(species_list_data_studies3[order(species_list_data_studies3$species),])
    colnames(species_list_data_studies4)<-c("species")
    species_list_data_studies5 <- as.data.frame(droplevels(species_list_data_studies4[!species_list_data_studies4$species == 'Sus scrofa domesticus',]))
    colnames(species_list_data_studies5)<-c("Species")
    species_list_data_studies5$Species<-as.character(species_list_data_studies5$Species)
    
    species_list_data_studies5$Nuevo_manual_ID<-1:276
    species_list_data_studies5[species_list_data_studies5 == "Ficus carica"] <- "Ficus carica L."
    #This is the list of species in the diet of the brown bear, but we modified manually the order and included additional information
    write.csv(file=paste0(folder_data_ouput,"/species_list_data_studies5.csv"), x=species_list_data_studies5, row.names = F)
    #This is the modified version:
    sheet<-read_excel(paste0(folder_data_input,"/Species_list_clean_binary.xlsx"), sheet = 1)
    write.xlsx(sheet, file=paste0(folder_data_ouput,"/Species_list_clean_binary.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   

    
