

############################################################################################################
#Readme:
############################################################################################################
#R code to calculate descriptive analysis for the SDMs of the brown bear food species
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
#13_Descriptive_Anlalysis_SDM_Food_Species
  #13.1  Plots SDMs food species

################################################################################################    
#13.1  Plots SDMs food species
################################################################################################   

  rm(list=ls()) 
  my_dir <-"writehereyourpath" # path of the folder where your data are stored
  folder_working<-paste0(my_dir,"/13_Descriptive_Anlalysis_SDM_Food_Species")
  setwd(folder_working)

  install.packages("cowplot")
  library(ggplot2)
  library(cowplot) 

  #We load a DF which contain the predicted values and can be used as a list to know which species have been modelled
  eval_DF_all<-read.xlsx(paste0(my_dir,"/4_SDM_Food_species/eval_DF_all3.xlsx"), sheetName = "Sheet1")
  #we add a column to classify the speces as modelled or not
  eval_DF_all$modelled_sps<-1
  merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF<-read.xlsx(paste0(my_dir,"/10_Descriptive_analysis_Brown_Bear_Food-web/merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF.xlsx"), sheetName = "Sheet1")
  merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF[214,17] <- 0
  merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF[276,17] <- 1
  table(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Human_origin, merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Diet_category2)
  table(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Human_origin)
  summary(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF)

  #We merge the data of taxonomy and energy of all species in the diet with the evaluation metrics of the SDMS
  merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Energy_sum[is.na(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Energy_sum)] <-0
  merged_category_energy_sps_eval<-merge(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF, eval_DF_all, by.x="i", by.y = "i",all.x = T)
  write.xlsx(merged_category_energy_sps_eval, file="merged_category_energy_sps_eval.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
  
  merged_category_energy_sps_eval_wild<-subset(merged_category_energy_sps_eval,Human_origin==0)
  write.xlsx(merged_category_energy_sps_eval_wild, file="merged_category_energy_sps_eval_wild.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   

  merged_category_energy_sps_eval_modeled<-subset(merged_category_energy_sps_eval,modelled_sps==1)
  str(merged_category_energy_sps_eval_modeled)#236 species modelled
  sum(merged_category_energy_sps_eval_modeled$Energy_sum)/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$Energy_sum)*100# Represent the 90.47% of the described energy at species level
  View(merged_category_energy_sps_eval_modeled)
  
  #We change the value of Human/Wild for the reindeer
  merged_category_energy_sps_eval_modeled[182,17] <- 0

  #For wild origin species  
    merged_category_energy_sps_eval_modeled_wild<-subset(merged_category_energy_sps_eval_modeled, Human_origin==0)  
    merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_wild<-subset(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF, Human_origin==0)  
    sum(merged_category_energy_sps_eval_modeled_wild$Energy_sum)/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_wild$Energy_sum)*100

    #Statistics for wild species  
    N_pixels_prensences_by_class_DF_wild<-tapply(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_wild$N_pixels_presence, merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_wild$class, FUN=mean)
    N_sps_by_class_DF_wild<-table(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_wild$class)
    GBIF_by_class_wild<-rbind(N_sps_by_class_DF_wild,N_pixels_prensences_by_class_DF_wild)
    #For the modelled species
    Eval_TSS_by_class_DF_wild<-tapply(merged_category_energy_sps_eval_modeled_wild$Testing.data, merged_category_energy_sps_eval_modeled_wild$class, FUN=mean)
    Eval_by_class_Sensitivity_wild<-tapply(merged_category_energy_sps_eval_modeled_wild$Sensitivity, merged_category_energy_sps_eval_modeled_wild$class, FUN=mean)
    Eval_by_class_Specificity_wild<-tapply(merged_category_energy_sps_eval_modeled_wild$Specificity, merged_category_energy_sps_eval_modeled_wild$class, FUN=mean)
    Eval_by_class_DF_wild<-rbind(Eval_TSS_by_class_DF_wild,Eval_by_class_Sensitivity_wild,Eval_by_class_Specificity_wild)

    DF_Agaricomycetes<-as.data.frame(c(NA,NA,NA))
    colnames(DF_Agaricomycetes)<-c("Agaricomycetes")
  
    Eval_by_class_DF2_wild<-cbind(DF_Agaricomycetes,Eval_by_class_DF_wild)
    GBIF_by_class2_wild<-rbind(GBIF_by_class_wild,Eval_by_class_DF2_wild)
    rownames(GBIF_by_class2_wild)<-c("N_species","Average number of pixels with presence","Average TSS", "Average Sensivity", "Average Specifity")
    write.xlsx(GBIF_by_class2_wild, file="GBIF_by_class2_wild.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   

    #Plot 3           
      X11()   
      #dev.off()
      summary(merged_category_energy_sps_eval_modeled_wild)
      str(merged_category_energy_sps_eval_modeled_wild$Diet_category2)
      summary(merged_category_energy_sps_eval_modeled_wild$Diet_category2)
      ord <- c("reproductive_plant_material", "vegetative_plant_material","unknown_plant_material_and_others" , "invertebrates","vertebrates")
      merged_category_energy_sps_eval_modeled_wild$Factor_Diet_category2 <- factor(merged_category_energy_sps_eval_modeled_wild$Diet_category2,levels=ord)
      par(mfrow=c(4,1))
      par(mar=c(3.5,3.5,0.5,0.5))
      par(xaxs="i",yaxs="i")
      par(mgp = c(2.3, 0.6, 0))
      boxplot((Testing.data)~Diet_category2,data=merged_category_energy_sps_eval_modeled_wild,#ylim=c(-1, 1),
              col="grey",
              border=gray(0.4),
              cex.axis=0.9,
              x.axes=F,
              y.axes=F,
              y.labels=T,
              yaxt='n',
              xaxt='n',
              outline=F, ylab="Index Use")
      abline((0),0,lty=5, col="black",lwd=1)
      axis(2, at=seq(-1,1.0,by=0.5), labels=c(seq(-1.0,1.0,by=0.5)),cex.axis=0.9,las=1)#rep('', 11),
      axis(1, at=0:10+0.5, labels=c("","","","","","","","","","",""),cex.axis=0.9)#rep('', 11),
 
  
   
  #X11()   
      #dev.off()
      par(mfrow=c(4,1))
      par(mar=c(3.5,3.5,0.5,0.5))
      par(xaxs="i",yaxs="i")
      par(mgp = c(2.3, 0.6, 0))
 
  X11()   
  # PDF 7 x 1.72 inches (aprox, later we used 6.7)
  #N presences    
  plots <- list()

  myplot<-( ggplot(merged_category_energy_sps_eval_modeled_wild, aes(y = log(N_pixels_presence),
  x = Factor_Diet_category2, fill =Factor_Diet_category2)) +
  geom_violin(colour = "grey50") +
  geom_jitter(shape = 16, position = position_jitter(0.3), size=0.8)+ theme(legend.position='none'))
  plots[[1]] <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))    

  #TSS# Change the size in photshop to 101.97%
  myplot<-(ggplot(merged_category_energy_sps_eval_modeled_wild, aes(y = Testing.data, x = Factor_Diet_category2, fill =Factor_Diet_category2)) +
  geom_violin(colour = "grey50") +
  geom_jitter(shape = 16, position = position_jitter(0.3), size=0.8)+ theme(legend.position='none'))
  plots[[2]] <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))    
   
  #Specifity# Change the size in photshop to 101.30%
  myplot<-(ggplot(merged_category_energy_sps_eval_modeled_wild, aes(y = Specificity, x = Factor_Diet_category2, fill =Factor_Diet_category2)) +
  geom_violin(colour = "grey50") +
  geom_jitter(shape = 16, position = position_jitter(0.3), size=0.8)+ theme(legend.position='none'))
  plots[[3]] <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))    

  #Sensivity # Change the size in photshop to 102.60%
  myplot<-(ggplot(merged_category_energy_sps_eval_modeled_wild, aes(y = Sensitivity, x = Factor_Diet_category2, fill =Factor_Diet_category2)) +
  geom_violin(colour = "grey50") +
  geom_jitter(shape = 16, position = position_jitter(0.3), size=0.8)+ theme(legend.position='none'))
  plots[[4]] <- myplot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))    
  
  # Create a grid of plots
  print(paste0("OUTPUT: Fig. 3 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  print(paste0("Fig. 3. Number of presences and validation of SDMs for wild food species according to category
    (i.e., reproductive plants, vegetative plants, unknown plant material, invertebrates and vertebrates).
    a, Number of presences from the GBIF used to fit a species distribution model (SDM) for each of the 205 wild food species
    out of the total 240 wild species recorded (85.4% of the wild species in the brown bear diet with enough data to fit a SDM).
    b, Predictive quality of SDMs using the true statistics skill (TSS), c, Specificity (percentage of absences correctly predicted)
    and d, Sensitivity (percentage of presences correctly predicted).~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  X11()
  print(plot_grid(plotlist = plots, nrow = 4, ncol=1)) #, unit="mm"

  pdf(file="Figure_3.pdf", width=16*89*0.0393701, height=16*89*0.0393701)#, pointsize = 12
  print(plot_grid(plotlist = plots, nrow = 4, ncol=1)) #, unit="mm"
  dev.off()

  #For human species  
    merged_category_energy_sps_eval_modeled_human<-subset(merged_category_energy_sps_eval_modeled, Human_origin==1)  
    merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_human<-subset(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF, Human_origin==1)  
    sum(merged_category_energy_sps_eval_modeled_human$Energy_sum)/sum(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF_human$Energy_sum)*100

  #For all species:  
      N_pixels_prensences_by_class_DF<-tapply(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$N_pixels_presence, merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$class, FUN=mean)
      N_sps_by_class_DF<-table(merged_sps_list_diet_category3b_taxonomy_Energy_sum_DF$class)
      GBIF_by_class<-rbind(N_sps_by_class_DF,N_pixels_prensences_by_class_DF)
      Eval_TSS_by_class_DF<-tapply(merged_category_energy_sps_eval_modeled$Testing.data, merged_category_energy_sps_eval_modeled$class, FUN=mean)
      Eval_by_class_Sensitivity<-tapply(merged_category_energy_sps_eval_modeled$Sensitivity, merged_category_energy_sps_eval_modeled$class, FUN=mean)
      Eval_by_class_Specificity<-tapply(merged_category_energy_sps_eval_modeled$Specificity, merged_category_energy_sps_eval_modeled$class, FUN=mean)
    Eval_by_class_DF<-rbind(Eval_TSS_by_class_DF,Eval_by_class_Sensitivity,Eval_by_class_Specificity)
    DF_Agaricomycetes<-as.data.frame(c(NA,NA,NA))
    colnames(DF_Agaricomycetes)<-c("Agaricomycetes")
    Eval_by_class_DF2<-cbind(DF_Agaricomycetes,Eval_by_class_DF)
    GBIF_by_class2<-rbind(GBIF_by_class,Eval_by_class_DF2)
    rownames(GBIF_by_class2)<-c("N_species","Average number of pixels with presence","Average TSS", "Average Sensivity", "Average Specifity")
    write.xlsx(GBIF_by_class2, file="GBIF_by_class2.xlsx", sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)   
 
  
  