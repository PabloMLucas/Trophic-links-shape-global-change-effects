

############################################################################################################
#Readme:
############################################################################################################
#R code to test whether biotic quantitative variables or biotic binary variables better explain the brown bear distribution
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run this script:
  #/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_model_all_5km.Rdata

#Data output to other scripts:
  #DF_AIC_univar_models_habitat_land_use.RData # Univariable models to select the best land use variables
  #DF_AIC_univar_models_habitat_energy2.RData # Univariable models to select the best biotic variables


##############################################################################################                
#Schema
##############################################################################################                
#8_Univarable_SDM_brown_bear_comparing_biotic_variables
  #8.1 We test the models for biotic variables # Supplementary Table 40
  #8.2 We test the models for land use variables # Supplementary Table 67

  ##############################################################################################
  ##############################################################################################
  #8.1 We test the models for biotic variables
  ##############################################################################################
  ##############################################################################################
  rm(list=ls())
  library(usdm)
  library(lme4)
  my_dir <-"writehereyourpath" 
  folder_working<-paste0(my_dir,"/8_Univarable_SDM_brown_bear_comparing_biotic_variables")
  setwd(folder_working)

  #We load the presences
  load(paste0(my_dir,"/12_Construction_Database_Brown_Bear_occurrences/RANDOM_presences_absences_for_model_all_5km.Rdata"))
  presences_absences_5km<-RANDOM_presences_absences_for_model_all_5km
  data.frame(colnames(presences_absences_5km))
  #A null model without predictors and using subpopulations as random factor:
  model_5km_null <- glmer(Europe_presences_0_1~  (1 | factor_subpop), data = presences_absences_5km, family = "binomial")

    presences_absences_5km_energy<-presences_absences_5km[,c(32:57)]
    #See https://stats.stackexchange.com/questions/89172/how-to-scale-new-observations-for-making-predictions-when-the-model-was-fitted-w
    scaled_presences_absences_5km_energy_for_univarmodels <- scale(presences_absences_5km_energy, scale=TRUE)
    save(scaled_presences_absences_5km_energy_for_univarmodels, file="scaled_presences_absences_5km_energy_for_univarmodels.RData")
    DF_scaled_presences_absences_5km_energy_for_univarmodels<-as.data.frame(cbind(scaled_presences_absences_5km_energy_for_univarmodels,presences_absences_5km$Europe_presences_0_1,presences_absences_5km$factor_subpop))
    colnames(DF_scaled_presences_absences_5km_energy_for_univarmodels)[27] <-c("Europe_presences_0_1")
    colnames(DF_scaled_presences_absences_5km_energy_for_univarmodels)[28] <-c("factor_subpop")
    DF_scaled_presences_absences_5km_energy_for_univarmodels$factor_subpop<-as.factor(DF_scaled_presences_absences_5km_energy_for_univarmodels$factor_subpop) 

    #Calculation of univariable models of habitat for variable selection  
    var_names_hab<-colnames(DF_scaled_presences_absences_5km_energy_for_univarmodels)
    names_var_hab<-var_names_hab[c(1:26)]
        DF_AIC_univar_models_habitat<-c()
        for (v in 1:26){
          var_in_loop_habitat<-names_var_hab[v]
          write_model_prefix_habitat<-c("(Europe_presences_0_1)~")
          fixed_term_habitat<-var_in_loop_habitat
          write_model_var_in_loop_habitat_1<-paste(var_in_loop_habitat,sep="")
          write_model_1_habitat<-paste(write_model_prefix_habitat,write_model_var_in_loop_habitat_1,"+ (1 | factor_subpop)")
          from_mod_1_habitat<-as.formula(write_model_1_habitat)
          mod_1_habitat<- glmer(formula = from_mod_1_habitat, data = DF_scaled_presences_absences_5km_energy_for_univarmodels, family = "binomial", na.action = "na.fail")  
          AIC_mod_1_habitat<-AIC(mod_1_habitat)
          DF_AIC_univar_models_habitat<-rbind(DF_AIC_univar_models_habitat,AIC_mod_1_habitat)
        }
        colnames(DF_AIC_univar_models_habitat)<-c("AICc_lin")
        rownames(DF_AIC_univar_models_habitat)<-names_var_hab
        
          
        DF_AIC_univar_models_habitat_energy<-as.data.frame(DF_AIC_univar_models_habitat)
        save(DF_AIC_univar_models_habitat_energy, file="DF_AIC_univar_models_habitat_energy.RData")
        write.csv(DF_AIC_univar_models_habitat_energy, file="DF_AIC_univar_models_habitat_energy.csv")        
       DF_AIC_univar_models_habitat_energy2<-as.data.frame(cbind(DF_AIC_univar_models_habitat_energy[8:13,],DF_AIC_univar_models_habitat_energy[21:26,]))
       DF_AIC_univar_models_habitat_energy2<- DF_AIC_univar_models_habitat_energy2[c(1,3,4,6,2,5),]
       colnames(DF_AIC_univar_models_habitat_energy2)<-c("AIC_Biotic_variables","AIC_Biotic_binary_variables")
       DF_AIC_univar_models_habitat_energy2$Group_of_species<-c("All_species","Reprod_plant","Veget_plant","Unknown_plant","Invertebrates","Vertebrates") 
       DF_AIC_univar_models_habitat_energy2$Biotic_variables<-c("Bio_All_species","Bio_Reprod_plant","Bio_Veget_plant","Bio_Unknown_plant","Bio_Invertebrates","Bio_Vertebrates") 
       DF_AIC_univar_models_habitat_energy2$Biotic_binary_variables<-c("Bio_binary_All_species","Bio_binary_Reprod_plant","Bio_binary_Veget_plant","Bio_binary_Unknown_plant","Bio_binary_Invertebrates","Bio_binary_Vertebrates") 
       DF_AIC_univar_models_habitat_energy2$AIC_Binary_less_Biotic<- DF_AIC_univar_models_habitat_energy2$AIC_Biotic_binary_variables - DF_AIC_univar_models_habitat_energy2$AIC_Biotic_variables
       row_null<-data.frame(AIC(model_5km_null),"-","Null","-","-","-")
       colnames(row_null)<-colnames(DF_AIC_univar_models_habitat_energy2)
       
       DF_AIC_univar_models_habitat_energy2<-rbind(DF_AIC_univar_models_habitat_energy2,row_null)
       DF_AIC_univar_models_habitat_energy2$Delta_AIC_Biotic_var<- DF_AIC_univar_models_habitat_energy2$AIC_Biotic_variables - (min(DF_AIC_univar_models_habitat_energy2$AIC_Biotic_variables))
       DF_AIC_univar_models_habitat_energy2<-DF_AIC_univar_models_habitat_energy2[,c(3,4,5,1,2,6,7)]
       
        save(DF_AIC_univar_models_habitat_energy2, file="DF_AIC_univar_models_habitat_energy2.RData")
        write.csv(DF_AIC_univar_models_habitat_energy2, file="DF_AIC_univar_models_habitat_energy2.csv")# Supplementary Table 40 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 40#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 40. Comparasion of AIC values for univarible models including the biotic variables and the biotic variables binary and the null model with only the intercept~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        DF_AIC_univar_models_habitat_energy2

      #Check the correlation of variables with VIF 
        #WE anlyse collinearity using spearman correlation
        cor_habitat_energy<-cor(scaled_presences_absences_5km_energy_for_univarmodels)#
        write.csv(cor_habitat_energy, file="cor_habitat_energy.csv")

        model_5km_four <- glmer(Europe_presences_0_1~  (1 | factor_subpop)+ene_Wild_rep_plant+ene_Wild_unk_plant_oth+ene_Wild_vertebrates+ene_Wild_invertebrates, data = DF_scaled_presences_absences_5km_energy_for_univarmodels, family = "binomial")
        AIC(model_5km_four)
        summary(model_5km_four)

  ##############################################################################################
  ##############################################################################################
  #8.2 We test the models for land use variables
  ##############################################################################################
  ##############################################################################################
    presences_absences_5km_land_use<-presences_absences_5km[,c(23:31,60,68)]
    #Calculation of univariable models of habitat for variable selection  
    var_names_hab<-colnames(presences_absences_5km_land_use)
    names_var_hab<-var_names_hab[c(1:9)]
        DF_AIC_univar_models_habitat<-c()
        for (v in 1:9){
          var_in_loop_habitat<-names_var_hab[v]
          write_model_prefix_habitat<-c("(Europe_presences_0_1)~")
          fixed_term_habitat<-var_in_loop_habitat
          write_model_var_in_loop_habitat_1<-paste(var_in_loop_habitat,sep="") 
          write_model_1_habitat<-paste(write_model_prefix_habitat,write_model_var_in_loop_habitat_1,"+ (1 | factor_subpop)")
          from_mod_1_habitat<-as.formula(write_model_1_habitat)
          mod_1_habitat<- glmer(formula = from_mod_1_habitat, data = presences_absences_5km_land_use, family = "binomial", na.action = "na.fail")  
          AIC_mod_1_habitat<-AIC(mod_1_habitat)
          DF_AIC_univar_models_habitat<-rbind(DF_AIC_univar_models_habitat,AIC_mod_1_habitat)
        }
        colnames(DF_AIC_univar_models_habitat)<-c("AIC_lin")
        rownames(DF_AIC_univar_models_habitat)<-names_var_hab
        DF_AIC_univar_models_habitat_land_use<-as.data.frame(DF_AIC_univar_models_habitat)
        DF_AIC_univar_models_habitat_land_use$Variable<-rownames(DF_AIC_univar_models_habitat_land_use)
        row_null_land_use<-data.frame(AIC(model_5km_null),"Null")
        colnames(row_null_land_use)<-colnames(DF_AIC_univar_models_habitat_land_use)
        DF_AIC_univar_models_habitat_land_use2<-rbind(DF_AIC_univar_models_habitat_land_use,row_null_land_use)
        rownames(DF_AIC_univar_models_habitat_land_use2)<-c()
        DF_AIC_univar_models_habitat_land_use2<-DF_AIC_univar_models_habitat_land_use2[,c(2,1)]
        save(DF_AIC_univar_models_habitat_land_use, file="DF_AIC_univar_models_habitat_land_use.RData")
        write.csv(DF_AIC_univar_models_habitat_land_use, file="DF_AIC_univar_models_habitat_land_use.csv")
        write.csv(DF_AIC_univar_models_habitat_land_use2, file="DF_AIC_univar_models_habitat_land_use2.csv") # Supplementary Table 68 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OUTPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        print(paste0("OUTPUT: Supplementary Table 68#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 68. Univariable models of land use for select the best land use variables~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        DF_AIC_univar_models_habitat_land_use2
        
        #WE anlyse collinearity using spearman correlation
        presences_absences_5km_land_use_var<-presences_absences_5km_land_use[,c(1:9)]
        cor_habitat_land_use<-cor(presences_absences_5km_land_use_var)#
        print(paste0("OUTPUT: Supplementary Table 67#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        print(paste0("Supplementary Table 67. Correlation among land use variables with current data of brown bear~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
        write.csv(cor_habitat_land_use, file="cor_habitat_land_use.csv")# Supplementary Table 67

        vif_result_habitat_land_use<-vifstep(presences_absences_5km_land_use_var, th=6) 
        vif_result_habitat_land_use
        save(vif_result_habitat_land_use, file="vif_result_habitat_land_use.RData")
        names_select_vif_habitat_land_use<-vif_result_habitat_land_use @ results $Variables   
        save(names_select_vif_habitat_land_use, file="names_select_vif_habitat_land_use.RData")
        
      #Check the correlation of variables
        presences_absences_5km_variables_for_model<-presences_absences_5km[,c(4:30,34:37,40:43)]
        str(presences_absences_5km_variables_for_model)
        correlation_variables_all<-cor(presences_absences_5km_variables_for_model)#
        write.csv(correlation_variables_all, file="correlation_variables_all.csv")
