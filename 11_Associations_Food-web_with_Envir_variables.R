

############################################################################################################
#Readme:
############################################################################################################
#R code to assess associations bwteen environmental variables and brown bear food web
#Authors: Pablo M. Lucas
#Last update: 13/03/2025

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run script:
  #/12_Construction_Database_Brown_Bear_occurrences/DF_Scen_current_all.RData
#Data output to other scripts:
  #No outputs to other scripts

##############################################################################################                
#Schema
##############################################################################################  
#11_Associations_Food-web_with_Envir_variables
  #11.1 Spatial analysis of environmental data for all study areas 
    #11.1.1 Spatial analysis for 6km radius scale 
    #11.1.2 Spatial analysis for 18km radius scale 
    #11.1.3 Spatial analysis for 40km radius scale 
  #11.2 Summary of energy by category for all study areas 
  #11.3 Calculations of univariable models
      #11.3.1 Calculations of univariable models for all studies
      #11.3.2 Calculations of univariable models for selected studies
      #11.3.3 Calculations of multivariable model for selected studies
      #11.3.4 Calculations of multivariable model of diversity for selected studies

################################################################################################################################################################
#11.1 Spatial analysis of environmental data for all study areas 
rm(list=ls()) 
setwd("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy")  
library(rgdal) 
library(raster) 
library(MuMIn)
#install.packages("usdm")
library(usdm)

#We load the environmental variables
load("H:/G/Project_Name/Brown_Bear_Presences/All_study_area/DF_Scen_current_all.RData")

  #11.1.1 Spatial analysis for 6km radius scale 
    print("########################### SCALE 6 KM #########################################")
    buffers_scale <- shapefile("H:/G/Project_Name/Database_diet/Final datasets/Meta_data3_EAEA_buf6km.shp")
    X11()
    plot(buffers_scale)
    df_buffers_scale<-as.data.frame(buffers_scale)  
    df_buffers_scale$code<-as.factor(df_buffers_scale$code)
    list_code<-df_buffers_scale$code
        #Raster of reference: 
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
        e2<- extent(-2800000,3100000,-800000,4800000) 
        #We are going to use Europe Albers Equal Area Conic  
        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
        #We defined the projection of our raster:
        projection(new_raster) <- newproj
        
        variables_in_study_all<-c()
        for (i in 1:47){#
          print(i)
        sub_buffer_in_loop <- subset(buffers_scale, code==list_code[i])
        sub_buffer_in_loop_raster<-rasterize(sub_buffer_in_loop, field=1,new_raster)
        vec_buffer<-extract(sub_buffer_in_loop_raster,e2)
        data_buffer<- cbind(DF_Scen_current_all,vec_buffer)
        #We select the areas in the buffer
        data_in_buffer<- subset(data_buffer,vec_buffer==1)
        #We select the areas with data of land use   
        data_in_buffer1<-subset(data_in_buffer,LC_1 >= 0 & LC_1 < 2)
        data_in_buffer1$factor_subpop<-as.factor(data_in_buffer1$vec_subpop_raster)
        #We select the areas without water 
        data_in_buffer12<-subset(data_in_buffer1,water == 0 )
        #variables_names<-   names(data_in_buffer12) [c(4:31,66:75)]
        data_in_buffer3<- data_in_buffer12 [,c(4:31,66:75)]
        variables_in_study<-colMeans(data_in_buffer3)
        variables_in_study_all<-rbind(variables_in_study_all,variables_in_study)
        }
    variables_in_study_all2<-as.data.frame(variables_in_study_all) 
    variables_in_study_all2$code<-list_code
    variables_in_study_all2<-variables_in_study_all2[,c(ncol(variables_in_study_all2),1:(ncol(variables_in_study_all2)-1))]
    rownames(variables_in_study_all2)<-c()
    variables_in_study_all_6km<-variables_in_study_all2
    save(variables_in_study_all_6km, file="variables_in_study_all_6km.RData")
    
  #11.1.2 Spatial analysis for 18km radius scale 
    print("########################### SCALE 18 KM ########################################")
    buffers_scale <- shapefile("H:/G/Project_Name/Database_diet/Final datasets/Meta_data3_EAEA_buf18km.shp")
      variables_in_study_all<-c()
        for (i in 1:47){#46
          print(i)
        sub_buffer_in_loop <- subset(buffers_scale, code==list_code[i])
        sub_buffer_in_loop_raster<-rasterize(sub_buffer_in_loop, field=1,new_raster)
        vec_buffer<-extract(sub_buffer_in_loop_raster,e2)
        data_buffer<- cbind(DF_Scen_current_all,vec_buffer)
        #We select the areas in the buffer
        data_in_buffer<- subset(data_buffer,vec_buffer==1)
        #We select the areas with data of land use   
        data_in_buffer1<-subset(data_in_buffer,LC_1 >= 0 & LC_1 < 2)
        data_in_buffer1$factor_subpop<-as.factor(data_in_buffer1$vec_subpop_raster)
        #We select the areas without water 
        data_in_buffer12<-subset(data_in_buffer1,water == 0 )
        #variables_names<-   names(data_in_buffer12) [c(4:31,66:75)]
        data_in_buffer3<- data_in_buffer12 [,c(4:31,66:75)]
        variables_in_study<-colMeans(data_in_buffer3)
        variables_in_study_all<-rbind(variables_in_study_all,variables_in_study)
        }
    variables_in_study_all2<-as.data.frame(variables_in_study_all) 
    variables_in_study_all2$code<-list_code
    variables_in_study_all2<-variables_in_study_all2[,c(ncol(variables_in_study_all2),1:(ncol(variables_in_study_all2)-1))]
    rownames(variables_in_study_all2)<-c()
    variables_in_study_all_18km<-variables_in_study_all2
    save(variables_in_study_all_18km, file="variables_in_study_all_18km.RData")
    
  #11.1.3 Spatial analysis for 40km radius scale 
    print("########################### SCALE 40 KM ########################################")
    #We read the original IUCN spatial distribution:
    buffers_scale <- shapefile("H:/G/Project_Name/Database_diet/Final datasets/Meta_data3_EAEA_buf40km.shp")
    variables_in_study_all<-c()
      for (i in 1:47){#46
          print(i)
        sub_buffer_in_loop <- subset(buffers_scale, code==list_code[i])
        sub_buffer_in_loop_raster<-rasterize(sub_buffer_in_loop, field=1,new_raster)
        vec_buffer<-extract(sub_buffer_in_loop_raster,e2)
        data_buffer<- cbind(DF_Scen_current_all,vec_buffer)
        #We select the areas in the buffer
        data_in_buffer<- subset(data_buffer,vec_buffer==1)
        #We select the areas with data of land use   
        data_in_buffer1<-subset(data_in_buffer,LC_1 >= 0 & LC_1 < 2)
        data_in_buffer1$factor_subpop<-as.factor(data_in_buffer1$vec_subpop_raster)
        #We select the areas without water 
        data_in_buffer12<-subset(data_in_buffer1,water == 0 )
        #variables_names<-   names(data_in_buffer12) [c(4:31,66:75)]
        data_in_buffer3<- data_in_buffer12 [,c(4:31,66:75)]
        variables_in_study<-colMeans(data_in_buffer3)
        variables_in_study_all<-rbind(variables_in_study_all,variables_in_study)
      }
    variables_in_study_all2<-as.data.frame(variables_in_study_all) 
    variables_in_study_all2$code<-list_code
    variables_in_study_all2<-variables_in_study_all2[,c(ncol(variables_in_study_all2),1:(ncol(variables_in_study_all2)-1))]
    rownames(variables_in_study_all2)<-c()
    variables_in_study_all_40km<-variables_in_study_all2
    save(variables_in_study_all_40km, file="variables_in_study_all_40km.RData")

#11.2 Summary of energy by category for all study areas 
    database_all_studies<-read.csv("database_original_with_energy_imputed2.csv")
    database_all_studies$Diet_category2<-database_all_studies$diet2
    levels(database_all_studies$Diet_category2) <- c(levels(database_all_studies$Diet_category2),"unknown_plant_material_and_others")
    database_all_studies$Diet_category2[database_all_studies$Diet_category2 == 'fungi_lichens_bryophytes_algae'] <- 'unknown_plant_material_and_others'
    database_all_studies$Diet_category2[database_all_studies$Diet_category2 == 'unknown_plant_material'] <- 'unknown_plant_material_and_others'
    database_all_studies$Diet_category2<-factor(database_all_studies$Diet_category2)
    table(database_all_studies$Diet_category2)
    list_categories<-data.frame(levels(database_all_studies$Diet_category2))
    colnames(list_categories)<-c("diet_categories")
    df_categories_all<-list_categories
      variables_in_study_all<-c()
      for (i in 1:47){
      print(i)
      sub_buffer_in_loop <- subset(database_all_studies, code==list_code[i])
      sub_buffer_in_loop$rEDEC_imputed[is.na(sub_buffer_in_loop$rEDEC_imputed)]<-0
      sum_REDEC<- aggregate(sub_buffer_in_loop$rEDEC_imputed, list(sub_buffer_in_loop$Diet_category2), FUN=sum) 
      colnames(sum_REDEC)<-c("diet_categories",list_code[i]) 
      df_categories_all<-merge(df_categories_all,sum_REDEC, by="diet_categories",all.x=T ) 
      }
    df_categories_all2<-df_categories_all[-1]
    df_categories_all_t<- t(df_categories_all2)
    df_categories_all_t[is.na(df_categories_all_t)]<-0
    colnames(df_categories_all_t) <-as.vector(df_categories_all[[1]])
    df_categories_all_t<-as.data.frame(df_categories_all_t) 
    df_categories_all_t$code<-list_code
    df_categories_all_t<-df_categories_all_t[,c(ncol(df_categories_all_t),1:(ncol(df_categories_all_t)-1))]
    rownames(df_categories_all_t)<-c()
    df_categories_all_rEDEC<-df_categories_all_t
    save(df_categories_all_rEDEC, file="df_categories_all_rEDEC.RData")

#11.3 Calculations of univariable models
    #11.3.1 Calculations of univariable models for all studies
      #We select the infor of latitutde, longitude code and subpopulation from all locations of studies 
      df_buffers_scale_sub<-df_buffers_scale[,c(2,3,19,23)]
      df_buffers_scale_sub$Subpopulation<-as.factor(df_buffers_scale_sub$Subpopulat)  
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      merged_variables_diet<-merge(df_buffers_scale_sub,df_categories_all_rEDEC,by="code",all=T)
      scale_data<-list(variables_in_study_all_6km,variables_in_study_all_18km,variables_in_study_all_40km)
      scale_names<-c("6km","18km","40km")
      results_lm_all<-c()
      for (s in 1:3){
        scale_selected <- scale_data[[s]]
        scale_names_selected<-scale_names[s]
        merged_variables_diet2<-merge(merged_variables_diet,scale_selected,by="code",all=T)
        variables_names<-names(merged_variables_diet2[,c(2,12:39)] ) 
        response_names<-names(merged_variables_diet2[,c(5,7:11)] )
        for (r in 1:(length(response_names))){
          response_selected <- response_names[r]
          for (v in 1:(length(variables_names))){
            variable_selected <- variables_names[v]
            print(variable_selected)
            response_selected <- response_names[r] 
            variable_response_selected<- paste0(c(response_selected,"~",variable_selected), collapse =" " )
            variable_response_selected
            formula_lm <- formula(variable_response_selected)
            model_lm <- lm(formula_lm,data=merged_variables_diet2)
            pvalue_model_lm<-summary(model_lm)$coefficients[2,4]  
            rsquared_model_lm<-summary(model_lm)$r.squared  
            AICc_model_lm<-AICc(model_lm)
            estimate_model_lm<-summary(model_lm)$ coefficients[2,1]
            Std_e_model_lm<-summary(model_lm)$ coefficients[2,2]
            results_lm<-cbind(estimate_model_lm,Std_e_model_lm,pvalue_model_lm,rsquared_model_lm,AICc_model_lm)
            results_lm<-as.data.frame(cbind(scale_names_selected,response_selected,variable_selected,results_lm))
            results_lm_all<-rbind(results_lm_all,results_lm)
          }
        }    
      }
      write.csv(file="results_lm_all.csv", x=results_lm_all, row.names = T)
      save(results_lm_all, file='results_lm_all.RData')#Supplementary Table 13# 
      print(paste0("OUTPUT: Supplemenry Table 13. #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 13. Results for the univariable lineal models, including all the 47 reviewed locations, n = 47 (See Supplementary Table 4), for explaining the rate of energy dietary content for each category of diet (Category of diet) as a function of climatic, land use and latitude (in decimal degrees; latitude_d) variables (Predictor; See Supplementary Tables 59-60 for the definition of predictors of climate and land use). Estimate and SE show the estimate and the standard error of the predictor variable used in the model. Models with significant p-values are in bold and the best models (based in AICc) for each diet category are shaded in grey.  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

    #11.3.2 Calculations of univariable models for selected studies
      #We select the infor of latitutde, longitude code and subpopulation from all locations of studies 
      df_buffers_scale_sub<-df_buffers_scale[,c(2,3,19,23,24)]
      df_buffers_scale_sub$Subpopulation<-as.factor(df_buffers_scale_sub$Subpopulat) 
      df_buffers_scale_sub$Included<-as.factor(df_buffers_scale_sub$Included)  
      df_buffers_scale_sub<-subset(df_buffers_scale_sub,Included=="yes")
      str(df_buffers_scale_sub)#31 studies selected
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      merged_variables_diet<-merge(df_buffers_scale_sub,df_categories_all_rEDEC,by="code",all.x=T)
      scale_data<-list(variables_in_study_all_6km,variables_in_study_all_18km,variables_in_study_all_40km)
      scale_names<-c("6km","18km","40km")
      results_lm_all<-c()
      for (s in 1:3){
        scale_selected <- scale_data[[s]]
        scale_names_selected<-scale_names[s]
        merged_variables_diet2<-merge(merged_variables_diet,scale_selected,by="code",all.x=T)
        variables_names<-names(merged_variables_diet2[,c(2,12:39)] ) 
        response_names<-names(merged_variables_diet2[,c(5,7,9:11)] )
        for (r in 1:(length(response_names))){
          response_selected <- response_names[r]
          for (v in 1:(length(variables_names))){
            variable_selected <- variables_names[v]
            print(variable_selected)
            response_selected <- response_names[r] 
            variable_response_selected<- paste0(c(response_selected,"~",variable_selected), collapse =" " )
            variable_response_selected
            formula_lm <- formula(variable_response_selected)
            model_lm <- lm(formula_lm,data=merged_variables_diet2)
            pvalue_model_lm<-summary(model_lm)$coefficients[2,4]  
            rsquared_model_lm<-summary(model_lm)$r.squared  
            AICc_model_lm<-AICc(model_lm)
            estimate_model_lm<-summary(model_lm)$ coefficients[2,1]
            Std_e_model_lm<-summary(model_lm)$ coefficients[2,2]
            results_lm<-cbind(estimate_model_lm,Std_e_model_lm,pvalue_model_lm,rsquared_model_lm,AICc_model_lm)
            results_lm<-as.data.frame(cbind(scale_names_selected,response_selected,variable_selected,results_lm))
            results_lm_all<-rbind(results_lm_all,results_lm)
          }
        }    
      }
      results_lm_all_selected_studies<-results_lm_all
      write.csv(file="results_lm_all_selected_studies.csv", x=results_lm_all, row.names = T)
      save(results_lm_all_selected_studies, file='results_lm_all_selected_studies.RData')#Supplementary Table 14# 
      print(paste0("OUTPUT: Supplemenry Table 14. #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
      print(paste0("Supplementary Table 14. Results for the univariable lineal models, including only the 31 selected locations, n = 31 (See Supplementary Table 4), for explaining the rate of energy dietary content for each category of diet (Category of diet) as a function of climatic, land use and latitude (in decimal degrees; latitude_d) variables (Predictor; See Supplementary Tables 59-60 for the definition of predictors of climate and land use). Estimate and SE show the estimate and the standard error of the predictor variable used in the model. Models with significant p-values are in bold and the best models (based in AICc) for each diet category are shaded in grey.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))

 
    #11.3.3 Calculations of multivariable model for selected studies
      library(rgdal) 
      library(raster) 
      library(MuMIn)
      library(usdm)
      load("df_categories_all_rEDEC.RData")
      load("variables_in_study_all_6km.RData")
      load("variables_in_study_all_18km.RData")
      load("variables_in_study_all_40km.RData")
      #We select columns which have data describing each study site (not scale dependen data, latitude, longitude, code subpopulation and Included/not included)
      buffers_scale <- shapefile("Meta_data3_EAEA_buf6km.shp")# No worries for the scale, we could have used other scale
      df_buffers_scale<-as.data.frame(buffers_scale)  
      df_buffers_scale$code<-as.factor(df_buffers_scale$code)
      df_buffers_scale_sub<-df_buffers_scale[,c(2,3,19,23,24)]
      str(df_buffers_scale_sub) 
      df_buffers_scale_sub$Subpopulation<-as.factor(df_buffers_scale_sub$Subpopulat) 
      df_buffers_scale_sub$Included<-as.factor(df_buffers_scale_sub$Included)  
      df_buffers_scale_sub<-subset(df_buffers_scale_sub,Included=="yes")
      str(df_buffers_scale_sub)#31 studies selected
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      merged_variables_diet<-merge(df_buffers_scale_sub,df_categories_all_rEDEC,by="code",all.x=T)
      scale_data<-list(variables_in_study_all_6km,variables_in_study_all_18km,variables_in_study_all_40km)
      scale_names<-c("6km","18km","40km")
      
      s=2
      for (s in 1:3){
        print(paste("############################### SCALE ",scale_names[s]," ###############################"))
        scale_selected <- scale_data[[s]]
        scale_names_selected<-scale_names[s]
        merged_variables_diet2<-merge(merged_variables_diet,scale_selected,by="code",all.x=T)
        variables_names<-names(merged_variables_diet2[,c(2,12:39)] ) #The environmental predictors
        response_names<-names(merged_variables_diet2[,c(5,7,9:11)] ) #The response variables, we will develop one model for each diet category
        #VIF calculation to remove high collinear variables
        environmental_var_by_scale<-merged_variables_diet2[,c(2,12:39)]
        vif_environmental_var_by_scale<-vifstep(environmental_var_by_scale, th=10) 
        save(vif_environmental_var_by_scale, file=paste0("vif_environmental_var_by_scale_",scale_names[s],".RData"))#
        vif_result_environmental_var_by_scale<-as.data.frame(vif_environmental_var_by_scale@ results)
        write.csv(vif_result_environmental_var_by_scale, file=paste("vif_result_environmental_var_by_scale",scale_names[s],".csv"))# Supplementary Table 15 Uncorrelated variables based on VIF with a threshold of 10 ##################################################################################output
        names_select_vif_environmental_var_by_scale<-vif_environmental_var_by_scale @ results $Variables
        save(names_select_vif_environmental_var_by_scale, file=paste0("names_select_vif_environmental_var_by_scale_",scale_names[s],".RData"))#
        write.csv(names_select_vif_environmental_var_by_scale, file=paste("names_select_vif_environmental_var_by_scale_",scale_names[s],".csv"))# 
        environmental_var_by_scale_selected<-   environmental_var_by_scale[,names_select_vif_environmental_var_by_scale]
        variables_selected <- colnames(environmental_var_by_scale_selected)
        print(variables_selected)
        #A loop for each diet category we fit a model
        for (r in 1:(length(response_names))){
          response_selected <- response_names[r]
          print(paste("##### DIET CATEGORY ",response_selected," ###############################"))
          variables_selected_list <- paste(c(variables_selected),collapse="+")
          variables_selected_list2 <- paste(c(response_selected,"~",variables_selected_list),collapse="")
          mod_lm <- formula(variables_selected_list2)
          #We write a full model which include all variables
          model_lm <- lm(mod_lm,data=merged_variables_diet2, na.action = "na.fail")
          #With dredge funtion we are going to obtain all possible models combining the variables included in the full model 
          dr_model_lm<- dredge(model_lm, rank= "AICc")
          save (dr_model_lm, file=paste("dr_model_lm_",scale_names[s],"_",response_selected,".RData"))
          sub_dr_model_lm<-subset(dr_model_lm,delta < 3)
          write.csv(sub_dr_model_lm, file=paste("sub_dr_model_lm_",scale_names[s],"_",response_selected,".csv"))# # Supplementary Table 16-20 Best multivariable models explaining the percentage of each diet category ##############################################################################################################output
          #We get the models subsetting the models which are equally good (delta<3)
          get_models_lm<-get.models(dr_model_lm, subset=delta < 3)
          save (get_models_lm, file=paste("get_models_lm_",scale_names[s],"_",response_selected,".RData"))
          print(paste("n models >3: ",length(get_models_lm)))
          if(length(get_models_lm)>1){
          #WE calculate the average model# S3 method for model.selection
          avg_model<-model.avg(get_models_lm)
          obj_sum<-summary(avg_model)
          coef_avg_model_subset<-obj_sum$coefmat.subset
          save (avg_model, file=paste("avg_model_",scale_names[s],"_",response_selected,".RData"))# Supplementary Table 
          save (coef_avg_model_subset, file=paste("coef_avg_model_subset_",scale_names[s],"_",response_selected,".RData"))# 
          write.csv(coef_avg_model_subset, file=paste("coef_avg_model_subset_",scale_names[s],"_",response_selected,".csv"))# Supplementary Tables 21-25 Averaged model using the best multivariable models###################################output
          }else{
          obj_sum<-summary(get_models_lm[[1]])
          coef_avg_model_subset<-obj_sum$coefficients
          save (coef_avg_model_subset, file=paste("coef_avg_model_subset_",scale_names[s],"_",response_selected,".RData"))#  
          write.csv(coef_avg_model_subset, file=paste("coef_avg_model_subset_",scale_names[s],"_",response_selected,".csv"))# Supplementary Tables 21-25 Averaged model using the best multivariable models############ 
          }
        }
      }          
      
 #Models for invertebrates
      #9.4.3 Plot of the regression line with a 95% confidence interval band around the regression line. These are the best 4 models explaining lagomorpha
      # Load the gridExtra package
      install.packages("gridExtra") 
      install.packages("ggplot2") 
      install.packages("ggrepel")
      install.packages("cowplot")
      library(cowplot)
      library(ggrepel)
      library(gridExtra) 
      library(ggplot2) 
      X11()  
      # Define a common theme for all plots to reduce margins
      common_theme <- theme(
        axis.title.x = element_text(size = 8),  # Reduce size of x-axis label
        axis.title.y = element_text(size = 8),  # Reduce size of y-axis label
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm")  # Reduce plot margins
      )
  data.frame(colnames(merged_variables_diet2))    
      # Create individual ggplot objects with the common theme
  
 # Reproductive_plant_material 
      p1 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = reproductive_plant_material)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p2 <- ggplot(data = merged_variables_diet2, aes(x = LC_7, y = reproductive_plant_material)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p3 <- ggplot(data = merged_variables_diet2, aes(x = LC_9, y = reproductive_plant_material)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme

# Unknown_plant_material_and_others           
      p4 <- ggplot(data = merged_variables_diet2, aes(x = LC_1, y = unknown_plant_material_and_others)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p5 <- ggplot(data = merged_variables_diet2, aes(x = LC_9, y = unknown_plant_material_and_others)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
## Invertebrates      
      p6 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = invertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p7 <- ggplot(data = merged_variables_diet2, aes(x = LC_9, y = invertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
## Vertebrates      
      
      p8 <- ggplot(data = merged_variables_diet2, aes(x = Clim_1, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p9 <- ggplot(data = merged_variables_diet2, aes(x = Clim_2, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p10 <- ggplot(data = merged_variables_diet2, aes(x = Clim_5, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p11 <- ggplot(data = merged_variables_diet2, aes(x = Clim_6, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p12 <- ggplot(data = merged_variables_diet2, aes(x = Clim_9, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p13 <- ggplot(data = merged_variables_diet2, aes(x = Clim_10, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p14 <- ggplot(data = merged_variables_diet2, aes(x = Clim_11, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p15 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      p16 <- ggplot(data = merged_variables_diet2, aes(x = LC_7, y = vertebrates)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme      
      
      # Create a multi-panel plot with reduced space between panels
      multi_panel_plot <- plot_grid(
        #p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
        p1, p6, p8, p9,#p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,
        nrow = 1,
        rel_widths = c(1, 1, 1, 1),  # Equal width for each panel
        align = "h",  # Align panels horizontally
        axis = "tb"  # Align only top and bottom axes
      )
      # Show the multi-panel plot
       multi_panel_plot
          #Save as PDF 4.85 x 1.48 inches 27/02/2025 in F:/G/Proyect_name/Database_diet/R_analysis/Imputation of dietary energy
      
    #Check the models
     #install.packages("performance")
     library("performance")
       
     #reproductive_plant_material  
     reproductive_plant_material_LC3 <-lm(merged_variables_diet2$reproductive_plant_material~merged_variables_diet2$LC_3)
     summary(reproductive_plant_material_LC3)
     X11()
     check_model(reproductive_plant_material_LC3)
     
     #invertebrates  
     invertebrates_LC3 <-lm(merged_variables_diet2$invertebrates~merged_variables_diet2$LC_3)
     summary(invertebrates_LC3)
     X11()
     check_model(invertebrates_LC3)

     #vertebrates  
     vertebrates_Clim_1 <-lm(merged_variables_diet2$vertebrates~merged_variables_diet2$Clim_1)
     summary(vertebrates_Clim_1)
     X11()
     check_model(vertebrates_Clim_1)
       
     vertebrates_Clim_2 <-lm(merged_variables_diet2$vertebrates~merged_variables_diet2$Clim_2)
     summary(vertebrates_Clim_2)
     X11()
     check_model(vertebrates_Clim_2)        
     
     vertebrates_LC3 <-lm(merged_variables_diet2$vertebrates~merged_variables_diet2$LC_3)
     summary(vertebrates_LC3)
     X11()
     check_model(vertebrates_LC3)
     
          vertebrates_LC5 <-lm(merged_variables_diet2$vertebrates~merged_variables_diet2$LC_5)
     summary(vertebrates_LC5)
     X11()
     check_model(vertebrates_LC5)
     
     
     vertebrates_LC7 <-lm(merged_variables_diet2$vertebrates~merged_variables_diet2$LC_7)
     summary(vertebrates_LC7)
     X11()
     check_model(vertebrates_LC7)
          

    #11.3.4 Calculations of multivariable model of diversity for selected studies
      rm(list=ls()) 
      setwd("H:/G/Project_Name/Database_diet/R_analysis/Imputation of dietary energy")  
      library(rgdal) 
      library(raster) 
      library(MuMIn)
      library(usdm)
      library(vegan)

      load("df_categories_all_rEDEC.RData")
      load("variables_in_study_all_6km.RData")
      load("variables_in_study_all_18km.RData")
      load("variables_in_study_all_40km.RData")
      
      #We calculate the diversity indexes
      BCI<-df_categories_all_rEDEC[,-c(1,3)] 
      div_simp <- diversity(BCI, "simpson")
      div_shannon <- diversity(BCI, "shannon")
      div_invsimpson <- diversity(BCI, "invsimpson")
      DF_index_div<-df_categories_all_rEDEC
      DF_index_div$div_simp<-div_simp
      DF_index_div$div_shannon<-div_shannon
      DF_index_div$div_invsimpson<-div_invsimpson
      DF_index_div<-DF_index_div[,c(1,9:11)]
      write.csv(DF_index_div, file=paste("DF_index_div.csv"))# Supplementary Table 26. Diversity indexes of broen bear diet (Simpson, Shannon and inverse Simpsons) for study sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
      #We select columns which have data describing each study site (not scale dependen data, latitude, longitude, code subpopulation and Included/not included)
      buffers_scale <- shapefile("H:/G/Project_Name/Database_diet/Final datasets/Meta_data3_EAEA_buf6km.shp")# 
      df_buffers_scale<-as.data.frame(buffers_scale)  
      df_buffers_scale$code<-as.factor(df_buffers_scale$code)
      df_buffers_scale_sub<-df_buffers_scale[,c(2,3,19,23,24)]
      str(df_buffers_scale_sub) 
      df_buffers_scale_sub$Subpopulation<-as.factor(df_buffers_scale_sub$Subpopulat) 
      df_buffers_scale_sub$Included<-as.factor(df_buffers_scale_sub$Included)  
      df_buffers_scale_sub<-subset(df_buffers_scale_sub,Included=="yes")
      str(df_buffers_scale_sub)#31 studies selected
      order(df_buffers_scale_sub)
  
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      df_buffers_scale_sub<-df_buffers_scale_sub[-4]
      merged_variables_diet<-merge(df_buffers_scale_sub,DF_index_div,by="code",all.x=T)
      newdata <- merged_variables_diet[order(merged_variables_diet$latitude_d),]
   
      scale_data<-list(variables_in_study_all_6km,variables_in_study_all_18km,variables_in_study_all_40km)
      scale_names<-c("6km","18km","40km")
      #data.frame(names(merged_variables_diet2))
      s=2
        for (s in 1:3){
          print(paste("############################### SCALE ",scale_names[s]," ###############################"))
          scale_selected <- scale_data[[s]]
          scale_names_selected<-scale_names[s]
          merged_variables_diet2<-merge(merged_variables_diet,scale_selected,by="code",all.x=T)
          variables_names<-names(merged_variables_diet2[,c(2,8:35)] ) #The environmental predictors
          response_names<-names(merged_variables_diet2[,c(5:7)] ) #The response variables, we will develop one model for each diet category
          #VIF calculation to remove high collinear variables
          environmental_var_by_scale<-merged_variables_diet2[,c(2,8:35)]
          vif_environmental_var_by_scale<-vifstep(environmental_var_by_scale, th=10) 
          save(vif_environmental_var_by_scale, file=paste0("vif_DIV_environmental_var_by_scale_",scale_names[s],".RData"))#
          vif_result_environmental_var_by_scale<-as.data.frame(vif_environmental_var_by_scale@ results)
          write.csv(vif_result_environmental_var_by_scale, file=paste("vif_DIV_result_environmental_var_by_scale",scale_names[s],".csv"))# 
          names_select_vif_environmental_var_by_scale<-vif_environmental_var_by_scale @ results $Variables
          save(names_select_vif_environmental_var_by_scale, file=paste0("names_select_DIV_vif_environmental_var_by_scale_",scale_names[s],".RData"))#
          write.csv(names_select_vif_environmental_var_by_scale, file=paste("names_select_DIV_vif_environmental_var_by_scale_",scale_names[s],".csv"))# 
          environmental_var_by_scale_selected<-   environmental_var_by_scale[,names_select_vif_environmental_var_by_scale]
          variables_selected <- colnames(environmental_var_by_scale_selected)
          print(variables_selected)
          #A loop for each diet category we fit a model
          for (r in 1:(length(response_names))){
            response_selected <- response_names[r]
            print(paste("##### DIVERSITY ",response_selected," ###############################"))
            variables_selected_list <- paste(c(variables_selected),collapse="+")
            variables_selected_list2 <- paste(c(response_selected,"~",variables_selected_list),collapse="")
            mod_lm <- formula(variables_selected_list2)
            #We write a full model which include all variables
            model_lm <- lm(mod_lm,data=merged_variables_diet2, na.action = "na.fail")
            #With dredge funtion we are going to obtain all possible models combining the variables included in the full model 
            dr_model_lm<- dredge(model_lm, rank= "AICc")
            save (dr_model_lm, file=paste("DIV_dr_model_lm_",scale_names[s],"_",response_selected,".RData"))
            sub_dr_model_lm<-subset(dr_model_lm,delta < 3)
            write.csv(sub_dr_model_lm, file=paste("DIV_sub_dr_model_lm_",scale_names[s],"_",response_selected,".csv"))# Supplementary Tables 27-29 ##################################################################################output
            #We get the models subsetting the models which are equally good (delta<3)
            get_models_lm<-get.models(dr_model_lm, subset=delta < 3)
            save (get_models_lm, file=paste("DIV_get_models_lm_",scale_names[s],"_",response_selected,".RData"))
            print(paste("n models >3: ",length(get_models_lm)))
            if(length(get_models_lm)>1){
            #WE calculate the average model# S3 method for model.selection
            avg_model<-model.avg(get_models_lm)
            obj_sum<-summary(avg_model)
            coef_avg_model_subset<-obj_sum$coefmat.subset
            save (avg_model, file=paste("DIV_avg_model_",scale_names[s],"_",response_selected,".RData"))# Supplementary Table 
            save (coef_avg_model_subset, file=paste("DIV_coef_avg_model_subset_",scale_names[s],"_",response_selected,".RData"))# 
            write.csv(coef_avg_model_subset, file=paste("DIV_coef_avg_model_subset_",scale_names[s],"_",response_selected,".csv"))# Supplementary Tables 30-32 #######################################################################################################output
            }else{
            obj_sum<-summary(get_models_lm[[1]])
            coef_avg_model_subset<-obj_sum$coefficients
            save (coef_avg_model_subset, file=paste("DIV_coef_avg_model_subset_",scale_names[s],"_",response_selected,".RData"))
            write.csv(coef_avg_model_subset, file=paste("DIV_coef_avg_model_subset_",scale_names[s],"_",response_selected,".csv"))# S# Supplementary Tables 30-32 ##################################################################
            }
          }
        }          
     
      
       # Define a common theme for all plots to reduce margins
      common_theme <- theme(
        axis.title.x = element_text(size = 8),  # Reduce size of x-axis label
        axis.title.y = element_text(size = 8),  # Reduce size of y-axis label
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm")  # Reduce plot margins
      )
  data.frame(colnames(merged_variables_diet2))    
      # Create individual ggplot objects with the common theme
  
  X11()
  
 # Clim_1 
      p01 <- ggplot(data = merged_variables_diet2, aes(x = Clim_1, y = div_simp)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p02 <- ggplot(data = merged_variables_diet2, aes(x = Clim_1, y = div_shannon)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p03 <- ggplot(data = merged_variables_diet2, aes(x = Clim_1, y = div_invsimpson)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
  
  
 # Clim_2 
      p1 <- ggplot(data = merged_variables_diet2, aes(x = Clim_1, y = div_simp)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p2 <- ggplot(data = merged_variables_diet2, aes(x = Clim_2, y = div_shannon)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p3 <- ggplot(data = merged_variables_diet2, aes(x = Clim_2, y = div_invsimpson)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
 
 # LC_3 
      p4 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = div_simp)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p5 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = div_shannon)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p6 <- ggplot(data = merged_variables_diet2, aes(x = LC_3, y = div_invsimpson)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
 # LC_8 
      p7 <- ggplot(data = merged_variables_diet2, aes(x = LC_8, y = div_simp)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p8 <- ggplot(data = merged_variables_diet2, aes(x = LC_8, y = div_shannon)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      p9 <- ggplot(data = merged_variables_diet2, aes(x = LC_8, y = div_invsimpson)) +
        geom_point(size = 1) +  # Add points
        geom_smooth(method = "lm") +  # Add regression line
        common_theme
      
      
      # Create a multi-panel plot with reduced space between panels
      multi_panel_plot <- plot_grid(
        #p01, p02, p03,p1, p2, p3, p4,p5,p6,p7,p8,p9,#p10,p11,p12,p13,p14,p15,p16,
        p01, p03,p4,p5,#p6,p7,p8,p9,#p10,p11,p12,p13,p14,p15,p16,
        nrow = 1,
        rel_widths = c(1, 1, 1, 1),  # Equal width for each panel
        align = "h",  # Align panels horizontally
        axis = "tb"  # Align only top and bottom axes
      )
      
      
     # Show the multi-panel plot
       multi_panel_plot
          #Save as PDF 4.85 x 1.48 inches 27/02/2025 in 
       

      
    #Check the models
     #install.packages("performance")
     library("performance")
       colnames(merged_variables_diet2)
       
       
 # Clim_1 
     div_simp_clim_1 <-lm(merged_variables_diet2$div_simp~merged_variables_diet2$Clim_1)
     summary(div_simp_clim_1)
     X11()
     check_model(div_simp_clim_1)      
      
 # Clim_1 
     div_shannon_clim_1 <-lm(merged_variables_diet2$div_shannon~merged_variables_diet2$Clim_1)
     summary(div_shannon_clim_1)
     X11()
     check_model(div_shannon_clim_1)      
      
 # Clim_1 
     div_invsimpson_clim_1 <-lm(merged_variables_diet2$div_invsimpson~merged_variables_diet2$Clim_1)
     summary(div_invsimpson_clim_1)
     X11()
     check_model(div_invsimpson_clim_1)      
     
 #################################################
     
 # Clim_2 
     div_simp_clim_2 <-lm(merged_variables_diet2$div_simp~merged_variables_diet2$Clim_2)
     summary(div_simp_clim_2)
     X11()
     check_model(div_simp_clim_2)      
      
 # Clim_2 
     div_shannon_clim_2 <-lm(merged_variables_diet2$div_shannon~merged_variables_diet2$Clim_2)
     summary(div_shannon_clim_2)
     X11()
     check_model(div_shannon_clim_2)      
      
 # Clim_2 
     div_invsimpson_clim_2 <-lm(merged_variables_diet2$div_invsimpson~merged_variables_diet2$Clim_2)
     summary(div_invsimpson_clim_2)
     X11()
     check_model(div_invsimpson_clim_2)      
     
#################################################
     
 # LC_3 
     div_simp_LC_3 <-lm(merged_variables_diet2$div_simp~merged_variables_diet2$LC_3)
     summary(div_simp_LC_3)
     X11()
     check_model(div_simp_LC_3)      
      
 # LC_3 
     div_shannon_LC_3 <-lm(merged_variables_diet2$div_shannon~merged_variables_diet2$LC_3)
     summary(div_shannon_LC_3)
     X11()
     check_model(div_shannon_LC_3)      
      
 # LC_3 
     div_invsimpson_LC_3 <-lm(merged_variables_diet2$div_invsimpson~merged_variables_diet2$LC_3)
     summary(div_invsimpson_LC_3)
     X11()
     check_model(div_invsimpson_LC_3)      
     
#################################################
     
 # LC_8 
     div_simp_LC_8 <-lm(merged_variables_diet2$div_simp~merged_variables_diet2$LC_8)
     summary(div_simp_LC_8)
     X11()
     check_model(div_simp_LC_8)      
      
 # LC_8 
     div_shannon_LC_8 <-lm(merged_variables_diet2$div_shannon~merged_variables_diet2$LC_8)
     summary(div_shannon_LC_8)
     X11()
     check_model(div_shannon_LC_8)      
      
 # LC_8 
     div_invsimpson_LC_8 <-lm(merged_variables_diet2$div_invsimpson~merged_variables_diet2$LC_8)
     summary(div_invsimpson_LC_8)
     X11()
     check_model(div_invsimpson_LC_8)      
       
      
