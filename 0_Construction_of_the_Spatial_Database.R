

############################################################################################################
#Readme:
############################################################################################################
#R code for construct the Spatial Database 
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
  #Data download from internet for CHELSA
  #Database GLOBIO4

#Data output:
  #/0_Construction_of_the_Spatial_Database/Scenario_current.RData
  #/0_Construction_of_the_Spatial_Database/data_matrix.RData #We load the data of coordinates and ID pixels
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP26.RData
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP1_RCP26.RData We load the vector of water
  #/0_Construction_of_the_Spatial_Database/vec_gaus_pr_sel_natural_SSP1_RCP26.Rdata # We load the vector of Frag           
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP60.RData
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP3_RCP70.RData
  #/0_Construction_of_the_Spatial_Database/vec_gaus_pr_sel_natural_SSP3_RCP70.Rdata
  #/0_Construction_of_the_Spatial_Database/Scenario_RCP85.RData
  #/0_Construction_of_the_Spatial_Database/vec_water_future_2050_SSP5_RCP85.RData

##############################################################################################                
#Schema
##############################################################################################  
#0_Construction_of_the_Spatial_Database

#0.1 Download/projection/creation of variables
    #0.1.1 CHELSA dataset  
      #0.1.1.1  CHELSA dataset CURRENT 
        #0.1.1.1.1 Data of temperature 
        #0.1.1.1.2 Data of precipitation   
        #0.1.1.1.3 Data of bioclimatic    
        #0.1.1.1.4 Data of tmax   
        #0.1.1.1.5 Data of tmin   
      #0.1.1.2  CHELSA dataset FUTURE 2040  
        #0.1.1.2.1 Data of temperature    
        #0.1.1.2.2 Data of precipitation   
        #0.1.1.2.3 Data of bioclimatic: 
        #0.1.1.2.4 Data of tmax    
        #0.1.1.2.5 Data of tmin   
      #0.1.1.3  CHELSA dataset FUTURE 2060  
        #0.1.1.3.1 Data of temperature:   
        #0.1.1.3.2 Data of precipitation:   
        #0.1.1.3.3 Data of bioclimatic 
        #0.1.1.3.4 Data of tmax    
        #0.1.1.3.5 Data of tmin:   
      #0.1.1.4 Read and transform Climatic database CHELSA  
    #0.1.2 LAND USE DATA
        #0.1.2.1 Scenario Current
            #0.1.3.1.1 For our category 1
            #0.1.3.1.2 For our category 2
            #0.1.3.1.3 For our category 3
            #0.1.3.1.4 For our category 4
            #0.1.3.1.5 For our category 5
            #0.1.3.1.6 For our category 6
            #0.1.3.1.7 For our category 7
            #0.1.3.1.8 For our category 8
            #0.1.3.1.9 For our category 9
        #0.1.3.2 Scenario RCP 26
            #0.1.3.2.1 For our category 1
            #0.1.3.2.2 For our category 2
            #0.1.3.2.3 For our category 3
            #0.1.3.2.4 For our category 4
            #0.1.3.2.5 For our category 5
            #0.1.3.2.6 For our category 6
            #0.1.3.2.7 For our category 7
            #0.1.3.2.8 For our category 8
            #0.1.3.2.9 For our category 9
        #0.1.3.3 Scenario RCP 60            
            #0.1.3.3.1 For our category 1
            #0.1.3.3.2 For our category 2
            #0.1.3.3.3 For our category 3
            #0.1.3.3.4 For our category 4
            #0.1.3.3.5 For our category 5
            #0.1.3.3.6 For our category 6
            #0.1.3.3.7 For our category 7
            #0.1.3.3.8 For our category 8
            #0.1.3.3.9 For our category 9
        #0.1.3.4 Scenario RCP 85
            #0.1.3.4.1 For our category 1
            #0.1.3.4.2 For our category 2
            #0.1.3.4.3 For our category 3
            #0.1.3.4.4 For our category 4
            #0.1.3.4.5 For our category 5
            #0.1.3.4.6 For our category 6
            #0.1.3.4.7 For our category 7
            #0.1.3.4.8 For our category 8
            #0.1.3.4.9 For our category 9
#0.2 Preparing the scenarios
    #0.2.1 Current land use
    #0.2.2 Current climate data 
    #0.2.3 Future climate data
        #0.2.3.1 We select the colums for create three DF one for each RCP
            #0.2.3.1.1 RCP26
            #0.2.3.1.2 RCP45
            #0.2.3.1.3 RCP60
                #0.2.3.1.3.1 RCP60 Part 1
                #0.2.3.1.3.2 RCP60 Part 2
            #0.2.3.1.4 RCP85
        #0.2.3.2 We load the created DF for future climate and we change the colums names for that they have the same names as current data (mandatory for do the predicitons for future with the models)   
    #0.2.4 Future land use data 
        #0.2.4.1 Future land use data SSP1_RCP26 
        #0.2.4.2 Future land use data SSP3_RCP70
        #0.2.4.3 Future land use data SSP5_RCP85
#0.3 We prepare the data by scenarios for BIOMOD analysis 
        #0.3.1 Scenario Current data
        #0.3.2 Future data
            #0.3.2.1 Scenario RCP26
            #0.3.2.2 Scenario RCP60
            #0.3.2.3 Scenario RCP85
        #0.3.3 Obtaintion of coordinates of our raster cells that we are going to use to calculate the buffer around presences and obtain pseudoabsences

#####################################################################################################
#####################################################################################################
#0.1 Download/projection/creation of variables
#####################################################################################################
#####################################################################################################

  #####################################################################################################
  #0.1.1 CHELSA dataset  
  #####################################################################################################
    rm(list=ls()) 
    setwd("D:/Project_Name/Spatial_Database/R_analysis")
    library(RCurl)
    library(XML)
    library(rlist)
    #0.1.1.1  CHELSA dataset CURRENT  ####
    #0.1.1.1.1 Data of temperature    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/temp/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model IPSL-CM5A-MR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*CHELSA_temp10*", clim) ]
       k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                  file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/temp/",file_to_download, sep=""))
                }
              
               
    #0.1.1.1.2 Data of precipitation   
        base_url<-"https://www.wsl.ch/lud/chelsa/data/climatologies/prec/"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*CHELSA_prec*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
               file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/climatologies/prec/",file_to_download, sep=""))
                }
    
        
    #0.1.1.1.3 Data of bioclimatic    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/bioclim/integer/"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*CHELSA_bio*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/bioclim/integer/",file_to_download, sep=""))
                }
     
           
        
    #0.1.1.1.4 Data of tmax   
        base_url<-"https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmax/"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*CHELSA_tmax10*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmax/",file_to_download, sep=""))
                }
          
      
           
    #0.1.1.1.5 Data of tmin    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmin/"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*CHELSA_tmin10*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmin/",file_to_download, sep=""))
                }
          
      
#0.1.1.2  CHELSA dataset FUTURE 2040  
    #0.1.1.2.1 Data of temperature    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/temp/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model IPSL-CM5A-MR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
       k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                  file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/temp/",file_to_download, sep=""))
                }
              
               
    #0.1.1.2.2 Data of precipitation    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/prec/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
               file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/prec/",file_to_download, sep=""))
                }
    
        
      #0.1.1.2.3 Data of bioclimatic:    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/bioclim/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/bioclim/",file_to_download, sep=""))
                }
     
           
        
      #0.1.1.2.4 Data of tmax    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/tmax/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/tmax/",file_to_download, sep=""))
                }
          
      
           
      #0.1.1.2.5 Data of tmin    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/tmin/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2041-2060/tmin/",file_to_download, sep=""))
                }
          
#0.1.1.3  CHELSA dataset FUTURE 2060  
    #0.1.1.3.1 Data of temperature:    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/temp/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model IPSL-CM5A-MR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
       k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                  file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/temp/",file_to_download, sep=""))
                }
              
               
    #0.1.1.3.2 Data of precipitation:    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/prec/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
               file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/prec/",file_to_download, sep=""))
                }
    
        
      #0.1.1.3.3 Data of bioclimatic    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/bioclim/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/bioclim/",file_to_download, sep=""))
                }
     
           
        
      #0.1.1.3.4 Data of tmax    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/tmax/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/tmax/",file_to_download, sep=""))
                }
          
      
           
      #0.1.1.3.5 Data of tmin:    
        base_url<-"https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/tmin/?C=D;O=A"
        base_html<-getURLContent(base_url)[[1]]
        links<-strsplit(base_html,"a href=")[[1]]
        
        get_data_url<-function(s) {
            u_split1<-strsplit(s,"/")[[1]][1]
            u_split2<-strsplit(u_split1,'\\"')[[1]][2]
            ifelse(grep("[[:upper:]]",u_split2)==1 & length(strsplit(u_split2,"#")[[1]])<2,return(u_split2),return(NA))
        }
        
        clim<-unlist(lapply(links,get_data_url))
        clim<-clim[which(is.na(clim)==FALSE)]
        
        #We are going to use the model MPI-ESM-LR (clim[1146:1181) and we are gonig to select the rcp45
        clim_sel<-clim[grep("*IPSL-CM5A-MR*", clim) ]
           k<-length(clim_sel)
         
             #Indica el n?mero de mapas que quieres ver cambiando 1:35 en el loop
                for (i in 1:k){
                file_to_download<-clim_sel[i]
                  shell.exec(paste("https://www.wsl.ch/lud/chelsa/data/cmip5/2061-2080/tmin/",file_to_download, sep=""))
                }         
           
##############################################################################################################################################          
#0.1.1.4 Read and transform Climatic database CHELSA  
##############################################################################################################################################          
    rm(list=ls()) 
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
   Original_folders<-c(        
   "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Original",        
   "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Precipitation/Original",        
   "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Temperature/Original",        
   "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Tmax/Original",        
   "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Tmin/Original",
   
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Precipitation/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Temperature/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Tmax/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Tmin/Original",  
   
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2060/Bioclimatic/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2060/Precipitation/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2060/Temperature/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2060/Tmax/Original",  
   "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2060/Tmin/Original")  

    DF_spatial_data<-c()
    data_all_folders<-c()
    files_list<-c()

    #We select the folder with "i"
    i=6
   
    
  #for (i in 1:15){#1:15, for some reason we do in groups of 5
    folder_in_loop<-Original_folders[i]
    
      list_files<-list.files(path=folder_in_loop)
    
      #######################################################################################################################################
      #######################################################################################################################################
      #THIS PART IS ONLY FOR WHEN WE WANT TO SPLIT A FOLDER TO RUN THE FIRST PART AND CLOSE R AND OPEN AGAIN R AND RUN THE SECOND PART
      #For folders with a lot of data (e.g. Bioclimatic data for future predictions) we divide the folder adding the constant k which
      #is the number of files that we are going to take in the fist loop
       
        #If we use "k" 
           k=50  #k is the number of files that we are going to take in the fist loop
          
              #First part with k:
                  #list_files<-list_files[1:k]
              #Second part with k: debemos anular la linea ID_original<-c(1:length_list_files) #This is 8 lines down
                  length_list_files_initial<-length(list_files)
                  list_files<-list_files[(k+1):(length_list_files_initial)]
                  ID_original<-c((k+1):(length_list_files_initial))
          #######################################################################################################################################
      #######################################################################################################################################
   

    length_list_files<-length(list_files)
    #ID_original<-c(1:length_list_files)# Debemos anular/quitar esta linea para Second part with k
    list_files_ID<-cbind(list_files, ID_original)
    DF_list_files_ID<-as.data.frame(list_files_ID)
    data_by_folder<-cbind(folder_in_loop,length_list_files)
    data_all_folders<-rbind(data_all_folders,data_by_folder)

     for (p in 1:length_list_files){ 
        file_in_loop<-list_files[p]
        r1 <- raster(paste(folder_in_loop,"/",file_in_loop,sep="")) 

        crs(r1) <- CRS('+init=EPSG:4326')

        #We are going to cut the map for reduce to Europe
        e <- extent(-20,105,15,75) 
           #sub_europe_lat<-subset(sub_year,decimalLongitude>-10&decimalLongitude<50)
            #            sub_europe_lon<-subset(sub_europe_lat,decimalLatitude>28&decimalLatitude<71.5)
                     
         m1_crop <- crop(r1,e) 
         #X11()
         #plot(m1_crop)
        Sys.sleep(5)  

        #Raster of reference: 
        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)

            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
        pr <- projectRaster(m1_crop, new_raster, method='bilinear')
         # X11()
         #plot(pr)
        
        #e2<- extent(-1800000,1800000,700000,4800000) 
        e2<- extent(-2800000,3100000,-800000,4800000) 
        
        m1_crope2 <- crop(pr,e2)
        #  X11()
         #plot(m1_crope2)
        #We save the raster:
          folder_root<-strsplit(folder_in_loop, "/Original")
          folder_to_save<-paste(folder_root,"/Projected/",sep="")
          
          files_list<-cbind(files_list,file_in_loop)
          #We are going to create a new folder if it does not exist
          if (file.exists(folder_to_save)){
           assign(file_in_loop,pr)
                save(list=file_in_loop, file=paste(folder_to_save, file_in_loop, ".Rdata",sep=""))         
            
            
} else {
    dir.create(file.path(folder_root, "Projected"))
 assign(file_in_loop,pr)
                save(list=file_in_loop, file=paste(folder_to_save, file_in_loop, ".Rdata",sep=""))       
}
                Sys.sleep(10)  
  
         vec_pr<-extract(pr,e2)
           
            DF_spatial_data<-cbind(DF_spatial_data,vec_pr) 
          
     }
  #}  
#1_5
colnames(DF_spatial_data)<-files_list
str(DF_spatial_data)
head(DF_spatial_data)
str(files_list)
DF_spatial_data6_k2<-DF_spatial_data
save(DF_spatial_data6_k2, file=paste(folder_to_save, "DF_spatial_data6_k2.Rdata",sep=""))


data_all_folders6_k2<-data_all_folders
save(data_all_folders6_k2, file=paste(folder_to_save, "data_all_folders6_k2.RData",sep=""))

files_list6_k2<-files_list
save(files_list6_k2, file=paste(folder_to_save, "files_list6_k2.RData",sep=""))


#The data that we are creating. Some saved files are corrected ("Corrected" written at the right) and other not corrected. 
#Corrected files are with the correcttion of the pixel dimensions. 
#For corrected the not corrected run the above script for the speceified folder with i values in the scrpt:
DF_spatial_data_1.RData# Corrected and saved in "folder_to_save" which means in the folder of current "G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Projected/"
DF_spatial_data_2.RData
DF_spatial_data_3.RData
DF_spatial_data_4.RData
DF_spatial_data_5.RData
DF_spatial_data_6_k1.RData# Corrected and saved in "folder_to_save" which means in the folder of 2040 projection "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/
DF_spatial_data_6_k2.RData# Corrected and saved in "folder_to_save" which means in the folder of 2040 projection "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/
DF_spatial_data_7.RData
DF_spatial_data_8.RData
DF_spatial_data_9.RData
DF_spatial_data_10.RData
DF_spatial_data_11_k1.RData
DF_spatial_data_11_k2.RData
DF_spatial_data_12.RData
DF_spatial_data_13.RData
DF_spatial_data_14.RData
DF_spatial_data_15.RData

data_all_folders_1.RData#Corrected
data_all_folders_2.RData
data_all_folders_3.RData
data_all_folders_4.RData
data_all_folders_5.RData
data_all_folders_6_k1.RData
data_all_folders_6_k2.RData
data_all_folders_7.RData
data_all_folders_8.RData
data_all_folders_9.RData
data_all_folders_10.RData
data_all_folders_11_k1.RData
data_all_folders_11_k2.RData
data_all_folders_12.RData
data_all_folders_13.RData
data_all_folders_14.RData
data_all_folders_15.RData

files_list_1.RData#Corrected
files_list_2.RData
files_list_3.RData
files_list_4.RData
files_list_5.RData
files_list_6_k1.RData
files_list_6_k2.RData
files_list_7.RData
files_list_8.RData
files_list_9.RData
files_list_10.RData
files_list_11_k1.RData
files_list_11_k2.RData
files_list_12.RData
files_list_13.RData
files_list_14.RData
files_list_15.RData
 #We exported the raster to ascii format. The we can convert it in Arcgis 9.4 using converts form Ascii to raster
#writeRaster(pr, 'pr.img',overwrite=TRUE) 

  ############################################################################################################################### 
  #0.1.2 LAND USE DATA
  ##############################################################################################################################################  
   #For our category 1
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
            e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
    # We are going to select a category for calculate the % of this category in bigger pixels
            #We extract the values of the raster to a vector:
            vector_val_select_1<-values(m1_crop)

             #We sustitute all values that are not our category by 0
              values <- c(0,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,200,201,202,210,220)
              for(i in values) vector_val_select_1[vector_val_select_1== i] <- 0
              values_category <- c(1,190)
              for(i in values_category) vector_val_select_1[vector_val_select_1== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_1<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_1)<-vector_val_select_1

            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_1_agreg<-aggregate(m1_crop_sel_1, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_1 <- projectRaster(m1_crop_sel_1_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_1)

        #We save the raster:
            save(pr_sel_1,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_1.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_1_current<-values(pr_sel_1)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_1_current, file='vec_LC_cat_1_current.RData')

        
 ##############################################################################################################################################  
   #For our category 2
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
            #We extract the values of the raster to a vector:
            vector_val_select_2<-values(m1_crop)

            #We sustitute all values that are not our category by 0
              values <- c(0,1,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_2[vector_val_select_2== i] <- 0
              values_category <- c(2,230,231,232,3)
              for(i in values_category) vector_val_select_2[vector_val_select_2== i] <- 1
              
            #We duplicate the raster with the exttension croped
            m1_crop_sel_2<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_2)<-vector_val_select_2
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_2_agreg<-aggregate(m1_crop_sel_2, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_2 <- projectRaster(m1_crop_sel_2_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_2)

        #We save the raster:
            save(pr_sel_2,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_2.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_2_current<-values(pr_sel_2)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_2_current, file='vec_LC_cat_2_current.RData')
        
          
##############################################################################################################################################  
   #For our category 3
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
             #We extract the values of the raster to a vector:
            vector_val_select_3<-values(m1_crop)
           
             #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_3[vector_val_select_3== i] <- 0
              values_category <- c(50,60,61,62)
              for(i in values_category) vector_val_select_3[vector_val_select_3== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_3<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_3 
            values(m1_crop_sel_3)<-vector_val_select_3
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_3_agreg<-aggregate(m1_crop_sel_3, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_3 <- projectRaster(m1_crop_sel_3_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_3)
        #We save the raster:
            save(pr_sel_3,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_3.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_3_current<-values(pr_sel_3)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_3_current, file='vec_LC_cat_3_current.RData')


 ##############################################################################################################################################  
   #For our category 4
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 

            #We extract the values of the raster to a vector:
            vector_val_select_4<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_4[vector_val_select_4== i] <- 0
              values_category <- c(70,71,72,80,81,82,90)
              for(i in values_category) vector_val_select_4[vector_val_select_4== i] <- 1
     
            #We duplicate the raster with the exttension croped
            m1_crop_sel_4<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_4 
            values(m1_crop_sel_4)<-vector_val_select_4
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_4_agreg<-aggregate(m1_crop_sel_4, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_4 <- projectRaster(m1_crop_sel_4_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_4)
            #We save the raster:
            save(pr_sel_4,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_4.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_4_current<-values(pr_sel_4)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_4_current, file='vec_LC_cat_4_current.RData')

##############################################################################################################################################  
   #For our category 5
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(100,110,120,121,122)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_5 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

            X11()
            plot(pr_sel_5)

            #We save the raster:
            save(pr_sel_5,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_5.Rdata")
    
            #We extract the values of the raster to a vector:
            vec_LC_cat_5_current<-values(pr_sel_5)
            #We save the vector of the category selected:
              save(vec_LC_cat_5_current, file='vec_LC_cat_5_current.RData')
    

##############################################################################################################################################  
   #For our category 6
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_6<-values(m1_crop)

                
                
              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,150,151,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_6[vector_val_select_6== i] <- 0
              values_category <- c(140)
              for(i in values_category) vector_val_select_6[vector_val_select_6== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_6<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_6 
            values(m1_crop_sel_6)<-vector_val_select_6
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_6_agreg<-aggregate(m1_crop_sel_6, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_6 <- projectRaster(m1_crop_sel_6_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_6)
            #We save the raster:
            save(pr_sel_6,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_6.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_6_current<-values(pr_sel_6)
            #We save the vector of the category selected:
              save(vec_LC_cat_6_current, file='vec_LC_cat_6_current.RData')
    
     ##############################################################################################################################################  
   #For our category 7
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_7<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_7[vector_val_select_7== i] <- 0
              values_category <- c(150,152,153)
              for(i in values_category) vector_val_select_7[vector_val_select_7== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_7<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_7)<-vector_val_select_7
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_7_agreg<-aggregate(m1_crop_sel_7, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_7 <- projectRaster(m1_crop_sel_7_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_7)
            #We save the raster:
            save(pr_sel_7,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_7.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_7_current<-values(pr_sel_7)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_7_current, file='vec_LC_cat_7_current.RData')
  
   ##############################################################################################################################################  
   #For our category 8
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_8<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_8[vector_val_select_8== i] <- 0
              values_category <- c(4)
              for(i in values_category) vector_val_select_8[vector_val_select_8== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_8<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_8)<-vector_val_select_8
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_8_agreg<-aggregate(m1_crop_sel_8, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_8 <- projectRaster(m1_crop_sel_8_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_8)
            #We save the raster:
            save(pr_sel_8,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_8.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_8_current<-values(pr_sel_8)
                     
            #We save the vector of the category selected:
            save(vec_LC_cat_8_current, file='vec_LC_cat_8_current.RData')
              
            ################################################################################################################################################  
            #0.1.3.1.9 For our category 9
                rm(list=ls()) 
                ##Establecemos el fichero de trabajo
                setwd("G:/Project_Name/Spatial_Database/R_analysis")
                library(rgdal) 
                library(raster) 
                #We select the original raster of land use that we want to transform:
                r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
                #We are going to cut/crop the map for reduce to Europe
                #We stablish an extension for crop the map:
                e <- extent(-20,105,15,75)#Extension que utlizamos de normal 
                #We use the extension to crop the original map    
                m1_crop <- crop(r1,e) 
                
                #We extract the values of the raster to a vector:
                vector_val_select_9<-values(m1_crop)
            
                #We sustitute all values that are not our category by 0
                values <- c(0,1,2,4,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,210,220)
                for(i in values) vector_val_select_9[vector_val_select_9== i] <- 0
                values_category <- c(200,201,202)
                for(i in values_category) vector_val_select_9[vector_val_select_9== i] <- 1
                        
                        #We duplicate the raster with the exttension croped
                        m1_crop_sel_9<-m1_crop
                        #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
                        values(m1_crop_sel_9)<-vector_val_select_9
                        #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
                        m1_crop_sel_9_agreg<-aggregate(m1_crop_sel_9, fact=3, fun=mean, expand=TRUE)
                       
                        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                        #We are going to use Europe Albers Equal Area Conic  
                        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                        #We defined the projection of our raster:
                        projection(new_raster) <- newproj
                        
                        #We project our data to the raster created with the new projection
                        pr_sel_9 <- projectRaster(m1_crop_sel_9_agreg, new_raster, method='bilinear')
                        X11()
                        plot(pr_sel_9)
                        #We save the raster:
                        save(pr_sel_9,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_9.Rdata")
            
                        #We extract the values of the raster to a vector:
                        vec_LC_cat_9_current<-values(pr_sel_9)
                                 
                        #We save the vector of the category selected:
                          save(vec_LC_cat_9_current, file='vec_LC_cat_9_current.RData')
             
            ################################################################################################################################################  
            #0.1.3.1.10 For our fragmentation variables
                rm(list=ls()) 
                ##Establecemos el fichero de trabajo
                setwd("G:/Project_Name/Spatial_Database/R_analysis")
                library(rgdal) 
                library(raster)
                library(spatialEco)
                library(readxl)
                #install.packages("spatialEco")
                
                #We select the original raster of land use that we want to transform:
                r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
                #We are going to cut/crop the map for reduce to Europe
                #We stablish an extension for crop the map:
                e <- extent(-20,105,15,75)#Extension que utlizamos de normal 
                #We use the extension to crop the original map    
                m1_crop <- crop(r1,e) 
                
                #We extract the values of the raster to a vector:
                vector_val_select_natural<-values(m1_crop)
            
                
                use_categories<-read_excel("H:/G/Project_Name/Spatial_Database/R_analysis/Land_Use_classification/land_use_categories.xlsx", sheet = 1)
                df_use_categories<-as.data.frame(use_categories)

                df_use_categories_natural<-subset(df_use_categories,Natural_Human_used == "Natural")
                list_natural_cat<-unique(df_use_categories_natural$ID_CCI_LC_classification)
             
                df_use_categories_no_natural<-subset(df_use_categories,Natural_Human_used != "Natural")
                list_no_natural_cat<-unique(df_use_categories_no_natural$ID_CCI_LC_classification)

                   
                #We sustitute all values that are not our category by 0
                values <- list_no_natural_cat
                for(i in values) vector_val_select_natural[vector_val_select_natural== i] <- 0
                values_category <- list_natural_cat
                for(i in values_category) vector_val_select_natural[vector_val_select_natural== i] <- 1
                        
                        #We duplicate the raster with the exttension croped
                        m1_crop_sel_natural<-m1_crop
                        #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
                        values(m1_crop_sel_natural)<-vector_val_select_natural
                        
                        
                        #Para distance to the border lo calcularemos sin agregar, es decir,
                        # exportamos m1_crop_sel_natural a idrisi y calculamos la distancia hacia adentro
                        
                        
                        
                        #Para el kernel lo calcularemos sin agregar
                        
                        
                        
                        
                        #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
                        m1_crop_sel_natural_agreg<-aggregate(m1_crop_sel_natural, fact=3, fun=mean, expand=TRUE)
                       
                        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                        #We are going to use Europe Albers Equal Area Conic  
                        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                        #We defined the projection of our raster:
                        projection(new_raster) <- newproj
                        
                        #We project our data to the raster created with the new projection
                        pr_sel_natural <- projectRaster(m1_crop_sel_natural_agreg, new_raster, method='bilinear')
                        X11()
                        plot(pr_sel_natural)
                        #We save the raster:
                        save(pr_sel_natural,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_natural.Rdata")
                        writeRaster(pr_sel_natural, 'G:/Project_Name/Spatial_Database/R_analysis/pr_sel_natural.img',overwrite=TRUE) 
                        #We extract the values of the raster to a vector:
                        vec_LC_cat_natural_current<-values(pr_sel_natural)
                                 
                        #We save the vector of the category selected:
                        save(vec_LC_cat_natural_current, file='vec_LC_cat_natural_current.RData')
                        gaus_pr_sel_natural<-raster.gaussian.smooth(pr_sel_natural, sigma = 2, n = 11, type = mean)
                        writeRaster(gaus_pr_sel_natural, 'G:/Project_Name/Spatial_Database/R_analysis/gaus_pr_sel_natural.img',overwrite=TRUE) 
                        
                        load("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/gaus_pr_sel_natural.Rdata")
                             #We extract the values of the raster to a vector:
                        vec_gaus_pr_sel_natural_current<-values(gaus_pr_sel_natural)
                        save(vec_gaus_pr_sel_natural_current,file="vec_gaus_pr_sel_natural_current.Rdata")

                        
##############################################################################################################################################  
   #For our category cober (Is the percentage of bushses and forests)
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
    #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
        #r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,130,140,150,151,152,153,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            map_cober <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

        #We save the raster:
           #save(map_cober,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/map_cober.Rdata")

        #We extract the values of the raster to a vector:
            vec_cober_current<-values(map_cober)
                     
        #We save the vector of the category selected:
          save(vec_cober_current, file='vec_cober_current.RData')
   

############################################################################################################
   #For category of water
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180,130,140,150,151,152,153,190,200,201,202,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(210)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1

            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #We reclass human to integer values to calculate water areas
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg >= 0.5 & m1_crop_sel_5_agreg <= 1 ] <- 1
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg > 0 & m1_crop_sel_5_agreg < 0.5 ] <- 0
            
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_water_current <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='ngb')


            #We save the raster:
            save(pr_water_current,file="pr_water_current.Rdata")
            #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            #writeRaster(pr_water, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_water_current<-values(pr_water_current)
            #We save the vector of the category selected:
            save(vec_water_current, file='vec_water_current.RData')
          

############################################################################################################################### 
#0.1.3.2 Scenario RCP 26
    #0.1.3.2.1 For our category 1
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
    r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
    #r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 

    #We are going to cut/crop the map for reduce to Europe
    #We stablish an extension for crop the map:
    e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

    #We use the extension to crop the original map    
    m1_crop <- crop(r1,e) 
  
    # We are going to select a category for calculate the % of this category in bigger pixels
    #We extract the values of the raster to a vector:
    vector_val_select_1<-values(m1_crop)

    #We sustitute all values that are not our category by 0
    values <- c(0,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,200,201,202,210,220)
              for(i in values) vector_val_select_1[vector_val_select_1== i] <- 0
              values_category <- c(1,190)
              for(i in values_category) vector_val_select_1[vector_val_select_1== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_1<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_1)<-vector_val_select_1

            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_1_agreg<-aggregate(m1_crop_sel_1, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_1 <- projectRaster(m1_crop_sel_1_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_1)

        #We save the raster:
            save(pr_sel_1,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_1.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_1_future_2050_SSP1_RCP26<-values(pr_sel_1)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_1_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_1_future_2050_SSP1_RCP26.RData')

        
 ##############################################################################################################################################  
   #0.1.4.3.2 For our category 2
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
            #We extract the values of the raster to a vector:
            vector_val_select_2<-values(m1_crop)

            #We sustitute all values that are not our category by 0
              values <- c(0,1,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_2[vector_val_select_2== i] <- 0
              values_category <- c(2,230,231,232,3)
              for(i in values_category) vector_val_select_2[vector_val_select_2== i] <- 1
              
            #We duplicate the raster with the exttension croped
            m1_crop_sel_2<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_2)<-vector_val_select_2
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_2_agreg<-aggregate(m1_crop_sel_2, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_2 <- projectRaster(m1_crop_sel_2_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_2)

        #We save the raster:
            save(pr_sel_2,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_2.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_2_future_2050_SSP1_RCP26<-values(pr_sel_2)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_2_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_2_future_2050_SSP1_RCP26.RData')
        
          
##############################################################################################################################################  
   #0.1.4.3.3 For our category 3
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
             #We extract the values of the raster to a vector:
            vector_val_select_3<-values(m1_crop)
           
             #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_3[vector_val_select_3== i] <- 0
              values_category <- c(50,60,61,62)
              for(i in values_category) vector_val_select_3[vector_val_select_3== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_3<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_3 
            values(m1_crop_sel_3)<-vector_val_select_3
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_3_agreg<-aggregate(m1_crop_sel_3, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_3 <- projectRaster(m1_crop_sel_3_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_3)
        #We save the raster:
            save(pr_sel_3,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_3.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_3_future_2050_SSP1_RCP26<-values(pr_sel_3)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_3_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_3_future_2050_SSP1_RCP26.RData')


 ##############################################################################################################################################  
   #0.1.4.3.4 For our category 4
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 

            #We extract the values of the raster to a vector:
            vector_val_select_4<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_4[vector_val_select_4== i] <- 0
              values_category <- c(70,71,72,80,81,82,90)
              for(i in values_category) vector_val_select_4[vector_val_select_4== i] <- 1
     
            #We duplicate the raster with the exttension croped
            m1_crop_sel_4<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_4 
            values(m1_crop_sel_4)<-vector_val_select_4
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_4_agreg<-aggregate(m1_crop_sel_4, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_4 <- projectRaster(m1_crop_sel_4_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_4)
            #We save the raster:
            save(pr_sel_4,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_4.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_4_future_2050_SSP1_RCP26<-values(pr_sel_4)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_4_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_4_future_2050_SSP1_RCP26.RData')

##############################################################################################################################################  
   #0.1.4.3.5 For our category 5
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(100,110,120,121,122)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_5 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

            X11()
            plot(pr_sel_5)

            #We save the raster:
            save(pr_sel_5,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_5.Rdata")
    
            #We extract the values of the raster to a vector:
            vec_LC_cat_5_future_2050_SSP1_RCP26<-values(pr_sel_5)
            #We save the vector of the category selected:
              save(vec_LC_cat_5_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_5_future_2050_SSP1_RCP26.RData')
    

##############################################################################################################################################  
   #0.1.4.3.6 For our category 6
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_6<-values(m1_crop)

                
                
              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,150,151,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_6[vector_val_select_6== i] <- 0
              values_category <- c(140)
              for(i in values_category) vector_val_select_6[vector_val_select_6== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_6<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_6 
            values(m1_crop_sel_6)<-vector_val_select_6
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_6_agreg<-aggregate(m1_crop_sel_6, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_6 <- projectRaster(m1_crop_sel_6_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_6)
            #We save the raster:
            save(pr_sel_6,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_6.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_6_future_2050_SSP1_RCP26<-values(pr_sel_6)
            #We save the vector of the category selected:
              save(vec_LC_cat_6_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_6_future_2050_SSP1_RCP26.RData')
    
     ##############################################################################################################################################  
   #0.1.4.3.7 For our category 7
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_7<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_7[vector_val_select_7== i] <- 0
              values_category <- c(150,152,153)
              for(i in values_category) vector_val_select_7[vector_val_select_7== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_7<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_7)<-vector_val_select_7
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_7_agreg<-aggregate(m1_crop_sel_7, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_7 <- projectRaster(m1_crop_sel_7_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_7)
            #We save the raster:
            save(pr_sel_7,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_7.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_7_future_2050_SSP1_RCP26<-values(pr_sel_7)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_7_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_7_future_2050_SSP1_RCP26.RData')
  
   ##############################################################################################################################################  
   #0.1.4.3.8 For our category 8
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_8<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_8[vector_val_select_8== i] <- 0
              values_category <- c(4)
              for(i in values_category) vector_val_select_8[vector_val_select_8== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_8<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_8)<-vector_val_select_8
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_8_agreg<-aggregate(m1_crop_sel_8, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_8 <- projectRaster(m1_crop_sel_8_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_8)
            #We save the raster:
            save(pr_sel_8,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_8.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_8_future_2050_SSP1_RCP26<-values(pr_sel_8)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_8_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_8_future_2050_SSP1_RCP26.RData')
              
   ##############################################################################################################################################  
   #0.1.4.3.9 For our category 9
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_9<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,4,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,210,220)
               for(i in values) vector_val_select_9[vector_val_select_9== i] <- 0
              values_category <- c(200,201,202)
              for(i in values_category) vector_val_select_9[vector_val_select_9== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_9<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_9)<-vector_val_select_9
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_9_agreg<-aggregate(m1_crop_sel_9, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_9 <- projectRaster(m1_crop_sel_9_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_9)
            #We save the raster:
            save(pr_sel_9,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_9.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_9_future_2050_SSP1_RCP26<-values(pr_sel_9)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_9_future_2050_SSP1_RCP26, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_9_future_2050_SSP1_RCP26.RData')

            ################################################################################################################################################  
            #0.1.3.1.10 For our fragmentation variables
                rm(list=ls()) 
                ##Establecemos el fichero de trabajo
                setwd("H:/G/Project_Name/Spatial_Database/R_analysis")
                library(rgdal) 
                library(raster)
                library(spatialEco)
                library(readxl)
                #install.packages("spatialEco")
                
                #We select the original raster of land use that we want to transform:
                r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
                #We are going to cut/crop the map for reduce to Europe
                #We stablish an extension for crop the map:
                e <- extent(-20,105,15,75)#Extension que utlizamos de normal 
                #We use the extension to crop the original map    
                m1_crop <- crop(r1,e) 
                
                #We extract the values of the raster to a vector:
                vector_val_select_natural<-values(m1_crop)
            
                
                use_categories<-read_excel("H:/G/Project_Name/Spatial_Database/R_analysis/Land_Use_classification/land_use_categories.xlsx", sheet = 1)
                df_use_categories<-as.data.frame(use_categories)

                df_use_categories_natural<-subset(df_use_categories,Natural_Human_used == "Natural")
                list_natural_cat<-unique(df_use_categories_natural$ID_CCI_LC_classification)
             
                df_use_categories_no_natural<-subset(df_use_categories,Natural_Human_used != "Natural")
                list_no_natural_cat<-unique(df_use_categories_no_natural$ID_CCI_LC_classification)

                   
                #We sustitute all values that are not our category by 0
                values <- list_no_natural_cat
                for(i in values) vector_val_select_natural[vector_val_select_natural== i] <- 0
                values_category <- list_natural_cat
                for(i in values_category) vector_val_select_natural[vector_val_select_natural== i] <- 1
                        
                        #We duplicate the raster with the exttension croped
                        m1_crop_sel_natural<-m1_crop
                        #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
                        values(m1_crop_sel_natural)<-vector_val_select_natural
                        #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
                        m1_crop_sel_natural_agreg<-aggregate(m1_crop_sel_natural, fact=3, fun=mean, expand=TRUE)
                       
                        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                        #We are going to use Europe Albers Equal Area Conic  
                        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                        #We defined the projection of our raster:
                        projection(new_raster) <- newproj
                        
                        #We project our data to the raster created with the new projection
                        pr_sel_natural <- projectRaster(m1_crop_sel_natural_agreg, new_raster, method='bilinear')
                        #X11()
                        #plot(pr_sel_natural)
                        #We save the raster:
                        save(pr_sel_natural,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/pr_sel_natural.Rdata")
                        #writeRaster(pr_sel_natural, 'H:/G/Project_Name/Spatial_Database/R_analysis/pr_sel_natural.img',overwrite=TRUE) 
                        #We extract the values of the raster to a vector:
                        vec_LC_cat_natural_SSP1_RCP26<-values(pr_sel_natural)
                                 
                        #We save the vector of the category selected:
                        save(vec_LC_cat_natural_SSP1_RCP26, file='H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_natural_SSP1_RCP26.RData')
                        gaus_pr_sel_natural<-raster.gaussian.smooth(pr_sel_natural, sigma = 2, n = 11, type = mean)
                        #writeRaster(gaus_pr_sel_natural, 'H:/G/Project_Name/Spatial_Database/R_analysis/gaus_pr_sel_natural.img',overwrite=TRUE) 
                        
                        #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/SSP1_RCP26/Projected/gaus_pr_sel_natural.Rdata")
                             #We extract the values of the raster to a vector:
                        vec_gaus_pr_sel_natural_SSP1_RCP26<-values(gaus_pr_sel_natural)
                        save(vec_gaus_pr_sel_natural_SSP1_RCP26,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_gaus_pr_sel_natural_SSP1_RCP26.Rdata")


##############################################################################################################################################  
   #For our category cober (Is the percentage of bushses and forests)
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
    #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
        #r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,130,140,150,151,152,153,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            map_cober <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

        #We save the raster:
          # save(map_cober,file="C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR/map_cober.Rdata")

        #We extract the values of the raster to a vector:
            vec_cober_future_2050_SSP1_RCP26<-values(map_cober)
                     
        #We save the vector of the category selected:
          save(vec_cober_future_2050_SSP1_RCP26, file='vec_cober_future_2050_SSP1_RCP26.RData')
          
          
############################################################################################################
   #For category of water
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180,130,140,150,151,152,153,190,200,201,202,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(210)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1

            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #We reclass human to integer values to calculate water areas
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg >= 0.5 & m1_crop_sel_5_agreg <= 1 ] <- 1
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg > 0 & m1_crop_sel_5_agreg < 0.5 ] <- 0
            
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_water_future_2050_SSP1_RCP26 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='ngb')


            #We save the raster:
            save(pr_water_future_2050_SSP1_RCP26,file="pr_water_future_2050_SSP1_RCP26.Rdata")
            #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            #writeRaster(pr_water, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_water_future_2050_SSP1_RCP26<-values(pr_water_future_2050_SSP1_RCP26)
            #We save the vector of the category selected:
            save(vec_water_future_2050_SSP1_RCP26, file='vec_water_future_2050_SSP1_RCP26.RData')
          
  
###############################################################################################         
#0.1.3.3 Scenario RCP 60
    #For our category 1
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        #r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
        #r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 

      #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
            e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
    # We are going to select a category for calculate the % of this category in bigger pixels
            #We extract the values of the raster to a vector:
            vector_val_select_1<-values(m1_crop)

             #We sustitute all values that are not our category by 0
              values <- c(0,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,200,201,202,210,220)
              for(i in values) vector_val_select_1[vector_val_select_1== i] <- 0
              values_category <- c(1,190)
              for(i in values_category) vector_val_select_1[vector_val_select_1== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_1<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_1)<-vector_val_select_1

            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_1_agreg<-aggregate(m1_crop_sel_1, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_1 <- projectRaster(m1_crop_sel_1_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_1)

        #We save the raster:
            save(pr_sel_1,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_1.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_1_future_2050_SSP3_RCP70<-values(pr_sel_1)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_1_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_1_future_2050_SSP3_RCP70.RData')

        
 ##############################################################################################################################################  
   #For our category 2
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
            #We extract the values of the raster to a vector:
            vector_val_select_2<-values(m1_crop)

            #We sustitute all values that are not our category by 0
              values <- c(0,1,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_2[vector_val_select_2== i] <- 0
              values_category <- c(2,230,231,232,3)
              for(i in values_category) vector_val_select_2[vector_val_select_2== i] <- 1
              
            #We duplicate the raster with the exttension croped
            m1_crop_sel_2<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_2)<-vector_val_select_2
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_2_agreg<-aggregate(m1_crop_sel_2, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_2 <- projectRaster(m1_crop_sel_2_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_2)

        #We save the raster:
            save(pr_sel_2,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_2.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_2_future_2050_SSP3_RCP70<-values(pr_sel_2)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_2_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_2_future_2050_SSP3_RCP70.RData')
        
          
##############################################################################################################################################  
   #For our category 3
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
             #We extract the values of the raster to a vector:
            vector_val_select_3<-values(m1_crop)
           
             #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_3[vector_val_select_3== i] <- 0
              values_category <- c(50,60,61,62)
              for(i in values_category) vector_val_select_3[vector_val_select_3== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_3<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_3 
            values(m1_crop_sel_3)<-vector_val_select_3
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_3_agreg<-aggregate(m1_crop_sel_3, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_3 <- projectRaster(m1_crop_sel_3_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_3)
        #We save the raster:
            save(pr_sel_3,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_3.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_3_future_2050_SSP3_RCP70<-values(pr_sel_3)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_3_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_3_future_2050_SSP3_RCP70.RData')


 ##############################################################################################################################################  
   #For our category 4
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 

            #We extract the values of the raster to a vector:
            vector_val_select_4<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_4[vector_val_select_4== i] <- 0
              values_category <- c(70,71,72,80,81,82,90)
              for(i in values_category) vector_val_select_4[vector_val_select_4== i] <- 1
     
            #We duplicate the raster with the exttension croped
            m1_crop_sel_4<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_4 
            values(m1_crop_sel_4)<-vector_val_select_4
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_4_agreg<-aggregate(m1_crop_sel_4, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_4 <- projectRaster(m1_crop_sel_4_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_4)
            #We save the raster:
            save(pr_sel_4,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_4.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_4_future_2050_SSP3_RCP70<-values(pr_sel_4)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_4_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_4_future_2050_SSP3_RCP70.RData')

##############################################################################################################################################  
   #For our category 5
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(100,110,120,121,122)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_5 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

            X11()
            plot(pr_sel_5)

            #We save the raster:
            save(pr_sel_5,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_5.Rdata")
    
            #We extract the values of the raster to a vector:
            vec_LC_cat_5_future_2050_SSP3_RCP70<-values(pr_sel_5)
            #We save the vector of the category selected:
              save(vec_LC_cat_5_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_5_future_2050_SSP3_RCP70.RData')
    

##############################################################################################################################################  
   #For our category 6
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_6<-values(m1_crop)

                
                
              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,150,151,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_6[vector_val_select_6== i] <- 0
              values_category <- c(140)
              for(i in values_category) vector_val_select_6[vector_val_select_6== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_6<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_6 
            values(m1_crop_sel_6)<-vector_val_select_6
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_6_agreg<-aggregate(m1_crop_sel_6, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_6 <- projectRaster(m1_crop_sel_6_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_6)
            #We save the raster:
            save(pr_sel_6,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_6.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_6_future_2050_SSP3_RCP70<-values(pr_sel_6)
            #We save the vector of the category selected:
              save(vec_LC_cat_6_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_6_future_2050_SSP3_RCP70.RData')
    
     ##############################################################################################################################################  
   #For our category 7
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_7<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_7[vector_val_select_7== i] <- 0
              values_category <- c(150,152,153)
              for(i in values_category) vector_val_select_7[vector_val_select_7== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_7<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_7)<-vector_val_select_7
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_7_agreg<-aggregate(m1_crop_sel_7, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_7 <- projectRaster(m1_crop_sel_7_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_7)
            #We save the raster:
            save(pr_sel_7,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_7.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_7_future_2050_SSP3_RCP70<-values(pr_sel_7)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_7_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_7_future_2050_SSP3_RCP70.RData')
  
   ##############################################################################################################################################  
   #For our category 8
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_8<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_8[vector_val_select_8== i] <- 0
              values_category <- c(4)
              for(i in values_category) vector_val_select_8[vector_val_select_8== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_8<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_8)<-vector_val_select_8
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_8_agreg<-aggregate(m1_crop_sel_8, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_8 <- projectRaster(m1_crop_sel_8_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_8)
            #We save the raster:
            save(pr_sel_8,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_8.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_8_future_2050_SSP3_RCP70<-values(pr_sel_8)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_8_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_8_future_2050_SSP3_RCP70.RData')
              
   ##############################################################################################################################################  
   #For our category 9
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_9<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,4,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,210,220)
               for(i in values) vector_val_select_9[vector_val_select_9== i] <- 0
              values_category <- c(200,201,202)
              for(i in values_category) vector_val_select_9[vector_val_select_9== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_9<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_9)<-vector_val_select_9
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_9_agreg<-aggregate(m1_crop_sel_9, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_9 <- projectRaster(m1_crop_sel_9_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_9)
            #We save the raster:
            save(pr_sel_9,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_9.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_9_future_2050_SSP3_RCP70<-values(pr_sel_9)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_9_future_2050_SSP3_RCP70, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_9_future_2050_SSP3_RCP70.RData')
 
              
            ################################################################################################################################################  
            #0.1.3.1.10 For our fragmentation variables
                rm(list=ls()) 
                ##Establecemos el fichero de trabajo
                setwd("H:/G/Project_Name/Spatial_Database/R_analysis")
                library(rgdal) 
                library(raster)
                library(spatialEco)
                library(readxl)
                #install.packages("spatialEco")
                
                #We select the original raster of land use that we want to transform:
                r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
                #We are going to cut/crop the map for reduce to Europe
                #We stablish an extension for crop the map:
                e <- extent(-20,105,15,75)#Extension que utlizamos de normal 
                #We use the extension to crop the original map    
                m1_crop <- crop(r1,e) 
                
                #We extract the values of the raster to a vector:
                vector_val_select_natural<-values(m1_crop)
            
                
                use_categories<-read_excel("H:/G/Project_Name/Spatial_Database/R_analysis/Land_Use_classification/land_use_categories.xlsx", sheet = 1)
                df_use_categories<-as.data.frame(use_categories)

                df_use_categories_natural<-subset(df_use_categories,Natural_Human_used == "Natural")
                list_natural_cat<-unique(df_use_categories_natural$ID_CCI_LC_classification)
             
                df_use_categories_no_natural<-subset(df_use_categories,Natural_Human_used != "Natural")
                list_no_natural_cat<-unique(df_use_categories_no_natural$ID_CCI_LC_classification)

                   
                #We sustitute all values that are not our category by 0
                values <- list_no_natural_cat
                for(i in values) vector_val_select_natural[vector_val_select_natural== i] <- 0
                values_category <- list_natural_cat
                for(i in values_category) vector_val_select_natural[vector_val_select_natural== i] <- 1
                        
                        #We duplicate the raster with the exttension croped
                        m1_crop_sel_natural<-m1_crop
                        #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
                        values(m1_crop_sel_natural)<-vector_val_select_natural
       
                        #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
                        m1_crop_sel_natural_agreg<-aggregate(m1_crop_sel_natural, fact=3, fun=mean, expand=TRUE)
                       
                        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                        #We are going to use Europe Albers Equal Area Conic  
                        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                        #We defined the projection of our raster:
                        projection(new_raster) <- newproj
                        
                        #We project our data to the raster created with the new projection
                        pr_sel_natural <- projectRaster(m1_crop_sel_natural_agreg, new_raster, method='bilinear')
                        #X11()
                        #plot(pr_sel_natural)
                        #We save the raster:
                        save(pr_sel_natural,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/pr_sel_natural.Rdata")
                        #writeRaster(pr_sel_natural, 'H:/G/Project_Name/Spatial_Database/R_analysis/pr_sel_natural.img',overwrite=TRUE) 
                        #We extract the values of the raster to a vector:
                        vec_LC_cat_natural_SSP3_RCP70<-values(pr_sel_natural)
                                 
                        #We save the vector of the category selected:
                        save(vec_LC_cat_natural_SSP3_RCP70, file='H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_natural_SSP3_RCP70.RData')
                        gaus_pr_sel_natural<-raster.gaussian.smooth(pr_sel_natural, sigma = 2, n = 11, type = mean)
                        #writeRaster(gaus_pr_sel_natural, 'H:/G/Project_Name/Spatial_Database/R_analysis/gaus_pr_sel_natural.img',overwrite=TRUE) 
                        
                        #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/SSP3_RCP70/Projected/gaus_pr_sel_natural.Rdata")
                             #We extract the values of the raster to a vector:
                        vec_gaus_pr_sel_natural_SSP3_RCP70<-values(gaus_pr_sel_natural)
                        save(vec_gaus_pr_sel_natural_SSP3_RCP70,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_gaus_pr_sel_natural_SSP3_RCP70.Rdata")

              
##############################################################################################################################################  
   #For our category cober (Is the percentage of bushses and forests)
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
    #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
        #r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,130,140,150,151,152,153,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            map_cober <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

        #We save the raster:
            #save(map_cober,file="H:/G/Project_Name/Spatial_Database/Future_2050/Projected/SSP3_RCP70/map_cober.Rdata")

        #We extract the values of the raster to a vector:
            vec_cober_future_2050_SSP3_RCP70<-values(map_cober)
                     
        #We save the vector of the category selected:
          save(vec_cober_future_2050_SSP3_RCP70, file='vec_cober_future_2050_SSP3_RCP70.RData')
              
          
############################################################################################################
   #For category of water
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP3_RCP70/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180,130,140,150,151,152,153,190,200,201,202,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(210)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1

            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #We reclass human to integer values to calculate water areas
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg >= 0.5 & m1_crop_sel_5_agreg <= 1 ] <- 1
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg > 0 & m1_crop_sel_5_agreg < 0.5 ] <- 0
            
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_water_future_2050_SSP3_RCP70 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='ngb')


            #We save the raster:
            save(pr_water_future_2050_SSP3_RCP70,file="pr_water_future_2050_SSP3_RCP70.Rdata")
            #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            #writeRaster(pr_water, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_water_future_2050_SSP3_RCP70<-values(pr_water_future_2050_SSP3_RCP70)
            #We save the vector of the category selected:
            save(vec_water_future_2050_SSP3_RCP70, file='vec_water_future_2050_SSP3_RCP70.RData')
          
          
          
          
          
 #0.1.3.4 Scenario RCP 85
                           
    #For our category 1
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        #r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP1_RCP26/Globio4_landuse_10sec_2050_cropint.tif") 
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
        #r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 

      #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
            e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
    # We are going to select a category for calculate the % of this category in bigger pixels
            #We extract the values of the raster to a vector:
            vector_val_select_1<-values(m1_crop)

             #We sustitute all values that are not our category by 0
              values <- c(0,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,200,201,202,210,220)
              for(i in values) vector_val_select_1[vector_val_select_1== i] <- 0
              values_category <- c(1,190)
              for(i in values_category) vector_val_select_1[vector_val_select_1== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_1<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_1)<-vector_val_select_1

            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_1_agreg<-aggregate(m1_crop_sel_1, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_1 <- projectRaster(m1_crop_sel_1_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_1)

        #We save the raster:
            save(pr_sel_1,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_1.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_1_future_2050_SSP5_RCP85<-values(pr_sel_1)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_1_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_1_future_2050_SSP5_RCP85.RData')

        
 ##############################################################################################################################################  
   #For our category 2
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
  
            #We extract the values of the raster to a vector:
            vector_val_select_2<-values(m1_crop)

            #We sustitute all values that are not our category by 0
              values <- c(0,1,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_2[vector_val_select_2== i] <- 0
              values_category <- c(2,230,231,232,3)
              for(i in values_category) vector_val_select_2[vector_val_select_2== i] <- 1
              
            #We duplicate the raster with the exttension croped
            m1_crop_sel_2<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_1 
            values(m1_crop_sel_2)<-vector_val_select_2
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_2_agreg<-aggregate(m1_crop_sel_2, fact=3, fun=mean, expand=TRUE)
                #X11()
                #plot(m1_crop_sel_agreg)

            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_2 <- projectRaster(m1_crop_sel_2_agreg, new_raster, method='bilinear')
            
       #For visualize the data, takes a lot of time:
              X11()
              plot(pr_sel_2)

        #We save the raster:
            save(pr_sel_2,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_2.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_2_future_2050_SSP5_RCP85<-values(pr_sel_2)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_2_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_2_future_2050_SSP5_RCP85.RData')
        
          
##############################################################################################################################################  
   #For our category 3
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
             #We extract the values of the raster to a vector:
            vector_val_select_3<-values(m1_crop)
           
             #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_3[vector_val_select_3== i] <- 0
              values_category <- c(50,60,61,62)
              for(i in values_category) vector_val_select_3[vector_val_select_3== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_3<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_3 
            values(m1_crop_sel_3)<-vector_val_select_3
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_3_agreg<-aggregate(m1_crop_sel_3, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_3 <- projectRaster(m1_crop_sel_3_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_3)
        #We save the raster:
            save(pr_sel_3,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_3.Rdata")

        #We extract the values of the raster to a vector:
            vec_LC_cat_3_future_2050_SSP5_RCP85<-values(pr_sel_3)
                     
        #We save the vector of the category selected:
          save(vec_LC_cat_3_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_3_future_2050_SSP5_RCP85.RData')


 ##############################################################################################################################################  
   #For our category 4
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 

    #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 

            #We extract the values of the raster to a vector:
            vector_val_select_4<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_4[vector_val_select_4== i] <- 0
              values_category <- c(70,71,72,80,81,82,90)
              for(i in values_category) vector_val_select_4[vector_val_select_4== i] <- 1
     
            #We duplicate the raster with the exttension croped
            m1_crop_sel_4<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_4 
            values(m1_crop_sel_4)<-vector_val_select_4
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_4_agreg<-aggregate(m1_crop_sel_4, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_4 <- projectRaster(m1_crop_sel_4_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_4)
            #We save the raster:
            save(pr_sel_4,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_4.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_4_future_2050_SSP5_RCP85<-values(pr_sel_4)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_4_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_4_future_2050_SSP5_RCP85.RData')

##############################################################################################################################################  
   #For our category 5
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(100,110,120,121,122)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_5 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

            X11()
            plot(pr_sel_5)

            #We save the raster:
            save(pr_sel_5,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_5.Rdata")
    
            #We extract the values of the raster to a vector:
            vec_LC_cat_5_future_2050_SSP5_RCP85<-values(pr_sel_5)
            #We save the vector of the category selected:
              save(vec_LC_cat_5_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_5_future_2050_SSP5_RCP85.RData')
    

##############################################################################################################################################  
   #For our category 6
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_6<-values(m1_crop)

                
                
              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,150,151,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_6[vector_val_select_6== i] <- 0
              values_category <- c(140)
              for(i in values_category) vector_val_select_6[vector_val_select_6== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_6<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_6 
            values(m1_crop_sel_6)<-vector_val_select_6
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_6_agreg<-aggregate(m1_crop_sel_6, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_6 <- projectRaster(m1_crop_sel_6_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_6)
            #We save the raster:
            save(pr_sel_6,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_6.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_6_future_2050_SSP5_RCP85<-values(pr_sel_6)
            #We save the vector of the category selected:
              save(vec_LC_cat_6_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_6_future_2050_SSP5_RCP85.RData')
    
     ##############################################################################################################################################  
   #For our category 7
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_7<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_7[vector_val_select_7== i] <- 0
              values_category <- c(150,152,153)
              for(i in values_category) vector_val_select_7[vector_val_select_7== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_7<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_7)<-vector_val_select_7
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_7_agreg<-aggregate(m1_crop_sel_7, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_7 <- projectRaster(m1_crop_sel_7_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_7)
            #We save the raster:
            save(pr_sel_7,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_7.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_7_future_2050_SSP5_RCP85<-values(pr_sel_7)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_7_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_7_future_2050_SSP5_RCP85.RData')
  
   ##############################################################################################################################################  
   #For our category 8
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_8<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,200,201,202,210,220)
               for(i in values) vector_val_select_8[vector_val_select_8== i] <- 0
              values_category <- c(4)
              for(i in values_category) vector_val_select_8[vector_val_select_8== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_8<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_8)<-vector_val_select_8
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_8_agreg<-aggregate(m1_crop_sel_8, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_8 <- projectRaster(m1_crop_sel_8_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_8)
            #We save the raster:
            save(pr_sel_8,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_8.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_8_future_2050_SSP5_RCP85<-values(pr_sel_8)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_8_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_8_future_2050_SSP5_RCP85.RData')
              
   ##############################################################################################################################################  
   #For our category 9
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_9<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,4,230,231,232,3,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,151,150,152,153,160,170,180,190,210,220)
               for(i in values) vector_val_select_9[vector_val_select_9== i] <- 0
              values_category <- c(200,201,202)
              for(i in values_category) vector_val_select_9[vector_val_select_9== i] <- 1
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_9<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
            values(m1_crop_sel_9)<-vector_val_select_9
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_9_agreg<-aggregate(m1_crop_sel_9, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_9 <- projectRaster(m1_crop_sel_9_agreg, new_raster, method='bilinear')
            X11()
            plot(pr_sel_9)
            #We save the raster:
            save(pr_sel_9,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_9.Rdata")

            #We extract the values of the raster to a vector:
            vec_LC_cat_9_future_2050_SSP5_RCP85<-values(pr_sel_9)
                     
            #We save the vector of the category selected:
              save(vec_LC_cat_9_future_2050_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_9_future_2050_SSP5_RCP85.RData')
 

                           
            ################################################################################################################################################  
            #Future For our fragmentation variables
                rm(list=ls()) 
                ##Establecemos el fichero de trabajo
                setwd("G:/Project_Name/Spatial_Database/R_analysis")
                library(rgdal) 
                library(raster)
                library(spatialEco)
                library(readxl)
                #install.packages("spatialEco")
                
                #We select the original raster of land use that we want to transform:
                r1 <- raster("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
                #We are going to cut/crop the map for reduce to Europe
                #We stablish an extension for crop the map:
                e <- extent(-20,105,15,75)#Extension que utlizamos de normal 
                #We use the extension to crop the original map    
                m1_crop <- crop(r1,e) 
                
                #We extract the values of the raster to a vector:
                vector_val_select_natural<-values(m1_crop)
            
                
                use_categories<-read_excel("G:/Project_Name/Spatial_Database/R_analysis/Land_Use_classification/land_use_categories.xlsx", sheet = 1)
                df_use_categories<-as.data.frame(use_categories)

                df_use_categories_natural<-subset(df_use_categories,Natural_Human_used == "Natural")
                list_natural_cat<-unique(df_use_categories_natural$ID_CCI_LC_classification)
             
                df_use_categories_no_natural<-subset(df_use_categories,Natural_Human_used != "Natural")
                list_no_natural_cat<-unique(df_use_categories_no_natural$ID_CCI_LC_classification)

                   
                #We sustitute all values that are not our category by 0
                values <- list_no_natural_cat
                for(i in values) vector_val_select_natural[vector_val_select_natural== i] <- 0
                values_category <- list_natural_cat
                for(i in values_category) vector_val_select_natural[vector_val_select_natural== i] <- 1
                        
                        #We duplicate the raster with the exttension croped
                        m1_crop_sel_natural<-m1_crop
                        #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_7 
                        values(m1_crop_sel_natural)<-vector_val_select_natural
                        
                        
                        #Para distance to the border lo calcularemos sin agregar, es decir,
                        # exportamos m1_crop_sel_natural a idrisi y calculamos la distancia hacia adentro
                        
                        
                        
                        #Para el kernel lo calcularemos sin agregar
                        
                        
                        
                        
                        #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
                        m1_crop_sel_natural_agreg<-aggregate(m1_crop_sel_natural, fact=3, fun=mean, expand=TRUE)
                       
                        #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                        new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                        #We are going to use Europe Albers Equal Area Conic  
                        newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                        #We defined the projection of our raster:
                        projection(new_raster) <- newproj
                        
                        #We project our data to the raster created with the new projection
                        pr_sel_natural <- projectRaster(m1_crop_sel_natural_agreg, new_raster, method='bilinear')
                        X11()
                        plot(pr_sel_natural)
                        #We save the raster:
                        save(pr_sel_natural,file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_natural.Rdata")
                        writeRaster(pr_sel_natural, 'G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/pr_sel_natural.img',overwrite=TRUE) 
                        #We extract the values of the raster to a vector:
                        vec_LC_cat_natural_current<-values(pr_sel_natural)
                        #We save the vector of the category selected:
                        save(vec_LC_cat_natural_current, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_natural_current.RData')
                        gaus_pr_sel_natural<-raster.gaussian.smooth(pr_sel_natural, sigma = 2, n = 11, type = mean)
                        writeRaster(gaus_pr_sel_natural, 'G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/gaus_pr_sel_natural.img',overwrite=TRUE) 
                        
                        #We extract the values of the raster to a vector:
                        vec_gaus_pr_sel_natural_SSP5_RCP85<-values(gaus_pr_sel_natural)

                        #We save the vector of the category selected:
                        save(vec_gaus_pr_sel_natural_SSP5_RCP85, file='G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_gaus_pr_sel_natural_SSP5_RCP85.Rdata')

                        
                        
              
##############################################################################################################################################  
   #For our category cober (Is the percentage of bushses and forests)
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
    #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 

    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,130,140,150,151,152,153,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            map_cober <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')

        #We save the raster:
           # save(map_cober,file="H:/G/Project_Name/Spatial_Database/Future_2050/Projected/SSP5_RCP85/map_cober.Rdata")

        #We extract the values of the raster to a vector:
            vec_cober_future_2050_SSP5_RCP85<-values(map_cober)
                     
        #We save the vector of the category selected:
          save(vec_cober_future_2050_SSP5_RCP85, file='vec_cober_future_2050_SSP5_RCP85.RData')
              
    
############################################################################################################
   #For category of water
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Original/SSP5_RCP85/Globio4_landuse_10sec_2050_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180,130,140,150,151,152,153,190,200,201,202,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(210)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1

            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #We reclass human to integer values to calculate water areas
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg >= 0.5 & m1_crop_sel_5_agreg <= 1 ] <- 1
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg > 0 & m1_crop_sel_5_agreg < 0.5 ] <- 0
            
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_water_future_2050_SSP5_RCP85 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='ngb')


            #We save the raster:
            save(pr_water_future_2050_SSP5_RCP85,file="pr_water_future_2050_SSP5_RCP85.Rdata")
            #load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            #writeRaster(pr_water, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_water_future_2050_SSP5_RCP85<-values(pr_water_future_2050_SSP5_RCP85)
            #We save the vector of the category selected:
            save(vec_water_future_2050_SSP5_RCP85, file='vec_water_future_2050_SSP5_RCP85.RData')
          
 
##########################################################################################              
##########################################################################################              
#0.2 Preparing the scenarios
##########################################################################################  
##########################################################################################              
            
    #0.2.1 Current land use
        rm(list=ls()) 
        setwd("G:/Project_Name/Spatial_Database/R_analysis")
        library(rgdal) 
        library(raster) 
        library(biomod2) 
        library(caret)
        #We load the vectors representing each variable of land use
        load("vec_LC_cat_1_current.RData")
        load("vec_LC_cat_2_current.RData")
        load("vec_LC_cat_3_current.RData")
        load("vec_LC_cat_4_current.RData")
        load("vec_LC_cat_5_current.RData")
        load("vec_LC_cat_7_current.RData")
        load("vec_LC_cat_8_current.RData")
        load("vec_LC_cat_9_current.RData")
        load("vec_gaus_pr_sel_natural_current.RData")
        #We bind the vecctors of land use as a DF
        DF_LC_current<-cbind(vec_LC_cat_1_current,vec_LC_cat_2_current,vec_LC_cat_3_current,vec_LC_cat_4_current,vec_LC_cat_5_current,vec_LC_cat_7_current,vec_LC_cat_8_current,vec_LC_cat_9_current,vec_gaus_pr_sel_natural_current)
        #We save the current land use data as a DF
        save(DF_LC_current, file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/DF_LC_current.RData")
   
    #0.2.2 Current climate data 
        rm(list=ls()) 
        setwd("G:/Project_Name/Spatial_Database/R_analysis")
        #We load the climate for current
        load("G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Projected/DF_spatial_data_1.RData")
        
        #We reorder the columns
        DF_spatial_data_2 <- DF_spatial_data_1[,c("CHELSA_bio10_1.tif","CHELSA_bio10_2.tif","CHELSA_bio10_3.tif",
                            "CHELSA_bio10_4.tif","CHELSA_bio10_5.tif","CHELSA_bio10_6.tif","CHELSA_bio10_7.tif",
                            "CHELSA_bio10_8.tif","CHELSA_bio10_9.tif","CHELSA_bio10_10.tif","CHELSA_bio10_11.tif",
                            "CHELSA_bio10_12.tif","CHELSA_bio10_13.tif","CHELSA_bio10_14.tif","CHELSA_bio10_15.tif",
                            "CHELSA_bio10_16.tif","CHELSA_bio10_17.tif","CHELSA_bio10_18.tif","CHELSA_bio10_19.tif")]
        DF_bioclim<-as.data.frame(DF_spatial_data_2)
        #We save as a DF
        save(DF_bioclim, file="G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Projected/DF_bioclim.RData")

    #0.2.3 Future climate data
        #0.2.3.1 We select the colums for create three DF one for each RCP
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_spatial_data6_k1.RData")
            DF_spatial_data6_k1<-as.data.frame(DF_spatial_data6_k1)
            #DF_spatial_data_6k2<-as.data.frame(DF_spatial_data6_k2)
            #DF_spatial_data6_k1<-cbind(DF_spatial_data6_k1,DF_spatial_data6_k2)  
            
            #0.2.3.1.1 RCP26
                #We reorder the columns
                DF_clim_fut_2040_rcp26 <- DF_spatial_data6_k1[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_1.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_2.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_3.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_4.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_5.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_6.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_7.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_8.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_9.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_10.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_11.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_12.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_13.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_14.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_15.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_16.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_17.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_18.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp26_r1i1p1_g025.nc_19.tif")]
         
                 save(DF_clim_fut_2040_rcp26, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp26.RData") 
     
            #0.2.3.1.2 RCP45
                #We reorder the columns
                DF_clim_fut_2040_rcp45 <- DF_spatial_data6_k1[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_1.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_2.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_3.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_4.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_5.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_6.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_7.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_8.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_9.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_10.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_11.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_12.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_13.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_14.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_15.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_16.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_17.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_18.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp45_r1i1p1_g025.nc_19.tif")]
                 
                 save(DF_clim_fut_2040_rcp45, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp45.RData") 
                 
            #0.2.3.1.3 RCP60
                #0.2.3.1.3.1 RCP60 Part 1
                    #We select the columns with data of RCP60 and create a DF for Part 1 of RCP60
                     DF_clim_fut_2040_rcp60_1 <- DF_spatial_data6_k1[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_1.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_2.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_10.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_11.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_12.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_13.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_14.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_15.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_16.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_17.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_18.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_19.tif")]
                     
                     save(DF_clim_fut_2040_rcp60_1, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60_1.RData") 
                     
                #0.2.3.1.3.2 RCP60 Part 2
                    #We select the columns with data of RCP60    
                    rm(list=ls()) 
                    setwd("G:/Project_Name/Spatial_Database/R_analysis")
                    load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_spatial_data6_k2.RData")
                    load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60_1.RData")
                    DF_spatial_data6_k2<-as.data.frame(DF_spatial_data6_k2)
            
                    #We select the columns with data of RCP60 and create a DF for Part 2 of RCP60
                    DF_clim_fut_2040_rcp60_2 <- DF_spatial_data6_k2[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_3.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_4.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_5.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_6.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_7.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_8.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_9.tif")]
                    
                    #We bind the two DFs of RCP60 and ordered the colums     
                    DF_clim_fut_2040_rcp60<-cbind(DF_clim_fut_2040_rcp60_1,DF_clim_fut_2040_rcp60_2)
                    DF_clim_fut_2040_rcp60<-DF_clim_fut_2040_rcp60[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_1.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_2.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_3.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_4.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_5.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_6.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_7.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_8.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_9.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_10.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_11.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_12.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_13.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_14.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_15.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_16.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_17.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_18.tif",
                                                                  "CHELSA_bio_mon_IPSL-CM5A-MR_rcp60_r1i1p1_g025.nc_19.tif")]
                    
                    #We save the data of RCP60 in a unique DF
                    save(DF_clim_fut_2040_rcp60, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60.RData") 
                    
                    
            #0.2.3.1.4 RCP85
                #We select the columns with data of RCP85 and create a DF for RCP85    
                 DF_clim_fut_2040_rcp85 <- DF_spatial_data6_k2[,c("CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_1.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_2.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_3.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_4.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_5.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_6.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_7.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_8.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_9.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_10.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_11.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_12.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_13.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_14.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_15.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_16.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_17.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_18.tif",
                                                              "CHELSA_bio_mon_IPSL-CM5A-MR_rcp85_r1i1p1_g025.nc_19.tif")]
                 
                 #We save the data of RCP85 in a unique DF
                 save(DF_clim_fut_2040_rcp85, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp85.RData") 


        #0.2.3.2 We load the created DF for future climate and we change the colums names for that they have the same names as current data (mandatory for do the predicitons for future with the models)   
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
            #We change the variable names for that they be in coincidence with current clim data:     
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp26.RData") 
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp45.RData") 
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60.RData") 
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp85.RData") 
            load("G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Projected/DF_bioclim.RData") 
            G:/Project_Name/Spatial_Database/Variables/Climate/Current/Bioclimatic/Projected/DF_bioclim.RData
            colnames(DF_clim_fut_2040_rcp26)<-c("CHELSA_bio10_1.tif",
                                                  "CHELSA_bio10_2.tif",
                                                  "CHELSA_bio10_3.tif",
                                                  "CHELSA_bio10_4.tif",
                                                  "CHELSA_bio10_5.tif",
                                                  "CHELSA_bio10_6.tif",
                                                  "CHELSA_bio10_7.tif",
                                                  "CHELSA_bio10_8.tif",
                                                  "CHELSA_bio10_9.tif",
                                                  "CHELSA_bio10_10.tif",
                                                  "CHELSA_bio10_11.tif",
                                                  "CHELSA_bio10_12.tif",
                                                  "CHELSA_bio10_13.tif",
                                                  "CHELSA_bio10_14.tif",
                                                  "CHELSA_bio10_15.tif",
                                                  "CHELSA_bio10_16.tif",
                                                  "CHELSA_bio10_17.tif",
                                                  "CHELSA_bio10_18.tif",
                                                  "CHELSA_bio10_19.tif")
               
            colnames(DF_clim_fut_2040_rcp45)<-c("CHELSA_bio10_1.tif",
                                                  "CHELSA_bio10_2.tif",
                                                  "CHELSA_bio10_3.tif",
                                                  "CHELSA_bio10_4.tif",
                                                  "CHELSA_bio10_5.tif",
                                                  "CHELSA_bio10_6.tif",
                                                  "CHELSA_bio10_7.tif",
                                                  "CHELSA_bio10_8.tif",
                                                  "CHELSA_bio10_9.tif",
                                                  "CHELSA_bio10_10.tif",
                                                  "CHELSA_bio10_11.tif",
                                                  "CHELSA_bio10_12.tif",
                                                  "CHELSA_bio10_13.tif",
                                                  "CHELSA_bio10_14.tif",
                                                  "CHELSA_bio10_15.tif",
                                                  "CHELSA_bio10_16.tif",
                                                  "CHELSA_bio10_17.tif",
                                                  "CHELSA_bio10_18.tif",
                                                  "CHELSA_bio10_19.tif") 
                
            colnames(DF_clim_fut_2040_rcp60)<-c("CHELSA_bio10_1.tif",
                                                  "CHELSA_bio10_2.tif",
                                                  "CHELSA_bio10_3.tif",
                                                  "CHELSA_bio10_4.tif",
                                                  "CHELSA_bio10_5.tif",
                                                  "CHELSA_bio10_6.tif",
                                                  "CHELSA_bio10_7.tif",
                                                  "CHELSA_bio10_8.tif",
                                                  "CHELSA_bio10_9.tif",
                                                  "CHELSA_bio10_10.tif",
                                                  "CHELSA_bio10_11.tif",
                                                  "CHELSA_bio10_12.tif",
                                                  "CHELSA_bio10_13.tif",
                                                  "CHELSA_bio10_14.tif",
                                                  "CHELSA_bio10_15.tif",
                                                  "CHELSA_bio10_16.tif",
                                                  "CHELSA_bio10_17.tif",
                                                  "CHELSA_bio10_18.tif",
                                                  "CHELSA_bio10_19.tif")
                  
            colnames(DF_clim_fut_2040_rcp85)<-c("CHELSA_bio10_1.tif",
                                                  "CHELSA_bio10_2.tif",
                                                  "CHELSA_bio10_3.tif",
                                                  "CHELSA_bio10_4.tif",
                                                  "CHELSA_bio10_5.tif",
                                                  "CHELSA_bio10_6.tif",
                                                  "CHELSA_bio10_7.tif",
                                                  "CHELSA_bio10_8.tif",
                                                  "CHELSA_bio10_9.tif",
                                                  "CHELSA_bio10_10.tif",
                                                  "CHELSA_bio10_11.tif",
                                                  "CHELSA_bio10_12.tif",
                                                  "CHELSA_bio10_13.tif",
                                                  "CHELSA_bio10_14.tif",
                                                  "CHELSA_bio10_15.tif",
                                                  "CHELSA_bio10_16.tif",
                                                  "CHELSA_bio10_17.tif",
                                                  "CHELSA_bio10_18.tif",
                                                  "CHELSA_bio10_19.tif")
                    
            save(DF_clim_fut_2040_rcp26, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp26.RData") 
            save(DF_clim_fut_2040_rcp45, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp45.RData") 
            save(DF_clim_fut_2040_rcp60, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60.RData") 
            save(DF_clim_fut_2040_rcp85, file= "G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp85.RData") 
    
    #0.2.4 Future land use data 
        #0.2.4.1 Future land use data SSP1_RCP26 
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
            #We load and merge in one DF the vectors representing land cover variables
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_1_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_2_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_3_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_4_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_5_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_7_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_8_future_2050_SSP1_RCP26.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/vec_LC_cat_9_future_2050_SSP1_RCP26.RData")
            DF_LC_future_2050_SSP1_RCP26<-cbind(vec_LC_cat_1_future_2050_SSP1_RCP26,vec_LC_cat_2_future_2050_SSP1_RCP26,vec_LC_cat_3_future_2050_SSP1_RCP26,vec_LC_cat_4_future_2050_SSP1_RCP26,vec_LC_cat_5_future_2050_SSP1_RCP26,vec_LC_cat_7_future_2050_SSP1_RCP26,vec_LC_cat_8_future_2050_SSP1_RCP26,vec_LC_cat_9_future_2050_SSP1_RCP26)
            save(DF_LC_future_2050_SSP1_RCP26, file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/DF_LC_future_2050_SSP1_RCP26.RData")
               
        #0.2.4.2 Future land use data SSP3_RCP70
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
            #We load and merge in one DF the vectors representing land cover variables
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_1_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_2_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_3_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_4_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_5_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_7_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_8_future_2050_SSP3_RCP70.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/vec_LC_cat_9_future_2050_SSP3_RCP70.RData")
            DF_LC_future_2050_SSP3_RCP70<-cbind(vec_LC_cat_1_future_2050_SSP3_RCP70,vec_LC_cat_2_future_2050_SSP3_RCP70,vec_LC_cat_3_future_2050_SSP3_RCP70,vec_LC_cat_4_future_2050_SSP3_RCP70,vec_LC_cat_5_future_2050_SSP3_RCP70,vec_LC_cat_7_future_2050_SSP3_RCP70,vec_LC_cat_8_future_2050_SSP3_RCP70,vec_LC_cat_9_future_2050_SSP3_RCP70)
            save(DF_LC_future_2050_SSP3_RCP70, file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/DF_LC_future_2050_SSP3_RCP70.RData")
               
        #0.2.4.3 Future land use data SSP5_RCP85
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
           #We load and merge in one DF the vectors representing land cover variables
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_1_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_2_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_3_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_4_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_5_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_7_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_8_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_LC_cat_9_future_2050_SSP5_RCP85.RData")
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/vec_gaus_pr_sel_natural_SSP5_RCP85.Rdata")
            DF_LC_future_2050_SSP5_RCP85<-cbind(vec_LC_cat_1_future_2050_SSP5_RCP85,vec_LC_cat_2_future_2050_SSP5_RCP85,vec_LC_cat_3_future_2050_SSP5_RCP85,vec_LC_cat_4_future_2050_SSP5_RCP85,vec_LC_cat_5_future_2050_SSP5_RCP85,vec_LC_cat_7_future_2050_SSP5_RCP85,vec_LC_cat_8_future_2050_SSP5_RCP85,vec_LC_cat_9_future_2050_SSP5_RCP85,vec_gaus_pr_sel_natural_SSP5_RCP85)
            save(DF_LC_future_2050_SSP5_RCP85, file="G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/DF_LC_future_2050_SSP5_RCP85.RData")
        #We save the current land use data as a DF
        save(DF_LC_current, file="G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/DF_LC_current.RData")

####################################################################################################
####################################################################################################
#0.3 We prepare the data by scenarios for BIOMOD analysis 
####################################################################################################
####################################################################################################
        
    #We have current data/scenario (called Scenario_current) and three different scenarios for future (RCP26, RCP60, RCP85)        

        #0.3.1 Scenario Current data
            rm(list=ls()) 
            setwd("G:/Project_Name/Spatial_Database/R_analysis")
            #Climate current. A data frame with columns representing each of the 19 bioclimatic variables and rows representing pixels:
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/DF_bioclim.RData")
            colnames(DF_bioclim)<-c("Clim_1","Clim_2","Clim_3","Clim_4","Clim_5","Clim_6","Clim_7","Clim_8","Clim_9","Clim_10","Clim_11","Clim_12","Clim_13","Clim_14","Clim_15","Clim_16","Clim_17","Clim_18","Clim_19")
            #Land use current. A data frame with columns representing each of the 8 land cover variables and rows representing pixels:
            load("G:/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/DF_LC_current.RData")
            colnames(DF_LC_current)<-c("LC_1","LC_2","LC_3","LC_4","LC_5","LC_7","LC_8","LC_9","Frag")
            Scenario_current<-as.data.frame(cbind(DF_bioclim,DF_LC_current))
            #We are going to specify the NA values (sea) using the DF_bioclim
            Scenario_current$LC_1[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_2[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_3[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_4[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_5[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_7[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_8[is.na(Scenario_current$Clim_1)]<-NA
            Scenario_current$LC_9[is.na(Scenario_current$Clim_1)]<-NA   
            Scenario_current$Frag[is.na(Scenario_current$Clim_1)]<-NA                             
            save(Scenario_current, file="G:/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_current.RData")
            
        #0.3.2 Future data
            #Land use future DFs. Three data frames. Each DF has columns representing each of the 8 land cover variables and rows representing pixels:
            #Climate future DFs. Three data frames. Each DF has columns representing each of the 19 bioclimatic variables and rows representing pixels:
            #0.3.2.1 Scenario RCP26
                rm(list=ls()) 
                setwd("G:/Project_Name/Spatial_Database/R_analysis")
                #Climate future RCP26
                load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp26.RData")
                colnames(DF_clim_fut_2040_rcp26)<-c("Clim_1","Clim_2","Clim_3","Clim_4","Clim_5","Clim_6","Clim_7","Clim_8","Clim_9","Clim_10","Clim_11","Clim_12","Clim_13","Clim_14","Clim_15","Clim_16","Clim_17","Clim_18","Clim_19")
                #Land use future RCP26
                load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP1_RCP26/DF_LC_future_2050_SSP1_RCP26.RData")
                colnames(DF_LC_future_2050_SSP1_RCP26)<-c("LC_1","LC_2","LC_3","LC_4","LC_5","LC_7","LC_8","LC_9")
                Scenario_RCP26<-as.data.frame(cbind(DF_clim_fut_2040_rcp26,DF_LC_future_2050_SSP1_RCP26))
                #We are going to specify the NA values (sea) using the DF_bioclim
                Scenario_RCP26$LC_1[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_2[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_3[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_4[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_5[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_7[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_8[is.na(Scenario_RCP26$Clim_1)]<-NA
                Scenario_RCP26$LC_9[is.na(Scenario_RCP26$Clim_1)]<-NA                
                save(Scenario_RCP26, file="G:/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP26.RData")

            #0.3.2.2 Scenario RCP60
                rm(list=ls()) 
                setwd("G:/Project_Name/Spatial_Database/R_analysis")                
                #Climate future RCP60
                load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp60.RData")
                colnames(DF_clim_fut_2040_rcp60)<-c("Clim_1","Clim_2","Clim_3","Clim_4","Clim_5","Clim_6","Clim_7","Clim_8","Clim_9","Clim_10","Clim_11","Clim_12","Clim_13","Clim_14","Clim_15","Clim_16","Clim_17","Clim_18","Clim_19")
                #Land use future RCP60 (There was a kind of error in the name of the original files of future scenarios of land use and instead of call it RCP70 they call RCP60, this is based in the table that Wilfred send me the 7 March 2018, 14:21)
                load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP3_RCP70/DF_LC_future_2050_SSP3_RCP70.RData")
                colnames(DF_LC_future_2050_SSP3_RCP70)<-c("LC_1","LC_2","LC_3","LC_4","LC_5","LC_7","LC_8","LC_9")
                Scenario_RCP60<-as.data.frame(cbind(DF_clim_fut_2040_rcp60,DF_LC_future_2050_SSP3_RCP70))
                #We are going to specify the NA values (sea) using the DF_bioclim
                Scenario_RCP60$LC_1[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_2[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_3[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_4[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_5[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_7[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_8[is.na(Scenario_RCP60$Clim_1)]<-NA
                Scenario_RCP60$LC_9[is.na(Scenario_RCP60$Clim_1)]<-NA
                save(Scenario_RCP60, file="G:/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP60.RData")
                
            #0.3.2.3 Scenario RCP85
                rm(list=ls()) 
                setwd("G:/Project_Name/Spatial_Database/R_analysis")                
                #Climate future RCP85
                load("G:/Project_Name/Spatial_Database/Variables/Climate/Future_2040/Bioclimatic/Projected/DF_clim_fut_2040_rcp85.RData")
                colnames(DF_clim_fut_2040_rcp85)<-c("Clim_1","Clim_2","Clim_3","Clim_4","Clim_5","Clim_6","Clim_7","Clim_8","Clim_9","Clim_10","Clim_11","Clim_12","Clim_13","Clim_14","Clim_15","Clim_16","Clim_17","Clim_18","Clim_19")
                #Land use future RCP85
                load("G:/Project_Name/Spatial_Database/Variables/Land_use/Future_2050/Projected/SSP5_RCP85/DF_LC_future_2050_SSP5_RCP85.RData")
                colnames(DF_LC_future_2050_SSP5_RCP85)<-c("LC_1","LC_2","LC_3","LC_4","LC_5","LC_7","LC_8","LC_9","Frag")
                Scenario_RCP85<-as.data.frame(cbind(DF_clim_fut_2040_rcp85,DF_LC_future_2050_SSP5_RCP85))
                #We are going to specify the NA values (sea) using the DF_bioclim
                Scenario_RCP85$LC_1[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_2[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_3[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_4[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_5[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_7[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_8[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$LC_9[is.na(Scenario_RCP85$Clim_1)]<-NA
                Scenario_RCP85$Frag[is.na(Scenario_RCP85$Clim_1)]<-NA                             
                save(Scenario_RCP85, file="G:/Project_Name/Spatial_Database/Variables/Final_scenarios/Scenario_RCP85.RData")

        #0.3.3 Obtaintion of coordinates of our raster cells that we are going to use to calculate the buffer around presences and obtain pseudoabsences
          #We create a new raster (with the same resolution of spatial extent as variables over Wurope)
          #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
          new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
          #We are going to use Europe Albers Equal Area Conic  
          newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
          #We defined the projection of our raster:
          projection(new_raster) <- newproj
          vec_values<-seq(1,33040000,1)
          values(new_raster)<-vec_values
          #We obtain a matrix with three columns (x, y and an ID for the pixel)
          data_matrix <- rasterToPoints(new_raster)
          colnames(data_matrix)<-c("x","y","ID_pixel")
          save(data_matrix, file="G:/Project_Name/Spatial_Database/Variables/data_matrix.RData")


          

####################################################################################################
####################################################################################################
#0.4 Creation of new variables 
####################################################################################################
####################################################################################################
          
    vec_LC_2_dis,
    cober,
    water,
           

#0.4.1 Creation of new variable distance to agriculture  
  rm(list=ls())
  library(rgdal) 
  library(raster) 
  library(readxl)  
  setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
  
  #We load the scenario current abiotic
  load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/Scenario_current.RData")

 #For agriculture 
   var_LC_2<-(Scenario_current$LC_2)
   mean(var_LC_2,na.rm = TRUE)
   median(var_LC_2,na.rm = TRUE)
   #We reclass human to integer values to calculate the variable distance to human 
   var_LC_2[var_LC_2 > 0 & var_LC_2 <= 1 ] <- 1
   mean(var_LC_2,na.rm = TRUE)
   median(var_LC_2,na.rm = TRUE)
   table(var_LC_2)
   #We are going to use Europe Albers Equal Area Conic  
   newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
   #Raster of reference: 
   map_LC_2 <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
   #We defined the projection of our raster:
   projection(map_LC_2) <- newproj
   values(map_LC_2)<-var_LC_2
   #save(Europe_presences_0_1_raster,file="Europe_presences_0_1_raster.RData")
   #We exported the raster
   writeRaster(map_LC_2, "map_LC_2.rst", datatype='INT4S', overwrite=TRUE) #
   #WE have calculated the distance variable in idrisi  map_LC_2_dis.rst
   raster_map_LC_2_dis<-raster("map_LC_2_dis.rst")
   vec_LC_2_dis<-values(raster_map_LC_2_dis)
   save(vec_LC_2_dis, file="vec_LC_2_dis.Rdata")  

 #For distance to agriculture and urban
   var_LC_1_and_2<-var_LC_2+var_LC_1
   table(var_LC_1_and_2)
   var_LC_1_and_2[var_LC_1_and_2 > 0 & var_LC_1_and_2 <= 2 ] <- 1
   table(var_LC_1_and_2)
   var_LC_1_and_2[is.na(var_LC_1_and_2)] <- 0
   
   #We are going to use Europe Albers Equal Area Conic  
   newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
   #Raster of reference: 
   map_LC_1_and_2 <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
   #We defined the projection of our raster:
   projection(map_LC_1_and_2) <- newproj
   values(map_LC_1_and_2)<-var_LC_1_and_2
              NAvalue(map_LC_1_and_2)
     map_LC_1_and_2[map_LC_1_and_2 == -Inf ] <- 0
   #save(Europe_presences_0_1_raster,file="Europe_presences_0_1_raster.RData")
   #We exported the raster
   writeRaster(map_LC_1_and_2, "map_LC_1_and_2.rst", datatype='INT4S', overwrite=TRUE) #
   #WE have calculated the distance variable in idrisi  map_LC_1_and_2_dis.rst
   raster_map_LC_1_and_2_dis<-raster("map_LC_1_and_2_dis.rst")
   vec_LC_1_and_2_dis<-values(raster_map_LC_1_and_2_dis)
   save(vec_LC_1_and_2_dis, file="vec_LC_1_and_2_dis.Rdata")  
   
   

   
############################################################################################################
#5C Creation of a new variable of bushes
############################################################################################################
   #For our category 10
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("H:/G/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,130,140,150,151,152,153,190,200,201,202,210,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(100,110,120,121,122,160,170,180)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1
     
            
            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_sel_10 <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='bilinear')


            #We save the raster:
            save(pr_sel_10,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_10.Rdata")
            load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_10.Rdata")
            writeRaster(pr_sel_10, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_sel_10.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_LC_cat_10_current<-values(pr_sel_10)
            #We save the vector of the category selected:
            save(vec_LC_cat_10_current, file='vec_LC_cat_10_current.RData')
            str(vec_LC_cat_10_current)
            
            
 #For kernels of % of cober areas
  setwd("C:/Users/First_Author_Name/Documents/Bear/Analisis_tochosR") 
  
  #We load the scenario current abiotic
  load("H:/D/Project_Name/Results_Biomod/Data_for_analysis/Scenario_current.RData")
            
   var_LC_3<-(Scenario_current$LC_3)
   var_LC_4<-(Scenario_current$LC_4)
   var_cober<-var_LC_3+var_LC_4+vec_LC_cat_10_current
   #We are going to use Europe Albers Equal Area Conic  
   save(var_cober, file='var_cober.RData')
 
   #We are going to use Europe Albers Equal Area Conic  
   newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
   #Raster of reference: 
   map_cober <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
   #We defined the projection of our raster:
   projection(map_cober) <- newproj
   values(map_cober)<-var_cober
              NAvalue(map_cober)
     map_cober[map_cober == -Inf ] <- 0
   #save(Europe_presences_0_1_raster,file="Europe_presences_0_1_raster.RData")
   #We exported the raster
   writeRaster(map_cober, "map_cober.rst", overwrite=TRUE) #
   
   #WE have calculated the kernels variable in idrisi 
   raster_map_cober_kernel_gauss_3<-raster("map_cober_kernel_gauss_3x3.rst")
   raster_map_cober_kernel_gauss_5<-raster("map_cober_kernel_gauss_5x5.rst")
   raster_map_cober_kernel_gauss_7<-raster("map_cober_kernel_gauss_7x7.rst")
   
   
   vec_cober_kernel_3<-values(raster_map_cober_kernel_gauss_3)
   vec_cober_kernel_5<-values(raster_map_cober_kernel_gauss_5)
   vec_cober_kernel_7<-values(raster_map_cober_kernel_gauss_7)
   
   save(vec_cober_kernel_3, file="vec_cober_kernel_3.Rdata")  
   save(vec_cober_kernel_5, file="vec_cober_kernel_5.Rdata")  
   save(vec_cober_kernel_7, file="vec_cober_kernel_7.Rdata")  


   
############################################################################################################
#5d Creation of a new variable of water
############################################################################################################
   #For our category 10
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("H:/G/Project_Name/Spatial_Database/R_analysis")
    library(rgdal) 
    library(raster) 
     #We select the original raster of land use that we want to transform:
        r1 <- raster("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Original/Globio4_landuse_10sec_2015_cropint.tif") 
   
    #We are going to cut/crop the map for reduce to Europe
        #We stablish an extension for crop the map:
        e <- extent(-20,105,15,75)#Extension que utlizamos de normal 

        #We use the extension to crop the original map    
            m1_crop <- crop(r1,e) 
    
            #We extract the values of the raster to a vector:
                vector_val_select_5<-values(m1_crop)

              #We sustitute all values that are not our category by 0
              values <- c(0,1,2,230,231,232,3,4,5,6,7,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,160,170,180,130,140,150,151,152,153,190,200,201,202,220)
              for(i in values) vector_val_select_5[vector_val_select_5== i] <- 0
              values_category <- c(210)
              for(i in values_category) vector_val_select_5[vector_val_select_5== i] <- 1

            #We duplicate the raster with the exttension croped
            m1_crop_sel_5<-m1_crop
            #We sustitute the data of the duplicated raster by the vector of data of our category with the sustitution values vector_val_select_5 
            values(m1_crop_sel_5)<-vector_val_select_5
            #Now we are going to aggregate the cells of the map. The term fact is the aggregation factor expressed as number of cells in each direction (horizontally and vertically)
            m1_crop_sel_5_agreg<-aggregate(m1_crop_sel_5, fact=3, fun=mean, expand=TRUE)
           
            #We reclass human to integer values to calculate water areas
            m1_crop_sel_5_agreg[m1_crop_sel_5_agreg >= 0.5 & m1_crop_sel_5_agreg <= 1 ] <- 1

            
            #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
            #We are going to use Europe Albers Equal Area Conic  
            newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            
            #We project our data to the raster created with the new projection
            pr_water <- projectRaster(m1_crop_sel_5_agreg, new_raster, method='ngb')


            #We save the raster:
            save(pr_water,file="H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            load("H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.Rdata")
            writeRaster(pr_water, "H:/G/Project_Name/Spatial_Database/Variables/Land_use/Current/Projected/pr_water.rst", overwrite=TRUE) #
            #We extract the values of the raster to a vector:
            vec_water<-values(pr_water)
            #We save the vector of the category selected:
            save(vec_water, file='vec_water.RData')
           
        
          
          