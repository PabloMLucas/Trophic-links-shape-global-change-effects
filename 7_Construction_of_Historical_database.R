

############################################################################################################
#Readme:
############################################################################################################
#R code to calculate the historical database. Using historical climate and historical range of brown bear
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
  #Files downloaded from CRU TS v. 4.02 download from https://crudata.uea.ac.uk/cru/data/hrg/ 
  #shapefile('Merged_digitized_hist_range_IUCN2.shp') # The constructed historical distribution of brown bear
#Data output:
  #presences_absences_eurasia_hist.RData

##############################################################################################                
#Schema
############################################################################################## 
#7_Construction_of_Historical_database
  #7.1 Download of files from CHELSA_TraCE21k V1
  #7.2 Project and aggregate the files
  #7.3 Rasterize historical distributions bear
  #7.4 Calculate historical climate for each historical distribution

########################################################################################################################################################
#7.1 Download of files from CHELSA_TraCE21k V1
########################################################################################################################################################

  rm(list=ls()) 
  setwd("D:/Project_name/Spatial_Database/R_analysis")
  #install.packages("RCurl")
  #install.packages("httr")
  #install.packages("XML")
  install.packages("sf")
  library(httr)
  library(RCurl)
  library(XML)
  library(raster)
  library(sf)
  
  #File example: "CHELSA_TraCE21k_bio1_-10_V1.0.tif"
  #vectortimes<-c(-10:20)
  vectortimes<-c(-1:20)#This is to download from yearID -3 becasue we have already downloaded the others
  length(vectortimes)
  
    for (t in 1:length(vectortimes)){
      timeID<-vectortimes[t]
      print(paste0("TimeID: ",timeID))
      for (b in 1:19){
        if (b<10){
          b<-paste0("0",b)
        }
        print(paste0("Bioclimatic_variable ",b)) 
          shell.exec(paste("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_trace/bio/CHELSA_TraCE21k_bio",b,"_",timeID,"_V1.0.tif", sep=""))
        
        #We will use this for that the code stop during 30 seconds avoiding problems with chrome when downloading data. Also notice that 
        # some servers can detect this code as a hacking atack, so it is better to stop the code with this to avoid it. Estimate how much time takes
          #your connection to download data and fit to it
        Sys.sleep(30)    
      }
    }
  
  
########################################################################################################################################################
#7.2 Project and aggregate the files
########################################################################################################################################################
      
  rm(list=ls()) 
  folder<-"D:/Project_name/Spatial_Database/R_analysis/"
  setwd(paste0(folder,"Historical/"))
  folder_data<-paste0(folder,"CHELSA-TraCE21k/")
  
  #We define the raster of reference at 50x50 km in CEA:
  raster_cea <- raster(xmn=-20000000, xmx=20000000, ymn=-6300000, ymx=6300000,  ncols=800, nrows=252)#Raster que utlizamos de normal 

  cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" 
  projection(raster_cea) <- cea_proj
  vectortimes<-c(-10:20)  

  for (y in 1:length(vectortimes)){ 
    timeID<-vectortimes[y]
    print(paste0("TimeID: ",timeID))
    for (b in 1:19){
      if (b<10){
        b<-paste0("0",b)
      }
      print(paste0("Bioclimatic_variable ",b)) 
      #We load the original download raster 
      var<-raster(paste0(folder_data,"CHELSA_TraCE21k_bio",b,"_",timeID,"_V1.0.tif", sep=""))
      #We aggregate the values to 66x66 cells to have moreless the same resolution that the CEA raster at 50x50km
      var_agreg<-aggregate(var, fact=50, fun=mean, expand=TRUE)
      #We defined the projection of our raster:
      raster_bio_cea <- projectRaster(var_agreg, raster_cea, method='bilinear')      
      vec_raster_bio_cea<-values(raster_bio_cea)
      assign(paste0("TraCE21k_bio_cea50km_",b,"_",timeID, sep=""),vec_raster_bio_cea)  
      save(list=paste0("TraCE21k_bio_cea50km_",b,"_",timeID, sep=""),file=paste0(folder_data,"TraCE21k_bio_cea50km_",b,"_",timeID,".RData", sep=""))       
      #file.remove(paste0(folder_data,"CHELSA_TraCE21k_bio",b,"_",timeID,"_V1.0.tif", sep=""))
      }
  }



########################################################################################################################################################
#7.3 Rasterize historical distributions bear
########################################################################################################################################################

  rm(list=ls()) 
  folder<-"D:/Project_name/Spatial_Database/R_analysis/"
  setwd(paste0(folder,"Historical/"))
  folder_data<-paste0(folder,"CHELSA-TraCE21k/")
  
  library(httr)
  library(RCurl)
  library(XML)
  library(raster)
  library(sf)
  
  #Raster of reference: 
    cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs" 
    #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data
    raster_cea <- raster(xmn=-20000000, xmx=20000000, ymn=-6300000, ymx=6300000,  ncols=800, nrows=252)
    projection(raster_cea) <- cea_proj

    shape_continents <- shapefile("F:/G/Project_name/Maps/continent shapefile/continent.shp")#Version 02/06/2
    shape_continents_cea<-  spTransform(shape_continents, cea_proj)
    #Rasterize the map
    raster_continents_cea<-rasterize(shape_continents_cea, raster_cea, field=1)

 
    #We obtain the coordinates of the cells:
    DF_xy<- as.data.frame(xyFromCell(raster_cea,1:201600))
    
  #File example: "CHELSA_TraCE21k_bio1_-10_V1.0.tif"
  vectortimes<-c(-1:20)#This is to download from yearID -3 becasue we have already downloaded the others
  length(vectortimes)

  #cea_proj <- "+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs" 
  folder_maps<-"D:/Project_name/Spatial_Database/R_analysis/Historical ranges/Ranges/"

  maps_bear<-c("Hist_range_Spain_varios_LAEA","Hist_range_Sp_Gal_Nav_LAEA","Hist_range_Spain_1800_LAEA",
              "Portugal_Bencatel_et_al_range_3035_v1.proj","Portugal_Alvares_and_Domingues_rangeConfirmado_3035_v1",
              "France_XV_XVI_range_3035_v1",
              "Italy_historical_range_3035_v1",
              "Yugoslavia_1800_range_3035",
              "Greece_II_century_range_3035_v2",
              "species_41688_dis2")
 
  area_bear<-c("Area_Spain_varios","Area_Gal_Nav","Area_Spain_1800",#1-3
              "Area_Portugal_Bencatel","Area_Portugal_Alvares",#4-5
              "Area_France",#6
              "Area_Italy",#7
              "Area_Yugoslavia",#8
              "Area_Greece",#9
              "Continent_sub")#10

  centuries_bear<-c("14th Century","18th Century","19th Century",#1-3
              "17th Century","16th-18th Centuries",#4-5
              "15-16th Centuries",#6
              "16-18th Centuries",#7
              "19th Century",#8
              "2nd Century",#9
              "16-20th Centuries")#10

  timeID_bear_list<-list(c(14),c(18),c(19),
            c(17),c(16,17,18),
            c(15,16),
            c(16,18),
            c(19),
            c(2),
            c(16,17,18,19,20))

env_var_bio_var_in_study_area_all_ranges<-c()  
for (m in 1:length(maps_bear)){#Loop for each distribution map
  print("############################################################")
  print(paste0("Map bear: ",maps_bear[m]))
  print(paste0("Area bear: ",area_bear[m]))
  print(paste0("TimeID: ",centuries_bear[m]))
   
  if (maps_bear[m]=="species_41688_dis2"){
  #Write the name of the files 
  historic_range_bear<- paste0(folder_maps,maps_bear[m],".shp") 
  historic_area<- paste0(folder_maps,area_bear[m],".shp") 

  #Convert the files to a SpatialPolygonsDataFrame to work in R
  shape_bear <- shapefile(historic_range_bear)
  shape_area <- shapefile(historic_area)   

  #We project the shapefile into CEA
  shape_bear_cea<-  spTransform(shape_bear, cea_proj)
  shape_area_cea<-  spTransform(shape_area, cea_proj)
  
  #Rasterize the maps
    raster_bear_cea<-rasterize(shape_bear_cea, raster_cea, field=1)
    raster_area_cea<-rasterize(shape_area_cea, raster_cea, field=1)
 
 }else{
  #Write the name of the files 
  historic_range_bear<- paste0(folder_maps,maps_bear[m],".shp") 
  historic_area<- paste0(folder_maps,area_bear[m],".shp") 
  #Convert the files to a SpatialPolygonsDataFrame to work in R
  shape_bear <- shapefile(historic_range_bear)
  shape_area <- shapefile(historic_area)#
  
  #We project the shapefile into CEA
  shape_bear_cea<-  spTransform(shape_bear, cea_proj)
  shape_area_cea<-  spTransform(shape_area, cea_proj)
  #Rasterize the maps
    raster_bear_cea0<-rasterize(shape_bear_cea, raster_cea, field=1)
      raster_bear_cea0[is.na(raster_bear_cea0)[]] <- 0 
    #rasterize_lines
      # convert SPDF to SLDF
      bound <- as(shape_bear_cea, 'SpatialLinesDataFrame')  
      raster_bear_cea_lines<-rasterize(bound, raster_cea, field=1)
      raster_bear_cea_lines[is.na(raster_bear_cea_lines)[]] <- 0 
      
      raster_bear_cea1 <- overlay(raster_bear_cea0, raster_bear_cea_lines, fun=sum)
      raster_bear_cea2<-reclassify(raster_bear_cea1,c(0.9, 2.1, 1))
      raster_bear_cea<-raster_bear_cea2
      raster_bear_cea[(raster_bear_cea==0[])] <- NA 
      
    raster_area_cea<-rasterize(shape_area_cea, raster_cea, field=1)  
  }
  #Apply a mask of continents  
  raster_bear_cea <- mask(raster_bear_cea, raster_continents_cea)
  raster_area_cea <- mask(raster_area_cea, raster_continents_cea)
  
  #Calculate the buffers
      #We calculate a raster with the buffer minimum   
      rasterOptions(todisk=TRUE)
      buf_max_bear_cea<-buffer(raster_bear_cea, width=3000000)
      buf_min_bear_cea<-buffer(raster_bear_cea, width=10) 
      rasterOptions(todisk=FALSE)  
      buf_max_bear_cea <- mask(buf_max_bear_cea, raster_continents_cea)
      buf_min_bear_cea <- mask(buf_min_bear_cea, raster_continents_cea)

  #We visuzlize the rater to check the results
  #Extract the coordinates from the vector to do a map smaller to better visualize
  x_min<-shape_area_cea@ bbox["x","min"]
  x_max<-shape_area_cea@ bbox["x","max"]
  y_min<-shape_area_cea@ bbox["y","min"]
  y_max<-shape_area_cea@ bbox["y","max"]

  if (maps_bear[m]=="species_41688_dis2"){
  x_min<--4000000
  y_min<-0
  } 
    
#To check the raster of the bear and area  
 x_min_for_clip_raster<-round((x_min-100000)/10000)*10000
 y_min_for_clip_raster<-round((y_min-100000)/10000)*10000
 x_max_for_clip_raster<-round((x_max+100000)/10000)*10000
 y_max_for_clip_raster<-round((y_max+100000)/10000)*10000
 
  e <- extent(x_min_for_clip_raster,x_max_for_clip_raster,y_min_for_clip_raster,y_max_for_clip_raster)
  extract_raster_bear_cea<-crop(raster_bear_cea, e)
  #extract_raster_bear_cea_lines<-crop(raster_bear_cea_lines, e)
  
  extract_raster_area_cea<-crop(raster_area_cea, e)
  shape_continents_cea1<-crop(shape_continents_cea,e)

#To check the raster of the buffers 
 x_min_for_clip_raster<-round((x_min-600000)/10000)*10000
 y_min_for_clip_raster<-round((y_min-600000)/10000)*10000
 x_max_for_clip_raster<-round((x_max+600000)/10000)*10000
 y_max_for_clip_raster<-round((y_max+600000)/10000)*10000
 
  e2 <- extent(x_min_for_clip_raster,x_max_for_clip_raster,y_min_for_clip_raster,y_max_for_clip_raster)
  crop_buf_max_bear_cea<-crop(buf_max_bear_cea, e2)
  crop_buf_min_bear_cea<-crop(buf_min_bear_cea, e2)
  shape_continents_cea2<-crop(shape_continents_cea,e2)
 
  #We create plots to check rasterization and area for pseudoabsences
  pdf(file=paste0("Plot_",maps_bear[m],".pdf"), width=180*0.0393701, height=180*0.0393701, pointsize = 5)
    par(mfrow = c(2, 2))
    par(mar=c(1.9,1.9,3,2))
    #Plot for rasterization of brown bear distribution
    plot(shape_continents_cea1, lwd=2, col='grey90', border=NA)
    plot(extract_raster_bear_cea, add=TRUE, col='grey60', border=NA) 
    plot(shape_bear_cea, add=TRUE, density=10, lwd=0.5, col=NA, border=NA)
    title(main = list("Historical range", cex=0.8))
    #Plot for rasterization of study area
    plot(shape_continents_cea1, lwd=0.01, col='grey90', border=NA)
    plot(extract_raster_area_cea, add=TRUE, col='grey75', border=NA) 
    plot(extract_raster_bear_cea, add=TRUE, col='grey60', border=NA) 
    plot(shape_area_cea, add=TRUE, density=5, lwd=0.5, col='black', border=NA)
    title(main = list("Study area", cex=0.8))
    #Plot for buffers for pseudoabsences
    plot(shape_continents_cea2, lwd=0.2, col='grey90', border=NA)
    plot(crop_buf_max_bear_cea, add=TRUE, col='yellow')
    plot(crop_buf_min_bear_cea, add=TRUE, col='red')  
    plot(shape_bear_cea, add=TRUE, density=10, lwd=2, col='black', border=NA)
    plot(shape_area_cea, add=TRUE, density=10, lwd=2, col='black', border=NA)
    title(main = list("Background for pseudoabsences", cex=0.8))
   mtext(paste0("Map bear: ",maps_bear[m]), side = 3, line = - 1, outer = TRUE)
   dev.off()     

#Extract the values  
         val_bear_cea<-   values(raster_bear_cea)
         val_area_cea<-   values(raster_area_cea)
         val_buf_max_bear_cea<-   values(buf_max_bear_cea)
         val_buf_min_bear_cea<-   values(buf_min_bear_cea)
         val_continent_cea<-   values(raster_continents_cea)  
  values_env_var<-data.frame("Range"=val_bear_cea, "Study_area"=val_area_cea, "Buf_max"=val_buf_max_bear_cea, "Buf_min"=val_buf_min_bear_cea, "Land"=val_continent_cea,"Coord_x"=DF_xy$x,"Coord_y"=DF_xy$y)

  
########################################################################################################################################################
#7.4 Calculate historical climate for each historical distribution
########################################################################################################################################################

  #Add the specific climate based on TimeID
  timeID_bear_list_in_loop<-timeID_bear_list[m]
  bio_var_all<-c()
  n_centuries<- length(timeID_bear_list_in_loop[[1]])
  print(paste0("Number of centuries ",n_centuries)) 
  
  if(length(timeID_bear_list_in_loop[[1]])==1){ #if the distribution is from one century, then we will just took the climate of that century:
    timeID<- timeID_bear_list_in_loop[[1]][1]
    print(paste0("TimeID ",timeID)) 
    for (b in 1:19){#loop for each bioclimatic variable
      if (b<10){#for bio less than 10 we will add a 0
        b<-paste0("0",b)
      }#for bio less than 10 we will add a 0
      print(paste0("Bioclimatic_variable ",b)) 
      #load the data
      load(paste0(folder_data,"TraCE21k_bio_cea50km_",b,"_",timeID,".RData", sep=""))
      bio_var<-get(paste0("TraCE21k_bio_cea50km_",b,"_",timeID, sep="")) 
      bio_var_all<-cbind(bio_var_all,bio_var)
    }#loop for each bioclimatic variable  
  }else{#if the distribution is from several centuries, then we will calculate the average climate of those centuries:
  
    for (b in 1:19){#loop for each bioclimatic variable
      if (b<10){#for bio less than 10 we will add a 0
        b<-paste0("0",b)
      }#for bio less than 10 we will add a 0
      print(paste0("Bioclimatic_variable ",b)) 
      bio_var_all_centuries<-c()
      for (c in 1:n_centuries){#Loop for each century
        timeID<- timeID_bear_list_in_loop[[1]][c] 
        print(paste0("TimeID ",timeID)) 
        #load the data
        load(paste0(folder_data,"TraCE21k_bio_cea50km_",b,"_",timeID,".RData", sep=""))
        bio_var<-get(paste0("TraCE21k_bio_cea50km_",b,"_",timeID, sep="")) 
        bio_var_all_centuries<-cbind(bio_var_all_centuries,bio_var)
      }#Loop for each century  
      bio_var_all_centuries<-data.frame(bio_var_all_centuries)
      bio_var_all_centuries_mean<-rowMeans(bio_var_all_centuries,na.rm=T)
      bio_var_all<-cbind(bio_var_all,bio_var_all_centuries_mean)
    }#loop for each bioclimatic variable
  
  }#if the distribution is from several centuries, then we will calculate the average climate of those centuries: 
   bio_var_all<-data.frame(bio_var_all)
   colnames(bio_var_all)<-c("bio01","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09","bio10",
                                    "bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
  env_var_bio_var<-cbind(values_env_var,bio_var_all) 
  env_var_bio_var$Original_map<-maps_bear[m]
  env_var_bio_var_in_study_area<-subset(env_var_bio_var,Study_area==1)
  env_var_bio_var_in_study_area_all_ranges<-rbind(env_var_bio_var_in_study_area_all_ranges,env_var_bio_var_in_study_area)
  }#Loop for each distribution map
  
  save(env_var_bio_var_in_study_area_all_ranges, file="env_var_bio_var_in_study_area_all_ranges.RData")#output
  summary(env_var_bio_var_in_study_area_all_ranges)
  
  env_var_bio_var_in_study_area_all_ranges2<-env_var_bio_var_in_study_area_all_ranges
  env_var_bio_var_in_study_area_all_ranges2$Buf_min[is.na(env_var_bio_var_in_study_area_all_ranges2$Buf_min)]<-0
  env_var_bio_var_in_study_area_all_ranges2$Buf_max[is.na(env_var_bio_var_in_study_area_all_ranges2$Buf_max)]<-0
  env_var_bio_var_in_study_area_all_ranges2$Range[is.na(env_var_bio_var_in_study_area_all_ranges2$Range)]<-0
  env_var_bio_var_in_study_area_all_ranges2$Study_area[is.na(env_var_bio_var_in_study_area_all_ranges2$Study_area)]<-0
  summary(env_var_bio_var_in_study_area_all_ranges2)
  nrow(env_var_bio_var_in_study_area_all_ranges2)
  
  #Subset within study areas (already subsetted)
  area_within_study_area<-subset(env_var_bio_var_in_study_area_all_ranges2,Study_area==1)
  nrow(area_within_study_area)#34552
  #Select the background area
    #Subset outside buffer min
    area_outside_buf_min<-subset(area_within_study_area,Buf_min==0)
    nrow(area_outside_buf_min)#22265
    #Subset inside buffer max
    background_area<-subset(area_outside_buf_min,Buf_max==1)
    nrow(background_area)#16924
  #Select the presences
    presences_range_scale<-subset(area_within_study_area,Range==1)
    nrow(presences_range_scale)#10831
  presences_absences_eurasia_hist<-data.frame(rbind(presences_range_scale,background_area))
    nrow(presences_absences_eurasia_hist)#27755
  presences_absences_eurasia_hist<-na.omit(presences_absences_eurasia_hist)
    nrow(presences_absences_eurasia_hist)#27644
  
  #Count pixels of presences and available in the background area:
    print(paste0("Number of pixels with presence ",nrow(subset(presences_absences_eurasia_hist,Range==1))))
    print(paste0("Number of pixels in background area ",nrow(subset(presences_absences_eurasia_hist,Range==0))))
  #Save the presences/pseudoabsences  
    save(presences_absences_eurasia_hist, file="presences_absences_eurasia_hist.RData")



