

############################################################################################################
#Readme:
############################################################################################################
#R code for download occurrences from GBIF for the species in the brown bear diet
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input:
  #/2_Representative_Food-web_and_Food_species_list/Species_list_clean_binary.xlsx # The list of species
#Data output:
  #df_datos_1_276_description_GBIF.RData

##############################################################################################                
#Schema
##############################################################################################  
#3_Download occurrences from GBIF
  #3.1 Download a satellite image from google maps for use as background in maps
  #3.2 Download of occurrences from GBIF
  #3.3 Join dataframes with results of GBIF downloads
  #3.4 Rerun the data from GBIF. ATENTION!!! This part is only for rerun a subset of species. Same code that before but for a ubset
  #3.5 Join the data of imitial run (joined in point 2 of this code) and the rerun info (point 3)
  #3.6 Create a derived dataset to reference the databases used in GBIF
    #3.6.1 Save a table with the downloaded occurrences from each dataset on GBIF (raw data without filtering)
    #3.6.2 Save a table with the filtered occurrences
    #3.6.3 Go to Zenodo (https://zenodo.org/) and upload the data
    #3.6.4 Go to GBIF and Register the Derived Dataset
    #3.6.5 HOW TO CITE/REFERENCE
  #3.7 Transformation of descriptive data of GBIF to eliminate NA
  #3.8 Descriptive analysis for data downloaded from GBIF
  #3.9 Creation of rasters in idrisi format with the species presences
  #3.10 Code for calculate buffers with macro command of Idrisi
    #3.10.1 Code the instructions for buffer 2000 for macro command of Idrisi
    #3.10.2 Code the instructions for buffer 5000 for macro command of Idrisi
    #3.10.3 Code the instructions for buffer 10000 for macro command of Idrisi

############################################################################################################################
#3.1 Download a satellite image from google maps for use as background in maps
############################################################################################################################
 
rm(list=ls()) 
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")


  # Required libraries ------------------------------------------------------------------
    #Install required packages 
    if(!requireNamespace("taxize", quietly=TRUE))
      install.packages("taxize", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("reshape2", quietly=TRUE))
      install.packages("reshape2", quiet=TRUE, dependencies=TRUE)   
    if(!requireNamespace("ggmap", quietly=TRUE))
      install.packages("ggmap", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("maps", quietly=TRUE))
      install.packages("maps", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("ggplot2", quietly=TRUE))
      install.packages("ggplot2", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("mapdata", quietly=TRUE))
      install.packages("mapdata", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("grid", quietly=TRUE))
      install.packages("grid", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("maptools", quietly=TRUE))
      install.packages("maptools", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("rgdal", quietly=TRUE))
      install.packages("rgdal", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("sp", quietly=TRUE))
      install.packages("sp", quiet=TRUE, dependencies=TRUE)    
    if(!requireNamespace("rgbif", quietly=TRUE))
      install.packages("rgbif", quiet=TRUE, dependencies=TRUE)    
    if(!requireNamespace("raster", quietly=TRUE))
      install.packages("raster", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("readxl", quietly=TRUE))
      install.packages("readxl", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("devtools", quietly=TRUE))
      install.packages("devtools", quiet=TRUE, dependencies=TRUE)    
    if(!requireNamespace("dplyr", quietly=TRUE))
      install.packages("dplyr", quiet=TRUE, dependencies=TRUE)    
    if(!requireNamespace("CoordinateCleaner", quietly=TRUE))
      install.packages("CoordinateCleaner", quiet=TRUE, dependencies=TRUE)    

    library(taxize)
    library(reshape2)
    library(ggmap)
    library(maps)
    library(ggplot2)
    library(mapdata)
    library(grid)  
    library(maptools)
    library(rgdal)
    library(sp)
    library(rgbif)     
    library(raster) 
    library(readxl) 
    library(devtools)
    library(dplyr)
    library(CoordinateCleaner)    

   #First we are going to download a satellite image from google maps for use as background image for create maps for each species
           #We are going to create maps with all the points
            #We have obtained an api key for google maps:
            key<-"writeyourgooglemapkey" 
            #We need to register it each time:
            register_google(key = key)
            ### Set a range for the map that we want to download
                lat <- c(30, 71)                
                lon <- c(-12, 55)  
            #We get the map from google maps with the hybrid option
                #We use get_map for download from google, it is not a map, it is an image
                Europe_map <- get_map(location = c(lon=mean(lon), lat=mean(lat)), zoom = 2,
                maptype = "hybrid", source = "google",api_key = key)
                save(Europe_map, file=paste0(folder_gbif,"/Europe_map.RData"))
                #We convert to a map with this process
                gmap<-Europe_map
                mgmap <- as.matrix(gmap)
                vgmap <- as.vector(mgmap)
                vgmaprgb <- col2rgb(vgmap)
                gmapr <- matrix(vgmaprgb[1, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
                gmapg <- matrix(vgmaprgb[2, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
                gmapb <- matrix(vgmaprgb[3, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
                rgmaprgb <- brick(raster(gmapr), raster(gmapg), raster(gmapb))
                rm(gmapr, gmapg, gmapb)
                rgmaprgbGM <- rgmaprgb
                #1.3 Assume output of get_map() is in ?Google PseudoMercator? (epsg:3857) POSSIBLE
                projection(rgmaprgbGM) <- CRS("+init=epsg:3857")
                #In order to set the extent (required to define a correct raster layer), we have to find out the extent in the units of the new projection
                #We create an SpPolDF and then use spTransform
                unlist(attr(gmap, which = "bb"))[c(2, 4, 1, 3)]
                rprobextSpDF <- as(extent(unlist(attr(gmap, which = "bb"))[c(2, 4, 1, 3)]),"SpatialPolygons")
                projection(rprobextSpDF) <- CRS("+init=epsg:4326")
                rprobextGM <- spTransform(rprobextSpDF, CRS("+init=epsg:3857"))
                #We assign the values of the extension to the map
                extent(rgmaprgbGM) <- c(rprobextGM@bbox[1, ], rprobextGM@bbox[2, ])
                #We have created a spatial point object with the projection  WGS84: mapa_sub_coord_uncert_564_df and we are going to define other projection and to project in this new one  
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                pr2 <- projectRaster(rgmaprgbGM, crs=newproj)
                e2<- extent(-2800000,3100000,-800000,4800000) 
                #We crop the satellite image
                pr2_crop <- crop(pr2,e2) 
                save(pr2_crop, file='G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData')

############################################################################################################################
#3.2 Download of occurrences from GBIF
############################################################################################################################
  
rm(list=ls()) 
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")

    if (!dir.exists(folder_data_ouput_Occ_search)){
      dir.create(folder_data_ouput_Occ_search)
      print("Folder created")
    }else{
      print("Folder exists")
    }
   if (!dir.exists(folder_data_ouput_Occ_search_Occurrences)){
      dir.create(folder_data_ouput_Occ_search_Occurrences)
      print("Folder created")
    }else{
      print("Folder exists")
    }
   if (!dir.exists(folder_data_ouput_Occ_search_Shapes)){
      dir.create(folder_data_ouput_Occ_search_Shapes)
      print("Folder created")
    }else{
      print("Folder exists")
    }
    if (!dir.exists(folder_data_ouput_Occ_search_Map_occurrences)){
      dir.create(folder_data_ouput_Occ_search_Map_occurrences)
      print("Folder created")
    }else{
      print("Folder exists")
    }    
    if (!dir.exists(folder_data_ouput_Occ_search_Rasters_in_vector)){
      dir.create(folder_data_ouput_Occ_search_Rasters_in_vector)
      print("Folder created")
    }else{
      print("Folder exists")
    }    
      
    
   #We load the satellite image for background maps created in part 3  
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
  #we load the list of all species
  sheet<-read_excel("H:/G/Project_Name/Database_diet/R_analysis/Descriptive_results/Species_list_clean_binary.xlsx", sheet = 1)
  data_taxo<-as.data.frame(sheet)
  head(data_taxo)
  str(data_taxo)

  #Check for non repetitons in the list of species
    str(data_taxo)#There are 276 lines which should be 276 species
    list_sps<-unique(data_taxo$Species)#We check that there are 276 unique species names

  #We have created a spatial point object with the projection  WGS84: mapa_sub_coord_uncert_564_df and we are going to define other projection and to project in this new one  
  #We are going to use Europe Albers Equal Area Conic  
  newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
  str(data_taxo)#276 speceis for the loop
  DF_spatial_data_sps<-c()
  DF_n_data_by_sps<-c()
  GBIF_data_by_sps<-c()
  list_sps_in_loop<-c()
  #check_acceptance_in_GBIF_df<-c()
  #Loop for download GBIF data
      X11()
      print("############################################################################################################################")
      print("############################################################################################################################")
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print(paste0("             ######## STARTING LOOP FOR DOWNLOAD DATA FROM GBIF AT ",Sys.time(),"########          "))
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print("############################################################################################################################")
      print("############################################################################################################################")
      for (i in 264:276){#276
        sps_in_loop<-data_taxo[i,]
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Starting Species ",sps_in_loop$Species,"  ### "))
        vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
        nam_sps<-as.data.frame(sps_in_loop$Species)
            key<-c(NA)
            num_occur<-c(NA)
            data_by_sps_in_loop<-c()
            Sys.sleep(2)    
            name_list<-c(NA)
            occ_search<-c(NA)
            sub_occ_search<-c(NA)
            sub_year<-c(NA)
            dat_clean<-c(NA)
            n_dat_clean<-c(NA)
        try({
        if (sps_in_loop$Rank_taxo=="species"){    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "species", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        } 
        if (sps_in_loop$Rank_taxo=="subspecies") {    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "subspecies", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        }  
        if (sps_in_loop$Rank_taxo=="form") {    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "form", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        }    
        })    
        Sys.sleep(2)    
        num_occur<-NA
        try({ 
            #Get number of occurrences records:  
            num_occur<-occ_count(taxonKey = key)#,basisOfRecord='HUMAN_OBSERVATION', georeferenced=TRUE
        }) 
        Sys.sleep(2) 
        n_dat_clean<-NA
        sum_flags_data<-t(as.data.frame(c(NA,NA,NA,NA,NA,NA,NA,NA)))
        colnames(sum_flags_data)<-c(".val", ".equ",".zer",".cap",".cen",".gbf",".inst",".summary")
        rownames(sum_flags_data)<-c("summary(flags)")
        if (num_occur>0){
        try({ 
            occ_search<-occ_search(taxonKey = key, return='data', limit = 50000, hasCoordinate=TRUE,decimalLatitude='15,75',decimalLongitude='-20,105') #,fields=c("species","coordinateUncertaintyInMeters","coordinatePrecision","year","decimalLongitude","decimalLatitude","basisOfRecord","countryCode")
            assign( paste( "occ_search_sps_", i, sep=""),  occ_search)
            save(list=paste("occ_search_sps_", i, sep=""), file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Occ_search/occ_search_sps_", i, ".Rdata",sep=""))  
            })
        data_by_sps_in_loop<-c(NA,NA,NA)
        try({
            sub_year<-subset(occ_search,year>1989)
            sub_europe_lat<-subset(sub_year,decimalLongitude>-20&decimalLongitude<105)
            sub_europe_lon<-subset(sub_europe_lat,decimalLatitude>15&decimalLatitude<75)
            #We subset by coordinate uncertainity https://terms.tdwg.org/wiki/dwc:coordinateUncertaintyInMeters    
            # or by coordinate precission https://terms.tdwg.org/wiki/dwc:coordinatePrecision
                #  1 decimal degree = 85.27990 km en latitud Madrid
                #  1 decimal degree = 85279.90 m
                #  564.1 m = 0.006636471 decimal degree
             attach(sub_europe_lon)
              if (exists("coordinatePrecision")){
              dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1 |  coordinatePrecision<0.006636471)
              } else {
              dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1)
              }
              detach(sub_europe_lon)
              #Basis of record https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
              sub_basis<-subset(dat_sub_uncert, basisOfRecord == "HUMAN_OBSERVATION" |  basisOfRecord == "LIVING_SPECIMEN" |  basisOfRecord == "MACHINE_OBSERVATION" |  basisOfRecord == "OBSERVATION")
            ####  COORDINATE CLEANER  ######
            #We are going to use coordinate cleaner, We have followed the turorial https://ropensci.github.io/CoordinateCleaner/articles/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
            dat<-sub_basis
            #plot data to get an overview
                #X11()
                #wm <- borders("world", colour="gray50", fill="gray50")
                #ggplot()+ coord_fixed()+ wm +
                #  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
                  #           colour = "darkred", size = 0.5)+
                  #theme_bw()
            #Option A) Using the clean_coordinates wrapper function
                #flag problems
                dat <- data.frame(dat)
        
                #funtion https://ropensci.github.io/CoordinateCleaner/reference/clean_coordinates.html
                flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
                                          countries = "countryCode", 
                                          species = "species",
                                          tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                                    "zeros")) # most test are on by default
                sum_flags_data<-t(as.data.frame(summary(flags)))
                #We take the clean data:
                dat_clean <- dat[flags$.summary,]
             }) 
              #We are going to create a shapefile of points for 564 uncert:
              n_dat_clean<-NA
                  try({ 
                  Y564<-dat_clean$decimalLatitude
                  X564<-dat_clean$decimalLongitude
                  mapa_dat_clean<-as.data.frame(cbind(X564,Y564))
                  class(mapa_dat_clean)
                  coordinates(mapa_dat_clean)=~X564+Y564
                  class(mapa_dat_clean)
                  proj4string(mapa_dat_clean) = CRS("+init=epsg:4326")
      
                  #We change the projection of the presences
                  mapa_dat_clean_df_EAEAC <- spTransform(mapa_dat_clean, newproj)
                  mapa_dat_clean_df_EAEAC_df=SpatialPointsDataFrame(mapa_dat_clean_df_EAEAC, data.frame(id=1:length(mapa_dat_clean_df_EAEAC)))
                  writeOGR(obj=mapa_dat_clean_df_EAEAC_df, "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Shapes",layer=paste("Occurrences_sps_", i, sep=""), "ESRI Shapefile",overwrite_layer=TRUE)
                  df_mapa_dat_clean_df_EAEAC<-as.data.frame(mapa_dat_clean_df_EAEAC)
                  n_dat_clean<-nrow(df_mapa_dat_clean_df_EAEAC)
                  }) 
            #Map for 1sq km data from GBIF
            try({
            plotRGB(pr2_crop) 
            })
            try({
            points(df_mapa_dat_clean_df_EAEAC$X564, df_mapa_dat_clean_df_EAEAC$Y564, col = "black", cex = 0.8,pch = 19) 
            })
            try({
            dev.copy(jpeg,file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Map_occurrences/Map_Occurrences_species_",i,".jpeg",sep=""))
            })
            dev.off()
            #We rasterize the data          
            #Raster of reference: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            #Rasterize the points of species selected
            try({
            sps_raster<-rasterize(mapa_dat_clean_df_EAEAC, new_raster, field=1)
            #pr <- projectRaster(m1_crop, new_raster, method='bilinear')
                e2<- extent(-2800000,3100000,-800000,4800000) 
                vec_sps<-extract(sps_raster,e2)
            #We save the raster:
                assign( paste( "Vector_sps_", i, sep=""),  vec_sps)
                save(list=paste("Vector_sps_", i, sep=""), file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters_in_vector/Vector_sps_", i, ".Rdata",sep=""))  
            })
            #We paste this code again here for species like Atadinus alpinus which we changed the name
            sps_in_loop<-data_taxo[i,]
            vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
            nam_sps<-as.data.frame(sps_in_loop$Species)
            #We count the number of presences:
            n_presences_sps<-NA
            try({
            no_na_vec<-na.omit(vec_sps)
            n_presences_sps<-length(no_na_vec)
            })
            DF_data_sps<-as.data.frame(cbind(i,vec_sps_name_in_loop,n_presences_sps,sum_flags_data))
            try({
            data_by_sps_in_loop<-cbind(key,num_occur,n_dat_clean)
            })

        } else {
        n_presences_sps<- 0  
        sum_flags_data<-t(as.data.frame(c(NA,NA,NA,NA,NA,NA,NA,NA)))
        colnames(sum_flags_data)<-c(".val", ".equ",".zer",".cap",".cen",".gbf",".inst",".summary")
        rownames(sum_flags_data)<-c("summary(flags)")
        n_dat_clean<-0
        DF_data_sps<-as.data.frame(cbind(i,vec_sps_name_in_loop,n_presences_sps,sum_flags_data))
          try({
              data_by_sps_in_loop<-cbind(key,num_occur,n_dat_clean)
          })
        }
        DF_n_data_by_sps<-rbind(DF_n_data_by_sps, DF_data_sps)
        GBIF_data_by_sps<-rbind(GBIF_data_by_sps,data_by_sps_in_loop)
        Sys.sleep(1) 
        list_sps_in_loop<-rbind(list_sps_in_loop,"Species"=vec_sps_name_in_loop)
                  print(paste0("FINISHED downloaded data for species ",sps_in_loop$Species))
                  print(paste0("             ######## Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
        rm(vec_sps,vec_sps_name_in_loop, data_by_sps_in_loop,no_na_vec, key,num_occur,n_dat_clean,df_mapa_dat_clean_df_EAEAC,Y564,X564,mapa_dat_clean,sum_flags_data,mapa_dat_clean_df_EAEAC,n_presences_sps)
    }
  GBIF_data_by_sps<-as.data.frame(GBIF_data_by_sps)
  df_datos<-cbind(DF_n_data_by_sps,GBIF_data_by_sps)
  df_datos_numeric<-df_datos_1_90
  df_datos_numeric$i <- as.numeric(as.character(df_datos_numeric$i))
  df_datos_numeric$n_presences_sps <- as.numeric(as.character(df_datos_numeric$n_presences_sps))
  df_datos_numeric$.val <- as.numeric(as.character(df_datos_numeric$.val))
  df_datos_numeric$.equ <- as.numeric(as.character(df_datos_numeric$.equ))
  df_datos_numeric$.zer <- as.numeric(as.character(df_datos_numeric$.zer))
  df_datos_numeric$.cap  <- as.numeric(as.character(df_datos_numeric$.cap))
  df_datos_numeric$.cen <- as.numeric(as.character(df_datos_numeric$.cen))
  df_datos_numeric$.gbf <- as.numeric(as.character(df_datos_numeric$.gbf))
  df_datos_numeric$.inst <- as.numeric(as.character(df_datos_numeric$.inst))
  df_datos_numeric$.summary <- as.numeric(as.character(df_datos_numeric$.summary))
  df_datos_numeric$key <- as.numeric(as.character(df_datos_numeric$key))
  df_datos_numeric$num_occur <- as.numeric(as.character(df_datos_numeric$num_occur))
  df_datos_numeric$n_dat_clean <- as.numeric(as.character(df_datos_numeric$n_dat_clean))
  rownames(df_datos_numeric)<-df_datos_numeric$i
  df_datos_numeric<-df_datos_numeric[,c(1,2,12,13,4,5,6,7,8,9,10,11,14,3)]
  #Description of data:
  #i: the ID for the species
  #Species. the species or taxonomic name (subspecies for wolf and dog)
  #GBIF_key: the key of the taxon in GBIF
  #Total_occur: the total occurrences in GBIF
  #"Coord_clean.val": the coordinates clean with var 
  #rtc
  #Filter_occur: the occurrences after filter by the area of study, spatial uncertainity and using coordinate cleaner
  #N_pixels_presence: the nember of pixeles with presence that we wil use to calculate the habitat
  colnames(df_datos_numeric)<-c("i","Species","GBIF_key","Total_occur","Coord_clean.val","Coord_clean.equ","Coord_clean.zer","Coord_clean.cap","Coord_clean.cen","Coord_clean.gbf","Coord_clean.inst","Coord_clean.summary","Filter_occur","N_pixels_presence")
  df_datos_numeric_1_90<-df_datos_numeric
  #From 1-90 save(df_datos_numeric_1_90, file='df_datos_numeric_1_90.RData')
  #From 91-150 save(df_datos_numeric_91_150, file='df_datos_numeric_91_150.RData')
  #From 151-210 save(df_datos_numeric_151_210, file='df_datos_numeric_151_210.RData')
  #From 211-263 save(df_datos_numeric_211_263, file='df_datos_numeric_211_263.RData')
  #From 264-276 save(df_datos_numeric_264_276, file='df_datos_numeric_264_276.RData')
  #load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/GBIF_data_by_sps2.RData")#This is the old data


############################################################################################################################
#3.3 Join dataframes with results of GBIF downloads
############################################################################################################################
  rm(list=ls()) 
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #install.packages("rgbif")   
  #install.packages("ggmap")
  #install.packages("maps")
  #install.packages("mapdata")
  #install.packages("maptools")
  library(ggmap)
  library(maps)
  library(ggplot2)
  library(mapdata)
  library(rgbif)     
  library(grid)  
  library(maptools)
  library(rgdal)
  library(sp)
  library(rgdal) 
  library(raster) 
  library(readxl) 
    #We load the all the dataframes with the number of presences/occurences downloaded from Gbif 
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_1_90.RData")
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_91_150.RData")
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_151_210.RData")
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_211_263.RData")
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_264_276.RData")
    df_datos_1_276_with_NA<-rbind(df_datos_numeric_1_90,df_datos_numeric_91_150,df_datos_numeric_151_210,df_datos_numeric_211_263,df_datos_numeric_264_276)  
    head(df_datos_1_276_with_NA,20)
    summary(df_datos_1_276_with_NA)  
    save(df_datos_1_276_with_NA, file='df_datos_1_276_with_NA.RData')
    sub_df_datos_1_276_with_NA<-subset(df_datos_1_276_with_NA,is.na(N_pixels_presence))  
    save(sub_df_datos_1_276_with_NA, file='sub_df_datos_1_276_with_NA.RData')
    write.csv(sub_df_datos_1_276_with_NA,file="sub_df_datos_1_276_with_NA.csv", row.names = F)


############################################################################################################################
#3.4 Rerun the data from GBIF. ATENTION!!! This part is only for rerun a subset of species. Same code that before but for a ubset
############################################################################################################################
  #We are adding some information for rerun the code for some selected species  
  rm(list=ls()) 
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #install.packages("rgbif")   
  #install.packages("ggmap")
  #install.packages("maps")
  #install.packages("mapdata")
  #install.packages("maptools")
  #install.packages("readxl")
  #install_github("ropensci/CoordinateCleaner")
  #install.packages("devtools")
  #install.packages("countrycode")
  library(ggmap)
  library(maps)
  library(ggplot2)
  library(mapdata)
  library(rgbif)     
  library(grid)  
  library(maptools)
  library(rgdal)
  library(sp)
  library(rgdal) 
  library(raster) 
  library(readxl) 
  library(devtools)
  library(dplyr)
  library(CoordinateCleaner)    
  #We load the satellite image for background maps created in part 3  
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Mapa_refencia/pr2_crop.RData")
  #we load the list of all species
  sheet<-read_excel("G:/Project_Name/Database_diet/R_analysis/Descriptive_results/Species_list_clean_binary.xlsx", sheet = 1)
  
  data_taxo<-as.data.frame(sheet)
  
  #Check for non repetitons in the list of species
    str(data_taxo)#There are 276 lines which should be 276 species
    list_sps<-unique(data_taxo$Species)#We check that there are 276 unique species names
    df_tab<-data.frame(table(data_taxo$Species))
    sub_df_tab<-subset(df_tab,Freq>1)#We check that any of the 276 speceis names are repeated
    str(list_sps)
    
  #############################  
  #ATENTION!!! This part is only for rerun a subset of species. This and 
    #the part "for (j in 1:11){#11" is the difference with the previous code   
    #We are going to rerun some species which we have NA values because 
    #we had this error when running the code "Error : 500 - Server error", probably due to some error in the servers or internet 
    #After run for first time all the species, we selected the speceis with NA in "N_pixels_presence", this can be due to a server error 
      #or due to they really do not have pixels with presences
    #We are going to load this subset of species with NA:
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/sub_df_datos_1_276_with_NA.RData")
    sub_df_datos_1_276_with_NA$Species  
    list_with_ID_i_species<-(sub_df_datos_1_276_with_NA$i)
  
  #############################  
       
  #We have created a spatial point object with the projection  WGS84: mapa_sub_coord_uncert_564_df and we are going to define other projection and to project in this new one  
  #We are going to use Europe Albers Equal Area Conic  
  newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 

  str(data_taxo)#276 speceis for the loop
  DF_spatial_data_sps<-c()
  DF_n_data_by_sps<-c()
  GBIF_data_by_sps<-c()
  list_sps_in_loop<-c()
  #check_acceptance_in_GBIF_df<-c()
  #Loop for download GBIF data

      X11()
      print("############################################################################################################################")
      print("############################################################################################################################")
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print(paste0("             ######## STARTING LOOP FOR DOWNLOAD DATA FROM GBIF AT ",Sys.time(),"########          "))
      print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
      print("############################################################################################################################")
      print("############################################################################################################################")
                  
      for (j in 1:11){#11
        i<-list_with_ID_i_species[j]#we have added this for select the subset of species with NA

        sps_in_loop<-data_taxo[i,]
            print("                       ")
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop j = ",j," ## Starting Species of the Subset of species with NA  ### "))
            print(paste0("Loop i = ",i," ## Starting Species ",sps_in_loop$Species,"  ### "))
        vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
        nam_sps<-as.data.frame(sps_in_loop$Species)
            key<-c(NA)
            num_occur<-c(NA)
            data_by_sps_in_loop<-c()
            Sys.sleep(2)    
            name_list<-c(NA)
            occ_search<-c(NA)
            sub_occ_search<-c(NA)
            sub_year<-c(NA)
            dat_clean<-c(NA)
            n_dat_clean<-c(NA)
        try({
        if (sps_in_loop$Rank_taxo=="species"){    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "species", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        } 
        if (sps_in_loop$Rank_taxo=="subspecies") {    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "subspecies", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        }  
        if (sps_in_loop$Rank_taxo=="form") {    
            name_list<-name_backbone(name=sps_in_loop$Species, rank = "form", kingdom = sps_in_loop$kingdom, phylum = sps_in_loop$phylum,class = sps_in_loop$class, order = sps_in_loop$order, family = sps_in_loop$family, genus = sps_in_loop$genus, verbose=TRUE, strict=TRUE)#sps_in_loop$Species
            key<-name_list$data$usageKey
        }    
        })    
        Sys.sleep(2)    
        num_occur<-NA
        try({ 
            #Get number of occurrences records:  
            num_occur<-occ_count(taxonKey = key)#,basisOfRecord='HUMAN_OBSERVATION', georeferenced=TRUE
        }) 
        
        Sys.sleep(2) 
        n_dat_clean<-NA
        sum_flags_data<-t(as.data.frame(c(NA,NA,NA,NA,NA,NA,NA,NA)))
        colnames(sum_flags_data)<-c(".val", ".equ",".zer",".cap",".cen",".gbf",".inst",".summary")
        rownames(sum_flags_data)<-c("summary(flags)")
        if (num_occur>0){
        try({ 
            occ_search<-occ_search(taxonKey = key, return='data', limit = 50000, hasCoordinate=TRUE,decimalLatitude='15,75',decimalLongitude='-20,105') #,fields=c("species","coordinateUncertaintyInMeters","coordinatePrecision","year","decimalLongitude","decimalLatitude","basisOfRecord","countryCode")
            assign( paste( "occ_search_sps_", i, sep=""),  occ_search)
            save(list=paste("occ_search_sps_", i, sep=""), file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Occ_search/occ_search_sps_", i, ".Rdata",sep=""))  
            })
        data_by_sps_in_loop<-c(NA,NA,NA)
        try({
            sub_year<-subset(occ_search,year>1989)
            sub_europe_lat<-subset(sub_year,decimalLongitude>-20&decimalLongitude<105)
            sub_europe_lon<-subset(sub_europe_lat,decimalLatitude>15&decimalLatitude<75)
            #We subset by coordinate uncertainity https://terms.tdwg.org/wiki/dwc:coordinateUncertaintyInMeters    
            # or by coordinate precission https://terms.tdwg.org/wiki/dwc:coordinatePrecision
                #  1 decimal degree = 85.27990 km en latitud Madrid
                #  1 decimal degree = 85279.90 m
                #  564.1 m = 0.006636471 decimal degree
             attach(sub_europe_lon)
              if (exists("coordinatePrecision")){
              dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1 |  coordinatePrecision<0.006636471)
              } else {
              dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1)
              }
              detach(sub_europe_lon)
              #Basis of record https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
              sub_basis<-subset(dat_sub_uncert, basisOfRecord == "HUMAN_OBSERVATION" |  basisOfRecord == "LIVING_SPECIMEN" |  basisOfRecord == "MACHINE_OBSERVATION" |  basisOfRecord == "OBSERVATION")
            ####  COORDINATE CLEANER  ######
            #We are going to use coordinate cleaner, We have followed the turorial https://ropensci.github.io/CoordinateCleaner/articles/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
            dat<-sub_basis
            #plot data to get an overview
                #X11()
                #wm <- borders("world", colour="gray50", fill="gray50")
                #ggplot()+ coord_fixed()+ wm +
                #  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
                  #           colour = "darkred", size = 0.5)+
                  #theme_bw()
            #Option A) Using the clean_coordinates wrapper function
                #flag problems
                dat <- data.frame(dat)
        
                #funtion https://ropensci.github.io/CoordinateCleaner/reference/clean_coordinates.html
                flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
                                          countries = "countryCode", 
                                          species = "species",
                                          tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                                    "zeros")) # most test are on by default
                sum_flags_data<-t(as.data.frame(summary(flags)))
                #We take the clean data:
                dat_clean <- dat[flags$.summary,]
             }) 
        
              #We are going to create a shapefile of points for 564 uncert:
              n_dat_clean<-NA
                  try({ 
                  Y564<-dat_clean$decimalLatitude
                  X564<-dat_clean$decimalLongitude
                  mapa_dat_clean<-as.data.frame(cbind(X564,Y564))
                  class(mapa_dat_clean)
                  coordinates(mapa_dat_clean)=~X564+Y564
                  class(mapa_dat_clean)
                  proj4string(mapa_dat_clean) = CRS("+init=epsg:4326")
      
                  #We change the projection of the presences
                  mapa_dat_clean_df_EAEAC <- spTransform(mapa_dat_clean, newproj)
                  mapa_dat_clean_df_EAEAC_df=SpatialPointsDataFrame(mapa_dat_clean_df_EAEAC, data.frame(id=1:length(mapa_dat_clean_df_EAEAC)))
                  writeOGR(obj=mapa_dat_clean_df_EAEAC_df, "G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Shapes",layer=paste("Occurrences_sps_", i, sep=""), "ESRI Shapefile",overwrite_layer=TRUE)
                  df_mapa_dat_clean_df_EAEAC<-as.data.frame(mapa_dat_clean_df_EAEAC)
                  n_dat_clean<-nrow(df_mapa_dat_clean_df_EAEAC)
                  }) 
            
            #Map for 1sq km data from GBIF
            try({
            plotRGB(pr2_crop) 
            })
            try({
            points(df_mapa_dat_clean_df_EAEAC$X564, df_mapa_dat_clean_df_EAEAC$Y564, col = "black", cex = 0.8,pch = 19) 
            })
            try({
            dev.copy(jpeg,file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Map_occurrences/Map_Occurrences_species_",i,".jpeg",sep=""))
            })
            dev.off()
        #We rasterize the data          
            #Raster of reference: 
            new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)
            #We defined the projection of our raster:
            projection(new_raster) <- newproj
            #Rasterize the points of species selected
            try({
            sps_raster<-rasterize(mapa_dat_clean_df_EAEAC, new_raster, field=1)
            #pr <- projectRaster(m1_crop, new_raster, method='bilinear')
                e2<- extent(-2800000,3100000,-800000,4800000) 
                vec_sps<-extract(sps_raster,e2)
            #We save the raster:
                assign( paste( "Vector_sps_", i, sep=""),  vec_sps)
                save(list=paste("Vector_sps_", i, sep=""), file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters_in_vector/Vector_sps_", i, ".Rdata",sep=""))  
            })
            #We paste this code again here for species like Atadinus alpinus which we changed the name
            sps_in_loop<-data_taxo[i,]
            vec_sps_name_in_loop<-as.vector(sps_in_loop$Species)
            nam_sps<-as.data.frame(sps_in_loop$Species)
        #We count the number of presences:
            n_presences_sps<-NA
            try({
            no_na_vec<-na.omit(vec_sps)
            n_presences_sps<-length(no_na_vec)
            })
        DF_data_sps<-as.data.frame(cbind(i,vec_sps_name_in_loop,n_presences_sps,sum_flags_data))
            try({
            data_by_sps_in_loop<-cbind(key,num_occur,n_dat_clean)
            })

        } else {
        n_presences_sps<- 0  
        sum_flags_data<-t(as.data.frame(c(NA,NA,NA,NA,NA,NA,NA,NA)))
        colnames(sum_flags_data)<-c(".val", ".equ",".zer",".cap",".cen",".gbf",".inst",".summary")
        rownames(sum_flags_data)<-c("summary(flags)")
        n_dat_clean<-0
        DF_data_sps<-as.data.frame(cbind(i,vec_sps_name_in_loop,n_presences_sps,sum_flags_data))
          try({
              data_by_sps_in_loop<-cbind(key,num_occur,n_dat_clean)
          })
        }
        DF_n_data_by_sps<-rbind(DF_n_data_by_sps, DF_data_sps)
        GBIF_data_by_sps<-rbind(GBIF_data_by_sps,data_by_sps_in_loop)
        Sys.sleep(1) 

        list_sps_in_loop<-rbind(list_sps_in_loop,"Species"=vec_sps_name_in_loop)
                  print(paste0("FINISHED downloaded data for species ",sps_in_loop$Species))
                  print(paste0("             ######## Finished at ",Sys.time(),"########          "))
                  print("############################################################################################################################")
                  print("############################################################################################################################")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   FINISHED    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")
                  print("=-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=   =-=-=-=-=-=-=-=    =-=-=-=-=-=-=-=")

                  rm(vec_sps,vec_sps_name_in_loop, data_by_sps_in_loop,no_na_vec, key,num_occur,n_dat_clean,df_mapa_dat_clean_df_EAEAC,Y564,X564,mapa_dat_clean,sum_flags_data,mapa_dat_clean_df_EAEAC,n_presences_sps)
    }
  GBIF_data_by_sps<-as.data.frame(GBIF_data_by_sps)
  df_datos<-cbind(DF_n_data_by_sps,GBIF_data_by_sps)
  df_datos_numeric<-df_datos
  df_datos_numeric$i <- as.numeric(as.character(df_datos_numeric$i))
  df_datos_numeric$n_presences_sps <- as.numeric(as.character(df_datos_numeric$n_presences_sps))
  df_datos_numeric$.val <- as.numeric(as.character(df_datos_numeric$.val))
  df_datos_numeric$.equ <- as.numeric(as.character(df_datos_numeric$.equ))
  df_datos_numeric$.zer <- as.numeric(as.character(df_datos_numeric$.zer))
  df_datos_numeric$.cap  <- as.numeric(as.character(df_datos_numeric$.cap))
  df_datos_numeric$.cen <- as.numeric(as.character(df_datos_numeric$.cen))
  df_datos_numeric$.gbf <- as.numeric(as.character(df_datos_numeric$.gbf))
  df_datos_numeric$.inst <- as.numeric(as.character(df_datos_numeric$.inst))
  df_datos_numeric$.summary <- as.numeric(as.character(df_datos_numeric$.summary))
  df_datos_numeric$key <- as.numeric(as.character(df_datos_numeric$key))
  df_datos_numeric$num_occur <- as.numeric(as.character(df_datos_numeric$num_occur))
  df_datos_numeric$n_dat_clean <- as.numeric(as.character(df_datos_numeric$n_dat_clean))
  rownames(df_datos_numeric)<-df_datos_numeric$i
  df_datos_numeric<-df_datos_numeric[,c(1,2,12,13,4,5,6,7,8,9,10,11,14,3)]
  #Description of data:
  #i: the ID for the species
  #Species. the species or taxonomic name (subspecies for wolf and dog)
  #GBIF_key: the key of the taxon in GBIF
  #Total_occur: the total occurrences in GBIF
  #"Coord_clean.val": the coordinates clean with var 
  #rtc
  #Filter_occur: the occurrences after filter by the area of study, spatial uncertainity and using coordinate cleaner
  #N_pixels_presence: the nember of pixeles with presence that we wil use to calculate the habitat
  colnames(df_datos_numeric)<-c("i","Species","GBIF_key","Total_occur","Coord_clean.val","Coord_clean.equ","Coord_clean.zer","Coord_clean.cap","Coord_clean.cen","Coord_clean.gbf","Coord_clean.inst","Coord_clean.summary","Filter_occur","N_pixels_presence")
  df_datos_numeric_rerun_of_NA_sps<-df_datos_numeric
  #From 1-90 save(df_datos_numeric_1_90, file='df_datos_numeric_1_90.RData')
  #From 91-150 save(df_datos_numeric_91_150, file='df_datos_numeric_91_150.RData')
  #From 151-210 save(df_datos_numeric_151_210, file='df_datos_numeric_151_210.RData')
  #From 211-263 save(df_datos_numeric_211_263, file='df_datos_numeric_211_263.RData')
  #From 264-276 save(df_datos_numeric_264_276, file='df_datos_numeric_264_276.RData')
  save(df_datos_numeric_rerun_of_NA_sps, file='df_datos_numeric_rerun_of_NA_sps.RData')#Rerun of the NA species

   
##########################################################################################################    
#3.5 Join the data of imitial run (joined in point 2 of this code) and the rerun info (point 3)
##########################################################################################################    
  rm(list=ls()) 
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #Info of the initial run:
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_with_NA.RData")
  df_datos_1_276_with_NA_2<-df_datos_1_276_with_NA

  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_numeric_rerun_of_NA_sps.RData")
  head(df_datos_numeric_rerun_of_NA_sps)
  str(df_datos_numeric_rerun_of_NA_sps)#We need to replace 11 species
  
  df_datos_numeric_rerun_of_NA_sps[1,]
  df_datos_1_276_with_NA_2[(13),]
  df_datos_1_276_with_NA_2[13,]<-df_datos_numeric_rerun_of_NA_sps[1,]
  df_datos_1_276_with_NA_2[(13),]
  
  df_datos_numeric_rerun_of_NA_sps[2,]
  df_datos_1_276_with_NA_2[(16),]
  df_datos_1_276_with_NA_2[16,]<-df_datos_numeric_rerun_of_NA_sps[2,]
  df_datos_1_276_with_NA_2[(16),]
  
  df_datos_numeric_rerun_of_NA_sps[3,]
  df_datos_1_276_with_NA_2[(104),]
  df_datos_1_276_with_NA_2[104,]<-df_datos_numeric_rerun_of_NA_sps[3,]
  df_datos_1_276_with_NA_2[(104),]
  
  df_datos_numeric_rerun_of_NA_sps[4,]
  df_datos_1_276_with_NA_2[(115),]
  df_datos_1_276_with_NA_2[115,]<-df_datos_numeric_rerun_of_NA_sps[4,]
  df_datos_1_276_with_NA_2[(115),]
   
  df_datos_numeric_rerun_of_NA_sps[5,]
  df_datos_1_276_with_NA_2[(156),]
  df_datos_1_276_with_NA_2[156,]<-df_datos_numeric_rerun_of_NA_sps[5,]
  df_datos_1_276_with_NA_2[(156),]
  
  df_datos_numeric_rerun_of_NA_sps[6,]    
  df_datos_1_276_with_NA_2[(204),]
  df_datos_1_276_with_NA_2[204,]<-df_datos_numeric_rerun_of_NA_sps[6,]
  df_datos_1_276_with_NA_2[(204),]
  
  df_datos_numeric_rerun_of_NA_sps[7,]    
  df_datos_1_276_with_NA_2[(214),]
  df_datos_1_276_with_NA_2[214,]<-df_datos_numeric_rerun_of_NA_sps[7,]
  df_datos_1_276_with_NA_2[(214),]
  
  df_datos_numeric_rerun_of_NA_sps[8,]    
  df_datos_1_276_with_NA_2[(215),]
  df_datos_1_276_with_NA_2[215,]<-df_datos_numeric_rerun_of_NA_sps[8,]
  df_datos_1_276_with_NA_2[(215),]
  
  df_datos_numeric_rerun_of_NA_sps[9,]    
  df_datos_1_276_with_NA_2[(259),]
  df_datos_1_276_with_NA_2[259,]<-df_datos_numeric_rerun_of_NA_sps[9,]
  df_datos_1_276_with_NA_2[(259),]
  
  df_datos_numeric_rerun_of_NA_sps[10,]    
  df_datos_1_276_with_NA_2[(266),]
  df_datos_1_276_with_NA_2[266,]<-df_datos_numeric_rerun_of_NA_sps[10,]
  df_datos_1_276_with_NA_2[(266),]
  
  df_datos_numeric_rerun_of_NA_sps[11,]    
  df_datos_1_276_with_NA_2[(267),]
  df_datos_1_276_with_NA_2[267,]<-df_datos_numeric_rerun_of_NA_sps[11,]
  df_datos_1_276_with_NA_2[(267),]
        
  summary(df_datos_1_276_with_NA_2)  #5 species with NA no data after filtering  
  str(df_datos_1_276_with_NA_2)#276 rows with 276 levels of species
  save(df_datos_1_276_with_NA_2, file='df_datos_1_276_with_NA_2.RData')
  
  
############################################################################################################################
#3.6 Create a derived dataset to reference the databases used in GBIF
############################################################################################################################
  
  #3.6.1 Save a table with the downloaded occurrences from each dataset on GBIF (raw data without filtering)
    #load data for each species with the datasetkey
    occ_search_sps_all<-c()
    for(i in 1:276){#276
      print("#####")
      print(i)
     try({load(paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Occ_search/occ_search_sps_", i, ".Rdata",sep=""))
      occ_search_sps_in_loop<-get(paste0("occ_search_sps_",i)) 
      occ_search_sps_in_loop_selec<-occ_search_sps_in_loop[,c("key","scientificName","kingdom","phylum","order","genus","family","species","decimalLatitude","datasetKey")]
      occ_search_sps_all<-rbind(occ_search_sps_all,occ_search_sps_in_loop_selec)
      rm(occ_search_sps_in_loop,occ_search_sps_in_loop_selec)
     })
      print(nrow(occ_search_sps_all))
    }
    derived_dataset_pre<-occ_search_sps_all
    save(derived_dataset_pre, file=paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset_pre.Rdata")  
    derived_dataset<- data.frame(table(occ_search_sps_all$datasetKey))
    colnames(derived_dataset)<-c("datasetKey","n")  
    save(derived_dataset, file="H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset.Rdata")  
    write.csv(derived_dataset, file = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset.csv")
    write.table(derived_dataset, "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset.txt", append = FALSE, sep="\t", dec = ".",
                row.names = F, col.names = T, quote = FALSE)

  #3.6.2 Save a table with the filtered occurrences
    install.packages("CoordinateCleaner")
    library(CoordinateCleaner)   
    dat_clean_selec_all<-c()        
    for(i in 1:276){#276
      print("#####")
      print(i)
      try({load(paste("H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/Occ_search/occ_search_sps_", i, ".Rdata",sep=""))
        occ_search_sps_in_loop<-get(paste0("occ_search_sps_",i)) 
      })      
      try({
        sub_year<-subset(occ_search_sps_in_loop,year>1989)
        sub_europe_lat<-subset(sub_year,decimalLongitude>-20&decimalLongitude<105)
        sub_europe_lon<-subset(sub_europe_lat,decimalLatitude>15&decimalLatitude<75)
        #We subset by coordinate uncertainity https://terms.tdwg.org/wiki/dwc:coordinateUncertaintyInMeters    
        # or by coordinate precission https://terms.tdwg.org/wiki/dwc:coordinatePrecision
          #  1 decimal degree = 85.27990 km en latitud Madrid
          #  1 decimal degree = 85279.90 m
          #  564.1 m = 0.006636471 decimal degree
        attach(sub_europe_lon)
        if (exists("coordinatePrecision")){
          dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1 |  coordinatePrecision<0.006636471)
        } else {
          dat_sub_uncert<-subset(sub_europe_lon,coordinateUncertaintyInMeters<564.1)
        }
        detach(sub_europe_lon)
        #Basis of record https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/BasisOfRecord.html
        sub_basis<-subset(dat_sub_uncert, basisOfRecord == "HUMAN_OBSERVATION" |  basisOfRecord == "LIVING_SPECIMEN" |  basisOfRecord == "MACHINE_OBSERVATION" |  basisOfRecord == "OBSERVATION")
        ####  COORDINATE CLEANER  ######
        #We are going to use coordinate cleaner, We have followed the turorial https://ropensci.github.io/CoordinateCleaner/articles/Tutorial_Cleaning_GBIF_data_with_CoordinateCleaner.html
        dat<-sub_basis
        #Option A) Using the clean_coordinates wrapper function
        #flag problems
        dat <- data.frame(dat)
        #funtion https://ropensci.github.io/CoordinateCleaner/reference/clean_coordinates.html
        flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
          countries = "countryCode", 
          species = "species",
          tests = c("capitals", "centroids", "equal","gbif", "institutions",
          "zeros")) # most test are on by default
        sum_flags_data<-t(as.data.frame(summary(flags)))
        #We take the clean data:
        dat_clean <- dat[flags$.summary,]
        dat_clean_selec<-dat_clean[,c("key","scientificName","kingdom","phylum","order","genus","family","species","decimalLatitude","decimalLongitude","datasetKey")]
        dat_clean_selec_all<-rbind(dat_clean_selec_all,dat_clean_selec)
        rm(dat_clean_selec)
      }) 
     } 
    save(dat_clean_selec_all, file="H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/dat_clean_selec_all.Rdata")  
    write.csv(dat_clean_selec_all, file = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/dat_clean_selec_all.csv")

    derived_dataset_filtered<- data.frame(table(dat_clean_selec_all$datasetKey))
    colnames(derived_dataset_filtered)<-c("datasetKey","n")  
    save(derived_dataset_filtered, file="H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset_filtered.Rdata")  
    write.csv(derived_dataset_filtered, file = "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset_filtered.csv")
    write.table(derived_dataset_filtered, "H:/G/Project_Name/Database_diet/R_analysis/GBIF_data2/derived_dataset_filtered.txt", append = FALSE, sep="\t", dec = ".",
                row.names = F, col.names = T, quote = FALSE)

  #3.6.3 Go to Zenodo (https://zenodo.org/) and upload the data
      #You should also upload your filtered GBIF dataset of occurrences to a public repository like Zenodo.
      #We are going to upload the data used

  #3.6.4 Go to GBIF and Register the Derived Dataset
    
  #3.6.5 HOW TO CITE/REFERENCE
    #Citation in text:
    #"y cross-plotting plant modern distributions from the Global Biodiversity Information Facility (GBIF.org 2018)"
    
    #In References:
    #GBIF.org, 2018. 27 June. GBIF Occurrence Download. https://doi.org/10.15468/dl.j6ob9o.


#################################################################################################################    
#3.7 Transformation of descriptive data of GBIF to eliminate NA
#################################################################################################################    
  rm(list=ls()) 
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #Info of the initial run:
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_with_NA_2.RData")
  summary(df_datos_1_276_with_NA_2)
  df_datos_1_276_with_NA_2[is.na(df_datos_1_276_with_NA_2)] <- 0
  df_datos_1_276_description_GBIF<-df_datos_1_276_with_NA_2
  save(df_datos_1_276_description_GBIF, file='df_datos_1_276_description_GBIF.RData')

#################################################################################################################    
#3.8 Descriptive analysis for data downloaded from GBIF
#################################################################################################################    
  rm(list=ls()) 
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #Info of the initial run:
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    #We load the data of selected species (based in selected studies)
    merged_species_selected_cat <- read.csv("G:/Project_Name/Database_diet/R_analysis/Descriptive_results/merged_species_selected_cat.csv") # header = TRUE, na.strings = "<NA>"
    #We filter the data of human origin
    merged_species_selected_cat_no_human<-subset(merged_species_selected_cat,Human_origin<1)
    #We exclude brown bear
    merged_species_selected_cat_no_human<-subset(merged_species_selected_cat_no_human, ID!=174)
    merged_species_selected_cat_no_human$Species<-factor(merged_species_selected_cat_no_human$Species)
    table(merged_species_selected_cat_no_human$Category)

    merged_species_selected_cat_GBIF_data<-merge(merged_species_selected_cat_no_human, df_datos_1_158, by.x="Species", by.y = "Species",all.x = T)
    write.csv(file="merged_species_selected_cat_GBIF_data.csv", x=merged_species_selected_cat_GBIF_data, row.names = F)
    save(merged_species_selected_cat_GBIF_data, file='merged_species_selected_cat_GBIF_data.RData')

    X11()
    par(mfrow=c(1,2), mar=c(5,5,1,1))#
    par(xaxs="i",yaxs="i")
    par(mgp = c(2.0, 0.50, 0))
    par("tck"=-0.025)
    
    hist(log(merged_species_selected_cat_GBIF_data$N_1km2_occur_GBIF),xlim=c(0,11),ylim=c(0,50),breaks=seq(0,11,1), main="Occurrences", xlab="Log (N occurrences)", ylab="Frecuency of species",pch=1, cex=1, cex.axis=1, cex.lab=1, cex.main=1,xaxt='n')
        abline(a=NULL, b= NULL, v=log(mean(merged_species_selected_cat_GBIF_data$N_1km2_occur_GBIF)),lty=2, col="grey50",lwd=3)
        abline(a=NULL, b= NULL, v=log(median(merged_species_selected_cat_GBIF_data$N_1km2_occur_GBIF)),lty=3, col="grey50",lwd=3)
        axis(1, at=seq(0,11,1), labels=seq(0,11,1),cex.axis=0.9)#rep('', 11),
    
    hist(log(merged_species_selected_cat_GBIF_data$N_pixeles_presence),xlim=c(0,11),ylim=c(0,50),breaks=seq(0,11,1), main="Pixels with presence", xlab="Log (N pixels)", ylab="Frecuency of species",pch=1, cex=1, cex.axis=1, cex.lab=1, cex.main=1,xaxt='n')
        abline(a=NULL, b= NULL, v=log(mean(merged_species_selected_cat_GBIF_data$N_pixeles_presence)),lty=2, col="grey50",lwd=3)
        abline(a=NULL, b= NULL, v=log(median(merged_species_selected_cat_GBIF_data$N_pixeles_presence)),lty=3, col="grey50",lwd=3)
        abline(a=NULL, b= NULL, v=log(150),lty=2, col="red",lwd=3)
        axis(1, at=seq(0,11,1), labels=seq(0,11,1),cex.axis=0.9)#rep('', 11),

        
#################################################################################################################    
#3.9 Creation of rasters in idrisi format with the species presences
#################################################################################################################    
  rm(list=ls()) 
  library(ggmap)
  library(maps)
  library(ggplot2)
  library(mapdata)
  library(rgbif)     
  library(grid)  
  library(maptools)
  library(rgdal)
  library(sp)
  library(rgdal) 
  library(raster) 
  library(readxl) 
  library(devtools)
  library(dplyr)
  library(CoordinateCleaner)    
  ##Establecemos el fichero de trabajo
  setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
  #Info of the initial run:
  load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
        for (i in 1:276){#The number of species to model Originally
          list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
          #Condition of a minimum number of presences below n=150 we consider that we have not enought data
          list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
          if ((list_sps_in_loop$N_pixels_presence)>=50){
            print("             ##############################################################          ")
            print("             ########           Starting a new species             ########          ")
            print("             ##############################################################          ")
            print(paste0("Loop i = ",i," ## Creating raster for species ",list_sps_in_loop$Species,"  ###   "))
            #We select the species that will be run in the loop. This is a column with 1 for presences and Na for the other areas of the map (including sea and pseudbasences)            
            #sp_1_europa<-DF_spatial_data_sps_2[,c(i)]
              file_load_sps=load(file=paste("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters_in_vector/Vector_sps_", list_sps_in_loop$i, ".Rdata",sep=""))       
              sp_1_europa = get(file_load_sps)
              sp_1_europa[is.na(sp_1_europa)] <- 0
             #We create a raster with the presences of the species
                #With this raster we will calculate the background area and pseudoabsences for models using buffers around each pixel with presences
                #Raster of reference, this is a raster with the spatial dimensions of the map in which you want to project your data: 
                new_raster <- raster(xmn=-2800000, xmx=3100000, ymn=-800000, ymx=4800000, ncols=5900, nrows=5600)#Raster que utlizamos de normal 
                #We are going to use Europe Albers Equal Area Conic  
                newproj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs" 
                #We defined the projection of our raster:
                projection(new_raster) <- newproj
                #Europe_presences_sum2[Europe_presences_sum2 > 1] <- 1
                values(new_raster)<-sp_1_europa 
                writeRaster(new_raster, filename=paste0("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters/Raster_sps_",i,".rst"), datatype='INT1U',overwrite=TRUE) 
            }
        }

       
#################################################################################################################    
#3.10 Code for calculate buffers with macro command of Idrisi
#################################################################################################################  
  
  #3.10.1 Code the instructions for buffer 2000 for macro command of Idrisi
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
    #Info of the initial run:
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    macroline_all<-c()
          for (i in 1:276){#The number of species to model Originally
            list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
            #Condition of a minimum number of presences below n=150 we consider that we have not enought data
            list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
            if ((list_sps_in_loop$N_pixels_presence)>=50){
                 macroline<-as.data.frame( paste0("buffer x G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters/Raster_sps_",i,".rst*G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Buffer_2000/Buffer_2000_sps_",i,".rst*1*1*0*2000"))
                 macroline_all<-rbind(macroline_all,macroline) 
              }
          }
      macroline_all_buffer_2000<-macroline_all
        save(macroline_all_buffer_2000,file="macroline_all_buffer_2000.RData")       
        head(macroline_all_buffer_2000)
        write.table(macroline_all_buffer_2000, file="macroline_all_buffer_2000.txt", append = FALSE, sep = "", quote=F,dec = ".", row.names = FALSE, col.names = FALSE)  

  #3.10.2 Code the instructions for buffer 5000 for macro command of Idrisi
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
    #Info of the initial run:
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    macroline_all<-c()
          for (i in 1:276){#The number of species to model Originally
            list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
            #Condition of a minimum number of presences below n=150 we consider that we have not enought data
            list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
            if ((list_sps_in_loop$N_pixels_presence)>=50){
   
                 macroline<-as.data.frame( paste0("buffer x G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters/Raster_sps_",i,".rst*G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Buffer_5000/Buffer_5000_sps_",i,".rst*1*1*0*5000"))
                 macroline_all<-rbind(macroline_all,macroline) 
              }
          }
        macroline_all_buffer_5000<-macroline_all
        save(macroline_all_buffer_5000,file="macroline_all_buffer_5000.RData")       
        head(macroline_all_buffer_5000)
        write.table(macroline_all_buffer_5000, file="macroline_all_buffer_5000.txt", append = FALSE, sep = "", quote=F,dec = ".", row.names = FALSE, col.names = FALSE)  

  #3.10.3 Code the instructions for buffer 10000 for macro command of Idrisi
    rm(list=ls()) 
    ##Establecemos el fichero de trabajo
    setwd("G:/Project_Name/Database_diet/R_analysis/GBIF_data2")
    #Info of the initial run:
    load("G:/Project_Name/Database_diet/R_analysis/GBIF_data2/df_datos_1_276_description_GBIF.RData")
    macroline_all<-c()
          for (i in 1:276){#The number of species to model Originally
            list_sps_in_loop<-df_datos_1_276_description_GBIF[c(i),]
            #Condition of a minimum number of presences below n=150 we consider that we have not enought data
            list_sps_in_loop$N_pixels_presence[is.na(list_sps_in_loop$N_pixels_presence)]<-0
            if ((list_sps_in_loop$N_pixels_presence)>=50){
   
                 macroline<-as.data.frame( paste0("buffer x G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Rasters/Raster_sps_",i,".rst*G:/Project_Name/Database_diet/R_analysis/GBIF_data2/Buffer_10000/Buffer_10000_sps_",i,".rst*1*1*0*10000"))
                 macroline_all<-rbind(macroline_all,macroline) 
              }
          }
        macroline_all_buffer_10000<-macroline_all
        save(macroline_all_buffer_10000,file="macroline_all_buffer_10000.RData")       
        head(macroline_all_buffer_10000)
        write.table(macroline_all_buffer_10000, file="macroline_all_buffer_10000.txt", append = FALSE, sep = "", quote=F,dec = ".", row.names = FALSE, col.names = FALSE)  
