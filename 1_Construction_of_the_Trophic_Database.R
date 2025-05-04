

############################################################################################################
#Readme:
############################################################################################################
#R code for construct the Trophic Database of brown bear in Europe
#Authors: Pablo M. Lucas
#Last update: 28/03/2023

##############################################################################################                
#Files on this folder
##############################################################################################                
#Data input to run script:
  #Example_taxo_diet_study.txt
  #Data_all_studies_in_a_table_csv6.csv

#Data output to other scripts:
  #/database_original_with_energy_imputed2.RData
  #/database_original_with_energy_imputed2_ALL_corrected_names2.csv
  #/database_original_with_energy_imputed2_ALL_corrected_names.xlsx

##############################################################################################                
#Schema
##############################################################################################  
#1_Construction_of_the_Trophic_Database
  #1.1 Retrieve taxonomic information for items in the diet
  #1.2 Imputation of energy
    #1.2.1 Data preparation
    #1.2.2 Imputation of relative volume
  #1.3 Manual modification of Supplementary Table 5 to create Supplementary Table 7
  #1.4 Manual modification of Supplementary Table 5 to create Supplementary Table 7
  
############################################################################################################################
#1.1 Retrieve taxonomic information for items in the diet
############################################################################################################################
    
#This code is to facilitate to add the taxonomic information for each item in the diet. It is designed to run for each
#study in the diet, not all studies at the same time
rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/1_Construction_of_the_Trophic_Database")
    setwd(folder_working)
 
  # Required libraries ------------------------------------------------------------------
    #Install required packages 
    if(!requireNamespace("taxize", quietly=TRUE))
      install.packages("taxize", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("reshape2", quietly=TRUE))
      install.packages("reshape2", quiet=TRUE, dependencies=TRUE)      
    library(taxize)
    library(reshape2)

    DF_vec_ranks<-c()
    ranks<-c()
    #We copy and paste the columns taxa and rank of the sheet that we want
    # read copied area in excel spreadsheet from clipboard into R (Macintosh version, code for windows is different)
    taxa <- read.table("Example_taxo_diet_study.txt",sep="\t", header = TRUE, na.strings = "<NA>")
    taxa[, 1] <- as.character(taxa[, 1])
    #Names of rank start with 
    vector_taxa<-c(taxa$taxon)
    # apply the function to the taxonlist
    # ATTEINTION: if several entries are available for a species, you need to select
    len_vector_taxa<-length(vector_taxa)
    taxa$ID<-c(1:len_vector_taxa)
    #taxa2<-taxa[1:5,] borrarlinea
    taxa$rank[taxa$rank=='NA'] <- NA
    sub_na_taxa <- taxa[is.na(taxa$rank),]
    sub_kingdom_taxa<-subset(taxa[complete.cases(taxa), ], rank == "kingdom")
    sub_vector_taxa<-subset(taxa[complete.cases(taxa), ], rank != "kingdom")
    
    #For data with NA in rank:
        len_sub_na_taxa<-length(sub_na_taxa$taxon)
        DF_vec_ranks<-c()
        vector_NA<-c()
        DF_vec_ranks_NA<-c()
        
        if (len_sub_na_taxa>0){
        for (i in 1:len_sub_na_taxa){
        vector_NA<-c(NA,NA,NA,NA,NA,NA,NA,NA)
              try({
            vector_NA[8]<-sub_na_taxa[i,"ID"]
            })
        DF_vec_ranks<-rbind(DF_vec_ranks,vector_NA)
      
        }
        DF_vec_ranks_NA<- DF_vec_ranks
        colnames(DF_vec_ranks_NA)<-c("kingdom", "phylum", "class", "order", "family", "genus", "species","ID")
        }
    
    #For data with kingdom in rank:
        len_sub_kingdom_taxa<-length(sub_kingdom_taxa$taxon)
        DF_vec_ranks<-c()
        vector_NA<-c()
        DF_vec_ranks_kingdom<-c()
        if (len_sub_kingdom_taxa>0){
         
          #For each item with kingdom level:
            for (i in 1:len_sub_kingdom_taxa){
              sub_kingdom_taxa_in_loop<-sub_kingdom_taxa[i,"taxon"]
              vector_NA<-c(NA,NA,NA,NA,NA,NA,NA,NA)
              vector_NA[8]<-sub_kingdom_taxa[i,"ID"]
            vector_NA[1]<-sub_kingdom_taxa_in_loop
            DF_vec_ranks_kingdom<-rbind(DF_vec_ranks_kingdom,vector_NA)
                 }
        #DF_vec_ranks_kingdom<-cbind(DF_vec_ranks,sub_kingdom_taxa$ID)
        DF_vec_ranks_kingdom<-as.data.frame(DF_vec_ranks_kingdom)
        colnames(DF_vec_ranks_kingdom)<-c("kingdom", "phylum", "class", "order", "family", "genus", "species","ID")
        }
        
    #For data other than NA or kingdom in rank:
    len_sub_vector_taxa<-length(sub_vector_taxa$taxon)   
    ranks<-c()
    DF_vec_ranks2<-c()

    ##CAUTION!!! Select only the loop and run only the loop, we need to select the taxon, if we select more than the loop will be an error   
     #Select from here    
        for (i in 1:len_sub_vector_taxa){#
          taxo_loop<-sub_vector_taxa[i,"taxon"]
          ####          taxo_loop<- "Helianthemum canum"
          try({
          #Interacts with a suite of web 'APIs' for taxonomic tasks:  
          pre_classes <- classification(taxo_loop, db = 'ncbi') # col: Catalogue of Life was used but changes 
          #in the API has made it unusable so we inlcude "ncbi" database instead to run the funtion; 
            #COL introduced rate limiting recently in 2019 - which has made the API essentially unusable;
            #CoL+ is coming soon and we'll incorporate it here when it's stable. See https://github.com/ropensci/colpluz
            #for the R implementation for CoL+
            sub_pre_classes<-pre_classes[taxo_loop]   
            #sub_pre_classes[c("name", "rank"),]
            df_sub_pre_classes<-as.data.frame(sub_pre_classes)
            colnames(df_sub_pre_classes)<-c("name","rank","id")
            m <- df_sub_pre_classes[,c("name","rank")]
            t_m<-t(m)
            DF_vec_ranks<-c()
            vector_NA<-c(NA,NA,NA,NA,NA,NA,NA,NA)
            DF_vec_ranks2<-rbind(DF_vec_ranks,vector_NA)
            colnames(t_m)<-t_m[2,]
            colnames(DF_vec_ranks2)<-c("kingdom", "phylum", "class", "order", "family", "genus", "species","ID")
            try({
              DF_vec_ranks2[1,"kingdom"]<-t_m[1,"kingdom"]
              DF_vec_ranks2[1,"phylum"]<-t_m[1,"phylum"]
              DF_vec_ranks2[1,"class"]<-t_m[1,"class"]
              DF_vec_ranks2[1,"order"]<-t_m[1,"order"]
              DF_vec_ranks2[1,"family"]<-t_m[1,"family"]
              DF_vec_ranks2[1,"genus"]<-t_m[1,"genus"]
              DF_vec_ranks2[1,"species"]<-t_m[1,"species"]
              DF_vec_ranks2[1,"ID"]<-sub_vector_taxa[i,"ID"]
            })
            ranks<-rbind(ranks, DF_vec_ranks2[1,])
          })
        }
      #To here (FINISH OF THE SELECTION)
    
    DF_ranks_normal<-as.data.frame(ranks)
    DF_ranks_normal$ID<-as.numeric(as.character(DF_ranks_normal$ID))
    DF_ranks_table<-DF_ranks_normal
    #We try to join the three created data frames with the taxonomic info. As we not allways will have three dataframes 
    # we are going to use try function:
        try({
          DF_ranks_table<-rbind(DF_ranks_table,DF_vec_ranks_kingdom)
        })
        
        try({
          DF_ranks_table<-rbind(DF_ranks_table,DF_vec_ranks_NA)
        }) 
    #We order with the ID
        DF_ranks_final2<-DF_ranks_table
        DF_ranks_final2$ID<-as.numeric(DF_ranks_final2$ID)
        str(DF_ranks_final2)
        DF_ranks_final3<-DF_ranks_final2[order(DF_ranks_final2$ID),]
        DF_ranks_final4<-DF_ranks_final3[,c(1:8)]
    #We save the file    
    write.csv(file="DF_ranks_final4.csv", x=DF_ranks_final4, row.names = F)

############################################################################################################################
#1.2 Imputation of energy
############################################################################################################################
  
rm(list=ls()) 
    my_dir <-"writehereyourpath" # path of the folder where your data are stored
    folder_working<-paste0(my_dir,"/1_Construction_of_the_Trophic_Database")
    setwd(folder_working)

  # Required libraries ------------------------------------------------------------------
    #Install required packages 
    if(!requireNamespace("readxl", quietly=TRUE))
      install.packages("readxl", quiet=TRUE, dependencies=TRUE)      
    if(!requireNamespace("reshape2", quietly=TRUE))
      install.packages("reshape2", quiet=TRUE, dependencies=TRUE)   
    if(!requireNamespace("sp", quietly=TRUE))
      install.packages("sp", quiet=TRUE, dependencies=TRUE)   
    if(!requireNamespace("raster", quietly=TRUE))
      install.packages("raster", quiet=TRUE, dependencies=TRUE)  
    if(!requireNamespace("lattice", quietly=TRUE))
      install.packages("lattice", quiet=TRUE, dependencies=TRUE)  
    if(!requireNamespace("MCMCglmm", quietly=TRUE))
      install.packages("MCMCglmm", quiet=TRUE, dependencies=TRUE)  

    library(readxl)
    library(reshape2)
    library(sp)
    library(raster)    
    library(lattice)
    library(MCMCglmm)    

    #1.2.1 Data preparation
    ################################################################################
    # Preparation of dataset for analysis of trophic interactions
    ################################################################################

    ################################################################################
    # set up geographic coordinate system
    ################################################################################
    WGS84.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    ################################################################################
    # load raw data
    ################################################################################
    if(! file.exists("trophicDatabase.RData")) {
        ##############################################################################
        ##############################################################################
        # load & process correction factor information
        ##############################################################################
        corrFact <- read.csv(paste0(folder_working,"/correction_factors.csv"), sep = ";")
        corrFact <- within(corrFact, {
          diet1 <- factor(diet1, levels = c("animals", "plants", "fungi_lichens_bryophytes_algae", "other"))
          diet2 <- factor(diet2, levels = c("vertebrates", "invertebrates", "unknown_animal_material", "vegetative_plant_material",
                                            "reproductive_plant_material", "unknown_plant_material", "fungi_lichens_bryophytes_algae", "other"))
          diet3 <- factor(diet3, levels = c("endotherm_vertebrates", "ectotherm_vertebrates", "fish_vertebrates", "unknown_vertebrates",
                                            "aquatic_invertebrates", "terrestrial_invertebrates", "unknown_invertebrates", "unknown_animal_material",
                                            "fruit", "seeds", "flowers_nectar_pollen", "unknown_reproductive_plant_material",
                                            "forbs_herbs_legumes", "grass", "roots_tubers_bulbs", "leaves_branches_bark", "unknown_vegetative_plant_material",
                                            "unknown_plant_material", "fungi_lichens_bryophytes_algae", "other"))
          diet4 <- as.character(diet3)
          diet4[diet4 %in% c("endotherm_vertebrates", "ectotherm_vertebrates", "fish_vertebrates", "unknown_vertebrates")] <- "vertebrates"
          diet4[diet4 %in% c("aquatic_invertebrates", "terrestrial_invertebrates", "unknown_invertebrates")] <- "invertebrates"
          diet4[diet4 %in% c("grass", "forbs_herbs_legumes", "roots_tubers_bulbs", "leaves_branches_bark",
                             "flowers_nectar_pollen", "unknown_vegetative_plant_material")] <- "vegetation"
          diet4[diet4 %in% c("fruit")] <- "fruit"
          diet4[diet4 %in% c("seeds")] <- "seeds"
          diet4[diet4 %in% c("fungi_lichens_bryophytes_algae", "unknown_animal_material", "unknown_plant_material",
                             "unknown_reproductive_plant_material", "other")] <- "other"
          diet4 <- factor(diet4, levels = c("vertebrates", "invertebrates", "vegetation", "fruit", "seeds", "other"))
        })
        corrFact <- subset(corrFact, select = c(diet1, diet2, diet3, diet4, cf1, cf2))
        str(corrFact)
        
      ##############################################################################
      # load & process metadata
      ##############################################################################
      sheet<-read_excel("meta_data_2b.xlsx", sheet = 1)
      metaData<-as.data.frame(sheet)
      sp::coordinates(metaData) <- ~ longitude_dd + latitude_dd
      raster::crs(metaData) <- raster::crs(WGS84.proj)
      metaData@data <- cbind(sp::coordinates(metaData), metaData@data)
      ##############################################################################
      # load & process raw trophic database
      ##############################################################################
      dataBase <- read.csv("Data_all_studies_in_a_table_csv6.csv") 
      dataBase<-dataBase[,c(2:21)]
      str(dataBase)#1356 obs of 20 variables
      dataBase$rF<-as.numeric(as.character(dataBase$rF))
      dataBase$rV<-as.numeric(as.character(dataBase$rV))
      dataBase$rIV<-as.numeric(as.character(dataBase$rIV))
      dataBase$rEDC<-as.numeric(as.character(dataBase$rEDC))
      dataBase$rEDEC<-as.numeric(as.character(dataBase$rEDEC))
        dataBase <- within(dataBase, {
          diet1 <- factor(diet1, levels = c("animals", "plants", "fungi_lichens_bryophytes_algae", "other"))
          diet2 <- factor(diet2, levels = c("vertebrates", "invertebrates", "unknown_animal_material", "vegetative_plant_material",
                                            "reproductive_plant_material", "unknown_plant_material", "fungi_lichens_bryophytes_algae", "other"))
          diet3 <- factor(diet3, levels = c("endotherm_vertebrates", "ectotherm_vertebrates", "fish_vertebrates", "unknown_vertebrates",
                                            "aquatic_invertebrates", "terrestrial_invertebrates", "unknown_invertebrates", "unknown_animal_material",
                                            "fruit", "seeds", "flowers_nectar_pollen", "unknown_reproductive_plant_material",
                                            "forbs_herbs_legumes", "grass", "roots_tubers_bulbs", "leaves_branches_bark", "unknown_vegetative_plant_material",
                                            "unknown_plant_material", "fungi_lichens_bryophytes_algae", "other"))
          diet4 <- as.character(diet1)
          diet4[diet4 %in% c("animals")] <- "animals"
          diet4[diet4 %in% c("plants", "fungi_lichens_bryophytes_algae")] <- "plants"
          diet4[diet4 %in% c("other")] <- "other"
          diet4 <- factor(diet4, levels = c("animals", "plants", "other"))
          diet5 <- as.character(diet3)
          diet5[diet5 %in% c("endotherm_vertebrates", "ectotherm_vertebrates", "fish_vertebrates", "unknown_vertebrates")] <- "vertebrates"
          diet5[diet5 %in% c("aquatic_invertebrates", "terrestrial_invertebrates", "unknown_invertebrates")] <- "invertebrates"
          diet5[diet5 %in% c("grass", "forbs_herbs_legumes", "roots_tubers_bulbs", "leaves_branches_bark",
                             "flowers_nectar_pollen", "unknown_vegetative_plant_material")] <- "vegetation"
          diet5[diet5 %in% c("fruit")] <- "fruit"
          diet5[diet5 %in% c("seeds")] <- "seeds"
          diet5[diet5 %in% c("fungi_lichens_bryophytes_algae", "unknown_animal_material",
                             "unknown_reproductive_plant_material", "unknown_plant_material", "other")] <- "other"
          diet5 <- factor(diet5, levels = c("vertebrates", "invertebrates", "vegetation", "fruit", "seeds", "other"))
        })
          print(str(dataBase))#1356 obs
          print(summary(dataBase))#1356 obs

        ##############################################################################
        # split database and process measures of interaction strength
        ##############################################################################
        x1 <- split(dataBase, dataBase$code)
        
        uni_code<-as.vector(unique(dataBase$code))
        uni_code_order<-order(uni_code)

        x1 <- lapply(1:length(x1), function(a) {x1[[a]]$id <- as.factor(sprintf("A%04.0f", a)); return(x1[[a]])})
        x1 <- lapply(x1, function(a) subset(a, !is.na(code)))
        x1 <- lapply(x1, function(a) if(any(is.na(a$rV))) {a$rV <- a$rF / a$rIV; a} else {a})
        x1 <- lapply(x1, function(a) if(any(is.na(a$rF))) {a} else {subset(a, rF > 0)})
        x1 <- lapply(x1, function(a) if(any(is.na(a$rV))) {a} else {subset(a, rV > 0)})
        x1 <- lapply(x1, function(a) if(any(is.na(a$rIV))) {a} else {subset(a, rIV > 0)})
        x1 <- lapply(x1, function(a) if(any(is.na(a$rEDC))) {a} else {subset(a, rEDC > 0)})
        x1 <- lapply(x1, function(a) if(any(is.na(a$rEDEC))) {a} else {subset(a, rEDEC > 0)})
        x1 <- lapply(x1, function(a) a <- within(a, {
          rF <- rF/sum(rF, na.rm = TRUE)
          rV <- rV/sum(rV, na.rm = TRUE)
          rEDC <- rEDC/sum(rEDC, na.rm = TRUE)
          rEDEC <- rEDEC/sum(rEDEC, na.rm = TRUE)
        }))
        x1 <- x1[sapply(x1, function(a) nrow(a) > 0)]
        ##############################################################################
        # bind all datasets together again
        ##############################################################################
        x2 <- do.call(rbind, x1)
        str(x2)
        fac_cod<-factor(x2$code)
        str(fac_cod)
        x2$cf1 <- corrFact$cf1[match(x2$diet3, corrFact$diet3)]
        x2$cf2 <- corrFact$cf2[match(x2$diet3, corrFact$diet3)]
        x2$uid <- 1:nrow(x2)
        x2 <- droplevels(x2)
        x2$rFP <- x2$rF
        x2$rVP <- x2$rV
        trophicData <- subset(x2, select = c(uid, id, code,
                                             diet1, diet2, diet3, diet4, diet5,
                                             rF, rFP, rV, rVP, rIV, rEDC, rEDEC, cf1, cf2,
                                             kingdom, phylum, class, order, family, genus, species, subspecies,
                                             taxon, rank, original_description))
        ##############################################################################
        # check of data structures
        ##############################################################################
        str(metaData)
        str(trophicData)
        str(corrFact)

        ##############################################################################
        # save all data to file
        ##############################################################################
        save(list = c("trophicData", "corrFact", "metaData"),
             file = paste0(folder_working,"/trophicDatabase.RData"))
    } else {
        ##############################################################################
        # alternatively load existing data
        ##############################################################################
        load(file = paste0(folder_working,"/trophicDatabase.RData"))
    }


  ################################################################################
  #1.2.2 Imputation of relative volume
    ################################################################################
    # fisher's z-transformation
    ################################################################################
    fisherZ <- function(x, type) {
      switch(type,
             fisherZ = 0.5 * log((1 + x) / (1 - x)),
             inverse = (exp(2 * x) - 1) / (exp(2 * x) + 1)
      ) 
    }
    ################################################################################
    # function to bootstrap correlation coefficients
    # and calculate weighted mean of r for each bootstrap sample
    ################################################################################
    bootCor <- function(dat, method = "weighted", nboot) {
      cors <- vector()
      for (i in 1:nboot){
        newDat <- dat[sample(1:nrow(dat), replace = TRUE),]
        cors[i] <- with(newDat, {
          switch(method,
                 weighted = fisherZ(weighted.mean(x = fisherZ(x = r, "fisherZ"), w = w), "inverse"),
                 unweighted = fisherZ(mean(x = fisherZ(x = r, "fisherZ")), "inverse")
          )
        }
        )
      }
      return(cors)
    }
    
    ################################################################################
    # plot relationships between relative frequency and relative volume
    # of food items
    ################################################################################
    str(trophicData)   
    corData <- subset(trophicData, !(is.na(rF) | is.na(rV)))
    corData <- droplevels(corData)
    
    corData$code<-factor(corData$code)
    cors <- sapply(levels(corData$code), function(x) {
      temp <- subset(corData, code == x, select = c(rF, rV))
      with(temp, cor(log(rF), log(rV)))
    })
    
    str(cors)
    
    bootData <- data.frame(r = cors,
                           n = metaData@data$n_rF[match(names(cors),
                                                        as.character(metaData$code))])
    bootData$n<-as.numeric(as.character(bootData$n))
    
    
    bootData$w <- with(bootData, n - 3) #In Ops.factor(n, 3) : ?-? not meaningful for factors
    cors1 <- bootCor(dat = bootData, method = "unweighted", nboot = 1e4)
    m1 <- round(median(cors1, na.rm = TRUE), 2)
    ci1 <- round(quantile(cors1, c(0.025,0.975), na.rm = TRUE), 2)
    m1; ci1
    
    ################################################################################
    # Supplementary_Figure_3. bootstrapped correlation between relative frequency and
    # relative volume
    ################################################################################
    print(paste0("OUTPUT: Supplementary Figure 3 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 3. bootstrapped correlation coefficient\nbetween rV and rF~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    pdf(paste0(folder_working,"/Supplementary_Figure_3.pdf"), width = 0.394 * 8, height = 0.394 * 8)
    par(mfrow = c(1, 1), cex = 0.7)
    hist(cors1, col = "black", border = FALSE, main = "bootstrapped correlation coefficient\nbetween rV and rF")
    dev.off()

    
    ################################################################################
    # Supplementary_Figure_4. Empirical relationships between relative frequency and relative volume
    ################################################################################
    print(paste0("OUTPUT: Supplementary Figure 4 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 4. Empirical relationships between relative frequency and relative volume ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    pdf(paste0(folder_working,"/Supplementary_Figure_4.pdf"), width = 0.394 * 16, height = 0.394 * 11)
    lattice.options(default.theme = standard.theme(color = FALSE))
    plot1 <- xyplot(rV ~ rF|diet5, trophicData, strip = TRUE,
                    par.strip.text = list(cex = 0.7),
                    layout = c(3, 2),
                    xlab = "Relative frequency",
                    ylab = "Relative volume",
                    scales = list(log = 10, at = c(1e-5, 1e-3, 1e-1),
                                  labels = c(expression(10^-5, 10^-3, 10^-1))),
                    panel = function(x, y,...) {
                      panel.xyplot(x, y, cex = 0.7, pch = 1, col = "black", ...)
                      panel.loess(x, y, col = "blue", lwd = 3, ...)
                      #panel.lmline(x, y, col = "blue", lwd = 1, ...)
                      panel.text(log10(0.5), log10(3e-4), adj = 1, cex = 0.7,
                                 label = bquote(italic(r)==.(round(cor(x, y, use = "complete.obs"), 2))))
                      panel.abline(0, 1, col = "red")
                      panel.text(log10(2e-4), log10(1e-3), "1:1", font = 2,
                                 col = "red", cex = 0.7, adj = 0)
                    })
    plot1
    dev.off()
    ################################################################################
    ################################################################################
    # fit MCMCglmm to impute missing data for relative volume
    ################################################################################
    trophicDataSub <- subset(trophicData, !(is.na(rF)), select = c(rF, rV, rFP, rVP, diet5, code, uid))
    trophicDataSub <- droplevels(trophicDataSub)
    prior1 = list(R = list(V = diag(1), nu = 0.002),
                  G = list(G1 = list(V = diag(1), nu = 0.002)))
    m0 <- MCMCglmm(log(rV) ~ log(rF) * diet5,
                   random = ~ code,
                   rcov = ~ units,
                   family = c("gaussian"),
                   nitt = 130000, thin = 100, burnin = 30000,
                   prior = prior1, data = trophicDataSub, pr = TRUE, pl = TRUE)

    ################################################################################
    # predict relative volume for studies with missing data based on model
    ################################################################################
    trophicDataSub$rVP <- colMeans(exp(m0$Liab))
    trophicDataSubList <- lapply(unique(trophicDataSub$code), function(x) {
      temp <- subset(trophicDataSub, code == x)
      temp <- within(temp, {
        rVP <- rVP/sum(rVP, na.rm = TRUE)
      })
      return(temp)
    })
    trophicDataSub <- do.call(rbind, trophicDataSubList)
    
    ################################################################################
    # Supplementary_Figure_5. observed and imputed data in one plot
    ################################################################################
    print(paste0("OUTPUT: Supplementary Figure 5 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Figure 5. Relationships between relative volume with relative frequency of observed and imputed data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    pdf(paste0(folder_working,"/Supplementary_Figure_5.pdf"), width = 0.394 * 10, height = 0.394 * 10)
    plot(rVP ~ rF, trophicDataSub, log = "xy", pch = 19,
         col = ifelse(is.na(trophicDataSub$rV), "red", "black"),
         main = "Observed (black) and imputed (red)\nrV data")
    abline(0, 1)
    dev.off()
 
    ################################################################################
    # calculate EDEC (estimated dietary energy content) based on relative Volume
    ################################################################################
    xtabs( ~ is.na(rV) + is.na(rVP), trophicData)
    trophicData$rVP[match(trophicDataSub$uid, trophicData$uid)] <- trophicDataSub$rVP
    trophicDataList <- lapply(unique(trophicData$code), function(x) {
      temp <- subset(trophicData, code == x)
      temp <- within(temp, {
        if(!any(is.na(rVP))) {rEDC <- rVP * cf1}
        rEDC <- rEDC/sum(rEDC, na.rm = TRUE)
        if(!any(is.na(rEDC))) {rEDEC <- rEDC * cf2}
        rEDEC <- rEDEC/sum(rEDEC, na.rm = TRUE)
      })
      return(temp)
    })
    trophicData <- do.call(rbind, trophicDataList)
    #trophicData <- subset(trophicData, rEDEC != 0)
    str(trophicData)
    save(list = c("trophicData", "corrFact", "metaData"),
         file = paste0(folder_working,"/trophicDatabaseImputed.RData"))

    trophicData$id_rows_num <- as.numeric(as.character(rownames(trophicData)))
    dataBase_original <- read.csv(paste0(folder_working,"/Data_all_studies_in_a_table_csv6.csv")) #
    database_original_with_energy_imputed<-merge(dataBase_original, trophicData, by.x="Item_ID", by.y = "id_rows_num",all.x = T)
    save(database_original_with_energy_imputed,  file = paste0(folder_working,"/database_original_with_energy_imputed.RData"))

    #We save selecting some fields and changing some names (Supplementary Table 5):
    database_original_with_energy_imputed2<-database_original_with_energy_imputed[,c("uid","id","Item_ID","code.x","original_description.x","taxon.x","rank.x","diet1.x","diet2.x","diet3.x","diet4","diet5","kingdom.x","phylum.x","class.x","order.x","family.x","genus.x","species.x","subspecies.x","rF.x","rV.x","rIV.x","rEDC.x","rEDEC.x","rF.y","rFP","rV.y","rVP","rIV.y","rEDC.y","rEDEC.y","cf1","cf2")]
    colnames(database_original_with_energy_imputed2)<-c("uid","id","Item_ID","code","original_description","taxon","rank","diet1","diet2","diet3","diet4","diet5","kingdom","phylum","class","order","family","genus","species","subspecies","rF_original","rV_original","rIV_original","rEDC_original","rEDEC_original","rF_imputed","rFP_imputed","rV_imputed","rVP_imputed","rIV_imputed","rEDC_imputed","rEDEC_imputed","cf1","cf2")
    write.csv(file=paste0(folder_working,"/database_original_with_energy_imputed2.csv"), x=database_original_with_energy_imputed2, row.names = T)
    save(database_original_with_energy_imputed2, file=paste0(folder_working,"/database_original_with_energy_imputed2.RData"))

    print(paste0("OUTPUT: Supplementary Table 5 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 5. Database of diet of brown bear ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    database_original_with_energy_imputed2


    #1.4 Manual modification of Supplementary Table 5 to create Supplementary Table 7
      # We select from Supplementary Table 5 ("database_original_with_energy_imputed2.csv")the items of the studies 
      #included in our study and we reassign its described rEDEC to improve the description at the species level: 
        #For each group containing several species, we assigned its described rEDECi to the species (s) 
        #that was most frequent in the bear diet, rEDECs = rEDECi, when explicitly stated in the paper. 
        #If several species belonging to the group were consumed by bears, but none explicitly 
        #mentioned as the most frequently eaten or occurring in scat, we divided equally the rEDECi 
        #among the species consumed by bears and present in the area, rEDECs = rEDECi,/nspecies. 
      #Supplementary Table 7: "database_original_with_energy_imputed2_ALL_corrected_names2.csv", we read
    #the manually modified file and export to the output folder:
    read.csv(file=paste0(folder_working,"/database_original_with_energy_imputed2_ALL_corrected_names2.csv"))
    write.csv(file=paste0(folder_working,"/database_original_with_energy_imputed2_ALL_corrected_names2.csv"), x=database_original_with_energy_imputed2_ALL_corrected_names2, row.names = T)
    print(paste0("OUTPUT: Supplementary Table 7 #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    print(paste0("Supplementary Table 7. Database of diet of brown bear ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
    database_original_with_energy_imputed2_ALL_corrected_names2

      
      