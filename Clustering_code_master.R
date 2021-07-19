library('cluster')
library("FactoMineR")
library("factoextra")
library("corrplot")
library(raster)
library(sp)
library(sf)
library(caTools)
library(factoextra)
library(stats)
library(NbClust)

##setwd("E:/Masters_Program/Semester 2/Spatial analysis/Class_project_GIS/Refined_Data/Final_dataset/R_data_inputs_transformed_factors/Cluster_Output_Junk")
rst=raster("E:/Masters_Program/Semester 2/Spatial analysis/Class_project_GIS/Refined_Data/Final_dataset/R_data_inputs_transformed_factors/PCA_Raster_outputs/PC123_merged.tif")

##rstDF <- values(rst)

rstDF <- as.data.frame(rst)

# Check NA's in the data
idx <- complete.cases(rstDF)

# Initiate the raster datasets that will hold all clustering solutions 
# from 2 groups/clusters up to 12
rstKM <- raster(rst[[1]])
rstCLARA <- raster(rst[[1]])


for(nClust in 2:30){
  
  cat("-> Clustering data for nClust =",nClust,"......")
  
  # Perform K-means clustering
  km <- kmeans(rstDF[idx,], centers = nClust, iter.max = 10)
  
  # Perform CLARA's clustering (using manhattan distance)
  cla <- clara(rstDF[idx, ], k = nClust, metric = "manhattan")
  
  # Create a temporary integer vector for holding cluster numbers
  kmClust <- vector(mode = "integer", length = ncell(rst))
  claClust <- vector(mode = "integer", length = ncell(rst))
  
  # Generate the temporary clustering vector for K-means (keeps track of NA's)
  kmClust[!idx] <- NA
  kmClust[idx] <- km$cluster
  
  # Generate the temporary clustering vector for CLARA (keeps track of NA's too ;-)
  claClust[!idx] <- NA
  claClust[idx] <- cla$clustering
  
  # Create a temporary raster for holding the new clustering solution
  # K-means
  tmpRstKM <- raster(rst[[1]])
  # CLARA
  tmpRstCLARA <- raster(rst[[1]])
  
  # Set raster values with the cluster vector
  # K-means
  values(tmpRstKM) <- kmClust
  # CLARA
  values(tmpRstCLARA) <- claClust
  
  # Stack the temporary rasters onto the final ones
  if(nClust==2){
    rstKM    <- tmpRstKM
    rstCLARA <- tmpRstCLARA
  }else{
    rstKM    <- stack(rstKM, tmpRstKM)
    rstCLARA <- stack(rstCLARA, tmpRstCLARA)
  }
  
  cat(" done!\n\n")
}

# Write the clustering solutions for each algorithm
##writeRaster(rstKM,"KMeans2_30.tif", overwrite=TRUE)
##writeRaster(rstCLARA,"CLARA2_30.tif", overwrite=TRUE)




# Start a data frame that will store all silhouette values
# for k-means and CLARA   
clustPerfSI <- data.frame(nClust = 2:30, SI_KM = NA, SI_CLARA = NA)


for(i in 1:nlayers(rstKM)){ # Iterate through each layer
  
  cat("-> Evaluating clustering performance for nClust =",(2:30)[i],"......")
  
  # Extract random cell samples stratified by cluster
  cellIdx_RstKM <- sampleStratified(rstKM[[i]], size = 2000)
  cellIdx_rstCLARA <- sampleStratified(rstCLARA[[i]], size = 2000)
  
  # Get cell values from the Stratified Random Sample from the raster 
  # data frame object (rstDF)
  rstDFStRS_KM <- rstDF[cellIdx_RstKM[,1], ]
  rstDFStRS_CLARA <- rstDF[cellIdx_rstCLARA[,1], ]
  
  # Make sure all columns are numeric (intCriteria function is picky on this)
  rstDFStRS_KM[] <- sapply(rstDFStRS_KM, as.numeric)
  rstDFStRS_CLARA[] <- sapply(rstDFStRS_CLARA, as.numeric)
  
  # Compute the sample-based Silhouette index for: 
  #    
  # K-means
 
  clCritKM <- intCriteria(traj = as.matrix(rstDFStRS_KM), 
                          part = as.integer(cellIdx_RstKM[,2]), 
                          crit = "Silhouette")
  # and CLARA
  clCritCLARA <- intCriteria(traj =  as.matrix(rstDFStRS_CLARA), 
                             part = as.integer(cellIdx_rstCLARA[,2]), 
                             crit = "Silhouette")
  
  # Write the silhouette index value to clustPerfSI data frame holding 
  # all results
  clustPerfSI[i, "SI_KM"]    <- clCritKM[[1]][1]
  clustPerfSI[i, "SI_CLARA"] <- clCritCLARA[[1]][1]
  
  cat(" done!\n\n")
  
}

write.csv(clustPerfSI, file = "clustPerfSI_combined_2_30.csv", row.names = FALSE)


