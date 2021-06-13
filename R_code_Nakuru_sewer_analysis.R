
library("FactoMineR")
library("factoextra")
library("corrplot")
library(raster)
library(sp)
library(caTools)
setwd("E:/Masters_Program/Semester 2/Spatial analysis/Class_project_GIS/Refined_Data/Final_dataset/R_data_inputs_transformed_factors")

ebu=raster('Final_transformed_Euclid_Built_up_factor.tif')
pd=raster('final_transformed_GHS_Population_factor.tif')
elv=raster('final_transformed_Elavation_factor.tif')
er=raster('Final_transformed_Euclid_River_factor.tif')
euc=raster('Final_transformed_Euclid_Urban_centers_factor.tif')
nl=raster('final_transformed_night_light_factor.tif')
rj=raster('final_transformed_road_junc_intensity_factor.tif')
slp=raster('final_transformed_slope_factor.tif')

ras_files=stack(ebu,pd,elv,er,euc,nl,rj,slp)
rast_df=as.data.frame(ras_files)
cleaned_df=na.omit(rast_df)

set.seed(101)
n=nrow(cleaned_df)
train_index=sample(1:n,size = round(0.8*n),replace=FALSE)
train=cleaned_df[train_index,]
test=train=cleaned_df[-train_index,]


pca_ana=PCA(train, scale.unit = FALSE, ncp = 5, graph = FALSE)

##to get amount of variance or information in each principle components or latent variables.
##eig.val <- get_eigenvalue(pca_ana)

##his function provides a list of matrices containing all the results for the active variables
var <- get_pca_var(pca_ana)

##Contributions of variables to PCs
head(var$contrib, 8)

##plot to highlight the most contributing variables for each dimension:
corrplot(var$contrib, is.corr=FALSE)

fviz_contrib(pca_ana, choice = "var", axes = 1:5, top = 10)
pdf("PCA_contr.pdf")
print(var$contrib.plot)
dev.off()

# Plot of variables
var.plot <- fviz_pca_var(pca_ana)

pdf("PCA.pdf") # Create a new pdf device
print(var.plot)
dev.off() # Close the pdf device
