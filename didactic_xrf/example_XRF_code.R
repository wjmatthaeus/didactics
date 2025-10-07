#example EDA, PCA, etc. code using several postgraduate datasets
#W.J. Matthaeus 7/10/2026
#this includes some data quality control on the way in
#and interpretation

#this is to load in required packages
pkgs <- c("readr", "ggplot2", "dplyr", "readxl",
          "tibble", "tidyr", "stringr",
          "PerformanceAnalytics", "psych",
          "FactoMineR", "factoextra","rstatix")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE)

#change the working directory. update the string "/Users.." to contain the 
#path to the directory containing this script and your data file
setwd("/Users/willmatthaeus/Dropbox/TCD Postdoc/didactic_xrf/")

#only change this here, then the rest of the code should work for *your* data
input_switch <-"Ant"
# input_switch <- "Sid"

if(input_switch=="Sid"){
#read in the data
xrf <- read_excel(path = "Mastersheet.xlsx")
#find the group ids from the sample name and make a new column
xrf$group<-NA
xrf[grep("CSD",xrf$SAMPLE),]$group <- "new"
xrf[is.na(xrf$group),]$group <- "old"
xrf$group <- as.factor(xrf$group)

#separate out the data columns and the error columns
elements<-xrf[seq(2,76,2)]
element_names <- colnames(elements)
errors<-xrf[seq(3,77,2)]
error_names <- colnames(errors)
}

if(input_switch=="Ant"){
  
  #read in the data
  xrf <- read_excel(path = "STANAST_XRF_GINK.xlsx")
  #find the group ids from the sample name and make a new column
  which(xrf$SAMPLE=="b4-51258-F")
  xrf[57,]$SAMPLE<-"b4-51258-F_2"
  which(xrf$SAMPLE=="b4-51408-F")
  xrf[77,]$SAMPLE<-"b4-51408-F_2"
  
  #grouping variables
  xrf$locality<-NA
  xrf[grep("AST",xrf$SAMPLE),]$locality <- "Astartekløft"
  xrf[is.na(xrf$locality),]$locality <- "South Tankrediakløft"
  xrf$locality <- as.factor(xrf$locality)
  which(is.na(xrf$locality))
  xrf$locality <- factor(xrf$locality)
  
  xrf$taxon <- factor(xrf$TAXON)
  
  #separate out the data columns and the error columns
  elements<-xrf[3:32]
  element_names <- colnames(elements)
  
}

##Boxplots

if(input_switch=="Ant"){
  #adjust the 'shape' of the data from a matrix of elements to a table where 
  #each row is a single element concentration with a other categorical variables
  #like taxon and locality
  xrf_long <- xrf%>%
    pivot_longer(cols = all_of(element_names),names_to = c("element"))%>%
    rename(ppm=value)
  
  #if you want to filter rows, keep certain ones or drop others
  # xrf_long <- filter(locality=="Astartekløft", 
  #                    taxon!=c("GINKGOITES MINUTA", "GINKGOITES")) 
  
  a_big_boxplot<-xrf_long %>%  
    ggplot()+#just to get things started
    geom_boxplot(aes(x=taxon, y=ppm, color=locality))+#make boxplot shapes, separate in space using taxon, and color using locality
    facet_wrap(element~., scales="free")+#separate element plots out into separate panels
    theme(axis.text.x=element_text(angle=45, hjust = 1))
  
  #some plot saving code
  #useful for making nice plots with high res at a particular size  
  #i just played with the size til it looked good
  #it saves into whatever directory you set at the the top
  plotName<-"Anto_element_boxplots.png"
  ggsave(plotName, plot = a_big_boxplot, device = "png", path = ".",
         scale = 1, height = 16, width = 14, units = c("in"),
         dpi = 300, limitsize = TRUE)
  
}

##remove columns with zero variance, which causes errors for subsequent analyses
#find zero variance columns
var_nz <- function(x) !is.na(var(x[x != 0]))
varying_elements_map<-apply(elements,2,var_nz)
#update data coulmns and names to only those with varianc
elements <- elements[,varying_elements_map]
element_names <- colnames(elements)


#individual tests of normality
elements %>% ungroup %>%  
  shapiro_test(element_names) %>% 
  arrange(variable)
#guide for interprting this and the next one:
#https://www.datanovia.com/en/lessons/normality-test-in-r/#shapiro-wilks-normality-test

#test of multivariate normality
elements %>% mshapiro_test
#the data are non-normal individually and as a multivariate group
#this is probably ok on it's own for PCA

#this is the important one!
chart.Correlation(elements)
ele_cor <- cor(elements) %>% abs 
# ele_cor_sig <- apply(ele_cor, c(1,2), function(x){any(abs(x)>0.3)})

#correlation matrix with significance and correlation coefficients in top 
#right triangle, histograms on diagonal, and biplots with best fit lines in red
#this is also ugly, but save it in a large format, zoom in and look at the 
#relationships between variables, are they linear?

#couple of different options for removing some columns
#things that don't have correlations about 0.3 with anything
#this is the more important thing for PCA

if(input_switch=="Sid"){
low_cor <- c("Ba","Sb","Rb")
#elements with normal-ish distributions based on the histograms
normal <- c("Bal","Mo","Nb","Cr","V","Ti","Ca","K")
}

if(input_switch=="Ant"){
  low_cor <- c("Ag","Mo")
  #elements with normal-ish distributions based on the histograms
  normal <- c("S","Cl","Si","P","K","Ca","Cr")
}

#normal elements
ele_norm <- elements%>%select(normal)#select(!low_cor) 
ele_norm_names <- colnames(ele_norm)
#individual tests of normality ... again
ele_norm %>% ungroup %>%  
  shapiro_test(ele_norm_names) %>% 
  arrange(variable)


#test of multivariate normality ... again
ele_norm %>% mshapiro_test
#still not normal

#to understand the linear relationship assumption let's break the correlation 
#matrix down into the normal and non-normal elements
chart.Correlation(ele_norm)
ele_notnorm <- elements%>%select(!normal)#select(!low_cor) 
chart.Correlation(ele_notnorm)

#filtering out the few elements with very low correlations
#this is just going to clarify things in the PCA
ele_cor <- elements%>%select(!low_cor) 
ele_cor_names <- colnames(ele_cor)
corr_mat_cor <- ele_cor %>% cor

#Conclusion: the outputs of this PCA are suspect because the underlying 
#relationships may not be linear, particularly for the non-normal elements
#you can still use PCA as an exploratory technique, but refer back to the
#correlation matrix when making element-level interpretations


ele_cor<-data.frame(ele_cor)
rownames(ele_cor) <- xrf$SAMPLE
# https://www.sthda.com/english/wiki/wiki.php?id_contents=7851
res.pca <- PCA(ele_cor, graph = FALSE)
eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
#the variance exlained in pc1 and 2 combined is less than 25%,
#add pc3 you get up to 35%
res.pca$var$contrib
#pc1 is mostly Cr, Mo, Ca, Nb, V; pc2 is Bal, Zr, Cu, Nb

fviz_pca_var(res.pca)
#the first two dimensions of the PCA don't single out any elements
#as driving variation
fviz_pca_var(res.pca, axes = c(1,3))
#the third two dimension (vertical there) 
#collapses a bit

fviz_pca_ind(res.pca, label="none", habillage = xrf$locality)
#no separation of old and new

fviz_pca_ind(res.pca, axes = c(1,3), label="none", habillage = xrf$locality)
#no separation of old and new
#but remember, the PCA (1) didn't work very well and (2) is not designed
#to test for differences between groups
