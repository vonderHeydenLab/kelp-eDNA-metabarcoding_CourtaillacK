##################################################################################
############  Refining eDNA metabarcoding in a kelp forest context  ##############
##################################################################################

#Kira-Lee Courtaillac, May 2024
#kiralee1506@gmail.com

#delete history 
rm(list =ls())

# set working directory (with data files) #ctrl+shift+H
setwd()

# Load required libraries
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(tidyverse)
library(labdsv)
library(readr)
library(ggordiplots)

# Load eDNA community matrix dataset (OTUs as columns, samples as rows) - OTU_data with ID column for sample names
# Load metadata file (variables = Site, Day, Transect, Depth, Repeat) - here called sample_metadata

sample_metadata <- read.csv("SampleInfoStandard.csv") #want to have simple categorical variables so can sync OTUs to metadata efficiently
OTUs_standard <- read_csv("OTUs_standard.csv")
#set row names
OTU_matrix <- OTUs_standard[, -1]
row.names(OTU_matrix) <- OTUs_standard$sample_id
rm(OTUs_standard)

# Set sample IDs as the common row name between the 2 documents
row.names(sample_metadata) <- sample_metadata$sample_id



############ Confirm sequencing depth sufficient to detect a-diversity in each sample (rarecurve) #####################

#rarefaction curves account for when have differing diversity among your samples - comparing communities
#rarefied is accounting for all the very low abundance OTUs (the rare things), which show up when pick up more species - this will help us make the curve with a good predictive power to check how much of the community is represented in each sample
#looking at the lowest number of samples you can include that will still be informative enough
raremax <- min(rowSums(OTU_matrix)) 
#returns the minimum reads of any one sequence within the matrix - if curves approach an asymptote here, all samples have likely achieved sufficient sequencing depth i.e.,  additional sequencing reads are unlikely to reveal many new species
#plot curve (ensure no rows / columns with all zero values):
rarecurve(OTU_matrix, step = 20, sample = raremax, col = "blue", cex = 0.6) # Takes some time
#all samples should plateau



############ DEVELOPING AN eDNA INDEX #########################
# Step 1: Calculate the relative abundance of each ASV within each sample - relative abundance of an OTU within its respective sample
relative_abundance <- OTU_matrix / rowSums(OTU_matrix)
# Step 2: Normalize the relative abundance of each taxa in each sample
max_abundance <- apply(OTU_matrix, 2, max) # maximum value for each column (i.e., each OTU) across all rows (samples).
eDNA_index <- relative_abundance / max_abundance
# Now eDNA_index contains the eDNA index values for each taxa in each sample

# Transform eDNA index into Qualitative data for Jaccard treatment: 
binary_data <- ifelse(eDNA_index > 0, 1, 0)



############ SPATIO-TEMPORAL COMMUNITY COMPARISONS ############
#SITE COMPARISONS: ----

#SEMI-QUANTITATIVE:
#1. create NMDS - use eDNA index and stipulate "bray" for method
nmdsBray<-metaMDS(eDNA_index, distance = "bray", trymax = 100) #stress <0.2 - okay to interpret the graphic & best solution found in 21 runs of real data
ordiplot_Site <- gg_ordiplot(nmdsBray, groups = sample_metadata$Site, pt.size=1)
Site_NMDS <- ordiplot_Site$plot #extract plot for editing
Site_NMDS + theme_classic() + labs(title = "NMDS Site") #use ggplot functions to edit to your preference
#2. PERMANOVA
#stipulate correct formular according to predictors - may need to use strat argument too (see https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/)
test1<-adonis2 (eDNA_index~Site, data = sample_metadata, permutations=9999, method = "bray") #modify model used according to your variables 
test1
    #F value is your test stat - larger means more likely for significance 
    #R2 shows that ~x% of variation of commun structure explained by Site
#3. PERMDISP (see https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permdisp/)
Braydistances<-vegdist(eDNA_index, method = "bray") #vegdist function creates bray-curtis distances by default btw sites
test2<-anova(betadisper(Braydistances,sample_metadata$Site)) #put distances into an anova with a betadispersal function - are there diffs in the dispersion of the data within these distance values we just created
test2
    #non-significant (p > 0.5), indicates a homogeneous dispersion of samples
#4. Indicator species analysis:
ind_val <- indval(eDNA_index, sample_metadata$Site)
summary(ind_val)
#OTUs with indicator_values close to 1 are strong indicators of those clusters and likely to be associated with that particular site (especially for low p-values)
#also investigate 'indicspecies’ v.1.7.15 package (De Cáceres, Jansen & Dell, 2024) - https://cran.r-universe.dev/indicspecies/doc/manual.html


#QUALITATIVE:
#1. create NMDS (change method - JACCARD):
nmdsJaccard<-metaMDS(binary_data, distance = "jaccard", trymax = 30)
ordiplot_Site_Jacc <- gg_ordiplot(nmdsJaccard, groups = sample_metadata$Site, pt.size=1)
Site_NMDS_Jacc <- ordiplot_Site_Jacc$plot
Site_NMDS_Jacc + theme_classic() + labs(title = "NMDS Site Jacc")
#2. PERMANOVA 
test3<-adonis2 (binary_data~Site, data = sample_metadata, permutations=9999, method = "jaccard") #change method
test3
#3. PERMDISP
Jaccdistances<-vegdist(binary_data, method = "jaccard") 
test4<-anova(betadisper(Jaccdistances,sample_metadata$Site)) 
test4
#4. Indicator species analysis:
ind_val <- indval(binary_data, sample_metadata$Site)
summary(ind_val)



#Analysis per site - Analysing individual site-level data in this way eliminates the portion of variance due to between-site differences, effectively amplifying the contributions of the remaining hierarchical sampling levels, including transect, depth and day

#make subsets:
#eDNA index:
remove_DH_rows <- 25:48
OTUs.MP <- eDNA_index[-remove_DH_rows, ]
Samples_MP <- sample_metadata[-remove_DH_rows, ]
row.names(OTUs.MP) <- Samples_MP$sample_id

remove_MP_rows <- 1:24
OTUs.DH <- eDNA_index[-remove_MP_rows, ]
Samples_DH <- sample_metadata[-remove_MP_rows, ]
row.names(OTUs.DH) <- Samples_DH$sample_id

#binary_data:
remove_DH_rows <- 25:48
OTUs.MP.binary <- binary_data[-remove_DH_rows, ]
remove_MP_rows <- 1:24
OTUs.DH.binary <- binary_data[-remove_MP_rows, ]


#SEMI-QUANTITATIVE PER SITE: ----

#Site 1, Miller's Point:
#nMDS:
    nmdsBrayMP<-metaMDS(OTUs.MP, distance = "bray", trymax = 100) #use subset and correct distance
    #Develop nMDS for each hierarchical level:
    #Day
    ordiplot_MPDays <- gg_ordiplot(nmdsBrayMP, groups = Samples_MP$Day, pt.size=1)
    MPDays <- ordiplot_MPDays$plot
    MPDays + theme_classic() + labs(title = "MP day 1 vs day 2") 
    #Transect
    ordiplot_MPTrans <- gg_ordiplot(nmdsBrayMP, groups = Samples_MP$Transect, pt.size=1)
    MPTrans <- ordiplot_MPTrans$plot
    MPTrans + theme_classic() + labs(title = "MP T1 vs T2") 
    #Depth
    ordiplot_MPDepth <- gg_ordiplot(nmdsBrayMP, groups = Samples_MP$Depth, pt.size=1)
    MPDepth <- ordiplot_MPDepth$plot
    MPDepth + theme_classic() + labs(title = "MP shallow vs deep")
#PERMANOVA:
    #use appropriate model for your variables, consider interactions and nesting, use correct method
    test5<-adonis2(OTUs.MP~Day+Transect+Depth+ Day:Transect + Day:Depth + Transect:Depth, data = Samples_MP, permutations=9999, method = "bray")
    test5
#PERMDISP:
    #perform for each variable
    #Day
    distances_MP<-vegdist(OTUs.MP, method = "bray") #changed method
    test6<-anova(betadisper(distances_MP,Samples_MP$Day)) 
    test6
    #Transect
    test7<-anova(betadisper(distances_MP,Samples_MP$Transect)) #change variable
    test7
    #Depth
    test8<-anova(betadisper(distances_MP,Samples_MP$Depth)) 
    test8
#Indicator species analysis:
    #perform for each variable, and decide indicator value cut-off and significance cut-off (e.g., >0.8 indicator value and <0.05 significance)
    #filter out absent OTUs for this site:
    absent_otus <- colSums(OTUs.MP) == 0
    filtered_MP <- OTUs.MP[, !absent_otus]#Day
    #Day
    ind_val <- indval(filtered_MP, Samples_MP$Day)
    summary(ind_val)
    #Transect
    ind_val <- indval(filtered_MP, Samples_MP$Transect)
    summary(ind_val)
    #Depth
    ind_val <- indval(filtered_MP, Samples_MP$Depth)
    summary(ind_val)

#REPEAT FOR EACH SITE
    


#QUALITATIVE PER SITE: ----

#Site 1, Miller's Point:
#nMDS:
    nmdsJaccMP<-metaMDS(OTUs.MP.binary, distance = "jaccard", trymax = 100) 
    #Day
    ordiplot_MPDays_Jacc <- gg_ordiplot(nmdsJaccMP, groups = Samples_MP$Day, pt.size=1)
    MPDays <- ordiplot_MPDays_Jacc$plot
    MPDays + theme_classic() + labs(title = "MP day 1 vs day 2 - Jaccard") 
    #Transect
    ordiplot_MPTrans_Jacc <- gg_ordiplot(nmdsJaccMP, groups = Samples_MP$Transect, pt.size=1)
    MPTrans <- ordiplot_MPTrans_Jacc$plot
    MPTrans + theme_classic() + labs(title = "MP T1 vs T2 - Jaccard") 
    #Depth
    ordiplot_MPDepth_Jacc <- gg_ordiplot(nmdsJaccMP, groups = Samples_MP$Depth, pt.size=1)
    MPDepth <- ordiplot_MPDepth_Jacc$plot
    MPDepth + theme_classic() + labs(title = "MP shallow vs deep - Jaccard") 
#PERMANOVA:
    #use appropriate model for your variables, consider interactions and nesting, use correct method
    test1<-adonis2(OTUs.MP.binary~Day+Transect+Depth + Day:Transect + Day:Depth + Transect:Depth, data = Samples_MP, permutations=9999, method = "jaccard")
    test1
#PERMDISP:
    distances_MP<-vegdist(OTUs.MP.binary, method = "jaccard") #changed method
    #Day
    test2<-anova(betadisper(distances_MP,Samples_MP$Day)) 
    test2
    #Transect
    test3<-anova(betadisper(distances_MP,Samples_MP$Transect)) 
    test3
    #Depth
    test4<-anova(betadisper(distances_MP,Samples_MP$Depth)) 
    test4
#Indicator species analysis (use correct binary dataset):
    absent_otus <- colSums(OTUs.MP.binary) == 0
    filtered_MP <- OTUs.MP.binary[, !absent_otus]
    #Day
    ind_val <- indval(filtered_MP, Samples_MP$Day)
    summary(ind_val)
    #Transect
    ind_val <- indval(filtered_MP, Samples_MP$Transect)
    summary(ind_val)
    #Depth
    ind_val <- indval(filtered_MP, Samples_MP$Depth)
    summary(ind_val)
    
#REPEAT FOR EACH SITE


############ INVESTIGATING REPLICATION ########################

#Load data for samples contributing to accumulation curve (biological replicates):
library(readr)
OTUs_A1_A10 <- read_csv("OTUs_A1_A10.csv", 
                        col_types = cols(...1 = col_skip()))

#### accumulation function ----
accumulation_curve_rep <- specaccum(OTUs_A1_A10, permutations = 1000)
plot(accumulation_curve_rep, main = "Species Accumulation Curve A1-A10", xlab = "Number of Samples", ylab = "Number of OTUs") 

#Multi-model fitting:
#1. AIC weights function
akaike.weights <- function(x)
{
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  rel.LL <- exp(-0.5 * delta.aic)
  sum.LL <- sum(rel.LL, na.rm = TRUE)
  weights.aic <- rel.LL/sum.LL
  return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
}


#2. Fit models (these are the ones used by me - could fit even more)
#some of these models may not work for your data, you can just omit them and even fit another one (see ?fitspecaccum for options)
lomo <- fitspecaccum(accumulation_curve_rep, "lomolino")
coef(lomo)
aic_lomo=AIC(lomo) 
plot(lomo, add = TRUE, col=2, lwd=2) #red

mm <- fitspecaccum(accumulation_curve_rep, "michaelis-menten")
coef(mm)
aic_mm=AIC(mm) #35.652
plot(mm, add = TRUE, col=3, lwd=2) #green

gom <- fitspecaccum(accumulation_curve_rep, "gompertz")
coef(gom)
aic_gom=AIC(gom) #27.136
plot(gom, add = TRUE, col=4, lwd=2) #blue

asy <- fitspecaccum(accumulation_curve_rep, "asymp")
coef(asy)
aic_asy=AIC(asy) #18.957
plot(asy, add = TRUE, col=5, lwd=2) #purple

gis <- fitspecaccum(accumulation_curve_rep, "logis")
coef(gis)
aic_gis=AIC(gis) #32.318
plot(gis, add = TRUE, col=6, lwd=2) #pink

wei <- fitspecaccum(accumulation_curve_rep, "weibull")
coef(wei)
aic_wei=AIC(wei) #-26.646
plot(wei, add = TRUE, col=11, lwd=2)


#3. Develop matrix containing each model's AIC, weights acc to lowest AIC and asymptotes:
res=matrix(NA,nrow=6,ncol=3)

rownames(res)=c("lomolino","michaelis-menten","gompertz","asymp","logis","weibull")
colnames(res)=c("AIC","Asymptote","Weigth")

res[,"AIC"]=c(aic_lomo,aic_mm,aic_gom,aic_asy,aic_gis,aic_wei)

res[,"Weigth"]=akaike.weights(c(aic_lomo,aic_mm,aic_gom,aic_asy,aic_gis,aic_wei))$weights

res[,"Asymptote"]=c(coef(lomo)[[1]],coef(mm)[[1]],coef(gom)[[1]],coef(asy)[[1]],coef(gis)[[1]],coef(wei)[[1]])


#4. Calcul Asymptote - mean of each models asymptote weighted in relevance by their AIC values
asymptote<-weighted.mean(res[,"Asymptote"], res[,"Weigth"])
# = average theoretical max no. of OTUs expected to be detected in your system
#By calculating the weighted mean of the asymptote values from different models, you're essentially trying to estimate a consensus asymptote value that best represents your data.


#5. Extrapolate using specslope and predict:
#Make a dataframe for the predictions - set upper limit here according to OTU richness
sat_dat<-data.frame(samps = seq(0,150))

#Run the specslope function using the best model from above - replace "wei" with best fit model as defined above)
sat_dat$slope<-specslope(wei, at = sat_dat$samps)

#Generate a factor variable that states whether a given slope is within the tolerance bound of 0.05.
#change 0.2 value to define which % of expected richness you want the sampling effort to cover
#(here, slope <= 0.2 which indicates extrapolation covers 80% of theoretical max of richness)
sat_dat$tolerance_OK<-ifelse(sat_dat$slope <= 0.2, "OK", "Not OK")

#Plot the results
ggplot(sat_dat, aes(x = samps, y = slope, shape = tolerance_OK)) +
  geom_point(colour = "black") +
  scale_shape_manual(values = c("Not OK" = 1, "OK" = 16))

#Extract the minimum number of samples at which sampling is associated with provided slope (here, slope <= 0.2 which indicates extrapolation covers 80% of theoretical max)
min(sat_dat[sat_dat$tolerance_OK == "OK",]$samps, na.rm = TRUE)
#29 samples will cover 80% of extrapolated richness


#6. Plot
#accumulation function
new_data <- data.frame(sites = 1:80)  # Define new sampling efforts (sites = samples) - increase this number depending on value returned for above function
predicted_curve <- predict(wei, newdata = new_data) #use best fit model for this

#change limits for plot according to max number of samples and richness:
#i.e., 29 samples covers 80% of richness, so I will plot up to 80 samples on x-axis, and OTU richness is <100 so I will plot to 100 on y axis:
plot(accumulation_curve_rep, xlab = "Number of Sites", ylab = "Species Richness", xlim = c(0,80), ylim = c(0,100), col = "black", lty = 1, lwd = 1.5) #set limits according to new_data samples and calculated asymptote
#add asymptote
abline(h=asymptote, col = "#36454F", lty = 1, lwd = 1.2)
#add predicted curve
lines(new_data$sites, predicted_curve, col = "#7CFC00", lty = 3, lwd = 1)
# Add vertical lines for sample indicators
abline(v = "29", col = "grey", lty = 2) #29 samples to detect 80% 
#add legend if you want
legend("bottomright", legend = c("Observed", "Predicted"), col = c("black", "#7CFC00"), lty = 1)


#7. Specpool can give indication of species observed vs species expected with different variants 
specpool(OTUs_A1_A10)
#     Species     chao  chao.se    jack1 jack1.se    jack2     boot  boot.se n
# All      76 84.12698 5.418178 90.22222 7.069322 93.13889 83.23884 4.409378 9




