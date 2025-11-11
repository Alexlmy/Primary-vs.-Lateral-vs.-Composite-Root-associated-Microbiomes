################ Raw data Processing ################
##########$$$ Make manifest for Local Desk $$$#########
# Set the path to where your sequence files are stored
## Note:when making the absolte-filepath, should be the same with this path!!!!!
path <- "/Users/mel270/Downloads/CLC-ITS2024/raw_data/"  
# Create a dataframe of the file locations
manifest<- data.frame(list.files(path), stringsAsFactors = F)
# Extract the sample names
sample_id <- sapply(strsplit(manifest$list.files.path., "_"), "[", 1)
# Make the file paths
absolute_filepath <- paste("/u2/mel270/CLC-ITS2024/raw_data/", manifest$list.files.path., sep = "")       
# Make the read directions (when there's 264 rows, direction should be 132)
direction <- rep(c("forward","reverse"), 132)           
# Make a new data frame including the required columns
manifest_ITS <- cbind.data.frame(sample_id, absolute_filepath, direction)
## Note: when showing "Error in data.frame(..., check.names = FALSE) : 
## arguments imply differing number of rows: 128, 20
## Check with the following code, making sure number of rows are the same.
NROW(absolute_filepath)
NROW(sample_id)
NROW(direction)
# Rename the columns
names(manifest_ITS) <- c("sample-id", "absolute-filepath", "direction")
View (manifest_ITS)
# Output the file
write.csv(manifest_ITS, "/Users/mel270/Downloads/CLC-ITS2024/manifest_ITS.csv", row.names = FALSE,quote=FALSE)
##########$$$##########$$$##########$$$##########$$$##########$$$##########$$$##########$$$##########$$$##########$$$










######################################################################################################
---
  title: "CLC-ITS"
author: "Mengying Liu"
date: "Dec 28, 2023"
output: pdf_document

######## library required ########
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("phyloseq"))
# install.packages('devtools')
# install.packages("usethis")
# devtools::install_github("gauravsk/ranacapa")
# install.packages("permute")
# install.packages("lattice")
# install.packages("pheatmap")
library(phyloseq)
library(usethis)
library(devtools)
library(ranacapa)
library(ggplot2) 
library(dplyr)
library(tidyr)
library(stringr)
library(lattice)
library(permute)
library(vegan)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(nlme)
library(permute)
library(lattice)
library(vegan)
library(bestNormalize)
library(emmeans)
library(multcompView)
library(mvtnorm)
library(survival)
library(TH.data)
library(MASS)
library(multcomp)
library(ALDEx2)
library(carData)
library(car)
library(CoDaSeq)
library(MASS)
library(phyloseq)
library(truncnorm)
library(survival)
library(NADA)
library(zCompositions)
library(emmeans)

##########
########## Load and manipulate data ##########
setwd("/Users/mel270/Downloads/Bobbi_DNAsequencing/CLC_BH/CLC-ITS2024/CLC-ITS-AnalysisR")
#setwd("/Users/mel270/Desktop/CARPITS_OldPipeline")
#feature<-read.table("ITS_feature-table-.tsv",sep='\t',row.names = 1, header = TRUE) # remove the "#" before "OTU"
feature<-read.table("ITS_feature-table.tsv",sep='\t',row.names = 1, header = TRUE) # remove the "#" before "OTU"
#View(feature)
metadata<-read.table("CLC-ITS-metadata.tsv",sep='\t',row.names = 1,header = TRUE) # change "-" to "."
#View(metadata)
names (metadata)
metadata$Plot<-factor(metadata$Plot)
metadata$SoilP<-factor(metadata$SoilP)
metadata$Tret<-factor(metadata$Tret)
metadata$Stage<-factor(metadata$Stage)
metadata$SampleType <-factor(metadata$SampleType)
metadata$RootType<-factor(metadata$RootType)
metadata$SampleLabel<-factor(metadata$SampleLabel)
metadata$PlantDensity<-factor(metadata$PlantDensity)
metadata$Yield<-factor(metadata$Yield)


ASV_1<-feature[,intersect(rownames(metadata),colnames(feature))]
#View(ASV_1)
setdiff (rownames(metadata),colnames(feature)) # Check whether metadata and ASV matches. Results are elements in rownames(metadata), but not in colnames(ASV), should be 0.
identical(colnames(ASV_1),rownames(metadata)) # TRUE
identical(ncol(ASV_1),nrow(metadata)) # TRUE
Taxon<-read.table("taxonomy.tsv",sep='\t',header = TRUE,row.names = 1,comment.char = "") # 1 indicates that the data file has row names. If we don't specify, the result will not have row names.
#View(Taxon)
Taxonomy<-Taxon[intersect(rownames(Taxon),rownames(feature)),] # The intersection of two sets is the material that they have in common
#View(Taxonomy)
setdiff (rownames(Taxon),rownames(feature)) # Check whether Taxon and ASV matches. Result should be 0.
rownames(ASV_1)[1:(dim(ASV_1)[1])]<-paste("ASV",1:dim(ASV_1)[1],sep='') # ????? dim: retrieve or set the dimension of an object.####### !!!!!!!!!!!!! ##### Change "feature" to "ASV" $$$%%^^&&*())*(&^$#^%^#&%$*$^^%#$@#$!@#~#~@%!^$^)
rownames(Taxonomy)<-rownames(ASV_1)
Ta<-separate(Taxonomy,1,c(LETTERS[1:7]),";") # seperate: turns a single character column into multiple columns. Here seperate into 7 columns.
#View(Ta)  ## 8114 ASV in total
Ta2<-sapply(Ta[,1:7],function(x) gsub("[a-z]__","",as.character(x))) 
# View(Ta2)
# sapply: take list, vector or data frame as input and gives output in vector or matrix.
# gsub("x","y",A) function in R is used for replacement operations (e.g.: replace x to y in the factor A). 
# as.character attempts to coerce its argument to character type
Tb<-as.data.frame(Ta2) # rownames: "1,2..."; columnames:"A,B,C..."
#View(Tb)
rownames(Tb)<-rownames(Taxonomy) # change rownames from "1,2,3..." to "ASV1,ASV2..."
#View(Tb)
colnames(Tb)<-c('Kingdom','Phylum','Class','Order','Family','Genus','Species') #change columnames from "A,B,C..." to "Kingdom,Phylum,Order..."
#View(Tb)
# Tb$Phylum<-factor(Tb$Phylum)
# levels(Tb$Phylum)
Taxonomy_0<-lapply(Tb,function(x)gsub("unidentified","Unclassified",x))
#View(Taxonomy_0) Note: at Phylum level, results like: "D_1__Proteobacteria"   "D_1__Proteobacteria"   "D_1__Proteobacteria"   "D_1__Actinobacteria"...
Taxonomy_cm<-as.data.frame(lapply(Taxonomy_0,function(x)replace_na(x,"Unclassified"))) # Replace the "Unclassified" with NA.
#View(Taxonomy_cm) Note:rownames are "1,2,3..."
rownames(Taxonomy_cm)<-rownames(Tb) #rownames change from "1,2,3..." to "ASV1,ASV2..."
#View(Taxonomy_cm)
#write.csv(Taxonomy_cm,"Fungi-allASV.csv") ### Change ASV names in Excel
# Remove unclassified at phylum level
Taxonomy_t1<-Taxonomy_cm[!(Taxonomy_cm$Phylum %in% c("Unclassified")),] #Remove "Unclassified" at phylum level.
# ## Confirm unclassified phylum have been removed, both results should be 0.
grep("Unclassified",Taxonomy_t1$Phylum)
Taxonomy_1<-Taxonomy_t1[!(Taxonomy_t1$Class %in% c("Unclassified")),] #Remove "Unclassified" at phylum level.
grep("Unclassified",Taxonomy_1$Class)
# View(Taxonomy_1)  ### 4211(note:4468 ASVs Old Pipeline)
write.csv(Taxonomy_1,"Fungi-ASV-NoUnclassified.csv")  ### Change ASV name in Excel
## Remove these ASVs for unclassified phylum in ASV table
# Taxonomy_1<-read.csv("Fungi-allASV.csv",row.names = 1,header = TRUE)
Taxonomy_1<-read.csv("Fungi-ASV-NoUnclassified.csv",row.names = 1,header = TRUE)
## Keep all ASVs, including the unclassified phylum, cuz only 4211 ASVs classified in total 7838 ASVs.
##Taxonomy_1<-read.csv("/Users/mel270/Desktop/CARPITS/Fungi-allASV.csv",row.names = 1,header = TRUE)
ASV<-ASV_1[intersect(rownames(ASV_1),rownames(Taxonomy_1)),] # The elements of intersect(x,y) are those elements in x and in y
#list (rownames(ASV_1))
#list (rownames(Taxonomy_1))

## Format data as phyloseq objects
OTU<-otu_table(as.matrix(ASV),taxa_are_rows = TRUE)
TAX<-phyloseq::tax_table(as.matrix(Taxonomy_1))
samples<-sample_data(metadata)
combined_plusctrl<-phyloseq(OTU,TAX,samples)
## Remove all control samples
combined <- subset_samples(combined_plusctrl, Ctrl!= "Y")
View (combined) ## in total:2800 ASV in 120 samples, Rt47 low reads filtered)
#levels (combined@sam_data$P_OW)
sub01 <-subset_samples(combined,Stage=="6wk") 
View (sub01)
######################################################################









# ######################################################
# ########## Venn Graph showing shared ASVs   ##########
# ## https://rdrr.io/github/Russel88/MicEco/man/ps_venn.html
# # install.packages("remotes")
# # remotes::install_github("Russel88/MicEco",force = TRUE)
# # library ("MicEco")
# ps_venn(combined,
#         "SampleType",
#         quantities = list(type=c("percent","counts"), font = 2), 
#         fraction = 0,
#         weight = FALSE,
#         relative = TRUE,
#         plot = TRUE,fill = c("#FFCCD2","#7ED180","#1B98E04D"))
# 
# ps_venn(combined_plusctrl,
#         "Controls",
#         quantities = list(type=c("percent","counts"), font = 2), 
#         fraction = 0,
#         weight = FALSE,
#         relative = TRUE,
#         plot = TRUE,fill = c("floralwhite","gray"))
# 
# ## Subset sample types
# sb1 <-subset_samples(combined,SampleType=="Rhizo")
# sb1 <-subset_samples(combined,SampleType=="Root")
# View(sb1)
# ## Year
# ps_venn(sb1,
#         "Year",
#         quantities = list(type=c("percent","counts"), font = 2), 
#         fraction = 0,
#         weight = FALSE,
#         relative = TRUE,
#         plot = TRUE,fill = c("lightblue","lightblue4","#1B98E04D"))
# ## Stage
# sample_data(sb1)$Stage <- factor(sample_data(sb1)$Stage,      
#                                  levels = c("L", "F"))
# ps_venn(sb1,
#         "Stage",
#         quantities = list(type=c("percent","counts"), font = 2), 
#         fraction = 0,
#         weight = FALSE,
#         relative = TRUE,
#         plot = TRUE,fill = c("yellowgreen","Yellow"))
# ######################################################




######## Rarefaction curve ########
## Curve removed contrls
sb1 <-subset_samples(combined,Stage=="6wk") 
curve<-ggrare(sb1, step = 250, se = FALSE)
## Curve doesn't remove contrls: curve<-ggrare(combined_plusctrl, step = 250, se = FALSE)

## SampleType
p<-curve+labs(title="Rarefaction_ITS")
p<-p+geom_line(aes(col=SampleType))+
  scale_color_manual(labels = c("Rhizo","Root"),
                     values = c("orange","darkgreen"))
p
ggsave("CLC_Rarefaction_ITS_SampleType_6wk_clean.png",width = 5,height=2.5)




#########################################
###### Relative abundance #####
##Determine which samples with lower reads need removing
#View(ASV)
#sort(colSums(ASV))
##Set the threshold to 12000
#ID_Sum<-as.data.frame(colSums(ASV))
#View(ID_Sum) #rownames: sample ID; columname: colSums(ASV)
#View(colSums(ASV)) #rownames: sample ID; no columname
#ID_Sum0<-cbind(rownames(ID_Sum),ID_Sum$`colSums(ASV)`) 
#Take a sequence of vector, matrix or data-frame arguments and combine by columns or rows, respectively.
#View(ID_Sum0) ##rownames:1,2,3....; no columname
#colnames(ID_Sum0)<-c("ID","abundance")
#ID_Sum1<-as.data.frame(ID_Sum0)
#View(ID_Sum1)
#ID_T<-ID_Sum1[as.numeric(as.character(ID_Sum1$abundance))>12000,]
#View(ID_T)
##select the samples used for downstream analysis
#sub_metadata1<-metadata[grep("N",metadata$Controls),] # Select only samples (remove controls).
#View(sub_metadata1)
#sub_metadata<-sub_metadata1[intersect(rownames(sub_metadata1),ID_T$ID),] # Select samples passed threshold of 12000
#View(sub_metadata)
#sub_ASV<-ASV[,intersect(rownames(sub_metadata),colnames(ASV))] ## Generate the the target sequencing sample ASVs
#View(sub_ASV)

## Format data as phyloseq objects
#OTU<-otu_table(as.matrix(sub_ASV),taxa_are_rows = TRUE)
#TAX<-tax_table(as.matrix(Taxonomy_1))
#samples<-sample_data(sub_metadata)
#sub_combind<-phyloseq(OTU,TAX,samples)
#View(sub_combind) ## 
View (combined)
## Resample an OTU table such that all samples have the same library size.
Sub_Rarefy<-rarefy_even_depth(combined,sample.size=2500,rngseed=999) 
## 1 samples removedbecause they contained fewer reads than `sample.size`:CARP2019SoilITSRhizo29CARP2020ITSBulk24CARP2019Root07


#########################################################
## Relative abundance in each Year at Class Level ###
#########################################################
## Resample an OTU table such that all samples have the same library size.
# library(pals)
# pal.bands(polychrome,show.names=TRUE)
View(Sub_Rarefy)
sample_data(Sub_Rarefy)$Stage_SampleType<-paste(sample_data(Sub_Rarefy)$Stage,sample_data(Sub_Rarefy)$SampleType,sep="-")
sample_data(Sub_Rarefy)$STPO<-paste(sample_data(Sub_Rarefy)$Stage_SampleType,sample_data(Sub_Rarefy)$Tret,sep="-")
sample_data(Sub_Rarefy)[]<-lapply(sample_data(Sub_Rarefy),factor) #Use lapply to apply the factor() function 
Sub_Rarefy_Class<-tax_glom(Sub_Rarefy, taxrank = "Class") # Merges species that have the same taxonomy at a certain taxaonomic rank (e.g.: "Family" here)
#View(Sub_Rarefy_Class) #43x284
## Order Class from high to low
A<-names(sort(rowSums(otu_table(Sub_Rarefy_Class)),decreasing=T))
otu_table(Sub_Rarefy_Class)<-otu_table(Sub_Rarefy_Class)[intersect(A,rownames(otu_table(Sub_Rarefy_Class))),]
tax_table(Sub_Rarefy_Class)<-tax_table(Sub_Rarefy_Class)[intersect(A,rownames(tax_table(Sub_Rarefy_Class))),]
View(A)
View(otu_table(Sub_Rarefy_Class))
View(tax_table(Sub_Rarefy_Class)) ## Check the number : 40

## Replace Class with low reads with "others" !!!!!!!!!!!!!!!!!!!!! low reads number?? how to define? 
tax_table(Sub_Rarefy_Class)[,3][29:40]<-rep("Others",12) ###### [,number of level][exhibit taxon:total taxon]######
Sub_Rarefy_Class<-tax_glom(Sub_Rarefy_Class, taxrank = "Class") # combine all the Family belonging to "others"
sample_data(Sub_Rarefy_Class)[]<-lapply(sample_data(Sub_Rarefy_Class),factor) 

Sub_Rarefy_merge<-merge_samples(Sub_Rarefy_Class,"STPO")
sample_data(Sub_Rarefy_merge)$SampleType<-levels(sample_data(Sub_Rarefy_Class)$SampleType)[get_variable(Sub_Rarefy_merge,"SampleType")] #### Need change to Bulk Soil, Rhizosphere Soil
sample_data(Sub_Rarefy_merge)$Stage<-levels(sample_data(Sub_Rarefy_Class)$Stage)[get_variable(Sub_Rarefy_merge,"Stage")]
sample_data(Sub_Rarefy_merge)$Tret<-levels(sample_data(Sub_Rarefy_Class)$Tret)[get_variable(Sub_Rarefy_merge,"Tret")]
sample_data(Sub_Rarefy_merge)[]<-lapply(sample_data(Sub_Rarefy_merge),factor) #this step is very important
# colourCount = length(unique(tax_table(Sub_Rarefy_merge)[,3])) ## number of levels (e.g.: Phylum is 2, Class is 3)
# getPalette = pal.bands(polychrome,show.names=FALSE) #Determine colors used in a visualization.
# otu_table(Sub_Rarefy_merge)<-otu_table(Sub_Rarefy_merge)/rowSums(otu_table(Sub_Rarefy_merge)) #Mean tax for samples
getPalette <- c("#a6cee3","#E1E7E9","#1f78b4",
                "#ff7f00","#de77ae", "#FFCC00","#999999","#EEE8AA","#5F9EA0","#99FFCC",
                "#FF6666", "dodgerblue2", "#CAB2D6", "green4", "#FB9A99", "deeppink1", "#FDBF6F", 
                "orchid1", "maroon", "#6A3D9A", "blue1", "steelblue4","darkturquoise", "green1", 
                "yellow4", "yellow3", "darkorange4", "brown","Strong_Blue")
#change the name of level which will be showed in figure
#transform the raw counts of reads to proportions within a sample
otu_table(Sub_Rarefy_merge)<-otu_table(Sub_Rarefy_merge)/rowSums(otu_table(Sub_Rarefy_merge)) #Mean tax for samples
#should only run one time when plotting!
SS<-levels(sample_data(Sub_Rarefy_merge)$SampleType)
SS<-gsub("Rhizo","Rhizosphere",SS)
levels(sample_data(Sub_Rarefy_merge)$SampleType)<-SS

# # Reordering group factor levels (default is alphabet,F before L. The below would reorder to L before F )
# sample_data(Sub_Rarefy_merge)$Stage <- factor(sample_data(Sub_Rarefy_merge)$Stage,      
#                                               levels = c("4wk", "6wk"))
# BB<-levels(sample_data(Sub_Rarefy_merge)$Stage)
# BB<-gsub("6wk","Flowering",BB)
# BB<-gsub("4wk","Vegetative",BB)
# levels(sample_data(Sub_Rarefy_merge)$Stage)<-BB

# Reorder stack barplot with: fill = fct_reorder(Genus, Abundance)    
p<-plot_bar(Sub_Rarefy_merge,"Tret", fill = "Class")+geom_bar(aes(fill=Class), stat="identity", position="stack")+scale_fill_manual(values = getPalette)+scale_y_continuous(labels=sprintf("%1.0f",seq(0,100,25)))
p+facet_grid(Stage~SampleType)+labs(y="Relative Abundance (%)",x="Tret")+guides(fill = guide_legend(nrow = 15))

ggsave("CLC_ITS_Relative Abund.png",width = 5,height=2.5)
##############################################################################:)







################################################################################
######### Alpha_diversity plot ######### Need to use rarefied data!!!!!! ######
################################################################################################################################################################################
################ Alpha Diversity-Boxplot with dots + Statistic Analysis ################################################################
################################################################################

## All alpha diversity index
View (combined)
View(Sub_Rarefy) ##
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy,SampleType=="Rhizo")
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy,SampleType=="Root")

Sub_Rarefy_RR<-subset_samples(Sub_Rarefy_Bulk,Stage=="4wk")
Sub_Rarefy_RR<-subset_samples(Sub_Rarefy_Bulk,Stage=="6wk")
View(Sub_Rarefy_RR)

# Sub_Rarefy_RR <- Sub_Rarefy_Bulk

## Calculate all alpha diversity indices
bac_alpha<-microbiome::alpha(Sub_Rarefy_RR,index="all")
# write.csv(bac_alpha,"CLC_ITS_Alpha_Root_6wk.csv")
# (1) Alpha - Observed richness ###############
observed<-subset(bac_alpha,select="observed")
mydata<-data.frame(sample_data(Sub_Rarefy_RR),observed)
# View (mydata)

## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=observed, fill=RootType)) + coord_cartesian(ylim = c(0, 100))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Observed richness")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))

## Rhizo 200

# library(viridis)
# mydata %>%
#   ggplot(aes(x=Plevel, y=evenness_pielou, fill=OW)) + coord_cartesian(ylim = c(0, 1))+ 
#   geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6) + geom_jitter(position=position_jitter(0.2))+
#   theme_bw() +labs(x="",y="Pielou's evenness")+ theme(axis.text.x = element_text(angle = 90),
#                                                       text = element_text(size = 12))
## Fit lm model
model<-lme(observed~SoilP*RootType,random=~1|Rep,data=mydata)
anova (model)
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/observed_anova_Root2019L.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/observed_anova_Root2019F.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/observed_anova_Root2020L.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/observed_anova_Root2020F.csv")

# Post Hoc Test- pairwase comparison for all sample types
model<-lme(observed~Tret_POW,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"Tret_POW")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/observed_PostHoc_Root2019L.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/observed_PostHoc_Root2019F.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/observed_PostHoc_Root2020L.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/observed_Plevel_PostHoc_Root2020F.csv")


# (2) Alpha - Shannon diversity ###############
diversity_shannon<-subset(bac_alpha,select="diversity_shannon")
mydata<-data.frame(sample_data(Sub_Rarefy_RR),diversity_shannon)
## Reorder for samples

## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=diversity_shannon, fill=RootType)) + coord_cartesian(ylim = c(0, 6))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Shannon diversity")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))

##
model<-lme(diversity_shannon~SoilP*RootType,random=~1|Rep,data=mydata)
anova(model)
# qqnorm(residuals(model))
# plot(residuals(model))
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/shannon_anova_Root2019L.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/shannon_anova_Root2019F.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/shannon_anova_Root2020L.csv")
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/shannon_anova_Root2020F.csv")

# Post Hoc Test- pairwase comparison for all sample types
model<-lme(diversity_shannon~Tret,random=~1|Rep,data=mydata)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
em<-emmeans(model,pairwise ~ Tret,adjust = "tukey") ##Estimated marginal means
em
write.csv(pair,file="Root_Alpha_shannon_PostHoc_Root_all.csv")

model<-lme(diversity_shannon~Plevel,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"Plevel")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/shannon_PostHoc_Root2019L.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/shannon_PostHoc_Root2019F.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/shannon_PostHoc_Root2020L.csv")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/shannon_Plevel_PostHoc_Root2020F.csv")

View (bac_alpha)

# (3) Alpha - Evenness ###############
evenness_pielou<-subset(bac_alpha,select="evenness_pielou")
mydata<-data.frame(sample_data(Sub_Rarefy_RR),evenness_pielou)
## Reorder for samples

## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=evenness_pielou, fill=RootType)) + coord_cartesian(ylim = c(0.5, 1.0))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Pielou's evenness")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))


## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=evenness_pielou, fill=SoilP)) + coord_cartesian(ylim = c(0.5, 1))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Pielou's evenness")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))

model<-lme(evenness_pielou~SoilP,random=~1|Rep,data=mydata)
anova(model)
##
model<-lme(evenness_pielou~SoilP*RootType,random=~1|Rep,data=mydata)
anova(model)
# qqnorm(residuals(model))
# plot(residuals(model))
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/evenness_anova_Root2019L.csv")



# Post Hoc Test- pairwase comparison for all sample types
model<-lme(evenness_pielou~Tret,random=~1|Rep,data=mydata)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_evenness_PostHoc_Root_6wk.csv")


model<-lme(evenness_pielou~RootType,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"RootType")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_evenness_PostHoc_Root_all_RootType.csv")

###########################################################################


















# ##################################################################
# View (combined) 
# sub_root00 <- subset_samples(combined, SampleType == "Root")
# View (sub_root00)
# sub_root01 <- subset_samples(sub_root00, YS == "2019L")
# sub_root01 <- subset_samples(sub_root00, YS == "2019F")
# sub_root01 <- subset_samples(sub_root00, YS == "2020L")
# sub_root01 <- subset_samples(sub_root00, YS == "2020F")
# View (sub_root01)
# 
# # GPr  = transform_sample_counts(sub_root01, function(x) x / sum(x) )
# # View (GPr)
# # GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)   ## Keep ASVs with a mean greater than 10^-5
# ## Order Class from high to low
# Sub_Rarefy_Class<-tax_glom(sub_root01, taxrank = "Class") # combine all the Family belonging to "others"
# A<-names(sort(rowSums(otu_table(Sub_Rarefy_Class)),decreasing=T))
# otu_table(Sub_Rarefy_Class)<-otu_table(Sub_Rarefy_Class)[intersect(A,rownames(otu_table(Sub_Rarefy_Class))),]
# tax_table(Sub_Rarefy_Class)<-tax_table(Sub_Rarefy_Class)[intersect(A,rownames(tax_table(Sub_Rarefy_Class))),]
# View(A)
# View(otu_table(Sub_Rarefy_Class))
# View(tax_table(Sub_Rarefy_Class)) ## Check the number : 43 for Rarefying at 12000
# View (Sub_Rarefy_Class)
# ## Replace Family with low reads with "others" !!!!!!!!!!!!!!!!!!!!! low reads number?? how to define? 
# PGroup <- transform_sample_counts(Sub_Rarefy_Class, function(x)100* x / sum(x))
# 
# p<-plot_bar(PGroup,"Tret_POW", fill = "Class")+geom_bar(aes(fill=Class), stat="identity", position="stack")+
#   scale_fill_manual(values = getPalette)+scale_y_continuous(labels=sprintf("%1.0f",seq(0,100,25)))
# 
# ########## 
# View (combined) 
# sub_root <- subset_samples(combined, SampleType == "Root")
# View (sub_root)
# sub_root01 <- subset_samples(sub_root, YS == "2019L")
# sub_root01 <- subset_samples(sub_root, YS == "2019F")
# sub_root01 <- subset_samples(sub_root, YS == "2020L")
# sub_root01 <- subset_samples(sub_root, YS == "2020F")
# View (sub_root01)
# 
# ##
# TGroup <- tax_glom(sub_root01, taxrank = "Class")
# # View (TGroup)
# PGroup <- transform_sample_counts(TGroup, function(x) log10(1 + x))
# View (PGroup)
# 
# ##
# ## Extract abundance matrix from the phyloseq object
# OTU1 = as(otu_table(PGroup), "matrix")
# # # transpose if necessary
# if(taxa_are_rows(PGroup)){OTU1 <- t(OTU1)}
# # # Coerce to data.frame
# OTUdf = as.data.frame(OTU1)
# # write.csv (OTUdf, "Relative abundance_percentage_ITS.csv")
# RA <- read.csv("Relative abundance_percentage_ITS.csv")
# # a33 <- read.csv ("metaRDA_16s.csv")
# # a3 <- subset (a33,SampleType=="Root")
# RA1 <- subset(a3, select= c(SampleID,Year,Rep,Tret,Tret_POW,Plevel,OW,Stage))
# # View (RA1)
# ## Create a table for correlation analysis
# RA2 <- merge(RA, RA1, by = 'SampleID', na.rm = TRUE)  ## Remove rows have NA
# View (RA2)
# ## Change first column to rownames
# RA3 <- data.frame(RA2, row.names = 1)
# View (RA3)
# names (RA3)
# model<-lm(ASV5~Tret_POW,data=RA3)
# anova (model)
# 
# # # Tret_POW = sample_data(PGroup)$Tret_POW
# # # table(Tret_POW)
# # mean_PGroup = sapply(levels(Tret_POW),function(i){
# #   rowMeans(otu_table(PGroup)[,Tret_POW==i])
# # })
# # phy = tax_table(PGroup)[rownames(mean_PGroup ),"Class"]
# # rownames(mean_PGroup) = phy
# # View (phy)
# # head(sort(rowMeans(mean_PGroup),decreasing=TRUE))
# # head (phy)
# # View (mean_PGroup)
# # test00 <- names(sort(rowMeans(mean_PGroup),decreasing=T))
# # 
# # ######
# # ps.rel = transform_sample_counts(TGroup, function(x)100* x / sum(x))
# # # agglomerate taxa
# # glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)
# # ps.melt <- psmelt(glom)
# # # change to character for easy-adjusted level
# # ps.melt$Phylum <- as.character(ps.melt$Phylum)
# # ps.melt <- ps.melt %>%
# #   group_by(Tret_POW, Phylum) %>%
# #   mutate(median=median(Abundance))
# # # select group median > 1
# # keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
# # ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
# # #to get the same rows together
# # ps.melt_sum <- ps.melt %>%
# #   group_by(Sample,Tret_POW,Phylum) %>%
# #   summarise(Abundance=sum(Abundance))
# # View (ps.melt_sum)
# # write.csv (ps.melt_sum,"Relative abundance_Phylum.csv")
# ################################################################################################################################################################################











################################################################################################################################################################################
############################# Alpha diversity Correlation ############################################################################################################################
## http://www.sthda.com/english/wiki/wiki.php?id_contents=7786
# install.packages("Hmisc")
library(Hmisc)
## Load table of Alpha diversity indices 
# bac_alpha<-microbiome::alpha(combined,index="all")
# write.csv(bac_alpha,file="Fungi-AlphaIndices-Redo2023.csv")
View (combined)
a1 <- read.csv("Fungi-AlphaIndices-Redo2023.csv")
View (a1)
names (a1)
## Select certain Alpha diversity indices
a2 <- a1[,c(1,2,6,10)]
# View (a2)

## Load table of Soil+Plant performance
a33 <- read.csv ("metaRDA_ITS.csv")
a03 <- subset (a33,SampleType=="Root")   ## Subset root samples
a3 <- subset (a03,Stage=="L")
a3 <- subset (a03,Stage=="F")

a3 <- subset (a3,Year=="2019")
a3 <- subset (a3,Year=="2020")
View (a3)
# names (a3)
## 
a4 <- subset(a3, select= -c(Year,Rep,Tret,Tret_POW,Plevel,OW,SampleType,Stage))
# View (a4)

## Create a table for correlation analysis
a5 <- merge(a2, a4, by = 'SampleID', na.rm = TRUE)  ## Remove rows have NA
View (a5)
## Change first column to rownames
df1 <- data.frame(a5, row.names = 1)
View (df1)

# # ## Calculate Pearson's Correlation  
# cor1 <- cor(df1[,1:3], df1[,4:16], method = c("pearson"))  ## , "kendall", "spearman"
# cor2 <- rcorr(as.matrix(df1), type=c("pearson","spearman"))
# 
# ## Visualize a correlation matrix using a correlogram
# # install.packages("corrplot")
# # library(corrplot)
# testRes = cor.mtest(cor1, conf.level = 0.95)
# corrplot(cor1, p.mat = testRes$p, type="upper",
#          addCoef.col = 'grey35', insig='blank',
#          number.cex = 0.7, tl.col="black", tl.srt=45)


######## 
## Pearson's for alpha diversity, Spearman for relative abundace.
### For compositional data, better use Spearman correlation instead
## https://doi.org/10.1038/ismej.2015.235
## https://stackoverflow.com/questions/70430772/generate-correlation-matrix-with-specific-columns-and-only-with-significant-valu
## 
library(corrplot)
c_df <- Hmisc::rcorr(as.matrix(df1), type='pearson')  ## spearman
# corrplot(corr=c_df$r[1:3, 4:16], p.mat=c_df$P[1:3, 4:16], sig.level=0.05,
#          method='color', diag=FALSE, addCoef.col=1, type='upper', insig='blank',
#          number.cex=.8)
corrplot(corr=c_df$r[1:3, 4:16], p.mat=c_df$P[1:3, 4:16], sig.level=0.05,
         addCoef.col=1,  insig='blank',
         number.cex=.8, tl.srt=45,col=brewer.pal(n=6, name="RdYlBu"))

## Check sig. levels
m <- c_df$P[1:3, 4:16] < .05
m[lower.tri(m, diag=TRUE)] <- ''
as.data.frame(replace(m, lower.tri(m, diag=TRUE), ''))

# ## e.g.: Visualize a correlation matrix using scatter plots
library(PerformanceAnalytics)
# # mydata <- df[, c(1,3,4,5,6,7)]
chart.Correlation(df1, histogram=TRUE, pch=19)
# 
# ## e.g.: Visualize correlation matrix using a heatmap
# # Get some colors
# col<- colorRampPalette(c("blue", "white", "red"))(20)
# heatmap(x = cor1, col = col, symm = TRUE)
# ################################################################################################################################################################################
# install.packages("corrr")
# library ("corrr")
# x <- correlate(df1)   ## Method: 'pearson'; Missing treated using: 'pairwise.complete.obs'
# # names (df1)
# cor_root <- focus(x, observed, diversity_shannon, evenness_pielou) # Focus on correlations of mpg and cyl with all other variables
# ##############################################################################################################################











# ##############################################################################################################################
# ######### Core Microbiome #########
# ############################################
# #https://web.stanford.edu/class/bios221/Pune/Labs/Lab_microbiome/Lab_microbiome.html#core_microbiota 
# #https://microbiome.github.io/tutorials/Core.html   for reference
# BiocManager::install("microbiome")
# install.packages("emmeans")
# install.packages("devtools")
# devtools::install_github("mdancho84/tidyquant")
# install.packages("lubridate")
# install.packages("zoo")
# install.packages("xts")
# install.packages("TTR")
# install.packages("quantmod")
# install.packages("PerformanceAnalytics")
# library(BiocManager)
# library(microbiome)
# library(emmeans)
# library(lubridate)
# library(zoo)
# library(xts)
# library(TTR)
# library(quantmod)
# library(PerformanceAnalytics)
# library(tidyquant)
# library(RColorBrewer)
# library(reshape2)
# library(knitr)
# #######################################################
# ## Format data as phyloseq objects for reads > 2000 ###
# sub_OTU<-otu_table(as.matrix(sub_ASV),taxa_are_rows = TRUE)
# sub_TAX<-tax_table(as.matrix(Taxonomy_1))
# sub_samples<-sample_data(sub_metadata)
# sub_combind<-phyloseq(sub_OTU,sub_TAX,sub_samples)
# View(sub_combind) ### 285 samples
# ## https://microbiome.github.io/tutorials/Core.html
# ## https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/core-microbiota.html
# install.packages("devtools")
# devtools::install_github("microsud/jeevanuDB")
# library(jeevanuDB)
# library(microbiome) # data analysis and visualisation
# library(phyloseq) # also the basis of data object. Data analysis and visualisation
# library(RColorBrewer) # nice color options
# library(ggpubr) # publication quality figures, based on ggplot2
# library(dplyr) # data handling  
# install.packages("devtools")
# devtools::install_github("microsud/microbiomeutilities")
# 
# ###############################
# ## Root Core Microbiome ##
# # Filter the data to include only Root smples 
# phyloseq_Root<-subset_samples(combined,SampleType=="Root") # 10335 taxa and 95 samples
# # keep only taxa with positive sums
# phyloseq_Root <- prune_taxa(taxa_sums(phyloseq_Root) > 0, phyloseq_Root) # 5032 taxa and 95 samples
# # Calculate compositional version of the data
# # (relative abundances)
# phyloseq_Root.rel <- microbiome::transform(phyloseq_Root, "compositional")
# taxa_names(phyloseq_Root.rel)[1:3] ##"feature1" "feature2" "feature3"
# 
# # A full phyloseq object of the core microbiota is obtained as follows:
# Root.core <- core(phyloseq_Root.rel, detection = 0.0001, prevalence = .5) # 164 taxa and 95 samples
# # list for the core 164 taxa
# core.taxa.standard <- core_members(phyloseq_Root.rel, detection = 0.0001, prevalence = 50/100)
# #Retrieving the associated taxa names from the phyloseq object:
# core.taxa <- taxa(Root.core)
# class(core.taxa) ## "character"
# 
# # get the taxonomy data
# tax.mat <- tax_table(Root.core)
# tax.df <- as.data.frame(tax.mat) # From Kingdom to species
# # add the OTus/ASV to last column
# tax.df$OTU <- rownames(tax.df)
# # select taxonomy of only those OTUs that are core memebers based on the thresholds that were used.
# core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
# knitr::kable(head(core.taxa.class))
# 
# ## Core Line Plot
# # Determine core microbiota across various abundance/prevalence thresholds with the blanket analysis (Salonen et al. CMI, 2012) based on various signal and prevalences.
# # With compositional (relative) abundances
# det <- c(0, 0.1, 0.5, 2, 5, 20)/100
# prevalences <- seq(.05, 1, .05)
# plot_core(phyloseq_Root.rel, prevalences = prevalences, 
#           detections = det, plot.type = "lineplot") + 
#   xlab("Relative Abundance (%)") + theme_bw()  #### ?????? Plot Interpretation ????
# 
# ## Core Heatmap
# # Core with compositionals:
# prevalences <- seq(.05, 1, .05)
# #detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# # Also define gray color palette
# Root.core <- plot_core(phyloseq_Root.rel, 
#                        plot.type = "heatmap", 
#                        colours = rev(brewer.pal(5, "Spectral")),
#                        prevalences = prevalences, 
#                        detections = detections, 
#                        min.prevalence = .5) + 
#   xlab("Detection Threshold (Relative Abundance (%))")
# print(Root.core) 
# ############ 
# ## Assign All Taxa levels to features ##
# # get the data used for plotting 
# View (Root.core)
# df <- Root.core$data
# # get the list of OTUs
# list <- df$Taxa 
# # check the OTU ids
# print(list) 
# # get the taxonomy data
# tax <- as.data.frame(tax_table(phyloseq_Root.rel))
# # add the ASVs to last column
# tax$ASV <- rownames(tax)
# # select taxonomy of only those OTUs that are used in the plot
# tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
# # head(tax2)
# # We will merege all the column into one except the Doamin as all is bacteria in this case
# tax.unit <- tidyr::unite(tax2, Taxa_level,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE)
# tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
# # add this new information into the plot data df
# df$Taxa <- tax.unit$Taxa_level
# # you can see now we have the taxonomic information
# knitr::kable(head(df))
# # replace the data in the plot object
# Root.core$data <- df
# plot(Root.core + theme(axis.text.y = element_text(face="italic")))
# print (Root.core)
# ###########
# ## Add Family level (after checking all taxa and find uncultured at Genus level, so choose Family level)
# phyloseq_Root.rel.fam <- aggregate_taxa(phyloseq_Root.rel, "Family")
# prevalences <- seq(.05, 1, .05)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# Root.core.fam <- plot_core(phyloseq_Root.rel.fam, 
#                            plot.type = "heatmap", 
#                            colours = rev(brewer.pal(5, "RdBu")),
#                            prevalences = prevalences, 
#                            detections = detections, min.prevalence = .5) +
#   xlab("Detection Threshold (Relative Abundance (%))") + theme_bw() + ylab("Family")
# print (Root.core.fam)
# 
# ###############################
# ## Bulk Soil Core Microbiome ##
# # Filter the data to include only Bulk Soil smples 
# phyloseq_Bulk<-subset_samples(combined,SampleType=="Bulk") # 10335 taxa and 95 samples
# # keep only taxa with positive sums
# phyloseq_Bulk <- prune_taxa(taxa_sums(phyloseq_Bulk) > 0, phyloseq_Bulk) # 9317 taxa and 95 samples
# # Calculate compositional version of the data
# # (relative abundances)
# phyloseq_Bulk.rel <- microbiome::transform(phyloseq_Bulk, "compositional")
# taxa_names(phyloseq_Bulk.rel)[1:3] ##"feature1" "feature2" "feature3"
# 
# # A full phyloseq object of the core microbiota is obtained as follows:
# Bulk.core <- core(phyloseq_Bulk.rel, detection = 0.0001, prevalence = .5) # 593 taxa and 95 samples
# # list for the core 164 taxa
# core.taxa.standard <- core_members(phyloseq_Bulk.rel, detection = 0.0001, prevalence = 50/100)
# #Retrieving the associated taxa names from the phyloseq object:
# core.taxa <- taxa(Bulk.core)
# class(core.taxa) ## "character"
# # get the taxonomy data
# tax.mat <- tax_table(Bulk.core)
# tax.df <- as.data.frame(tax.mat) # From Kingdom to species
# # add the OTus/ASV to last column
# tax.df$OTU <- rownames(tax.df)
# # select taxonomy of only those OTUs that are core memebers based on the thresholds that were used.
# core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
# knitr::kable(head(core.taxa.class))
# 
# ## Core Heatmap for Bulk Soil
# # Core with compositionals:
# prevalences <- seq(.05, 1, .05)
# #detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# # Also define gray color palette
# Bulk.core <- plot_core(phyloseq_Bulk.rel, 
#                        plot.type = "heatmap", 
#                        colours = rev(brewer.pal(5, "Spectral")),
#                        prevalences = prevalences, 
#                        detections = detections, 
#                        min.prevalence = .5) + 
#   xlab("Detection Threshold (Relative Abundance (%))")
# print(Bulk.core) 
# ############ 
# ## Assign All Taxa levels to features ##
# # get the data used for plotting 
# View (Bulk.core)
# df <- Bulk.core$data
# # get the list of OTUs
# list <- df$Taxa 
# # check the OTU ids
# print(list) 
# # get the taxonomy data
# tax <- as.data.frame(tax_table(phyloseq_Bulk.rel))
# # add the ASVs to last column
# tax$ASV <- rownames(tax)
# # select taxonomy of only those OTUs that are used in the plot
# tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
# # head(tax2)
# # We will merege all the column into one except the Doamin as all is bacteria in this case
# tax.unit <- tidyr::unite(tax2, Taxa_level,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE)
# tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
# # add this new information into the plot data df
# df$Taxa <- tax.unit$Taxa_level
# # you can see now we have the taxonomic information
# knitr::kable(head(df))
# # replace the data in the plot object
# Bulk.core$data <- df
# plot(Bulk.core + theme(axis.text.y = element_text(face="italic")))
# print (Bulk.core)
# ###########
# ## Add Family level (after checking all taxa and find uncultured at Genus level, so choose Family level)
# phyloseq_Bulk.rel.fam <- aggregate_taxa(phyloseq_Bulk.rel, "Family")
# prevalences <- seq(.05, 1, .05)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# Bulk.core.fam <- plot_core(phyloseq_Bulk.rel.fam, 
#                            plot.type = "heatmap", 
#                            colours = rev(brewer.pal(5, "RdBu")),
#                            prevalences = prevalences, 
#                            detections = detections, min.prevalence = .5) +
#   xlab("Detection Threshold (Relative Abundance (%))") + theme_bw() + ylab("Family")
# print (Bulk.core.fam)
# 
# 
# ###############################
# ## Rhizo Soil Core Microbiome ##
# # Filter the data to include only Rhizo Soil smples 
# phyloseq_Rhizo<-subset_samples(combined,SampleType=="Rhizo") # 10335 taxa and 95 samples
# # keep only taxa with positive sums
# phyloseq_Rhizo <- prune_taxa(taxa_sums(phyloseq_Rhizo) > 0, phyloseq_Rhizo) # 9052 taxa and 95 samples
# # Calculate compositional version of the data
# # (relative abundances)
# phyloseq_Rhizo.rel <- microbiome::transform(phyloseq_Rhizo, "compositional")
# taxa_names(phyloseq_Rhizo.rel)[1:3] ##"feature1" "feature2" "feature3"
# # A full phyloseq object of the core microbiota is obtained as follows:
# Rhizo.core <- core(phyloseq_Rhizo.rel, detection = 0.0001, prevalence = .5) # 303 taxa and 95 samples
# # list for the core 303 taxa
# core.taxa.standard <- core_members(phyloseq_Rhizo.rel, detection = 0.0001, prevalence = 50/100)
# #Retrieving the associated taxa names from the phyloseq object:
# core.taxa <- taxa(Rhizo.core)
# class(core.taxa) ## "character"
# # get the taxonomy data
# tax.mat <- tax_table(Rhizo.core)
# tax.df <- as.data.frame(tax.mat) # From Kingdom to species
# # add the OTus/ASV to last column
# tax.df$OTU <- rownames(tax.df)
# # select taxonomy of only those OTUs that are core memebers based on the thresholds that were used.
# core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
# knitr::kable(head(core.taxa.class))
# 
# ## Core Heatmap for Rhizo Soil
# # Core with compositionals:
# prevalences <- seq(.05, 1, .05)
# #detections <- 10^seq(log10(1e-3), log10(.2), length = 10)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# # Also define gray color palette
# Rhizo.core <- plot_core(phyloseq_Rhizo.rel, 
#                         plot.type = "heatmap", 
#                         colours = rev(brewer.pal(5, "Spectral")),
#                         prevalences = prevalences, 
#                         detections = detections, 
#                         min.prevalence = .5) + 
#   xlab("Detection Threshold (Relative Abundance (%))")
# print(Rhizo.core) 
# ############ 
# ## Assign All Taxa levels to features ##
# # get the data used for plotting 
# View (Rhizo.core)
# df <- Rhizo.core$data
# # get the list of OTUs
# list <- df$Taxa 
# # check the OTU ids
# print(list) 
# # get the taxonomy data
# tax <- as.data.frame(tax_table(phyloseq_Rhizo.rel))
# # add the ASVs to last column
# tax$ASV <- rownames(tax)
# # select taxonomy of only those OTUs that are used in the plot
# tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
# # head(tax2)
# # We will merege all the column into one except the Doamin as all is bacteria in this case
# tax.unit <- tidyr::unite(tax2, Taxa_level,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE)
# tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
# # add this new information into the plot data df
# df$Taxa <- tax.unit$Taxa_level
# # you can see now we have the taxonomic information
# knitr::kable(head(df))
# # replace the data in the plot object
# Rhizo.core$data <- df
# plot(Rhizo.core + theme(axis.text.y = element_text(face="italic")))
# print (Rhizo.core)
# ###########
# ## Add Family level (after checking all taxa and find uncultured at Genus level, so choose Family level)
# phyloseq_Rhizo.rel.fam <- aggregate_taxa(phyloseq_Rhizo.rel, "Family")
# prevalences <- seq(.05, 1, .05)
# detections <- round(10^seq(log10(1e-3), log10(.2), length = 10), 3)
# Rhizo.core.fam <- plot_core(phyloseq_Rhizo.rel.fam, 
#                             plot.type = "heatmap", 
#                             colours = rev(brewer.pal(5, "RdBu")),
#                             prevalences = prevalences, 
#                             detections = detections, min.prevalence = .5) +
#   xlab("Detection Threshold (Relative Abundance (%))") + theme_bw() + ylab("Family")
# print (Rhizo.core.fam)
# 
# 
# 
# 
# ###################################
# ######### Heatmap #########
# ###################################
# #Write code for saving pheatmap figure
# save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
#   png(filename, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()}
# heatmap_Root<-subset_samples(sub_combind,RootSoil=="Root")
# 
# hp_tax<-tax_table(heatmap_Root)
# hp_tax<-as.data.frame(hp_tax)
# #remove the unclassified genus in taxanomy data
# hp_tax1<-hp_tax[!(hp_tax$Genus %in% c("Unclassified")),]
# 
# #remove the corresponded ASVs in ASV table
# ASV_Root<-as.data.frame(otu_table(heatmap_Root))
# ASV_Root1<-ASV_Root[intersect(rownames(ASV_Root),rownames(hp_tax1)),]
# #
# AVS_Root2<-cbind(hp_tax1$Genus,ASV_Root1)
# colnames(AVS_Root2)[1]<-"Genus"
# AVS_Root2<-AVS_Root2%>%group_by(Genus)%>%summarize_all(sum)
# AVS_Root2<-data.frame(AVS_Root2,row.names = 1)
# AVS_Root3<-t(AVS_Root2)
# Meta_root<-sample_data(heatmap_Root)
# identical(rownames(Meta_root),rownames(AVS_Root3))
# ASV_Meta<-cbind.data.frame(Meta_root$Variety,AVS_Root3)
# colnames(ASV_Meta)[1]<-"Variety"
# ASV_Meta1<-ASV_Meta%>%group_by(Variety)%>%summarize_all(mean)
# ASV_Meta1<-data.frame(ASV_Meta1,row.names = 1)
# p<-pheatmap(ASV_Meta1,cluster_cols = F,labels_col = "Genus")
# save_pheatmap_png(p,"Root_Genus.png",width=700,height=500)









###############################################################
##################*** Statistical Analysis ***##################
##### Alpha_diversity - Peter's Method (rarefied data)#########
install.packages("devtools")
devtools::install_github("petersonR/bestNormalize")
install.packages("emmeans")
install.packages("multcompView")
library(nlme)
library(permute)
library(lattice)
library(vegan)
library(bestNormalize)
library(emmeans)
library(multcompView)
library(mvtnorm)
library(survival)
library(TH.data)
library(MASS)
library(multcomp)

######## Rarefy at Same Depth ######
## Resample an OTU table such that all samples have the same library size.
## (2) https://microbiome.github.io/tutorials/Alphadiversity.html
library(microbiome)
library(knitr)
library(tidyverse)
# ## If only interested in richness: 
# rich<- richness(Sub_Rarefy) ##e.g.: observed, chao1
# ## If only interested in evenness: 
# even<-evenness(Sub_Rarefy) ## e.g.: pielou, simpson, etc
# # ## Include all index:
bac_alpha<-microbiome::alpha(Sub_Rarefy,index="all")
# kable(head(bac_alpha))

### First Rarefy to 12000 ######
Sub_Rarefy<-rarefy_even_depth(combined,sample.size=12000,rngseed=999) ##3582 taxa and 284 samples
########################################
##### Alpha diversity -- Observed #####
View(Sub_Rarefy)  ## 3392 x 284
bac_alpha<-microbiome::alpha(Sub_Rarefy,index="all")
View(bac_alpha)
names(bac_alpha)
observed<-subset(bac_alpha,select="observed")
mydata<-data.frame(sample_data(Sub_Rarefy),observed)
View(mydata)
# ## model for All SampleTypes
# model<-lme(observed~SampleType,random=~1|Rep,data=mydata)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(residuals(model))
# SampleType<-anova(model)
# ## Post hoc test 
# model<-lme(observed~SampleType,random=~1|Rep,data=mydata)
# #a factorial (two-way) repeated measures ANOVA with lme(), 
# #and the F-test has revealed a main effect of at least one of the factors, and possibly an interaction.
# em<-emmeans(model,pairwise ~ SampleType) ##Estimated marginal means
# # All pairwise comparisons
# contrast(em)
# # pairs(em,adjust = "tukey") # Not what I want, cuz compare (Bulk - Rhizo) - (Bulk - Root)
# library(multcompView)
# library(multcomp)
# Multcomp<-cld(em$emmeans,alpha=0.05, Letters=letters)### Correct!!!!!!
# plot(Multcomp)
######
model<-lme(observed~SampleType*Year*Stage,random=~1|Rep,data=mydata)
anova(model)
write.csv(anova(model),file="ITS_observed_anova.csv")
#Anova(model,type="III")
##
model<-lme(observed~SampleType,random=~1|Rep,data=mydata)
em<-emmeans(model,"SampleType")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="ITS_observed_post hoc.csv")
#Anova(model,type="III")
# anova(model)
# em<-emmeans(model,pairwise ~ SampleType)
# contrast(em)

## (1) Bulk Soil
mydata_Bulk<-subset(mydata,SampleType=="Bulk")
model<-lme(observed~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Bulk)
anova(model)
#Anova(model,type="III")
qqnorm(residuals(model))
plot(residuals(model))
Bulk_observed<-anova(model)
write.csv(Bulk_observed,file="ITS-bulk_observed.csv")
# All pairwise comparisons
# em<-emmeans(model,pairwise ~ SampleType) ##Estimated marginal means
#contrast(em)
# library(multcompView)
# library(multcomp)
# #Multcomp<-cld(em$emmeans,alpha=0.05, Letters=letters)
# #plot(Multcomp)
# # Since sig difference between Year and no interaction of it, so analysis 2019 and 2020 seperately.
# # subset bulk2019 samples
# mydata_Bulk2019<-subset(mydata_Bulk,Year=="2019")
# View(mydata_Bulk2019) # 48 samples
# model<-lme(observed~Stage,random=~1|Rep,data=mydata_Bulk2019)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# Bulk2019_observed <- anova(model)
# write.csv(Bulk2019_observed,file="Fungi-Bulk2019_observed_anova.csv")
# # subset bulk2020 samples
# mydata_Bulk2020<-subset(mydata_Bulk,Year=="2020")
# View(mydata_Bulk2020) # 47 samples
# model<-lme(observed~Stage,random=~1|Rep,data=mydata_Bulk2020)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# Bulk2020_observed <- anova(model)
# write.csv(Bulk2020_observed,file="Fungi-Bulk2020_observed_anova.csv")

## (2) Rhizosphere Soil
mydata_Rhizo<-subset(mydata,SampleType=="Rhizo")
model<-lme(observed~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Rhizo)
anova(model)
qqnorm(residuals(model))
plot(residuals(model))
Rhizo_observed <- anova(model)
write.csv(Rhizo_observed,file="ITS-Rhizo_observed.csv")
# # subset Rhizo Vegetative stage samples
# mydata_RhizoL<-subset(mydata_Rhizo,Stage=="L")
# View(mydata_RhizoL) # 48 samples
# model<-lme(observed~Year*Plevel,random=~1|Rep,data=mydata_RhizoL)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoL_observed <- anova(model)
# write.csv(RhizoL_observed,file="Fungi-Rhizo_Vegetative_observed_anova.csv")
# ##Post Hoc Test
# model<-lme(observed~Year*Plevel,random=~1|Rep,data=mydata_RhizoL)
# em<-emmeans(model,pairwise ~ Year*Plevel) ##Estimated marginal means
# em
# # subset Rhizo Flowering stage samples
# mydata_RhizoF<-subset(mydata_Rhizo,Stage=="F")
# View(mydata_RhizoF) # 47 samples
# model<-lme(observed~Year+Plevel,random=~1|Rep,data=mydata_RhizoF)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoF_observed <- anova(model)
# write.csv(RhizoF_observed,file="Fungi-Rhizo_Flowering_observed_anova.csv")
# ##Post Hoc Test
# model<-lme(observed~Year+Plevel,random=~1|Rep,data=mydata_RhizoF)
# anova(model)
# em<-emmeans(model,pairwise ~Year+Plevel) ##Estimated marginal means
# subset Rhizo Flowering stage samples into different Years
# ## (a) Rhizo Flowering stage 2019
# mydata_RhizoF2019<-subset(mydata_RhizoF,Year=="2019")
# View(mydata_RhizoF2019) # 23 samples
# model<-lme(observed~Plevel*OW,random=~1|Rep,data=mydata_RhizoF2019)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# ## (b) Rhizo Flowering stage 2020
# mydata_RhizoF2020<-subset(mydata_RhizoF,Year=="2020")
# View(mydata_RhizoF2020) # 24 samples
# model<-lme(observed~Plevel*OW,random=~1|Rep,data=mydata_RhizoF2020)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# ## (a) Rhizo Flowering stage 2019
# mydata_Rhizo2019<-subset(mydata_Rhizo,Year=="2019")
# View(mydata_Rhizo2019) # 47 samples
# model<-lme(observed~Stage*Plevel*OW,random=~1|Rep,data=mydata_Rhizo2019)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# ## (b) Rhizo Flowering stage 2020
# mydata_Rhizo2020<-subset(mydata_Rhizo,Year=="2020")
# View(mydata_Rhizo2020) # 24 samples
# model<-lme(observed~Stage*Plevel*OW,random=~1|Rep,data=mydata_Rhizo2020)
# #anova(model)
# Anova(model,type="III")
# ##Post Hoc Test
# model<-lme(observed~Plevel,random=~1|Rep,data=mydata_RhizoF2020)
# em<-emmeans(model,pairwise ~ Plevel) ##Estimated marginal means

## (3) Root
mydata_Root<-subset(mydata,SampleType=="Root")
model<-lme(observed~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Root)
anova(model)
#Anova(model,type="III")
qqnorm(residuals(model))
plot(residuals(model))
Root_observed<-anova(model)
write.csv(Root_observed,file="ITS-Root_observed.csv")
# ##Post Hoc Test
# model<-lme(observed~Year*Stage*OW,random=~1|Rep,data=mydata_Root)
# em<-emmeans(model,pairwise ~ Year*Stage*OW) ##Estimated marginal means
# # All pairwise comparisons
# contrast(em)
# subset Root Vegetative stage samples
mydata_RootL<-subset(mydata_Root,Stage=="L")
mydata_RootL<-subset(mydata_Root,Stage=="F")
##
mydata_RootL<-subset(mydata_RootL,Year=="2019")
mydata_RootL<-subset(mydata_RootL,Year=="2020")

View(mydata_RootL) # 47 samples
model<-lme(observed~OW,random=~1|Rep,data=mydata_RootL)
#anova(model)
Anova(model,type="III")
qqnorm(residuals(model))
plot(model)
plot(residuals(model))
RootL_observed <- anova(model)
write.csv(RootL_observed,file="Fungi-Root_Vegetative_observed_anova.csv")
##Post Hoc Test
model<-lme(observed~P_OW,random=~1|Rep,data=mydata_RootL)
em<-emmeans(model,pairwise ~P_OW) ##Estimated marginal means
em
write.csv(em,file="Fungi-Root2019L_observed_pairwise.csv")

# subset Root Flowering stage samples
mydata_RootF<-subset(mydata_Root,Stage=="F")
View(mydata_RootF) # 48 samples
model<-lme(observed~Year+OW,random=~1|Rep,data=mydata_RootF)
#anova(model)
Anova(model,type="III")
qqnorm(residuals(model))
plot(model)
plot(residuals(model))
RootF_observed <- anova(model)
write.csv(RootF_observed,file="Fungi-Root_Flowering_observed_anova.csv")
##Post Hoc Test
model<-lme(observed~Year+Plevel,random=~1|Rep,data=mydata_RootF)
#anova(model)
Anova(model,type="III")
em<-emmeans(model,pairwise ~Year+Plevel) ##Estimated marginal means





# ##### Alpha diversity -- Chao1 #####
# View(bac_alpha)
# write.csv(bac_alpha,file="Fungi-AlphaAllIndex.csv")
# names(bac_alpha)
# chao1<-subset(bac_alpha,select="chao1")
# mydata<-data.frame(sample_data(Sub_Rarefy),chao1)
# View(mydata)
# ## Year_SampleType_Stage
# model<-lme(chao1~SampleType*Year*Stage,random=~1|Rep,data=mydata)
# anova(model)
# qqnorm(residuals(model))
# plot(residuals(model))
# ## Post hoc test
# model<-lme(chao1~SampleType,random=~1|Rep,data=mydata)
# #anova(model)
# Anova(model,type="III")
# em<-emmeans(model,pairwise ~ SampleType) ##Estimated marginal means
# contrast(em)
# library(multcompView)
# library(multcomp)
# Multcomp<-cld(em$emmeans,alpha=0.05, Letters=letters)
# plot(Multcomp)
# #####
# # model<-lme(chao1~SampleType*Year*Stage,random=~1|Rep,data=mydata)
# # em<-emmeans(model,pairwise ~ SampleType+Year*Stage)
# # contrast(em)
# 
# ## (1) Bulk Soil
# mydata_Bulk<-subset(mydata,SampleType=="Bulk")
# model<-lme(chao1~Year*Stage*Plevel*OW,random=~1|Rep,data=mydata_Bulk)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(residuals(model))
# Bulk_chao1 <- anova(model)
# write.csv(Bulk_chao1,file="Fungi-Bulk_chao1_anova.csv")
# ## Since sig difference between Year and no interaction of Year, so analysis 2019 and 2020 seperately.
# # subset bulk2019 samples
# mydata_Bulk2019<-subset(mydata_Bulk,Year=="2019")
# View(mydata_Bulk2019) # 48 samples
# model<-lme(chao1~Stage,random=~1|Rep,data=mydata_Bulk2019)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# Bulk2019_chao1 <- anova(model)
# write.csv(Bulk2019_chao1,file="Fungi-Bulk2019_chao1_anova.csv")
# # subset bulk2020 samples
# mydata_Bulk2020<-subset(mydata_Bulk,Year=="2020")
# View(mydata_Bulk2020) # 46 samples
# model<-lme(chao1~Stage,random=~1|Rep,data=mydata_Bulk2020)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# Bulk2020_chao1 <- anova(model)
# write.csv(Bulk2020_chao1,file="Fungi-Bulk2020_chao1_anova.csv")
# 
# ## (2) Rhizosphere Soil
# mydata_Rhizo<-subset(mydata,SampleType=="Rhizo")
# model<-lme(chao1~Year*Stage*Plevel*OW,random=~1|Rep,data=mydata_Rhizo)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(residuals(model))
# Rhizo_chao1 <- anova(model)
# write.csv(Rhizo_chao1,file="Fungi-Rhizo_chao1_anova.csv")
# # subset Rhizo Vegetative stage samples
# mydata_RhizoL<-subset(mydata_Rhizo,Stage=="L")
# View(mydata_RhizoL) # 48 samples
# model<-lme(chao1~Year*Plevel,random=~1|Rep,data=mydata_RhizoL)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoL_chao1 <- anova(model)
# write.csv(RhizoL_chao1,file="Fungi-Rhizo_Vegetative_chao1_anova.csv")
# # subset Rhizo Flowering stage samples
# mydata_RhizoF<-subset(mydata_Rhizo,Stage=="F")
# View(mydata_RhizoF) # 47 samples
# model<-lme(chao1~Year*Plevel,random=~1|Rep,data=mydata_RhizoF)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoF_chao1 <- anova(model)
# write.csv(RhizoF_chao1,file="Fungi-Rhizo_Flowering_chao1_anova.csv")
# 
# ## (3) Root
# mydata_Root<-subset(mydata,SampleType=="Root")
# model<-lme(chao1~Year*Stage*Plevel*OW,random=~1|Rep,data=mydata_Root)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(residuals(model))
# Root_chao<-anova(model)
# write.csv(Root_chao,file="Fungi-Root_chao1_anova.csv")
# ##Post Hoc
# model<-lme(chao1~Year*Stage*OW,random=~1|Rep,data=mydata_Root)
# em<-emmeans(model,pairwise ~ Year*Stage*OW) ##Estimated marginal means
# write.csv(em$contrasts,file="Fungi-Root_chao1_post.csv")
# # subset Root Vegetative stage samples
# mydata_RootL<-subset(mydata_Root,Stage=="L")
# View(mydata_RootL) # 47 samples
# model<-lme(chao1~Year,random=~1|Rep,data=mydata_RootL)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RootL_chao1 <- anova(model)
# write.csv(RootL_chao1,file="Fungi-Root_Vegetative_chao1_anova.csv")
# # subset Root Flowering stage samples
# mydata_RootF<-subset(mydata_Root,Stage=="F")
# View(mydata_RootF) # 48 samples
# model<-lme(chao1~Year,random=~1|Rep,data=mydata_RootF)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RootF_chao1 <- anova(model)
# write.csv(RootF_chao1,file="Fungi-Root_Flowering_chao1_anova.csv")


## Double check for Shannon Diversity ####### Correct ## Same with Peter's method #####
# shannon<-diversity(otu_table(Sub_Rarefy),"shannon")
# View(shannons)
# mydata<-data.frame(sample_data(Sub_Rarefy),shannons)
# mydata_Root<-subset(mydata,SampleType=="Root")
# View(mydata_Root)
# ##Year_SampleType
# model<-lme(diversity_shannon~SampleType*Year*Stage,random=~1|Rep,data=mydata)
# anova(model)
#@@@@@@@@@@@@@ %%%%%%%%%%%%%%%% &&&&&&&&&&& ################
#### Alpha diversity- shannon -Peter's method ####
View (Sub_Rarefy) ## 284 samples
# x<-t(otu_table(Sub_Rarefy))
View (otu_table(Sub_Rarefy))
shannon<-subset(bac_alpha,select="diversity_shannon")
View(shannon)
mydata<-data.frame(sample_data(Sub_Rarefy),shannon)
View(mydata)
# All sample types 
model<-lme(diversity_shannon~SampleType*Year*Stage,random=~1|Rep,data=mydata)
anova(model)
summary(model)
qqnorm(residuals(model))
qqnorm(residuals(model));abline(0,1)
plot(residuals(model))
write.csv(anova(model),file="ITS-shannon-anova.csv")
##Post doc test
model<-lme(diversity_shannon~SampleType,random=~1|Rep,data=mydata)
##Estimated marginal means
em<-emmeans(model,"SampleType")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="ITS_shannon_post hoc.csv")
#####
# model<-lme(shannon~SampleType+Year*Stage,random=~1|Rep,data=mydata)
# Anova(model,type="III") ### ANOVA Type III (the sequence of factor doesn't matter)
# em<-emmeans(model,pairwise ~ SampleType+Year*Stage)
# contrast(em)
# library(multcomp)
# summary(glht(model, linfct = mcp(Group = "Tukey")), test = adjusted("holm"))

### (1)Bulk Soil
mydata_Bulk<-subset(mydata,SampleType=="Bulk")
model<-lme(diversity_shannon~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Bulk)
anova(model)
#Anova(model,type="III")
qqnorm(residuals(model))
plot(model)
plot(residuals(model))
Bulk_shannon<-anova(model)
write.csv(Bulk_shannon,file="ITS-Bulk_shannon.csv")
View(mydata_Bulk) # 94 samples
# ## Post Hoc Test
# model<-lme(shannon~Year*Stage,random=~1|Rep,data=mydata_Bulk)
# em<-emmeans(model,pairwise~Year*Stage)

### (2)Rhizosphere
mydata_Rhizo<-subset(mydata,SampleType=="Rhizo")
model<-lme(diversity_shannon~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Rhizo)
anova(model)
#Anova(model,type="III")
qqnorm(residuals(model))
plot(residuals(model))
Rhizo<-anova(model)
write.csv(Rhizo,file="ITS-Rhizo_shannon.csv")
# subset Rhizo Vegetative Stage samples
# mydata_RhizoL<-subset(mydata_Rhizo,Stage=="L")
# View(mydata_RhizoL) # 48 samples
# model<-lme(shannon~Year*Plevel,random=~1|Rep,data=mydata_RhizoL)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoL_shannon <- anova(model)
# write.csv(RhizoL_shannon,file="Fungi-Rhizo_Vegetative_shannon_anova.csv")
# 
# # subset Rhizo Flowering Stage samples
# mydata_RhizoF<-subset(mydata_Rhizo,Stage=="F")
# View(mydata_RhizoF) # 47 samples
# model<-lme(shannon~Year*Plevel,random=~1|Rep,data=mydata_RhizoF)
# #anova(model)
# Anova(model,type="III")
# qqnorm(residuals(model))
# plot(model)
# plot(residuals(model))
# RhizoF_shannon <- anova(model)
# write.csv(RhizoF_shannon,file="Fungi-Rhizo_Flowering_shannon_anova.csv")

### (3)Root
mydata_Root<-subset(mydata,SampleType=="Root")
model<-lme(diversity_shannon~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Root)
anova(model)
qqnorm(residuals(model))
plot(residuals(model))
Root_shannon<-anova(model)
write.csv(Root_shannon,file="ITS-Root_shannon.csv")
# #Post Hoc Test
# model<-lme(shannon~Year*Stage*OW,random=~1|Rep,data=mydata_Root)
# em<-emmeans(model,pairwise~Year*Stage*OW)


###### Alpha diversity -- Pielou evenness ######
evenness<-subset(bac_alpha,select="evenness_pielou")
mydata<-data.frame(sample_data(Sub_Rarefy),evenness)
View(mydata)
model<-lme(evenness_pielou~SampleType*Year*Stage,random=~1|Rep,data=mydata)
anova(model)
qqnorm(residuals(model))
plot(residuals(model))
write.csv(anova(model),file="ITS_evenness.csv")
# Post Hoc Test- pairwase comparison for all sample types
model<-lme(evenness_pielou~SampleType,random=~1|Rep,data=mydata)
em<-emmeans(model,"SampleType")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="ITS_evenness_post hoc.csv")

## Bulk Soil
mydata_Bulk<-subset(mydata,SampleType=="Bulk")
model<-lme(evenness_pielou~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Bulk)
anova(model)
write.csv(anova(model),file="ITS_Bulk_evenness.csv")
## Rhizo Soil
mydata_Rhizo<-subset(mydata,SampleType=="Rhizo")
model<-lme(evenness_pielou~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Rhizo)
anova(model)
write.csv(anova(model),file="ITS_Rhizo_evenness.csv")
## Root
mydata_Root<-subset(mydata,SampleType=="Root")
model<-lme(evenness_pielou~Year*Stage*Tret_POW,random=~1|Rep,data=mydata_Root)
anova(model)
write.csv(anova(model),file="ITS_Root_evenness.csv")













##################################################################################
##################### Beta_diversity #################### 
## sub combind -- delete reads lower than 2000
# sub_OTU<-otu_table(as.matrix(sub_ASV),taxa_are_rows = TRUE)
# sub_TAX<-tax_table(as.matrix(Taxonomy_1))
# sub_samples<-sample_data(sub_metadata)
# sub_combind<-phyloseq(sub_OTU,sub_TAX,sub_samples)
###########################################################################
############## Normolization Method - PCA Plot ##############################
################### Peter adapted from Navid ####################################
install.packages("MASS")
install.packages("NADA")
install.packages("truncnorm")
install.packages("survival")
install.packages("zCompositions")
library(MASS)
library(phyloseq)
library(truncnorm)
library(survival)
library(NADA)
library(zCompositions)
library(emmeans)
#install.packages('devtools')
#devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
#if (!requireNamespace("BiocManager",quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager :: install("ALDEx2")
library(ALDEx2)
library(carData)
library(car)
library(CoDaSeq)
## https://ourcodingclub.github.io/tutorials/ordination/ 
View(combined) 
sub1 <-subset_samples(combined,Stage=="6wk")
ASV_Beta<-phyloseq::filter_taxa(sub1,function(x) mean(x)>1,TRUE)

# ## Include all samples, Remove OTUs with a mean read count across all samples less than or equal to 1.
# ASV_Beta<-phyloseq::filter_taxa(combined,function(x) mean(x)>1,TRUE)
#data transformation
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta)),method = "BL", output = "p-counts") ### No. corrected values: 853
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE) ## "clr": gives the centered log ratio transform; clrInv gives closed compositions with the given clr-transform
meta<-data.frame(sample_data(ASV_Beta))
taxa<-data.frame(tax_table(ASV_Beta))
#generate a new phyloseq with transformed data
OTU_clr<-otu_table(as.matrix(t(ASV_clr)),taxa_are_rows = TRUE) #format data as phyloseq objects
View (OTU_clr) ## mean>1: 706 x 288  ; mean>3: 430x 288 (PCA 11.4%+8.8%) ; mean>5: 351x288 (PCA 11.7%+9.6%)
TAX1<-tax_table(as.matrix(taxa))
samples1<-sample_data(meta)
combined_clr<-phyloseq(OTU_clr,TAX1,samples1)
View(combined_clr) ## 591 x 120
###### Use principal components analysis method for PCA plot ###### ))))))))))))))))))))))))))))))))))
# RDA: performs redundancy analysis, or optionally principal components analysis
ord_clr <- ordinate(combined_clr, "RDA")   ### Perform An Ordination On Phyloseq Data
#Examine eigenvalues and % prop. variance explained for each PC
head(ord_clr$CA)
## (1)Show proportion explained by samples
summary(ord_clr)
plot_ordination(combined_clr, ord_clr, type="samples", color="SampleType")+ geom_point(size=2)+labs(color="SampleType") +stat_ellipse(aes(group = SampleType), linetype = 2)+scale_x_continuous(breaks = seq(-6, 6, by = 2))+scale_y_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                                                                                               panel.background = element_blank())
## (2)Show proportion explained by taxa
# combined_clr_Phylum<-tax_glom(combined_clr, taxrank = "Phylum") 
# sample_data(combined_clr)$SampleType<-levels(sample_data(combined_clr_Phylum)$SampleType)[get_variable(combined_clr,"SampleType")]
# plot_ordination(combined_clr, ord_clr, type="taxa", color="Phylum",shape="SampleType")+ geom_point(size=3)+labs(shape="SampleType",color="Phylum")+theme_bw()
# #Add circle in plot: +stat_ellipse(aes(group = SampleType), linetype = 2)
# #plot_ordination(combined_clr, ord_clr, type="samples", color="YS",shape="SampleType")+ geom_point(size=3)+labs(shape="SampleType",color="YS")
# #plot_ordination(combined_clr, ord_clr, type="samples", color="Plevel",shape="OW")+ geom_point(size=3)+labs(shape="OW",color="Plevel")
# dev.off()
# print(plot(1))
###########

###############
##Rhizosphere Soil PCA
rhizo_clr<-subset_samples(combined_clr,SampleType=="Rhizo")
rhizo_clr<-subset_samples(combined_clr,SampleType=="Root")

# ord_clr_rhizo <- ordinate(rhizo_clr, "RDA")
# plot_ordination(rhizo_clr, ord_clr_rhizo, type="samples", color="Stage",shape="SoilP")+ geom_point(size=2)+labs(color="Stage",shape="SoilP")+theme_bw()+scale_x_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+
#   scale_y_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())

ord_clr_rhizo <- ordinate(rhizo_clr, "RDA")
plot_ordination(rhizo_clr, ord_clr_rhizo, type="samples", color="RootType",shape="SoilP")+ geom_point(size=2)+labs(color="RootType",shape="SoilP")+theme_bw()+scale_x_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+
  scale_y_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())

plot_ordination(rhizo_clr, ord_clr_rhizo, type="samples", color="SoilP",shape="RootType")+ geom_point(size=2)+labs(color="SoilP",shape="RootType")+theme_bw()+scale_x_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+
  scale_y_continuous(breaks = seq(-6, 6, by = 2),limits=c(-6, 6))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())


#ggsave("PCA_Site_Compartment.png")
##### separate rhizo soil into 2 Stages ##########
# (1) Rhizo 4wk
rhizo2019_clr<-subset_samples(rhizo_clr,Stage=="4wk")
rhizo2019_clr<-subset_samples(rhizo_clr,Stage=="6wk")
View (rhizo2019_clr)
# sample_data(rhizo2019_clr)$POW<-paste(sample_data(rhizo2019_clr)$Plevel,sample_data(rhizo2019_clr)$OW)
ord_clr_rhizo2019 <- ordinate(rhizo2019_clr, "RDA")
plot_ordination(rhizo2019_clr, ord_clr_rhizo2019, type="samples", color="SoilP",shape="RootType")+ geom_point(size=3)+labs(color="SoilP",shape="RootType")+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                             panel.background = element_blank())

plot_ordination(rhizo2019_clr, ord_clr_rhizo2019, type="samples", color="Tret",shape="Rep")+ geom_point(size=3)+labs(color="Tret",shape="Rep")+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                             panel.background = element_blank())




############################################################################
################### PREMERMOVA ########################??????????????????????????????????????
## https://ourcodingclub.github.io/tutorials/ordination/ 
View(combined)
sub1 <-subset_samples(combined,Stage=="6wk")
ASV_Beta<-phyloseq::filter_taxa(sub1,function(x) mean(x)>1,TRUE)

# # ## Include all samples, Remove OTUs with a mean read count across all samples less than or equal to 1.
# ASV_Beta<-filter_taxa(combined,function(x) mean(x)>1,TRUE)
# data transformation -- clr-transform
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta)),method = "SQ", output = "p-counts") # "cmultRepl": Bayesian-Multiplicative Replacement Of Count Zeros
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE)
sample_data(ASV_Beta)[]<-lapply(sample_data(ASV_Beta),factor)
metadata_Beta<-data.frame(sample_data(ASV_Beta))
# View(metadata_Beta)
# View(ASV_clr)
## BL Method
## ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta)),method = "BL", output = "p-counts") ### No. corrected values: 31503
## ord_clr <- ordinate(combined_clr, method="NMDS")   ### Perform An Ordination On Phyloseq Data
## sample_data(combined_clr)$YS<-paste(sample_data(combined_clr)$Year,sample_data(combined_clr)$Stage)
## plot_ordination(combined_clr, ord_clr, type="samples", color="SampleType",shape="YS")+ geom_point(size=3)+labs(shape="Year-Stage",color="SampleType")
## plot_ordination(combined_clr, ord_clr, type="samples", color="Year",shape="Stage")+ geom_point(size=3)+labs(shape="Stage",color="Year")
#ggsave("PCA_Site_Compartment.png")

#set Blocks for PREMERMOVA
set.seed(2021)
perm <- how(nperm = 999)
setBlocks(perm) <- with (metadata_Beta,Rep)
# check the significance of Stage, Year, and SampleType on bacterial community structure
fit1<-adonis2(ASV_clr~SampleType*SoilP*RootType,data=metadata_Beta,method = "euclidean",by="terms",permutations = 999)
write.csv(fit1,file="CLC-ITS-PERMANO-All_6wk.csv")



##Rhizosphere Soil + Root Beta Diversity #########
##  Remove OTUs with a mean read count across all samples less than or equal to 1.
#generate Rhizosphere Soi phyloseq
combined <- sub1

sub_rhizo<-subset_samples(sub1,SampleType=="Rhizo") 
sub_rhizo<-subset_samples(sub1,SampleType=="Root") 
ASV_Beta_rhizo<-filter_taxa(sub_rhizo,function(x) mean(x)>1,TRUE)
#data transformation
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta_rhizo)),method = "BL", output = "p-counts") ### No. corrected values: 1716
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE) ## "clr": gives the centered log ratio transform; clrInv gives closed compositions with the given clr-transform
meta<-data.frame(sample_data(ASV_Beta_rhizo))
taxa<-data.frame(tax_table(ASV_Beta_rhizo))
#generate a new phyloseq with transformed data
OTU_clr<-otu_table(as.matrix(t(ASV_clr)),taxa_are_rows = TRUE) #format data as phyloseq objects
# View (OTU_clr) 
TAX1<-tax_table(as.matrix(taxa))
samples1<-sample_data(meta)
rhizo_clr<-phyloseq(OTU_clr,TAX1,samples1)
#make ordinate for PCA plot
ord_clr <- ordinate(rhizo_clr, "RDA")   ### Perform An Ordination On Phyloseq Data
plot_ordination(rhizo_clr, ord_clr, type="samples", color="Stage",shape="Tret")+ geom_point(size=3)+labs(color="Stage",shape="Tret")+theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                       panel.background = element_blank())


## PERMANOVA for Rhizosphere soil
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta_rhizo)),method = "BL", output = "p-counts") ##No. corrected values 73
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE)
sample_data(ASV_Beta_rhizo)[]<-lapply(sample_data(ASV_Beta_rhizo),factor)
metadata_Beta_rhizo<-data.frame(sample_data(ASV_Beta_rhizo))
# View(metadata_Beta_rhizo)
#set Blocks
set.seed(2021)
perm <- how(nperm = 999)
setBlocks(perm) <- with (metadata_Beta_rhizo,Rep)
# check the significance of Stage, Year on bacterial community structure
fit1<-adonis2(ASV_clr~SoilP*RootType,data=metadata_Beta_rhizo,method = "euclidean",by="terms",permutations = perm)
# write.csv(fit1,file="CLC-ITS_Perman_Rhizo_6wk.csv")
write.csv(fit1,file="CLC-ITS_Perman_Root_6wk.csv")

# fit1<-adonis2(ASV_clr~Stage*SoilP*RootType,data=metadata_Beta_rhizo,method = "euclidean",by="terms",permutations = perm)
# write.csv(fit1,file="CLC-ITS_Perman_Root.csv")
########################################################

library(devtools) 
library(cluster)
library(pairwiseAdonis) # Citation:Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4
fit1<-adonis2(ASV_clr~SoilP*RootType,data=metadata_Beta_rhizo,method = "euclidean",by="terms",permutations = perm)
pair <- pairwise.adonis(ASV_clr,metadata_Beta_rhizo$RootType,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
write.csv (pair,"Pairwise PERMANO_Root_6wk.csv")



## Root Beta Diversity #########
##  Remove OTUs with a mean read count across all samples less than or equal to 1.
#generate Bulk Soi phyloseq 
View (combined)
sub_root<-subset_samples(combined,SampleType=="Root")  ##4468 taxa and 96 samples
# sub_root<-subset_samples(sub_root,Year=="2019")  ##4468 taxa and 96 samples
# sub_root<-subset_samples(sub_root,Year=="2020")  ##4468 taxa and 96 samples
# sub_root<-subset_samples(sub_root,Stage=="L")  ##4468 taxa and 96 samples
# sub_root<-subset_samples(sub_root,Stage=="F")  ##4468 taxa and 96 samples
sub_root <- subset_taxa(sub_root, Class == "Agaricomycetes")  ## 208 ASVs
sub_root <- subset_taxa(sub_root, Class == "Sordariomycetes")  ## 208 ASVs
sub_root <- subset_taxa(sub_root, Class == "Dothideomycetes")  ## 208 ASVs
sub_root <- subset_taxa(sub_root, Class == "Mortierellomycetes")  ## 208 ASVs
sub_root <- subset_taxa(sub_root, Class == "Olpidiomycetes")  ## 208 ASVs
View(sub_root)
ASV_Beta_root<-filter_taxa(sub_root,function(x) mean(x)>1,TRUE)
#data transformation
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta_root)),method = "BL", output = "p-counts") ### No. corrected values: 1716
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE) ## "clr": gives the centered log ratio transform; clrInv gives closed compositions with the given clr-transform
meta<-data.frame(sample_data(ASV_Beta_root))
taxa<-data.frame(tax_table(ASV_Beta_root))
#generate a new phyloseq with transformed data
OTU_clr<-otu_table(as.matrix(t(ASV_clr)),taxa_are_rows = TRUE) #format data as phyloseq objects
View (OTU_clr) ## 440 x 96
TAX1<-tax_table(as.matrix(taxa))
samples1<-sample_data(meta)
root_clr<-phyloseq(OTU_clr,TAX1,samples1)
#make ordinate for PCA plot
ord_clr <- ordinate(root_clr, "RDA")   ### Perform An Ordination On Phyloseq Data
# sample_data(root_clr)$YS<-paste(sample_data(root_clr)$Year,sample_data(root_clr)$Stage)
# ## Reorder legend
# sample_data(root_clr)$Stage <- factor(sample_data(root_clr)$Stage,      
#                                               levels = c("L", "F"))
# # plot_ordination(root_clr, ord_clr, type="samples", color="Year",shape="Stage")+ geom_point(size=3)+labs(shape="Stage",color="Year")
# plot_ordination(root_clr, ord_clr, type="samples", color="Tret_POW",shape="Stage")+ geom_point(size=3)+labs(shape="Stage",color="P(OW)")+
#   theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                     panel.background = element_blank())
plot_ordination(root_clr, ord_clr, type="samples", color="Plevel",shape="OW")+ geom_point(size=2)+labs(shape="OW",color="P rate")+
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank())

View (root_clr)
## PERMANOVA for Root
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta_root)),method = "BL", output = "p-counts")
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE)
sample_data(ASV_Beta_root)[]<-lapply(sample_data(ASV_Beta_root),factor)
metadata_Beta_root<-data.frame(sample_data(ASV_Beta_root))
View (metadata_Beta_root)
#set Blocks
set.seed(2021)
perm <- how(nperm = 999)
setBlocks(perm) <- with (metadata_Beta_root,Rep)
# check the significance of Stage, Year on bacterial community structure
# fit1<-adonis2(ASV_clr~Year*Stage*Tret,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
fit1<-adonis2(ASV_clr~Plevel*OW,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)  # Non sig. for different P levels in 2019F, 2020F, F  
# fit1<-adonis2(ASV_clr~Plevel,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
# fit1<-adonis2(ASV_clr~Year*Stage*Tret+Rep,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
write.csv(fit1,file="PERMANOVA_root_Sordariomycetes_2020F.csv")

library(devtools) 
library(cluster)
library(pairwiseAdonis) # Citation:Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4
pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$Plevel,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$OW,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$Tret,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$Rep,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
#write.csv(fit1,file="Bac_pairwise_PERMA_root2019FPlevel.csv")
write.csv(fit1,file="root_Agaricomycetes_Plevelpairwise_2019L.csv")


# # ## Subset root samples into different sampling years, different growth stages
# # ASV_Beta_root1<-subset_samples(ASV_Beta_root,Year=="2019")
# # ASV_Beta_root1<-subset_samples(ASV_Beta_root,Year=="2020")
# ASV_Beta_root2<-subset_samples(ASV_Beta_root1,Stage=="L")
# ASV_Beta_root2<-subset_samples(ASV_Beta_root,Stage=="F")
# #
# View (ASV_Beta_root2)
# View (ASV_Beta_root)
# 
# ## PERMANOVA for Root
# ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta_root2)),method = "BL", output = "p-counts")
# ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE)
# sample_data(ASV_Beta_root)[]<-lapply(sample_data(ASV_Beta_root2),factor)
# metadata_Beta_root<-data.frame(sample_data(ASV_Beta_root2))
# #set Blocks
# set.seed(2021)
# perm <- how(nperm = 999)
# setBlocks(perm) <- with (metadata_Beta_root,Rep)
# # check the significance of Stage, Year on bacterial community structure
# # fit1<-adonis2(ASV_clr~Stage*Tret,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
# # fit1<-adonis2(ASV_clr~Year*Tret,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
# fit1<-adonis2(ASV_clr~Tret*Rep,data=metadata_Beta_root,method = "euclidean",by="terms",permutations = perm)
# # # write.csv(fit1,file="Fungi_PERMA_root2019_StageP.csv")
# # # write.csv(fit1,file="Fungi_PERMA_root2020_StageP.csv")
# # write.csv(fit1,file="Fungi_PERMA_root2020FP.csv")
# # 
# # ## Pairwise comparison
# # pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$Rep,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
# # pair <- pairwise.adonis(ASV_clr,metadata_Beta_root$Plevel,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
# # # write.csv(pair,"Fungi_PERMA_root2019_StageP-pairwise.csv")
# # # write.csv(pair,"Fungi_PERMA_root2020_StageP-pairwise.csv")
# # View (metadata_Beta_root)



##############################################################################
sub_OTU<-otu_table(as.matrix(sub_ASV),taxa_are_rows = TRUE)
sub_TAX<-tax_table(as.matrix(Taxonomy_1))
sub_samples<-sample_data(sub_metadata)
sub_combind<-phyloseq(sub_OTU,sub_TAX,sub_samples)

########### original phyloseq (delete low reads <2000)###############
otus<-otu_table(as.matrix(sub_ASV),taxa_are_rows = TRUE)
tax<-tax_table(as.matrix(Taxonomy_1))
map<-sample_data(sub_metadata)
sub_combind<-phyloseq(sub_OTU,sub_TAX,sub_samples)

##### log-transformed phyloseq ######
OTU_clr<-otu_table(as.matrix(t(ASV_clr)),taxa_are_rows = TRUE) #format data as phyloseq objects
TAX1<-tax_table(as.matrix(taxa))
samples1<-sample_data(meta)
combined_clr<-phyloseq(OTU_clr,TAX1,samples1)








######################################################################
######### Alpha_diversity plot Combined boxplot and dotplot_20220620 #########
#####################################################
### 16s Bacteria--Alpha Diversity-Boxplot with dots
Sub_Rarefy<-rarefy_even_depth(combined,sample.size=12000,rngseed=999)##Sub_Rarefy<-rarefy_even_depth(sub_combind,sample.size=median(sort(sample_sums(sub_combind))),rngseed=999)
View(Sub_Rarefy)
##all alpha diversity index
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy,SampleType=="Root")
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy_Bulk,Year=="2019")
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy_Bulk,Year=="2020")
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy_Bulk,Stage=="L")
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy_Bulk,Stage=="F")
View(Sub_Rarefy_Bulk)
write.csv(sample_data(Sub_Rarefy_Bulk),file="/Users/mel270/Desktop/CARPITS/Fungi-metadata-rarefied-Root2020F.csv")
fun_alpha<-microbiome::alpha(Sub_Rarefy_Bulk,index="all")
write.csv(fun_alpha,file="Fun-alpha_allIndex_Root2020F.csv")  ### Add SampleID in this csv
alpha_diversity <- read.csv ("/Users/mel270/Desktop/CARPITS/Fun-alpha_allIndex_Root2020F.csv") 
names(alpha_diversity)
View(alpha_diversity)
metadata<-read.csv("/Users/mel270/Desktop/CARPITS/Fungi-metadata-rarefied-Root2020F.csv")
View(metadata)
metadata_alpha <- inner_join(metadata, alpha_diversity,
                             by=c('SampleID'))

# ## By P_OW + Stage
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=observed, fill=P_OW)) + coord_cartesian(ylim = c(0, 1800))+ 
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19,aes(color=as.factor(Stage)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Growth stage",y="Observed richness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                        panel.background = element_blank())
# 
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=diversity_shannon, fill=P_OW)) + coord_cartesian(ylim = c(0, 7))+ 
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19,aes(color=as.factor(Stage)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Growth stage",x="P_OW",y="Shannon diversity")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                                 panel.background = element_blank())
# 
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=evenness_pielou, fill=P_OW)) + coord_cartesian(ylim = c(0, 1))+ 
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25,size=1, shape=19,aes(color=as.factor(Stage)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Growth stage",x="P_OW",y="Pielou's evenness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                                 panel.background = element_blank())

# ###### By P_OW + Year
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=observed, fill=P_OW)) + coord_cartesian(ylim = c(0, 1800))+
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19,aes(color=as.factor(Year)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Year",y="Observed richness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                        panel.background = element_blank())
# 
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=diversity_shannon, fill=P_OW)) + coord_cartesian(ylim = c(0, 7))+
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19,aes(color=as.factor(Year)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Year",x="P_OW",y="Shannon diversity")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                                 panel.background = element_blank())
# 
# metadata_alpha %>%
#   ggplot(aes(x=P_OW, y=evenness_pielou, fill=P_OW)) + coord_cartesian(ylim = c(0, 1))+
#   geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
#   geom_jitter(show.legend=TRUE, width=0.25,size=1, shape=19,aes(color=as.factor(Year)))+
#   stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
#                color="black", width=0.6, size=0.5)+labs(color="Year",x="P_OW",y="Pielou's evenness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                                                 panel.background = element_blank())

## By P_OW + Year + Stage
metadata_alpha %>%
  ggplot(aes(x=P_OW, y=observed, fill=P_OW)) + coord_cartesian(ylim = c(0, 1800))+ 
  geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
  geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19)+
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
               color="black", width=0.6, size=0.5)+labs(y="Observed richness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                  panel.background = element_blank())

metadata_alpha %>%
  ggplot(aes(x=P_OW, y=diversity_shannon, fill=P_OW)) + coord_cartesian(ylim = c(0, 7))+ 
  geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
  geom_jitter(show.legend=TRUE, width=0.25, size=1, shape=19)+
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
               color="black", width=0.6, size=0.5)+labs(x="P_OW",y="Shannon diversity")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank())

metadata_alpha %>%
  ggplot(aes(x=P_OW, y=evenness_pielou, fill=P_OW)) + coord_cartesian(ylim = c(0, 1))+ 
  geom_boxplot(outlier.shape = NA, alpha=0.5, width=0.6)+
  geom_jitter(show.legend=TRUE, width=0.25,size=1, shape=19)+
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar",
               color="black", width=0.6, size=0.5)+labs(x="P_OW",y="Pielou's evenness")+ theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                           panel.background = element_blank())

####### Alpha diversity Statistic analysis #######
## Observed richness
observed<-subset(alpha_diversity,select="observed")
mydata<-data.frame(sample_data(Sub_Rarefy_Bulk),observed)
View(mydata)
####
model<-lme(observed~P_OW,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"P_OW")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="Fun_observed_Root2020F_posthoc.csv")

## Shannon Diversity
diversity_shannon<-subset(alpha_diversity,select="diversity_shannon")
mydata<-data.frame(sample_data(Sub_Rarefy_Bulk),diversity_shannon)
View(mydata)
####
model<-lme(diversity_shannon~P_OW,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"P_OW")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="Fun_shannon_Root2020F_posthoc.csv")

## Evenness
evenness_pielou<-subset(alpha_diversity,select="evenness_pielou")
mydata<-data.frame(sample_data(Sub_Rarefy_Bulk),evenness_pielou)
View(mydata)
####
model<-lme(evenness_pielou~P_OW,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"P_OW")
pair<-pairs(em,adjust = "tukey")
write.csv(pair,file="Fun_evenness_Root2019L_post hoc.csv")





######################### Differential Abundance Analysis ###########################
############# ANCOM-BC ############# ############# ############# #############
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC",force = TRUE)
library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
                    '#000', 'color': '#fff'});","}")))
?ancombc2

## Check phyloseq
View(combined)

# ## Filter low reads data
# sub<- phyloseq::filter_taxa(combined,function(x) mean(x)>1,TRUE)

## Subset samples
sub0 <-subset_samples(combined,SampleType=="Root")
# sub1 <-subset_samples(sub0,Stage=="L")
# sub1 <-subset_samples(sub0,Stage=="F")
sub1 <-subset_samples(sub0,Year=="2019")
sub1 <-subset_samples(sub0,Year=="2020")
# sub1 <-subset_samples(sub0,YS=="2019L")
# sub1 <-subset_samples(sub0,YS=="2019F")
# sub1 <-subset_samples(sub0,YS=="2020L")
# sub1 <-subset_samples(sub0,YS=="2020F")

sub1 <- sub0
View (sub1)
#####

# ## Aggregate to phylum level
# phylum_data = aggregate_taxa(sub1, "ASV") ##67 taxa and 26 samples
# ## The taxonomy table
# tax_mat = as(tax_table(phylum_data), "matrix")
## ANCOMBC
out = ancombc(phyloseq = sub1, formula = "Tret_POW",tax_level = "Species",
              p_adj_method = "holm",
              group = "Tret_POW", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.001, global = TRUE)

## ANCOMBC - Results
ancom = out$res
res_global = out$res_global
write.csv(res_global,"CARP_ITS_DA_Root_Tret.csv")
#write.csv(res_global,"CARP_ITS_DA_Root_L_Tret.csv")
#write.csv(res_global,"CARP_ITS_DA_Root_F_Tret.csv")
write.csv(res_global,"CARP_ITS_DA_Root2019_Tret.csv")
write.csv(res_global,"CARP_ITS_DA_Root2020_Tret.csv")
# write.csv(res_global,"CARP_ITS_DA_Root_2019L_Tret.csv")
# write.csv(res_global,"CARP_ITS_DA_Root_2019F_Tret.csv")
# write.csv(res_global,"CARP_ITS_DA_Root_2020L_Tret.csv")
# write.csv(res_global,"CARP_ITS_DA_Root_2020F_Tret.csv")


## Extract the "TRUE" ASVs
ancom <- read.csv("CARP_ITS_DA_Root_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_L_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_F_Tret.csv")
ancom <- read.csv("CARP_ITS_DA_Root2019_Tret.csv")
ancom <- read.csv("CARP_ITS_DA_Root2020_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_2019L_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_2019F_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_2020L_Tret.csv")
# ancom <- read.csv("CARP_ITS_DA_Root_2020F_Tret.csv")

attach (ancom)
#View (ancom)
ancom_TRUE<-ancom[ancom$diff_abn=="TRUE",] ## 12 genus different for P0 vs P35, which is the same as P0 vs P65
write.csv(ancom_TRUE,"CARP_ITS_DA_Root_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_L_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_F_Tret_TRUE.csv")
write.csv(ancom_TRUE,"CARP_ITS_DA_Root2019_Tret_TRUE.csv")
write.csv(ancom_TRUE,"CARP_ITS_DA_Root2020_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_2019L_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_2019F_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_2020L_Tret_TRUE.csv")
# write.csv(ancom_TRUE,"CARP_ITS_DA_Root_2020F_Tret_TRUE.csv")
detach (ancom)
############################################################################################
#############################################################################################
####### Construct differential table that includes all taxon levels ##########################
df01 <- read.csv("Fungi-allASV copy.csv")
View (df01)
# df00 <- read.csv("CARP_ITS_DA_Root_2019L_Tret_TRUE.csv")
# df00 <- read.csv("CARP_ITS_DA_Root_2019F_Tret_TRUE.csv")
# df00 <- read.csv("CARP_ITS_DA_Root_2020L_Tret_TRUE.csv")
# df00 <- read.csv("CARP_ITS_DA_Root_2020F_Tret_TRUE.csv")
# df00 <- read.csv("CARP_ITS_DA_Root2019_Tret_TRUE.csv")
# df00 <- read.csv("CARP_ITS_DA_Root2020_Tret_TRUE.csv")
df00 <- read.csv("CARP_ITS_DA_Root_L_Tret_TRUE.csv")
df00 <- read.csv("CARP_ITS_DA_Root_F_Tret_TRUE.csv")
df00 <- read.csv("CARP_ITS_DA_Root_Tret_TRUE.csv")
View (df00)
df02 <- merge(df00, df01, by = 'X')  ## Remove rows have NA
View (df02)
# write.csv (df02,"CARP_ITS_DA_Root_2019L_Tret_AllTaxon.csv")
# write.csv (df02,"CARP_ITS_DA_Root_2019F_Tret_AllTaxon.csv")
# write.csv (df02,"CARP_ITS_DA_Root_2020L_Tret_AllTaxon.csv")
# write.csv (df02,"CARP_ITS_DA_Root_2020F_Tret_AllTaxon.csv")
# write.csv (df02,"CARP_ITS_DA_Root2019_Tret_AllTaxon.csv")
# write.csv (df02,"CARP_ITS_DA_Root2020_Tret_AllTaxon.csv")
write.csv (df02,"CARP_ITS_DA_Root_L_Tret_AllTaxon.csv")
write.csv (df02,"CARP_ITS_DA_Root_F_Tret_AllTaxon.csv")
write.csv (df02,"CARP_ITS_DA_Root_Tret_AllTaxon.csv")








View (combined)  # 10335 x 288

##################################################################
############ ANCOM-BC Heatmap ###################################
################################################################
View (metadata)
# t1 <- read.csv("CARP_ITS_DA_Root_Tret_TRUE.csv") ## mean >1
t1 <- read.csv("CARP_ITS_DA_Root_L_Tret_TRUE.csv") ## mean >1
t1 <- read.csv("CARP_ITS_DA_Root_F_Tret_TRUE.csv") ## mean >1

View (t1)
ASV_test<-ASV_1[intersect(t1$X,rownames(Taxonomy_1)),] # The elements of intersect(x,y) are those elements in x and in y
# View (ASV_test)
#ASV_test<-ASV[,intersect(rownames(metadata),colnames(ASV))]
# phyloseq
OTU0<-otu_table(as.matrix(ASV_test),taxa_are_rows = TRUE)
TAX0<-phyloseq::tax_table(as.matrix(Taxonomy_1))
samples0<-sample_data(metadata)
combined_l<-phyloseq(OTU0,TAX0,samples0)
#View(combined_l) 
combined_lol <- subset_samples (combined_l,Controls=="N")
#View (combined_lol)
## Bulk soil
#combined_666 <- subset_samples(combined_lol, SampleType == "Bulk")
## Rhizo soil
#combined_666 <- subset_samples(combined_lol, SampleType == "Rhizo")
combined_333 <- subset_samples(combined_lol, SampleType == "Root")
combined_666 <- subset_samples(combined_333, Stage == "L")
combined_666 <- subset_samples(combined_333, Stage == "F")

# combined_666 <-subset_samples (combined_333,NAMline=="NoPlant")  # (1) Bulk 
# combined_666 <-subset_samples (combined_333,NAMline=="NoPlantTime0")  # (2) Bulk_Time0
# combined_666 <-subset_samples (combined_333,NAMline=="NAM23")  # (1) Rhizo_NAM23 
# combined_666 <-subset_samples (combined_333,NAMline=="NAM37")  # (2) Rhizo_NAM37
combined_666 <- combined_333
View (combined_666)

########
sample_data(combined_666)$P_OW<-paste(sample_data(combined_666)$P_OW)
sample_data(combined_666)[]<-lapply(sample_data(combined_666),factor) #Use lapply to apply the factor() function

## Merges species that have the same taxonomy at a certain taxaonomic rank (e.g.: "Phylum" here)
Sub_Rarefy_Phylum<-tax_glom(combined_666, taxrank = "Species")   ## !!! Changed from ASV to Species cuz updated ANCOMBC!!!!
View (Sub_Rarefy_Phylum)
## Order phyla at abundance from high to low
A<-names(sort(rowSums(otu_table(Sub_Rarefy_Phylum)),decreasing=T))
otu_table(Sub_Rarefy_Phylum)<-otu_table(Sub_Rarefy_Phylum)[intersect(A,rownames(otu_table(Sub_Rarefy_Phylum))),]
tax_table(Sub_Rarefy_Phylum)<-phyloseq::tax_table(Sub_Rarefy_Phylum)[intersect(A,rownames(phyloseq::tax_table(Sub_Rarefy_Phylum))),]
#View (Sub_Rarefy_Phylum)
## Merge sample for plotting [combine taxonomy info and 'Tret']
Sub_Rarefy_merge<-merge_samples(Sub_Rarefy_Phylum,"Tret_POW")
sample_data(Sub_Rarefy_merge)$Tret_POW<-levels(sample_data(Sub_Rarefy_Phylum)$Tret_POW)[get_variable(Sub_Rarefy_merge,"Tret_POW")] #### Need change to Bulk Soil, Rhizosphere Soil
sample_data(Sub_Rarefy_merge)[]<-lapply(sample_data(Sub_Rarefy_merge),factor)
#View (Sub_Rarefy_merge)
# Transform abundance data for plotting
# RA_filter666 <-transform_sample_counts(Sub_Rarefy_merge, function(x) x/sum(x))
RA_filter666 <-transform_sample_counts(Sub_Rarefy_merge, function(x) log10(1 + x))

ASV_Beta666 <- phyloseq::filter_taxa(RA_filter666, function(x) mean(x) > 0, TRUE)
#View (RA_filter666)

##################
## Reorder Tret
# sample_data(ASV_Beta666)$Tret_POW <- factor(sample_data(ASV_Beta666)$Tret_POW,
#                                         levels = c("P0(1)","P35(1)","P65(1)","P0(4)","P35(4)","P65(4)"))

sample_data(ASV_Beta666)$Tret_POW <- factor(sample_data(ASV_Beta666)$Tret_POW,
                                            levels = c("P0(1)","P0(4)","P35(1)","P35(4)","P65(1)","P65(4)"))

# ## Order level
# plot_heatmap(ASV_Beta666,sample.label="Tret_POW",taxa.label="ASV",taxa.order="Order",
#              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Tret_POW,scales="free_x")+labs(fill = "Relative\nabundance")
# plot_heatmap(ASV_Beta666,sample.label="Tret_POW",taxa.label="Order",taxa.order="Order",
#              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Tret_POW,scales="free_x")+labs(fill = "Relative\nabundance")

## Class level
plot_heatmap(ASV_Beta666,sample.label="Tret_POW",taxa.label="Species",taxa.order="Class",
             low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Tret_POW,scales="free_x")+labs(fill = "Log10\nabundance")+
  theme(axis.text.y = element_text(size = 8)) ##+ scale_fill_continuous(limits=c(0, 0.3), breaks=seq(0,0.3,by=0.01))

plot_heatmap(ASV_Beta666,sample.label="Tret_POW",taxa.label="Class",taxa.order="Class",
             low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Tret_POW,scales="free_x")+labs(fill = "Log10\nabundance")+
  theme(axis.text.y = element_text(size = 8))
########################################################################################################################################
# df01 <- read.csv("Fungi-ASV-NoUnclassified.csv")
# df00 <- read.csv("CARP_RootITS_DA_ASVs.csv")
# View (df00)
# df02 <- merge(df00, df01, by = 'ASV')  ## Remove rows have NA
# View (df02)
# write.csv (df02,"CARPITS_RootDA_ASVs_Taxon.csv")








######################## ########### ###########  ###################### 
# # Can combine Year or Stage according to PERMANOVA, mean > 5
# ######### Compare Root between Year_Stage #########
# # Comparison among All 6 Treatments 
# ASV_data = combinedBulk
# View (ASV_data)
# ###### Loop Start ######
# # The taxonomy table
# tax_mat = as(tax_table(ASV_data), "matrix")
# # Run ancombc function on "Treatment (Plevel * OW)"
# out = ancombc(phyloseq = ASV_data, formula = "YS",
#               p_adj_method = "holm", lib_cut = 1000,
#               group = "YS", struc_zero = TRUE, neg_lb = FALSE,
#               tol = 1e-5, max_iter = 100, conserve = TRUE,
#               alpha = 0.05, global = TRUE)
# #ancom = data.frame(out$res,check.names = F)
# ancom = out$res      ## pair comparison 
# res_global = out$res_global     ## multiple comparison
# ###### Loop End #####
# # Can combine Year or Stage according to PERMANOVA, mean > 10
# write.csv(res_global,"/Users/mel270/Desktop/CARPITS/RootITS_YSmean5.csv")
# ancom_root_P0OW <- read.csv("RootITS_YSmean5.csv")
# attach (ancom_root_P0OW)
# ancom_root_P0OW_TRUE<-ancom_root_P0OW[diff_abn=="TRUE",] ## 12 genus different for P0 vs P35, which is the same as P0 vs P65
# write.csv(ancom_root_P0OW_TRUE,"RootITS_YSmean5_TRUE.csv")
# detach (ancom_root_P0OW)
# ######### Compare Root between Year_Stage #########
# 
# 
# ## Plot the ANCOM_BC heatmap
# View (metadata)
# ####### Root differece between Years #######
# t1 <- read.csv("RootITS_YSmean5_TRUE.csv")
# ASV_test<-ASV_1[intersect(t1$ASV,rownames(Taxonomy_1)),]
# View (t1)
# ASV_test<-ASV_1[intersect(t1$X,rownames(Taxonomy_1)),]
# # phyloseq
# OTU0<-otu_table(as.matrix(ASV_test),taxa_are_rows = TRUE)
# TAX0<-tax_table(as.matrix(Taxonomy_1))
# samples0<-sample_data(metadata)
# combined_lol<-phyloseq(OTU0,TAX0,samples0)
# View(combined_lol) 
# combined_666 <- subset_samples(combined_lol, SampleType == "Root")
# combined_666 <- subset_samples(combined_666, Controls == "N")
# View (combined_666)
# 
# ## Merges species that have the same taxonomy at a certain taxaonomic rank (e.g.: "Phylum" here)
# Sub_Rarefy_Phylum<-tax_glom(combined_666, taxrank = "ASV")
# 
# ## Order phyla at abundance from high to low
# A<-names(sort(rowSums(otu_table(Sub_Rarefy_Phylum)),decreasing=T))
# otu_table(Sub_Rarefy_Phylum)<-otu_table(Sub_Rarefy_Phylum)[intersect(A,rownames(otu_table(Sub_Rarefy_Phylum))),]
# tax_table(Sub_Rarefy_Phylum)<-tax_table(Sub_Rarefy_Phylum)[intersect(A,rownames(tax_table(Sub_Rarefy_Phylum))),]
# 
# ## Merge sample for plotting [combine taxonomy info and 'Tret']
# Sub_Rarefy_merge<-merge_samples(Sub_Rarefy_Phylum,"YS")
# sample_data(Sub_Rarefy_merge)$YS<-levels(sample_data(Sub_Rarefy_Phylum)$YS)[get_variable(Sub_Rarefy_merge,"YS")] #### Need change to Bulk Soil, Rhizosphere Soil
# sample_data(Sub_Rarefy_merge)[]<-lapply(sample_data(Sub_Rarefy_merge),factor)
# 
# # Transform abundance data for plotting
# RA_filter666 <-transform_sample_counts(Sub_Rarefy_merge, function(x) x/sum(x))
# ASV_Beta666 <- phyloseq::filter_taxa(RA_filter666, function(x) mean(x) > 0, TRUE)
# View (Sub_Rarefy_merge)
# View (ASV_Beta666)
# sample_data(ASV_Beta666)$YS <- factor(sample_data(ASV_Beta666)$YS,      
#                                       levels = c("2019L", "2019F","2020L","2020F"))
# plot_heatmap(ASV_Beta666,sample.label="YS",taxa.label="ASV",taxa.order="Class",
#              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~YS,scales="free_x")+labs(fill = "Relative\nabundance")
# 
# plot_heatmap(ASV_Beta666,sample.label="YS",taxa.label="Class",taxa.order="Class",
#              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~YS,scales="free_x")+labs(fill = "Relative\nabundance")
# 
# 
# 
# 



# ######################## ########### ###########  ###################### 
# # Can combine Year or Stage according to PERMANOVA, mean > 5
# ######### Compare Root between Year_Stage #########
# # Comparison among All 6 Treatments 
# ASV_data = combinedBulk
# View (ASV_data)
# ###### Loop Start ######
# # The taxonomy table
# tax_mat = as(tax_table(ASV_data), "matrix")
# # Run ancombc function on "Treatment (Plevel * OW)"
# out = ancombc(phyloseq = ASV_data, formula = "YS",
#               p_adj_method = "holm", lib_cut = 1000,
#               group = "YS", struc_zero = TRUE, neg_lb = FALSE,
#               tol = 1e-5, max_iter = 100, conserve = TRUE,
#               alpha = 0.05, global = TRUE)
# #ancom = data.frame(out$res,check.names = F)
# ancom = out$res      ## pair comparison 
# res_global = out$res_global     ## multiple comparison
# ###### Loop End #####
# # Can combine Year or Stage according to PERMANOVA, mean > 10
# write.csv(res_global,"/Users/mel270/Desktop/CARPITS/RootITS_YSmean5.csv")
# ancom_root_P0OW <- read.csv("RootITS_YSmean5.csv")
# attach (ancom_root_P0OW)
# ancom_root_P0OW_TRUE<-ancom_root_P0OW[diff_abn=="TRUE",] ## 12 genus different for P0 vs P35, which is the same as P0 vs P65
# write.csv(ancom_root_P0OW_TRUE,"RootITS_YSmean5_TRUE.csv")
# detach (ancom_root_P0OW)
# ######### Compare Root between Year_Stage #########
# ## Plot the ANCOM_BC heatmap
# View (metadata)
# ####### Root differece between Years #######
# t1 <- read.csv("RootITS_Yearmean5_TRUE.csv")
# ASV_test<-ASV_1[intersect(t1$ASV,rownames(Taxonomy_1)),]
# View (t1)
# ASV_test<-ASV_1[intersect(t1$X,rownames(Taxonomy_1)),]
# # phyloseq
# OTU0<-otu_table(as.matrix(ASV_test),taxa_are_rows = TRUE)
# TAX0<-tax_table(as.matrix(Taxonomy_1))
# samples0<-sample_data(metadata)
# combined_lol<-phyloseq(OTU0,TAX0,samples0)
# View(combined_lol) 
# combined_666 <- subset_samples(combined_lol, SampleType == "Root")
# combined_666 <- subset_samples(combined_666, Controls == "N")
# View (combined_666)
# 
# ## Merges species that have the same taxonomy at a certain taxaonomic rank (e.g.: "Phylum" here)
# Sub_Rarefy_Phylum<-tax_glom(combined_666, taxrank = "ASV")
# 
# ## Order phyla at abundance from high to low
# A<-names(sort(rowSums(otu_table(Sub_Rarefy_Phylum)),decreasing=T))
# otu_table(Sub_Rarefy_Phylum)<-otu_table(Sub_Rarefy_Phylum)[intersect(A,rownames(otu_table(Sub_Rarefy_Phylum))),]
# tax_table(Sub_Rarefy_Phylum)<-tax_table(Sub_Rarefy_Phylum)[intersect(A,rownames(tax_table(Sub_Rarefy_Phylum))),]
# 
# ## Merge sample for plotting [combine taxonomy info and 'Tret']
# Sub_Rarefy_merge<-merge_samples(Sub_Rarefy_Phylum,"Year")
# sample_data(Sub_Rarefy_merge)$Year<-levels(sample_data(Sub_Rarefy_Phylum)$Year)[get_variable(Sub_Rarefy_merge,"Year")] #### Need change to Bulk Soil, Rhizosphere Soil
# sample_data(Sub_Rarefy_merge)[]<-lapply(sample_data(Sub_Rarefy_merge),factor)
# 
# # Transform abundance data for plotting
# RA_filter666 <-transform_sample_counts(Sub_Rarefy_merge, function(x) x/sum(x))
# ASV_Beta666 <- phyloseq::filter_taxa(RA_filter666, function(x) mean(x) > 0, TRUE)
# View (RA_filter666)
# View (ASV_Beta666)
# # plot_heatmap(ASV_Beta666,sample.label="Year",taxa.label="ASV",taxa.order="Class",
# #              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Year,scales="free_x")+labs(fill = "Relative\nabundance")
# # 
# # plot_heatmap(ASV_Beta666,sample.label="Year",taxa.label="Class",taxa.order="Class",
# #              low="#165CAA",high="orange",na.value="#BFC2C5")+facet_grid(~Year,scales="free_x")+labs(fill = "Relative\nabundance")
# 
# plot_bar(ASV_Beta666, "ASV", y="Abundance", fill="Class",
#          title=NULL, facet_grid=~Year)+ theme(panel.background = element_blank())







############################################################################################
############### Differential Abundant ASVs _ Correlation ##############################
########################################################################################
library(corrplot)
View (combined)
sub_root01 <- subset_samples(combined, SampleType == "Root")
View (sub_root)

# sub_root01 <- subset_samples(sub_root, Stage == "L")
# sub_root01 <- subset_samples(sub_root, Stage == "F")
View (sub_root01)

## Change ASV column name to OTU
# DA_ASV <- read.csv("CARP_ITS_DA_Root_L_Tret_TRUE_copy.csv",row.names = 1) ## mean >1
# DA_ASV <- read.csv("CARP_ITS_DA_Root_F_Tret_TRUE_copy.csv",row.names = 1) ## mean >1
DA_ASV <- read.csv("CARP_ITS_DA_Root_Tret_TRUE.csv",row.names = 1) ## mean >1

View (DA_ASV)

OTU = as(otu_table(sub_root01), "matrix")
# if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
root_otu = as.data.frame(OTU)
View (root_otu)
# write.csv (root_otu, "Root_otu_L.csv")
# df_root <- read.csv ("Root_otu_L.csv")
# write.csv (root_otu, "Root_otu_F.csv")
# df_root <- read.csv ("Root_otu_F.csv")  ## Add ASV column name as OTU
# write.csv (root_otu, "Root_otu.csv")
df_root <- read.csv ("Root_otu.csv")  ## Add ASV column name as OTU
View (df_root) 



t11 <- merge(df_root, DA_ASV,by='X', na.rm = TRUE)  ## Remove rows have NA
View (t11)
names (t11)
t55 <- t11[,-1]
rownames(t55) <- t11[,1]  
View (t55)
# t66 <- t55[,-c(49,50,51,52)]   # Separate into different stages
t66 <- t55[,-c(97:100)]   ## Include all samples
View (t66)

r_meta = as(sample_data(sub_root01), "matrix")
# if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
metadf = as.data.frame(r_meta)
View (metadf)
# write.csv (metadf, "root_metadata_L.csv")
# metadf <- read.csv ("root_metadata_L.csv",row.names = 1)
# write.csv (metadf, "root_metadata_F.csv")   ## Clean metadata
# metadf <- read.csv ("root_metadata_F.csv",row.names = 1)

# write.csv (metadf, "root_metadata.csv")  ## Clean metadata
metadf <- read.csv ("root_metadata.csv",row.names = 1)
View (metadf)
## Confirm 'metadata' and 'abundance table' matches
setdiff (colnames(t66),rownames(metadf))  # character(0)
identical(colnames(t66),rownames(metadf)) # TRUE
identical(ncol(t66),nrow(metadf))  # TRUE

Taxonomy_1<-read.csv("Fungi-ASV-NoUnclassified.csv",row.names = 1,header = TRUE)
# grep("Unclassified",Taxonomy_1$Phylum) ## Should be 0
# grep("Chloroplast",Taxonomy_1$Order)  ## Should be 0
# grep("Mitochondria",Taxonomy_1$Family)  ## Should be 0

## Remove these ASVs for chloroplast, mitochondria, unclassified phylum in ASV table
ASV<-ASV_1[intersect(rownames(t66),rownames(Taxonomy_1)),] 

## Format data as phyloseq objects
OTU<-otu_table(as.matrix(ASV),taxa_are_rows = TRUE)
TAX<-phyloseq::tax_table(as.matrix(Taxonomy_1))
samples<-sample_data(metadf)
combined_root_RA <-phyloseq(OTU,TAX,samples)

View (combined_root_RA)
View (Taxonomy_1)

library(MicrobiotaProcess)
data1 <- get_alltaxadf(
  combined_root_RA,
  method = "count",
  type = "species")    ## for relative abundance use "method = NULL"
# write.csv (data1,"RootL_all taxon.csv")
# write.csv (data1,"RootF_all taxon.csv")  
# write.csv (data1,"Root_all taxon.csv")  
View (data1)
## Subset a dataframe to certain taxonomy level, and change taxon names

## Log transform abundance data
# lol <- read.csv ("RootL_taxon_Order.csv",row.names = 1)
# lol <- read.csv ("RootF_taxon_Order.csv",row.names = 1)

# lol <- read.csv ("RootL_taxon_Genus.csv",row.names = 1)
# lol <- read.csv ("RootF_taxon_Genus.csv",row.names = 1)
View (Taxonomy_1)

# lol <- read.csv ("Root_All_Genus.csv",row.names = 1)   ## Taxon order by alphabet 
lol <- read.csv ("Root_taxon_Genus_copy.csv",row.names = 1)  ## Taxon reorder by Class
# haha1 <-lol[intersect(colnames(lol),rownames(Taxonomy_1)),] 
View (haha1)

View (lol)
log_r <- log (lol+1) 
View (log_r)

thaha <- merge(log_r, metadf,by='row.names', na.rm = TRUE)  ## Remove rows have NA
View (thaha)
tha0 <- data.frame(thaha[,-1], row.names=thaha[,1])
View (tha0)
colnames(tha0)
# tha <- tha0[, c(1:31,48:60)]  ## Vegetative stage & Order
# tha <- tha0[, c(1:30,47:59)]  ## Flowering stage & Order
# tha <- tha0[, c(1:66,83:95)]  ## Flowering stage & Genus

tha <- tha0[, c(1:34,51:63)]  ## Flowering stage & Genus

View (tha)
c_df <- Hmisc::rcorr(as.matrix(tha), type='spearman')  ## pearson
View (c_df)
# write.csv (c_df$r,"CARP_RootITS_DA_Correlation_Coefficient.csv")
# write.csv (c_df$P,"CARP_RootITS_DA_Correlation_P-value.csv")

corrplot(corr=c_df$r[35:47,1:34], p.mat=c_df$P[35:47,1:34], sig.level = 0.05,
         addCoef.col=1,  insig='blank', number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

## Sig. level 0.01
corrplot(corr=c_df$r[35:47,1:34], p.mat=c_df$P[35:47,1:34], sig.level=0.01,
         addCoef.col=1,  insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))
## Sig. level 0.001
corrplot(corr=c_df$r[35:47,1:34], p.mat=c_df$P[35:47,1:34], sig.level=0.001,
         addCoef.col=1,  insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

################### Plot method1    ########################
# ## For less variables 
# chart.Correlation(tha, histogram=TRUE, pch=19)

##################### Plot method2  #####################
# cor1 <- cor(tha[,1:45], tha[,46:62], method = c("spearman"))  ## , "kendall", "spearman"
# View (cor1)
# write.csv (cor1,"Root_Order_Cor.csv")
# testRes = cor.mtest(cor1, conf.level = 0.95)
# corrplot(cor1, p.mat = testRes$p, type="upper",
#          addCoef.col = 'grey35', insig='blank',
#          number.cex = 0.7, tl.col="black", tl.srt=45)

#################### Plot method3 #########################################
library (corrplot)
c_df <- Hmisc::rcorr(as.matrix(tha), type='spearman')  ## pearson
View (c_df$r)
# ## (1) Vegetative stage
corrplot(corr=c_df$r[32:44, 1:31], p.mat=c_df$P[32:44, 1:31], sig.level=0.05,
         addCoef.col=1,  insig='blank',
         number.cex=.8, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))
# 
# ## Check sig. levels
# m <- c_df$P[32:44, 1:31] < .05
# m[lower.tri(m, diag=TRUE)] <- ''
# as.data.frame(replace(m, lower.tri(m, diag=TRUE), ''))

## (2) Flowering stage
corrplot(corr=c_df$r[31:43, 1:30], p.mat=c_df$P[31:43, 1:30], sig.level=0.05,
         addCoef.col=1,  insig='blank',
         number.cex=.8, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

## Check sig. levels
m <- c_df$P[31:43, 1:30] < .05
m[lower.tri(m, diag=TRUE)] <- ''
as.data.frame(replace(m, lower.tri(m, diag=TRUE), ''))
##############################################################






#############################################################################
############# Plant Performance -- Soil Property -- Correlation #############
df00 <- read.csv ("/Users/mel270/Desktop/CARPNutrient_sum/CARP_PlantSoil_Correlation.csv")
View (df00)
## Subset samples
df01 <- subset (df00, Year == "2019")
df01 <- subset (df00, Year == "2020")

df02 <- subset (df01, Stage == "L")
df02 <- subset (df01, Stage == "F")
df03 <- df02[,c(9:24)]
View (df03)

## Correlation
library (corrplot)
c_df <- Hmisc::rcorr(as.matrix(df03), type='spearman')  ## pearson
View (c_df$r)
# ## (1) Vegetative stage
corrplot(corr=c_df$r[1:6,7:16], p.mat=c_df$P[1:6,7:16], sig.level=0.01,
         addCoef.col=1,  insig='blank',
         number.cex=.8, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))
#####################################################################











############################################################################################################################
############ Network Analysis ################################################################################################
############################################################################################################################
## Load packages
library (SPRING)
library (SpiecEasi)
library (NetCoMi)
library (limma)
View (combined)
sub001 <- subset_samples(combined, SampleType == "Root")
View (sub001)
sub002 <- filter_taxa(sub001, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
View (sub002)


# sub001 <- subset_samples(combined, SampleType == "Root")
# View (sub001)
# Rarefy_root <- filter_taxa(sub001, function(x) sum(x > 3) > (0.3*length(x)), TRUE)

## Differential abundant ASVs
t1 <- read.csv("CARP_ITS_DA_Root_Tret_TRUE.csv") ## mean >1
View (t1)
ASV_test<-ASV_1[intersect(t1$X,rownames(Taxonomy_1)),]
# phyloseq
OTU0<-otu_table(as.matrix(ASV_test),taxa_are_rows = TRUE)
TAX0<-tax_table(as.matrix(Taxonomy_1))
samples0<-sample_data(metadata)
combined_lol<-phyloseq(OTU0,TAX0,samples0)
# View(combined_lol) 
combined_333 <- subset_samples(combined_lol, SampleType == "Root")
sub002 <- subset_samples(combined_333, Controls == "N")
View (combined_666) ## 67 x 96

sub001 <- subset_samples(combined, SampleType == "Root")

sub002 <- filter_taxa(sub001, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
View (sub002)





## Assign ASV to certain Phylum
taxtab0 <- as(phyloseq::tax_table(sub002), "matrix")
phyla0 <- as.factor(taxtab0[, "Class"])
names(phyla0) <- taxtab0[, "Species"]
# levels (phyla0)
# table (phyla0)

## Subset groups
# soil_warm_yes <- phyloseq::subset_samples(sub002, P_OW == "P65_1")
# soil_warm_no  <- phyloseq::subset_samples(sub002, P_OW == "P0_1")
soil_warm_yes <- phyloseq::subset_samples(sub002, Plevel == "P65")
# soil_warm_yes <- phyloseq::subset_samples(sub002, Plevel == "P35")
soil_warm_no  <- phyloseq::subset_samples(sub002, Plevel == "P0")
View (soil_warm_yes)
View (soil_warm_no)
## p adjust method: lfdr (local false discovery rate correction)
# filtTax = "highestVar",
# filtTaxPar = list(highestVar = 50),
# filtSamp = "totalReads",
# filtSampPar = list(totalReads = 1000),
net_seas_p <- netConstruct(soil_warm_no,soil_warm_yes,
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           zeroMethod = "pseudo",
                           normMethod = "clr",
                           measure = "spearman",
                           verbose = 2,
                           sparsMethod = "threshold", 
                           thresh = 0.3,
                           alpha = 0.01,
                           seed = 123456)

View (net_seas_p)
netprops1 <- netAnalyze(net_seas_p, centrLCC = FALSE,
                        avDissIgnoreInf = TRUE,
                        sPathNorm = FALSE,
                        clustMethod = "cluster_fast_greedy",
                        hubPar = c("degree", "eigenvector"),
                        hubQuant = 0.9,
                        lnormFit = TRUE,
                        normDeg = FALSE,
                        normBetw = FALSE,
                        normClose = FALSE,
                        normEigen = FALSE)

# summary (netprops1)
phylcol0 <- c( "#FF6666", "dodgerblue2", "#CAB2D6",  "#6A3D9A","#FB9A99", "#FDBF6F", 
               "orchid1", "maroon", "blue1", "steelblue4","darkturquoise", 
               "#FFCC00","#999999","#99FFCC","yellow4", "yellow3", "darkorange4", "brown",
               "#a6cee3","#E1E7E9","#1f78b4","#ff7f00","#de77ae","#EEE8AA","#5F9EA0","green4","green1","deeppink1")
## "#a6cee3","#E1E7E9","#1f78b4","#ff7f00","#de77ae","#EEE8AA","#5F9EA0",

## Bacteria color legends
# phylcol0 <- c( "#a6cee3","#E1E7E9","#1f78b4","#ff7f00","#de77ae","#EEE8AA","#5F9EA0")

# View (netprops1) 
## nodeFilterPar = 50  -- Only plot 50 nodes with highest degree
plot(netprops1, 
     sameLayout = TRUE, 
     repulsion = 0.95,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "clr", 
     nodeColor = "feature",
     colorVec = phylcol0,
     featVecCol = phyla0, 
     labelScale = FALSE,
     cexNodes = 1, 
     cexLabels = 0.8,
     cexHubLabels = 0.8,
     cexTitle = 1.5,
     groupNames = c("P0", "P65"),
     hubBorderCol  = "gray")
## c("P0", "P65")

# print (names(featVecCol))
# ls ()
# Colors used in the legend should be equally transparent as in the plot
phylcol_transp <- colToTransp(phylcol0, 60)
legend(-0.13, 1.2, cex = 0.8, pt.cex = 0.8, title = "Phylum:", 
       legend=levels(phyla0), col = phylcol_transp, bty = "n", pch = 16)

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.01, cex = 0.5, lty = 0.5, lwd = 0.5, 
       bty = "n", horiz = TRUE)


## Quantitative network comparison
comp_season <- netCompare(netprops1, 
                          permTest = TRUE, 
                          verbose = FALSE,
                          seed = 123456,
                          testRand = TRUE)
summary(comp_season, pAdjust = TRUE,
        groupNames = c("P35", "P65"),
        showCentr = c("all"), 
        numbNodes = 10)
## c("P0", "P65")
## "degree", "between", "closeness"
#############################################################################################################################














############################################################### ######
###### Relative abundace comparison between P treatments
############################################################### ######
View(combined) 
View (sub01)
## Subset root samples
sub_root <- subset_samples(sub01, SampleType == "Rhizo")
sub_root <- subset_samples(sub01, SampleType == "Root")

View (sub_root)
sub_root01 <- sub_root
View (sub_root01)
cyanos_glom <- tax_glom(sub_root01, taxrank = "Class")
Tax <-names(sort(rowSums(otu_table(cyanos_glom)),decreasing=T))
otu_table(cyanos_glom)<-otu_table(cyanos_glom)[intersect(Tax,rownames(otu_table(cyanos_glom))),]
tax_table(cyanos_glom)<-phyloseq::tax_table(cyanos_glom)[intersect(Tax,rownames(phyloseq::tax_table(cyanos_glom))),]
# View (cyanos_glom)
TT <- transform_sample_counts(cyanos_glom, function(x) log(x+1) ) ## relative abundance
phylum.sum = tapply(taxa_sums(TT), phyloseq::tax_table(TT)[, "Class"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:12]
GP1 = prune_taxa((phyloseq::tax_table(TT)[, "Class"] %in% top5phyla),TT) ## 67 x 96
View (GP1)

cyanos_df <- psmelt(GP1)   ## melting phyloseq to dataframe
cyanos_df$Class <- as.factor(cyanos_df$Class)
# cyanos_df$Species <- as.factor(cyanos_df$Species)
View (cyanos_df)
# write.csv (cyanos_df, "Bulk_RA_Pairwise.csv")
levels (cyanos_df$Class)
unique (cyanos_df$Phylum)

## Fungal color legend
# library(RColorBrewer)
# display.brewer.pal(n = 12, name = 'Paired') ## Display color
# brewer.pal(n = 12, name = "Paired")  ## Get color names
# ## "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
# ## "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

# ## "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
# ## "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"
speciesPalette <- c(Ascomycota = "#A6CEE3", Basidiomycota = "#1F78B4", Chytridiomycota = "#B2DF8A", Glomeromycota = "#33A02C",
                    Olpidiomycota = "#FB9A99", Mortierellomycota= "#E31A1C",  Rozellomycota = "#FDBF6F",
                    Mucoromycota = "#FF7F00", Basidiobolomycota = "#CAB2D6", Calcarisporiellomycota = "#FFFF99", Entomophthoromycota = "#6A3D9A", Kickxellomycota = "#B15928", 
                    Unclassified = "#6A3D9A", Zoopagomycota = "#6A3D9A",Entorrhizomycota = "#6A3D9A")
# ## "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
# ## "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"


# ############ Fungi
# speciesPalette <- c(Ascomycota = "#FF7F00", Basidiomycota = "#B15928", Mortierellomycota = "#B2DF8A", Olpidiomycota = "#33A02C", Unclassified = "#6A3D9A",
#                     Chytridiomycota = "#6A3D9A", Zoopagomycota = "#6A3D9A", Rozellomycota = "#6A3D9A",Glomeromycota = "#6A3D9A", Mucoromycota = "#6A3D9A",
#                     Basidiobolomycota = "#6A3D9A", Calcarisporiellomycota = "#6A3D9A", Entomophthoromycota = "#6A3D9A", Kickxellomycota = "#6A3D9A", Entorrhizomycota = "#6A3D9A")
# phyloseq::psmelt(GP1) %>%
ggplot(data = cyanos_df, aes(x=factor(Tret, level=c('NP_LR', 'P_LR', 'NP_PR', 'P_PR', 'NP_WR', 'P_WR')), y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Fungal Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+ 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor = element_blank(),panel.background = element_blank())
##+scale_colour_manual( values = speciesPalette)

## Compare RootType only if no interaction was detected    [Used for Fungi]
ggplot(data = cyanos_df, aes(x = RootType, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Fungal Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())
##+scale_colour_manual( values = speciesPalette)

##### Compare SoilP only 
ggplot(data = cyanos_df, aes(x = SoilP, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Fungal Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank())



#### Fungi - Pairwise comparison among the Tret
levels (cyanos_df$Class)
# [1] "Atractiellomycetes"  "Cystobasidiomycetes"
# [3] "Dothideomycetes"     "Eurotiomycetes"     
# [5] "Laboulbeniomycetes"  "Leotiomycetes"      
# [7] "Microbotryomycetes"  "Mortierellomycetes" 
# [9] "Olpidiomycetes"      "Orbiliomycetes"     
# [11] "Sordariomycetes"     "Tremellomycetes"
ta1 <- subset (cyanos_df, Class=="Agaricomycetes")   ## Basidiomycota
ta1 <- subset (cyanos_df, Class=="Atractiellomycetes")   ## Basidiomycota
# ta1 <- subset (cyanos_df, Class=="Cystobasidiomycetes")   ## Basidiomycota
ta1 <- subset (cyanos_df, Class=="Dothideomycetes")    ## Ascomycota
ta1 <- subset (cyanos_df, Class=="Eurotiomycetes")   ## Ascomycota
# ta1 <- subset (cyanos_df, Class=="Laboulbeniomycetes")   ## Ascomycota
ta1 <- subset (cyanos_df, Class=="Leotiomycetes")   ## Ascomycota
ta1 <- subset (cyanos_df, Class=="Microbotryomycetes")   ## Basidiomycota
ta1 <- subset (cyanos_df, Class=="Mortierellomycetes")   ## Mortierellomycota
ta1 <- subset (cyanos_df, Class=="Olpidiomycetes")   ## Olpidiomycota
ta1 <- subset (cyanos_df, Class=="Orbiliomycetes")   ## Ascomycota
# ta1 <- subset (cyanos_df, Class=="Rhizophlyctidomycetes")   ## Olpidiomycota
ta1 <- subset (cyanos_df, Class=="Saccharomycetes")   ## Ascomycota
ta1 <- subset (cyanos_df, Class=="Sordariomycetes")   ## Ascomycota
ta1 <- subset (cyanos_df, Class=="Tremellomycetes")   ## Basidiomycota
# ta1 <- subset (cyanos_df, Class=="Unclassified")   ## Unclassified

ta1
#####
model<-lme(Abundance~SoilP*RootType,random=~1|Rep,data=ta1)
anova (model)
# write.csv (anova (model),"ANOVA_LogRA_Tret_Orbiliomycetes_RootClass_6wk.csv")


##### Separate RootTypes to check "soil P" effect ##### 
ta2 <- subset (ta1, RootType =="Lateral")
ta2 <- subset (ta1, RootType =="Primary")
ta2 <- subset (ta1, RootType =="Whole")
model<-lme(Abundance~SoilP,random=~1|Rep,data=ta2)
anova (model)
# write.csv (anova (model),"ANOVA_LogRA_Tret_Eurotiomycetes_RhizoClass_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret__RootClass_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidobacteriae_RootClass_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Saccharimonadia_RootClass_6wk.csv")

model<-lme(Abundance~RootType,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"RootType")
pair<-pairs(em,adjust = "tukey")
# write.csv (pair,"ANOVA_LogRA_Tret_Saccharimonadia_RootClass_6wk_Pairwise.csv")

model<-lme(Abundance~Tret,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
# write.csv (pair,"ANOVA_LogRA_Tret_Orbiliomycetes_RootClass_6wk_Pairwise.csv")


model<-lme(Abundance~SoilP,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"SoilP")
pair<-pairs(em,adjust = "tukey")