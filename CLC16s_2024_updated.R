################ Raw data Processing ################
##########$$$ Make manifest for Local Desk $$$#########
# Set the path to where your sequence files are stored
## Note:when making the absolte-filepath, should be the same with this path!!!!!
path <- "/Users/mel270/Downloads/Test_Qiime2/raw_data/"  
# Create a dataframe of the file locations
manifest<- data.frame(list.files(path), stringsAsFactors = F)
# Extract the sample names
sample_id <- sapply(strsplit(manifest$list.files.path., "_"), "[", 1)
# Make the file paths
absolute_filepath <- paste("/u2/mel270/TestQiime2_16s/raw_data/", manifest$list.files.path., sep = "")       
# Make the read directions (when there's 440 rows, direction should be 220)
direction <- rep(c("forward","reverse"), 3)           
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
write.csv(manifest_ITS, "/Users/mel270/Downloads/manifest_TestQiime2.csv", row.names = FALSE,quote=FALSE)
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
library(nloptr)
library(dplyr)
library(pldist)
library(ANCOMBC)

########## Load and manipulate data ##########
setwd("/Users/mel270/Downloads/Bobbi_DNAsequencing/CLC_BH/CLC16s-2024/CLC-16s-AnalysisR/DADA2-Output")
#setwd("/Users/apple/Downloads/CLC16s/")
#feature<-read.table("ITS_feature-table-.tsv",sep='\t',row.names = 1, header = TRUE) # remove the "#" before "OTU"
feature<-read.table("feature-table_updated.tsv",sep='\t',row.names = 1, header = TRUE) # remove the "#" before "OTU"
#View(feature)
metadata<-read.table("CLC-16s-metadata_updated.tsv",sep='\t',row.names = 1,header = TRUE) # change "-" to "."
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
Ta2<-sapply(Ta[,1:7],function(x)gsub("[a-z]__","",as.character(x))) 
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
## Remove chloroplast, mitochondria,and unclassified phylum
Taxonomy_u<-Taxonomy_cm[!(Taxonomy_cm$Phylum %in% c("Unclassified")),] #Remove "Unclassified" at phylum level.
#View(Taxonomy_u)
Taxonomy_0 <- Taxonomy_u[!grepl("Mitochondria", Taxonomy_u$Genus),]
#View (Taxonomy_0)
Taxonomy_1 <- Taxonomy_0[!grepl("Chloroplast", Taxonomy_0$Genus),]
# View(Taxonomy_1)
Taxonomy_1$Phylum <- factor (Taxonomy_1$Phylum)
levels(Taxonomy_1$Phylum)
# write.csv(Taxonomy_1,"Bacteria-allASV-clean.csv")### Change ASV names in Excel
# Taxonomy_1<-read.csv("Bacteria-allASV-deleteD0.csv",row.names = 1,header = TRUE)
Taxonomy_1<-read.csv("Bacteria-allASV-clean.csv",row.names = 1,header = TRUE)
Taxonomy_1<-read.csv("Bacteria-allASV-DA.csv",row.names = 1,header = TRUE)  ### Exchange ASV and Species name for DA analysis
# View(Taxonomy_1)
levels (factor(Taxonomy_1$Class))
## Confirm chloroplast, Mitochondria and unclassified phylum have been removed at certain level, both results should be 0.
grep("Unclassified",Taxonomy_1$Phylum) ## Should be 0
## grep("Chloroplast|Mitochondria",Taxonomy_1$Genus) ########## **********Remove from Genus level********* 
grep("Chloroplast",Taxonomy_1$Order)  ## Should be 0
grep("Mitochondria",Taxonomy_1$Family)  ## Should be 0

## Remove these ASVs for chloroplast, mitochondria, unclassified phylum in ASV table
ASV<-ASV_1[intersect(rownames(ASV_1),rownames(Taxonomy_1)),] # The elements of intersect(x,y) are those elements in x and in y
#View (ASV)
#list (rownames(ASV_1))
#list (rownames(Taxonomy_1))
#############################################
## Format data as phyloseq objects
OTU<-otu_table(as.matrix(ASV),taxa_are_rows = TRUE)
TAX<-phyloseq::tax_table(as.matrix(Taxonomy_1))
samples<-sample_data(metadata)
combined_plusctrl<-phyloseq(OTU,TAX,samples)
## Remove all control samples ### Note 
combined <- subset_samples(combined_plusctrl, Ctrl!= "Y")
View (combined)
sub01 <- subset_samples(combined, Stage!= "4wk")
View (sub01) ##8021 otu in 96 samples
#############################################



######## Rarefaction curve ########
## Curve removed contrls
curve<-ggrare(sub01, step = 250, se = FALSE)
## Curve doesn't remove contrls: curve<-ggrare(combined_plusctrl, step = 250, se = FALSE)

## SampleType
p<-curve+labs(title="Rarefaction_16S")
p<-p+geom_line(aes(col=SampleType))+
  scale_color_manual(labels = c("Rhizo","Root"),
                     values = c("orange","darkgreen"))
p
ggsave("CLC_Rarefaction_16s_SampleType_6wk_clean.png",width = 5,height=2.5)
########################################################################










########################################################################
################## Relative abundance ###############################
########################################################################
View (combined)
## Resample an OTU table such that all samples have the same library size.
Sub_Rarefy<-rarefy_even_depth(sub01,sample.size=2500,rngseed=999) 

#########################################################
## Relative abundance in each Year at Class Level ###
#########################################################
## Resample an OTU table such that all samples have the same library size.
# library(pals)
# pal.bands(polychrome,show.names=TRUE)
View(Sub_Rarefy)
sample_data(Sub_Rarefy)$STPO<-paste(sample_data(Sub_Rarefy)$SampleType,sample_data(Sub_Rarefy)$Tret,sep="-")
sample_data(Sub_Rarefy)[]<-lapply(sample_data(Sub_Rarefy),factor) #Use lapply to apply the factor() function 
Sub_Rarefy_Class<-tax_glom(Sub_Rarefy, taxrank = "Phylum") # Merges species that have the same taxonomy at a certain taxaonomic rank (e.g.: "Family" here)
#View(Sub_Rarefy_Class) #43x284
## Order Class from high to low
A<-names(sort(rowSums(otu_table(Sub_Rarefy_Class)),decreasing=T))
otu_table(Sub_Rarefy_Class)<-otu_table(Sub_Rarefy_Class)[intersect(A,rownames(otu_table(Sub_Rarefy_Class))),]
tax_table(Sub_Rarefy_Class)<-tax_table(Sub_Rarefy_Class)[intersect(A,rownames(tax_table(Sub_Rarefy_Class))),]
View(A)
View(otu_table(Sub_Rarefy_Class))
View(tax_table(Sub_Rarefy_Class)) ## Check the number : 40

## Replace Class with low reads with "others" !!!!!!!!!!!!!!!!!!!!! low reads number?? how to define? 
# tax_table(Sub_Rarefy_Class)[,2][26:27]<-rep("Others",2) ###### [,number of level][exhibit taxon:total taxon]######
# Sub_Rarefy_Class<-tax_glom(Sub_Rarefy_Class, taxrank = "Phylum") # combine all the Family belonging to "others"
# sample_data(Sub_Rarefy_Class)[]<-lapply(sample_data(Sub_Rarefy_Class),factor) 

Sub_Rarefy_merge<-merge_samples(Sub_Rarefy_Class,"STPO")
sample_data(Sub_Rarefy_merge)$SampleType<-levels(sample_data(Sub_Rarefy_Class)$SampleType)[get_variable(Sub_Rarefy_merge,"SampleType")]
sample_data(Sub_Rarefy_merge)$Tret<-levels(sample_data(Sub_Rarefy_Class)$Tret)[get_variable(Sub_Rarefy_merge,"Tret")]
sample_data(Sub_Rarefy_merge)[]<-lapply(sample_data(Sub_Rarefy_merge),factor) #this step is very important
# colourCount = length(unique(tax_table(Sub_Rarefy_merge)[,3])) ## number of levels (e.g.: Phylum is 2, Class is 3)
# getPalette = pal.bands(polychrome,show.names=FALSE) #Determine colors used in a visualization.
# # otu_table(Sub_Rarefy_merge)<-otu_table(Sub_Rarefy_merge)/rowSums(otu_table(Sub_Rarefy_merge)) #Mean tax for samples
# getPalette <- c("#a6cee3","#E1E7E9","#1f78b4",
#                 "#ff7f00","#de77ae", "#FFCC00","#999999","#EEE8AA","#5F9EA0","#99FFCC",
#                 "#FF6666", "dodgerblue2", "#CAB2D6", "green4", "#FB9A99", "deeppink1", "#FDBF6F", 
#                 "orchid1", "maroon", "#6A3D9A", "blue1", "steelblue4","darkturquoise", "green1", 
#                 "yellow4", "yellow3", "darkorange4", "brown")

getPalette <- c("#ff7f00", "#EEE8AA", "#1f78b4", "#FFCC00","#de77ae",
                "#a6cee3","#999999","#5F9EA0","#99FFCC","#FF6666",
                "#9ACD32","darkblue","#006400","#00FA9A","chocolate4", 
                "#20B2AA","#6a3d9a","#E1E7E9","#9999FF" , "#3399FF","pink", "blue1", "steelblue4","darkturquoise", "green1", 
                "yellow4", "yellow3", "darkorange4", "brown","red","black")




#change the name of level which will be showed in figure
#transform the raw counts of reads to proportions within a sample
otu_table(Sub_Rarefy_merge)<-otu_table(Sub_Rarefy_merge)/rowSums(otu_table(Sub_Rarefy_merge)) #Mean tax for samples
#should only run one time when plotting!
SS<-levels(sample_data(Sub_Rarefy_merge)$SampleType)
SS<-gsub("Rhizo","Rhizosphere",SS)
levels(sample_data(Sub_Rarefy_merge)$SampleType)<-SS

# Reorder stack barplot with: fill = fct_reorder(Genus, Abundance)    
p<-plot_bar(Sub_Rarefy_merge,"Tret", fill = "Phylum")+geom_bar(aes(fill=Phylum), stat="identity", position="stack")+scale_fill_manual(values = getPalette)+scale_y_continuous(labels=sprintf("%1.0f",seq(0,100,25)))
p+facet_grid(Stage~SampleType)+labs(y="Relative Abundance (%)",x="Tret")+guides(fill = guide_legend(nrow = 15))
########################################################################################################################















################################################################################
######### Alpha_diversity plot ######### Need to use rarefied data!!!!!! ######
################################################################################################################################################################################
################ Alpha Diversity-Boxplot with dots + Statistic Analysis ################################################################
################################################################################

## All alpha diversity index
View (sub01) ## 71 x 13
Sub_Rarefy<-rarefy_even_depth(sub01,sample.size=3000,rngseed=999) 
View(Sub_Rarefy) ##
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy,SampleType=="Rhizo")

Sub_Rarefy<-rarefy_even_depth(sub01,sample.size=1500,rngseed=999) 
Sub_Rarefy_Bulk<-subset_samples(Sub_Rarefy,SampleType=="Root")

# Sub_Rarefy_RR<-subset_samples(Sub_Rarefy_Bulk,Stage=="4wk")
# Sub_Rarefy_RR<-subset_samples(Sub_Rarefy_Bulk,Stage=="6wk")


Sub_Rarefy_RR <- Sub_Rarefy_Bulk
View(Sub_Rarefy_RR)
## Calculate all alpha diversity indices
bac_alpha<-microbiome::alpha(Sub_Rarefy_RR,index="all")
# write.csv(bac_alpha,"CLC_16s_Alpha_Rhizo_6wk.csv")
# write.csv(bac_alpha,"CLC_16s_Alpha_Root_6wk.csv")
# (1) Alpha - Observed richness ###############
observed<-subset(bac_alpha,select="observed")
mydata<-data.frame(sample_data(Sub_Rarefy_RR),observed)
# View (mydata)

## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=observed, fill=RootType)) + coord_cartesian(ylim = c(0, 200))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Observed richness")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))

## Rhizo 500

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


# Post Hoc Test- pairwase comparison for all sample types
model<-lme(observed~Tret,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_PostHoc_6wk.csv")


# (2) Alpha - Shannon diversity ###############
diversity_shannon<-subset(bac_alpha,select="diversity_shannon")
mydata<-data.frame(sample_data(Sub_Rarefy_RR),diversity_shannon)
## Reorder for samples

## Plot diversity
mydata %>%
  ggplot(aes(x=SoilP, y=diversity_shannon, fill=RootType)) + coord_cartesian(ylim = c(2, 5))+ 
  geom_boxplot(outlier.shape = 20, alpha=0.5, width=0.6)+ 
  theme_bw() +labs(x="",y="Shannon diversity")+ theme(axis.text.x = element_text(angle = 90),
                                                      text = element_text(size = 12))

## Rhizo 6
model<-lme(diversity_shannon~SoilP*RootType,random=~1|Rep,data=mydata)
anova(model)
# qqnorm(residuals(model))
# plot(residuals(model))
# write.csv(anova(model),file="Root_Alpha_Redo_Mar04/shannon_anova_Root2020F.csv")

# Post Hoc Test- pairwase comparison for all sample types
model<-lme(diversity_shannon~Tret,random=~1|Rep,data=mydata)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
em<-emmeans(model,pairwise ~ Tret,adjust = "tukey") ##Estimated marginal means
em
write.csv(pair,file="Root_Alpha_shannon_PostHoc_Root_6wk.csv")

model<-lme(diversity_shannon~Plevel,random=~1|Rep,data=mydata)
anova(model)
em<-emmeans(model,"Plevel")
pair<-pairs(em,adjust = "tukey")
# write.csv(pair,file="Root_Alpha_Redo_Mar04/shannon_PostHoc_Root2019L.csv")

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










###########################################################################
############### PCA  #############################################
###########################################################################
View(combined) 
View (sub01)
ASV_Beta<-phyloseq::filter_taxa(sub01,function(x) mean(x)>1,TRUE)
#data transformation
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta)),method = "BL", output = "p-counts") ### No. corrected values: 853
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE) ## "clr": gives the centered log ratio transform; clrInv gives closed compositions with the given clr-transform
meta<-data.frame(sample_data(ASV_Beta))
taxa<-data.frame(tax_table(ASV_Beta))
#generate a new phyloseq with transformed data
OTU_clr<-otu_table(as.matrix(t(ASV_clr)),taxa_are_rows = TRUE) #format data as phyloseq objects
View (OTU_clr) 
TAX1<-tax_table(as.matrix(taxa))
samples1<-sample_data(meta)
combined_clr<-phyloseq(OTU_clr,TAX1,samples1)
View(combined_clr) 
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
### Rhizo & Root  
# rhizo2019_clr<-subset_samples(rhizo_clr,Stage=="4wk")
# rhizo2019_clr<-subset_samples(rhizo_clr,Stage=="6wk")
rhizo2019_clr <- rhizo_clr
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
View (sub01)
# ## Include all samples, Remove OTUs with a mean read count across all samples less than or equal to 1.
ASV_Beta<-filter_taxa(sub01,function(x) mean(x)>1,TRUE)
# data transformation -- clr-transform
ASV_ZR<-cmultRepl(t(otu_table(ASV_Beta)),method = "SQ", output = "p-counts") # "cmultRepl": Bayesian-Multiplicative Replacement Of Count Zeros
ASV_clr <- codaSeq.clr(ASV_ZR, samples.by.row = TRUE)
sample_data(ASV_Beta)[]<-lapply(sample_data(ASV_Beta),factor)
metadata_Beta<-data.frame(sample_data(ASV_Beta))

#set Blocks for PREMERMOVA
set.seed(2021)
perm <- how(nperm = 999)
setBlocks(perm) <- with (metadata_Beta,Rep)
# check the significance of Stage, Year, and SampleType on bacterial community structure
fit1<-adonis2(ASV_clr~SampleType*SoilP*RootType,data=metadata_Beta,method = "euclidean",by="terms",permutations = 999)
write.csv(fit1,file="CLC-16s-PERMANO-All_6wk.csv")



##Rhizosphere Soil + Root Beta Diversity #########
##  Remove OTUs with a mean read count across all samples less than or equal to 1.
#generate Rhizosphere Soi phyloseq

sub_rhizo<-subset_samples(sub01,SampleType=="Rhizo") 
sub_rhizo<-subset_samples(sub01,SampleType=="Root") 
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
# write.csv(fit1,file="CLC-16s_Perman_Rhizo_6wk.csv")
write.csv(fit1,file="CLC-16s_Perman_Root_6wk.csv")

# fit1<-adonis2(ASV_clr~Stage*SoilP*RootType,data=metadata_Beta_rhizo,method = "euclidean",by="terms",permutations = perm)
# write.csv(fit1,file="CLC-16s_Perman_Root.csv")
########################################################

library(devtools) 
library(cluster)
library(pairwiseAdonis) # Citation:Martinez Arbizu, P. (2020). pairwiseAdonis: Pairwise multilevel comparison using adonis. R package version 0.4
fit1<-adonis2(ASV_clr~SoilP*RootType,data=metadata_Beta_rhizo,method = "euclidean",by="terms",permutations = perm)
pair <- pairwise.adonis(ASV_clr,metadata_Beta_rhizo$Tret,sim.function='vegdist',sim.method='euclidian',p.adjust.m='holm')
write.csv (pair,"Pairwise PERMANO_Root_6wk.csv")
#############################################################################



############## Visualize Prevalence Distributions for Taxa ############################
library (microbiome)
p <- plot_taxa_prevalence(sub0, 'Phylum', detection = 5)
print(p)
####################################################################################
######################### Differential Abundance Analysis ###########################
############# ANCOM-BC ############# ############# ############# #############
View (combined)
View (sub01)
## Subset Root samples
sub0 <-subset_samples(sub01,SampleType=="Rhizo")
sub0 <-subset_samples(sub01,SampleType=="Root")
View (sub0)
## Filter according to mean value
# sub0 <- phyloseq::filter_taxa(sub,function(x) mean(x)>4,TRUE)
## This function will remove taxa with low prevalence (e.g.: 5% of samples), where prevalence is the fraction of total samples in which an OTU is observed. 
## https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_filter_prevalence.html#google_vignette
# install.packages("remotes")
# remotes::install_github("vmikk/metagMisc")
# library ("metagMisc")
# sub1 <- phyloseq_filter_prevalence(
#   sub0,
#   prev.trh = 0.05,
#   abund.trh = NULL,
#   threshold_condition = "OR",
#   abund.type = "total")

###### Filter rare taxa: (1) either prevalence of read count > 10 in more than 50% of samples, OR mean relative abundance of >0.1% and prevalence of read count > 10 in more than 10% of samples.
# dat_pr_high = filter_taxa(dat_pr_OBIT,  function(x){(sum(x > 10) > nsamples*0.5) | ((sum(x > 10) > (nsamples*0.1)) & (mean(x/sample_sums(dat_pr_OBIT)) > 0.001))}, prune = T)
## (2) Remove taxa not seen more than 3 times in at least 20% of the samples.
# filter_taxa(GlobalPatterns, function(x) sum(x > 3) > (0.2*length(x)), TRUE)



# sub0 <-subset_samples(sub,SampleType=="Root")
# sub1 <-subset_samples(sub0,Stage=="L")
# sub1 <-subset_samples(sub0,Stage=="F")

##### Filter 20% prevalence
# sub1 <- filter_taxa(sub0, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

sub1 <- sub0
View (sub1)
## Include all samples !!!!!!!!
# tax_mat = as(phyloseq::tax_table(sub1), "matrix")
# out = ancombc(phyloseq = sub1, formula = "Tret",
#               p_adj_method = "holm",tax_level = "Species",
#               group = "Tret", struc_zero = TRUE, neg_lb = FALSE,
#               tol = 1e-5, max_iter = 100, conserve = TRUE,
#               alpha = 0.001, global = TRUE)
# ## ANCOMBC - Results
# ancom = out$res
# res_global = out$res_global
# # write.csv(res_global,"CLC_16S_DA_Rhizo_Tret_6wk_3_20%.csv")
# # write.csv(res_global,"CLC_16S_DA_Root_Tret_6wk_3_20%.csv")


########## Compare RootTypes ##########
tax_mat = as(phyloseq::tax_table(sub1), "matrix")
out = ancombc(phyloseq = sub1, formula = "RootType",
              p_adj_method = "holm",tax_level = "Species",
              group = "RootType", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.001, global = TRUE)
ancom = out$res
res_global = out$res_global
# write.csv(res_global,"CLC_16S_DA_RootType_Rhizo_6wk.csv")
# write.csv(res_global,"CLC_16S_DA_RootType_Root_6wk.csv")


##### Compare SoilP
tax_mat = as(phyloseq::tax_table(sub1), "matrix")
out = ancombc(phyloseq = sub1, formula = "SoilP",
              p_adj_method = "holm",tax_level = "Species",
              group = "SoilP", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.001, global = TRUE)
ancom = out$res     ## pairwise comparison
# res_global = out$res_global  ## multiple comparison
res_global = out$res_pair    
# write.csv(ancom,"CLC_16S_DA_Rhizo_SoilP_6wk.csv")
# write.csv(ancom,"CLC_16S_DA_Root_SoilP_6wk.csv")




#####################
## Extract the "TRUE" ASVs
# ancom <- read.csv("CLC_16S_DA_RootType_Rhizo_6wk.csv")
# ancom <- read.csv("CLC_16S_DA_Rhizo_SoilP_6wk.csv")
ancom <- read.csv("CLC_16S_DA_RootType_Root_6wk.csv")
ancom <- read.csv("CLC_16S_DA_Root_SoilP_6wk.csv")

attach (ancom)
#View (ancom)
ancom_TRUE<-ancom[ancom$diff_abn=="TRUE",] ## 12 genus different for P0 vs P35, which is the same as P0 vs P65
##### Compare SoilP  
ancom_TRUE<-ancom[ancom$diff_abn.SoilPP=="TRUE",] ## 12 genus different for P0 vs P35, which is the same as P0 vs P65
# write.csv(ancom_TRUE,"CLC_16S_DA_RootType_Rhizo_6wk_TRUE.csv") diff_abn.SoilPP
# write.csv(ancom_TRUE,"CLC_16S_DA_Rhizo_SoilP_6wk_TRUE.csv")
# write.csv(ancom_TRUE,"CLC_16S_DA_RootType_Root_6wk_TRUE.csv")
# write.csv(ancom_TRUE,"CLC_16S_DA_Root_SoilP_6wk_TRUE.csv")
detach (ancom)


####### Construct differential table that includes all taxon levels ##########################
df01 <- read.csv("Bacteria-allASV-DA.csv")
View (df01)

# df00 <- read.csv("CLC_16S_DA_RootType_Rhizo_6wk_TRUE.csv")
# df00 <- read.csv("CLC_16S_DA_Rhizo_SoilP_6wk_TRUE.csv")
# df00 <- read.csv("CLC_16S_DA_RootType_Root_6wk_TRUE.csv")
df00 <- read.csv("CLC_16S_DA_Root_SoilP_6wk_TRUE.csv")

View (df00)
df02 <- merge(df00, df01, by = 'X')  ## Remove rows have NA
View (df02)

# write.csv (df02,"CLC_16S_DA_RootType_Rhizo_6wk_TRUE_AllTaxon.csv")
# write.csv (df02,"CLC_16S_DA_Rhizo_SoilP_6wk_TRUE_AllTaxon.csv")
# write.csv (df02,"CLC_16S_DA_RootType_Root_6wk_TRUE_AllTaxon.csv")
write.csv (df02,"CLC_16S_DA_Root_SoilP_6wk_TRUE_AllTaxon.csv")
###################################################################








############################################################################################
############### Differential Abundant ASVs _ Correlation ##############################
########################################################################################
View (combined)
View (sub01)
sub_root01 <- subset_samples(sub01, SampleType == "Rhizo")
sub_root01 <- subset_samples(sub01, SampleType == "Root")

View (sub_root01)

# sub_root01 <- subset_samples(sub_root, Stage == "L")
# sub_root01 <- subset_samples(sub_root, Stage == "F")
# View (sub_root01)

## Change ASV column name to OTU
DA_ASV <- read.csv("CLC_16S_DA_Rhizo_Tret_6wk_TRUE.csv",row.names = 1) ## mean >1
DA_ASV <- read.csv("CLC_16S_DA_Root_Tret_6wk_TRUE.csv",row.names = 1) ## mean >1
View (DA_ASV)

OTU = as(otu_table(sub_root01), "matrix")
# if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
root_otu = as.data.frame(OTU)
View (root_otu)
# write.csv (root_otu, "Rhizo_otu_Tret_6wk.csv")
# write.csv (root_otu, "Root_otu_Tret_6wk.csv")

# df_root <- read.csv ("Root_otu_L.csv")
# write.csv (root_otu, "Root_otu_F.csv")
# df_root <- read.csv ("Root_otu_F.csv")  ## Add ASV column name as OTU
# write.csv (root_otu, "Root_otu.csv")
df_root <- read.csv ("Rhizo_otu_Tret_6wk.csv")
df_root <- read.csv ("Root_otu_Tret_6wk.csv")
View (df_root) 

t11 <- merge(df_root, DA_ASV,by='X', na.rm = TRUE)  ## Remove rows have NA
View (t11)
names (t11)
t55 <- t11[,-1]
rownames(t55) <- t11[,1]  
View (t55)
# t66 <- t55[,-c(49,50,51,52)]   # Separate into different stages
# t66 <- t55[,-c(37:40)]   ## Rhizo
t66 <- t55[,-c(36:39)]   ## Root

View (t66)

r_meta = as(sample_data(sub_root01), "matrix")
# if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
metadf = as.data.frame(r_meta)
View (metadf)
# write.csv (metadf, "root_metadata_L.csv")
# metadf <- read.csv ("root_metadata_L.csv",row.names = 1)
# write.csv (metadf, "root_metadata_F.csv") ## Clean metadata


## confirm identity
identical(colnames(t66),rownames(metadf)) # TRUE
identical(ncol(t66),nrow(metadf))  # TRUE

Taxonomy_1<-read.csv("Bacteria-allASV-DA.csv",row.names = 1,header = TRUE)
# grep("Unclassified",Taxonomy_1$Phylum) ## Should be 0
# ## grep("Chloroplast|Mitochondria",Taxonomy_1$Genus) ########## **********Remove from Genus level********* 
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
# write.csv (data1,"Rhizo_DA_Tret_6wk_all taxon.csv")
# write.csv (data1,"Root_DA_Tret_6wk_all taxon.csv")

# write.csv (data1,"RootF_all taxon.csv")  
# write.csv (data1,"Root_all taxon.csv")  
View (data1)
## Subset a dataframe to certain taxonomy level, and change taxon names

lol <- read.csv ("Rhizo_DA_Tret_6wk_Order.csv",row.names = 1)
lol <- read.csv ("Root_DA_Tret_6wk_Family.csv",row.names = 1)

View (lol)
log_r <- log (lol+1) 
View (log_r)

thaha <- merge(log_r, metadf,by='row.names', na.rm = TRUE)  ## Remove rows have NA
View (thaha)
tha0 <- data.frame(thaha[,-1], row.names=thaha[,1])
View (tha0)
colnames(tha0)

# tha <- tha0[, c(1:52,64:65)]  ## Rhizo
tha <- tha0[, c(1:40,50:51)]  ## Root

View (tha)
colnames(tha)
View (Taxonomy_1$Genus)

c_df <- Hmisc::rcorr(as.matrix(tha), type='spearman')  ## pearson
# write.csv (c_df$r,"CARP_Root16s_DA_Correlation_Coefficient.csv")
# write.csv (c_df$P,"CARP_Root16s_DA_Correlation_P-value.csv")

## Root
corrplot(corr=c_df$r[41:42,1:40], p.mat=c_df$P[41:42,1:40], sig.level=0.05,
         addCoef.col=1,  insig='blank', number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

corrplot(corr=c_df$r[1:40,1:40], p.mat=c_df$P[1:40,1:40], sig.level=0.05,
         addCoef.col=1,  insig='blank', number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

##### Rhizo
corrplot(corr=c_df$r[53:54,1:52], p.mat=c_df$P[53:54,1:52], sig.level=0.05,
         addCoef.col=1,  insig='blank', number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

corrplot(corr=c_df$r[53:54,1:52], p.mat=c_df$P[53:54,1:52], sig.level=0.01,
         addCoef.col=1,  insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

corrplot(corr=c_df$r[53:54,1:52], p.mat=c_df$P[53:54,1:52], sig.level=0.001,
         addCoef.col=1, insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))
################################################################################












#######################################################################################################################
########## Relative abundance between Treatments ###################################################################
##################################################################################################################
##### Use clr transformed data
# unique(phyloseq::tax_table(combined)[,"Phylum"])
View(combined) #10335 x 288
View (sub01)
## Subset root samples
sub_root <- subset_samples(sub01, SampleType == "Rhizo")
sub_root <- subset_samples(sub01, SampleType == "Root")
View (sub_root)
# unique(tax_table(sub_root)[,"Phylum"])
# # ## Subset samples for specific sampling year and growth stage
# sub_root01 <- subset_samples(sub_root, YS == "2019L")
# sub_root01 <- subset_samples(sub_root, YS == "2019F")
sub_root01 <- sub_root
View (sub_root01)
cyanos_glom <- tax_glom(sub_root01, taxrank = "Class")
Tax <-names(sort(rowSums(otu_table(cyanos_glom)),decreasing=T))
otu_table(cyanos_glom)<-otu_table(cyanos_glom)[intersect(Tax,rownames(otu_table(cyanos_glom))),]
tax_table(cyanos_glom)<-phyloseq::tax_table(cyanos_glom)[intersect(Tax,rownames(phyloseq::tax_table(cyanos_glom))),]
TT <- transform_sample_counts(cyanos_glom, function(x) log(x+1) ) ## relative abundance

phylum.sum = tapply(taxa_sums(TT), phyloseq::tax_table(TT)[, "Class"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:12]
GP1 = prune_taxa((phyloseq::tax_table(TT)[, "Class"] %in% top5phyla),TT) ## 67 x 96
View (GP1)
########
# # ## Family level (top12)
# cyanos_glom <- tax_glom(sub_root01, taxrank = "Family")
# Tax <-names(sort(rowSums(otu_table(cyanos_glom)),decreasing=T))
# otu_table(cyanos_glom)<-otu_table(cyanos_glom)[intersect(Tax,rownames(otu_table(cyanos_glom))),]
# tax_table(cyanos_glom)<-phyloseq::tax_table(cyanos_glom)[intersect(Tax,rownames(phyloseq::tax_table(cyanos_glom))),]
# TT <- transform_sample_counts(cyanos_glom, function(x) log(x+1) ) ## relative abundance
# 
# phylum.sum = tapply(taxa_sums(TT), phyloseq::tax_table(TT)[, "Family"], sum, na.rm=TRUE)
# top5phyla = names(sort(phylum.sum, TRUE))[1:12]
# GP1 = prune_taxa((phyloseq::tax_table(TT)[, "Family"] %in% top5phyla),TT) ## 67 x 96
# View (GP1)

######## 
cyanos_df <- psmelt(GP1)   ## melting phyloseq to dataframe
cyanos_df$Class <- as.factor(cyanos_df$Class)
# cyanos_df$Family <- as.factor(cyanos_df$Family)
View (cyanos_df)
# write.csv (cyanos_df, "Bulk_RA_Pairwise.csv")
levels (cyanos_df$Class)
levels (cyanos_df$Family)
# unique (cyanos_df$Phylum)

## Rhizo Class: "Acidimicrobiia","Acidobacteriae","Actinobacteria"     
# "Alphaproteobacteria", "Anaerolineae", "Blastocatellia"     
# "Chloroflexia","Clostridia","KD4-96"             
#  "Ktedonobacteria","Phycisphaerae","Saccharimonadia" 

## Root Class: "Bacteroidia" Gammaproteobacteria" 

######### Bacteria
speciesPalette <- c(Acidimicrobiia = "#8E2043", Acidobacteriae = "#59C7EB", Actinobacteria = "#FEA090",
                    Alphaproteobacteria ="#3E5496", Anaerolineae = "#9AA0A7",Blastocatellia = "#EFDC60",
                    Chloroflexia = "#E0607E", Clostridia = "orange", "KD4-96" = "black",
                    Ktedonobacteria = "gray", Phycisphaerae = "yellow", Saccharimonadia = "#0A9086")

##                    
## Planctomycetota,Patescibacteria,Gemmatimonadetes,Firmicutes,Verrucomicrobia,Patescibacteria
## Bacteroidota = "#EFDC60",Gammaproteobacteria= "",Patescibacteria="orange"
ggplot(data = cyanos_df, aes(x=factor(Tret, level=c('NP_LR', 'P_LR', 'NP_PR', 'P_PR', 'NP_WR', 'P_WR')), y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Bacterial Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+ 
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.minor = element_blank(),panel.background = element_blank()) ##+scale_colour_manual( values = speciesPalette)

# ggplot(data = cyanos_df, aes(x = Tret, y = Abundance)) +
#   geom_boxplot(outlier.shape  = NA) +
#   geom_jitter(aes(color = Class), height = 0, width = .2) +
#   labs(x = "", y = "Bacterial Log10 Abundance") +
#   facet_wrap(~ Family, scales = "free")+ theme_bw()+ 
#   theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         panel.grid.minor = element_blank(),panel.background = element_blank()) +scale_colour_manual( values = speciesPalette)


ggplot(data = cyanos_df, aes(x = RootType, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Bacterial Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())##+scale_colour_manual( values = speciesPalette)

## Compare P only if no interaction was detected    [Used for Bacteria]
ggplot(data = cyanos_df, aes(x = SoilP, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Phylum), height = 0, width = .2) +
  labs(x = "", y = "Bacterial Log10 Abundance") +
  facet_wrap(~ Class, scales = "free")+ theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())
# +scale_colour_manual( values = speciesPalette)

# View (cyanos_df)
# names (cyanos_df)
# levels (cyanos_df$Phylum)
levels (cyanos_df$Class)

## Rhizo/Root Class: Acidimicrobiia,Acidobacteriae, Actinobacteria, Alphaproteobacteria,
## Anaerolineae, Blastocatellia, Chloroflexia, Clostridia, KD4-96, Ktedonobacteria,
## Phycisphaerae, Saccharimonadia

View (cyanos_df)
ta1 <- subset (cyanos_df, Class =="Acidimicrobiia")  ## Acidobacteria
ta1 <- subset (cyanos_df, Class=="Acidobacteriae")    ## Actinobacteria
ta1 <- subset (cyanos_df, Class=="Actinobacteria")   ## Proteobacteria
ta1 <- subset (cyanos_df, Class=="Alphaproteobacteria")    ## Firmicutes
# ta1 <- subset (cyanos_df, Class=="Anaerolineae")  
ta1 <- subset (cyanos_df, Class=="Bacteroidia")  ## Bactero idetes
ta1 <- subset (cyanos_df, Class=="Blastocatellia")  ## Bactero idetes
ta1 <- subset (cyanos_df, Class=="Chloroflexia")  ## Firmicutes Chloroflexia
ta1 <- subset (cyanos_df, Class=="Clostridia")  ## Firmicutes Chloroflexia
ta1 <- subset (cyanos_df, Class=="Gammaproteobacteria")  ## Proteobacteria
ta1 <- subset (cyanos_df, Class=="KD4-96")  ## Gemmatimonadetes
ta1 <- subset (cyanos_df, Class=="Ktedonobacteria")  ## Chloroflexi
ta1 <- subset (cyanos_df, Class=="Phycisphaerae")  ## Verrucomicrobia
ta1 <- subset (cyanos_df, Class=="Saccharimonadia")  ## Patescibacteria
###### Verrucomicrobia
##### All RootTypes together ##### 
ta1
model<-lme(Abundance~SoilP*RootType,random=~1|Rep,data=ta1)
anova (model)

##### Separate RootTypes to check "soil P" effect ##### 
ta2 <- subset (ta1, RootType =="Lateral")
ta2 <- subset (ta1, RootType =="Primary")
ta2 <- subset (ta1, RootType =="Whole")
model<-lme(Abundance~SoilP,random=~1|Rep,data=ta2)
anova (model)
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidimicrobiia_RootClass_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidobacteriae_RootClass_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Saccharimonadia_RootClass_6wk.csv")

model<-lme(Abundance~RootType,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"RootType")
pair<-pairs(em,adjust = "tukey")
# write.csv (pair,"ANOVA_LogRA_Tret_Ktedonobacteria_RhizoClass_6wk_Pairwise.csv")
#####
model<-lme(Abundance~Tret,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"Tret")
pair<-pairs(em,adjust = "tukey")
write.csv (pair,"ANOVA_LogRA_Tret_Saccharimonadia_RootClass_6wk_Pairwise.csv")

model<-lme(Abundance~SoilP,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"SoilP")
pair<-pairs(em,adjust = "tukey")

######
###### Rhizo Family: 
## [1] "Blastocatellaceae","KD4-96","Microbacteriaceae","Micrococcaceae"      
# [5] "Nocardioidaceae"      "Propionibacteriaceae" "Rhizobiaceae"         "Sphingomonadaceae"   
# [9] "Streptomycetaceae"    "Unclassified"         "uncultured"           "Xanthobacteraceae" 
# ta1 <- subset (cyanos_df, Family == "Blastocatellaceae")  ## Acidobacteria
# ta1 <- subset (cyanos_df, Family=="KD4-96")    ## Actinobacteria
# ta1 <- subset (cyanos_df, Family=="Microbacteriaceae")   ## Proteobacteria
# ta1 <- subset (cyanos_df, Family=="Micrococcaceae")    ## Firmicutes
# ta1 <- subset (cyanos_df, Family=="Nocardioidaceae")  ## Bactero idetes
# ta1 <- subset (cyanos_df, Family=="Propionibacteriaceae")  ## Bactero idetes
# ta1 <- subset (cyanos_df, Family=="Rhizobiaceae")  ## Firmicutes Chloroflexia
# ta1 <- subset (cyanos_df, Family=="Sphingomonadaceae")  ## Proteobacteria
# ta1 <- subset (cyanos_df, Family=="Streptomycetaceae")  ## Gemmatimonadetes
# ta1 <- subset (cyanos_df, Family=="Unclassified")  ## Chloroflexi
# ta1 <- subset (cyanos_df, Family=="uncultured")  ## Verrucomicrobia
# ta1 <- subset (cyanos_df, Family=="Xanthobacteraceae")  ## Patescibacteria


ta1
model<-lme(Abundance~SoilP*RootType,random=~1|Rep,data=ta1)
anova (model)
qqnorm (model)
plot (model)
# pvalue.ttest <- t.test(Log10_Abundance ~ Group, data = df)$p.value  ## normal distribution
# pvalue.wilcoxon <- wilcox.test(Abundance ~ Tret, data = ta1)  ## un-normal distribution

# write.csv (anova (model),"ANOVA_LogRA_Tret_KD4-96_Rhizo_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Rhizobiaceae_Rhizo_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk.csv")
# write.csv (anova (model),"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk.csv")


model<-lme(Abundance~RootType,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"RootType")
pair<-pairs(em,adjust = "tukey")
# write.csv (pair,"ANOVA_LogRA_Tret_KD4-96_Rhizo_6wk_Pairwise.csv")
# write.csv (pair,"ANOVA_LogRA_Tret_Rhizobiaceae_Rhizo_6wk_Pairwise.csv")
write.csv (pair,"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk_Pairwise.csv")
write.csv (pair,"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk_Pairwise.csv")
write.csv (pair,"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk_Pairwise.csv")

write.csv (pair,"ANOVA_LogRA_Tret_Acidimicrobiia_Rhizo_6wk_Pairwise.csv")



model<-lme(Abundance~SoilP*RootType,random=~1|Rep,data=ta1)
anova (model)
em<-emmeans(model,"SoilP")
###################################################################################























############################################################################################
############### Differential Abundant ASVs _ Correlation ##############################
########################################################################################
View (combined)
sub_root01 <- subset_samples(combined, SampleType == "Root")
View (sub_root01)

# sub_root01 <- subset_samples(sub_root, Stage == "L")
# sub_root01 <- subset_samples(sub_root, Stage == "F")
# View (sub_root01)

## Change ASV column name to OTU
# DA_ASV <- read.csv("CARP_16S_DA_Root_Tret_L_TRUE_copy.csv",row.names = 1) ## mean >1
# DA_ASV <- read.csv("CARP_16S_DA_Root_Tret_F_TRUE_copy.csv",row.names = 1) ## mean >1
DA_ASV <- read.csv("CARP_16S_DA_Root_Tret_TRUE_copy.csv",row.names = 1) ## mean >1
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
View (df_root) 
# write.csv (root_otu, "Root_otu.csv")
df_root <- read.csv ("Root_otu.csv")


t11 <- merge(df_root, DA_ASV,by='OTU', na.rm = TRUE)  ## Remove rows have NA
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
# write.csv (metadf, "root_metadata_F.csv") ## Clean metadata
# metadf <- read.csv ("root_metadata_F.csv",row.names = 1)
# write.csv (metadf, "root_metadata.csv") ## Clean metadata
metadf <- read.csv ("root_metadata.csv",row.names = 1)

## confirm identity
identical(colnames(t66),rownames(metadf)) # TRUE
identical(ncol(t66),nrow(metadf))  # TRUE

Taxonomy_1<-read.csv("Bacteria-allASV-deleteD0.csv",row.names = 1,header = TRUE)
# grep("Unclassified",Taxonomy_1$Phylum) ## Should be 0
# ## grep("Chloroplast|Mitochondria",Taxonomy_1$Genus) ########## **********Remove from Genus level********* 
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

lol <- read.csv ("Root_taxon_Genus_copy.csv",row.names = 1)
View (lol)
log_r <- log (lol+1) 
View (log_r)

thaha <- merge(log_r, metadf,by='row.names', na.rm = TRUE)  ## Remove rows have NA
View (thaha)
tha0 <- data.frame(thaha[,-1], row.names=thaha[,1])
View (tha0)
colnames(tha0)

tha <- tha0[, c(1:47,66:78)]  ## All samples & Genus
View (tha)
colnames(tha)
View (Taxonomy_1$Genus)

c_df <- Hmisc::rcorr(as.matrix(tha), type='spearman')  ## pearson
# write.csv (c_df$r,"CARP_Root16s_DA_Correlation_Coefficient.csv")
# write.csv (c_df$P,"CARP_Root16s_DA_Correlation_P-value.csv")

corrplot(corr=c_df$r[48:60,1:47], p.mat=c_df$P[48:60,1:47], sig.level=0.05,
         addCoef.col=1,  insig='blank', number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

corrplot(corr=c_df$r[48:60,1:47], p.mat=c_df$P[48:60,1:47], sig.level=0.01,
         addCoef.col=1,  insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))

corrplot(corr=c_df$r[48:60,1:47], p.mat=c_df$P[48:60,1:47], sig.level=0.001,
         addCoef.col=1, insig='blank',number.digits = 1, tl.col="black",tl.cex = 0.8,
         number.cex=.6, tl.srt=45,col=brewer.pal(n=10, name="RdYlBu"))


