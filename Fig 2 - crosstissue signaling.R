#Load required packages
library(WGCNA)
library(ggplot2)
library(colormaps)
library(dplyr)
library(ggrepel)
library(Rmisc)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(ggrepel)
allowWGCNAThreads()
#Load in the environment for full and filtered GTEx data
load('GTEx NA included env.RData')

#Here we will use the filtered dataset, where indidvudals must share > 1.2e6 gene-tissue combinations - listed as working dataset 
working_dataset=GTEx_subfiltered

#set row.names and transpose for downstream analyses
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))

#read in Secreted Proteins list
Secreted_proteins <- read.delim("uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab", header = T, check.names = F)

#Read in the annotation table conaining sex information
sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)
new_trts = sex_table[sex_table$GTEx_ID %in% row.names(working_dataset),]
table(new_trts$sexMF)
#  F   M 
#327 653

males = new_trts[new_trts$sexMF=='M',]
females = new_trts[!new_trts$sexMF=='M',]

#Subset two datasets for expression based on sex
working_datasetM = working_dataset[row.names(working_dataset) %in% males$GTEx_ID,]
working_datasetF = working_dataset[row.names(working_dataset) %in% females$GTEx_ID,]

#From these data, we will construct correlations between all myokines and target tissue gene expression, then bind results together.  A more efficent approach would be to directly correlated all myokines and a matrix of all tissue_gene combinations, then subset; however, given that the local memory required to create such a matrix in one-step exceed most individual workstations/servers, we performed individually
#create a matrix of myokines to correlate with all genes in subcutaneous adipose tissue
#First in males
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Adipose - Subcutaneous', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))

# Calculate cross-tissue correlation coefficents (midweight bicorrelation)
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'male_bicor_musc_adipSQ')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

#We create a matrix with which to bind following results
Ssec_bicor_matrix = ttp

#Repeat the process in females, using the working_datasetF
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Adipose - Subcutaneous', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))

tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

#Add the female coefficents to previous dataframe
Ssec_bicor_matrix$female_bicor_musc_adipSQ = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#These scripts are repeated for individual target tissues until line 252.  This can be tailored to any tissue of interest, where many more are provided in working_dataset
#Male Liver
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Liver', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_liver = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Female Liver
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Liver', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_liver = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Male Pancreas
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Pancreas', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_pancreas = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Feale Pancreas
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Pancreas', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_pancreas = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Male Visceral Adipose
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Adipose - Visceral', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_viscAdip = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Female Visceral Adipose
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Adipose - Visceral', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_viscAdip = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Male Hypothalamus
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Brain - Hypothalamus', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_hypothalamus = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Female Hypothalamus
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Brain - Hypothalamus', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_hypothalamus = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Male Heart
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Heart - Left Ventricle', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_heart = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Female Heart
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Heart - Left Ventricle', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_heart = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Male Small Intestine
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetM[,grepl('Small Intestine - Terminal Ileum', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$male_bicor_musc_intestine = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#Female Small Intestine
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue1 = tissue1[,colnames(tissue1) %in% Secreted_proteins$`Gene names  (primary )`]
tissue2 <- working_datasetF[,grepl('Small Intestine - Terminal Ileum', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue.tissue.p = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$bicor
ttp1 = tissue.tissue.p
ttp = melt(ttp1)
colnames(ttp) = c('gene1', 'gene2', 'bicor')
ttp$gene_gene = paste0(ttp$gene1, '_', ttp$gene2)

Ssec_bicor_matrix$female_bicor_musc_intestine = ttp$bicor[match(Ssec_bicor_matrix$gene_gene, ttp$gene_gene)]

#All of the coefficient results are stored in Ssec_bicor_matrix, where the columns are sex_tissue and rows are gene symbols
#The male vs female coefficent results are plotted for each tissue here:
ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_adipSQ, y=female_bicor_musc_adipSQ)) + geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> adipSQ Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_liver, y=female_bicor_musc_liver)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> Liver Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_pancreas, y=female_bicor_musc_pancreas)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> Pancreas Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_viscAdip, y=female_bicor_musc_viscAdip)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> viscAdip Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_hypothalamus, y=female_bicor_musc_hypothalamus)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> hypothalamus Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_heart, y=female_bicor_musc_heart)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> heart Correlations (by Sex)') + theme_minimal()

ggplot(Ssec_bicor_matrix, aes(x=male_bicor_musc_intestine, y=female_bicor_musc_intestine)) + 
  geom_hex(bins = 80) +
  scale_fill_continuous(type = "viridis")+ xlab('male bicor') + ylab('female bicor') + ggtitle('Muscle -> intestine Correlations (by Sex)') + theme_minimal()

#Next the dataframe will be melted so that specific combinations can be extracted directly.  The new dataframe (nn_melt1) will be used extensively in the remainder of analyses.      

new_nor = Ssec_bicor_matrix
row.names(new_nor) = new_nor$gene_gene
nn_melt = melt(as.matrix(new_nor))

#add column names for ID
colnames(nn_melt) = c('gene_gene', 'ID', 'bicor')
nn_melt1$bicor = as.numeric(nn_melt1$bicor)
nn_melt1 =nn_melt[!is.na(nn_melt$bicor),]

#This object is quite large since it is a melted table.  To limit space we will drop other large objects no longer being used:
nn_melt=NULL
GTEx_full=NULL
GTEx_subfiltered=NULL
tissue.tissue.p=NULL
new_nor = NULL
Ssec_bicor_matrix=NULL
ttp=NULL
ttp1=NULL

nn_melt1$sex = ifelse(grepl('female', nn_melt1$ID), 'F', 'M')
nn_melt1$target_tissue = gsub(".*musc_","", nn_melt1$ID)
nn_melt1$target_gene = gsub(".*_","", nn_melt1$gene_gene)
nn_melt1$musc_gene = gsub("\\_.*","",nn_melt1$gene_gene)
nn_melt1$bicor = as.numeric(nn_melt1$bicor)
#Next, sex-specificity, horomone receptor correlations are assessed for significant individual myokine ~ target tissue gene correlations.  The level of significance is calculated based on regression coefficent beyond 2 stadard deviations of the mean for all myokine ~ target tissue genes.  First we calculate the mean and SD for each target tissue here:  
summary_table = nn_melt1 %>% dplyr::select(gene_gene, sex, target_tissue, bicor) %>% dplyr::group_by(sex, target_tissue) %>% dplyr::summarise(mean=mean(bicor), sd=sd(bicor))

#Initially, we will calulate the pvalue (students pvalue from midweight bicorrelation) between all myokines and ESR1 or AR in skeletal muscle across both sexes or within each sex
#First calulate for both sexes
tissue1 <- working_dataset[,colnames(working_dataset)=='ESR1_Muscle - Skeletal' | colnames(working_dataset)=='AR_Muscle - Skeletal']
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_dataset[,grepl('Muscle - Skeletal', colnames(working_dataset)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
myobinning = melt(bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$p)
myobinning = dcast(myobinning, Var2 ~ Var1, value.var = 'value', fun.aggregate = mean)
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1<0.05, 'Both', 'non-hormone')
myobinning$sig = ifelse(myobinning$AR>0.05 & myobinning$ESR1<0.05, 'ESR1', paste0(myobinning$sig))
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1>0.05, 'AR', paste0(myobinning$sig))

#Make a table to merge all results
myobins_table = myobinning %>% mutate(sex =paste0('Both'))

#Correlate all myokines with hormone receptors in males
tissue1 <- working_datasetM[,colnames(working_datasetM)=='ESR1_Muscle - Skeletal' | colnames(working_datasetM)=='AR_Muscle - Skeletal']
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
myobinning = melt(bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$p)
myobinning = dcast(myobinning, Var2 ~ Var1, value.var = 'value', fun.aggregate = mean)
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1<0.05, 'Both', 'non-hormone')
myobinning$sig = ifelse(myobinning$AR>0.05 & myobinning$ESR1<0.05, 'ESR1', paste0(myobinning$sig))
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1>0.05, 'AR', paste0(myobinning$sig))

#Bind the male results to previous table
myobinning$sex = paste0('Male')
myobins_table = as.data.frame(rbind(myobins_table, myobinning))

#Correlate all myokines with hormone receptors in females
tissue1 <- working_datasetF[,colnames(working_datasetF)=='ESR1_Muscle - Skeletal' | colnames(working_datasetF)=='AR_Muscle - Skeletal']
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
myobinning = melt(bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')$p)
myobinning = dcast(myobinning, Var2 ~ Var1, value.var = 'value', fun.aggregate = mean)
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1<0.05, 'Both', 'non-hormone')
myobinning$sig = ifelse(myobinning$AR>0.05 & myobinning$ESR1<0.05, 'ESR1', paste0(myobinning$sig))
myobinning$sig = ifelse(myobinning$AR<0.05 & myobinning$ESR1>0.05, 'AR', paste0(myobinning$sig))

#Bind all of the results together for myokines ~ AR or ESR1 in females, males or both sexes
myobinning$sex = paste0('Female')
myobins_table = as.data.frame(rbind(myobins_table, myobinning))
myobins_table$gene_sex = paste0(myobins_table$Var2, myobins_table$sex)

#Next the top myokine ~ target tissue gene correlations will be extracted individually.  Note: the mean and sd bicor coefficents are entered from the table produced above (summary_table) containing these values for all tissues
#First Subcutaneous adipose tissue
one_tissue = nn_melt1[nn_melt1$target_tissue=='adipSQ',]
#set cutoffs for 2SD
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.019327166)+(2*0.1334038), abs(one_tissue$bicor) > (0.018629905)+(2*0.1013276))
table(one_tissue$sd_CO)

#subset significant correlations
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
#create a new table indicating significance per sex
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
head(new_table)
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))

#Now for the significant myokine ~ target tissue genes, bin by sex-specific or shared (both)
colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of SubQ Adipose')

#For the significant correlations, bin myokines by correlation with AR, ESR1, both or neither
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
#Plot
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Adipose Subq')



#Repeat for other tissues (heart)
one_tissue = nn_melt1[nn_melt1$target_tissue=='heart',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.025865997)+(2*0.1592068), abs(one_tissue$bicor) > (0.02017142)+(2*0.1133469))
table(one_tissue$sd_CO)
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
head(new_table)
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Heart')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr:: summarise(n = n()) %>%
  dplyr:: mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of Heart')


#Hypothalamus
one_tissue = nn_melt1[nn_melt1$target_tissue=='hypothalamus',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.029130417)+(2*0.2019200), abs(one_tissue$bicor) > (0.037425897)+(2*0.1119655))
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Hypothalamus')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of Hypothalamus')


#Intestine
one_tissue = nn_melt1[nn_melt1$target_tissue=='intestine',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.013006558)+(2*0.1982339), abs(one_tissue$bicor) > (0.011318279)+(2*0.1527021))
table(one_tissue$sd_CO)
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Intestine')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of Intestine')

#Liver
one_tissue = nn_melt1[nn_melt1$target_tissue=='liver',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.007574209)+(2*0.1975969), abs(one_tissue$bicor) > (0.021257840)+(2*0.1342818))
table(one_tissue$sd_CO)
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Liver')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of Liver')

#Pancreas
one_tissue = nn_melt1[nn_melt1$target_tissue=='pancreas',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.003544995)+(2*0.1624668), abs(one_tissue$bicor) > (0.008705634)+(2*0.1394398))
table(one_tissue$sd_CO)
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Pancreas')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of pancreas')


#viscAdip
one_tissue = nn_melt1[nn_melt1$target_tissue=='viscAdip',]
one_tissue$sd_CO = ifelse(one_tissue$sex=='F', abs(one_tissue$bicor) > (0.020451051)+(2*0.1437336), abs(one_tissue$bicor) > (0.026364824)+(2*0.1103125))
table(one_tissue$sd_CO)
one_tissue1 = one_tissue[!one_tissue$sd_CO=="FALSE",]
new_table = dcast(one_tissue1, gene_gene ~ sex, value.var = 'sd_CO')
new_table$sig_count = rowSums(is.na(new_table[,c(-1)]))
new_table$sig = ifelse(new_table$sig_count==0, 'Both', paste0(new_table$M))
new_table$sig = ifelse(new_table$sig=='TRUE', 'Male', paste0(new_table$sig))
new_table$sig = ifelse(new_table$sig=='NA', 'Female', paste0(new_table$sig))
new_table$myokine = sub("_.*", "", new_table$gene_gene)
new_table$gene_sex = paste0(new_table$myokine, new_table$sig)
myobinning1 = myobins_table[myobins_table$gene_sex %in% new_table$gene_sex,]
colors = c('palevioletred1', 'darkturquoise', 'forestgreen', 'purple1')
names(colors) = c('AR', 'ESR1', 'Both', "non-hormone")
binned_sig_prots= myobinning1 %>%
  dplyr::group_by(sex, sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
binned_sig_prots = na.omit(binned_sig_prots)
ggplot(data=binned_sig_prots, aes(x=" ", y=freq, group=sig, colour=sig, fill=sig)) +  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_grid(.~ sex) +theme_void() + ggtitle('hormone receptor cors Visceral adipose')

colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('Female', 'Male', 'Both')
binned_sig_prots= new_table %>%
  dplyr::group_by(sig) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$sig, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$sig,col = binned_sig_prots$cols,  main = '2SD cut-off sex enirchment of viscAdip')


#Now the class of myokine or crosstissue correlations will be subsetted to infer relative quantities of myokine regulation (ESR1) or crosstissue regression (sex-specific or not).  We use a bicor cutoff of 0.3, as that typically corresponds to P<1e-4
nn_new = nn_melt1[abs(nn_melt1$bicor) >0.3,]
nn_melt1 = NULL
nn_new = nn_new[!is.na(nn_new$bicor),]
class(nn_new$bicor)
nn_new$new_sex = ifelse(duplicated(nn_new$gene_gene), 'Both', paste0(nn_new$sex))
nn_new$gene_sex = paste0(nn_new$musc_gene, '_', nn_new$new_sex)
head(nn_new)

#Read in the ESR1 DEGs, filter for significant DE intersect with list of mouse-human orthologs
degs_female = read.csv('Muscle ERa KO Female DEG.csv')
degs_female$gene_sex = paste0(degs_female$Gene.Name, '_', 'F')
degs_male = read.csv('Muscle ERa KO Male DEG.csv')
degs_male$gene_sex = paste0(degs_male$Gene.Name, '_', 'M')
degs_male$X=NULL

full_degs = as.data.frame(rbind(degs_female, degs_male))
full_degs$gene_sex = ifelse(duplicated(full_degs$Gene.Name), paste0(full_degs$Gene.Name, '_', 'Both'), paste0(full_degs$gene_sex))
table(full_degs$gene_sex)
full_degs = full_degs[full_degs$pvalue<0.05,]
deg_orths = read.delim('Mouse Gene info with Human Orthologues.txt')

#filter orthologs list for significant DEGS
full_degs$human_orth_sex = deg_orths$human_orth[match(full_degs$Gene.Name, deg_orths$Symbol)]
full_degs$human_orth_sex = paste0(full_degs$human_orth_sex, '_', gsub(".*_","", full_degs$gene_sex))

#Now that nn_new has significant myokine-target gene combinations, categories of regulation (ESR1-drived) or if crosstissue signaling is sex-specific is added
nn_new$ESR1_driven = ifelse(nn_new$gene_sex %in% full_degs$human_orth_sex, 'ESR1_driven', 'Non-ESR1')
nn_new$source_count = ifelse(nn_new$ESR1_driven=='ESR1_driven' & nn_new$new_sex=='Both', 'ESR1_driven', 'none')
nn_new$source_count = ifelse(nn_new$ESR1_driven=='ESR1_driven' & nn_new$new_sex=='F', 'ESR1_driven_and_SS', paste0(nn_new$source_count))
nn_new$source_count = ifelse(nn_new$ESR1_driven=='ESR1_driven' & nn_new$new_sex=='M', 'ESR1_driven_and_SS', paste0(nn_new$source_count))
nn_new$source_count = ifelse(nn_new$ESR1_driven=='Non-ESR1' & nn_new$new_sex=='M', 'SS', paste0(nn_new$source_count))
nn_new$source_count = ifelse(nn_new$ESR1_driven=='Non-ESR1' & nn_new$new_sex=='F', 'SS', paste0(nn_new$source_count))

#see the total numbers of significant corelations
table(nn_new$source_count)
#ESR1_driven ESR1_driven_and_SS               none                 SS 
#1471860             497063           19837067           43961752

#summarzize these as number of total significant crossitssue correlations per myokine, per condition 
quant_tables = nn_new %>% dplyr::group_by(target_tissue, source_count) %>% dplyr::summarise(n=n(), n_myo=length(unique((musc_gene))))
quant_tables$normalized_signal = quant_tables$n/quant_tables$n_myo

#plot the result
plot_data = quant_tables
plot_data$source_count = factor(plot_data$source_count, levels = c('none', 'ESR1_driven', 'SS', 'ESR1_driven_and_SS'))
ggplot(plot_data, aes(x=target_tissue, y=normalized_signal, fill=source_count)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c( 'burlywood3', 'darkgoldenrod1', 'brown2',  'darkslategray')) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('regulation per myokine') + xlab('')






