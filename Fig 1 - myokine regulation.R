#Load required packages
library(WGCNA)
library(ggplot2)
library(colormap)
library(dplyr)
library(ggrepel)
library(Rmisc)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
#with numbers
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

#DE of muscle genes in sex.  First we will create an index of sex variables then perform DE using DESeq2
sample_table = as.data.frame(row.names(working_dataset))
colnames(sample_table) = 'GTEx_ID'
row.names(sample_table) = sample_table$GTEx_ID
sample_table$Sex = sex_table$SEX[match(sample_table$GTEx_ID, sex_table$GTEx_ID)]

#Filter for myokines
myokine_table = working_dataset[,grepl('Muscle - Skeletal', colnames(working_dataset))]
colnames(myokine_table) = gsub('_Muscle - Skeletal', '', colnames(myokine_table))
mm1 = as.data.frame(t(myokine_table))
mm1 = mm1[,colSums(mm1 != 0,na.rm = TRUE) > 0]
mm1[is.na(mm1)] = 0

#Since DESeq2 requires integers, the class is applied to our dataframe
mm2 = as.data.frame(lapply(mm1, as.integer))
colnames(mm2) = colnames(mm1)
row.names(mm2) = row.names(mm1)

#Make a sample table used to contrast 
sample_table = sample_table[sample_table$GTEx_ID %in% colnames(mm2),]
sample_table$Sex = ifelse(sample_table$Sex==1, 'M', 'F')

#Filter for Myokines only
mm2 = mm2[row.names(mm2) %in% Secreted_proteins$`Gene names  (primary )`,]

#create a DESeq object
dds = DESeqDataSetFromMatrix(mm2, sample_table, design=~Sex)

#Filter for genes where >50 individauls have a count >10 in each myokine gene
keep = rowSums(counts(dds) >= 10) >= 50
dds <- dds[keep,]

#perform DE
dd1 = DESeq(dds)
res = results(dd1)
res1 = as.data.frame(res[order(res$padj, decreasing = F),])
res1$gene_symbol = row.names(res1)
res1 = na.omit(res1)

#The DE results are stored in res1.  The full table can be written:
write.csv(res1, file = 'Suppl Table - DE of all myokines.csv', row.names = F)

#The following parameters are used to create the volcano plot of myokines changing with sex
label_key = res1$gene_symbol[res1$padj<0.001]
res1$label2 = match(row.names(res1), label_key, nomatch = 0)
res1$label2 = ifelse(res1$label2>0, paste0(row.names(res1)), '')
res1$label_col1 = ifelse(res1$label2>0, 'blueviolet', '')
res1$label_col = ifelse(res1$log2FoldChange > 0, 'darkgoldenrod4', 'cyan3')
res1$label_col2 = ifelse(-log10(res1$pvalue) < 3, 'gray74', paste0(res1$label_col))
#Number of genes which will be labelled
table(res1$label_col2)
#cyan3 darkgoldenrod4         gray74 
#52             67           1672

#Volcano plot
ggplot(res1, aes(x=log2FoldChange, y=-log10(pvalue))) + theme_classic() +
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue)), color=res1$label_col2) +
  geom_label_repel(aes(x=log2FoldChange, y=-log10(pvalue), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = 400, alpha = .6, segment.color = 'grey50') +   ggtitle('Volcano plot sex-specific myokines status - female lower male higher')


#Subset genes based on DE category: male-specific, female-specific or non-sex
res1$cat1 = ifelse(res1$log2FoldChange>0, 'male-specific', 'female-specific')
res1$cat = ifelse(res1$pvalue<0.05, paste0(res1$cat1), 'non-sex regulated')
#plot the proportions of myokines in each DE category
ggplot(res1, aes(x = 1, fill = factor(cat))) +
  geom_bar(position="fill")


#Look at pathway enrichment corresponding to the muscle genes which show strongest 

#For each DE category, correlate all myokines with all skeletal muscle genes
#first in males
tissue1 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_datasetM[,grepl('Muscle - Skeletal', colnames(working_datasetM)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% Secreted_proteins$`Gene names  (primary )`]
tissue.tissue.all = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')

#melt the bicor coef and ad the pvalue
set1 = melt(tissue.tissue.all$bicor)
colnames(set1) = c('gene_symbol', 'myokine', 'male_bicor')
set2 = melt(tissue.tissue.all$p)
set1$male_pvalue = set2$value
set1$myo_gene = paste0(set1$myokine, '_', set1$gene_symbol)
sex_pathways = set1

#now repeat for females
tissue1 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_datasetF[,grepl('Muscle - Skeletal', colnames(working_datasetF)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% Secreted_proteins$`Gene names  (primary )`]

tissue.tissue.all = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')
set1 = melt(tissue.tissue.all$bicor)
colnames(set1) = c('gene_symbol', 'myokine', 'female_bicor')
set2 = melt(tissue.tissue.all$p)
set1$female_pvalue = set2$value
set1$myo_gene = paste0(set1$myokine, '_', set1$gene_symbol)

#add the female-specific cors to all pathways
sex_pathways$female_bicor = set1$female_bicor[match(sex_pathways$myo_gene, set1$myo_gene)]
sex_pathways$female_pvalue = set1$female_pvalue[match(sex_pathways$myo_gene, set1$myo_gene)]

#now repeat for both males and females
tissue1 <- working_dataset[,grepl('Muscle - Skeletal', colnames(working_dataset)),]
colnames(tissue1) = gsub("\\_.*","",colnames(tissue1))
tissue2 <- working_dataset[,grepl('Muscle - Skeletal', colnames(working_dataset)),]
colnames(tissue2) = gsub("\\_.*","",colnames(tissue2))
tissue2 = tissue2[,colnames(tissue2) %in% Secreted_proteins$`Gene names  (primary )`]
tissue.tissue.all = bicorAndPvalue(tissue1, tissue2, use='pairwise.complete.obs')

set1 = melt(tissue.tissue.all$bicor)
colnames(set1) = c('gene_symbol', 'myokine', 'bicor')
set2 = melt(tissue.tissue.all$p)
set1$pvalue = set2$value
set1$myo_gene = paste0(set1$myokine, set1$gene_symbol)

#Add the both sex data to table
sex_pathways$both_bicor = set1$bicor[match(sex_pathways$myo_gene, set1$myo_gene)]
sex_pathways$both_pvalue = set1$pvalue[match(sex_pathways$myo_gene, set1$myo_gene)]

#For each myokine add a category based on DE
sex_pathways$de_category = res1$cat[match(sex_pathways$myokine, res1$gene_symbol)]

#Now the dataframe sex_pathways contains all myokine ~ muscle gene correlations in each DE category
#subset by DE category then
myokine_reg_paths = sex_pathways[sex_pathways$de_category=='non-sex regulated',]

#group by skeletal muscle gene and average the -log10pvalue across all myokines to summarize strength of correlations
new_paths = myokine_reg_paths %>% dplyr::group_by(gene_symbol) %>% dplyr::summarise(avg_log = mean(-log10(both_pvalue)))
new_paths = new_paths[!is.infinite(new_paths$avg_log),]
new_paths = new_paths[!is.na(new_paths$avg_log),]
#order by significance across all myokines (-log10pvalue)
new_paths = new_paths[order(new_paths$avg_log, decreasing = T),]
#write file to be used for pathway enrichment
write.csv(new_paths, file = 'non-sex-specific muscle pathway enrichment input.csv', row.names = F)

#Repeat for other 2 DE categories
#female-specific DE:
myokine_reg_paths = sex_pathways[sex_pathways$de_category=='female-specific',]
new_paths = myokine_reg_paths %>% dplyr::group_by(gene_symbol) %>% dplyr::summarise(avg_log = mean(-log10(female_pvalue)))
new_paths = new_paths[!is.infinite(new_paths$avg_log),]
new_paths = new_paths[!is.na(new_paths$avg_log),]
new_paths = new_paths[order(new_paths$avg_log, decreasing = T),]
write.csv(new_paths, file = 'female-specific muscle pathway enrichment input.csv', row.names = F)

#male-specific DE
myokine_reg_paths = sex_pathways[sex_pathways$de_category=='male-specific',]
new_paths = myokine_reg_paths %>% dplyr::group_by(gene_symbol) %>% dplyr::summarise(avg_log = mean(-log10(male_pvalue)))
new_paths = new_paths[!is.infinite(new_paths$avg_log),]
new_paths = new_paths[!is.na(new_paths$avg_log),]
new_paths = new_paths[order(new_paths$avg_log, decreasing = T),]
write.csv(new_paths, file = 'male-specific muscle pathway enrichment input.csv', row.names = F)

#Next, for each DE category, look at enrichment of ESR1, AR, both or neither.  For this we use gene correlations across both sexes, but subset distributions of hormone receptor cors based on myokines in each DE category
#start by retrieving the coef and pvalue from dataframe above (both sexes)
melt_musc_cors = as.data.frame(melt(tissue.tissue.all$bicor))
melt_musc_corsp = as.data.frame(melt(tissue.tissue.all$p))
colnames(melt_musc_cors) = c('gene_1', 'myokine', 'bicor')
melt_musc_cors$pvalue = melt_musc_corsp$value

#exract all myokines ~ muscle ESR1
esr1_cors = melt_musc_cors[melt_musc_cors$gene_1=='ESR1', ]
esr1_cors = na.omit(esr1_cors)

#Summarize average statistics for each myokine
esr1_cors1 = esr1_cors %>% dplyr::group_by(myokine) %>% dplyr::summarise(pval = mean(pvalue), bicor = mean(bicor))

#Name categories for each hormone
esr1_cors1$myo_sig = ifelse(esr1_cors1$pval<0.05, 'ER', 'non-hormone')

#Add ESR1 cors to DE table from above
res1$esr1_cor = esr1_cors1$myo_sig[match(res1$gene_symbol, esr1_cors1$myokine)]

#Repeat for AR
esr1_cors = melt_musc_cors[melt_musc_cors$gene_1=='AR' , ]
esr1_cors = na.omit(esr1_cors)
esr1_cors1 = esr1_cors %>% dplyr::group_by(myokine) %>% dplyr::summarise(pval = mean(pvalue), bicor = mean(bicor))
esr1_cors1$myo_sig = ifelse(esr1_cors1$pval<0.05, 'AR', 'non-hormone')
res1$AR_cor = esr1_cors1$myo_sig[match(res1$gene_symbol, esr1_cors1$myokine)]

#Label categories for hormone receptor correlations
res1$hormone_cor = ifelse(res1$esr1_cor =='non-hormone' & res1$AR_cor == 'non-hormone', 'non-hormone', paste0(res1$esr1_cor))
res1$hormone_cor = ifelse(res1$esr1_cor =='ER' & res1$AR_cor == 'AR', 'Both ER and AR', paste0(res1$hormone_cor))
res1$hormone_cor = ifelse(res1$esr1_cor =='non-hormone' & res1$AR_cor == 'AR', 'AR', paste0(res1$hormone_cor))

#plot the proportions of homone receptor correlations within each category 
ggplot(res1, aes(x = factor(cat), fill = factor(hormone_cor))) +
  geom_bar(position="fill")

#Next, we will integrate ESR1KO DEGS
#Read in the DESeq2 results from ESR1KO vs WT in males or females
male_deg = read.csv('Muscle ERa KO Male DEG.csv')
female_deg = read.csv('Muscle ERa KO Female DEG.csv')

#Read in all human-mouse orthologs, downloaded from; http://www.informatics.jax.org/homology.shtml
orths = read.delim('Mouse Gene info with Human Orthologues.txt')
#Load mouse secreted protein list, also accesessed via uniprot
ms_secs = read.delim('secreted_proteins.txt')

#Filter DEGs for pvalue cutoffs
#male
deg_table = male_deg[male_deg$pvalue<0.05,]
deg_table$sex = paste0('M')
deg_table = deg_table[!(duplicated(deg_table$Gene.Name)),]
deg_table$X=NULL
#female
deg_table1 = female_deg[female_deg$pvalue<0.05,]
deg_table1$sex = paste0('F')
deg_table1 = deg_table1[!(duplicated(deg_table1$Gene.Name)),]

#Join lists of both sets of DEGs
deg_table  = as.data.frame(rbind(deg_table, deg_table1))

#look for significanlty DE in bothed sexes
deg_table$sex = ifelse(duplicated(deg_table$Gene.Name), 'Both', paste0(deg_table$sex))

#Annotated secreted proteins
deg_table$secreted = ifelse(match(deg_table$Gene.Name, ms_secs$Gene.names...primary.., nomatch = 0) >0 , 'Secreted', 'Non-secreted')
table(deg_table$secreted)
#Non-secreted     Secreted 
#1133          138 

#Add human orthologs to mosue gene symbols
deg_table$human_ort = orths$human_orth[match(deg_table$Gene.Name, orths$Symbol)]
deg_orths = deg_table[!is.na(deg_table$human_ort),]

#Join mouse-human orthologs with DE results from humans
er_cor = res1$gene_symbol[res1$esr1_cor=='ER']
new_deg = deg_table[deg_table$human_ort %in% er_cor,]
new_deg$gene_sex = paste0(new_deg$human_ort, '_', new_deg$sex)
table(res1$cat)

#Join by gene symbol - sex 
res1$new_sex = ifelse(res1$cat=='female-specific', 'F', 'Both')
res1$new_sex = ifelse(res1$cat=='male-specific', 'M', paste0(res1$new_sex))
res1$gene_sex = paste0(res1$gene_symbol, '_', res1$new_sex)

#Subset human myokines by ESR1-driven (gene correlates with muslce ESR1 in humans AND DE in WT vs KO Esr1 mice) or not
res1$ESR1_driven = ifelse(res1$gene_sex %in% new_deg$gene_sex, 'ESR1_driven', 'Non-ESR1')

#Pie charts for proportions of sex-specific DE in either ESR1 or non-ESR1-driven categories 
set1 = res1[res1$ESR1_driven=='ESR1_driven',]
colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('F', 'M', 'Both')
binned_sig_prots= set1 %>%
  dplyr::group_by(new_sex) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$new_sex, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$new_sex,col = binned_sig_prots$cols,  main = 'ESR1-driven myokines')

#distribution of myokines in non-ESR1-driven category
set1 = res1[!res1$ESR1_driven=='ESR1_driven',]
colors = c('orangered1', 'seagreen', 'mediumblue')
names(colors) = c('F', 'M', 'Both')
binned_sig_prots= set1 %>%
  dplyr::group_by(new_sex) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n))
binned_sig_prots$cols = colors[match(binned_sig_prots$new_sex, names(colors))]
pie(binned_sig_prots$freq, labels = binned_sig_prots$new_sex,col = binned_sig_prots$cols,  main = 'Non ESR1-driven myokines')

#Highlight MSTN1 as an example.  Add log2FC and hormone receptor correlations for MSTN
select_table = as.data.frame(res1$log2FoldChange[row.names(res1)=='MSTN'])
colnames(select_table) = 'human_log2FC'

esr1_myocor = sex_pathways[grepl('ESR1', sex_pathways$gene_symbol),]
esr1_myocor = esr1_myocor[grepl('MSTN', esr1_myocor$myokine),]
select_table$Female_ESR1cor = esr1_myocor$female_bicor
select_table$Male_ESR1cor = esr1_myocor$male_bicor

ar_myocor = sex_pathways[grepl('AR', sex_pathways$gene_symbol),]
ar_myocor = ar_myocor[grepl('MSTN', ar_myocor$myokine),]
ar_myocor = ar_myocor[ar_myocor$gene_symbol=='AR',]
select_table$Female_ARcor = ar_myocor$female_bicor
select_table$Male_ARcor = ar_myocor$male_bicor

#View MSTN ~ ESR1 or MSTN ~ AR correlations in muscle
breaksList = seq(-.7, .7, by = .1)
pheatmap(select_table[,2:5], fontsize_number = 30,  number_color = "black", color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), breaks = breaksList, main='Hormone receptor core - MSTN', fontsize_row = 5, fontsize_col = 5, cluster_rows = F, cluster_cols = F)

#subset log2FC of MSTN in both humans (sex) and mice (Esr1 KO vs WT)
mm1 = deg_orths[deg_orths$sex=='M',]
select_table$Male_MERK0log2FC = 2^(mm1$log2FoldChange[mm1$human_ort=='MSTN'])
mm1 = deg_orths[deg_orths$sex=='F',]
select_table$Female_MERK0log2FC = 1.02
new_table = cbind(select_table$human_log2FC, select_table$Male_MERK0log2FC, select_table$Female_MERK0log2FC)

#plot log2FC
breaksList = seq(0.1, 1.6, by = .01)
pheatmap(new_table, number_color = "black", color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(length(breaksList)), breaks = breaksList, main='log2FC human merko', fontsize_row = 5, fontsize_col = 5, cluster_rows = F, cluster_cols = F)

#write a csv file corresponding to the MSTN correlations in either males or females
targeted_pathways = sex_pathways[sex_pathways$myokine=='MSTN',]
write.csv(targeted_pathways, file = 'MSTN pathways.csv', row.names = F)
