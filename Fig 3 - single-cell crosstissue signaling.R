#Load required packages
library(reticulate)
library(Rcpp)
library(Seurat)
library(reshape2)
library(colormap)
library(patchwork)
library(dplyr)
library(limma)
library(devtools)
library(ADAPTS)
library(preprocessCore)
library(pheatmap)
library(ggplot2)
library(WGCNA)
library(mclust)
library(pheatmap)

doParallel::registerDoParallel(cores = parallel::detectCores())


#Read in the single-cell RNA_seq datasets.  Note, these are 4 individuals which are read in one-by-one then combined in Seurat using integration (below)
#Read in first counts matrix
dataset1 = read.csv('GSM3746212_Muscle_1_Counts.csv')
row.names(dataset1) = dataset1$X
dataset1$X=NULL
colnames(dataset1) = gsub('_', '-', colnames(dataset1), fixed = T)
row.names(dataset1) = gsub('_', '-', row.names(dataset1), fixed = T)

#create table to merge results
musc1 = CreateSeuratObject(counts = dataset1, project = "skl_musc_decon", min.cells = 3, min.features = 200)

#read in second individual
dataset1 = read.csv('GSM3746213_Muscle_2_Counts.csv')
row.names(dataset1) = dataset1$X
dataset1$X=NULL
colnames(dataset1) = gsub('_', '-', colnames(dataset1), fixed = T)
row.names(dataset1) = gsub('_', '-', row.names(dataset1), fixed = T)

musc2 = CreateSeuratObject(counts = dataset1, project = "skl_musc_decon", min.cells = 3, min.features = 200)

#Third individual
dataset1 = read.csv('GSM3746214_Muscle_3_Counts.csv')
row.names(dataset1) = dataset1$X
dataset1$X=NULL
colnames(dataset1) = gsub('_', '-', colnames(dataset1), fixed = T)
row.names(dataset1) = gsub('_', '-', row.names(dataset1), fixed = T)

musc3 = CreateSeuratObject(counts = dataset1, project = "skl_musc_decon", min.cells = 3, min.features = 200)

#Fourth individual
dataset1 = read.csv('GSM3746215_Muscle_4_Counts.csv')
row.names(dataset1) = dataset1$X
dataset1$X=NULL
colnames(dataset1) = gsub('_', '-', colnames(dataset1), fixed = T)
row.names(dataset1) = gsub('_', '-', row.names(dataset1), fixed = T)

musc4 = CreateSeuratObject(counts = dataset1, project = "skl_musc_decon", min.cells = 3, min.features = 200)

#Next, we will use Seurat to merge datasets into seurat objects.  These were adopted from the guided tutorial available at: https://satijalab.org/seurat/archive/v3.1/merge_vignette.html

#list muscle sc-seq datasets
set_test = c(musc1, musc2, musc3, musc4)

ifnb.list <- lapply(X = set_test, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- JackStraw(immune.combined, num.replicate = 100)
immune.combined <- ScoreJackStraw(immune.combined, dims = 1:20)
JackStrawPlot(immune.combined, dims = 1:20)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:18)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:18)
#immune.combined <- FindClusters(immune.combined, resolution = 0.5) #12 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.3) #11 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.2) #11 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.1) #9 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.4) #12 clusters
immune.combined <- FindClusters(immune.combined, resolution = 0.6) #14 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.7) #11 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.8) #11 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 0.9) #11 clusters
#immune.combined <- FindClusters(immune.combined, resolution = 1.1) #11 clusters

DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)

#looks like 12 cell types
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindMarkers(immune.combined, ident.1 = 0, verbose = FALSE)
#nk.markers2 <- limma()
#looks like myotubes
nk.markers$cluster = paste0('0')
marker_set = nk.markers

#Next, individual gene enrichment per cluster are written to a merged table
#Cluster 0
nk.markers <- FindMarkers(immune.combined, ident.1 = 0, verbose = FALSE)
nk.markers$cluster = paste0('0')
marker_set = nk.markers
#Cluster 1
nk.markers <- FindMarkers(immune.combined, ident.1 = 1, verbose = FALSE)
nk.markers$cluster = paste0('1')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 2
nk.markers <- FindMarkers(immune.combined, ident.1 = 2, verbose = FALSE)
nk.markers$cluster = paste0('2')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 3
nk.markers <- FindMarkers(immune.combined, ident.1 = 3, verbose = FALSE)
nk.markers$cluster = paste0('3')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 4
nk.markers <- FindMarkers(immune.combined, ident.1 = 4, verbose = FALSE)
nk.markers$cluster = paste0('4')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 5
nk.markers <- FindMarkers(immune.combined, ident.1 = 5, verbose = FALSE)
nk.markers$cluster = paste0('5')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 6
nk.markers <- FindMarkers(immune.combined, ident.1 = 6, verbose = FALSE)
nk.markers$cluster = paste0('6')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 7
nk.markers <- FindMarkers(immune.combined, ident.1 = 7, verbose = FALSE)
nk.markers$cluster = paste0('7')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 8
nk.markers <- FindMarkers(immune.combined, ident.1 = 8, verbose = FALSE)
nk.markers$cluster = paste0('8')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 9
nk.markers <- FindMarkers(immune.combined, ident.1 = 9, verbose = FALSE)
nk.markers$cluster = paste0('9')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 10
nk.markers <- FindMarkers(immune.combined, ident.1 = 10, verbose = FALSE)
nk.markers$cluster = paste0('10')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 11
nk.markers <- FindMarkers(immune.combined, ident.1 = 11, verbose = FALSE)
nk.markers$cluster = paste0('11')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 12
nk.markers <- FindMarkers(immune.combined, ident.1 = 12, verbose = FALSE)
nk.markers$cluster = paste0('12')
marker_set = as.data.frame(rbind(marker_set, nk.markers))
#Cluster 13
nk.markers <- FindMarkers(immune.combined, ident.1 = 13, verbose = FALSE)
nk.markers$cluster = paste0('13')
marker_set = as.data.frame(rbind(marker_set, nk.markers))

#The combined table contains all marker genes and enrichments per cell cluster
head(marker_set)
new_seurat_object = immune.combined

#From the enrichment of top 30 genes (by pvalue of cluster enrichment) cell types annotations were assigned manually.  These are read in as a .csv file:
marker_genes = read.csv('cluster_annotations.csv')
new.cluster.ids <- marker_genes$Cell.type.Call
names(new.cluster.ids) <- as.numeric(as.character(gsub('cluster ', '', marker_genes$Cluster)))

#Now the cell clusters are renamed in the Seurat object
new_seurat_object <- RenameIdents(new_seurat_object, new.cluster.ids)

#Plot the renamed and revised UMAP
DimPlot(new_seurat_object, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

#Raw counts with renamed cell types is extracted to input into deconvolution 
raw_counts =as.data.frame(new_seurat_object@assays$RNA@counts)
annot_link = as.data.frame(marker_genes)
annot_link$cluster_number = as.numeric(as.character(gsub('cluster ', '', marker_genes$Cluster)))
cell_ID = new_seurat_object@meta.data
cell_ID$cell_call = annot_link$Cell.type.Call[match(cell_ID$seurat_clusters, annot_link$cluster_number)]
colnames(raw_counts) = cell_ID$cell_call[match(colnames(raw_counts), row.names(cell_ID))]

#The raw_counts dataframe now contains names matched to cell types
#Filtered bulk data is now read in to deconvolute.  Note: additional filtering was applied to skeletal muscle here as deconvolution tools can be subjected to influence of genes with high variance.  Here, muscle data was further filtered for genes present in >50% of individuals.  Further, individuals were removed with <50% genes present, resulting in 37,413 skeletal muscle genes in  298 individuals 
#read in bulk data
bulk_data = read.csv('filtered muscle bulk data for adapts.csv')
bulk_data = bulk_data[!is.na(bulk_data$Gene),]
row.names(bulk_data) = bulk_data$Gene
bulk_data$Gene=NULL

#named objects and feed into ADAPTS
musc_bulk = as.data.frame(bulk_data)
musc_sc = raw_counts
musc_sc = musc_sc[!duplicated(row.names(musc_sc)),]
colnames(musc_sc) = paste0('Celltype', '_', colnames(musc_sc))
musc_bulk1 = musc_bulk[, colSums(musc_bulk != 0) > 0]
new_sc = as.data.frame(sapply(unique(colnames(musc_sc)), function(x) rowMeans(musc_sc[,colnames(musc_sc) == x])))

#Calculate cell proportions using the 'proportionsInAdmixture' approch.  Other approaches (nnls and DCQ) failed to capture as many numbers of cell types
musc_proportions = as.data.frame(estCellPercent.proportionsInAdmixture(new_sc, musc_bulk1))
musc_proportions$cell_type_mean = rowMeans(musc_proportions)
musc_proportions$cell_type = row.names(musc_proportions)

#Proportions are estimated and compared between sexes
musc_proportions = musc_proportions[!row.names(musc_proportions)=='others',]
mm1 = as.data.frame(t(musc_proportions))
colnames(mm1) = gsub('Celltype_', '', colnames(mm1))
row.names(mm1) = gsub('Celltype_', '', row.names(mm1))

#Read in table containing biologic sex info
sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
sex_table$xID = gsub(".*-","",sex_table$SUBJID)

proportions_melt = reshape2::melt(as.matrix(mm1))
colnames(proportions_melt) =c('GTEx_ID', 'Musc_cell', 'proportion')
proportions_melt$GTEx_ID = substring(proportions_melt$GTEx_ID, 2)

proportions_melt$sex = sex_table$sexMF[match(proportions_melt$GTEx_ID, sex_table$xID)]
proportions_melt$proportion = as.numeric(proportions_melt$proportion)
proportions_melt = proportions_melt[!is.na(proportions_melt$proportion),]
proportions_melt = proportions_melt[!is.na(proportions_melt$sex),]

binned_sig_prots= proportions_melt %>%
  group_by(sex, Musc_cell) %>%
  summarise(mean = mean(proportion, na.rm=T)) 

#plot the cell compositions between sexes
ggplot(binned_sig_prots, aes(x=sex, y=mean, fill=Musc_cell)) + geom_bar(position="fill", stat="identity") + theme_classic() + ggtitle('estimated cell proportions')

#Next, cell compositions will be integrated across tissues, so the full GTEx dataset is read in
load('GTEx NA included env.RData')

#Extract the dataset which has been filtered
working_dataset=GTEx_subfiltered
row.names(working_dataset) = working_dataset$gene_tissue
working_dataset$gene_tissue=NULL
working_dataset = as.data.frame(t(working_dataset))

#read in sex annotations and filter cell proportions by sex
sex_table = read.delim('GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
sex_table$GTEx_ID = gsub('GTEX-', '', sex_table$SUBJID)
sex_table$sexMF = ifelse(sex_table$SEX==1, 'M', 'F')
table(sex_table$sexMF)
sex_table$xID = paste0('X', sex_table$GTEx_ID)
fems = sex_table[sex_table$sexMF=='F',]

#Filter cell proportions dataframe (mm1) for females only then assess correlation between cell compositions
m2 = mm1[row.names(mm1) %in% fems$xID,]
cc1 = bicorAndPvalue(m2, m2, use = 'p')
cc2 = cc1$bicor
cc2[is.na(cc2)] = 0
tt4 = cc1$p
tt4[is.na(tt4)] = 1
tt3 = ifelse(tt4 < 0.001,"*","")

#plot composition correlation structure
pheatmap(cc2, fontsize_number = 20, display_numbers = tt3, number_color = "black", color = colormap(colormap = colormaps$viridis, nshades = 50), main='Decon cell proportions Females', fontsize_row = 5, fontsize_col = 5)

#Filter for males only then do the same
m2 = mm1[!row.names(mm1) %in% fems$xID,]
cc1 = bicorAndPvalue(m2, m2, use = 'p')
cc2 = cc1$bicor
cc2[is.na(cc2)] = 0
tt4 = cc1$p
tt4[is.na(tt4)] = 1
tt3 = ifelse(tt4 < 0.001,"*","")

#plot heatmap
pheatmap(cc2, fontsize_number = 20, display_numbers = tt3, number_color = "black", color = colormap(colormap = colormaps$viridis, nshades = 50), main='Decon cell proportions Males', fontsize_row = 5, fontsize_col = 5)

#Next the cell proportions will be separated by sex, then correlated with the same sex across tissues
fem_decon = mm1[row.names(mm1) %in% fems$xID,]
fem_decon$GTEx_ID = substring(row.names(fem_decon), 2)
male_decon = mm1[!row.names(mm1) %in% fems$xID,]
male_decon$GTEx_ID = substring(row.names(male_decon), 2)

#A more efficient approach compared to Fig 2 (since there are few cell proportions) is to select all tissues first then correlate accross
males = sex_table[sex_table$sexMF=='M',]
females = sex_table[!sex_table$sexMF=='M',]
working_dataset1 = working_dataset[ ,grepl('Adipose - Subcutaneous', colnames(working_dataset)) | grepl('Liver', colnames(working_dataset)) | grepl('Pancreas', colnames(working_dataset)) | grepl('Adipose - Visceral', colnames(working_dataset)) | grepl('Brain - Hypothalamus', colnames(working_dataset)) | grepl('Adipose - Visceral', colnames(working_dataset)) | grepl('Heart - Left Ventricle', colnames(working_dataset)) | grepl('Small Intestine - Terminal Ileum', colnames(working_dataset)) | grepl('Muscle - Skeletal', colnames(working_dataset))]

#remove muscle
working_dataset2 = working_dataset1[,!grepl('Muscle - Skeletal', colnames(working_dataset1))]

#separate by sex
working_datasetM = working_dataset1[row.names(working_dataset1) %in% males$GTEx_ID,]
working_datasetF = working_dataset1[row.names(working_dataset1) %in% females$GTEx_ID,]

#Since the cell proportions contains less individuals, filter for those across tissues
working_datasetF = working_datasetF[row.names(working_datasetF) %in% fem_decon$GTEx_ID,]
fem_decon = fem_decon[fem_decon$GTEx_ID %in% row.names(working_datasetF),]
working_datasetF = working_datasetF[order(row.names(working_datasetF)), ]
fem_decon = fem_decon[order(fem_decon$GTEx_ID),]
row.names(fem_decon) = fem_decon$GTEx_ID
fem_decon$GTEx_ID=NULL

#Correlate female muscle proportions across tissues
fems_cors = bicorAndPvalue(fem_decon, working_datasetF, use = 'p')
cell_cors = melt(fems_cors$bicor)
colnames(cell_cors) = c('muscle_cell_type', 'gene_tissue', 'bicor')
ccc1 = melt(fems_cors$p)
cell_cors$pvalue = ccc1$value
cell_cors$sex = paste0('F')
full_cell_cors = cell_cors

#Filter the males similarly and correlate across tissues
working_datasetM = working_datasetM[row.names(working_datasetM) %in% male_decon$GTEx_ID,]
male_decon = male_decon[male_decon$GTEx_ID %in% row.names(working_datasetM),]
working_datasetM = working_datasetM[order(row.names(working_datasetM)), ]
male_decon = male_decon[order(male_decon$GTEx_ID),]
row.names(male_decon) = male_decon$GTEx_ID
male_decon$GTEx_ID=NULL
fems_cors = bicorAndPvalue(male_decon, working_datasetM, use = 'p')
cell_cors = melt(fems_cors$bicor)
colnames(cell_cors) = c('muscle_cell_type', 'gene_tissue', 'bicor')
ccc1 = melt(fems_cors$p)
cell_cors$pvalue = ccc1$value
cell_cors$sex = paste0('M')

#create a full table for crosstissue correlations of cell proportions in each sex
full_cell_cors = as.data.frame(rbind(full_cell_cors, cell_cors))

#Subset significant crosstissue correlations and plot frequency per target tissue
sig_cell_cors = full_cell_cors
sig_cell_cors$gene_symbol = gsub("\\_.*","",sig_cell_cors$gene_tissue)
sig_cell_cors$tissue = sub(".*?_", "", sig_cell_cors$gene_tissue)
sig_cell_cors=na.omit(sig_cell_cors)
sig_cell_cors$musccell_sex = paste0(sig_cell_cors$muscle_cell_type, '_', sig_cell_cors$sex)

#remove muscle since these are crosstissue and filter for pvalue<0.001
sig_cell_cors1 = sig_cell_cors[!sig_cell_cors$tissue=='Muscle - Skeletal',]
sig_cell_cors1 =sig_cell_cors1[sig_cell_cors1$pvalue<0.001,]

#plot
ggplot(sig_cell_cors1, aes(x=musccell_sex, fill = tissue)) +
  geom_bar(position = "identity") + theme_classic() + ggtitle('Cell composition x sig cors p less0.01') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Next look at skeletal muscle hormone receptor correlations with cell composition
binned_bicor_table = sig_cell_cors[sig_cell_cors$tissue=='Muscle - Skeletal',]
binned_bicor_table$cell_sex = paste0(binned_bicor_table$muscle_cell_type, '_', binned_bicor_table$sex)
binned_bicor_table$gene_sex = paste0(binned_bicor_table$gene_symbol, '_', binned_bicor_table$sex)
binned_bicor_table$gee_tissue_sex = paste0(binned_bicor_table$gene_sex, '_', binned_bicor_table$tissue)

#Next correlate skeletal muscle hormone receptor expression with composition metrics
#start with females
esr_cell_cors = binned_bicor_table[binned_bicor_table$gene_symbol=='ESR1' | binned_bicor_table$gene_symbol=='AR',]
esr_heat = esr_cell_cors[esr_cell_cors$sex=='F',]
esr_heat = dcast(esr_heat,  muscle_cell_type ~ gene_tissue, value.var = 'bicor',fun.aggregate = mean)
row.names(esr_heat) = esr_heat$muscle_cell_type
esr_heat$muscle_cell_type=NULL
esr_heat = as.matrix(esr_heat)
esr_heat1 = esr_cell_cors[esr_cell_cors$sex=='F',]
esr_heat1 = dcast(esr_heat1,  muscle_cell_type ~ gene_tissue, value.var = 'pvalue',fun.aggregate = mean)
row.names(esr_heat1) = esr_heat1$muscle_cell_type
esr_heat1$muscle_cell_type=NULL
esr_heat1 = as.matrix(esr_heat1)
labs_pval = ifelse(esr_heat1 < 0.05, '*', '')

#plot
pheatmap(esr_heat, fontsize_number = 20, display_numbers = labs_pval, number_color = "black", color = colormap(colormap = colormaps$plasma, nshades = 50), main='Decon cell proportions AR and ER Females', fontsize_row = 5, fontsize_col = 5)

#The same in males
esr_cell_cors = binned_bicor_table[binned_bicor_table$gene_symbol=='ESR1' | binned_bicor_table$gene_symbol=='AR',]
esr_heat = esr_cell_cors[esr_cell_cors$sex=='M',]
esr_heat = dcast(esr_heat,  muscle_cell_type ~ gene_tissue, value.var = 'bicor',fun.aggregate = mean)
row.names(esr_heat) = esr_heat$muscle_cell_type
esr_heat$muscle_cell_type=NULL
esr_heat = as.matrix(esr_heat)
esr_heat1 = esr_cell_cors[esr_cell_cors$sex=='M',]
esr_heat1 = dcast(esr_heat1,  muscle_cell_type ~ gene_tissue, value.var = 'pvalue',fun.aggregate = mean)
row.names(esr_heat1) = esr_heat1$muscle_cell_type
esr_heat1$muscle_cell_type=NULL
esr_heat1 = as.matrix(esr_heat1)

labs_pval = as.matrix(ifelse(esr_heat1 < 0.05, '*', ''))


pheatmap(esr_heat, fontsize_number = 20, display_numbers = labs_pval, number_color = "black", color = colormap(colormap = colormaps$plasma, nshades = 50),  main='Decon cell proportions AR and ER Males', fontsize_row = 5, fontsize_col = 5)

#Next find the strongest-correlated myokines with each cell proporition in each sex to plot
musc_genes = binned_bicor_table[binned_bicor_table$tissue=='Muscle - Skeletal',]
musc_genes = musc_genes[order(musc_genes$pvalue, decreasing = F),]
sec_prots = read.delim('uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab')
musc_genes = musc_genes[musc_genes$gene_symbol %in% sec_prots$Gene.names...primary..,]
musc_genes$cell_sex = paste0(musc_genes$muscle_cell_type, '_', musc_genes$sex)
musc_genes1 = musc_genes[!duplicated(musc_genes$cell_sex),]
musc_genes1$x_order = factor(musc_genes1$cell_sex[order(musc_genes1$cell_sex)])
musc_genes1$x_order
musc_genes1$logp11 = -log10(musc_genes1$pvalue[order(musc_genes1$cell_sex)])
musc_genes1$dir_col = ifelse(musc_genes1$bicor>0, 'firebrick3', 'dodgerblue3')
musc_genes1 = musc_genes1[order(musc_genes1$cell_sex),]

#plot the top-ranked myokines for each cell proportion and sex
ggplot(musc_genes1, aes(x=cell_sex, y=logp11)) +
  geom_point(size = 3, colour = musc_genes1$dir_col, aes(color=))   +geom_segment(aes(xend = cell_sex, yend = 0), size = 1.2, colour = musc_genes1$dir_col)+
  geom_label(aes(cell_sex, logp11+1.5, label = gene_symbol, colour = dir_col,  size = 0.02))+
  labs(y= "-log10(pvalue)", x="cell_sex") + theme_classic()+  theme(axis.text.x = element_text(angle =45, vjust = 0, hjust=0.1))


#Finally, find example correlations between cell comp and tissue pathways adjusted for myokine.  Here we will focus on glycolytic fast-twitch fiber types
#subset crosstissue glyolytic fiber correlations
glyco_pathways = full_cell_cors[full_cell_cors$muscle_cell_type=='fast_twitch_glycolytic_fiber',]

#order by pvalue, which are used for pathway enrichemnts
glyco_pathways = glyco_pathways[order(glyco_pathways$pvalue, decreasing = T),]
write.csv(glyco_pathways, file = 'fast-twitch cell comp crosstissue no pval filter.csv', row.names = F)

#The enrichments suggested female-specific enrichments with pancreatic protein synthesis and male-specific enrichments with liver immune response and exocytosis.  We took the top 5 genes for each pathways for this analysis
#Select specific genes corresponding to pathways in pancreas and liver
panc_genes = c('RPL10', 'RPL10A', 'RPL11', 'RPL13A', 'RPL14')
panc_genes = paste0( panc_genes, '_Pancreas')
liv_immune_genes = c('CD68', 'ALOX5', 'GRN', 'NFKB1', 'TNFRSF1B')
liv_immune_genes = paste0(liv_immune_genes, '_Liver')
liv_exo_genes = c('STX11', 'CRP', 'EXOC3L1', 'VPS26A', 'SAA1')
liv_exo_genes = paste0(liv_exo_genes, '_Liver')

#make a key to isolate gene specifically
gene_key = c(panc_genes, liv_immune_genes, liv_exo_genes)

#focus correlations on these genes only
male_pathway_genes = working_datasetM[,colnames(working_datasetM) %in% gene_key]
female_pathay_genes = working_datasetF[,colnames(working_datasetF) %in% gene_key]
mm1 = male_decon[row.names(male_decon) %in% row.names(male_pathway_genes),]
mm2 = male_pathway_genes[row.names(male_pathway_genes) %in% row.names(mm1),]

#Now correlate cell compositions in each sex with the top genes
raw_cors = bicorAndPvalue(mm2, mm1$fast_twitch_glycolytic_fiber, use = 'p')

#Make an new table to bind results from adjusted regressions
adj_table = raw_cors$p
ad_table = as.data.frame(adj_table[order(ordered(row.names(adj_table), levels = gene_key)),])
colnames(ad_table) = 'male_fasttwitch'
heat_table = ad_table

#focus on myokines now - gpx3
gpx_male = as.data.frame(working_dataset$`GPX3_Muscle - Skeletal`)
row.names(gpx_male) = row.names(working_dataset)
colnames(gpx_male) = 'GPX3'
gp3 = gpx_male[row.names(gpx_male) %in% row.names(mm1),]

#calculate residuals from glycolytic fiber composition ~ skeletal muscle GPX3
resids = summary(lm(mm1$fast_twitch_glycolytic_fiber ~ gp3))$residuals
raw_cors = bicorAndPvalue(mm2, resids, use = 'p')
adj_table = raw_cors$p
ad_table = as.data.frame(adj_table[order(ordered(row.names(adj_table), levels = gene_key)),])
colnames(ad_table) = 'male_fasttwitch_gpx3_adj'

#Bind these results together
heat_table$male_fasttwitch_gpxadj = ad_table$male_fasttwitch_gpx3_adj[match(row.names(heat_table), row.names(ad_table))]

#Repeat the analysis for females
mm1 = fem_decon[row.names(fem_decon) %in% row.names(female_pathay_genes),]
mm2 = female_pathay_genes[row.names(female_pathay_genes) %in% row.names(mm1),]
raw_cors = bicorAndPvalue(mm2, mm1$fast_twitch_glycolytic_fiber, use = 'p')
adj_table = raw_cors$p
ad_table = as.data.frame(adj_table[order(ordered(row.names(adj_table), levels = gene_key)),])
colnames(ad_table) = 'female_fasttwitch'
heat_table$fem_fastwticch = ad_table$female_fasttwitch[match(row.names(heat_table), row.names(ad_table))]

gpx_male = as.data.frame(working_dataset$`CES4A_Muscle - Skeletal`)
row.names(gpx_male) = row.names(working_dataset)
colnames(gpx_male) = 'GPX3'
gp3 = gpx_male[row.names(gpx_male) %in% row.names(mm1),]

resids = summary(lm(mm1$fast_twitch_glycolytic_fiber ~ gp3))$residuals
raw_cors = bicorAndPvalue(mm2, resids, use = 'p')
adj_table = raw_cors$p
ad_table = as.data.frame(adj_table[order(ordered(row.names(adj_table), levels = gene_key)),])
colnames(ad_table) = 'male_fasttwitch_gpx3_adj'
heat_table$female_fasttwitch_ces4a = ad_table$male_fasttwitch_gpx3_adj[match(row.names(heat_table), row.names(ad_table))]

#Now all of the results from fast-twitch fiber composition correlated to th top genes in both sexes and +/- residuals for myokines are present in the same matix
heat_table1 = as.matrix(heat_table)
#set a pvalue threshold
heat_annot = ifelse(heat_table1<1e-6, '*', '')
#plot
pheatmap(heat_table, fontsize_number = 20, display_numbers = heat_annot, number_color = "white", color = colormap(colormap = colormaps$freesurface_red, nshades = 50), cluster_rows = F, cluster_cols = F,  main='Decon cell proportions AR and ER Males', fontsize_row = 5, fontsize_col = 5)

