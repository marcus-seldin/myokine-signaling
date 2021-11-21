# Human myokine signaling study overview
This repository provides all scripts and a detailed walk-through to repeat analyses shown in Velez, Van et al. 2021

# Analyses are broken down into 3 sections:
1: Myokine differential expression by sex and muslce-specific enrichment (Fig 1)

2: Cross-tissue signaling for myokines in the context of sex and hormones (Fig 2)

3: Generation of pseudo-single-cell muscle maps and cross-tissue regressions (Fig 3)

# Datasets (pre-processed) used in this study are available here:
https://drive.google.com/drive/folders/1YKT8lkGzGVFk74CqS5FmT6lIH0ZVfgfJ?usp=sharing

## Publicly available data used here were acquired from:

GTEx V8: https://gtexportal.org/home/datasets

Uniprot annotations for secreted proteins in humans: https://www.uniprot.org/uniprot/?query=locations%3A%28location%3A%22Secreted+%5BSL-0243%5D%22+type%3Acomponent%29+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score

Uniprot annotations for secreted proteins in mice: https://www.uniprot.org/uniprot/?query=locations:(location:%22Secreted%20[SL-0243]%22%20type:component)&fil=organism%3A%22Mus+musculus+%28Mouse%29+%5B10090%5D%22&sort=score

MGI list of mouse-human orthologs: http://www.informatics.jax.org/homology.shtml

Human skeletal muslce sc-seq data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130977

Mouse muscle-specific Esr1 KO (MERKO) and WT skeletal muslce RNA-Seq were generated as part of this study and deposited in: 
