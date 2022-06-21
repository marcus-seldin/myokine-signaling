# Human myokine signaling study overview
This repository provides all scripts and a detailed walk-through to repeat analyses shown in Velez, Van et al. 2022 (https://elifesciences.org/articles/76887)

# Analyses are broken down into 3 sections:
Myokine differential expression by sex and muscle-specific enrichment (Fig 1)

Cross-tissue signaling for myokines in the context of sex and hormones (Fig 2)

Generation of pseudo-single-cell muscle maps and cross-tissue regressions (Fig 3)

Each analysis should operate independently and all use the datasets provided below.

# All datasets (pre-processed) used in this study are available here:
https://drive.google.com/drive/folders/1YKT8lkGzGVFk74CqS5FmT6lIH0ZVfgfJ?usp=sharing

## Publicly available data used (above) were acquired from:

GTEx V8: https://gtexportal.org/home/datasets (the raw data were filtered where individuals were required to show counts > 0 in 1.2e6 gene_tissue combinations across all data.  Given that our goal was to look across tissues at enrichments, this was done to limit spurious influence of genes only expressed in specific tissues.  Post-filtering consists of 310 individuals and 1.8e7 gene_tissue combinations). This results in 52% NA values in the final matrix used. 

Uniprot annotations for secreted proteins in humans: https://www.uniprot.org/uniprot/?query=locations%3A%28location%3A%22Secreted+%5BSL-0243%5D%22+type%3Acomponent%29+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score

Uniprot annotations for secreted proteins in mice: https://www.uniprot.org/uniprot/?query=locations:(location:%22Secreted%20[SL-0243]%22%20type:component)&fil=organism%3A%22Mus+musculus+%28Mouse%29+%5B10090%5D%22&sort=score

MGI list of mouse-human orthologs: http://www.informatics.jax.org/homology.shtml

Human skeletal muslce sc-seq data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130977

Mouse muscle-specific Esr1 KO (MERKO) and WT skeletal muslce RNA-Seq were generated as part of this study and deposited in Sequence Read Archive (SRA) under the accession: PRJNA785746

# Additional Considerations
While the tissues shown in this study capture notable metabolic processes, all available tissues are included in the filtered dataset to enable exapnsion of the analysis.  Specifically, these can be identified as columns in 'working_dataset'.  The scripts cna be easily repurposed to incorporate other tissues (for example, coronary artery)

We apologize for bulky and recursive portions of the scripts.  Realizing that there and more efficent mechanisms of executing these analyses, we felt that breaking down step-by-step is helpful to facilitate learning.  

All feedback is welcome, please email: lmvelez@uci.edu & mseldin@uci.edu

