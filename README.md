# Current strategy:

Aim: We want to use spatial scRNA to understand the development (driving genes and mechanisms) of the cancer. 

Background: It is believed that the cancer starts by a collection of unfortunate mutations in key cell mechanisms (for example they start multiplying way more, they don’t want to do apoptosis (kill themselves) when instructed by other cells and so on..). We are interested if and how are these mutations (copy number variants - CNVs) related to gene expression and mainly different biological mechanisms. To do that, we use spatial scRNA data. Main difference between scRNA and spatial is that we know the exact position of the spots, which is great, but instead of single cells, we measure RNA in spots which usually consist of 5-20 cells depending on the method.

Plan: 
1)	We start with the hest (https://github.com/mahmoodlab/HEST) database and focus on prostate cancer. This is a database that collects spatial scRNA data from multiple other papers. Pros: It is a lot of data which are structure in the same way. Cons: We only have basic information about the samples (e.g. we are missing study specific annotation (if the cell is malignant or not, what is a stage of cancer, what kind of cell type)) \
Scripts: \
python/HEST.ipynb to download data from HEST database \
python/h5ad_to_counts.ipynb to save the data in format that can be loaded by Scevan \



2)	We use scevan (https://www.nature.com/articles/s41467-023-36790-9) to find out which spots are normal and which are malignant. In the latter case, we cluster them based on the CNV profile using the same method. 
Scripts: \
R/scevan.R to load data saved in python/h5ad_to_counts.ipynb and run scevan  \

 
3)	We use stlearn (https://www.nature.com/articles/s41467-023-43120-6) to compute the pseudotime trajectories from normal spots, through different stages (CNV clusters) of cancer, to the most malignant spots and we find the genes that are the most associated with these trajectories.
 

4)	We use gseapy (https://academic.oup.com/bioinformatics/article/39/1/btac757/6847088) to infer biological mechanisms based on the 
 

5)	The next steps would be to do this on multiple datasets, find some shared biological mechanisms and try to validate it in other datasets.








# Old strategy:


# Spatial_cell_cell_communication in Cancers - The mighty gang

Try to identify the CCC change in cancer spatial transcriptomic data. First, try two methods in the lung cancer data (Delineating the dynamic evolution from preneoplasia to invasive lung adenocarcinoma by integrating single-cell RNA sequencing and spatial transcriptomics https://www.nature.com/articles/s12276-022-00896-9#Sec2):

Scrabian - Comparative analysis of cell–cell communication at single-cell resolution(https://www.nature.com/articles/s41587-023-01782-z). Designed for scRNA - Martin has the code run nicely, this one might move quicker
COMMOT - Screening cell–cell communication in spatial transcriptomics via collective optimal transport (https://www.nature.com/articles/s41592-022-01728-4)
Ideally, we would like have it as a quick project. Find something in the lung cancer paper using the CCC method, then replicate it in multiple different cancer ST datasets

=================== 

Below is Mikael's theory: Here, we hypothesized that the same principle would apply to premalignant cells. We constructed a spatially resolved network of all cells in a tumour. In this network each cell got a number depicting how many cells of the same type were in its vicinity. In other words, a central premalignant cell in a large cluster of premalignant cells got a high number, while a central premalignant cell in a small cluster got a medium number. A single premalignant cell surrounded by stroma and immune cell got a low number. We found that premalignant cells with high numbers were more homogenous and similar to malignant cells than those with low numbers. The high number cells only had stimulatory interactions, while low number cells had multiple mixed interactions. We ranked the stimulatory interactions and their URs. We showed that the top-scoring URs were associated with survival in TCGA, CPTAC and UKBB. Blocking a top UR killed premalignant clusters and tumor development in mice. In summary, the street gang principle applies to premalignant cells: Those with many similar cells in the vicinity successfully outcompete other premalignant cells and become malignant because of friendly support from its neighbours.
