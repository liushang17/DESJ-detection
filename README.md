# DSJ-detection:Differential splicing junction Detection between two groups at single cell level

## About DSJ-Detection


A computation framework designed to detect differential junctions usage between two groups at single cell level. The software have three highlights:
1. iterative K-means to discern cells failing to reflect the real junctions count distribution of each gene due to low expression and technical noise. 
2. a new normalization method, which normalize the read count per junction of each gene using the reads count aligned to all junctions of the gene rather than uniquely mapped reads in each cell. 
3. DSJ-detection can detect differential junction in all positions of genes, so it can discover any pattern of alternative splicing

