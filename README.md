# DSJ-detection:Differential splicing junction Detection between two groups at single cell level

About DSJ-Detection


A computation framework designed to detect differential junctions usage between two groups at single cell level. The software applies iterative K-means to discern cells failing to reflect the real junctions count distribution of each gene due to low expression and technical noise. Then, the second problem is solved with a new normalization method, which normalize the read count per junction of each gene using the reads count aligned to all junctions of the gene rather than uniquely mapped reads in each cell. Meanwhile, DSJ-detection can detect differential junction in all positions of genes, so it can discover any pattern of alternative splicing
