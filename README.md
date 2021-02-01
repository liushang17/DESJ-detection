# DESJ-detection:Differential splicing junction Detection between two groups at single cell level

## About DESJ-Detection


A computation framework designed to detect differential junctions usage between two groups at single cell level. The software have three highlights:
1. iterative K-means to discern cells failing to reflect the real junctions count distribution of each gene due to low expression and technical noise. 
2. a new normalization method, which normalize the read count per junction of each gene using the reads count aligned to all junctions of the gene rather than uniquely mapped reads in each cell. 
3. DESJ-detection can detect differential junction in all positions of genes, so it can discover any pattern of alternative splicing.

## Requirements
optparse
ggplot2
parallel
pheatmap
edgeR
limma

## Release

### Release v2.0.5
Initial release of DESJ-detection

### Usage
Just download it. To check that it can run perfectly, try the command with the help option (-h), Rscript DESJ-detection.v0.1.r  -h. The output should look like this:

Usage: DESJ-detection.v0.1.r [options]

Options:


        --junction_matrix=CHARACTER
                The junction_cell count matrix

        --col_ann=CHARACTER
                The text file, containing the clustering info, possess two columns, including cell name and clusters name

        --min_pct=NUMERIC
                The minimum scale in either of the two groups of cells, is able to detect the junctions of any genes

        --cellinfo=CHARACTER
                The text file, containing the uniquely mapped reads number, possess two colums, including cell name and
                reads number

        --row_ann=CHARACTER
                The text file, containing the junctions info, including three columns, including junction id, gene id
                and gene symbol

        --clus=CHARACTER
                The name of the two clusters, will be detected for differential junctions. for example: clusterA,clusterB

        --outdir=CHARACTER
                The pathway of the outdir

        --mincell=NUMERIC
                The minimum cell number of the cluster determined by iterative k-means would be filtered

        --maxsd=NUMERIC
                The maximum standard deviation of the expression of the junctions of the cells in the
                filtered cluster determined by iterative k-means

        --maxmean=NUMERIC
                The maximum mean of the expression of the junctions of the cells in the
                filtered cluster determined by iterative k-means

        --cpu=NUMERIC
                The cpu number

        -h, --help
                Show this help message and exit

## About Plot.cell.filter.r
This would output a heatmap file, showing the expression distribution of all junction of the selected genes across cell clusters determined by iterative K-means. 

### Usage
To check that it can run perfectly, try the command with the help option (-h), Rscript Plot.cell.filter.r  -h. The output should look like this:

Usage: Plot.cell.filter.r [options]
Options:
        --junction_matrix=CHARACTER
                The junction_cell count matrix

        --col_ann=CHARACTER
                The text file, containing the clustering info, possess two columns, including cell name and clusters name

        --min_pct=NUMERIC
                The minimum scale in either of the two groups of cells, is able to detect the junctions of any genes

        --cellinfo=CHARACTER
                The text file, containing the uniquely mapped reads number, possess two colums, including cell name and
                reads number

        --row_ann=CHARACTER
                The text file, containing the junctions info, including three columns, including junction id, gene id
                and gene symbol

        --clus=CHARACTER
                The name of the two clusters, will be detected for differential junctions. for example: clusterA,clusterB

        --outdir=CHARACTER
                The pathway of the outdir

        --mincell=NUMERIC
                The minimum cell number of the cluster determined by iterative k-means would be filtered

        --maxsd=NUMERIC
                The maximum standard deviation of the expression of the junctions of the cells in the
                filtered cluster determined by iterative k-means

        --maxmean=NUMERIC
                The maximum mean of the expression of the junctions of the cells in the
                filtered cluster determined by iterative k-means

        --cpu=NUMERIC
                The cpu number

        --genename=CHARACTER
                The name of gene is used to show how to filter cells

        -h, --help
                Show this help message and exit

