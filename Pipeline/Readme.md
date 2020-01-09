# DESJ-pipeline

This pipeline is developed for the preparation of input files of DESJ-detection.

## merge.sj.r

### Description
  This script input path of the splicing junction information file(SJ.out.tab) outputed by STAR and the threshold of read count and cell number. It would output splicing junction list.
 
### Usage

Usage: merge.sj.r [options]
Options:

        --sj_dir=CHARACTER
                The path of splicing junction outputed by STAR

        --outdir=CHARACTER
                The pathway of outdir

        --min_read=NUMERIC
                The read count threshold

        --min_cell=NUMERIC
                The cell count threshold

        --cpu=NUMERIC
                The cpu number

        -h, --help
                Show this help message and exit

### Output: The junction list

chr21_33547906_33549256_1

chr2_151692363_151694322_2

chr3_180985997_180986997_2

chr3_37832939_37834678_2

chr11_125653848_125655224_1

chr3_123740002_123752330_2

## Junction.ann.pl
### Descript
This script needs four inputs: the junction list (as output by merge.sj.r), the path of directory containing SJ.out.tab file, the genome annotation file (gtf), the output directory.
This will output four files: two Junction annotation files (Alljunction.filter.list.ann.xls, Alljunction.filter.list.ann.onegene.xls), a junction cell count matrix (merge.count.txt)

### Usage
perl Junction.ann.pl

Describe:
        Annotation junction with gene names, filtering junctions, creating juntion cell count matrix.
Author:
        liushang@genomics.cn
Usage:
        perl  Junction.ann.pl  <in.gtf>  <Alljunction.filter.list.xls>  <sj.dir>  <Outdir>

### Output
Alljunction.filter.list.ann.xls:
  Junction  GeneID    GeneName
  
  chr2_151692363_151694322_2      ENSG00000183091.19      NEB
  
  chr3_180985997_180986997_2      ENSG00000205981.6       DNAJC19
  
  chr3_37832939_37834678_2        ENSG00000235257.8       ITGA9-AS1
  
  chr11_125653848_125655224_1     ENSG00000149554.12      CHEK1
  
  chr3_123740002_123752330_2      ENSG00000065534.18      MYLK

Alljunction.filter.list.ann.onegene.xls:
  Junction  GeneID    GeneName
  
  chr2_151692363_151694322_2      ENSG00000183091.19      NEB
  
  chr3_180985997_180986997_2      ENSG00000205981.6       DNAJC19
  
  chr3_37832939_37834678_2        ENSG00000235257.8       ITGA9-AS1

merge.count.txt:
cell    chr2_151692363_151694322_2      chr3_180985997_180986997_2      chr3_37832939_37834678_2

P08_T_0976      0       0       0

P15_A_0065      0       0       0

P09_T_0669      0       0       0

P09_T_1398      0       2       0
