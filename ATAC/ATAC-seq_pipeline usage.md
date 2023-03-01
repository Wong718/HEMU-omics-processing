# ATAC-seq_pipeline usage

## Overview of the pipeline

### Language

Python 3, Python 2.7

### Original input files:

ATAC-seq data formatting in .sra, .fastq or .fastq.gz

### Functions of the pipeline

The packages used in the pipline should be installed initially through `conda install`.

| Function | Description | Packages | Environment |
| --- | --- | --- | --- |
| mvfile | move files with specific filetype into a directory |  | Python 3 |
| extsra | transfer the .sra files into the .fastq files | fastq-dump | Python 3 |
| fastp | remove adaptor and filter the .fastq or .fastq.gz files | fastp | Python 3 |
| renm | remove the ‘_Filtered’ label of the filtered files |  | Python 3 |
| bowtie | align the ATAC-seq data to the reference genome | bowtie2
samtools | Python 3 |
| rmdup | remove the repeated sequencing and transfer the .bam files into .bw files | bedtools
samtools
deeptools | Python 3 |
| TSS | TSS enrichment and visulization | deeptools | Python 3 |
| Peak | Peak calling | macs2 | Python 2 |
| tn5 | Tn5 filtering | deeptools | Python3, R |

### Usage

1. Process the following command.
    
    ```bash
    python ATAC-seq_pipeline.1.0.0.py [para1] <function> [para2] [para3] ...
    ```
    
2. Check the output line and process the output script.
    
    ```bash
    nohup sh Command*.sh &
    ```
    

## Function: `mvfile`

### Input

specific formating files in the working directory

### Output

specific formating files in the target directory

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] mvfile [file format]
```

### Parameters

- target_directory: The target directory for specific formating files
- file format: The ending of the specific files, e.g. ‘.sra’

## Function: `extsra`

### Input

.sra files

### Output

.fastq.gz files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] extsra [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .sra files
- parallel_task_number: the number of parallel tasks on the server

## Function: `fastp`

### Input

.fastq or fastq.gz files

### Output

.fastq_Filtered.gz, _Report.html, _Report.json files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] fastp [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .fastq or fastq.gz files
- parallel_task_number: the number of parallel tasks on the server

## Function: `renm`

### Input

.fastq_Filtered.gz files

### Output

.fastq.gz files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] renm
```

### Parameters

- target_directory: The target directory for .fastq_Filtered.gz files

## Function: `bowtie`

### Input

.fastq.gz files

### Output

.sort.bam and .sort.bam.bai files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] bowtie [parallel_task_number] [ref]
```

### Parameters

- target_directory: The target directory for .fastq.gz files
- parallel_task_number: the number of parallel tasks on the server
- ref: the path of the reference genome after bowtie2 indexing
    
    ```python
    ref_genome_zm_V4_path = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/B73V4"
    ref_genome_sor_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sor_ref/sbi"
    
    if sys.argv[4] == "zmv4":
        ref_genome_path = ref_genome_zm_V4_path
    elif sys.argv[4] == "sor":
        ref_genome_path = ref_genome_sor_path
    else:
        ref_genome_path = sys.argv[4]
    ```
    

### Attention

- For comprehensive alignment, the reference genome and relative annotation should be downloaded from *ensemble plant*.
- If your target species are not included above, the reference genome and bowtie2 indexing should be prepared and add your path as the parameter.
    
    ```python
    # Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz downloaded from ensemblplants
    gunzip Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz
    nohup bowtie2-build B73.fa B73
    ```
    
- For effective downstream analysis, the chromosome name of the .fa and the .gff3 files should be the same. Use `sed -i` to rename.

## Function: `rmdup`

### Input

.sort.bam files

### Output

.sort.rmdup.bam, .sort.rmdup.bed and .sort.rmdup.bw files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] rmdup [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .sort.bam files
- parallel_task_number: the number of parallel tasks on the server

## Function: `TSS`

### Input

.bw files

### Output

_TSS.gz, _Profile.png files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] TSS [parallel_task_number] [ref] [up_stream] [down_stream]
```

### Parameters

- target_directory: The target directory for .fastq or fastq.gz files
- parallel_task_number: the number of parallel tasks on the server
- ref: the path of the reference genome TSS extraction txt file
    
    ```python
    ref_genome_zm_V4_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/ZmV4/zmv4_refseq_genes_TSS.txt"
    ref_genome_sor_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/sorghum/Sbi_refseq_genes_TSS.txt"
    
    if sys.argv[4] == "zmv4":
        ref_genome_path = ref_genome_zm_V4_path
    elif sys.argv[4] == "sor":
        ref_genome_path = ref_genome_sor_path
    else:
    		ref_genome_path = sys.argv[4]
    ```
    
    - If your target species are not included above, the reference genome TSS extraction txt file should be prepared and add your path as the parameter. The above txt file can be generated by `gfftobed` tools.
    
    ```python
    # step1. download and install the gfftobed tools. link: https://github.com/jacobbierstedt/gfftobed
    # step2. process with the following bash script
    nohup ../gfftobed-main/gfftobed -t [reference genome gff path] | cut -f 1,2,3 | grep chr > ./[reference genome gff path]_refseq_genes_TSS.txt &
    ```
    
    - the above script may go wrong for the chromosome names in .gff file are not started with chr
- up_stream/down_stream: the range of the TSS region and 3000 is recommended.

## Function: `Peak`

### Input

.bed files

### Output

.narrowPeak, .log files

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] peak [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .sort.bam files
- parallel_task_number: the number of parallel tasks on the server

### Attention

conda environment should change from python3 to python2, which has installed macs2

## Function: `tn5`

### Input

.rmdup.bed, .narrowPeak, .fa.fai, files

### Output

.narrowPeak files (after tn5 filtering, locating in a new ‘filter_OCR’ file

### Usage

```bash
python ATAC-seq_pipeline.1.0.0.py [target_directory] tn5 [parallel_task_number] [genome_fai] [fdr]
```

### Parameters

- target_directory: The target directory for .sort.bam files
- genome_fai: the .fa.fai file of the reference genome

```bash
ref_genome_zm_V4_path = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa.fai"
ref_genome_sor_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sor_ref/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai"

if sys.argv[4] == "zmv4":
    ref_genome_path = ref_genome_zm_V4_path
elif sys.argv[4] == "sor":
    ref_genome_path = ref_genome_sor_path
else:
    ref_genome_path = sys.argv[4]
```

- fdr: the fdr setting, recommending 0.01

### Attention

- R base is needed for this function and it can be easily install in the conda environment by `conda install -c conda-forge r-base`
- Two script should add to the **target directory**
    - filter_ocr_by_tn5_density.sh
        
        ```bash
        #!/bin/bash
        ##usage example: bash filter_ocr_by_tn5_density.sh B73v4.fa.fai 0.01
        ls *sort.rmdup.bed |while read bed
        do
        sample_name=$(echo $bed |cut -d "." -f 1)
        input_narrowpeak=${sample_name}.sort.rmdup_peaks.narrowPeak
        genome_fai=$1
        fdr=$2
        ##提取reads的Tn5 site
        awk 'BEGIN {OFS = "\t"} ; {print $1,($2+$3)/2,($2+$3)/2+1}' $bed > $bed.Tn5
        ##计算每个ocr的tn5 site densify（tn5 count/ ocr length)
        bedtools coverage -a $input_narrowpeak -b $bed.Tn5 | awk '{print $1,$2,$3,$4,$11}' > input_${sample_name}_peaks_and_tn5count.bed
        ##bedtools shuffle specifically excluding ACRs from the randomized selection space
        bedtools shuffle -i $input_narrowpeak -g $genome_fai |awk '{print $1,$2,$3}' | tr ' ' '\t' > ${sample_name}.shuffle_randomized.bed
        ##计算randomized区域的tn5 density
        bedtools coverage -a ${sample_name}.shuffle_randomized.bed -b $bed.Tn5 | awk '{print $1,$2,$3,$4}' > input_${sample_name}_randomized_and_tn5count.bed #四列 chr start end tn5count
        
        Rscript ./filter_eFDR.R input_${sample_name}_peaks_and_tn5count.bed input_${sample_name}_randomized_and_tn5count.bed $fdr ${sample_name}_peaks_filted_by_tn5density.fdr$fdr.bed
        
        done
        
        ##move filter OCR to another dir
        # ls *_filted_by_tn5density.fdr0.01.bed |while read id
        # do
        # name=$(echo $id |cut -d "_" -f 1)
        # cp $id filter_OCR/${name}_peaks.narrowPeak
        # done
        
        ###将pvalue等信息加上
        ls *_filted_by_tn5density.fdr0.01.bed |while read id
        do
        name=$(echo $id |cut -d "_" -f 1)
        awk '{print $4}' $id >$id.fdr0.01.idex
        grep -w -f  $id.fdr0.01.idex ${name}_peaks.narrowPeak >filter_OCR/${name}_peaks.narrowPeak
        done%cd 
        ```
        
    - filter_eFDR.R
        
        ```r
        # load data
        args <- commandArgs(T)
        true.acrs <- as.character(args[1])
        fake.acrs <- as.character(args[2])
        fdr <- as.numeric(args[3])
        output <- as.character(args[4])
        
        # import
        a <- read.table(true.acrs)
        b <- read.table(fake.acrs)
        
        # estimate Tn5 density per 1000bp
        a$tn5_density <- a[,ncol(a)]/((a$V3-a$V2)/1000)
        b$tn5_density <- b[,ncol(b)]/((b$V3-b$V2)/1000)
        
        # find emprical threshold
        threshold <- quantile(b$tn5_density, 1-fdr)
        
        # filter ACRs
        a.filtered <- subset(a, a$tn5_density > threshold)
        a.filtered$tn5_density <- NULL
        a.filtered[,c(1:(ncol(a.filtered)-1))]
        write.table(a.filtered, file=output, quote=F, row.names=F, col.names=F, sep="\t")
        
        ##usage:    >Rscript filter_eFDR.R <input_ACRs.bed> <controls.bed> <FDR> <output_filename.bed>
        ```