# miRNA-seq_pipeline usage

## Overview of the pipeline

### Language

Python 3

### Original input files:

miRNA-seq data formatting in .sra, .fastq or .fastq.gz

### Functions of the pipeline

The packages used in the pipline should be installed initially through `conda install`.

| Function | Description | Packages | Environment |
| --- | --- | --- | --- |
| mvfile | move files with specific filetype into a directory |  | Python 3 |
| extsra | transfer the .sra files into the .fastq files | fastq-dump | Python 3 |
| fastp | remove adaptor and filter the .fastq or .fastq.gz files | fastp | Python 3 |
| renm | remove the ‘_Filtered’ label of the filtered files and unzip the .gz files |  | Python 3 |
| merge | merge the double end .fastq file into a single file | pandaseq | Python 3 |
| premiRNA | preprocess the sequencing data before quantification | bowtie
mirdeep2 | Python 3 |
| qualify | align miRNA-seq to the miRNA database and qualify their expression | mirdeep2 | Python 3 |
| expression | integrate the .csv expression data into an expression matrix | pandas | Python 3 |

### Usage

1. Process the following command.
    
    ```bash
    python miRNA_pipeline.1.0.0.py [para1] <function> [para2] [para3] ...
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
python miRNA_pipeline.1.0.0.py [target_directory] mvfile [file format]
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
python miRNA_pipeline.1.0.0.py [target_directory] extsra [parallel_task_number]
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
python miRNA_pipeline.1.0.0.py [target_directory] fastp [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .fastq or fastq.gz files
- parallel_task_number: the number of parallel tasks on the server

## Function: `renm`

### Input

.fastq_Filtered.gz files

### Output

.fastq files

### Usage

```bash
python miRNA_pipeline.1.0.0.py [target_directory] renm [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .fastq_Filtered.gz files

## Function: `merge`

### Input

_1.fastq, _2.fastq files

### Output

.fastq files

### Usage

```bash
python miRNA_pipeline.1.0.0.py [target_directory] merge [parallel_task_number]
```

### Parameters

- target_directory: The target directory for _1.fastq, _2.fastq files
- parallel_task_number: the number of parallel tasks on the server

## Function: `premiRNA`

### Input

.fastq files

### Output

_collapsed.fa, _vs_refdb.arf

### Usage

```bash
python miRNA_pipeline.1.0.0.py [target_directory] premiRNA [parallel_task_number] [ref]
```

### Parameters

- target_directory: The target directory for .fastq files
- parallel_task_number: the number of parallel tasks on the server
- ref: the path of the reference genome and the bowtie indexes
    
    ```python
    ref_genome_zm_V4_path_new = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/B73V4"
    ref_genome_zm_V5_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Zmsubsp/S"
    ref_genome_coixlac_path = "/data5/zhuyuzhi/Reference_Genomes/Coix_New/Coix_lacryma/CoixLac"
    ref_genome_MF_path = "/data5/RNA_Seq_Database/wzj_reference_genome/MF_ref/Mlu_HiC_tran"
    ref_genome_sac_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sac_ref/sac_tran"
    ref_genome_sor_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/sorghum/Sorghum_bowtie.fa"
    ref_genome_Chy_ser_path = "/data5/Andropogoneae_LTR_Project/Step0_Genome_Download/GCA_015844335.1/GCA_015844335.1.fa"
    ref_genome_ssp_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Ssp/Saccharum_spo_bowtie.fa"
    ref_genome_shc_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Shc/Shc_bowtie.fa"
    ref_genome_mfl_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Miscanthus/Mfl/Mfl_bowtie.fa"
    ref_genome_mlu_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Miscanthus/Mlu/Mlu_bowtie.fa"
    ref_genome_msa_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/Miscanthus/Msa/Msa_bowtie.fa"
    
    if sys.argv[4] == "zmv4":
        ref_genome_path = ref_genome_zm_V4_path_new
    elif sys.argv[4] == "zmsub":
        ref_genome_path = ref_genome_zm_V5_path
    elif sys.argv[4] == "coix":
        ref_genome_path = ref_genome_coixlac_path
    elif sys.argv[4] == "MF":
        ref_genome_path = ref_genome_MF_path
    elif sys.argv[4] == "sac":
        ref_genome_path = ref_genome_sac_path
    elif sys.argv[4] == "sor":
        ref_genome_path = ref_genome_sor_path
    elif sys.argv[4] == "Chy_ser":
        ref_genome_path = ref_genome_Chy_ser_path
    elif sys.argv[4] == "ssp":
        ref_genome_path = ref_genome_ssp_path
    elif sys.argv[4] == 'shc':
        ref_genome_path = ref_genome_shc_path
    elif sys.argv[4] == 'Mfl':
        ref_genome_path = ref_genome_mfl_path
    elif sys.argv[4] == 'Mlu':
        ref_genome_path = ref_genome_mlu_path
    elif sys.argv[4] == 'Msa':
        ref_genome_path = ref_genome_msa_path
    else:
        ref_genome_path = sys.argv[4]
    ```
    

### Attention

The reference genome shoube be indexed by bowtie-build (not bowtie2 !!)

```bash
bowtie-build refdb.fa refdb.fa
```

## Function: `qualify`

### Input

.fa files

### Output

.csv files

### Usage

```bash
python miRNA_pipeline.1.0.0.py [target_directory] premiRNA [parallel_task_number] [ref]
```

### Parameters

- target_directory: The target directory for .fa files
- parallel_task_number: the number of parallel tasks on the server
- ref: the species of the miRNA-seq; only ‘sbi’ (sorghum) and ‘zma’ (maize) is avaliable in the pipeline
    - the microRNA of the specific species should be extrated initially following the steps
    1. Dowload the microRNA reference sequence from miRbase ([https://www.mirbase.org/ftp.shtml](https://www.mirbase.org/ftp.shtml))
    2. Extract the sequence of the target species and place them in ‘ref’ directory under the working directory
    
    ```python
    # extract mature
    extract_miRNAs.pl mature.fa.gz zma > mature_ref.fa
    # extract hairpin
    extract_miRNAs.pl hairpin.fa.gz zma > hairpin_ref.fa
    ```
    

## Function: `expression`

### Input

.csv files

### Output

Expression_miRNA.csv

### Usage

```bash
python miRNA_pipeline.1.0.0.py [target_directory] expression 
```

### Parameters

- target_directory: The target directory for .csv files

### Attention

The pandas packages should be installed in the Python3 environment initially.