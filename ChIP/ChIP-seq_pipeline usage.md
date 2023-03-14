# ChIP-seq_pipeline usage

## Overview of the pipeline

### Language

Python 3

### Original input files:

ChIP-seq data formatting in .sra, .fastq or .fastq.gz

### Functions of the pipeline

The packages used in the pipline should be installed initially through `conda install`.

| Function | Description | Packages | Environment |
| --- | --- | --- | --- |
| mvfile | move files with specific filetype into a directory |  | Python 3 |
| extsra | transfer the .sra files into the .fastq files | fastq-dump | Python 3 |
| fastp | remove adaptor and filter the .fastq or .fastq.gz files | fastp | Python 3 |
| renm | remove the ‘_Filtered’ label of the filtered files and unzip the .gz files |  | Python 3 |
| bowtie | align ChIP-seq to the reference genome | bowtie
samtools | Python 3 |
| picard | .bam files dedulplicate | picard | Python 3 |
| Peak | call peak from the dedulplicated .bam files | macs2 | Python 2 |

### Usage

1. Process the following command.
    
    ```bash
    python ChIP_pipeline.1.0.0.py [para1] <function> [para2] [para3] ...
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
python ChIP_pipeline.1.0.0.py [target_directory] mvfile [file format]
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
python ChIP_pipeline.1.0.0.py [target_directory] extsra [parallel_task_number]
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
python ChIP_pipeline.1.0.0.py [target_directory] fastp [parallel_task_number]
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
python ChIP_pipeline.1.0.0.py [target_directory] renm [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .fastq_Filtered.gz files

## Function: `bowtie`

### Input

fastq or fastq.gz files

### Output

.bam files

### Usage

```bash
python ChIP_pipeline.1.0.0.py [target_directory] bowtie [parallel_task_number] [ref]
```

### Parameters

- target_directory: The target directory for fastq or fastq.gz files
- parallel_task_number: the number of parallel tasks on the server
- ref: path of the reference genome and the bowtie2 genome index
    
    ```python
    ref_genome_zm_V4_path = "/data21/wongzj/Other_Seq/miRNA-seq/bowtie/ZmV4/ZmV4_bowtie.fa"
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
        ref_genome_path = ref_genome_zm_V4_path
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

The index of the reference genome should be generated initially by `bowtie2-build`

## Function: `picard`

### Input

.bam files

### Output

.bam_deduplicate_bam files

### Usage

```bash
python ChIP_pipeline.1.0.0.py [target_directory] picard [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .bam files
- parallel_task_number: the number of parallel tasks on the server

## Function: `Peak`

### Input

.bam_deduplicate_bam files

### Output

_model.r, _peaks.narrowPeak, _peak.xls, _summits.bed, *treat*pileup.bdg

### Usage

```bash
python ChIP_pipeline.1.0.0.py [target_directory] Peak [parallel_task_number]
```

### Parameters

- target_directory: The target directory for .fa files
- parallel_task_number: the number of parallel tasks on the server

### Attention

The commands generated by the scripts would change the environment to a Python 2 environment called ‘`CHIP_analysis`'.

Please create the conda environment in Python 2 initially.

```bash
conda create -n CHIP_analysis python==2
```