import os, sys, shutil, random

command_list = []


def rename_file(directory):
    file_list = []
    target_fq_list = []
    target_fqgz_list = []
    target_fq_single_list = []
    target_fqgz_single_list = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for file_dir in file_list:
        if file_dir.endswith("_1_Filtered.fastq") or file_dir.endswith("_2_Filtered.fastq"):
            target_fq_list.append(file_dir)
            continue
        if file_dir.endswith("_1_Filtered.fastq.gz") or file_dir.endswith("_2_Filtered.fastq.gz"):
            target_fqgz_list.append(file_dir)
            continue
        if file_dir.endswith("_Filtered.fastq"):
            target_fq_single_list.append(file_dir)
            continue
        if file_dir.endswith("_Filtered.fastq.gz"):
            target_fqgz_single_list.append(file_dir)
            continue

    for file_dir in target_fq_list:
        new_name = os.path.split(file_dir)[0] + "/" + os.path.split(file_dir)[1].split("_")[0] + "_" + \
                   os.path.split(file_dir)[1].split("_")[1] + ".fastq"
        print("%s -> %s" % (file_dir, new_name))
        os.rename(file_dir, new_name)

    for file_dir in target_fqgz_list:
        new_name = os.path.split(file_dir)[0] + "/" + os.path.split(file_dir)[1].split("_")[0] + "_" + \
                   os.path.split(file_dir)[1].split("_")[1] + ".fastq.gz"
        print("%s -> %s" % (file_dir, new_name))
        os.rename(file_dir, new_name)

    for file_dir in target_fq_single_list:
        new_name = os.path.split(file_dir)[0] + "/" + os.path.split(file_dir)[1].split("_")[0] + ".fastq"
        print("%s -> %s" % (file_dir, new_name))
        os.rename(file_dir, new_name)

    for file_dir in target_fqgz_single_list:
        new_name = os.path.split(file_dir)[0] + "/" + os.path.split(file_dir)[1].split("_")[0] + ".fastq.gz"
        print("%s -> %s" % (file_dir, new_name))
        os.rename(file_dir, new_name)


def copy_to_wd(directory, file_type):
    file_list = []
    target_file_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(str(file_type)):
            target_file_list.append(filedir)

    for file_dir in target_file_list:
        shutil.move(file_dir, os.getcwd())


def extract_sra(directory):
    global filedir
    file_list = []
    sra_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        # if filedir.endswith(".sra"):
        #     sra_list.append(filedir)
        if filedir.split('/')[1].startswith("SRR"):
            sra_list.append(filedir)

    outdir = os.getcwd() + "/FASTQs_BGI_s/"
    for sra_file in sra_list:
        tmpinst = "fastq-dump --gzip -outdir %s %s &" % (outdir, sra_file)
        # os.system(tmpinst)
        # time.sleep(sra_extract_interval)
        command_list.append(tmpinst)


def fastp_zip(directory):
    file_list = []
    fastq_list = []

    double_edge_seq_fastq_comp = []
    single_edge_seq_fastq_comp = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".fastq.gz"):
            fastq_list.append(filedir)

    for fastq_file in fastq_list:
        if fastq_file.endswith("_1.fastq.gz") or fastq_file.endswith(
                "_2.fastq.gz"):  # Identify Double-Edge Sequenced File

            fastq_file = os.path.split(fastq_file)[0] + "/" + os.path.split(fastq_file)[1].split("_")[0]

            if fastq_file not in double_edge_seq_fastq_comp:
                double_edge_seq_fastq_comp.append(fastq_file)
        else:
            single_edge_seq_fastq_comp.append(fastq_file)

    for fastq_file in single_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        output_fq = fastq_file.rstrip(".fastq.gz") + "_Filtered.fastq.gz"
        output_json = fastq_file.rstrip(".fastq.gz") + "_Report.json"
        output_html = fastq_file.rstrip(".fastq.gz") + "_Report.html"

        tmpinst = 'fastp -i %s -o %s -j %s -h %s &' % (
            fastq_file, output_fq, output_json, output_html)
        # os.system(tmpinst)
        # time.sleep(fastp_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq.gz"
        input_fastq2 = fastq_file + "_2.fastq.gz"
        output_fastq1 = fastq_file + "_1_Filtered.fastq.gz"
        output_fastq2 = fastq_file + "_2_Filtered.fastq.gz"
        output_json = fastq_file + "_Report.json"
        output_html = fastq_file + "_Report.html"

        tmpinst = 'fastp -i %s -I %s -o %s -O %s -j %s -h %s &' % (
            input_fastq1, input_fastq2, output_fastq1, output_fastq2, output_json, output_html)
        # os.system(tmpinst)
        # time.sleep(fastp_interval)
        command_list.append(tmpinst)
        # print(tmpinst)


def fastp_unzip(directory):
    file_list = []
    fastq_list = []

    double_edge_seq_fastq_comp = []
    single_edge_seq_fastq_comp = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".fastq"):
            fastq_list.append(filedir)

    for fastq_file in fastq_list:
        if fastq_file.endswith("_1.fastq") or fastq_file.endswith(
                "_2.fastq"):  # Identify Double-Edge Sequenced File

            fastq_file = os.path.split(fastq_file)[0] + "/" + os.path.split(fastq_file)[1].split("_")[0]

            if fastq_file not in double_edge_seq_fastq_comp:
                double_edge_seq_fastq_comp.append(fastq_file)
        else:
            single_edge_seq_fastq_comp.append(fastq_file)

    for fastq_file in single_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        output_fq = fastq_file.rstrip(".fastq") + "_Filtered.fastq.gz"
        output_json = fastq_file.rstrip(".fastq") + "_Report.json"
        output_html = fastq_file.rstrip(".fastq") + "_Report.html"

        tmpinst = 'fastp -i %s -o %s -j %s -h %s &' % (
            fastq_file, output_fq, output_json, output_html)
        # os.system(tmpinst)
        # time.sleep(fastp_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq"
        input_fastq2 = fastq_file + "_2.fastq"
        output_fastq1 = fastq_file + "_1_Filtered.fastq.gz"
        output_fastq2 = fastq_file + "_2_Filtered.fastq.gz"
        output_json = fastq_file + "_Report.json"
        output_html = fastq_file + "_Report.html"

        tmpinst = 'fastp -i %s -I %s -o %s -O %s -j %s -h %s &' % (
            input_fastq1, input_fastq2, output_fastq1, output_fastq2, output_json, output_html)
        # os.system(tmpinst)
        # time.sleep(fastp_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def picard(directory):
    file_list = []
    bam_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".bam"):
            bam_list.append(filedir.rstrip(".bam"))

    for bam in bam_list:
        tmp_imput = bam + ".bam"
        tmp_output_gtf = bam + ".bam_deduplicate_bam"
        tmp_output_log = bam + ".log"
        tmpinst = "picard MarkDuplicates REMOVE_DUPLICATES=true I= %s O= %s M= %s &" % (
            tmp_imput, tmp_output_gtf, tmp_output_log)
        # os.system(tmpinst)
        # time.sleep(stringtie_interval)
        command_list.append(tmpinst)

def bowtie_sort_zip(directory, ref_genome_path):
    file_list = []
    fastq_list = []
    bam_list = []  # Used for filtering out Pre-Analyzed bams

    double_edge_seq_fastq_comp = []
    single_edge_seq_fastq_comp = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".fastq.gz"):
            fastq_list.append(filedir)
        if filedir.endswith(".bam"):
            bam_list.append(os.path.split(filedir)[1].rstrip(".bam"))

    for fastq_file in fastq_list:
        if fastq_file.endswith("_1.fastq.gz") or fastq_file.endswith(
                "_2.fastq.gz"):  # Identify Double-Edge Sequenced File
            fastq_file = os.path.split(fastq_file)[0] + "/" + os.path.split(fastq_file)[1].split("_")[0]
            if os.path.split(fastq_file)[1].split("_")[0] not in bam_list:
                if fastq_file not in double_edge_seq_fastq_comp:
                    double_edge_seq_fastq_comp.append(fastq_file)
        else:
            if os.path.split(fastq_file)[1].rstrip(".fastq.gz") not in bam_list:
                single_edge_seq_fastq_comp.append(fastq_file)

    for fastq_file in single_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        tmp_output_sam = fastq_file.rstrip(".fastq.gz") + "V4.sam"
        tmp_output_bam = fastq_file.rstrip(".fastq.gz") + "V4.bam"

        tmpinst = 'bowtie -q -v 2 -l 10 -k 15 %s %s -S %s && samtools sort -@ 4 -o %s %s && rm %s & '% (
            ref_genome_path, fastq_file, tmp_output_sam, tmp_output_bam, tmp_output_sam, tmp_output_sam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq.gz"
        input_fastq2 = fastq_file + "_2.fastq.gz"
        tmp_output_sam = fastq_file + "V4.sam"
        tmp_output_bam = fastq_file + "V4.bam"

        tmpinst = 'bowtie -q -v 2 -l 10 -k 15 %s -1 %s -2 %s -S %s && samtools sort -@ 4 -o %s %s && rm %s &' % (
            ref_genome_path, input_fastq1, input_fastq2, tmp_output_sam, tmp_output_bam, tmp_output_sam, tmp_output_sam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def bowtie_sort_unzip(directory,ref_genome_path):
    file_list = []
    fastq_list = []
    bam_list = []

    double_edge_seq_fastq_comp = []
    single_edge_seq_fastq_comp = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".fastq"):
            fastq_list.append(filedir)
        if filedir.endswith(".bam"):
            bam_list.append(os.path.split(filedir)[1].rstrip(".bam"))

    for fastq_file in fastq_list:
        if fastq_file.endswith("_1.fastq") or fastq_file.endswith(
                "_2.fastq"):  # Identify Double-Edge Sequenced File
            fastq_file = os.path.split(fastq_file)[0] + "/" + os.path.split(fastq_file)[1].split("_")[0]

            if os.path.split(fastq_file)[1].split("_")[0] not in bam_list:
                if fastq_file not in double_edge_seq_fastq_comp:
                    double_edge_seq_fastq_comp.append(fastq_file)
        else:
            if os.path.split(fastq_file)[1].rstrip(".fastq") not in bam_list:
                single_edge_seq_fastq_comp.append(fastq_file)

    for fastq_file in single_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        tmp_output_sam = fastq_file.rstrip(".fastq") + "V4.sam"
        tmp_output_bam = fastq_file.rstrip(".fastq") + "V4.bam"

        tmpinst = 'bowtie -q -v 2 -l 10 -k 15 %s %s -S %s && samtools sort -@ 4 -o %s %s && rm %s & ' % (
            ref_genome_path, fastq_file, tmp_output_sam, tmp_output_bam, tmp_output_sam, tmp_output_sam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq"
        input_fastq2 = fastq_file + "_2.fastq"
        tmp_output_sam = fastq_file + "V4.sam"
        tmp_output_bam = fastq_file + "V4.bam"

        tmpinst = 'bowtie -q -v 2 -l 10 -k 15 %s -1 %s -2 %s -S %s && samtools sort -@ 4 -o %s %s && rm %s &' % (
            ref_genome_path, input_fastq1, input_fastq2, tmp_output_sam, tmp_output_bam, tmp_output_sam, tmp_output_sam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def Peakcalling(directory):
    file_list = []
    bed_list = []
    command_list.append('conda activate CHIP_analysis')
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith("bed"):
            bed_list.append(filedir.rstrip(".bed"))

    for bed in bed_list:
        tmp_imput = bed + ".bed"
        tmp_log = bed + ".log"

        tmpinst = "macs2 callpeak -f BAM -t %s -n %s -g 2.1e9 --outdir ./peak/ --bdg -q 0.05 >%s 2>&1" % (
                      tmp_imput, tmp_imput, tmp_log)
        # os.system(tmpinst)
        # time.sleep(stringtie_interval)
        command_list.append(tmpinst)

def main():
    # USAGE
    # python Automaton.py [Target Directory] [ActionCommand] [SubsidiaryCommand/ParallelTaskNumber] [maize/coix/MF]
    global parallel_task_number
    if sys.argv[2] == "mvfile":  # [SubsidiaryCMD = Filetype(Example: .sra )]
        dir_main = sys.argv[1]
        file_type = sys.argv[3]
        copy_to_wd(dir_main, file_type)
    elif sys.argv[2] == "extsra":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        extract_sra(dir_main)
    elif sys.argv[2] == "fastp":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        fastp_zip(dir_main)
        fastp_unzip(dir_main)
    elif sys.argv[2] == "renm":  # Remove "_Filtered"
        dir_main = sys.argv[1]
        rename_file(dir_main)

    elif sys.argv[2] == "picard":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        picard(dir_main)

    elif sys.argv[2] == "bowtie":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]

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

        print("Using Ref %s" % ref_genome_path)

        bowtie_sort_zip(dir_main, ref_genome_path)
        bowtie_sort_unzip(dir_main, ref_genome_path)

    elif sys.argv[2] == "Peak":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        Peakcalling(dir_main)

    else:
        print("No Match ActionCommand.")

    try:
        sh_name = ("Command_List%s.sh" % random.randint(10000,99999))
        with open(sh_name, "w") as cmd_file:
            tmp_task_count = 0
            for inst in command_list:
                if str(tmp_task_count) == str(parallel_task_number):
                    cmd_file.write("wait\n")
                    tmp_task_count = 0
                cmd_file.write(inst + "\n")
                tmp_task_count += 1
        print("CMD Ready: Execute 'nohup sh %s &'" % sh_name)

    except IOError:
        print("CMD Write Failed: IO Error")


if __name__ == "__main__":
    main()
