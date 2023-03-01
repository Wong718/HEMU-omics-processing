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
        tmp_output_sam = fastq_file.rstrip(".fastq.gz") + ".sort.sam"
        tmp_output_bam = fastq_file.rstrip(".fastq.gz") + ".sort.bam"

        tmpinst = 'bowtie2 -p 8 --very-sensitive -X 2000 -x %s %s | samtools sort -O bam -@ 5 -o %s - && samtools index %s & '% (
            ref_genome_path, fastq_file, tmp_output_bam, tmp_output_bam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval) sor
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq.gz"
        input_fastq2 = fastq_file + "_2.fastq.gz"
        tmp_output_sam = fastq_file + ".sort.sam"
        tmp_output_bam = fastq_file + ".sort.bam"

        tmpinst = 'bowtie2 -p 8 --very-sensitive -X 2000 -x %s -1 %s -2 %s | samtools sort -O bam -@ 5 -o %s - && samtools index %s' % (
            ref_genome_path, input_fastq1, input_fastq2, tmp_output_bam, tmp_output_bam)
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
        tmp_output_sam = fastq_file.rstrip(".fastq") + "sort.sam"
        tmp_output_bam = fastq_file.rstrip(".fastq") + "sort.bam"

        tmpinst = 'bowtie2 -p 8 --very-sensitive -X 2000 -x %s %s | samtools sort -O bam -@ 5 -o %s - && samtools index %s & '% (
            ref_genome_path, fastq_file, tmp_output_bam, tmp_output_bam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

    for fastq_file in double_edge_seq_fastq_comp:
        input_fastq1 = fastq_file + "_1.fastq"
        input_fastq2 = fastq_file + "_2.fastq"
        tmp_output_sam = fastq_file + ".sort.sam"
        tmp_output_bam = fastq_file + ".sort.bam"

        tmpinst = 'bowtie2 -p 8 --very-sensitive -X 2000 -x %s -1 %s -2 %s | samtools sort -O bam -@ 5 -o %s - && samtools index %s &' % (
            ref_genome_path, input_fastq1, input_fastq2, tmp_output_bam, tmp_output_bam)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def rmdup(directory):
    file_list = []
    bam_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith("sort.bam"):
            bam_list.append(filedir.rstrip(".sort.bam"))

    for bam in bam_list:
        tmp_imput = bam + ".sort.bam"
        tmp_output_bam = bam + ".sort.rmdup.bam"
        tmp_output_bed = bam + ".sort.rmdup.bed"
        tmp_output_bw = bam +'.sort.rmdup.bw'

        tmpinst = "samtools rmdup %s %s &&" \
                  "samtools index %s &&" \
                  "bedtools bamtobed -i %s > %s &&" \
                  "bamCoverage -o %s -b %s -p 8 --normalizeUsing RPGC --effectiveGenomeSize 2100000000 --binSize 10 &" % (
            tmp_imput, tmp_output_bam, tmp_output_bam, tmp_output_bam, tmp_output_bed, tmp_output_bw, tmp_output_bam )
        # os.system(tmpinst)
        # time.sleep(stringtie_interval)
        command_list.append(tmpinst)

def TSS(directory, ref, up_stream, down_stream):
    file_list = []
    bw_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith("bw"):
            bw_list.append(filedir.rstrip(".bw"))

    for bw in bw_list:
        tmp_imput = bw + ".bw"
        tmp_output_gz = bw + "._3kb_TSS.gz"
        tmp_output_png = bw + "._TSS_3kb_Profile.png"

        tmpinst = "computeMatrix reference-point --referencePoint TSS -p 30 -b %s -a %s -R %s -S %s --skipZeros -o %s &&" \
                  "plotProfile -m %s -out %s &" % (
                      up_stream, down_stream, ref, tmp_imput, tmp_output_gz,tmp_output_gz,tmp_output_png)
        # os.system(tmpinst)
        # time.sleep(stringtie_interval)
        command_list.append(tmpinst)

def Peakcalling(directory):
    file_list = []
    bed_list = []
    file_line = 'cd %s'%(directory)
    command_list.append(file_line)
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(dir))
        for file in files:
            file_list.append(str(file))

    for filedir in file_list:
        if filedir.endswith("bed"):
            bed_list.append(filedir.rstrip(".bed"))

    for bed in bed_list:
        tmp_imput = bed + ".bed"
        tmp_log = bed + ".log"
        tmpinst = "macs2 callpeak -t %s -g 2.1e9 -q 0.01 --nomodel --shift -100 --extsize 200 -n %s --outdir peak/ >%s 2>&1" % (
                      tmp_imput, bed, tmp_log)
        # os.system(tmpinst)
        # time.sleep(stringtie_interval)
        command_list.append(tmpinst)

def Tn5filtering(directory, genome_fai, fdr):
    file_line = 'cd %s' % (directory)

    command_list.append(file_line)
    command_line = 'bash filter_ocr_by_tn5_density.sh %s %s' % (genome_fai,fdr)
    command_list.append(command_line)



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

    elif sys.argv[2] == "rmdup":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        rmdup(dir_main)

    elif sys.argv[2] == "bowtie":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]

        ref_genome_zm_V4_path = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/B73V4"
        ref_genome_sor_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sor_ref/sbi"

        if sys.argv[4] == "zmv4":
            ref_genome_path = ref_genome_zm_V4_path
        elif sys.argv[4] == "sor":
            ref_genome_path = ref_genome_sor_path
        else:
            ref_genome_path = sys.argv[4]

        print("Using Ref %s" % ref_genome_path)

        bowtie_sort_zip(dir_main, ref_genome_path)
        bowtie_sort_unzip(dir_main, ref_genome_path)

    elif sys.argv[2] == "TSS":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        up_stream = sys.argv[5]
        down_stream = sys.argv[6]

        ref_genome_zm_V4_path = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/Zmv4_refseq_genes_TSS.txt"
        ref_genome_sor_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sor_ref/Sbi_refseq_genes_TSS.txt"


        if sys.argv[4] == "zmv4":
            ref_genome_path = ref_genome_zm_V4_path
        elif sys.argv[4] == "sor":
            ref_genome_path = ref_genome_sor_path
        else:
            ref_genome_path = sys.argv[4]


        print("Using Ref %s" % ref_genome_path)

        TSS(dir_main, ref_genome_path, up_stream, down_stream)


    elif sys.argv[2] == "peak":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        Peakcalling(dir_main)

    elif sys.argv[2] == "tn5":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]

        ref_genome_zm_V4_path = "/data5/zhuyuzhi/Reference_Genomes/Mays_New/B73V4/Zm-B73-REFERENCE-GRAMENE-4.0.fa.fai"
        ref_genome_sor_path = "/data5/RNA_Seq_Database/wzj_reference_genome/sor_ref/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.fai"

        if sys.argv[4] == "zmv4":
            ref_genome_path = ref_genome_zm_V4_path
        elif sys.argv[4] == "sor":
            ref_genome_path = ref_genome_sor_path
        else:
            ref_genome_path = sys.argv[4]
        fdr = sys.argv[5]

        Tn5filtering(dir_main,ref_genome_path,fdr)


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
