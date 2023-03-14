import os, sys, shutil, random, re
import pandas as pd

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

    file_list = []

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".fastq.gz"):
            fastq_file = os.path.split(filedir)[0] + "/" + filedir
            tmpinst = 'gzip -df %s &' % (
                fastq_file)
            command_list.append(tmpinst)
        # print(tmpinst)


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


def merge_double(directory):
    file_list = []
    fastq_list = []


    double_edge_seq_fastq_comp = []


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

    for fastq_file in double_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        input_fastq1 = fastq_file + "_1.fastq"
        input_fastq2 = fastq_file + "_2.fastq"
        output_fastq = fastq_file + ".fastq"
        output_tab = fastq_file + ".tab"
        # pandaseq -f forward.fastq -r reverse.fastq -k 5 -o 5 -F > out.fastq 2> out.tab
        tmpinst = 'pandaseq -f %s -r %s -k 5 -o 5 -F > %s 2> %s && rm %s %s &' % (
            input_fastq1, input_fastq2, output_fastq, output_tab, input_fastq1, input_fastq2)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def premiRNA(directory,ref_genome_path):
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

    for fastq_file in single_edge_seq_fastq_comp:  # Process Single-Edge Seq Files
        reads_collapsed = fastq_file.rstrip(".fastq") + "_collapsed.fa"
        reads_vs_refdb = fastq_file.rstrip(".fastq") + "_vs_refdb.arf"
# # mapper.pl example_small_rna_file.fastq -e -h -i -j -k TGGAATTC -l 18 -m -p refdb.fa -s reads_collapsed.fa -t reads_vs_refdb.arf -v -o 4
        tmpinst = 'mapper.pl %s -e -h -i -j -k TGGAATTC -l 18 -m -p %s -s %s -t %s -v -o 4 &' % (
            fastq_file, ref_genome_path, reads_collapsed, reads_vs_refdb)
        # os.system(tmpinst)
        # time.sleep(hisort_interval)
        command_list.append(tmpinst)
        # print(tmpinst)

def predictmiRNA(directory,species):
    file_list = []
    collapsed_fa_list = []

    mature_zma = "./ref/mature_zma_ref.fa"
    mature_sbi = "./ref/mature_sbi_ref.fa"
    hairpin_zma = "./ref/hairpin_zma_ref.fa"
    hairpin_sbi = "./ref/hairpin_sbi_ref.fa"

    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith("_collapsed.fa"):
            collapsed_fa_list.append(filedir)

# quantifier.pl -p ./ref/hairpin_sbi_ref.fa -m ./ref/mature_sbi_ref.fa -r Sorghum_filtered/SRR7628851_collapsed.fa -t sbi -y 16_19
    for fa_file in collapsed_fa_list:
        reads_collapsed = fa_file
        reads_name = fa_file.split('/')[-1].rstrip('_collapsed.fa')
        reads_vs_refdb = fa_file.rstrip("_collapsed.fa") + "_vs_refdb.arf"
        if species == "zma":
            reference_genome = './bowtie/ZmV4/GCF_000005005.2.fa'
            tmpinst = 'quantifier.pl -p %s -m %s -r %s  -y %s &' % (
                hairpin_zma, mature_zma, reads_collapsed, reads_name)
            # os.system(tmpinst)
            # time.sleep(hisort_interval)
            command_list.append(tmpinst)
            # print(tmpinst)
        elif species == "sbi":
            reference_genome = './bowtie/sorghum/GCF_000003195.3.fa'
            tmpinst = 'quantifier.pl -p %s -m %s -r %s  -y %s &' % (
                hairpin_sbi, mature_sbi, reads_collapsed, reads_name)
            # os.system(tmpinst)
            # time.sleep(hisort_interval)
            command_list.append(tmpinst)
            # print(tmpinst)

def miRNA_expression(directory):
    file_list = []
    csv_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))

    for filedir in file_list:
        if filedir.endswith(".csv"):
            csv_list.append(filedir)

    file_merge = pd.DataFrame()
    for file in csv_list:
        file1 = pd.read_csv(file, index_col=0, sep='\t')
        file1_read = file1.loc[:, ['read_count', 'seq(norm)']]
        file1_sequence = file.split('_')[-1].rstrip('.csv')
        file1_header = [file1_sequence + '_counts', file1_sequence + '_normalized seq']
        file1_read.columns = file1_header
        file_merge = pd.concat([file_merge, file1_read], axis=1)

    file_merge.to_csv('Expression_miRNA.csv')




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

    elif sys.argv[2] == "premiRNA":

        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]

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


        print("Using Ref %s" % ref_genome_path)

        premiRNA(dir_main, ref_genome_path)

    elif sys.argv[2] == "qualify":

        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        species = sys.argv[4]

        predictmiRNA(dir_main, species)

    elif sys.argv[2] == "expression":
        dir_main = sys.argv[1]
        miRNA_expression(dir_main)

    elif sys.argv[2] == "merge":
        dir_main = sys.argv[1]
        parallel_task_number = sys.argv[3]
        merge_double(dir_main)


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
