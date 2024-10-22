import subprocess
import os

# Define the input and output directories
input_dir1 = "/home/analysis"
output_dir1 = "/home/analysis/result"

# Define the reference genome and known sites
ref_genome = "/home/n9261834/phd/data/gatk/Homo_sapiens_assembly38.fasta"
known_sites = "/home/n9261834/phd/data/gatk/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

# Define the resources for VariantRecalibrator
hapmap = os.path.join(input_dir1, "hapmap_3.3.hg38.vcf.gz")
omni = os.path.join(input_dir1, "1000G_omni2.5.hg38.vcf.gz")
thousand_genomes = os.path.join(input_dir1, "1000G_phase1.snps.high_confidence.hg38.vcf.gz")
dbsnp = os.path.join(input_dir1, "Homo_sapiens_assembly38.dbsnp138.vcf.gz")
mills = os.path.join(input_dir1, "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
axiomPoly = os.path.join(input_dir1, "Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz")
blacklist_bed = os.path.join(input_dir1, "blacklist.bed")

def run_command(command):
   process = subprocess.run(command, stderr=subprocess.PIPE)
   if process.returncode != 0:
       print(f"Error running {' '.join(command)}: {process.stderr.decode()}")
       return process.returncode

# Loop through the input files and run bwa mem and samtools view
for input1 in os.listdir(input_dir1):
    if input1.endswith(".fq"):
        input_base1 = os.path.basename(input1).replace(".fq", "")
        output1 = os.path.join(output_dir1, f"{input_base1}.hg38.bam")

        bwa_command = [
            "bwa", "mem", "-M", "-t", "2", "-R", f"@RG\\tID:{input_base1}\\tSM:{input_base1}\\tLB:library\\tPL:IONTORRENT",
            ref_genome, os.path.join(input_dir1, input1)
        ]
        samtools_command = [
            "samtools", "view", "-b", "-h", "-o", output1
        ]

        with subprocess.Popen(bwa_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as bwa_process:
            with subprocess.Popen(samtools_command, stdin=bwa_process.stdout, stderr=subprocess.PIPE) as samtools_process:
                bwa_process.stdout.close()  # Allow bwa_process to receive a SIGPIPE if samtools_process exits.
                bwa_stderr = bwa_process.stderr.read()
                samtools_stderr = samtools_process.stderr.read()

        if bwa_process.returncode != 0:
            print(f"Error running bwa mem: {bwa_stderr.decode()}")
            continue

        if samtools_process.returncode != 0:
            print(f"Error running samtools view: {samtools_stderr.decode()}")
            continue

# Loop through the input files again and run picard SortSam
for input2 in os.listdir(output_dir1):
    if input2.endswith(".hg38.bam"):
        input_base2 = os.path.basename(input2).replace(".bam", "")
        output2 = os.path.join(output_dir1, f"{input_base2}.sorted.bam")

        picard_command = [
            "picard", "-Xmx7g", "SortSam", "I=" + os.path.join(output_dir1, input2),
            "O=" + output2, "VALIDATION_STRINGENCY=LENIENT", "SORT_ORDER=coordinate",
            "MAX_RECORDS_IN_RAM=3000000", "CREATE_INDEX=True"
        ]

        run_command(picard_command)

# Loop through the input files again and run picard MarkDuplicates
for input3 in os.listdir(output_dir1):
    if input3.endswith(".sorted.bam"):
        input_base3 = os.path.basename(input3).replace(".bam", "")
        output3 = os.path.join(output_dir1, f"{input_base3}.dup.bam")
        metrics_file = f"/home/n9261834/phd/data/data-hg38/joseph/{input_base3}.txt"

        picard_command = [
            "picard", "-Xmx7g", "MarkDuplicates", "I=" + os.path.join(output_dir1, input3),
            "O=" + output3, "METRICS_FILE=" + metrics_file
        ]

        run_command(picard_command)

# Index the duplicate-marked BAM files
for file in os.listdir(output_dir1):
    if file.endswith(".dup.bam"):
        index_command = ["samtools", "index", os.path.join(output_dir1, file)]
        run_command(index_command)

# Run gatk BaseRecalibrator on the duplicate-marked BAM files
for input4 in os.listdir(output_dir1):
    if input4.endswith(".dup.bam"):
        input_base4 = os.path.basename(input4).replace(".bam", "")
        output4 = os.path.join(output_dir1, f"{input_base4}.recal_data.table")

        gatk_command = [
            "gatk", "--java-options", "-Xmx7g", "BaseRecalibrator",
            "-I", os.path.join(output_dir1, input4),
            "-R", ref_genome,
            "--known-sites", known_sites,
            "-O", output4
        ]

        run_command(gatk_command)

# Apply BQSR to the duplicate-marked BAM files
for input5 in os.listdir(output_dir1):
    if input5.endswith(".dup.bam"):
        input_base5 = os.path.basename(input5).replace(".bam", "")
        recal_table = os.path.join(output_dir1, f"{input_base5}.recal_data.table")
        output5 = os.path.join(output_dir1, f"{input_base5}.bqsr.bam")

        gatk_command = [
            "gatk", "--java-options", "-Xmx7g", "ApplyBQSR",
            "-I", os.path.join(output_dir1, input5),
            "-R", ref_genome,
            "--bqsr-recal-file", recal_table,
            "-O", output5
        ]

        run_command(gatk_command)

# Index the recalibrated BAM files
for file in os.listdir(output_dir1):
    if file.endswith(".bqsr.bam"):
        index_command = ["samtools", "index", os.path.join(output_dir1, file)]
        run_command(index_command)

# Run gatk HaplotypeCaller on the recalibrated BAM files
for input6 in os.listdir(output_dir1):
    if input6.endswith(".bqsr.bam"):
        input_base6 = os.path.basename(input6).replace(".bam", "")
        output6 = os.path.join(output_dir1, f"{input_base6}.g.vcf.gz")

        gatk_command = [
                "gatk", "--java-options", "-Xmx7g", "HaplotypeCaller",
                "-I", os.path.join(output_dir1, input6),
                "-R", ref_genome,
                "-ERC", "GVCF",
                "-O", output6
            ]
        run_command(gatk_command)

# Create a sample map for GenomicsDBImport
gvcfs = [os.path.join(output_dir1, file) for file in os.listdir(output_dir1) if file.endswith(".g.vcf.gz")]
sample_map_path = os.path.join(output_dir1, "sample_map.txt")

with open(sample_map_path, 'w') as sample_map_file:
    for gvcf in gvcfs:
        sample_name = os.path.basename(gvcf).replace(".g.vcf.gz", "")
        sample_map_file.write(f"{sample_name}\t{gvcf}\n")

# Create tmp directory for GenomicsDBImport
tmp_dir = os.path.join(output_dir1, "tmp")
os.makedirs(tmp_dir, exist_ok=True)

# Run GenomicsDBImport with appropriate parameters for large datasets
genomicsdbimport_command = [
        "gatk", "--java-options", "-Xmx7g", "GenomicsDBImport",
        "--genomicsdb-workspace-path", os.path.join(output_dir1, "genomicsdb"),
        "--sample-name-map", sample_map_path,
        "--reader-threads", "5",
        "--tmp-dir", os.path.join(output_dir1, "tmp"),  # Ensure sufficient space in tmp directory
        "--batch-size", "50",  # Adjust batch size as needed
        "--consolidate", "true",
        ]

intervals = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
for interval in intervals:
    genomicsdbimport_command.extend(["--intervals", interval])

run_command(genomicsdbimport_command)

# Genotype the gVCF files from GenomicsDB using GATK GenotypeGVCFs
genotyped_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.vcf.gz")

gatk_genotype_command = [
        "gatk", "--java-options", "-Xmx7g", "GenotypeGVCFs",
        "-R", ref_genome,
        "-V", "gendb://" + os.path.join(output_dir1, "genomicsdb"),
        "-O", genotyped_vcf
        ]

run_command(gatk_genotype_command)

# Annotate the genotyped VCF file
annotated_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.vcf.gz")
bcftools_annotate_command = [
        "bcftools", "annotate", "-c", "ID",
        "-a", known_sites,
        "-o", annotated_vcf,
        "-O", "z",
        genotyped_vcf
        ]
run_command(bcftools_annotate_command)

# Index the annotated VCF file using tabix
tabix_index_command = ["tabix", "-p", "vcf", annotated_vcf]
run_command(tabix_index_command)

# Generate SNP VCF file
snp_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.snp.vcf.gz")
gatk_select_snps_command = [
        "gatk", "--java-options", "-Xmx7g", "SelectVariants",
        "-R", ref_genome,
        "-V", annotated_vcf,
        "--select-type", "SNP",
        "-O", snp_vcf
        ]
run_command(gatk_select_snps_command)

# Generate Indel VCF file
indel_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.indel.vcf.gz")
gatk_select_indels_command = [
        "gatk", "--java-options", "-Xmx7g", "SelectVariants",
        "-R", ref_genome,
        "-V", annotated_vcf,
        "--select-type", "INDEL",
        "-O", indel_vcf
        ]
run_command(gatk_select_indels_command)

# Hard filter SNPs
snp_filtered_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.snp.hardFilter.vcf.gz")

gatk_hard_filter_snp_command = [
        "gatk", "--java-options", "-Xmx7g", "VariantFiltration",
        "-R", ref_genome,
        "-V", snp_vcf,
        "--filter-expression", "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
        "--filter-name", "hardfilter_snp",
        "-O", snp_filtered_vcf
        ]
run_command(gatk_hard_filter_snp_command)

# Hard filter Indels
indel_filtered_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.indel.hardFilter.vcf.gz")

gatk_hard_filter_indel_command = [
        "gatk", "--java-options", "-Xmx7g", "VariantFiltration",
        "-R", ref_genome,
        "-V", indel_vcf,
        "--filter-expression", "QD < 2.0 || QUAL < 30.0 || FS > 200.0 || ReadPosRankSum < -20.0",
        "--filter-name", "hardfilter_indel",
        "-O", indel_filtered_vcf
        ]
run_command(gatk_hard_filter_indel_command)

# Run VariantRecalibrator for SNPs
recal_snp = os.path.join(output_dir1, "snp.recal")
tranches_snp = os.path.join(output_dir1, "snp.tranches")

gatk_variantrecalibrator_snp_command = [
        "gatk", "--java-options", "-Xmx7g", "VariantRecalibrator",
        "-R", ref_genome,
        "-V", snp_filtered_vcf,
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0", hapmap,
        "--resource:omni,known=false,training=true,truth=true,prior=12.0", omni,
        "--resource:1000G,known=false,training=true,truth=false,prior=10.0", thousand_genomes,
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0", dbsnp,
        "-an", "QD", "-an", "MQ", "-an", "MQRankSum", "-an", "ReadPosRankSum",
        "-an", "FS", "-an", "SOR", "-an", "DP",
        "-mode", "SNP",
        "-O", recal_snp,
        "--tranches-file", tranches_snp,
        ]
run_command(gatk_variantrecalibrator_snp_command)

# Run VariantRecalibrator for indels
recal_indel = os.path.join(output_dir1, "indel.recal")
tranches_indel = os.path.join(output_dir1, "indel.tranches")

gatk_variantrecalibrator_indel_command = [
        "gatk", "--java-options", "-Xmx7g", "VariantRecalibrator",
        "-R", ref_genome,
        "-V", indel_filtered_vcf,
        "--resource:mills,known=false,training=true,truth=true,prior=12.0", mills,
        "--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0", axiomPoly,
        "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0", dbsnp,
        "-an", "QD", "-an", "DP", "-an", "FS", "-an", "SOR", "-an", "ReadPosRankSum", "-an", "MQRankSum",
        "-mode", "INDEL",
        "-O", recal_indel,
        "--tranches-file", tranches_indel,
        ]
run_command(gatk_variantrecalibrator_indel_command)

# Apply VQSR to SNPs
snp_vqsr_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.snp.hardFilter.vqsr.vcf.gz")

gatk_apply_vqsr_snp_command = [
        "gatk", "--java-options", "-Xmx7g", "ApplyVQSR",
        "-R", ref_genome,
        "-V", snp_filtered_vcf,
        "--recal-file", recal_snp,
        "--tranches-file", tranches_snp,
        "--mode", "SNP",
        "--truth-sensitivity-filter-level", "99.0",
        "-O", snp_vqsr_vcf
        ]

run_command(gatk_apply_vqsr_snp_command)

# Apply VQSR to indels
indel_vqsr_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.indel.hardFilter.vqsr.vcf.gz")

gatk_apply_vqsr_indel_command = [
        "gatk", "--java-options", "-Xmx7g", "ApplyVQSR",
        "-R", ref_genome,
        "-V", indel_filtered_vcf,
        "--recal-file", recal_indel,
        "--tranches-file", tranches_indel,
        "--mode", "INDEL",
        "--truth-sensitivity-filter-level", "99.0",
        "-O", indel_vqsr_vcf
        ]

run_command(gatk_apply_vqsr_indel_command)

# Extract only PASS variants for SNPs
snp_pass_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.snp.hardFilter.vqsr.pass.vcf.gz")
bcftools_extract_pass_snp_command = [
        "bcftools", "view", "-f", "PASS", "-O", "z", "-o", snp_pass_vcf, snp_vqsr_vcf
        ]

run_command(bcftools_extract_pass_snp_command)

# Extract only PASS variants for Indels
indel_pass_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.indel.hardFilter.vqsr.pass.vcf.gz")

bcftools_extract_pass_indel_command = [
        "bcftools", "view", "-f", "PASS", "-O", "z", "-o", indel_pass_vcf, indel_vqsr_vcf
        ]

run_command(bcftools_extract_pass_indel_command)

# Index the final PASS VCF files
snp_pass_vcf_idx = f"{snp_pass_vcf}.tbi"
indel_pass_vcf_idx = f"{indel_pass_vcf}.tbi"

tabix_index_snp_command = ["tabix", "-p", "vcf", snp_pass_vcf]
tabix_index_indel_command = ["tabix", "-p", "vcf", indel_pass_vcf]

run_command(tabix_index_snp_command)
run_command(tabix_index_indel_command)


# Exclude variants in blacklist regions from the final SNPs and Indels files
snp_no_blacklist_vcf_temp = os.path.join(output_dir1, "snps.no_blacklist.temp.vcf")
indel_no_blacklist_vcf_temp = os.path.join(output_dir1, "indels.no_blacklist.temp.vcf")

bedtools_intersect_snp_command = [
        "bedtools", "intersect", "-v", "-a", snp_pass_vcf, "-b", blacklist_bed, "-f", "0.5", "-header"
        ]

bedtools_intersect_indel_command = [
        "bedtools", "intersect", "-v", "-a", indel_pass_vcf, "-b", blacklist_bed, "-f", "0.5", "-header"
        ]

with open(snp_no_blacklist_vcf_temp, 'w') as temp_snp, open(indel_no_blacklist_vcf_temp, 'w') as temp_indel:
    subprocess.run(bedtools_intersect_snp_command, stdout=temp_snp)
    subprocess.run(bedtools_intersect_indel_command, stdout=temp_indel)

# Compress and index the final VCF files without blacklist regions
snp_no_blacklist_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.snp.hardFilter.vqsr.pass.blacklist.snp.final.vcf.gz")
indel_no_blacklist_vcf = os.path.join(output_dir1, "BWA.HaplotypeCaller.rsID.indel.hardFilter.vqsr.pass.blacklist.indel.final.vcf.gz")

bgzip_command_snp = ["bgzip", "-c", snp_no_blacklist_vcf_temp]
bgzip_command_indel = ["bgzip", "-c", indel_no_blacklist_vcf_temp]

with open(snp_no_blacklist_vcf, 'wb') as f_out_snp, open(indel_no_blacklist_vcf, 'wb') as f_out_indel:
    subprocess.run(bgzip_command_snp, stdout=f_out_snp)
    subprocess.run(bgzip_command_indel, stdout=f_out_indel)

tabix_index_snp_no_blacklist_command = ["tabix", "-p", "vcf", snp_no_blacklist_vcf]
tabix_index_indel_no_blacklist_command = ["tabix", "-p", "vcf", indel_no_blacklist_vcf]

run_command(tabix_index_snp_no_blacklist_command)
run_command(tabix_index_indel_no_blacklist_command)
