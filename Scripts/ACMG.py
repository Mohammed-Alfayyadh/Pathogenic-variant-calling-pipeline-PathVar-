import os
import subprocess
import glob

# Define the directories and paths
# Define the directories and paths
input_dir = "/home/analysis/result"
output_dir = "/home/analysis/result"

output_dir_snp_vep = os.path.join(output_dir, "snp_VEP/")
output_dir_indel_vep = os.path.join(output_dir, "indel_VEP/")
output_dir_snp_annovar = os.path.join(output_dir, "snp_annovar/")
output_dir_indel_annovar = os.path.join(output_dir, "indel_annovar/")

# Create output directories
os.makedirs(output_dir_snp_vep, exist_ok=True)
os.makedirs(output_dir_indel_vep, exist_ok=True)
os.makedirs(output_dir_snp_annovar, exist_ok=True)
os.makedirs(output_dir_indel_annovar, exist_ok=True)

tapes_path = "/home/n9261834/tapes/tapes.py"

# Iterate over input files
for input_file in glob.glob(os.path.join(input_dir, "*final.vcf.gz")):
    input_base = os.path.basename(input_file).replace(".vcf.gz", "")

    output_annovar_path = os.path.join(output_dir, f"{input_base}.annovar.annotated.vcf")

    command = [
        "python",
        tapes_path,
        "annotate",
        "-i", input_file,
        "-o", output_annovar_path,
        "-a", "hg38",
        "--ref_anno", "ensGene"

        ]

    # Execute the command
    subprocess.run(command, check=True)
# Step 2: Sort the annotated VCF files
def sort_vcf_files(annotated, input_dir, output_dir):
    for input_file in glob.glob(os.path.join(input_dir, annotated)):

        subprocess.run([
            "python",
            tapes_path,
            "sort",
            "-i", input_file,
            "-o", output_dir,
            "--tab",
            "-a", "hg38"
        ], check=True)

# Sort files for Annovar and VEP
sort_vcf_files("*snp.final.annovar.annotated.hg38_multianno.vcf", output_dir, output_dir_snp_annovar)
sort_vcf_files("*indel.final.annovar.annotated.hg38_multianno.vcf", output_dir, output_dir_indel_annovar)
sort_vcf_files("*snp.final.vep.annotated.vcf", output_dir, output_dir_snp_vep)
sort_vcf_files("*indel.final.vep.annotated.vcf", output_dir, output_dir_indel_vep)
