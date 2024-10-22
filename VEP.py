import os
import subprocess
import glob

# Define the directories and paths
input_dir1 = "/home/analysis/result"
output_dir1 = "/home/analysis/result"

# VEP executable path
vep_path = "/home/n9261834/ensembl-vep/./vep"

# Plugin paths
clinvar_path = "/home/n9261834/.vep/Plugins/clinvar.vcf.gz"
hgmd_path = "/home/n9261834/.vep/Plugins/HGMD_Pro_2024.1_hg38.vcf.gz"
revel_path = "/home/n9261834/.vep/Plugins/new_tabbed_revel_grch38.tsv.gz"
dbscsnv_path = "/home/n9261834/.vep/Plugins/dbscSNV1.1_GRCh38.txt.gz"
gnomad_path = "/home/n9261834/.vep/Plugins/gnomad.ch.genomesv3.tabbed.tsv.gz"
alpha_missense_path = "/home/n9261834/.vep/Plugins/AlphaMissense_hg38.tsv.gz"
dbnsfp_path = "/home/n9261834/.vep/Plugins/dbNSFP4.3a_grch38.gz"

# Iterate over input files
for input1 in glob.glob(os.path.join(input_dir1, "*final.vcf.gz")):
    input_base1 = os.path.basename(input1).replace(".vcf.gz", "")
    output1 = os.path.join(output_dir1, f"{input_base1}.vep.annotated.vcf")

    # Define the VEP command
    command = [
        vep_path,
        "-i", input1,
        "-o", output1,
        "--vcf",
        "--hgvs",
        "--cache",
        "--offline",
        "--assembly", "GRCh38",
        "--custom", f"file={clinvar_path},short_name=clinvar,format=vcf,type=exact,coords=0,fields=ALLELEID%CLNDN%CLNDISDB%CLNREVSTAT%CLNSIG",
        "--custom", f"file={hgmd_path},short_name=HGMD,format=vcf,type=exact,coords=0,fields=%CLASS%STRAND%DB%PHEN%RANKSCORE",
        "--plugin", f"REVEL,{revel_path}",
        "--plugin", f"dbscSNV,{dbscsnv_path}",
        "--plugin", f"gnomADc,{gnomad_path}",
        "--plugin", f"AlphaMissense,file={alpha_missense_path}",
        "--plugin", f"dbNSFP,{dbnsfp_path},rs_dbSNP,SIFT_score,SIFT4G_score,LRT_score,MutationTaster_score,MutationAssessor_score,FATHMM_score,fathmm-MKL_coding_score,PROVEAN_score,MetaSVM_score,MetaLR_score,M-CAP_score,GenoCanyon_score,GERP++_NR,GERP++_RS,Polyphen2_HDIV_score,Polyphen2_HVAR_score,DANN_score,Interpro_domain,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ExAC_AC,ExAC_AF,gnomAD_exomes_flag,gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_genomes_flag,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF"
        ]

    # Execute the command
    subprocess.run(command, check=True)
