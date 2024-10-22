import pandas as pd
import numpy as np
# Read the text files into pandas DataFrames
file1 = '/home/analysis/result/snp_annovar/snp_annovar.txt'
file2 = '/home/analysis/result/snp_VEP/snp_VEP.txt'
file3 = '/home/analysis/result/indel_annovar/indel_annovar.txt'
file4 = '/home/analysis/result/indel_VEP/indel_VEP.txt'

df1 = pd.read_csv(file1, sep='\t', low_memory=False)  # Assuming tab-separated values, adjust sep if necessary
df2 = pd.read_csv(file2, sep='\t', low_memory=False)

# Identify the first three columns to use for the merge
merge_columns = df1.columns[:3].tolist()

# Merge the DataFrames on the first three columns, keeping columns from df1
VEP_ANNOVAR = pd.merge(df1, df2, on=merge_columns, suffixes=('_ANN', '_VEP'))

# Drop the duplicate columns from df2
for col in merge_columns:
    if col + '_VEP' in VEP_ANNOVAR.columns:
        VEP_ANNOVAR.drop(col + '_VEP', axis=1, inplace=True)

# Columns to keep based on criteria
column_indices_to_keep = [1, 2, 3, 4, 5, 6, 7, 9, 14, 15, 16, 17, 18, 32, 37, 39, 42, 45, 48, 54, 57, 60, 80, 100, 108, 111, 130, 131, 148, 156]
columns_to_keep = [VEP_ANNOVAR.columns[i-1] for i in column_indices_to_keep if i-1 < len(VEP_ANNOVAR.columns)]

# Columns from 195 to column that ends with 'Hom count_ANN'
hom_count_ANN_column = [col for col in VEP_ANNOVAR.columns if col.endswith('Hom count_ANN')]
if hom_count_ANN_column:
    hom_count_ANN_index = VEP_ANNOVAR.columns.get_loc(hom_count_ANN_column[0])
    if len(VEP_ANNOVAR.columns) > 194:
        columns_to_keep.extend(VEP_ANNOVAR.columns[194:hom_count_ANN_index + 1])

# 4- These columns (column names end with ensGene_VEP, EXON, Feature, Amino_acids, am_class, 1000Gp3_AF, ExAC_AF, exomes_AF, genomes_AF, rs_dbSNP, clinvar_clnsig, clinvar_CLNDN, HGMD, HGMD_CLASS, HGMD_PHEN, tapes_VEP)
additional_columns = ['ensGene_VEP', 'EXON', 'Feature', 'HGVSp', 'Amino_acids', 'am_class', '1000Gp3_AF', 'ExAC_AF', 'exomes_AF', 'genomes_AF', 'rs_dbSNP', 'clinvar', 'clinvar_CLNDN', 'HGMD', 'HGMD_CLASS', 'HGMD_PHEN', 'tapes_ANN', 'tapes_VEP']
columns_to_keep.extend([col for col in VEP_ANNOVAR.columns if any(col.endswith(suffix) for suffix in additional_columns)])

# Filter the DataFrame to keep only the selected columns
VEP_ANNOVAR_filtered = VEP_ANNOVAR[columns_to_keep]

# Fix the prediction categories, merge benign and benign auto #
def replace_categories(column):
        return column.apply(lambda x: 'Benign' if 'Benign' in x and 'Likely' not in x else x)

VEP_ANNOVAR_filtered = VEP_ANNOVAR_filtered.copy()

VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_ANN'] = replace_categories(VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_ANN'])
VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_VEP'] = replace_categories(VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_VEP'])

# Find the column ending with "WT count_ANN"
ending_column = [col for col in VEP_ANNOVAR_filtered.columns if col.endswith("WT count_ANN")][0]
ending_index = VEP_ANNOVAR_filtered.columns.get_loc(ending_column)

# Select columns from the 30th to the one immediately before the ending column
columns_to_check = VEP_ANNOVAR_filtered.columns[30:ending_index]

VEP_ANNOVAR_filtered.loc[:, columns_to_check] = VEP_ANNOVAR_filtered.loc[:, columns_to_check].replace('./.', '0/0')
# Calculate Het.count2 and Hom.count2
VEP_ANNOVAR_filtered['Het.count2'] = VEP_ANNOVAR_filtered[columns_to_check].apply(
            lambda row: sum((row != "0/0") & (row != "0|0") & (row != "1/1") & (row != "1|1") & (row != "2/2") & (row != "2|2")), axis=1
            )
VEP_ANNOVAR_filtered['Hom.count2'] = VEP_ANNOVAR_filtered[columns_to_check].apply(
            lambda row: sum((row == "1/1") & (row == "1|1") & (row == "2/2") & (row == "2|2")), axis=1
            )

# Calculate IDs
VEP_ANNOVAR_filtered['IDs'] = VEP_ANNOVAR_filtered.apply(
            lambda row: row['avsnp150'] if row['avsnp150'] != '.' else row['rs_dbSNP'], axis=1
            )

# Calculate AF
VEP_ANNOVAR_filtered['AF'] = VEP_ANNOVAR_filtered.apply(
            lambda row: row['gnomAD_genome_ALL'] if row['gnomAD_genome_ALL'] != '.' else (
                        row['1000Gp3_AF'] if row['1000Gp3_AF'] != '.' else row['ExAC_ALL']
                            ), axis=1
            )

# Calculate ClinvarID
VEP_ANNOVAR_filtered['ClinvarID'] = np.where(
            VEP_ANNOVAR_filtered['CLNSIG'] == '.', VEP_ANNOVAR_filtered['clinvar'], VEP_ANNOVAR_filtered['CLNSIG']
            )

# Calculate Function
VEP_ANNOVAR_filtered['Function'] = np.where(
            VEP_ANNOVAR_filtered['ExonicFunc.ensGene_VEP'] == '.', VEP_ANNOVAR_filtered['ExonicFunc.ensGene_ANN'], VEP_ANNOVAR_filtered['ExonicFunc.ensGene_VEP']
            )
VEP_ANNOVAR_filtered['Function2'] = np.where(
            VEP_ANNOVAR_filtered['ExonicFunc.ensGene_ANN'] == '.', VEP_ANNOVAR_filtered['ExonicFunc.ensGene_VEP'], VEP_ANNOVAR_filtered['ExonicFunc.ensGene_ANN']
            )


# Define the list of clinical significance categories for concordance
clin = ['Pathogenic', 'Likely Pathogenic']

# Filter the DataFrame for concordant entries
Concordant_list = VEP_ANNOVAR_filtered[
        (VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_ANN'].isin(clin)) &
        (VEP_ANNOVAR_filtered['Prediction_ACMG_tapes_VEP'].isin(clin))

        ]


# Save the merged DataFrame to a new text file
VEP_ANNOVAR_filtered.to_csv('/home/analysis/result/SNP/snps_all_merged_file.csv', index=False)
Concordant_list.to_csv('/home/analysis/result/SNP/Pathogenic_snps_ACMG_list.csv', index=False)


# Indels
df3 = pd.read_csv(file3, sep='\t', low_memory=False)
df4 = pd.read_csv(file4, sep='\t', low_memory=False)

merge_columns_indels = df3.columns[:3].tolist()


VEP_ANNOVAR_indels = pd.merge(df3, df4, on=merge_columns_indels, suffixes=('_ANN', '_VEP'))


for col in merge_columns_indels:
    if col + '_VEP' in VEP_ANNOVAR_indels.columns:
        VEP_ANNOVAR_indels.drop(col + '_VEP', axis=1, inplace=True)


column_indices_to_keep = [1, 2, 3, 4, 5, 6, 7, 9, 14, 15, 16, 17, 18, 32, 37, 39, 42, 45, 48, 54, 57, 60, 80, 100, 108, 111, 130, 131, 148, 156]
columns_to_keep_indels = [VEP_ANNOVAR_indels.columns[i-1] for i in column_indices_to_keep if i-1 < len(VEP_ANNOVAR_indels.columns)]


hom_count_ANN_column_indels = [col for col in VEP_ANNOVAR_indels.columns if col.endswith('Hom count_ANN')]
if hom_count_ANN_column_indels:
    hom_count_ANN_index_indels = VEP_ANNOVAR_indels.columns.get_loc(hom_count_ANN_column_indels[0])
    if len(VEP_ANNOVAR_indels.columns) > 194:
        columns_to_keep_indels.extend(VEP_ANNOVAR_indels.columns[194:hom_count_ANN_index_indels + 1])


additional_columns_indels = ['ensGene_VEP', 'EXON', 'Feature', 'HGVSp', 'Amino_acids', 'am_class', '1000Gp3_AF', 'ExAC_AF', 'exomes_AF', 'genomes_AF', 'rs_dbSNP', 'clinvar', 'clinvar_CLNDN', 'HGMD', 'HGMD_CLASS', 'HGMD_PHEN', 'tapes_ANN', 'tapes_VEP']
columns_to_keep_indels.extend([col for col in VEP_ANNOVAR_indels.columns if any(col.endswith(suffix) for suffix in additional_columns_indels)])

VEP_ANNOVAR_filtered_indels = VEP_ANNOVAR_indels[columns_to_keep_indels]


def replace_categories(column):
        return column.apply(lambda x: 'Benign' if 'Benign' in x and 'Likely' not in x else x)

VEP_ANNOVAR_filtered_indels = VEP_ANNOVAR_filtered_indels.copy()

VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_ANN'] = replace_categories(VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_ANN'])
VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_VEP'] = replace_categories(VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_VEP'])


ending_column_indels = [col for col in VEP_ANNOVAR_filtered_indels.columns if col.endswith("WT count_ANN")][0]
ending_index_indels = VEP_ANNOVAR_filtered_indels.columns.get_loc(ending_column_indels)


columns_to_check_indels = VEP_ANNOVAR_filtered_indels.columns[30:ending_index_indels]
VEP_ANNOVAR_filtered_indels.loc[:, columns_to_check_indels] = VEP_ANNOVAR_filtered_indels.loc[:, columns_to_check_indels].replace('./.', '0/0')

VEP_ANNOVAR_filtered_indels['Het.count2'] = VEP_ANNOVAR_filtered_indels[columns_to_check_indels].apply(
            lambda row: sum((row != "0/0") & (row != "0|0") & (row != "1/1") & (row != "1|1") & (row != "2/2") & (row != "2|2")), axis=1
            )
VEP_ANNOVAR_filtered_indels['Hom.count2'] = VEP_ANNOVAR_filtered_indels[columns_to_check_indels].apply(
            lambda row: sum((row == "1/1") & (row == "1|1") & (row == "2/2") & (row == "2|2")), axis=1
            )
VEP_ANNOVAR_filtered_indels['IDs'] = VEP_ANNOVAR_filtered_indels.apply(
            lambda row: row['avsnp150'] if row['avsnp150'] != '.' else row['rs_dbSNP'], axis=1
            )
VEP_ANNOVAR_filtered_indels['AF'] = VEP_ANNOVAR_filtered_indels.apply(
            lambda row: row['gnomAD_genome_ALL'] if row['gnomAD_genome_ALL'] != '.' else (
                        row['1000Gp3_AF'] if row['1000Gp3_AF'] != '.' else row['ExAC_ALL']
                            ), axis=1
            )
VEP_ANNOVAR_filtered_indels['ClinvarID'] = np.where(
            VEP_ANNOVAR_filtered_indels['CLNSIG'] == '.', VEP_ANNOVAR_filtered_indels['clinvar'], VEP_ANNOVAR_filtered_indels['CLNSIG']
            )
VEP_ANNOVAR_filtered_indels['Function'] = np.where(
            VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_VEP'] == '.', VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_ANN'], VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_VEP']
            )
VEP_ANNOVAR_filtered_indels['Function2'] = np.where(
            VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_ANN'] == '.', VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_VEP'], VEP_ANNOVAR_filtered_indels['ExonicFunc.ensGene_ANN']
            )

clin = ['Pathogenic', 'Likely Pathogenic']

Concordant_list_indels = VEP_ANNOVAR_filtered_indels[
        (VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_ANN'].isin(clin)) &
        (VEP_ANNOVAR_filtered_indels['Prediction_ACMG_tapes_VEP'].isin(clin))

        ]

VEP_ANNOVAR_filtered_indels.to_csv('/home/analysis/result/INDEL/indels_all_merged_file.csv', index=False)
Concordant_list_indels.to_csv('/home/analysis/result/INDEL/Pathogenic_indels_ACMG_list.csv', index=False)
