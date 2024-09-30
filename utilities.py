import pandas as pd
import numpy as np
import csv
import re
from snps import SNPs

def find_column_name(df, substrings):
    if isinstance(substrings, str):
        substrings = [substrings]
    
    # Iterate over columns and substrings
    for col in df.columns:
        for substring in substrings:
            if substring in col:
                return col
    return None

def to_snake_case(s):
    if s.isupper():                     
        return s.lower()
    elif re.match(r'.+[A-Z]{2,}$', s):  
        return s.lower()
    else:
        s = re.sub(r'(?<=[a-z0-9])(?=[A-Z])', '_', s).lower()
        return s

def rsid_find(file_path):
    with open(file_path, mode='r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for line_num, row in enumerate(reader, start=1):
            if any("rsid" in cell.lower() for cell in row):
                return line_num - 1

    return None 

def is_ambiguous(a1, a2):
    return {a1, a2} == {'A', 'T'} or \
            {a1, a2} == {'C', 'G'}

# Function to count effect alleles
def count_effect_alleles(row):
    genotype = row['genotype']
    effect_allele = row["eff"]
    other_allele = row["oth"]

    # Exclude genotypes with missing data or unusual length
    if pd.isna(genotype) or len(genotype) != 2:
        return np.nan

    genotype = genotype.upper()
    effect_allele = effect_allele.upper()
    other_allele = other_allele.upper()

    alleles = list(genotype)

    # Check if alleles match the effect and other alleles
    if set(alleles).issubset({effect_allele, other_allele}):
        return alleles.count(effect_allele)
    else:
        # Attempt strand flip
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        flipped_alleles = [complement.get(a, 'N') for a in alleles]
        flipped_effect_allele = complement.get(effect_allele, 'N')
        flipped_other_allele = complement.get(other_allele, 'N')

        # Check for invalid nucleotides
        if 'N' in flipped_alleles or 'N' in [flipped_effect_allele, flipped_other_allele]:
            return np.nan

        if set(flipped_alleles).issubset({flipped_effect_allele, flipped_other_allele}):
            return flipped_alleles.count(flipped_effect_allele)
        else:
            return np.nan

