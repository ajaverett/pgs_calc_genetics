import pandas as pd
import numpy as np
import csv
import re
from snps import SNPs
from utilities import find_column_name, to_snake_case, rsid_find, is_ambiguous, count_effect_alleles

def calculate_polygenic_score(
        genotype_file, 
        pgs_weights_file):
    
    df_genotype = SNPs(genotype_file)\
        .snps\
        .reset_index()\
        .assign(
            chr=lambda df: df[find_column_name(df, 'chr')].astype(str),
            pos=lambda df: df[find_column_name(df, 'pos')].astype(int))

    df_pgs = pd\
        .read_csv(pgs_weights_file, sep='\t', skiprows=rsid_find(pgs_weights_file))\
        .rename(columns=to_snake_case)\
        .assign(
            chr=lambda df: df[find_column_name(df, 'chr')].astype(str),
            pos=lambda df: df[find_column_name(df, 'pos')].astype(int),
            eff=lambda df: df[find_column_name(df, 'eff')].astype(str),
            oth=lambda df: df[find_column_name(df, 'oth')].astype(str),
            weight=lambda df: df[find_column_name(df, ['weight', 'beta'])].astype(float))\
        .merge(df_genotype, on=['chr', 'pos'])\
        .assign(ambiguous=lambda df: df.apply(
            lambda x: is_ambiguous(x["eff"], x["oth"]), axis=1))\
        .query('ambiguous == False')\
        .assign(effect_allele_count=lambda x: x.apply(count_effect_alleles, axis=1))\
        .dropna(subset=['effect_allele_count'])\
        .assign(effect_allele_count=lambda x: x['effect_allele_count'].astype(int))\
        .assign(weighted_effect=lambda x: x['effect_allele_count'] * x['weight'])

    PGS_total = df_pgs['weighted_effect'].sum()

    return PGS_total

# usage
calculate_polygenic_score(
    genotype_file=r"my_genotype_file.txt",
    pgs_weights_file=r"my_pgs_weights_file.txt"
)