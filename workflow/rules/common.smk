import pandas as pd
import sys
import os

# set configfile
configfile: "config/config.yaml"

# configfile parameters
reference = config["reference"]

# read sample data
if os.path.exists(config["samples"]):
    sample_data = pd.read_csv(config["samples"]).set_index("id", drop=False)
else:
    sys.exit(f"Error: samples.csv file '{config['samples']}' does not exist")

def get_fastq(wildcards):
    fwd = sample_data.loc[wildcards.sample, "forward"]
    rev = sample_data.loc[wildcards.sample, "reverse"]
    return [fwd, rev]

# dictionary of groups to sample ids
dict_groups = sample_data.groupby('group')['id'].apply(list).to_dict()

def get_technical_replicates_forward(wildcards):
    return [ f'results/fastp/{x}_R1.fastq' for x in dict_groups[wildcards.group]]

def get_technical_replicates_reverse(wildcards):
    return [ f'results/fastp/{x}_R2.fastq' for x in dict_groups[wildcards.group]]