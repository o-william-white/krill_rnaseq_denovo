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