import pandas as pd
import matplotlib.pyplot as plt
import argparse

# argparse
parser = argparse.ArgumentParser(description="Plot sequence length and GC content distributions.")
parser.add_argument("-i", "--input", required=True, help="Input TSV file with sequence data.")
parser.add_argument("-o", "--output", required=True, help="Output image file for the plots.")
args = parser.parse_args()

# read input data
dat = pd.read_csv(args.input, sep="\t", header=None, names=['name', 'length', 'gc_content'])

fig, axs = plt.subplots(1, 3, figsize=(12, 4))

axs[0].hist(dat['length'], bins=50, color="#E69F00")
axs[0].set_title('Sequence Length Distribution')
axs[0].set_xlabel('Length')
axs[0].set_ylabel('Frequency')

axs[1].hist(dat['gc_content'], bins=50, color="#56B4E9")
axs[1].set_title('GC Content Distribution')
axs[1].set_xlabel('GC Content (%)')
axs[1].set_ylabel('Frequency')  

axs[2].scatter(dat['gc_content'], dat['length'], alpha=0.2, color="#009E73")
axs[2].set_title('Length vs GC Content')
axs[2].set_xlabel('GC Content (%)')
axs[2].set_ylabel('Length') 

plt.tight_layout()

# save the figure
plt.savefig(args.output, dpi=300)