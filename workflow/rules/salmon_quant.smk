rule salmon_quant:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1="results/fastp/{sample}_R1.fastq",
        r2="results/fastp/{sample}_R2.fastq",
        index="results/salmon_index",
    output:
        quant="results/salmon/{sample}/quant.sf",
        lib="results/salmon/{sample}/lib_format_counts.json",
        json="results/salmon/{sample}/aux_info/meta_info.json",
    log:
        "logs/salmon/{sample}.log",
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 2
    wrapper:
        "v6.0.0/bio/salmon/quant"
