rule salmon_quant:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1=get_technical_replicates_forward,
        r2=get_technical_replicates_reverse,
        index="results/salmon_index",
    output:
        quant="results/salmon/{group}/quant.sf",
        lib="results/salmon/{group}/lib_format_counts.json",
        json="results/salmon/{group}/aux_info/meta_info.json",
    log:
        "logs/salmon/{group}.log",
    params:
        # optional parameters
        libtype="A",
        extra="",
    threads: 8
    wrapper:
        "v6.0.0/bio/salmon/quant"
