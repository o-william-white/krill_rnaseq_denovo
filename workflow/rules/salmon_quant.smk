rule salmon_quant:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r1=get_technical_replicates_forward,
        r2=get_technical_replicates_reverse,
        index=config["salmon_index"],
    output:
        quant="results/salmon/{group}/quant.sf",
        lib="results/salmon/{group}/lib_format_counts.json",
        json="results/salmon/{group}/aux_info/meta_info.json",
    log:
        "logs/salmon/{group}.log",
    conda:
        "../envs/salmon.yaml"
    threads: 8
    shell:
        """
        salmon quant \
            -i {input.index} \
            -l A \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            -o results/salmon/{wildcards.group} >{log} 2>&1
        """
