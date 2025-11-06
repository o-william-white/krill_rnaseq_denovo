rule fastqc_fwd:
    input:
        "results/fastp/{sample}_R1.fastq"
    output:
        html="results/fastqc/{sample}_R1.html",
        zip="results/fastqc/{sample}_R1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/{sample}.log"
    threads: 8
    resources:
        mem_mb = 1024
    wrapper:
        "v7.6.0/bio/fastqc"

rule fastqc_rev:
    input:
        "results/fastp/{sample}_R2.fastq"
    output:
        html="results/fastqc/{sample}_R2.html",
        zip="results/fastqc/{sample}_R2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        extra = "--quiet",
        mem_overhead_factor=0.1,
    log:
        "logs/fastqc/{sample}.log"
    threads: 8
    resources:
        mem_mb = 1024
    wrapper:
        "v7.6.0/bio/fastqc"