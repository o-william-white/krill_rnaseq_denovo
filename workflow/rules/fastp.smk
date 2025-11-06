rule fastp:
    input:
        sample=get_fastq,
    output:
        trimmed=["results/fastp/{sample}_R1.fastq", "results/fastp/{sample}_R2.fastq"],
        unpaired1="results/fastp/{sample}_U1.fastq",
        unpaired2="results/fastp/{sample}_U2.fastq",
        failed="results/fastp/{sample}_failed.fastq",
        html="results/fastp/{sample}_fastp.html",
        json="results/fastp/{sample}_fastp.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --trim_poly_g --trim_poly_x",
    threads: 8
    wrapper:
        "v7.1.0/bio/fastp"