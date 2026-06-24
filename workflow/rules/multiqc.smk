rule multiqc:
    input:
        expand(
            "results/fastp/{sample}_fastp.json",
            sample=sample_data.index.tolist(),
        ),
        expand(
            "results/fastqc/{sample}_R1.html",
            sample=sample_data.index.tolist(),
        ),
        expand(
            "results/fastqc/{sample}_R2.html",
            sample=sample_data.index.tolist(),
        ),
        expand(
            "results/salmon/{group}/aux_info/meta_info.json",
            group=dict_groups.keys(),
        ),
        config="config/config_multiqc.yaml",
    output:
        "results/multiqc.html",
        directory("results/multiqc_data"),
    log:
        "logs/multiqc.log",
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    wrapper:
        "v7.5.0/bio/multiqc"
