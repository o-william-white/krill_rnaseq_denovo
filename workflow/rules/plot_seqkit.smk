rule plot_seqkit:
    input:
        fx2tab = "results/seqkit/fx2tab.txt",
    output:
        png = "results/seqkit/fx2tab_mqc.png",
    log: 
        "logs/plot_seqkit.txt",
    conda:
        "../envs/npm.yaml"
    shell:
        """
        python workflow/scripts/plot_transcript_length_gc.py --input {input.fx2tab} --output {output.png}
        """