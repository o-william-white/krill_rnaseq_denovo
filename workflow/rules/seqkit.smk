rule seqkit_stats:
    input:
        fastx=reference,
    output:
        stats="results/seqkit/stats_mqc.txt",
    log:
        "logs/seqkit/stats.log",
    params:
        command="stats",
        extra="--all --tabular",
    threads: 2
    wrapper:
        "v7.6.0/bio/seqkit"

rule seqkit_fx2tab:
    input:
        fastx=reference,
    output:
        tsv="results/seqkit/fx2tab_mqc.txt",
    log:
        "logs/seqkit/fx2tab.log",
    params:
        command="fx2tab",
        extra="--name --length --gc",
    threads: 8
    wrapper:
        "v7.6.0/bio/seqkit"

