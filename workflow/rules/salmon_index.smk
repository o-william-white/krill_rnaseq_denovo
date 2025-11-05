rule salmon_index:
    input:
        sequences=reference,
    output:
        multiext(
            "results/salmon_index/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        directory("results/salmon_index"), # need to specifiy output direcotry here so it can be picked up by salmon quant
    log:
        "logs/salmon_index.log",
    threads: 2
    params:
        # optional parameters
        extra="",
    wrapper:
        "v3.5.3/bio/salmon/index"
