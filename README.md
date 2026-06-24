
Conda env setup
```bash
conda env create -n rnaseq_denovo_env -f workflow/envs/conda_env.yaml

conda activate rnaseq_denovo_env
```

Dry run 
```bash 
snakemake --profile workflow/profiles/ --dry-run
```

Build rule graph 
```bash 
snakemake --profile workflow/profiles/ --rulegraph | dot -Tsvg > rulegraph.svg
```

Run pipeline
```
run_rnaseq_denovo_Ecry.sh
```
