# HCMVmut
Identifying low frequency mutations in the HCMV genome


Download the github repository

Navigate into the directory

Install the environment
```
conda env create -f environment.yml
```

Run the snakemake pipeline
```
snakemake -s HCMVmut.smk --cores 4
```
