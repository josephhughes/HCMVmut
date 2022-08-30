# HCMVmut
Identifying low frequency mutations in the HCMV genome


 - Download the github repository
```
git clone https://github.com/josephhughes/HCMVmut.git
```

 - Navigate into the directory
```
cd HCMVmut
```


 - Install the environment
```
conda env create -f environment.yml
```
 - Change paths in config.yml for you input files

 - Run the snakemake pipeline
```
snakemake -s HCMVmut.smk --cores 4
```
