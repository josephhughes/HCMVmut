# HCMVmut

The genomes for each patient were assembled using a combination of SPAdes and manual curation based on 
visualisation in Tablet. The detailed method is described in:

Suárez, N. M., Wilkie, G. S., Hage, E., Camiolo, S., Holton, M., Hughes, J., Maabar, M., Vattipally, S. B., Dhingra, A., Gompels, U. A., Wilkinson, G. W. G., Baldanti, F., Furione, M., Lilleri, D., Arossa, A., Ganzenmueller, T., Gerna, G., Hubáček, P., Schulz, T. F., … Davison, A. J. (2019a). Human Cytomegalovirus Genomes Sequenced Directly From Clinical Material: Variation, Multiple-Strain Infection, Recombination, and Gene Loss. The Journal of Infectious Diseases, 220(5), 781–791. [https://doi.org/10.1093/infdis/jiz208](https://doi.org/10.1093/infdis/jiz208)

These genomes are then used for identifying low frequency mutations in each of the HCMV patients with the HCMVmut pipeline.


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

 - Activate the conda environment
```
conda activate HCMVmut
```

 - Change paths in config.yml for you input files

 - Run the snakemake pipeline
```
snakemake -s HCMVmut.smk --cores 4
```

 - When finished, deactivate the conda environment
```
conda deactivate
```