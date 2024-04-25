Analysis of 16S data from Oxford Nanopore

The workflow analysis 16S sequencing data from Oxford Nanopore. The final output is a report of reads across taxonomies including figures.

R
```
snakemake --config rundir=/Users/admin1/Documents/KMA_AUH/Nanopore16s_workflow/test_data  samplesheet=/Users/admin1/Documents/KMA_AUH/Nanopore16s_workflow/test_data/test_sample.xlsx  live_analysis=False --cores=1  --use-conda --conda-frontend mamba
}
```


# 16S_Nanopore
Input: 
- Raw nanopore fastq files as output from MinKnow. The directory must contain a folder titled "fastq_pass". This matches the default output format.   
- Samplesheet containing sample ids and barcodes. Barcodes must be in the format "NB01", "NB02" etc.  
Example:
![image](https://user-images.githubusercontent.com/90172976/157239868-b8989c11-0dce-4d4e-b886-5e89bc3bab1a.png)


Output:  
- A report containing QC and taksonomy results for each sample.

## Installation. 
Requires conda, openpyxl and snakemake to run.  
https://docs.conda.io/en/latest/miniconda.html. 
https://openpyxl.readthedocs.io/en/stable/.  
https://snakemake.readthedocs.io/en/stable/. 



### Install with git
```
git clone https://github.com/KMA-Aarhus/16S_Nanopore.git
```
## How to run
The workflow is made for running on a local workstation as it is able to produce fast results which can assign taxonomies before the analysis is done running. Input can be specificied in the config.yaml file or alternatively as part of the snakemake cmd. The workflow requires as input:
* A run directory (rundir)
* An excel samplesheet (samplesheet)
* An emu database that should cover all relevant taxa
* A human reference to filter out host reads
To run:
* Navigate to 16S_Nanopore directory where the snakefile is.  
* Make sure your rundir exists and only contain sequencing data from one run as the scripts does not allow overlaps in barcodes.
To run locally:
```
snakemake --config rundir=<path_to_rundir>  samplesheet=<path_to_samplesheet>/  live_analysis=False --cores=1  --use-conda --conda-frontend mamba
```
To utilise live_analysis which will produce results as soon as a fastq pass file is generated for a barcode, create a conda environment from the live_analysis.yaml file. Then add live_analysis=True to you command.
```
snakemake --config rundir=<path_to_rundir>  samplesheet=<path_to_samplesheet>/  live_analysis=True --cores=1  --use-conda --conda-frontend mamba
```

If sequencing data is already generated it is possible to run on GenomeDK. To do so Set up an alias for snakemake to run on slurm:
```
alias snakeslurm='mkdir -p logs/old; mv logs/*.{err,out} logs/old 2> /dev/null; snakemake --profile configs/slurm --use-conda --conda-frontend mamba'
```


### Output:  
- A report containing QC and taksonomy results for each sample.

