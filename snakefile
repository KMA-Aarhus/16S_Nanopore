# snakemake --profile configs/slurm


print("/*")
__author__ = "Tine Ebsen" # Please add your name here if you make changes.
__version__ = "0.2"

import sys
import os
from os import listdir
from os.path import isfile, isdir, join, exists, expanduser
import yaml
import pandas as pd
import numpy as np
from datetime import datetime
import glob
import time
import atexit
import datetime

#TO DO
# Merge reads
# Filtering with seqkit retaining read between 1300–1950bp
# Masking of repeats  using the TANTAN program?
# Map to human reference genome
# For each unmatched read, a minimap2 search with 5850 representative bacterial genome sequences stored in the GenomeSync database was performed.
# For each read, the species showing the highest minimap2 score were assigned to the query sequence. 
# When more than one species showed the same similarity score, the reads were classified at any higher taxonomic rank covering all the identified species. 
# Taxa were determined based on the NCBI taxonomy database [48]. Low-abundance taxa with less than 0.01% of total reads were discarded from the analysis.

configfile: "config.yaml"

# Actually given on the command line:
samplesheet = config["samplesheet"]
rundir = config["rundir"]
human_reference = config["human_reference"]
emu_db = config["emu_db"]
live_analysis = config["live_analysis"]




tab = "\t"
nl = "\n"


# Check that input was given.
if config["samplesheet"] == "NA":
    raise Exception("No samplesheet file was given. Please specify a samplesheet by appending --config samplesheet=\"path/to/samplesheet/\" to the command line call.")
if config["rundir"] == "NA":
    raise Exception("No rundir path was given. Please specify a rundir by appending --config rundir=\"path/to/rundir/\" to the command line call.")
# TODO: Implement additional input validation, like checking that the objects given are file and dir respectively.

print("   _   _                                       __   __  _____  KMA,AUH")
print("  | \\ | |                                     /_ | / / / ____|        ")
print("  |  \\| | __ _ _ __   ___  _ __   ___  _ __ ___| |/ /_| (___          ")
print("  | . ` |/ _` | '_ \\ / _ \\| '_ \\ / _ \\| '__/ _ \\ | '_ \\___ \\         ")
print("  | |\\  | (_| | | | | (_) | |_) | (_) | | |  __/ | (_) |___) |        ")
print("  |_| \\_|\\__,_|_| |_|\\___/| .__/ \\___/|_|  \\___|_|\\___/_____/         ")
print("                          | |                                         ")
print("                          |_|                                         ")


print(f"These are the parameters given:")
print(f"  samplesheet: {samplesheet}")
print(f"  rundir: {rundir}")
print(f"  human_reference: {human_reference}")
print()


#########################
# Parse the samplesheet #
#########################

samplesheet_extension = samplesheet.split(".")[-1]
print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

if samplesheet_extension == "xlsx":
    # Uses openpyxl
    df = pd.read_excel(samplesheet, dtype = str)

elif samplesheet_extension == "xls":
    df = pd.read_excel(samplesheet)

else:
    raise Exception(f"The spreadsheet file extension {samplesheet_extension} is not yet implemented.")

# Clean up the spreadsheet
print("Cleaning sample sheet ...                              ", end = "", flush = True)
df.columns = map(str.lower, df.columns) # Lowercase
df.columns = map(str.strip, df.columns) # Remove edge-spaces
df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
df = df.dropna(subset = ["sample_id"])# remove rows not containing a sample ID
print("✓")

# Check that the spreadsheet complies
print("Checking that the necessary columns exist ...          ", end = "", flush = True)
for i in ["barcode", "sample_id"]:
    if not i in df.columns:
        raise Exception(f"The sample sheet is missing a necessary column. The sample sheet must contain the column {i}, but it only contains {df.columns.tolist()}")
print("✓")

acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)]

print("Checking that the barcodes are correctly formatted and that the barcodes are unique... ", end = "", flush = True)
for i in df["barcode"]:
    #print("checking", i)
    if not i in acceptable_barcodes: 
        raise Exception(f"The given barcode \'{i}\' is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
if not len(df["barcode"]) == len(set(df["barcode"])):
    raise Exception(f"One or more barcodes are duplicated. Each barcode may only be used once")
print("✓")


print("Checking that the sample id's are unique ...           ", end = "", flush = True)
if not len(df["sample_id"]) == len(set(df["sample_id"])):
    raise Exception(f"One or more sample_id's are duplicated. Each sample_id may only be used once")
print("✓")

print()
print("These are the samples from the samplesheet you have given:")
print(df.to_string())
print("//")
print()

###################
# Validate rundir #
###################

# Wait for the rundir to occur in the specified path. 
# If it doesn't occur after a specified waiting time, then stop the p

if rundir[-1] == "/":
    print("Removing trailing slash from rundir")
    rundir = rundir[0:-1]


print("Checking that the rundir exists ...                    ", end = "", flush = True)
if not os.path.isdir(rundir):
    raise Exception(f"The rundir does not exist.")
print("✓")

print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)

for i in range(200):
    print("  Looking ... ", end = "", flush = True)
    fastq_pass_bases = glob.glob(rundir + "/**/fastq_pass", recursive = True) # Find any occurrence of the wanted path
    if len(fastq_pass_bases) == 0:
        print("nothing found yet, waiting 10 secs ...")
        time.sleep(10) # Wait 10 seconds.
    elif(i == 10):
        print() # clean newline
        raise Exception("nothing found after 10 tries. Aborting.")
    else: 
        print(f"Found                                    ✓")
        break


if not len(fastq_pass_bases) == 1:
    raise Exception(f"There seems to be more than one fastq_pass sub-directory beneath the given rundir. These paths were found:{nl} {str(nl + ' ').join(fastq_pass_bases)}{nl}Please specify a more specific rundir.")


fastq_pass_base = fastq_pass_bases[0]
del fastq_pass_bases

# base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
print(f"This is the batch base directory:{nl}  {base_dir}")



out_base = os.path.join(base_dir, "16s_output") # out_base is the directory where the pipeline will write its output to.


# And here is the code from the rule wait_for_minknow

sample_sheet_given_file = f"{fastq_pass_base}/../sample_sheet_given.tsv"
print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
df.to_csv(sample_sheet_given_file, sep = "\t")
print("✓")

print()

disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})

disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])



print("Continuing with the following barcodes:")

# the workflow_table is the table that contains the records where the barcode could be found on the disk.
workflow_table = disk_barcodes_df.merge(df, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
workflow_table = workflow_table.dropna(subset = ["sample_id"])

#print(workflow_table[["barcode", "sample_id", "type"]].to_string(index = False))
print(workflow_table)
print("//")
print()


####### Live analysis
if config["live_analysis"]:
# and exists(expanduser("/Users/admin1/Documents/KMA_AUH/Nanopore16s_workflow/LiveAnalysis.flag")):

    # Now we have all resources to start monitoring in the background
    #os.system("rm ~/Nanopore16S/LiveAnalysis.flag") 
    print(f"Starting analysis in the background ... ")
    live_analysis_sh = open("start_analysis.sh", "w")
    command = f"#!/bin/bash{nl}cd /Users/admin1/Documents/KMA_AUH/Nanopore16s_workflow/LiveAnalysis{nl}mamba activate live_analysis {nl}python live_analysis.py {samplesheet} {rundir} {emu_db}"
    print(command)
    live_analysis_sh.write(command)
    live_analysis_sh.close()
    # Start the sequence monitoring in a new terminal
    os.system("osascript -e \'tell app \"Terminal\" to do script \" source /Users/admin1/Documents/KMA_AUH/Nanopore16s_workflow/start_analysis.sh\"\'")
    #os.system("gnome-terminal --tab -- bash start_analysis.sh")

#############

# And here is the code from the rule wait_for_minknow
minutes_wait = 10
print("Checking that the sequencing_summary_*.txt-file has been written to disk ...")
while True:
    sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
    if len(sequencing_summary_file) == 0:
        print(f"  Still sequencing/basecalling; waiting {minutes_wait} minutes ...")
        time.sleep(60*minutes_wait)
    else:
        break

sequencing_summary_file = sequencing_summary_file[0]
print("  The sequencing summary has been found                ✓")
#print(f"  This is the sequencing_summary_*.txt-file: \"{sequencing_summary_file.split('/')[-1]}\"")
print(f"  This is the sequencing_summary_*.txt-file (full): \"{sequencing_summary_file}\"")

###############################
# Now we can run the analysis #
###############################


rule all:
    input:
        expand(["{out_base}/{sample_id}/merged_reads/{sample_id}.fastq.gz", \
                "{out_base}/{sample_id}/merged_reads/{sample_id}_nanostat", \
                "{out_base}/{sample_id}/merged_reads/{sample_id}_clean.fastq.gz", \
                "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz", \
                "{out_base}/{sample_id}/emu_abundance/{sample_id}_rel-abundance.tsv", \
                "{out_base}/family_abundance_table.csv", \
                "{out_base}/final_report.html" \
                ], \
                out_base = out_base, human_reference = human_reference, sample_id = df["sample_id"])





###########################
# Setup for data analysis #
###########################

rule merge_reads:
    input:
        barcode_dir = directory(lambda wildcards: workflow_table[workflow_table["sample_id"] == wildcards.sample_id]["barcode_path"].values[0])
    output: "{out_base}/{sample_id}/merged_reads/{sample_id}.fastq.gz"
    threads: 1
    shell: """

    cat {input.barcode_dir}/* > {output}

    """

rule nanostat:
    input:
        "{out_base}/{sample_id}/merged_reads/{sample_id}.fastq.gz"
    output: 
        "{out_base}/{sample_id}/merged_reads/{sample_id}_nanostat"
    conda: "configs/qc.yaml"
    threads: 1
    shell: """

    NanoStat --fastq {input} -o {out_base}/{wildcards.sample_id}/merged_reads/ -n {wildcards.sample_id}_nanostat

    """

rule remove_human:
    input:
        "{out_base}/{sample_id}/merged_reads/{sample_id}.fastq.gz"
    output: 
        clean_fastq = "{out_base}/{sample_id}/merged_reads/{sample_id}_clean.fastq.gz"
    conda: "configs/minimap2.yaml"
    params:
        human_reference = human_reference
    threads: 1
    shell: """
    
    minimap2 -ax map-ont {params.human_reference} {input} > {wildcards.sample_id}.contam.sam
    samtools fastq -n -f 4 {wildcards.sample_id}.contam.sam | gzip > {output.clean_fastq}
    rm {wildcards.sample_id}.contam.sam


    """


# Trim adapters
rule trim_adapt:
    input: 
        "{out_base}/{sample_id}/merged_reads/{sample_id}_clean.fastq.gz"
    output: 
        "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz"
    conda: "configs/porechop.yaml"
    threads: 4
    shell: """

    mkdir -p {out_base}/{wildcards.sample_id}/trimmed

    porechop -i {input} --format fastq.gz -t 4 -o {output}

    """



rule emu_abundance:
    input: 
        "{out_base}/{sample_id}/trimmed/{sample_id}_trimmed.fastq.gz"
    output: 
        rel_abundance = "{out_base}/{sample_id}/emu_abundance/{sample_id}_rel-abundance.tsv",
        alignments = "{out_base}/{sample_id}/emu_abundance/{sample_id}_emu_alignments.sam"
    conda: "configs/emu.yaml"
    params:
        emu_db = emu_db
    threads: 4
    shell: """
    mkdir -p {out_base}/{wildcards.sample_id}/emu_abundance
    emu abundance --type map-ont {input} --threads {threads} --keep-files --keep-counts --output-unclassified --db {params.emu_db} --output-dir {out_base}/{wildcards.sample_id}/emu_abundance --output-basename {wildcards.sample_id}

    """


rule aggregate_compare:
    input: 
        expand("{out_base}/{sample_id}/emu_abundance/{sample_id}_rel-abundance.tsv", out_base = out_base, sample_id = df["sample_id"])
    output: 
        full = "{out_base}/full_abundance_table.csv",
        species = "{out_base}/species_abundance_table.csv",
        genus = "{out_base}/genus_abundance_table.csv",
        family = "{out_base}/family_abundance_table.csv"
    threads: 1
    shell: """
    mkdir -p {out_base}/emu_abundance_tables
    find {out_base} -type f -name "*_rel-abundance.tsv" -exec cp {{}} {out_base}/emu_abundance_tables \\;
    Rscript scripts/aggregate_compare.R {out_base}/emu_abundance_tables {out_base}

    """

rule final_report:
    input:
        full = "{out_base}/full_abundance_table.csv",
        nanostat = expand("{out_base}/{sample_id}/merged_reads/{sample_id}_nanostat", out_base = out_base, sample_id = df["sample_id"])
    output:
        "{out_base}/final_report.html"
    conda: "configs/report.yaml"
    threads: 1
    shell: """
        Rscript -e \"rmarkdown::render(input = \'scripts/final_report.Rmd\',params = list(full_abundance_table = \'{input.full}\', samplesheet = \'{samplesheet}\', path = \'{out_base}/\'))\"
        mv scripts/final_report.html {out_base}

        """

