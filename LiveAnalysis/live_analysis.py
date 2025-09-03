__author__ = "Tine Sneibjerg Ebsen"
__version__ = "0.2"

import sys
import os
from os import listdir
from os.path import isfile, isdir, join, exists
import pandas as pd
import numpy as np
import subprocess
import glob
import time

# new imports
from pathlib import Path



def start_covermon():
    tab = "\t"
    nl = "\n"

    # Defines some helper functions

    def samplesheet_parse(samplesheet):
        samplesheet_extension = samplesheet.split(".")[-1]
        print(f"Reading .{samplesheet_extension}-type sample sheet \"{samplesheet}\"")

        if samplesheet_extension == "xlsx":
                # Uses openpyxl
            df = pd.read_excel(samplesheet, dtype = str)
        elif samplesheet_extension == "xls":
            df = pd.read_excel(samplesheet)
        else:
            raise Exception(f"The spreadsheet must be excel formatted (.xlsx or .xls)")

        # Clean up the spreadsheet
        print("Cleaning sample sheet ...                              ", end = "", flush = True)
        df.columns = map(str.lower, df.columns) # Lowercase
        df.columns = map(str.strip, df.columns) # Remove edge-spaces
        df.columns = map(lambda x: str(x).replace(" ", "_"), df.columns) # Replace spaces with underscore
        df["barcode"] = df["barcode"].apply(np.vectorize(lambda x: str(x).strip().replace(" ", ""))) # Because we are later going to join using this column, it is necessary to strip it for spaces.
        df = df.dropna(subset = ["sample_id"])# remove rows not containing a barcode
        print("✓")
        print(df)
        return df

    def validate_samplesheet(df):
        # Check that the barcodes look correct
        acceptable_barcodes = [f"NB{i:02d}" for i in range(1,97)] + [f"RB{i:02d}" for i in range(1,97)]

        print("Checking that the barcodes are correctly formatted ... ", end = "", flush = True)

        for i in df["barcode"]:
            if not i in acceptable_barcodes: 
                raise Exception(f"The given barcode ", i, " is not an acceptable barcode. Here is a list of acceptable barcodes for inspiration:{nl} {' '.join(acceptable_barcodes)}")
        print("✓")


        print("Checking that the barcodes are unique ...              ", end = "", flush = True)
        if not len(df["barcode"]) == len(set(df["barcode"])):
            bc_counts = pd.DataFrame(df['barcode'].value_counts())
            bc_counts.columns = ["count"]
            bc_counts = bc_counts[bc_counts["count"] > 1]
            #print(nl, bc_counts)
            raise Exception(f"{nl}One or more barcodes are duplicated. Each barcode may only be used once:{nl}{bc_counts}")
        print("✓")

        print()
        print("These are the samples from the samplesheet you have given:")
        print(df.to_string())
        print("//")
        print()


    def validate_rundir(rundir):
        if rundir[-1] == "/":
            print("Removing trailing slash from rundir")
            rundir = rundir[0:-1]


        print("Checking that the rundir exists ...                    ", end = "", flush = True)
        if not os.path.isdir(rundir):
            raise Exception(f"The rundir does not exist.")
        print("✓")

        print(f"Looking for MinKNOW-characteristic output:") #, end = "", flush = True)
        # Wait for the rundir to occur in the specified path. 
        # If it doesn't occur after a specified waiting time, then stop the p
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
        print(f"Found the following fastq_pass base which will be given to CoverMon: {nl}  {fastq_pass_base}{nl}")


        # base_dir is the place where fastq_pass, fast5_pass and the sequencing summary resides.
        base_dir = os.path.dirname(fastq_pass_base) # This only works because there is NOT a trailing slash on the fastq_pass_base
        print(f"This is the batch base directory:{nl}  {base_dir}")

        return base_dir, fastq_pass_base

    def create_workflow_table(df):
        disk_barcodes_list  = sorted(glob.glob(fastq_pass_base + "/barcode*")) # Find all fastq_pass/barcode* directories
        disk_barcodes_df = pd.DataFrame({'barcode_path': disk_barcodes_list})



        disk_barcodes_df = disk_barcodes_df.assign(barcode_basename = [i.split("/")[-1] for i in disk_barcodes_df["barcode_path"]])
        if "RB" in df["barcode"][0]:
            disk_barcodes_df = disk_barcodes_df.assign(barcode = ["RB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
        elif "NB" in df["barcode"][0]:
            disk_barcodes_df = disk_barcodes_df.assign(barcode = ["NB" + i[-2:] for i in disk_barcodes_df["barcode_path"]])
        else:
            raise Exception(f"Barcodes in samplesheet are not acceptable")

        print("Continuing with the following barcodes:")

        # the workflow_table is the table that contains the records where the barcode could be found on the disk.
        workflow_table = disk_barcodes_df.merge(df, how='left', on='barcode') # left join (merge) the present barcodes onto the df table.
        workflow_table = workflow_table.dropna(subset = ["sample_id"])

        print(workflow_table)
        print("//")
        print()
        return workflow_table


    def update_plot(workflow_table, samplesheet,open_report):
        print("Scanning for new fastq files...")
        for index, row in workflow_table.iterrows():
            emu_out = out_base+row['sample_id']+"_rel-abundance.tsv"
            barcode_path = row['barcode_path']
            sample_id = row['sample_id']
            fastq_passed = listdir(row['barcode_path'])
            print(row['barcode_path'])
            print(fastq_passed)
            if not exists(emu_out) and fastq_passed != "":
                cat_cmd = f'cp {barcode_path}/{fastq_passed[0]} {out_base}{sample_id}.fastq.gz'
                porechop_cmd = f'porechop -i {out_base}{sample_id}.fastq.gz -o {out_base}{sample_id}_trimmed.fastq.gz'
                #Getting the below line to run was a pain. Hence, it is presented as a list.

                print(cat_cmd.split())
                subprocess.run(cat_cmd.split())
                print(porechop_cmd.split())
                subprocess.run(porechop_cmd.split())

                emu_cmd = f"emu abundance --type map-ont {out_base}{sample_id}_trimmed.fastq.gz --threads 1 --keep-files --keep-counts --output-unclassified --db {emu_db} --output-dir {out_base} --output-basename {sample_id}"
                print(emu_cmd.split())
                subprocess.run(emu_cmd.split())

                aggregate_cmd = f"Rscript scripts/aggregate_compare.R {out_base} {out_base}"
                print(aggregate_cmd.split())
                subprocess.run(aggregate_cmd.split())


                create_live_report = ["Rscript", "-e", "\"rmarkdown::render(input = ", "\'scripts/live_report.Rmd\',", "params = list(path = ", "\'"+out_base+"\', samplesheet = ", "\'"+samplesheet+"\'))\"" ]
                    
                print(create_live_report)
                subprocess.run(" ".join(create_live_report), shell=True)
                subprocess.run(['mv', 'scripts/live_report.html', out_base])

                # Starts browser-sync in a new terminal if the report is not open. This will not work on windows or macOS.
                print(f"Updated plot with sample {sample_id}")
            if not exists(out_base+"/live_report.html") and fastq_passed != "":
                aggregate_cmd = f"Rscript scripts/aggregate_compare.R {out_base} {out_base}"
                print(aggregate_cmd.split())
                subprocess.run(aggregate_cmd.split())
                create_live_report = ["Rscript", "-e", "\"rmarkdown::render(input = ", "\'scripts/live_report.Rmd\',", "params = list(path = ", "\'"+out_base+"\', samplesheet = ", "\'"+samplesheet+"\'))\"" ]
                    
                subprocess.run(" ".join(create_live_report), shell=True)
                subprocess.run(['mv', 'scripts/live_report.html', out_base])
            if open_report == False:
                print("open_report = ", open_report, ". Opening report")
                subprocess.run("gnome-terminal --tab -- browser-sync start -w --no-notify -s \"" + out_base +"\" --host 127.0.0.1 --port 9000 --index \"live_report.html\"", shell=True)        
                open_report = True
        return open_report

    #####################
    # Start the monitor #
    #####################

    # Parse and check arguments
    if len(sys.argv) < 2:
        raise Exception(f"Missing arguments. The script must contain (1) samplesheet, (2) path to run directory, (3) path to emu database.")

    samplesheet = sys.argv[1]
    rundir = sys.argv[2]
    emu_db = sys.argv[3]

    print(f"These are the parameters given:")
    print("This is the samplesheet: ", samplesheet)
    print("This is the run directory: ", rundir)


    #########################
    # Parse the samplesheet #
    #########################

    df = samplesheet_parse(samplesheet)

    validate_samplesheet(df)


    ###################
    # Validate rundir #
    ###################

    base_dir, fastq_pass_base = validate_rundir(rundir)

    ###########################
    # Create output directory #
    ###########################
    out_base = os.path.join(base_dir, "LiveAnalysis/") # out_base is the directory where the pipeline will write its output to.
    print("Creating output directory ", out_base,"...")
    subprocess.run(["mkdir","-p", out_base])
    print()

    sample_sheet_out = f"{out_base}/sample_sheet_given.tsv"
    print(f"Backing up the original sample sheet ...               ", end = "", flush = True)
    df.to_csv(sample_sheet_out, sep = "\t", index=False, na_rep='NA')
    print("✓")

    #########################
    # Create workflow table #
    #########################

    workflow_table = create_workflow_table(df)

    ##################
    # Start CoverMon #
    ##################

    # Set an open_report state to stop opening multiple reports
    print("Setting open_report to false")
    open_report = False

    # Keep track of processed files to avoid starting from scratch if script is terminated
    if exists(f"{out_base}/live_report.html") and open_report == False:
        # Starts browser-sync in a new terminal. This will not work on windows or macOS.
        subprocess.run("gnome-terminal --tab -- browser-sync start -w --no-notify -s \"" + out_base +"\" --host 127.0.0.1 --port 9000 --index \"live_report.html\"", shell=True)                
        open_report = True


    # When sequencing, we will check for new files every 60 seconds
    seconds_wait = 60

    # Set a sequencing state to recognise when to stop looking for new files
    still_sequencing = True

    while still_sequencing:
        # Scans for new files and updates the plot if any are found
        open_report = update_plot(workflow_table, sample_sheet_out,open_report)

        # Continue the monitor as long as the sequence summary does not exist. Wait <seconds_wait> between scans.
        sequencing_summary_file = glob.glob(base_dir + "/sequencing_summary_*.txt")
        if len(sequencing_summary_file) == 0:
            print(f"  Still sequencing/basecalling; waiting {seconds_wait} seconds before next scan ...")
            time.sleep(seconds_wait)
        else:
            still_sequencing = False


    sequencing_summary_file = sequencing_summary_file[0]
    print("  The sequencing summary has been found. Run complete    ✓")

if __name__ == "__main__":
    start_covermon()

