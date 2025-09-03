#!/bin/bash


# You may source this bashrc-extension, if you wish to have some shortcuts which can be nice for development and deployment.



# Neat aliases
alias mconda='mamba'
alias survey='watch "sensors; nvidia-smi"'
alias citament='git add -u && git commit -m "amend" &&  git pull && git push && echo OK'



start_16s () {

	example="\n\texample:\n\t start_16s path/to/my_samplesheet.xlsx path/to/my_rundir\n"

if [ -z "$1" ]; then
    echo "Input error: The samplesheet argument is empty. Please specify a samplesheet."
    echo -e $example
     
elif [ -z "$2" ]; then
    echo "Input error: The rundir argument is empty. Please specify a rundir."
    echo -e $example
else

    #clear
    cd ~/16S_Nanopore&& conda activate 16S_Nanopore && snakemake --use-conda --conda-frontend mamba --config samplesheet="${1}" rundir="${2}" live_analysis=True ${3} --cores 4 && echo && cowsay -f fish "16S_Nanopore finished successfully, you may now close this window."
fi
}



