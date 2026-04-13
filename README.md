# nf-varcall
Nextflow pipeline for calling variants on population genomics datasets


## Nextflow

Nextflow is installed as a module on dardel, and will by default be loaded inside the run_on_slurm.sh script.

If not on dardel or just want your own installation, install nextflow as per instructions https://www.nextflow.io/docs/latest/install.html, and remove the ```module load Nextflow``` line from run_on_slurm.sh.
    
    

## Creating environment

All software dependencies of this pipeline have been packaged in a conda environment. To generate the conda environment and install all the necessary software, run:

    conda env create -f environment.yml

    # once this is done, nextflow need the path to the conda environment at the top of the nextflow.config file, can be set with this command:

    sed -i "s|CONDA_ENVIRONMENT_PATH|$(conda env list | awk ' $1 == "varcall_env" {print $2} ')|" nextflow.config

