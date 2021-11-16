# sWGS_Pipeline

## How to run
### Initial setup
1. Navigate to a folder and `git clone git@github.ugent.be:DePreterLab/sWGS_Pipeline.git`. This is the pipeline directory.
1. Change line 5 in `runsWGSPipeline.sh` so that it reflects your current configuration.

### Initiating the pipeline
1. Navigate to a new folder. This is the analysis directory. Add your reads to a subfolder called "demultiplexed_reads".
1. Make sure your reads are named `{samplename}_R1.fastq.gz`. The {samplename} can be anything but the filenames should end on `_R1.fastq.gz`
1. Move the example `config.yaml` file from the pipeline directory to the analysis directory and edit so that it reflects your configuration.
1. (optional but highly recommended) open a new screen session with `screen`.
1. Start the pipeline by navigating to the analysis directory and run `bash {pipeline_directory}/runsWGSPipeline.sh`. Jobs are constantly being submitted.
1. If you keep the screen session open and your connection is terminated (timeout or if you close your laptop) the job submission stops. To avoid this, close the screen session with ctrl+A followed by ctrl+D and reattach it with `screen -ls ; screen -r {screen_ID}` to inspect the status of the pipeline. 

## Troubleshooting
### Insufficient walltime, out-of-memory errors
Adjust the walltime or ppn (in case of OOM errors) in the snakefile.
### I get an error, but don't know whats happening
Snakemake will report the external jobid of the error. You can inspect it with `cat *e{jobid}`/`cat *o{jobid}` e.g. `cat *e9549216*`. 
### I cant find the external jobid of the error
In the analysis directory, there is a hidden `.snakemake/logs` folder. Find the most recent log with `ls -laht | head` and use `grep -B 5 -A 5 "Error" logfile.log`.
### Nothing happens after starting the script
1. The filenames are wrong, or the files are not in a folder called "demultiplexed_reads".
1. The pipeline was initiated earlier and finished without errors. If you want to force a restart, it's better to remove the multiqc report file and folder, along with all the (except demultiplexed_reads) folders.
