ml purge
ml snakemake/5.2.4-foss-2018b-Python-3.6.6
ml pygraphviz/1.5-foss-2018b-Python-3.6.6

snakefile="/data/gent/vo/000/gvo00027/vsc42220/sWGS_Pipeline/sWGS_pipeline_SE.snakefile"

if [[ $2 = "pe" ]];
  then
    echo "starting sWGS - PE pipeline "
    snakefile=snakefile="/data/gent/vo/000/gvo00027/vsc42220/sWGS_Pipeline/sWGS_pipeline_PE.snakefile"
  else
if [[ $2 = "se" ]];
  then
    echo "starting sWGS - SE pipeline"
  fi
fi

#source activate snakemake
snakemake -s $snakefile --cluster "qsub -V -l nodes=1:ppn={params.ppn} -l walltime={params.walltime} " --jobs 200 -R multiqc
current_time=$(date "+%Y.%m.%d-%H.%M.%S") &&
report_fileName=${current_time}-report.html &&
snakemake -s $snakefile --report $report_fileName
