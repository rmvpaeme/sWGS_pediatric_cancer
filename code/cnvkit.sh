#!/bin/bash
#PBS -N cnvkit
#PBS -l nodes=1:ppn=5
#PBS -l walltime=4:50:00

cd $PBS_O_WORKDIR
ml Python/3.8.2-GCCcore-9.3.0
ml R/4.0.0-foss-2020a
NORMAL_FOLDER=/data/gent/vo/000/gvo00027/vsc42220/2020_WES_pediatric_cancer_Curie/germline_WES/nodups/demultiplexed_reads/

# NK003
cnvkit.py batch  G28R3*bam --drop-low-coverage -m hybrid -t SureSelect_Clinical_Research_Exome_Regions_hg38.bed  \
			  -f hg38.fa \
			  -n $NORMAL_FOLDER/G28R4*bam  \
			  --access access-excludes.hg38.bed --output-reference flatreg.cnn -d test

#NK006
cnvkit.py batch  B192-B194E33*bam --drop-low-coverage -m hybrid -t SeqCap_EZ_Exome_v3_capture_hg38.bed  \
                          -f hg38.fa \
			    -n $NORMAL_FOLDER/B224E14*bam $NORMAL_FOLDER/B192-B194E19*bam \
                          --access access-excludes.hg38.bed --output-reference flatreg.cnn -d test

# NK025
cnvkit.py batch  B251E10*bam --drop-low-coverage -m hybrid -t SeqCap_EZ_Exome_v3_capture_hg38.bed  \
                          -f hg38.fa \
			   -n $NORMAL_FOLDER/B224E14*bam $NORMAL_FOLDER/B192-B194E19*bam \
                          --access access-excludes.hg38.bed --output-reference flatreg.cnn -d test

#UGS_352
cnvkit.py batch  B210E33*bam --drop-low-coverage -m hybrid -t MedExome_hg38_capture_targets.bed  \
                          -f hg38.fa \
			   -n $NORMAL_FOLDER/D275E05*bam $NORMAL_FOLDER/B202-B203E128*bam \
                          --access access-excludes.hg38.bed --output-reference flatreg.cnn -d test

# 33_190
cnvkit.py batch  D275E22*bam --drop-low-coverage -m hybrid -t MedExome_hg38_capture_targets.bed  \
                          -f hg38.fa \
			  -n $NORMAL_FOLDER/D275E05*bam $NORMAL_FOLDER/B202-B203E128*bam \
                          --access access-excludes.hg38.bed --output-reference flatreg.cnn -d test
