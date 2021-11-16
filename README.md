# RMarkdown
Data-analysis can be found at [link](./code/CNV.html).

# Data preprocessing
## sWGS
The fastq files were preprocessed with `./code/preprocessing/sWGS_pipeline_SE.snakefile`, which runs bwa-mem, picard markduplicates and WisecondorX/ichorCNA. The output (`bins.bed`, `segments.bed`) from WisecondorX were further used in the downstream processing. The reference for WisecondorX was generated from an in-house dataset of healthy volunteers. More information on how to obtain a reference dataset for WisecondorX is available at https://github.com/CenterForMedicalGeneticsGhent/WisecondorX. The PoN and supporting files for ichorCNA were obtained from the ichorCNA repository (https://github.com/broadinstitute/ichorCNA).

## WES
The fastq files were preprocessed with `./code/preprocessing/sWGS_pipeline_PE.snakefile`, which runs bwa-mem, picard markduplicates and WisecondorX/ichorCNA. 

CNVkit was run with the parameters described in `cnvkit.sh` on the deduplicated bam files with the paired germline sample. Target region bed files (e.g. SureSelect) was downloaded from the Agilent website. The link between the ID from germline WES and ID from tumor WES is also in the `cnvkit.sh` file. 

## Illumina HumanCytoSNP
Illumina Genomestudio 2.0 (https://www.illumina.com/techniques/microarrays/array-data-analysis-experimental-design/genomestudio.html) was used to obtain the log2ratio (Robs/Rexp) per bin from the IDAT files. The SampleSheets are available in `./resources/`. The other files (`.egt` and `.bpm` are too large to host here and are available on the Illumina website (e.g. https://emea.support.illumina.com/downloads/humancytosnp-12v2-1_product_files-ns.html)

## Agilent 180K array
Raw data for the Agilent 180K arrays was not available. The processed data (`sample_raw.csv`) was used for downstream processing. The data was processed according to the methods in https://pubmed.ncbi.nlm.nih.gov/23308108/. 

## Integrating the different platforms
`./code/convertRaw.snakefile`: converts data from obtained from WES, sWGS, Illumina HumanCytoSNP-12 or 180K array (Agilent) into a structure that allows the comparison of these different data-types, and harmonises the bin size between all data types. Also uses DNAcopy to obtain segmentations. If necessary (i.e. array data), the files are liftover to GRCh38.

Example input files are present in the `./examples/` folder. 

## Other scripts
`./code/makeBins.sh` obtain bins of choice starting from `.fa.fai` file, is a dependency of `./code/convertRaw.snakefile`

`./code/DNAcopy_segment.R`: runs DNAcopy, is a dependency of `./code/convertRaw.snakefile`
