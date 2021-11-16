# Import config and make report
configfile: "config.yaml"
report: "report/workflow.rst"

IDS, = glob_wildcards("{sample}_R1.fastq.gz")

rule all:
    input:
        "multiqc_report.html"

rule mapping_bwa:
    input:
        fq = "{sample}_R1.fastq.gz"
    output:
        "mapped_reads/{sample}_R1.bam"
    threads: 8
    params:
        walltime="20:00:00",
        ppn = 10,
        genome = config["humangenome"],
        outdir = "mapped_reads"
    shell:
        "ml purge; ml FastQC/0.11.9-Java-11 ; "
        "fastqc {input.fq} --outdir {params.outdir} ; "
        "ml purge; ml BWA/0.7.17-intel-2018a; ml SAMtools/1.8-intel-2018a; "
        "bwa mem -t {threads} {params.genome} {input.fq} | samtools sort -O bam -o {output}"

rule deduplicate:
    input:
        "mapped_reads/{sample}_R1.bam"
    output:
        bam = "nodups/{sample}_nodups.bam",
        metrics = "logs/{sample}_dup_metrics.txt"
    params:
        walltime="20:45:00",
        ppn = 14,
    shell:
        "ml purge && ml picard/2.21.6-Java-11; "
        "java -jar $EBROOTPICARD/picard.jar MarkDuplicates "
        "I={input} "
        "O={output.bam} "
        "M={output.metrics}  "
        "VALIDATION_STRINGENCY=SILENT "
        "REMOVE_DUPLICATES=true; "
        "ml purge && ml SAMtools/1.9-intel-2018b; "
        "samtools index {output.bam}"

rule convertnpz:
    input:
        "nodups/{sample}_nodups.bam"
    output:
        "convertnpz/{sample}.npz"
    params:
        walltime="1:20:00",
        ppn = 1,
        binsize = config["wisexkb"] * 1000,
        wiseXref = config["wisexref"]
    shell:
        "ml purge && ml WisecondorX/1.1.6-foss-2020a-Python-3.8.2;  "
        "WisecondorX convert --binsize {params.binsize} {input} {output}"

rule plotnpz:
    input:
        "convertnpz/{sample}.npz"
    output:
        bed = "plotnpz/{sample}_bins.bed"
    params:
        walltime="1:20:00",
        ppn = 1,
        binsize = config["wisexkb"] * 1000,
        wiseXref = config["wisexref"],
        blacklist = config["blacklist"],
        outputid = "plotnpz/{sample}",
    shell:
        "ml purge && ml WisecondorX/1.1.6-foss-2020a-Python-3.8.2; "
        "WisecondorX predict {input} {params.wiseXref} {params.outputid} --blacklist {params.blacklist} --plot --bed --cairo --ylim [-2,2] "

rule hmmcopy:
    input:
        "nodups/{sample}_nodups.bam"
    output:
        "wigfiles/{sample}.wig"
    params:
        walltime="0:20:00",
        ppn = 1,
        hmmcopy = config["hmmcopy"],
        binsize = config["ichorkb"] * 1000
    shell:
        "{params.hmmcopy} --window {params.binsize} --quality 20 "
        "--chromosome 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' "
        "{input} > {output}"

rule ichorCNA_high:
    input:
        tum = "wigfiles/{sample}.wig"
    output:
        directory("ichorCNA_high/{sample}/")
    params:
        walltime="1:40:00",
        ppn = 1,
        runichorCNA = config["runichorCNA"],
        ichorCNAfolder = config["ichorCNAfolder"],
        PoNichorCNA = config["PoNichorCNA"],
        ichorgc_wig = config["gcwig"],
        ichormap_wig = config["mapwig"],
        kb = config["ichorkb"],
        id = "{sample}"
    shell:
        "ml purge && ml ichorCNA/0.3.2-20191219-foss-2020a;   "
        "mkdir -p {output} ; "
        "Rscript {params.runichorCNA} --id {params.id} "
        "--WIG {input.tum} "
        "--gcWig {params.ichorgc_wig} "
        "--mapWig {params.ichormap_wig} "
        "--normalPanel {params.PoNichorCNA} "
        "--chrs 'c(1:22)' "
        "--chrTrain 'c(1:22)' "
        "--outDir {output} "
        "--ploidy 'c(2:4)' "
        "--scStates 'c()' --estimateScPrevalence FALSE --estimatePloidy TRUE "
        "--txnE 0.99 --txnStrength 100 "
        "--normal 'c(.01,.25,.5,.75,.99)' --maxCN 5 "

rule ichorCNA_low:
    input:
        tum = "wigfiles/{sample}.wig"
    output:
        directory("ichorCNA_low/{sample}/")
    params:
        walltime="1:40:00",
        ppn = 1,
        runichorCNA = config["runichorCNA"],
        ichorCNAfolder = config["ichorCNAfolder"],
        PoNichorCNA = config["PoNichorCNA"],
        ichorgc_wig = config["gcwig"],
        ichormap_wig = config["mapwig"],
        kb = config["ichorkb"],
        id = "{sample}"
    shell:
        "ml purge && ml ichorCNA/0.3.2-20191219-foss-2020a;   "
        "mkdir -p {output} ; "
        "Rscript {params.runichorCNA} --id {params.id} "
        "--WIG {input.tum} "
        "--gcWig {params.ichorgc_wig} "
        "--mapWig {params.ichormap_wig} "
        "--normalPanel {params.PoNichorCNA} "
        "--chrs 'c(1:22)' "
        "--chrTrain 'c(1:22)' "
        "--outDir {output} "
        "--ploidy 'c(2)' "
        "--scStates 'c()' --estimateScPrevalence FALSE --estimatePloidy TRUE "
        "--txnE 0.9999 --txnStrength 10000 "
        "--normal 'c(.95,.99,.995,.999)' --maxCN 3 "

rule collectLogs:
    input:
        "nodups/{sample}_nodups.bam"
    output:
        "logs/{sample}_samtools.stats"
    params:
        walltime="1:40:00",
        ppn = 1
    shell:
        "ml purge && ml SAMtools/1.10-iccifort-2019.5.281;"
        "samtools stats {input} > {output};"

rule multiqc:
    input:
        expand("logs/{sample}_samtools.stats", sample = IDS),
        expand(directory("ichorCNA_low/{sample}/"), sample = IDS),
        expand(directory("ichorCNA_high/{sample}/"), sample = IDS),
        expand("plotnpz/{sample}_bins.bed", sample = IDS)
    output:
        report("multiqc_report.html")
    params:
        walltime="1:20:00",
        ppn = 1
    shell:
        "ml purge && ml Python/3.6.6-intel-2018b; "
        "multiqc -f ."
