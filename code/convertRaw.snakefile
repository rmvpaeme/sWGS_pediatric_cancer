configfile: "config_CNV.yaml"
# Sample list
if config["input"] == "SNParray":
    baseIDS, = glob_wildcards("{sample}_snpArray.tsv")
elif config["input"] == "array":
    baseIDS, = glob_wildcards("{sample}_raw.csv")
elif config["input"] == "sWGS":
    baseIDS, = glob_wildcards("{sample}_bins.bed")
elif config["input"] == "WES":
    baseIDS, = glob_wildcards("{sample}.cnr")
# Import config and make report
#configfile: "config_RVP1908_vcf.yaml"

if config["input"] == "SNParray":
    rule all:
        input:
            expand("{sample}.done", sample = baseIDS)

    rule convert:
        input:
            "{sample}_snpArray.tsv"
        output:
            "{sample}.done"
        params:
            liftOver = "/Users/rmvpaeme/liftOver",
            chain = "/Users/rmvpaeme/GRCh37_to_GRCh38.chain.gz",
            binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/200kb_GRCh38.bed",
            #binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/15kbBins_GRCh38.bed",
            #binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/50kbBins_GRCh38.bed",
            DNAcopy = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/DNAcopy_segment.R",
            centromeres = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/centromere_repeats.hg38.txt",
            binsize = "200kb"
        shell:
            """
            tail -n +2 {input} > {wildcards.sample}_noH.tsv.tmp
            cut -f 4,5,5,15 {wildcards.sample}_noH.tsv.tmp > {wildcards.sample}_noH_small.tsv.tmp
            awk 'BEGIN{{FS=OFS="\t"}} {{$2 = $2 - 1 OFS $2}} 1' {wildcards.sample}_noH_small.tsv.tmp > {wildcards.sample}_noH_small.bed.tmp
            sort -k 1,1 -k2,2n {wildcards.sample}_noH_small.bed.tmp > {wildcards.sample}_noH_small_sort.bed.tmp
            {params.liftOver} {wildcards.sample}_noH_small_sort.bed.tmp {params.chain} {wildcards.sample}_noH_small_GRCh38.tsv.tmp  unmapped.tsv.tmp
            bedtools subtract -a {wildcards.sample}_noH_small_GRCh38.tsv.tmp -b {params.centromeres} > {wildcards.sample}_noH_small_GRCh38.mask.tsv.tmp
            bedtools intersect -a {wildcards.sample}_noH_small_GRCh38.mask.tsv.tmp  -b {params.binfile}  -wa -wb > {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp
            bedtools groupby -g 5,6,7 -c 4 -o mean -i {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp > {wildcards.sample}_{params.binsize}.tsv
            Rscript {params.DNAcopy} {wildcards.sample}_{params.binsize}.tsv
            bedtools intersect -a {wildcards.sample}_{params.binsize}_segments.tsv -b {params.binfile}  -wa -wb | awk '{{ print $5 "\t" $6 "\t" $7 "\t" $4 "\t"}}'  >{wildcards.sample}_segments_per_{params.binsize}_mask.tsv
            touch {output}
            """
elif config["input"] == "array":
    rule all:
        input:
            expand("{sample}.done", sample = baseIDS)

    rule convert:
        input:
            "{sample}_raw.csv"
        output:
            "{sample}.done"
        params:
            liftOver = "/Users/rmvpaeme/liftOver",
            chain = "/Users/rmvpaeme/GRCh37_to_GRCh38.chain.gz ",
            binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/200kb_GRCh38.bed",
            #binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/15kbBins_GRCh38.bed",
            #binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/50kbBins_GRCh38.bed",
            DNAcopy = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/DNAcopy_segment.R",
            centromeres = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/centromere_repeats.hg38.txt",
            binsize = "200kb"
        shell:
            """
            tr ',' '\t' < {input} > {wildcards.sample}.tsv
            tail -n +2 {wildcards.sample}.tsv > {wildcards.sample}_noH.tsv.tmp
            cut -f 1,2,3,4,5,6,7 {wildcards.sample}_noH.tsv.tmp > {wildcards.sample}_noH_small.tsv.tmp
            {params.liftOver} {wildcards.sample}_noH_small.tsv.tmp {params.chain} {wildcards.sample}_noH_small_GRCh38.tsv.tmp  unmapped.tsv.tmp
            bedtools subtract -a {wildcards.sample}_noH_small_GRCh38.tsv.tmp -b {params.centromeres} > {wildcards.sample}_noH_small_GRCh38.mask.tsv.tmp
            bedtools intersect -a {wildcards.sample}_noH_small_GRCh38.tsv.tmp  -b {params.binfile}  -wa -wb > {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp
            bedtools groupby -g 8,9,10 -c 4 -o mean -i {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp > {wildcards.sample}_{params.binsize}.tsv
            Rscript {params.DNAcopy} {wildcards.sample}_{params.binsize}.tsv
            bedtools intersect -a {wildcards.sample}_{params.binsize}_segments.tsv -b {params.binfile}  -wa -wb | awk '{{ print $5 "\t" $6 "\t" $7 "\t" $4 "\t"}}'  >{wildcards.sample}_segments_per_{params.binsize}_mask.tsv
            touch {output}
            """
elif config["input"] == "sWGS":
    rule all:
        input:
            "CPA_all.csv"

    rule convert:
        input:
            segs = "{sample}_segments.bed",
            bins = "{sample}_bins.bed",
        output:
            "{sample}_cpa.csv"
        params:
            binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/200kb_GRCh38.bed",
            binsize = "200kb",
            centromeres = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/centromere_repeats.hg38.txt",
            calcCPA = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/calcCPA.R"
        shell:
            """
            bedtools subtract -a {input.segs} -b {params.centromeres} > {wildcards.sample}_segments_mask.bed
            bedtools subtract -a {input.bins} -b {params.centromeres} > {wildcards.sample}_bins_mask.bed
            bedtools intersect -a {input.segs}  -b {params.binfile}  -wa -wb | awk '{{ print $6 "\t" $7 "\t" $8 "\t" $4 "\t"}}'  >{wildcards.sample}_segments_per_{params.binsize}_mask.tsv
            Rscript {params.calcCPA} {input.segs}
            touch {output}
            """

    rule paste:
        input:
            cpa = expand("{sample}_cpa.csv", sample = baseIDS)
        output:
            "CPA_all.csv"
        shell:
            """
            echo "SampleID,CPA,CPAm" | cat - {input.cpa} > {output}
            """

elif config["input"] == "WES":
    rule all:
        input:
            expand("{sample}.done", sample = baseIDS)

    rule convert:
        input:
            "{sample}.cnr"
        output:
            "{sample}.done"
        params:
            binfile = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/data/200kb_GRCh38.bed",
            DNAcopy = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/DNAcopy_segment.R",
            centromeres = "/Users/rmvpaeme/Repos/RVPCVP2012_sWGS/code/centromere_repeats.hg38.txt",
            binsize = "200kb"
        shell:
            """
            tail -n +2 {input} > {wildcards.sample}_noH.tsv.tmp
            cat {wildcards.sample}_noH.tsv.tmp | sed 's/^chr//' > {wildcards.sample}_noH.tsv_nochr.tmp
            awk '$6>-15' {wildcards.sample}_noH.tsv_nochr.tmp > {wildcards.sample}_noH.tsv_nochr_filt.tmp
            bedtools subtract -a {wildcards.sample}_noH.tsv_nochr_filt.tmp -b {params.centromeres} > {wildcards.sample}_noH_small_GRCh38.mask.tsv.tmp
            bedtools intersect -a {wildcards.sample}_noH_small_GRCh38.mask.tsv.tmp -b {params.binfile}  -wa -wb > {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp
            bedtools groupby -g 8,9,10 -c 6 -o mean -i {wildcards.sample}_noH_small_GRCh38_intersect.tsv.tmp > {wildcards.sample}_{params.binsize}.tsv
            Rscript {params.DNAcopy} {wildcards.sample}_{params.binsize}.tsv
            bedtools intersect -a {wildcards.sample}_{params.binsize}_segments.tsv -b {params.binfile}  -wa -wb | awk '{{ print $5 "\t" $6 "\t" $7 "\t" $4 "\t"}}'  >{wildcards.sample}_segments_per_{params.binsize}_mask.tsv
            touch {output}
            """
