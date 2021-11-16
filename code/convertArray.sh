
filename=D1209065
tr ',' '\t' < ${filename}_raw.csv > ${filename}_raw.tsv
tail -n +2 ${filename}_raw.tsv > ${filename}_raw_noH.tsv
cut -f 1,2,3,4,5,6,7 ${filename}_raw_noH.tsv > ${filename}_raw_noH_small.tsv
/Users/rmvpaeme/liftOver  ${filename}_raw_noH_small.tsv /Users/rmvpaeme/GRCh37_to_GRCh38.chain.gz ${filename}_raw_noH_GRCh38.tsv unmapped.tsv
bedtools intersect -a ${filename}_raw_noH_GRCh38.tsv -b ../200kb_GRCh38.bed  -wa -wb > ${filename}_raw_intersect.tsv
bedtools groupby -g 8,9,10 -c 4 -o mean -i ${filename}_raw_intersect.tsv > ${filename}_raw_200kb.tsv
Rscript ../../code/DNAcopy_segment.R ${filename}_raw_200kb.tsv

CFD1800942 vs D1209065
