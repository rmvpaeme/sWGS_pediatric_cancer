awk 'FS=OFS="\t"{print $1, 0, $2}' GRCh38.fa.fai | bedtools  makewindows -b - -w 50000 > 50kbBins_GRCh38.bed
