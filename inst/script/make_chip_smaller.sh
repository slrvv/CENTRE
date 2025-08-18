##script to reduce the example files

wd=/project/CRUP_scores/CENTRE/inst/extdata

for hist in H3K4me3 H3K4me1 H3K27ac input; do
  bedtools intersect -abam $wd/example/HeLa_$hist.REF_chr19.bam \
  -b $wd/regions.bed > $wd/example/HeLa_$hist.REF_chr19_reduced.bam
done

for hist in H3K4me3 H3K4me1 H3K27ac input; do
  samtools index $wd/example/HeLa_$hist.REF_chr19_reduced.bam
done

