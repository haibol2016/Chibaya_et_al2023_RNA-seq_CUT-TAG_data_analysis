## adding group labels to peakID

for i in  *.narrowPeak; do 
   base=`basename $i .q005.Genrich.narrowPeak`; awk 'BEGIN{FS=OFS="\t"}{$4=base"."$4; print}' base=$base $i > ${i}.bed ;

 done

## merge peaks using bedtools merge
sort -k1,1 -k2,2n *.bed | bedtools merge  -c 4 -o collapse,count  -i - > q005.merged.peaks.bed

## convert bed to gtf
perl  ~/michael_lodato/Narendra/scripts/016.bed.gtf.pl q005.merged.peaks.bed  > ../q005.all.samples.merged.peaks.gtf

## format header lines
grep -v '^#' q005.all.samples.merged.peaks.Marcus.ATACseq.count.table |perl -p -e  's{results/009.samtools.merge.out/}{}g; s{_S0_L001.filted.name.sort.merge.bam}{}g; s{.20049D-07}{}g; s{-}{_}g' > q005.all.samples.merged.peaks.refmt.count.table

awk 'BEGIN{FS=OFS="\t"} {print $1, $2+1, $3, $4, $5}' q005.merged.peaks.bed | awk 'BEGIN{FS=OFS="\t"} NR ==FNR {a[$1"_"$2"_"$3] = $4"\t"$5; next} {if ($2"_"$3"_"$4 in a) {print a[$2"_"$3"_"$4], $0} else if (/^Geneid/) {print "MergepeakID\tPeakNumber", $0}}' -   q005.all.samples.merged.peaks.refmt.count.table > q005.all.samples.merged.peaks.final.refmt.count.table
