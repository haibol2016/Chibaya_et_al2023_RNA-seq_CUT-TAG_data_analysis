samtools view -H NW01-1A.sort.bam |grep '@SQ' |grep -v 'chrM' | perl -p -e 's{.+?SN:([^\t]+).+}{$1}' > ../../docs/human.gencodeV34.nuc.chr.list
