#!/bin/bash
query="peta"
db="spombe_genes.fa"
start_reads="read/start_reads.txt"
start_reads_gene="read/start_reads_gene.txt"
blastdb="../ncbi-blast-2.2.25+/db/"
similarity=90
cp read/$query.fa $blastdb
cp read/$db $blastdb

makeblastdb -in $blastdb/$db -out $blastdb/$db -dbtype nucl

echo "Blasting $query.fa to $blastdb/$db: -perc_identity $similarity -evalue 0.001 ..." 1>&2

#blastn -task blastn -query $blastdb$query.fa -db $blastdb$db -out $blastdb$query.blastn -outfmt "7 qacc sseqid pident length mismatch gaps sstart send evalue qstart qend" -num_threads 8 -perc_identity $similarity -evalue 0.001

echo
cp $blastdb$query.blastn read/

#while read LINE
#do
#    nos=$nos" "$LINE
#done < $start_reads

#echo "Getting target transcripts..." 1>&2

#echo "-------------------------------------------------------------------"
#echo "Getting the starting reads..."
#echo "./show_gene_no.pl $nos"
#./show_gene_no.pl $nos
#echo "-------------------------------------------------------------------"

echo "Evaluating..." 1>&2
#eva read/ass_contigs.blastn read/result.txt graph/gene_10_lengths.txt graph/gene_10_edges.txt read/ass_contigs.fa read/tx.fa read/start_reads_gene.txt
./peta eva -o 0 -m read/$query.blastn -r read/$query.result.txt -e graph/gene_10_lengths.txt -g graph/gene_10_edges.txt -c read/$query.fa -t read/$db -s $start_reads_gene
echo "See result in $query.result.txt" 1>&2

echo
echo "Preparing the sam file..." 1>&2
bwa bwasw -t 8 read/$db read/$query.fa > read/$query.sam
echo "Done. Check sorted bam file $query.sam" 1>&2
