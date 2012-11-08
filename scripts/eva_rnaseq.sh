#!/bin/bash
query="peta"
db="spombe.broad.tx.fa"
blastdb="../ncbi-blast-2.2.25+/db/"
similarity=90
query_dir=SRR097897/
cp $query_dir/$query.fa $blastdb
cp $query_dir/$db $blastdb 

makeblastdb -in $blastdb/$db -out $blastdb/$db -dbtype nucl

echo "Blasting $query.fa to $blastdb/$db: -perc_identity $similarity -evalue 0.001 ..." 1>&2

blastn -task blastn -query $blastdb$query.fa -db $blastdb$db -out $blastdb$query.blastn -outfmt "7 qacc sseqid pident length mismatch gaps sstart send evalue qstart qend" -num_threads 8 -perc_identity $similarity -evalue 0.001

echo
cp $blastdb$query.blastn $query_dir

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

echo "Done"
echo "Evaluating..." 1>&2
#eva read/ass_contigs.blastn read/result.txt graph/gene_10_lengths.txt graph/gene_10_edges.txt read/ass_contigs.fa read/tx.fa read/start_reads_gene.txt
echo "./peta eva -o 0 -m $query_dir/$query.blastn -r $query_dir/$query.result.txt -e graph/gene_10_lengths.txt -g graph/gene_10_edges.txt -c $query_dir/$query.fa -t $query_dir/$db -s $start_reads_gene" 1>&2
./peta eva -o 0 -m $query_dir/$query.blastn -r $query_dir/$query.result.txt -c $query_dir/$query.fa -t $query_dir/$db
echo "See result in $query_dir/$query.result.txt" 1>&2

#echo
#echo "Preparing the sam file..." 1>&2
#bwa bwasw -t 8 $query_dir/$db $query_dir/$query.fa > $query_dir/$query.sam
#echo "Done. Check sorted bam file $query.sam" 1>&2

