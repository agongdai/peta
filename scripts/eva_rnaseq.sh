#!/bin/bash

query="peta.1"
db_dir="../rnaseq/Spombe/genome/"
db="spombe.broad.tx.fasta"
blastdb="../../ncbi-blast-2.2.26+/db/"
similarity=90
query_dir=../SRR097897_out/

if [ $1 = "new" ]; 
	then
	cp $query_dir/$query.fa $blastdb
	cp $db_dir$db $blastdb 
	makeblastdb -in $blastdb/$db -out $blastdb/$db -dbtype nucl

	echo "=======================================" 1>&2
	echo "blastn -task blastn -query $blastdb$query.fa -db $blastdb$db -out $blastdb$query.blastn -outfmt \"7 qacc sseqid pident length mismatch gaps sstart send evalue qstart qend\" -num_threads 8 -perc_identity $similarity -evalue 0.001" 1>&2
	blastn -task blastn -query $blastdb$query.fa -db $blastdb$db -out $blastdb$query.blastn -outfmt "7 qacc sseqid pident length mismatch gaps sstart send evalue qstart qend" -num_threads 8 -perc_identity $similarity -evalue 0.001
	echo "=======================================" 1>&2

	echo
	cp $blastdb$query.blastn $query_dir
fi

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
#
echo "Evaluating..." 1>&2
echo "=======================================" 1>&2
echo "python eva.py blast -t $db_dir$db -c $query_dir$query.fa -b $query_dir$query.blastn -o $query_dir"
python eva.py blast -t $db_dir$db -c $query_dir$query.fa -b $query_dir$query.blastn -o $query_dir
echo "=======================================" 1>&2
if [ $1 = "bwa" ]; 
	then
	echo "bwa bwasw -t 2 $db_dir$db $query_dir$query.fa > $query_dir$query.tx.sam"
	bwa bwasw -t 6 $db_dir$db $query_dir$query.fa > $query_dir$query.tx.sam
fi
echo "=======================================" 1>&2
#echo "python eva.py bwa -t $db_dir$db -c $query_dir$query.fa -s $query_dir$query.tx.sam -o $query_dir"
#python eva.py bwa -t $db_dir$db -c $query_dir$query.fa -s $query_dir$query.tx.sam -o $query_dir
#eva read/ass_contigs.blastn read/result.txt graph/gene_10_lengths.txt graph/gene_10_edges.txt read/ass_contigs.fa read/tx.fa read/start_reads_gene.txt
#echo "../peta eva -o 0 -m $query_dir/$query.blastn -r $query_dir -c $query_dir/$query.fa -t $db_dir/$db" 1>&2
#../src/peta eva -o 0 -m $query_dir/$query.blastn -r $query_dir -c $query_dir/$query.fa -t $db_dir/$db
#echo "=======================================" 1>&2
#echo "See result in $query_dir/result.txt" 1>&2
echo "Done"

#echo
#echo "Preparing the sam file..." 1>&2
#bwa bwasw -t 8 $query_dir/$db $query_dir/$query.fa > $query_dir/$query.sam
#echo "Done. Check sorted bam file $query.sam" 1>&2

