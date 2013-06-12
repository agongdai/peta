#!/bin/bash

query="single"
root_dir="/home/carl/Projects/"
eva_exe="$root_dir/peta_dev/scripts/eva.py"
db_dir="/home/carl/Projects/peta/rnaseq/Spombe/genome/"
db="spombe.broad.tx.fasta.rev"
db="tx.630.oracle.fa"
blastdb="$root_dir/ncbi-blast-2.2.26+/db/"
similarity="98"
query_dir="$root_dir/peta_dev/SRR097897_out/"
#query_dir="$root_dir/peta_dev/spombe_630/"

blastn_exe="blastn"
blat_exe="$root_dir/blat/blat"
occ="$root_dir/blat/11.ooc"
bwa_exe="bwa"

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 2 ] || die "2 arguments required, $# provided. run: sh eva.sh [aligner] [new/other]"

if [ $1 = "blastn" ]; 
	then
	if [ $2 = "new" ];
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
	echo "Evaluating..." 1>&2
	python eva.py blast -t $db_dir$db -c $query_dir$query.fa -b $query_dir$query.blastn -o $query_dir
fi

if [ $1 = "blat" ];
	then
	if [ $2 = "new" ];
		then
		echo "=======================================" 1>&2
		echo "$blat_exe $db_dir$db $query_dir$query.fa -ooc=$occ $query_dir$query.fa.psl"
		$blat_exe $db_dir$db $query_dir$query.fa -ooc=$occ $query_dir$query.fa.psl
		echo "=======================================" 1>&2
	fi
	echo "Evaluating..." 1>&2
	echo "python $eva_exe blat -t $db_dir$db -c $query_dir$query.fa -p $query_dir$query.fa.psl -o $query_dir -s 0.$similarity"
	python $eva_exe blat -t $db_dir$db -c $query_dir$query.fa -p $query_dir$query.fa.psl -o $query_dir -s 0.$similarity
fi

if [ $1 = "bwa" ];
	then
	if [ $2 = "new" ];
		then
		echo "bwa bwasw -t 6 $db_dir$db $query_dir$query.fa > $query_dir$query.tx.sam"
        bwa bwasw -t 6 $db_dir$db $query_dir$query.fa > $query_dir$query.tx.sam
	fi
fi

echo "=======================================" 1>&2
echo "Done"


