query="ass_contigs"
db="cufflinks.fa"
ref="spombe/spombe.fa"
blastdb="../ncbi-blast-2.2.25+/db/"
reads_used="reads_used.fa"
cp read/$query.fa $blastdb
cp read/$db $blastdb

makeblastdb -in $blastdb/$db -out $blastdb/$db -dbtype nucl

echo "Blasting -perc_identity 96 -evalue 0.001 ..."

blastn -task blastn -query $blastdb$query.fa -db $blastdb$db -out $blastdb$query.blastn -num_threads 8 -perc_identity 90 -evalue 0.001

cp $blastdb$query.blastn read/

echo "Mapping contigs to assembled contigs..." 

bwa index read/$query.fa
bwa aln -t 8 read/$query.fa read/$reads_used > read/$reads_used.sai
bwa samse read/$query.fa read/$reads_used.sai read/$reads_used > read/$reads_used.sam

echo

echo "Mapping contigs to target transcripts..." 
#bwa aln -t 8 read/$db read/$query.fa > read/$query.sai
#bwa samse read/$db read/$query.sai read/$query.fa > read/$query.sam
bwa bwasw read/$db read/$query.fa > read/$query.sam
#samtools view -bS read/$query.sam > read/$query.bam
#samtools sort read/$query.bam read/$query.sorted
#samtools index read/$query.sorted.bam
echo "Done. Check sam file $query.sam"
