err=0.02
cov=40
base="./"
petfile="read/rnapet.fa"
seqfile="read/rnaseq_200.fa"
echo 
./peta simu -c $cov -g read/gene_10.txt -e $err -m 200 -d 20 -b 1 -t read/tx.fa -p $base$petfile -l 100 -a 100 -s 0 read/chr10.fa
#echo 
#./peta simu -c $cov -g read/gene_10.txt -e $err -m 500 -s 50 -t read/tx.fa -p $base$petfile read/chr21.fa
#echo 
#./peta simu -c $cov -g read/gene_10.txt -e $err -m 1000 -s 100 -t read/tx.fa -p $base$petfile read/chr21.fa
