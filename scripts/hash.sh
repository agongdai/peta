query_dir=/home/ariyaratnep/shaojiang/peta/rnaseq/Spombe/SRR097897/
query_dir=/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/

#../src/peta hash -k 11 -l 45 -i 2 -b 5 -s 4 /home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876_corrected.fa
../src/peta k_hash -k 11 -l 68 -i 2 -b 6 -s 9 $query_dir/simu.fa
#./peta hash -k 7 -l 67 -i 2 -b 12 -s 2 read/rnapet.fa
#./peta hash -k 11 -l 45 -i 2 -b 5 -s 5 rnaseq/SRX011545.fasta read/rnapet.fa
#./peta hash -k 11 -l 45 -i 2 -b 5 -s 5 rnaseq/SRR027876.fasta read/rnapet.fa
#../src/peta hash /home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa 25
#../src/peta hash /home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa 25
#../src/peta hash /home/carl/Projects/peta/rnaseq/Spombe/genome/simu.fa 25
