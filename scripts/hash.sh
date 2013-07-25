query_dir=/home/ariyaratnep/shaojiang/peta/rnaseq/Spombe/SRR097897/
query_dir=/home/carl/Projects/peta/rnaseq/Spombe/genome/
#query_dir=/home/carl/Projects/peta/rnaseq/hg19/


#../src/peta k_hash -k 11 -l 68 -i 2 -b 6 -s 9 $query_dir/simu.pe.fa

../src/peta k_hash -k 11 -l 45 -i 2 -b 6 -s 4 $query_dir/simu.pe.fa

../src/peta group -t 4 $query_dir/simu.pe.fa 

