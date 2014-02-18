ulimit -c unlimited
query_dir=/home/hoiqz/shaojiang/peta/rnaseq/Spombe/SRR097897/
#query_dir=/home/carl/Projects/peta/rnaseq/Spombe/genome

../src/peta ass -k 25 -m 202 -s 62 -t 1 -o ../SRR097897_branch/ $query_dir/SRR097897.fa

query_dir=/home/carl/Projects/peta/rnaseq/hg19/SRX011545/
#../src/peta ass -k 25 -m 197 -s 62 -t 1 -o ../SRR027876_branch/ $query_dir/SRR027876.fa $query_dir/SRR097897.solid
