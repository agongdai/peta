#query_dir=../rnaseq/Spombe/SRR097897
#../src/peta pair -p 0 -k 33 -m 197 -s 62 -t 6 -d 1 -b /home/carl/Projects/blat/blat -o ../SRR097897_out/ $query_dir/SRR097897.fa $query_dir/SRR097897.solid 
#query_dir=../rnaseq/hg19/SRX011545/
#../src/peta pair -p 0 -k 27 -m 174 -s 80 -t 5 -d 1 -b /home/carl/Projects/blat/blat -o ../SRR027876_out/ $query_dir/SRR027876.trim2.fa $query_dir/SRR027876.trim2.solid 
query_dir=../rnaseq/hg19/SRX011545/
../src/peta pair -p 0 -k 27 -m 174 -s 80 -t 5 -d 1 -b /home/carl/Projects/blat/blat -o ../SRR027876_trim/ $query_dir/SRR027876.trim2.fa $query_dir/SRR027876.trim2.solid 

