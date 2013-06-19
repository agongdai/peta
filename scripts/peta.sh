#query_dir=../rnaseq/Spombe/SRR097897
#./peta pair -p 0 -k 39 -m 197 -d 62 -t 1 -o ../SRR097897_out/ $query_dir/SRR097897.fa $query_dir/SRR097897.solid 
query_dir=../rnaseq/hg19/SRX011545/

../src/peta clean -s 0.15 -k 15 -l $query_dir/SRR027878 $query_dir/SRR027878.fa

../src/peta pair -p 0 -k 27 -m 174 -d 300 -t 1 -o ../SRR027878_out/ $query_dir/SRR027878.fa $query_dir/SRR027878.solid 

