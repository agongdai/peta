query_dir=/home/carl/Projects/peta/rnaseq/Spombe/SRR097897
#../src/peta pair -p 0 -k 25 -m 197 -s 62 -t 4 -d 1 -b /home/carl/Projects/blat/blat -o ../SRR097897_out/ $query_dir/SRR097897.fa $query_dir/SRR097897.solid ../SRR097897_out/kmer.freq 
../src/peta ass -k 25 -m 197 -s 62 -t 4 -o ../SRR097897_part/ $query_dir/simu.fa $query_dir/SRR097897.solid
query_dir=/home/carl/Projects/peta/rnaseq/hg19/SRX011545/
#../src/peta ass -k 25 -m 197 -s 62 -t 4 -o ../SRR027876_out/ $query_dir/SRR027876.fa $query_dir/SRR027876.solid
#../src/peta pair -p 0 -k 33 -m 174 -s 60 -t 4 -d 1 -b /home/carl/Projects/blat/blat -o ../SRR027876_out/ $query_dir/SRR027876.fa $query_dir/SRR027876.solid ../SRR027876_out/kmer.freq 
#../src/peta kmer -k 25 -m 174 -s 60 -t 4 -o ../SRR027876_out/ $query_dir/SRR027876_corrected.fa $query_dir/SRR027876.solid ../SRR027876_out/kmer.cor.freq
