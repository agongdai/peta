ulimit -c unlimited
query_dir=/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/
#query_dir=/home/carl/Projects/peta/rnaseq/Spombe/genome
#query_dir=/home/carl/Projects/peta/rnaseq/hg19/

#../src/peta ass -k 29 -m 197 -s 62 -t 4 -o ../simu_out/ $query_dir/simu.pe.fa $query_dir/SRR097897.solid
#../src/peta ass -k 29 -m 197 -s 62 -t 4 -o ../simu_out/ $query_dir/simu.pe.fa $query_dir/SRR097897.solid

../src/peta ass -k 25 -m 197 -s 62 -t 4 -o ../SRR097897_part/ $query_dir/SRR097897.part.fa $query_dir/SRR097897.solid
