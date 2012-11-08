#rm *.o
#rm peta
make
echo ------------------------------------------------
query_dir=rnaseq/Spombe/SRR097897
./peta ass -n 2 -o 29 -r 67 -a 120 -m 120 -d 48 -s 0 -p 1 -o ../SRR097897_out -b $query_dir/SRR097897.solid $query_dir/SRR097897.fa 
#./peta ass -n 2 -o 29 -r 45 -a 201 -m 201 -d 20 -s 0 -p 1 -b rnaseq/SRR027876.solid rnaseq/SRR027876.fasta read/rnapet.fa
