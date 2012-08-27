#rm *.o
#rm peta
make
echo ------------------------------------------------
./peta ass -n 2 -o 29 -r 67 -a 326 -m 326 -d 78 -s 1 -p 1 -b read/SRR097897.solid read/SRR097897.fa read/rnapet.fa
#./peta ass -n 2 -o 29 -r 45 -a 201 -m 201 -d 20 -s 0 -p 1 -b rnaseq/SRR027876.solid rnaseq/SRR027876.fasta read/rnapet.fa
