#rm *.o
#rm peta
make
echo ------------------------------------------------
./peta ass -n 2 -o 29 -r 67 -a 300 -m 300 -d 60 -s 1 -p 1 read/SRR097897.fa read/rnapet.fa
