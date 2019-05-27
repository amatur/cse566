./exclude/dsk-v2.3.0-bin-Darwin/bin/dsk  -file exclude/staph31/list_reads.unitigs.fa  -kmer-size 31 -abundance-min 1 
./exclude/dsk-v2.3.0-bin-Darwin/bin/dsk  -file stitchedUnitigs.txt  -kmer-size 31 -abundance-min 1 

./exclude/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file list_reads.unitigs.h5 -out output-bcalm.txt
./exclude/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file stitchedUnitigs.h5 -out output-my.txt


sort output-my.txt > a.txt
sort output-bcalm.txt > b.txt

diff a.txt b.txt
