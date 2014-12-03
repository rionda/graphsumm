#!/bin/bash
set -x
time ./summ -p10000 -k 10 data/gen1.txt 2> /dev/null
time ./summ -p -k 2 data/gr1.txt 2> /dev/null
time ./summ -p -k 3 data/gr2.txt 2> /dev/null
time ./summ -p -k 4 data/gr3.txt 2> /dev/null
time ./summ -p -k 5 data/gr4.txt 2> /dev/null
time ./summ -p -a -k 500 data/ca-HepPh.txt  2> /dev/null
time ./summ -p -a -k 1200 -d 800 data/ca-AstroPh.txt 2> /dev/null
time ./summ -p -a -k 400 data/wiki-Vote.txt  2> /dev/null
time ./summ -a -n -k100 -d300 data/TheMarkerAnonymized.csv 2> /dev/null
time ./summ -p -a -k2000 data/iref/edges.txt 2> /dev/null
time ./summ -p200000 -a -k 3000 ~/data/email-Enron.txt 
