#! /usr/bin/env bash

rm -f trackDb.txt 

bash bigwig/makeBigWigDb.sh > bigwig/trackDb.txt

dbfiles=$(find . -maxdepth 2 -name 'trackDb.txt')
cat $dbfiles >> trackDb.txt
