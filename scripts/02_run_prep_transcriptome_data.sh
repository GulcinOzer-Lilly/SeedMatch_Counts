#!/bin/bash


# create  log folder
mkdir -p ../log

# create output folder
OUT=../data_prelim/
mkdir -p $OUT


#10090   mouse   Mus musculus

while read TX NAME TXNAME
do
   echo $TX
   echo $NAME
   echo $TXNAME

   qsub -v "TX=$TX,NAME=$NAME,OUT=$OUT" \
        -N "prep-tr-$NAME" -o "../log/prep-tr-$NAME"  prep_transcriptome_data.sge

done <  ../data/genomes.txt


