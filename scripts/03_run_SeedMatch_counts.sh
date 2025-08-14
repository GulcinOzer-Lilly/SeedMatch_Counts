#!/bin/bash


# create  log folder
mkdir -p ../log

# create output folder
OUT=../data/
mkdir -p $OUT



while read TR SEED OUTFILE
do
   echo $TR
   echo $SEED
   echo $OUTFILE

   qsub -v "TR=$TR,SEED=$SEED,OUT=$OUT/$OUTFILE" \
        -N "SeedMatch_counts-$OUTFILE" -o "../log/SeedMatch_counts-$OUTFILE"  SeedMatch_counts.sge

done <  input_SeedMatch.txt


