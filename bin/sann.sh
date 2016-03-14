#!/bin/bash

# Where the NCBI programs have been installed
NCBIDIR=$HOME/ncbi

# Where the sequence database (uniref90filt)
nr90=$NCBIDIR/uniref/uniref90filt

# Where the NNDB_HOME and SANN_HOME
#NNDB_HOME=$HOME/database/nndb
#SANN_HOME=$HOME/sann

# Get basename
basename=$1

fasta=$basename.fa
if [ ! -e $fasta ]; then
    echo "$fasta is not exist! check it."
    exit
fi

NCPUS=8

# run psi-blast with standard parameters
echo "Running PSI-BLAST with sequence" $basename "..."

$NCBIDIR/bin/blastpgp -i $basename.fa -o $basename.bla -j 3 -e 1.000000e-03 -h 1.000000e-03 -v 1000 -b 1000 -a 4 -d $nr90 -C $basename.chk -Q $basename.prf

# transfomation
mkprf=$SANN_HOME/bin/mkchk2.py
$mkprf $basename


echo "Running SANN prediction..."
$SANN_HOME/bin/sann -i $basename -np $NCPUS


