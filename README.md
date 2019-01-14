# SANN: 
A Protein solvent accessibility prediction based on nearest neighbor method

## Pre-requisite:
    - Python, numpy
    - Psi-blast (blastpgp 2.2.23 [Feb-03-2010])
    - uniref90 sequence databse

## Installation:

*Download SANN package and build
~~~
  $ tar xvzf sann.tar.gz
  $ cd sann/src
  $ make
  $ make install
~~~

*Download database
~~~
  nndb: http://lee.kias.re.kr/insilico/data/nndb.tar.gz
  seqdb: http://lee.kias.re.kr/insilico/data/seqdb.tar.gz
~~~

  
*Set environment variable in your shell (ex: .bashrc)
~~~
  export NNDB_HOME=$HOME/database/nndb
  export SANN_HOME=$HOME/sann
  export PATH=$PATH:$SANN_HOME/bin
~~~

*Edit NCBIDIR and nr90 variable in the script “sann.sh”

## Update History:

- Current       2016-02-01
- First release 2015-01-21

## Run a example

input file: test.fa (a protein sequence in fasta format)
~~~
   $ cd example
   $ sann.sh test
~~~

## References

[1] Keehyoung Joo, Sung Jong Lee and Jooyoung Lee, Sann: Solvent accessibility prediction of proteins by nearest neighbor method, PROTEINS-STRUCTURE FUNCTION AND BIOINFORMATICS, 80, 1791-1797 (2012)
