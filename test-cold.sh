#!/bin/sh
if [ $# -lt 1 ]; then cat <<- EOF
Usage:

  [ z=80000 .. ] $0 <distance>

This is a benchmarking harness for 'suggest' similar to test-suggest.sh that
takes as primary arguments the max distance and max query distance and creates
an estimated number of cold cache random disk accesses file.
EOF
    exit
fi
: ${freqs:="/tmp/freqs"}
: ${z:=80000}
: ${d:=1}       #Number of deletes in makeTypos -d {1..4}
: ${qm:=5}      #Matches query -m {1..10}

set -e; rm -rf /tmp/[0-9][0-9]*

dir=/tmp/$z
mkdir -p $dir/typos
cd $dir
head -n $z $freqs > freqs
suggest makeTypos -d$d -s10000 -n1 -p freqs
suggest update -pp -d$1 -i freqs
suggest query -pp -vd$1 -m$qm `cat typos.0`        | #-v ->instrumented
  grep '^ '                                        | #Non-summary lines
  awk '{ print 3*$7 + 2*(($9 < 188) ? $9 : 188) }' | #Estimate rand accesses
  sort -n > ../randAcc.d$1                           #Sort and store
