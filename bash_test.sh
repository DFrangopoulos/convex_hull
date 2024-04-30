#!/bin/bash

#ensure smooth exit
trap ctrl_c INT

function ctrl_c() {
        echo "Cleaning up..."
        rm  a.txt b.txt c.txt d.txt
        echo "Exiting..."
        exit 0
}

#check args
if [ -z "$1" ] || [ ! -f "$1" ]
then
  echo "Usage: $0 <path_to_test_build_binary>"
  exit 1
fi

#stress test
while true
do
    #testP1 
    echo "Testing P1"
    "$1" generate a.txt 35 || break;
    uniq a.txt > b.txt
    "$1" convex b.txt || break;

    echo "Testing P2"
    "$1" generate a.txt 35 || break;
    "$1" generate c.txt 1 || break;
    uniq a.txt > b.txt
    arg=$(head -n1 c.txt)
    "$1" check b.txt $arg || break;

    echo "Testing P3"
    "$1" generate a.txt 15 || break;
    uniq a.txt > b.txt
    "$1" generate c.txt 13 || break;
    uniq c.txt > d.txt
    "$1" intersect b.txt d.txt || break;

done
