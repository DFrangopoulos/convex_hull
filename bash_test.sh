#!/bin/bash
#dfrangopoulos

gcc main.c -lm -DUNITTEST=true -o main

while true
do

    #testP1 
    echo "Testing P1"
    ./main generate a.txt 35
    uniq a.txt > b.txt
    ./main convex b.txt

    echo "Testing P2"
    ./main generate a.txt 35
    ./main generate c.txt 1
    uniq a.txt > b.txt
    arg=$(head -n1 c.txt)
    ./main check b.txt $arg

    echo "Testing P3"
    ./main generate a.txt 15
    uniq a.txt > b.txt
    ./main generate c.txt 13
    uniq c.txt > d.txt
    ./main intersect b.txt d.txt


done
