#!/bin/bash

while true
do

    ./main generate a.txt 35
    uniq a.txt > b.txt
    ./main convex b.txt


done
