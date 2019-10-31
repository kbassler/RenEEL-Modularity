#!/bin/bash
awk '{if ($1<$2) print $2" "$1; else print $1" "$2}' $1 | sort -n -k 1,1 -k 2,2 | awk '{ if (NF==1) {A=$1; B=$2; print A" "B; deg[A]+=1; deg[B]+=1;} else { if ((A!=$1)||(B!=$2)) {A=$1;B=$2;print A" "B;deg[A]+=1;deg[B]+=1;} }} END{  for (i=1;i<=$1;i++) print deg[i] > "degree.txt" }' | awk '{if ($1!=$2) print $1" "$2}' > clean.txt



let edges=`awk 'END{print NR}' clean.txt`
echo edges=$edges
let nodes=`tail -n 1 clean.txt | cut -d ' ' -f 1`
echo nodes=$nodes


awk '{if (NF>0) {print $1} else {print 0} }' degree.txt > degree1.txt
rm degree.txt
mv degree1.txt degree.txt

echo $nodes $edges > info.txt
