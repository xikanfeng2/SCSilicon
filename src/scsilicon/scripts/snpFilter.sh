#!/bin/bash

echo 

function usage() {
	echo -e "$0 <snpFile> <outFile>\n"
}

if [ $# -ne 2 ]; then
	echo "input arguments are required!"
	echo -n "Usage: " 
	usage
	exit 0
fi

snpFile=$1
outFile=$2

[ ! -e $snpFile ] && echo "Error: $snpFile does not exist!" && exit 0

cat $snpFile | grep -v "^#" | grep "single" | grep "exact" | awk '{print $5 "\t" $2 "\t" $3+1 "\t" $10 "\t" $7 "\t" $8}' > $outFile
snp_num=$(cat $outFile | wc -l)

echo -e "$snp_num SNPs were saved into file $outFile.\n"

exit 0
