#!/bin/bash
#
# Run the dg code for all the p=0,1,... for all the dx in those
# directories
#
# Usage: rundx
#
DG='../../../../../../code/dg1d -d deck.inp'
WORKDIR=$(pwd)

# Clean the directory trees
rm -rf $(find p* -type d);

# Make the directories and input deck for 
# p=0 to 4 and dx=0 to 4
flux=('llf' 'ncf' 'roe');
model=('invgamma' 'gammamod');
limiter='hrl';
for i in {0..2}; do
    for j in {3..3}; do
	for k in {0..0}; do
	    for p in {0..1}; do
		dir='p'$i'/dx'$j'/'${flux[k]}'/'${model[p]};
		echo 'Creating directory '$dir;
		mkdir -p $dir;
		sed -e 's|ORDER|'$i'|g' \
		    -e 's|DX|'$j'|g' \
		    -e 's|MODEL|'${model[p]}'|g' \
		    -e 's|LIMITER|'$limiter'|g' \
		    -e 's|FLUX|'${flux[k]}'|g' <sample.inp >$dir'/deck.inp'; 
	    done
	done
    done
done


for DIR in $(find p*/dx*/*/*/ -type d); do
    echo '' ; 
    echo '------ Entering' $DIR '------'; 
    echo '' ; 
    cd $DIR;
    $DG
    cd $WORKDIR;
done