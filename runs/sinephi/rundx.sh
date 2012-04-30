#!/bin/bash
#
# Run the dg code for all the p=0,1,... for all the dx in those
# directories
#
# Usage: rundx
#
DG='../../../../code/dg1d -d deck.inp'

# Clean the directory trees
rm -rf $(find p* -type d);

# Make the directories and input deck for 
# p=0 to 4 and dx=0 to 4
flux='llf';
limiter='none';
for i in {0..4}; do
    for j in {0..4}; do
	echo 'Creating directory p'$i'/dx'$j;
	mkdir -p 'p'$i'/dx'$j;
	sed -e 's|ORDER|'$i'|g' \
	    -e 's|DX|'$j'|g' \
	    -e 's|LIMITER|'$limiter'|g' \
	    -e 's|FLUX|'$flux'|g' <sample.inp >'p'$i'/dx'$j'/deck.inp'; 
    done
done


for DIR in $(find p*/dx*/ -type d); do
    echo '' ; 
    echo '------ Entering' $DIR '------'; 
    echo '' ; 
    cd $DIR;
    $DG
    cd ../..;
done



