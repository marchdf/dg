#!/bin/bash
#
# THIS IS A SAMPLE SCRIPT. YOU NEED TO CHANGE THE DIRECTORIES FOR LOOP
# TO LOOP ON THE ACTUAL THINGS YOU ARE INTERESTED IN.
# 
# Converts the 1D .pos gmsh files into a .txt file
# to be read by a python script
# Requires pos2txt_sample.geo
#
# Usage: pos2txt.sh
#
GMSH='gmsh - pos2txt.geo'

field=('rho' 'ux' 'p' 'g' 'sensor');
for dir in ./shuoshe_p*; do
    for n in {0..4}; do # loop on fields
	file=$dir'/'${field[n]}'.pos';
	echo 'Converting the following files to .txt '$file;

	# get the order from the directory name
        order=$(echo $dir | grep -o -E '[0-9]+' | head -1)

	# output gmsh file
	sed -e 's|DIR|'$dir'|g' \
	    -e 's|ORDER|'$order'|g' \
	    -e 's|FIELD|'${field[n]}'|g' <pos2txt_sample.geo >'pos2txt.geo'; 

	# run gmsh command
	$GMSH
    done
done

rm pos2txt.geo