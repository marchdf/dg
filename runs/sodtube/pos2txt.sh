#!/bin/bash
#
# Converts the .pos gmsh files into a .txt file
# to be read by a python script
# Requires pos2txt_sample.geo
#
# Usage: pos2txt.sh
#
GMSH='gmsh pos2txt.geo'

# Make the directories and input deck for 
# p=0 to 4 and dx=0 to 4
flux=('llf' 'ncf' 'roe');
model=('invgamma' 'gammamod');
field=('rho' 'ux' 'et' 'p' 'g');
for i in {1..2}; do # loop on dg order
    for j in 4; do #loop on delta x
	for k in {0..2}; do # loop on fluxes
	    for p in {0..1}; do  # loop on models
		for n in {0..4}; do # loop on fields
		    file='p'$i'/dx'$j'/'${flux[k]}'/'${model[p]}'/'${field[n]}'.pos';
		    echo 'Converting the following files to .txt '$file;
		    sed -e 's|ORDER|'$i'|g' \
			-e 's|DX|'$j'|g' \
			-e 's|MODEL|'${model[p]}'|g' \
			-e 's|FIELD|'${field[n]}'|g' \
			-e 's|FLUX|'${flux[k]}'|g' <pos2txt_sample.geo >'pos2txt.geo'; 
		    $GMSH
		done
	    done
	done
    done
done

rm pos2txt.geo;