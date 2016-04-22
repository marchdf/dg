#!/bin/bash
#
#
# Generate some square meshes useful for profiling the code. 
#
# This script uses sample_square_quad.geo and sample_square_quad.geo
#GMSH='/Applications/Gmsh.app/Contents/MacOS/gmsh'

# create the mesh directory
wrkdir=mesh

# different resolutions
res=(16 32 64 128 256 512 1024 2048 ) #(16 32 64 128 256 512 1024 2048);
#res=(2048)
# number of processors 
nprocs=(1 2 4 8 16); 

# Loop in resolutions
for ((i=0; i<${#res[@]}; i++)); do #C style loop with the array length
    
    # Make the .geo file with the right resolution
    sed -e 's|RESOLUTION|'${res[i]}'|g'<sample_square_qua.geo >'my_square_qua.geo';
    sed -e 's|RESOLUTION|'${res[i]}'|g'<sample_square_tri.geo >'my_square_tri.geo';

    # Loop on number of processors
    for ((k=0; k<${#nprocs[@]}; k++)); do #C style loop with the array length
        # Loop on p
	for ((p=1; p <= 5; p++)); do 

	    echo ''
	    echo 'Creating mesh with with p='$p' and N='${res[i]}' and '${nprocs[k]}' partitions';
	    echo ''
	    gmsh -2 -order $p -part ${nprocs[k]} my_square_qua.geo -o ./$wrkdir/square_qua_${res[i]}_${p}_np${nprocs[k]}.msh;
	    gmsh -2 -order $p -part ${nprocs[k]} my_square_tri.geo -o ./$wrkdir/square_tri_${res[i]}_${p}_np${nprocs[k]}.msh;

	    $GMSH -2 -order $p -part ${nprocs[k]} my_square_qua.geo -o ./$wrkdir/square_qua_${res[i]}_${p}_np${nprocs[k]}.msh;
	    $GMSH -2 -order $p -part ${nprocs[k]} my_square_tri.geo -o ./$wrkdir/square_tri_${res[i]}_${p}_np${nprocs[k]}.msh;

	done
    done
    rm -f my_square*.geo;    
done
