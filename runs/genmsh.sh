#!/bin/bash 
#
# Generate the meshes for varying spatial resolutions
#
#

# Clean the directory trees                                                                                                                                                            
rm -rf $(find msh* -type d);    

#loop on the spatial resolutions
#num=(11 21 41 81 101); # create an array
num=(9 17 33 65 129 257 513); # create an array
for ((i=0; i<${#num[@]}; i++)); do #C style loop with the array length
    echo 'Creating temporary file line_dx'$i'.geo with '${num[i]}' intervals';
    sed 's|NUM|'${num[i]}'|g' <sample_line.geo >'line_dx'$i'.geo';   
done

# loop on the DG orders (p=0..4)
# and create the meshes
for ((i=0; i<= 4; i++)); do
    echo 'Creating directory msh'$i;
    mkdir -p 'msh'$i;
    for ((j=0; j<${#num[@]}; j++)); do
	echo '   meshing for dx'$j;
	gmsh -1 -order $i -o ./msh$i/my_line_$i\_dx$j.msh line_dx$j.geo;
    done    
done

# Clean up the temporary files created
rm -rf $(find line_dx*.geo);
