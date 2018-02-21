#!/bin/bash

#echo "bash script looking for vis_SCHEME = 3 to change to vis_SCHEME = 4"
echo "Running mesh refinement study; code should be hit with run.py -f BEFORE running this script, with p=1 and Nx=A resolution."

#The initial run:
./run.py
cp HiOWvtxProject/ErrorPhil.csv ErrorPhil_p1_NxA.csv

#Change the mesh name in run.py:
sed -i -e 's/MESHFILE=MESHPFX+"1_p1.msh"/MESHFILE=MESHPFX+"2_p1.msh"/g' run.py
./run.py
cp HiOWvtxProject/ErrorPhil.csv ErrorPhil_p1_NxB.csv

#Change the mesh name in run.py:
sed -i -e 's/MESHFILE=MESHPFX+"2_p1.msh"/MESHFILE=MESHPFX+"3_p1.msh"/g' run.py
./run.py
cp HiOWvtxProject/ErrorPhil.csv ErrorPhil_p1_NxC.csv
