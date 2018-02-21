#!/bin/bash

#echo "bash script looking for vis_SCHEME = 3 to change to vis_SCHEME = 4"
echo "Running mesh refinement study; code should be hit with run.py -f BEFORE running this script, with p=1 and Nx=A resolution."

#The initial run:
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxA.csv

cd mesh
sed -i -e 's/Nx = 24+1;/Nx = 32+1;/g' PERI_setup.geo
gmsh -1 -order 1 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxB.csv

cd mesh
sed -i -e 's/Nx = 32+1;/Nx = 48+1;/g' PERI_setup.geo
gmsh -1 -order 1 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxC.csv

cd mesh
sed -i -e 's/Nx = 48+1;/Nx = 64+1;/g' PERI_setup.geo
gmsh -1 -order 1 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxD.csv

cd mesh
sed -i -e 's/Nx = 64+1;/Nx = 96+1;/g' PERI_setup.geo
gmsh -1 -order 1 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxE.csv

cd mesh
sed -i -e 's/Nx = 96+1;/Nx = 128+1;/g' PERI_setup.geo
gmsh -1 -order 1 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p1_NxF.csv

echo "Concluded p=1 refinement study."
