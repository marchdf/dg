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

echo "Concluded p=1 refinement study, moving to p=2"

sed -i -e 's/PDG="1"/PDG="2"/g' run.py

cd mesh
sed -i -e 's/Nx = 128+1;/Nx = 16+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py -f
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxA.csv

cd mesh
sed -i -e 's/Nx = 16+1;/Nx = 21+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxB.csv

cd mesh
sed -i -e 's/Nx = 21+1;/Nx = 32+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxC.csv

cd mesh
sed -i -e 's/Nx = 32+1;/Nx = 43+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxD.csv

cd mesh
sed -i -e 's/Nx = 43+1;/Nx = 64+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxE.csv

cd mesh
sed -i -e 's/Nx = 64+1;/Nx = 85+1;/g' PERI_setup.geo
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2_NxF.csv

echo "Concluded p=2 refinement study, moving to p=3"

sed -i -e 's/PDG="2"/PDG="3"/g' run.py

cd mesh
sed -i -e 's/Nx = 85+1;/Nx = 12+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py -f
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxA.csv

cd mesh
sed -i -e 's/Nx = 12+1;/Nx = 16+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxB.csv

cd mesh
sed -i -e 's/Nx = 16+1;/Nx = 24+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxC.csv

cd mesh
sed -i -e 's/Nx = 24+1;/Nx = 32+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxD.csv

cd mesh
sed -i -e 's/Nx = 32+1;/Nx = 48+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxE.csv

cd mesh
sed -i -e 's/Nx = 48+1;/Nx = 64+1;/g' PERI_setup.geo
gmsh -1 -order 3 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p3_NxF.csv

echo "Concluded p=3 refinement study, moving to p=4"

sed -i -e 's/PDG="3"/PDG="4"/g' run.py

cd mesh
sed -i -e 's/Nx = 64+1;/Nx = 9+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py -f
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxA.csv

cd mesh
sed -i -e 's/Nx = 9+1;/Nx = 12+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxB.csv

cd mesh
sed -i -e 's/Nx = 12+1;/Nx = 20+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxC.csv

cd mesh
sed -i -e 's/Nx = 20+1;/Nx = 25+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxD.csv

cd mesh
sed -i -e 's/Nx = 25+1;/Nx = 39+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxE.csv

cd mesh
sed -i -e 's/Nx = 39+1;/Nx = 51+1;/g' PERI_setup.geo
gmsh -1 -order 4 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p4_NxF.csv

echo "Concluded p=4 refinement study, moving to p=5"

sed -i -e 's/PDG="4"/PDG="5"/g' run.py

cd mesh
sed -i -e 's/Nx = 51+1;/Nx = 8+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py -f
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxA.csv

cd mesh
sed -i -e 's/Nx = 8+1;/Nx = 11+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxB.csv

cd mesh
sed -i -e 's/Nx = 11+1;/Nx = 16+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxC.csv

cd mesh
sed -i -e 's/Nx = 16+1;/Nx = 21+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxD.csv

cd mesh
sed -i -e 's/Nx = 21+1;/Nx = 32+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxE.csv

cd mesh
sed -i -e 's/Nx = 32+1;/Nx = 43+1;/g' PERI_setup.geo
gmsh -1 -order 5 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py
cp sinphilProject/ErrorPhil.csv ErrorPhil_p5_NxF.csv

