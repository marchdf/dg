#!/bin/bash


echo "Concluded p=1 refinement study, moving to p=2"

sed -i -e 's/PDG="1"/PDG="2"/g' run.py
#Yes, this works. Need to come back and make it work p1 to p5
cd mesh
gmsh -1 -order 2 -o sample_mesh_1.msh PERI_setup.geo
cd ..
./run.py -f
cp sinphilProject/ErrorPhil.csv ErrorPhil_p2.csv

