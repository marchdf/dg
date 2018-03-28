#!/usr/bin/env python
#
# Given an image, generate a mesh and for each node in a mesh, 
#     - get its (x,y) coordinate, 
#     - interpolate that to a jpeg pixel value and grab the RGB value
#     - output a file with 3 columns: node id, RBG values
# 
#
__author__ = 'marchdf'
def usage():
    print '\nUsage: {0:s} [options ...]\nOptions:\n -i, --image\timage filename\n -n, --nx\tnumber of cells in x direction (width of image)\n -h, --help\tshow this message and exit\n'.format(sys.argv[0]) 

#================================================================================
#
# Imports
#
#================================================================================
# Ignore deprecation warnings
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import os,sys,getopt
import numpy as np
import re
from PIL import Image
from scipy import interpolate
import subprocess as sp

#================================================================================
#
# Define some functions
#
#================================================================================
def get_num_nodes(fname):
    # Open file and return the number of nodes in the mesh
    with open(fname) as f:
        for line in f:
            if '$Nodes' in line:
                return int(next(f))

#================================================================================
#
# Default arguments
#
#================================================================================
# image file name
iname = 'image.jpg'

# how many cells in the x-direction
# approximate for triangles, exact for quads
num_cells = '16'

#================================================================================
#
# Parse arguments
#
#================================================================================
try:
    myopts, args = getopt.getopt(sys.argv[1:],"hn:i:",["help","nx","image"])
except getopt.GetoptError as e:
    print (str(e))
    usage()
    sys.exit(2)

for o, a in myopts:
    if o in ('-h', '--help'):
        usage()
        sys.exit()
    elif o in ('-n', '--nx'):
        num_cells = a
    elif o in ('-i', '--image'):
        iname = a


#================================================================================
#
# Setup
#
#================================================================================

# element order (p)
order = 1

# location of gmsh sample scripts
sdir = '/home/marchdf/dg/code/scripts/'

# sample script to mesh the jpeg (triangles or quads)
#sname = sdir+'jpeg_mesh_sample_tri.geo'
sname = sdir+'jpeg_mesh_sample_qua.geo'

# script to mesh the jpeg
gname = 'jpeg_mesh.geo'

# mesh file name
fname = 'jpeg_mesh_'+str(order)+'.msh'

# output file name
ofile = 'jpeg_data.dat'

#================================================================================
#
# Get pixel RGB values from image and some other info
#
#================================================================================
print 'Mesh the jpeg'

# Replace the generic IMAGENAME in the gmsh script by the correct one
# Do the same for the number of cells
# sort of taken from http://stackoverflow.com/questions/13089234/replacing-text-in-a-file-with-python
replacements = {'IMAGENAME':iname, 'NUM_CELLS':num_cells}
with open(sname, "r") as f: 
    lines = f.readlines()
with open(gname, "w") as f:
    for line in lines:
        for src, tgt in replacements.iteritems():
            line = line.replace(src, tgt)
        f.write(line)

# Launch gmsh to mesh the jpeg
cmd = 'gmsh -2 -order '+str(order)+' '+gname+' -o '+fname
proc = sp.Popen(cmd, shell=True, stderr=sp.PIPE)
out,err = proc.communicate()

#================================================================================
#
# Get pixel RGB values from image and some other info
#
#================================================================================
print 'Parse image'

# Get image data
im = Image.open(iname)
pix = im.load()
w,h = im.size 

# # OTHER WAY OF ACCOMPLISHING THE SAME THING. When interpolating,
# # instead of using pix[x,y], use fr(x,y),fg(x,y),fb(x,y). Probably
# # more accurate but also slower?

# # Create grid of (x,y) values for R, G, B matrices
# # to create interpolation functions
# x = range(w)
# y = range(h)
# RGB = np.asarray(im)
# fr = interpolate.interp2d(x, y, RGB[:,:,0], kind='cubic')
# fg = interpolate.interp2d(x, y, RGB[:,:,1], kind='cubic')
# fb = interpolate.interp2d(x, y, RGB[:,:,2], kind='cubic')


#================================================================================
#
# Link node ID with the RGB value in the image
#
#================================================================================
print 'Link node ID and RGB image value'
start = False
cnt = 0

# Initialize variables
num_nodes = get_num_nodes(fname)
nodes = np.zeros((num_nodes,4))

# loop through lines in mesh
with open(fname) as f:
    for line in f:


        # start parsing the nodes after reaching this
        if '$Nodes' in line:
            start = True
            next(f)

        # stop parsing the mesh when you reach the end of the nodes
        elif '$EndNodes' in line:
            break

        elif start:
            line = line.split(' ')
            id = int(line[0])
            x  = float(line[1])
            y  = float(line[2])

            # first method
            nodes[cnt,0] = id-1 # to start counting at 0 like the c++ code
            rgb = pix[x,h-y-1]
            nodes[cnt,1] = rgb[0]
            nodes[cnt,2] = rgb[1]
            nodes[cnt,3] = rgb[2]
            
            # # second method more accurate. Uncomment above interpolation functions
            # nodes[cnt,0] = id
            # nodes[cnt,1] = fr(x,y)[0]
            # nodes[cnt,2] = fg(x,y)[0]
            # nodes[cnt,3] = fb(x,y)[0]
            cnt = cnt+1
           

# Now output the nodes array to a file which will then be read by our
# C++ initial condition function. It will use the RGB values to decide
# what the different fields should be
np.savetxt(ofile,nodes,fmt='%.0f',header=str(num_nodes), comments='')
