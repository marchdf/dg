######################################################################
#
# Function to read a 2D .txt file from a .pos file generated by gmsh
# and do some manipulations based on this data
#
######################################################################

import numpy as np
import pandas as pd

#================================================================================
def read2dpos_nodal_qua(fname):
    # This works for 2D quads for any order and assuming all output is
    # at the same time step

    #
    # Read the first row of data
    #
    df = pd.read_csv(fname, delim_whitespace = True, header=None,nrows=1)
    
    #
    # Rename the columns appropriately (based on this first row data)
    #
    # each chunck of data is formed by a set of 8 values
    # ['time.step','time','element.type','element.num',x,y,z,f]
    ncols = len(df.columns)
    nchunks = int(ncols/8)
    new_names = ['time.step','time.value','element.type','element.num','x0','y0','z0','f0']
    for chunk in range(1,nchunks):
        new_names.extend(['-','-','-','-','x{0:d}'.format(chunk),'y{0:d}'.format(chunk),'z{0:d}'.format(chunk),'f{0:d}'.format(chunk)])
    df.columns = new_names
    timestep = df['time.step'].iloc[0]
    time     = df['time.value'].iloc[0]
    
    #
    # Only read in the necessary columns (i.e. the ones not containing
    # the string '-', '.' or 'z')
    #
    cols_to_keep = [col for col in df.columns if (('-' not in col) and ('z' not in col) and ('.' not in col)) ]

    #
    # Read in only the data we want to keep
    # 
    df = pd.read_csv(fname, delim_whitespace = True, header=None, names=new_names, usecols=cols_to_keep)
    
    # Return the dataframe
    return timestep,time,df


#================================================================================
def df_qua_cellcenter(df):
    # From a dataframe generated by read2dpos_nodal_qua, add columns
    # which will be [x,y,z,favg] the coordinates of the cell center
    # and function value at the cell center

    # 3*(order+1)^2 = ncols => order = sqrt(ncols/3) - 1
    ncols = len(df.columns)
    order = int(np.sqrt(ncols/3) - 1)

    # if you are first order, calculate the average (= cell center value)
    if order == 1:
        df['xc'] = df[['x0','x1','x2','x3']].mean(axis=1)
        df['yc'] = df[['y0','y1','y2','y3']].mean(axis=1)
        df['fc'] = df[['f0','f1','f2','f3']].mean(axis=1)

    # if you are second order, the cell center value is the middle value
    elif order == 2:
        df['xc'] = df['x8']
        df['yc'] = df['y8']
        df['fc'] = df['f8']

    else:
        print("You are trying to calculate the cell center value for an element order > 2.")
        print("This has not been implemented yet.")
        print("You have 2 options: ")
        print("\t- you can code an interpolation scheme to get the right cell center value")
        print("\t- you can code an integration scheme to get the right cell average value")        

    return df


#================================================================================
def df_qua_cellcenter_to_numpy(df):
    # This converts the Pandas dataframe cell centers to a nice numpy
    # meshgrid we can use to plot with matplotlib

    # Because we are dealing with (x,y) coordinates as floats,
    # sometimes nodes that should have the same (x,y) position do
    # not. This next procedure rounds to 12 decimal points. The hope
    # is that nodes that differ by less than that amount will now
    # actually share the same coordinate (and the pivot operation will
    # work)
    df['xc'] = df['xc'].round(8)
    df['yc'] = df['yc'].round(8)

    # Pivot the table to organize by the cell center coordinates
    df=df.pivot('yc', 'xc', 'fc')

    # Get the unique x and y coordinates as well as the function
    # evaluated at those points
    x = df.columns.values
    y = df.index.values
    Z = df.values

    # Turn these coordinates into a meshgrid
    X,Y=np.meshgrid(x, y)

    return X,Y,Z
    

# # If you feel like playing with these functions, here's some stuff you can do.
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import matplotlib.axis as axis
# pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])

# # Field name and time step to load
# fieldname = 'ux'
# timestep  = 8
# fname = fieldname+'{0:010d}_p2.txt'.format(timestep)

# # Load the data
# print('Load data')
# timestep,time,df = read2dpos_nodal_qua(fname)
# print('Data loaded')

# # Get centroid data
# df = df_qua_cellcenter(df)

# # Get a numpy meshgrid of the centroid data
# X,Y,Z = df_qua_cellcenter_to_numpy(df)

# # Plot
# plt.pcolormesh(X, Y, Z)
# plt.colorbar()
# plt.axis('equal')
# plt.show()
