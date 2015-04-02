######################################################################
#
# Function to read a 1D .pos file generated by gmsh
# and integrate it over the domain as a function of time
#
######################################################################

from numpy import *
from gaussxw import *

def integratepos(fname):

    dat = loadtxt(fname)
    nt  = dat[:,0]        # extract the time step list (number)
    dt  = dat[:,1]        # extract time step list (value)
    el  = dat[:,3]        # extract the element list
    N_t = nt.max(0)+1     # number of time steps + 1 to account for initial
    N_E = el.max(0)+1     # number of elements
    N_x = dat.shape[1]/8  # number of x in an element 8 columns per element
    order = N_x-1
 
    I = zeros((N_t,N_E))

    x       = zeros(N_x)
    y       = zeros(N_x)
    w       = zeros(N_x)
    wtmp    = zeros(N_x)
    xgauss  = zeros(N_x)
    ygauss  = zeros(N_x)
    indx    = range(4,dat.shape[1],8)       # index for x position
    indy    = range(7,dat.shape[1],8)       # index for y data
    time = dt[ ::N_E];

    for cnt,T in enumerate(time):
        ind = nonzero(dt==T)                    # get the indexes at the time we want to plot
        indsort = argsort(dat[ind[0][0],indx])  # index for sorting x in ascending order
        print cnt, T

        for i in range(0,int(N_E)):
            x = dat[ind[0][i],indx][indsort]
            y = dat[ind[0][i],indy][indsort]
            xtmp,wtmp = gaussxwab(N_x,x[0],x[-1])
            xgauss = xtmp[ ::-1] # fliplr
            w = wtmp[ ::-1]
            ygauss = polyval(polyfit(x,y,order),xgauss)
            
            for k in range(0,N_x):
                I[cnt,i] = I[cnt,i]+w[k]*ygauss[k]

    return I,time
