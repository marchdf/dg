#!/usr/bin/env python2
import sys
from numpy import *
import matplotlib.axis as axis
from pylab import*

rc('text', usetex=True)
rc('font', family='sans-serif')

# Load the files
pdir  = ['p0','p1','p2','p3','p4']
dxdir = ['dx0','dx1','dx2','dx3','dx4']
fdir  = ['llf','ncf','roe']
mdir  = ['invgamma','gammamod']
fieldFile=['rho', 'ux', 'p', 'g']
fieldLatex = ['$\\rho$', '$u$', '$p$', '$\\gamma$']

pd  = pdir[1]
dxd = dxdir[4]
fd  = fdir[2]
field = fieldFile[3]
md = mdir[1]
T = 0.2
cmap = ['r','g','b','m','c']

for j in range(0,len(mdir)): # loop on the models 
    md  = mdir[j]
    for k in range(0,len(fieldFile)): # loop on fields
        field = fieldFile[k]
        print "Reading file: ", pd+'/'+dxd+'/'+fd+'/'+md+'/'+field+'.txt'
        dat = loadtxt(pd+'/'+dxd+'/'+fd+'/'+md+'/'+field+'.txt')
        nt  = dat[:,0]        # extract the time step list
        dt  = dat[:,1]        # extract time step list
        el  = dat[:,3]        # extract the element list
        N_t = nt.max(0)+1     # number of time steps + 1 to account for initial
        N_E = el.max(0)+1     # number of elements
        N_x = dat.shape[1]/8  # number of x in an element 8 columns per element
    
        ind = nonzero(dt==T)  # get the indexes at the time we want to plot
        
        x       = zeros(N_E*N_x)      
        y       = zeros(N_E*N_x)      
        indx    = range(4,dat.shape[1],8)       # index for x position
        indy    = range(7,dat.shape[1],8)       # index for y data
        indsort = argsort(dat[ind[0][0],indx])  # index for sorting x in ascending order
        refine  = 40                            # refine x resolution
        xinterp = zeros(N_E*refine)
        yinterp = zeros(N_E*refine)      
        
        for i in range(0,int(N_E)):
            x[i*N_x:(i*N_x+N_x)] = dat[ind[0][i],indx][indsort]
            y[i*N_x:(i*N_x+N_x)] = dat[ind[0][i],indy][indsort]
            xinterp[i*refine:(i*refine+refine)] = linspace(x[i*N_x], x[i*N_x+N_x-1],refine);
            yinterp[i*refine:(i*refine+refine)] = polyval(polyfit(x[i*N_x:(i*N_x+N_x)],y[i*N_x:(i*N_x+N_x)],N_x-1),xinterp[i*refine:(i*refine+refine)])
        figure(k)
        plot(x,y,'.',markerfacecolor=cmap[j],markeredgecolor=cmap[j])
    #plot(xinterp,yinterp,linewidth=2,color='b')
        xlabel('$x$',fontsize=22,fontweight='bold')
        ylabel(fieldLatex[k],fontsize=22,fontweight='bold')
        setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
        setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
        if k < 3 : # for rho, u, p
            setp(gca(),ylim=[0,1.1])
        else:  # for gamma
            setp(gca(),ylim=[1,2])
        if k==2:
            # this is an inset axes over the main axes
            a = axes([.52, .5, .35, .35])
            plot(x,y,'.',markerfacecolor=cmap[j],markeredgecolor=cmap[j]);
            plot(x,y,'-',color=cmap[j]);
            title('Zoom')
            setp(a,xlim=[0,0.4],ylim=[0.30,0.315], xticks=[], yticks=[])
        savefig(field+'_'+fd+'.png',format='png')
        savefig(field+'_'+fd+'.eps',format='eps')