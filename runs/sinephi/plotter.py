#!/usr/bin/env python2
import sys
from numpy import *
import matplotlib.axis as axis
from pylab import*

rc('text', usetex=True)
rc('font', family='sans-serif')

#load the files
pdir  = ['p0/','p1/','p2/','p3/','p4/']
dxdir = ['dx0/','dx1/','dx2/','dx3/','dx4/']
dxFile  = 'error.dat'

dxs    = zeros((size(pdir),size(dxdir)))
errors = zeros((size(pdir),size(dxdir),5))

# Load all the files
pdcnt = 0
for pd in pdir:
    dxdcnt = 0
    for dxd in dxdir:
        dat    = loadtxt(pd+dxd+dxFile)
        print 'Loading', pd,dxd, 'data:', dat
        dxs[pdcnt,dxdcnt] = dat[0,0]
        errors[pdcnt,dxdcnt,0] = dat[0,1]
        errors[pdcnt,dxdcnt,1] = dat[0,2]
        errors[pdcnt,dxdcnt,2] = dat[0,3]
        errors[pdcnt,dxdcnt,3] = dat[0,4]
        errors[pdcnt,dxdcnt,4] = dat[0,5]
        dxdcnt = dxdcnt+1
    pdcnt = pdcnt+1

# Theoretical error
err_th = zeros((size(pdir),size(dxdir)))
x  = array([[1.,2.,3.],
            [4.,5.,6.]])
for i in range(0,len(pdir)):
    err_th[i,:] = dxs[i,:]**(2*i+1)*errors[i,0,0]/dxs[i,0]**(2*i+1)

# plot these
for i in range(0,len(pdir)):
    figure(i)
    clf()
    # error for rho
    loglog(dxs[i,:],errors[i,:,0],color='r',linewidth=2)
    loglog(dxs[i,:],errors[i,:,0],'ks',markerfacecolor='r',markersize=12)
    # error for u
    loglog(dxs[i,:],errors[i,:,1],color='g',linewidth=2)
    loglog(dxs[i,:],errors[i,:,1],'ks',markerfacecolor='g',markersize=12)
    # error for E
    loglog(dxs[i,:],errors[i,:,2],color='b',linewidth=2)
    loglog(dxs[i,:],errors[i,:,2],'ks',markerfacecolor='b',markersize=12)
    # error for phic
    loglog(dxs[i,:],errors[i,:,3],color='m',linewidth=2)
    loglog(dxs[i,:],errors[i,:,3],'ks',markerfacecolor='m',markersize=12)
    # error for phinc
    loglog(dxs[i,:],errors[i,:,4],color='c',linewidth=2)
    loglog(dxs[i,:],errors[i,:,4],'ks',markerfacecolor='c',markersize=12)
    # theoretical: p+1 slope
    loglog(dxs[i,:],err_th[i,:],color='k',linewidth=2,linestyle='dashed')
    xlabel('$\Delta x$',fontsize=22,fontweight='bold')
    ylabel('$E$',fontsize=22,fontweight='bold')
    setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
    setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
    savefig('p'+str(i)+'.png',format='png')
# setp(gca(),xlim=[0,10000])
# setp(gca(),ylim=[0,60])


# make a latex table
print '\\begin{table}[!htb]'
print '  \\footnotesize'
print '  \\begin{minipage}[b]{0.45\\textwidth}'
print '  \\begin{center}'
print '    \\begin{tabular}{ccc}'
for i in range(0,len(pdir)):
    print '      \\toprule'
    print '      \\multicolumn{3}{c}{$p=',i,'$} \\\\'
    print '      \\midrule'
    print '      $\\Delta x$ &  $||E||_2$       & rate\\\\'
    print '      {0:12.5f} & {1:12.3e} & -\\\\'.format(dxs[i,0], errors[i,0,3])
    for j in range(1,len(dxdir)):
        rate = log(errors[i,j,3]/errors[i,j-1,3])/log(dxs[i,j]/dxs[i,j-1])
        print '      {0:12.5f} & {1:12.3e} & {2:12.4f}\\\\'.format(dxs[i,j], errors[i,j,3], rate)
print '      \\bottomrule'
print '    \\end{tabular}'
print '    \\caption{Convergence for different $p$ for conservative $\phi$.}'
print '    \\label{tab:conv_phic}'
print '  \\end{center}'
print '  \\end{minipage}'
print '  \\hspace{1cm}'
print '  \\begin{minipage}[b]{0.45\\textwidth}'
print '  \\begin{center}'
print '    \\begin{tabular}{ccc}'
for i in range(0,len(pdir)):
    print '      \\toprule'
    print '      \\multicolumn{3}{c}{$p=',i,'$} \\\\'
    print '      \\midrule'
    print '      $\\Delta x$ &  $||E||_2$       & rate\\\\'
    print '      {0:12.5f} & {1:12.3e} & -\\\\'.format(dxs[i,0], errors[i,0,4])
    for j in range(1,len(dxdir)):
        rate = log(errors[i,j,4]/errors[i,j-1,4])/log(dxs[i,j]/dxs[i,j-1])
        print '      {0:12.5f} & {1:12.3e} & {2:12.4f}\\\\'.format(dxs[i,j], errors[i,j,4], rate)
print '      \\bottomrule'
print '    \\end{tabular}'
print '    \\caption{Convergence for different $p$ for non-conservative $\phi$.}'
print '    \\label{tab:conv_phinc}'
print '  \\end{center}'
print '  \\end{minipage}'
print '\\end{table}'
