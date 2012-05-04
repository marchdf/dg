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
markertype = ['ks','kp','kh','kd','ko']
for i in range(0,len(pdir)):
    figure(i)
    clf()
    # error for rho
    loglog(dxs[i,:],errors[i,:,0],color='r',linewidth=2)
    loglog(dxs[i,:],errors[i,:,0],markertype[0],markerfacecolor='r',markersize=12)
    # error for u
    loglog(dxs[i,:],errors[i,:,1],color='g',linewidth=2)
    loglog(dxs[i,:],errors[i,:,1],markertype[1],markerfacecolor='g',markersize=12)
    # error for E
    loglog(dxs[i,:],errors[i,:,2],color='b',linewidth=2)
    loglog(dxs[i,:],errors[i,:,2],markertype[2],markerfacecolor='b',markersize=12)
    # error for phic
    loglog(dxs[i,:],errors[i,:,3],color='m',linewidth=2)
    loglog(dxs[i,:],errors[i,:,3],markertype[3],markerfacecolor='m',markersize=12)
    # error for phinc
    loglog(dxs[i,:],errors[i,:,4],color='c',linewidth=2)
    loglog(dxs[i,:],errors[i,:,4],markertype[4],markerfacecolor='c',markersize=12)
    # theoretical: p+1 slope
    loglog(dxs[i,:],err_th[i,:],color='k',linewidth=2,linestyle='dashed')
    xlabel('$\Delta x$',fontsize=22,fontweight='bold')
    ylabel('$L_2$ error',fontsize=22,fontweight='bold')
    setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
    setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
    setp(gca(),ylim=[1e-16,1e-1])
    savefig('p'+str(i)+'.png',format='png')
    savefig('p'+str(i)+'.eps',format='eps')


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

# make another latex table
for i in range(3,5): # loop on the different phi formulations
    print '\\begin{table}[!htb]'
    print '  \\footnotesize'
    print '  \\begin{center}'
    print '    \\begin{tabular}{cccp{0.1cm}ccp{0.1cm}cc}'
    print '      \\toprule'
    print '      & \\multicolumn{2}{c}{$p=',1,'$} & &  \\multicolumn{2}{c}{$p=',2,'$} & &  \\multicolumn{2}{c}{$p=',3,'$} \\\\'
    print '      \\cmidrule{2-3} \\cmidrule{5-6} \\cmidrule{8-9}'
    print '      $\Delta x$ & $L_2$ error &rate & & $L_2$ error &rate & & $L_2$ error &rate  \\\\'
    print '      $\\nicefrac{{1}}{{{0:2.0f}}}$ & {1:12.2e} & - & & {2:12.2e} & -  & & {2:12.2e} & - \\\\'.format(1/dxs[1,0], errors[1,0,i], errors[2,0,i], errors[3,0,i])
    for j in range(1,len(dxdir)):
        rate1 = log(errors[1,j,i]/errors[1,j-1,i])/log(dxs[1,j]/dxs[1,j-1])
        rate2 = log(errors[2,j,i]/errors[2,j-1,i])/log(dxs[1,j]/dxs[1,j-1])
        rate3 = log(errors[3,j,i]/errors[3,j-1,i])/log(dxs[1,j]/dxs[1,j-1])
        print '      $\\nicefrac{{1}}{{{0:2.0f}}}$ & {1:12.2e} & {2:12.2f} & & {3:12.2e} & {4:12.2f} & & {5:12.2e} & {6:12.2f}\\\\'.format(1/dxs[1,j], errors[1,j,i], rate1, errors[2,j,i], rate2, errors[3,j,i], rate3)
    print '      \\bottomrule'
    print '    \\end{tabular}'
    print '    \\caption{$L_2$ cell average error of $\phi_{nc}$ and convergence rates for different polynomials $p$. The conservative formulation is $2p+1$ super-convergent.}'
    print '    \\label{tab:conv_phinc}'
    print '  \\end{center}'
    print '\\end{table}'
        
