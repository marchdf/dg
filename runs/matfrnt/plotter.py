#!/usr/bin/env python2
import sys
from numpy import *
import matplotlib.axis as axis
from pylab import*

rc('text', usetex=True)
rc('font', family='sans-serif')

#load the files
pdir  = ['p1','p2']
dxdir = ['dx0','dx1','dx2','dx3']#,'dx4']
fdir  = ['llf','ncf','roe']
mdir  = ['invgamma','gammamod']
dxFile  = 'error.dat'

dxs    = zeros((size(pdir),size(dxdir)))
errors = zeros((size(pdir),size(dxdir),size(fdir),size(mdir),4))

# Load all the files
pdcnt = 0
for pd in pdir:
    dxdcnt = 0
    for dxd in dxdir:
        fdcnt = 0
        for fd in fdir:
            mdcnt=0
            for md in mdir:
                dat    = loadtxt(pd+'/'+dxd+'/'+fd+'/'+md+'/'+dxFile)
                print 'Loading', pd,dxd,fd,md,'data:'
                print dat
                dxs   [pdcnt,dxdcnt]               = dat[4,0]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,0] = dat[4,1]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,1] = dat[4,2]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,2] = dat[4,3]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,3] = dat[4,4]
                mdcnt = mdcnt+1
            fdcnt = fdcnt+1
        dxdcnt = dxdcnt+1
    pdcnt = pdcnt+1

# Theoretical error
err_th = zeros((size(pdir),size(dxdir)))
for i in range(0,len(pdir)):
    err_th[i,:] = dxs[i,:]**(2*i+1)*errors[i,0,0,0,3]/dxs[i,0]**(2*i+1)

# plot these
eps = 1e-16
cmap = ['r','g','b','m','c']
markertype = ['s','*','o']
field = ['den','vel','pre','gam']
fill = ['w','']
cnt = 0
for i in range(0,len(pdir)): # loop on element orders
    for j in range(0,4):     # loop on error for rho, u, p, gamma
        figure(cnt)
        clf()
        for k in range(0,len(fdir),2):#len(fdir)): # loop on fluxes
            loglog(dxs[i,:]+eps,errors[i,:,k,0,j]+eps,color='r',linewidth=2)
            loglog(dxs[i,:]+eps,errors[i,:,k,0,j]+eps,markertype[k],markerfacecolor='w',markeredgecolor='r',markersize=12)
            loglog(dxs[i,:]+eps,errors[i,:,k,1,j]+eps,color='g',linewidth=2)
            loglog(dxs[i,:]+eps,errors[i,:,k,1,j]+eps,markertype[k],markerfacecolor='g',markeredgecolor='g',markersize=12)
        # theoretical: 2p+1 slope
        # loglog(dxs[i,:],err_th[i,:],color='k',linewidth=2,linestyle='dashed')
        xlabel('$\Delta x$',fontsize=22,fontweight='bold')
        ylabel('$E$',fontsize=22,fontweight='bold')
        setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
        setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
        setp(gca(),ylim=[1e-16,1e-1])
        savefig(pdir[i]+field[j]+'.png',format='png')
        savefig(pdir[i]+field[j]+'.eps',format='eps')
        cnt = cnt+1

# # make a latex table
# print '\\begin{table}[!htb]'
# print '  \\footnotesize'
# print '  \\begin{minipage}[b]{0.45\\textwidth}'
# print '  \\begin{center}'
# print '    \\begin{tabular}{ccc}'
# for i in range(0,len(pdir)):
#     print '      \\toprule'
#     print '      \\multicolumn{3}{c}{$p=',i,'$} \\\\'
#     print '      \\midrule'
#     print '      $\\Delta x$ &  $||E||_2$       & rate\\\\'
#     print '      {0:12.5f} & {1:12.3e} & -\\\\'.format(dxs[i,0], errors[i,0,3])
#     for j in range(1,len(dxdir)):
#         rate = log(errors[i,j,3]/errors[i,j-1,3])/log(dxs[i,j]/dxs[i,j-1])
#         print '      {0:12.5f} & {1:12.3e} & {2:12.4f}\\\\'.format(dxs[i,j], errors[i,j,3], rate)
# print '      \\bottomrule'
# print '    \\end{tabular}'
# print '    \\caption{Convergence for different $p$ for conservative $\phi$.}'
# print '    \\label{tab:conv_phic}'
# print '  \\end{center}'
# print '  \\end{minipage}'
# print '  \\hspace{1cm}'
# print '  \\begin{minipage}[b]{0.45\\textwidth}'
# print '  \\begin{center}'
# print '    \\begin{tabular}{ccc}'
# for i in range(0,len(pdir)):
#     print '      \\toprule'
#     print '      \\multicolumn{3}{c}{$p=',i,'$} \\\\'
#     print '      \\midrule'
#     print '      $\\Delta x$ &  $||E||_2$       & rate\\\\'
#     print '      {0:12.5f} & {1:12.3e} & -\\\\'.format(dxs[i,0], errors[i,0,4])
#     for j in range(1,len(dxdir)):
#         rate = log(errors[i,j,4]/errors[i,j-1,4])/log(dxs[i,j]/dxs[i,j-1])
#         print '      {0:12.5f} & {1:12.3e} & {2:12.4f}\\\\'.format(dxs[i,j], errors[i,j,4], rate)
# print '      \\bottomrule'
# print '    \\end{tabular}'
# print '    \\caption{Convergence for different $p$ for non-conservative $\phi$.}'
# print '    \\label{tab:conv_phinc}'
# print '  \\end{center}'
# print '  \\end{minipage}'
# print '\\end{table}'
