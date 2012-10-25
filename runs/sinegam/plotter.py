#!/usr/bin/env python2
import sys
from numpy import *
import matplotlib.axis as axis
from pylab import*

rc('text', usetex=True)
#rc('font', family='sans-serif')
rc('font', family='serif', serif='Times')

#load the files
#pdir  = ['p0','p1','p2','p3','p4']
pdir  = ['p0','p1','p2','p3','p4']
dxdir = ['dx0','dx1','dx2','dx3','dx4','dx5','dx6']
fdir  = ['roe']#,'ncf','roe']
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
                dxs   [pdcnt,dxdcnt]               = dat[2,0]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,0] = dat[2,1]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,1] = dat[2,2]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,2] = dat[2,3]
                errors[pdcnt,dxdcnt,fdcnt,mdcnt,3] = dat[2,4]
                mdcnt = mdcnt+1
            fdcnt = fdcnt+1
        dxdcnt = dxdcnt+1
    pdcnt = pdcnt+1

# Theoretical error
err_th = zeros((size(pdir),size(dxdir)))
for i in range(0,len(pdir)):
    order = int(pdir[i][-1])
    err_th[i,:] = dxs[i,:]**(2*order+1)*errors[i,0,0,0,3]/dxs[i,0]**(2*order+1)

# plot these
eps = 1e-16
cmap = ['r','g','b','m','c']
markertype = ['s','d','o']
field = ['den','vel','pre','gam']
fill = ['w','']
cnt = 0
for i in range(0,len(pdir)): # loop on element orders
    if i==0:
        pltdx = range(1,len(dxdir))
    else:
        pltdx = range(0,len(dxdir))
    for j in range(0,4):     # loop on error for rho, u, p, gamma
        figure(cnt)
        clf()
        for k in range(0,1):#len(fdir)): # loop on fluxes
            loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,color='r',linewidth=2)
            loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,markertype[k],markerfacecolor='w',markeredgecolor='r',markersize=12)
            loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,color='g',linewidth=2)
            loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,markertype[k],markerfacecolor='g',markeredgecolor='g',markersize=12)
        # theoretical: 2p+1 slope
        loglog(dxs[i,pltdx],err_th[i,pltdx],color='k',linewidth=2,linestyle='dashed')
        xlabel(r"$\Delta x$",fontsize=22)
        ylabel(r"$L_\infty$ \textit{error}",fontsize=22)
        setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
        setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
        setp(gca(),xlim=[0.5*dxs[i,:].min(),2*dxs[i,:].max()],ylim=[1e-16,1e-0])
        savefig(pdir[i]+field[j]+'.png',format='png')
        savefig(pdir[i]+field[j]+'.pdf',format='pdf')
        cnt = cnt+1

cnt = 0
for i in range(0,len(pdir)): # loop on element orders
    if i==0:
        pltdx = range(1,len(dxdir))
    else:
        pltdx = range(0,len(dxdir))
    figure(i)
    clf()
    k = 0 # flux
    # pressure
    j = 2 # field
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,color='r',linewidth=2)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,markertype[cnt],markerfacecolor='w',markeredgecolor='r',markersize=12)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,color='g',linewidth=2)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,markertype[cnt],markerfacecolor='g',markeredgecolor='g',markersize=12)
    cnt = cnt + 1
    # gamma
    j = 3 # field
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,color='r',linewidth=2)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,0,j]+eps,markertype[cnt],markerfacecolor='w',markeredgecolor='r',markersize=12)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,color='g',linewidth=2)
    loglog(dxs[i,pltdx]+eps,errors[i,pltdx,k,1,j]+eps,markertype[cnt],markerfacecolor='g',markeredgecolor='g',markersize=12)
    cnt = 0
    # theoretical: 2p+1 slope
    loglog(dxs[i,pltdx],err_th[i,pltdx],color='k',linewidth=2,linestyle='dashed')
    xlabel(r"$\Delta x$",fontsize=22)
    ylabel(r"$L_\infty$ \textit{error}",fontsize=22)
    setp(gca().get_ymajorticklabels(),fontsize=18,fontweight='bold');
    setp(gca().get_xmajorticklabels(),fontsize=18,fontweight='bold');    
    setp(gca(),xlim=[0.5*dxs[i,:].min(),2*dxs[i,:].max()],ylim=[1e-16,1e-0])
    savefig(pdir[i]+fdir[k]+'_p_g.png',format='png')
    savefig(pdir[i]+fdir[k]+'_p_g.pdf',format='pdf')


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
