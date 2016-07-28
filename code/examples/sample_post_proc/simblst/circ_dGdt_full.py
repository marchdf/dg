#!/usr/bin/env python
#
# Using the .txt files, get the circulation and time derivatives of
# circulation in the full domain (do this for the Vazsonyi equation too)
# 
#
__author__ = 'marchdf'

#================================================================================
#
# Imports
#
#================================================================================
import sys, os
import numpy as np
import read2dpos_nodal as r2d
from scipy.integrate import simps


#================================================================================
#
# Functions
#
#================================================================================


#================================================================================
#
# Setup
#
#================================================================================
fnames = [filename for filename in os.listdir('.') if filename.startswith("rho") and filename.endswith(".txt")]
timesteps = range(len(fnames))


#================================================================================
#
# Initialize
#
#================================================================================
t                  = np.zeros(len(timesteps))
G                  = np.zeros(len(timesteps))
dGdt_advective     = np.zeros(len(timesteps))
dGdt_compressible  = np.zeros(len(timesteps))
dGdt_baroclinic    = np.zeros(len(timesteps))
VG                 = np.zeros(len(timesteps))
dVGdt_advective    = np.zeros(len(timesteps))
dVGdt_compressible = np.zeros(len(timesteps))
dVGdt_baroclinic   = np.zeros(len(timesteps))

#================================================================================
#
# Loop on time steps
#
#================================================================================
for k,timestep in enumerate(timesteps): 

    # fields to load
    rhoname = 'rho{0:010d}.txt'.format(timestep)
    uxname  = 'ux{0:010d}.txt'.format(timestep)
    uyname  = 'uy{0:010d}.txt'.format(timestep)
    pname   = 'p{0:010d}.txt'.format(timestep)

    # Load the data
    timestep,time,rho_df = r2d.read2dpos_nodal_qua(rhoname)
    timestep,time,ux_df  = r2d.read2dpos_nodal_qua(uxname)
    timestep,time,uy_df  = r2d.read2dpos_nodal_qua(uyname)
    timestep,time,p_df   = r2d.read2dpos_nodal_qua(pname)

    # Get centroid data
    rho_df = r2d.df_qua_cellcenter(rho_df)
    ux_df  = r2d.df_qua_cellcenter(ux_df)
    uy_df  = r2d.df_qua_cellcenter(uy_df)
    p_df   = r2d.df_qua_cellcenter(p_df)

    # Get a numpy meshgrid of the centroid data
    X,Y,RHO = r2d.df_qua_cellcenter_to_numpy(rho_df)
    X,Y,UX  = r2d.df_qua_cellcenter_to_numpy(ux_df)
    X,Y,UY  = r2d.df_qua_cellcenter_to_numpy(uy_df)
    X,Y,P   = r2d.df_qua_cellcenter_to_numpy(p_df)

    # Calculate the gradients
    dx = np.gradient(X)[1]
    dy = np.gradient(Y)[0]
    drhody, drhodx = np.gradient(RHO)
    dudy, dudx     = np.gradient(UX)
    dvdy, dvdx     = np.gradient(UY)
    dpdy, dpdx     = np.gradient(P)
    drhodx = drhodx/dx; drhody = drhody/dy
    dudx = dudx/dx;     dudy = dudy/dy
    dvdx = dvdx/dx;     dvdy = dvdy/dy
    dpdx = dpdx/dx;     dpdy = dpdy/dy

    # Calculate the vorticity (curl of velocity) and it's gradients
    w = dvdx - dudy
    dwdy, dwdx = np.gradient(w)
    dwdx = dwdx/dx;     dwdy = dwdy/dy

    #
    # Vorticity equation terms
    # 
    
    # Advective term: -(\mbf{u}\cdot \nabla) \omega
    advective = -(UX*dwdx + UY*dwdy)

    # Compressible term: -\omega (\nabla \cdot \mbf{u})
    compressible = -w*(dudx + dvdy)
    
    # Baroclinic term: \frac{\nabla \rho \times \nabla p}{\rho^2}
    baroclinic = (drhodx*dpdy - drhody*dpdx)/(RHO**2)

    #
    # Vazsonyi equation terms
    #
    w_rho = w/RHO
    dw_rhody, dw_rhodx = np.gradient(w_rho)
    dw_rhodx = dw_rhodx/dx; dw_rhody = dw_rhody/dy

    # Vazsonyi advective term: -(\mbf{u}\cdot \nabla) \frac{\omega}{\rho}
    Vadvective = -(UX*dw_rhodx + UY*dw_rhody) 

    # Vazsonyi compressible term: -\frac{\omega}{\rho} (\nabla \cdot \mbf{u})
    Vcompressible = -w_rho*(dudx + dvdy)
    
    # Vazsonyi baroclinic term: \frac{\nabla \rho \times \nabla p}{\rho^3}
    Vbaroclinic = (drhodx*dpdy - drhody*dpdx)/(RHO**3)
    
    #
    # Integrate all these terms to get what we need
    #

    # Get circulation in the full domain
    t[k]  = time
    G[k]  = simps(simps(w,X[0,:]),Y[:,0])
    VG[k] = simps(simps(w_rho,X[0,:]),Y[:,0]) 

    # Get the dG/dt for the various terms in the vorticity equation
    dGdt_advective[k]    = simps(simps(advective,X[0,:]),Y[:,0])
    dGdt_compressible[k] = simps(simps(compressible,X[0,:]),Y[:,0])
    dGdt_baroclinic[k]   = simps(simps(baroclinic,X[0,:]),Y[:,0])

    dVGdt_advective[k]    = simps(simps(Vadvective,X[0,:]),Y[:,0])
    dVGdt_compressible[k] = simps(simps(Vcompressible,X[0,:]),Y[:,0])
    dVGdt_baroclinic[k]   = simps(simps(Vbaroclinic,X[0,:]),Y[:,0])


#================================================================================
#
# Save data to a file 
#
#================================================================================
oname = 'circ_dGdt.dat'
np.savetxt(oname, np.vstack([t,
                             G,
                             dGdt_advective,
                             dGdt_compressible,
                             dGdt_baroclinic,
                             dGdt_advective+dGdt_compressible+dGdt_baroclinic,
                             VG,
                             dVGdt_advective,
                             dVGdt_compressible,
                             dVGdt_baroclinic,
                             dVGdt_advective+dVGdt_compressible+dVGdt_baroclinic]).transpose(),
           delimiter=',',
           header='time, G, dGdt_advective, dGdt_compressible, dGdt_baroclinic, dGdt_total, VG, dVGdt_advective, dVGdt_compressible, dVGdt_baroclinic, dVGdt_total')
