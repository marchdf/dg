/*!
  \file init_cond.h  
  \brief Header file for initial condition functions
  \copyright Copyright (C) 2012-2015, Regents of the University of Michigan
  \license This project is released under the GNU Public License. See LICENSE.
  \author Marc T. Henry de Frahan <marchdf@umich.edu>, Computational Flow Physics Laboratory, University of Michigan
*/
#ifndef INIT_COND_H
#define INIT_COND_H
#include "fullMatrix.h"
#include "scalar_def.h"
#include "constants.h"
#include "misc.h"
#include "simpleMesh.h"
#include <vector>

// Define the different initial condition functions
void buildLRstates_multifluid(scalar rhoL, scalar uL, scalar EtL, scalar gammaL, scalar rhoR, scalar uR, scalar EtR, scalar gammaR, const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
scalar rtaylor_integrate_density(scalar A, scalar B, scalar rho01, scalar rho02, scalar yint, scalar H);
void init_dg_shallow(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_tranvtx_singlefluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs);
void init_dg_tranvtx_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs);
void init_dg_velpert_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs);
void init_dg_simplew_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sodtube_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sodmono_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_contact_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_rhotact_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_matfrnt_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_sinegam_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_expogam_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_shckint_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const double Tf);
void init_dg_shuoshe_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_multint_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_blast1d_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_simblst_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_shckrar_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_rarecon_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_sodcirc_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_rminstb_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_rmmulti_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_rtaylor_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_khdrake_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_khuramp_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_khinstb_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_khblast_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_khpertu_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_blastrm_multifluid(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_sinephi_passive(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sodmono_passive(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_stffrnt_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_stfshck_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_stfbubl_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_shckdrp_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_drpwall_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U);
void init_dg_jetcrss_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_prsrflw_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_injectr_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_bblwedg_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_cfplrun_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const simpleMesh &m, int typeElement, const std::vector<double> &ic_inputs = std::vector<double>());
void init_dg_rmawave_stiffened(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const fullMatrix<scalar> &XYZCen, fullMatrix<scalar> &U, const std::vector<double> &ic_inputs = std::vector<double>());
/* void init_dg_euler1D_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_euler2D_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_sodtube_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_explode_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_ovortex_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_mhdroto_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U); */
/* void init_dg_brio_wu_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
/* void init_dg_alfvenw_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
/* void init_dg_fastshk_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
/* void init_dg_slowshk_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
/* void init_dg_fastrar_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
/* void init_dg_slowrar_mhd(const int N_s, const int N_E, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U); */
#endif

