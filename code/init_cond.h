#ifndef INIT_COND_H
#define INIT_COND_H
#include "fullMatrix.h"
#include <scalar_def.h>

// Define the different initial condition functions

void init_dg_shallow(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_simplew_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sodtube_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_contact_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_rhotact_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_matfrnt_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sinegam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_expogam_multifluid(const int N_s, const int N_E, const int N_F, const int D, const int model, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_sinephi_passive(const int N_s, const int N_E, const int N_F, const int D, scalar &gamma, const fullMatrix<scalar> &XYZNodes, fullMatrix<scalar> &U);
void init_dg_euler1D_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_euler2D_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_sodtube_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_explode_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_ovortex_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_mhdroto_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, const scalar gamma, fullMatrix<scalar> &U);
void init_dg_brio_wu_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
void init_dg_alfvenw_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
void init_dg_fastshk_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
void init_dg_slowshk_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
void init_dg_fastrar_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
void init_dg_slowrar_mhd(const int N_s, const int N_E, const int N_F, const int D, const fullMatrix<scalar> &XYZNodes, scalar &gamma, fullMatrix<scalar> &U);
#endif
