/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_Fippel.h"

void Fippel_Stop_Pow_correction(Hadron *hadron, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_result, 64);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    if(v_density[v] >= 0.9) v_result[v] = 1.0123 - 3.386e-5 * hadron->v_T[v]/UMeV + 0.291*(1+pow(hadron->v_T[v]/UMeV, -0.3421)) * (pow(v_density[v], -0.7) - 1);
    else if(v_density[v] > 0.0012 && v_density[v] <= 0.26) v_result[v] = ((0.9925 - 0.8815) / (0.26 - 0.0012))*(v_density[v] - 0.0012) + 0.8815;  // lung
    else if(v_density[v] > 0.26 && v_density[v] < 0.9){
      v_result[v] = 1.0123 - 3.386e-5 * hadron->v_T[v]/UMeV + 0.291*(1+pow(hadron->v_T[v]/UMeV, -0.3421)) * 0.0765;
      v_result[v] = ((v_result[v] - 0.9925) / (0.9 - 0.26))*(v_density[v] - 0.26) + 0.9925;
    }
    else v_result[v] = 0.8815;  // for rho <= 0.0012 // air
  }

  return;
}



void Compute_dE2_Fippel(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_s, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  // Valeur précalculées

  ALIGNED_(64) VAR_COMPUTE v_StpCorr[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_L[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_tau1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_e1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_C[VLENGTH];


  Fippel_Stop_Pow_correction(hadron, v_density, v_StpCorr);
  Compute_L(hadron, v_N_el, v_density, material, Te_min, v_material_label, v_L);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_L[v] = v_L[v] * v_StpCorr[v];
    v_dE1[v] = v_L[v] * v_s[v];
    v_tau1[v] = hadron->v_T[v] / MC2_PRO;
    v_e1[v] = v_dE1[v] / hadron->v_T[v];
    v_C[v] = v_L[v] * hadron->v_beta2[v];
  }

  // Calcul numérique de la dérivée de C(E)

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    tmp.v_T[v] = hadron->v_T[v] * CONST_DERIV;
  }

  Update_Hadron(&tmp);

  ALIGNED_(64) VAR_COMPUTE v_StpCorr2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_L2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_C2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_deriv_C[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_b[VLENGTH];


  Fippel_Stop_Pow_correction(&tmp, v_density, v_StpCorr2);
  Compute_L(&tmp, v_N_el, v_density, material, Te_min, v_material_label, v_L2);

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_L2[v] = v_L2[v] * v_StpCorr2[v];
    v_C2[v] = v_L2[v] * tmp.v_beta2[v];
    v_deriv_C[v] = (v_C2[v] - v_C[v]) / (tmp.v_T[v] - hadron->v_T[v]);
    v_b[v] = hadron->v_T[v] * v_deriv_C[v] / v_C[v];

    // Calcul de dE2
    v_result[v] = v_dE1[v] * (	1 
					+ (v_e1[v] / ((1+v_tau1[v]) * (2+v_tau1[v]))) 
					+ (	v_e1[v]*v_e1[v] 
						* (2+2*v_tau1[v]+v_tau1[v]*v_tau1[v]) 
						/ ((1+v_tau1[v])*(1+v_tau1[v])*(2+v_tau1[v])*(2+v_tau1[v]))
					  ) 
					- (v_b[v] * v_e1[v] * (0.5 + 2*v_e1[v]/(3*(1+v_tau1[v])*(2+v_tau1[v])) + (1-v_b[v]) * v_e1[v]/6)) 
			 	);
  }

  return;
}


