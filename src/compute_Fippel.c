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


  if(v_density[vALL] >= 0.9){
    v_result[vALL] = 1.0123 - 3.386e-5 * hadron->v_T[vALL]/UMeV + 0.291*(1+pow(hadron->v_T[vALL]/UMeV, -0.3421)) * (pow(v_density[vALL], -0.7) - 1);
  }

  // lung
  else if(v_density[vALL] > 0.0012 && v_density[vALL] <= 0.26){
    v_result[vALL] = ((0.9925 - 0.8815) / (0.26 - 0.0012))*(v_density[vALL] - 0.0012) + 0.8815;
  }

  
  else if(v_density[vALL] > 0.26 && v_density[vALL] < 0.9){
    v_result[vALL] = 1.0123 - 3.386e-5 * hadron->v_T[vALL]/UMeV + 0.291*(1+pow(hadron->v_T[vALL]/UMeV, -0.3421)) * 0.0765;
    v_result[vALL] = ((v_result[vALL] - 0.9925) / (0.9 - 0.26))*(v_density[vALL] - 0.26) + 0.9925;
  }

  // air
  else{
    v_result[vALL] = 0.8815;  // for rho <= 0.0012
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
  Fippel_Stop_Pow_correction(hadron, v_density, v_StpCorr);

  ALIGNED_(64) VAR_COMPUTE v_L[VLENGTH];
  Compute_L(hadron, v_N_el, v_density, material, Te_min, v_material_label, v_L);
  v_L[vALL] = v_L[vALL] * v_StpCorr[vALL];

  ALIGNED_(64) VAR_COMPUTE v_dE1[VLENGTH];
  v_dE1[vALL] = v_L[vALL] * v_s[vALL];

  ALIGNED_(64) VAR_COMPUTE v_tau1[VLENGTH];
  v_tau1[vALL] = hadron->v_T[vALL] / MC2_PRO;

  ALIGNED_(64) VAR_COMPUTE v_e1[VLENGTH];
  v_e1[vALL] = v_dE1[vALL] / hadron->v_T[vALL];

  ALIGNED_(64) VAR_COMPUTE v_C[VLENGTH];
  v_C[vALL] = v_L[vALL] * hadron->v_beta2[vALL];


  // Calcul numérique de la dérivée de C(E)

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);
  tmp.v_T[vALL] = hadron->v_T[vALL] * CONST_DERIV;
  Update_Hadron(&tmp);

  ALIGNED_(64) VAR_COMPUTE v_StpCorr2[VLENGTH];
  Fippel_Stop_Pow_correction(&tmp, v_density, v_StpCorr2);

  ALIGNED_(64) VAR_COMPUTE v_L2[VLENGTH];
  Compute_L(&tmp, v_N_el, v_density, material, Te_min, v_material_label, v_L2);
  v_L2[vALL] = v_L2[vALL] * v_StpCorr2[vALL];

  ALIGNED_(64) VAR_COMPUTE v_C2[VLENGTH];
  v_C2[vALL] = v_L2[vALL] * tmp.v_beta2[vALL];

  ALIGNED_(64) VAR_COMPUTE v_deriv_C[VLENGTH];
  v_deriv_C[vALL] = (v_C2[vALL] - v_C[vALL]) / (tmp.v_T[vALL] - hadron->v_T[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_b[VLENGTH];
  v_b[vALL] = hadron->v_T[vALL] * v_deriv_C[vALL] / v_C[vALL];


  // Calcul de dE2

  v_result[vALL] = v_dE1[vALL] * (	1 
					+ (v_e1[vALL] / ((1+v_tau1[vALL]) * (2+v_tau1[vALL]))) 
					+ (	v_e1[vALL]*v_e1[vALL] 
						* (2+2*v_tau1[vALL]+v_tau1[vALL]*v_tau1[vALL]) 
						/ ((1+v_tau1[vALL])*(1+v_tau1[vALL])*(2+v_tau1[vALL])*(2+v_tau1[vALL]))
					  ) 
					- (v_b[vALL] * v_e1[vALL] * (0.5 + 2*v_e1[vALL]/(3*(1+v_tau1[vALL])*(2+v_tau1[vALL])) + (1-v_b[vALL]) * v_e1[vALL]/6)) 
			 	);

  return;
}
