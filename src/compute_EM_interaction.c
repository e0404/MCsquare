/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_EM_interaction.h"

void Total_Stop_Pow(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_stop_pow){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_stop_pow, 64);

  int i;

  ALIGNED_(64) VAR_COMPUTE v_scaled_T[VLENGTH];
  v_scaled_T[vALL] = hadron->v_T[vALL]/hadron->v_mass[vALL];


  ALIGNED_(64) int v_index[VLENGTH];
//  v_index[vALL] = (int)(v_scaled_T[vALL] / (UMeV*PSTAR_BIN));
//  v_index[vALL] = (int)floor(v_scaled_T[vALL] / (UMeV*PSTAR_BIN));
  ALIGNED_(64) VAR_COMPUTE v_scaled_T2[VLENGTH];
  v_scaled_T2[vALL] = v_scaled_T[vALL] / (UMeV*PSTAR_BIN);
  v_index[vALL] = (int)floor(v_scaled_T2[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_data_Energy1[VLENGTH];
  v_data_Energy1[vALL] = v_index[vALL] * UMeV * PSTAR_BIN;

  ALIGNED_(64) VAR_COMPUTE v_data_Energy2[VLENGTH];
//  v_data_Energy2[vALL] = v_data_Energy1[vALL] + UMeV * PSTAR_BIN;
  v_data_Energy2[vALL] = (v_index[vALL]+1) * UMeV * PSTAR_BIN;


  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow1[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Stop_Pow1[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow2[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Stop_Pow2[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_index[i]+1];
  }

  vec_Linear_Interpolation(v_scaled_T, v_data_Energy1, v_data_Energy2, v_Stop_Pow1, v_Stop_Pow2, v_stop_pow);

  return;
}


void Total_Hard_Cross_Section(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, DATA_config *config, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_dE_max, 64);
  __assume_aligned(v_result, 64);

  ALIGNED_(64) VAR_COMPUTE v_cross_section[VLENGTH];
  cross_section_ionization(hadron, v_N_el, Te_min, v_cross_section);

  ALIGNED_(64) VAR_COMPUTE v_tmp_result[VLENGTH];
  if(config->Simulate_Nuclear_Interactions == 1){
    total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_tmp_result);
    v_cross_section[vALL] += v_tmp_result[vALL];
  }


  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);
  tmp.v_T[vALL] = hadron->v_T[vALL] - v_dE_max[vALL];
  if(tmp.v_T[vALL] <= 0) tmp.v_T[vALL] = hadron->v_T[vALL];
  Update_Hadron(&tmp);
    
  ALIGNED_(64) VAR_COMPUTE v_cross_section2[VLENGTH];
  cross_section_ionization(&tmp, v_N_el, Te_min, v_cross_section2);

  if(config->Simulate_Nuclear_Interactions == 1){
    total_Nuclear_cross_section(&tmp, material, v_material_label, v_density, v_tmp_result);
    v_cross_section2[vALL] += v_tmp_result[vALL];
  }

  v_result[vALL] = fmax(v_cross_section[vALL], v_cross_section2[vALL]);

  return;
}


void get_interaction_type(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, VAR_COMPUTE *v_tot_section, VSLStreamStatePtr RNG_Stream, DATA_config *config, int *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_dE_max, 64);
  __assume_aligned(v_tot_section, 64);
  __assume_aligned(v_result, 64);


  v_result[vALL] = 0;

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_ionization_section[VLENGTH];
  cross_section_ionization(hadron, v_N_el, Te_min, v_ionization_section);
  v_ionization_section[vALL] = v_ionization_section[vALL] / v_tot_section[vALL];

  if(v_rnd[vALL] <= v_ionization_section[vALL]){
    v_result[vALL] = 1;
  }

  if(config->Simulate_Nuclear_Interactions == 1){
    ALIGNED_(64) VAR_COMPUTE v_nuclear_section[VLENGTH];
    total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_nuclear_section);
    v_nuclear_section[vALL] = (v_nuclear_section[vALL] / v_tot_section[vALL]) + v_ionization_section[vALL];

    if(v_rnd[vALL] <= v_ionization_section[vALL]){
      v_result[vALL] = 1;
    }
    else if(v_rnd[vALL] <= v_nuclear_section[vALL]){
      v_result[vALL] = 2;
    }
  }

  else{
    if(v_rnd[vALL] <= v_ionization_section[vALL]){
      v_result[vALL] = 1;
    }
  }

  return;
}

void cross_section_ionization(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_result, 64);



  ALIGNED_(64) VAR_COMPUTE v_log_result[VLENGTH];
  v_log_result[vALL] = hadron->v_Te_max[vALL]/Te_min;
  v_log_result[vALL] = log(v_log_result[vALL]);

  v_result[vALL] =  	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] 
			* (	((1.0/Te_min) - (1.0/hadron->v_Te_max[vALL])) 
				- (hadron->v_beta2[vALL]/hadron->v_Te_max[vALL]) * v_log_result[vALL] 
				+ (hadron->v_Te_max[vALL]-Te_min) / (2*hadron->v_E[vALL]*hadron->v_E[vALL])
			  )
			/ (hadron->v_beta2[vALL]);

  if(hadron->v_Te_max[vALL] <= Te_min) v_result[vALL] = 0.0;

  return;
}




void Compute_L(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  ALIGNED_(64) VAR_COMPUTE v_log_result[VLENGTH];
  v_log_result[vALL] = hadron->v_Te_max[vALL]/Te_min;
  v_log_result[vALL] = log(v_log_result[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_M[VLENGTH];

  v_M[vALL] =	(2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] / hadron->v_beta2[vALL]) 
		* ( 	v_log_result[vALL]
			- (hadron->v_Te_max[vALL] - Te_min) * hadron->v_beta2[vALL] / hadron->v_Te_max[vALL]
			+ (hadron->v_Te_max[vALL]*hadron->v_Te_max[vALL] - Te_min*Te_min) / (4*hadron->v_E[vALL]*hadron->v_E[vALL])
		);

  if(hadron->v_Te_max[vALL] <= Te_min) v_M[vALL] = 0;

  Total_Stop_Pow(hadron, material, v_material_label, v_result);
  v_result[vALL] = v_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_result[vALL] - v_M[vALL];	// Pouvoir d'arrêt restreint en eV / cm

  return;
}


void Compute_dE2(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result){

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

  ALIGNED_(64) VAR_COMPUTE v_L[VLENGTH];
  Compute_L(hadron, v_N_el, v_density, material, Te_min, v_material_label, v_L);

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

  ALIGNED_(64) VAR_COMPUTE v_L2[VLENGTH];
  Compute_L(&tmp, v_N_el, v_density, material, Te_min, v_material_label, v_L2);

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


void Compute_Energy_straggling(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result){

  __assume_aligned(v_N_el, 64);
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

  v_result[vALL] = 	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_s[vALL] 
			* fmin(Te_min, hadron->v_Te_max[vALL]) 
			* (1 - 0.5*hadron->v_beta2[vALL]) / hadron->v_beta2[vALL];

  return;
}


void Compute_MS_Fippel(Hadron *hadron, VAR_COMPUTE *v_s, VAR_COMPUTE *v_X0, VAR_COMPUTE *v_result){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_X0, 64);
  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);
  
  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  v_result[vALL] = (CONST_MS_Fippel*UMeV * hadron->v_charge[vALL] / (hadron->v_beta2[vALL]*hadron->v_gamma[vALL]*MC2_PRO)) * sqrt(v_s[vALL]/v_X0[vALL]);

  return;
}


void Compute_Ionization_Energy(Hadron *hadron, VAR_COMPUTE Te_min, VSLStreamStatePtr RNG_Stream, VAR_COMPUTE *v_result){

  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_Te[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_g[VLENGTH];

  ALIGNED_(64) VAR_COMPUTE v_mask[VLENGTH];
  v_mask[vALL] = 1.0;

  if(hadron->v_Te_max[vALL] < Te_min){
    v_result[vALL] = 0.0;
    v_mask[vALL] = 0.0;
  }

  int run = __sec_reduce_add(v_mask[vALL]);

  while(run != 0.0){
    rand_uniform(RNG_Stream, v_rnd);
    v_Te[vALL] = ( Te_min * hadron->v_Te_max[vALL]) / ((1-v_rnd[vALL]) * hadron->v_Te_max[vALL] + v_rnd[vALL] * Te_min);
    v_g[vALL] = 1.0 - hadron->v_beta2[vALL] * (v_Te[vALL]/hadron->v_Te_max[vALL]) + v_Te[vALL]*v_Te[vALL]/(2*hadron->v_E[vALL]*hadron->v_E[vALL]);
    rand_uniform(RNG_Stream, v_rnd);

    if(v_rnd[vALL] <= v_g[vALL] && v_mask[vALL] == 1.0){
      v_result[vALL] = v_Te[vALL];
      v_mask[vALL] = 0.0;
    }

  run = __sec_reduce_add(v_mask[vALL]);
  }

  return;
}



