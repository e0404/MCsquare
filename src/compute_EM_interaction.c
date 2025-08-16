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

  ALIGNED_(64) VAR_COMPUTE v_scaled_T[VLENGTH];
  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_scaled_T2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_data_Energy1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_data_Energy2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_Stop_Pow2[VLENGTH];

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_scaled_T[v] = hadron->v_T[v]/hadron->v_mass[v];
    v_scaled_T2[v] = v_scaled_T[v] / (UMeV*PSTAR_BIN);
    v_index[v] = (int)floor(v_scaled_T2[v]);
    v_data_Energy1[v] = v_index[v] * UMeV * PSTAR_BIN;
    v_data_Energy2[v] = (v_index[v]+1) * UMeV * PSTAR_BIN;
  }
  
  #pragma omp simd
  for (v = 0; v < VLENGTH; v++) {
    v_Stop_Pow1[v] = (VAR_COMPUTE)material[v_material_label[v]].Stop_Pow[v_index[v]];
    v_Stop_Pow2[v] = (VAR_COMPUTE)material[v_material_label[v]].Stop_Pow[v_index[v]+1];
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

  int v; // simd counter

  ALIGNED_(64) VAR_COMPUTE v_tmp_result[VLENGTH];
  if(config->Simulate_Nuclear_Interactions == 1){
    total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_tmp_result);

    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      v_cross_section[v] += v_tmp_result[v];
    }
  }


  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    tmp.v_T[v] = hadron->v_T[v] - v_dE_max[v];
    if(tmp.v_T[v] <= 0) tmp.v_T[v] = hadron->v_T[v];
  }
  Update_Hadron(&tmp);
    
  ALIGNED_(64) VAR_COMPUTE v_cross_section2[VLENGTH];
  cross_section_ionization(&tmp, v_N_el, Te_min, v_cross_section2);

  if(config->Simulate_Nuclear_Interactions == 1){
    total_Nuclear_cross_section(&tmp, material, v_material_label, v_density, v_tmp_result);
    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      v_cross_section2[v] += v_tmp_result[v];
    }
  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] = fmax(v_cross_section[v], v_cross_section2[v]);
  }

  return;
}


void get_interaction_type(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, VAR_COMPUTE *v_tot_section, VAR_RND_SEED RNG_Stream, DATA_config *config, int *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_N_el, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_dE_max, 64);
  __assume_aligned(v_tot_section, 64);
  __assume_aligned(v_result, 64);



  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_ionization_section[VLENGTH];
  cross_section_ionization(hadron, v_N_el, Te_min, v_ionization_section);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] = 0;
    v_ionization_section[v] = v_ionization_section[v] / v_tot_section[v];
    if(v_rnd[v] <= v_ionization_section[v]) v_result[v] = 1;
  }

  if(config->Simulate_Nuclear_Interactions == 1){
    ALIGNED_(64) VAR_COMPUTE v_nuclear_section[VLENGTH];
    total_Nuclear_cross_section(hadron, material, v_material_label, v_density, v_nuclear_section);

    #pragma omp simd
    for(v = 0; v<VLENGTH; v++)
      v_nuclear_section[v] = (v_nuclear_section[v] / v_tot_section[v]) + v_ionization_section[v];
    
    #pragma omp simd
    for (v = 0; v < VLENGTH; v++) {
    if(v_rnd[v] <= v_ionization_section[v]) v_result[v] = 1;
      else if(v_rnd[v] <= v_nuclear_section[v]) v_result[v] = 2;
    }

  }

  else{
    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      if(v_rnd[v] <= v_ionization_section[v]) v_result[v] = 1;
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


  // simd loops separated due to weird MSVC bug
  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_log_result[v] = hadron->v_Te_max[v]/Te_min;
    v_log_result[v] = log(v_log_result[v]);
  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] =  	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[v] * hadron->v_charge[v]*hadron->v_charge[v] 
			* (	((1.0/Te_min) - (1.0/hadron->v_Te_max[v])) 
				- (hadron->v_beta2[v]/hadron->v_Te_max[v]) * v_log_result[v] 
				+ (hadron->v_Te_max[v]-Te_min) / (2*hadron->v_E[v]*hadron->v_E[v])
			  )
			/ (hadron->v_beta2[v]);
  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++)
    if(hadron->v_Te_max[v] <= Te_min) v_result[v] = 0.0;
  

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
  ALIGNED_(64) VAR_COMPUTE v_M[VLENGTH];
  
  Total_Stop_Pow(hadron, material, v_material_label, v_result);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_log_result[v] = hadron->v_Te_max[v]/Te_min;
    v_log_result[v] = log(v_log_result[v]);

    v_M[v] =	(2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[v] * hadron->v_charge[v]*hadron->v_charge[v] / hadron->v_beta2[v]) 
		* ( 	v_log_result[v]
			- (hadron->v_Te_max[v] - Te_min) * hadron->v_beta2[v] / hadron->v_Te_max[v]
			+ (hadron->v_Te_max[v]*hadron->v_Te_max[v] - Te_min*Te_min) / (4*hadron->v_E[v]*hadron->v_E[v])
		);

    // Moved to separate simd loops due to MSVC
    //if(hadron->v_Te_max[v] <= Te_min) v_M[v] = 0;
    //v_result[v] = v_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_result[v] - v_M[v];	// Pouvoir d'arrêt restreint en eV / cm
  }

  #pragma omp simd
  for (v = 0; v < VLENGTH; v++) 
    if (hadron->v_Te_max[v] <= Te_min) v_M[v] = 0;

  #pragma omp simd
  for (v = 0; v < VLENGTH; v++) 
      v_result[v] = v_density[v] * hadron->v_charge[v] * hadron->v_charge[v] * v_result[v] - v_M[v];


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
  ALIGNED_(64) VAR_COMPUTE v_dE1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_tau1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_e1[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_C[VLENGTH];

  Compute_L(hadron, v_N_el, v_density, material, Te_min, v_material_label, v_L);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
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

  ALIGNED_(64) VAR_COMPUTE v_L2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_C2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_deriv_C[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_b[VLENGTH];

  Compute_L(&tmp, v_N_el, v_density, material, Te_min, v_material_label, v_L2);

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
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

  int v;
  // simd loops separated due to MSVC
  #pragma omp simd
  for (v = 0; v < VLENGTH; v++)
    Te_min = fmin(Te_min, hadron->v_Te_max[v]);

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] = 	2*M_PI*R_ELEC*R_ELEC*MC2_ELEC * v_N_el[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_s[v] 
			* Te_min
			* (1 - 0.5*hadron->v_beta2[v]) / hadron->v_beta2[v];
  }

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

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] = (CONST_MS_Fippel*UMeV * hadron->v_charge[v] / (hadron->v_beta2[v]*hadron->v_gamma[v]*MC2_PRO)) * sqrt(v_s[v]/v_X0[v]);
  }

  return;
}


void Compute_Ionization_Energy(Hadron *hadron, VAR_COMPUTE Te_min, VAR_RND_SEED RNG_Stream, VAR_COMPUTE *v_result){

  __assume_aligned(v_result, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_Te[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_g[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_mask[VLENGTH];

  int run = 0;
  int v;

  #pragma omp simd reduction(+:run)
  for(v = 0; v<VLENGTH; v++){
    v_mask[v] = 1.0;
    if(hadron->v_Te_max[v] < Te_min){
      v_result[v] = 0.0;
      v_mask[v] = 0.0;
    }
    run += v_mask[v];
  }

  while(run != 0.0){
    rand_uniform(RNG_Stream, v_rnd);

    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      v_Te[v] = ( Te_min * hadron->v_Te_max[v]) / ((1-v_rnd[v]) * hadron->v_Te_max[v] + v_rnd[v] * Te_min);
      v_g[v] = 1.0 - hadron->v_beta2[v] * (v_Te[v]/hadron->v_Te_max[v]) + v_Te[v]*v_Te[v]/(2*hadron->v_E[v]*hadron->v_E[v]);
    }

    rand_uniform(RNG_Stream, v_rnd);

    run = 0;

    #pragma omp simd reduction(+:run)
    for(v = 0; v<VLENGTH; v++){
      if(v_rnd[v] <= v_g[v] && v_mask[v] == 1.0){
        v_result[v] = v_Te[v];
        v_mask[v] = 0.0;
      }
      run += v_mask[v];
    }
  }

  return;
}



