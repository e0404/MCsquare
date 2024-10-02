/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_semi_infinite_slab.h"

void SemiInfiniteSlab_step(Hadron *hadron, Materials *material, Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Hadron_ID, int *Nbr_hadrons, VAR_COMPUTE *RS_exit_position, DATA_config *config, machine_parameters *machine, VAR_RND_SEED RNG_Stream){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64);

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(&hadron->v_type, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  Update_Hadron(hadron);


  int i,j,r;

  // Compute physical quantities
  ALIGNED_(64) int v_material_label[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_init_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_N_el[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    if(hadron->v_type[i] == Unknown){
      v_material_label[i] = 0;
      v_init_density[i] = 1;
      v_N_el[i] = 1;
      v_X0[i] = 1;
    }
    else{
      r = field_data[Hadron_ID[i]]->RS_num;
      v_material_label[i] = machine->RS_Material[r];
      v_init_density[i] = machine->RS_Density[r];
      v_N_el[i] = material[machine->RS_Material[r]].N_el * v_init_density[i];
      v_X0[i] = material[machine->RS_Material[r]].X0 / v_init_density[i];
    }
  }

  // Compute total cross section
  ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step_max[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_max[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_section[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_mean_dE[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_straggling[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_MS[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  ALIGNED_(64) int v_interaction_type[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_hard[VLENGTH];


  Update_Hadron(hadron);

  // Compute physical quantities
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_material_label[v] = machine->RS_Material;
    v_init_density[v] = machine->RS_Density;
    v_N_el[v] = material[machine->RS_Material].N_el * v_init_density[v];
    v_Dist_Interface[v] = hadron->v_z[v] - RS_exit_position[v] + 1e-4;
    if(v_Dist_Interface[v] < 0) v_Dist_Interface[v] = 0;
  }

  // Compute total cross section
  Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
  
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_stop_pow[v] = v_init_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_stop_pow[v];
    v_step_max[v] = fmin(fmin(v_Dist_Interface[v], config->D_Max), (config->Epsilon_Max * hadron->v_T[v] / v_stop_pow[v]));
    v_dE_max[v] = v_step_max[v] * v_stop_pow[v];
  }

  Total_Hard_Cross_Section(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, config, v_section);
 
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_section[v] += 1e-10;
    v_section[v] *= 1.017;
  }

  // Compute step length
  rand_uniform(RNG_Stream, v_rnd);

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_step[v] = -log(v_rnd[v])/v_section[v];
    if(v_step[v] > v_step_max[v]) v_step[v] = v_step_max[v];  // stop at step_max
  }

  // Compute CSDA + MS
  Compute_dE2(hadron, v_N_el, v_init_density, material, (config->Te_Min*UMeV), v_material_label, v_step, v_mean_dE);  
  Compute_Energy_straggling(hadron, v_N_el, (config->Te_Min*UMeV), v_step, v_straggling);

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_straggling[v] = sqrt(v_straggling[v]);
    v_X0[v] = material[machine->RS_Material].X0 / v_init_density[v];
  }

  rand_normal(RNG_Stream, v_dE, v_mean_dE, v_straggling);
  Compute_MS_Fippel(hadron, v_step, v_X0, v_MS);
  rand_normal_zero(RNG_Stream, v_theta, v_MS);
  rand_uniform(RNG_Stream, v_phi);

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_phi[v] = 2*M_PI*v_phi[v];
  
    // Lose energy
    if(hadron->v_type[v] == Unknown) v_dE[v] = 0;
      hadron->v_T[v] = hadron->v_T[v] - v_dE[v];
    if(hadron->v_type[vALL] != Unknown && hadron->v_T[v] <= (config->Ecut_Pro * UMeV)){
      hadron->v_type[v] = Unknown;
      hadron_list[Hadron_ID[v]].type = Unknown;
    }
  }

  Update_Hadron(hadron);

  // Update position and direction
  Update_position(hadron, v_step);
  Update_direction(hadron, v_theta, v_phi);


  // Hard interaction


  get_interaction_type(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, v_section, RNG_Stream, config, v_interaction_type);
  Compute_Ionization_Energy(hadron, (config->Te_Min*UMeV), RNG_Stream, v_dE_hard);

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    if(v_step[v] == v_step_max[v]) v_interaction_type[v] = 0; // force ficitious interaction if step >= step_max
    if(hadron->v_type[v] == Unknown) v_dE_hard[v] = 0;  // Ionization
    if(v_interaction_type[v] == 1) hadron->v_T[v] = hadron->v_T[v] - v_dE_hard[v];
  }

  // Nuclear interaction
  DATA_Scoring tmp;
  int previous_Nbr_hadrons;

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    if(hadron->v_type[v] != Unknown && v_interaction_type[v] == 2){
      previous_Nbr_hadrons = *Nbr_hadrons;
      Compute_Nuclear_interaction(i, hadron, material, v_material_label[v], hadron_list, Nbr_hadrons, &tmp, RNG_Stream, config);
      if(hadron->v_type[v] == Unknown) hadron_list[Hadron_ID[i]].type = Unknown;
      for(j=previous_Nbr_hadrons; j<*Nbr_hadrons; j++){
         layer_data[j] = layer_data[Hadron_ID[v]];
         field_data[j] = field_data[Hadron_ID[v]];
      }
    }

    if(hadron->v_T[v] <= (config->Ecut_Pro * UMeV)){
      hadron->v_type[v] = Unknown;
      hadron_list[Hadron_ID[v]].type = Unknown;
    }
  }

}
