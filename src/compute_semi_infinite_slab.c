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

void SemiInfiniteSlab_step(Hadron *hadron, Materials *material, Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Hadron_ID, int *Nbr_hadrons, VAR_COMPUTE *RS_exit_position, DATA_config *config, machine_parameters *machine, VSLStreamStatePtr RNG_Stream){

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

  // Compute physical quantities
  ALIGNED_(64) int v_material_label[VLENGTH];
  v_material_label[vALL] = machine->RS_Material;

  ALIGNED_(64) VAR_COMPUTE v_init_density[VLENGTH];
  v_init_density[vALL] = machine->RS_Density;

  ALIGNED_(64) VAR_COMPUTE v_N_el[VLENGTH];
  v_N_el[vALL] = material[machine->RS_Material].N_el * v_init_density[vALL];


  // Compute total cross section
  ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
  v_Dist_Interface[vALL] = hadron->v_z[vALL] - RS_exit_position[vALL] + 1e-4;
  if(v_Dist_Interface[vALL] < 0) v_Dist_Interface[vALL] = 0;

  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
  v_stop_pow[vALL] = v_init_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_stop_pow[vALL];

  ALIGNED_(64) VAR_COMPUTE v_step_max[VLENGTH];
  v_step_max[vALL] = fmin(fmin(v_Dist_Interface[vALL], config->D_Max), (config->Epsilon_Max * hadron->v_T[vALL] / v_stop_pow[vALL]));

  ALIGNED_(64) VAR_COMPUTE v_dE_max[VLENGTH];
  v_dE_max[vALL] = v_step_max[vALL] * v_stop_pow[vALL];

  ALIGNED_(64) VAR_COMPUTE v_section[VLENGTH];
  Total_Hard_Cross_Section(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, config, v_section); 
  v_section[vALL] += 1e-10;
  v_section[vALL] *= 1.017;


  // Compute step length
  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  v_step[vALL] = -log(v_rnd[vALL])/v_section[vALL];
  if(v_step[vALL] > v_step_max[vALL]) v_step[vALL] = v_step_max[vALL];  // stop at step_max


  // Compute CSDA + MS
  ALIGNED_(64) VAR_COMPUTE v_mean_dE[VLENGTH];
  Compute_dE2(hadron, v_N_el, v_init_density, material, (config->Te_Min*UMeV), v_material_label, v_step, v_mean_dE);  

  ALIGNED_(64) VAR_COMPUTE v_straggling[VLENGTH];
  Compute_Energy_straggling(hadron, v_N_el, (config->Te_Min*UMeV), v_step, v_straggling);
  v_straggling[vALL] = sqrt(v_straggling[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_dE[VLENGTH];
  rand_normal(RNG_Stream, v_dE, v_mean_dE, v_straggling);

  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
  v_X0[vALL] = material[machine->RS_Material].X0 / v_init_density[vALL];

  ALIGNED_(64) VAR_COMPUTE v_MS[VLENGTH];
  Compute_MS_Fippel(hadron, v_step, v_X0, v_MS);

  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  rand_normal_zero(RNG_Stream, v_theta, v_MS);

  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  rand_uniform(RNG_Stream, v_phi);
  v_phi[vALL] = 2*M_PI*v_phi[vALL];


  // Lose energy
  if(hadron->v_type[vALL] == Unknown) v_dE[vALL] = 0;
  hadron->v_T[vALL] = hadron->v_T[vALL] - v_dE[vALL];
  if(hadron->v_T[vALL] <= (config->Ecut_Pro * UMeV)){
    hadron->v_type[vALL] = Unknown;
    hadron_list[Hadron_ID[vALL]].type = Unknown;
  }
  Update_Hadron(hadron);

  // Update position and direction
  Update_position(hadron, v_step);
  Update_direction(hadron, v_theta, v_phi);


  // Hard interaction

  ALIGNED_(64) int v_interaction_type[VLENGTH];

  get_interaction_type(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, v_section, RNG_Stream, config, v_interaction_type);
  if(v_step[vALL] == v_step_max[vALL]) v_interaction_type[vALL] = 0; // force ficitious interaction if step >= step_max

  // Ionization
  ALIGNED_(64) VAR_COMPUTE v_dE_hard[VLENGTH];
  Compute_Ionization_Energy(hadron, (config->Te_Min*UMeV), RNG_Stream, v_dE_hard);
  if(hadron->v_type[vALL] == Unknown) v_dE_hard[vALL] = 0;
  if(v_interaction_type[vALL] == 1) hadron->v_T[vALL] = hadron->v_T[vALL] - v_dE_hard[vALL];

  // Nuclear interaction
  DATA_Scoring tmp;
  int i,j;
  int previous_Nbr_hadrons;
  for(i=0; i<VLENGTH; i++){
    if(hadron->v_type[i] != Unknown && v_interaction_type[i] == 2){
      previous_Nbr_hadrons = *Nbr_hadrons;
      Compute_Nuclear_interaction(i, hadron, material, machine->RS_Material, hadron_list, Nbr_hadrons, 0, &tmp, RNG_Stream, config);
      if(hadron->v_type[i] == Unknown) hadron_list[Hadron_ID[i]].type = Unknown;
      for(j=previous_Nbr_hadrons; j<*Nbr_hadrons; j++){
         layer_data[j] = layer_data[Hadron_ID[i]];
         field_data[j] = field_data[Hadron_ID[i]];
      }
    }
  }

  if(hadron->v_T[vALL] <= (config->Ecut_Pro * UMeV)){
    hadron->v_type[vALL] = Unknown;
    hadron_list[Hadron_ID[vALL]].type = Unknown;
  }

}
