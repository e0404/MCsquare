/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_hadron.h"


void hadron_step(Hadron *hadron, DATA_Scoring *scoring, Materials *material, DATA_CT *ct, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, VSLStreamStatePtr RNG_Stream, DATA_config *config){
  
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

  // Init variables
  int i;

  // Compute CT index and remove particles out of geometry
  ALIGNED_(64) int v_index[VLENGTH];
  get_CT_Offset(hadron, ct, v_index);

  ALIGNED_(64) int v_material_label[VLENGTH];
  v_material_label[vALL] = ct->material[v_index[vALL]];

  ALIGNED_(64) VAR_COMPUTE v_init_density[VLENGTH];
  v_init_density[vALL] = ct->density[v_index[vALL]];

  ALIGNED_(64) VAR_COMPUTE v_N_el[VLENGTH];
//v_N_el[vALL] = material[v_material_label[vALL]].N_el * v_init_density[vALL];
  for(i=0; i<VLENGTH; i++){
    v_N_el[i] = material[v_material_label[i]].N_el * v_init_density[i];
  }



  // calcul de la section efficace

  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step_max[VLENGTH];

  #if EM_Method==EM_FIPPEL
    ALIGNED_(64) int v_water_label[VLENGTH];
    v_water_label[vALL] = WATER_LABEL;

    ALIGNED_(64) VAR_COMPUTE v_N_el_water[VLENGTH];
    v_N_el_water[vALL] = material[v_water_label[vALL]].N_el * v_init_density[vALL];

    ALIGNED_(64) VAR_COMPUTE v_StpCorr[VLENGTH];
    Fippel_Stop_Pow_correction(hadron, v_init_density, v_StpCorr);

    Total_Stop_Pow(hadron, material, v_water_label, v_stop_pow);
    v_stop_pow[vALL] = v_init_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_StpCorr[vALL] * v_stop_pow[vALL];

  #else	// full PSTAR
    Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
    v_stop_pow[vALL] = v_init_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_stop_pow[vALL];
  #endif


  #if InterfaceCrossing==FictitiousInteraction
    ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
    Dist_To_Material_Interface(hadron, ct, config->D_Max, v_index, v_init_density, v_Dist_Interface);

    v_step_max[vALL] = fmin(fmin(v_Dist_Interface[vALL], config->D_Max), (config->Epsilon_Max * hadron->v_T[vALL] / v_stop_pow[vALL]));

  #elif InterfaceCrossing==VoxelInterface
    ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
    Dist_To_Interface(hadron, ct, v_Dist_Interface);

    v_step_max[vALL] = fmin(fmin(v_Dist_Interface[vALL], config->D_Max), (config->Epsilon_Max * hadron->v_T[vALL] / v_stop_pow[vALL]));

  #else	// No interface or Random Hinge or Fippel Transport
    v_step_max[vALL] = fmin(config->D_Max, (config->Epsilon_Max * hadron->v_T[vALL] / v_stop_pow[vALL]));
  #endif

  ALIGNED_(64) VAR_COMPUTE v_dE_max[VLENGTH];
  v_dE_max[vALL] = v_step_max[vALL] * v_stop_pow[vALL];

  ALIGNED_(64) VAR_COMPUTE v_section[VLENGTH];
  Total_Hard_Cross_Section(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, config, v_section); 
  v_section[vALL] += 1e-10;
  v_section[vALL] *= 1.017;  // facteur pour palier l'approximation.


  // calcul du SPR pour la conversion dose to water
  ALIGNED_(64) VAR_COMPUTE v_SPR[VLENGTH];
  if(config->DoseToWater == 2){
    ALIGNED_(64) int v_water_ID[VLENGTH];
    v_water_ID[vALL] = config->Water_Material_ID;
    Total_Stop_Pow(hadron, material, v_water_ID, v_SPR);
    v_SPR[vALL] = v_stop_pow[vALL] / (v_init_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_SPR[vALL]);
    if(hadron->v_type[vALL] == Unknown) v_SPR[vALL] = 1.0;
  }
  else v_SPR[vALL] = 1.0;

  // calcul de la distance pour arriver au prochain step

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  v_step[vALL] = -log(v_rnd[vALL])/v_section[vALL];

  if(v_step[vALL] > v_step_max[vALL]) v_step[vALL] = v_step_max[vALL];  // on se limite à une distance step_max

  
  // Compute CSDA + MS

  ALIGNED_(64) VAR_COMPUTE v_mean_dE[VLENGTH];					// mean energy loss
  #if EM_Method==EM_FIPPEL
    Compute_dE2_Fippel(hadron, v_N_el_water, v_init_density, material, (config->Te_Min*UMeV), v_water_label, v_step, v_mean_dE);  
  #else
    Compute_dE2(hadron, v_N_el, v_init_density, material, (config->Te_Min*UMeV), v_material_label, v_step, v_mean_dE);  
  #endif

  ALIGNED_(64) VAR_COMPUTE v_straggling[VLENGTH];
  Compute_Energy_straggling(hadron, v_N_el, (config->Te_Min*UMeV), v_step, v_straggling);		// energy straggling
  v_straggling[vALL] = sqrt(v_straggling[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_dE[VLENGTH];
  rand_normal(RNG_Stream, v_dE, v_mean_dE, v_straggling);					// energie perdue

  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
//v_X0[vALL] = material[v_material_label[vALL]].X0 / v_init_density[vALL];			// longueur de radiation
for(i=0; i<VLENGTH; i++){
v_X0[i] = material[v_material_label[i]].X0 / v_init_density[i];
}

  ALIGNED_(64) VAR_COMPUTE v_MS[VLENGTH];
  Compute_MS_Fippel(hadron, v_step, v_X0, v_MS);						// MS
//v_MS[vALL] = (config->MCS_const*UMeV * hadron->v_charge[vALL] / (hadron->v_beta2[vALL]*hadron->v_gamma[vALL]*MC2_PRO)) * sqrt(v_step[vALL]/v_X0[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  rand_normal_zero(RNG_Stream, v_theta, v_MS);							// deviation angle (theta)

  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];					// déviation angle (phi)
  rand_uniform(RNG_Stream, v_phi);
  v_phi[vALL] = 2*M_PI*v_phi[vALL];


  // Dépot de l'énergie à un point choisi aléatoirement dans le step

  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE v_tau[VLENGTH];
  v_tau[vALL] = v_rnd[vALL] * v_step[vALL];

  ALIGNED_(64) int v_hinge_index[VLENGTH];

  #if InterfaceCrossing==RandomHinge
	// Gestion des interfaces par la méthode du Random Hinge
    ALIGNED_(64) VAR_COMPUTE v_mask[VLENGTH];
    v_mask[vALL] = 1.0;
    CT_Transport_Random_Hinge(hadron, ct, v_step, v_tau, v_index, v_hinge_index, v_init_density, v_mask);
  #elif InterfaceCrossing==FictitiousInteraction
	// Gestion des interfaces par interaction fictives
    CT_Transport(hadron, ct, v_step, v_tau, v_index, v_hinge_index, v_init_density);
  #elif InterfaceCrossing==VoxelInterface
    v_hinge_index[vALL] = v_index[vALL];
    Update_position(hadron, v_step);
  #else
        // NoInterface
    CT_Transport_SPR(hadron, ct, material, v_step, v_tau, v_index, v_hinge_index, v_init_density);
  #endif


  // scoring de la perte d'énergie

  if(hadron->v_type[vALL] == Unknown) v_dE[vALL] = 0;
  #if InterfaceCrossing==RandomHinge
    if(v_hinge_index[vALL] == -1){	// si on croise une interface, on ne continue pas le step.
      v_theta[vALL] = 0.0;
      v_dE[vALL] = 0;
    }
  #endif
  if((hadron->v_T[vALL] - v_dE[vALL]) <= (config->Ecut_Pro * UMeV)){
    v_dE[vALL] = hadron->v_T[vALL];
    hadron->v_type[vALL] = Unknown;
  }
  hadron->v_T[vALL] = hadron->v_T[vALL] - v_dE[vALL];

  for(i=0; i<VLENGTH; i++) Energy_Scoring(scoring, v_hinge_index[i], hadron->v_M[i], v_dE[i], v_SPR[i]);

/*
  for(i=0; i<VLENGTH; i++){
    if(hadron->v_type[i] == Unknown) continue;
    #if InterfaceCrossing==RandomHinge
      if(v_hinge_index[i] == -1){	// si on croise une interface, on ne continue pas le step.
	v_theta[i] = 0.0;
	continue;
      }
    #endif

    if((hadron->v_T[i] - v_dE[i]) <= (config->Ecut_Pro * UMeV)){
      v_dE[i] = hadron->v_T[i];
      hadron->v_type[i] = Unknown;
    }
    else{
      hadron->v_T[i] = hadron->v_T[i] - v_dE[i];
    }
    Energy_Scoring(scoring, v_hinge_index[i], hadron->v_M[i], v_dE[i], v_SPR[i]);
  }
*/

  Update_Hadron(hadron);


  // Compute CT index and remove particles out of geometry
  get_CT_Offset(hadron, ct, v_index);

  // Changement de direction (Multiple scattering)
  Update_direction(hadron, v_theta, v_phi);



  // Interaction HARD
  ALIGNED_(64) VAR_COMPUTE v_dE_hard[VLENGTH];
  v_dE_hard[vALL] = 0.0;

  ALIGNED_(64) int v_interaction_type[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_tmp[VLENGTH];

  get_interaction_type(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, v_section, RNG_Stream, config, v_interaction_type);
  if(v_step[vALL] == v_step_max[vALL]) v_interaction_type[vALL] = 0;
  // Si step > step_max : step = step_max et force interaction fictive
  // Sinon, on détermine aléatoirement le type d'interaction

  // interaction discrète d'ionisation
  Compute_Ionization_Energy(hadron, (config->Te_Min*UMeV), RNG_Stream, v_dE_tmp);
  if(v_interaction_type[vALL] == 1) v_dE_hard[vALL] = v_dE_tmp[vALL];
  hadron->v_T[vALL] = hadron->v_T[vALL] - v_dE_hard[vALL];

  ALIGNED_(64) VAR_COMPUTE v_stop_pow2[VLENGTH];
  if(config->Score_LET == 1 && config->LET_Calculation_Method == 1){
    Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow2);
    v_stop_pow[vALL] = 0.5 * (v_stop_pow[vALL] + v_init_density[vALL] * hadron->v_charge[vALL]*hadron->v_charge[vALL] * v_stop_pow2[vALL]);
  }

  // interaction delta
  // if(v_interaction_type[vALL] == 0){
    		// on fait rien
  // }

  // Interaction Nucléaire et Scoring de la perte d'énergie
  for(i=0; i<VLENGTH; i++){
    #if InterfaceCrossing==RandomHinge
      if(v_mask[i] == 0.0) continue;	// si on croise une interface, on ne continue pas le step.
    #endif

    if(hadron->v_type[i] == Unknown) v_dE_hard[i] = 0;
    else{

      // Nuclear interactions
      if(v_interaction_type[i] == 2){
        v_dE_hard[i] = Compute_Nuclear_interaction(i, hadron, material, v_material_label[i], secondary_hadron, Nbr_secondaries, v_index[i], scoring, RNG_Stream, config);

        if(hadron->v_T[i] <= (config->Ecut_Pro * UMeV)){
	  hadron->v_type[i] = Unknown;
	  v_dE_hard[i] += hadron->v_T[i];
        }

	Energy_Scoring(scoring, v_index[i], hadron->v_M[i], v_dE_hard[i], 1.0);
      }

      // Scoring for EM interactions
      else if(v_interaction_type[i] == 1 || v_interaction_type[i] == 0){
        if(hadron->v_T[i] <= (config->Ecut_Pro * UMeV)){
          v_dE_hard[i] += hadron->v_T[i];
          hadron->v_type[i] = Unknown;
        }

        if(config->Score_LET == 1) LET_Scoring(scoring, v_hinge_index[i], hadron->v_M[i], v_dE[i]+v_dE_hard[i], v_step[i], v_stop_pow[i], config);

        Energy_Scoring(scoring, v_index[i], hadron->v_M[i], v_dE_hard[i], v_SPR[i]);
      }

    }

  }

//  Update_Hadron(hadron);


  return;
}








