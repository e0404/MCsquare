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


void hadron_step(Hadron *hadron, DATA_Scoring *scoring, Materials *material, DATA_CT *ct, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, VAR_RND_SEED RNG_Stream, DATA_config *config){
  
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
  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) int v_material_label[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_init_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_N_el[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step_max[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_max[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_section[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_SPR[VLENGTH];
  ALIGNED_(64) int v_water_ID[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_mean_dE[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_straggling[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_X0[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_MS[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_tau[VLENGTH];
  ALIGNED_(64) int v_hinge_index[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_hard[VLENGTH];
  ALIGNED_(64) int v_interaction_type[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_dE_tmp[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_stop_pow2[VLENGTH];

  // Compute CT index and remove particles out of geometry
  get_CT_Offset(hadron, ct, v_index);
  
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_material_label[v] = ct->material[v_index[v]];
    v_init_density[v] = ct->density[v_index[v]];
    v_N_el[v] = material[v_material_label[v]].N_el * v_init_density[v];
  }

  // calcul de la section efficace

  #if EM_Method==EM_FIPPEL
    ALIGNED_(64) int v_water_label[VLENGTH];
    ALIGNED_(64) VAR_COMPUTE v_N_el_water[VLENGTH];
    ALIGNED_(64) VAR_COMPUTE v_StpCorr[VLENGTH];

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_water_label[v] = WATER_LABEL;
      v_N_el_water[v] = material[v_water_label[v]].N_el * v_init_density[v];
    }

    Fippel_Stop_Pow_correction(hadron, v_init_density, v_StpCorr);
    Total_Stop_Pow(hadron, material, v_water_label, v_stop_pow);

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_stop_pow[v] = v_init_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_StpCorr[v] * v_stop_pow[v];
    }

  #else	// full PSTAR
    Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_stop_pow[v] = v_init_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_stop_pow[v];
    }
  #endif


  #if InterfaceCrossing==FictitiousInteraction
    ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
    Dist_To_Material_Interface(hadron, ct, config->D_Max, v_index, v_init_density, v_Dist_Interface);

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_step_max[v] = fmin(fmin(v_Dist_Interface[v], config->D_Max), (config->Epsilon_Max * hadron->v_T[v] / v_stop_pow[v]));
    }

  #elif InterfaceCrossing==VoxelInterface
    ALIGNED_(64) VAR_COMPUTE v_Dist_Interface[VLENGTH];
    Dist_To_Interface(hadron, ct, v_Dist_Interface);

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_step_max[v] = fmin(fmin(v_Dist_Interface[v], config->D_Max), (config->Epsilon_Max * hadron->v_T[v] / v_stop_pow[v]));
    }

  #else	// No interface or Random Hinge or Fippel Transport

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_step_max[v] = fmin(config->D_Max, (config->Epsilon_Max * hadron->v_T[v] / v_stop_pow[v]));
    }
  #endif

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_dE_max[v] = v_step_max[v] * v_stop_pow[v];
  }

  Total_Hard_Cross_Section(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, config, v_section); 
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_section[v] += 1e-10;
    v_section[v] *= 1.017;  // facteur pour palier l'approximation.
  }

  // calcul du SPR pour la conversion dose to water
  if(config->DoseToWater == 2){
    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_water_ID[v] = config->Water_Material_ID;
    }

    Total_Stop_Pow(hadron, material, v_water_ID, v_SPR);
    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_SPR[v] = v_stop_pow[v] / (v_init_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_SPR[v]);
      if(hadron->v_type[v] == Unknown) v_SPR[v] = 1.0;
    }
  }
  else{
    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_SPR[v] = 1.0;
    }
  }
  // calcul de la distance pour arriver au prochain step

  rand_uniform(RNG_Stream, v_rnd);

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_step[v] = -log(v_rnd[v])/v_section[v];
    if(v_step[v] > v_step_max[v]) v_step[v] = v_step_max[v];  // on se limite à une distance step_max
  }
  
  // Compute CSDA + MS
  // mean energy loss
  #if EM_Method==EM_FIPPEL
    Compute_dE2_Fippel(hadron, v_N_el_water, v_init_density, material, (config->Te_Min*UMeV), v_water_label, v_step, v_mean_dE);  
  #else
    Compute_dE2(hadron, v_N_el, v_init_density, material, (config->Te_Min*UMeV), v_material_label, v_step, v_mean_dE);  
  #endif


  Compute_Energy_straggling(hadron, v_N_el, (config->Te_Min*UMeV), v_step, v_straggling);		// energy straggling
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_straggling[v] = sqrt(v_straggling[v]);
  }
  rand_normal(RNG_Stream, v_dE, v_mean_dE, v_straggling);					// energie perdue
			
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_X0[v] = material[v_material_label[v]].X0 / v_init_density[v];  // longueur de radiation
  }

  Compute_MS_Fippel(hadron, v_step, v_X0, v_MS);	// MS
  rand_normal_zero(RNG_Stream, v_theta, v_MS);		// deviation angle (theta)  
  rand_uniform(RNG_Stream, v_phi);			// déviation angle (phi)

  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_phi[v] = 2*M_PI*v_phi[v];
  }

  // Dépot de l'énergie à un point choisi aléatoirement dans le step
  rand_uniform(RNG_Stream, v_rnd);

  ALIGNED_(64) VAR_COMPUTE scoring_x[VLENGTH], scoring_y[VLENGTH], scoring_z[VLENGTH];
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_tau[v] = v_rnd[v] * v_step[v];
    scoring_x[v] = hadron->v_x[v] + v_tau[v] * hadron->v_u[v];
    scoring_y[v] = hadron->v_y[v] + v_tau[v] * hadron->v_v[v];
    scoring_z[v] = hadron->v_z[v] + v_tau[v] * hadron->v_w[v];
  }

  #if InterfaceCrossing==RandomHinge // Gestion des interfaces par la méthode du Random Hinge
    ALIGNED_(64) VAR_COMPUTE v_mask[VLENGTH];

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_mask[v] = 1.0;
    }
    CT_Transport_Random_Hinge(hadron, ct, v_step, v_tau, v_index, v_hinge_index, v_init_density, v_mask);
  #elif InterfaceCrossing==FictitiousInteraction // Gestion des interfaces par interaction fictives
    CT_Transport(hadron, ct, v_step, v_tau, v_index, v_hinge_index, v_init_density);
  #elif InterfaceCrossing==VoxelInterface

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_hinge_index[v] = v_index[v];
    }
    Update_position(hadron, v_step);
  #else // NoInterface
    CT_Transport_SPR(hadron, ct, material, v_step, v_tau, v_index, v_hinge_index, v_init_density);
  #endif


  // scoring de la perte d'énergie
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    if(hadron->v_type[v] == Unknown) v_dE[v] = 0;

    #if InterfaceCrossing==RandomHinge
      if(v_hinge_index[v] == -1){	// si on croise une interface, on ne continue pas le step.
        v_theta[v] = 0.0;
        v_dE[v] = 0;
      }
    #endif

    if((hadron->v_T[v] - v_dE[v]) <= (config->Ecut_Pro * UMeV)){
      v_dE[v] = hadron->v_T[v];
      hadron->v_type[v] = Unknown;
    }
    hadron->v_T[v] = hadron->v_T[v] - v_dE[v];
  }
    //Energy_Scoring(scoring, v_hinge_index[v], hadron->v_M[v], v_dE[v], v_SPR[v]);
  if(config->Independent_scoring_grid == 0) Energy_Scoring_from_index(scoring, v_hinge_index, hadron->v_M, v_dE, v_init_density, v_SPR, config);
  else Energy_Scoring_from_coordinates(scoring, scoring_x, scoring_y, scoring_z, hadron->v_M, v_dE, v_init_density, v_SPR, config);


  Update_Hadron(hadron);

  // Compute CT index and remove particles out of geometry
  get_CT_Offset(hadron, ct, v_index);

  // Changement de direction (Multiple scattering)
  Update_direction(hadron, v_theta, v_phi);



  // Interaction HARD
  get_interaction_type(hadron, material, v_material_label, v_N_el, v_init_density, (config->Te_Min*UMeV), v_dE_max, v_section, RNG_Stream, config, v_interaction_type);
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    v_dE_hard[v] = 0.0;
    if(v_step[v] == v_step_max[v]) v_interaction_type[v] = 0;
    // Si step > step_max : step = step_max et force interaction fictive
    // Sinon, on détermine aléatoirement le type d'interaction
  }

  // interaction discrète d'ionisation
  Compute_Ionization_Energy(hadron, (config->Te_Min*UMeV), RNG_Stream, v_dE_tmp);
  #pragma omp simd
  for(int v = 0; v<VLENGTH; v++){
    if(v_interaction_type[v] == 1) v_dE_hard[v] = v_dE_tmp[v];
    hadron->v_T[v] = hadron->v_T[v] - v_dE_hard[v];
  }

  if(config->Score_LET == 1 && config->LET_Calculation_Method == 1){
    Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow2);
    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_stop_pow[v] = 0.5 * (v_stop_pow[v] + v_init_density[v] * hadron->v_charge[v]*hadron->v_charge[v] * v_stop_pow2[v]);
    }
  }

  // interaction delta
  // if(v_interaction_type[vALL] == 0){
    		// on fait rien
  // }

  // Interaction Nucléaire et Scoring de la perte d'énergie
  for(int v = 0; v<VLENGTH; v++){
    #if InterfaceCrossing==RandomHinge
      if(v_mask[v] == 0.0) continue;	// si on croise une interface, on ne continue pas le step.
    #endif

    if(hadron->v_type[v] == Unknown) v_dE_hard[v] = 0;
    else{

      // Nuclear interactions
      if(v_interaction_type[v] == 2){
        v_dE_hard[v] = Compute_Nuclear_interaction(v, hadron, material, v_material_label[v], secondary_hadron, Nbr_secondaries, scoring, RNG_Stream, config);

        if(hadron->v_T[v] <= (config->Ecut_Pro * UMeV)){
	        hadron->v_type[v] = Unknown;
	        v_dE_hard[v] += hadron->v_T[v];
        }

	      //Energy_Scoring(scoring, v_index[v], hadron->v_M[v], v_dE_hard[v], 1.0);
        v_SPR[v] = 1.0; //SPR not use for dose-to-water conversion of nuclear interaction
      }

      // Scoring for EM interactions
      else if(v_interaction_type[v] == 1 || v_interaction_type[v] == 0){
        if(hadron->v_T[v] <= (config->Ecut_Pro * UMeV)){
          v_dE_hard[v] += hadron->v_T[v];
          hadron->v_type[v] = Unknown;
        }

        //if(config->Score_LET == 1) LET_Scoring(scoring, v_hinge_index[v], hadron->v_M[v], v_dE[v]+v_dE_hard[v], v_step[v], v_stop_pow[v], config);
        if(config->Score_LET == 1) LET_Scoring(scoring, scoring_x[v], scoring_y[v], scoring_z[v], hadron->v_M[v], v_dE[v]+v_dE_hard[v], v_step[v], v_stop_pow[v], config);


        //Energy_Scoring(scoring, v_index[v], hadron->v_M[v], v_dE_hard[v], v_SPR[v]);
      }

    }

  }
    
  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  #pragma omp simd
  for(int v=0; v<VLENGTH; v++){
    v_density[v] = ct->density[v_index[v]];
  }
  
  if(config->Independent_scoring_grid == 0) Energy_Scoring_from_index(scoring, v_index, hadron->v_M, v_dE_hard, v_density, v_SPR, config);
  else Energy_Scoring_from_coordinates(scoring, hadron->v_x, hadron->v_y, hadron->v_z, hadron->v_M, v_dE_hard, v_density, v_SPR, config);

  return;
}








