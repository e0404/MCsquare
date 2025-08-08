/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_geometry.h"

void verif_position(Hadron *hadron, DATA_CT *ct){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);  
  __assume_aligned(&hadron->v_type, 64);

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    if(	hadron->v_x[v] < 0 || hadron->v_y[v] < 0 || hadron->v_z[v] < 0 || 
	hadron->v_x[v] >= ct->Length[0] || hadron->v_y[v] >= ct->Length[1] || hadron->v_z[v] >= ct->Length[2]) hadron->v_type[v] = Unknown;
  }

  return;
}


void get_CT_Offset(Hadron *hadron, DATA_CT *ct, int *v_index){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64);
  __assume_aligned(&hadron->v_type, 64);  
  __assume_aligned(v_index, 64);

  // Calcul de l'offset : Offset = x + ct->Nx * y + ct->Nx * ct->Ny * z
  // L'axe x du repère simulation ne correspond pas à l'axe x du repère CT : x_simu = -x_ct + Lx
  // Conversion position -> index CT : index = floor(x/dx)

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_index[v] = 	(int)floor( (-hadron->v_x[v] + ct->Length[0]) / ct->VoxelLength[0] ) 
			+ ct->GridSize[0] * (int)floor( hadron->v_y[v] / ct->VoxelLength[1] ) 
			+ ct->GridSize[0] * ct->GridSize[1] * (int)floor( hadron->v_z[v] / ct->VoxelLength[2] );

    if(	hadron->v_x[v] < 0 || hadron->v_y[v] < 0 || hadron->v_z[v] < 0 || 
	hadron->v_x[v] >= ct->Length[0] || hadron->v_y[v] >= ct->Length[1] || hadron->v_z[v] >= ct->Length[2]) hadron->v_type[v] = Unknown;

    if(v_index[v] < 0 || v_index[v] > ct->Nbr_voxels) hadron->v_type[v] = Unknown;

    if(hadron->v_type[v] == Unknown) v_index[v] = 0;
  }

  return;
}


void Dist_To_Material_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE dist_max, int *v_init_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 
  __assume_aligned(&hadron->v_type, 64);  

  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_result, 64);


  ALIGNED_(64) VAR_COMPUTE v_mass_distance[VLENGTH];
  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) int v_index2[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);

  int run = 0;
  int v;

  #pragma omp simd reduction(+:run)
  for(v = 0; v<VLENGTH; v++){
    v_mass_distance[v] = 0;
    v_index[v] = v_init_index[v];
    v_index2[v] = v_init_index[v];
    v_run[v] = 1.0;
    if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
    run += v_run[v];
  }

  while(run != 0){

    Dist_To_Interface(&tmp, ct, v_step);
    
    run = 0;

    #pragma omp simd reduction(+:run)
    for(v = 0; v<VLENGTH; v++){
      if((v_mass_distance[v] / v_init_density[v]) < dist_max && ct->material[v_index[v]] == ct->material[v_index2[v]]){
        v_index[v] = v_index2[v];
        v_mass_distance[v] += v_step[v] * ct->density[v_index[v]];
      }
      else{
        v_step[v] = 0.0;
        v_run[v] = 0.0;
      }

      run += v_run[v];
    }

    Update_position(&tmp, v_step);
    get_CT_Offset(&tmp, ct, v_index2);
  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_result[v] = v_mass_distance[v] / v_init_density[v];
  }

  return;
}


void Dist_To_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_result){
  
  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 
  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 
  __assume_aligned(v_result, 64); 


  ALIGNED_(64) VAR_COMPUTE v_DistX[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_DistY[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_DistZ[VLENGTH];

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_DistX[v] = fabs(((floor(hadron->v_x[v]/ct->VoxelLength[0]) + (hadron->v_u[v] > 0)) * ct->VoxelLength[0] - hadron->v_x[v])/hadron->v_u[v]);
    v_DistY[v] = fabs(((floor(hadron->v_y[v]/ct->VoxelLength[1]) + (hadron->v_v[v] > 0)) * ct->VoxelLength[1] - hadron->v_y[v])/hadron->v_v[v]);
    v_DistZ[v] = fabs(((floor(hadron->v_z[v]/ct->VoxelLength[2]) + (hadron->v_w[v] > 0)) * ct->VoxelLength[2] - hadron->v_z[v])/hadron->v_w[v]);
  }
  // Separated simd loops to avoid MSVC issues with the fmin functions
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    // Add safety increment to compensate for rounding errors and to be sure to pass the interface (2e-4 for float, 1.5e-8 for double);
    #if VAR_COMPUTE_PRECISION==1
      v_result[v] = fmin(v_DistX[v], fmin(v_DistY[v], v_DistZ[v])); 
      if(v_result[v] < 1e-3) v_result[v] += 2e-4;
      else v_result[v] += 5e-5;
    #else
      v_result[v] = fmin(v_DistX[v], fmin(v_DistY[v], v_DistZ[v])) + 1.5e-8; 
    #endif
  }


  return;
}


void Update_position(Hadron *hadron, VAR_COMPUTE *v_step){

  __assume_aligned(&hadron->v_x, 64);
  __assume_aligned(&hadron->v_y, 64);
  __assume_aligned(&hadron->v_z, 64); 
  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 
  __assume_aligned(v_step, 64); 

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    hadron->v_x[v] += v_step[v] * hadron->v_u[v];
    hadron->v_y[v] += v_step[v] * hadron->v_v[v];
    hadron->v_z[v] += v_step[v] * hadron->v_w[v];
  }

  return;
}


void CT_Transport(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(&hadron->v_type, 64); 

  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];

  int run = 0;
  int v;

  Dist_To_Interface(hadron, ct, v_step);

  #pragma omp simd reduction(+:run)
  for(v = 0; v<VLENGTH; v++){
    v_MassDistance[v] = v_s[v] * v_init_density[v];
    v_index[v] = v_init_index[v];
    v_density[v] = ct->density[v_index[v]];
    v_run[v] = 1.0;
    if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
    v_hinge_index[v] = -1;
    v_HingeDistance[v] = v_MassDistance[v] - v_tau[v] * v_init_density[v];
    run += v_run[v];
  }

  while(run != 0){

    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      if(v_MassDistance[v] > v_step[v] * v_density[v]) v_MassDistance[v] -= v_step[v] * v_density[v];
      else{
        v_run[v] = 0.0;
        v_step[v] = 0.0;
      }

      if(v_MassDistance[v] < v_HingeDistance[v] && v_hinge_index[v] == -1){
        v_hinge_index[v] = v_index[v];
      }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    get_CT_Offset(hadron, ct, v_index);
    Dist_To_Interface(hadron, ct, v_step);

    run = 0;

    #pragma omp simd reduction(+:run)
    for(v = 0; v<VLENGTH; v++){
      if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
      v_density[v] = ct->density[v_index[v]];
      if(v_index[v] > ct->Nbr_voxels || v_index[v] < 0){
        v_run[v] = 0.0;
        v_hinge_index[v] = 0;
      }
      run +=v_run[v];
    }
  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_step[v] = v_MassDistance[v] / v_density[v];
    if(v_hinge_index[v] == -1) v_hinge_index[v] = v_index[v];
  }

  Update_position(hadron, v_step);

  return;
}


void CT_Transport_SPR(Hadron *hadron, DATA_CT *ct, Materials *material, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(&hadron->v_type, 64); 


  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) int v_material_label[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) int v_data_index[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];

  Dist_To_Interface(hadron, ct, v_step);

  int run = 0;
  int v;

  #pragma omp simd reduction(+:run)
  for(v = 0; v<VLENGTH; v++){
    v_index[v] = v_init_index[v];
    v_material_label[v] = ct->material[v_index[v]];
    v_density[v] = v_init_density[v];
    v_data_index[v] = (int)ceil( hadron->v_T[v] / (UMeV*PSTAR_BIN * hadron->v_mass[v]));
    v_stop_pow[v] = (VAR_COMPUTE)material[v_material_label[v]].Stop_Pow[v_data_index[v]];
    v_MassDistance[v] = v_s[v] * v_init_density[v] * v_stop_pow[v];
    v_HingeDistance[v] = v_MassDistance[v] - v_tau[v] * v_init_density[v] * v_stop_pow[v];
    v_hinge_index[v] = -1;
    if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
    else v_run[v] = 1.0;

    run += v_run[v];
  }

  while(run != 0){

    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      if(v_MassDistance[v] > v_step[v] * v_density[v] * v_stop_pow[v]){ 
        v_MassDistance[v] -= v_step[v] * v_density[v] * v_stop_pow[v];
      }
      else{
        v_run[v] = 0.0;
        v_step[v] = 0.0;
      }

      if(v_MassDistance[v] < v_HingeDistance[v] && v_hinge_index[v] == -1){
        v_hinge_index[v] = v_index[v];
      }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    get_CT_Offset(hadron, ct, v_index);
    Dist_To_Interface(hadron, ct, v_step);

    run = 0;

    #pragma omp simd reduction(+:run)
    for(v = 0; v<VLENGTH; v++){
      if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
      v_density[v] = ct->density[v_index[v]];
      if(v_material_label[v] != ct->material[v_index[v]]){
	v_material_label[v] = ct->material[v_index[v]];
  	v_stop_pow[v] = (VAR_COMPUTE)material[v_material_label[v]].Stop_Pow[v_data_index[v]];
      }

      if(v_index[v] > ct->Nbr_voxels || v_index[v] < 0){
        v_run[v] = 0.0;
        v_hinge_index[v] = 0;
      }

      run += v_run[v];
    }

  }

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_step[v] = v_MassDistance[v] / (v_density[v] * v_stop_pow[v]);
    if(v_hinge_index[v] == -1) v_hinge_index[v] = v_index[v];
  }

  Update_position(hadron, v_step);

  return;
}


void CT_Transport_Random_Hinge(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_mask){

  __assume_aligned(v_s, 64);
  __assume_aligned(v_tau, 64);
  __assume_aligned(v_init_index, 64);
  __assume_aligned(v_init_density, 64);
  __assume_aligned(v_hinge_index, 64);
  __assume_aligned(v_mask, 64);
  __assume_aligned(&hadron->v_type, 64); 

  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) int v_index[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  ALIGNED_(64) unsigned short int v_material[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];

  Dist_To_Interface(hadron, ct, v_step);

  int run = 0;
  int v;

  #pragma omp simd reduction(+:run)
  for(v = 0; v<VLENGTH; v++){
    v_MassDistance[v] = v_s[v] * v_init_density[v];
    v_index[v] = v_init_index[v];
    v_density[v] = ct->density[v_index[v]];
    v_run[v] = v_mask[v];
    if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
    v_material[v] = ct->material[v_index[v]];
    v_hinge_index[v] = -1;
    v_HingeDistance[v] = v_MassDistance[v] - v_tau[v] * v_init_density[v];
    run += v_run[v];
  }

  while(run != 0){

    int v;
    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      if(v_MassDistance[v] > v_step[v] * v_density[v] && v_run[v] != 0.0){
        v_MassDistance[v] -= v_step[v] * v_density[v];
      }
      else{
        v_run[v] = 0.0;
        v_step[v] = 0.0;
      }

      if(v_MassDistance[v] < v_HingeDistance[v] && v_hinge_index[v] == -1){
        v_hinge_index[v] = v_index[v];
      }
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    get_CT_Offset(hadron, ct, v_index);
    Dist_To_Interface(hadron, ct, v_step);

    run = 0;

    #pragma omp simd reduction(+:run)
    for(v = 0; v<VLENGTH; v++){
      if(hadron->v_type[v] == Unknown) v_run[v] = 0.0;
      v_density[v] = ct->density[v_index[v]];
      if(v_index[v] > ct->Nbr_voxels || v_index[v] < 0) v_run[v] = 0.0;
      else if(v_material[v] != ct->material[v_index[v]]){
        v_run[v] = 0.0;
        v_mask[v] = 0.0;
      }

      run += v_run[v];
    }
  }  

  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    if(v_mask[v] == 0.0) v_MassDistance[v] = 0.0;
    v_step[v] = v_MassDistance[v] / v_density[v];
    if(v_hinge_index[v] == -1 && v_mask[v] != 0.0) v_hinge_index[v] = v_index[v];
  }

  Update_position(hadron, v_step);

  return;
}


void Update_direction(Hadron *hadron, VAR_COMPUTE *v_theta, VAR_COMPUTE *v_phi){

  __assume_aligned(v_theta, 64); 
  __assume_aligned(v_phi, 64); 

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 

  ALIGNED_(64) VAR_COMPUTE v_cosT[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_sinT[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_cosP[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_sinP[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_prev_u[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_norme[VLENGTH];

  int v;
  #pragma omp simd
  for(v = 0; v<VLENGTH; v++){
    v_cosT[v] = cos(v_theta[v]);
    v_sinT[v] = sin(v_theta[v]);
    v_cosP[v] = cos(v_phi[v]);
    v_sinP[v] = sin(v_phi[v]);
    v_prev_u[v] = hadron->v_u[v];

    if(hadron->v_w[v] < 0.999999 && hadron->v_w[v] > -0.999999){	// Si direction non parallèle à l'axe z

      hadron->v_u[v] = (v_prev_u[v]*v_cosT[v] 
			+ (v_sinT[v] / (sqrt(1.0-hadron->v_w[v]*hadron->v_w[v]))) * (v_prev_u[v]*hadron->v_w[v]*v_cosP[v] - hadron->v_v[v]*v_sinP[v]));

      hadron->v_v[v] = (hadron->v_v[v]*v_cosT[v] 
			+ (v_sinT[v] / (sqrt(1.0-hadron->v_w[v]*hadron->v_w[v]))) * (hadron->v_v[v]*hadron->v_w[v]*v_cosP[v] + v_prev_u[v]*v_sinP[v]));

      hadron->v_w[v] = (hadron->v_w[v]*v_cosT[v] - sqrt(1.0-hadron->v_w[v]*hadron->v_w[v])*v_sinT[v]*v_cosP[v]);

    }
    else{
      hadron->v_v[v] = v_sinT[v] * v_sinP[v];

      if(hadron->v_w[v] > 0){		// Si direction parallère à l'axe z
        hadron->v_u[v] = v_sinT[v] * v_cosP[v];
        hadron->v_w[v] = v_cosT[v];
      }
      else{				// Si direction antiparallère à l'axe z
        hadron->v_u[v] = -v_sinT[v] * v_cosP[v];
        hadron->v_w[v] = -v_cosT[v];
      }
    }
  }

  // for loop separated to avoid MSVC issues with the sqrt function
  #pragma omp simd
  for (v = 0; v < VLENGTH; v++) {
    // Si la norme dévie trop de 1, on renormalise
    v_norme[v] = sqrt(hadron->v_u[v]*hadron->v_u[v] + hadron->v_v[v]*hadron->v_v[v] + hadron->v_w[v]*hadron->v_w[v]);
    hadron->v_u[v] = hadron->v_u[v] / v_norme[v];
    hadron->v_v[v] = hadron->v_v[v] / v_norme[v];
    hadron->v_w[v] = hadron->v_w[v] / v_norme[v];
  }

  return;
}


void Update_buffer_direction(Hadron_buffer *secondary_hadron, VAR_COMPUTE theta, VAR_COMPUTE phi){
  VAR_COMPUTE cosT = cos(theta);
  VAR_COMPUTE sinT = sin(theta);
  VAR_COMPUTE cosP = cos(phi);
  VAR_COMPUTE sinP = sin(phi);

  if(fabs(secondary_hadron->w) < 0.999999){	// Si direction non parallèle à l'axe z

    VAR_COMPUTE Prev_u = secondary_hadron->u;	// valeur initiale de u
					// nécessaire car particule->u est mis à jour
    secondary_hadron->u = Prev_u*cosT + (sinT/(sqrt(1-secondary_hadron->w*secondary_hadron->w)))*(Prev_u*secondary_hadron->w*cosP - secondary_hadron->v*sinP);
    secondary_hadron->v = secondary_hadron->v*cosT + (sinT/(sqrt(1-secondary_hadron->w*secondary_hadron->w)))*(secondary_hadron->v*secondary_hadron->w*cosP + Prev_u*sinP);
    secondary_hadron->w = secondary_hadron->w*cosT - sqrt(1-secondary_hadron->w*secondary_hadron->w)*sinT*cosP;
  }
  else{
    secondary_hadron->v = sinT*sinP;

    if(secondary_hadron->w > 0){		// Si direction parallère à l'axe z
      secondary_hadron->u = sinT*cosP;
      secondary_hadron->w = cosT;
    }
    else{				// Si direction antiparallère à l'axe z
      secondary_hadron->u = -sinT*cosP;
      secondary_hadron->w = -cosT;
    }
  }

  // Si la norme dévie trop de 1, on renormalise
  VAR_COMPUTE norme = sqrt(secondary_hadron->u*secondary_hadron->u + secondary_hadron->v*secondary_hadron->v + secondary_hadron->w*secondary_hadron->w);
//  if(fabs(norme-1) > 1e-14){	

    secondary_hadron->u = secondary_hadron->u / norme;
    secondary_hadron->v = secondary_hadron->v / norme;
    secondary_hadron->w = secondary_hadron->w / norme;
//  }
}
