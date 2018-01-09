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

  if(	hadron->v_x[vALL] < 0 || hadron->v_y[vALL] < 0 || hadron->v_z[vALL] < 0 || 
	hadron->v_x[vALL] >= ct->Length[0] || hadron->v_y[vALL] >= ct->Length[1] || hadron->v_z[vALL] >= ct->Length[2]) hadron->v_type[vALL] = Unknown;

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

  v_index[vALL] = 	(int)floor( (-hadron->v_x[vALL] + ct->Length[0]) / ct->VoxelLength[0] ) 
			+ ct->GridSize[0] * (int)floor( hadron->v_y[vALL] / ct->VoxelLength[1] ) 
			+ ct->GridSize[0] * ct->GridSize[1] * (int)floor( hadron->v_z[vALL] / ct->VoxelLength[2] );

  if(	hadron->v_x[vALL] < 0 || hadron->v_y[vALL] < 0 || hadron->v_z[vALL] < 0 || 
	hadron->v_x[vALL] >= ct->Length[0] || hadron->v_y[vALL] >= ct->Length[1] || hadron->v_z[vALL] >= ct->Length[2]) hadron->v_type[vALL] = Unknown;

  if(v_index[vALL] < 0 || v_index[vALL] > ct->Nbr_voxels) hadron->v_type[vALL] = Unknown;

  if(hadron->v_type[vALL] == Unknown) v_index[vALL] = 0;

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
  v_mass_distance[vALL] = 0;

  Hadron tmp;
  Copy_Hadron_struct(&tmp, hadron);

  ALIGNED_(64) int v_index[VLENGTH];
  v_index[vALL] = v_init_index[vALL];

  ALIGNED_(64) int v_index2[VLENGTH];
  v_index2[vALL] = v_init_index[vALL];

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  v_run[vALL] = 1.0;
  if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;

  int run = __sec_reduce_add(v_run[vALL]);

  while(run != 0){

    if((v_mass_distance[vALL] / v_init_density[vALL]) < dist_max && ct->material[v_index[vALL]] == ct->material[v_index2[vALL]]){

      Dist_To_Interface(&tmp, ct, v_step);
      Update_position(&tmp, v_step);

      v_index[vALL] = v_index2[vALL];
      get_CT_Offset(&tmp, ct, v_index2);

      v_mass_distance[vALL] += v_step[vALL] * ct->density[v_index[vALL]];
    }
    else{
      v_run[vALL] = 0.0;
    }

  run = __sec_reduce_add(v_run[vALL]);
  }

  v_result[vALL] = v_mass_distance[vALL] / v_init_density[vALL];

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
  v_DistX[vALL] = fabs(((floor(hadron->v_x[vALL]/ct->VoxelLength[0]) + (hadron->v_u[vALL] > 0)) * ct->VoxelLength[0] - hadron->v_x[vALL])/hadron->v_u[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_DistY[VLENGTH];
  v_DistY[vALL] = fabs(((floor(hadron->v_y[vALL]/ct->VoxelLength[1]) + (hadron->v_v[vALL] > 0)) * ct->VoxelLength[1] - hadron->v_y[vALL])/hadron->v_v[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_DistZ[VLENGTH];
  v_DistZ[vALL] = fabs(((floor(hadron->v_z[vALL]/ct->VoxelLength[2]) + (hadron->v_w[vALL] > 0)) * ct->VoxelLength[2] - hadron->v_z[vALL])/hadron->v_w[vALL]);

  // Add safety increment to compensate for rounding errors and to be sure to pass the interface (2e-4 for float, 1.5e-8 for double);
  #if VAR_COMPUTE_PRECISION==1
    v_result[vALL] = fmin(v_DistX[vALL], fmin(v_DistY[vALL], v_DistZ[vALL])); 
    if(v_result[vALL] < 1e-3) v_result[vALL] += 2e-4;
    else v_result[vALL] += 5e-5;
  #else
    v_result[vALL] = fmin(v_DistX[vALL], fmin(v_DistY[vALL], v_DistZ[vALL])) + 1.5e-8; 
  #endif
  


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

  hadron->v_x[vALL] += v_step[vALL] * hadron->v_u[vALL];
  hadron->v_y[vALL] += v_step[vALL] * hadron->v_v[vALL];
  hadron->v_z[vALL] += v_step[vALL] * hadron->v_w[vALL];

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
  v_MassDistance[vALL] = v_s[vALL] * v_init_density[vALL];

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_index[VLENGTH];
  v_index[vALL] = v_init_index[vALL];

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  v_density[vALL] = ct->density[v_index[vALL]];

  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  v_run[vALL] = 1.0;
  if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;

  v_hinge_index[vALL] = -1;

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
  v_HingeDistance[vALL] = v_MassDistance[vALL] - v_tau[vALL] * v_init_density[vALL];

  int run = __sec_reduce_add(v_run[vALL]);

  while(run != 0){

    if(v_MassDistance[vALL] > v_step[vALL] * v_density[vALL]){ 
      v_MassDistance[vALL] -= v_step[vALL] * v_density[vALL];
    }
    else{
      v_run[vALL] = 0.0;
      v_step[vALL] = 0.0;
    }

    if(v_MassDistance[vALL] < v_HingeDistance[vALL] && v_hinge_index[vALL] == -1){
      v_hinge_index[vALL] = v_index[vALL];
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;
    get_CT_Offset(hadron, ct, v_index);
    v_density[vALL] = ct->density[v_index[vALL]];
    Dist_To_Interface(hadron, ct, v_step);

    if(v_index[vALL] > ct->Nbr_voxels || v_index[vALL] < 0){
      v_run[vALL] = 0.0;
      v_hinge_index[vALL] = 0;
    }

  run = __sec_reduce_add(v_run[vALL]);
  }

  v_step[vALL] = v_MassDistance[vALL] / v_density[vALL];
  Update_position(hadron, v_step);

  if(v_hinge_index[vALL] == -1){
    v_hinge_index[vALL] = v_index[vALL];
  }



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
  v_index[vALL] = v_init_index[vALL];

  ALIGNED_(64) int v_material_label[VLENGTH];
  v_material_label[vALL] = ct->material[v_index[vALL]];

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  v_density[vALL] = v_init_density[vALL];

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_data_index[VLENGTH];
  v_data_index[vALL] = (int)ceil( hadron->v_T[vALL] / (UMeV*PSTAR_BIN * hadron->v_mass[vALL]));

  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  int i;
  for(i=0; i<VLENGTH; i++){
    v_stop_pow[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_data_index[i]];
  }
/*
  ALIGNED_(64) VAR_COMPUTE v_stop_pow[VLENGTH];
  Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
*/


  ALIGNED_(64) VAR_COMPUTE v_MassDistance[VLENGTH];
  v_MassDistance[vALL] = v_s[vALL] * v_init_density[vALL] * v_stop_pow[vALL];

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
  v_HingeDistance[vALL] = v_MassDistance[vALL] - v_tau[vALL] * v_init_density[vALL] * v_stop_pow[vALL];


  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  v_run[vALL] = 1.0;
  if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;

  v_hinge_index[vALL] = -1;

  int run = __sec_reduce_add(v_run[vALL]);

  while(run != 0){

    if(v_MassDistance[vALL] > v_step[vALL] * v_density[vALL] * v_stop_pow[vALL]){ 
      v_MassDistance[vALL] -= v_step[vALL] * v_density[vALL] * v_stop_pow[vALL];
    }
    else{
      v_run[vALL] = 0.0;
      v_step[vALL] = 0.0;
    }

    if(v_MassDistance[vALL] < v_HingeDistance[vALL] && v_hinge_index[vALL] == -1){
      v_hinge_index[vALL] = v_index[vALL];
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;
    get_CT_Offset(hadron, ct, v_index);
    v_density[vALL] = ct->density[v_index[vALL]];
    Dist_To_Interface(hadron, ct, v_step);

    if(v_material_label[vALL] != ct->material[v_index[vALL]]){
	v_material_label[vALL] = ct->material[v_index[vALL]];
//	Total_Stop_Pow(hadron, material, v_material_label, v_stop_pow);
  	for(i=0; i<VLENGTH; i++){
  	  v_stop_pow[i] = (VAR_COMPUTE)material[v_material_label[i]].Stop_Pow[v_data_index[i]];
  	}
    }

    if(v_index[vALL] > ct->Nbr_voxels || v_index[vALL] < 0){
      v_run[vALL] = 0.0;
      v_hinge_index[vALL] = 0;
    }

  run = __sec_reduce_add(v_run[vALL]);
  }

  v_step[vALL] = v_MassDistance[vALL] / (v_density[vALL] * v_stop_pow[vALL]);
  Update_position(hadron, v_step);

  if(v_hinge_index[vALL] == -1){
    v_hinge_index[vALL] = v_index[vALL];
  }



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
  v_MassDistance[vALL] = v_s[vALL] * v_init_density[vALL];

  ALIGNED_(64) VAR_COMPUTE v_step[VLENGTH];
  Dist_To_Interface(hadron, ct, v_step);

  ALIGNED_(64) int v_index[VLENGTH];
  v_index[vALL] = v_init_index[vALL];

  ALIGNED_(64) VAR_COMPUTE v_density[VLENGTH];
  v_density[vALL] = ct->density[v_index[vALL]];

  ALIGNED_(64) VAR_COMPUTE v_run[VLENGTH];
  v_run[vALL] = v_mask[vALL];
  if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;

  ALIGNED_(64) unsigned short int v_material[VLENGTH];
  v_material[vALL] = ct->material[v_index[vALL]];

  v_hinge_index[vALL] = -1;

  ALIGNED_(64) VAR_COMPUTE v_HingeDistance[VLENGTH];
  v_HingeDistance[vALL] = v_MassDistance[vALL] - v_tau[vALL] * v_init_density[vALL];

  int run = __sec_reduce_add(v_run[vALL]);

  while(run != 0){

    if(v_MassDistance[vALL] > v_step[vALL] * v_density[vALL] && v_run[vALL] != 0.0){
      v_MassDistance[vALL] -= v_step[vALL] * v_density[vALL];
    }
    else{
      v_run[vALL] = 0.0;
      v_step[vALL] = 0.0;
    }

    if(v_MassDistance[vALL] < v_HingeDistance[vALL] && v_hinge_index[vALL] == -1){
      v_hinge_index[vALL] = v_index[vALL];
    }

    Update_position(hadron, v_step);
    verif_position(hadron, ct);
    if(hadron->v_type[vALL] == Unknown) v_run[vALL] = 0.0;
    get_CT_Offset(hadron, ct, v_index);
    v_density[vALL] = ct->density[v_index[vALL]];
    Dist_To_Interface(hadron, ct, v_step);

    if(v_index[vALL] > ct->Nbr_voxels || v_index[vALL] < 0) v_run[vALL] = 0.0;

    else if(v_material[vALL] != ct->material[v_index[vALL]]){
      v_run[vALL] = 0.0;
      v_mask[vALL] = 0.0;
    }

  run = __sec_reduce_add(v_run[vALL]);
  }  

  if(v_mask[vALL] == 0.0){
    v_MassDistance[vALL] = 0.0;
  }

  v_step[vALL] = v_MassDistance[vALL] / v_density[vALL];
  Update_position(hadron, v_step);

  if(v_hinge_index[vALL] == -1 && v_mask[vALL] != 0.0){
    v_hinge_index[vALL] = v_index[vALL];
  }

  return;
}


void Update_direction(Hadron *hadron, VAR_COMPUTE *v_theta, VAR_COMPUTE *v_phi){

  __assume_aligned(v_theta, 64); 
  __assume_aligned(v_phi, 64); 

  __assume_aligned(&hadron->v_u, 64);
  __assume_aligned(&hadron->v_v, 64);
  __assume_aligned(&hadron->v_w, 64); 

  ALIGNED_(64) VAR_COMPUTE v_cosT[VLENGTH];
  v_cosT[vALL] = cos(v_theta[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_sinT[VLENGTH];
  v_sinT[vALL] = sin(v_theta[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_cosP[VLENGTH];
  v_cosP[vALL] = cos(v_phi[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_sinP[VLENGTH];
  v_sinP[vALL] = sin(v_phi[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_prev_u[VLENGTH];
  v_prev_u[vALL] = hadron->v_u[vALL];


//  if(fabs(1.0 - hadron->v_w[vALL]) > 1e-10){	// Si direction non parallèle à l'axe z
  if(hadron->v_w[vALL] < 0.999999 && hadron->v_w[vALL] > -0.999999){	// Si direction non parallèle à l'axe z

    hadron->v_u[vALL] = (v_prev_u[vALL]*v_cosT[vALL] 
			+ (v_sinT[vALL] / (sqrt(1.0-hadron->v_w[vALL]*hadron->v_w[vALL]))) * (v_prev_u[vALL]*hadron->v_w[vALL]*v_cosP[vALL] - hadron->v_v[vALL]*v_sinP[vALL]));

    hadron->v_v[vALL] = (hadron->v_v[vALL]*v_cosT[vALL] 
			+ (v_sinT[vALL] / (sqrt(1.0-hadron->v_w[vALL]*hadron->v_w[vALL]))) * (hadron->v_v[vALL]*hadron->v_w[vALL]*v_cosP[vALL] + v_prev_u[vALL]*v_sinP[vALL]));

    hadron->v_w[vALL] = (hadron->v_w[vALL]*v_cosT[vALL] - sqrt(1.0-hadron->v_w[vALL]*hadron->v_w[vALL])*v_sinT[vALL]*v_cosP[vALL]);

  }
  else{
    hadron->v_v[vALL] = v_sinT[vALL] * v_sinP[vALL];

    if(hadron->v_w[vALL] > 0){		// Si direction parallère à l'axe z
      hadron->v_u[vALL] = v_sinT[vALL] * v_cosP[vALL];
      hadron->v_w[vALL] = v_cosT[vALL];
    }
    else{				// Si direction antiparallère à l'axe z
      hadron->v_u[vALL] = -v_sinT[vALL] * v_cosP[vALL];
      hadron->v_w[vALL] = -v_cosT[vALL];
    }
  }



  // Si la norme dévie trop de 1, on renormalise
  ALIGNED_(64) VAR_COMPUTE v_norme[VLENGTH];
  v_norme[vALL] = sqrt(hadron->v_u[vALL]*hadron->v_u[vALL] + hadron->v_v[vALL]*hadron->v_v[vALL] + hadron->v_w[vALL]*hadron->v_w[vALL]);

//  if(fabs(v_norme[vALL]-1) > 1e-10){	

    hadron->v_u[vALL] = hadron->v_u[vALL] / v_norme[vALL];
    hadron->v_v[vALL] = hadron->v_v[vALL] / v_norme[vALL];
    hadron->v_w[vALL] = hadron->v_w[vALL] / v_norme[vALL];
//  }  


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
