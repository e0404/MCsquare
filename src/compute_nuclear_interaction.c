/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_nuclear_interaction.h"

void proton_proton_cross_section(Hadron *hadron, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(v_density, 64);


  if(hadron->v_T[vALL] > 10*UMeV){
    v_result[vALL] = v_density[vALL] * (0.315*pow(hadron->v_T[vALL]/UMeV, -1.126) + 3.78e-6 * hadron->v_T[vALL]/UMeV) / 0.1119;
  }
  else{
    v_result[vALL] = 0.0;
  }
 
  return;
}


void total_Nuclear_cross_section(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result){

  __assume_aligned(&hadron->v_T, 64);

  __assume_aligned(v_material_label, 64);
  __assume_aligned(v_density, 64);
  __assume_aligned(v_result, 64);

  int i;

  ALIGNED_(64) int v_index[VLENGTH];
//  v_index[vALL] = (int)floor(hadron->v_T[vALL]/(UMeV*INTERP_BIN));
  ALIGNED_(64) VAR_COMPUTE v_T[VLENGTH];
  v_T[vALL] = hadron->v_T[vALL] / (UMeV*INTERP_BIN);
  v_index[vALL] = (int)floor(v_T[vALL]);

  ALIGNED_(64) VAR_COMPUTE v_T1[VLENGTH];
  v_T1[vALL] = v_index[vALL] * UMeV * INTERP_BIN;

  ALIGNED_(64) VAR_COMPUTE v_T2[VLENGTH];
  v_T2[vALL] = (v_index[vALL]+1) * UMeV * INTERP_BIN;

  ALIGNED_(64) VAR_COMPUTE v_Cross_section1[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Cross_section1[i] = (VAR_COMPUTE)material[v_material_label[i]].Interp_Total_Nuclear_Cross_Section[v_index[i]];
  }

  ALIGNED_(64) VAR_COMPUTE v_Cross_section2[VLENGTH];
  for(i=0; i<VLENGTH; i++){
    v_Cross_section2[i] = (VAR_COMPUTE)material[v_material_label[i]].Interp_Total_Nuclear_Cross_Section[v_index[i]+1];
  }

  vec_Linear_Interpolation(hadron->v_T, v_T1, v_T2, v_Cross_section1, v_Cross_section2, v_result);

  v_result[vALL] = v_result[vALL] * v_density[vALL];


  return;
}


VAR_COMPUTE Compute_Nuclear_interaction(int hadron_index, Hadron *hadron, Materials *material, int material_label, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, int scoring_index, DATA_Scoring *scoring, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  VAR_COMPUTE rnd, dE;
  int index;

  VAR_COMPUTE Energy = hadron->v_T[hadron_index]/UMeV;

  // Material = Hydrogen
  if(material[material_label].Nuclear_data_type == Proton_Proton){  
    if(Energy > 10.0) return Compute_Elastic_PP(hadron_index, hadron, secondary_hadron, Nbr_secondaries, RNG_Stream, config);
    else return 0.0;
  }

  // Pure material (no mixture)
  else if(material[material_label].Nuclear_data_type == ICRU){  

    if(Energy > 7.0){
      rnd = single_rand_uniform(RNG_Stream);
 
      index = (int)floor(Energy / INTERP_BIN);
      VAR_COMPUTE Total_cross_section = Linear_Interpolation(	Energy, 
								(VAR_COMPUTE)index * INTERP_BIN, 
								(VAR_COMPUTE)(index+1) * INTERP_BIN, 
								(VAR_COMPUTE)material[material_label].Interp_Total_Nuclear_Cross_Section[index], 
								(VAR_COMPUTE)material[material_label].Interp_Total_Nuclear_Cross_Section[index+1]);

      index = Binary_Search(Energy, material[material_label].Elastic_Energy_List, material[material_label].Nbr_Elastic_Energy);
      VAR_COMPUTE Elastic_cross_section = Linear_Interpolation(	Energy, 
								material[material_label].Elastic_Energy_List[index], 
								material[material_label].Elastic_Energy_List[index+1], 
								material[material_label].Nuclear_Elastic[index].Cross_section, 
								material[material_label].Nuclear_Elastic[index+1].Cross_section);
      Elastic_cross_section = Elastic_cross_section / Total_cross_section;


      // Nuclear elastic interaction
      if(rnd <= Elastic_cross_section){
	return Compute_Elastic_ICRU(hadron_index, hadron, &material[material_label], RNG_Stream);
      }

      // Nuclear inelastic interaction
      else{
	dE = 0.0;
	index = Binary_Search(Energy, material[material_label].Inelastic_Energy_List, material[material_label].Nbr_Inelastic_Energy);
	dE += Compute_Nuclear_Inelastic_proton(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[material_label], index, RNG_Stream, config);
	dE += Compute_Nuclear_Inelastic_deuteron(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[material_label], index, RNG_Stream, config);
	dE += Compute_Nuclear_Inelastic_alpha(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[material_label], index, RNG_Stream, config);
	dE += Compute_Nuclear_Inelastic_recoils(hadron->v_T[hadron_index], &material[material_label], index);
	if(config->Score_PromptGammas == 1){
	  Compute_PromptGamma(hadron_index, hadron, &material[material_label], scoring_index, scoring, RNG_Stream, config);
	}
	hadron->v_type[hadron_index] = Unknown;
	return dE; 
      }

    }

    // Energy < 7 MeV (no interaction)
    else return 0.0;
  }

  // Material = mixture of several component
  else if(material[material_label].Nuclear_data_type == Mixture){  

    int i, target_label;

    index = (int)floor(Energy / INTERP_BIN);
    VAR_COMPUTE Total_cross_section = Linear_Interpolation(	Energy, 
								(VAR_COMPUTE)index * INTERP_BIN, 
								(VAR_COMPUTE)(index+1) * INTERP_BIN, 
								(VAR_COMPUTE)material[material_label].Interp_Total_Nuclear_Cross_Section[index], 
								(VAR_COMPUTE)material[material_label].Interp_Total_Nuclear_Cross_Section[index+1]);

    VAR_COMPUTE cross_section = 0.0;
    rnd = single_rand_uniform(RNG_Stream) * Total_cross_section;

    for(i=0; i < material[material_label].NbrComponents; i++){
      target_label = material[material_label].Mixture_Components_label[i];

      if(material[target_label].Nuclear_data_type == Proton_Proton){
	if(Energy > 10.0) cross_section += material[material_label].Mixture_Components_fraction[i] * (0.315*pow(Energy, -1.126) + 3.78e-6 * Energy)/0.1119;
	if(rnd <= cross_section){
	  return Compute_Elastic_PP(hadron_index, hadron, secondary_hadron, Nbr_secondaries, RNG_Stream, config);
	}
      }

      else if(material[target_label].Nuclear_data_type == ICRU){
	if(Energy > 7.0){
	  index = Binary_Search(Energy, material[target_label].Elastic_Energy_List, material[target_label].Nbr_Elastic_Energy);
	  cross_section += material[material_label].Mixture_Components_fraction[i] * Linear_Interpolation(	Energy, 
														material[target_label].Elastic_Energy_List[index], 
														material[target_label].Elastic_Energy_List[index+1], 
														material[target_label].Nuclear_Elastic[index].Cross_section, 
														material[target_label].Nuclear_Elastic[index+1].Cross_section);
	  if(rnd <= cross_section){
	    return Compute_Elastic_ICRU(hadron_index, hadron, &material[target_label], RNG_Stream);
	  }

      	  index = Binary_Search(Energy, material[target_label].Inelastic_Energy_List, material[target_label].Nbr_Inelastic_Energy);
      	  cross_section += material[material_label].Mixture_Components_fraction[i] * Linear_Interpolation(	Energy, 
														material[target_label].Inelastic_Energy_List[index], 
														material[target_label].Inelastic_Energy_List[index+1], 
														material[target_label].Nuclear_Inelastic[index].Cross_section, 
														material[target_label].Nuclear_Inelastic[index+1].Cross_section);
	  if(rnd <= cross_section){
	    dE = 0.0;
	    dE += Compute_Nuclear_Inelastic_proton(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[target_label], index, RNG_Stream, config);
	    dE += Compute_Nuclear_Inelastic_deuteron(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[target_label], index, RNG_Stream, config);
	    dE += Compute_Nuclear_Inelastic_alpha(hadron_index, hadron, secondary_hadron, Nbr_secondaries, &material[target_label], index, RNG_Stream, config);
	    dE += Compute_Nuclear_Inelastic_recoils(hadron->v_T[hadron_index], &material[target_label], index);
	    if(config->Score_PromptGammas == 1){
	      Compute_PromptGamma(hadron_index, hadron, &material[target_label], scoring_index, scoring, RNG_Stream, config);
	    }
	    hadron->v_type[hadron_index] = Unknown;

	    return dE;
	  }

	} 	// end if 7 MeV
      } 	// end if ICRU
    } 		// end for

    return 0.0;
  }		// end if mixture
 
  else return 0.0;

}


VAR_COMPUTE Compute_Elastic_PP(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  VAR_COMPUTE cos_theta_CM = 2.0*single_rand_uniform(RNG_Stream) - 1.0;
  VAR_COMPUTE dE = hadron->v_T[hadron_index] * (1-cos_theta_CM) / 2.0;

  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  v_theta[vALL] = 0.0;
  v_theta[hadron_index] = acos( (cos_theta_CM+1.0)/sqrt((cos_theta_CM+1.0)*(cos_theta_CM+1.0) + (1 - hadron->v_T[hadron_index]/(hadron->v_T[hadron_index]+2*MC2_PRO))*(1.0-cos_theta_CM*cos_theta_CM)) );

  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  v_phi[vALL] = 0.0;
  v_phi[hadron_index] = 2*M_PI*single_rand_uniform(RNG_Stream);
      	  
  // Secondary proton :
  secondary_hadron[*Nbr_secondaries].T = dE;
  secondary_hadron[*Nbr_secondaries].x = hadron->v_x[hadron_index];
  secondary_hadron[*Nbr_secondaries].y = hadron->v_y[hadron_index];
  secondary_hadron[*Nbr_secondaries].z = hadron->v_z[hadron_index];
  secondary_hadron[*Nbr_secondaries].u = hadron->v_u[hadron_index];
  secondary_hadron[*Nbr_secondaries].v = hadron->v_v[hadron_index];
  secondary_hadron[*Nbr_secondaries].w = hadron->v_w[hadron_index];

  secondary_hadron[*Nbr_secondaries].M = hadron->v_M[hadron_index];
  secondary_hadron[*Nbr_secondaries].charge = 1;
  secondary_hadron[*Nbr_secondaries].mass = 1;
  secondary_hadron[*Nbr_secondaries].type = Proton;

  VAR_COMPUTE secondary_theta = acos( (1.0-cos_theta_CM)/sqrt((1.0-cos_theta_CM)*(1.0-cos_theta_CM) + (1 - hadron->v_T[hadron_index]/(hadron->v_T[hadron_index]+2*MC2_PRO))*(1.0-cos_theta_CM*cos_theta_CM)) );
  VAR_COMPUTE secondary_phi;
  if(v_phi[hadron_index] > M_PI) secondary_phi = v_phi[hadron_index] - M_PI;
  else secondary_phi = v_phi[hadron_index] + M_PI;

  Update_direction(hadron, v_theta, v_phi);
  Update_buffer_direction(&secondary_hadron[*Nbr_secondaries], secondary_theta, secondary_phi);

  hadron->v_T[hadron_index] = hadron->v_T[hadron_index] - dE;

  if(secondary_hadron[*Nbr_secondaries].T < config->Ecut_Pro * UMeV) return secondary_hadron[*Nbr_secondaries].T;

  *Nbr_secondaries += 1;

  return 0.0;
}


VAR_COMPUTE Compute_Elastic_ICRU(int hadron_index, Hadron *hadron, Materials *material, VSLStreamStatePtr RNG_Stream){

  int index, i;
  VAR_COMPUTE cos_theta_CM, rnd;
  VAR_DATA diff_cross_section[36];    

  index = Binary_Search(hadron->v_T[hadron_index]/UMeV, material->Elastic_Energy_List, material->Nbr_Elastic_Energy);

  for(i=0; i<36; i++){
    diff_cross_section[i] = (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Elastic_Energy_List[index], 
								material->Elastic_Energy_List[index+1], 
								material->Nuclear_Elastic[index].Diff_Cross_section[i], 
								material->Nuclear_Elastic[index+1].Diff_Cross_section[i]);
  }

  // échantillonnage de l'angle
  // la section efficace a déjà été convertie en densité de proba cumulative (data_nuclear.c)
  rnd = single_rand_uniform(RNG_Stream) * diff_cross_section[35];
  i = Sequential_Search(rnd, diff_cross_section, 36) + 1;

  if(i == 0) cos_theta_CM = cos(5*M_PI/180) + single_rand_uniform(RNG_Stream) * (cos(7.5*M_PI/180)-cos(5*M_PI/180));
  else if(i == 35) cos_theta_CM = cos(175.5*M_PI/180) + single_rand_uniform(RNG_Stream) * (cos(180*M_PI/180)-cos(177.5*M_PI/180));
  else cos_theta_CM = cos(((i+1)*5 - 2.5)*M_PI/180) + single_rand_uniform(RNG_Stream) * (cos(((i+1)*5 + 2.5)*M_PI/180)-cos(((i+1)*5 - 2.5)*M_PI/180));
    
  VAR_COMPUTE beta2_CM = (hadron->v_T[hadron_index]*(hadron->v_T[hadron_index]+2*MC2_PRO)) / ((hadron->v_T[hadron_index] + MC2_PRO + material->A*Uamu)*(hadron->v_T[hadron_index] + MC2_PRO + material->A*Uamu));
  VAR_COMPUTE gamma2_CM = 1.0/(1-beta2_CM);
  VAR_COMPUTE tau = sqrt( (MC2_PRO/(material->A*Uamu)) * (MC2_PRO/(material->A*Uamu)) * (1-beta2_CM) + beta2_CM );

  ALIGNED_(64) VAR_COMPUTE v_theta[VLENGTH];
  v_theta[vALL] = 0.0;
  v_theta[hadron_index] = acos( (cos_theta_CM+tau)/sqrt( (cos_theta_CM+tau)*(cos_theta_CM+tau) + (1 - cos_theta_CM*cos_theta_CM)/gamma2_CM ) );

  ALIGNED_(64) VAR_COMPUTE v_phi[VLENGTH];
  v_phi[vALL] = 0.0;
  v_phi[hadron_index] = 2*M_PI*single_rand_uniform(RNG_Stream);

  Update_direction(hadron, v_theta, v_phi);

  VAR_COMPUTE dE = gamma2_CM*beta2_CM*material->A*Uamu * (1-cos_theta_CM);
  hadron->v_T[hadron_index] -= dE;

  return dE;
}


VAR_COMPUTE Compute_Nuclear_Inelastic_recoils(VAR_COMPUTE Hadron_T, Materials *material, int index){

  return Hadron_T * Linear_Interpolation(	Hadron_T/UMeV, 
						material->Inelastic_Energy_List[index], 
						material->Inelastic_Energy_List[index+1], 
						material->Nuclear_Inelastic[index].Energy_fraction_recoils, 
						material->Nuclear_Inelastic[index+1].Energy_fraction_recoils);

}


VAR_COMPUTE Compute_Nuclear_Inelastic_proton(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  if(material->Nuclear_Inelastic[index].Proton_Mult == 0.0 || material->Nuclear_Inelastic[index+1].Proton_Mult == 0.0) return 0.0;


  // Interpolation de la multiplicité en fonction de l'énergie du proton incident
  secondary_hadron[*Nbr_secondaries].M = hadron->v_M[hadron_index] * Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
												material->Inelastic_Energy_List[index], 
												material->Inelastic_Energy_List[index+1], 
												material->Nuclear_Inelastic[index].Proton_Mult, 
												material->Nuclear_Inelastic[index+1].Proton_Mult);

  // Interpolation de la section eff diff (intégrée sur l'angle) en fonction de l'énergie du proton incident
  // + conversion en proba cumulative

  VAR_DATA *diff_cross_section = (VAR_DATA*) malloc(material->Nuclear_Inelastic[index].P_Nbr_data * sizeof(VAR_DATA));

  diff_cross_section[0] = (VAR_DATA) Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index], 
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].P_D_Cross_section[0], 
								material->Nuclear_Inelastic[index+1].P_D_Cross_section[0]);
  int i;
  for(i=1; i<material->Nuclear_Inelastic[index].P_Nbr_data; i++){
    diff_cross_section[i] = diff_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
											material->Inelastic_Energy_List[index], 
											material->Inelastic_Energy_List[index+1], 
											material->Nuclear_Inelastic[index].P_D_Cross_section[i], 
											material->Nuclear_Inelastic[index+1].P_D_Cross_section[i]);
  }
    

  // échantillonnage de l'énergie de la particule secondaire
  VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream) * diff_cross_section[material->Nuclear_Inelastic[index].P_Nbr_data - 1];
  int secondary_index = Binary_Search(rnd, diff_cross_section, material->Nuclear_Inelastic[index].P_Nbr_data);

  secondary_hadron[*Nbr_secondaries].T =  Linear_Interpolation(	rnd, 
								diff_cross_section[secondary_index], 
								diff_cross_section[secondary_index+1],
								material->Nuclear_Inelastic[index].P_Energy_list[secondary_index], 
								material->Nuclear_Inelastic[index].P_Energy_list[secondary_index+1]);
  free(diff_cross_section);

  // scaling de l'énergie
  secondary_hadron[*Nbr_secondaries].T = ((hadron->v_T[hadron_index]/UMeV) / material->Inelastic_Energy_List[index]) * UMeV * secondary_hadron[*Nbr_secondaries].T;
  if(secondary_hadron[*Nbr_secondaries].T < config->Ecut_Pro * UMeV) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];  
  if(config->Simulate_Secondary_Protons == 0) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];  // divise par M car dE sera remultiplié par M plus tard

  secondary_hadron[*Nbr_secondaries].type = Proton;
  secondary_hadron[*Nbr_secondaries].charge = 1;
  secondary_hadron[*Nbr_secondaries].mass = 1;
  secondary_hadron[*Nbr_secondaries].x = hadron->v_x[hadron_index];
  secondary_hadron[*Nbr_secondaries].y = hadron->v_y[hadron_index];
  secondary_hadron[*Nbr_secondaries].z = hadron->v_z[hadron_index];
  secondary_hadron[*Nbr_secondaries].u = hadron->v_u[hadron_index];
  secondary_hadron[*Nbr_secondaries].v = hadron->v_v[hadron_index];
  secondary_hadron[*Nbr_secondaries].w = hadron->v_w[hadron_index];

  VAR_COMPUTE phi = 2*M_PI*single_rand_uniform(RNG_Stream);


  // interpolation de la section eff double diff en fonction de l'énergie du proton incident
  // cette interpolation est effectuée pour les deux index voisin de l'énergie échantillonnée pour la particule secondaire

  VAR_DATA dd_cross_section1[13], dd_cross_section2[13];

  for(i=0; i<13; i++){
    dd_cross_section1[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].P_DD_Cross_section[13*secondary_index+i], 
								material->Nuclear_Inelastic[index+1].P_DD_Cross_section[13*secondary_index+i]);

    dd_cross_section2[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].P_DD_Cross_section[13*(secondary_index+1)+i], 
								material->Nuclear_Inelastic[index+1].P_DD_Cross_section[13*(secondary_index+1)+i]);

  }


  // interpolation de la section eff double diff en fonction de l'énergie de la particule secondaire
  // + conversion en fonction de proba cumulative

  VAR_DATA dd_cross_section[13];

  dd_cross_section[0] = (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
							material->Nuclear_Inelastic[index].P_Energy_list[secondary_index],
							material->Nuclear_Inelastic[index].P_Energy_list[secondary_index+1], 
							dd_cross_section1[0],
							dd_cross_section2[0]);

  for(i=1; i<13; i++){
    dd_cross_section[i] = dd_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
											material->Nuclear_Inelastic[index].P_Energy_list[secondary_index],
											material->Nuclear_Inelastic[index].P_Energy_list[secondary_index+1], 
											dd_cross_section1[i],
											dd_cross_section2[i]);
  }

  // Echantillonage de l'angle theta d'émission de la particule secondaire
  rnd = single_rand_uniform(RNG_Stream) * dd_cross_section[12];
  int angle_index = Binary_Search(rnd, dd_cross_section, 13)+1;

  static const double ICRU_angles[13] = { 0, 10, 20, 30, 40, 50, 60, 70, 90, 110, 130, 150, 180 };

  VAR_COMPUTE theta;
  if(angle_index == 0) theta = acos(1 + single_rand_uniform(RNG_Stream) * (cos(5*M_PI/180) - 1) );
  else if(angle_index == 12) theta = acos( cos(165*M_PI/180) + single_rand_uniform(RNG_Stream) * (-1.0 - cos(165*M_PI/180)));
  else theta = acos( cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180)) + single_rand_uniform(RNG_Stream) * (cos((ICRU_angles[angle_index]+ICRU_angles[angle_index+1])*M_PI/(2*180)) - cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180))) );

  Update_buffer_direction(&secondary_hadron[*Nbr_secondaries], theta, phi);

  *Nbr_secondaries += 1;

  return 0.0;
}


VAR_COMPUTE Compute_Nuclear_Inelastic_deuteron(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  if(material->Nuclear_Inelastic[index].Deuteron_Mult == 0.0 || material->Nuclear_Inelastic[index+1].Deuteron_Mult == 0.0) return 0.0;


  // Interpolation de la multiplicité en fonction de l'énergie du proton incident
  secondary_hadron[*Nbr_secondaries].M = hadron->v_M[hadron_index] * Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
												material->Inelastic_Energy_List[index], 
												material->Inelastic_Energy_List[index+1], 
												material->Nuclear_Inelastic[index].Deuteron_Mult, 
												material->Nuclear_Inelastic[index+1].Deuteron_Mult);

  // Interpolation de la section eff diff (intégrée sur l'angle) en fonction de l'énergie du proton incident
  // + conversion en proba cumulative

  VAR_DATA *diff_cross_section = (VAR_DATA*) malloc(material->Nuclear_Inelastic[index].D_Nbr_data * sizeof(VAR_DATA));

  diff_cross_section[0] = (VAR_DATA) Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index], 
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].D_D_Cross_section[0], 
								material->Nuclear_Inelastic[index+1].D_D_Cross_section[0]);
  int i;
  for(i=1; i<material->Nuclear_Inelastic[index].D_Nbr_data; i++){
    diff_cross_section[i] = diff_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
											material->Inelastic_Energy_List[index], 
											material->Inelastic_Energy_List[index+1], 
											material->Nuclear_Inelastic[index].D_D_Cross_section[i], 
											material->Nuclear_Inelastic[index+1].D_D_Cross_section[i]);
  }
    

  // échantillonnage de l'énergie de la particule secondaire
  VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream) * diff_cross_section[material->Nuclear_Inelastic[index].D_Nbr_data - 1];
  int secondary_index = Binary_Search(rnd, diff_cross_section, material->Nuclear_Inelastic[index].D_Nbr_data);

  secondary_hadron[*Nbr_secondaries].T =  Linear_Interpolation(	rnd, 
								diff_cross_section[secondary_index], 
								diff_cross_section[secondary_index+1],
								material->Nuclear_Inelastic[index].D_Energy_list[secondary_index], 
								material->Nuclear_Inelastic[index].D_Energy_list[secondary_index+1]);
  free(diff_cross_section);

  // scaling de l'énergie
  secondary_hadron[*Nbr_secondaries].T = ((hadron->v_T[hadron_index]/UMeV) / material->Inelastic_Energy_List[index]) * UMeV * secondary_hadron[*Nbr_secondaries].T;
  if(secondary_hadron[*Nbr_secondaries].T < config->Ecut_Pro * UMeV) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];
  if(config->Simulate_Secondary_Deuterons == 0) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];  // divise par M car dE sera remultiplié par M plus tard

  secondary_hadron[*Nbr_secondaries].type = Deuteron;
  secondary_hadron[*Nbr_secondaries].charge = 1;
  secondary_hadron[*Nbr_secondaries].mass = 2;
  secondary_hadron[*Nbr_secondaries].x = hadron->v_x[hadron_index];
  secondary_hadron[*Nbr_secondaries].y = hadron->v_y[hadron_index];
  secondary_hadron[*Nbr_secondaries].z = hadron->v_z[hadron_index];
  secondary_hadron[*Nbr_secondaries].u = hadron->v_u[hadron_index];
  secondary_hadron[*Nbr_secondaries].v = hadron->v_v[hadron_index];
  secondary_hadron[*Nbr_secondaries].w = hadron->v_w[hadron_index];

  VAR_COMPUTE phi = 2*M_PI*single_rand_uniform(RNG_Stream);



  // interpolation de la section eff double diff en fonction de l'énergie du proton incident
  // cette interpolation est effectuée pour les deux index voisin de l'énergie échantillonnée pour la particule secondaire

  VAR_DATA dd_cross_section1[13], dd_cross_section2[13];

  for(i=0; i<13; i++){
    dd_cross_section1[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].D_DD_Cross_section[13*secondary_index+i], 
								material->Nuclear_Inelastic[index+1].D_DD_Cross_section[13*secondary_index+i]);

    dd_cross_section2[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].D_DD_Cross_section[13*(secondary_index+1)+i], 
								material->Nuclear_Inelastic[index+1].D_DD_Cross_section[13*(secondary_index+1)+i]);

  }


  // interpolation de la section eff double diff en fonction de l'énergie de la particule secondaire
  // + conversion en fonction de proba cumulative

  VAR_DATA dd_cross_section[13];

  dd_cross_section[0] = (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
							material->Nuclear_Inelastic[index].D_Energy_list[secondary_index],
							material->Nuclear_Inelastic[index].D_Energy_list[secondary_index+1], 
							dd_cross_section1[0],
							dd_cross_section2[0]);

  for(i=1; i<13; i++){
    dd_cross_section[i] = dd_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
											material->Nuclear_Inelastic[index].D_Energy_list[secondary_index],
											material->Nuclear_Inelastic[index].D_Energy_list[secondary_index+1], 
											dd_cross_section1[i],
											dd_cross_section2[i]);
  }

  // Echantillonage de l'angle theta d'émission de la particule secondaire
  rnd = single_rand_uniform(RNG_Stream) * dd_cross_section[12];
  int angle_index = Binary_Search(rnd, dd_cross_section, 13)+1;

  static const double ICRU_angles[13] = { 0, 10, 20, 30, 40, 50, 60, 70, 90, 110, 130, 150, 180 };

  VAR_COMPUTE theta;
  if(angle_index == 0) theta = acos(1 + single_rand_uniform(RNG_Stream) * (cos(5*M_PI/180) - 1) );
  else if(angle_index == 12) theta = acos( cos(165*M_PI/180) + single_rand_uniform(RNG_Stream) * (-1.0 - cos(165*M_PI/180)));
  else theta = acos( cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180)) + single_rand_uniform(RNG_Stream) * (cos((ICRU_angles[angle_index]+ICRU_angles[angle_index+1])*M_PI/(2*180)) - cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180))) );

  Update_buffer_direction(&secondary_hadron[*Nbr_secondaries], theta, phi);


  *Nbr_secondaries += 1;

  return 0.0;
}


VAR_COMPUTE Compute_Nuclear_Inelastic_alpha(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  if(material->Nuclear_Inelastic[index].Alpha_Mult == 0.0 || material->Nuclear_Inelastic[index+1].Alpha_Mult == 0.0) return 0.0;


  // Interpolation de la multiplicité en fonction de l'énergie du proton incident
  secondary_hadron[*Nbr_secondaries].M = hadron->v_M[hadron_index] * Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
												material->Inelastic_Energy_List[index], 
												material->Inelastic_Energy_List[index+1], 
												material->Nuclear_Inelastic[index].Alpha_Mult, 
												material->Nuclear_Inelastic[index+1].Alpha_Mult);

  // Interpolation de la section eff diff (intégrée sur l'angle) en fonction de l'énergie du proton incident
  // + conversion en proba cumulative

  VAR_DATA *diff_cross_section = (VAR_DATA*) malloc(material->Nuclear_Inelastic[index].A_Nbr_data * sizeof(VAR_DATA));

  diff_cross_section[0] = (VAR_DATA) Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index], 
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].A_D_Cross_section[0], 
								material->Nuclear_Inelastic[index+1].A_D_Cross_section[0]);
  int i;
  for(i=1; i<material->Nuclear_Inelastic[index].A_Nbr_data; i++){
    diff_cross_section[i] = diff_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
											material->Inelastic_Energy_List[index], 
											material->Inelastic_Energy_List[index+1], 
											material->Nuclear_Inelastic[index].A_D_Cross_section[i], 
											material->Nuclear_Inelastic[index+1].A_D_Cross_section[i]);
  }
    

  // échantillonnage de l'énergie de la particule secondaire
  VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream) * diff_cross_section[material->Nuclear_Inelastic[index].A_Nbr_data - 1];
  int secondary_index = Binary_Search(rnd, diff_cross_section, material->Nuclear_Inelastic[index].A_Nbr_data);

  secondary_hadron[*Nbr_secondaries].T =  Linear_Interpolation(	rnd, 
								diff_cross_section[secondary_index], 
								diff_cross_section[secondary_index+1],
								material->Nuclear_Inelastic[index].A_Energy_list[secondary_index], 
								material->Nuclear_Inelastic[index].A_Energy_list[secondary_index+1]);
  free(diff_cross_section);

  // scaling de l'énergie
  secondary_hadron[*Nbr_secondaries].T = ((hadron->v_T[hadron_index]/UMeV) / material->Inelastic_Energy_List[index]) * UMeV * secondary_hadron[*Nbr_secondaries].T;
  if(secondary_hadron[*Nbr_secondaries].T < config->Ecut_Pro * UMeV) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];
  if(config->Simulate_Secondary_Alphas == 0) return secondary_hadron[*Nbr_secondaries].M * secondary_hadron[*Nbr_secondaries].T / hadron->v_M[hadron_index];  // divise par M car dE sera remultiplié par M plus tard

  secondary_hadron[*Nbr_secondaries].type = Alpha;
  secondary_hadron[*Nbr_secondaries].charge = 2;
  secondary_hadron[*Nbr_secondaries].mass = 4;
  secondary_hadron[*Nbr_secondaries].x = hadron->v_x[hadron_index];
  secondary_hadron[*Nbr_secondaries].y = hadron->v_y[hadron_index];
  secondary_hadron[*Nbr_secondaries].z = hadron->v_z[hadron_index];
  secondary_hadron[*Nbr_secondaries].u = hadron->v_u[hadron_index];
  secondary_hadron[*Nbr_secondaries].v = hadron->v_v[hadron_index];
  secondary_hadron[*Nbr_secondaries].w = hadron->v_w[hadron_index];

  VAR_COMPUTE phi = 2*M_PI*single_rand_uniform(RNG_Stream);

  // interpolation de la section eff double diff en fonction de l'énergie du proton incident
  // cette interpolation est effectuée pour les deux index voisin de l'énergie échantillonnée pour la particule secondaire

  VAR_DATA dd_cross_section1[13], dd_cross_section2[13];

  for(i=0; i<13; i++){
    dd_cross_section1[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].A_DD_Cross_section[13*secondary_index+i], 
								material->Nuclear_Inelastic[index+1].A_DD_Cross_section[13*secondary_index+i]);

    dd_cross_section2[i] =  (VAR_DATA)Linear_Interpolation(	hadron->v_T[hadron_index]/UMeV, 
								material->Inelastic_Energy_List[index],
								material->Inelastic_Energy_List[index+1], 
								material->Nuclear_Inelastic[index].A_DD_Cross_section[13*(secondary_index+1)+i], 
								material->Nuclear_Inelastic[index+1].A_DD_Cross_section[13*(secondary_index+1)+i]);

  }


  // interpolation de la section eff double diff en fonction de l'énergie de la particule secondaire
  // + conversion en fonction de proba cumulative

  VAR_DATA dd_cross_section[13];

  dd_cross_section[0] = (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
							material->Nuclear_Inelastic[index].A_Energy_list[secondary_index],
							material->Nuclear_Inelastic[index].A_Energy_list[secondary_index+1], 
							dd_cross_section1[0],
							dd_cross_section2[0]);

  for(i=1; i<13; i++){
    dd_cross_section[i] = dd_cross_section[i-1] + (VAR_DATA)Linear_Interpolation(	secondary_hadron[*Nbr_secondaries].T, 
											material->Nuclear_Inelastic[index].A_Energy_list[secondary_index],
											material->Nuclear_Inelastic[index].A_Energy_list[secondary_index+1], 
											dd_cross_section1[i],
											dd_cross_section2[i]);
  }

  // Echantillonage de l'angle theta d'émission de la particule secondaire
  rnd = single_rand_uniform(RNG_Stream) * dd_cross_section[12];
  int angle_index = Binary_Search(rnd, dd_cross_section, 13)+1;

  static const double ICRU_angles[13] = { 0, 10, 20, 30, 40, 50, 60, 70, 90, 110, 130, 150, 180 };

  VAR_COMPUTE theta;
  if(angle_index == 0) theta = acos(1 + single_rand_uniform(RNG_Stream) * (cos(5*M_PI/180) - 1) );
  else if(angle_index == 12) theta = acos( cos(165*M_PI/180) + single_rand_uniform(RNG_Stream) * (-1.0 - cos(165*M_PI/180)));
  else theta = acos( cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180)) + single_rand_uniform(RNG_Stream) * (cos((ICRU_angles[angle_index]+ICRU_angles[angle_index+1])*M_PI/(2*180)) - cos((ICRU_angles[angle_index-1]+ICRU_angles[angle_index])*M_PI/(2*180))) );

  Update_buffer_direction(&secondary_hadron[*Nbr_secondaries], theta, phi);

  *Nbr_secondaries += 1;

  return 0.0;
}


void Compute_PromptGamma(int hadron_index, Hadron *hadron, Materials *material, int scoring_index, DATA_Scoring *scoring, VSLStreamStatePtr RNG_Stream, DATA_config *config){

  VAR_COMPUTE Energy = hadron->v_T[hadron_index]/UMeV;
  int index, interp_index;
  index = Binary_Search(Energy, material->PG_Energy_List, material->Nbr_PG_Energy);

  if( index < 0 || (index+1) > material->Nbr_PG_Energy ) return;
  if(material->PG_data[index].Multiplicity == 0.0 || material->PG_data[index+1].Multiplicity == 0.0) return;


  // Interpolation de la multiplicité en fonction de l'énergie du proton incident
  VAR_COMPUTE M = hadron->v_M[hadron_index] * Linear_Interpolation(	Energy, 
									material->PG_Energy_List[index], 
									material->PG_Energy_List[index+1], 
									material->PG_data[index].Multiplicity, 
									material->PG_data[index+1].Multiplicity);



  // Nearest-neighbor interpolation of diferential cross section
  if((Energy - material->PG_Energy_List[index]) < (material->PG_Energy_List[index+1] - Energy)) interp_index = index;
  else interp_index = index + 1;

  VAR_DATA *diff_cross_section = (VAR_DATA*) malloc(material->PG_data[interp_index].Nbr_data * sizeof(VAR_DATA));
  
  diff_cross_section[0] = material->PG_data[interp_index].Diff_Cross_section[0];
  int i;
  for(i=1; i<material->PG_data[interp_index].Nbr_data; i++){
    diff_cross_section[i] = diff_cross_section[i-1] + material->PG_data[interp_index].Diff_Cross_section[i];
  }
  
  // échantillonnage de l'énergie de la particule secondaire
  VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream) * diff_cross_section[material->PG_data[interp_index].Nbr_data - 1];
  int secondary_index = Binary_Search(rnd, diff_cross_section, material->PG_data[interp_index].Nbr_data);

  VAR_COMPUTE PG_energy =  Linear_Interpolation(rnd, 
						diff_cross_section[secondary_index], 
						diff_cross_section[secondary_index+1],
						material->PG_data[interp_index].Energy_list[secondary_index], 
						material->PG_data[interp_index].Energy_list[secondary_index+1]);
  free(diff_cross_section);


  if(PG_energy >= config->PG_LowEnergyCut && PG_energy <= config->PG_HighEnergyCut){
    scoring->PG_particles[scoring_index] += M;
    if(PG_energy >= config->PG_Spectrum_Binning*config->PG_Spectrum_NumBin) scoring->PG_spectrum[config->PG_Spectrum_NumBin-1] += M;
    else scoring->PG_spectrum[(int)floor(PG_energy / config->PG_Spectrum_Binning)] += M;
  }

  return;
}
