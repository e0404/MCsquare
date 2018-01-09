/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_scoring.h"

DATA_Scoring Init_Scoring(DATA_config *config, int Nbr_voxels, int init_dose){

  DATA_Scoring scoring;
  scoring.energy = (VAR_SCORING*)calloc(Nbr_voxels, sizeof(VAR_SCORING));

  if(init_dose > 0) scoring.dose = (VAR_SCORING*)calloc(Nbr_voxels, sizeof(VAR_SCORING));
  else scoring.dose = NULL;

  if(config->Score_PromptGammas == 1){
    scoring.PG_particles = (VAR_SCORING*)calloc(Nbr_voxels, sizeof(VAR_SCORING));
    scoring.PG_spectrum = (VAR_SCORING*)calloc(config->PG_Spectrum_NumBin, sizeof(VAR_SCORING));
  }
  else{
    scoring.PG_particles = NULL;
    scoring.PG_spectrum = NULL;
  }

  if(config->Score_LET == 1){
    scoring.LET = (VAR_SCORING*)calloc(Nbr_voxels, sizeof(VAR_SCORING));
    scoring.LET_denominator = (VAR_SCORING*)calloc(Nbr_voxels, sizeof(VAR_SCORING));
  }
  else{
    scoring.LET = NULL;
    scoring.LET_denominator = NULL;
  }


  return scoring;

}


void Energy_Scoring(DATA_Scoring *scoring, int index, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE SPR){

  scoring->energy[index] += multiplicity * dE / SPR;

}


void LET_Scoring(DATA_Scoring *scoring, int index, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE step, VAR_COMPUTE stop_pow, DATA_config *config){

  if(config->LET_Calculation_Method == 0){			// 0 = DepositedEnergy
    scoring->LET[index] += multiplicity * dE * dE / step;
  }
  else{				
    scoring->LET[index] += multiplicity * dE * stop_pow;
  }

  scoring->LET_denominator[index] += multiplicity * dE;

}


void PostProcess_Scoring(DATA_Scoring *scoring, DATA_CT *ct, Materials *material, VAR_COMPUTE normalization, unsigned long Nbr_simulated_primaries, DATA_config *config){

  double voxel_volume = ct->VoxelLength[0]*ct->VoxelLength[1]*ct->VoxelLength[2];
  int ii,j,k,index=0;


  if(config->DoseToWater == 1){
    for(ii=0; ii<ct->Nbr_voxels; ii++){
      scoring->energy[ii] /= material[ct->material[ii]].SPR;
    }
  }

  for(ii=0; ii<ct->Nbr_voxels; ii++){
    scoring->energy[ii] = scoring->energy[ii] * normalization / Nbr_simulated_primaries;
    scoring->energy[ii] *= (scoring->energy[ii] > 0);
    scoring->dose[ii] = scoring->energy[ii] / (voxel_volume*ct->density[ii]);
  }

  if(config->Dose_Segmentation != 0){
    for(ii=0; ii<ct->Nbr_voxels; ii++){
      scoring->dose[ii] *= (ct->density[ii] > config->Segmentation_Density_Threshold);
    }
  }

  if(config->Score_PromptGammas == 1){
    for(ii=0; ii<ct->Nbr_voxels; ii++){
      scoring->PG_particles[ii] = scoring->PG_particles[ii] * normalization / Nbr_simulated_primaries;
    }
  }

  if(config->Score_LET == 1){
    for(ii=0; ii<ct->Nbr_voxels; ii++){
      scoring->LET[ii] = (scoring->LET[ii] > 0) * scoring->LET[ii] / ((scoring->LET_denominator[ii] * 1e7) + FLT_EPSILON);
    }
  }

}


void Free_Scoring(DATA_Scoring *scoring){

  if(scoring->energy != NULL) free(scoring->energy);
  if(scoring->dose != NULL) free(scoring->dose);

  if(scoring->PG_particles != NULL) free(scoring->PG_particles);
  if(scoring->PG_spectrum != NULL) free(scoring->PG_spectrum);

  if(scoring->LET != NULL) free(scoring->LET);
  if(scoring->LET_denominator != NULL) free(scoring->LET_denominator);

}
