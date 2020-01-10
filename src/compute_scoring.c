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

DATA_Scoring Init_Scoring(DATA_config *config, DATA_CT *ct, int init_dose_squared){

  DATA_Scoring scoring;

  if(config->Independent_scoring_grid == 0){
    scoring.Nbr_voxels = ct->Nbr_voxels;
    scoring.Origin[0] = ct->Origin[0];
    scoring.Origin[1] = ct->Origin[1];
    scoring.Origin[2] = ct->Origin[2];
    scoring.Offset[0] = 0.0;
    scoring.Offset[1] = 0.0;
    scoring.Offset[2] = 0.0;
    scoring.Length[0] = ct->Length[0];
    scoring.Length[1] = ct->Length[1];
    scoring.Length[2] = ct->Length[2];
    scoring.Grid_end[0] = ct->Length[0];
    scoring.Grid_end[1] = ct->Length[1];
    scoring.Grid_end[2] = ct->Length[2];
    scoring.GridSize[0] = ct->GridSize[0];
    scoring.GridSize[1] = ct->GridSize[1];
    scoring.GridSize[2] = ct->GridSize[2];
    scoring.VoxelLength[0] = ct->VoxelLength[0];
    scoring.VoxelLength[1] = ct->VoxelLength[1];
    scoring.VoxelLength[2] = ct->VoxelLength[2];
  }
  else{
    scoring.Nbr_voxels = config->Scoring_grid_size[0] * config->Scoring_grid_size[1] * config->Scoring_grid_size[2];
    scoring.Origin[0] = config->Scoring_origin[0];
    scoring.Origin[1] = config->Scoring_origin[1];
    scoring.Origin[2] = config->Scoring_origin[2];
    scoring.Offset[0] = ct->Origin[0] - config->Scoring_origin[0];
    scoring.Offset[1] = ct->Origin[1] - config->Scoring_origin[1];
    scoring.Offset[2] = ct->Origin[2] - config->Scoring_origin[2];
    scoring.Length[0] = config->Scoring_grid_size[0] * config->Scoring_voxel_spacing[0];
    scoring.Length[1] = config->Scoring_grid_size[1] * config->Scoring_voxel_spacing[1];
    scoring.Length[2] = config->Scoring_grid_size[2] * config->Scoring_voxel_spacing[2];
    scoring.Grid_end[0] = scoring.Offset[0] + scoring.Length[0];
    scoring.Grid_end[1] = scoring.Offset[1] + scoring.Length[1];
    scoring.Grid_end[2] = scoring.Offset[2] + scoring.Length[2];
    scoring.GridSize[0] = config->Scoring_grid_size[0];
    scoring.GridSize[1] = config->Scoring_grid_size[1];
    scoring.GridSize[2] = config->Scoring_grid_size[2];
    scoring.VoxelLength[0] = config->Scoring_voxel_spacing[0];
    scoring.VoxelLength[1] = config->Scoring_voxel_spacing[1];
    scoring.VoxelLength[2] = config->Scoring_voxel_spacing[2];
  }

  scoring.dose = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));

  if(config->Score_Energy == 1) scoring.energy = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));
  else scoring.energy = NULL;

  if(init_dose_squared > 0) scoring.dose_squared = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));
  else scoring.dose_squared = NULL;

  if(config->Score_PromptGammas == 1){
    scoring.PG_particles = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));
    scoring.PG_spectrum = (VAR_SCORING*)calloc(config->PG_Spectrum_NumBin, sizeof(VAR_SCORING));
  }
  else{
    scoring.PG_particles = NULL;
    scoring.PG_spectrum = NULL;
  }

  if(config->Score_LET == 1){
    scoring.LET = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));
    scoring.LET_denominator = (VAR_SCORING*)calloc(scoring.Nbr_voxels, sizeof(VAR_SCORING));
  }
  else{
    scoring.LET = NULL;
    scoring.LET_denominator = NULL;
  }


  return scoring;

}


int get_scoring_index(DATA_Scoring *scoring, VAR_COMPUTE position_x, VAR_COMPUTE position_y, VAR_COMPUTE position_z){
  
  if(position_x < scoring->Offset[0] || position_y < scoring->Offset[1] || position_z < scoring->Offset[2] || position_x > scoring->Grid_end[0] || position_y > scoring->Grid_end[1] || position_z > scoring->Grid_end[2]) return -1;

   int index = 	(int)floor( (-position_x + scoring->Length[0]) / scoring->VoxelLength[0] ) 
			+ scoring->GridSize[0] * (int)floor( position_y / scoring->VoxelLength[1] ) 
			+ scoring->GridSize[0] * scoring->GridSize[1] * (int)floor( position_z / scoring->VoxelLength[2] );

  if(index > scoring->Nbr_voxels) return -1;

  return index;
  
}


void Energy_Scoring(DATA_Scoring *scoring, VAR_COMPUTE position_x, VAR_COMPUTE position_y, VAR_COMPUTE position_z, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE density, VAR_COMPUTE SPR, DATA_config *config){

  int index = get_scoring_index(scoring, position_x, position_y, position_z);
  if(index < 0) return;

  
  if(config->Dose_weighting_algorithm == 0) scoring->dose[index] += multiplicity * dE / (density * SPR); // Volume weighting
  else scoring->dose[index] += multiplicity * dE / SPR;  // Mass weighting
  // SPR is used when online dose to water conversion is enabled

  if(config->Score_Energy == 1) scoring->energy[index] += multiplicity * dE;

}


void LET_Scoring(DATA_Scoring *scoring, VAR_COMPUTE position_x, VAR_COMPUTE position_y, VAR_COMPUTE position_z, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE step, VAR_COMPUTE stop_pow, DATA_config *config){

  int index = get_scoring_index(scoring, position_x, position_y, position_z);
  if(index < 0) return;

  if(config->LET_Calculation_Method == 0){			// 0 = DepositedEnergy
    scoring->LET[index] += multiplicity * dE * dE / step;
  }
  else{				
    scoring->LET[index] += multiplicity * dE * stop_pow;
  }

  scoring->LET_denominator[index] += multiplicity * dE;

}


void PG_Scoring(DATA_Scoring *scoring, VAR_COMPUTE position_x, VAR_COMPUTE position_y, VAR_COMPUTE position_z, VAR_COMPUTE multiplicity, VAR_COMPUTE PG_energy, DATA_config *config){

  int index = get_scoring_index(scoring, position_x, position_y, position_z);
  if(index < 0) return;

  if(PG_energy >= config->PG_LowEnergyCut && PG_energy <= config->PG_HighEnergyCut){
    scoring->PG_particles[index] += multiplicity;
    if(PG_energy >= config->PG_Spectrum_Binning*config->PG_Spectrum_NumBin) scoring->PG_spectrum[config->PG_Spectrum_NumBin-1] += multiplicity;
    else scoring->PG_spectrum[(int)floor(PG_energy / config->PG_Spectrum_Binning)] += multiplicity;
  }

}


void PostProcess_Scoring(DATA_Scoring *scoring, DATA_CT *ct, Materials *material, VAR_COMPUTE normalization, unsigned long Nbr_simulated_primaries, DATA_config *config){

  double voxel_volume = scoring->VoxelLength[0]*scoring->VoxelLength[1]*scoring->VoxelLength[2];
  int ii,j,k,index=0;


  if(config->DoseToWater == 1){ // dose to water conversion by post-processing
    for(ii=0; ii<scoring->Nbr_voxels; ii++){
      scoring->dose[ii] /= material[ct->material[ii]].SPR;
    }
  }

  for(ii=0; ii<scoring->Nbr_voxels; ii++){
    scoring->dose[ii] = scoring->dose[ii] * normalization / Nbr_simulated_primaries;
    scoring->dose[ii] *= (scoring->dose[ii] > 0);
    scoring->dose[ii] = scoring->dose[ii] / voxel_volume; // Volume weighting
  }

  if(config->Dose_Segmentation != 0){
    for(ii=0; ii<scoring->Nbr_voxels; ii++){
      scoring->dose[ii] *= (ct->density[ii] > config->Segmentation_Density_Threshold);
    }
  }

  if(config->Score_Energy == 1){
    for(ii=0; ii<scoring->Nbr_voxels; ii++){
      scoring->energy[ii] = scoring->energy[ii] * normalization / Nbr_simulated_primaries;
      scoring->energy[ii] *= (scoring->energy[ii] > 0);
    }
  }

  if(config->Score_PromptGammas == 1){
    for(ii=0; ii<scoring->Nbr_voxels; ii++){
      scoring->PG_particles[ii] = scoring->PG_particles[ii] * normalization / Nbr_simulated_primaries;
    }
  }

  if(config->Score_LET == 1){
    for(ii=0; ii<scoring->Nbr_voxels; ii++){
      scoring->LET[ii] = (scoring->LET[ii] > 0) * scoring->LET[ii] / ((scoring->LET_denominator[ii] * 1e7) + FLT_EPSILON);
    }
  }

}


VAR_SCORING Process_batch(DATA_Scoring *Tot_scoring, DATA_Scoring *batch, Materials *material, DATA_CT *ct, int Num_batch, DATA_config *config){

  //double voxel_volume = Tot_scoring->VoxelLength[0]*Tot_scoring->VoxelLength[1]*Tot_scoring->VoxelLength[2];
  VAR_SCORING tmp, sigma=0, max_dose=0;
  int count=0;

  #pragma omp parallel for reduction(max: max_dose)
  for(int j=0; j<Tot_scoring->Nbr_voxels; j++){
    //batch->dose[j] = batch->dose[j] / voxel_volume; // Volume weighting
    Tot_scoring->dose[j] += batch->dose[j];
    Tot_scoring->dose_squared[j] += batch->dose[j] * batch->dose[j];
    if(ct->density[j] > 0.1 && max_dose < Tot_scoring->dose[j]) max_dose = Tot_scoring->dose[j];
  }

  if(config->Score_Energy == 1){
    #pragma omp parallel for
    for(int j=0; j<Tot_scoring->Nbr_voxels; j++){
      Tot_scoring->energy[j] += batch->energy[j];
    }
  }

  if(config->Score_LET == 1){
    #pragma omp parallel for
    for(int j=0; j<Tot_scoring->Nbr_voxels; j++){
      Tot_scoring->LET[j] += batch->LET[j];
      Tot_scoring->LET_denominator[j] += batch->LET_denominator[j];
    }
  }

  if(config->Score_PromptGammas == 1){
    #pragma omp parallel for
    for(int j=0; j<Tot_scoring->Nbr_voxels; j++){
      Tot_scoring->PG_particles[j] += batch->PG_particles[j];
    }

    #pragma omp parallel for
    for(int j=0; j<config->PG_Spectrum_NumBin; j++) Tot_scoring->PG_spectrum[j] += batch->PG_spectrum[j];
  }

  count = 0;
  for(int j=0; j<Tot_scoring->Nbr_voxels; j++){
    if(Tot_scoring->dose[j] > 0.5*max_dose){
      //sigma += sqrt(Tot_scoring->dose_squared[j] - Tot_scoring->dose[j]*Tot_scoring->dose[j]/Num_batch) / (Tot_scoring->dose[j]/Num_batch);
      if(config->Ignore_low_density_voxels == 0 || ct->density[j] > 0.1){
        sigma += sqrt(Num_batch * (Tot_scoring->dose_squared[j]*Num_batch/(Tot_scoring->dose[j]*Tot_scoring->dose[j]) - 1.0) );
        count++;
      }
    }
  }

  sigma = sigma / ((long)count*(long)Num_batch);



  VAR_SCORING *batch_dose = NULL;

  if(config->Export_batch_dose == 1){
    batch_dose = (VAR_SCORING*)calloc(Tot_scoring->Nbr_voxels, sizeof(VAR_SCORING));

    double voxel_volume = Tot_scoring->VoxelLength[0]*Tot_scoring->VoxelLength[1]*Tot_scoring->VoxelLength[2];

    int ii;
    for(ii=0; ii<Tot_scoring->Nbr_voxels; ii++){
      batch_dose[ii] = Tot_scoring->dose[ii] / (Num_batch*config->Num_Primaries/MIN_NUM_BATCH);
      if(config->DoseToWater == 1) batch_dose[ii] /= material[ct->material[ii]].SPR; // dose to water conversion by post-processing
      batch_dose[ii] *= (batch_dose[ii] > 0);
      batch_dose[ii] = batch_dose[ii] / voxel_volume; // Volume weighting
      if(config->Dose_Segmentation != 0) batch_dose[ii] *= (ct->density[ii] > config->Segmentation_Density_Threshold);
    }

    char file_path[100];
    strcpy(file_path, config->Output_Directory);
    strcat(file_path, "Batch_Dose");
    strcat(file_path, config->output_robustness_suffix);
    strcat(file_path, config->output_beamlet_suffix);
    strcat(file_path, config->output_4D_suffix);
    strcat(file_path, config->output_beams_suffix);
    strcat(file_path, ".mhd");
    export_MHD_image(file_path, Tot_scoring->GridSize, Tot_scoring->VoxelLength, Tot_scoring->Origin, batch_dose);

    if(batch_dose != NULL) free(batch_dose);

    file_path[strlen(file_path)-4] = '\0';
    strcat(file_path, "_stat.txt");
    FILE *file = NULL;
    file = fopen(file_path, "w");
    fprintf(file, "Num_simulated_batches = %d\n", Num_batch);
    fprintf(file, "Num_simulated_primaries = %ld\n", (Num_batch*config->Num_Primaries/MIN_NUM_BATCH));
    fprintf(file, "Estimated_mean_uncertainty = %.3f %%\n", 100*sigma);
    fclose(file);
  }

  return sigma;
}


void Free_Scoring(DATA_Scoring *scoring){

  if(scoring->energy != NULL) free(scoring->energy);
  if(scoring->dose != NULL) free(scoring->dose);
  if(scoring->dose_squared != NULL) free(scoring->dose_squared);

  if(scoring->PG_particles != NULL) free(scoring->PG_particles);
  if(scoring->PG_spectrum != NULL) free(scoring->PG_spectrum);

  if(scoring->LET != NULL) free(scoring->LET);
  if(scoring->LET_denominator != NULL) free(scoring->LET_denominator);

}
