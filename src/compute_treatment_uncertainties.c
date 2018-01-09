/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_treatment_uncertainties.h"

void Translation_uncertainty(Hadron_buffer *hadron, DATA_config *config, VSLStreamStatePtr RNG_Stream){

  int i;

/*
  if(config->Current_Random_setup[0] != 0.0) hadron->x = hadron->x - single_rand_normal(RNG_Stream, 0, config->Current_Random_setup[0]) - config->Current_Systematic_setup[0];
  else hadron->x = hadron->x - config->Current_Systematic_setup[0];

  if(config->Current_Random_setup[1] != 0.0) hadron->y = hadron->y - single_rand_normal(RNG_Stream, 0, config->Current_Random_setup[1]) - config->Current_Systematic_setup[1];
  else hadron->y = hadron->y - config->Current_Systematic_setup[1];

  if(config->Current_Random_setup[2] != 0.0) hadron->z = hadron->z - single_rand_normal(RNG_Stream, 0, config->Current_Random_setup[2]) - config->Current_Systematic_setup[2];
  else hadron->z = hadron->z - config->Current_Systematic_setup[2];
*/

hadron->x = hadron->x - config->Current_Systematic_setup[0] - config->Current_Random_setup[0];
hadron->y = hadron->y - config->Current_Systematic_setup[1] - config->Current_Random_setup[1];
hadron->z = hadron->z - config->Current_Systematic_setup[2] - config->Current_Random_setup[2];

	// Sampling of sigma before the Gaussian sampling (test for Ana)
  	//VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream);
	//hadron->x = hadron->x - single_rand_normal(RNG_Stream, 0, config->Current_Random_error[0]*rnd) - config->Current_Systematic_error[0];
	//hadron->y = hadron->y - single_rand_normal(RNG_Stream, 0, config->Current_Random_error[1]*rnd) - config->Current_Systematic_error[1];
	//hadron->z = hadron->z - single_rand_normal(RNG_Stream, 0, config->Current_Random_error[2]*rnd) - config->Current_Systematic_error[2];


}



void Density_scaling(VAR_DATA *Nominal_density, VAR_DATA *Scaled_density, int Num_voxels, VAR_DATA scaling_factor){

  int i;
  for(i=0; i<Num_voxels; i++){
    Scaled_density[i] = scaling_factor * Nominal_density[i];
  }
}



void Breathing_amplitude_variation(DATA_config *config, DATA_CT *ct, DATA_CT **CT_phases, DATA_4D_Fields *Fields){

  int phaseID, i;
  VAR_DATA *tmp = NULL;

  double time_init = omp_get_wtime();
  printf("Generating new 4DCT with %f %% motion amplitude ", 100*config->Current_Breathing_amplitude);

  for(phaseID=0; phaseID < config->Num_4DCT_phases; phaseID++){

    tmp = Field_multiplication(config->Current_Breathing_amplitude, Fields->Ref2Phase_log[phaseID], Fields->GridSize);

    if(Fields->Phase2Ref[phaseID] != NULL) free(Fields->Phase2Ref[phaseID]);
    Fields->Phase2Ref[phaseID] = Field_exponentiation(tmp, Fields->GridSize, Fields->Spacing, Fields->Origin, 1);

    if(Fields->Ref2Phase[phaseID] != NULL) free(Fields->Ref2Phase[phaseID]);
    Fields->Ref2Phase[phaseID] = Field_exponentiation(tmp, Fields->GridSize, Fields->Spacing, Fields->Origin, 0);

    free(tmp);

    if(CT_phases[phaseID]->density == CT_phases[phaseID]->Scaled_density){
      if(CT_phases[phaseID]->density != NULL) free(CT_phases[phaseID]->density);
      CT_phases[phaseID]->density = Image_deformation(ct->density, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Ref2Phase[phaseID], Fields->GridSize, Fields->Spacing, Fields->Origin);
      CT_phases[phaseID]->Scaled_density = CT_phases[phaseID]->density;
    }
    else{
      if(CT_phases[phaseID]->density != NULL) free(CT_phases[phaseID]->density);
      CT_phases[phaseID]->density = Image_deformation(ct->density, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Ref2Phase[phaseID], Fields->GridSize, Fields->Spacing, Fields->Origin);
      CT_phases[phaseID]->Nominal_density = CT_phases[phaseID]->density;
    }

    #pragma omp parallel for private(i)
    for(i=0; i<CT_phases[phaseID]->Nbr_voxels; i++){
      CT_phases[phaseID]->material[i] = (unsigned short int)Density_to_Material_convertion(CT_phases[phaseID]->density[i], CT_phases[phaseID]);
    }
  }

  printf("(%f s) \n", omp_get_wtime() - time_init);

}



plan_parameters* Spot_Sorting(DATA_config *config, int phase, plan_parameters *Plan){

  double time_init = omp_get_wtime();
  printf("Generating partial plan for phase %d ", phase+1);

  VAR_DATA SpotPhase, cumulative_weight=0;
  VAR_DATA PhaseDuration = config->Current_Breathing_period * 1000 / config->Num_4DCT_phases;
  VAR_DATA InitTime;

  plan_parameters *partial_plan = (plan_parameters*)malloc(sizeof(plan_parameters));
  strcpy(partial_plan->PlanName, "Partial plan"); 
  partial_plan->NumberOfFractions = Plan->NumberOfFractions;
  partial_plan->FractionID = Plan->FractionID;
  partial_plan->NumberOfFields = Plan->NumberOfFields;
  partial_plan->TotalMetersetWeightOfAllFields = 0.0;
  partial_plan->fields = (field_parameters*)malloc(Plan->NumberOfFields * sizeof(field_parameters));
  partial_plan->Fields_cumulative_PDF = (VAR_DATA*)malloc(Plan->NumberOfFields * sizeof(VAR_DATA));
  partial_plan->FieldsID = (int*)malloc(Plan->NumberOfFields * sizeof(int));

  int i,j,k;
  for(i=0; i<Plan->NumberOfFields; i++){
    if(config->Export_Beam_dose == 1) InitTime = config->Current_init_delivery_points[config->Current_Beam] * config->Current_Breathing_period * 1000;
    else InitTime = config->Current_init_delivery_points[i] * config->Current_Breathing_period * 1000;

    partial_plan->FieldsID[i] = Plan->FieldsID[i];
    partial_plan->fields[i].FieldID = Plan->fields[i].FieldID;
    partial_plan->fields[i].FinalCumulativeMeterSetWeight = 0.0;
    partial_plan->fields[i].GantryAngle = Plan->fields[i].GantryAngle;
    partial_plan->fields[i].PatientSupportAngle = Plan->fields[i].PatientSupportAngle;
    partial_plan->fields[i].IsocenterPositionX = Plan->fields[i].IsocenterPositionX;
    partial_plan->fields[i].IsocenterPositionY = Plan->fields[i].IsocenterPositionY;
    partial_plan->fields[i].IsocenterPositionZ = Plan->fields[i].IsocenterPositionZ;
    partial_plan->fields[i].NumberOfControlPoints = Plan->fields[i].NumberOfControlPoints;
    partial_plan->fields[i].RS_Type = Plan->fields[i].RS_Type;
    partial_plan->fields[i].ControlPoints = (ControlPoint_parameters*)malloc(Plan->fields[i].NumberOfControlPoints * sizeof(ControlPoint_parameters));
    partial_plan->fields[i].ControlPoints_cumulative_PDF = (VAR_DATA*)malloc(Plan->fields[i].NumberOfControlPoints * sizeof(VAR_DATA));
  
    for(j=0; j<Plan->fields[i].NumberOfControlPoints; j++){
      partial_plan->fields[i].ControlPoints[j].ControlPointIndex = Plan->fields[i].ControlPoints[j].ControlPointIndex;
      partial_plan->fields[i].ControlPoints[j].SpotTunnedID = Plan->fields[i].ControlPoints[j].SpotTunnedID;
      partial_plan->fields[i].ControlPoints[j].CumulativeMetersetWeight = 0.0;
      partial_plan->fields[i].ControlPoints[j].Energy = Plan->fields[i].ControlPoints[j].Energy;
      partial_plan->fields[i].ControlPoints[j].NbOfScannedSpots = Plan->fields[i].ControlPoints[j].NbOfScannedSpots;
      partial_plan->fields[i].ControlPoints[j].RS_setting = Plan->fields[i].ControlPoints[j].RS_setting;
      partial_plan->fields[i].ControlPoints[j].RS_IsocenterDist = Plan->fields[i].ControlPoints[j].RS_IsocenterDist;
      partial_plan->fields[i].ControlPoints[j].RS_WET = Plan->fields[i].ControlPoints[j].RS_WET;
      partial_plan->fields[i].ControlPoints[j].RS_Thickness = Plan->fields[i].ControlPoints[j].RS_Thickness;
      partial_plan->fields[i].ControlPoints[j].spots = (spot_parameters*)malloc(Plan->fields[i].ControlPoints[j].NbOfScannedSpots * sizeof(spot_parameters));
      partial_plan->fields[i].ControlPoints[j].Spots_cumulative_PDF = (VAR_DATA*)malloc(Plan->fields[i].ControlPoints[j].NbOfScannedSpots * sizeof(VAR_DATA));

      for(k=0; k<Plan->fields[i].ControlPoints[j].NbOfScannedSpots; k++){
        partial_plan->fields[i].ControlPoints[j].spots[k].Spot_X = Plan->fields[i].ControlPoints[j].spots[k].Spot_X;
        partial_plan->fields[i].ControlPoints[j].spots[k].Spot_Y = Plan->fields[i].ControlPoints[j].spots[k].Spot_Y;
        partial_plan->fields[i].ControlPoints[j].spots[k].Spot_Time = Plan->fields[i].ControlPoints[j].spots[k].Spot_Time;

        SpotPhase = fmod(floor((Plan->fields[i].ControlPoints[j].spots[k].Spot_Time + InitTime) / PhaseDuration), config->Num_4DCT_phases);
        if(SpotPhase == phase) partial_plan->fields[i].ControlPoints[j].spots[k].Spot_Weight = Plan->fields[i].ControlPoints[j].spots[k].Spot_Weight;
        else partial_plan->fields[i].ControlPoints[j].spots[k].Spot_Weight = 0.0;

        cumulative_weight += partial_plan->fields[i].ControlPoints[j].spots[k].Spot_Weight;
        partial_plan->fields[i].ControlPoints[j].Spots_cumulative_PDF[k] = cumulative_weight;
      }
      partial_plan->fields[i].ControlPoints_cumulative_PDF[j] = cumulative_weight;
      partial_plan->fields[i].ControlPoints[j].CumulativeMetersetWeight = cumulative_weight;
    }

    partial_plan->Fields_cumulative_PDF[i] = cumulative_weight;
    partial_plan->fields[i].FinalCumulativeMeterSetWeight = cumulative_weight;
  }

  partial_plan->TotalMetersetWeightOfAllFields = cumulative_weight;
  partial_plan->cumulative_weight = cumulative_weight;

  partial_plan->normalization_factor = Plan->normalization_factor * partial_plan->Fields_cumulative_PDF[Plan->NumberOfFields-1] / Plan->Fields_cumulative_PDF[Plan->NumberOfFields-1];
  if(config->Dose_4D_Accumulation == 1) partial_plan->normalization_factor *= config->Num_4DCT_phases;

  printf("(%f s) \n", omp_get_wtime() - time_init);

  return partial_plan;
}



