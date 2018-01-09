/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/



#include "include/define.h"
#include "include/struct.h"
#include "include/data_config.h"
#include "include/data_materials.h"
#include "include/data_ct.h"
#include "include/data_4D.h"
#include "include/data_ct_penelope.h"
#include "include/data_mhd.h"
#include "include/data_ct_mhd.h"
#include "include/compute_4D.h"
#include "include/compute_simulation.h"
#include "include/compute_simulation_beamlet.h"
#include "include/compute_treatment_uncertainties.h"
#include "include/compute_random.h"
#include "include/data_beam_model.h"
#include "include/compute_beam_model.h"
#include "include/compute_math.h"
#include "include/compute_range_shifter.h"
#include "include/compute_scenario_simulation.h"

//#include <unistd.h>

int main(int argc, char *argv[]){

  //////////////////////
  // Initialisation
  //////////////////////

  //Parse configuration file
  DATA_config config;
  int error;
  if(argc < 2) error = Parse_Config(&config, "config.txt");
  else error = Parse_Config(&config, argv[1]);

  if(error != 0) return 1;
  //display_config(&config);

  config.timestamp = time(NULL);

  if(config.Num_Threads == 0) config.Num_Threads = omp_get_num_procs();
  else if(config.Num_Threads > omp_get_num_procs()){
    printf("\n Warning: Num_Threads is set to %u but only %u processing units are available\n\n", config.Num_Threads, omp_get_num_procs());
  }

  omp_set_num_threads(config.Num_Threads);
  double time_start = omp_get_wtime();

  if(config.Robustness_Mode == 1 && config.Simu_4D_Mode == 1 && (config.Systematic_Amplitude_Error != 0.0 || config.Random_Amplitude_Error != 0.0) && config.Field_type != 1){
    printf("\nError: Velocity fields are required to simulate variation of breathing amplitude!\n");
    return 1;
  }
  if(config.Robustness_Mode == 1 && config.Simu_4D_Mode == 1 && (config.Systematic_Amplitude_Error != 0.0 || config.Random_Amplitude_Error != 0.0) && config.Create_4DCT_from_Ref == 0){
    printf("\nError: Generating 4D-CT from reference image is required to simulate variation of breathing amplitude!\n");
    return 1;
  }

  // Import deformation fields
  DATA_4D_Fields *Fields = NULL;
  if(config.Simu_4D_Mode == 1 && (config.Dose_4D_Accumulation == 1 || config.Create_Ref_from_4DCT == 1 || config.Create_4DCT_from_Ref == 1)){
    Fields = Import_4D_Fields(&config);
    if(Fields == NULL) return 1;
  }
  else Fields = NULL;

  // Import CT image(s)
  DATA_CT **CT_phases = NULL;
  DATA_CT *ct = NULL;
  if(config.Simu_4D_Mode == 1){  // 4D mode

    if(config.Create_Ref_from_4DCT == 0){	
      if(config.Create_4DCT_from_Ref == 1){
        ct = Read_CT_MHD(&config);		// Import reference image
        if(ct == NULL){
          Free_4D_Fields(Fields);
          return 1;
        }
      }
    }
    else{					// Generate reference image from 4DCT
      CT_phases = Import_4DCT(&config);
      if(CT_phases == NULL){
        Free_4D_Fields(Fields);
        return 1;
      }
      ct = Create_Ref_from_4DCT(&config, Fields, CT_phases);
    }

    if(config.Create_4DCT_from_Ref == 0){	// Import 4DCT
      if(config.Create_Ref_from_4DCT == 0){
        CT_phases = Import_4DCT(&config);
        if(CT_phases == NULL){
          Free_4D_Fields(Fields);
          return 1;
        }
        ct = CT_phases[0];
      }
    }
    else{					// Generate 4DCT from reference image
      CT_phases = Create_4DCT_from_Ref(&config, Fields, ct);
      if(CT_phases == NULL){
        Free_4D_Fields(Fields);
        return 1;
      }
    }

    config.Num_Primaries = (unsigned long)config.Num_Primaries / config.Num_4DCT_phases;
  }

  else{						// 3D mode
    //ct = Read_PENCT(config.CT_File);
    ct = Read_CT_MHD(&config);
    if(ct == NULL) return 1;
    config.Num_4DCT_phases = 0;
  }

  //Display_Density_conversion_data(ct);
  //Display_Material_conversion_data(ct);


  // Import material database
  Materials *material = Init_materials(&config.Num_Materials, &config.Num_Components, &config);
  if(material == NULL){
    if(config.Simu_4D_Mode == 0) Free_CT_DATA(ct);
    else{
	Free_4DCT(CT_phases, config.Num_4DCT_phases);
	if(config.Dose_4D_Accumulation == 1) Free_4D_Fields(Fields);
    }
    Free_Materials_DATA(material, config.Num_Materials);
    return 1;
  }

  machine_parameters machine;
  error = read_machine_parameters(config.BDL_machine, &machine);
  if(error != 0) return 1;
//  display_machine_parameters(&machine);

  // Import PBS plan
  plan_parameters *plan = read_plan_parameters(config.BDL_plan, &config);
  if(plan == NULL){
    if(config.Simu_4D_Mode == 0) Free_CT_DATA(ct);
    else{
	Free_4DCT(CT_phases, config.Num_4DCT_phases);
	if(config.Dose_4D_Accumulation == 1) Free_4D_Fields(Fields);
    }
    Free_Materials_DATA(material, config.Num_Materials);
    Free_Plan_Parameters(plan);
    return 1;
  }
  if(config.Export_Beam_dose == 1) config.Num_Primaries = (unsigned long)config.Num_Primaries / plan->NumberOfFields;
//  display_plan_parameters(plan);

  error = Init_RangeShifter_Data(plan, &machine, material, &config);
  if(error == 1){
    if(config.Simu_4D_Mode == 0) Free_CT_DATA(ct);
    else{
	Free_4DCT(CT_phases, config.Num_4DCT_phases);
	if(config.Dose_4D_Accumulation == 1) Free_4D_Fields(Fields);
    }
    Free_Materials_DATA(material, config.Num_Materials);
    Free_Plan_Parameters(plan);
    return 1;
  }
  Display_RangeShifter_Data(plan, &machine, material);


  // Export density map
  char file_path[200];
  int a;
  if(config.Densities_Output == 1){
    if(config.Simu_4D_Mode == 0){
      sprintf(file_path, "%sDensities.mhd", config.Output_Directory);
      export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, ct->density);
    }
    else{
      sprintf(file_path, "%sDensities_ref.mhd", config.Output_Directory);
      export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, ct->density);
      for(a=0; a <config.Num_4DCT_phases; a++){
        sprintf(file_path, "%sDensities_%d.mhd", config.Output_Directory,a);
      	export_MHD_image(file_path, CT_phases[a]->GridSize, CT_phases[a]->VoxelLength, CT_phases[a]->density);
      }
    }
  }

  // Export material map
  if(config.Materials_Output == 1){
    int o;
    VAR_SCORING *mat = (VAR_SCORING*)malloc(ct->Nbr_voxels * sizeof(VAR_SCORING));
    if(config.Simu_4D_Mode == 0){
      for(o=0; o<ct->Nbr_voxels; o++){ mat[o] = (VAR_SCORING)ct->material[o]; }
      sprintf(file_path, "%sMaterials_out.mhd", config.Output_Directory);
      export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, mat);
    }
    else{
      for(a=0; a <config.Num_4DCT_phases; a++){
	for(o=0; o<CT_phases[a]->Nbr_voxels; o++){ mat[o] = (VAR_SCORING)CT_phases[a]->material[o]; }
        sprintf(file_path, "%sMaterials_out_%d.mhd", config.Output_Directory,a);
      	export_MHD_image(file_path, CT_phases[a]->GridSize, CT_phases[a]->VoxelLength, mat);
      }
    }
    free(mat);
  }

  //////////////////////
  // Start computation
  //////////////////////

  config.Current_fraction = plan->NumberOfFractions;
  config.Fraction_accumulation = 0;

  if(config.Robustness_Mode == 0){	// Not Robustness_Mode

    double time_init = omp_get_wtime();
    printf("\nInitialization time: %f s \n\n", (time_init-time_start));

    config.Current_scenario_type = Regular;
    config.Current_Systematic_setup[0] = 0.0;
    config.Current_Systematic_setup[1] = 0.0;
    config.Current_Systematic_setup[2] = 0.0;
    config.Current_Random_setup[0] = 0.0;
    config.Current_Random_setup[1] = 0.0;
    config.Current_Random_setup[2] = 0.0;
    config.Current_Range_error = 0.0;
    config.Current_Systematic_amplitude = 0.0;
    config.Current_Random_amplitude = 0.0;
    config.Current_Breathing_amplitude = 1.0;
    strcpy(config.output_robustness_suffix, "");
    config.Current_Systematic_period = 0.0;
    config.Current_Random_period = 0.0;
    config.Current_Breathing_period = config.Breathing_period;
    for(a=0;a<<plan->NumberOfFields;a++) config.Current_init_delivery_points[a] = 0.0;

    Scenario_simulation(&config, material, ct, CT_phases, plan, &machine, Fields);

    double time_end = omp_get_wtime();
    printf("Total computation time: %f s \n", (time_end-time_start));
  }


  else{ // Robustness Mode

    ct->Nominal_density = ct->density;
    ct->Scaled_density = (VAR_DATA*)malloc(ct->Nbr_voxels * sizeof(VAR_DATA));
    if(config.Simu_4D_Mode == 1){
      for(a=0; a <config.Num_4DCT_phases; a++){
        CT_phases[a]->Nominal_density = CT_phases[a]->density;
        CT_phases[a]->Scaled_density = (VAR_DATA*)malloc(CT_phases[a]->Nbr_voxels * sizeof(VAR_DATA));
      }
    }

    strcpy(file_path, config.Output_Directory);
    strcat(file_path, "Robustness_scenarios.txt");

    FILE *file_hdl = NULL;
    file_hdl = fopen(file_path, "w");
    fprintf(file_hdl, "Robustness Parameters:\n");
    fprintf(file_hdl, "----------------------\n");
    fprintf(file_hdl, "Systematic Setup Error = %.3f %.3f %.3f cm\n", config.Systematic_Setup_Error[0], config.Systematic_Setup_Error[1], config.Systematic_Setup_Error[2]);
    fprintf(file_hdl, "Random Setup Error = %.3f %.3f %.3f cm\n", config.Random_Setup_Error[0], config.Random_Setup_Error[1], config.Random_Setup_Error[2]);
    fprintf(file_hdl, "Systematic Range Error = %.2f %%\n", config.Systematic_Range_Error);
    if(config.Simu_4D_Mode == 1){
      fprintf(file_hdl, "Systematic motion amplitude error = %.2f %%\n", config.Systematic_Amplitude_Error);
      fprintf(file_hdl, "Random motion amplitude error = %.2f %%\n", config.Random_Amplitude_Error);
    }
    if(config.Dynamic_delivery == 1){
      fprintf(file_hdl, "Systematic motion period error = %.2f %%\n", config.Systematic_Period_Error);
      fprintf(file_hdl, "Random motion period error = %.2f %%\n", config.Random_Period_Error);
    }
    if(config.Scenario_selection == 1) fprintf(file_hdl, "Scenario selection: random sampling\n");
    else fprintf(file_hdl, "Scenario selection: all combinations\n");
    fprintf(file_hdl, "\n");
    fprintf(file_hdl, "Uncertainty scenarios:\n");
    fprintf(file_hdl, "----------------------\n");
    fclose(file_hdl);

    double time_init = omp_get_wtime();
    printf("\nInitialization time: %f s \n\n", (time_init-time_start));


    // Nominal plan:
    if(config.Simulate_nominal_plan == 1){
	config.Current_scenario_type = Nominal;
    	config.Current_Systematic_setup[0] = 0.0;
    	config.Current_Systematic_setup[1] = 0.0;
    	config.Current_Systematic_setup[2] = 0.0;
    	config.Current_Random_setup[0] = 0.0;
    	config.Current_Random_setup[1] = 0.0;
    	config.Current_Random_setup[2] = 0.0;
    	config.Current_Range_error = 0.0;
    	config.Current_Systematic_amplitude = 0.0;
    	config.Current_Random_amplitude = 0.0;
    	config.Current_Breathing_amplitude = 1.0;
    	config.Current_Systematic_period = 0.0;
    	config.Current_Random_period = 0.0;
        config.Current_Breathing_period = config.Breathing_period;
        for(a=0;a<<plan->NumberOfFields;a++) config.Current_init_delivery_points[a] = 0.0;

    	strcpy(config.output_beamlet_suffix, "");
    	strcpy(config.output_robustness_suffix, "_Nominal");

	file_hdl = fopen(file_path, "a");
	fprintf(file_hdl, "Nominal: ");
	fprintf(file_hdl, "Systematic_Setup(%.2f %.2f %.2f mm) ", 10*config.Current_Systematic_setup[0], 10*config.Current_Systematic_setup[1], 10*config.Current_Systematic_setup[2]);
	fprintf(file_hdl, "Random_Setup(%.2f %.2f %.2f mm) ", 10*config.Current_Random_setup[0], 10*config.Current_Random_setup[1], 10*config.Current_Random_setup[2]);
	fprintf(file_hdl, "Systematic_Range(%+.2f %%) ", config.Current_Range_error);
	if(config.Simu_4D_Mode == 1) fprintf(file_hdl, "Motion_amplitude(%.1f %%) ", 100*config.Current_Breathing_amplitude);
        if(config.Dynamic_delivery == 1){
          fprintf(file_hdl, "Motion_period(%.1f %%) ", 100*config.Current_Breathing_period);
          fprintf(file_hdl, "Start_delivery(%.1f %% period) ", 100*config.Current_init_delivery_points[0]);
        }
	fprintf(file_hdl, "\n");

	fclose(file_hdl);

	Scenario_simulation(&config, material, ct, CT_phases, plan, &machine, Fields);
    }

    // Robustness scenarios:
    if(config.Scenario_selection == 1) Scenarios_selection_random(&config, material, ct, CT_phases, plan, &machine, Fields, file_path);
    else Scenarios_selection_all(&config, material, ct, CT_phases, plan, &machine, Fields, file_path);


    // Free dynamic variables
    if(config.Simu_4D_Mode == 0){ 	// 3D mode
      if(ct->Nominal_density != NULL) free(ct->Nominal_density);
      if(ct->Scaled_density != NULL) free(ct->Scaled_density);
      ct->density = NULL;
    }
    else{
      for(a=0; a <config.Num_4DCT_phases; a++){
        if(CT_phases[a]->Nominal_density != NULL) free(CT_phases[a]->Nominal_density);
        if(CT_phases[a]->Scaled_density != NULL) free(CT_phases[a]->Scaled_density);
        CT_phases[a]->density = NULL;
      }
    }


  }



  //////////////////////
  // Clean Memory
  //////////////////////


  if(config.Simu_4D_Mode == 0) Free_CT_DATA(ct);
  else{
	Free_4DCT(CT_phases, config.Num_4DCT_phases);
	if(config.Dose_4D_Accumulation == 1 || config.Create_4DCT_from_Ref == 1) Free_4D_Fields(Fields);
  }
  Free_Materials_DATA(material, config.Num_Materials);
  Free_Plan_Parameters(plan);
  Free_Machine_Parameters(&machine);

  return 0;
}
