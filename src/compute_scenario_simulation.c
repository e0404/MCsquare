/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_scenario_simulation.h"


void Scenarios_selection_all(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path, enum Scenario_selection_shape Selection_shape){

  int a, x, y, z, r;
  FILE *file_hdl = NULL;
  VAR_COMPUTE setup_norm, range_norm;

  config->TotalNumScenarios = 81;
  config->Current_scenario = 0;
  config->Current_scenario_type = Uncertainty;

  if(config->Systematic_Setup_Error[0] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios / 3;
  if(config->Systematic_Setup_Error[1] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios / 3;
  if(config->Systematic_Setup_Error[2] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios / 3;
  if(config->Systematic_Range_Error == 0.0) config->TotalNumScenarios = config->TotalNumScenarios / 3;

  for(r=-1; r<=1; r++){
    for(x=-1; x<=1; x++){
      for(y=-1; y<=1; y++){
        for(z=-1; z<=1; z++){
        
          // skip simulation of scenarios when error is null
          if(config->Systematic_Range_Error == 0.0 && r != 0) continue; 
          if(config->Systematic_Setup_Error[0] == 0.0 && x != 0) continue;
          if(config->Systematic_Setup_Error[2] == 0.0 && y != 0) continue;
          if(config->Systematic_Setup_Error[1] == 0.0 && z != 0) continue;
          
	      config->Current_scenario += 1;
	      
	      // scenarios selection shape in the error space
	      if(Selection_shape == HyperCube){
	        setup_norm = 1.0;
	        range_norm = 1.0;
	      }
	      else if(Selection_shape == HyperCylinder){
	        setup_norm = sqrt((double)x*x + (double)y*y + (double)z*z) + FLT_EPSILON;
	        range_norm = 1.0;
	      }
	      else{ // HyperSphere
	        setup_norm = sqrt((double)x*x + (double)y*y + (double)z*z + (double)r*r) + FLT_EPSILON;
	        range_norm = sqrt((double)x*x + (double)y*y + (double)z*z + (double)r*r) + FLT_EPSILON;
	      }
          
          // Systematic range error
          config->Current_Range_error = r * config->Systematic_Range_Error / range_norm;
          if(config->Simu_4D_Mode == 0){
            Density_scaling(ct->Nominal_density, ct->Scaled_density, ct->Nbr_voxels, 1.0 + (config->Current_Range_error * 0.01));
            ct->density = ct->Scaled_density;
          }
          else{ 
            for(a=0; a <config->Num_4DCT_phases; a++){
              Density_scaling(CT_phases[a]->Nominal_density, CT_phases[a]->Scaled_density, ct->Nbr_voxels, 1.0 + (config->Current_Range_error * 0.01));
              CT_phases[a]->density = CT_phases[a]->Scaled_density;
            }
          }

          // Systematic setup error
          config->Current_Systematic_setup[0] = x * config->Systematic_Setup_Error[0] / setup_norm;
          config->Current_Systematic_setup[1] = y * config->Systematic_Setup_Error[1] / setup_norm;
          config->Current_Systematic_setup[2] = z * config->Systematic_Setup_Error[2] / setup_norm;

          // Random setup error
          config->Current_Random_setup[0] = config->Random_Setup_Error[0];
          config->Current_Random_setup[1] = config->Random_Setup_Error[1];
          config->Current_Random_setup[2] = config->Random_Setup_Error[2];

          // log
          sprintf(config->output_robustness_suffix, "_Scenario_%d-%d", config->Current_scenario, config->TotalNumScenarios);
          file_hdl = fopen(file_path, "a");
          fprintf(file_hdl, "Scenario (%d/%d): Systematic_Setup(%.3f %.3f %.3f mm) Random_Setup(%.3f %.3f %.3f mm) Systematic_Range(%.2f %%)\n", config->Current_scenario, config->TotalNumScenarios, 10*config->Current_Systematic_setup[0], 10*config->Current_Systematic_setup[1], 10*config->Current_Systematic_setup[2], 10*config->Current_Random_setup[0], 10*config->Current_Random_setup[1], 10*config->Current_Random_setup[2], config->Current_Range_error);
          fclose(file_hdl);

          // run simulation
          Scenario_simulation(config, material, ct, CT_phases, plan, machine, Fields);
          
        }
      }
    }
  }

}


void Scenarios_selection_reduced(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path, enum Scenario_selection_shape Selection_shape){

  int a, setup_direction, setup_shift, r;
  FILE *file_hdl = NULL;
  VAR_COMPUTE norm;

  config->TotalNumScenarios = 21;
  config->Current_scenario = 0;
  config->Current_scenario_type = Uncertainty;

  if(config->Systematic_Setup_Error[0] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios - 6;
  if(config->Systematic_Setup_Error[1] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios - 6;
  if(config->Systematic_Setup_Error[2] == 0.0) config->TotalNumScenarios = config->TotalNumScenarios - 6;
  if(config->Systematic_Range_Error == 0.0) config->TotalNumScenarios = config->TotalNumScenarios / 3;

  for(r=-1; r<=1; r++){
    for (setup_direction=0; setup_direction<=2; setup_direction++){
      for (setup_shift=-1; setup_shift<=1; setup_shift++){
        
        // skip simulation of scenarios when error is null
        if(config->Systematic_Range_Error == 0.0 && r != 0) continue; 
        if(config->Systematic_Setup_Error[setup_direction] == 0.0 && setup_shift != 0) continue;
        if(setup_shift == 0 && setup_direction != 0) continue;
          
	    config->Current_scenario += 1;
	      
	    // scenarios selection shape in the error space
	    if(Selection_shape == HyperCube) norm = 1.0;
	    else if(Selection_shape == HyperCylinder) norm = 1.0;
	    else norm = sqrt((double)setup_shift*setup_shift + (double)r*r) + FLT_EPSILON; // HyperSphere
          
        // Systematic range error
        config->Current_Range_error = r * config->Systematic_Range_Error / norm;
        if(config->Simu_4D_Mode == 0){
          Density_scaling(ct->Nominal_density, ct->Scaled_density, ct->Nbr_voxels, 1.0 + (config->Current_Range_error * 0.01));
          ct->density = ct->Scaled_density;
        }
        else{ 
          for(a=0; a <config->Num_4DCT_phases; a++){
            Density_scaling(CT_phases[a]->Nominal_density, CT_phases[a]->Scaled_density, ct->Nbr_voxels, 1.0 + (config->Current_Range_error * 0.01));
            CT_phases[a]->density = CT_phases[a]->Scaled_density;
          }
        }

        // Systematic setup error
        config->Current_Systematic_setup[0] = 0.0;
	    config->Current_Systematic_setup[1] = 0.0;
	    config->Current_Systematic_setup[2] = 0.0;
	    config->Current_Systematic_setup[setup_direction] = setup_shift * config->Systematic_Setup_Error[setup_direction] / norm;

        // Random setup error
        config->Current_Random_setup[0] = config->Random_Setup_Error[0];
        config->Current_Random_setup[1] = config->Random_Setup_Error[1];
        config->Current_Random_setup[2] = config->Random_Setup_Error[2];

        // log
        sprintf(config->output_robustness_suffix, "_Scenario_%d-%d", config->Current_scenario, config->TotalNumScenarios);
        file_hdl = fopen(file_path, "a");
        fprintf(file_hdl, "Scenario (%d/%d): Systematic_Setup(%.3f %.3f %.3f mm) Random_Setup(%.3f %.3f %.3f mm) Systematic_Range(%.2f %%)\n", config->Current_scenario, config->TotalNumScenarios, 10*config->Current_Systematic_setup[0], 10*config->Current_Systematic_setup[1], 10*config->Current_Systematic_setup[2], 10*config->Current_Random_setup[0], 10*config->Current_Random_setup[1], 10*config->Current_Random_setup[2], config->Current_Range_error);
        fclose(file_hdl);

        // run simulation
        Scenario_simulation(config, material, ct, CT_phases, plan, machine, Fields);
          
      }
    }
  }

}





void Scenarios_selection_random(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, char *file_path){

  int i, j, a, f;
  FILE *file_hdl = NULL;

  // Init RNG
  VSLStreamStatePtr RNDstream;					// un stream de RNG
  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];			// vecteur de nbr aleatoires
  if(config->RNG_Seed == 0){
    vslNewStream(&RNDstream, VSL_BRNG_MCG59, time(NULL));	// initialisation du stream du RNG avec le seed (time+thread_id)
  }
  else{
    vslNewStream(&RNDstream, VSL_BRNG_MCG59, config->RNG_Seed);
  }
  rand_uniform(RNDstream, v_rnd);				// on genere une première fois un set de nbr car les premiers semblent mal distribués


  config->TotalNumScenarios = config->Num_random_scenarios;
  config->Current_scenario = 0;
  config->Current_scenario_type = Uncertainty;
  config->Fraction_accumulation = 1;
  config->Num_Primaries = (unsigned long)config->Num_Primaries / plan->NumberOfFractions;

  for(i=0; i<config->TotalNumScenarios; i++){

    config->Current_scenario += 1;

    // Systematic setup error
    config->Current_Systematic_setup[0] = single_rand_normal(RNDstream, 0, config->Systematic_Setup_Error[0]);
    config->Current_Systematic_setup[1] = single_rand_normal(RNDstream, 0, config->Systematic_Setup_Error[1]);
    config->Current_Systematic_setup[2] = single_rand_normal(RNDstream, 0, config->Systematic_Setup_Error[2]);

    // Systematic Range error
    config->Current_Range_error = single_rand_normal(RNDstream, 0, config->Systematic_Range_Error);
    Density_scaling(ct->Nominal_density, ct->Scaled_density, ct->Nbr_voxels, 1.0 - (config->Current_Range_error * 0.01));
    ct->density = ct->Scaled_density;
    if(config->Simu_4D_Mode == 1){ 
      for(a=0; a <config->Num_4DCT_phases; a++){
	Density_scaling(CT_phases[a]->Nominal_density, CT_phases[a]->Scaled_density, ct->Nbr_voxels, 1.0 - (config->Current_Range_error * 0.01));
	CT_phases[a]->density = CT_phases[a]->Scaled_density;
      }
    }

    // Systematic breathing amplitude variation
    if(config->Simu_4D_Mode != 0 && config->Systematic_Amplitude_Error != 0.0) config->Current_Systematic_amplitude = single_rand_normal(RNDstream, 0, config->Systematic_Amplitude_Error);
    else config->Current_Systematic_amplitude = 0.0;

    // Systematic breathing period variation
    if(config->Simu_4D_Mode != 0 && config->Dynamic_delivery == 1 && config->Systematic_Period_Error != 0.0) config->Current_Systematic_period = single_rand_normal(RNDstream, 0, config->Systematic_Period_Error);
    else config->Current_Systematic_period = 0.0;
    
    // 4D not accumulated mode: randomly sample a phase (ITV like evaluation)
    if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0){
      config->Current_4D_phase = floor(config->Num_4DCT_phases * single_rand_uniform(RNDstream));
      if(config->Current_4D_phase >= config->Num_4DCT_phases) config->Current_4D_phase = config->Num_4DCT_phases - 1;
    }

    config->Current_fraction = 0;

    for(j=0; j<plan->NumberOfFractions; j++){

      config->Current_fraction += 1;
    
      // Random setup error
      config->Current_Random_setup[0] = single_rand_normal(RNDstream, 0, config->Random_Setup_Error[0]);
      config->Current_Random_setup[1] = single_rand_normal(RNDstream, 0, config->Random_Setup_Error[1]);
      config->Current_Random_setup[2] = single_rand_normal(RNDstream, 0, config->Random_Setup_Error[2]);

      // Random breathing amplitude variation
      if(config->Simu_4D_Mode != 0 && config->Random_Amplitude_Error != 0.0) config->Current_Random_amplitude = single_rand_normal(RNDstream, 0, config->Random_Amplitude_Error);
      else config->Current_Random_amplitude = 0.0;
      config->Current_Breathing_amplitude = 1.0 + 0.01 * (config->Current_Systematic_amplitude + config->Current_Random_amplitude);
      if(config->Simu_4D_Mode != 0 && config->Current_Breathing_amplitude != 1.0) Breathing_amplitude_variation(config, ct, CT_phases, Fields);

      // Systematic breathing period variation
      if(config->Simu_4D_Mode != 0 && config->Dynamic_delivery == 1 && config->Random_Period_Error != 0.0) config->Current_Random_period = single_rand_normal(RNDstream, 0, config->Random_Period_Error);
      else config->Current_Random_period = 0.0;
      config->Current_Breathing_period = config->Breathing_period * (1.0 + 0.01 * (config->Current_Systematic_amplitude + config->Current_Random_amplitude));

      // Starting delivery point (for each beam)
      if(config->Simu_4D_Mode != 0 && config->Dynamic_delivery == 1){
	for(f=0;f<plan->NumberOfFields;f++) config->Current_init_delivery_points[f] = single_rand_uniform(RNDstream);
      }
      else{
	for(f=0;f<plan->NumberOfFields;f++) config->Current_init_delivery_points[f] = 0.0;
      }

    
      // Simulation
      sprintf(config->output_robustness_suffix, "_Scenario_%d-%d", config->Current_scenario, config->TotalNumScenarios);
	    
      file_hdl = fopen(file_path, "a");
      fprintf(file_hdl, "Scenario (%d/%d): ", config->Current_scenario, config->TotalNumScenarios);
      fprintf(file_hdl, "Systematic_Setup(%.2f %.2f %.2f mm) ", 10*config->Current_Systematic_setup[0], 10*config->Current_Systematic_setup[1], 10*config->Current_Systematic_setup[2]);
      fprintf(file_hdl, "Random_Setup(%.2f %.2f %.2f mm) ", 10*config->Current_Random_setup[0], 10*config->Current_Random_setup[1], 10*config->Current_Random_setup[2]);
      fprintf(file_hdl, "Systematic_Range(%+.2f %%) ", config->Current_Range_error);
      if(config->Simu_4D_Mode == 1) fprintf(file_hdl, "Motion_amplitude(%.1f %%) ", 100*config->Current_Breathing_amplitude);
      if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0) fprintf(file_hdl, "4D_phase(%.1d) ", config->Current_4D_phase+1);
      if(config->Dynamic_delivery == 1){
        fprintf(file_hdl, "Motion_period(%.1f s) ", config->Current_Breathing_period);
        fprintf(file_hdl, "Start_delivery(");
        for(f=0;f<plan->NumberOfFields;f++) fprintf(file_hdl, "%.1f%% ", 100*config->Current_init_delivery_points[f]);
        fprintf(file_hdl, "period) ");
      }
      fprintf(file_hdl, "\n");
      fclose(file_hdl);

      Scenario_simulation(config, material, ct, CT_phases, plan, machine, Fields);

  }}


  // reset to nominal data
  if(config->Simu_4D_Mode != 0 && (config->Systematic_Amplitude_Error != 0.0 || config->Random_Amplitude_Error != 0.0)){
    config->Current_Breathing_amplitude = 1.0;
    Breathing_amplitude_variation(config, ct, CT_phases, Fields);
    config->Num_Primaries = (unsigned long)config->Num_Primaries * plan->NumberOfFractions;
  }

}




void Scenario_simulation(DATA_config *config, Materials *material, DATA_CT *ct, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields){

  int a,f;
  int tmp_current_phase, tmp_num_phases;
  VAR_DATA *tmp_field_ptr;

  strcpy(config->output_beams_suffix, "");

  if(config->Beamlet_Mode == 1 && config->Beamlet_Parallelization == 1){	// Beamlet mode
    if(config->Simu_4D_Mode == 0) Run_simulation_beamlet(config, material, &ct, plan, machine, Fields);
    else Run_simulation_beamlet(config, material, CT_phases, plan, machine, Fields);
  }
  else if(config->Beamlet_Mode == 1 && config->Beamlet_Parallelization == 0){

    int b,c,d, current_spot;

    // Create a new plan containing only the current spot
    plan_parameters *Beamlet = Init_single_spot_plan(plan);

    current_spot = 0;
    for(b=0; b < plan->NumberOfFields; b++){
      for(c=0; c < plan->fields[b].NumberOfControlPoints; c++){
        for(d=0; d < plan->fields[b].ControlPoints[c].NbOfScannedSpots; d++){
          current_spot++;
          Select_spot(plan, Beamlet, b, c, d);
          sprintf(config->output_beamlet_suffix, "_Beamlet_%d_%d_%d", b, c, d);

          if(config->Simu_4D_Mode == 0){ 	// 3D mode
            config->Current_4D_phase = 0;
            if(config->Current_scenario_type == Nominal)
              printf("\nRobustness simulation (Nominal - Beamlet %d/%d) ", current_spot, config->TotalNbrSpots);
            else if(config->Current_scenario_type == Uncertainty)
              printf("\nRobustness simulation (scenario %d/%d - Beamlet %d/%d) ", config->Current_scenario, config->TotalNumScenarios, current_spot, config->TotalNbrSpots);
            else
              printf("\nBeamlet %d / %d \n", current_spot, config->TotalNbrSpots);

            sprintf(config->output_4D_suffix, "");
            Run_simulation(config, material, ct, Beamlet, machine, Fields);
          }
          else{ 	// 4D mode
            for(a=0; a <config->Num_4DCT_phases; a++){
              if(config->Current_scenario_type == Nominal)
                printf("\nRobustness simulation (Nominal - Beamlet %d/%d - phase %d/%d) \n", current_spot, config->TotalNbrSpots, a+1, config->Num_4DCT_phases);
              else if(config->Current_scenario_type == Uncertainty)
                printf("\nRobustness simulation (scenario %d/%d - Beamlet %d/%d - phase %d/%d) \n", config->Current_scenario, config->TotalNumScenarios, current_spot, config->TotalNbrSpots, a+1, config->Num_4DCT_phases);
              else
                printf("\nBeamlet %d / %d  (phase %d) \n", current_spot, config->TotalNbrSpots, a+1);

              sprintf(config->output_4D_suffix, "_Phase%d", a+1);
              config->Current_4D_phase = a;
              Run_simulation(config, material, CT_phases[a], Beamlet, machine, Fields);
            }
          }
        }
      }
    }

    Free_Plan_Parameters(Beamlet);

  }
  else{	// Full simulation mode (not beamlet)
  
    strcpy(config->output_beamlet_suffix, "");

    plan_parameters *full_plan;
    int SubPlans = 1;

    if(config->Export_Beam_dose == 1){
      SubPlans = plan->NumberOfFields;
      full_plan = plan;
    }

    for(f=0; f<SubPlans; f++){
      if(config->Export_Beam_dose == 1){
        sprintf(config->output_beams_suffix, "_Beam%d", f+1);
        plan = Select_beam(full_plan, f);
        config->Current_Beam = f;
      }

      if(config->Simu_4D_Mode == 0){ 	// 3D mode
        config->Current_4D_phase = 0;
        sprintf(config->output_4D_suffix, "");
        if(config->Current_scenario_type == Nominal) printf("\nRobustness simulation (Nominal)\n");
        else if(config->Current_scenario_type == Uncertainty) printf("\nRobustness simulation (scenario %d/%d - fraction %d/%d)\n", config->Current_scenario, config->TotalNumScenarios, config->Current_fraction, plan->NumberOfFractions);
      
        Run_simulation(config, material, ct, plan, machine, Fields);
      }
      
      else if(config->Robustness_Mode == 1 && config->Current_scenario_type != Nominal && config->Dose_4D_Accumulation == 0){ // 4D robustness test without phase accumulation (ITV like)
        if(config->Current_scenario_type == Nominal) printf("\nRobustness simulation (Nominal)\n");
        else if(config->Current_scenario_type == Uncertainty) printf("\nRobustness simulation (scenario %d/%d - fraction %d/%d)\n", config->Current_scenario, config->TotalNumScenarios, config->Current_fraction, plan->NumberOfFractions);
        
        // force accumulation of 1 single phase to enable deformation to reference phase
        config->Dose_4D_Accumulation = 1;
        tmp_current_phase = config->Current_4D_phase;
        config->Current_4D_phase = 0;
        tmp_num_phases = config->Num_4DCT_phases;
        config->Num_4DCT_phases = 1;
        tmp_field_ptr = Fields->Phase2Ref[0];
        Fields->Phase2Ref[0] = Fields->Phase2Ref[tmp_current_phase];
        
        Run_simulation(config, material, CT_phases[tmp_current_phase], plan, machine, Fields);
        
        // restore initial config
        config->Dose_4D_Accumulation = 0;
        config->Current_4D_phase = tmp_current_phase;
        config->Num_4DCT_phases = tmp_num_phases;
        Fields->Phase2Ref[0] = tmp_field_ptr;
        
        
      } // 4D mode condition
      
      else{ 	// regular 4D mode
      
        for(a=0; a <config->Num_4DCT_phases; a++){
          if(config->Current_scenario_type == Nominal) printf("\nRobustness simulation (Nominal - phase %d/%d) \n", a+1, config->Num_4DCT_phases);
          else if(config->Current_scenario_type == Uncertainty) printf("\nRobustness simulation (scenario %d/%d - fraction %d/%d - phase %d/%d) \n", config->Current_scenario, config->TotalNumScenarios, config->Current_fraction, plan->NumberOfFractions, a+1, config->Num_4DCT_phases);
          else printf("\n4D simulation (phase %d / %d) \n", a+1, config->Num_4DCT_phases);
        
          sprintf(config->output_4D_suffix, "_Phase%d", a+1);
          config->Current_4D_phase = a;

          if(config->Dynamic_delivery == 1 && config->Current_scenario_type != Nominal){	// Interplay simulation (dynamic delivery)
            plan_parameters *partial_plan = Spot_Sorting(config, a, plan);
            Run_simulation(config, material, CT_phases[a], partial_plan, machine, Fields);
            Free_Plan_Parameters(partial_plan);
          }
          else{
            Run_simulation(config, material, CT_phases[a], plan, machine, Fields);	// No interplay simulation (breathing motion only)
          }

        } // 4D loop
      }   // 4D mode condition

      if(config->Export_Beam_dose == 1) Free_Plan_Parameters(plan);
    } 	  // beam loop

    if(config->Export_Beam_dose == 1) plan = full_plan;
  } 	  // beamlet mode condition
} 	  // end function


