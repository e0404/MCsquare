/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_simulation.h"


void Run_simulation(DATA_config *config, Materials *material, DATA_CT *ct, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields){


  
  double time_init, time_MC, time_end;
  time_init = omp_get_wtime();

  char file_path[100];
  unsigned long Num_simulated_primaries = 0;
  VAR_SCORING stat_uncertainty = 1.0;
  int Num_batch = MIN_NUM_BATCH;

  DATA_Scoring Batch_scoring, Tot_scoring;
  Tot_scoring = Init_Scoring(config, ct->Nbr_voxels, 1);

  if(config->Compute_stat_uncertainty == 1 && (config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 0) && config->Fraction_accumulation == 0){

    int batch = 1;
    while(batch<=Num_batch){
      Batch_scoring = Init_Scoring(config, ct->Nbr_voxels, 1);
      Num_simulated_primaries += Simulation_loop(config, material, ct, plan, machine, Fields, &Batch_scoring, (unsigned long)config->Num_Primaries/MIN_NUM_BATCH);
      stat_uncertainty = Process_batch(&Tot_scoring, &Batch_scoring, ct, batch, config);
      Free_Scoring(&Batch_scoring);

      if(config->Stat_uncertainty == 0.0){
        if(batch < 5) printf(" %.1f %% \n", batch*(100.0/MIN_NUM_BATCH));
        else printf(" %.1f %% (stat uncertainty: %.2f %%) \n", batch*(100.0/MIN_NUM_BATCH), stat_uncertainty*100);
      }
      else{
        if(batch < 5) printf("batch %d completed \n", batch);
        else printf("batch %d completed (stat uncertainty: %.2f %%) \n", batch, stat_uncertainty*100);

	if(batch == Num_batch && stat_uncertainty*100 > config->Stat_uncertainty)  Num_batch++;
      }

      batch++;
    }

  }

  else{  // no uncertainty evaluation (1 batch)

    Num_simulated_primaries += Simulation_loop(config, material, ct, plan, machine, Fields, &Tot_scoring, (unsigned long)config->Num_Primaries);

  }

  PostProcess_Scoring(&Tot_scoring, ct, material, plan->normalization_factor, Num_simulated_primaries, config); // convert energy to dose per proton, etc...
  
  time_MC = omp_get_wtime();

  // Accumulation of multiple 4D phases and fractions (if required):

  int export_results = 0;
  static VAR_SCORING *energy_accumulation = NULL;
  static VAR_SCORING *dose_accumulation = NULL;
  static VAR_SCORING *PG_accumulation = NULL;
  static VAR_SCORING *PG_Spectrum_accumulation = NULL;
  static VAR_SCORING *LET_accumulation = NULL;
  VAR_COMPUTE *deformed;

  int ii;

  if((config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 0) && config->Fraction_accumulation == 0) export_results = 1;
  else{

    strcpy(config->output_4D_suffix, "");


    if(config->Energy_ASCII_Output == 1 || config->Energy_MHD_Output == 1 || config->Energy_Sparse_Output == 1){

      if(config->Current_4D_phase == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) energy_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

      if(config->Simu_4D_Mode == 0){
	for(ii=0; ii<ct->Nbr_voxels; ii++) energy_accumulation[ii] += Tot_scoring.energy[ii] / plan->NumberOfFractions;
      }
      else{
	deformed = Image_deformation(Tot_scoring.energy, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[config->Current_4D_phase], Fields->GridSize, Fields->Spacing, Fields->Origin);
	for(ii=0; ii<ct->Nbr_voxels; ii++) energy_accumulation[ii] += deformed[ii] / (config->Num_4DCT_phases*plan->NumberOfFractions);
      	free(deformed);
      }
      
      if(config->Current_4D_phase == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	export_results = 1;
	if(Tot_scoring.energy != NULL) free(Tot_scoring.energy);
        Tot_scoring.energy = energy_accumulation;
      }
    }


    if(config->Compute_DVH == 1 || config->Dose_ASCII_Output == 1 || config->Dose_MHD_Output == 1 || config->Dose_Sparse_Output == 1){

      if(config->Current_4D_phase == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) dose_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

      if(config->Simu_4D_Mode == 0){
	for(ii=0; ii<ct->Nbr_voxels; ii++) dose_accumulation[ii] += Tot_scoring.dose[ii] / plan->NumberOfFractions;
      }
      else{
	deformed = Image_deformation(Tot_scoring.dose, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[config->Current_4D_phase], Fields->GridSize, Fields->Spacing, Fields->Origin);
	for(ii=0; ii<ct->Nbr_voxels; ii++) dose_accumulation[ii] += deformed[ii] / (config->Num_4DCT_phases*plan->NumberOfFractions);
      	free(deformed);
      }
      
      if(config->Current_4D_phase == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	export_results = 1;
	if(Tot_scoring.dose != NULL) free(Tot_scoring.dose);
        Tot_scoring.dose = dose_accumulation;
      }
    }


    if(config->Score_PromptGammas == 1){

      if(config->Current_4D_phase == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)){
	PG_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));
	PG_Spectrum_accumulation = (VAR_SCORING*)calloc(config->PG_Spectrum_NumBin, sizeof(VAR_SCORING));
      }

      if(config->Simu_4D_Mode == 0){
	for(ii=0; ii<ct->Nbr_voxels; ii++) PG_accumulation[ii] += Tot_scoring.PG_particles[ii] / plan->NumberOfFractions;
      }
      else{
	deformed = Image_deformation(Tot_scoring.PG_particles, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[config->Current_4D_phase], Fields->GridSize, Fields->Spacing, Fields->Origin);
	for(ii=0; ii<ct->Nbr_voxels; ii++) PG_accumulation[ii] += deformed[ii] / (config->Num_4DCT_phases*plan->NumberOfFractions);
      	free(deformed);
      }

      for(ii=0; ii<config->PG_Spectrum_NumBin; ii++) PG_Spectrum_accumulation[ii] += Tot_scoring.PG_spectrum[ii];
      
      if(config->Current_4D_phase == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	export_results = 1;
	if(Tot_scoring.PG_particles != NULL) free(Tot_scoring.PG_particles);
	if(Tot_scoring.PG_spectrum != NULL) free(Tot_scoring.PG_spectrum);
        Tot_scoring.PG_particles = PG_accumulation;
	Tot_scoring.PG_spectrum = PG_Spectrum_accumulation;
      }
    }


    if(config->Score_LET == 1){

      if(config->Current_4D_phase == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) LET_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

      if(config->Simu_4D_Mode == 0){
	for(ii=0; ii<ct->Nbr_voxels; ii++) LET_accumulation[ii] += Tot_scoring.LET[ii] / plan->NumberOfFractions;
      }
      else{
	deformed = Image_deformation(Tot_scoring.LET, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[config->Current_4D_phase], Fields->GridSize, Fields->Spacing, Fields->Origin);
      	for(ii=0; ii<ct->Nbr_voxels; ii++) LET_accumulation[ii] += deformed[ii] / (config->Num_4DCT_phases*plan->NumberOfFractions);
      	free(deformed);
      }
      
      if(config->Current_4D_phase == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	export_results = 1;
	if(Tot_scoring.LET != NULL) free(Tot_scoring.LET);
        Tot_scoring.LET = LET_accumulation;
      }
    }

  }


  // Export results

  // Convert the dose in Gray units for DVH calculation
  // 1 eV = 1.602176e-19 J
  VAR_SCORING DoseScaling = 1.602176e-19 * 1000 * config->Num_delivered_protons * plan->NumberOfFractions;

  if(export_results == 1 && config->Compute_DVH == 1) compute_all_DVH(config, Tot_scoring.dose, DoseScaling);

  if(export_results == 1 && config->Energy_ASCII_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Energy");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".dat");
	export_dose_ascii(file_path, ct->GridSize, Tot_scoring.energy);
  }
  if(export_results == 1 && config->Energy_MHD_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Energy");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".mhd");
	export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.energy);
  }
  if(export_results == 1 && config->Energy_Sparse_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Sparse_Energy");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".txt");
	export_Sparse_image(file_path, config, ct, plan, Tot_scoring.energy, config->Energy_Sparse_Threshold);
  }
  if(export_results == 1 && config->Dose_ASCII_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Dose");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".dat");
	export_dose_ascii(file_path, ct->GridSize, Tot_scoring.dose);
  }
  if(export_results == 1 && config->Dose_MHD_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Dose");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".mhd");
	export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.dose);
  }
  if(export_results == 1 && config->Dose_Sparse_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Sparse_Dose");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".txt");
	export_Sparse_image(file_path, config, ct, plan, Tot_scoring.dose, config->Dose_Sparse_Threshold);
  }
  if(export_results == 1 && config->LET_ASCII_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "LET");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".dat");
	export_dose_ascii(file_path, ct->GridSize, Tot_scoring.LET);
  }
  if(export_results == 1 && config->LET_MHD_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "LET");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_beamlet_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".mhd");
	export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.LET);
  }
  if(export_results == 1 && config->LET_Sparse_Output == 1){
	strcpy(file_path, config->Output_Directory);
	strcat(file_path, "Sparse_LET");
	strcat(file_path, config->output_robustness_suffix);
	strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
	strcat(file_path, ".txt");
	export_Sparse_image(file_path, config, ct, plan, Tot_scoring.LET, config->LET_Sparse_Threshold);
  }


  if(export_results == 1 && config->Score_PromptGammas == 1){
    
    strcpy(file_path, config->Output_Directory);
    strcat(file_path, "PromptGamma");
    strcat(file_path, config->output_robustness_suffix);
    strcat(file_path, config->output_beamlet_suffix);
    strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
    strcat(file_path, ".dat");
    export_PG_ascii(file_path, ct->GridSize, Tot_scoring.PG_particles);
    
    strcpy(file_path, config->Output_Directory);
    strcat(file_path, "PromptGamma_spectrum");
    strcat(file_path, config->output_robustness_suffix);
    strcat(file_path, config->output_beamlet_suffix);
    strcat(file_path, config->output_4D_suffix);
	strcat(file_path, config->output_beams_suffix);
    strcat(file_path, ".dat");
    export_PG_spectrum_ascii(file_path, config->PG_Spectrum_NumBin, config->PG_Spectrum_Binning, Tot_scoring.PG_spectrum);
  }


  time_end = omp_get_wtime();

  if(config->Robustness_Mode == 0 && config->Beamlet_Mode == 0){
    if(config->Particle_Generated_outside != 0) printf("\nNbr primaries simulated: %lu (%u generated outside the geometry) \n", Num_simulated_primaries, config->Particle_Generated_outside);
    else printf("\nNbr primaries simulated: %lu \n", Num_simulated_primaries);
  }

  printf("MC computation time: %f s \n", (time_MC-time_init));
  printf("Output computation time: %f s \n\n", (time_end-time_MC));

  // Delete dynamic variables
  Free_Scoring(&Tot_scoring);

}


unsigned long Simulation_loop(DATA_config *config, Materials *material, DATA_CT *ct, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields, DATA_Scoring *Tot_scoring, unsigned long Num_primaries){

  static int Num_call = 0;
  Num_call++;

  unsigned long Num_simulated_primaries = 0;

  VAR_SCORING *ptr_energy_scoring[config->Num_Threads];
  VAR_SCORING *ptr_PG_scoring[config->Num_Threads];
  VAR_SCORING *ptr_PG_spectrum[config->Num_Threads];
  VAR_SCORING *ptr_LET_scoring[config->Num_Threads];
  VAR_SCORING *ptr_LET_denominator[config->Num_Threads];

  // Parallelisation
  #pragma omp parallel shared(config, material, ct, plan, machine, Fields, Num_simulated_primaries, Tot_scoring, ptr_energy_scoring, ptr_PG_scoring, ptr_PG_spectrum, ptr_LET_scoring, ptr_LET_denominator)
  {
    int tid = omp_get_thread_num();

    // Init RNG
    VSLStreamStatePtr RNDstream;				// variable for the random number generator
    ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];			// variable that contains random numbers
    if(config->RNG_Seed == 0){
      vslNewStream(&RNDstream, VSL_BRNG_MCG59, time(NULL)+tid*1e4+Num_call*1e5);	// initialize the RNG for each thread individually with a seed = time+thread_id*10000
    }
    else{
      vslNewStream(&RNDstream, VSL_BRNG_MCG59, config->RNG_Seed+tid*1e4+Num_call*1e5);
    }
    rand_uniform(RNDstream, v_rnd);				// the RNG is called here because random numbers seems not well distributed the first time.

    // Init scoring
    DATA_Scoring scoring = Init_Scoring(config, ct->Nbr_voxels, 0);


    // Init particle stacks
    Hadron hadron;
    Init_particles(&hadron);
    Hadron_buffer HadronToSimulate[100];
    int Nbr_HadronToSimulate = 0;

    // variables
    int i, count, stop = 0;
    float progress_binning = 10.0;
    float progress_interval = progress_binning * Num_primaries / 100;
    float progress_next = progress_interval;

    int display_progress;
    if(config->Compute_stat_uncertainty == 0 && (config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 0) && config->Robustness_Mode == 0 && config->Beamlet_Mode == 0) display_progress = 1;
    else display_progress = 0;


    // Compute simulation
    while(stop == 0){
      for(i=0; i<VLENGTH; i++){
	if(hadron.v_type[i] == Unknown){

	  if(Nbr_HadronToSimulate > 0){
	    Nbr_HadronToSimulate -= 1;
	    Insert_particle(&hadron, i, &HadronToSimulate[Nbr_HadronToSimulate]);
	  }

	  else if(Num_simulated_primaries < Num_primaries){
	    #pragma omp atomic
	    Num_simulated_primaries += VLENGTH;

	    //Generate_particle(&hadron, i, BeamPOSx, BeamPOSy, BeamPOSz, PEnergy*UMeV);
	    Generate_PBS_particle(HadronToSimulate, &Nbr_HadronToSimulate, ct->Length, plan, machine, RNDstream, config, material);
	    if(Nbr_HadronToSimulate > 0){
	      Nbr_HadronToSimulate -= 1;
	      Insert_particle(&hadron, i, &HadronToSimulate[Nbr_HadronToSimulate]);
	    }
	  }

	  else{
	    count = __sec_reduce_add(hadron.v_type[vALL]);
	    if(count == 0) stop = 1;
	  }

	  if(tid == 0 && display_progress == 1 && Num_simulated_primaries > progress_next){
	    printf(" %.1f %% \n", floor(Num_simulated_primaries/progress_interval)*progress_binning);
	    fflush(stdout);
	    progress_next += progress_interval;
	  }

/*
	  if(tid == 0){
	    if(config->Robustness_Mode == 0 && config->Beamlet_Mode == 0 && Num_simulated_primaries > progress_next){
	      printf(" %.1f %% \n", floor(Num_simulated_primaries/progress_interval)*progress_binning);
	      fflush(stdout);
	      progress_next += progress_interval;
	    }
	  }
*/
	}
      }

      hadron_step(&hadron, &scoring, material, ct, HadronToSimulate, &Nbr_HadronToSimulate, RNDstream, config);
    }

    if(tid == 0 && display_progress == 1) printf(" 100.0 %% \n");


//    if(tid == 0 && config->Robustness_Mode == 0 && config->Beamlet_Mode == 0) printf(" 100.0 %% \n");

    // Agreggate results

/*
    #pragma omp critical
    {
      int l;
      for(l=0; l<ct->Nbr_voxels; l++){
        Tot_scoring.energy[l] += scoring.energy[l];
        if(config->Score_PromptGammas == 1){
          Tot_scoring.PG_particles[l] += scoring.PG_particles[l];
	}
      }
      if(config->Score_PromptGammas == 1){
        for(l=0; l<config->PG_Spectrum_NumBin; l++){
          Tot_scoring.PG_spectrum[l] += scoring.PG_spectrum[l];
        }
      }
    }
*/

    ptr_energy_scoring[tid] = scoring.energy;
    if(config->Score_PromptGammas == 1){
      ptr_PG_scoring[tid] = scoring.PG_particles;
      ptr_PG_spectrum[tid] = scoring.PG_spectrum;
    }
    if(config->Score_LET == 1){
      ptr_LET_scoring[tid] = scoring.LET;
      ptr_LET_denominator[tid] = scoring.LET_denominator;
    }

    int k,p;

    #pragma omp barrier

    #pragma omp for
    for(k=0; k<ct->Nbr_voxels; ++k){
      for(p=0; p<config->Num_Threads; ++p){
         Tot_scoring->energy[k] += ptr_energy_scoring[p][k]; 
	 if(config->Score_PromptGammas == 1){
           Tot_scoring->PG_particles[k] += ptr_PG_scoring[p][k];   
	 }
      }
    }

    if(config->Score_PromptGammas == 1){
      #pragma omp for
      for(k=0; k<config->PG_Spectrum_NumBin; ++k){
        for(p=0; p<config->Num_Threads; ++p){
          Tot_scoring->PG_spectrum[k] += ptr_PG_spectrum[p][k]; 
        }
      }
    }

    if(config->Score_LET == 1){
      #pragma omp for
      for(k=0; k<ct->Nbr_voxels; ++k){
        for(p=0; p<config->Num_Threads; ++p){
	 Tot_scoring->LET[k] += ptr_LET_scoring[p][k];
	 Tot_scoring->LET_denominator[k] += ptr_LET_denominator[p][k];
        }
      }
    }


    // Delete dynamic variables
    Free_Scoring(&scoring);

  }  // end of Parallelization

  return Num_simulated_primaries;

}
