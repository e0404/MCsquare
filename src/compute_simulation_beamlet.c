/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_simulation_beamlet.h"


void Run_simulation_beamlet(DATA_config *config, Materials *material, DATA_CT **CT_phases, plan_parameters *plan, machine_parameters *machine, DATA_4D_Fields *Fields){

  int spotID;

  int tot_phases;
  if(config->Simu_4D_Mode == 0) tot_phases = 1;
  else tot_phases = config->Num_4DCT_phases;

  // Parallelisation  
  #pragma omp parallel for shared(config, material, CT_phases, plan, machine, Fields, tot_phases) schedule(dynamic,1)
//  #pragma omp parallel for shared(config, material, CT_phases, plan, machine, Fields, tot_phases) ordered schedule(dynamic,1)
  for(spotID=0; spotID<config->TotalNbrSpots; spotID++){

    double time_init, time_MC, time_end;
    char file_path[200], output_beamlet_suffix[200], output_4D_suffix[200];

    int tid = omp_get_thread_num();

    VAR_SCORING *energy_accumulation;
    VAR_SCORING *dose_accumulation;
    VAR_SCORING *PG_accumulation;
    VAR_SCORING *PG_Spectrum_accumulation;
    VAR_SCORING *LET_accumulation;
    VAR_COMPUTE *deformed;
    VAR_COMPUTE norm_factor = 1;

    DATA_CT *ct = NULL;

    int a,b,c,d;

    // Create a new plan containing only the current spot
    plan_parameters *Beamlet = Init_single_spot_plan(plan);

    int current_spot = 0;
    int spot_selected = 0;
    for(b=0; b < plan->NumberOfFields; b++){
      for(c=0; c < plan->fields[b].NumberOfControlPoints; c++){
        for(d=0; d < plan->fields[b].ControlPoints[c].NbOfScannedSpots; d++){
	  if(current_spot == spotID){
	    Select_spot(plan, Beamlet, b, c, d);
	    spot_selected = 1;
	    break;
	  }
	  current_spot++;
	}
	if(spot_selected == 1) break;
      }
      if(spot_selected == 1) break;
    }

    for(a=0; a <tot_phases; a++){

      time_init = omp_get_wtime();

      if(config->Simu_4D_Mode == 0) ct = CT_phases[0];
      else ct = CT_phases[a];


      unsigned long Nbr_simulated_primaries = 0;

      // Init RNG
      VSLStreamStatePtr RNDstream;				// un stream de RNG par thread
      ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];			// vecteur de nbr aleatoires
      if(config->RNG_Seed == 0){
        vslNewStream(&RNDstream, VSL_BRNG_MCG59, time(NULL)+tid);	// initialisation du stream du RNG avec le seed (time+thread_id)
      }
      else{
        vslNewStream(&RNDstream, VSL_BRNG_MCG59, config->RNG_Seed+tid);
      }
      rand_uniform(RNDstream, v_rnd);				// on genere une première fois un set de nbr car les premiers semblent mal distribués

      // Init scoring    
      DATA_Scoring Tot_scoring = Init_Scoring(config, ct->Nbr_voxels, 1);

      // Init particle stacks
      Hadron hadron;
      Init_particles(&hadron);
      Hadron_buffer HadronToSimulate[100];
      int Nbr_HadronToSimulate = 0;

      // variables
      int i, count, stop = 0;

      // Compute simulation
      while(stop == 0){
        for(i=0; i<VLENGTH; i++){
	  if(hadron.v_type[i] == Unknown){

	    if(Nbr_HadronToSimulate > 0){
	      Nbr_HadronToSimulate -= 1;
	      Insert_particle(&hadron, i, &HadronToSimulate[Nbr_HadronToSimulate]);
	    }

	    else if(Nbr_simulated_primaries < config->Num_Primaries){

	      Nbr_simulated_primaries += VLENGTH;

	      //Generate_particle(&hadron, i, BeamPOSx, BeamPOSy, BeamPOSz, PEnergy*UMeV);
	      Generate_PBS_particle(HadronToSimulate, &Nbr_HadronToSimulate, ct->Length, Beamlet, machine, RNDstream, config, material);
	      if(Nbr_HadronToSimulate > 0){
	        Nbr_HadronToSimulate -= 1;
	        Insert_particle(&hadron, i, &HadronToSimulate[Nbr_HadronToSimulate]);
	      }
	    }

	    else{
	      count = __sec_reduce_add(hadron.v_type[vALL]);
	      if(count == 0) stop = 1;
	    }
	  }
        }

        hadron_step(&hadron, &Tot_scoring, material, ct, HadronToSimulate, &Nbr_HadronToSimulate, RNDstream, config);
      }

      time_MC = omp_get_wtime();


      // Scoring post processing: (convert energy to dose per proton, etc...)
      PostProcess_Scoring(&Tot_scoring, ct, material, plan->normalization_factor, Nbr_simulated_primaries, config);


      // 4D Accumulation:

      int export_results = 0;
      int ii;

      if((config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 0) && config->Fraction_accumulation == 0) export_results = 1;
      else{

        strcpy(config->output_4D_suffix, "");

        if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 1) norm_factor /= config->Num_4DCT_phases;
        if(config->Fraction_accumulation == 1) norm_factor /= plan->NumberOfFractions;

        if(config->Energy_ASCII_Output == 1 || config->Energy_MHD_Output == 1 || config->Energy_Sparse_Output == 1){

          if(a == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) energy_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

          if(config->Simu_4D_Mode == 0){
            for(ii=0; ii<ct->Nbr_voxels; ii++) energy_accumulation[ii] += Tot_scoring.energy[ii] * norm_factor;
          }
          else{
            deformed = Image_deformation(Tot_scoring.energy, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[a], Fields->GridSize, Fields->Spacing, Fields->Origin);
            for(ii=0; ii<ct->Nbr_voxels; ii++) energy_accumulation[ii] += deformed[ii] * norm_factor;
            free(deformed);
          }
      
          if(a == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	    export_results = 1;
	    if(Tot_scoring.energy != NULL) free(Tot_scoring.energy);
            Tot_scoring.energy = energy_accumulation;
          }
        }


        if(config->Compute_DVH == 1 || config->Dose_ASCII_Output == 1 || config->Dose_MHD_Output == 1 || config->Dose_Sparse_Output == 1){

          if(a == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) dose_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

          if(config->Simu_4D_Mode == 0){
            for(ii=0; ii<ct->Nbr_voxels; ii++) dose_accumulation[ii] += Tot_scoring.dose[ii] * norm_factor;
          }
          else{
            deformed = Image_deformation(Tot_scoring.dose, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[a], Fields->GridSize, Fields->Spacing, Fields->Origin);
            for(ii=0; ii<ct->Nbr_voxels; ii++) dose_accumulation[ii] += deformed[ii] * norm_factor;
            free(deformed);
          }
      
          if(a == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	    export_results = 1;
	    if(Tot_scoring.dose != NULL) free(Tot_scoring.dose);
            Tot_scoring.dose = dose_accumulation;
          }
        }


        if(config->Score_PromptGammas == 1){

          if(a == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)){
	    PG_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));
	    PG_Spectrum_accumulation = (VAR_SCORING*)calloc(config->PG_Spectrum_NumBin, sizeof(VAR_SCORING));
          }

	  if(config->Simu_4D_Mode == 0){
	    for(ii=0; ii<ct->Nbr_voxels; ii++) PG_accumulation[ii] += Tot_scoring.PG_particles[ii] * norm_factor;
	  }
	  else{
	    deformed = Image_deformation(Tot_scoring.PG_particles, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[a], Fields->GridSize, Fields->Spacing, Fields->Origin);
	    for(ii=0; ii<ct->Nbr_voxels; ii++) PG_accumulation[ii] += deformed[ii] * norm_factor;
	    free(deformed);
	  }

          for(ii=0; ii<config->PG_Spectrum_NumBin; ii++) PG_Spectrum_accumulation[ii] += Tot_scoring.PG_spectrum[ii];
      
          if(a == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	    export_results = 1;
	    if(Tot_scoring.PG_particles != NULL) free(Tot_scoring.PG_particles);
	    if(Tot_scoring.PG_spectrum != NULL) free(Tot_scoring.PG_spectrum);
            Tot_scoring.PG_particles = PG_accumulation;
	    Tot_scoring.PG_spectrum = PG_Spectrum_accumulation;
          }
        }


        if(config->Score_LET == 1){

          if(a == 0 && (config->Current_fraction == 1 || config->Fraction_accumulation == 0)) LET_accumulation = (VAR_SCORING*)calloc(ct->Nbr_voxels, sizeof(VAR_SCORING));

          if(config->Simu_4D_Mode == 0){
            for(ii=0; ii<ct->Nbr_voxels; ii++) LET_accumulation[ii] += Tot_scoring.LET[ii] * norm_factor;
          }
          else{
            deformed = Image_deformation(Tot_scoring.LET, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Phase2Ref[a], Fields->GridSize, Fields->Spacing, Fields->Origin);
            for(ii=0; ii<ct->Nbr_voxels; ii++) LET_accumulation[ii] += deformed[ii] * norm_factor;
            free(deformed);
          }
      
          if(a == (config->Num_4DCT_phases-1) && config->Current_fraction == plan->NumberOfFractions){
	    export_results = 1;
	    if(Tot_scoring.LET != NULL) free(Tot_scoring.LET);
            Tot_scoring.LET = LET_accumulation;
          }
        }

      }


      // Export results

      // Convert the dose in Gray units for DVH calculation
      // 1 eV = 1.602176e-19 J
      VAR_SCORING DoseScaling = 1.602176e-19 * 1000 * Beamlet->cumulative_weight * Beamlet->NumberOfFractions;

      if(export_results == 1){

          sprintf(output_beamlet_suffix, "_Beamlet_%d_%d_%d", b, c, d);
          if(config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 1){
	    sprintf(output_4D_suffix, "");
          }
          else{ 	// 4D mode
	    sprintf(output_4D_suffix, "_Phase%d", a+1);
          }

          if(config->Compute_DVH == 1){
	    #pragma omp critical (Outputs)
            {
	      sprintf(config->output_beamlet_suffix, output_beamlet_suffix);
	      sprintf(config->output_4D_suffix, output_4D_suffix);
	      config->Current_4D_phase = a;
	      compute_all_DVH(config, Tot_scoring.dose, DoseScaling);
	    }
	  }

          if(config->Energy_ASCII_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "Energy");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".dat");
	    export_dose_ascii(file_path, ct->GridSize, Tot_scoring.energy);
          }
          if(config->Energy_MHD_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "Energy");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".mhd");
	    export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.energy);
          }
          if(config->Dose_ASCII_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "Dose");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".dat");
	    export_dose_ascii(file_path, ct->GridSize, Tot_scoring.dose);
          }
          if(config->Dose_MHD_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "Dose");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".mhd");
	    export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.dose);
          }
  	  if(export_results == 1 && config->LET_ASCII_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "LET");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".dat");
	    export_dose_ascii(file_path, ct->GridSize, Tot_scoring.LET);
  	  }
  	  if(export_results == 1 && config->LET_MHD_Output == 1){
	    strcpy(file_path, config->Output_Directory);
	    strcat(file_path, "LET");
	    strcat(file_path, config->output_robustness_suffix);
	    strcat(file_path, output_beamlet_suffix);
	    strcat(file_path, output_4D_suffix);
	    strcat(file_path, ".mhd");
	    export_MHD_image(file_path, ct->GridSize, ct->VoxelLength, Tot_scoring.LET);
  	  }
          


          if(config->Score_PromptGammas == 1){
            strcpy(file_path, config->Output_Directory);
            strcat(file_path, "PromptGamma");
            strcat(file_path, config->output_robustness_suffix);
            strcat(file_path, output_beamlet_suffix);
            strcat(file_path, output_4D_suffix);
            strcat(file_path, ".dat");
            export_PG_ascii(file_path, ct->GridSize, Tot_scoring.PG_particles);
        
            strcpy(file_path, config->Output_Directory);
            strcat(file_path, "PromptGamma_spectrum");
            strcat(file_path, config->output_robustness_suffix);
            strcat(file_path, output_beamlet_suffix);
            strcat(file_path, output_4D_suffix);
            strcat(file_path, ".dat");
            export_PG_spectrum_ascii(file_path, config->PG_Spectrum_NumBin, config->PG_Spectrum_Binning, Tot_scoring.PG_spectrum);
          }


	  if(export_results == 1 && config->Energy_Sparse_Output == 1){
	      strcpy(file_path, config->Output_Directory);
	      strcat(file_path, "tmp/");
	      CreateDir(file_path);
	      sprintf(file_path, "%sBeamlet_%d/", file_path, current_spot+1);
	      CreateDir(file_path);
	      strcat(file_path, "Sparse_Energy");
	      strcat(file_path, config->output_robustness_suffix);
	      strcat(file_path, output_4D_suffix);
	      strcat(file_path, ".txt");
	      if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0){
		#pragma omp critical (Outputs)
                {
	          sprintf(config->output_beamlet_suffix, output_beamlet_suffix);
	          sprintf(config->output_4D_suffix, output_4D_suffix);
	          config->Current_4D_phase = a;
	          export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.energy, config->Energy_Sparse_Threshold);
	        }
	      }
	      else export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.energy, config->Energy_Sparse_Threshold);
	      
          }
	  if(export_results == 1 && config->Dose_Sparse_Output == 1){
	      strcpy(file_path, config->Output_Directory);
	      strcat(file_path, "tmp/");
	      CreateDir(file_path);
	      sprintf(file_path, "%sBeamlet_%d/", file_path, current_spot+1);
	      CreateDir(file_path);
	      strcat(file_path, "Sparse_Dose");
	      strcat(file_path, config->output_robustness_suffix);
	      strcat(file_path, output_4D_suffix);
	      strcat(file_path, ".txt");
	      if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0){
		#pragma omp critical (Outputs)
                {
	          sprintf(config->output_beamlet_suffix, output_beamlet_suffix);
	          sprintf(config->output_4D_suffix, output_4D_suffix);
	          config->Current_4D_phase = a;
	          export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.dose, config->Dose_Sparse_Threshold);
	        }
	      }
	      else export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.dose, config->Dose_Sparse_Threshold);
	      
          }
	  if(export_results == 1 && config->LET_Sparse_Output == 1){
	      strcpy(file_path, config->Output_Directory);
	      strcat(file_path, "tmp/");
	      CreateDir(file_path);
	      sprintf(file_path, "%sBeamlet_%d/", file_path, current_spot+1);
	      CreateDir(file_path);
	      strcat(file_path, "Sparse_LET");
	      strcat(file_path, config->output_robustness_suffix);
	      strcat(file_path, output_4D_suffix);
	      strcat(file_path, ".txt");
	      if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0){
		#pragma omp critical (Outputs)
                {
	          sprintf(config->output_beamlet_suffix, output_beamlet_suffix);
	          sprintf(config->output_4D_suffix, output_4D_suffix);
	          config->Current_4D_phase = a;
	          export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.LET, config->LET_Sparse_Threshold);
	        }
	      }
	      else export_Sparse_image(file_path, config, ct, Beamlet, Tot_scoring.LET, config->LET_Sparse_Threshold);
	      
          }

      } // end export_results == 1

      #pragma omp critical (Display)
      {

        time_end = omp_get_wtime();

        if(config->Simu_4D_Mode == 0) printf("\nSimulation of beamlet %d/%d %s \n", current_spot+1, config->TotalNbrSpots, config->output_robustness_suffix);
        else printf("\nSimulation of beamlet %d/%d Phase%d %s \n", current_spot+1, config->TotalNbrSpots, a+1, config->output_robustness_suffix);

        printf("MC computation time: %f s \n", (time_MC-time_init));
        printf("Output computation time: %f s \n", (time_end-time_MC));
 
      }
 
      // Delete dynamic variables
      Free_Scoring(&Tot_scoring);

    } // end 4D phases loop
    Free_Plan_Parameters(Beamlet);
  }  // end Parallelisation

  #pragma omp barrier

  if(config->Energy_Sparse_Output == 1 || config->Dose_Sparse_Output == 1){

    char InPath[200], InFile[200], OutPath[200];

    if(config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 1) tot_phases = 1;
    else tot_phases = config->Num_4DCT_phases;

    int aa;
    for(aa=0; aa <tot_phases; aa++){

      if(config->Simu_4D_Mode == 0 || config->Dose_4D_Accumulation == 1){
        sprintf(config->output_4D_suffix, "");
      }
      else{ 	// 4D mode
        sprintf(config->output_4D_suffix, "_Phase%d", aa+1);
        config->Current_4D_phase = aa;
      }

      if(config->Energy_Sparse_Output == 1){
        sprintf(InPath, "%stmp/Beamlet_", config->Output_Directory);
        sprintf(InFile, "Sparse_Energy%s%s.txt", config->output_robustness_suffix, config->output_4D_suffix);
        sprintf(OutPath, "%sSparse_Energy%s%s.txt", config->Output_Directory, config->output_robustness_suffix, config->output_4D_suffix);
	Merge_Sparse_Files(InPath, InFile, config->TotalNbrSpots, OutPath);
      }
      if(config->Dose_Sparse_Output == 1){
        sprintf(InPath, "%stmp/Beamlet_", config->Output_Directory);
        sprintf(InFile, "Sparse_Dose%s%s.txt", config->output_robustness_suffix, config->output_4D_suffix);
        sprintf(OutPath, "%sSparse_Dose%s%s.txt", config->Output_Directory, config->output_robustness_suffix, config->output_4D_suffix);
	Merge_Sparse_Files(InPath, InFile, config->TotalNbrSpots, OutPath);
      }
      if(config->LET_Sparse_Output == 1){
        sprintf(InPath, "%stmp/Beamlet_", config->Output_Directory);
        sprintf(InFile, "Sparse_LET%s%s.txt", config->output_robustness_suffix, config->output_4D_suffix);
        sprintf(OutPath, "%sSparse_LET%s%s.txt", config->Output_Directory, config->output_robustness_suffix, config->output_4D_suffix);
	Merge_Sparse_Files(InPath, InFile, config->TotalNbrSpots, OutPath);
      }
    }
  
    
  }
}
