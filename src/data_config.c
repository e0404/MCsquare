/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_config.h"


DATA_config_dictionary *Init_Config(DATA_config *config){

  unsigned int Num_Config_Tags = 63;

  DATA_config_dictionary *config_dictionary = (DATA_config_dictionary*) malloc(Num_Config_Tags * sizeof(DATA_config_dictionary));

  Add_uint_Config_element("_Internal_Num_Config_Tags", &config_dictionary[0], &config->Num_Config_Tags, 1, Num_Config_Tags, 0, UINT_MAX);
  Add_uint_Config_element("Num_Threads", &config_dictionary[1], &config->Num_Threads, 1, 0, 0, UINT_MAX);
  Add_ulong_Config_element("Num_Primaries", &config_dictionary[2], &config->Num_Primaries, 1, 10000000, 1, ULONG_MAX);
  Add_string_Config_element("CT_File", &config_dictionary[3], config->CT_File, 1, "CT.mhd", 1, 200);
  Add_ureal_Config_element("E_Cut_Pro", &config_dictionary[4], &config->Ecut_Pro, 1, 0.5, 0.001, 200);
  Add_ureal_Config_element("D_Max", &config_dictionary[5], &config->D_Max, 1, 0.2, 0.001, 100);
  Add_ureal_Config_element("Epsilon_Max", &config_dictionary[6], &config->Epsilon_Max, 1, 0.25, 0.001, 1.0);
  Add_ureal_Config_element("Te_Min", &config_dictionary[7], &config->Te_Min, 1, 0.05, 0.001, 10.0);
  Add_bool_Config_element("Simulate_Secondary_Protons", &config_dictionary[8], &config->Simulate_Secondary_Protons, 1, 1);
  Add_bool_Config_element("Simulate_Secondary_Deuterons", &config_dictionary[9], &config->Simulate_Secondary_Deuterons, 1, 1);
  Add_bool_Config_element("Simulate_Secondary_Alphas", &config_dictionary[10], &config->Simulate_Secondary_Alphas, 1, 1);
  Add_bool_Config_element("Score_PromptGammas", &config_dictionary[11], &config->Score_PromptGammas, 1, 0);
  Add_ureal_Config_element("PG_LowEnergyCut", &config_dictionary[12], &config->PG_LowEnergyCut, 1, 0.0, 0.0, 1000.0);
  Add_ureal_Config_element("PG_HighEnergyCut", &config_dictionary[13], &config->PG_HighEnergyCut, 1, 50.0, 0.001, 1000.0);
  Add_uint_Config_element("PG_Spectrum_NumBin", &config_dictionary[14], &config->PG_Spectrum_NumBin, 1, 150, 1, UINT_MAX);
  Add_ureal_Config_element("PG_Spectrum_Binning", &config_dictionary[15], &config->PG_Spectrum_Binning, 1, 0.1, 0.001, 1000.0);
  Add_string_Config_element("BDL_Machine_Parameter_File", &config_dictionary[16], config->BDL_machine, 1, "BDL.txt", 1, 200);
  Add_string_Config_element("BDL_Plan_File", &config_dictionary[17], config->BDL_plan, 1, "Plan.txt", 1, 200);
  Add_bool_Config_element("Energy_ASCII_Output", &config_dictionary[18], &config->Energy_ASCII_Output, 1, 0);
  Add_bool_Config_element("Energy_MHD_Output", &config_dictionary[19], &config->Energy_MHD_Output, 1, 0);
  Add_bool_Config_element("Energy_Sparse_Output", &config_dictionary[20], &config->Energy_Sparse_Output, 1, 0);
  Add_bool_Config_element("Dose_ASCII_Output", &config_dictionary[21], &config->Dose_ASCII_Output, 1, 0);
  Add_bool_Config_element("Dose_MHD_Output", &config_dictionary[22], &config->Dose_MHD_Output, 1, 1);
  Add_bool_Config_element("Dose_Sparse_Output", &config_dictionary[23], &config->Dose_Sparse_Output, 1, 0);
  Add_bool_Config_element("Densities_Output", &config_dictionary[24], &config->Densities_Output, 1, 0);
  Add_bool_Config_element("Materials_Output", &config_dictionary[25], &config->Materials_Output, 1, 0);
  Add_string_Config_element("HU_Density_Conversion_File", &config_dictionary[26], config->HU_Density_File, 1, "HU_Density_Conversion.txt", 1, 200);
  Add_string_Config_element("HU_Material_Conversion_File", &config_dictionary[27], config->HU_Material_File, 1, "HU_Material_Conversion.txt", 1, 200);
  Add_string_Config_element("Output_Directory", &config_dictionary[28], config->Output_Directory, 1, "Outputs", 1, 200);
  Add_bool_Config_element("Simulate_Nuclear_Interactions", &config_dictionary[29], &config->Simulate_Nuclear_Interactions, 1, 1);
  Add_bool_Config_element("Dose_Segmentation", &config_dictionary[30], &config->Dose_Segmentation, 1, 0);
  Add_ureal_Config_element("Density_Threshold_for_Segmentation", &config_dictionary[31], &config->Segmentation_Density_Threshold, 1, 0.01, 0.00001, 20);
  Add_uint_Config_element("RNG_Seed", &config_dictionary[32], &config->RNG_Seed, 1, 0, 0, UINT_MAX);
  Add_bool_Config_element("4D_Mode", &config_dictionary[33], &config->Simu_4D_Mode, 1, 0);
  Add_bool_Config_element("Robustness_Mode", &config_dictionary[34], &config->Robustness_Mode, 1, 0);
  Add_bool_Config_element("Beamlet_Mode", &config_dictionary[35], &config->Beamlet_Mode, 1, 0);
  Add_ureal_Config_element("Dose_Sparse_Threshold", &config_dictionary[36], &config->Dose_Sparse_Threshold, 1, 0.0, 0.0, 1e10);
  Add_ureal_Config_element("Energy_Sparse_Threshold", &config_dictionary[37], &config->Energy_Sparse_Threshold, 1, 0.0, 0.0, 1e10);
  Add_bool_Config_element("Compute_DVH", &config_dictionary[38], &config->Compute_DVH, 1, 0);
  Add_vec_ureal_Config_element("Systematic_Setup_Error", &config_dictionary[39], config->Systematic_Setup_Error, 1, 0.25, 0.25, 0.25, 0.0, 10.0);
  Add_vec_ureal_Config_element("Random_Setup_Error", &config_dictionary[40], config->Random_Setup_Error, 1, 0.1, 0.1, 0.1, 0.0, 10.0);
  Add_ureal_Config_element("Systematic_Range_Error", &config_dictionary[41], &config->Systematic_Range_Error, 1, 3.0, 0.0, 99.9);
  Add_bool_Config_element("Simulate_nominal_plan", &config_dictionary[42], &config->Simulate_nominal_plan, 1, 1);
  Add_bool_Config_element("4D_Dose_Accumulation", &config_dictionary[43], &config->Dose_4D_Accumulation, 1, 0);
  Add_bool_Config_element("LET_ASCII_Output", &config_dictionary[44], &config->LET_ASCII_Output, 1, 0);
  Add_bool_Config_element("LET_MHD_Output", &config_dictionary[45], &config->LET_MHD_Output, 1, 0);
  Add_bool_Config_element("LET_Sparse_Output", &config_dictionary[46], &config->LET_Sparse_Output, 1, 0);
  Add_ureal_Config_element("LET_Sparse_Threshold", &config_dictionary[47], &config->LET_Sparse_Threshold, 1, 0.0, 0.0, 1e10);
  Add_Enum_Config_element("LET_Calculation_Method", &config_dictionary[48], &config->LET_Calculation_Method, 1, 1, "DepositedEnergy;StopPow");
  Add_bool_Config_element("Beamlet_Parallelization", &config_dictionary[49], &config->Beamlet_Parallelization, 1, 0);
  Add_Enum_Config_element("Dose_to_Water_conversion", &config_dictionary[50], &config->DoseToWater, 1, 0, "Disabled;PostProcessing;OnlineSPR");
  Add_ureal_Config_element("MCS_const", &config_dictionary[51], &config->MCS_const, 1, 20.3, 0.001, 200);
  Add_Enum_Config_element("Scenario_selection", &config_dictionary[52], &config->Scenario_selection, 1, 0, "All;Random");
  Add_Enum_Config_element("Field_type", &config_dictionary[53], &config->Field_type, 1, 1, "Displacement;Velocity");
  Add_bool_Config_element("Create_4DCT_from_Ref", &config_dictionary[54], &config->Create_4DCT_from_Ref, 1, 0);
  Add_bool_Config_element("Create_Ref_from_4DCT", &config_dictionary[55], &config->Create_Ref_from_4DCT, 1, 0);
  Add_bool_Config_element("Dynamic_delivery", &config_dictionary[56], &config->Dynamic_delivery, 1, 0);
  Add_ureal_Config_element("Breathing_period", &config_dictionary[57], &config->Breathing_period, 1, 7.0, 0.01, 100);
  Add_ureal_Config_element("Systematic_Period_Error", &config_dictionary[58], &config->Systematic_Period_Error, 1, 5, 0.0, 200);
  Add_ureal_Config_element("Random_Period_Error", &config_dictionary[59], &config->Random_Period_Error, 1, 5, 0.0, 200);
  Add_ureal_Config_element("Systematic_Amplitude_Error", &config_dictionary[60], &config->Systematic_Amplitude_Error, 1, 5, 0.0, 200);
  Add_ureal_Config_element("Random_Amplitude_Error", &config_dictionary[61], &config->Random_Amplitude_Error, 1, 5, 0.0, 200);
  Add_bool_Config_element("Export_Beam_dose", &config_dictionary[62], &config->Export_Beam_dose, 1, 0);


  return config_dictionary;
}



void Add_bool_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned int *config, int use_default, unsigned int default_value){

  strcpy(config_dictionary->type, "bool");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    *config = default_value;
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.uint_adr = config;

  return;
}


void Add_string_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, char *config, int use_default, char *default_value, unsigned int min_size,  unsigned int max_size){

  strcpy(config_dictionary->type, "string");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    strcpy(config, default_value);
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.string_adr = config;
  config_dictionary->min_value.uint_value = min_size;
  config_dictionary->max_value.uint_value = max_size;

  return;
}


void Add_uint_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned int *config, int use_default, unsigned int default_value, unsigned int min_value,  unsigned int max_value){

  strcpy(config_dictionary->type, "uint");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    *config = default_value;
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.uint_adr = config;
  config_dictionary->min_value.uint_value = min_value;
  config_dictionary->max_value.uint_value = max_value;

  return;
}


void Add_ulong_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned long *config, int use_default, unsigned long default_value, unsigned long min_value,  unsigned long max_value){

  strcpy(config_dictionary->type, "ulong");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    *config = default_value;
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.ulong_adr = config;
  config_dictionary->min_value.ulong_value = min_value;
  config_dictionary->max_value.ulong_value = max_value;

  return;
}


void Add_ureal_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, VAR_DATA *config, int use_default, VAR_DATA default_value, VAR_DATA min_value,  VAR_DATA max_value){

  strcpy(config_dictionary->type, "ureal");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    *config = default_value;
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.real_adr = config;
  config_dictionary->min_value.real_value = min_value;
  config_dictionary->max_value.real_value = max_value;

  return;
}


void Add_vec_ureal_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, VAR_DATA *config, int use_default, VAR_DATA default_value_X, VAR_DATA default_value_Y, VAR_DATA default_value_Z, VAR_DATA min_value,  VAR_DATA max_value){

  strcpy(config_dictionary->type, "Vureal");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    config[0] = default_value_X;
    config[1] = default_value_Y;
    config[2] = default_value_Z;
  }

  strcpy(config_dictionary->Tag, Tag);
  config_dictionary->adress.real_adr = config;
  config_dictionary->min_value.real_value = min_value;
  config_dictionary->max_value.real_value = max_value;

  return;
}


void Add_Enum_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, int *config, int use_default, int default_value, char *List){

  strcpy(config_dictionary->type, "enum");
  config_dictionary->is_defined = 0;

  if(use_default == 0) config_dictionary->is_default = 0;
  else{
    config_dictionary->is_default = 1;
    *config = default_value;
  }

  strcpy(config_dictionary->Tag, Tag);
  strcpy(config_dictionary->List, List);
  config_dictionary->adress.int_adr = config;

  return;
}



int Parse_Config(DATA_config *config, char *file_name){

  FILE *file = NULL;
  char read[500], list[100], *read_token, *read_list, *save_token, *save_list;
  int i, j, k, parsed;

  DATA_config_dictionary *config_dictionary = Init_Config(config);
  
  file = fopen(file_name, "r");
  if(file == NULL){
    printf("\n\n ERROR: Configuration file \"%s\" not found! \n\n", file_name);
    return 1;
  }

  while (fgets(read, 500, file) != NULL){
    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok_r(read, "#", &save_token);

    read_token = strtok_r(read, " \t\r\n", &save_token);
    if(read_token == NULL) continue;

    parsed = 0;

    for(i = 0; i < config->Num_Config_Tags; i++){
      if(strcmp(read_token, config_dictionary[i].Tag) == 0){

	read_token = strtok_r(NULL, " \t\r\n", &save_token);

	///// Bool /////
        if(strcmp(config_dictionary[i].type, "bool") == 0){
	  if(!isBoolean(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"\n\n", read_token, config_dictionary[i].Tag, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else{
	    *config_dictionary[i].adress.uint_adr = getBoolean(read_token);
	    config_dictionary[i].is_defined = 1;
	  }
        }

	///// String /////
        else if(strcmp(config_dictionary[i].type, "string") == 0){
	  if(strlen(read_token) < config_dictionary[i].min_value.uint_value || strlen(read_token) > config_dictionary[i].max_value.uint_value){
	    printf("\n Error: %s value must contain between %u and %u characters in \"%s\"\n\n", config_dictionary[i].Tag, config_dictionary[i].min_value.uint_value, config_dictionary[i].max_value.uint_value, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else{
	    strcpy(config_dictionary[i].adress.string_adr, read_token);
	    config_dictionary[i].is_defined = 1;
	  }
        }

	///// UINT /////
        else if(strcmp(config_dictionary[i].type, "uint") == 0){
	  if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"\n\n", read_token, config_dictionary[i].Tag, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else if((unsigned)atoi(read_token) < config_dictionary[i].min_value.uint_value || (unsigned)atoi(read_token) > config_dictionary[i].max_value.uint_value){
	    printf("\n Error: %s value must be between %u and %u in \"%s\"\n\n", config_dictionary[i].Tag, config_dictionary[i].min_value.uint_value, config_dictionary[i].max_value.uint_value, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else{
	    *config_dictionary[i].adress.uint_adr = (unsigned)atoi(read_token);
	    config_dictionary[i].is_defined = 1;
	  }
        }

	///// ULONG /////
        else if(strcmp(config_dictionary[i].type, "ulong") == 0){
	  if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"\n\n", read_token, config_dictionary[i].Tag, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else if((unsigned long)atof(read_token) < config_dictionary[i].min_value.ulong_value || (unsigned long)atof(read_token) > config_dictionary[i].max_value.ulong_value){
	    printf("\n Error: %s value must be between %lu and %lu in \"%s\"\n\n", config_dictionary[i].Tag, config_dictionary[i].min_value.ulong_value, config_dictionary[i].max_value.ulong_value, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else{
	    *config_dictionary[i].adress.ulong_adr = (unsigned long)atof(read_token);
	    config_dictionary[i].is_defined = 1;
	  }
        }

	///// UReal /////
        else if(strcmp(config_dictionary[i].type, "ureal") == 0){
	  if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"\n\n", read_token, config_dictionary[i].Tag, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else if((VAR_DATA)atof(read_token) < config_dictionary[i].min_value.real_value || (VAR_DATA)atof(read_token) > config_dictionary[i].max_value.real_value){
	    printf("\n Error: %s value must be between %f and %f in \"%s\"\n\n", config_dictionary[i].Tag, config_dictionary[i].min_value.real_value, config_dictionary[i].max_value.real_value, file_name);
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
	  else{
	    *config_dictionary[i].adress.real_adr = (VAR_DATA)atof(read_token);
	    config_dictionary[i].is_defined = 1;
	  }
        }

	///// Vector UReal /////
        else if(strcmp(config_dictionary[i].type, "Vureal") == 0){
	  for(j=0; j<3; j++){
	    if(!isUnsignedFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"\n\n", read_token, config_dictionary[i].Tag, file_name);
	      fclose(file);
  	      free(config_dictionary);
	      return 1;
	    }
	    else if((VAR_DATA)atof(read_token) < config_dictionary[i].min_value.real_value || (VAR_DATA)atof(read_token) > config_dictionary[i].max_value.real_value){
	      printf("\n Error: %s value must be between %f and %f in \"%s\"\n\n", config_dictionary[i].Tag, config_dictionary[i].min_value.real_value, config_dictionary[i].max_value.real_value, file_name);
	      fclose(file);
  	      free(config_dictionary);
	      return 1;
	    }
	    else{
	      config_dictionary[i].adress.real_adr[j] = (VAR_DATA)atof(read_token);
	      if(j==2) config_dictionary[i].is_defined = 1;
	    }
	    read_token = strtok_r(NULL, " \t", &save_token);
	  }
        }

	///// Enum /////
        else if(strcmp(config_dictionary[i].type, "enum") == 0){
	  strcpy(list, config_dictionary[i].List);
	  read_list = strtok_r(list, ";", &save_list);
          j = -1;
	  k = -1;
	  while(read_list != NULL){
	    j += 1;
	    if(strcmp(read_list, read_token) == 0){
	      k = j;
	      *config_dictionary[i].adress.int_adr = k;
	      config_dictionary[i].is_defined = 1;
	      break;
	    }
	    read_list = strtok_r(NULL, ";", &save_list);
	  }
	  if(k == -1){
	    printf("\n Error: \"%s\" is not a valid value for %s Tag in \"%s\"", read_token, config_dictionary[i].Tag, file_name);
	    printf("\n Possible values are: \n");
	    strcpy(list, config_dictionary[i].List);
	    read_list = strtok_r(list, ";", &save_list);
	    while(read_list != NULL){
	      printf("\t%s\n", read_list);
	      read_list = strtok_r(NULL, ";", &save_list);
	    }
	    printf("\n\n");
	    fclose(file);
  	    free(config_dictionary);
	    return 1;
	  }
        }

	else{
	  printf("\n Error: Variable type \"%s\" is not supported for %s Tag in \"%s\"\n\n", config_dictionary[i].type, config_dictionary[i].Tag, file_name);
	  fclose(file);
  	  free(config_dictionary);
	  return 1;
	}
	
	parsed = 1;
	break;
      }
    } // end for

    if(parsed == 0){
	printf("\n Warning: Unknown tag \"%s\" in \"%s\"\n\n", read_token, file_name);
    }

  } // end while


  fclose(file);


  for(i = 0; i < config->Num_Config_Tags; i++){
    if(config_dictionary[i].is_default == 0 && config_dictionary[i].is_defined == 0){
      printf("\n Error: %s value not found.  It must be defined in \"%s\"\n\n", config_dictionary[i].Tag, file_name);
      free(config_dictionary);
      return 1;
    }
  }

  if(config->Output_Directory[strlen(config->Output_Directory)-1] != '/') strcat(config->Output_Directory, "/");

  CreateDir(config->Output_Directory);

  config->Particle_Generated_outside = 0;
  config->RangeShifter_enabled = 0;
  config->Water_Material_ID = 17;
  config->Current_Breathing_amplitude = 1.0;

  if(config->LET_ASCII_Output == 1 || config->LET_MHD_Output == 1 || config->LET_Sparse_Output == 1) config->Score_LET = 1;
  else config->Score_LET = 0;

  free(config_dictionary);

  return 0;

}

void display_config(DATA_config *config){

if(config->Num_Threads == 0) printf("\n\nNum threads = auto (%u)\n", omp_get_num_procs());
else printf("\n\nNum threads = %u \n", config->Num_Threads);
printf("Num primaries = %lu \n", config->Num_Primaries);
printf("RNG_Seed = %u \n", config->RNG_Seed);
printf("Ecut_Pro = %f \n", config->Ecut_Pro);
printf("D_Max = %f \n", config->D_Max);
printf("Epsilon_Max = %f \n", config->Epsilon_Max);
printf("Te_Min = %f \n\n", config->Te_Min);

printf("CT File = %s \n", config->CT_File);
printf("HU_Density_Conversion_File = %s \n", config->HU_Density_File);
printf("HU_Material_Conversion_File = %s \n", config->HU_Material_File);
printf("BDL machine parameter file = %s \n", config->BDL_machine);
printf("BDL plan file = %s \n\n", config->BDL_plan);

printf("Simulate_Nuclear_Interactions = %u \n", config->Simulate_Nuclear_Interactions);
printf("Simulate_Secondary_Protons = %u \n", config->Simulate_Secondary_Protons);
printf("Simulate_Secondary_Deuterons = %u \n", config->Simulate_Secondary_Deuterons);
printf("Simulate_Secondary_Alphas = %u \n\n", config->Simulate_Secondary_Alphas);

printf("4D_Mode = %u \n", config->Simu_4D_Mode);
printf("Robustness_Mode = %u \n", config->Robustness_Mode);
printf("Beamlet_Mode = %u \n\n", config->Beamlet_Mode);

printf("4D_Dose_Accumulation = %u \n\n", config->Dose_4D_Accumulation);
printf("Field_type = %d \n", config->Field_type);
printf("Create_4DCT_from_Ref = %u \n", config->Create_4DCT_from_Ref);
printf("Create_Ref_from_4DCT = %u \n", config->Create_Ref_from_4DCT);
printf("Dynamic_delivery = %u \n", config->Dynamic_delivery);
printf("Breathing_period = %f \n\n", config->Breathing_period);

printf("Simulate_nominal_plan = %u \n", config->Simulate_nominal_plan);
printf("Scenario_selection = %d \n", config->Scenario_selection);
printf("Systematic_Setup_Error = %f %f %f\n", config->Systematic_Setup_Error[0], config->Systematic_Setup_Error[1], config->Systematic_Setup_Error[2]);
printf("Random_Setup_Error = %f %f %f\n", config->Random_Setup_Error[0], config->Random_Setup_Error[1], config->Random_Setup_Error[2]);
printf("Systematic_Range_Error = %f \n", config->Systematic_Range_Error);
printf("Systematic_Amplitude_Error = %f \n", config->Systematic_Amplitude_Error);
printf("Random_Amplitude_Error = %f \n", config->Random_Amplitude_Error);
printf("Systematic_Period_Error = %f \n", config->Systematic_Period_Error);
printf("Random_Period_Error = %f \n\n", config->Random_Period_Error);

printf("Beamlet_Parallelization = %u \n", config->Beamlet_Parallelization);

printf("Output_Directory = %s \n\n", config->Output_Directory);

printf("Energy_ASCII_Output = %u \n", config->Energy_ASCII_Output);
printf("Energy_MHD_Output = %u \n", config->Energy_MHD_Output);
printf("Energy_Sparse_Output = %u \n", config->Energy_Sparse_Output);
printf("Dose_ASCII_Output = %u \n", config->Dose_ASCII_Output);
printf("Dose_MHD_Output = %u \n", config->Dose_MHD_Output);
printf("Dose_Sparse_Output = %u \n", config->Energy_Sparse_Output);
printf("LET_ASCII_Output = %u \n", config->LET_ASCII_Output);
printf("LET_MHD_Output = %u \n", config->LET_MHD_Output);
printf("LET_Sparse_Output = %u \n\n", config->LET_Sparse_Output);

printf("Densities_Output = %u \n", config->Densities_Output);
printf("Materials_Output = %u \n\n", config->Materials_Output);

printf("Dose_Sparse_Threshold = %f \n", config->Dose_Sparse_Threshold);
printf("Energy_Sparse_Threshold = %f \n", config->Energy_Sparse_Threshold);
printf("LET_Sparse_Threshold = %f \n\n", config->LET_Sparse_Threshold);

printf("Compute_DVH = %u \n\n", config->Compute_DVH);

printf("Score_PromptGammas = %u \n", config->Score_PromptGammas);
printf("PG_LowEnergyCut = %f \n", config->PG_LowEnergyCut);
printf("PG_HighEnergyCut = %f \n", config->PG_HighEnergyCut);
printf("PG_Spectrum_NumBin = %u \n", config->PG_Spectrum_NumBin);
printf("PG_Spectrum_Binning = %f \n\n", config->PG_Spectrum_Binning);

printf("Score_LET = %u \n", config->Score_LET);
printf("LET_Calculation_Method = %d \n\n", config->LET_Calculation_Method);

printf("Export_Beam_dose = %u \n", config->Export_Beam_dose);
printf("Dose_to_Water_conversion = %d \n\n", config->DoseToWater);

printf("Dose_Segmentation = %u \n", config->Dose_Segmentation);
printf("Segmentation_Density_Threshold = %f \n\n", config->Segmentation_Density_Threshold);
}



