/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_nuclear.h"

int Init_Nuclear_Data(Materials *material, int Nbr_Materials, DATA_config *config){

  int i, j, fail;

  for(i=0; i<Nbr_Materials; i++){
    if(material[i].Nuclear_data_type == Mixture){
      for(j=0; j < material[i].NbrComponents; j++){
        if(material[material[i].Mixture_Components_label[j]].Nuclear_data_type != Proton_Proton && material[material[i].Mixture_Components_label[j]].Nuclear_data_type != ICRU){
          printf("\n Error: The Nuclear_Data value for material %s must be \"ICRU\" or \"proton-proton\" because it is used in Mixture %s.\n\n", material[material[i].Mixture_Components_label[j]].Name, material[i].Name);
	  return -1;
        }
      }
    }

    else if(material[i].Nuclear_data_type == ICRU){

      fail = read_Nuclear_Elastic_ICRU(&material[i], config);
      if(fail != 0) return -1;

      fail = read_Nuclear_Inelastic_ICRU(&material[i], config);
      if(fail != 0) return -1;

  if(config->Score_PromptGammas == 1){
      fail = read_PromptGamma_ICRU(&material[i], config);
      if(fail != 0) return -1;
  }

    }
  }

  Interp_Nuclear_Cross_section(material, Nbr_Materials);

  return 0;

}



int read_Nuclear_Elastic_ICRU(Materials *material, DATA_config *config){
  FILE *file = NULL;
  char read[500], file_name[100], *read_token;
  int i, Nbr_Energy = 0;

  strcpy(file_name, config->Materials_Dir);
  strcat(file_name, material->Name);
  strcat(file_name, "/ICRU_Nuclear_elastic.dat");
  file = fopen(file_name, "r");

  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
    return -1;
  }

  // Count the number of Energy in the file
  while (fgets(read, 500, file) != NULL){
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      Nbr_Energy++;
    }
  }

  material->Nbr_Elastic_Energy = Nbr_Energy;

  material->Nuclear_Elastic = (DATA_Nuclear_Elastic*) malloc(Nbr_Energy * sizeof(DATA_Nuclear_Elastic));
  material->Elastic_Energy_List = (VAR_DATA*) malloc(Nbr_Energy * sizeof(VAR_DATA));

  // Go to the beginning of the file
  rewind(file);
  Nbr_Energy = 0;

  // Data Processing
  while (fgets(read, 500, file) != NULL){
    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Energy in \"%s\"\n\n", read_token, file_name);
	return -1;
      }
      Nbr_Energy++;
      material->Elastic_Energy_List[Nbr_Energy-1] = atof(read_token);
    }

    else if(strcmp(read_token, "Cross_section") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Cross_section\" in \"%s\"\n\n", file_name);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Cross_section in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Elastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      material->Nuclear_Elastic[Nbr_Energy-1].Cross_section = N_AVO * atof(read_token) * 1e-24 / material->A ;
    }

    else if(strcmp(read_token, "Differential_cross_section") == 0){
      i = 0;
      while (fgets(read, 500, file) != NULL && i < 36){
	// on ignore les commentaires
    	if(read[0] == '#') continue;
    	strtok(read, "#");
    	read_token = strtok(read, " \t\r\n");
    	if(read_token == NULL) continue;

	if(atof(read_token) != (i+1)*5.0){
	  printf("\n\n Error: \"%s\" is not a valid angle for Differential_cross_section in \"%s\" for Energy = %f\n", read_token, file_name, material->Elastic_Energy_List[Nbr_Energy-1]);
	  printf("The angle must be a multible of 5 between 5 and 180 degrees\n\n");
	  return -1;
	}
	
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Differential_cross_section in \"%s\" for Energy = %f & Angle = %f\n\n", read_token, file_name, material->Elastic_Energy_List[Nbr_Energy-1], (i+1)*5.0);
	  return -1;
        }
	if(i == 0) material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i] = atof(read_token) * (SolidAngle(7.5) - SolidAngle(5));
	else if(i > 0 && i < 35) material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i] = material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i-1] + atof(read_token) * (SolidAngle((i+1)*5 + 2.5) - SolidAngle((i+1)*5 - 2.5));
	else if(i == 35) material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i] = material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i-1] + atof(read_token) * (4*M_PI - SolidAngle(177.5));
	i++;
      }
      
      // Normalisation de la fonction de probabilité cumulative
      for(i=0; i<36; i++){
	material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i] = material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[i] / material->Nuclear_Elastic[Nbr_Energy-1].Diff_Cross_section[35];
      }
    }
  }
  
  return 0;
}




int read_Nuclear_Inelastic_ICRU(Materials *material, DATA_config *config){

  static const double ICRU_angles[13] = { 0, 10, 20, 30, 40, 50, 60, 70, 90, 110, 130, 150, 180 };

  FILE *file = NULL;
  char read[500], file_name[100], *read_token, secondary[25];
  int Nbr_Energy = 0, Nbr_rows, i, j;
  VAR_DATA *pt_Energy, *pt_cross_section, *pt_diff_cross_section;

  strcpy(file_name, config->Materials_Dir);
  strcat(file_name, material->Name);
  strcat(file_name, "/ICRU_Nuclear_inelastic.dat");
  file = fopen(file_name, "r");

  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
    return -1;
  }

  // Count the number of Energy in the file
  while (fgets(read, 500, file) != NULL){
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      Nbr_Energy++;
    }
  }

  material->Nbr_Inelastic_Energy = Nbr_Energy;

  material->Nuclear_Inelastic = (DATA_Nuclear_Inelastic*) malloc(Nbr_Energy * sizeof(DATA_Nuclear_Inelastic));
  material->Inelastic_Energy_List = (VAR_DATA*) malloc(Nbr_Energy * sizeof(VAR_DATA));

  // Go to the beginning of the file
  rewind(file);
  Nbr_Energy = 0;

  while (fgets(read, 500, file) != NULL){
    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Energy in \"%s\"\n\n", read_token, file_name);
	return -1;
      }
      Nbr_Energy++;
      material->Inelastic_Energy_List[Nbr_Energy-1] = atof(read_token);

      material->Nuclear_Inelastic[Nbr_Energy-1].P_Nbr_data = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].Proton_Mult = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].P_Energy_list = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].P_D_Cross_section = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].P_DD_Cross_section = NULL;

      material->Nuclear_Inelastic[Nbr_Energy-1].A_Nbr_data = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].Alpha_Mult = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].A_Energy_list = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].A_D_Cross_section = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].A_DD_Cross_section = NULL;

      material->Nuclear_Inelastic[Nbr_Energy-1].D_Nbr_data = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].Deuteron_Mult = 0;
      material->Nuclear_Inelastic[Nbr_Energy-1].D_Energy_list = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].D_D_Cross_section = NULL;
      material->Nuclear_Inelastic[Nbr_Energy-1].D_DD_Cross_section = NULL;
    }

    else if(strcmp(read_token, "Cross_section") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Cross_section\" in \"%s\"\n\n", file_name);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Cross_section in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      material->Nuclear_Inelastic[Nbr_Energy-1].Cross_section = N_AVO * atof(read_token) * 1e-24 / material->A ;
    }

    else if(strcmp(read_token, "Multiplicity") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Multiplicity\" in \"%s\"\n\n", file_name);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(strcmp(read_token, "proton") == 0){
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Multiplicity proton in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	material->Nuclear_Inelastic[Nbr_Energy-1].Proton_Mult = atof(read_token);
      }
      else if(strcmp(read_token, "alpha") == 0){
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Multiplicity alpha in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	material->Nuclear_Inelastic[Nbr_Energy-1].Alpha_Mult = atof(read_token);
      }
      else if(strcmp(read_token, "deuteron") == 0){
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Multiplicity deuteron in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	material->Nuclear_Inelastic[Nbr_Energy-1].Deuteron_Mult = atof(read_token);
      }
      else{
	printf("\n\n Error: \"%s\" is not a valid value for Multiplicity in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }
    }

    else if(strcmp(read_token, "Energy_fraction_recoils") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Energy_fraction_recoils\" in \"%s\"\n\n", file_name);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Energy_fraction_recoils in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      material->Nuclear_Inelastic[Nbr_Energy-1].Energy_fraction_recoils = atof(read_token);
    }

    else if(strcmp(read_token, "Differential_cross_section") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Differential_cross_section\" in \"%s\"\n\n", file_name);
	return -1;
      }

      read_token = strtok(NULL, " \t\r\n");
      strcpy (secondary, read_token);

      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedInt(read_token)){
	printf("\n\n Error: \"%s\" is not a valid \"Number of rows\" for \"Differential_cross_section\" in \"%s\" for Energy = %f\n\n", read_token, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      Nbr_rows = atoi(read_token);

      if(strcmp(secondary, "proton") == 0){
	material->Nuclear_Inelastic[Nbr_Energy-1].P_Nbr_data = Nbr_rows;
	material->Nuclear_Inelastic[Nbr_Energy-1].P_Energy_list = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].P_D_Cross_section = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].P_DD_Cross_section = (VAR_DATA*) malloc(Nbr_rows*13 * sizeof(VAR_DATA));

	pt_Energy = material->Nuclear_Inelastic[Nbr_Energy-1].P_Energy_list;
	pt_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].P_D_Cross_section;
	pt_diff_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].P_DD_Cross_section;
      }
      else if(strcmp(secondary, "alpha") == 0){
	material->Nuclear_Inelastic[Nbr_Energy-1].A_Nbr_data = Nbr_rows;
	material->Nuclear_Inelastic[Nbr_Energy-1].A_Energy_list = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].A_D_Cross_section = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].A_DD_Cross_section = (VAR_DATA*) malloc(Nbr_rows*13 * sizeof(VAR_DATA));

	pt_Energy = material->Nuclear_Inelastic[Nbr_Energy-1].A_Energy_list;
	pt_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].A_D_Cross_section;
	pt_diff_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].A_DD_Cross_section;
      }
      else if(strcmp(secondary, "deuteron") == 0){
	material->Nuclear_Inelastic[Nbr_Energy-1].D_Nbr_data = Nbr_rows;
	material->Nuclear_Inelastic[Nbr_Energy-1].D_Energy_list = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].D_D_Cross_section = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->Nuclear_Inelastic[Nbr_Energy-1].D_DD_Cross_section = (VAR_DATA*) malloc(Nbr_rows*13 * sizeof(VAR_DATA));

	pt_Energy = material->Nuclear_Inelastic[Nbr_Energy-1].D_Energy_list;
	pt_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].D_D_Cross_section;
	pt_diff_cross_section = material->Nuclear_Inelastic[Nbr_Energy-1].D_DD_Cross_section;
      }
      else{
	printf("\n\n Error: \"%s\" is not a valid value for Differential_cross_section in \"%s\" for Energy = %f\n\n", secondary, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	return -1;
      }

      i = 0;
      while (fgets(read, 500, file) != NULL && i < Nbr_rows){

	// on ignore les commentaires
    	if(read[0] == '#') continue;
    	strtok(read, "#");
    	read_token = strtok(read, " \t\r\n");
    	if(read_token == NULL) continue;
	
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid Energy for the %s Differential cross section in \"%s\" for Energy = %f\n\n", read_token, secondary, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	pt_Energy[i] = atof(read_token);

	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for the %s Integrated Differential cross section in \"%s\" for Energy = %f\n\n", read_token, secondary, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	pt_cross_section[i] = atof(read_token);

	for(j=0; j<13; j++){
	  read_token = strtok(NULL, " \t\r\n");

	  if(!isUnsignedFloat(read_token)){
	    printf("\n\n %d-%d Error: \"%s\" is not a valid value for the %s Differential cross section in \"%s\" for Energy = %f\n\n", i, j, read_token, secondary, file_name, material->Inelastic_Energy_List[Nbr_Energy-1]);
	    return -1;
          }

          if(j == 0) pt_diff_cross_section[i*13+j] = atof(read_token) * SolidAngle(5);
          else if(j == 12) pt_diff_cross_section[i*13+j] = atof(read_token) * (4*M_PI - SolidAngle(165));
          else pt_diff_cross_section[i*13+j] = atof(read_token) * (SolidAngle((ICRU_angles[j]+ICRU_angles[j+1])/2) - SolidAngle((ICRU_angles[j-1]+ICRU_angles[j])/2));

	}
	i++;
      }
    }

  }

  return 0;
}


int read_PromptGamma_ICRU(Materials *material, DATA_config *config){
  FILE *file = NULL;
  char read[500], file_name[100], *read_token, secondary[25];
  int i, Nbr_Energy = 0, Nbr_rows;
  VAR_DATA *pt_Energy, *pt_diff_cross_section;

  strcpy(file_name, config->Materials_Dir);
  strcat(file_name, material->Name);
  strcat(file_name, "/ICRU_PromptGamma.dat");
  file = fopen(file_name, "r");

  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
    return -1;
  }

  // Count the number of Energy in the file
  while (fgets(read, 500, file) != NULL){
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      Nbr_Energy++;
    }
  }

  material->Nbr_PG_Energy = Nbr_Energy;

  material->PG_data = (DATA_PG*) malloc(Nbr_Energy * sizeof(DATA_PG));
  material->PG_Energy_List = (VAR_DATA*) malloc(Nbr_Energy * sizeof(VAR_DATA));

  // Go to the beginning of the file
  rewind(file);
  Nbr_Energy = 0;

  // Data Processing
  while (fgets(read, 500, file) != NULL){
    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "Energy") == 0){
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Energy in \"%s\"\n\n", read_token, file_name);
	return -1;
      }
      Nbr_Energy++;
      material->PG_Energy_List[Nbr_Energy-1] = atof(read_token);
      material->PG_data[Nbr_Energy-1].Energy_list = NULL;
      material->PG_data[Nbr_Energy-1].Diff_Cross_section = NULL;
    }

    else if(strcmp(read_token, "Multiplicity") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Multiplicity\" in \"%s\"\n\n", file_name);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(strcmp(read_token, "PromptGamma") != 0){
	printf("\n\n Error: \"%s\" is not a valid secondary particle for Multiplicity in \"%s\" for Energy = %f\n\n", read_token, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedFloat(read_token)){
	printf("\n\n Error: \"%s\" is not a valid value for Multiplicity in \"%s\" for Energy = %f\n\n", read_token, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      material->PG_data[Nbr_Energy-1].Multiplicity = atof(read_token);
    }

    else if(strcmp(read_token, "Differential_cross_section") == 0){
      if(Nbr_Energy == 0){
	printf("\n\n Error: \"Energy\" must be defined before \"Differential_cross_section\" in \"%s\"\n\n", file_name);
	return -1;
      }

      read_token = strtok(NULL, " \t\r\n");
      strcpy (secondary, read_token);

      read_token = strtok(NULL, " \t\r\n");
      if(!isUnsignedInt(read_token)){
	printf("\n\n Error: \"%s\" is not a valid \"Number of rows\" for \"Differential_cross_section\" in \"%s\" for Energy = %f\n\n", read_token, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	return -1;
      }
      Nbr_rows = atoi(read_token);

      if(strcmp(secondary, "PromptGamma") == 0){
	material->PG_data[Nbr_Energy-1].Nbr_data = Nbr_rows;
	material->PG_data[Nbr_Energy-1].Energy_list = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));
	material->PG_data[Nbr_Energy-1].Diff_Cross_section = (VAR_DATA*) malloc(Nbr_rows * sizeof(VAR_DATA));

	pt_Energy = material->PG_data[Nbr_Energy-1].Energy_list;
	pt_diff_cross_section = material->PG_data[Nbr_Energy-1].Diff_Cross_section;
      }
      else{
	printf("\n\n Error: \"%s\" is not a valid value for Differential_cross_section in \"%s\" for Energy = %f\n\n", secondary, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	return -1;
      }

      i = 0;
      while (fgets(read, 500, file) != NULL && i < Nbr_rows){

	// on ignore les commentaires
    	if(read[0] == '#') continue;
    	strtok(read, "#");
    	read_token = strtok(read, " \t\r\n");
    	if(read_token == NULL) continue;
	
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid Energy for the %s Differential cross section in \"%s\" for Energy = %f\n\n", read_token, secondary, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	pt_Energy[i] = atof(read_token);

	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for the %s Integrated Differential cross section in \"%s\" for Energy = %f\n\n", read_token, secondary, file_name, material->PG_Energy_List[Nbr_Energy-1]);
	  return -1;
        }
	pt_diff_cross_section[i] = atof(read_token);

	i++;
      }
    }
  }
  
  return 0;
}


void Interp_Nuclear_Cross_section(Materials *material, int Nbr_Materials){

  int Nbr_bin = ceil(250.0 / INTERP_BIN);
  VAR_DATA tmp;

  int i, j, k, index, label;
  for(i=0; i<Nbr_Materials; i++){

    material[i].Interp_Total_Nuclear_Cross_Section = (VAR_DATA*) malloc(Nbr_bin * sizeof(VAR_DATA));

    if(material[i].Nuclear_data_type == Proton_Proton){
      for(j=0; j<Nbr_bin; j++){
	if(j*INTERP_BIN < 10.0) material[i].Interp_Total_Nuclear_Cross_Section[j] = 0.0;
	else material[i].Interp_Total_Nuclear_Cross_Section[j] = (0.315*pow((j*INTERP_BIN), -1.126) + 3.78e-6 * (j*INTERP_BIN))/0.1119;
      }
    }

    else if(material[i].Nuclear_data_type == ICRU){
      for(j=0; j<Nbr_bin; j++){

	material[i].Interp_Total_Nuclear_Cross_Section[j] = 0;

        index = Binary_Search((j*INTERP_BIN), material[i].Elastic_Energy_List, material[i].Nbr_Elastic_Energy);
	if(index < 0) index = 0;
	else if(index >= material[i].Nbr_Elastic_Energy-1) index = material[i].Nbr_Elastic_Energy - 2;
	
	tmp = (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
						material[i].Elastic_Energy_List[index], 
						material[i].Elastic_Energy_List[index+1], 
						material[i].Nuclear_Elastic[index].Cross_section, 
						material[i].Nuclear_Elastic[index+1].Cross_section);

	if(tmp > 0) material[i].Interp_Total_Nuclear_Cross_Section[j] += tmp;

	index = Binary_Search((j*INTERP_BIN), material[i].Inelastic_Energy_List, material[i].Nbr_Inelastic_Energy);
	if(index < 0) index = 0;
	else if(index >= material[i].Nbr_Inelastic_Energy-1) index = material[i].Nbr_Inelastic_Energy - 2;

	tmp = (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
						material[i].Inelastic_Energy_List[index], 
						material[i].Inelastic_Energy_List[index+1], 
						material[i].Nuclear_Inelastic[index].Cross_section, 
						material[i].Nuclear_Inelastic[index+1].Cross_section);

	material[i].Interp_Total_Nuclear_Cross_Section[j] += tmp;


/*
	if(j*INTERP_BIN <= material[i].Elastic_Energy_List[0] || j*INTERP_BIN <= material[i].Inelastic_Energy_List[0]) material[i].Interp_Total_Nuclear_Cross_Section[j] = 0.0;
        else{
	  index = Binary_Search((j*INTERP_BIN), material[i].Elastic_Energy_List, material[i].Nbr_Elastic_Energy);
	  material[i].Interp_Total_Nuclear_Cross_Section[j] = (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
											material[i].Elastic_Energy_List[index], 
											material[i].Elastic_Energy_List[index+1], 
											material[i].Nuclear_Elastic[index].Cross_section, 
											material[i].Nuclear_Elastic[index+1].Cross_section);

	  index = Binary_Search((j*INTERP_BIN), material[i].Inelastic_Energy_List, material[i].Nbr_Inelastic_Energy);
	  material[i].Interp_Total_Nuclear_Cross_Section[j] += (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
											material[i].Inelastic_Energy_List[index], 
											material[i].Inelastic_Energy_List[index+1], 
											material[i].Nuclear_Inelastic[index].Cross_section, 
											material[i].Nuclear_Inelastic[index+1].Cross_section);
	}
*/
      }
    }

    else if(material[i].Nuclear_data_type == Mixture){
      for(j=0; j<Nbr_bin; j++){
        material[i].Interp_Total_Nuclear_Cross_Section[j] = 0.0;

        for(k=0; k < material[i].NbrComponents; k++){
	  label = material[i].Mixture_Components_label[k];
	
	  if(material[label].Nuclear_data_type == Proton_Proton){
	    if(j*INTERP_BIN >= 10.0) 
	      material[i].Interp_Total_Nuclear_Cross_Section[j] += 	material[i].Mixture_Components_fraction[k] 
									* (0.315*pow((j*INTERP_BIN), -1.126) + 3.78e-6 * (j*INTERP_BIN))/0.1119;
	  }

	  else if(material[label].Nuclear_data_type == ICRU){
	    if(j*INTERP_BIN > 7){
	      index = Binary_Search((j*INTERP_BIN), material[label].Elastic_Energy_List, material[label].Nbr_Elastic_Energy);
	      material[i].Interp_Total_Nuclear_Cross_Section[j] += 	material[i].Mixture_Components_fraction[k] 
									* (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
												material[label].Elastic_Energy_List[index], 
												material[label].Elastic_Energy_List[index+1], 
												material[label].Nuclear_Elastic[index].Cross_section, 
												material[label].Nuclear_Elastic[index+1].Cross_section);

	      index = Binary_Search((j*INTERP_BIN), material[label].Inelastic_Energy_List, material[label].Nbr_Inelastic_Energy);
	      material[i].Interp_Total_Nuclear_Cross_Section[j] += 	material[i].Mixture_Components_fraction[k] 
									* (VAR_DATA)Linear_Interpolation(	(j*INTERP_BIN), 
												material[label].Inelastic_Energy_List[index], 
												material[label].Inelastic_Energy_List[index+1], 
												material[label].Nuclear_Inelastic[index].Cross_section, 
												material[label].Nuclear_Inelastic[index+1].Cross_section);
	    }
	  }

        }

      }
    }

  }

  return;
}


void Free_DATA_Nuclear_Inelastic(DATA_Nuclear_Inelastic *data, int Nbr_data){
  if (data == NULL) return;

  int i;
  for(i=0; i<Nbr_data; i++){
    if(data[i].P_Energy_list != NULL) free(data[i].P_Energy_list);
    if(data[i].P_D_Cross_section != NULL) free(data[i].P_D_Cross_section);
    if(data[i].P_DD_Cross_section != NULL) free(data[i].P_DD_Cross_section);

    if(data[i].A_Energy_list != NULL) free(data[i].A_Energy_list);
    if(data[i].A_D_Cross_section != NULL) free(data[i].A_D_Cross_section);
    if(data[i].A_DD_Cross_section != NULL) free(data[i].A_DD_Cross_section);

    if(data[i].D_Energy_list != NULL) free(data[i].D_Energy_list);
    if(data[i].D_D_Cross_section != NULL) free(data[i].D_D_Cross_section);
    if(data[i].D_DD_Cross_section != NULL) free(data[i].D_DD_Cross_section);
  }

  free(data);
}


void Free_DATA_PG(DATA_PG *data, int Nbr_data){
  if (data == NULL) return;

  int i;
  for(i=0; i<Nbr_data; i++){
    if(data[i].Energy_list != NULL) free(data[i].Energy_list);
    if(data[i].Diff_Cross_section != NULL) free(data[i].Diff_Cross_section);
  }

  free(data);
}



