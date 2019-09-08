/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_materials.h"

Materials *Init_materials(unsigned int *Nbr_Materials, unsigned int *Max_Components, DATA_config *config){

  FILE *file = NULL;
  char read[500], file_name[200], *GetEnvValue, *read_token;
  int current_component, fail;

  GetEnvValue = NULL;
  GetEnvValue = getenv("MCsquare_Materials_Dir");

  if(File_exists("Materials/list.dat") == 1){
    strcpy(config->Materials_Dir, "Materials/");
  }
  else if(GetEnvValue == NULL){
      printf("\n\n ERROR: Material database not found! \n\n");
      return NULL;
  }
  else{
    strcpy(config->Materials_Dir, GetEnvValue);
    strcat(config->Materials_Dir, "/list.dat");
    if(File_exists(config->Materials_Dir) == 1){
      strcpy(config->Materials_Dir, GetEnvValue);
      strcat(config->Materials_Dir, "/");
    }
    else{
      printf("\n\n ERROR: Material database not found! \n\n");
      return NULL;
    }
  }

  Materials *material = List_materials(Nbr_Materials, config);
  *Max_Components = 0;

  if(material == NULL) return NULL;


  int i = 0;
  for(i=1; i<*Nbr_Materials; i++){
    if(strlen(material[i].Name) == 0){
      printf("\n\n Warning : Material label %d is not defined \n\n", i);
      continue;
    }

    strcpy(file_name, config->Materials_Dir);
    strcat(file_name, material[i].Name);
    strcat(file_name, "/Material_Properties.dat");
    file = fopen(file_name, "r");

    current_component = 0;

    if(file == NULL){
      printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    while (fgets(read, 500, file) != NULL){
      // on ignore les commentaires
      if(read[0] == '#') continue;
      strtok(read, "#");

      read_token = strtok(read, " \t\r\n");
      if(read_token == NULL) continue;

      if(strcmp(read_token, "Name") == 0){
	
      }

      else if(strcmp(read_token, "Density") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Density in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
        }
	material[i].Density = atof(read_token);
      }

      else if(strcmp(read_token, "Atomic_Weight") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Atomic_Weight in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
        }
	material[i].A = atof(read_token);
      }

      else if(strcmp(read_token, "Molecular_Weight") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(!isUnsignedFloat(read_token)){
	  printf("\n\n Error: \"%s\" is not a valid value for Molecular_Weight in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
        }
	material[i].A = atof(read_token);
      }

      else if(strcmp(read_token, "Electron_Density") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(!isUnsignedFloat(read_token)){
	  printf("\n Error: \"%s\" is not a valid value for Electron_Density in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
        }
	material[i].N_el = atof(read_token);
      }

      else if(strcmp(read_token, "Radiation_Length") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(!isUnsignedFloat(read_token)){
	  printf("\n Error: \"%s\" is not a valid value for Radiation_Length in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
        }
	material[i].X0 = atof(read_token);
      }

      else if(strcmp(read_token, "Nuclear_Data") == 0){
	read_token = strtok(NULL, " \t\r\n");
        if(strcmp(read_token, "ICRU") == 0){
	  material[i].Nuclear_data_type = ICRU;
	  material[i].NbrComponents = 1;
        }
        else if(strcmp(read_token, "Mixture") == 0){
	  material[i].Nuclear_data_type = Mixture;
	  read_token = strtok(NULL, " \t\r\n");
	  if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for Nuclear_Data Mixture in \"%s\"\n\n", read_token, file_name);
	    Free_Materials_DATA(material, *Nbr_Materials);
	    return NULL;
	  }
          else{
	    material[i].NbrComponents = atoi(read_token);
	    if(*Max_Components < material[i].NbrComponents) *Max_Components = material[i].NbrComponents;
	    material[i].Mixture_Components_label = (int*)malloc(material[i].NbrComponents * sizeof(int));
	    material[i].Mixture_Components_fraction = (VAR_DATA*)malloc(material[i].NbrComponents * sizeof(VAR_DATA));
	  }
        }
        else if(strcmp(read_token, "proton-proton") == 0){
	  material[i].Nuclear_data_type = Proton_Proton;
	  material[i].NbrComponents = 1;
        }
	else{
	  printf("\n Error: \"%s\" is not a valid value for Nuclear_Data in \"%s\"\n\n", read_token, file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
	}
      }

      else if(strcmp(read_token, "Mixture_Component") == 0){
	if(material[i].Nuclear_data_type != Mixture){
	  printf("\n Error: \"Nuclear_Data Mixture\" must be declared before \"Mixture_Component\" in \"%s\"\n\n", file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
	}
	else if(current_component >= material[i].NbrComponents){
	  printf("\n Error: too much \"Mixture_Component\"!  You have to modify the \"Nuclear_Data Mixture\" value in \"%s\"\n\n", file_name);
	  Free_Materials_DATA(material, *Nbr_Materials);
	  return NULL;
	}
	else{
	  read_token = strtok(NULL, " \t\r\n");
	  if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Mixture_Component label in \"%s\"\n\n", read_token, file_name);
	    Free_Materials_DATA(material, *Nbr_Materials);
	    return NULL;
	  }
	  material[i].Mixture_Components_label[current_component] = atoi(read_token);
	  read_token = strtok(NULL, " \t\r\n");
	  read_token = strtok(NULL, " \t\r\n");
	  if(!isUnsignedFloat(read_token) || atof(read_token) < 0.0 || atof(read_token) > 100.0){
	    printf("\n Error: \"%s\" is not a valid value for Mixture_Component fraction in \"%s\"\n\n", read_token, file_name);
	    Free_Materials_DATA(material, *Nbr_Materials);
	    return NULL;
	  }
	  material[i].Mixture_Components_fraction[current_component] = atof(read_token)/100.0;
	  current_component++;
        }
      }

      else{
	printf("\n Warning: Unknown parameter \"%s\" in \"%s\"\n\n", read_token, file_name);
      }

    }
    
    fclose(file);

    if(material[i].Density == 0){
      printf("\n Error: Parameter \"Density\" is missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    if(material[i].N_el == 0){
      printf("\n Error: Parameter \"Electron_Density\" is missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    if(material[i].X0 == 0){
      printf("\n Error: Parameter \"Radiation_Length\" is missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    if(material[i].Nuclear_data_type == None){
      printf("\n Error: Parameter \"Nuclear_Data\" is missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    if((material[i].Nuclear_data_type == ICRU || material[i].Nuclear_data_type == Proton_Proton) && material[i].A == 0){
      printf("\n Error: Parameter \"Atomic_Weight\" is missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    if(material[i].Nuclear_data_type == Mixture && current_component != material[i].NbrComponents){
      printf("\n Error: Some \"Mixture_Component\" definition are missing in \"%s\"\n\n", file_name);
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

    // On divise la densité électronique par la densité massique pour pouvoir multiplier cette valeur par la densité locale dans la simulation.
    material[i].N_el = material[i].N_el / material[i].Density;

    fail = Read_Stop_Pow(material[i].Name, &material[i], config);
    if(fail !=  0){
      Free_Materials_DATA(material, *Nbr_Materials);
      return NULL;
    }

  }



  // Initialisation des données des interactions nucléaires
  fail = Init_Nuclear_Data(material, *Nbr_Materials, config);

  // Si l'initialisation nucléaire échoue, on stop le programme
  if(fail != 0){
    Free_Materials_DATA(material, *Nbr_Materials);
    return NULL;
  }

  // copy water material in label 0
    strcpy(material[0].Name, "water");
    material[0].A = material[config->Water_Material_ID].A;
    material[0].Density = material[config->Water_Material_ID].Density;
    material[0].N_el =  material[config->Water_Material_ID].N_el;
    material[0].X0 = material[config->Water_Material_ID].X0;
    material[0].SPR = material[config->Water_Material_ID].SPR;
    material[0].Nuclear_data_type = material[config->Water_Material_ID].Nuclear_data_type;
    material[0].NbrComponents = material[config->Water_Material_ID].NbrComponents;
    material[0].Mixture_Components_label = material[config->Water_Material_ID].Mixture_Components_label;
    material[0].Mixture_Components_fraction = material[config->Water_Material_ID].Mixture_Components_fraction;
    material[0].Nuclear_Elastic = material[config->Water_Material_ID].Nuclear_Elastic;
    material[0].Elastic_Energy_List = material[config->Water_Material_ID].Elastic_Energy_List;
    material[0].Nuclear_Inelastic = material[config->Water_Material_ID].Nuclear_Inelastic;
    material[0].Inelastic_Energy_List = material[config->Water_Material_ID].Inelastic_Energy_List;
    material[0].PG_data = material[config->Water_Material_ID].PG_data;
    material[0].PG_Energy_List = material[config->Water_Material_ID].PG_Energy_List;
    material[0].SP_Energy = material[config->Water_Material_ID].SP_Energy;
    material[0].Stop_Pow = material[config->Water_Material_ID].Stop_Pow;
    material[0].Interp_Total_Nuclear_Cross_Section = material[config->Water_Material_ID].Interp_Total_Nuclear_Cross_Section;

  
  // Compute SPR:
  const int ID_Energy = (int)floor(100/PSTAR_BIN);
  for(i=1; i<*Nbr_Materials; i++){
    if(strlen(material[i].Name) == 0) continue;
    if(i == config->Water_Material_ID) continue;
    material[i].SPR = material[i].Stop_Pow[ID_Energy] / material[config->Water_Material_ID].Stop_Pow[ID_Energy];
  }

  return material;
}



Materials *List_materials(unsigned int *Nbr_Materials, DATA_config *config){
  FILE *file = NULL;
  char read[500], file_name[200], *read_token;
  int label, max_label = 0;

  strcpy(file_name, config->Materials_Dir);
  strcat(file_name, "list.dat");

  file = fopen(file_name, "r");
  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
    return NULL;
  }

  // Get the max material label
  while (fgets(read, 500, file) != NULL){

    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");	

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL || !isUnsignedInt(read_token)) continue;

    label = atoi(read_token);
    if(label > max_label) max_label = label;
  }

  *Nbr_Materials = max_label+1;

  // Go to the beginning of the file
  rewind(file);

  // Generate the Materials structure
  Materials *material = (Materials*)malloc((max_label+1) * sizeof(Materials));
  for(label=0; label<=max_label; label++){
    strcpy(material[label].Name, "");
    material[label].A = 0;
    material[label].Density = 0;
    material[label].N_el = 0;
    material[label].X0 = 0;
    material[label].SPR = 1;
    material[label].Nuclear_data_type = None;
    material[label].NbrComponents = 0;
    material[label].Mixture_Components_label = NULL;
    material[label].Mixture_Components_fraction = NULL;
    material[label].Nuclear_Elastic = NULL;
    material[label].Elastic_Energy_List = NULL;
    material[label].Nuclear_Inelastic = NULL;
    material[label].Inelastic_Energy_List = NULL;
    material[label].PG_data = NULL;
    material[label].PG_Energy_List = NULL;
    material[label].SP_Energy = NULL;
    material[label].Stop_Pow = NULL;
    material[label].Interp_Total_Nuclear_Cross_Section = NULL;
  }
  config->Water_Material_ID = 17;

  while (fgets(read, 500, file) != NULL){

    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");	

    read_token = strtok(read, " \t\r\n");
    if(read_token == NULL || !isUnsignedInt(read_token)) continue;

    label = atoi(read_token);
    
    read_token = strtok(NULL, " \t\r\n");
    if(read_token == NULL){
      printf("\n\n ERROR: No material name for label %d \n\n", label);
      free(material);
      return NULL;
    }

    if(strlen(material[label].Name) > 0){
      printf("\n\n ERROR: Material label %d is already defined : %s \n\n", label, material[label].Name);
      free(material);
      return NULL;
    }

    strcpy(material[label].Name, read_token);
    if(strcmp(read_token, "Water") == 0 || strcmp(read_token, "water") == 0) config->Water_Material_ID = label;
  }

  fclose(file);

  return material;
}


void Free_Materials_DATA(Materials *material, unsigned int Nbr_Materials){
  if (material == NULL) return;

  int i;
  for(i=1; i<Nbr_Materials; i++){
    if(material[i].SP_Energy != NULL) free(material[i].SP_Energy);
    if(material[i].Stop_Pow != NULL) free(material[i].Stop_Pow);
    if(material[i].Mixture_Components_label != NULL) free(material[i].Mixture_Components_label);
    if(material[i].Mixture_Components_fraction != NULL) free(material[i].Mixture_Components_fraction);
    if(material[i].Interp_Total_Nuclear_Cross_Section != NULL) free(material[i].Interp_Total_Nuclear_Cross_Section);
    if(material[i].Nuclear_Elastic != NULL) free(material[i].Nuclear_Elastic);
    if(material[i].Elastic_Energy_List != NULL) free(material[i].Elastic_Energy_List);
    if(material[i].Nuclear_Inelastic != NULL) Free_DATA_Nuclear_Inelastic(material[i].Nuclear_Inelastic, material[i].Nbr_Inelastic_Energy);
    if(material[i].Inelastic_Energy_List != NULL) free(material[i].Inelastic_Energy_List);
    #if Simulate_Secondary_PromptGammas==1
      if(material[i].PG_data != NULL) Free_DATA_PG(material[i].PG_data, material[i].Nbr_PG_Energy);
      if(material[i].PG_Energy_List != NULL) free(material[i].PG_Energy_List);
    #endif
  }

  free(material);
}

