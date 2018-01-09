/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_DVH.h"


void compute_all_DVH(DATA_config *config, VAR_SCORING *Dose, VAR_SCORING DoseScaling){

  char DVH_file_path[200];
  char struct_file_path[200];
  char *file_extension;

  int ListSize, NbrStructs = 0;
  int *Struct;

  DIR *directory;
  struct dirent *struct_file;
  
  if((directory = opendir("./structs")) != NULL){
    
    while((struct_file = readdir(directory)) != NULL){

      file_extension = strrchr(struct_file->d_name, '.');
      if(file_extension==NULL) continue;
      if(strcmp(file_extension, ".txt")!=0 && strcmp(file_extension, ".TXT")!=0 && strcmp(file_extension, ".mhd")!=0 && strcmp(file_extension, ".MHD")!=0) continue;

      //printf ("%s\n", struct_file->d_name);
      strcpy(struct_file_path, "./structs/");
      strcat(struct_file_path, struct_file->d_name);

      Struct = import_mask(struct_file_path, &ListSize);
      if(Struct == NULL) continue;

      NbrStructs++;

      strcpy(DVH_file_path, config->Output_Directory);
      strcat(DVH_file_path, "DVH_");
      strcat(DVH_file_path, struct_file->d_name);
      DVH_file_path[strlen(DVH_file_path)-4] = '\0';
      strcat(DVH_file_path, config->output_robustness_suffix);
      strcat(DVH_file_path, config->output_beamlet_suffix);
      strcat(DVH_file_path, config->output_4D_suffix);
      strcat(DVH_file_path, config->output_beams_suffix);
      strcat(DVH_file_path, ".txt");

      //printf ("%s\n", DVH_file_path);

      compute_DVH(Struct, ListSize, Dose, DoseScaling, DVH_file_path);
      free(Struct);
    }

    closedir(directory);

    if(NbrStructs == 0) printf("\n Warning: No RT-struct found in \"./structs\"\n\n");
  }
  else printf("\n Warning: Unable to open directory \"./structs\"\n\n");

  return;
}


int *import_mask(char *file_name, int *ListSize){

  FILE *header_file;
  char read[500], *read_token;
  int isMHD = 0, isSparse = 0;

  header_file = fopen(file_name,"r");
  if(header_file == NULL) return NULL;

  while(fgets(read, 500, header_file) != NULL){
    if(read[0] == '#') continue;
    strtok(read, "#");
    read_token = strtok(read, " \t=\r\n");
    if(read_token == NULL) continue;
    if(strcmp(read_token, "ElementDataFile") == 0) isMHD = 1;
    if(strcmp(read_token, "BinaryFile") == 0) isSparse = 1;
  }

  fclose(header_file);

  int GridSize[3];
  VAR_DATA VoxelLength[3], Origin[3];
  VAR_DATA *mask = NULL;

  if(isMHD == 1) mask = import_MHD_image(file_name, GridSize, VoxelLength, Origin);
  else if(isSparse == 1) mask = import_Sparse_image(file_name, GridSize, VoxelLength, Origin);
  else return NULL;

  if(mask == NULL) return NULL;

  int NbrVoxels = GridSize[0]*GridSize[1]*GridSize[2];

  int i, ListCount = 0;
  for(i=0; i<NbrVoxels; i++){
    if(mask[i] != 0) ListCount++;
  }

  *ListSize = ListCount;
  int *ListIndex = (int*)malloc(ListCount * sizeof(int));

  ListCount = 0;
  for(i=0; i<NbrVoxels; i++){
    if(mask[i] != 0){
      ListIndex[ListCount] = i;
      ListCount++;
    }
  }

  if(mask != NULL) free(mask);

  return ListIndex;
}


void compute_DVH(int *ListIndex, int ListSize, VAR_SCORING *Dose, VAR_SCORING DoseScaling, char *Output_file){

  int NbrBins = 4096;
  float interval[2] = {0, 100};

  float bin_size = (interval[1]-interval[0])/NbrBins;
  float *DVH = (float*)calloc(NbrBins, sizeof(float));

  int i, bin;
  for(i=0; i<ListSize; i++){
    bin = floor(Dose[ListIndex[i]] * DoseScaling / bin_size);
    if(bin >= NbrBins) bin = NbrBins-1;
    DVH[bin] += 1;
  }

  for(i=NbrBins-2; i>=0; i--){
    DVH[i] += DVH[i+1];
  }

  for(i=NbrBins-1; i>=0; i--){
    DVH[i] = 100*DVH[i]/DVH[0];
  }


  FILE *file = NULL;
  file = fopen(Output_file, "w");
  for(i=0; i<NbrBins; i++){
    fprintf(file, "%f\t%lf\n", (i+1)*bin_size, DVH[i]);
  }
  fclose(file);

  if(DVH != NULL) free(DVH);

  return;
}
