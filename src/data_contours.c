/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_contours.h"


DATA_StructList *load_all_structs(){

  int NbrStructs = 0;
  char struct_file_path[200];
  char *file_extension;

  DIR *directory;
  struct dirent *struct_file;
  
  if((directory = opendir("./structs")) != NULL){
    
    // start by counting the number of structures to load
    while((struct_file = readdir(directory)) != NULL){
      file_extension = strrchr(struct_file->d_name, '.');
      if(file_extension==NULL) continue;
      if(strcmp(file_extension, ".txt")!=0 && strcmp(file_extension, ".TXT")!=0 && strcmp(file_extension, ".mhd")!=0 && strcmp(file_extension, ".MHD")!=0) continue;
      NbrStructs++;
    }

    if(NbrStructs == 0){
      printf("\n Warning: No RT-struct found in \"./structs\"\n\n");
      closedir(directory);
      return NULL;
    }
    else{

      DATA_StructList *StructList = (DATA_StructList*) malloc(sizeof(DATA_StructList));
      StructList->Nbr_Structs = NbrStructs;
      StructList->Structs = (DATA_Struct*) malloc(NbrStructs * sizeof(DATA_Struct));

      NbrStructs = 0;
      rewinddir(directory);

      // import the structures
      while((struct_file = readdir(directory)) != NULL){
        file_extension = strrchr(struct_file->d_name, '.');

        if(file_extension==NULL) continue;
        if(strcmp(file_extension, ".txt")!=0 && strcmp(file_extension, ".TXT")!=0 && strcmp(file_extension, ".mhd")!=0 && strcmp(file_extension, ".MHD")!=0) continue;

	strcpy(StructList->Structs[NbrStructs].Name, struct_file->d_name);
        StructList->Structs[NbrStructs].Name[strlen(StructList->Structs[NbrStructs].Name)-4] = '\0';

	strcpy(struct_file_path, "./structs/");
	strcat(struct_file_path, struct_file->d_name);
	import_struct(&StructList->Structs[NbrStructs], struct_file_path);
	if(StructList->Structs[NbrStructs].Mask == NULL) continue;

	NbrStructs++;
      }
  
      if(NbrStructs < StructList->Nbr_Structs) StructList->Nbr_Structs = NbrStructs;

      closedir(directory);

      return StructList;
    }
  }
  else{
    printf("\n Warning: Unable to open directory \"./structs\"\n\n");
    return NULL;
  }

}


void import_struct(DATA_Struct *Struct, char *file_name){

  FILE *header_file;
  char read[500], *read_token;
  int isMHD = 0, isSparse = 0;

  Struct->Mask = NULL;

  header_file = fopen(file_name,"r");
  if(header_file == NULL) return;

  while(fgets(read, 500, header_file) != NULL){
    if(read[0] == '#') continue;
    strtok(read, "#");
    read_token = strtok(read, " \t=\r\n");
    if(read_token == NULL) continue;
    if(strcmp(read_token, "ElementDataFile") == 0) isMHD = 1;
    if(strcmp(read_token, "BinaryFile") == 0) isSparse = 1;
  }

  fclose(header_file);

  if(isMHD == 1) Struct->Mask = import_MHD_image(file_name, Struct->GridSize, Struct->VoxelLength, Struct->Origin);
  else if(isSparse == 1) Struct->Mask = import_Sparse_image(file_name, Struct->GridSize, Struct->VoxelLength, Struct->Origin);
  else return;

  if(Struct->Mask == NULL) return;

  int NbrVoxels = Struct->GridSize[0] * Struct->GridSize[1] * Struct->GridSize[2];

  int i;
  for(i=0; i<NbrVoxels; i++){
    if(Struct->Mask[i] != 0){
      Struct->Mask[i] = 1;
      Struct->N_Index++;
    }
  }

  Struct->IndexList = (int*)malloc(Struct->N_Index * sizeof(int));

  Struct->N_Index = 0;
  for(i=0; i<NbrVoxels; i++){
    if(Struct->Mask[i] != 0){
      Struct->IndexList[Struct->N_Index] = i;
      Struct->N_Index++;
    }
  }
}



void display_structs_information(DATA_StructList *StructList){

  
  if(StructList == NULL){ 
    printf("\nStructures not loaded\n\n");
    return;
  }
  else printf("\nLoaded structures:\n\n");

  int i;
  for(i=0; i<StructList->Nbr_Structs; i++){
    printf("Struct %d:\n", i);
    printf("Name: %s\n", StructList->Structs[i].Name);
    printf("Mask GridSize: %d %d %d\n", StructList->Structs[i].GridSize[0], StructList->Structs[i].GridSize[1], StructList->Structs[i].GridSize[2]);
    printf("Mask VoxelLength: %f %f %f (cm)\n", StructList->Structs[i].VoxelLength[0], StructList->Structs[i].VoxelLength[1], StructList->Structs[i].VoxelLength[2]);
    printf("Mask Origin: %f %f %f (cm)\n", StructList->Structs[i].Origin[0], StructList->Structs[i].Origin[1], StructList->Structs[i].Origin[2]);
    printf("Number voxels inside mask: %d\n", StructList->Structs[i].N_Index);
    printf("\n");
  }

  printf("\n");

}



void Free_all_structs(DATA_StructList **StructList){

  if(*StructList == NULL) return;

  DATA_StructList *ptr = *StructList;
  
  int i;
  for(i=0; i<ptr->Nbr_Structs; i++){
    if(ptr->Structs[i].Mask != NULL) free(ptr->Structs[i].Mask);
    if(ptr->Structs[i].IndexList != NULL) free(ptr->Structs[i].IndexList);
  }

  if(ptr->Structs != NULL) free(ptr->Structs);

  if(*StructList != NULL) free(*StructList);
  *StructList = NULL;
}




