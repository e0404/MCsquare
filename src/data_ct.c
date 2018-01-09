/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_ct.h"

void Init_CT_DATA(DATA_CT *ct, int Nx, int Ny, int Nz, VAR_DATA Lx, VAR_DATA Ly, VAR_DATA Lz){

  ct->GridSize[0] = Nx;
  ct->GridSize[1] = Ny;
  ct->GridSize[2] = Nz;
  ct->Nbr_voxels = Nx*Ny*Nz;
  ct->Length[0] = Lx;
  ct->Length[1] = Ly;
  ct->Length[2] = Lz;
  ct->VoxelLength[0] = Lx/Nx;
  ct->VoxelLength[1] = Ly/Ny;
  ct->VoxelLength[2] = Lz/Nz;

  ct->material = (unsigned short int*)malloc(Nx*Ny*Nz * sizeof(unsigned short int));
  ct->density = (VAR_DATA*)malloc(Nx*Ny*Nz * sizeof(VAR_DATA));

  ct->Conversion_HU_Density = NULL;
  ct->Conversion_Densities = NULL;
  ct->Conversion_HU_Material = NULL;
  ct->Conversion_Density_Material = NULL;
  ct->Conversion_Material_labels = NULL;
}


int Read_Density_conversion_data(char *file_name, DATA_CT *ct){

  FILE *file;
  char tmp[256];
  char *data;
  int Num_data=0, i=0;
  VAR_DATA previous = -9999;

  file = fopen(file_name,"r");
  if(file == NULL){
	printf("Error: Unable to open \"%s\". File required for the HU to density conversion.\n", file_name);
	return 1;
  }

  // Firstly, we pass through the file to count the number of data points
  while (fgets(tmp, 256, file) != NULL){

    // on ignore les commentaires
    if(tmp[0] == '#') continue;
    strtok(tmp, "#");	

    data = strtok(tmp, " \t\r\n");
    if(data == NULL) continue;

    Num_data++;
  }

  ct->Num_Density_Data = Num_data;
  ct->Conversion_HU_Density = (VAR_DATA*)malloc(Num_data * sizeof(VAR_DATA));
  ct->Conversion_Densities = (VAR_DATA*)malloc(Num_data * sizeof(VAR_DATA));

  // Go to the beginning of the file
  rewind(file);

  // Secondly, data are stored in arrays
  while (fgets(tmp, 256, file) != NULL){

    // on ignore les commentaires
    if(tmp[0] == '#') continue;
    strtok(tmp, "#");	

    data = strtok(tmp, " \t\r\n");
    if(data == NULL) continue;

    ct->Conversion_HU_Density[i] = atof(data);
    if(previous >= ct->Conversion_HU_Density[i]){
	printf("Warning: HU to density conversion data are not sorted in ascending order in %s.  This may lead to conversion errors\n", file_name);
    }
    previous = ct->Conversion_HU_Density[i];

    data = strtok(NULL, " \t\r\n");
    ct->Conversion_Densities[i] = atof(data);

    i++;
  }

  fclose(file);

  return 0;
}


void Display_Density_conversion_data(DATA_CT *ct){

  printf("\n\nHU to density conversion data:\n");
  printf("HU\t\tDensity (g/cm3)\n");

  int i;
  for(i=0; i<ct->Num_Density_Data; i++){
    printf("%.1f\t\t%.4f\n", ct->Conversion_HU_Density[i], ct->Conversion_Densities[i]);
  }

  printf("\n\n");

}


VAR_DATA HU_to_Density_convertion(VAR_DATA HU, DATA_CT *ct){

  int index = Binary_Search(HU, ct->Conversion_HU_Density, ct->Num_Density_Data);
  if(index < 0) index = 0;
  else if(index > ct->Num_Density_Data - 1) index = ct->Num_Density_Data - 1;

  VAR_DATA density = (VAR_DATA)Linear_Interpolation(HU, ct->Conversion_HU_Density[index], ct->Conversion_HU_Density[index+1], ct->Conversion_Densities[index], ct->Conversion_Densities[index+1]);

  if(density <= 0) density = 1e-6;

  return density;

}


int Read_Material_conversion_data(char *file_name, DATA_CT *ct){

  FILE *file;
  char tmp[256];
  char *data;
  int Num_data=0, i=0;
  VAR_DATA previous = -9999;

  if(ct->Conversion_HU_Density == NULL){
	printf("Error: HU to density conversion data are required before reading the Material conversion data in %s.\n", file_name);
	return 1;
  }

  file = fopen(file_name,"r");
  if(file == NULL){
	printf("Error: Unable to open \"%s\". File required for the HU to material conversion.\n", file_name);
	return 1;
  }

  // Firstly, we pass through the file to count the number of data points
  while (fgets(tmp, 256, file) != NULL){

    // on ignore les commentaires
    if(tmp[0] == '#') continue;
    strtok(tmp, "#");	

    data = strtok(tmp, " \t\r\n");
    if(data == NULL) continue;

    Num_data++;
  }

  ct->Num_Materials_Data = Num_data;
  ct->Conversion_HU_Material = (VAR_DATA*)malloc(Num_data * sizeof(VAR_DATA));
  ct->Conversion_Density_Material = (VAR_DATA*)malloc(Num_data * sizeof(VAR_DATA));
  ct->Conversion_Material_labels = (unsigned short int*)malloc(Num_data * sizeof(unsigned short int));

  // Go to the beginning of the file
  rewind(file);

  // Secondly, data are stored in arrays
  while (fgets(tmp, 256, file) != NULL){

    // on ignore les commentaires
    if(tmp[0] == '#') continue;
    strtok(tmp, "#");	

    data = strtok(tmp, " \t\r\n");
    if(data == NULL) continue;

    ct->Conversion_HU_Material[i] = atof(data);
    if(previous >= ct->Conversion_HU_Material[i]){
	printf("Warning: HU to material conversion data are not sorted in ascending order in %s.  This may lead to conversion errors\n", file_name);
    }
    previous = ct->Conversion_HU_Material[i];
    ct->Conversion_Density_Material[i] = HU_to_Density_convertion(ct->Conversion_HU_Material[i], ct);

    data = strtok(NULL, " \t\r\n");
    ct->Conversion_Material_labels[i] = (unsigned short int)atoi(data);

    i++;
  }

  fclose(file);

  return 0;
}



void Display_Material_conversion_data(DATA_CT *ct){

  printf("\n\nHU to material conversion data:\n");
  printf("HU\t\tDensity (g/cm3)\tMaterial label\n");

  int i;
  for(i=0; i<ct->Num_Materials_Data; i++){
    printf("%.1f\t\t%.4f\t\t%d\n", ct->Conversion_HU_Material[i], ct->Conversion_Density_Material[i], ct->Conversion_Material_labels[i]);
  }

  printf("\n\n");

}


unsigned short int HU_to_Material_convertion(VAR_DATA HU, DATA_CT *ct){

  int index = Binary_Search(HU, ct->Conversion_HU_Material, ct->Num_Materials_Data);
  if(index < 0) index = 0;
  else if(index > ct->Num_Materials_Data - 1) index = ct->Num_Materials_Data - 1;

  return ct->Conversion_Material_labels[index];

}


unsigned short int Density_to_Material_convertion(VAR_DATA Density, DATA_CT *ct){

  int index = Binary_Search(Density, ct->Conversion_Density_Material, ct->Num_Materials_Data);
  if(index < 0) index = 0;
  else if(index > ct->Num_Materials_Data - 1) index = ct->Num_Materials_Data - 1;

  return ct->Conversion_Material_labels[index];

}


void Free_CT_DATA(DATA_CT *ct){
  if(ct == NULL) return;

  if(ct->material != NULL) free(ct->material);
  if(ct->density != NULL) free(ct->density);

  if(ct->Conversion_HU_Density != NULL) free(ct->Conversion_HU_Density);
  if(ct->Conversion_Densities != NULL) free(ct->Conversion_Densities);
  if(ct->Conversion_HU_Material != NULL) free(ct->Conversion_HU_Material);
  if(ct->Conversion_Density_Material != NULL) free(ct->Conversion_Density_Material);
  if(ct->Conversion_Material_labels != NULL) free(ct->Conversion_Material_labels);

  free(ct);
}
