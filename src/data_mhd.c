/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_mhd.h"

void export_MHD_image(char *file_name, int GridSize[3], VAR_DATA VoxelLength[3], VAR_SCORING *data){

  char file_path[100];
  char file_mhd_name[100];
  char file_mhd_path[100];
  char file_raw_name[100];
  char file_raw_path[100];

  char *path_ptr = strrchr(file_name, '/');
  if(path_ptr==NULL){
    strcpy(file_path, "./");
    strcpy(file_mhd_name, file_name);
  }
  else{
    strncpy(file_path, file_name, strlen(file_name)-strlen(path_ptr)+1);
    file_path[strlen(file_name)-strlen(path_ptr)+1] = '\0';
    strcpy(file_mhd_name, &path_ptr[1]);
  }


  char *file_extension = strrchr(file_mhd_name, '.');
  if(file_extension==NULL || strcmp(file_extension, ".mhd")!=0){
    strcpy(file_raw_name, file_mhd_name);
    strcat(file_raw_name, ".raw");
    strcat(file_mhd_name, ".mhd");
  }
  else{
    strncpy(file_raw_name, file_mhd_name, strlen(file_mhd_name)-4);
    file_raw_name[strlen(file_mhd_name)-4] = '\0';
    strcat(file_raw_name, ".raw");
  }

  strcpy(file_mhd_path, file_path);
  strcat(file_mhd_path, file_mhd_name);
  strcpy(file_raw_path, file_path);
  strcat(file_raw_path, file_raw_name);

  FILE *file_mhd = NULL;	// Header file (ASCII)
  FILE *file_raw = NULL;	// Data file (Binary)

  file_mhd = fopen(file_mhd_path, "w");
  fprintf(file_mhd, "ObjectType = Image\n");
  fprintf(file_mhd, "NDims = 3\n");
  fprintf(file_mhd, "DimSize = %d %d %d\n", GridSize[0], GridSize[1], GridSize[2]);
  #if VAR_DATA_PRECISION==1
    fprintf(file_mhd, "ElementSpacing = %f %f %f\n", 10*VoxelLength[0], 10*VoxelLength[1], 10*VoxelLength[2]);
  #else
    fprintf(file_mhd, "ElementSpacing = %lf %lf %lf\n", 10*VoxelLength[0], 10*VoxelLength[1], 10*VoxelLength[2]);
  #endif
  #if VAR_SCORING_PRECISION==1
    fprintf(file_mhd, "ElementType = MET_FLOAT\n");
  #else
    fprintf(file_mhd, "ElementType = MET_DOUBLE\n");
  #endif
  fprintf(file_mhd, "ElementByteOrderMSB = False\n");
  fprintf(file_mhd, "ElementDataFile = %s\n", file_raw_name);
  fclose(file_mhd);

  file_raw = fopen(file_raw_path, "wb");
  fwrite(data, sizeof(VAR_SCORING), GridSize[0]*GridSize[1]*GridSize[2], file_raw);
  fclose(file_raw);

}


int Parse_MHD_header(char *file_name, MHD_header *header){


  // Initialize Header informations
  header->NDims = 0;
  header->DimSize[0] = 0;
  header->DimSize[1] = 0;
  header->DimSize[2] = 0;
  header->ElementNumberOfChannels = 1;
  header->ElementSpacing[0] = 1.0;
  header->ElementSpacing[1] = 1.0;
  header->ElementSpacing[2] = 1.0;
  header->Offset[0] = 0.0;
  header->Offset[1] = 0.0;
  header->Offset[2] = 0.0;
  header->ElementType = NotDefined;
  header->ElementByteOrderMSB = 0;
  strcpy(header->ElementDataFile, "");


  char read[500], *read_token;

  FILE *file_mhd = fopen(file_name,"r");
  if(file_mhd == NULL){
	printf("Error: Unable to open \"%s\".\n", file_name);
	return 0;
  }


  while (fgets(read, 500, file_mhd) != NULL){

    // remove comments
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t=\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "NDims") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for NDims in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->NDims = atoi(read_token);
    }

    else if(strcmp(read_token, "ElementNumberOfChannels") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for ElementNumberOfChannels in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->ElementNumberOfChannels = atoi(read_token);
    }

    else if(strcmp(read_token, "DimSize") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for DimSize in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->DimSize[0] = atoi(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for DimSize in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->DimSize[1] = atoi(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for DimSize in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->DimSize[2] = atoi(read_token);
    }

    else if(strcmp(read_token, "ElementSpacing") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for ElementSpacing in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->ElementSpacing[0] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for ElementSpacing in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->ElementSpacing[1] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for ElementSpacing in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->ElementSpacing[2] = atof(read_token);
    }

    else if(strcmp(read_token, "Offset") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->Offset[0] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->Offset[1] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return 0;
	}
	header->Offset[2] = atof(read_token);
    }

    else if(strcmp(read_token, "ElementType") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(strcmp(read_token, "MET_FLOAT") == 0) header->ElementType = MET_FLOAT;
	else if(strcmp(read_token, "MET_DOUBLE") == 0) header->ElementType = MET_DOUBLE;
	else{
	    printf("\n Error: \"%s\" is not a valid value for ElementType in \"%s\".  Only MET_FLOAT and MET_DOUBLE formats are supported currently.\n\n", read_token, file_name);
	    return 0;
	}
    }

    else if(strcmp(read_token, "ElementByteOrderMSB") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(strcmp(read_token, "True") == 0 || strcmp(read_token, "true") == 0) header->ElementByteOrderMSB = 1;
	else if(strcmp(read_token, "False") == 0 || strcmp(read_token, "false") == 0) header->ElementByteOrderMSB = 0;
	else{
	    printf("\n Error: \"%s\" is not a valid value for ElementByteOrderMSB in \"%s\".  The value must be True or False.\n\n", read_token, file_name);
	    return 0;
	}
    }

    else if(strcmp(read_token, "ElementDataFile") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(read_token == NULL){
	    printf("\n Error: \"%s\" is not a valid value for ElementDataFile in \"%s\".\n\n", read_token, file_name);
	    return 0;
	}
	strcpy(header->ElementDataFile, read_token);
    }

  }

  fclose(file_mhd);

  return 1;

}


VAR_DATA *import_MHD_image(char *file_name, int *GridSize, VAR_DATA *VoxelLength, VAR_DATA *Origin){

  // Header informations
  MHD_header header;
  if(Parse_MHD_header(file_name, &header) != 1){
    return NULL;
  }

  // verification of saved parameters
  if(header.NDims != 3){
    printf("\n Error: NDims value must be 3 in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(header.DimSize[0] <= 0 || header.DimSize[1] <= 0 || header.DimSize[2] <= 0){
    printf("\n Error: the three DimSize values must > 0 in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(header.ElementSpacing[0] <= 0.0 || header.ElementSpacing[1] <= 0.0 || header.ElementSpacing[2] <= 0.0){
    printf("\n Error: the three ElementSpacing values must > 0.0 in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(header.ElementType == NotDefined){
    printf("\n Error: ElementType is not defined in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(strcmp(header.ElementDataFile, "") == 0){
    printf("\n Error: ElementDataFile is not defined in \"%s\".\n\n", file_name);
    return NULL;
  }


  // Read binary data
  char *path_ptr = strrchr(file_name, '/');
  char file_path[200];
  char file_raw_path[200];
  if(path_ptr!=NULL){
    strncpy(file_path, file_name, strlen(file_name)-strlen(path_ptr)+1);
    file_path[strlen(file_name)-strlen(path_ptr)+1] = '\0';
    strcpy(file_raw_path, file_path);
    strcat(file_raw_path, header.ElementDataFile);
  }
  else{
    strcpy(file_raw_path, header.ElementDataFile);
  }

  FILE *file_raw = fopen(file_raw_path,"rb");
  if(file_raw == NULL){
	printf("Error: Unable to open \"%s\".\n", file_raw_path);
	return NULL;
  }

  GridSize[0] = header.DimSize[0];
  GridSize[1] = header.DimSize[1];
  GridSize[2] = header.DimSize[2];
  
  VoxelLength[0] = (VAR_DATA)header.ElementSpacing[0]/10;
  VoxelLength[1] = (VAR_DATA)header.ElementSpacing[1]/10;
  VoxelLength[2] = (VAR_DATA)header.ElementSpacing[2]/10;
  
  Origin[0] = (VAR_DATA)header.Offset[0]/10;
  Origin[1] = (VAR_DATA)header.Offset[1]/10;
  Origin[2] = (VAR_DATA)header.Offset[2]/10;

  VAR_DATA *data = (VAR_DATA*)malloc(header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2] * sizeof(VAR_DATA));

  if(header.ElementType == MET_FLOAT && VAR_DATA_PRECISION==1){
    fread(data, sizeof(float), header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2], file_raw);
  }
  else if(header.ElementType == MET_FLOAT && VAR_DATA_PRECISION==2){
    float *buffer = (float*)malloc(header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2] * sizeof(float));
    fread(buffer, sizeof(float), header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2], file_raw);
    int i;
    for(i=0; i<header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2]; i++){
	data[i] = (VAR_DATA)buffer[i];
    }
    free(buffer);
  }
  else if(header.ElementType == MET_DOUBLE && VAR_DATA_PRECISION==1){
    double *buffer = (double*)malloc(header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2] * sizeof(double));
    fread(buffer, sizeof(double), header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2], file_raw);
    int i;
    for(i=0; i<header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2]; i++){
	data[i] = (VAR_DATA)buffer[i];
    }
    free(buffer);
  }
  else if(header.ElementType == MET_DOUBLE && VAR_DATA_PRECISION==2){
    fread(data, sizeof(double), header.ElementNumberOfChannels*GridSize[0]*GridSize[1]*GridSize[2], file_raw);
  }
  else{
	printf("Error: Unable to read data in \"%s\".\n", file_raw_path);
	if(data != NULL) free(data);
	return NULL;
  }

  fclose(file_raw);

  return data;
}


