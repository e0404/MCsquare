/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_sparse.h"

void export_Sparse_image(char *file_name, DATA_config *config, DATA_CT *ct, plan_parameters *plan, VAR_SCORING *data, VAR_SCORING threshold){

  char file_path[200];
  char file_header_name[200];
  char file_header_path[200];
  char file_bin_name[200];
  char file_bin_path[200];

  char *path_ptr = strrchr(file_name, '/');
  if(path_ptr==NULL){
    strcpy(file_path, "./");
    strcpy(file_header_name, file_name);
  }
  else{
    strncpy(file_path, file_name, strlen(file_name)-strlen(path_ptr)+1);
    file_path[strlen(file_name)-strlen(path_ptr)+1] = '\0';
    strcpy(file_header_name, &path_ptr[1]);
  }


  char *file_extension = strrchr(file_header_name, '.');
  if(file_extension==NULL || strcmp(file_extension, ".txt")!=0){
    strcpy(file_bin_name, file_header_name);
    strcat(file_bin_name, ".bin");
    strcat(file_header_name, ".txt");
  }
  else{
    strncpy(file_bin_name, file_header_name, strlen(file_header_name)-4);
    file_bin_name[strlen(file_header_name)-4] = '\0';
    strcat(file_bin_name, ".bin");
  }

  strcpy(file_header_path, file_path);
  strcat(file_header_path, file_header_name);
  strcpy(file_bin_path, file_path);
  strcat(file_bin_path, file_bin_name);

  FILE *file_header = NULL;	// Header file (ASCII)
  FILE *file_bin = NULL;	// Data file (Binary)

  int FileAlreadyCreated = 0;
  DATA_Sparse_Header Header;

  if(File_exists(file_header_path) == 1){
    Header = Read_Sparse_Header(file_header_path);
    if(Header.SimulationTimestamp == config->timestamp && strcmp(Header.PlanName, plan->PlanName) == 0) FileAlreadyCreated = 1;
  }
    
  if(FileAlreadyCreated == 0){
    file_header = fopen(file_header_path, "w");
    fprintf(file_header, "# MCsquare sparse matrix format\n");
    struct tm *date = localtime(&config->timestamp);
    fprintf(file_header, "SimulationDate = %d/%d/%d %d:%d:%d\n", date->tm_year+1900, date->tm_mon+1, date->tm_mday, date->tm_hour, date->tm_min, date->tm_sec);
    fprintf(file_header, "PlanName = %s\n", plan->PlanName);
    if(config->Beamlet_Mode == 1){ 
	fprintf(file_header, "SimulationMode = Beamlet\n");
	fprintf(file_header, "NbrSpots = %u\n", config->TotalNbrSpots);
    }
    if(config->Robustness_Mode == 1){ 
	fprintf(file_header, "SimulationMode = Robustness\n");
	fprintf(file_header, "RobustParam_SystematicSetupError = %f %f %f\n", config->Systematic_Setup_Error[0], config->Systematic_Setup_Error[1], config->Systematic_Setup_Error[2]);
	fprintf(file_header, "RobustParam_RandomSetupError = %f %f %f\n", config->Random_Setup_Error[0], config->Random_Setup_Error[1], config->Random_Setup_Error[2]);
	fprintf(file_header, "RobustParam_SystematicRangeError = %f\n", config->Systematic_Range_Error);
	fprintf(file_header, "Scenario_SystematicSetupError = %f %f %f\n", config->Current_Systematic_setup[0], config->Current_Systematic_setup[1], config->Current_Systematic_setup[2]);
	fprintf(file_header, "Scenario_RandomSetupError = %f %f %f\n", config->Current_Random_setup[0], config->Current_Random_setup[1], config->Current_Random_setup[2]);
	fprintf(file_header, "Scenario_SystematicRangeError = %f\n", config->Current_Range_error);
    }
    if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 0){
	fprintf(file_header, "SimulationMode = 4D\n");
	fprintf(file_header, "Dose_Accumulation = disabled\n");
	fprintf(file_header, "4D_Phase = %d\n", config->Current_4D_phase);
    }
    else if(config->Simu_4D_Mode == 1 && config->Dose_4D_Accumulation == 1){
	fprintf(file_header, "SimulationMode = 4D\n");
	fprintf(file_header, "Dose_Accumulation = enabled\n");
    }
    fprintf(file_header, "ImageSize = %d %d %d\n", ct->GridSize[0], ct->GridSize[1], ct->GridSize[2]);
    #if VAR_DATA_PRECISION==1
      fprintf(file_header, "VoxelSpacing = %f %f %f\n", 10*ct->VoxelLength[0], 10*ct->VoxelLength[1], 10*ct->VoxelLength[2]);
    #else
      fprintf(file_header, "VoxelSpacing = %lf %lf %lf\n", 10*ct->VoxelLength[0], 10*ct->VoxelLength[1], 10*ct->VoxelLength[2]);
    #endif
    fprintf(file_header, "BinaryFile = %s\n", file_bin_name);
    fclose(file_header);

    file_bin = fopen(file_bin_path, "wb");
  }
  else{
    strcpy(file_bin_path, file_path);
    strcat(file_bin_path, Header.BinaryFile);
    file_bin = fopen(file_bin_path, "ab");
  }

  uint32_t i,j,k, index=0, NbrSelectedVoxel=0;
  for(i=1; i<=ct->GridSize[2]; i++){
    for(j=1; j<=ct->GridSize[1]; j++){
      for(k=1; k<=ct->GridSize[0]; k++){
        if(data[index] > threshold) NbrSelectedVoxel++;
        index++;
      }
    }
  }


/*
// Standard sparse format
  if(config->Robustness_Mode == 0 && config->Beamlet_Mode == 0) printf("\nSparse matrix (threshold=%f): %u / %u voxels stored\n", threshold, NbrSelectedVoxel, index);

  fwrite(&NbrSelectedVoxel, sizeof(uint32_t), 1, file_bin); // Number of nonzero voxels

  if(config->Beamlet_Mode == 1){
    uint32_t FieldID = (uint32_t)plan->fields[0].FieldID;
    fwrite(&FieldID, sizeof(uint32_t), 1, file_bin); // FieldID
    uint32_t ControlPointIndex = (uint32_t)plan->fields[0].ControlPoints[0].ControlPointIndex;
    fwrite(&ControlPointIndex, sizeof(uint32_t), 1, file_bin); // ControlPointID
    float Spot_X = (float)plan->fields[0].ControlPoints[0].spots[0].Spot_X;
    fwrite(&Spot_X, sizeof(float), 1, file_bin); // Spot_X
    float Spot_Y = (float)plan->fields[0].ControlPoints[0].spots[0].Spot_Y;
    fwrite(&Spot_Y, sizeof(float), 1, file_bin); // Spot_Y
  }

  index = 0;
  float value;

  for(i=1; i<=ct->GridSize[2]; i++){
    for(j=1; j<=ct->GridSize[1]; j++){
      for(k=1; k<=ct->GridSize[0]; k++){
        if(data[index] > threshold){
	  fwrite(&index, sizeof(uint32_t), 1, file_bin); // 1D-index
	  value = (float)data[index];
	  fwrite(&value, sizeof(float), 1, file_bin); 	 // voxel value
	}
        index++;
      }
    }
  }
  
  fclose(file_bin);
*/



// Continuous sparse format:
  fwrite(&NbrSelectedVoxel, sizeof(uint32_t), 1, file_bin); // Number of nonzero voxels

  if(config->Beamlet_Mode == 1){
    uint32_t FieldID = (uint32_t)plan->fields[0].FieldID;
    fwrite(&FieldID, sizeof(uint32_t), 1, file_bin); // FieldID
    uint32_t ControlPointIndex = (uint32_t)plan->fields[0].ControlPoints[0].ControlPointIndex;
    fwrite(&ControlPointIndex, sizeof(uint32_t), 1, file_bin); // ControlPointID
    float Spot_X = (float)plan->fields[0].ControlPoints[0].spots[0].Spot_X;
    fwrite(&Spot_X, sizeof(float), 1, file_bin); // Spot_X
    float Spot_Y = (float)plan->fields[0].ControlPoints[0].spots[0].Spot_Y;
    fwrite(&Spot_Y, sizeof(float), 1, file_bin); // Spot_Y
  }

  index = 0;
  float *value_array = (float*)malloc(ct->GridSize[0]*ct->GridSize[1]*ct->GridSize[2] * sizeof(float));
  uint32_t FirstIndex, NumContinuousValues = 0;

  for(i=1; i<=ct->GridSize[2]; i++){
    for(j=1; j<=ct->GridSize[1]; j++){
      for(k=1; k<=ct->GridSize[0]; k++){
        if(data[index] > threshold){
	  value_array[NumContinuousValues] = (float)data[index];
	  if(NumContinuousValues == 0) FirstIndex = index;
	  NumContinuousValues += 1;
	}
	else if(NumContinuousValues != 0){
	  fwrite(&NumContinuousValues, sizeof(uint32_t), 1, file_bin); // Number of continuous non-zero voxels
	  fwrite(&FirstIndex, sizeof(uint32_t), 1, file_bin); // First 1D-index of the list of continuous voxels
	  fwrite(value_array, sizeof(float)*NumContinuousValues, 1, file_bin); 	 // voxel values
	  NumContinuousValues = 0;
	}
        index++;
      }
    }
  }

  if(NumContinuousValues != 0){
	  fwrite(&NumContinuousValues, sizeof(uint32_t), 1, file_bin); // Number of continuous non-zero voxels
	  fwrite(&FirstIndex, sizeof(uint32_t), 1, file_bin); // First 1D-index of the list of continuous voxels
	  fwrite(value_array, sizeof(float)*NumContinuousValues, 1, file_bin); 	 // voxel values
  }

  free(value_array);

  fclose(file_bin);

  return;  
}


VAR_DATA *import_Sparse_image(char *file_name, int *GridSize, VAR_DATA *VoxelLength, VAR_DATA *Origin){


  // Read Header data
  DATA_Sparse_Header Header = Read_Sparse_Header(file_name);

  if(Header.SimulationTimestamp == 0) return NULL;

  if(Header.Beamlet_mode != 0){
    printf("Error: Unable to import Sparse file \"%s\" in Beamlet mode.\n", file_name);
    return NULL;
  }
  if(Header.ImageSize[0] <= 0 || Header.ImageSize[1] <= 0 || Header.ImageSize[2] <= 0){
    printf("\n Error: the three ImageSize values must > 0 in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(Header.VoxelSpacing[0] <= 0.0 || Header.VoxelSpacing[1] <= 0.0 || Header.VoxelSpacing[2] <= 0.0){
    printf("\n Error: the three VoxelSpacing values must > 0.0 in \"%s\".\n\n", file_name);
    return NULL;
  }
  if(strcmp(Header.BinaryFile, "") == 0){
    printf("\n Error: BinaryFile is not defined in \"%s\".\n\n", file_name);
    return NULL;
  }

  
  // Read binary data
  char file_path_bin[200];
  char file_path[200];

  char *path_ptr = strrchr(file_name, '/');
  if(path_ptr==NULL) strcpy(file_path, "./");
  else{
    strncpy(file_path, file_name, strlen(file_name)-strlen(path_ptr)+1);
    file_path[strlen(file_name)-strlen(path_ptr)+1] = '\0';
  }
  
  strcpy(file_path_bin, file_path);
  strcat(file_path_bin, Header.BinaryFile);

  FILE *fid = fopen(file_path_bin,"rb");
  if(fid == NULL){
	printf("Error: Unable to open \"%s\".\n", file_path_bin);
	return NULL;
  }

  uint32_t NonZeroVoxels, NbrContinuousValues, ReadVoxels, FirstIndex;
  int NbrVoxels = Header.ImageSize[0]*Header.ImageSize[1]*Header.ImageSize[2];

  fread(&NonZeroVoxels, sizeof(uint32_t), 1, fid);
  if(NonZeroVoxels > NbrVoxels){
	printf("Error: Number of non-zero values must be lower than number of voxels in \"%s\".\n", file_path_bin);
	fclose(fid);
	return NULL;
  }

  VAR_DATA *data = (VAR_DATA*)calloc(NbrVoxels, sizeof(VAR_DATA));

  ReadVoxels = 0;

  while(1){
    fread(&NbrContinuousValues, sizeof(uint32_t), 1, fid);
    ReadVoxels += NbrContinuousValues;

    fread(&FirstIndex, sizeof(uint32_t), 1, fid);
    if((FirstIndex+NbrContinuousValues) >= NbrVoxels){
	printf("Error: Exceeds the maximum index in \"%s\".\n", file_path_bin);
	fclose(fid);
	free(data);
	return NULL;
    }

    fread(&data[FirstIndex], sizeof(float), NbrContinuousValues, fid);

    if(ReadVoxels >= NonZeroVoxels) break;
  }

  fclose(fid);

  GridSize[0] = Header.ImageSize[0];
  GridSize[1] = Header.ImageSize[1];
  GridSize[2] = Header.ImageSize[2];
  
  VoxelLength[0] = (VAR_DATA)Header.VoxelSpacing[0]/10;
  VoxelLength[1] = (VAR_DATA)Header.VoxelSpacing[1]/10;
  VoxelLength[2] = (VAR_DATA)Header.VoxelSpacing[2]/10;
  
  Origin[0] = (VAR_DATA)Header.Offset[0]/10;
  Origin[1] = (VAR_DATA)Header.Offset[1]/10;
  Origin[2] = (VAR_DATA)Header.Offset[2]/10;

  return data;
}


DATA_Sparse_Header Read_Sparse_Header(char *file_name){

  DATA_Sparse_Header Header = Init_Sparse_Header();

  char read[500], *read_token;

  FILE *file = fopen(file_name,"r");
  if(file == NULL){
	printf("Error: Unable to open \"%s\".\n", file_name);
	return Init_Sparse_Header();
  }


  while (fgets(read, 500, file) != NULL){

    // on ignore les commentaires
    if(read[0] == '#') continue;
    strtok(read, "#");

    read_token = strtok(read, " \t=\r\n");
    if(read_token == NULL) continue;

    if(strcmp(read_token, "SimulationDate") == 0){
	time_t timestamp = time(NULL);
	struct tm *tmp = localtime(&timestamp);
	struct tm date;
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_year = atoi(read_token) - 1900;
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_mon = atoi(read_token) - 1;
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_mday = atoi(read_token);
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_hour = atoi(read_token);
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_min = atoi(read_token);
	read_token = strtok(NULL, " \t=\r\n/:");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SimulationDate in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
	date.tm_sec = atoi(read_token);
	date.tm_isdst = tmp->tm_isdst;
	Header.SimulationTimestamp = mktime(&date);
    }

    else if(strcmp(read_token, "PlanName") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	strcpy(Header.PlanName, read_token);
    }

    else if(strcmp(read_token, "ImageSize") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for ImageSize in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
        Header.ImageSize[0] = atoi(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for ImageSize in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.ImageSize[1] = atoi(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for ImageSize in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.ImageSize[2] = atoi(read_token);
    }

    else if(strcmp(read_token, "VoxelSpacing") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for VoxelSpacing in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.VoxelSpacing[0] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for VoxelSpacing in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.VoxelSpacing[1] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for VoxelSpacing in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.VoxelSpacing[2] = atof(read_token);
    }

    else if(strcmp(read_token, "Offset") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.Offset[0] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.Offset[1] = atof(read_token);
	read_token = strtok(NULL, " \t\r\n");
	if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Offset in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
	}
	Header.Offset[2] = atof(read_token);
    }

    else if(strcmp(read_token, "BinaryFile") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	strcpy(Header.BinaryFile, read_token);
    }

    else if(strcmp(read_token, "SimulationMode") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(strcmp(read_token, "Beamlet") == 0) Header.Beamlet_mode = 1;
    }

    else if(strcmp(read_token, "NbrSpots") == 0){
	read_token = strtok(NULL, " \t=\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for NbrSpots in \"%s\"\n\n", read_token, file_name);
	    return Init_Sparse_Header();
        }
        Header.NbrSpots = atoi(read_token);
    }

  }

  fclose(file);

  //Display_Sparse_Header(&Header);

  return Header;
}


void Display_Sparse_Header(DATA_Sparse_Header *Header){

  printf("\n");
  printf("SimulationTimestamp = %d\n", (int)Header->SimulationTimestamp);	
  printf("PlanName = %s\n", Header->PlanName);
  printf("Image Size = %d %d %d\n", Header->ImageSize[0], Header->ImageSize[1], Header->ImageSize[2]);
  printf("Voxel Size = %f %f %f\n", Header->VoxelSpacing[0], Header->VoxelSpacing[1], Header->VoxelSpacing[2]);	
  printf("BinaryFile = %s\n", Header->BinaryFile);
  printf("\n");

  return;

}


DATA_Sparse_Header Init_Sparse_Header(){

  DATA_Sparse_Header Header;
  Header.SimulationTimestamp = 0;
  strcpy(Header.PlanName, "");	
  Header.ImageSize[0] = 0;	
  Header.ImageSize[1] = 0;	
  Header.ImageSize[2] = 0;
  Header.VoxelSpacing[0] = 1.0;
  Header.VoxelSpacing[1] = 1.0;
  Header.VoxelSpacing[2] = 1.0;
  strcpy(Header.BinaryFile, "");
  Header.Beamlet_mode = 0;
  Header.NbrSpots = 0;

  return Header;
}


int Merge_Sparse_Files(char *InputPath, char *FileName, int NbrDirectories, char *OutputFile){

  char file_header_path[200];
  char file_bin_path[200];
  char out_header_path[200];
  char out_bin_path[200];

  char *file_extension = strrchr(FileName, '.');
  if(file_extension==NULL || strcmp(file_extension, ".txt")!=0){
    strcpy(file_bin_path, FileName);
    strcat(file_bin_path, ".bin");
    strcpy(file_header_path, FileName);
    strcat(file_header_path, ".txt");
  }
  else{
    strncpy(file_bin_path, FileName, strlen(FileName)-4);
    file_bin_path[strlen(FileName)-4] = '\0';
    strcpy(file_header_path, file_bin_path);
    strcat(file_bin_path, ".bin");
    strcat(file_header_path, ".txt");
  }

  file_extension = strrchr(OutputFile, '.');
  if(file_extension==NULL || strcmp(file_extension, ".txt")!=0){
    strcpy(out_bin_path, OutputFile);
    strcat(out_bin_path, ".bin");
    strcpy(out_header_path, OutputFile);
    strcat(out_header_path, ".txt");
  }
  else{
    strncpy(out_bin_path, OutputFile, strlen(OutputFile)-4);
    out_bin_path[strlen(OutputFile)-4] = '\0';
    strcpy(out_header_path, out_bin_path);
    strcat(out_bin_path, ".bin");
    strcat(out_header_path, ".txt");
  }
  
  char from[200], ID[10];
  sprintf(from, "%s1/%s", InputPath, file_header_path);
  CopyFile(out_header_path, from);

  FILE *fd_to, *fd_from;
  char buf[4096];
  ssize_t nread;

  fd_to = fopen(out_bin_path, "wb");
  if(fd_to == NULL) return -1;

  int i;
  for(i=0; i<NbrDirectories; i++){
    sprintf(from, "%s%d/%s", InputPath, i+1, file_bin_path);

    fd_from = fopen(from, "rb");
    if(fd_from == NULL){
      printf("\nError: unable to open temporary file %s \n", from);
      fclose(fd_to);
      return -1;
    }

    while (nread = fread(buf, 1, 4096, fd_from), nread > 0){
      char *out_ptr = buf;
      ssize_t nwritten;

      do{
        nwritten = fwrite(out_ptr, 1, nread, fd_to);
        if (nwritten >= 0){
          nread -= nwritten;
          out_ptr += nwritten;
        }
        else{
	  printf("\nError: unable to write merged sparse file %s \n", out_bin_path);
	  fclose(fd_to);
	  fclose(fd_from);
	  return -1;
        }
      } while(nread > 0);
    }


    fclose(fd_from);

    // Remove Bin and header files
    remove(from);
    strcpy(from, file_header_path);
    str_replace("{Num}" , ID , from);
    remove(from);

  }

  fclose(fd_to);
  return 0;
}



