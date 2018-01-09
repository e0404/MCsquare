/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_4D.h"

DATA_CT **Import_4DCT(DATA_config *config){

  char phase_file_path[200];
  config->Num_4DCT_phases = 0;
  int i;
  DATA_CT **CT =  NULL;

  // count number of phases.
  sprintf(phase_file_path, "./4DCT/CT_1.mhd");
  FILE *file_mhd = fopen(phase_file_path,"r");
  while(file_mhd != NULL){
    config->Num_4DCT_phases += 1;
    fclose(file_mhd);
    sprintf(phase_file_path, "./4DCT/CT_%d.mhd", config->Num_4DCT_phases+1);
    file_mhd = fopen(phase_file_path,"r");
  }

  if(config->Num_4DCT_phases == 0){
    printf("\n Error: 4D-CT phases not found in directory \"./4DCT\"\n\n");
    return NULL;
  }

  printf("\n4D data: %d phases found in directory \"./4DCT\"\n", config->Num_4DCT_phases);

  CT = (DATA_CT**)malloc(config->Num_4DCT_phases * sizeof(DATA_CT*));
  for(i=0; i < config->Num_4DCT_phases; i++) CT[i] = NULL;



  // import HU conversion data
  CT[0] = (DATA_CT*) malloc(sizeof(DATA_CT));
  CT[0]->Conversion_HU_Density = NULL;
  CT[0]->Conversion_Densities = NULL;
  CT[0]->Conversion_HU_Material = NULL;
  CT[0]->Conversion_Density_Material = NULL;
  CT[0]->Conversion_Material_labels = NULL;

  if(Read_Density_conversion_data(config->HU_Density_File, CT[0]) != 0){
    Free_4DCT(CT, config->Num_4DCT_phases);
    return NULL;
  }
//  Display_Density_conversion_data(CT);

  if(Read_Material_conversion_data(config->HU_Material_File, CT[0]) != 0){
    Free_4DCT(CT, config->Num_4DCT_phases);
    return NULL;
  }
//  Display_Material_conversion_data(CT);




  // import all phases
  int GridSize[3];
  VAR_DATA VoxelLength[3], Origin[3];
  int phaseID = 0;
  VAR_DATA *hu;


  for(phaseID=0; phaseID<config->Num_4DCT_phases; phaseID++){

    printf(" Loading phase %d\n", phaseID+1);

    if(phaseID != 0){
      CT[phaseID] = (DATA_CT*) malloc(sizeof(DATA_CT));
      CT[phaseID]->Conversion_HU_Density = CT[0]->Conversion_HU_Density;
      CT[phaseID]->Conversion_Densities = CT[0]->Conversion_Densities;
      CT[phaseID]->Conversion_HU_Material = CT[0]->Conversion_HU_Material;
      CT[phaseID]->Conversion_Density_Material = CT[0]->Conversion_Density_Material;
      CT[phaseID]->Conversion_Material_labels = CT[0]->Conversion_Material_labels;
      CT[phaseID]->Num_Density_Data = CT[0]->Num_Density_Data;
      CT[phaseID]->Num_Materials_Data = CT[0]->Num_Materials_Data;	
    }

    sprintf(phase_file_path, "./4DCT/CT_%d.mhd", phaseID+1);
    hu = import_MHD_image(phase_file_path, GridSize, VoxelLength, Origin);
    if(hu == NULL){
      Free_4DCT(CT, config->Num_4DCT_phases);
      return NULL;
    }

    CT[phaseID]->GridSize[0] = GridSize[0];
    CT[phaseID]->GridSize[1] = GridSize[1];
    CT[phaseID]->GridSize[2] = GridSize[2];

    CT[phaseID]->Nbr_voxels = GridSize[0]*GridSize[1]*GridSize[2];
  
    CT[phaseID]->Length[0] = VoxelLength[0]*GridSize[0];
    CT[phaseID]->Length[1] = VoxelLength[1]*GridSize[1];
    CT[phaseID]->Length[2] = VoxelLength[2]*GridSize[2];

    CT[phaseID]->VoxelLength[0] = VoxelLength[0];
    CT[phaseID]->VoxelLength[1] = VoxelLength[1];
    CT[phaseID]->VoxelLength[2] = VoxelLength[2];

    CT[phaseID]->Origin[0] = Origin[0];
    CT[phaseID]->Origin[1] = Origin[1];
    CT[phaseID]->Origin[2] = Origin[2];
 
    CT[phaseID]->density = (VAR_DATA*)malloc(CT[phaseID]->Nbr_voxels * sizeof(VAR_DATA));
    CT[phaseID]->material = (unsigned short int*)malloc(CT[phaseID]->Nbr_voxels * sizeof(unsigned short int));

    #pragma omp parallel for private(i)
    for(i=0; i<CT[phaseID]->Nbr_voxels; i++){
      CT[phaseID]->density[i] = (VAR_DATA)HU_to_Density_convertion(hu[i], CT[phaseID]);
      CT[phaseID]->material[i] = (unsigned short int)Density_to_Material_convertion(CT[phaseID]->density[i], CT[phaseID]);
    }
    
    free(hu);
  }


  // Verify that all CT images have the same properties
  for(phaseID=0; phaseID<config->Num_4DCT_phases; phaseID++){
    if(CT[phaseID]->GridSize[0] != GridSize[0] || CT[phaseID]->GridSize[1] != GridSize[1] || CT[phaseID]->GridSize[2] != GridSize[2]){
      printf("\n Error: all phases of the imported 4DCT doesn't have the same grid size\n\n");
      Free_4DCT(CT, config->Num_4DCT_phases);
      return NULL;
    }
    if(CT[phaseID]->VoxelLength[0] != VoxelLength[0] || CT[phaseID]->VoxelLength[1] != VoxelLength[1] || CT[phaseID]->VoxelLength[2] != VoxelLength[2]){
      printf("\n Error: all phases of the imported 4DCT doesn't have the same voxel size\n\n");
      Free_4DCT(CT, config->Num_4DCT_phases);
      return NULL;
    }
    if(CT[phaseID]->Origin[0] != Origin[0] || CT[phaseID]->Origin[1] != Origin[1] || CT[phaseID]->Origin[2] != Origin[2]){
      printf("\n Error: all phases of the imported 4DCT doesn't have the same origin\n\n");
      Free_4DCT(CT, config->Num_4DCT_phases);
      return NULL;
    }
  }

  printf("\n");

  return CT;
}


VAR_DATA *Import_Def_Field(char *file_path, int *GridSize, VAR_DATA *Spacing, VAR_DATA *Origin){

  int i, j, k, l, m;
  
  VAR_DATA *Field_tmp = import_MHD_image(file_path, GridSize, Spacing, Origin);
  if(Field_tmp == NULL) return NULL;

  int NumVoxels = GridSize[0]*GridSize[1]*GridSize[2];
  VAR_DATA *Field = (VAR_DATA*)malloc(3*NumVoxels * sizeof(VAR_DATA));

  m = 0;
  for(j=0; j<GridSize[2]; j++){
    for(k=0; k<GridSize[1]; k++){
      for(l=0; l<GridSize[0]; l++){

	Field[m] = Field_tmp[3*m+1]/(10*Spacing[1]);
	Field[m+NumVoxels] = Field_tmp[3*m+0]/(10*Spacing[0]);
	Field[m+2*NumVoxels] = Field_tmp[3*m+2]/(10*Spacing[2]);

        m += 1;
  }}}

  free(Field_tmp);
 
  return Field;
}


DATA_4D_Fields *Import_4D_Fields(DATA_config *config){

  int i;
  char file_path[200];
  DATA_4D_Fields *Fields =  NULL;
  int Fields_GridSize[4];
  VAR_DATA Fields_Spacing[3], Fields_Origin[3];

  config->Num_4DCT_phases = 0;

  // count number of phases.
  sprintf(file_path, "./Fields/Field_Ref_to_phase1.mhd");
  FILE *file_mhd = fopen(file_path,"r");
  while(file_mhd != NULL){
    config->Num_4DCT_phases += 1;
    fclose(file_mhd);
    sprintf(file_path, "./Fields/Field_Ref_to_phase%d.mhd", config->Num_4DCT_phases+1);
    file_mhd = fopen(file_path,"r");
  }

  if(config->Num_4DCT_phases < 1){
    printf("\n Error: The number of fields to be imported must be >= 1.\n\n");
    return NULL;
  }

  printf("\n4D data: %d phases found in directory \"./Fields\"\n", config->Num_4DCT_phases);

  Fields = (DATA_4D_Fields*) malloc(sizeof(DATA_4D_Fields));
  Fields->Nbr_Fields = config->Num_4DCT_phases;
  Fields->Phase2Ref = (VAR_DATA**)malloc(config->Num_4DCT_phases * sizeof(VAR_DATA*));
  Fields->Ref2Phase = (VAR_DATA**)malloc(config->Num_4DCT_phases * sizeof(VAR_DATA*));
  for(i=0; i < config->Num_4DCT_phases; i++){
    Fields->Phase2Ref[i] = NULL;
    Fields->Ref2Phase[i] = NULL;
  }
  
  if(config->Field_type == 1){ // Field type == Velocity

    Fields->Ref2Phase_log = (VAR_DATA**)malloc(config->Num_4DCT_phases * sizeof(VAR_DATA*));

    for(i=0; i<config->Num_4DCT_phases; i++){
      printf(" Loading deformation field %d\n", i+1);

      sprintf(file_path, "./Fields/Field_Ref_to_phase%d.mhd", i+1);
      Fields->Ref2Phase_log[i] = Import_Def_Field(file_path, Fields_GridSize, Fields_Spacing, Fields_Origin);

      if(Fields->Ref2Phase_log[i] == NULL){
        Free_4D_Fields(Fields);
        return NULL;
      }
    }

    Fields->GridSize[0] = 3;
    Fields->GridSize[1] = Fields_GridSize[0];
    Fields->GridSize[2] = Fields_GridSize[1];
    Fields->GridSize[3] = Fields_GridSize[2];

    Fields->Spacing[0] = Fields_Spacing[0];
    Fields->Spacing[1] = Fields_Spacing[1];
    Fields->Spacing[2] = Fields_Spacing[2];

    Fields->Origin[0] = Fields_Origin[0];
    Fields->Origin[1] = Fields_Origin[1];
    Fields->Origin[2] = Fields_Origin[2];

    printf(" Fields exponentiation\n");
    
    for(i=0; i<config->Num_4DCT_phases; i++){
      Fields->Phase2Ref[i] = Field_exponentiation(Fields->Ref2Phase_log[i], Fields->GridSize, Fields->Spacing, Fields->Origin, 1);
      Fields->Ref2Phase[i] = Field_exponentiation(Fields->Ref2Phase_log[i], Fields->GridSize, Fields->Spacing, Fields->Origin, 0);
    }

  }
  else{  // Field type == Displacement

    Fields->Ref2Phase_log = NULL;

    for(i=0; i<config->Num_4DCT_phases; i++){
      printf(" Loading deformation field %d\n", i+1);

      sprintf(file_path, "./Fields/Field_phase%d_to_Ref.mhd", i+1);
      Fields->Phase2Ref[i] = Import_Def_Field(file_path, Fields_GridSize, Fields_Spacing, Fields_Origin);

      if(Fields->Phase2Ref[i] == NULL){
        Free_4D_Fields(Fields);
        return NULL;
      }

      sprintf(file_path, "./Fields/Field_Ref_to_phase%d.mhd", i+1);
      Fields->Ref2Phase[i] = Import_Def_Field(file_path, Fields_GridSize, Fields_Spacing, Fields_Origin);

      if(Fields->Ref2Phase[i] == NULL){
        Free_4D_Fields(Fields);
        return NULL;
      }
    }

    Fields->GridSize[0] = 3;
    Fields->GridSize[1] = Fields_GridSize[0];
    Fields->GridSize[2] = Fields_GridSize[1];
    Fields->GridSize[3] = Fields_GridSize[2];

    Fields->Spacing[0] = Fields_Spacing[0];
    Fields->Spacing[1] = Fields_Spacing[1];
    Fields->Spacing[2] = Fields_Spacing[2];

    Fields->Origin[0] = Fields_Origin[0];
    Fields->Origin[1] = Fields_Origin[1];
    Fields->Origin[2] = Fields_Origin[2];

  }

  printf("\n");

  return Fields;
}


void Free_4D_Fields(DATA_4D_Fields *Fields){

  if(Fields == NULL) return;

  int i;
  for(i=0; i<Fields->Nbr_Fields; i++){
    if(Fields->Phase2Ref[i] != NULL) free(Fields->Phase2Ref[i]);
    if(Fields->Ref2Phase[i] != NULL) free(Fields->Ref2Phase[i]);
  }

  if(Fields->Ref2Phase_log != NULL){
    for(i=0; i<Fields->Nbr_Fields; i++){
      if(Fields->Ref2Phase_log[i] != NULL) free(Fields->Ref2Phase_log[i]);
    }
  }

  if(Fields->Phase2Ref != NULL) free(Fields->Phase2Ref);
  if(Fields->Ref2Phase != NULL) free(Fields->Ref2Phase);
  if(Fields->Ref2Phase_log != NULL) free(Fields->Ref2Phase_log);

  return;
}


void Free_4DCT(DATA_CT **CT, int Nbr_phases){

  if(CT == NULL) return;

  int i;
  for(i=0; i<Nbr_phases; i++){
    if(i != 0){
      CT[i]->Conversion_HU_Density = NULL;
      CT[i]->Conversion_Densities = NULL;
      CT[i]->Conversion_HU_Material = NULL;
      CT[i]->Conversion_Density_Material = NULL;
      CT[i]->Conversion_Material_labels = NULL;
    }

    if(CT[i] != NULL) Free_CT_DATA(CT[i]);
  }

  free(CT);

  return;
}
