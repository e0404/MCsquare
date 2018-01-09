/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_4D.h"
#include "include/data_mhd.h"


VAR_COMPUTE *Image_deformation(VAR_DATA *image, int *image_size, VAR_DATA *image_spacing, VAR_DATA *image_origin, VAR_DATA *field, int *field_size, VAR_DATA *field_spacing, VAR_DATA *field_origin){

  int j,k,l,m;

  int ImageVoxels = image_size[0]*image_size[1]*image_size[2];
  int FieldVoxels = field_size[1]*field_size[2]*field_size[3];

  VAR_COMPUTE *deformed = (VAR_COMPUTE*)malloc(ImageVoxels * sizeof(VAR_COMPUTE));
  VAR_COMPUTE Deformation[3];

  VAR_COMPUTE Spacing_Ratio[3];
  Spacing_Ratio[0] = image_spacing[0] / field_spacing[0];
  Spacing_Ratio[1] = image_spacing[1] / field_spacing[1];
  Spacing_Ratio[2] = image_spacing[2] / field_spacing[2];

  VAR_COMPUTE Origin[3];
  Origin[0] = (image_origin[0] - field_origin[0]) / field_spacing[0];
  Origin[1] = -(image_size[1]*image_spacing[1] - field_size[2]*field_spacing[1] - field_origin[1] + image_origin[1]) / field_spacing[1];
  Origin[2] = -(image_size[2]*image_spacing[2] - field_size[3]*field_spacing[2] - field_origin[2] + image_origin[2]) / field_spacing[2];

  
  if(Spacing_Ratio[0] == 1.0 && Spacing_Ratio[1] == 1.0 && Spacing_Ratio[2] == 1.0 && Origin[0] == 0.0 && Origin[1] == 0.0 && Origin[2] == 0.0){

    #pragma omp parallel for private(j,k,l,m, Deformation)
    for(j=0; j<image_size[2]; j++){
      for(k=0; k<image_size[1]; k++){
        for(l=0; l<image_size[0]; l++){
	  m = l + k*image_size[0] + j*image_size[0]*image_size[1];
	  Deformation[0] = l + field[m];
	  Deformation[1] = k + field[m+FieldVoxels];
	  Deformation[2] = j + field[m+2*FieldVoxels];
          deformed[m] = (VAR_COMPUTE)Trilinear_Interpolation(Deformation, image, image_size);
    }}}
  }
  else{

    VAR_COMPUTE FieldCoordinates[3];

    #pragma omp parallel for private(j,k,l,m, Deformation, FieldCoordinates)
    for(j=0; j<image_size[2]; j++){
      for(k=0; k<image_size[1]; k++){
        for(l=0; l<image_size[0]; l++){

	  m = l + k*image_size[0] + j*image_size[0]*image_size[1];

	  FieldCoordinates[0] = Origin[0] + l*Spacing_Ratio[0];
	  FieldCoordinates[1] = Origin[1] + k*Spacing_Ratio[1];
	  FieldCoordinates[2] = Origin[2] + j*Spacing_Ratio[2];

	  Deformation[0] = l + (VAR_COMPUTE)Trilinear_Interpolation(FieldCoordinates, field, &field_size[1]) / Spacing_Ratio[0];
	  Deformation[1] = k + (VAR_COMPUTE)Trilinear_Interpolation(FieldCoordinates, &field[FieldVoxels], &field_size[1]) / Spacing_Ratio[1];
	  Deformation[2] = j + (VAR_COMPUTE)Trilinear_Interpolation(FieldCoordinates, &field[2*FieldVoxels], &field_size[1]) / Spacing_Ratio[2];
          deformed[m] = (VAR_COMPUTE)Trilinear_Interpolation(Deformation, image, image_size);
    }}}

  }

  return deformed;
}



VAR_DATA *Field_multiplication(VAR_DATA Scalar, VAR_DATA *InitField, int *GridSize){

  int m;

  int Elements = GridSize[0]*GridSize[1]*GridSize[2]*GridSize[3];
  
  VAR_DATA *Field = (VAR_DATA*)malloc(Elements * sizeof(VAR_DATA));
  
  #pragma omp parallel for private(m)
  for(m=0; m<Elements; m++) Field[m] = Scalar*InitField[m];
  
  return Field;
}


void Field_addition(VAR_DATA *Field1, VAR_DATA *Field2, int *GridSize){

  int m;

  int Elements = GridSize[0]*GridSize[1]*GridSize[2]*GridSize[3];
  
  
  #pragma omp parallel for private(m)
  for(m=0; m<Elements; m++) Field1[m] += Field2[m];
  
  return;
}


VAR_DATA *Field_exponentiation(VAR_DATA *InitField, int *GridSize, VAR_DATA *Spacing, VAR_DATA *Origin, int inverse){

  int r;
  int N = 6;
  int NumVoxels = GridSize[1]*GridSize[2]*GridSize[3];

  int ComponentSize[4];
  ComponentSize[0] = 1;
  ComponentSize[1] = GridSize[1];
  ComponentSize[2] = GridSize[2];
  ComponentSize[3] = GridSize[3];
  
  VAR_DATA Factor = pow(2, -N);
  if(inverse != 0) Factor *= -1.0;
  VAR_DATA *Field = Field_multiplication(Factor, InitField, GridSize);

  VAR_DATA *Field_x, *Field_y, *Field_z;
  
  for(r=1; r<=N; r++){
    Field_x = Image_deformation(Field, &GridSize[1], Spacing, Origin, Field, GridSize, Spacing, Origin);
    Field_y = Image_deformation(&Field[NumVoxels], &GridSize[1], Spacing, Origin, Field, GridSize, Spacing, Origin);
    Field_z = Image_deformation(&Field[2*NumVoxels], &GridSize[1], Spacing, Origin, Field, GridSize, Spacing, Origin);

    Field_addition(Field, Field_x, ComponentSize);
    free(Field_x);

    Field_addition(&Field[NumVoxels], Field_y, ComponentSize);
    free(Field_y);

    Field_addition(&Field[2*NumVoxels], Field_z, ComponentSize);
    free(Field_z);
  }
  
  return Field;
}


DATA_CT **Create_4DCT_from_Ref(DATA_config *config, DATA_4D_Fields *Fields, DATA_CT *ct){

  printf("\nCreating 4DCT images from Reference phase\n");

  DATA_CT **CT =  (DATA_CT**)malloc(config->Num_4DCT_phases * sizeof(DATA_CT*));

  int phaseID, i;
  for(phaseID=0; phaseID < config->Num_4DCT_phases; phaseID++){

    printf(" Phase %d\n", phaseID+1);

    CT[phaseID] = (DATA_CT*) malloc(sizeof(DATA_CT));

    CT[phaseID]->Conversion_HU_Density = ct->Conversion_HU_Density;
    CT[phaseID]->Conversion_Densities = ct->Conversion_Densities;
    CT[phaseID]->Conversion_HU_Material = ct->Conversion_HU_Material;
    CT[phaseID]->Conversion_Density_Material = ct->Conversion_Density_Material;
    CT[phaseID]->Conversion_Material_labels = ct->Conversion_Material_labels;
    CT[phaseID]->Num_Density_Data = ct->Num_Density_Data;
    CT[phaseID]->Num_Materials_Data = ct->Num_Materials_Data;	

    CT[phaseID]->GridSize[0] = ct->GridSize[0];
    CT[phaseID]->GridSize[1] = ct->GridSize[1];
    CT[phaseID]->GridSize[2] = ct->GridSize[2];

    CT[phaseID]->Nbr_voxels = ct->Nbr_voxels;

    CT[phaseID]->Length[0] = ct->Length[0];
    CT[phaseID]->Length[1] = ct->Length[1];
    CT[phaseID]->Length[2] = ct->Length[2];

    CT[phaseID]->VoxelLength[0] = ct->VoxelLength[0];
    CT[phaseID]->VoxelLength[1] = ct->VoxelLength[1];
    CT[phaseID]->VoxelLength[2] = ct->VoxelLength[2];

    CT[phaseID]->Origin[0] = ct->Origin[0];
    CT[phaseID]->Origin[1] = ct->Origin[1];
    CT[phaseID]->Origin[2] = ct->Origin[2];
 
    CT[phaseID]->density = Image_deformation(ct->density, ct->GridSize, ct->VoxelLength, ct->Origin, Fields->Ref2Phase[phaseID], Fields->GridSize, Fields->Spacing, Fields->Origin);
    CT[phaseID]->material = (unsigned short int*)malloc(CT[phaseID]->Nbr_voxels * sizeof(unsigned short int));

    #pragma omp parallel for private(i)
    for(i=0; i<CT[phaseID]->Nbr_voxels; i++){
      CT[phaseID]->material[i] = (unsigned short int)Density_to_Material_convertion(CT[phaseID]->density[i], CT[phaseID]);
    }
  }

  return CT;
  
}


DATA_CT *Create_Ref_from_4DCT(DATA_config *config, DATA_4D_Fields *Fields, DATA_CT **CT_phases){

  printf("\nCreating reference phase image from 4DCT\n");

  DATA_CT *ct = (DATA_CT*) malloc(sizeof(DATA_CT));

  ct->Conversion_HU_Density = CT_phases[0]->Conversion_HU_Density;
  ct->Conversion_Densities = CT_phases[0]->Conversion_Densities;
  ct->Conversion_HU_Material = CT_phases[0]->Conversion_HU_Material;
  ct->Conversion_Density_Material = CT_phases[0]->Conversion_Density_Material;
  ct->Conversion_Material_labels = CT_phases[0]->Conversion_Material_labels;
  ct->Num_Density_Data = CT_phases[0]->Num_Density_Data;
  ct->Num_Materials_Data = CT_phases[0]->Num_Materials_Data;	

  ct->GridSize[0] = CT_phases[0]->GridSize[0];
  ct->GridSize[1] = CT_phases[0]->GridSize[1];
  ct->GridSize[2] = CT_phases[0]->GridSize[2];

  ct->Nbr_voxels = CT_phases[0]->Nbr_voxels;

  ct->Length[0] = CT_phases[0]->Length[0];
  ct->Length[1] = CT_phases[0]->Length[1];
  ct->Length[2] = CT_phases[0]->Length[2];

  ct->VoxelLength[0] = CT_phases[0]->VoxelLength[0];
  ct->VoxelLength[1] = CT_phases[0]->VoxelLength[1];
  ct->VoxelLength[2] = CT_phases[0]->VoxelLength[2];

  ct->Origin[0] = CT_phases[0]->Origin[0];
  ct->Origin[1] = CT_phases[0]->Origin[1];
  ct->Origin[2] = CT_phases[0]->Origin[2];

  ct->density = (VAR_COMPUTE*)malloc(ct->Nbr_voxels * sizeof(VAR_COMPUTE));
  ct->material = (unsigned short int*)malloc(ct->Nbr_voxels * sizeof(unsigned short int));


  VAR_COMPUTE **deformed = (VAR_COMPUTE**)malloc(config->Num_4DCT_phases * sizeof(VAR_COMPUTE*));

  int phaseID;
  for(phaseID=0; phaseID < config->Num_4DCT_phases; phaseID++){
    deformed[phaseID] = Image_deformation(CT_phases[phaseID]->density, CT_phases[phaseID]->GridSize, CT_phases[phaseID]->VoxelLength, CT_phases[phaseID]->Origin, Fields->Phase2Ref[phaseID], Fields->GridSize, Fields->Spacing, Fields->Origin);
  }

  int i;
  #pragma omp parallel for private(i, phaseID)
  for(i=0; i<ct->Nbr_voxels; i++){

    VAR_COMPUTE tmp[20];
    for(phaseID=0; phaseID < config->Num_4DCT_phases; phaseID++) tmp[phaseID] = deformed[phaseID][i];

    ct->density[i] = compute_median(tmp, config->Num_4DCT_phases);
    ct->material[i] = (unsigned short int)Density_to_Material_convertion(ct->density[i], ct);
  }

  // clean memory
  for(phaseID=0; phaseID < config->Num_4DCT_phases; phaseID++) free(deformed[phaseID]);
  free(deformed);

  return ct;
}




