/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_ct_mhd.h"

DATA_CT *Read_CT_MHD(DATA_config *config){

  int GridSize[3];
  VAR_DATA VoxelLength[3], Origin[3], *hu;

  DATA_CT *CT = (DATA_CT*) malloc(sizeof(DATA_CT));

//  CT->density = import_MHD_image(config->CT_File, GridSize, VoxelLength);
//  if(CT->density == NULL) return NULL;

  hu = import_MHD_image(config->CT_File, GridSize, VoxelLength, Origin);
  if(hu == NULL) return NULL;

  CT->GridSize[0] = GridSize[0];
  CT->GridSize[1] = GridSize[1];
  CT->GridSize[2] = GridSize[2];

  CT->Nbr_voxels = GridSize[0]*GridSize[1]*GridSize[2];
  
  CT->Length[0] = VoxelLength[0]*GridSize[0];
  CT->Length[1] = VoxelLength[1]*GridSize[1];
  CT->Length[2] = VoxelLength[2]*GridSize[2];

  CT->VoxelLength[0] = VoxelLength[0];
  CT->VoxelLength[1] = VoxelLength[1];
  CT->VoxelLength[2] = VoxelLength[2];

  CT->Origin[0] = Origin[0];
  CT->Origin[1] = Origin[1];
  CT->Origin[2] = Origin[2];

  CT->Conversion_HU_Density = NULL;
  CT->Conversion_Densities = NULL;
  CT->Conversion_HU_Material = NULL;
  CT->Conversion_Density_Material = NULL;
  CT->Conversion_Material_labels = NULL;

  if(Read_Density_conversion_data(config->HU_Density_File, CT) != 0){
    Free_CT_DATA(CT);
    free(hu);
    return NULL;
  }
//  Display_Density_conversion_data(CT);

  if(Read_Material_conversion_data(config->HU_Material_File, CT) != 0){
    Free_CT_DATA(CT);
    free(hu);
    return NULL;
  }
//  Display_Material_conversion_data(CT);


  CT->density = (VAR_DATA*)malloc(CT->Nbr_voxels * sizeof(VAR_DATA));
  CT->material = (unsigned short int*)malloc(CT->Nbr_voxels * sizeof(unsigned short int));

  int i;

  #pragma omp parallel for private(i)
  for(i=0; i<CT->Nbr_voxels; i++){
    CT->density[i] = (VAR_DATA)HU_to_Density_convertion(hu[i], CT);
    //CT->material[i] = (unsigned short int)Density_to_Material_convertion(CT->density[i], CT);
    CT->material[i] = (unsigned short int)HU_to_Material_convertion(hu[i], CT);
  }

  free(hu);

  return CT;
}
