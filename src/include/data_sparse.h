/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_sparse
#define H_data_sparse

#include "define.h"
#include "struct.h"
#include "File_process.h"
#include "data_beam_model.h"


typedef struct DATA_Sparse_Header DATA_Sparse_Header;
struct DATA_Sparse_Header{
	time_t SimulationTimestamp;	
	char PlanName[100];
	unsigned int ImageSize[3];
	float VoxelSpacing[3];
	float Offset[3];
	char BinaryFile[100];
	int Beamlet_mode;
	int NbrSpots;
};


void export_Sparse_image(char *file_name, DATA_config *config, DATA_CT *ct, plan_parameters *plan, VAR_SCORING *data, VAR_SCORING threshold);
VAR_DATA *import_Sparse_image(char *file_name, int *GridSize, VAR_DATA *VoxelLength, VAR_DATA *Origin);
DATA_Sparse_Header Read_Sparse_Header(char *file_name);
void Display_Sparse_Header(DATA_Sparse_Header *Header);
DATA_Sparse_Header Init_Sparse_Header();
int Merge_Sparse_Files(char *InputPath, char *FileName, int NbrDirectories, char *OutputFile);

#endif
