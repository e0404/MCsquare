/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_4D
#define H_data_4D

#include "define.h"
#include "struct.h"
#include "data_mhd.h"
#include "data_ct.h"
#include "compute_math.h"
#include "compute_4D.h"


DATA_CT **Import_4DCT(DATA_config *config);
VAR_DATA *Import_Def_Field(char *file_path, int *GridSize, VAR_DATA *Spacing, VAR_DATA *Origin);
DATA_4D_Fields *Import_4D_Fields(DATA_config *config);
void Free_4D_Fields(DATA_4D_Fields *Fields);
void Free_4DCT(DATA_CT **CT, int Nbr_phases);

#endif
