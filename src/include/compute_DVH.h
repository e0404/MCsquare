/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_DVH
#define H_compute_DVH

#include "define.h"
#include "struct.h"
#include "data_mhd.h"
#include "data_sparse.h"

#if defined(_MSC_VER)
  #include "lib/win_dirent.h"
#else
  #include <dirent.h>
#endif


void compute_all_DVH(DATA_config *config, VAR_SCORING *Dose, VAR_SCORING DoseScaling);
int *import_mask(char *file_name, int *ListSize);
void compute_DVH(int *ListIndex, int ListSize, VAR_SCORING *Dose, VAR_SCORING DoseScaling, char *Output_file);

#endif
