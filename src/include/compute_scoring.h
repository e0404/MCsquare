/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_scoring
#define H_compute_scoring

#include "define.h"
#include "struct.h"
#include "compute_EM_interaction.h"

DATA_Scoring Init_Scoring(DATA_config *config, int Nbr_voxels, int init_dose);
void Energy_Scoring(DATA_Scoring *scoring, int index, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE SPR);
void LET_Scoring(DATA_Scoring *scoring, int index, VAR_COMPUTE multiplicity, VAR_COMPUTE dE, VAR_COMPUTE step, VAR_COMPUTE stop_pow, DATA_config *config);
void PostProcess_Scoring(DATA_Scoring *scoring, DATA_CT *ct, Materials *material, VAR_COMPUTE normalization, unsigned long Nbr_simulated_primaries, DATA_config *config);
void Free_Scoring(DATA_Scoring *scoring);

#endif
