/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_treatment_uncertainties
#define H_compute_treatment_uncertainties

#include "define.h"
#include "struct.h"
#include "compute_random.h"
#include "data_beam_model.h"
#include "compute_beam_model.h"
#include "compute_4D.h"
#include "compute_math.h"
#include "data_ct.h"

void Translation_uncertainty(Hadron_buffer *hadron, DATA_config *config, VSLStreamStatePtr RNG_Stream);
void Density_scaling(VAR_DATA *Nominal_density, VAR_DATA *Scaled_density, int Num_voxels, VAR_DATA scaling_factor);
void Breathing_amplitude_variation(DATA_config *config, DATA_CT *ct, DATA_CT **CT_phases, DATA_4D_Fields *Fields);
plan_parameters* Spot_Sorting(DATA_config *config, int phase, plan_parameters *Plan);

#endif
