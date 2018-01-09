/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_semi_infinite_slab
#define H_compute_semi_infinite_slab

#include "define.h"
#include "struct.h"
#include "compute_geometry.h"
#include "compute_EM_interaction.h"
#include "compute_nuclear_interaction.h"
#include "compute_random.h"
#include "data_beam_model.h"

void SemiInfiniteSlab_step(Hadron *hadron, Materials *material, Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Hadron_ID, int *Nbr_hadrons, VAR_COMPUTE *RS_exit_position, DATA_config *config, machine_parameters *machine, VSLStreamStatePtr RNG_Stream);

#endif
