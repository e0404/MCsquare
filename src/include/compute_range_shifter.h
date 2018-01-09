/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_range_shifter
#define H_compute_range_shifter

#include "define.h"
#include "struct.h"
#include "data_beam_model.h"
#include "compute_random.h"
#include "compute_semi_infinite_slab.h"

int Init_RangeShifter_Data(plan_parameters *plan, machine_parameters *machine, Materials *material, DATA_config *config);
void Display_RangeShifter_Data(plan_parameters *plan, machine_parameters *machine, Materials *material);
void Simulate_RangeShifter(Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Nbr_hadrons, DATA_config *config, machine_parameters *machine, Materials *material, VSLStreamStatePtr RNG_Stream);

#endif
