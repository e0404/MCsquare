/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_geometry
#define H_compute_geometry

#include "define.h"
#include "struct.h"
#include "compute_math.h"
#include "compute_EM_interaction.h"

void verif_position(Hadron *hadron, DATA_CT *ct);
void get_CT_Offset(Hadron *hadron, DATA_CT *ct, int *v_index);
void Dist_To_Material_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE dist_max, int *v_init_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_result);
void Dist_To_Interface(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_step);
void Update_position(Hadron *hadron, VAR_COMPUTE *v_step);
void Update_direction(Hadron *hadron, VAR_COMPUTE *v_theta, VAR_COMPUTE *v_phi);
void Update_buffer_direction(Hadron_buffer *secondary_hadron, VAR_COMPUTE theta, VAR_COMPUTE phi);
void CT_Transport(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density);
void CT_Transport_SPR(Hadron *hadron, DATA_CT *ct, Materials *material, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density);
void CT_Transport_Random_Hinge(Hadron *hadron, DATA_CT *ct, VAR_COMPUTE *v_s, VAR_COMPUTE *v_tau, int *v_init_index, int *v_hinge_index, VAR_COMPUTE *v_init_density, VAR_COMPUTE *v_mask);

#endif
