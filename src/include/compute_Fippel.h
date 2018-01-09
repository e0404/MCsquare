/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_Fippel
#define H_compute_Fippel

#include "define.h"
#include "struct.h"
#include "compute_math.h"
#include "compute_random.h"
#include "compute_EM_interaction.h"
#include "compute_geometry.h"

#define CONST_DERIV 0.99 	// Constante pour le calcul des dérivées


void Fippel_Stop_Pow_correction(Hadron *hadron, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result);
void Compute_dE2_Fippel(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result);

#endif
