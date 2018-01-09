/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_nuclear_interaction
#define H_compute_nuclear_interaction

#include "define.h"
#include "struct.h"
#include "compute_math.h"
#include "compute_random.h"
#include "compute_geometry.h"

void proton_proton_cross_section(Hadron *hadron, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result);
void total_Nuclear_cross_section(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_density, VAR_COMPUTE *v_result);
VAR_COMPUTE Compute_Nuclear_interaction(int hadron_index, Hadron *hadron, Materials *material, int material_label, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, int scoring_index, DATA_Scoring *scoring, VSLStreamStatePtr RNG_Stream, DATA_config *config);
VAR_COMPUTE Compute_Elastic_PP(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, VSLStreamStatePtr RNG_Stream, DATA_config *config);
VAR_COMPUTE Compute_Elastic_ICRU(int hadron_index, Hadron *hadron, Materials *material, VSLStreamStatePtr RNG_Stream);
VAR_COMPUTE Compute_Nuclear_Inelastic_recoils(VAR_COMPUTE Hadron_T, Materials *material, int index);
VAR_COMPUTE Compute_Nuclear_Inelastic_proton(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config);
VAR_COMPUTE Compute_Nuclear_Inelastic_deuteron(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config);
VAR_COMPUTE Compute_Nuclear_Inelastic_alpha(int hadron_index, Hadron *hadron, Hadron_buffer *secondary_hadron, int *Nbr_secondaries, Materials *material, int index, VSLStreamStatePtr RNG_Stream, DATA_config *config);
void Compute_PromptGamma(int hadron_index, Hadron *hadron, Materials *material, int scoring_index, DATA_Scoring *scoring, VSLStreamStatePtr RNG_Stream, DATA_config *config);

#endif
