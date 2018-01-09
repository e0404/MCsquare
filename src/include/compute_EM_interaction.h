/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_EM_interaction
#define H_compute_EM_interaction

#include "define.h"
#include "struct.h"
#include "compute_math.h"
#include "compute_random.h"
#include "compute_nuclear_interaction.h"

#define CONST_DERIV 0.99 	// Constante pour le calcul des dérivées

#define CONST_MS_Fippel 18.3 	// Constante pour le modèle du Multiple Scattering de Fippel
				// 12.4 MeV pour GEANT4 v5.1 (article Fippel)
				// 12.0 MeV pour GEANT4 v5.2 (article Fippel)
				// 12.9 MeV pour FLUKA2002.4 (article Fippel)
				// 20.3 MeV pour PENELOPE (ce code)
				// 16.5 MeV pour GEANT4 (ce code)
				// 18.3 MeV tuned by Sheng Huang at UPenn facility (October 2016)

void Total_Stop_Pow(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_stop_pow);
void Total_Hard_Cross_Section(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, DATA_config *config, VAR_COMPUTE *v_result);
void get_interaction_type(Hadron *hadron, Materials *material, int *v_material_label, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, VAR_COMPUTE Te_min, VAR_COMPUTE *v_dE_max, VAR_COMPUTE *v_tot_section, VSLStreamStatePtr RNG_Stream, DATA_config *config, int *v_result);
void cross_section_ionization(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_result);
void Compute_L(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_result);
void Compute_dE2(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE *v_density, Materials *material, VAR_COMPUTE Te_min, int *v_material_label, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result);
void Compute_Energy_straggling(Hadron *hadron, VAR_COMPUTE *v_N_el, VAR_COMPUTE Te_min, VAR_COMPUTE *v_s, VAR_COMPUTE *v_result);
void Compute_MS_Fippel(Hadron *hadron, VAR_COMPUTE *v_s, VAR_COMPUTE *v_X0, VAR_COMPUTE *v_result);
void Compute_Ionization_Energy(Hadron *hadron, VAR_COMPUTE Te_min, VSLStreamStatePtr RNG_Stream, VAR_COMPUTE *v_result);

#endif
