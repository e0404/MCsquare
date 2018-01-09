/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_Random
#define H_Random

#include "define.h"

void rand_uniform(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd);
VAR_COMPUTE single_rand_uniform(VSLStreamStatePtr stream);
void rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_mu, VAR_COMPUTE *v_sigma);
void rand_normal_zero(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_sigma);
VAR_COMPUTE single_rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE mu, VAR_COMPUTE sigma);

#endif
