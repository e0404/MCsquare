/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_math
#define H_compute_math

#include "define.h"
#include <stdint.h>

void vector_floor(VAR_COMPUTE *v_vec, VAR_COMPUTE *v_result);
int sign(VAR_COMPUTE x);
void vec_sign(VAR_COMPUTE *v_vec, VAR_COMPUTE *v_sign);
VAR_COMPUTE SolidAngle(VAR_COMPUTE theta);
VAR_COMPUTE Linear_Interpolation(VAR_COMPUTE x, VAR_COMPUTE x1, VAR_COMPUTE x2, VAR_COMPUTE y1, VAR_COMPUTE y2);
void vec_Linear_Interpolation(VAR_COMPUTE *v_x, VAR_COMPUTE *v_x1, VAR_COMPUTE *v_x2, VAR_COMPUTE *v_y1, VAR_COMPUTE *v_y2, VAR_COMPUTE *v_result);
VAR_COMPUTE Trilinear_Interpolation(VAR_DATA *x, VAR_DATA *image, int *GridSize);
int Sequential_Search(VAR_COMPUTE value, VAR_DATA *list, int Nbr_Values);
int Binary_Search(VAR_COMPUTE value, VAR_DATA *list, int Nbr_Values);
void my_log(VAR_COMPUTE *v_data, VAR_COMPUTE *v_result);
VAR_COMPUTE find_min(VAR_COMPUTE *data, int NumElements);
VAR_COMPUTE find_max(VAR_COMPUTE *data, int NumElements);
VAR_COMPUTE compute_mean(VAR_COMPUTE *data, int NumElements);
VAR_COMPUTE compute_median(VAR_COMPUTE *data, int NumElements);

#endif
