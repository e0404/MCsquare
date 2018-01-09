/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_4D
#define H_compute_4D

#include "define.h"
#include "struct.h"
#include "compute_math.h"
#include "data_ct.h"

VAR_COMPUTE *Image_deformation(VAR_DATA *image, int *image_size, VAR_DATA *image_spacing, VAR_DATA *image_origin, VAR_DATA *field, int *field_size, VAR_DATA *field_spacing, VAR_DATA *field_origin);
VAR_DATA *Field_multiplication(VAR_DATA Scalar, VAR_DATA *InitField, int *GridSize);
void Field_addition(VAR_DATA *Field1, VAR_DATA *Field2, int *GridSize);
VAR_DATA *Field_exponentiation(VAR_DATA *InitField, int *GridSize, VAR_DATA *Spacing, VAR_DATA *Origin, int inverse);
DATA_CT **Create_4DCT_from_Ref(DATA_config *config, DATA_4D_Fields *Fields, DATA_CT *ct);
DATA_CT *Create_Ref_from_4DCT(DATA_config *config, DATA_4D_Fields *Fields, DATA_CT **CT_phases);

#endif
