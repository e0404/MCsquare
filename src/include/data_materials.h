/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_materials
#define H_data_materials

#include "define.h"
#include "struct.h"
#include "File_process.h"
#include "data_Stop_Pow.h"
#include "data_nuclear.h"


Materials *Init_materials(unsigned int *Nbr_Materials, unsigned int *Max_Components, DATA_config *config);
Materials *List_materials(unsigned int *Nbr_Materials, DATA_config *config);
void Free_Materials_DATA(Materials *material, unsigned int Nbr_Materials);

#endif
