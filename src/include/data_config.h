/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_config
#define H_data_config

#include "define.h"
#include "struct.h"
#include "File_process.h"
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

union Config_Value
{
    int int_value;
    unsigned int uint_value;
    unsigned long ulong_value;
    VAR_DATA real_value;
};

union Config_Adress
{
    int *int_adr;
    unsigned int *uint_adr;
    unsigned long *ulong_adr;
    VAR_DATA *real_adr;
    char *string_adr;
};

typedef struct DATA_config_dictionary DATA_config_dictionary;
struct DATA_config_dictionary{
	char type[10];
	int is_defined;
	int is_default;
	char Tag[40];
	char List[100];
	union Config_Adress adress;
	union Config_Value min_value;
	union Config_Value max_value;
};


DATA_config_dictionary *Init_Config(DATA_config *config);
void Add_bool_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned int *config, int use_default, unsigned int default_value);
void Add_string_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, char *config, int use_default, char *default_value, unsigned int min_size,  unsigned int max_size);
void Add_uint_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned int *config, int use_default, unsigned int default_value, unsigned int min_value,  unsigned int max_value);
void Add_ulong_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, unsigned long *config, int use_default, unsigned long default_value, unsigned long min_value,  unsigned long max_value);
void Add_ureal_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, VAR_DATA *config, int use_default, VAR_DATA default_value, VAR_DATA min_value,  VAR_DATA max_value);
void Add_vec_ureal_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, VAR_DATA *config, int use_default, VAR_DATA default_value_X, VAR_DATA default_value_Y, VAR_DATA default_value_Z, VAR_DATA min_value,  VAR_DATA max_value);
void Add_Enum_Config_element(char *Tag, DATA_config_dictionary *config_dictionary, int *config, int use_default, int default_value, char *List);
int Parse_Config(DATA_config *config, char *file_name);
void display_config(DATA_config *config);

#endif
