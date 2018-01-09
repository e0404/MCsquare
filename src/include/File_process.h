/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_File_Process
#define H_File_Process

#include "define.h"
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>


int File_exists(const char *FileName);
double get_units(char *str);
int isUnsignedInt(char *str);
int isUnsignedFloat(char *str);
int isFloat(char *str);
int isBoolean(char *str);
unsigned int getBoolean(char *str);
void CreateDir(char *DirName);
int CopyFile(char *to, char *from);
void str_replace(char *search , char *replace , char *subject);

#endif
