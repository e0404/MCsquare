/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_ct_penelope.h"

DATA_CT *Read_PENCT(char *file_name){
  FILE *file = NULL;
  char *read_line, *read_data;
  int line = 0, CT_dim = 0;
  int CT_Nx = -1, CT_Ny = -1, CT_Nz = -1;
  VAR_DATA LengthX = -1, LengthY = -1, LengthZ = -1;
  DATA_CT *ct_struct = NULL;

  file = fopen(file_name, "r");
  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", file_name);
    return NULL;
  }

  read_line = (char*) malloc(1000 * sizeof(char));

  while (fgets(read_line, 1000, file) != NULL){

    // on ignore les commentaires
    if(read_line[0] == '#') continue;
    strtok(read_line, "#");	

    read_data = strtok(read_line, " \t\r\n");
    if(read_data == NULL) continue;

    line++;

    if(line == 1){
      LengthX = atof(read_data);
      read_data = strtok(NULL, " \t\r\n");
      LengthY = atof(read_data);
      read_data = strtok(NULL, " \t\r\n");
      LengthZ = atof(read_data);
      continue;
    }
    else if(line == 2){
      CT_Nx = atoi(read_data);
      read_data = strtok(NULL, " \t\r\n");
      CT_Ny = atoi(read_data);
      read_data = strtok(NULL, " \t\r\n");
      CT_Nz = atoi(read_data);
      continue;
    }
    else if(line == 3){
      if(CT_Nx == -1 || CT_Ny == -1 || CT_Nz == -1 || LengthX == -1 || LengthY == -1 || LengthZ == -1){
        printf("\n\nError: Invalid CT header!\n\n");
        return NULL;
      }
      ct_struct = (DATA_CT*) malloc(sizeof(DATA_CT));
      Init_CT_DATA(ct_struct, CT_Nx, CT_Ny, CT_Nz, LengthX, LengthY, LengthZ);
      CT_dim = CT_Nx * CT_Ny * CT_Nz;
    }

    if(line >= 3){
      if(CT_dim < line-2) break;
      ct_struct->material[line-3] = atoi(read_data);
      read_data = strtok(NULL, " \t\r\n");
      ct_struct->density[line-3] = atof(read_data);
    }
  }

  fclose(file);
  return ct_struct;
}

