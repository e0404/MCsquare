/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_ascii.h"

void export_dose_ascii(char *file_name, int GridSize[3], VAR_SCORING *data){

  FILE *file = NULL;
  file = fopen(file_name, "w");

  int i,j,k,index=0;
  for(i=1; i<=GridSize[2]; i++){
    for(j=1; j<=GridSize[1]; j++){
      for(k=1; k<=GridSize[0]; k++){
        fprintf(file, "%d\t%d\t%d\t%lf\n", k, j, i, data[index]);
        index++;
      }
    }
  }

  fclose(file);

  return;
}


void export_PG_ascii(char *file_name, int GridSize[3], VAR_SCORING *data){

  FILE *file = NULL;
  file = fopen(file_name, "w");

  int i,j,k,index=0;
  for(i=1; i<=GridSize[2]; i++){
    for(j=1; j<=GridSize[1]; j++){
      for(k=1; k<=GridSize[0]; k++){
        fprintf(file, "%d\t%d\t%d\t%.6e\n", k, j, i, data[index]);
        index++;
      }
    }
  }

  fclose(file);

  return;
}


void export_PG_spectrum_ascii(char *file_name, int NumBin, VAR_DATA Binning, VAR_SCORING *data){

  FILE *file = NULL;
  file = fopen(file_name, "w");

  int i;
  for(i=0; i<=NumBin; i++){
    fprintf(file, "%f\t%f\n", i*Binning, data[i]);
  }

  fclose(file);

  return;
}
