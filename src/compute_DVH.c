/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_DVH.h"


void compute_all_DVH(DATA_config *config, VAR_SCORING *Dose, VAR_SCORING DoseScaling){

  char DVH_file_path[200];

  if(config->StructList == NULL) config->StructList = load_all_structs();

  int i;
  for(i=0; i<config->StructList->Nbr_Structs; i++){
    
    strcpy(DVH_file_path, config->Output_Directory);
    strcat(DVH_file_path, "DVH_");
    strcat(DVH_file_path, config->StructList->Structs[i].Name);
    strcat(DVH_file_path, config->output_robustness_suffix);
    strcat(DVH_file_path, config->output_beamlet_suffix);
    strcat(DVH_file_path, config->output_4D_suffix);
    strcat(DVH_file_path, config->output_beams_suffix);
    strcat(DVH_file_path, ".txt");

    compute_DVH(config->StructList->Structs[i].IndexList, config->StructList->Structs[i].N_Index, Dose, DoseScaling, DVH_file_path);
  }
    
}


void compute_DVH(int *ListIndex, int ListSize, VAR_SCORING *Dose, VAR_SCORING DoseScaling, char *Output_file){

  int NbrBins = 4096;
  float interval[2] = {0, 100};

  float bin_size = (interval[1]-interval[0])/NbrBins;
  float *DVH = (float*)calloc(NbrBins, sizeof(float));

  int i, bin;
  for(i=0; i<ListSize; i++){
    bin = floor(Dose[ListIndex[i]] * DoseScaling / bin_size);
    if(bin >= NbrBins) bin = NbrBins-1;
    DVH[bin] += 1;
  }

  for(i=NbrBins-2; i>=0; i--){
    DVH[i] += DVH[i+1];
  }

  for(i=NbrBins-1; i>=0; i--){
    DVH[i] = 100*DVH[i]/DVH[0];
  }


  FILE *file = NULL;
  file = fopen(Output_file, "w");
  for(i=0; i<NbrBins; i++){
    fprintf(file, "%f\t%lf\n", (i+1)*bin_size, DVH[i]);
  }
  fclose(file);

  if(DVH != NULL) free(DVH);

  return;
}
