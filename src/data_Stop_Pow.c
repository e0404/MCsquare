/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_Stop_Pow.h"

int Read_Stop_Pow(char *MaterialName, Materials *material, DATA_config *config){

  FILE *file = NULL;
  VAR_DATA tmp1, tmp2;
  int i = 0, nbr_lines = 0;

  char FileName[200];

  strcpy(FileName, config->Materials_Dir);
  strcat(FileName, MaterialName);

#if DB_STOP_POW==SP_PSTAR
  strcat(FileName, "/PSTAR_Stop_Pow.dat");
#else
  strcat(FileName, "/G4_Stop_Pow.dat");
#endif

  file = fopen(FileName, "r");

  if(file == NULL){
    printf("\n\n ERROR: File \"%s\" is missing! \n\n", FileName);
    return 1;
  }

#if VAR_DATA_PRECISION==1
  while(fscanf(file, "%f\t%f", &tmp1, &tmp2) > 0){
#else
  while(fscanf(file, "%lf\t%lf", &tmp1, &tmp2) > 0){
#endif

    nbr_lines++;					// compte le nombre de ligne pour pouvoir allouer la mémoire dynamiquement
  }
  //printf("\n\n Nbr lignes: %d\n\n", nbr_lines);

  rewind(file);	// retour au début du fichier

  material->SP_Energy = (VAR_DATA*) malloc(nbr_lines * sizeof(VAR_DATA));
  material->Stop_Pow = (VAR_DATA*) malloc(nbr_lines * sizeof(VAR_DATA));

#if VAR_DATA_PRECISION==1
  while(fscanf(file, "%f\t%f", &material->SP_Energy[i], &material->Stop_Pow[i]) > 0){
#else
  while(fscanf(file, "%lf\t%lf", &material->SP_Energy[i], &material->Stop_Pow[i]) > 0){
#endif

    material->SP_Energy[i] = material->SP_Energy[i] * UMeV;	// Energie cinétique en eV
    material->Stop_Pow[i] = material->Stop_Pow[i] * UMeV;	// Pouvoir d'arrêt en eV cm² / g
    i++;
  }

  fclose(file);
  return 0;
}

