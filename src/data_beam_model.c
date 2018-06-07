/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/data_beam_model.h"


int read_machine_parameters(char* machine_name, machine_parameters *mac){

    char tmp [256];
    FILE *machine;
    machine=fopen(machine_name,"r");

    if (machine==NULL)
    {
	printf("unable to open machine parameters\n");
	return 1;
    }

    fgets(tmp,256,machine);
    fclose(machine);

    mac->RS_Type = none;
    mac->RS_Density = 0.0;
    mac->RS_Material = 17;

    int error;

    if(strcmp(tmp, "--UPenn beam model (double gaussian)--\n") == 0){
	//printf("\n Upenn beam model\n");
	mac->Beam_Model = UPenn;
	error = read_UPenn_BDL(machine_name, mac);
	if(error != 0) return 1;
    }
    else{
	mac->Beam_Model = Grevillot;
	error = read_Grevillot_BDL(machine_name, mac);
	if(error != 0) return 1;
    }

    return 0;

}


int read_UPenn_BDL(char* machine_name, machine_parameters *mac){

  char read[500], *read_token;

  FILE *file = fopen(machine_name,"r");
  if(file == NULL){
	printf("Error: Unable to open \"%s\".\n", machine_name);
	return 1;
  }

  int i;


  while (fgets(read, 500, file) != NULL){

    // remove comments
    if(read[0] == '#' || read[0] == '\n') continue;

    if(strcmp(read, "Nozzle exit to Isocenter distance\n") == 0){
	fgets(read,500,file);
	read_token = strtok(read, " \t\r\n");
	if(read_token == NULL || !isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Nozzle exit to Isocenter distance in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	}
	mac->mDistanceSourcePatient = atof(read_token);
    }
    else if(strcmp(read, "SMX to Isocenter distance\n") == 0){
	fgets(read,500,file);
	read_token = strtok(read, " \t\r\n");
	if(read_token == NULL || !isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SMX to Isocenter distance in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	}
	mac->mDistanceSMXToIsocenter = atof(read_token);
    }
    else if(strcmp(read, "SMY to Isocenter distance\n") == 0){
	fgets(read,500,file);
	read_token = strtok(read, " \t\r\n");
	if(read_token == NULL || !isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SMY to Isocenter distance in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	}
	mac->mDistanceSMYToIsocenter = atof(read_token);
    }
    else if(strcmp(read, "Range Shifter parameters\n") == 0){
	while(1){
	  fgets(read,500,file);
	  read_token = strtok(read, "= \t#\r\n");
	  if(read_token == NULL) break;
	  else if(strcmp(read_token, "RS_type") == 0){
	    read_token = strtok(NULL, "= \t#\r\n");
	    if(strcmp(read_token, "none") == 0) mac->RS_Type = none;
	    else if(strcmp(read_token, "binary") == 0) mac->RS_Type = binary;
	    else if(strcmp(read_token, "analog") == 0) mac->RS_Type = analog;
	    else{
		printf("\n Error: \"%s\" is not a valid value for RS_type in \"%s\"\n\n", read_token, machine_name);
		return 1;
	    }
	  }
/*
	  else if(strcmp(read_token, "RS_position") == 0){
	    read_token = strtok(NULL, "= \t#\r\n");
	    if(read_token == NULL || !isUnsignedFloat(read_token)){
		printf("\n Error: \"%s\" is not a valid value for RS_position in \"%s\"\n\n", read_token, machine_name);
		return 1;
	    }
	    mac->RS_Position = atof(read_token);
	  }
	  else if(strcmp(read_token, "RS_thickness") == 0){
	    read_token = strtok(NULL, "= \t#\r\n");
	    if(read_token == NULL || !isUnsignedFloat(read_token)){
		printf("\n Error: \"%s\" is not a valid value for RS_thickness in \"%s\"\n\n", read_token, machine_name);
		return 1;
	    }
	    mac->RS_Thickness = atof(read_token);
	  }
*/
	  else if(strcmp(read_token, "RS_material") == 0){
	    read_token = strtok(NULL, "= \t#\r\n");
	    if(read_token == NULL || !isUnsignedInt(read_token)){
		printf("\n Error: \"%s\" is not a valid value for RS_material in \"%s\"\n\n", read_token, machine_name);
		return 1;
	    }
	    mac->RS_Material = atoi(read_token);
	  }
	  else if(strcmp(read_token, "RS_density") == 0){
	    read_token = strtok(NULL, "= \t#\r\n");
	    if(read_token == NULL || !isUnsignedFloat(read_token)){
		printf("\n Error: \"%s\" is not a valid value for RS_density in \"%s\"\n\n", read_token, machine_name);
		return 1;
	    }
	    mac->RS_Density = atof(read_token);
	  }
	  else break;
	}
    }
    else if(strcmp(read, "Beam parameters\n") == 0){

	fgets(read,500,file);
	read_token = strtok(read, " \t\r\n");
	if(!isUnsignedInt(read_token) || atoi(read_token) == 0){
	    printf("\n Error: \"%s\" is not a valid value for Number of beam parameter energies in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	}
	mac->Number_Energies = atoi(read_token);

	mac->Nominal_Energies = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Mean_Energies = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Energy_Spread = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Proton_Per_MU = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Weight1 = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Weight2 = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));

	mac->SpotSize1x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Divergence1x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Correlation1x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));

	mac->SpotSize1y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Divergence1y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Correlation1y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));

	mac->SpotSize2x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Divergence2x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Correlation2x = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));

	mac->SpotSize2y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Divergence2y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));
	mac->Correlation2y = (VAR_DATA*)malloc(mac->Number_Energies * sizeof(VAR_DATA));

	fgets(read,500,file);
	fgets(read,500,file);

	for(i=0; i<mac->Number_Energies; i++){
	  fgets(read,500,file);

	  read_token = strtok(read, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for NominalEnergy in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Nominal_Energies[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for MeanEnergy in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Mean_Energies[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for EnergySpread in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Energy_Spread[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for ProtonsMU in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Proton_Per_MU[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Weight1 in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Weight1[i] = atof(read_token);


	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SpotSize1x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->SpotSize1x[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Divergence1x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Divergence1x[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Correlation1x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Correlation1x[i] = atof(read_token);


	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SpotSize1y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->SpotSize1y[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Divergence1y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Divergence1y[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Correlation1y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Correlation1y[i] = atof(read_token);


	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Weight2 in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Weight2[i] = atof(read_token);


	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SpotSize2x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->SpotSize2x[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Divergence2x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Divergence2x[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Correlation2x in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Correlation2x[i] = atof(read_token);


	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for SpotSize2y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->SpotSize2y[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Divergence2y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Divergence2y[i] = atof(read_token);

	  read_token = strtok(NULL, " \t\r\n");
	  if(!isFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for Correlation2y in \"%s\"\n\n", read_token, machine_name);
	    return 1;
	  }
	  mac->Correlation2y[i] = atof(read_token);
	}
    }

  }

  fclose(file);

  return 0;

}


int read_Grevillot_BDL(char* machine_name, machine_parameters *mac)
{


    char dummy [256];
    FILE *machine;
    machine=fopen(machine_name,"r");
    int i=0;

    if (machine==NULL)
    {
		printf("unable to open machine parameters\n");
		return 1;
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%lf",&mac->mDistanceSourcePatient);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%lf",&mac->mDistanceSMXToIsocenter);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%lf",&mac->mDistanceSMYToIsocenter);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->mEnergy_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->mEnergy_order;i++)
    {
        fscanf(machine,"%lf",&mac->mEnergy_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->sEnergy_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->sEnergy_order;i++)
    {
        fscanf(machine,"%lf",&mac->sEnergy_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->mX_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->mX_order;i++)
    {
        fscanf(machine,"%lf",&mac->mX_poly[i]);
    //    printf("%lf\n",mac->mX_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->mTheta_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->mTheta_order;i++)
    {
        fscanf(machine,"%lf",&mac->mTheta_poly[i]);

        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->mY_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->mY_order;i++)
    {
        fscanf(machine,"%lf",&mac->mY_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->mPhi_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->mPhi_order;i++)
    {
        fscanf(machine,"%lf",&mac->mPhi_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->eXTheta_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    for (i=0;i<=mac->eXTheta_order;i++)
    {
        fscanf(machine,"%lf",&mac->eXTheta_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
    fscanf(machine,"%d",&mac->eYPhi_order);
//    printf("%d\n",mac->eYPhi_order);

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);
//    printf("%s",dummy);
    for (i=0;i<=mac->eYPhi_order;i++)
    {
        fscanf(machine,"%lf",&mac->eYPhi_poly[i]);
//        printf("%lf\n",mac->eYPhi_poly[i]);
        fgets(dummy,256,machine);
    }

    fgets(dummy,256,machine);
    fgets(dummy,256,machine);

    fclose(machine);

    return 0;
}


void display_machine_parameters (machine_parameters *machine){

  int i;

  printf("\n\n");
  printf("mDistanceSourcePatient: %lf \n", machine->mDistanceSourcePatient);
  printf("mDistanceSMXToIsocenter: %lf \n", machine->mDistanceSMXToIsocenter);
  printf("mDistanceSMYToIsocenter: %lf \n\n", machine->mDistanceSMYToIsocenter);

  printf("Range Shifter:\n");
  printf("RS_type: %d \n", machine->RS_Type);
//  printf("RS_position: %lf \n", machine->RS_Position);
//  printf("RS_thickness: %lf \n", machine->RS_Thickness);
  printf("RS_material: %d \n", machine->RS_Material);
  printf("RS_density: %lf \n\n", machine->RS_Density);

  if(machine->Beam_Model == UPenn){
    printf("Beam model type: UPenn\n\n");

    printf("Nominal Energies\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Nominal_Energies[i]);

    printf("\n\nMean Energies\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Mean_Energies[i]);

    printf("\n\nEnergy Spread\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Energy_Spread[i]);

    printf("\n\nProtons per MU\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %e ", machine->Proton_Per_MU[i]);

    printf("\n\nWeight 1\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Weight1[i]);


    printf("\n\nSpotSize 1 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->SpotSize1x[i]);

    printf("\n\nDivergence 1 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Divergence1x[i]);

    printf("\n\nCorrelation 1 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Correlation1x[i]);


    printf("\n\nSpotSize 1 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->SpotSize1y[i]);

    printf("\n\nDivergence 1 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Divergence1y[i]);

    printf("\n\nCorrelation 1 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Correlation1y[i]);


    printf("\n\nWeight 2\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Weight2[i]);


    printf("\n\nSpotSize 2 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->SpotSize2x[i]);

    printf("\n\nDivergence 2 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Divergence2x[i]);

    printf("\n\nCorrelation 2 X\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Correlation2x[i]);


    printf("\n\nSpotSize 2 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->SpotSize2y[i]);

    printf("\n\nDivergence 2 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Divergence2y[i]);

    printf("\n\nCorrelation 2 Y\n");
    for(i=0; i<machine->Number_Energies; i++) printf(" %f ", machine->Correlation2y[i]);
  }

  else if(machine->Beam_Model == Grevillot){
    printf("Beam model type: Grevillot\n\n");

    //printf("mEnergy_order: %d \n", machine->mEnergy_order);
    printf("mEnergy_poly:");
    for(i=0; i<machine->mEnergy_order; i++) printf(" %e E^%d +", machine->mEnergy_poly[i], i);
    printf(" %e E^%d \n", machine->mEnergy_poly[i], i);

    //printf("sEnergy_order: %d \n", machine->sEnergy_order);
    printf("sEnergy_poly:");
    for(i=0; i<machine->sEnergy_order; i++) printf(" %e E^%d +", machine->sEnergy_poly[i], i);
    printf(" %e E^%d \n", machine->sEnergy_poly[i], i);

    //printf("mX_order: %d \n", machine->mX_order);
    printf("mX_poly:");
    for(i=0; i<machine->mX_order; i++) printf(" %e E^%d +", machine->mX_poly[i], i);
    printf(" %e E^%d \n", machine->mX_poly[i], i);

    //printf("mTheta_order: %d \n", machine->mTheta_order);
    printf("mTheta_poly:");
    for(i=0; i<machine->mTheta_order; i++) printf(" %e E^%d +", machine->mTheta_poly[i], i);
    printf(" %e E^%d \n", machine->mTheta_poly[i], i);

    //printf("mY_order: %d \n", machine->mY_order);
    printf("mY_poly:");
    for(i=0; i<machine->mY_order; i++) printf(" %e E^%d +", machine->mY_poly[i], i);
    printf(" %e E^%d \n", machine->mY_poly[i], i);

    //printf("mPhi_order: %d \n", machine->mPhi_order);
    printf("mPhi_poly:");
    for(i=0; i<machine->mPhi_order; i++) printf(" %e E^%d +", machine->mPhi_poly[i], i);
    printf(" %e E^%d \n", machine->mPhi_poly[i], i);

    //printf("eXTheta_order: %d \n", machine->eXTheta_order);
    printf("eXTheta_poly:");
    for(i=0; i<machine->eXTheta_order; i++) printf(" %e E^%d +", machine->eXTheta_poly[i], i);
    printf(" %e E^%d \n", machine->eXTheta_poly[i], i);

    //printf("eYPhi_order: %d \n", machine->eYPhi_order);
    printf("eYPhi_poly:");
    for(i=0; i<machine->eYPhi_order; i++) printf(" %e E^%d +", machine->eYPhi_poly[i], i);
    printf(" %e E^%d \n", machine->eYPhi_poly[i], i);

  }

  printf("\n\n");

}


plan_parameters* read_plan_parameters(char* plan_name, DATA_config *config, machine_parameters *machine)
 {

    char read[500], *read_token;
    int i, j, k, l;
    double cumulative_weight = 0;

    // variables for Proton_Per_MU interpolation
    int EnergyID;
    VAR_COMPUTE Energy1, Energy2;

    FILE *plan_file;
    plan_file=fopen(plan_name,"r");

    if (plan_file==NULL)
    {
		printf("unable to open plan parameters\n");
		return NULL;
    }

    config->TotalNbrSpots = 0;

    plan_parameters *plan = (plan_parameters*)malloc(sizeof(plan_parameters));
    plan->fields = NULL;
    plan->Fields_cumulative_PDF = NULL;
    plan->FieldsID = NULL;

    // init plan
    strcpy(plan->PlanName, "");
    plan->NumberOfFractions = 0;
    plan->FractionID = 0;
    plan->NumberOfFields = 0;
    plan->TotalMetersetWeightOfAllFields = 0;

    while (fgets(read, 500, plan_file) != NULL){
      read_token = strtok(read, " \t\r\n");
      if(read_token == NULL) continue;

      if(strcmp(read_token, "#PlanName") == 0){
	fgets(read, 500, plan_file);
	read_token = strtok(read, "\t\r\n");
	if(read_token != NULL) strcpy(plan->PlanName, read_token);
      }
      else if(strcmp(read_token, "#NumberOfFractions") == 0){
	fgets(read, 500, plan_file);
	read_token = strtok(read, " \t\r\n");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for NumberOfFractions in \"%s\"\n\n", read_token, plan_name);
	}
	else plan->NumberOfFractions = atoi(read_token);
      }
      else if(strcmp(read_token, "##FractionID") == 0){
	fgets(read, 500, plan_file);
	read_token = strtok(read, " \t\r\n");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for FractionID in \"%s\"\n\n", read_token, plan_name);
	}
	else plan->FractionID = atoi(read_token);
      }
      else if(strcmp(read_token, "##NumberOfFields") == 0){
	fgets(read, 500, plan_file);
	read_token = strtok(read, " \t\r\n");
	if(!isUnsignedInt(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for NumberOfFields in \"%s\"\n\n", read_token, plan_name);
	}
	else{
	  plan->NumberOfFields = atoi(read_token);
          plan->fields = (field_parameters*)malloc(plan->NumberOfFields * sizeof(field_parameters));
          plan->Fields_cumulative_PDF = (VAR_DATA*)malloc(plan->NumberOfFields * sizeof(VAR_DATA));
          plan->FieldsID = (int*)malloc(plan->NumberOfFields * sizeof(int));
	  l = 0;
	  for (i=0;i<plan->NumberOfFields;i++){
            fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(read_token == NULL || strcmp(read_token, "###FieldsID") != 0){
	      printf("\n Warning: FieldsID list is interrupted in \"%s\"\n\n", plan_name);
	      break;
	    }
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isUnsignedInt(read_token)){
	      printf("\n Warning: \"%s\" is not a valid value for FieldsID in \"%s\"\n\n", read_token, plan_name);
	      break;
	    }
	    else plan->FieldsID[i] = atoi(read_token);
	  }
        }
      }
      else if(strcmp(read_token, "#TotalMetersetWeightOfAllFields") == 0){
	fgets(read, 500, plan_file);
	read_token = strtok(read, " \t\r\n");
	if(!isUnsignedFloat(read_token)){
	    printf("\n Error: \"%s\" is not a valid value for TotalMetersetWeightOfAllFields in \"%s\"\n\n", read_token, plan_name);
	}
	else plan->TotalMetersetWeightOfAllFields = atof(read_token);
      }
      else if(strcmp(read_token, "#FIELD-DESCRIPTION") == 0){
	plan->fields[l].FieldID = 0;
	plan->fields[l].FinalCumulativeMeterSetWeight = 0.0;
	plan->fields[l].GantryAngle = 0.0;
	plan->fields[l].PatientSupportAngle = 0.0;
	plan->fields[l].IsocenterPositionX = 0.0;
	plan->fields[l].IsocenterPositionY = 0.0;
	plan->fields[l].IsocenterPositionZ = 0.0;
	plan->fields[l].NumberOfControlPoints = 0;
	plan->fields[l].RS_Type = none;

	while(fgets(read, 500, plan_file) != NULL){
	  read_token = strtok(read, " \t\r\n");
	  if(read_token == NULL) continue;

	  if(strcmp(read_token, "#FIELD-DESCRIPTION") == 0){

	    l += 1;
	    plan->fields[l].FieldID = 0;
	    plan->fields[l].FinalCumulativeMeterSetWeight = 0.0;
	    plan->fields[l].GantryAngle = 0.0;
	    plan->fields[l].PatientSupportAngle = 0.0;
	    plan->fields[l].IsocenterPositionX = 0.0;
	    plan->fields[l].IsocenterPositionY = 0.0;
	    plan->fields[l].IsocenterPositionZ = 0.0;
	    plan->fields[l].NumberOfControlPoints = 0;
	    plan->fields[l].RS_Type = none;
	  }

	  if(strcmp(read_token, "###FieldID") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isUnsignedInt(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for FieldID in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].FieldID = atoi(read_token);
	  }
	  else if(strcmp(read_token, "###FinalCumulativeMeterSetWeight") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isUnsignedFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for FinalCumulativeMeterSetWeight in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].FinalCumulativeMeterSetWeight = atof(read_token);
	  }
	  else if(strcmp(read_token, "###GantryAngle") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for GantryAngle in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].GantryAngle = atof(read_token) * M_PI / 180;
	  }
	  else if(strcmp(read_token, "###PatientSupportAngle") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for PatientSupportAngle in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].PatientSupportAngle = atof(read_token) * M_PI / 180;
	  }
	  else if(strcmp(read_token, "###IsocenterPosition") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for IsocenterPosition in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].IsocenterPositionX = atof(read_token);
	    read_token = strtok(NULL, " \t\r\n");
	    if(!isFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for IsocenterPosition in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].IsocenterPositionY = atof(read_token);
	    read_token = strtok(NULL, " \t\r\n");
	    if(!isFloat(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for IsocenterPosition in \"%s\"\n\n", read_token, plan_name);
	    }
	    else plan->fields[l].IsocenterPositionZ = atof(read_token);
	  }
	  else if(strcmp(read_token, "###RangeShifterType") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(strcmp(read_token, "none") == 0) plan->fields[l].RS_Type = none;
	    else if(strcmp(read_token, "binary") == 0) plan->fields[l].RS_Type = binary;
	    else if(strcmp(read_token, "analog") == 0) plan->fields[l].RS_Type = analog;
	    else{
	      printf("\n Error: \"%s\" is not a valid value for RangeShifterType in \"%s\"\n\n", read_token, plan_name);
	    }
	  }
	  else if(strcmp(read_token, "###NumberOfControlPoints") == 0){
	    fgets(read, 500, plan_file);
	    read_token = strtok(read, " \t\r\n");
	    if(!isUnsignedInt(read_token)){
	      printf("\n Error: \"%s\" is not a valid value for NumberOfControlPoints in \"%s\"\n\n", read_token, plan_name);
	    }
	    else{
	      plan->fields[l].NumberOfControlPoints = atoi(read_token);
	      plan->fields[l].ControlPoints = (ControlPoint_parameters*)malloc(plan->fields[l].NumberOfControlPoints * sizeof(ControlPoint_parameters));
	      plan->fields[l].ControlPoints_cumulative_PDF = (VAR_DATA*)malloc(plan->fields[l].NumberOfControlPoints * sizeof(VAR_DATA));
	      j = 0;
	    }
	  }
	  else if(strcmp(read_token, "#SPOTS-DESCRIPTION") == 0){

	    while(j < plan->fields[l].NumberOfControlPoints){
		fgets(read, 500, plan_file);
	    	read_token = strtok(read, " \t\r\n");
      		if(read_token == NULL) continue;

		if(strcmp(read_token, "####ControlPointIndex") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedInt(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for ControlPointIndex in \"%s\"\n\n", read_token, plan_name);
	    	  }
		  else{
	    	    plan->fields[l].ControlPoints[j].ControlPointIndex = atoi(read_token);
	    	    plan->fields[l].ControlPoints[j].SpotTunnedID = 0;
	    	    plan->fields[l].ControlPoints[j].CumulativeMetersetWeight = 0.0;
	    	    plan->fields[l].ControlPoints[j].Energy = 0.0;
		    plan->fields[l].ControlPoints[j].NbOfScannedSpots = 0;
		    plan->fields[l].ControlPoints[j].RS_setting = OUT;
		    plan->fields[l].ControlPoints[j].RS_IsocenterDist = 40.0;
		    plan->fields[l].ControlPoints[j].RS_WET = 0.0;
		  }
	  	}
		else if(strcmp(read_token, "####SpotTunnedID") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedInt(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for SpotTunnedID in \"%s\"\n\n", read_token, plan_name);
	    	  }
	    	  else plan->fields[l].ControlPoints[j].SpotTunnedID = atoi(read_token);
	  	}
		else if(strcmp(read_token, "####CumulativeMetersetWeight") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedFloat(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for CumulativeMetersetWeight in \"%s\"\n\n", read_token, plan_name);
	    	  }
	    	  else plan->fields[l].ControlPoints[j].CumulativeMetersetWeight = atof(read_token);
	  	}
		else if(strcmp(read_token, "####Energy") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedFloat(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for Energy in \"%s\"\n\n", read_token, plan_name);
	    	  }
	    	  else plan->fields[l].ControlPoints[j].Energy = atof(read_token);
	  	}
		else if(strcmp(read_token, "####RangeShifterSetting") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(strcmp(read_token, "OUT") == 0) plan->fields[l].ControlPoints[j].RS_setting = OUT;
	    	  else if(strcmp(read_token, "IN") == 0) plan->fields[l].ControlPoints[j].RS_setting = IN;
		  else{
	      	    printf("\n Error: \"%s\" is not a valid value for RangeShifterSetting in \"%s\"\n\n", read_token, plan_name);
	    	  }
	  	}
		else if(strcmp(read_token, "####IsocenterToRangeShifterDistance") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedFloat(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for IsocenterToRangeShifterDistance in \"%s\"\n\n", read_token, plan_name);
	    	  }
	    	  else plan->fields[l].ControlPoints[j].RS_IsocenterDist = atof(read_token) / 10.0;
	  	}
		else if(strcmp(read_token, "####RangeShifterWaterEquivalentThickness") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedFloat(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for RangeShifterWaterEquivalentThickness in \"%s\"\n\n", read_token, plan_name);
	    	  }
	    	  else plan->fields[l].ControlPoints[j].RS_WET = atof(read_token) / 10.0;
	  	}
		else if(strcmp(read_token, "####NbOfScannedSpots") == 0){
	    	  fgets(read, 500, plan_file);
	    	  read_token = strtok(read, " \t\r\n");
	    	  if(!isUnsignedInt(read_token)){
	      	    printf("\n Error: \"%s\" is not a valid value for NbOfScannedSpots in \"%s\"\n\n", read_token, plan_name);
	    	  }
		  else{
	    	    plan->fields[l].ControlPoints[j].NbOfScannedSpots = atoi(read_token);
		    plan->fields[l].ControlPoints[j].spots = (spot_parameters*)malloc(plan->fields[l].ControlPoints[j].NbOfScannedSpots * sizeof(spot_parameters));
	    	    plan->fields[l].ControlPoints[j].Spots_cumulative_PDF = (VAR_DATA*)malloc(plan->fields[l].ControlPoints[j].NbOfScannedSpots * sizeof(VAR_DATA));
		  }
	  	}
		else if(strcmp(read_token, "####X") == 0){
		  for (k=0; k<plan->fields[l].ControlPoints[j].NbOfScannedSpots; k++){
			fgets(read, 500, plan_file);
			read_token = strtok(read, " \t\r\n");
			if(!isFloat(read_token)){
	      	    	  printf("\n Error: \"%s\" is not a valid value for Spot_X in \"%s\"\n\n", read_token, plan_name);
	    	  	}
			else plan->fields[l].ControlPoints[j].spots[k].Spot_X = atof(read_token);
			read_token = strtok(NULL, " \t\r\n");
			if(!isFloat(read_token)){
	      	    	  printf("\n Error: \"%s\" is not a valid value for Spot_Y in \"%s\"\n\n", read_token, plan_name);
	    	  	}
			else plan->fields[l].ControlPoints[j].spots[k].Spot_Y = atof(read_token);
			read_token = strtok(NULL, " \t\r\n");
			if(!isUnsignedFloat(read_token)){
	      	    	  printf("\n Error: \"%s\" is not a valid value for Spot_Weight in \"%s\"\n\n", read_token, plan_name);
	    	  	}
			else{
			  plan->fields[l].ControlPoints[j].spots[k].Spot_Weight = atof(read_token);
        // if beam is UPenn, Proton_Per_MU must be interpolated from beam model parameters
        if(machine->Beam_Model == UPenn){
          EnergyID = floor(plan->fields[l].ControlPoints[j].Energy / Machine_Param_BIN) - floor(machine->Nominal_Energies[0] / Machine_Param_BIN);
          Energy1 = machine->Nominal_Energies[EnergyID];
          Energy2 = Energy1 + Machine_Param_BIN;
          plan->fields[l].ControlPoints[j].spots[k].Spot_Weight = Linear_Interpolation(plan->fields[l].ControlPoints[j].Energy, Energy1, Energy2, machine->Proton_Per_MU[EnergyID], machine->Proton_Per_MU[EnergyID+1]);
          }
        else{
          plan->fields[l].ControlPoints[j].spots[k].Spot_Weight = ConvertMuToProtons(plan->fields[l].ControlPoints[j].spots[k].Spot_Weight,
          plan->fields[l].ControlPoints[j].Energy);
        }
			  cumulative_weight += plan->fields[l].ControlPoints[j].spots[k].Spot_Weight;
			  plan->fields[l].ControlPoints[j].Spots_cumulative_PDF[k] = cumulative_weight;
			  config->TotalNbrSpots++;
			}
			read_token = strtok(NULL, " \t\r\n");
			if(read_token != NULL){
			  if(!isFloat(read_token)){
	      	    	    printf("\n Error: \"%s\" is not a valid value for Spot_Time in \"%s\"\n\n", read_token, plan_name);
			    plan->fields[l].ControlPoints[j].spots[k].Spot_Time = 0;
	    	  	  }
			  else plan->fields[l].ControlPoints[j].spots[k].Spot_Time = atof(read_token);
			}
			else if(config->Dynamic_delivery == 1){
			  printf("\n Error: Spot_Time is expected in \"%s\" for dynamic delivery simulation\n\n", plan_name);
			  return NULL;
			}
			else plan->fields[l].ControlPoints[j].spots[k].Spot_Time = 0;
		  }
		  plan->fields[l].ControlPoints_cumulative_PDF[j] = cumulative_weight;
		  j += 1;
		}// end spot parameters
	    }// end layer loop
	    plan->Fields_cumulative_PDF[l] = cumulative_weight;
	  }// end SPOT DESCRIPTION
	}// end loop field parameters
      }// end FIELD DESCRIPTION
    }// end plan loop

    plan->normalization_factor = 1.0;
    plan->cumulative_weight = cumulative_weight;
    config->Num_delivered_protons = cumulative_weight;

    fclose(plan_file);

    return plan;
}


void display_plan_parameters(plan_parameters *plan){

  int i, j, k;

  printf("\n\n");
  printf("PlanName: %s \n", plan->PlanName);
  printf("NumberOfFractions: %d \n", plan->NumberOfFractions);
  //printf("FractionID: %d \n", plan->FractionID);
  printf("NumberOfFields: %d \n", plan->NumberOfFields);
  //printf("FieldsID: %d \n", plan->FieldsID[0]);
  //printf("TotalMetersetWeightOfAllFields: %d \n", plan->TotalMetersetWeightOfAllFields);
  printf("CumulativeWeight: %.4e \n", plan->cumulative_weight);

  for(i=0; i<plan->NumberOfFields; i++){
    printf("\n");
    printf("\tFieldID: %d \n", plan->fields[i].FieldID);
    //printf("\tFinalCumulativeMeterSetWeight: %d \n", plan->fields[0].FinalCumulativeMeterSetWeight);
    printf("\tCumulativeWeight: %.4e \n", plan->Fields_cumulative_PDF[i]);
    printf("\tGantryAngle: %.2f \n", plan->fields[i].GantryAngle);
    printf("\tPatientSupportAngle: %.2f \n", plan->fields[i].PatientSupportAngle);
    printf("\tIsocenterPosition: (%.2f ; %.2f ; %.2f) \n", plan->fields[i].IsocenterPositionX, plan->fields[i].IsocenterPositionY, plan->fields[i].IsocenterPositionZ);
    printf("\tRangeShifter Type: %d \n", plan->fields[i].RS_Type);
    printf("\tNumberOfControlPoints: %d \n", plan->fields[i].NumberOfControlPoints);

    for(j=0; j<plan->fields[i].NumberOfControlPoints; j++){
      printf("\n");
      printf("\t\tControlPointIndex: %d \n", plan->fields[i].ControlPoints[j].ControlPointIndex);
      printf("\t\tSpotTunnedID: %d \n", plan->fields[i].ControlPoints[j].SpotTunnedID);
      //printf("\t\tCumulativeMetersetWeight: %d \n", plan->fields[i].ControlPoints[j].CumulativeMetersetWeight);
      printf("\t\tCumulativeWeight: %.4e \n", plan->fields[i].ControlPoints_cumulative_PDF[j]);
      printf("\t\tEnergy: %.2f \n", plan->fields[i].ControlPoints[j].Energy);
      printf("\t\tRangeShifter Setting: %d \n", plan->fields[i].ControlPoints[j].RS_setting);
      printf("\t\tRangeShifter Position: %lf \n", plan->fields[i].ControlPoints[j].RS_IsocenterDist);
      printf("\t\tRangeShifter WET: %lf \n", plan->fields[i].ControlPoints[j].RS_WET);
      printf("\t\tNbOfScannedSpots: %d \n", plan->fields[i].ControlPoints[j].NbOfScannedSpots);
      printf("\n");

      for(k=0; k<plan->fields[i].ControlPoints[j].NbOfScannedSpots; k++){
	if(plan->fields[i].ControlPoints[j].spots[k].Spot_Time == 0) printf("\t\t\tSpot %d: (%.2f ; %.2f) -> %.4e - %.4e\n", k, plan->fields[i].ControlPoints[j].spots[k].Spot_X, plan->fields[i].ControlPoints[j].spots[k].Spot_Y, plan->fields[i].ControlPoints[j].spots[k].Spot_Weight, plan->fields[i].ControlPoints[j].Spots_cumulative_PDF[k]);
	else printf("\t\t\tSpot %d: (%.2f ; %.2f) at %.2f ms -> %.4e - %.4e\n", k, plan->fields[i].ControlPoints[j].spots[k].Spot_X, plan->fields[i].ControlPoints[j].spots[k].Spot_Y, plan->fields[i].ControlPoints[j].spots[k].Spot_Time, plan->fields[i].ControlPoints[j].spots[k].Spot_Weight, plan->fields[i].ControlPoints[j].Spots_cumulative_PDF[k]);
      }
    }
  }

  printf("\n\n");
}




void Free_Plan_Parameters(plan_parameters *plan){

  if (plan == NULL) return;

  int i, j;
  for(i=0; i < plan->NumberOfFields; i++){

    for(j=0; j < plan->fields[i].NumberOfControlPoints; j++){

      if(plan->fields[i].ControlPoints[j].spots != NULL) free(plan->fields[i].ControlPoints[j].spots);
      if(plan->fields[i].ControlPoints[j].Spots_cumulative_PDF != NULL) free(plan->fields[i].ControlPoints[j].Spots_cumulative_PDF);
    }
    if(plan->fields[i].ControlPoints != NULL) free(plan->fields[i].ControlPoints);
    if(plan->fields[i].ControlPoints_cumulative_PDF != NULL) free(plan->fields[i].ControlPoints_cumulative_PDF);
  }

  if(plan->fields != NULL) free(plan->fields);

  if(plan->Fields_cumulative_PDF != NULL) free(plan->Fields_cumulative_PDF);

  if(plan->FieldsID != NULL) free(plan->FieldsID);

  free(plan);

}


void Free_Machine_Parameters(machine_parameters *mac){

  if(mac->Beam_Model == UPenn){
	if(mac->Nominal_Energies != NULL) free(mac->Nominal_Energies);
	if(mac->Mean_Energies != NULL) free(mac->Mean_Energies);
	if(mac->Energy_Spread != NULL) free(mac->Energy_Spread);
	if(mac->Proton_Per_MU != NULL) free(mac->Proton_Per_MU);
	if(mac->Weight1 != NULL) free(mac->Weight1);
	if(mac->Weight2 != NULL) free(mac->Weight2);

	if(mac->SpotSize1x != NULL) free(mac->SpotSize1x);
	if(mac->Divergence1x != NULL) free(mac->Divergence1x);
	if(mac->Correlation1x != NULL) free(mac->Correlation1x);

	if(mac->SpotSize1y != NULL) free(mac->SpotSize1y);
	if(mac->Divergence1y != NULL) free(mac->Divergence1y);
	if(mac->Correlation1y != NULL) free(mac->Correlation1y);

	if(mac->SpotSize2x != NULL) free(mac->SpotSize2x);
	if(mac->Divergence2x != NULL) free(mac->Divergence2x);
	if(mac->Correlation2x != NULL) free(mac->Correlation2x);

	if(mac->SpotSize2y != NULL) free(mac->SpotSize2y);
	if(mac->Divergence2y != NULL) free(mac->Divergence2y);
	if(mac->Correlation2y != NULL) free(mac->Correlation2y);
  }

}
