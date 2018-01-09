/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_range_shifter.h"

int Init_RangeShifter_Data(plan_parameters *plan, machine_parameters *machine, Materials *material, DATA_config *config){

  int i, j, RS_beam_enabled = 0, RS_layer_enabled = 0;

  double SPR = 1.0;
  if(machine->RS_Type != none) SPR = machine->RS_Density * material[machine->RS_Material].Stop_Pow[(int)floor(100/PSTAR_BIN)] / material[config->Water_Material_ID].Stop_Pow[(int)floor(100/PSTAR_BIN)];

  for(i=0; i<plan->NumberOfFields; i++){
    if(plan->fields[i].RS_Type != none){
	// Check RS type
	if(plan->fields[i].RS_Type != machine->RS_Type){
	  printf("\nERROR: The range shifter type for beam %d does not match the BDL description\n", i);
	  return 1;
	}
	RS_beam_enabled = 1;
	for(j=0; j<plan->fields[i].NumberOfControlPoints; j++){
	  if(plan->fields[i].ControlPoints[j].RS_setting == IN && plan->fields[i].ControlPoints[j].RS_WET > 0){
	    RS_layer_enabled = 1;
	    plan->fields[i].ControlPoints[j].RS_Thickness = plan->fields[i].ControlPoints[j].RS_WET / SPR;
	  }
	  else{
	    plan->fields[i].ControlPoints[j].RS_setting == OUT;
	    plan->fields[i].ControlPoints[j].RS_WET == 0.0;
	    plan->fields[i].ControlPoints[j].RS_Thickness == 0.0;
	  }
	}
    }
    else{
      for(j=0; j<plan->fields[i].NumberOfControlPoints; j++){
	plan->fields[i].ControlPoints[j].RS_setting == OUT;
      }
    }
  }

  if(RS_beam_enabled == 1 && RS_layer_enabled == 1) config->RangeShifter_enabled = 1;

  return 0;
}


void Display_RangeShifter_Data(plan_parameters *plan, machine_parameters *machine, Materials *material){

  int i, j, AlwaysIN, AlwaysOUT, FixedPosition, FixedWET;
  double FirstPosition, FirstWET, FirstThickness;

  for(i=0; i<plan->NumberOfFields; i++){
    if(plan->fields[i].RS_Type != none){
      FixedPosition = 1;
      FixedWET = 1;
      AlwaysIN = 1;
      AlwaysOUT = 1;
      FirstPosition = -1;
      for(j=0; j<plan->fields[i].NumberOfControlPoints; j++){
	if(FirstPosition == -1 && plan->fields[i].ControlPoints[j].RS_setting == IN){
	  FirstPosition = plan->fields[i].ControlPoints[j].RS_IsocenterDist;
	  FirstWET = plan->fields[i].ControlPoints[j].RS_WET;
	  FirstThickness = plan->fields[i].ControlPoints[j].RS_Thickness;
	}
	if(plan->fields[i].ControlPoints[j].RS_setting == IN){
	  AlwaysOUT = 0;
	  if(FirstPosition != plan->fields[i].ControlPoints[j].RS_IsocenterDist) FixedPosition = 0;
	  if(FirstWET != plan->fields[i].ControlPoints[j].RS_WET) FixedWET = 0;
	}
	if(plan->fields[i].ControlPoints[j].RS_setting == OUT) AlwaysIN = 0;
      }
      printf("\nRange shifter initialized for beam %d:\n", i);
      if(plan->fields[i].RS_Type == binary) printf("\ttype=binary\n");
      else if(plan->fields[i].RS_Type == analog) printf("\ttype=analog\n");
      printf("\tMaterial: %s (ID %d)\n", material[machine->RS_Material].Name, machine->RS_Material);
      printf("\tDensity: %.2lf g/cm3\n", machine->RS_Density);
      if(AlwaysIN == 1) printf("\tEnabled for all layers\n");
      if(AlwaysOUT == 1) printf("\tDisabled for all layers\n");
      if(FixedPosition == 1) printf("\tStatic position: %.2lf cm from isocenter\n", FirstPosition);
      if(FixedWET == 1) printf("\tStatic thickness: %.2lf cm (WET=%.2lf)\n", FirstThickness, FirstWET);
    }
  }

}


void Simulate_RangeShifter(Hadron_buffer *hadron_list, ControlPoint_parameters **layer_data, field_parameters **field_data, int *Nbr_hadrons, DATA_config *config, machine_parameters *machine, Materials *material, VSLStreamStatePtr RNG_Stream){

  Hadron hadron;
  Init_particles(&hadron);

  int Nbr_HadronSimulated = 0;
  int i, j, count;

  ALIGNED_(64) VAR_COMPUTE RS_exit_position[VLENGTH];
  ALIGNED_(64) int Hadron_ID[VLENGTH];

  while(1){
    for(i=0; i<VLENGTH; i++){
      if(hadron.v_type[i] != Unknown && hadron.v_z[i] <= RS_exit_position[i]){
	Extract_particle(&hadron_list[Hadron_ID[i]], i, &hadron);
	hadron.v_type[i] = Unknown;
      }

      if(hadron.v_type[i] == Unknown){
	for(j=Nbr_HadronSimulated; j<*Nbr_hadrons; j++){
	  if(layer_data[j]->RS_setting == OUT){
	    Nbr_HadronSimulated += 1;
	    continue;
	  }
	  else{
	    Insert_particle(&hadron, i, &hadron_list[Nbr_HadronSimulated]);
	    RS_exit_position[i] = layer_data[j]->RS_IsocenterDist;
	    Hadron_ID[i] = j;
	    Nbr_HadronSimulated += 1;
	    break;
	  }
	} // for loop HadronSimulated
      } // if unknown
    } // for loop VLENGTH

    count = __sec_reduce_add(hadron.v_type[vALL]);
    if(count == 0){
      break;
    }

    SemiInfiniteSlab_step(&hadron, material, hadron_list, layer_data, field_data, Hadron_ID, Nbr_hadrons, RS_exit_position, config, machine, RNG_Stream);
  } // end while

}


