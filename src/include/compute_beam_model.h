/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_compute_beam_model
#define H_compute_beam_model

#include "define.h"
#include "struct.h"
#include "data_beam_model.h"
#include "compute_random.h"
#include "compute_math.h"
#include "compute_treatment_uncertainties.h"
#include "compute_range_shifter.h"


double ConvertMuToProtons(double weight, double energy);
void deviates (double U[2][2], double sigmas[2], double R[2], VSLStreamStatePtr RNG_Stream);
void diagonalize (double A[2][2],double B[2], double T[2][2]);
double compute_mEnergy (machine_parameters *mac, double energy);
double compute_sEnergy (machine_parameters *mac, double energy);
double compute_mX (machine_parameters *mac, double energy);
double compute_mY (machine_parameters *mac, double energy);
double compute_mTheta (machine_parameters *mac, double energy);
double compute_mPhi (machine_parameters *mac, double energy);
double compute_eXTheta (machine_parameters *mac, double energy);
double compute_eYPhi (machine_parameters *mac, double energy);
void rotateX (double angle, double vector[3]);
void rotateY (double angle, double vector[3]);
void rotateZ (double angle, double vector[3]);
void Sample_particle (Hadron_buffer *hadron, VAR_DATA CT_Length[3], machine_parameters *mac, ControlPoint_parameters *ControlPoint, spot_parameters *spot, VSLStreamStatePtr RNG_Stream);
void BEV_to_CT_frame(Hadron_buffer *hadron, machine_parameters *mac, field_parameters *field);
void Transport_to_CT(Hadron_buffer *hadron, VAR_DATA CT_Length[3]);
void Transport_to_RangeShifter(Hadron_buffer *hadron, ControlPoint_parameters **layer_data, int Nbr_hadrons);
void Generate_PBS_particle(Hadron_buffer *hadron, int *Nbr_hadrons, VAR_DATA CT_Length[3], plan_parameters *plan, machine_parameters *machine, VSLStreamStatePtr RNG_Stream, DATA_config *config, Materials *material);
plan_parameters* Init_single_spot_plan(plan_parameters *Plan);
void Select_spot(plan_parameters *Plan, plan_parameters *Beamlet, int FieldID, int ControlPointID, int SpotID);
plan_parameters* Select_beam(plan_parameters *Plan, int Beam);

#endif
