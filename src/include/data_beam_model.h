/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_data_beam_model
#define H_data_beam_model

#include "define.h"
#include "struct.h"
#include "File_process.h"

#define poly_size 10

typedef struct spot_parameters spot_parameters;
struct spot_parameters
{
    double Spot_X;
    double Spot_Y;
    double Spot_Weight;
    double Spot_Time;
};

enum RangeShifter_setting{
	OUT,
	IN
} ;

typedef struct ControlPoint_parameters ControlPoint_parameters;
struct ControlPoint_parameters
{
    int ControlPointIndex;
    int SpotTunnedID;
    double CumulativeMetersetWeight;
    double Energy;
    int NbOfScannedSpots;
    spot_parameters *spots;
    VAR_DATA *Spots_cumulative_PDF;
    enum RangeShifter_setting RS_setting;
    double RS_IsocenterDist;
    double RS_WET;
    double RS_Thickness;
};

enum RangeShifter_type{
	none,
	binary,
	analog
} ;

typedef struct field_parameters field_parameters;
struct field_parameters
{
    int FieldID;
    double FinalCumulativeMeterSetWeight;
    double GantryAngle;
    double PatientSupportAngle;
    double IsocenterPositionX;
    double IsocenterPositionY;
    double IsocenterPositionZ;
    int NumberOfControlPoints;
    ControlPoint_parameters *ControlPoints;
    VAR_DATA *ControlPoints_cumulative_PDF;
    enum RangeShifter_type RS_Type;
};

typedef struct plan_parameters plan_parameters;
struct plan_parameters
{
    char PlanName[100];
    int NumberOfFractions;
    int FractionID;
    int NumberOfFields;
    int *FieldsID;
    double TotalMetersetWeightOfAllFields;
    VAR_DATA cumulative_weight;
    VAR_DATA normalization_factor;
    field_parameters *fields;
    VAR_DATA *Fields_cumulative_PDF;
};

enum Beam_Model_type{
	Grevillot,
	UPenn
} ;

typedef struct machine_parameters machine_parameters;
struct machine_parameters
{
    enum Beam_Model_type Beam_Model;

    double mDistanceSourcePatient;
    double mDistanceSMXToIsocenter;
    double mDistanceSMYToIsocenter;

    // Range shifter parameters
    int RS_defined;
    char RS_ID[100];
    enum RangeShifter_type RS_Type;
//    double RS_Position;
//    double RS_Thickness;
    double RS_Density;
    int RS_Material;
    double RS_WET;

    // Loic Grevillot beam model
    int mEnergy_order;
    double mEnergy_poly[poly_size];
    int sEnergy_order;
    double sEnergy_poly[poly_size];
    int mX_order;
    double mX_poly[poly_size];
    int mTheta_order;
    double mTheta_poly[poly_size];
    int mY_order;
    double mY_poly[poly_size];
    int mPhi_order;
    double mPhi_poly[poly_size];
    int eXTheta_order;
    double eXTheta_poly[poly_size];
    int eYPhi_order;
    double eYPhi_poly[poly_size];

    // UPenn beam model
    int Number_Energies;
    VAR_DATA *Nominal_Energies;
    VAR_DATA *Mean_Energies;
    VAR_DATA *Energy_Spread;
    VAR_DATA *Proton_Per_MU;
    VAR_DATA *Weight1;
    VAR_DATA *Weight2;

    VAR_DATA *SpotSize1x;
    VAR_DATA *Divergence1x;
    VAR_DATA *Correlation1x;

    VAR_DATA *SpotSize1y;
    VAR_DATA *Divergence1y;
    VAR_DATA *Correlation1y;

    VAR_DATA *SpotSize2x;
    VAR_DATA *Divergence2x;
    VAR_DATA *Correlation2x;

    VAR_DATA *SpotSize2y;
    VAR_DATA *Divergence2y;
    VAR_DATA *Correlation2y;
};



#include "compute_beam_model.h"


int read_machine_parameters(char* machine_name, machine_parameters *mac);
int read_UPenn_BDL(char* machine_name, machine_parameters *mac);
int read_Grevillot_BDL(char* machine_name, machine_parameters *mac);
void display_machine_parameters(machine_parameters *machine);
plan_parameters* read_plan_parameters(char* plan_name, DATA_config *config, machine_parameters *machine);
void display_plan_parameters(plan_parameters *plan);
void Free_Plan_Parameters(plan_parameters *plan);
void Free_Machine_Parameters(machine_parameters *mac);

#endif
