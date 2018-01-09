/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_struct
#define H_struct

#include "define.h"

/*
typedef struct DATA_Stop_Pow DATA_Stop_Pow;
struct DATA_Stop_Pow{
	VAR_DATA Energy;	// Energie cinétique
	VAR_DATA Stop_Pow;	// Pouvoir d'arrêt

};
*/

typedef struct DATA_Nuclear_Elastic DATA_Nuclear_Elastic;
struct DATA_Nuclear_Elastic{
//	VAR_DATA Energy;			// Energie
	VAR_DATA Cross_section;			// Section efficace totale à cette énergie
	VAR_DATA Diff_Cross_section[36];	// Sections efficaces différentielles pour chaque angle
};


typedef struct DATA_Nuclear_Inelastic DATA_Nuclear_Inelastic;
struct DATA_Nuclear_Inelastic{
//	VAR_DATA Energy;			// Energie
	VAR_DATA Cross_section;			// Section efficace totale à cette énergie
	VAR_DATA Energy_fraction_recoils;	// Fraction de l'énergie qui sera transmise au reste du noyau cible

	VAR_DATA Proton_Mult;			// Multiplicité des protons à cette énergie
	int P_Nbr_data;				// Nbr de lignes (énergies) pour la double section efficace differentielle
	VAR_DATA *P_Energy_list;		// liste des énergies définies dans les données
	VAR_DATA *P_D_Cross_section;		// Section efficace différentielle intégrée sur l'angle
	VAR_DATA *P_DD_Cross_section;		// Section efficace double différentielle

	VAR_DATA Alpha_Mult;			// Multiplicité des alphas à cette énergie
	int A_Nbr_data;				// Nbr de lignes (énergies) pour la double section efficace differentielle
	VAR_DATA *A_Energy_list;		// liste des énergies définies dans les données
	VAR_DATA *A_D_Cross_section;		// Section efficace différentielle intégrée sur l'angle
	VAR_DATA *A_DD_Cross_section;		// Section efficace double différentielle

	VAR_DATA Deuteron_Mult;			// Multiplicité des deuterons à cette énergie
	int D_Nbr_data;				// Nbr de lignes (énergies) pour la double section efficace differentielle
	VAR_DATA *D_Energy_list;		// liste des énergies définies dans les données
	VAR_DATA *D_D_Cross_section;		// Section efficace différentielle intégrée sur l'angle
	VAR_DATA *D_DD_Cross_section;		// Section efficace double différentielle

};


typedef struct DATA_PG DATA_PG;
struct DATA_PG{
//	VAR_DATA Energy;			// Energie
	VAR_DATA Multiplicity;			// Multiplicité à cette énergie
	int Nbr_data;				// Nbr de lignes (énergies) pour la section efficace differentielle
	VAR_DATA *Energy_list;			// liste des énergies définies dans les données
	VAR_DATA *Diff_Cross_section;		// Sections efficaces différentielles pour chaque angle
};


typedef struct Materials Materials;
struct Materials{
	char Name[50];		// Nom du matériau
	VAR_DATA A;		// Atomic or Molecular weight
	VAR_DATA Density;	// Densité (g/cm3)
	VAR_DATA X0;		// Longueur de radiation
	VAR_DATA N_el;		// Densité électronique
//	DATA_Stop_Pow *SP_list;// Pouvoir d'arret
	VAR_DATA *SP_Energy;	// energies corresponding to the SP table
	VAR_DATA *Stop_Pow;	// Stopping power table

	VAR_DATA SPR;

        // Nuclear data :
	enum { 
	  None,
	  ICRU,
	  Proton_Proton,
	  Mixture
	} Nuclear_data_type;

	int NbrComponents;
	int *Mixture_Components_label;
	VAR_DATA *Mixture_Components_fraction;
	VAR_DATA *Interp_Total_Nuclear_Cross_Section;

	// Nuclear Elastic
	int Nbr_Elastic_Energy;
	VAR_DATA *Elastic_Energy_List;
	DATA_Nuclear_Elastic *Nuclear_Elastic;

	// Nuclear Inelastic
	int Nbr_Inelastic_Energy;
	VAR_DATA *Inelastic_Energy_List;
	DATA_Nuclear_Inelastic *Nuclear_Inelastic;

	// Gamma Prompt
	int Nbr_PG_Energy;
	VAR_DATA *PG_Energy_List;
	DATA_PG *PG_data;
};


typedef struct DATA_CT DATA_CT;
struct DATA_CT{
	int GridSize[3];		// grid size (nbr voxel par axe)
	int Nbr_voxels;			// nbr total voxel
	VAR_DATA Length[3];		// image lengths (cm)
	VAR_DATA VoxelLength[3];	// voxel lengths (Lx/Nx, Ly/Ny, Lz/Nz)
	VAR_DATA Origin[3];		// voxel lengths (Lx/Nx, Ly/Ny, Lz/Nz)
	unsigned short int *material;	// material label array
	VAR_DATA *density;		// density array

	// HU conversion data:
	int Num_Density_Data;			// Number of data points for HU to density conversion
	int Num_Materials_Data;			// Number of data points for HU to material conversion
	VAR_DATA *Conversion_HU_Density;	// List of HU values for HU to density conversion
	VAR_DATA *Conversion_Densities;		// List of density values for HU to density conversion
	VAR_DATA *Conversion_HU_Material;	// List of HU values for HU to material conversion
	VAR_DATA *Conversion_Density_Material;	// List of density values for HU to material conversion
	unsigned short int *Conversion_Material_labels;	// List of material labels for HU to material conversion

	// Robustness test density scaling:
	VAR_DATA *Nominal_density;
	VAR_DATA *Scaled_density;
};


typedef struct DATA_4D_Fields DATA_4D_Fields;
struct DATA_4D_Fields{
	int Nbr_Fields;
	int GridSize[4];
	VAR_DATA Spacing[3];
	VAR_DATA Origin[3];
	VAR_DATA **Phase2Ref;
	VAR_DATA **Ref2Phase;
	VAR_DATA **Ref2Phase_log;
};


typedef struct DATA_Scoring DATA_Scoring;
struct DATA_Scoring{
	VAR_SCORING *energy;
	VAR_SCORING *dose;
	VAR_SCORING *PG_particles;
	VAR_SCORING *PG_spectrum;
	VAR_SCORING *LET;
	VAR_SCORING *LET_denominator;
};


enum Hadron_type{ 
	Unknown,
	Proton,
	Deuteron,
	Triton,
	Alpha
} ;

typedef struct Hadron Hadron;
struct Hadron{
	// Position
	ALIGNED_(64) VAR_COMPUTE v_x[VLENGTH];
	ALIGNED_(64) VAR_COMPUTE v_y[VLENGTH];
	ALIGNED_(64) VAR_COMPUTE v_z[VLENGTH];

	// Direction
	ALIGNED_(64) VAR_COMPUTE v_u[VLENGTH];
	ALIGNED_(64) VAR_COMPUTE v_v[VLENGTH];
	ALIGNED_(64) VAR_COMPUTE v_w[VLENGTH];

	ALIGNED_(64) VAR_COMPUTE v_T[VLENGTH];		// Energie cinétique (en eV)
	ALIGNED_(64) VAR_COMPUTE v_M[VLENGTH];		// Multiplicité (par défaut = 1)
	ALIGNED_(64) VAR_COMPUTE v_charge[VLENGTH];	// Charge de la particule
	ALIGNED_(64) VAR_COMPUTE v_mass[VLENGTH];	// Nbr de masses de l'ion

	ALIGNED_(64) enum Hadron_type v_type[VLENGTH];	// type de particule

	ALIGNED_(64) VAR_COMPUTE v_E[VLENGTH];		// Energie totale
	ALIGNED_(64) VAR_COMPUTE v_gamma[VLENGTH];	// Paramètre relativiste gamma
	ALIGNED_(64) VAR_COMPUTE v_beta2[VLENGTH];	// Betta au carré (v/c)^2
	ALIGNED_(64) VAR_COMPUTE v_Te_max[VLENGTH];	// Energie maximum transferable à l'e-
};


typedef struct Hadron_buffer Hadron_buffer;
struct Hadron_buffer{
	// Position
	VAR_COMPUTE x;
	VAR_COMPUTE y;
	VAR_COMPUTE z;

	// Direction
	VAR_COMPUTE u;
	VAR_COMPUTE v;
	VAR_COMPUTE w;

	VAR_COMPUTE T;		// Energie cinétique (en eV)
	VAR_COMPUTE M;		// Multiplicité (par défaut = 1)
	VAR_COMPUTE charge;	// Charge de la particule
	VAR_COMPUTE mass;	// Nbr de masses de l'ion

	enum Hadron_type type;	// type de particule

	VAR_COMPUTE E;		// Energie totale
	VAR_COMPUTE gamma;	// Paramètre relativiste gamma
	VAR_COMPUTE beta2;	// Betta au carré (v/c)^2
	VAR_COMPUTE Te_max;	// Energie maximum transferable à l'e-
};

enum Scenario_type{ 
	Regular,
	Nominal,
	Uncertainty
} ;

typedef struct DATA_config DATA_config;
struct DATA_config{

	unsigned int Num_Config_Tags;

	// Simulation parameters
	unsigned int Num_Threads;
	unsigned int RNG_Seed;
	unsigned long Num_Primaries;
	VAR_DATA Ecut_Pro;
	VAR_DATA D_Max;
	VAR_DATA Epsilon_Max;
	VAR_DATA Te_Min;

	// Input files
	char CT_File[200];
	char HU_Density_File[200];
	char HU_Material_File[200];
	char BDL_machine[200];
	char BDL_plan[200];

	// Physical parameters
	unsigned int Simulate_Nuclear_Interactions;
	unsigned int Simulate_Secondary_Protons;
	unsigned int Simulate_Secondary_Deuterons;
	unsigned int Simulate_Secondary_Alphas;

	// 4D simulation
	unsigned int Simu_4D_Mode;
	unsigned int Dose_4D_Accumulation;
	int Field_type;
	unsigned int Create_Ref_from_4DCT;
	unsigned int Create_4DCT_from_Ref;
	unsigned int Dynamic_delivery;
	VAR_DATA  Breathing_period;
	
	// Robustness simulation
	unsigned int Robustness_Mode;
	int Scenario_selection;
	unsigned int Simulate_nominal_plan;
	VAR_DATA  Systematic_Setup_Error[3];
	VAR_DATA  Random_Setup_Error[3];
	VAR_DATA  Systematic_Range_Error;
	VAR_DATA  Systematic_Amplitude_Error;
	VAR_DATA  Random_Amplitude_Error;
	VAR_DATA  Systematic_Period_Error;
	VAR_DATA  Random_Period_Error;

	// Beamlet simulation
	unsigned int Beamlet_Mode;
	unsigned int Beamlet_Parallelization;

	// Output parameters
	char Output_Directory[200];
	unsigned int Energy_ASCII_Output;
	unsigned int Energy_MHD_Output;
	unsigned int Energy_Sparse_Output;
	unsigned int Dose_ASCII_Output;
	unsigned int Dose_MHD_Output;
	unsigned int Dose_Sparse_Output;
	unsigned int LET_ASCII_Output;
	unsigned int LET_MHD_Output;
	unsigned int LET_Sparse_Output;
	unsigned int Densities_Output;
	unsigned int Materials_Output;
	unsigned int Compute_DVH;
	VAR_DATA Dose_Sparse_Threshold;
	VAR_DATA Energy_Sparse_Threshold;
	VAR_DATA LET_Sparse_Threshold;
	unsigned int Score_PromptGammas;
	VAR_DATA PG_LowEnergyCut;
	VAR_DATA PG_HighEnergyCut;
	unsigned int PG_Spectrum_NumBin;
	VAR_DATA PG_Spectrum_Binning;
	int LET_Calculation_Method;
	unsigned int Export_Beam_dose;
	int DoseToWater;
	unsigned int Dose_Segmentation;
	VAR_DATA Segmentation_Density_Threshold;

	// Internal variables
	unsigned int Particle_Generated_outside;
	VAR_DATA Num_delivered_protons;
	unsigned int Num_4DCT_phases;
	unsigned int TotalNbrSpots;
	int TotalNumScenarios;
	unsigned int Num_Materials;
	unsigned int Water_Material_ID;
	unsigned int Num_Components;
	char Materials_Dir[200];
	time_t timestamp;
	unsigned int RangeShifter_enabled;
	char output_robustness_suffix[100];
	char output_beamlet_suffix[100];
	char output_4D_suffix[100];
	char output_beams_suffix[100];
	VAR_DATA Current_Systematic_setup[3];
	VAR_DATA Current_Random_setup[3];
	VAR_DATA Current_Range_error;
	VAR_DATA Current_Systematic_amplitude;
	VAR_DATA Current_Random_amplitude;
	VAR_DATA Current_Breathing_amplitude;
	VAR_DATA Current_Systematic_period;
	VAR_DATA Current_Random_period;
	VAR_DATA Current_Breathing_period;
	VAR_DATA Current_init_delivery_points[10];
	int Current_4D_phase;
	int Current_Beam;
	int Current_scenario;
	int Current_fraction;
	int Fraction_accumulation;
	enum Scenario_type Current_scenario_type;
	VAR_DATA MCS_const;
	unsigned int Score_LET;

};

void Init_particles(Hadron *hadron);
void Insert_particle(Hadron *destination, int index, Hadron_buffer *origin);
void Extract_particle(Hadron_buffer *destination, int index, Hadron *origin);
void Copy_particle_buffer(Hadron_buffer *destination, Hadron_buffer *origin);
void Generate_particle(Hadron *destination, int index, VAR_COMPUTE beamx, VAR_COMPUTE beamy, VAR_COMPUTE beamz, VAR_COMPUTE Energy);
void Update_Hadron(Hadron *hadron);
void Copy_Hadron_struct(Hadron *destination, Hadron *origin);

#endif
