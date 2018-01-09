/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#ifndef H_define
#define H_define

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include <omp.h> 
#include <mkl_vsl.h>
#include <mkl.h>

#define VAR_DATA_PRECISION 1	// 1 =  float, 2 = double
#define VAR_SCORING_PRECISION 1	// 1 =  float, 2 = double
#define VAR_COMPUTE_PRECISION 1	// 1 =  float, 2 = double

#define INTERP_BIN 1.0 // MeV
#define PSTAR_BIN 0.5 // MeV

#define BeamPOSx 2.0		// Position du point de départ des protons (cm)
#define BeamPOSy 2.0
#define BeamPOSz 0
#define PEnergy 200		// Energie des protons (MeV)

#define FictitiousInteraction 1
#define RandomHinge 2
#define VoxelInterface 3
#define NoInterface 4
#define InterfaceCrossing VoxelInterface

#define EM_PSTAR 1
#define EM_FIPPEL 2
#define EM_Method EM_PSTAR

#define SP_PSTAR 1
#define SP_GEANT4 2
#define DB_STOP_POW SP_GEANT4

#define R_ELEC 2.8179403e-13	// Rayon classique de l'électron en cm
#define MC2_ELEC 0.5109989e6	// Energie de masse de l'électron en eV
#define MC2_PRO 938.272046e6	// Energie de masse du proton en eV
#define N_AVO 6.0221415e23 	// Nombre d'avogadro en 1/mol

// Définition des unités:
#define Umm 0.1	// 1 mm = 0.1 cm
#define Ucm 1	// 1 cm = 1 cm
#define Udm 10	// 1 dm = 10 cm
#define Um 100	// 1 m = 100 cm

#define UeV 1		// 1 eV = 1 eV
#define UkeV 1e3	// 1 keV = 1e3 eV
#define UMeV 1e6	// 1 MeV = 1e6 eV

#define Uamu 931.46e6	// 1 unité de masse atomique = 931.46 MeV

#define WATER_LABEL 17

#if VAR_DATA_PRECISION==1
  #define VAR_DATA float
#else
  #define VAR_DATA double
#endif

#if VAR_SCORING_PRECISION==1
  #define VAR_SCORING float
#else
  #define VAR_SCORING double
#endif

#if VAR_COMPUTE_PRECISION==1
  #define VAR_COMPUTE float
  #define VLENGTH 16
#else
  #define VAR_COMPUTE double
  #define VLENGTH 8
#endif

#define vALL 0:VLENGTH

// Cross platform compatibility
#if defined(_MSC_VER)
  #define ALIGNED_(n) __declspec(align(n))
  #define M_PI 3.14159265359
  #include <mathimf.h>
  #include <BaseTsd.h>
  typedef SSIZE_T ssize_t;
  #define strtok_r strtok_s
#else
  #define ALIGNED_(n) __attribute__((aligned(n)))
  #include <math.h>
#endif

#ifndef __INTEL_COMPILER
  #pragma GCC diagnostic ignored "-Wunused-value"
  #define __assume_aligned __builtin_assume_aligned
#endif

#endif
