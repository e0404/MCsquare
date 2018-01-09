/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/struct.h"


void Init_particles(Hadron *hadron){

  hadron->v_x[vALL] = 0.0;
  hadron->v_y[vALL] = 0.0;
  hadron->v_z[vALL] = 0.0;

  hadron->v_u[vALL] = 0.0;
  hadron->v_v[vALL] = 0.0;
  hadron->v_w[vALL] = 0.1;

  hadron->v_T[vALL] = 0.0;

  hadron->v_M[vALL] = 1.0;
  hadron->v_type[vALL] = Unknown;
  hadron->v_charge[vALL] = 1.0;
  hadron->v_mass[vALL] = 1.0;

  return;
}


void Insert_particle(Hadron *destination, int index, Hadron_buffer *origin){

  destination->v_x[index] = origin->x;
  destination->v_y[index] = origin->y;
  destination->v_z[index] = origin->z;

  destination->v_u[index] = origin->u;
  destination->v_v[index] = origin->v;
  destination->v_w[index] = origin->w;

  destination->v_T[index] = origin->T;

  destination->v_M[index] = origin->M;
  destination->v_type[index] = origin->type;
  destination->v_charge[index] = origin->charge;
  destination->v_mass[index] = origin->mass;

  return;
}


void Extract_particle(Hadron_buffer *destination, int index, Hadron *origin){

  destination->x = origin->v_x[index];
  destination->y = origin->v_y[index];
  destination->z = origin->v_z[index];

  destination->u = origin->v_u[index];
  destination->v = origin->v_v[index];
  destination->w = origin->v_w[index];

  destination->T = origin->v_T[index];

  destination->M = origin->v_M[index];
  destination->type = origin->v_type[index];
  destination->charge = origin->v_charge[index];
  destination->mass = origin->v_mass[index];

  return;
}

void Copy_particle_buffer(Hadron_buffer *destination, Hadron_buffer *origin){

  destination->x = origin->x;
  destination->y = origin->y;
  destination->z = origin->z;

  destination->u = origin->u;
  destination->v = origin->v;
  destination->w = origin->w;

  destination->T = origin->T;

  destination->M = origin->M;
  destination->type = origin->type;
  destination->charge = origin->charge;
  destination->mass = origin->mass;

  return;
}


void Generate_particle(Hadron *destination, int index, VAR_COMPUTE beamx, VAR_COMPUTE beamy, VAR_COMPUTE beamz, VAR_COMPUTE Energy){

  destination->v_x[index] = beamx;
  destination->v_y[index] = beamy;
  destination->v_z[index] = beamz;

  destination->v_u[index] = 0.0;
  destination->v_v[index] = 0.0;
  destination->v_w[index] = 1.0;

  destination->v_T[index] = Energy;

  destination->v_M[index] = 1.0;
  destination->v_charge[index] = 1.0;
  destination->v_mass[index] = 1.0;
  destination->v_type[index] = Proton;
  
  return;
}


void Update_Hadron(Hadron *hadron){

  __assume_aligned(&hadron->v_T, 64);
  __assume_aligned(&hadron->v_M, 64);
  __assume_aligned(&hadron->v_charge, 64);
  __assume_aligned(&hadron->v_mass, 64);

  __assume_aligned(&hadron->v_E, 64);
  __assume_aligned(&hadron->v_gamma, 64);
  __assume_aligned(&hadron->v_beta2, 64);
  __assume_aligned(&hadron->v_Te_max, 64);


  hadron->v_E[vALL] = hadron->v_T[vALL] + hadron->v_mass[vALL] * MC2_PRO;	// Energie totale du proton
  hadron->v_gamma[vALL] = hadron->v_E[vALL] / (hadron->v_mass[vALL] * MC2_PRO);	// Paramètre relativiste gamma du proton
  hadron->v_beta2[vALL] = 1 - (1/(hadron->v_gamma[vALL]*hadron->v_gamma[vALL]));	// Betta au carré (v/c)^2

  // Energie maximum transferable à l'e-
  //hadron->v_Te_max[vALL] = 	(2*MC2_ELEC * (hadron->v_gamma[vALL]*hadron->v_gamma[vALL] - 1)) / (1 + 2*hadron->v_gamma[vALL]*(MC2_ELEC/(hadron->v_mass[vALL]*MC2_PRO)) 
  //				+ (MC2_ELEC/(hadron->v_mass[vALL]*MC2_PRO)) * (MC2_ELEC/(hadron->v_mass[vALL]*MC2_PRO)));
  hadron->v_Te_max[vALL] = 	(2*MC2_ELEC * (hadron->v_mass[vALL]*MC2_PRO)*(hadron->v_mass[vALL]*MC2_PRO) * (hadron->v_gamma[vALL]*hadron->v_gamma[vALL] - 1)) / 
				((hadron->v_mass[vALL]*MC2_PRO)*(hadron->v_mass[vALL]*MC2_PRO) + 2*MC2_ELEC*hadron->v_gamma[vALL]*(hadron->v_mass[vALL]*MC2_PRO) 
				+ MC2_ELEC*MC2_ELEC);

  return;
}


void Copy_Hadron_struct(Hadron *destination, Hadron *origin){
  __assume_aligned(&destination->v_x, 64);
  __assume_aligned(&destination->v_y, 64);
  __assume_aligned(&destination->v_z, 64);

  __assume_aligned(&destination->v_u, 64);
  __assume_aligned(&destination->v_v, 64);
  __assume_aligned(&destination->v_w, 64);

  __assume_aligned(&destination->v_T, 64);
  __assume_aligned(&destination->v_M, 64);
  __assume_aligned(&destination->v_charge, 64);
  __assume_aligned(&destination->v_mass, 64);

  __assume_aligned(&destination->v_type, 64);

  __assume_aligned(&destination->v_E, 64);
  __assume_aligned(&destination->v_gamma, 64);
  __assume_aligned(&destination->v_beta2, 64);
  __assume_aligned(&destination->v_Te_max, 64);

  __assume_aligned(&origin->v_x, 64);
  __assume_aligned(&origin->v_y, 64);
  __assume_aligned(&origin->v_z, 64);

  __assume_aligned(&origin->v_u, 64);
  __assume_aligned(&origin->v_v, 64);
  __assume_aligned(&origin->v_w, 64);

  __assume_aligned(&origin->v_T, 64);
  __assume_aligned(&origin->v_M, 64);
  __assume_aligned(&origin->v_charge, 64);
  __assume_aligned(&origin->v_mass, 64);

  __assume_aligned(&origin->v_type, 64);

  __assume_aligned(&origin->v_E, 64);
  __assume_aligned(&origin->v_gamma, 64);
  __assume_aligned(&origin->v_beta2, 64);
  __assume_aligned(&origin->v_Te_max, 64);


  destination->v_x[:] = origin->v_x[:];
  destination->v_y[:] = origin->v_y[:];
  destination->v_z[:] = origin->v_z[:];

  destination->v_u[:] = origin->v_u[:];
  destination->v_v[:] = origin->v_v[:];
  destination->v_w[:] = origin->v_w[:];

  destination->v_T[:] = origin->v_T[:];
  destination->v_M[:] = origin->v_M[:];
  destination->v_charge[:] = origin->v_charge[:];
  destination->v_mass[:] = origin->v_mass[:];

  destination->v_type[:] = origin->v_type[:];

  destination->v_E[:] = origin->v_E[:];
  destination->v_gamma[:] = origin->v_gamma[:];
  destination->v_beta2[:] = origin->v_beta2[:];
  destination->v_Te_max[:] = origin->v_Te_max[:];
}
