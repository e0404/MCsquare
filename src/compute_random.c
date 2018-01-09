/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_random.h"

// float encode entre 1E-38 et 1E38
// float encode entre 1E-308 et 1E308


void rand_uniform(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd){

  __assume_aligned(v_rnd, 64);

  #if VAR_COMPUTE_PRECISION==1
    vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, VLENGTH, v_rnd, FLT_EPSILON , (1.0-FLT_EPSILON) );
  #else
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, VLENGTH, v_rnd, DBL_EPSILON, (1.0-DBL_EPSILON) );
  #endif

  return;
}


VAR_COMPUTE single_rand_uniform(VSLStreamStatePtr stream){

  VAR_COMPUTE rnd;

  #if VAR_COMPUTE_PRECISION==1
    vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &rnd, FLT_EPSILON, (1.0-FLT_EPSILON) );
  #else
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &rnd, DBL_EPSILON, (1.0-DBL_EPSILON) );
  #endif

  return rnd;
}


void rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_mu, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_mu, 64);
  __assume_aligned(v_sigma, 64);

  #if VAR_COMPUTE_PRECISION==1										// Methods :
    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, VLENGTH, v_rnd, 0.0, 1.0 );		// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
  #else													// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, VLENGTH, v_rnd, 0.0, 1.0);		// VSL_RNG_METHOD_GAUSSIAN_ICDF
  #endif

  v_rnd[vALL] = v_sigma[vALL] * v_rnd[vALL] + v_mu[vALL];

  return;
}


void rand_normal_zero(VSLStreamStatePtr stream, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_sigma, 64);

  #if VAR_COMPUTE_PRECISION==1										// Methods :
    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, VLENGTH, v_rnd, 0.0, 1.0 );		// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
  #else													// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, VLENGTH, v_rnd, 0.0, 1.0);		// VSL_RNG_METHOD_GAUSSIAN_ICDF
  #endif

  v_rnd[vALL] = v_sigma[vALL] * v_rnd[vALL];

  return;
}


VAR_COMPUTE single_rand_normal(VSLStreamStatePtr stream, VAR_COMPUTE mu, VAR_COMPUTE sigma){

  VAR_COMPUTE rnd;

  #if VAR_COMPUTE_PRECISION==1
    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, 1, &rnd, 0.0, 1.0 );
  #else
   vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, 1, &rnd, 0.0, 1.0);
  #endif

  rnd = sigma * rnd + mu;

  return rnd;
}
