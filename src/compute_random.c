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

void Init_RND(DATA_config *config, VAR_RND_SEED RNDstream, int offset){

  #if USE_MKL_LIB==1

    if(config->RNG_Seed == 0){
      vslNewStream(&RNDstream, VSL_BRNG_MCG59, time(NULL)+offset);	// initialisation du stream du RNG avec le seed (time+thread_id)
    }
    else{
      vslNewStream(&RNDstream, VSL_BRNG_MCG59, config->RNG_Seed+offset);
    }

  #else

    if(config->RNG_Seed == 0){
      *RNDstream = time(NULL)+offset;	// initialisation du stream du RNG avec le seed (time+thread_id)
    }
    else{
      *RNDstream = config->RNG_Seed+offset;
    }

  #endif

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNDstream, v_rnd);				// on genere une première fois un set de nbr car les premiers semblent mal distribués

}


void rand_uniform(VAR_RND_SEED seedp, VAR_COMPUTE *v_rnd){

  __assume_aligned(v_rnd, 64);

  #if USE_MKL_LIB==1

    #if VAR_COMPUTE_PRECISION==1
      vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, seedp, VLENGTH, v_rnd, FLT_EPSILON , (1.0-FLT_EPSILON) );
    #else
      vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, seedp, VLENGTH, v_rnd, DBL_EPSILON, (1.0-DBL_EPSILON) );
    #endif

  #else

    int i;
    for(i=0; i<VLENGTH; i++){
      v_rnd[i] = single_rand_uniform(seedp);
    }

  #endif

  return;
}


VAR_COMPUTE single_rand_uniform(VAR_RND_SEED seedp){

  VAR_COMPUTE rnd;

  #if USE_MKL_LIB==1

    #if VAR_COMPUTE_PRECISION==1
      vsRngUniform( VSL_RNG_METHOD_UNIFORM_STD, seedp, 1, &rnd, FLT_EPSILON, (1.0-FLT_EPSILON) );
    #else
      vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, seedp, 1, &rnd, DBL_EPSILON, (1.0-DBL_EPSILON) );
    #endif

  #else

    #if VAR_COMPUTE_PRECISION==1
      rnd = (1.0-FLT_EPSILON) * ((VAR_COMPUTE)rand_r(seedp) / RAND_MAX) + FLT_EPSILON;
    #else
      rnd = (1.0-DBL_EPSILON) * ((VAR_COMPUTE)rand_r(seedp) / RAND_MAX) + DBL_EPSILON;
    #endif

  #endif

  return rnd;
}


void rand_normal(VAR_RND_SEED seedp, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_mu, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_mu, 64);
  __assume_aligned(v_sigma, 64);

  #if USE_MKL_LIB==1

    #if VAR_COMPUTE_PRECISION==1										// Methods :
      vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, VLENGTH, v_rnd, 0.0, 1.0 );		// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    #else													// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
      vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, VLENGTH, v_rnd, 0.0, 1.0);		// VSL_RNG_METHOD_GAUSSIAN_ICDF
    #endif

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_rnd[v] = v_sigma[v] * v_rnd[v] + v_mu[v];
    }

  #else

    int i;
    for(i=0; i<VLENGTH; i++){
      v_rnd[i] = single_rand_normal(seedp, v_mu[i], v_sigma[i]);
    }

  #endif

  return;
}


void rand_normal_zero(VAR_RND_SEED seedp, VAR_COMPUTE *v_rnd, VAR_COMPUTE *v_sigma){

  __assume_aligned(v_rnd, 64);
  __assume_aligned(v_sigma, 64);

  #if USE_MKL_LIB==1

    #if VAR_COMPUTE_PRECISION==1										// Methods :
      vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, VLENGTH, v_rnd, 0.0, 1.0 );		// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    #else													// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
      vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, VLENGTH, v_rnd, 0.0, 1.0);		// VSL_RNG_METHOD_GAUSSIAN_ICDF
    #endif

    #pragma omp simd
    for(int v = 0; v<VLENGTH; v++){
      v_rnd[v] = v_sigma[v] * v_rnd[v];
    }

  #else

    int i;
    for(i=0; i<VLENGTH; i++){
      v_rnd[i] = single_rand_normal(seedp, 0.0, v_sigma[i]);
    }

  #endif

  return;
}


VAR_COMPUTE single_rand_normal(VAR_RND_SEED seedp, VAR_COMPUTE mu, VAR_COMPUTE sigma){

  #if USE_MKL_LIB==1

    VAR_COMPUTE rnd;

    #if VAR_COMPUTE_PRECISION==1
      vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, 1, &rnd, 0.0, 1.0 );
    #else
     vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, seedp, 1, &rnd, 0.0, 1.0);
    #endif

    rnd = sigma * rnd + mu;
    return rnd;

  #else

    VAR_COMPUTE rnd1 = single_rand_uniform(seedp); 
    VAR_COMPUTE rnd2 = single_rand_uniform(seedp);

    return sigma * sqrt(-2*log(rnd1))*cos(2*M_PI*rnd2) + mu;
    // 2nd number: sqrt(-2*log(u1))*sin(2*M_PI*u2)

  #endif

}
