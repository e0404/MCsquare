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


void Init_RND(DATA_config *config, VAR_RND_SEED* RNDstream, int offset){

  #if USE_MKL_LIB==1

    if(config->RNG_Seed == 0){
      vslNewStream(RNDstream, VSL_BRNG_MCG59, time(NULL)+offset);	// initialize the RNG for each thread individually with a seed = time+thread_id*10000
    }
    else{
      vslNewStream(RNDstream, VSL_BRNG_MCG59, config->RNG_Seed+offset);
    }

  #else

    if(config->RNG_Seed == 0){
      //**RNDstream = (VAR_RND_SEED) time(NULL)+offset;	// initialize the RNG for each thread individually with a seed = time+offset
      pcg32_init((VAR_RND_SEED_TYPE) (time(NULL) + offset), *RNDstream);
    }
    else{
      //**RNDstream = (VAR_RND_SEED) config->RNG_Seed+offset;
      pcg32_init((VAR_RND_SEED_TYPE) (config->RNG_Seed + offset), *RNDstream);
    }

  #endif

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(*RNDstream, v_rnd);				// the RNG is called here because random numbers seems not well distributed the first time.

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
      rnd = (1.0 - FLT_EPSILON) * ((VAR_COMPUTE) pcg32(seedp) / UINT32_MAX) + FLT_EPSILON;
    #else
      rnd = (1.0 - DBL_EPSILON) * ((VAR_COMPUTE) pcg32(seedp) / UINT32_MAX) + DBL_EPSILON;
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

    int v;
    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      v_rnd[v] = v_sigma[v] * v_rnd[v] + v_mu[v];
    }

  #else

  int i;
  for (i = 0; i < VLENGTH; i += 2) {
      //v_rnd[i] = single_rand_normal(seedp, 0.0, v_sigma[i]);
      box_muller_rand_normal(seedp, &v_rnd[i], &v_rnd[i + 1], v_mu[i], v_mu[i+1], v_sigma[i], v_sigma[i + 1]);
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

    int v;
    #pragma omp simd
    for(v = 0; v<VLENGTH; v++){
      v_rnd[v] = v_sigma[v] * v_rnd[v];
    }

  #else

    int i;
    for(i=0; i<VLENGTH; i+=2){
        //v_rnd[i] = single_rand_normal(seedp, 0.0, v_sigma[i]);
        box_muller_rand_normal(seedp, &v_rnd[i],&v_rnd[i+1], 0.0,0.0, v_sigma[i], v_sigma[i+1]);
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

    VAR_COMPUTE rnd1, rnd2;
    box_muller_rand_normal(seedp, &rnd1, &rnd2, mu, 0.0, sigma, 1.0);
    return rnd1;

  #endif

}

//Direct call to box muller transform.
void box_muller_rand_normal(VAR_RND_SEED seedp, VAR_COMPUTE* rnd1, VAR_COMPUTE* rnd2, const VAR_COMPUTE mu1, const  VAR_COMPUTE mu2, const VAR_COMPUTE sigma1, const  VAR_COMPUTE sigma2)
{
    VAR_COMPUTE rnd1_uni = single_rand_uniform(seedp);
    VAR_COMPUTE rnd2_uni = single_rand_uniform(seedp);

    VAR_COMPUTE r = sqrt(-2.0 * log(rnd1_uni));
    VAR_COMPUTE phi = 2 * M_PI * rnd2_uni;

    *rnd1 = sigma1 * r * cos(phi) + mu1;
    *rnd2 = sigma2 * r * sin(phi) + mu2;
}

//PCG-XSH-RR
uint32_t rotr32(uint32_t x, unsigned r)
{
    return x >> r | x << (-r & 31);
}

uint32_t pcg32(uint64_t* seedp)
{
    uint64_t const multiplier = 6364136223846793005u;
    uint64_t const increment = 1442695040888963407u;

    uint64_t x = *seedp;
    unsigned count = (unsigned)(x >> 59);
    *seedp = x * multiplier + increment;
    x ^= x >> 18;

    return rotr32((uint32_t)(x >> 27), count);
}

void pcg32_init(uint64_t seed, uint64_t* seedp)
{
    uint64_t const increment = 1442695040888963407u;
    *seedp = seed + increment;
    (void)pcg32(seedp); //initial run
}



