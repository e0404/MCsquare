/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_math.h"


void vector_floor(VAR_COMPUTE *v_vec, VAR_COMPUTE *v_result){

  __assume_aligned(v_vec, 64);
  __assume_aligned(v_result, 64);

  #if VAR_COMPUTE_PRECISION==1
    vsFloor(VLENGTH, v_vec, v_result);
  #else
    vdFloor(VLENGTH, v_vec, v_result);
  #endif

  return;
}


int sign(VAR_COMPUTE x){
/*
  if(x > 0) return 1;
  else if(x < 0) return -1;
  else return 0;
*/
 
  return (x > 0.0) - (x < 0.0);
 
}


void vec_sign(VAR_COMPUTE *v_vec, VAR_COMPUTE *v_sign){

  __assume_aligned(v_vec, 64);
  __assume_aligned(v_sign, 64);

/*
  if(v_vec[vALL] > 0) v_sign[vALL] = 1;
  else if(v_vec[vALL] < 0) v_sign[vALL] = -1;
  else v_sign[vALL] = 0;
*/

  v_sign[vALL] = (v_vec[vALL] > 0.0) - (v_vec[vALL] < 0.0);

  return;
}


VAR_COMPUTE SolidAngle(VAR_COMPUTE theta){
  return 2*M_PI*(1-cos(theta*M_PI/180));
}


VAR_COMPUTE Linear_Interpolation(VAR_COMPUTE x, VAR_COMPUTE x1, VAR_COMPUTE x2, VAR_COMPUTE y1, VAR_COMPUTE y2){
  // Calcul de la pente p = (y2 - y1) / (x2 - x1)
  // Calcul de la valeur interpolée y = p*(x-x1) + y1

  return ((y2 - y1) / (x2 - x1))*(x - x1) + y1;
}


void vec_Linear_Interpolation(VAR_COMPUTE *v_x, VAR_COMPUTE *v_x1, VAR_COMPUTE *v_x2, VAR_COMPUTE *v_y1, VAR_COMPUTE *v_y2, VAR_COMPUTE *v_result){

  __assume_aligned(v_x, 64);
  __assume_aligned(v_x1, 64);
  __assume_aligned(v_x2, 64);
  __assume_aligned(v_y1, 64);
  __assume_aligned(v_y2, 64);
  __assume_aligned(v_result, 64);

  // Calcul de la pente p = (y2 - y1) / (x2 - x1)
  // Calcul de la valeur interpolée y = p*(x-x1) + y1

  v_result[vALL] = ((v_y2[vALL] - v_y1[vALL]) / (v_x2[vALL] - v_x1[vALL]))*(v_x[vALL] - v_x1[vALL]) + v_y1[vALL];

  return;
}


VAR_COMPUTE Trilinear_Interpolation(VAR_DATA *x, VAR_DATA *image, int *GridSize){


  // We assume that the deformation vector x is given in voxel units with the same spacing than the image.

  int Id_x = floor(x[0]);
  VAR_COMPUTE X = x[0];
  if(Id_x > GridSize[0]-2){
    Id_x = GridSize[0]-2;
    X = GridSize[0]-1;
  }
  else if(Id_x < 0){
    Id_x = 0;
    X = 0;
  }

  int Id_y = floor(x[1]);
  VAR_COMPUTE Y = x[1];
  if(Id_y > GridSize[1]-2){
    Id_y = GridSize[1]-2;
    Y = GridSize[1]-1;
  }
  else if(Id_y < 0){
    Id_y = 0;
    Y = 0;
  }

  int Id_z = floor(x[2]);
  VAR_COMPUTE Z = x[2];
  if(Id_z > GridSize[2]-2){
    Id_z = GridSize[2]-2;
    Z = GridSize[2]-1;
  }
  else if(Id_z < 0){
    Id_z = 0;
    Z = 0;
  }

  int Id1 = (Id_x) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z);
  int Id2 = (Id_x+1) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z);
  VAR_COMPUTE C00 = (image[Id2] - image[Id1]) * (X - Id_x) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z);
  VAR_COMPUTE C10 = (image[Id2] - image[Id1]) * (X - Id_x) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z+1);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z+1);
  VAR_COMPUTE C01 = (image[Id2] - image[Id1]) * (X - Id_x) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z+1);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z+1);
  VAR_COMPUTE C11 = (image[Id2] - image[Id1]) * (X - Id_x) + image[Id1];

  VAR_COMPUTE C0 = (C10 - C00) * (Y - Id_y) + C00;
  VAR_COMPUTE C1 = (C11 - C01) * (Y - Id_y) + C01;

  return (C1 - C0) * (Z - Id_z) + C0;

/*
  // We assume that the deformation vector x is given in distance units with the same spacing than the image.

  int Id_x = floor(x[0] / Spacing[0]);
  if(Id_x > GridSize[0]-2) Id_x = GridSize[0]-2;
  else if(Id_x < 0) Id_x = 0;

  int Id_y = floor(x[1] / Spacing[1]);
  if(Id_y > GridSize[1]-2) Id_y = GridSize[1]-2;
  else if(Id_y < 0) Id_y = 0;

  int Id_z = floor(x[2] / Spacing[2]);
  if(Id_z > GridSize[2]-2) Id_z = GridSize[2]-2;
  else if(Id_z < 0) Id_z = 0;

  int Id1 = (Id_x) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z);
  int Id2 = (Id_x+1) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z);
  VAR_COMPUTE C00 = ((image[Id2] - image[Id1]) / Spacing[0]) * (x[0] - Id_x*Spacing[0]) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z);
  VAR_COMPUTE C10 = ((image[Id2] - image[Id1]) / Spacing[0]) * (x[0] - Id_x*Spacing[0]) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z+1);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y) + GridSize[0]*GridSize[1]*(Id_z+1);
  VAR_COMPUTE C01 = ((image[Id2] - image[Id1]) / Spacing[0]) * (x[0] - Id_x*Spacing[0]) + image[Id1];

  Id1 = (Id_x) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z+1);
  Id2 = (Id_x+1) + GridSize[0]*(Id_y+1) + GridSize[0]*GridSize[1]*(Id_z+1);
  VAR_COMPUTE C11 = ((image[Id2] - image[Id1]) / Spacing[0]) * (x[0] - Id_x*Spacing[0]) + image[Id1];

  VAR_COMPUTE C0 = ((C10 - C00) / Spacing[1]) * (x[1] - Id_y*Spacing[1]) + C00;
  VAR_COMPUTE C1 = ((C11 - C01) / Spacing[1]) * (x[1] - Id_y*Spacing[1]) + C01;

  return ((C1 - C0) / Spacing[2]) * (x[2] - Id_z*Spacing[2]) + C0;
*/

}


int Sequential_Search(VAR_COMPUTE value, VAR_DATA *list, int Nbr_Values){
  int index = 0;

  if(Nbr_Values <= 1) return -1;

  while(index < Nbr_Values){
    if(value <= list[index]) return index - 1;
    else index++;
  }
  return index - 1;
}


int Binary_Search(VAR_COMPUTE value, VAR_DATA *list, int Nbr_Values){
  int k, i=-1, j=Nbr_Values;

  if(Nbr_Values <= 1) return -1;

  while((j-i) > 1){
    k = (i+j)/2;
    if(list[k] < value) i=k;
    else j=k;
  }
  
  return i;
}


inline void my_log(VAR_COMPUTE *v_data, VAR_COMPUTE *v_result){

  __assume_aligned(v_data, 64);
  __assume_aligned(v_result, 64);

//  v_result[vALL] = log(v_data[vALL]);

//  vsLog10( VLENGTH, v_data, v_result );

  union { 
	ALIGNED_(64) float f[VLENGTH];
	ALIGNED_(64) uint32_t i[VLENGTH]; 
	} vx;

  vx.f[vALL] = (float)v_data[vALL];

  union { 
	ALIGNED_(64) uint32_t i[VLENGTH]; 
	ALIGNED_(64) float f[VLENGTH];
	} mx;

  mx.i[vALL] = (vx.i[vALL] & 0x007FFFFF) | 0x3f000000;

  v_result[vALL] = vx.i[vALL];
  v_result[vALL] *= 1.1920928955078125e-7f;
  v_result[vALL] = v_result[vALL] - 124.22551499f - 1.498030302f * mx.f[vALL] - 1.72587999f / (0.3520887068f + mx.f[vALL]);
  v_result[vALL] *= 0.69314718f;

  return;
}


VAR_COMPUTE find_min(VAR_COMPUTE *data, int NumElements){

  VAR_COMPUTE minimum = FLT_MAX; // 1E+37
  int i;

  for(i=0; i<NumElements; i++){
    if(data[i] < minimum) minimum = data[i];
  }

  return minimum;
}


VAR_COMPUTE find_max(VAR_COMPUTE *data, int NumElements){

  VAR_COMPUTE maximum = -FLT_MAX; // -1E+37
  int i;

  for(i=0; i<NumElements; i++){
    if(data[i] > maximum) maximum = data[i];
  }

  return maximum;
}


VAR_COMPUTE compute_mean(VAR_COMPUTE *data, int NumElements){

  VAR_COMPUTE sum = 0;
  int i;

  for(i=0; i<NumElements; i++) sum += data[i];

  return sum / NumElements;
}


VAR_COMPUTE compute_median(VAR_COMPUTE *data, int NumElements){

  VAR_COMPUTE tmp;
  int i, j;

  // sort the array in ascending order
  for(i=0; i<NumElements-1; i++) {
    for(j=i+1; j<NumElements; j++) {
      if(data[j] < data[i]) {
        tmp = data[i];
        data[i] = data[j];  // swap elements
        data[j] = tmp;
      }
    }
  }

  if(NumElements%2 == 0) return ((data[NumElements/2] + data[NumElements/2 - 1]) / 2.0);
  else return data[NumElements/2];

}


