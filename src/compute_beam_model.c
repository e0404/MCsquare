/*
This file is part of the MCsquare software
Copyright © 2016-2017 Université catholique de Louvain (UCL)
All rights reserved.

The MCsquare software has been developed by Kevin Souris from UCL in the context of a collaboration with IBA s.a.
Each use of this software must be attributed to Université catholique de Louvain (UCL, Louvain-la-Neuve). Any other additional authorizations may be asked to LTTO@uclouvain.be.
The MCsquare software is released under the terms of the open-source Apache 2.0 license. Anyone can use or modify the code provided that the Apache 2.0 license conditions are met. See the Apache 2.0 license for more details https://www.apache.org/licenses/LICENSE-2.0
The MCsquare software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*/


#include "include/compute_beam_model.h"


double ConvertMuToProtons(double weight, double energy)
{

/*  http://thread.gmane.org/gmane.comp.science.opengate.user/1564/focus=1572

A MU is defined as a number of nC collected in the Ionization Chamber (IC), which is filled with air.
SP corresponds to proton stopping power in air and it is based on a fit from ICRU data.
K is a constant which depends on the mean energy loss (W) to create an electron/hole pair.
PTP are the temperature and pression corrections.
Using all these parameters correctly allows for MU to absolute number of protons conversion.

*/

// Constant which depends on the mean energy loss (W) to create an electron/hole pair
double K = 35.87; // in eV (other value 34.23 ?)

// Air stopping power (fit ICRU) multiplied by air density
double SP = (9.6139e-9*pow(energy,4) - 7.0508e-6*pow(energy,3) + 2.0028e-3*pow(energy,2) - 2.7615e-1*energy + 2.0082e1) * 1.20479E-3 * 1E6; // in eV / cm

// Temp & Pressure correction
double PTP = 1.0; 

// MU calibration (1 MU = 3 nC/cm)
// 1cm de gap effectif
double C = 3.0E-9; // in C / cm

// Gain: 1eV = 1.602176E-19 J
double Gain = (C*K) / (SP*PTP*1.602176E-19);

return weight*Gain;


/*
// Loic's formula (not correct ?)

    double K=37.60933;
    double SP=9.6139E-09*pow(energy,4)-7.0508E-06*pow(energy,3)+2.0028E-03*pow(energy,2)-2.7615E-01*pow(energy,1)+2.0082E+01*pow(energy,0);
    double PTP=1;
    double Gain=3./(K*SP*PTP*1.602176E-10);
    return (weight*Gain);
*/
}


void deviates (double U[2][2], double sigmas[2], double R[2], VSLStreamStatePtr RNG_Stream)
{
    // Returns vector of gaussian randoms based on sigmas, rotated by U,
    // with means of 0.

// NE FONCTIONNE PAS  (car les nombres gaussiens générés doivent probablement être compris entre -1 et 1)?
//		      (ou bien la distribution gaussienne 2D =! de 2 distribution gaussiennes 1D) ?
/*
    VAR_COMPUTE rnd[2];
    
  #if VAR_COMPUTE_PRECISION==1									// Methods :
    vsRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, RNG_Stream, 2, rnd, 0.0, 1.0 );		// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
  #else												// VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2
    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, RNG_Stream, 2, rnd, 0.0, 1.0);		// VSL_RNG_METHOD_GAUSSIAN_ICDF
  #endif
    
    R[0] = U[0][0]*rnd[0] + U[0][1]*rnd[1];
    R[1] = U[1][0]*rnd[0] + U[1][1]*rnd[1];
*/

    int n = 2;
    double V[2];  // The vector to be returned is R now
    
    double r,v1,v2,fac;
    
    int i = 1;
    
    while ( i <= n ) {
        do {
            v1 = 2.0 * single_rand_uniform(RNG_Stream) - 1.0;
            v2 = 2.0 * single_rand_uniform(RNG_Stream) - 1.0;
            r = v1*v1 + v2*v2;
        } while ( r > 1.0 );
        fac = sqrt(-2.0*log(r)/r);
        V[i-1] = v1*fac;
        i++;
        if ( i <= n ) {
            V[i-1] = v2*fac;
            i++;
        } 
    } 
    
    for ( i = 0; i < n; i++ ) {
        V[i] *= sigmas[i];
    }
    
    R[0]=U[0][0]*V[0]+U[0][1]*V[1];
    R[1]=U[1][0]*V[0]+U[1][1]*V[1];

    
}



void diagonalize (double A[2][2],double B[2], double T[2][2])
{
    double b=-A[0][0]-A[1][1];
    double c=A[1][1]*A[0][0]-A[0][1]*A[1][0];
    B[0]=0.5*(-b+sqrt(b*b-4*c));
    B[1]=0.5*(-b-sqrt(b*b-4*c));
    double r_1=sqrt(1+(A[1][0]/(A[1][1]-B[0]))*(A[1][0]/(A[1][1]-B[0])));
    double r_2=sqrt(1+(A[1][0]/(A[1][1]-B[1]))*(A[1][0]/(A[1][1]-B[1])));
    T[0][0]=1/r_1;
    T[0][1]=1/r_2;
    T[1][0]=1/r_1*(-A[1][0]/(A[1][1]-B[0]));
    T[1][1]=1/r_2*(-A[1][0]/(A[1][1]-B[1]));
    
//    printf("%f %f %f %f %f %f\n",r_1, r_2,T[1][1],T[1][0],B[0],B[1]);
}


double compute_mEnergy (machine_parameters *mac, double energy)
{
	// Ce qui compte c'est le range des protons.  
	// Or, l'énergie configurée sur la machine pour atteindre un certain range n'est pas forcément la même que l'énergie Monte Carlo.
	// Il faut donc convertir l'énergie du plan à l'aide d'un fit.  Les paramètres du fit doivent être recalculé pour les différentes machines et différents modèles MC.

    int i=0;
    double val=0;
    for (i=0;i<mac->mEnergy_order;i++)
    {
        val+=mac->mEnergy_poly[i]*pow(energy,mac->mEnergy_order-i);        
    }
    val+=mac->mEnergy_poly[mac->mEnergy_order];

    return val;
}


double compute_sEnergy (machine_parameters *mac, double energy)
{
    int i=0;
    double val=0;
    for (i=0;i<mac->sEnergy_order;i++)
    {
        val+=mac->sEnergy_poly[i]*pow(energy,mac->sEnergy_order-i);
    }
    val+=mac->sEnergy_poly[mac->sEnergy_order];
    
    return val;
}


double compute_mX (machine_parameters *mac, double energy)
{
    int i=0;
    double val=0;
    for (i=0;i<mac->mX_order;i++)
    {
        val+=mac->mX_poly[i]*pow(energy,mac->mX_order-i);
    }
    val+=mac->mX_poly[mac->mX_order];
    
    return val;
}


double compute_mY (machine_parameters *mac, double energy)
{
    int i=0;
    double val=0;
    for (i=0;i<mac->mY_order;i++)
    {
        val+=mac->mY_poly[i]*pow(energy,mac->mY_order-i);
    }
    val+=mac->mY_poly[mac->mY_order];
    
    return val;
}

double compute_mTheta (machine_parameters *mac, double energy)
{
    int i=0;
    double val=0;
    for (i=0;i<mac->mTheta_order;i++)
    {
        val+=mac->mTheta_poly[i]*pow(energy,mac->mTheta_order-i);
    }
    val+=mac->mTheta_poly[mac->mTheta_order];
    
    return val;
}

double compute_mPhi (machine_parameters *mac, double energy)
{
    int i=0;
    double val=0;
    for (i=0;i<mac->mPhi_order;i++)
    {
        val+=mac->mPhi_poly[i]*pow(energy,mac->mPhi_order-i);
    }
    val+=mac->mPhi_poly[mac->mPhi_order];
    
    return val;
}

double compute_eXTheta (machine_parameters *mac, double energy)
{
/*
    int i=0;
    double val=0;
    for (i=0;i<mac->eXTheta_order;i++)
    {
        val+=mac->eXTheta_poly[i]*pow(energy,mac->eXTheta_order-i);
    }
    val+=mac->eXTheta_poly[mac->eXTheta_order];
    
    return val;
*/
    
    return 0.5*M_PI * compute_mX(mac, energy) * compute_mTheta(mac, energy);
}

double compute_eYPhi (machine_parameters *mac, double energy)
{
/*
    int i=0;
    double val=0;
    for (i=0;i<mac->eYPhi_order;i++)
    {
        val+=mac->eYPhi_poly[i]*pow(energy,mac->eYPhi_order-i);
    }
    val+=mac->eYPhi_poly[mac->eYPhi_order];
    
    return val;
*/

    return 0.5*M_PI * compute_mY(mac, energy) * compute_mPhi(mac, energy);
}

void rotateX (double angle, double vector[3])
{
    double cx=cos(angle);
    double sx=sin(angle);
    double tempX=vector [0];
    double tempY=vector [1];
    double tempZ=vector [2];
    
    vector[0]= tempX;
    vector[1]=tempY*cx - tempZ*sx;
    vector[2]=tempY*sx + tempZ*cx;
/*
    double norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    vector[0]/=norm;
    vector[1]/=norm;
    vector[2]/=norm;
*/
}

void rotateY (double angle, double vector[3])
{
    double cx=cos(angle);
    double sx=sin(angle);
    double tempX=vector [0];
    double tempY=vector [1];
    double tempZ=vector [2];
    
    vector[0]= tempX*cx + tempZ*sx;
    vector[1]= tempY;
    vector[2]=-tempX*sx + tempZ*cx;
/*
    double norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    vector[0]/=norm;
    vector[1]/=norm;
    vector[2]/=norm;
*/
}

void rotateZ (double angle, double vector[3])
{
    double cx=cos(angle);
    double sx=sin(angle);
    double tempX=vector [0];
    double tempY=vector [1];
    double tempZ=vector [2];
    
    vector[0]= tempX*cx - tempY*sx;
    vector[1]= tempX*sx + tempY*cx;
    vector[2]= tempZ;
/*
    double norm=sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    vector[0]/=norm;
    vector[1]/=norm;
    vector[2]/=norm;
*/
}


void Sample_particle (Hadron_buffer *hadron, VAR_DATA CT_Length[3], machine_parameters *mac, ControlPoint_parameters *ControlPoint, spot_parameters *spot, VSLStreamStatePtr RNG_Stream)
{
   //here we need to sample particle parameters    

    double E;
    double A[2][2];
    double T[2][2];
    double XTheta[2];
    double YPhi[2];
    double sigmas[2];


    if(mac->Beam_Model == UPenn){
	
	int EnergyID =  Sequential_Search(ControlPoint->Energy, mac->Nominal_Energies, mac->Number_Energies);
	if(EnergyID < 0) EnergyID = 0;
	if(EnergyID > (mac->Number_Energies - 2)) EnergyID = mac->Number_Energies - 2;
	VAR_COMPUTE Energy1 = mac->Nominal_Energies[EnergyID];
	VAR_COMPUTE Energy2 = mac->Nominal_Energies[EnergyID+1];

	E = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Mean_Energies[EnergyID], mac->Mean_Energies[EnergyID+1]);
	VAR_COMPUTE sE = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Energy_Spread[EnergyID], mac->Energy_Spread[EnergyID+1]);
	E = single_rand_normal(RNG_Stream, E, sE*ControlPoint->Energy/100);

	VAR_COMPUTE SpotSizeX, DivergenceX, CorrelationX, SpotSizeY, DivergenceY, CorrelationY;

	VAR_COMPUTE weight1 = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Weight1[EnergyID], mac->Weight1[EnergyID+1]);
	VAR_COMPUTE weight2 = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Weight2[EnergyID], mac->Weight2[EnergyID+1]);
	VAR_COMPUTE rnd = single_rand_uniform(RNG_Stream);
    	rnd = rnd * (weight1 + weight2);
	if(rnd < weight1){
	  SpotSizeX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->SpotSize1x[EnergyID], mac->SpotSize1x[EnergyID+1]);
	  DivergenceX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Divergence1x[EnergyID], mac->Divergence1x[EnergyID+1]);
	  CorrelationX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Correlation1x[EnergyID], mac->Correlation1x[EnergyID+1]);
	  SpotSizeY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->SpotSize1y[EnergyID], mac->SpotSize1y[EnergyID+1]);
	  DivergenceY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Divergence1y[EnergyID], mac->Divergence1y[EnergyID+1]);
	  CorrelationY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Correlation1y[EnergyID], mac->Correlation1y[EnergyID+1]);
	}
	else{
	  SpotSizeX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->SpotSize2x[EnergyID], mac->SpotSize2x[EnergyID+1]);
	  DivergenceX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Divergence2x[EnergyID], mac->Divergence2x[EnergyID+1]);
	  CorrelationX = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Correlation2x[EnergyID], mac->Correlation2x[EnergyID+1]);
	  SpotSizeY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->SpotSize2y[EnergyID], mac->SpotSize2y[EnergyID+1]);
	  DivergenceY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Divergence2y[EnergyID], mac->Divergence2y[EnergyID+1]);
	  CorrelationY = Linear_Interpolation(ControlPoint->Energy, Energy1, Energy2, mac->Correlation2y[EnergyID], mac->Correlation2y[EnergyID+1]);
	}

	// Sample X direction
	A[0][0] = SpotSizeX * SpotSizeX;
	A[1][1] = DivergenceX * DivergenceX;
	A[0][1] = CorrelationX * SpotSizeX * DivergenceX;
	A[1][0] = A[0][1];

	diagonalize(A,sigmas,T);
	sigmas[0]=sqrt(sigmas[0]);
	sigmas[1]=sqrt(sigmas[1]);
	deviates(T,sigmas,XTheta, RNG_Stream);

	// Sample Y direction
	A[0][0] = SpotSizeY * SpotSizeY;
	A[1][1] = DivergenceY * DivergenceY;
	A[0][1] = CorrelationY * SpotSizeY * DivergenceY;
	A[1][0] = A[0][1];

	diagonalize(A,sigmas,T);
	sigmas[0]=sqrt(sigmas[0]);
	sigmas[1]=sqrt(sigmas[1]);
	deviates(T,sigmas,YPhi, RNG_Stream);

    }

    else{

	E = compute_mEnergy(mac,ControlPoint->Energy);

	double E_corr = ControlPoint->Energy;
	if(ControlPoint->Energy < 70.0){
	  E_corr = 70.0;
	}

	double sE = compute_sEnergy (mac,E_corr);
	double sX=compute_mX(mac,E_corr);
	double sY=compute_mY(mac,E_corr);
	double sTheta=compute_mTheta(mac,E_corr);
	double sPhi=compute_mPhi(mac,E_corr);
	double epsilonXTheta=compute_eXTheta(mac,E_corr);
	double epsilonYPhi=compute_eYPhi(mac,E_corr);

	// sampling of beam characteristics
	E = single_rand_normal(RNG_Stream, E, sE*ControlPoint->Energy/100);

	// Sample X direction
	epsilonXTheta/=M_PI;
	double beta=sX*sX/epsilonXTheta;
	double gamma=sTheta*sTheta/epsilonXTheta;
	double alpha=-sqrt(beta*gamma-1.);
    
	A[0][0]=sX*sX;
	A[1][1]=sTheta*sTheta;
	A[0][1]=-alpha*epsilonXTheta;
	A[1][0]=A[0][1];    
    
	diagonalize(A,sigmas,T);
	sigmas[0]=sqrt(sigmas[0]);
	sigmas[1]=sqrt(sigmas[1]);
	deviates(T,sigmas,XTheta, RNG_Stream);

	// Sample Y direction
	epsilonYPhi/=M_PI;
	beta=sY*sY/epsilonYPhi;
	gamma=sPhi*sPhi/epsilonYPhi;
	alpha=-sqrt(beta*gamma-1.);
    
	A[0][0]=sY*sY;
	A[1][1]=sPhi*sPhi;
	A[0][1]=-alpha*epsilonYPhi;
	A[1][0]=A[0][1];
    
	diagonalize(A,sigmas,T);
	sigmas[0]=sqrt(sigmas[0]);
	sigmas[1]=sqrt(sigmas[1]);
	deviates(T,sigmas,YPhi, RNG_Stream);

    }

    // Spot eye view to Beam eye view:
    //---- at this point, we have all information of particle characteristics
    //---- we perform all preliminary calculations according to BEV reference frame.
    
    double rotation[3];
    rotation[0]= M_PI+atan(spot->Spot_Y/mac->mDistanceSMYToIsocenter);
    rotation[1]= -atan(spot->Spot_X/mac->mDistanceSMXToIsocenter);
    rotation[2]=0.0;



    double particle_position[3];
    particle_position[0]=XTheta[0];
    particle_position[1]=YPhi[0];
    particle_position[2]=0.0;
//printf("\nSample: init particle_position: [%.3f ; %.3f ; %.3f]", particle_position[0], particle_position[1], particle_position[2]);

    rotateX(rotation[0],particle_position);
    rotateY(rotation[1],particle_position);
    rotateZ(rotation[2],particle_position);

    particle_position[0] += spot->Spot_X*(mac->mDistanceSMXToIsocenter - mac->mDistanceSourcePatient)/mac->mDistanceSMXToIsocenter;
    particle_position[1] += spot->Spot_Y*(mac->mDistanceSMYToIsocenter - mac->mDistanceSourcePatient)/mac->mDistanceSMYToIsocenter;
    particle_position[2] += mac->mDistanceSourcePatient;
//printf("\nSample: add source distance: [%.3f ; %.3f ; %.3f]", particle_position[0], particle_position[1], particle_position[2]);

    
    double particle_direction[3];
    particle_direction[0]=tan(XTheta[1]);
    particle_direction[1]=tan(YPhi[1]);
    particle_direction[2]=1.0;

    double norm=sqrt(particle_direction[0]*particle_direction[0] + particle_direction[1]*particle_direction[1] + particle_direction[2]*particle_direction[2]);

    particle_direction[0]/=norm;
    particle_direction[1]/=norm;
    particle_direction[2]/=norm;
    
    rotateX(rotation[0],particle_direction);
    rotateY(rotation[1],particle_direction);
    rotateZ(rotation[2],particle_direction);
    


    hadron->x = particle_position[0] / 10.0;
    hadron->y = particle_position[1] / 10.0;
    hadron->z = particle_position[2] / 10.0;
//printf("\nSample: convert to cm: [%.3f ; %.3f ; %.3f]", hadron->x, hadron->y, hadron->z);

    hadron->u = particle_direction[0];
    hadron->v = particle_direction[1];
    hadron->w = particle_direction[2];

    hadron->T = E * 1e6;	// energy in eV

    hadron->M = 1.0;
    hadron->charge = 1.0;
    hadron->mass = 1.0;
    hadron->type = Proton;

}



void BEV_to_CT_frame(Hadron_buffer *hadron, machine_parameters *mac, field_parameters *field){
    
    //---- we need now to express all particles properties to values consistent with MCsquare
    //---- after that, we express all variables in CT reference frame, which is used in MCsquare

    double particle_position[3];
    particle_position[0] = hadron->x;
    particle_position[1] = hadron->y;
    particle_position[2] = hadron->z;

    rotateY(field->GantryAngle, particle_position);
    rotateZ(-field->PatientSupportAngle, particle_position);

    // Beam model: XYZ (mm)
    // MCsquare: XZY (cm)
    hadron->x = particle_position[0] + (field->IsocenterPositionX / 10.0);
    hadron->y = particle_position[2] + (field->IsocenterPositionY / 10.0);
    hadron->z = particle_position[1] + (field->IsocenterPositionZ / 10.0);



    double particle_direction[3];
    particle_direction[0] = hadron->u;
    particle_direction[1] = hadron->v;
    particle_direction[2] = hadron->w;

    rotateY(field->GantryAngle, particle_direction);
    rotateZ(-field->PatientSupportAngle, particle_direction);

    // Beam model: XYZ (mm)
    // MCsquare: XZY (cm)
    hadron->u = particle_direction[0];
    hadron->v = particle_direction[2];
    hadron->w = particle_direction[1];

}



void Transport_to_CT(Hadron_buffer *hadron, VAR_DATA CT_Length[3]){

  if(hadron->x >= 0 && hadron->y >= 0 && hadron->z >= 0 && hadron->x <= CT_Length[0] && hadron->y <= CT_Length[1] && hadron->z <= CT_Length[2]) return;

  // Translate the particle to the CT image
  double Translation[3], new_position[3];

  if(hadron->u > 0)	Translation[0] = (0 - hadron->x) / hadron->u;
  else 			Translation[0] = (CT_Length[0] - hadron->x) / hadron->u;
  if(hadron->v > 0)	Translation[1] = (0 - hadron->y) / hadron->v;
  else 			Translation[1] = (CT_Length[1] - hadron->y) / hadron->v;
  if(hadron->w > 0)	Translation[2] = (0 - hadron->z) / hadron->w;
  else 			Translation[2] = (CT_Length[2] - hadron->z) / hadron->w;

  int i;
  for(i=0; i<3; i++){

    if(Translation[i] < 0) continue;


    Translation[i] += 1e-4;
    new_position[0] = hadron->x + Translation[i] * hadron->u;
    new_position[1] = hadron->y + Translation[i] * hadron->v;
    new_position[2] = hadron->z + Translation[i] * hadron->w;



    if(new_position[0] > 0.0 && new_position[1] > 0.0 && new_position[2] > 0.0 && new_position[0] < CT_Length[0] &&  new_position[1] < CT_Length[1] &&  new_position[2] < CT_Length[2] && !isnan(new_position[0]) && !isnan(new_position[1]) && !isnan(new_position[2])) break;
  }

  hadron->x = new_position[0];
  hadron->y = new_position[1];
  hadron->z = new_position[2];

  double energy = hadron->T / (UMeV*hadron->mass);
  double SP_air = hadron->charge * (9.6139e-9*pow(energy,4) - 7.0508e-6*pow(energy,3) + 2.0028e-3*pow(energy,2) - 2.7615e-1*energy + 2.0082e1) * 1.20479E-3 * 1E6;
  double dE = SP_air * Translation[i];
  
  hadron->T = hadron->T - dE;
  if(hadron->T < 0) hadron->type = Unknown;

}


void Transport_to_RangeShifter(Hadron_buffer *hadron, ControlPoint_parameters **layer_data, int Nbr_hadrons){

  int i;
  double IsocenterDistance;

  for(i=0; i<Nbr_hadrons; i++){
    if(layer_data[i]->RS_setting == OUT || layer_data[i]->RS_Thickness <= 0.0) continue;
    IsocenterDistance = layer_data[i]->RS_IsocenterDist + layer_data[i]->RS_Thickness;
//printf("\nTransportRS: IsoDist: [%.3f]", IsocenterDistance);
//printf("\nTransportRS: init position: [%.3f ; %.3f ; %.3f]", hadron[i].x, hadron[i].y, hadron[i].z);
//printf("\nTransportRS: direction: [%.3f ; %.3f ; %.3f]", hadron[i].u, hadron[i].v, hadron[i].w);

    hadron[i].x += hadron[i].u * (hadron[i].z - IsocenterDistance) / fabs(hadron[i].w);
    hadron[i].y += hadron[i].v * (hadron[i].z - IsocenterDistance) / fabs(hadron[i].w);
    hadron[i].z = IsocenterDistance;
//printf("\nTransportRS: post: [%.3f ; %.3f ; %.3f]", hadron[i].x, hadron[i].y, hadron[i].z);
  }

}


void Generate_PBS_particle(Hadron_buffer *hadron, int *Nbr_hadrons, VAR_DATA CT_Length[3], plan_parameters *plan, machine_parameters *machine, VSLStreamStatePtr RNG_Stream, DATA_config *config, Materials *material){

  ALIGNED_(64) VAR_COMPUTE v_rnd[VLENGTH];
  rand_uniform(RNG_Stream, v_rnd);
  v_rnd[vALL] = v_rnd[vALL] * plan->cumulative_weight;

  ALIGNED_(64) int v_field_index[VLENGTH];
  ALIGNED_(64) int v_ControlPoint_index[VLENGTH];
  ALIGNED_(64) int v_spot_index[VLENGTH];

  Hadron_buffer New_hadrons[50];
  field_parameters *New_hadrons_field[50];
  ControlPoint_parameters *New_hadrons_layer[50];
  int Nbr_New_hadrons = VLENGTH;

  int use_RS = 0;

  // Generate new particles according to PBS plan
  int i;
  for(i=0; i<Nbr_New_hadrons; i++){
    v_field_index[i] = Sequential_Search(v_rnd[i], plan->Fields_cumulative_PDF, plan->NumberOfFields) + 1;
    v_ControlPoint_index[i] = Binary_Search(v_rnd[i], plan->fields[v_field_index[i]].ControlPoints_cumulative_PDF, plan->fields[v_field_index[i]].NumberOfControlPoints) + 1;
    v_spot_index[i] = Binary_Search(v_rnd[i], plan->fields[v_field_index[i]].ControlPoints[v_ControlPoint_index[i]].Spots_cumulative_PDF, plan->fields[v_field_index[i]].ControlPoints[v_ControlPoint_index[i]].NbOfScannedSpots) + 1;
    
    Sample_particle (&New_hadrons[i], CT_Length, machine, &plan->fields[v_field_index[i]].ControlPoints[v_ControlPoint_index[i]], &plan->fields[v_field_index[i]].ControlPoints[v_ControlPoint_index[i]].spots[v_spot_index[i]], RNG_Stream);

    New_hadrons_field[i] = &plan->fields[v_field_index[i]];
    New_hadrons_layer[i] = &plan->fields[v_field_index[i]].ControlPoints[v_ControlPoint_index[i]];
    if(New_hadrons_layer[i]->RS_setting == IN && New_hadrons_layer[i]->RS_WET > 0) use_RS = 1;
  }

  // Range shifter simulation
  if(use_RS == 1){
    Transport_to_RangeShifter(New_hadrons, New_hadrons_layer, Nbr_New_hadrons);
    Simulate_RangeShifter(New_hadrons, New_hadrons_layer, New_hadrons_field, &Nbr_New_hadrons, config, machine, material, RNG_Stream);
  }

  // Convert BEV to CT reference frame, simulate setup uncertainties, and translate all particles to CT
  double position[3];
  for(i=0; i<Nbr_New_hadrons; i++){

    if(New_hadrons[i].type == Unknown) continue;

//printf("\nAfterSim: init: [%.3f ; %.3f ; %.3f]", New_hadrons[i].x, New_hadrons[i].y, New_hadrons[i].z);
//printf("\nAfterSim: direction BEV: [%.3f ; %.3f ; %.3f]", New_hadrons[i].u, New_hadrons[i].v, New_hadrons[i].w);
    BEV_to_CT_frame(&New_hadrons[i], machine, New_hadrons_field[i]);
//printf("\nAfterSim: BEV to CT: [%.3f ; %.3f ; %.3f]", New_hadrons[i].x, New_hadrons[i].y, New_hadrons[i].z);
//printf("\nAfterSim: direction CT: [%.3f ; %.3f ; %.3f]", New_hadrons[i].u, New_hadrons[i].v, New_hadrons[i].w);

    Translation_uncertainty(&New_hadrons[i], config, RNG_Stream);
//printf("\nAfterSim: translation uncertainty: [%.3f ; %.3f ; %.3f]", New_hadrons[i].x, New_hadrons[i].y, New_hadrons[i].z);

    Transport_to_CT(&New_hadrons[i], CT_Length);
//printf("\nAfterSim: transport CT: [%.3f ; %.3f ; %.3f]", New_hadrons[i].x, New_hadrons[i].y, New_hadrons[i].z);

    if(New_hadrons[i].x < 0 || New_hadrons[i].y < 0 || New_hadrons[i].z < 0 || New_hadrons[i].x > CT_Length[0] || New_hadrons[i].y > CT_Length[1] || New_hadrons[i].z > CT_Length[2] || isnan(New_hadrons[i].x) || isnan(New_hadrons[i].y) || isnan(New_hadrons[i].z)){

      #pragma omp atomic
      config->Particle_Generated_outside += 1;
/*
      printf("\nWarning: Particle generated outside the geometry: Field:%d - ControlPoint:%d - Spot:%d \n", v_field_index[i], v_ControlPoint_index[i], v_spot_index[i]);
      printf("Position: (%f;%f;%f) - Direction: (%f;%f;%f) \n", New_hadrons[i].x, New_hadrons[i].y, New_hadrons[i].z, New_hadrons[i].u, New_hadrons[i].v, New_hadrons[i].w);
	if(New_hadrons[i].x < 0) printf("x: %f < 0\n", New_hadrons[i].x);
	if(New_hadrons[i].y < 0) printf("y: %f < 0\n", New_hadrons[i].y);
	if(New_hadrons[i].z < 0) printf("z: %f < 0\n", New_hadrons[i].z);
	if(New_hadrons[i].x > CT_Length[0]) printf("x: %f > %f\n", New_hadrons[i].x, CT_Length[0]);
	if(New_hadrons[i].y > CT_Length[1]) printf("y: %f > %f\n", New_hadrons[i].y, CT_Length[1]);
	if(New_hadrons[i].z > CT_Length[2]) printf("z: %f > %f\n", New_hadrons[i].z, CT_Length[2]);
*/
    }
    else if(New_hadrons[i].type != Unknown){
      Copy_particle_buffer(&hadron[*Nbr_hadrons], &New_hadrons[i]);
      *Nbr_hadrons = *Nbr_hadrons + 1;
    }
  }
}


plan_parameters* Init_single_spot_plan(plan_parameters *Plan){

  plan_parameters *Beamlet = (plan_parameters*)malloc(sizeof(plan_parameters));

  strcpy(Beamlet->PlanName, Plan->PlanName);
  Beamlet->NumberOfFractions = Plan->NumberOfFractions;
  Beamlet->FractionID = Plan->FractionID;
  Beamlet->NumberOfFields = 1;
  Beamlet->FieldsID = (int*)malloc(sizeof(int));
  Beamlet->fields = (field_parameters*)malloc(sizeof(field_parameters));
  Beamlet->Fields_cumulative_PDF = (VAR_DATA*)malloc(sizeof(VAR_DATA));
  Beamlet->fields[0].NumberOfControlPoints = 1;
  Beamlet->fields[0].ControlPoints = (ControlPoint_parameters*)malloc(sizeof(ControlPoint_parameters));
  Beamlet->fields[0].ControlPoints_cumulative_PDF = (VAR_DATA*)malloc(sizeof(VAR_DATA));
  Beamlet->fields[0].ControlPoints[0].NbOfScannedSpots = 1;
  Beamlet->fields[0].ControlPoints[0].spots = (spot_parameters*)malloc(sizeof(spot_parameters));
  Beamlet->fields[0].ControlPoints[0].Spots_cumulative_PDF = (VAR_DATA*)malloc(sizeof(VAR_DATA));
  
  return Beamlet;
}


void Select_spot(plan_parameters *Plan, plan_parameters *Beamlet, int FieldID, int ControlPointID, int SpotID){

  Beamlet->FieldsID[0] = Plan->fields[FieldID].FieldID;
  Beamlet->TotalMetersetWeightOfAllFields = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->cumulative_weight = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->Fields_cumulative_PDF[0] = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->normalization_factor = Plan->normalization_factor;
  Beamlet->fields[0].FieldID = Plan->fields[FieldID].FieldID;
  Beamlet->fields[0].FinalCumulativeMeterSetWeight = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->fields[0].GantryAngle = Plan->fields[FieldID].GantryAngle;
  Beamlet->fields[0].PatientSupportAngle = Plan->fields[FieldID].PatientSupportAngle;
  Beamlet->fields[0].IsocenterPositionX = Plan->fields[FieldID].IsocenterPositionX;
  Beamlet->fields[0].IsocenterPositionY = Plan->fields[FieldID].IsocenterPositionY;
  Beamlet->fields[0].IsocenterPositionZ = Plan->fields[FieldID].IsocenterPositionZ;
  Beamlet->fields[0].ControlPoints_cumulative_PDF[0] = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->fields[0].RS_Type = Plan->fields[FieldID].RS_Type;
  Beamlet->fields[0].ControlPoints[0].ControlPointIndex = Plan->fields[FieldID].ControlPoints[ControlPointID].ControlPointIndex;
  Beamlet->fields[0].ControlPoints[0].SpotTunnedID = Plan->fields[FieldID].ControlPoints[ControlPointID].SpotTunnedID;
  Beamlet->fields[0].ControlPoints[0].CumulativeMetersetWeight = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->fields[0].ControlPoints[0].Energy = Plan->fields[FieldID].ControlPoints[ControlPointID].Energy;
  Beamlet->fields[0].ControlPoints[0].Spots_cumulative_PDF[0] = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->fields[0].ControlPoints[0].RS_setting = Plan->fields[FieldID].ControlPoints[ControlPointID].RS_setting;
  Beamlet->fields[0].ControlPoints[0].RS_IsocenterDist = Plan->fields[FieldID].ControlPoints[ControlPointID].RS_IsocenterDist;
  Beamlet->fields[0].ControlPoints[0].RS_WET = Plan->fields[FieldID].ControlPoints[ControlPointID].RS_WET;
  Beamlet->fields[0].ControlPoints[0].RS_Thickness = Plan->fields[FieldID].ControlPoints[ControlPointID].RS_Thickness;
  Beamlet->fields[0].ControlPoints[0].spots[0].Spot_X = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_X;
  Beamlet->fields[0].ControlPoints[0].spots[0].Spot_Y = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Y;
  Beamlet->fields[0].ControlPoints[0].spots[0].Spot_Weight = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Weight;
  Beamlet->fields[0].ControlPoints[0].spots[0].Spot_Time = Plan->fields[FieldID].ControlPoints[ControlPointID].spots[SpotID].Spot_Time;

  return;
}


plan_parameters* Select_beam(plan_parameters *Plan, int Beam){

  printf("\nGenerating sub-plan for beam %d \n\n", Beam+1);

  VAR_DATA cumulative_weight = 0;

  plan_parameters *beam_plan = (plan_parameters*)malloc(sizeof(plan_parameters));
  strcpy(beam_plan->PlanName, "Beam plan"); 
  beam_plan->NumberOfFractions = Plan->NumberOfFractions;
  beam_plan->FractionID = Plan->FractionID;
  beam_plan->NumberOfFields = 1;
  beam_plan->TotalMetersetWeightOfAllFields = 0.0;
  beam_plan->fields = (field_parameters*)malloc(1 * sizeof(field_parameters));
  beam_plan->Fields_cumulative_PDF = (VAR_DATA*)malloc(1 * sizeof(VAR_DATA));
  beam_plan->FieldsID = (int*)malloc(1 * sizeof(int));
  beam_plan->FieldsID[0] = Plan->FieldsID[Beam];
  beam_plan->fields[0].FieldID = Plan->fields[Beam].FieldID;
  beam_plan->fields[0].FinalCumulativeMeterSetWeight = 0.0;
  beam_plan->fields[0].GantryAngle = Plan->fields[Beam].GantryAngle;
  beam_plan->fields[0].PatientSupportAngle = Plan->fields[Beam].PatientSupportAngle;
  beam_plan->fields[0].IsocenterPositionX = Plan->fields[Beam].IsocenterPositionX;
  beam_plan->fields[0].IsocenterPositionY = Plan->fields[Beam].IsocenterPositionY;
  beam_plan->fields[0].IsocenterPositionZ = Plan->fields[Beam].IsocenterPositionZ;
  beam_plan->fields[0].NumberOfControlPoints = Plan->fields[Beam].NumberOfControlPoints;
  beam_plan->fields[0].RS_Type = Plan->fields[Beam].RS_Type;
  beam_plan->fields[0].ControlPoints = (ControlPoint_parameters*)malloc(Plan->fields[Beam].NumberOfControlPoints * sizeof(ControlPoint_parameters));
  beam_plan->fields[0].ControlPoints_cumulative_PDF = (VAR_DATA*)malloc(Plan->fields[Beam].NumberOfControlPoints * sizeof(VAR_DATA));

  int j,k;
  for(j=0; j<Plan->fields[Beam].NumberOfControlPoints; j++){
    beam_plan->fields[0].ControlPoints[j].ControlPointIndex = Plan->fields[Beam].ControlPoints[j].ControlPointIndex;
    beam_plan->fields[0].ControlPoints[j].SpotTunnedID = Plan->fields[Beam].ControlPoints[j].SpotTunnedID;
    beam_plan->fields[0].ControlPoints[j].CumulativeMetersetWeight = 0.0;
    beam_plan->fields[0].ControlPoints[j].Energy = Plan->fields[Beam].ControlPoints[j].Energy;
    beam_plan->fields[0].ControlPoints[j].NbOfScannedSpots = Plan->fields[Beam].ControlPoints[j].NbOfScannedSpots;
    beam_plan->fields[0].ControlPoints[j].RS_setting = Plan->fields[Beam].ControlPoints[j].RS_setting;
    beam_plan->fields[0].ControlPoints[j].RS_IsocenterDist = Plan->fields[Beam].ControlPoints[j].RS_IsocenterDist;
    beam_plan->fields[0].ControlPoints[j].RS_WET = Plan->fields[Beam].ControlPoints[j].RS_WET;
    beam_plan->fields[0].ControlPoints[j].RS_Thickness = Plan->fields[Beam].ControlPoints[j].RS_Thickness;
    beam_plan->fields[0].ControlPoints[j].spots = (spot_parameters*)malloc(Plan->fields[Beam].ControlPoints[j].NbOfScannedSpots * sizeof(spot_parameters));
    beam_plan->fields[0].ControlPoints[j].Spots_cumulative_PDF = (VAR_DATA*)malloc(Plan->fields[Beam].ControlPoints[j].NbOfScannedSpots * sizeof(VAR_DATA));

    for(k=0; k<Plan->fields[Beam].ControlPoints[j].NbOfScannedSpots; k++){
      beam_plan->fields[0].ControlPoints[j].spots[k].Spot_X = Plan->fields[Beam].ControlPoints[j].spots[k].Spot_X;
      beam_plan->fields[0].ControlPoints[j].spots[k].Spot_Y = Plan->fields[Beam].ControlPoints[j].spots[k].Spot_Y;
      beam_plan->fields[0].ControlPoints[j].spots[k].Spot_Time = Plan->fields[Beam].ControlPoints[j].spots[k].Spot_Time;
      beam_plan->fields[0].ControlPoints[j].spots[k].Spot_Weight = Plan->fields[Beam].ControlPoints[j].spots[k].Spot_Weight;

      cumulative_weight += beam_plan->fields[0].ControlPoints[j].spots[k].Spot_Weight;
      beam_plan->fields[0].ControlPoints[j].Spots_cumulative_PDF[k] = cumulative_weight;
    }
    beam_plan->fields[0].ControlPoints_cumulative_PDF[j] = cumulative_weight;
    beam_plan->fields[0].ControlPoints[j].CumulativeMetersetWeight = cumulative_weight;
  }

  beam_plan->Fields_cumulative_PDF[0] = cumulative_weight;
  beam_plan->fields[0].FinalCumulativeMeterSetWeight = cumulative_weight;
  beam_plan->TotalMetersetWeightOfAllFields = cumulative_weight;
  beam_plan->cumulative_weight = cumulative_weight;

  beam_plan->normalization_factor = Plan->normalization_factor * beam_plan->Fields_cumulative_PDF[0] / Plan->Fields_cumulative_PDF[Plan->NumberOfFields-1];

  return beam_plan;

}


