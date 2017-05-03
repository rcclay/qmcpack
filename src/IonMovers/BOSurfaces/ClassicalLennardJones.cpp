//
//  classical Lennard Jones potential and its force
//  Form of potential U(R)= epsilon*((sigma/R)^12-(sigma/R)^6)
//
//  Created by Dominik Domin on 4/27/17.

#include "classicalLJ.h"


double LennardJones::lennard_jones_energy(double * sqrt_epsilon , double * half_sigma, double * Rnuc, int nnuclei,int ndim)
{
    int i, j, k;
    double energy = 0.0;
    double epsilon = 0.0;
    double sigma = 0.0;
    double d = 0.0;
    
    for (i=0; i<nnuclei; i++) {
        for (j=0; j<i; j++) {
            sigma =half_sigma[i]+half_sigma[j];
            epsilon = sqrt_epsilon[i]*sqrt_epsilon[j];
            d = 0.0;
            for (k=0; k<ndim; k++) {
                d+=pow((Rnuc[(ndim*i)+k]-Rnuc[(ndim*j)+k]), 2.0);
            }
            d = sqrt(d);
            energy += epsilon*(pow(sigma/d,12.0)-pow(sigma/d,6.0));
        }
    }
    return energy;
}

void LennardJones::lennard_jones_force(double * sqrt_epsilon , double * half_sigma, double * Rnuc, int nnuclei,int ndim, double * force_array){
    double f[ndim];
    int i, j, k;
    double epsilon = 0.0;
    double sigma = 0.0;
    double d = 0.0;
    
    for (i=0; i<nnuclei; i++) {
        for (k=0; k<ndim; k++){
            f[k]=0.0;
            force_array[(ndim*i)+k]=0.0;
        }
        for (j=0; j<nnuclei; j++) {
            if (i!=j){
                sigma =half_sigma[i]+half_sigma[j];
                epsilon = sqrt_epsilon[i]*sqrt_epsilon[j];
                d = 0.0;
                for (k=0; k<ndim; k++){
                    d+=pow((Rnuc[(ndim*i)+k]-Rnuc[(ndim*j)+k]), 2.0);
                }
                d=sqrt(d);
                for (k=0; k<ndim; k++){
                    f[k]-=(1.0-2.0*pow(sigma/d,6.0))*24.0*epsilon*pow(sigma,6.0)*(Rnuc[(ndim*i)+k]-Rnuc[(ndim*j)+k])/pow(d,8.0);
                }
            }
        }
        for (k=0;k<ndim ; k++) {
            force_array[(ndim*i)+k]=f[k];
        }
    }
}

bool LennardJones::E(const ParticlePos_t& R1, double& dE, double& err){
//d= ParticleSet R1.distance(i,j);
//dummy functions
return 0;
}

bool LennardJones::dE(const ParticlePos_t& R0, const ParticlePos_t& R1, double& dE, double& err){
//dummy functions
return 1;
}

bool LennardJones::dE(const ParticlePos_t& R0, const ParticlePos_t& R1, int iat, double& dE, double& err){
//dummy functions
return 1;
}

bool LennardJones::F(const ParticlePos_t& R1, ParticlePos_t& F){
//dummy functions
return 1;
}
