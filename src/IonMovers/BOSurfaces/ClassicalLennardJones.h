//
//  classical Lennard Jones potential
//
//  Created by Dominik Domin on 4/27/17.
//

#ifndef __classicalLJ__
#define __classicalLJ__

//#include <iostream>
//#include <stdio.h>
#include <stdlib.h>
#include <math.h>
class LennardJones : public BOSurfaceBase {

private:
   int ndim;
   /* dimensionality of the simulation (e.g. 1D, 2D, 3D) */

   int nnuclei;
   /* number of nuclei (atoms) */

   double * sqrt_epsilon;
   /* array of square root of epsilon of the varius nuclei since private functions assume epsilon_ij = sqrt(epsilon_ii*epsilon_jj) */
   
   double * half_sigma;
   /* half of sigma values since functions assume sigma_ij = (sigma_ii+sigma_jj)/2  */
   /* ordering of the array is the same as for the Rnuc array and sqrt_epsilon */

   double * force_array;
   /* Nnuclei*NDIM-dimensional array containing the forces on all the nuclei */

   double * Rnuc;
   /* Nnuclei*NDIM-dimensional array containing the positions of the nuclei in NDIM dimensions */

   double lennard_jones_energy(double * sqrt_epsilon, double * half_sigma, double * Rnuc, int nnuclei,int ndim); 
   /* computes the potential energy based on the Lennard Jones potential */

   void lennard_jones_forces(double * sqrt_epsilon, double * half_sigma, double * Rnuc, int nnuclei, int ndim, double * force_array);
   /* computes the forces based on the Lennard Jones potential */

public:
   LennardJones(double * epsilon, double * sigma, int nnuclei, int ndim);
   /* constructor to initialize class */

   ~LennardJones();
   /* destructor to free memory */

//*************************************************************************************/
   //DO I NEED TO PLACE declaration of the overloaded BOSurfaceBase functions HERE?

}
#endif /* defined(__classicalLJ__) */
