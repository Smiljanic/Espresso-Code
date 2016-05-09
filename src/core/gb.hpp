/*
  Copyright (C) 2010,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _GB_HPP
#define _GB_HPP

/** \file gb.hpp
 *  Routines to calculate the Gay-Berne energy and force 
 *  for a pair of particles.
 *  \ref forces.cpp
*/

#include "utils.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "mol_cut.hpp"
#include "grid.hpp"

#ifdef GAY_BERNE

//zuerst nur mit aehnlichen Ellipsoide
int gay_berne_set_params(int part_type_a, int part_type_b,
			 double eps1, double eps2, double eps3, 
                         double sig1, double sig2, double sig3, 
                         double cut,
			 double mu, double nu);


//double connecting_vector[3];
//void inline calculate_vector(double *vector, double *a, double *b)
//{
//  vector[0]=b[0]-a[0];
//  vector[1]=b[1]-a[1];
//  vector[2]=b[2]-a[2];
//  return;
//}

//particle orientation/quat is taken with the particle
inline void add_gb_pair_force(const Particle * const p1, const Particle * const p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3], 
                              double torque1[3], double torque2[3])

{  
  if (!CUTOFF_CHECK(dist < ia_params->GB_cut))   
    return;
  
// evaluate current particle orientation
 double orient_1[3], orient_2[3];
 orient_1[0]=p1->r.quatu[0];
 orient_1[1]=p1->r.quatu[1];
 orient_1[2]=p1->r.quatu[2];
 orient_2[0]=p2->r.quatu[0];
 orient_2[1]=p2->r.quatu[1];
 orient_2[2]=p2->r.quatu[2];

// evaluate current particle position
 double pos_1[3], pos_2[3];
 pos_1[0]=p1->r.p[0];
 pos_1[1]=p1->r.p[1];
 pos_1[2]=p1->r.p[2];
 pos_2[0]=p2->r.p[0];
 pos_2[1]=p2->r.p[1];
 pos_2[2]=p2->r.p[2];
 printf("particle 1 at the position %f %f %f\n", pos_1[0], pos_1[1], pos_1[2]);
 
 double r_vector[3];
 vecsub(pos_1, pos_2, r_vector);
 double r_moduo;
 r_moduo=normr(r_vector);
 double r_unit[3];
 unit_vector(r_vector, r_unit);

 double sig1=ia_params->GB_sig1;
 double sig2=ia_params->GB_sig2;
 double sig3=ia_params->GB_sig3;

 //implementation according Michael P. Allen - EXPRESSIONS FOR FORCES AND TORQUES IN MOLECULAR SIMULATIONS USING RIGID BODIES 
 //double sigma_min=dmin(dmin(ia_params->GB_sig1,ia_params->GB_sig2),ia_params->GB_sig3);
 double sigma_min=dmin(dmin(sig1,sig2),sig3);
 double sigma[3][3];
 
//identity matrix
 double I[3][3];
 for (int i=0;i<3;i++){
   for (int j=0; j<3;j++) {
     if (i==j) I[i][j]=1;
     else I[i][j]=0;
   }
 }

//orientation matrices a_hat, b_hat


// double a_hat[3][3], b_hat[3][3];
 double a_hat[3][3]={{0,0,0},{0,0,0},{0,0,0}};
 double b_hat[3][3]={{0,0,0},{0,0,0},{0,0,0}};
 for (int i=0; i<3; i++){
   for (int j=0; j<3; j++){
     if (i==j) { 
       a_hat[i][j]=orient_1[i];
       b_hat[i][j]=orient_1[i];
     }
   }
 }

//shape matrix S1, S2
double S1[3][3]={{sig1,0,0},{0,sig2,0},{0,0,sig3}};

//A, B
 


 
// a_hat=scalar(orient_1, I);    check utils.hpp --> scalar(a[3][3],b[3][3])
// b_hat=scalar(orient_2, I); 

// double qu=(r_moduo-sigma+sigma_min)/sigma_min;
// double dU_dfi = 24*ia_params->GB_eps;

//for (int i = 0; i<3; i++)
//   for (int j = 0; j<3; j++)
//     result[j][i] = matrix[i][j];

}

inline double gb_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
			       double d[3], double dist)
{
  if (!CUTOFF_CHECK(dist < ia_params->GB_cut))   
    return 0.0;
  
  double a,b,c, X, Xcut,
    Brack,BrackCut,
    u1x, u1y, u1z,
    u2x, u2y, u2z,	  
    E,E1,E2, Sigma,
    Plus1, Minus1,
    Plus2, Minus2;
	
    
    u1x = p1->r.quatu[0]; u1y = p1->r.quatu[1]; u1z = p1->r.quatu[2];
    u2x = p2->r.quatu[0]; u2y = p2->r.quatu[1]; u2z = p2->r.quatu[2]; 

    a = d[0]*u1x + d[1]*u1y + d[2]*u1z;
    b = d[0]*u2x + d[1]*u2y + d[2]*u2z;
    c =  u1x*u2x +  u1y*u2y +  u1z*u2z;

    Plus1 = (a+b)/(1+ia_params->GB_chi1*c);
    Plus2 = (a+b)/(1+ia_params->GB_chi2*c);
    Minus1 = (a-b)/(1-ia_params->GB_chi1*c);
    Minus2 = (a-b)/(1-ia_params->GB_chi2*c);
    E1 = 1/sqrt(1-ia_params->GB_chi1*ia_params->GB_chi1*c*c);
    E2 = 1-0.5*(ia_params->GB_chi2/dist/dist)*(Plus2*(a+b) + Minus2*(a-b));
    E = 4*ia_params->GB_eps1*pow(E1,ia_params->GB_nu)*pow(E2,ia_params->GB_mu);  
    Sigma = ia_params->GB_sig1/sqrt(1-0.5*(ia_params->GB_chi1/dist/dist)*(Plus1*(a+b) + Minus1*(a-b)));
        
    X = 1/(dist - Sigma + ia_params->GB_sig1);
    Xcut = 1/(ia_params->GB_cut - Sigma + ia_params->GB_sig1);

    Brack = X*X*X;
    BrackCut = Xcut*Xcut*Xcut;
    Brack = Brack*Brack;
    BrackCut = BrackCut*BrackCut;
    Brack = Brack*(Brack-1);
    BrackCut = BrackCut*(BrackCut-1);

    return E*(Brack-BrackCut);
}

#endif
#endif
