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


inline double matrix_product(double a[3][3], double b[3][3], double product[3][3])
{  
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
       product[i][j]=a[i][j]*b[i][j];
    }
  } 
}

/**multiply vectors (3x1)*(1x3) */
inline double multiply_3x1_1x3(double v1[3], double v2[3], double res[3][3])
{
 for (int i=0;i<3;i++)
 { 
   for (int j=0;j<3;j++)
   {
    res[i][j]=v1[i]*v2[j];
   } 
 }
}

/**calculate the inverse of matrix 3x3 */
inline double matrix_3x3_inverse(double matrix[3][3], double inverse[3][3])
{
   double a11=matrix[0][0];  
   double a12=matrix[0][1];  
   double a13=matrix[0][2];
   double a21=matrix[1][0];  
   double a22=matrix[1][1];  
   double a23=matrix[1][2];
   double a31=matrix[2][0];  
   double a32=matrix[2][1];  
   double a33=matrix[2][2];
  
   double determinant=a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a32*a23*a11 - a33*a21*a12;
   if (!(determinant==0)) 
     determinant=1/determinant;

   inverse[0][0]=determinant*a22*a33-a32*a23;
   inverse[0][1]=determinant*a13*a32-a33*a12;
   inverse[0][2]=determinant*a12*a23-a22*a13;
   inverse[1][0]=determinant*a23*a31-a33*a21; 
   inverse[1][1]=determinant*a11*a32-a33*a12;
   inverse[1][2]=determinant*a13*a21-a23*a11;
   inverse[2][0]=determinant*a21*a32-a31*a22;
   inverse[2][1]=determinant*a12*a31-a32*a11;
   inverse[2][2]=determinant*a11*a22-a12*a21;
}

inline double determinant_of_3x3_matrix(double matrix[3][3], double determinant)
{
  determinant = matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[2][1]*matrix[1][2])-
                matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[2][0]*matrix[1][2])+
                matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[2][0]*matrix[1][1]); 
}

/** orthogonal rotation matrix (ORM) is here filled in with random numbers, however its rows should represent current semiaxes of ellipsoid, which means that they have to be recalculated in each integration step*/
inline double get_mi_random_oriented_ellipsoid(int control_parameter, double (&ORM)[3][3])
{
   for (int i=0;i<3;i++){
     for (int j=0;j<3;j++){
       ORM[i][j]=(d_random()*2-1)+i+j;
     }
   }
/*

//ellipsoid 'a' projections of its semiaxes in x,y and z direction
double semiaxa1[3]={siga1_x,siga1_y,siga1_z};
double semiaxa2[3]={siga2_x,siga2_y,siga2_z};
double semiaxa3[3]={siga3_x,siga3_y,siga3_z};

//ellipsoid 'b' projections of its semiaxes in x,y and z direction
double semiaxb1[3]={sigb1_x,sigb1_y,sigb1_z};
double semiaxb2[3]={sigb2_x,sigb2_y,sigb2_z};
double semiaxb3[3]={sigb3_x,sigb3_y,sigb3_z};
*/
}


inline double solve_quadratic_form(double matrix[3][3], double vector[3], double result) 
{
  double transposed_vector_dot_matrix[3];
  transposed_vector_dot_matrix[0] = vector[0]*matrix[0][0] + vector[1]*matrix[1][0] + vector[2]*matrix[2][0];
  transposed_vector_dot_matrix[1] = vector[0]*matrix[0][1] + vector[1]*matrix[1][1] + vector[2]*matrix[2][1];
  transposed_vector_dot_matrix[2] = vector[0]*matrix[0][2] + vector[1]*matrix[1][2] + vector[2]*matrix[2][2];
  
  result = transposed_vector_dot_matrix[0]*vector[0] + transposed_vector_dot_matrix[1]*vector[1] + transposed_vector_dot_matrix[2]*vector[2];
}

inline double multiply_matrix_3x3_by_vector_3x1(double matrix[3][3], double vector[3], double result[3]) 
{
  result[0]=matrix[0][0]*vector[0]+matrix[0][1]*vector[1]+matrix[0][2]*vector[2];
  result[1]=matrix[1][0]*vector[0]+matrix[1][1]*vector[1]+matrix[1][2]*vector[2];
  result[0]=matrix[2][0]*vector[0]+matrix[2][1]*vector[1]+matrix[2][2]*vector[2];
}

inline double multiply_vector_1x3_by_matrix_3x3(double matrix[3][3], double vector[3], double result[3]) 
{
  result[0]=vector[0]*matrix[0][0]+vector[1]*matrix[1][0]+vector[2]*matrix[2][0];
  result[1]=vector[0]*matrix[0][1]+vector[1]*matrix[1][1]+vector[2]*matrix[2][1];
  result[2]=vector[0]*matrix[0][2]+vector[1]*matrix[1][2]+vector[2]*matrix[2][2];
}


//* calculate Gay-Berne potential derivatives*/
inline void add_gb_pair_force(const Particle * const p1, const Particle * const p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3], 
                              double torque1[3], double torque2[3])
{  
  if (!CUTOFF_CHECK(dist < ia_params->GB_cut))   
    return;
  
// evaluate current particle orientation
 double orient_1[3], orient_2[3];
 memcpy(orient_1, p1->r.quatu,  3*sizeof(double));
 memcpy(orient_2, p2->r.quatu,  3*sizeof(double));

/*
 orient_1[0]=p1->r.quatu[0];
 orient_1[1]=p1->r.quatu[1];
 orient_1[2]=p1->r.quatu[2];
 orient_2[0]=p2->r.quatu[0];
 orient_2[1]=p2->r.quatu[1];
 orient_2[2]=p2->r.quatu[2];
*/
// evaluate current particle position
 double pos_1[3], pos_2[3];
 pos_1[0]=p1->r.p[0];
 pos_1[1]=p1->r.p[1];
 pos_1[2]=p1->r.p[2];
 pos_2[0]=p2->r.p[0];
 pos_2[1]=p2->r.p[1];
 pos_2[2]=p2->r.p[2];
 printf("particle at the position %f %f %f\n", pos_1[0], pos_1[1], pos_1[2]);
 printf("particle at the position %f %f %f\n", p2->r.p[0],p2->r.p[1],p2->r.p[2]);
 
 double r_vector[3];
 vecsub(pos_1, pos_2, r_vector);
 double r_moduo;
 r_moduo=normr(r_vector);
 double r_unit[3];
 unit_vector(r_vector, r_unit);

 double sig1=ia_params->GB_sig1;
 double sig2=ia_params->GB_sig2;
 double sig3=ia_params->GB_sig3;

 double eps1=ia_params->GB_eps1;
 double eps2=ia_params->GB_eps2;
 double eps3=ia_params->GB_eps3;

 double mu=ia_params->GB_mu;
 double nu=ia_params->GB_nu;
// implementation according Michael P. Allen - EXPRESSIONS FOR FORCES AND TORQUES IN MOLECULAR SIMULATIONS USING RIGID BODIES 
 double sigma_min=dmin(dmin(sig1,sig2),sig3);
 double sigma[3][3];
 
// identity matrix
 double I[3][3];
 for (int i=0;i<3;i++){
   for (int j=0; j<3;j++) {
     if (i==j) I[i][j]=1;
     else I[i][j]=0;
   }
 }

// orientation matrices a_hat, b_hat
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

//shape matrix S1, S2 - matrices 3x3 with ellipsoids semiradii on the main diagonals
double S1[3][3]={{sig1,0,0},{0,sig2,0},{0,0,sig3}};
double S2[3][3]={{sig1,0,0},{0,sig2,0},{0,0,sig3}};
double a_hat_inverse[3][3];
double b_hat_inverse[3][3];

matrix_3x3_inverse(a_hat, a_hat_inverse);
matrix_3x3_inverse(b_hat, b_hat_inverse);

double A[3][3], B[3][3], A_1[3][3], B_1[3][3], A_SS[3][3], B_SS[3][3];
matrix_product(S1,S1,A_SS);
matrix_product(S2,S2,B_SS);
matrix_product(a_hat_inverse,A_SS,A_1);
matrix_product(b_hat_inverse,B_SS,B_1);
matrix_product(A_1,a_hat,A);
matrix_product(B_1,b_hat,B);

printf("matrix A=[%f %f %f %f %f %f %f %f %f\n]", A[0,0], A[0,1], A[0,2], A[1,0], A[1,1], A[1,2], A[2,0], A[2,1], A[2,2]);
printf("matrix B=[%f %f %f %f %f %f %f %f %f\n]", B[0,0], B[0,1], B[0,2], B[1,0], B[1,1], B[1,2], B[2,0], B[2,1], B[2,2]);
double H[3][3];
for (int i=0;i<3;i++){
  for (int j=0;j<3;j++){
    H[i][j]=A[i][j]+B[i][j];
  }
}

double epsilon0=0;
double eps_AB;
double determinant;
eps_AB=pow((2*sigma_min*sigma_min/determinant_of_3x3_matrix(H,determinant)),0.5);
double eps_AB_nu;
eps_AB_nu=pow(eps_AB,nu);

/* matrix of inversial of potential depths */
double B_AB[3][3];
for (int i=0;i<3;i++){
 for (int j=0;j<3;j++){
   if (!(i==j)) 
   {
    B_AB[i][j]=0;
   }
 }
}  
B_AB[0][0]=pow((epsilon0/eps1),1./mu);
B_AB[1][1]=pow((epsilon0/eps2),1./mu);
B_AB[2][2]=pow((epsilon0/eps3),1./mu);

double inv_B_AB[3][3];
matrix_3x3_inverse(B_AB,inv_B_AB);
double eps_AB_prim;
solve_quadratic_form(inv_B_AB, r_vector, eps_AB_prim);
double eps_Ab_prim_mu=pow(eps_AB_prim, mu);

//derivative with respect to r
double d_eps_AB_prim_dr[3];
d_eps_AB_prim_dr[0]=4*(inv_B_AB[0][0]*r_vector[0]+inv_B_AB[0][1]*r_vector[1]+inv_B_AB[0][2]*r_vector[2]);
d_eps_AB_prim_dr[1]=4*(inv_B_AB[1][0]*r_vector[0]+inv_B_AB[1][1]*r_vector[1]+inv_B_AB[1][2]*r_vector[2]);
d_eps_AB_prim_dr[2]=4*(inv_B_AB[2][0]*r_vector[0]+inv_B_AB[2][1]*r_vector[1]+inv_B_AB[2][2]*r_vector[2]);

//sigma_AB
double sigma_AB;
double inv_H[3][3];
matrix_3x3_inverse(H, inv_H);
double intermediate;
solve_quadratic_form(inv_H,r_unit,intermediate);
sigma_AB=pow(2*intermediate,-0.5);

//derivative with respect to r 
double temporarly;
solve_quadratic_form(inv_H, r_unit,temporarly);
double matrvec[3];
multiply_matrix_3x3_by_vector_3x1(inv_H, r_unit, matrvec);
double temp=-0.5*pow(temporarly,-3./2)*4;
double d_sigma_AB_dr[3];
d_sigma_AB_dr[0]=temp*matrvec[0];
d_sigma_AB_dr[1]=temp*matrvec[1];
d_sigma_AB_dr[2]=temp*matrvec[2];

double prederivative=4*epsilon0*eps_AB_nu;
//aplied chain rule --> left and right derivative
double left_derivative; double right_derivative;
double some[3]; 
multiply_matrix_3x3_by_vector_3x1(inv_B_AB, r_vector, some);
double somesth[3];
multiply_matrix_3x3_by_vector_3x1(inv_H, r_unit, somesth);

//sigma_c/(r_12-sigma_AB+sigma_c)
double q=sigma_min/(r_moduo-sigma_AB+sigma_min);

left_derivative=mu*pow(eps_AB_prim,(mu-1))*4*(pow(q,12)-pow(q,6));
double left[3];
for (int i=0;i<3;i++)
  {
  left[i]=left_derivative*some[i];
  }
right_derivative=pow(eps_AB_prim,mu)*6*(-4*sigma_min)*(2*pow(q,11)*(1+0.5*pow(sigma_AB,3))-pow(q,5)*1+0.5*pow(sigma_AB,3));
double right[3];
for (int i=0;i<3;i++)
  {
  right[i]=right_derivative*somesth[i];
  }

double dU_dr[3];
for (int i=0;i<3;i++)
  {
  dU_dr[i]=prederivative*(left[i]-right[i]);
  }
printf("potential derivative is: %f %f %f\n", dU_dr[0], dU_dr[1],dU_dr[2]);

}


inline double gb_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
			       double d[3], double dist)
{
  if (!CUTOFF_CHECK(dist < ia_params->GB_cut))   
    return 0.0;

// evaluate current particle orientation
 double orient_1[3], orient_2[3];
 memcpy(orient_1, p1->r.quatu,  3*sizeof(double));
 memcpy(orient_2, p2->r.quatu,  3*sizeof(double));

// evaluate current particle position
 double pos_1[3], pos_2[3];
 pos_1[0]=p1->r.p[0];
 pos_1[1]=p1->r.p[1];
 pos_1[2]=p1->r.p[2];
 pos_2[0]=p2->r.p[0];
 pos_2[1]=p2->r.p[1];
 pos_2[2]=p2->r.p[2];
 printf("particle at the position %f %f %f\n", pos_1[0], pos_1[1], pos_1[2]);
 printf("particle at the position %f %f %f\n", p2->r.p[0],p2->r.p[1],p2->r.p[2]);
 
 double r_vector[3];
 vecsub(pos_1, pos_2, r_vector);
 double r_moduo;
 r_moduo=normr(r_vector);
 double r_unit[3];
 unit_vector(r_vector, r_unit);

 double sig1=ia_params->GB_sig1;
 double sig2=ia_params->GB_sig2;
 double sig3=ia_params->GB_sig3;

 double eps1=ia_params->GB_eps1;
 double eps2=ia_params->GB_eps2;
 double eps3=ia_params->GB_eps3;

 double mu=ia_params->GB_mu;
 double nu=ia_params->GB_nu;
// implementation according Michael P. Allen - EXPRESSIONS FOR FORCES AND TORQUES IN MOLECULAR SIMULATIONS USING RIGID BODIES 
 double sigma_min=dmin(dmin(sig1,sig2),sig3);
 double sigma[3][3];
 
// identity matrix
 double I[3][3];
 for (int i=0;i<3;i++){
   for (int j=0; j<3;j++) {
     if (i==j) I[i][j]=1;
     else I[i][j]=0;
   }
 }

// orientation matrices a_hat, b_hat
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

//shape matrix S1, S2 - matrices 3x3 with ellipsoids semiradii on the main diagonals
double S1[3][3]={{sig1,0,0},{0,sig2,0},{0,0,sig3}};
double S2[3][3]={{sig1,0,0},{0,sig2,0},{0,0,sig3}};
double a_hat_inverse[3][3];
double b_hat_inverse[3][3];

matrix_3x3_inverse(a_hat, a_hat_inverse);
matrix_3x3_inverse(b_hat, b_hat_inverse);

double A[3][3], B[3][3], A_1[3][3], B_1[3][3], A_SS[3][3], B_SS[3][3];
matrix_product(S1,S1,A_SS);
matrix_product(S2,S2,B_SS);
matrix_product(a_hat_inverse,A_SS,A_1);
matrix_product(b_hat_inverse,B_SS,B_1);
matrix_product(A_1,a_hat,A);
matrix_product(B_1,b_hat,B);

double H[3][3];
for (int i=0;i<3;i++){
  for (int j=0;j<3;j++){
    H[i][j]=A[i][j]+B[i][j];
  }
}

double epsilon0=0;
double eps_AB;
double determinant;
eps_AB=pow((2*sigma_min*sigma_min/determinant_of_3x3_matrix(H,determinant)),0.5);
double eps_AB_nu;
eps_AB_nu=pow(eps_AB,nu);

/* matrix of inversial of potential depths */
double B_AB[3][3];
for (int i=0;i<3;i++){
 for (int j=0;j<3;j++){
   if (!(i==j)) 
   {
    B_AB[i][j]=0;
   }
 }
}  
B_AB[0][0]=pow((epsilon0/eps1),1./mu);
B_AB[1][1]=pow((epsilon0/eps2),1./mu);
B_AB[2][2]=pow((epsilon0/eps3),1./mu);

double inv_B_AB[3][3];
matrix_3x3_inverse(B_AB,inv_B_AB);
double eps_AB_prim;
solve_quadratic_form(inv_B_AB, r_vector, eps_AB_prim);
double eps_Ab_prim_mu=pow(eps_AB_prim, mu);

//derivative with respect to r
double d_eps_AB_prim_dr[3];
d_eps_AB_prim_dr[0]=4*(inv_B_AB[0][0]*r_vector[0]+inv_B_AB[0][1]*r_vector[1]+inv_B_AB[0][2]*r_vector[2]);
d_eps_AB_prim_dr[1]=4*(inv_B_AB[1][0]*r_vector[0]+inv_B_AB[1][1]*r_vector[1]+inv_B_AB[1][2]*r_vector[2]);
d_eps_AB_prim_dr[2]=4*(inv_B_AB[2][0]*r_vector[0]+inv_B_AB[2][1]*r_vector[1]+inv_B_AB[2][2]*r_vector[2]);

//sigma_AB
double sigma_AB;
double inv_H[3][3];
matrix_3x3_inverse(H, inv_H);
printf("matrix H=[ %f %f %f %f %f %f %f %f %f\n]",H[0,0],H[0,1],H[0,2],H[1,0],H[1,1],H[1,2],H[2,0],H[2,1],H[2,2]);
printf("matrix inv_H=[ %f %f %f %f %f %f %f %f %f\n]",inv_H[0,0],inv_H[0,1],inv_H[0,2],inv_H[1,0],inv_H[1,1],inv_H[1,2],inv_H[2,0],inv_H[2,1],inv_H[2,2]);
double intermediate;
solve_quadratic_form(inv_H,r_unit,intermediate);
sigma_AB=pow(2*intermediate,-0.5);

//derivative with respect to r 
double temporarly;
solve_quadratic_form(inv_H, r_unit,temporarly);
double matrvec[3];
multiply_matrix_3x3_by_vector_3x1(inv_H, r_unit, matrvec);
double temp=-0.5*pow(temporarly,-3./2)*4;
double d_sigma_AB_dr[3];
d_sigma_AB_dr[0]=temp*matrvec[0];
d_sigma_AB_dr[1]=temp*matrvec[1];
d_sigma_AB_dr[2]=temp*matrvec[2];

double prederivative=4*epsilon0*eps_AB_nu;
//aplied chain rule --> left and right derivative
double left_derivative; double right_derivative;
double some[3]; 
multiply_matrix_3x3_by_vector_3x1(inv_B_AB, r_vector, some);
double somesth[3];
multiply_matrix_3x3_by_vector_3x1(inv_H, r_unit, somesth);

//sigma_c/(r_12-sigma_AB+sigma_c)
double q=sigma_min/(r_moduo-sigma_AB+sigma_min);

printf("r moduo= %f\n", r_moduo);
printf("sigma AB= %f\n", sigma_AB);
printf("sigma min= %f\n", sigma_min);
printf("q= %f\n", q);
left_derivative=mu*pow(eps_AB_prim,(mu-1))*4*(pow(q,12)-pow(q,6));
printf("left derivative= %f\n", left_derivative);
double left[3];
for (int i=0;i<3;i++)
  {
  left[i]=left_derivative*some[i];
  }
right_derivative=pow(eps_AB_prim,mu)*6*(-4*sigma_min)*(2*pow(q,11)*(1+0.5*pow(sigma_AB,3))-pow(q,5)*1+0.5*pow(sigma_AB,3));
//printf("right derivative is %f\n", right_derivative);
double right[3];
for (int i=0;i<3;i++)
  {
  right[i]=right_derivative*somesth[i];
  }

double dU_dr[3];
for (int i=0;i<3;i++)
  {
  dU_dr[i]=prederivative*(left[i]-right[i]);
  }
//until here everything is the same as in force derivative calculation

//* auxilary vector k=inv(H)\scalar r  */

double invH[3][3];
matrix_3x3_inverse(H, invH);
double k[3];
multiply_matrix_3x3_by_vector_3x1(invH, r_vector, k);

//as force=-d_U/d_fi \times k, -d_U/d_fi=force k (-1)
double f_k[3][3];
multiply_3x1_1x3(dU_dr,r_vector,f_k); //here f_k became a 3x3 array

//here define torque tau_A=-k \times A * force
double kA[3];
multiply_vector_1x3_by_matrix_3x3(A,k,kA);
double tau_A[3];
multiply_vector_1x3_by_matrix_3x3(f_k,kA,tau_A);


double kB[3];
multiply_vector_1x3_by_matrix_3x3(B,k,kB);
double tau_B[3];
multiply_vector_1x3_by_matrix_3x3(f_k,kB,tau_B);

printf("Increment of the Gay-Berne energy is: %f %f %f\n", tau_A[0], tau_A[1], tau_A[2]);
}

#endif
#endif
