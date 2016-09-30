/*****************************************************************************/
/*                                                                           */
/*                   Copyright 08/2006 by Dr. Andres Bruhn,                  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef OF_HORN_SCHUNCK_INCLUDED
#define OF_HORN_SCHUNCK_INCLUDED

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc_mem_linear.c"
#include "alloc_mem_linear_mult.c"
#include "funct_lib.c"


/* ------------------------------------------------------------------------- */
/*
 Computes one Jacobi iteration
*/
void horn_schunck_jacobi
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **u,              /* in+out : x-component of displacement field        */
float **v,              /* in+out : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float alpha             /* in     : smoothness weight                        */
                        /*****************************************************/
)



{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float hx_2,hy_2;        /* time saver variables                              */
float xp,xm,yp,ym;      /* neighbourhood weights                             */
float xpp, xmm, ypp, ymm;
float sum;              /* central weight                                    */      
float **u_old;          /* x-component from old time step                    */
float **v_old;          /* y-component from old time step                    */
                        /*****************************************************/
float omega = 1.96;
int choice;

/* define time saver variables */
hx_2=alpha/(hx*hx);
hy_2=alpha/(hy*hy);

/* alloc memory */
ALLOC_MATRIX(2, nx+2*bx, ny+2*by, 
	     &u_old,
	     &v_old);

/* copy previous result */
copy_matrix_2d(u,u_old,nx,ny,bx,by);
copy_matrix_2d(v,v_old,nx,ny,bx,by);

/* set boundaries zero */
set_bounds_2d(u_old,nx,ny,bx,by,0.0);
set_bounds_2d(v_old,nx,ny,bx,by,0.0);

choice = 2;
/* Jacobi iteration */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute weights */
	
        /* TODO
          ----- fill in your code for computing the neighbourhood weights
                here - do not forget the correct boundary conditions  ----    
	    */
	    xp = (i <= nx+bx-2)*(alpha/(hx*hx));
	    xm = (i >= bx+1)*(alpha/(hx*hx));
	    yp = (j <= ny+by-2)*(alpha/(hy*hy));
	    ym = (j >= by+1)*(alpha/(hx*hx));
	    
	    xpp = (i <= nx+bx-3)*(alpha/2*(hx*hx));
	    xmm = (i >= bx+2)*(alpha/2*(hx*hx));
	    ypp = (j <= ny+by-3)*(alpha/2*(hy*hy));
	    ymm = (j >= by+2)*(alpha/2*(hx*hx));

		/* compute the sum of weights */
		//sum = xp + xm + yp + ym;
		
		switch (choice)
		{
          case 1:
             //printf("Erste Konsistenzordnung ... ");
             sum = xp + xm + yp + ym;
             u[i][j] = (1.0-omega)*u[i][j] +  
						    omega*((-J_13 [i][j]-J_12 [i][j]*v[i][j]+xp*u[i+1][j]+xm*u[i-1][j]+yp*u[i][j+1]+ym*u[i][j-1])
													  /(J_11[i][j]+sum));
			 v[i][j] = (1.0-omega)*v[i][j] +  
							omega*((-J_23 [i][j]-J_12 [i][j]*u[i][j]+xp*v[i+1][j]+xm*v[i-1][j]+yp*v[i][j+1]+ym*v[i][j-1])
													  /(J_22[i][j]+sum));	
             break;
          case 2:
             //printf("Zweite Konsistenzordnung ... ");
             sum = xpp + xmm + ypp + ymm;
             u[i][j] = (1.0-omega)*u[i][j] +  
						    omega*((-J_13 [i][j]-J_12 [i][j]*v[i][j]+xpp*u[i+2][j]+xmm*u[i-2][j]+ypp*u[i][j+2]+ymm*u[i][j-2])
													  /(J_11[i][j]+sum));
			 v[i][j] = (1.0-omega)*v[i][j] +  
							omega*((-J_23 [i][j]-J_12 [i][j]*u[i][j]+xpp*v[i+2][j]+xmm*v[i-2][j]+ypp*v[i][j+2]+ymm*v[i][j-2])
													  /(J_22[i][j]+sum));
             break;
          case 3:
             //printf("Dritte Konsistenzordnung ... ");
             break;
          case 4:
             //printf("Vierte Konsistenzordnung ... ");
             break;
          default:
             //printf("Default ");
             break;
		}
	}

/* free memory */
FREE_MATRIX(2, nx+2*bx, ny+2*by, u_old, v_old);
}


/* ------------------------------------------------------------------------- */

void compute_motion_tensor
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float **J_11,           /* out    : entry 11 of the motion tensor            */
float **J_22,           /* out    : entry 22 of the motion tensor            */
float **J_33,           /* out    : entry 33 of the motion tensor            */
float **J_12,           /* out    : entry 12 of the motion tensor            */
float **J_13,           /* out    : entry 13 of the motion tensor            */
float **J_23            /* out    : entry 23 of the motion tensor            */
                        /*****************************************************/
)

/*
 Computes the motion tensor entries from the given image pair
*/

{
                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float fx,fy,ft;         /* image derivatives                                 */
float hx_1,hy_1;        /* time saver variables                              */
                        /*****************************************************/


/* define time saver variables */
hx_1=1.0/(2.0*hx);
hy_1=1.0/(2.0*hy);

/* mirror boundaries */
mirror_bounds_2d(f1,nx,ny,bx,by);
mirror_bounds_2d(f2,nx,ny,bx,by);

/* compute motion tensor entries entries */        
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	/* compute derivatives */

        /* TODO
         ----- fill in your code for computing the derivatives here ----
        */
		fx = 0.5*((f2 [i+1] [j]- f2[i-1][j]) + (f1 [i+1] [j]- f1[i-1][j]))*hx_1;
		fy = 0.5*((f2 [i] [j+1]- f2[i][j-1]) + (f1 [i] [j+1]- f1[i][j-1]))*hy_1;
		ft = f2 [i][j] - f1 [i][j];
		/* ------------------------------------------------------------ */

	/* set up motion tensor */

        /* TODO
         ----- fill in your code for the motion tensor entries here  ----
        */
		J_11 [i][j] = fx*fx;
		J_22 [i][j] = fy*fy;
		J_33 [i][j] = ft*ft;
		J_12 [i][j] = fx*fy;
		J_13 [i][j] = fx*ft;
		J_23 [i][j] = fy*ft;
		/* ------------------------------------------------------------ */
	
    }

}


/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
float **u,              /* out    : x-component of displacement field        */
float **v,              /* out    : y-component of displacement field        */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float m_alpha,          /* in     : smoothness weight                        */
int   n_iter            /* in     : number of iterations                     */
                        /*****************************************************/
)

/* computes optic flow with Horn/Schunck */

{

                        /*****************************************************/
int   i,j;              /* loop variables                                    */
float **J_11;           /* entry 11 of the motion tensor                     */
float **J_22;           /* entry 22 of the motion tensor                     */
float **J_33;           /* entry 33 of the motion tensor                     */
float **J_12;           /* entry 12 of the motion tensor                     */
float **J_13;           /* entry 13 of the motion tensor                     */
float **J_23;           /* entry 23 of the motion tensor                     */
                        /*****************************************************/
  


/* ---- alloc memory ---- */
ALLOC_MATRIX (6, nx+2*bx,  ny+2*by, 
	      &J_11,
	      &J_22,
	      &J_33,
	      &J_12,
	      &J_13,
	      &J_23);


/* ---- initialise displacement field with zero ---- */
for (i=bx; i<bx+nx; i++)
    for (j=by; j<by+ny; j++)	 
    {
	u[i][j]=0;
	v[i][j]=0;
    }


/* ---- compute motion tensor ---- */
compute_motion_tensor(f1, f2, nx, ny, bx, by, hx, hy,
		      J_11, J_22, J_33, J_12, J_13, J_23);



/* ---- perform Jacobi iterations ---- */
for(i=1;i<=n_iter;i++)
 {
     horn_schunck_jacobi(J_11, J_22, J_33, J_12, J_13, J_23, 
			 u, v, nx, ny, bx, by, hx, hy, m_alpha);
 }



/* ---- free memory ---- */
FREE_MATRIX (6, nx+2*bx,  ny+2*by, 
	     J_11,
	     J_22,
	     J_33,
	     J_12,
	     J_13,
	     J_23);


}
/* ------------------------------------------------------------------------- */

#endif
