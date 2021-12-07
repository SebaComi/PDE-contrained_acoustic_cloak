/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch> 
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

// [Crows, Ccols, Ccoefx, Ccoefy] = ADR_assembler_C_omp(MESH.dim, MESH.elements, FE_SPACE.numElemDof,FE_SPACE.quad_weights, ...
//                                                      MESH.invjac(index_subd,:,:), MESH.jac(index_subd), FE_SPACE.phi, FE_SPACE.dphi_ref);

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=8) {
        mexErrMsgTxt("15 inputs are required.");
    } else if(nlhs>4) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[1]);
    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[1]);
    int nln2    = nln*nln;
    
    /**/
    plhs[0] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
       
    
    double* myCrows     = mxGetPr(plhs[0]);
    double* myCcols     = mxGetPr(plhs[1]);
    double* myCcoef_x   = mxGetPr(plhs[2]);
    double* myCcoef_y   = mxGetPr(plhs[3]);
    
    int k,l,q;
    
    int NumQuadPoints     = mxGetN(prhs[3]);
    
    double* w   = mxGetPr(prhs[3]);
    double* invjac = mxGetPr(prhs[4]);
    double* detjac = mxGetPr(prhs[5]);
    double* phi = mxGetPr(prhs[6]);
    double* gradrefphi = mxGetPr(prhs[7]);
    
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(invjac,detjac,elements, myCcols,myCrows, myCcoef_x,myCcoef_y) private(gradphi,ie,k,l,q) firstprivate(phi,gradrefphi, w, numRowsElements, nln2, nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int d1, d2;
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }
        
        int iii = 0;
        int ii = 0;
        int a, b;
    
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                double aloc_x = 0;
                double aloc_y = 0; 
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    aloc_x = aloc_x + gradphi[0][b][q] * phi[a+q*nln] * w[q];   // Sommo i trasporti lungo x
                    aloc_y = aloc_y + gradphi[1][b][q] * phi[a+q*nln] * w[q];   // Sommo i trasporti lungo y
                }
 
                myCrows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myCcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myCcoef_x[ie*nln2+iii] = aloc_x * detjac[ie];
                myCcoef_y[ie*nln2+iii] = aloc_y * detjac[ie];
                
                iii = iii + 1;
            }

        }
        
    }
            
}

