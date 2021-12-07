/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch> 
 *  Sebastiano Cominelli
 */

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
// #define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
// #define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]


void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=8) {
        mexErrMsgTxt("8 inputs are required.");
    } else if(nlhs>3) {
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
       
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);

    
    int k,l,q;
    int NumQuadPoints     = mxGetN(prhs[4]);
    

    double* si   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[4]);
    double* invjac = mxGetPr(prhs[5]);
    double* detjac = mxGetPr(prhs[6]);
    double* phi = mxGetPr(prhs[7]);
//     double* gradrefphi = mxGetPr(prhs[8]);
    
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(invjac,si,detjac,elements, myAcols, myArows, myAcoef)    private(ie,k,l,q)   firstprivate(phi, w, numRowsElements, nln2, nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int iii = 0;
        int ii = 0;
        int a, b;
    
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            for (b = 0; b < nln; b = b + 1 )
            {
                double aloc = 0;                
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    aloc = aloc + si[ie+q*noe] * phi[b+q*nln] * phi[a+q*nln] * w[q];
                }
 
                myArows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myAcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myAcoef[ie*nln2+iii] = aloc*detjac[ie];
                
                iii = iii + 1;
            }
            
            ii = ii + 1;
        }
        
    }
            
}

