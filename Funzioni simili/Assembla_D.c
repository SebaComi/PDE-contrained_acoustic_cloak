
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#include "matrix.h"

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Check for proper number of arguments. */
    if(nrhs!=7) {
        mexErrMsgTxt("7 inputs are required.");
    } else if(nlhs>2) {
        mexErrMsgTxt("Too many output arguments.");
    }

    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[1]);
    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[1]);
    int nln2    = nln*nln;
    

    plhs[0] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
//     plhs[1] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe,1, mxCOMPLEX);
           

    double* myRrows         = mxGetPr(plhs[0]);
//     double* myRcoef         = mxGetPr(plhs[1]);
    mxComplexDouble* myRcoef    = mxGetComplexDoubles(plhs[1]);
    
    int q;
    int NumQuadPoints     = mxGetN(prhs[4]);
    

//     double* f  = mxGetPr(prhs[3]);
    mxComplexDouble* f  = mxGetComplexDoubles(prhs[3]);
    double* w       = mxGetPr(prhs[4]);
    double* detjac = mxGetPr(prhs[5]);
    double* phi = mxGetPr(prhs[6]);

    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
	
    #pragma omp parallel for shared(f,detjac,elements, myRrows, myRcoef) private(ie,q) firstprivate(phi, w, numRowsElements, nln2, nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int ii = 0;
        int a;
        
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            double floc_real = 0;
            double floc_imag = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
//                  floc_real = floc_real + ( phi[a+q*nln] * f[ie+q*noe] ) * w[q];
                floc_real = floc_real + ( phi[a+q*nln] * f[ie+q*noe].real ) * w[q];
                floc_imag = floc_imag + ( phi[a+q*nln] * f[ie+q*noe].imag ) * w[q];
            }
            
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
//             myRcoef[ie*nln+ii] = floc_real * detjac[ie];
            myRcoef[ie*nln+ii].real = floc_real * detjac[ie];
            myRcoef[ie*nln+ii].imag = floc_imag * detjac[ie];
    
            ii = ii + 1;
        }
    }
            
}

