
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#include "matrix.h"
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Check for proper number of arguments. */
    if(nrhs != 10) {
        mexErrMsgTxt("10 inputs are required.");
    } else if(nlhs > 2) {
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
    
    int q,k;
    int NumQuadPoints     = mxGetN(prhs[5]);
    

//     double* f  = mxGetPr(prhs[3]);
    mxComplexDouble* f_x  = mxGetComplexDoubles(prhs[3]);
    mxComplexDouble* f_y  = mxGetComplexDoubles(prhs[4]);

    double* w       = mxGetPr(prhs[5]);
    double* invjac = mxGetPr(prhs[6]);
    double* detjac = mxGetPr(prhs[7]);
    double* phi = mxGetPr(prhs[8]);
    double* gradrefphi = mxGetPr(prhs[9]);
    
    double* elements  = mxGetPr(prhs[1]);
    double gradphi[dim][nln][NumQuadPoints];

    
    /* Assembly: loop over the elements */
    int ie;
	
    #pragma omp parallel for shared(invjac,f_x,f_y,detjac,elements, myRrows, myRcoef) private(gradphi,ie,k,q) firstprivate(phi,gradrefphi, w, numRowsElements, nln2, nln)
    
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
        
        
        int ii = 0;
        int a;
        
        /* a tes, b trial */
        for (a = 0; a < nln; a = a + 1 )
        {
            
            double valore_real = 0;
            double valore_imag = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                double floc_real = 0;
                double floc_imag = 0;                
                
                floc_real = floc_real + gradphi[0][a][q] * f_x[ie+q*noe].real;
                floc_imag = floc_imag + gradphi[0][a][q] * f_x[ie+q*noe].imag;
                
                floc_real = floc_real + gradphi[1][a][q] * f_y[ie+q*noe].real;
                floc_imag = floc_imag + gradphi[1][a][q] * f_y[ie+q*noe].imag;

                valore_real = valore_real + floc_real * w[q];
                valore_imag = valore_imag + floc_imag * w[q];
            }

            
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
//             myRcoef[ie*nln+ii] = floc_real * detjac[ie];
            myRcoef[ie*nln+ii].real = valore_real * detjac[ie];
            myRcoef[ie*nln+ii].imag = valore_imag * detjac[ie];
    
            ii = ii + 1;
        }
    }
            
}

