
#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>
#include "blas.h"
#include <string.h>
#define IsNonZero(d) ((d)!=0.0)
#define INVJAC(i,j,k) invjac[i+(j+k*dim)*noe]
#define GRADREFPHI(i,j,k) gradrefphi[i+(j+k*NumQuadPoints)*nln]

// Questa funzione si usa così:
// [Brows, Bcols, Balte, Bcoef] = Assembla_B( MESH.dim, MESH.elements, FE_SPACE.numElemDof, ...
//                                            FE_SPACE.quad_weigths, MESH.invjac, MESH.jac, FE_SPACE.phi )
// Assumiamo che la matrice sparsa abbia nzmax = 1% dell'intera matrice, facilmente si scende anche sotto lo 0.1% !!
void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /* Check for proper number of arguments. */
    if(nrhs!=7) {
        mexErrMsgTxt("8 inputs are required.");
    } else if(nlhs>4) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    // Definiamo le cose per la matrice sparsa
                        // Declare variable
    /*{    mwSize m,n;
        mwSize nzmax;
        mwIndex *irs,*jcs;
        double *pr,*sr;
        double percent_sparse;
        
        // Get the size and pointers to input data
        m  = mxGetM(prhs[0]);
        n  = mxGetN(prhs[0]);
        pr = mxGetPr(prhs[0]);

//         Allocate space for sparse matrix
//          * NOTE:  Assume at most 1% of the data is sparse.  Use ceil
//          * to cause it to round up.
        

        percent_sparse = 0.01;
        nzmax=(mwSize)ceil((double)m*(double)n*percent_sparse);

        plhs[3] = mxCreateSparse(m,n,nzmax,mxREAL);
        sr  = mxGetPr(plhs[0]);
        irs = mxGetIr(plhs[0]);
        jcs = mxGetJc(plhs[0]);
    */
    // Queste sono le cose vecchie
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[1]);
    double* nln_ptr = mxGetPr(prhs[2]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[1]);
    int nln3    = nln*nln*nln;
    
    /**/
//     mxArray* plhs_finto[3]
    plhs[0] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln3*noe,1, mxREAL);
    
    double* Brows   = mxGetPr(plhs[0]);
    double* Bcols   = mxGetPr(plhs[1]);
    double* Balte   = mxGetPr(plhs[2]);
    double* Bcoef   = mxGetPr(plhs[3]);
    
    int k,l;
    int q;
    int NumQuadPoints = mxGetN(prhs[3]);
    
    double* w   = mxGetPr(prhs[3]);
    double* invjac = mxGetPr(prhs[4]);
    double* detjac = mxGetPr(prhs[5]);
    double* phi = mxGetPr(prhs[6]);
    
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[1]);

    /* Assembly: loop over the elements */
    int ie;
            
    #pragma omp parallel for shared(invjac,detjac,elements, Bcols,Brows,Bcoef) private(gradphi,ie,k,l,q) firstprivate(phi,gradrefphi, w, numRowsElements, nln3, nln)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        int iii = 0;
        int ii = 0;
        int a, b, c;
        
        for(c = 0; c < nln; c = c + 1)
        {
            /* a tes, b trial */
            for (a = 0; a < nln; a = a + 1 )
            {
                for (b = 0; b < nln; b = b + 1 )
                {
                    double aloc = 0;                
                    for (q = 0; q < NumQuadPoints; q = q + 1 )  // Questo somma sui punti di quadratura q
                    {
                        aloc = aloc + phi[b+q*nln] * phi[a+q*nln] * phi[c+q*nln] * w[q]; // Solo il termine di reazione
                    }

                    Brows[ie*nln3+iii] = elements[a+ie*numRowsElements];
                    Bcols[ie*nln3+iii] = elements[b+ie*numRowsElements];
                    Balte[ie*nln3+iii] = elements[c+ie*numRowsElements]; // Questa è la terza dimensioine di B
                    Bcoef[ie*nln3+iii] = aloc*detjac[ie];

                    iii = iii + 1;
                } // for di b
                
                ii = ii + 1;
            } // for di a
            
//             A questo punto i vettori Brows, Bcols e Balte richiamano più volte i punti
            
            

/*  Questa parte sarebbe per la matrice sparsa
//         Copy nonzeros
        mxIndex jj, kk;
        kk = 0;
        for (jj = 0; (jj < n); jj++) {
            mwSize i;
            jcs[jj] = kk;
            for (i=0; (i<m ); i++) {
                if (IsNonZero(pr[i])) {

                    /* Check to see if non-zero element will fit in
                     * allocated output array.  If not, increase percent_sparse
                     * by 10%, recalculate nzmax, and augment the sparse array
                     *
                    if (kk >= nzmax) {
                        mwSize oldnzmax = nzmax;
                        percent_sparse += 0.1;
                        nzmax = (mwSize)ceil((double)m*(double)n*percent_sparse);

//                      make sure nzmax increases at least by 1
                        if (oldnzmax >= nzmax) {
                            nzmax = oldnzmax+1;
                        }

                        mxSetNzmax(plhs[0], nzmax);
                        mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
                        mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

                        sr  = mxGetPr(plhs[0]);
                        irs = mxGetIr(plhs[0]);
                        
                    }
                    sr[kk] = pr[i];
                    irs[kk] = i;
                    kk++;
                }
            }
            pr += m;
        }
        jcs[n] = kk;
*/
            
            
            
            
        } // for di c
    }
            
}

