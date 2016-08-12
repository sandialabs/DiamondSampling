#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <mex.h>

#include "cs_mex.h" 
    
/* functions implemented in sampling.c */
cs * sample_dense (const double *A, const double *B, csi M, csi K, csi N, csi numSamples, bool diamond, bool optimized, double *timers, double *col_norms_A, double *cumsums_A);
cs * sample_sparse (const cs *A, const cs *B, csi numSamples, bool diamond, bool optimized, double *timers, double *col_norms_A, double *cumsums_A);

/* sample_ab_mex: approximate sparse matrix multiply A'*B by diamond sampling */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    time_t tic, ttic = clock();
    bool sparse, diamond, optimized;
    double *timers;
    csi M, P, N, numsamples, numtimers;
    cs Amatrix, Bmatrix, *As, *Bs, *Cs, *Ds;
    double *A, *B, *col_norms_A, *cumsums_A;
        
    /* output timers */
    numtimers = 12;
    pargout[1] = mxCreateDoubleMatrix(numtimers, 1, mxREAL);
    timers = mxGetPr(pargout[1]);
    char const * labels[numtimers];
    labels[0]  = "transpose B";
    labels[1]  = "colnorms/cumsums";
    labels[2]  = "center weights";
    labels[3]  = "random number sort";
    labels[4]  = "center sample";
    labels[5]  = "center transposes";
    labels[6]  = "left sample";
    labels[7]  = "right sample";
    labels[8]  = "compress 4th edges";
    labels[9]  = "set output entry";
    labels[10] = "compress/sort output";
    labels[11] = "total";
    pargout[2] = mxCreateCharMatrixFromStrings(numtimers, (const char **) labels);
            
    if ((nargin != 6 && nargin != 8) || (nargout != 1 && nargout != 4))
    {
        mexErrMsgTxt ("Usage: [C,(wtimers,labels,successes)] = sample_ab_mex(A,B,'sparse' or 'dense','wedge' or 'diamond','optimized' or 'naive',numsamples,(col_norms_A,cumsums_A))\n"       
                "\t(parameters in parentheses are optional)\n") ;
    }
    
    /* convert text parameters from Matlab to bools */
    if ( mxIsChar(pargin[2]) == 1 || mxIsChar(pargin[3]) == 1 ) 
    {
        if (!strcmp(mxArrayToString(pargin[2]),"sparse")) sparse = true;
        else if (!strcmp(mxArrayToString(pargin[2]),"dense")) sparse = false;
        else mexErrMsgTxt ("Third input must be 'sparse' or 'dense'");
        if (!strcmp(mxArrayToString(pargin[3]),"diamond")) diamond = true;
        else if (!strcmp(mxArrayToString(pargin[3]),"wedge")) diamond = false;
        else mexErrMsgTxt("Fourth input must be 'wedge' or 'diamond'");
        if (!strcmp(mxArrayToString(pargin[4]),"optimized")) optimized = true;
        else if (!strcmp(mxArrayToString(pargin[4]),"naive")) optimized = false;
        else mexErrMsgTxt("Fifth input must be 'optimized' or 'naive'");
    } 
    else 
        mexErrMsgTxt("Third, fourth, and fifth inputs must be strings");

    /* get input matrices A and B from Matlab */
    if( sparse )
    {
        /* get cs matrices */
        As  = cs_mex_get_sparse (&Amatrix, 0, 1, pargin [0]) ;    
        Bs  = cs_mex_get_sparse (&Bmatrix, 0, 1, pargin [1]) ;
        if (As->n != Bs->m)
            mexErrMsgTxt("Inner dimensions of A*B do not match");
    } else {
        /* get dense matrices and their dimensions */
        A = mxGetPr(pargin[0]); 
        B = mxGetPr(pargin[1]); 
        M = (csi) mxGetM(pargin[0]);  
        P = (csi) mxGetN(pargin[0]);
        N = (csi) mxGetN(pargin[1]);
        if (P != (csi) mxGetM(pargin[1])) 
            mexErrMsgTxt("Inner dimensions of A*B do not match");
    }
    
    /* get numsamples from Matlab */
    if ( mxIsScalar(pargin[5]) == 1 )
       numsamples = (csi) mxGetScalar(pargin[5]);
    else
       mexErrMsgTxt ("Sixth input must be scalar") ;
    
    /* get precomputed col_norms and cumsums of A (if provided) */
    if (nargin == 8) {
        col_norms_A = mxGetPr(pargin[6]);
        cumsums_A = mxGetPr(pargin[7]);
        if ((sparse ? As->n : P) != (csi) mxGetM(pargin[6]))
        {
            mexPrintf("length of col_norms: %d; number of columns: %d\n", mxGetM(pargin[6]), As->n);
            mexErrMsgTxt("Length of col_norms_A does not match number of columns of A");
        }
        if ((sparse ? As->p[As->n] : M*P)+1 != (csi) mxGetM(pargin[7]))
        {
            mexPrintf("length of cumsums: %d; number of nonzeros: %d\n", mxGetM(pargin[7]), (sparse ? As->p[As->n] : M*P));
            mexErrMsgTxt("Length of cumsums_A does not match number of nonzeros of A (plus one)");
        }
    } else {
        col_norms_A = NULL;
        cumsums_A   = NULL;
    }
                   
    /* run code to sample entries of A*B */
    if( sparse ) 
        Ds = sample_sparse(As,Bs,numsamples,diamond,optimized,timers,col_norms_A,cumsums_A);
    else 
        Ds = sample_dense(A,B,M,P,N,numsamples,diamond,optimized,timers,col_norms_A,cumsums_A);
        
    /* output number of successes */
    pargout[3] = mxCreateDoubleScalar(Ds->nz);
            
    /* compress, dedup, and transpose back and forth (to sort) */
    tic = clock();
    Cs = cs_compress( Ds );
    cs_dupl( Cs );
    Ds = cs_transpose( Cs, 1 );
    Cs = cs_transpose( Ds, 1 );
    timers[10] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* return C in Matlab format */
    pargout [0] = cs_mex_put_sparse (&Cs) ;
        
    /* free up temporary memory */
    cs_spfree (Ds) ;
    
    timers[11] = (double)(clock()-ttic)/CLOCKS_PER_SEC;

}

