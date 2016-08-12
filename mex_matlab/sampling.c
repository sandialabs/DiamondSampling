#include <stdlib.h> 
#include <time.h>
#include "cs_mex.h" 

    
#define SIGN(x) ((x>0) ? 1 : -1)

/* compare function for qsort */
int fcompare (const void * a, const void * b)
{
  if (*(double*) a > *(double*) b)
    return(1);
  return(-1);
}

/* fills in array s with pre-sorted uniform samples between 0 and max 
       ** thanks to Cliff Anderson-Bergman for this trick  **         */
void generate_sorted_uniform_samples(double *s, csi n, double max)
{
    csi t;
    s[0] = ( 1 - exp( log((double)rand()/(double)RAND_MAX) / n) ) * max;
    for(t = 1; t < n; ++t)
        s[t] = (1.0 - exp( log( (double)rand() / (double)RAND_MAX ) / (n - t) ) ) * (max - s[t-1]) + s[t-1];
}

/* finds largest index of sorted (ascending) array a whose value is less than s */
csi bin_search(double *a, csi lb, csi ub, double s) 
{
    csi m;
    
    if (s < a[lb] || s > a[ub])
        mexErrMsgTxt("bin search: s out of bounds\n");
    
    while (lb < ub-1) 
    {
        m = (lb + ub )/2;
        if (s < a[m])
            ub = m;
        else
            lb = m;
    }
    return (lb);
}

/* finds index of sorted (ascending) array a whose value is s
   and returns -1 if s is not in array */
csi bin_search_locate(csi *a, csi lb, csi ub, csi s) 
{
    csi m;
    
    /* look for quick returns, ensure a[lb] < s < a[ub] */
    if (s < a[lb] || s > a[ub])
        return -1;
    if (s == a[lb])
        return lb;
    if (s == a[ub])
        return ub;
    
    /* find i such that a[i] < s < a[i+1] or return success */
    while (lb < ub-1) 
    {
        m = (lb + ub )/2;
        if (s == a[m])
            return m;
        else if (s < a[m])
            ub = m;
        else
            lb = m;
    }
    
    /* return failure */
    return -1;
}

/* return value of (i,j) entry of A */
double entry(const cs *A, csi i, csi j)
{
    csi k;
    
    /* use binary search (assumes each column's row indices are sorted) */
    k = bin_search_locate(A->i,A->p[j],A->p[j+1]-1,i);
    return ( k != -1 ? A->x[k] : 0.0 );
    
    /* use linear search 
    for (k = A->p[j]; k < A->p[j+1]; k++)
        if (A->i[k] == i)
            return A->x[k];
    return 0.0; */
    
}

/* compute cumulative sums and column norms of sparse matrix 
     length of colnorms needs to be number of columns of A   
     length of cumsums needs to be nnz of A plus 1           */
void cumulative_sums( const cs *A, double *colnorms, double *cumsums )
{
    csi i,l;
    cumsums[0] = 0.0;
    for (l = 0; l < A->n; l++)
    {   
        colnorms[l] = 0.0;
        for (i = A->p[l]; i < A->p[l+1]; i++)
        {
            /* accumulate 1-norm of lth column of A */
            colnorms[l] += fabs( A->x[i] );
            
            /* compute cumsum for nonzero values of A */
            cumsums[i+1] = cumsums[i] + fabs( A->x[i] );
        } 
    }
}

/* compute cumulative sums and column norms of dense MxN matrix 
     length of colnorms needs to be number of columns of A   
     length of cumsums needs to be nnz of A plus 1           */
void cumulative_sums_dense( const double *A, csi M, csi N, double *colnorms, double *cumsums )
{
    csi i,l;
    cumsums[0] = 0.0;
    for (l = 0; l < N; l++)
    {   
        colnorms[l] = 0.0;
        for (i = l*M; i < (l+1)*M; i++)
        {
            /* accumulate 1-norm of lth column of A */
            colnorms[l] += fabs( A[i] );
            
            /* compute cumsum for nonzero values of A */
            cumsums[i+1] = cumsums[i] + fabs( A[i] );
        } 
    }
}

/* computes approximation of A*B using 3-path/diamond or wedge samples 
   assumes A (M x K) and B (K x N) are sparse and stored in CSC        */
cs * sample_sparse (const cs *A, const cs *B, csi numSamples, bool diamond, bool optimized, double* timers, double* col_norms_A, double* cumsums_A)
{
    time_t tic;
    csi M, K, N, KN, NK, i, j, k, l, t, e, ind1, ind2, ii, jj, kk, ll, dd;
    csi *centers, *lefts, *rights, *jjs, *Ep, *bins, *B_cols;
    double x, ce_wt_tmp, tot_wt, wt_val, l_wt, val;
    double *S, *col_norms_B, *cumsums_B, *W;
    cs *C, *Bt, *D, *E, *Et, *F;
    bool preprocessed;
    
    /* A is MxK and B is KxN */
    M = A->m;
    K = A->n;
    N = B->n;
    if (B->m != K)
        mexErrMsgTxt ("sample_sparse: inner dimensions don't match") ;
    
    /* check if A has been preprocessed */
    if (col_norms_A == NULL && cumsums_A == NULL)
        preprocessed = false;
    else if (col_norms_A && cumsums_A)
        preprocessed = true;
    else
        mexErrMsgTxt ("sample_sparse: preprocessed arrays aren't consistent") ;
    
    /* set N and K switches because of wedge/diamond discrepancy */
    KN = (diamond ? K : N);
    NK = (diamond ? N : K);
                
    /* allocate result (with values, in triplet form) */
    C = cs_spalloc (M, N, numSamples, 1, 1) ;
   
    /* allocate temporary memory */
    if (!preprocessed)
        col_norms_A = (double *) mxMalloc( K * sizeof(double) );
    col_norms_B = (double *) mxMalloc( NK * sizeof(double) );
    if (!preprocessed)
        cumsums_A   = (double *) mxMalloc( (A->p[K]+1) * sizeof(double) );  
    cumsums_B   = (double *) mxMalloc( (B->p[N]+1) * sizeof(double) );
    lefts       = (csi *)    mxMalloc( numSamples * sizeof(csi) );
    rights      = (csi *)    mxMalloc( numSamples * sizeof(csi) );
    if (diamond) {
        centers = (csi *)    mxMalloc( numSamples * sizeof(csi) );
        if (optimized) {
            S   = (double *) mxMalloc( numSamples * sizeof(double) );
            jjs = (csi *)    mxMalloc( numSamples * sizeof(csi) );
            Ep  = (csi *)    mxMalloc( N * sizeof(csi) );
        } else {
            W      = (double *) mxMalloc( (B->p[N]+1) * sizeof(double) );
            B_cols = (csi *)    mxMalloc( B->p[N] * sizeof(csi) );
        }
    } else {
        W = (double *) mxMalloc( (K+1) * sizeof(double) );
        if (optimized)
            bins = (csi *) mxCalloc( K, sizeof(csi) );
        else
            centers = (csi *) mxMalloc( numSamples * sizeof(csi) );
    }
    
    /* allocate matrix for possible output and central edges */
    if (diamond && optimized)
    {
        /* possible 4th edges with values, in triplet form */
        D = cs_spalloc (M, K, numSamples, 1, 1) ;
        /* central edges without values, in compressed form */
        E = cs_spalloc (K, N, numSamples, 0, 0) ;
    }

    /* convert B to CSR for wedge case */   
    if (!diamond)
    {   
        tic = clock();
        Bt = cs_transpose (B, 1);
        timers[0] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    } else {
        timers[0] = 0;
    }
    
    /* precompute column norms and cumulative sums of edge weights of A and B/Bt */
    tic = clock();
    if (!preprocessed)
        cumulative_sums( A, col_norms_A, cumsums_A );
    if (diamond)
        cumulative_sums( B, col_norms_B, cumsums_B );
    else
        cumulative_sums( Bt, col_norms_B, cumsums_B );
    timers[1] = (double)(clock()-tic)/CLOCKS_PER_SEC;
        
    /* compute total weight of 3-paths or wedges 
       (and store cumsums if not optimized diamond) */
    tic = clock();
    tot_wt = 0.0;
    if (diamond) {
        if (!optimized) 
            W[0] = 0.0;
        for (j = 0; j < N; ++j)
            for (e = B->p[j]; e < B->p[j+1]; ++e)
            {
                if (optimized) 
                    tot_wt += col_norms_A[B->i[e]] * fabs(B->x[e]) * col_norms_B[j];
                else
                {
                    W[e+1] = W[e] + col_norms_A[B->i[e]] * fabs(B->x[e]) * col_norms_B[j];
                    B_cols[e] = j;
                }
            }
        if (!optimized)
            tot_wt = W[B->p[N]];
    } else {
        W[0] = 0.0;
        for (l = 0; l < K; l++) 
            W[l+1] = W[l] + col_norms_A[l] * col_norms_B[l];
        tot_wt = W[K];
    }
    timers[2] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* generate samples and sort them */ 
    tic = clock();
    if (diamond && optimized)
    {
/*         for(t = 0; t < numSamples; t++)
             S[t] = ((double)rand() / (double)RAND_MAX) * tot_wt;
         qsort(S,numSamples,sizeof(double),fcompare); */
        generate_sorted_uniform_samples( S, numSamples, tot_wt);
    }
    timers[3] = (double)(clock()-tic)/CLOCKS_PER_SEC;

    /* choose center edges/vertices */
    tic = clock();
    if (diamond)
    {
        if (optimized) {
            /* consider random numbers in order and walk through edges of B;
               sparse matrix E will store central edges (with duplication)   */
            t = 0;
            wt_val = 0.0;
            for (j = 0; j < N; ++j)
            {
                E->p[j] = t;
                for (k = B->p[j]; k < B->p[j+1]; ++k)
                {
                    wt_val += col_norms_A[B->i[k]] * fabs(B->x[k]) * col_norms_B[j];
                    while( S[t] < wt_val && t < numSamples)
                    {
                        E->i[t] = B->i[k];
                        centers[t] = k;
                        jjs[t] = j;
                        t++;
                    }
                }
            }
            E->p[N] = numSamples;
        } else {
            /* perform binary search for each sample (store center edge index) */
            for( t = 0; t < numSamples; t++) 
            {
                x = ((double)rand() / (double)RAND_MAX) * tot_wt;
                centers[t] = bin_search(W, 0, B->p[N], x);
            }
        }
    } else {
        /* perform binary search for each sample */
        for( t = 0; t < numSamples; t++) 
        {
            x = ((double)rand() / (double)RAND_MAX) * tot_wt;
            e = bin_search(W, 0, K, x);
            if (optimized)
                /* increment count for center vertex */
                bins[e]++;
            else
                /* store center vertex */
                centers[t] = e;
        }
    }
    timers[4] = (double)(clock()-tic)/CLOCKS_PER_SEC;

    /* choose left neighbor of center edge/vertex */
    if (diamond) {
        if (optimized) {
            /* transpose E (without values) to group left endpoints */
            tic = clock();
            Et = cs_transpose(E,0);
            timers[5] = (double)(clock()-tic)/CLOCKS_PER_SEC;

            /* store "tranposed" left neighbors temporarily in rights array */
            tic = clock();
            for( kk = 0; kk < K; ++kk)
                for (e = Et->p[kk]; e < Et->p[kk+1]; ++e)
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[A->p[kk]];
                    rights[e] = bin_search(cumsums_A, A->p[kk], A->p[kk+1], x);
                }
            timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;

            /* transpose values back to match ordering of right neighbors */
            /* Ep is copy of E->p that can be updated as lefts gets filled in */
            tic = clock();
            for( j = 0; j < N; j++ )
                Ep[j] = E->p[j];
            /* "transpose" rights back into lefts */
            for( l = 0; l < K; ++l )
                for ( t = Et->p[l]; t < Et->p[l+1]; t++ )
                    lefts[Ep[Et->i[t]]++] = rights[t] ;
            timers[5] += (double)(clock()-tic)/CLOCKS_PER_SEC;
        } else {
            /* corresponds to center transposes timer */
            timers[5] = 0.0;
            
            tic = clock();
            for( t = 0; t < numSamples; ++t )
            {
                kk = B->i[centers[t]];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[A->p[kk]];
                lefts[t] = bin_search(cumsums_A, A->p[kk], A->p[kk+1], x);
            }
            timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;
        }
    } else {
        /* corresponds to center transposes timer */
        timers[5] = 0.0;
        
        tic = clock();
        if (optimized) {
            t = 0;
            for( kk = 0; kk < K; ++kk ) 
                for( l = 0; l < bins[kk]; ++l )
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[A->p[kk]];
                    lefts[t++] = bin_search(cumsums_A, A->p[kk], A->p[kk+1], x);
                }
        } else {
            for( t = 0; t < numSamples; ++t )
            {
                kk = centers[t];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[A->p[kk]];
                lefts[t] = bin_search(cumsums_A, A->p[kk], A->p[kk+1], x);
            }
        }
        timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    }
   
    /* choose right neighbor of center edge/vertex */
    tic = clock();
    if (diamond) {
        if (optimized) {
            for( jj = 0; jj < N; ++jj )
                for (t = E->p[jj]; t < E->p[jj+1]; ++t )
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[B->p[jj]];
                    rights[t] = bin_search(cumsums_B, B->p[jj], B->p[jj+1], x);
                }
        } else {
            for( t = 0; t < numSamples; ++t )
            {
                jj = B_cols[centers[t]];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[B->p[jj]];
                rights[t] = bin_search(cumsums_B, B->p[jj], B->p[jj+1], x);
            }
        }
    } else {
        if (optimized) {
            t = 0;
            for( jj = 0; jj < K; ++jj ) 
                for( l = 0; l < bins[jj]; ++l )
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[Bt->p[jj]];
                    rights[t++] = bin_search(cumsums_B, Bt->p[jj], Bt->p[jj+1], x);
                }
        } else {
            for( t = 0; t < numSamples; ++t )
            {
                jj = centers[t];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[Bt->p[jj]];
                rights[t] = bin_search(cumsums_B, Bt->p[jj], Bt->p[jj+1], x);
            }
        }
    }
    timers[7] = (double)(clock()-tic)/CLOCKS_PER_SEC;
  
    /* set output entry */
    if (diamond) {
        if (optimized) {
            /* dereference lefts and rights to store actual ii and ll  
               indices temporarily in output triplet matrix (with values),
               then compress into CSC to order by column                   */
            tic = clock();
            for( t = 0; t < numSamples; t++ )
            {
                D->i[t] = A->i[lefts[t]];
                D->p[t] = B->i[rights[t]];
                D->x[t] = (double) t;
            }
            D->nz = numSamples;
            F = cs_compress( D );
            timers[8] = (double)(clock()-tic)/CLOCKS_PER_SEC;

            /* look up A(ii,ll) vals and if 4th edge is found, set output entry */
            tic = clock();
            C->nz = 0;
            for( ll = 0; ll < K; ll++ )
                for(t = F->p[ll]; t < F->p[ll+1]; t++ )
                {
                    ii = F->i[t];
                    dd = (csi) F->x[t];
                    val = entry(A,ii,ll);
                    if (val)
                    {
                        C->i[C->nz] = ii;        
                        C->p[C->nz] = jjs[dd];        
                        C->x[C->nz] = val * SIGN(A->x[lefts[dd]]) * SIGN(B->x[rights[dd]]) * SIGN(B->x[centers[dd]]);       
                        C->nz++;
                    }
                }
            timers[9] = (double)(clock()-tic)/CLOCKS_PER_SEC;
        } else {
            /* corresponds to compress 4th edges time */
            timers[8] = 0;
            
            tic = clock();
            for( t = 0; t < numSamples; ++t )
            {
                /* dereference lefts and rights to get actual ii and ll values */
                ii = A->i[lefts[t]];
                ll = B->i[rights[t]];

                /* see if (ii,ll) exists in A */
                val = entry(A,ii,ll);

                /* set output entry if diamond found */
                if (val)
                {
                    C->i[C->nz] = ii;        
                    C->p[C->nz] = B_cols[centers[t]];        
                    C->x[C->nz] = val * SIGN(A->x[lefts[t]]) * SIGN(B->x[rights[t]]) * SIGN(B->x[centers[t]]);       
                    C->nz++;
                }
            }
            timers[9] = (double)(clock()-tic)/CLOCKS_PER_SEC;
        }
    } else {

        /* corresponds to compress 4th edges time */
        timers[8] = 0;

        /* set output entry by dereferencing lefts and rights */
        tic = clock();
        for(t = 0; t < numSamples; t++) 
        { 
            C->i[C->nz] = A->i[lefts[t]]; 
            C->p[C->nz] = Bt->i[rights[t]];
            C->x[C->nz] = SIGN(A->x[lefts[t]]) * SIGN(Bt->x[rights[t]]);
            C->nz++;
        }
        timers[9] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    }  

    mxFree(lefts);
    mxFree(rights);
    mxFree(cumsums_B);
    if (!preprocessed)
        mxFree(cumsums_A);
    mxFree(col_norms_B);
    if (!preprocessed)
        mxFree(col_norms_A);
    if (diamond) {
        mxFree(centers);
        if (optimized) {
            mxFree(S);
            mxFree(jjs);
            mxFree(Ep);
            cs_spfree(D);
            cs_spfree(E); 
            cs_spfree(Et);
            cs_spfree(F);
        } else {
            mxFree(W);
            mxFree(B_cols);
        }
    } else {
        cs_spfree (Bt);
        mxFree(W);
        if (optimized)
            mxFree(bins);
        else
            mxFree(centers);
    }
        
    return C;
}

/* computes approximation of A*B using 3-path/diamond or wedge samples 
   assumes A (M x K) and B (K x N) are dense and stored in col-major */
cs * sample_dense (const double *A, const double *B, csi M, csi K, csi N, csi numSamples, bool diamond, bool optimized, double* timers, double *col_norms_A, double *cumsums_A)
{
    time_t tic;
    double timer;
    csi i, l, j, t, e, ind1, ind2, ii, jj, kk, ll, KN, NK;
    csi *B_col, *centers, *lefts, *rights, *bins, *Ep;
    double *Bt, x, ce_wt_tmp, tot_wt, wt_val, l_wt, val;
    double *S, *col_norms_B, *cumsums_B, *W;
    cs *C, *D, *E, *Et;
    bool preprocessed;
        
    /* allocate result (with values, in triplet form) */
    C = cs_spalloc (M, N, numSamples, 1, 1) ;
    
    /* check if A has been preprocessed */
    if (col_norms_A == NULL && cumsums_A == NULL)
        preprocessed = false;
    else if (col_norms_A && cumsums_A)
        preprocessed = true;
    else
        mexErrMsgTxt ("sample_sparse: preprocessed arrays aren't consistent") ;
    
    /* set N and K switches because of wedge/diamond discrepancy */
    KN = (diamond ? K : N);
    NK = (diamond ? N : K);
    
    /* allocate temporary memory */
    if (!preprocessed)
        col_norms_A = (double *) mxMalloc( K * sizeof(double) );
    col_norms_B = (double *) mxMalloc( NK * sizeof(double) );
    if (!preprocessed)
        cumsums_A   = (double *) mxMalloc( (M*K+1) * sizeof(double) );  
    cumsums_B   = (double *) mxMalloc( (K*N+1) * sizeof(double) );
    lefts       = (csi *)    mxMalloc( numSamples * sizeof(csi) );
    rights      = (csi *)    mxMalloc( numSamples * sizeof(csi) );
    if (diamond) {
        centers = (csi *) mxMalloc( numSamples * sizeof(csi) );
        if (optimized) {
            S  = (double *) mxMalloc( numSamples * sizeof(double) );
            Ep = (csi *)    mxMalloc( N * sizeof(csi) );
        } else {
            W = (double *) mxMalloc( (K*N+1) * sizeof(double) );
        }
    } else {
        W = (double *) mxMalloc( (K+1) * sizeof(double) );
        if (optimized) {
            bins = (csi *) mxCalloc( K, sizeof(csi) );
        } else {
            centers = (csi *) mxMalloc( numSamples * sizeof(csi) );
        }
    }
    
    /* allocate matrix for possible output and central edges */
    if (diamond && optimized)
    {
        /* possible 4th edges with values, in triplet form */
        D = cs_spalloc (M, K, numSamples, 1, 1) ;
        /* central edges without values, in compressed form */
        E = cs_spalloc (K, N, numSamples, 0, 0) ;
    }
    
    /* convert B to row-major ordering */
    tic = clock();
    if (!diamond)
    {
        Bt = (double *) mxMalloc( K*N * sizeof(double) );
        for (j=0; j<N; ++j)
            for (l=0; l<K; ++l)
                Bt[j+l*N] = B[l+j*K];
    } 
    timers[0] = (double)(clock()-tic)/CLOCKS_PER_SEC;
        
    /* precompute column norms and cumulative sums */
    tic = clock();
    if (!preprocessed)
        cumulative_sums_dense( A, M, K, col_norms_A, cumsums_A );
    if (diamond)
        cumulative_sums_dense( B, K, N, col_norms_B, cumsums_B );
    else
        cumulative_sums_dense( Bt, N, K, col_norms_B, cumsums_B );
    timers[1] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* compute total weight of 3-paths or wedges */
    tic = clock();
    if (diamond) {
        if (optimized)
            tot_wt = 0.0;
        else
            W[0] = 0.0;
        for (j = 0; j < N; ++j)
        {
            for (l = 0; l < K; ++l)
            {
                if (optimized) {
                    tot_wt += col_norms_A[l] * fabs(B[l+j*K]) * col_norms_B[j];
                } else {
                    W[l+j*K+1] = W[l+j*K] + col_norms_A[l] * fabs(B[l+j*K]) * col_norms_B[j];
                }
            }
        }
        if (!optimized)
            tot_wt = W[K*N];
    } else {
        W[0] = 0.0;
        for (l = 0; l < K; l++)
            W[l+1] = W[l] + col_norms_A[l] * col_norms_B[l];
        tot_wt = W[K];
    }
    timers[2] = (double)(clock()-tic)/CLOCKS_PER_SEC;
               
    /* generate samples and sort them */ 
    tic = clock();
    if (diamond && optimized)
    {
/*         for(t = 0; t < numSamples; ++t)
             S[t] = ((double)rand() / (double)RAND_MAX) * tot_wt;
         qsort(S,numSamples,sizeof(double),fcompare); */
        generate_sorted_uniform_samples( S, numSamples, tot_wt);
    }
    timers[3] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* choose centers */
    tic = clock();
    if (diamond) {
        if (optimized) {
            /* consider random numbers in order and walk through edges of B;
               sparse matrix E will store central edges (with duplication)   */
            t = 0;
            wt_val = 0.0;
            for (j = 0; j < N; ++j)
            {
                E->p[j] = t;
                for (l = 0; l < K; ++l)
                {
                    wt_val += col_norms_A[l] * fabs(B[l+j*K]) * col_norms_B[j];
                    while( t < numSamples && S[t] < wt_val )
                    {
                        E->i[t] = l;
                        centers[t] = l+j*K;
                        t++;
                    }
                }
            }
            E->p[N] = numSamples;
        } else {
            /* perform binary search for each sample (store center edge index) */
            for( t = 0; t < numSamples; t++) 
            {
                x = ((double)rand() / (double)RAND_MAX) * tot_wt;
                centers[t] = bin_search(W, 0, K*N, x);
            }
        }
    } else {
        /* perform binary search for each sample (store center vertex index) */
        for( t = 0; t < numSamples; ++t) 
        {
            x = ((double)rand() / (double)RAND_MAX) * tot_wt;
            e = bin_search(W, 0, K, x);
            if (optimized) {
                /* increment count for center vertex */
                bins[e]++;
            } else {
                /* store center vertex */
                centers[t] = e;
            }
        }
    }
    timers[4] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* choose lefts */
    if (diamond) {
        if (optimized) {
            /* transpose E (without values) to group left endpoints */
            tic = clock();
            Et = cs_transpose(E,0);
            timers[5] = (double)(clock()-tic)/CLOCKS_PER_SEC;

            /* store "tranposed" left neighbors temporarily in rights array */
            tic = clock();
            for( kk = 0; kk < K; ++kk)
            {
                for (e = Et->p[kk]; e < Et->p[kk+1]; ++e)
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[kk*M];
                    rights[e] = bin_search(cumsums_A, kk*M, (kk+1)*M, x);
                }
            }
            timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;

            /* transpose values back to match ordering of right neighbors */
            /* Ep is copy of E->p that can be updated as lefts gets filled in */
            tic = clock();
            for( j = 0; j < N; j++ )
                Ep[j] = E->p[j];
            /* "transpose" rights back into lefts */
            for( l = 0; l < K; ++l )
                for ( t = Et->p[l]; t < Et->p[l+1]; t++ )
                    lefts[Ep[Et->i[t]]++] = rights[t] ;
            timers[5] += (double)(clock()-tic)/CLOCKS_PER_SEC;
        } else {
            /* corresponds to center transposes */
            timers[5] = 0;
            
            tic = clock();
            for( t = 0; t < numSamples; ++t )
            {
                kk = centers[t] % K;
                x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[kk*M];
                lefts[t] = bin_search(cumsums_A, kk*M, (kk+1)*M, x);
            }
            timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;
            
        }
    } else {
        /* corresponds to center transposes */
        timers[5] = 0;
        
        tic = clock();
        if (optimized) {
            t = 0;
            for( kk = 0; kk < K; ++kk ) 
                for( l = 0; l < bins[kk]; ++l )
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[kk*M];
                    lefts[t] = bin_search(cumsums_A, kk*M, (kk+1)*M, x);
                    t++;
                }
        } else {
            for( t = 0; t < numSamples; ++t )
            {
                kk = centers[t];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_A[kk] + cumsums_A[kk*M];
                lefts[t] = bin_search(cumsums_A, kk*M, (kk+1)*M, x);
            }
        }
        timers[6] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    }
    
    /* choose rights */
    tic = clock();
    if (diamond) {
        for( t = 0; t < numSamples; ++t )
        {
            jj = centers[t] / K;
            x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[jj*KN];
            rights[t] = bin_search(cumsums_B, jj*KN, (jj+1)*KN, x);
        }
    } else {
        if (optimized) {
            t = 0;
            for( jj = 0; jj < K; ++jj ) 
                for( l = 0; l < bins[jj]; ++l )
                {
                    x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[jj*KN];
                    rights[t] = bin_search(cumsums_B, jj*KN, (jj+1)*KN, x);
                    t++;
                }
        } else {
            for( t = 0; t < numSamples; ++t )
            {
                jj = centers[t];
                x = ((double)rand() / (double)RAND_MAX) * col_norms_B[jj] + cumsums_B[jj*KN];
                rights[t] = bin_search(cumsums_B, jj*KN, (jj+1)*KN, x);
            }
        }
    }
    timers[7] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* corresponds to compress 4th edges */
    timers[8] = 0;
    
    /* set output entry */
    tic = clock();
    for( t = 0; t < numSamples; ++t )
    {
        /* set output entry */
        ii = lefts[t] % M;
        ll = rights[t] % KN;
        if (diamond) {
            jj = centers[t] / K;
            C->i[C->nz] = ii;        
            C->p[C->nz] = jj;        
            C->x[C->nz] = A[ii+ll*M] * SIGN(A[lefts[t]]) * SIGN(B[rights[t]]) * SIGN(B[centers[t]]);
        } else {
            C->i[C->nz] = ii;                      
            C->p[C->nz] = ll;                      
            C->x[C->nz] = SIGN(A[lefts[t]]) * SIGN(Bt[rights[t]]);
        }
        C->nz++;
    }
    timers[9] = (double)(clock()-tic)/CLOCKS_PER_SEC;
    
    /* free memory */
    mxFree(cumsums_B);
    if (!preprocessed)
        mxFree(cumsums_A);
    mxFree(col_norms_B);
    if (!preprocessed)
        mxFree(col_norms_A);
    mxFree(lefts);
    mxFree(rights);
    if (diamond) {
        mxFree(centers);
        if (optimized) {
            mxFree(S);
            mxFree(Ep);
            cs_spfree(D);
            cs_spfree(E);
            cs_spfree(Et);
        } else {
            mxFree(W);
        }
    } else {
        mxFree(Bt);
        mxFree(W);
        if (optimized) {
            mxFree(bins);
        } else {
            mxFree(centers);
        }
    }
    
    return C;
}
