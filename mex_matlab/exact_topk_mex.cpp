extern "C" 
{
    #include "cs_mex.h"
    #include "blas.h"
}
#include <cstring>
#include <queue>

using namespace std;

struct triplet {
	csi row;
	csi col;
	double val;
};

class compare {
public:
    bool operator()(const triplet& t1, const triplet& t2)
    {
       return fabs(t1.val) > fabs(t2.val);
    }
};

/* triplet comparison function (by abs val, inverted)  */
//bool operator<(const triplet& a, const triplet& b) {
//    return fabs(a.val) > fabs(b.val);
//}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
void scatter (const cs *A, csi j, double beta, csi *w, double *x, csi mark, csi *Ci, csi &nz)
{
    csi i, p;
    for (p = A->p[j] ; p < A->p[j+1] ; p++)
    {
        i = A->i[p] ;                            /* A(i,j) is nonzero */
        if (w [i] < mark)
        {
            w [i] = mark ;                       /* i is new entry in column j */
            Ci [nz++] = i ;                      /* add i to pattern of C(:,j) */
            if (x) x [i] = beta * A->x[p] ;      /* x(i) = beta*A(i,j) */
        }
        else if (x) x [i] += beta * A->x[p] ;    /* i exists in C(:,j) already */
    }
}

/* C = A*B, return top k values in triplet format
   (edited code from CSparse library)                               */
cs * cs_multiply_topk (const cs *A, const cs *B, csi k, bool upper = false)
{
    csi p, q, j, nz, anz, *Cp, *Bp, m, n, bnz, *w, *Ci, *Bi ;
    double *x, *Bx, *Cx, val ;
    cs *C ;
    
    triplet t;
    priority_queue<triplet, vector<triplet>, compare> pq;
    
    /* check inputs */
    if (!CS_CSC (A) || !CS_CSC (B)) return (NULL) ;   
    if (A->n != B->m) return (NULL) ;
    
    /* define shortcuts */
    m = A->m ; anz = A->p [A->n] ;
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
    
    /* get workspace */
    w  = (csi *) cs_calloc (m, sizeof (csi)) ;   
    Ci = (csi *) cs_calloc (m, sizeof (csi)) ;
    x  = (double *) cs_malloc (m, sizeof (double)) ; 
    
    /* work column by column */
    for (j = 0 ; j < n ; j++)
    {
        /* work over nonzeros in jth column of B, scattering
           scaled columns of A into temporary workspace       */
        nz = 0;
        for (p = Bp [j] ; p < Bp [j+1] ; p++)
        {
            scatter (A, Bi [p], Bx [p], w, x, j+1, Ci, nz) ;
        }
              
        /* try to push triplets into priority queue */
        for (q = 0 ; q < nz ; q++)
        {   
            t.row = Ci[q];
            t.col = j;
            t.val = x[ Ci[q] ];

            /* push if pq is not full or val is in top-k so far */
            if (pq.size() < k) {
                /* if upper is set, entry must be in strict upper triangle */
                if (upper) {
                    if (t.row < t.col)
                        pq.push(t);
                } else {
                    pq.push(t);
                }
            } else if (fabs(t.val) > fabs(pq.top().val)) {
                /* if upper is set, entry must be in strict upper triangle */
                if (upper) {
                    if (t.row < t.col)
                    {
                        pq.pop();
                        pq.push(t);
                    }
                } else {      
                    pq.pop();
                    pq.push(t);
                }
            }     
        }
        
    }
    
    nz = pq.size();
    C = cs_spalloc (m, n, nz, 1, 1) ;
    while (! pq.empty()) {
       t = pq.top();
       C->i[C->nz] = t.row;
       C->p[C->nz] = t.col;
       C->x[C->nz] = t.val;
       C->nz++;
       pq.pop();
    }

    cs_free( Ci );
    return (cs_done (C, w, x, 1)) ;     /* success; free workspace, return C */
}

cs * cs_multiply_topk_dense (double *A, double *B, csi k, csi M, csi P, csi N)
{
    csi i, j, jj, cols, nz;
    char NoTrans = 'N';
    double *Ctmp, one = 1.0, zero = 0.0;
    cs *C;
    triplet t;
    priority_queue<triplet, vector<triplet>, compare> pq;
    
    /* allocate temporary memory to store pieces of output */
    /* (allow as much space as A requires)                  */
    Ctmp = (double *) mxMalloc( M*P * sizeof(double) );
    
    /* loop over column panels of output */
    for ( jj = 0; jj < N; jj += P )
    {
        /* call BLAS */
        cols = min(P,N-jj);
        dgemm(&NoTrans,&NoTrans,&M,&cols,&P,&one,A,&M,B+P*jj,&P,&zero,Ctmp,&M);
        
        /* try to push computed entries into priority queue */
        for( j = 0; j < cols; ++j)
            for( i = 0; i < M; ++i)
            {
                /* push if pq is not full or val is in top-k so far */
                if (pq.size() < k) {
                    t.row = i;
                    t.col = jj+j;
                    t.val = Ctmp[i+j*M];
                	pq.push(t);
                } else if (fabs(Ctmp[i+j*M]) > fabs(pq.top().val)) {
                    t.row = i;
                    t.col = jj+j;
                    t.val = Ctmp[i+j*M];
                    pq.pop();
                    pq.push(t);
                }
            }
    }
    
    /* empty pq into a sparse matrix */
    nz = pq.size();
    C = cs_spalloc (M, N, nz, 1, 1) ;
    while (! pq.empty()) {
       t = pq.top();
       C->i[C->nz] = t.row;
       C->p[C->nz] = t.col;
       C->x[C->nz] = t.val;
       C->nz++;
       pq.pop();
    }

    mxFree(Ctmp);
    
    return C;  

}

/* exact_topk_mex: compute A*B, returning only top k values */
void mexFunction
(
    int nargout,
    mxArray *pargout [ ],
    int nargin,
    const mxArray *pargin [ ]
)
{
    cs Amatrix, Bmatrix, *As, *Bs, *Cs, *Ds, *Es, *Fs ;
    double *A, *B;
    csi k, M, P, N;
    bool sparse, upper = false;
    
    if (nargout > 1 || nargin < 4 || nargin > 5)
    {
        mexErrMsgTxt ("Usage: C = exact_topk_mex(A,B,k,'sparse' or 'dense',['upper'])") ;
    }
    
    /* convert text parameters from Matlab to bools */
    if ( mxIsChar(pargin[3]) == 1 ) 
    {
        if (!strcmp(mxArrayToString(pargin[3]),"sparse")) sparse = true;
        else if (!strcmp(mxArrayToString(pargin[3]),"dense")) sparse = false;
        else mexErrMsgTxt ("Fourth input must be 'sparse' or 'dense'");
    } 
    else 
        mexErrMsgTxt("Fourth input must be 'sparse' or 'dense'");
    
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
    
    /* get top k parameter */
    k = mxGetScalar(pargin[2]);
    
    /* if 5th argument provided, compute only strict upper triangle */
    if (nargin == 5 && !strcmp(mxArrayToString(pargin[4]),"upper"))
        upper = true;
       
    /* compute A*B, returning only top k output values */
    if (sparse)
        Fs  = cs_multiply_topk (As, Bs, k, upper) ;
    else
        Fs  = cs_multiply_topk_dense (A, B, k, M, P, N);
        
    /* compress triplet to CSC */
    Es = cs_compress( Fs );
        
    /* transpose back and forth (to sort columns) */
    Ds = cs_transpose( Es, 1 );
    Cs = cs_transpose( Ds, 1 );
    
    /* free up temporary memory */
    cs_spfree (Fs) ;
    cs_spfree (Es) ;
    cs_spfree (Ds) ;
    
    /* return C in Matlab format */
    pargout [0] = cs_mex_put_sparse (&Cs) ;
}



