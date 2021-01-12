//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_wamf: Weighted approximate minimum fill ordering
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: This function is an implementation of the weighted approximate minimum
 * fill ordering. In essence, we modify the minimum degree and minimum fill-in algorithms
 * in order to account for the bit lengths of entries as well. At each iteration, this 
 * algorithm attempts to select a node which reduces a combination of fill in and 
 * entry size criterion.
 * 
 * This function has the following input arguments:
 * 
 *  order: an int which determines how the matrix will be analyzed. It can take the
 *  following values:
 *      0: The matrix will be analyzed directly with no modification. This is appropriate
 *         for SPD instances.
 * 
 *      1: The graph of A'+A will be analyzed. This is also appropriate for SPD instances
 *         and some unsymmetric instances
 * 
 *      2: Analyze the pattern of A'*A. This is appropriate for most unsymmetric matrices.
 *  
 *  A: The users input matrix
 * 
 *  option: Indicates the way in which the weights of each node are selected. It can take 
 *  the following values:
 *      0: w[k] = sum A(i,k)*i
 *      1: w[k] = sum A(:,k)
 *      2: w[k] = mean (A(:,k))
 *      3: w[k] = min A(:,k)
 *      4: w[k] = max A(:,k)
 *      5: w[k] = A(k,k), if A(k,k) = 0 sum.
 * 
 *  alpha: A parameter which is used in the ordering to help determine pivotal nodes. WAMF
 *  selects the node to eliminate based on a convex combination of sparsity and entry size
 *  criterion. alpha is the weight assigned to sparsity. It is recommended that this weight
 *  be larger than 0.5
 * 
 *  aggressive: A boolean which determines whether to be aggressive with dense columns. If
 *  aggressive = true, then dense columns are eliminated from the graph prior to ordering.
 *  If aggressive = false, then dense columns are not eliminated.
 * 
 */

int64_t* SPEX_WAMF_wamf
(
    int64_t order,   
    SPEX_matrix *A,   
    int64_t option,
    double alpha,
    bool aggressive
)
{
    SPEX_matrix *C, *A2, *AT ;
    int64_t *Cp, *Ci, *last, *W, *len, *nv, *next, *P, *head, *elen, *degree, *w,
        *hhead, *ATp, *ATi, d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
        k2, k3, jlast, ln, dense, nzmax, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
        ok, cnz, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q, n, m, t ;
    int64_t h ;
    
    //------------------------------------------------------------------------------
    // Check all inputs
    //------------------------------------------------------------------------------
    
    if (order < 0 || order > 3 || option > 5 || option < 0) 
        return NULL;
    if (!A || !A->p || !A->x.fp64 || !A->i)
        return NULL;
    
    if (alpha < 0 || alpha > 1)
        alpha = SPEX_WAMF_DEFAULT_ALPHA;
    
    
    //------------------------------------------------------------------------------
    // Construct the workspace matrix C. C is a workspace matrix which is either
    // a copy of A, A + A' or A'*A 
    //------------------------------------------------------------------------------
    
    SPEX_transpose (&AT, A) ;              /* compute A' */
    if (!AT) return (NULL) ;
    m = A->m ; n = A->n ;
    
    // Compute a dense threshold. The dense threshold is used to automatically exclude
    // columns whose nnz > dense if aggressive = true. By default, we utilize 10 * sqrt(n)    
    dense = SPEX_MAX (16, 10 * (int64_t) sqrt ((double) n)) ;   /* find dense threshold */
    dense = SPEX_MIN (n-2, dense) ;
    
    // order = 0, set C to be a copy of A    
    if (order == 0)     // Matrix is already SPD, make a copy
    {
        SPEX_matrix_allocate(&C, SPEX_CSC, SPEX_FP64, A->m, A->n, A->nzmax,
        false, true, NULL);
        for (i = 0; i <=n; i++)
        {
            C->p[i] = A->p[i];
        }
        for (i = 0; i < A->nzmax; i++)
        {
            C->x.fp64[i] = A->x.fp64[i];
            C->i[i] = A->i[i];
        }
        C->nz = A->nz;
    }
    
    // Order = 1, set C to be A+A'
    else if (order == 1 && n == m)
    {
        C = SPEX_WAMF_symbolic_add(A, AT);
    }
    
    // Order = 2, set C = A'*A
    else if (order == 2)
    {
        C = SPEX_WAMF_symbolic_multiply (AT, A);
    }
    
    // Free AT
    SPEX_matrix_free(&AT, NULL);
   
    // Ensure that C was properly allocated
    if (!C) return (NULL) ;
    
    // Get the weights associated with each node
    double* weight = SPEX_WAMF_get_weights(C, option);    
    
    // Drop the diagonal entries of C from the graph
    SPEX_WAMF_fkeep (C, &SPEX_WAMF_diag, NULL) ;
    
    // Create workspace, result, and expand C
    Cp = C->p ;
    cnz = Cp [n] ;
    P = (int64_t*) SPEX_malloc( (n+1)*sizeof(int64_t));       /* allocate result */
    W = (int64_t*) SPEX_malloc( (8*(n+1)) * sizeof(int64_t));   /* get workspace */
    t = cnz + cnz/5 + 2*n ;                 /* add elbow room to C */
   
    SPEX_WAMF_sparse_realloc(C);
    
    if (!P || !W || !C) return NULL;
    
    // Create various workspace vectors
    len  = W           ; nv     = W +   (n+1) ; next   = W + 2*(n+1) ;
    head = W + 3*(n+1) ; elen   = W + 4*(n+1) ; degree = W + 5*(n+1) ;
    w    = W + 6*(n+1) ; hhead  = W + 7*(n+1) ;
    last = P ;                              /* use P as workspace for last */
    
    //------------------------------------------------------------------------------
    // Initialize the quotient graph of the matrix. This entails the following steps:
    //
    //      1) Get the size of each column of C, len(k) = nnz(C(:,k)). We then
    //         initialize the degree of each node as deg(k) = len(k)
    //
    //      2) Initialize all linked lists. 
    //
    //      3) Create and initialize two additional arrays for AMF. score is the 
    //         computed AMF score of each node. cliques(i) is the size of the most 
    //         recent clique node i was in.
    //
    //      4) Create additional arrays for WAMF which are:
    //          a) updated: Used when updating the size of each node when simulating
    //             IPGE
    //          b) sizes: Used when approximating the size of each pivot when 
    //             simulating IPGE. Notice that this parameter is different than
    //             the w array because its based on the size when a pivot is selected
    //          c) history: contains the pivot which was the last time an entry was 
    //             updated. Used to simulate history updates.
    //          d) w_score: Computes the WAMF score at each iteration. This is different
    //             from the AMF score array discussed in item 3.
    //------------------------------------------------------------------------------
    
    // Initialize degree list of each node
    for (k = 0 ; k < n ; k++)
    {
        len [k] = Cp [k+1] - Cp [k] ;
    }
    len [n] = 0 ;
    nzmax = C->nzmax ;
    Ci = C->i ;
    
    // Initialize all linked lists
    for (i = 0 ; i <= n ; i++)
    {
        head [i] = -1 ;                     /* degree list i is empty */
        last [i] = -1 ;
        next [i] = -1 ;
        hhead [i] = -1 ;                    /* hash list i is empty */
        nv [i] = 1 ;                        /* node i is just one node */
        w [i] = 1 ;                         /* node i is alive */
        elen [i] = 0 ;                      /* Ek of node i is empty */
        degree [i] = len [i] ;              /* degree of node i */
    }
    mark = SPEX_WAMF_wclear (0, 0, w, n) ;         /* clear w */
    elen [n] = -2 ;                         /* n is a dead element */
    Cp [n] = -1 ;                           /* n is a root of assembly tree */
    w [n] = 0 ;                             /* n is a dead element */
    
        
    // Initialize AMF score array
    double* score = SPEX_calloc(n, sizeof(double));
    // Initialize valid array. valid[k] = 0 if it can be selected as a pivot
    int64_t* valid = SPEX_calloc(n, sizeof(int64_t));
    // Initialize updated array. updated[k] = 0 if it has not been updated
    int64_t* updated = SPEX_calloc(n, sizeof(int64_t));
    // Initialize cliques array. cliques[k] is the size of the most recent clique
    // which node k appears in
    double* cliques = SPEX_calloc(n, sizeof(double));
    // Initialize sizes array. size[k] is the size of eliminating node k
    double* sizes = SPEX_calloc(n, sizeof(double));
    // initialize history. history[k] is the last pivot node k was updated WRT
    double* history = SPEX_calloc(n, sizeof(double));
    // initialize w_score which is the WAMF score of node k
    double* w_score = SPEX_calloc(n, sizeof(double));
        
    // Extra doubles for WAMF
    double minscore, prev_piv = 0, mins = 0, min, maxweigh = 0, maxscore = 0;
    
        
    //------------------------------------------------------------------------------
    // First pass and initialization of the degree lists. At this step, we begin to
    // analyze the graph of A. If aggressive = true, we will automatically eliminate 
    // any node of degree 0 or degree greater than dense. 
    // Then, the linked lists are initialized for each node in the graph.
    //------------------------------------------------------------------------------
    for (i = 0 ; i < n ; i++)
    {
        d = degree [i] ;
         if (d == 0 && aggressive == true)     /* node i is empty */
         {
             valid[i] = -1;
             elen [i] = -2 ;                 /* element i is dead */
             nel++ ;
             Cp [i] = -1 ;                   /* i is a root of assembly tree */
             w [i] = 0 ;
         }
         else if (d > dense && aggressive == true)  /* node i is dense */
         {
             valid[i] = -1;
             nv [i] = 0 ;                    /* absorb i into element n */
             elen [i] = -1 ;                 /* node i is dead */
             nel++ ;
             Cp [i] = SPEX_FLIP (n) ;
             nv [n]++ ;
         }
         else
         {
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;           /* put node i in degree list d */
            head [d] = i ;
         }
    }

    //------------------------------------------------------------------------------
    // Begin WAMF ordering
    //------------------------------------------------------------------------------
    while (nel < n)                         /* while (selecting pivots) do */
    {
        //--------------------------------------------------------------------------
        // Compute the AMF score and sizes of each node. The AMF score is done as
        // described by Rothberg and Eisenstat. The size of each node is computed
        // by performing a history update to each node.
        //--------------------------------------------------------------------------
        for (k = 0; k < n; k++)
        {
            // score[k] = d[k]^2 - d[k]
            score[k] = degree[k]*degree[k]-degree[k];
            // score[k] = score[k] - ( c[k]^2 - c[k])
            score[k] -= cliques[k]*cliques[k] - cliques[k];
            // score[k] = score[k] / 2
            score[k] = score[k]/2;
            
            // score[k] must be >= 0
            if (score[k] < 0)
                score[k] = 0;
                    
            ASSERT(score[k] >= 0);
            // Compute the size of each node. The size of each node is given by
            // performing a history update to each node.
            sizes[k] = weight[k] + prev_piv - history[k];
            
            // Sizes must be >=0
            ASSERT(sizes[k] >= 0);
        }
        
        
        //--------------------------------------------------------------------------
        // At this point, each node has a size and score. However, these parameters 
        // are not comparable as they can be an order of magnitude apart. In order
        // to fix this, we standardize the measurements by finding the maximum values
        // of score and sizes.
        //--------------------------------------------------------------------------
        for (k = 0; k < n; k++)
        {
            if (valid[k] == 0)
            {
                // Obtain maximum size
                if (sizes[k] > maxweigh)
                {
                    maxweigh = sizes[k];
                }
                // Obtain maximum score
                if (score[k] > maxscore)
                {
                    maxscore = score[k];
                }
            }
        }
        
        //--------------------------------------------------------------------------
        // We now have a score and size associated with each node. We also have the 
        // maximum weight and maximum score of each node. Note that for all k,
        // 0 <= score[k] / maxscore <= 1, and 0 <= sizes[k] / maxweigh <= 1. We will
        // use this to compute the w_score of each node. The w_score is a convex 
        // combination of these standardized values.
        //
        // Note that the parameter alpha is applied to the score (sparsity) criterion.
        // If alpha is within the recommended range, this means more weight is given 
        // towards the sparsity criterion.
        //
        // Once we compute the w_score associated with each node, we eliminate the 
        // node of lowest w_score from the graph.
        //--------------------------------------------------------------------------
        for (k  = 0; k < n; k++)
        {
            w_score[k] = alpha * score[k] / maxscore + (1-alpha) * sizes[k] / maxweigh;           
        }
        
        // Find a valid pivot element
        int64_t flag = 0;
        j = 0;
        while (flag == 0)
        {
            if (valid[j] == 0)
            {
                k = j;
                min = w_score[k];   // Current min score
                flag = 1;
            }
            else
                j+=1;            
        }
        
        // Now that we have a valid pivot, find the one of min score
        for ( ; j < n; j++)
        {
            if ( w_score[j] < min && valid[j] == 0)
            {
                min = w_score[j];
                k = j;
            }
        }
        
        // Node k is no longer valid
        valid[k] = -1;
        
        // mindegree is the degree of node k (for linked lists)
        mindeg = degree[k];
        
        //--------------------------------------------------------------------------
        // Select the node from linked lists
        //--------------------------------------------------------------------------
        
        if (next [k] != -1) last [next [k]] = -1 ;
        head [mindeg] = next [k] ;          /* remove k from degree list */
        elenk = elen [k] ;                  /* elenk = |Ek| */
        nvk = nv [k] ;                      /* # of nodes k represents */
        nel += nvk ;                        /* nv[k] nodes of A eliminated */
        //--------------------------------------------------------------------------
        // Garbage collection
        //--------------------------------------------------------------------------
        if (elenk > 0 && cnz + mindeg >= nzmax)
        {
            for (j = 0 ; j < n ; j++)
            {
                if ((p = Cp [j]) >= 0)      /* j is a live node or element */
                {
                    Cp [j] = Ci [p] ;       /* save first entry of object */
                    Ci [p] = SPEX_FLIP (j) ;  /* first entry is now CS_FLIP(j) */
                }
            }
            for (q = 0, p = 0 ; p < cnz ; ) /* scan all of memory */
            {
                if ((j = SPEX_FLIP (Ci [p++])) >= 0)  /* found object j */
                {
                    Ci [q] = Cp [j] ;       /* restore first entry of object */
                    Cp [j] = q++ ;          /* new pointer to object j */
                    for (k3 = 0 ; k3 < len [j]-1 ; k3++) Ci [q++] = Ci [p++] ;
                }
            }
            cnz = q ;                       /* Ci [cnz...nzmax-1] now free */
        }
        //--------------------------------------------------------------------------
        // Construct new element
        //--------------------------------------------------------------------------
        dk = 0 ;
        nv [k] = -nvk ;                     /* flag k as in Lk */
        p = Cp [k] ;
        pk1 = (elenk == 0) ? p : cnz ;      /* do in place if elen[k] == 0 */
        pk2 = pk1 ;
        for (k1 = 1 ; k1 <= elenk + 1 ; k1++)
        {
            if (k1 > elenk)
            {
                e = k ;                     /* search the nodes in k */
                pj = p ;                    /* list of nodes starts at Ci[pj]*/
                ln = len [k] - elenk ;      /* length of list of nodes in k */
            }
            else
            {
                e = Ci [p++] ;              /* search the nodes in e */
                pj = Cp [e] ;
                ln = len [e] ;              /* length of list of nodes in e */
            }
            for (k2 = 1 ; k2 <= ln ; k2++)
            {
                i = Ci [pj++] ;
                if ((nvi = nv [i]) <= 0) continue ; /* node i dead, or seen */
                dk += nvi ;                 /* degree[Lk] += size of node i */
                nv [i] = -nvi ;             /* negate nv[i] to denote i in Lk*/
                Ci [pk2++] = i ;            /* place i in Lk */
                if (next [i] != -1) last [next [i]] = last [i] ;
                if (last [i] != -1)         /* remove i from degree list */
                {
                    next [last [i]] = next [i] ;
                }
                else
                {
                    head [degree [i]] = next [i] ;
                }
            }
            if (e != k)
            {
                Cp [e] = SPEX_FLIP (k) ;      /* absorb e into k */
                w [e] = 0 ;                 /* e is now a dead element */
            }
        }
        if (elenk != 0) cnz = pk2 ;         /* Ci [cnz...nzmax] is free */
        degree [k] = dk ;                   /* external degree of k - |Lk\i| */
        Cp [k] = pk1 ;                      /* element k is in Ci[pk1..pk2-1] */
        len [k] = pk2 - pk1 ;
        elen [k] = -2 ;                     /* k is now an element */
        //--------------------------------------------------------------------------
        // Compute the set differences
        //--------------------------------------------------------------------------
        mark = SPEX_WAMF_wclear (mark, lemax, w, n) ;  /* clear w if necessary */
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan 1: find |Le\Lk| */
        {
            i = Ci [pk] ;
            if ((eln = elen [i]) <= 0) continue ;/* skip if elen[i] empty */
            nvi = -nv [i] ;                      /* nv [i] was negated */
            wnvi = mark - nvi ;
            for (p = Cp [i] ; p <= Cp [i] + eln - 1 ; p++)  /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] >= mark)
                {
                    w [e] -= nvi ;          /* decrement |Le\Lk| */
                }
                else if (w [e] != 0)        /* ensure e is a live element */
                {
                    w [e] = degree [e] + wnvi ; /* 1st time e seen in scan 1 */
                }
            }
        }
        //--------------------------------------------------------------------------
        // Update the degree of neighboring nodes
        //--------------------------------------------------------------------------
        for (pk = pk1 ; pk < pk2 ; pk++)    /* scan2: degree update */
        {
            i = Ci [pk] ;                   /* consider node i in Lk */
            p1 = Cp [i] ;
            p2 = p1 + elen [i] - 1 ;
            pn = p1 ;
            for (h = 0, d = 0, p = p1 ; p <= p2 ; p++)    /* scan Ei */
            {
                e = Ci [p] ;
                if (w [e] != 0)             /* e is an unabsorbed element */
                {
                    dext = w [e] - mark ;   /* dext = |Le\Lk| */
                    if (dext > 0)
                    {
                        d += dext ;         /* sum up the set differences */
                        Ci [pn++] = e ;     /* keep e in Ei */
                        h += e ;            /* compute the hash of node i */
                    }
                    else
                    {
                        Cp [e] = SPEX_FLIP (k) ;  /* aggressive absorb. e->k */
                        w [e] = 0 ;             /* e is a dead element */
                    }
                }
            }
            elen [i] = pn - p1 + 1 ;        /* elen[i] = |Ei| */
            p3 = pn ;
            p4 = p1 + len [i] ;
            for (p = p2 + 1 ; p < p4 ; p++) /* prune edges in Ai */
            {
                j = Ci [p] ;
                if ((nvj = nv [j]) <= 0) continue ; /* node j dead or in Lk */
                d += nvj ;                  /* degree(i) += |j| */
                Ci [pn++] = j ;             /* place j in node list of i */
                h += j ;                    /* compute hash for node i */
            }
            degree [i] = SPEX_MIN (degree [i], d) ;   /* update degree(i) */
            Ci [pn] = Ci [p3] ;         /* move first node to end */
            Ci [p3] = Ci [p1] ;         /* move 1st el. to end of Ei */
            Ci [p1] = k ;               /* add k as 1st element in of Ei */
            len [i] = pn - p1 + 1 ;     /* new len of adj. list of node i */
            h = ((h<0) ? (-h):h) % n ;  /* finalize hash of i */
            next [i] = hhead [h] ;      /* place i in hash bucket */
            hhead [h] = i ;
            last [i] = h ;              /* save hash of i in last[i] */
        }                                   /* scan2 is done */
        degree [k] = dk ;                   /* finalize |Lk| */
        lemax = SPEX_MAX (lemax, dk) ;
        mark = SPEX_WAMF_wclear (mark+lemax, lemax, w, n) ;    /* clear w */
        //--------------------------------------------------------------------------
        // Finalize the new element
        //--------------------------------------------------------------------------
        for (p = pk1, pk = pk1 ; pk < pk2 ; pk++)   /* finalize Lk */
        {
            i = Ci [pk] ;
            if ((nvi = -nv [i]) <= 0) continue ;/* skip if i is dead */
            nv [i] = nvi ;                      /* restore nv[i] */
            d = degree [i] + dk - nvi ;         /* compute external degree(i) */
            d = SPEX_MIN (d, n - nel - nvi) ;
            if (head [d] != -1) last [head [d]] = i ;
            next [i] = head [d] ;               /* put i back in degree list */
            last [i] = -1 ;
            head [d] = i ;
            mindeg = SPEX_MIN (mindeg, d) ;       /* find new minimum degree */
            degree [i] = d ;
            Ci [p++] = i ;                      /* place i in Lk */
            /* Set the size of clique [i] = mindeg - 1 */
            cliques[i] = SPEX_MAX(mindeg, 1);
            
            //----------------------------------------------------------------------
            // We now update the weight of each node, i, adjacent to the eliminated node. 
            // This is done via a symbolic version of SLIP LU's triangular solve.
            //
            // First, we perform a history update to node i by additing the size of
            // the previous pivot and subtracting the size of its history (which mimics
            // the multiplication by p^k-1 and division by p^hi.
            //
            // Next, we simulate the IPGE update on A. The for an entry aij, the general
            // IPGE update is given as:  ( p^k * aij - aik * ajk) / p^k-1.
            // In this algorithm, we are not storing the value of aik and ajk, thus we
            // cannot perform an exact approximation of this operation. Instead, we 
            // approximate the size of aik and ajk to be the same size of the pivot element.
            // Thus, in essence we are computing the numerator as (p^k*aij-p^k*p^k).
            // Now, note that the number of bits in the numerator is upper bounded as 
            // the maximum number of bits in the left term or right term. Thus, 
            // to simulate the IPGE update on node i, we take the max of w[i]+w[k] and 
            // 2*w[k]. Then, we subtract the value of the previous pivot. 
            //
            // Finally, we update the two book keeping updated and history arrays.
            //----------------------------------------------------------------------
            if (updated[i] == 0 && i!= k)
            {
                // Simulate IPGE update on i
                weight[i] += prev_piv - history[i];
                weight[i] = SPEX_MAX(weight[i]+weight[k], 2*weight[k]);
                weight[i] -= prev_piv;
                
                // node i has been updated
                updated[i] = -1;
                // node i's history is w[k]
                history[i] = weight[k];
            }
            
        }
        nv [k] = nvk ;                      /* # nodes absorbed into k */
        if ((len [k] = p-pk1) == 0)         /* length of adj list of element k*/
        {
            Cp [k] = -1 ;                   /* k is a root of the tree */
            w [k] = 0 ;                     /* k is now a dead element */
        }
        if (elenk != 0) cnz = p ;           /* free unused space in Lk */
        // Set the previous pivot to be node k
        prev_piv = weight[k];
        // Reset updated
        for (j = 0; j < n; j++) 
            updated[j] = 0;
    }
    /* --- Postordering ----------------------------------------------------- */
    for (i = 0 ; i < n ; i++) Cp [i] = SPEX_FLIP (Cp [i]) ;/* fix assembly tree */
    for (j = 0 ; j <= n ; j++) head [j] = -1 ;
    for (j = n ; j >= 0 ; j--)              /* place unordered nodes in lists */
    {
        if (nv [j] > 0) continue ;          /* skip if j is an element */
        next [j] = head [Cp [j]] ;          /* place j in list of its parent */
        head [Cp [j]] = j ;
    }
    for (e = n ; e >= 0 ; e--)              /* place elements in lists */
    {
        if (nv [e] <= 0) continue ;         /* skip unless e is an element */
        if (Cp [e] != -1)
        {
            next [e] = head [Cp [e]] ;      /* place e in list of its parent */
            head [Cp [e]] = e ;
        }
    }
    for (k = 0, i = 0 ; i <= n ; i++)       /* postorder the assembly tree */
    {
        if (Cp [i] == -1) k = SPEX_WAMF_tdfs (i, k, head, next, P, w) ;
    }
    SPEX_matrix_free(&C, NULL);
    SPEX_FREE(W);
    SPEX_FREE(weight);
    SPEX_FREE(score);
    SPEX_FREE(valid);
    SPEX_FREE(updated);
    SPEX_FREE(cliques);
    SPEX_FREE(sizes);
    SPEX_FREE(history);
    SPEX_FREE(w_score);
    return P;
}

