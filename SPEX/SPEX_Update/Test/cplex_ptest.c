//------------------------------------------------------------------------------
// SPEX_Update/Test/cplex_ptest.c: performance test for SPEX_Update library
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/*
 * performance test using CPLEX library
 */

#include "cplex_test.h"

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    //char *prob_name = "lp_25fv47";//5.5018458883E+03
    char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize for SPEX libray
    //------------------------------------------------------------------

    SPEX_initialize () ;

    int CPLEX_status = 0;
    CPXENVptr env = NULL;
    CPXLPptr LP = NULL;
    SPEX_options* option = NULL;
    SPEX_matrix *MA = NULL, *Mb = NULL, *Mc = NULL;
    /*SPEX_matrix *Prob_A = NULL, *Prob_c = NULL, *tempA = NULL;
    SPEX_matrix *L1 = NULL, *U1 = NULL, *rhos = NULL, *A1 = NULL, *x1 = NULL;
    SPEX_matrix *rhos2 = NULL, *rhos3 = NULL;
    SPEX_matrix *L2 = NULL, *U2 = NULL, *A2 = NULL;
    SPEX_matrix *b = NULL, *c = NULL,  *x2 = NULL, *y = NULL;
    SPEX_matrix *L3 = NULL, *U3 = NULL, *A3 = NULL;
    SPEX_vector *tmpv, *vk = NULL, *vk3 = NULL;
    mpz_t z, minz, tmpz;
    mpq_t minq, tmpq, obj_mpq;
    SPEX_LU_analysis* analysis = NULL;
    int64_t *P1_inv = NULL, *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
    int64_t *P3 = NULL, *P3_inv = NULL, *Q3 = NULL, *Q3_inv = NULL;
    int64_t *basis = NULL, *h = NULL, *used_as_basis = NULL;
    double *glpk_rhs = NULL, *col_val = NULL;
    int *col_ind = NULL, *row_ind = NULL;
    clock_t start_llu, start_luu, end_llu, end_luu, start1, end1, start2, end2,
            start_luu3 = 0, end_luu3 = 0;
    FILE *out_file = NULL;*/
    char file_name[1000];
    OK(SPEX_create_default_options(&option));

    env = CPXopenCPLEX (&CPLEX_status);
    if (env == NULL)
    {
        printf ("Could not open CPLEX environment.\n");
        return 0;
    }
    sprintf(file_name, "TestMats/%s/%s",prob_name,prob_name);

    OK(SPEX_get_CPLEX_LP(env, &LP, &CPLEX_status, &MA, &Mb, &Mc, file_name,
        option));
    printf("CPLEX status=%d\n", CPLEX_status);

   CPLEX_status = CPXsetintparam (env, CPXPARAM_Simplex_Display, 2);
   if ( CPLEX_status ) {
      fprintf (stderr,"Failed to turn up simplex display level.\n");
   }
    int nvars = CPXgetnumcols(env, LP);
    int neqs = CPXgetnumrows(env, LP);
    //int nnz = CPXgetnumnz(env,LP);
    double objval;

    CPLEX_status = CPXsetintparam (env, CPXPARAM_LPMethod, CPX_ALG_DUAL);
    printf("CPLEX status=%d\n", CPLEX_status);
    CPLEX_status = CPXlpopt(env, LP);
    printf("CPLEX status=%d\n", CPLEX_status);
   CPLEX_status = CPXgetobjval (env, LP, &objval);
    printf("CPLEX status=%d\n", CPLEX_status);

    CPLEX_status = CPXpresolve(env, LP, CPX_ALG_PRIMAL);
    printf("presolve CPLEX status=%d\n", CPLEX_status);

    CPXLPptr redLP = NULL;
    CPLEX_status = CPXgetredlp(env, LP, (CPXCLPptr*)&redLP);
    printf("get presolve status=%d\n", CPLEX_status);


    nvars = CPXgetnumcols(env, LP);
    neqs = CPXgetnumrows(env, LP);
    printf("size is %dx%d\n",neqs,nvars);

    int *cstat = SPEX_malloc(nvars*sizeof(int));
    int *rstat = SPEX_malloc(neqs *sizeof(int));
    int *head = (int*) SPEX_malloc(neqs*sizeof(int));
    double *x = (double*) SPEX_malloc(neqs*sizeof(double));

    int count[4] = {0, 0, 0,0};
    CPLEX_status = CPXgetbase(env, LP, cstat, rstat);
    printf("get basis status=%d\n", CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        double rhs_i = 0;
        CPLEX_status = CPXgetrhs(env, LP, &rhs_i, i, i);
        if (CPLEX_status) printf("get rhs error =%d\n",CPLEX_status);
        printf("r[%d]=%d, rhs=%f\n",i, rstat[i],rhs_i);
        count[rstat[i]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;
    for (int j = 0 ; j < nvars; j++)
    {
        double lb = 0, ub = 0;
        CPLEX_status = CPXgetlb(env, LP, &lb, j, j);
        if (CPLEX_status) printf("get lb error =%d\n",CPLEX_status);
        CPLEX_status = CPXgetub(env, LP, &ub, j, j);
        if (CPLEX_status) printf("get ub error =%d\n",CPLEX_status);
        printf("c[%d]=%d, [%f %f]\n", j,cstat[j],lb,ub);
        count[cstat[j]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;

    nvars = CPXgetnumcols(env, redLP);
    neqs = CPXgetnumrows(env, redLP);
    printf("size of redLP is %dx%d\n",neqs,nvars);
    CPLEX_status = CPXgetbase(env, redLP, cstat, rstat);
    printf("get basis status=%d\n", CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        double rhs_i = 0;
        CPLEX_status = CPXgetrhs(env, redLP, &rhs_i, i, i);
        if (CPLEX_status) printf("get rhs error =%d\n",CPLEX_status);
        printf("r[%d]=%d, rhs=%f\n",i, rstat[i],rhs_i);
        count[rstat[i]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;
    for (int j = 0 ; j < nvars; j++)
    {
        double lb = 0, ub = 0;
        CPLEX_status = CPXgetlb(env, redLP, &lb, j, j);
        if (CPLEX_status) printf("get lb error =%d\n",CPLEX_status);
        CPLEX_status = CPXgetub(env, redLP, &ub, j, j);
        if (CPLEX_status) printf("get ub error =%d\n",CPLEX_status);
        printf("c[%d]=%d, [%f %f]\n", j,cstat[j],lb,ub);
        count[cstat[j]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);

    CPLEX_status = CPXgetbhead(env, LP, head, x);
    if (CPLEX_status) printf("get bhead error =%d\n",CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        printf("%d ",head[i]);
    }
    printf("\n");

    //--------------------------------------------------------------------------
    CPLEX_status = CPXprimopt(env, LP);
    printf("simplex CPLEX status=%d\n", CPLEX_status);

    CPLEX_status = CPXgetobjval (env, LP, &objval);
    printf("presplve status=%d\n", CPLEX_status);
    printf ("After network optimization, objective is %.10g\n", objval);

    CPLEX_status = CPXprimopt(env, redLP);
    printf("simplex CPLEX status=%d\n", CPLEX_status);

    CPLEX_status = CPXgetobjval (env, redLP, &objval);
    printf("presplve status=%d\n", CPLEX_status);
    printf ("After network optimization, objective is %.10g\n", objval);


    nvars = CPXgetnumcols(env, LP);
    neqs = CPXgetnumrows(env, LP);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;
    CPLEX_status = CPXgetbase(env, LP, cstat, rstat);
    printf("get basis status=%d\n", CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        double rhs_i = 0;
        CPLEX_status = CPXgetrhs(env, LP, &rhs_i, i, i);
        if (CPLEX_status) printf("get rhs error =%d\n",CPLEX_status);
        printf("r[%d]=%d, rhs=%f\n",i, rstat[i],rhs_i);
        count[rstat[i]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;
    for (int j = 0 ; j < nvars; j++)
    {
        double lb = 0, ub = 0;
        CPLEX_status = CPXgetlb(env, LP, &lb, j, j);
        if (CPLEX_status) printf("get lb error =%d\n",CPLEX_status);
        CPLEX_status = CPXgetub(env, LP, &ub, j, j);
        if (CPLEX_status) printf("get ub error =%d\n",CPLEX_status);
        printf("c[%d]=%d, [%f %f]\n", j,cstat[j],lb,ub);
        count[cstat[j]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;

    CPLEX_status = CPXgetbhead(env, LP, head, x);
    if (CPLEX_status) printf("get bhead error =%d\n",CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        printf("%d ",head[i]);
    }
    printf("\n");

    nvars = CPXgetnumcols(env, redLP);
    neqs = CPXgetnumrows(env, redLP);
    printf("size of redLP is %dx%d\n",neqs,nvars);
    CPLEX_status = CPXgetbase(env, redLP, cstat, rstat);
    printf("get basis status=%d\n", CPLEX_status);
    for (int i = 0; i < neqs; i++)
    {
        double rhs_i = 0;
        CPLEX_status = CPXgetrhs(env, redLP, &rhs_i, i, i);
        if (CPLEX_status) printf("get rhs error =%d\n",CPLEX_status);
        printf("r[%d]=%d, rhs=%f\n",i, rstat[i],rhs_i);
        count[rstat[i]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    count[0]=0;count[1]=0;count[2]=0;count[3]=0;
    for (int j = 0 ; j < nvars; j++)
    {
        double lb = 0, ub = 0;
        CPLEX_status = CPXgetlb(env, redLP, &lb, j, j);
        if (CPLEX_status) printf("get lb error =%d\n",CPLEX_status);
        CPLEX_status = CPXgetub(env, redLP, &ub, j, j);
        if (CPLEX_status) printf("get ub error =%d\n",CPLEX_status);
        printf("c[%d]=%d, [%f %f]\n", j,cstat[j],lb,ub);
        count[cstat[j]]++; 
    }
    printf("count0=%d count1=%d count2=%d count3=%d\n",count[0],count[1],count[2],count[3]);
    return 0;
}

