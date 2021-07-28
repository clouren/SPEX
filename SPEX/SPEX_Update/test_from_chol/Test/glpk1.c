#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
int main(void)
{
    glp_prob *P;
    glp_smcp parm;
    int m, n;
    P = glp_create_prob();
    glp_read_mps(P, GLP_MPS_DECK, NULL, "25fv47.mps");
    m = glp_get_num_rows(P);
    n = glp_get_num_cols(P);
    printf("nrows=%d ncols=%d\n",m,n);
    glp_adv_basis(P, 0);
    glp_init_smcp(&parm);
    parm.it_lim=1;
    glp_factorize(P);
    printf("nrows=%d ncols=%d\n",m,n);
    int *basis = (int*) malloc(m * sizeof(int));
    if (!basis)  {        return 0;    }
    for (int i =1; i <=m;i++)
    {
        basis[i-1] = glp_get_bhead(P, i);
        if (glp_get_mat_row(P,i,NULL, NULL)==0)
        {
            printf("%d\n",i);
        }
    }
    printf("\n\n\n");
    parm.it_lim=1;
    while (glp_get_status(P) != GLP_OPT)
    {
        if(glp_get_row_stat(P,1)!=1) return 0;
        int num[5]={0,0,0,0,0};
        for (int i = 1; i <= m; i++)
        {
            int stat=glp_get_row_stat(P,i);
            num[stat-1]++;
            if (glp_get_mat_row(P,i,NULL, NULL)==0)
            {
                printf("%d\n",i);
            }
        }
        for (int i = 1; i <= n; i++)
        {
            int stat=glp_get_col_stat(P,i);
            num[stat-1]++;
        }
        printf("n:%d m:%d bs:%d nl:%d nu:%d nf:%d ns:%d\n",n,m,num[0],num[1],num[2],num[3],num[4]);
        glp_simplex(P, &parm);
        //glp_simplex(P, NULL);
        int k = -1;
        for (int i = 1; i <= m; i++)
        {
            int j = glp_get_bhead(P, i);
            if (basis[i-1] != j)
            {
                if (k == -1)
                {
                    k = 1;
                    printf("replace %d-th basis (%d) with %d\n",i,basis[i-1],j);
                    basis[i-1]=j;
                }
                else
                {
                    printf("more than 1 col replacement %d -> %d\n",basis[i-1],j);
                    return 0;
                }
            }
        }
        if (k == -1)
        {
            printf("no column replcement\n");
                    return 0;
        }
        glp_print_sol(P, "25fv47.txt");
    }
    glp_delete_prob(P);
    return 0;
}
