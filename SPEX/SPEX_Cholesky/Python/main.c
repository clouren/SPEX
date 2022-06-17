#include <stdio.h>
#include "SPEX_Chol_connect.h"
//TODO rewrite this to make it into a demo for the python interfase
//TODO I still need the makefile :(

int main(){
 
     double* x=NULL;                 
     int64_t* col_pointers=NULL;    
     int64_t* row_index=NULL;        
     double* data=NULL;             
     double* rhs=NULL;            
     int n;                  
     int nz;
     
     n=3;
     nz=9;
     int64_t cp[4]={0, 3, 6, 9};
     int64_t ri[9]={0, 1, 2, 0, 1, 2, 0, 1, 2};
     double d[9]={4,  12, -16,  12,  37, -43, -16, -43,  98};
     double r[3]={1,1,1};
     
     col_pointers=cp;
     row_index=ri;
     data=d;
     rhs=r;
     

     int res;

     //res=python_backslash(&x,col_pointers,row_index,data,rhs,n,nz);
     //TOASK: freeing memory?? double free or corruption (out)

     //printf("%f,%f,%f\n",x[0],x[1],x[2]);
     //free(x);
     //free(col_pointers);
     //free(row_index);
     //free(data);
     //free(rhs);
/*
     int size=4;
     double* input=NULL;
     input=(double *)malloc(size * sizeof(double));
*/
     /*read(input,size);
     for (int i = 0; i < size; ++i)
     {
          printf("%f,", input[i]);
     }
     printf("\n");*/
    /*
     double* x2=NULL;
     x2=(double *)malloc(n*sizeof(double));
     for (int i = 0; i < n; ++i)
    {
        x2[i]=0;
    }
     python_backslash_double(x2,col_pointers,row_index,data,rhs,n,nz,2); //passing by reference on python is weird and complicated
     printf("%f,%f,%f\n",x2[0],x2[1],x2[2]);
*/
    /*
    mpz_t var;
    mpz_init(var);
    mpz_set_ui(var, 40005442);
    char* varS;
    varS=(char *)malloc(n*sizeof(char));
    varS=mpz_get_str(NULL,10,var);
    printf("%s\n",varS);
    */
    /*char** x2=NULL;
    x2=(char **)malloc(n*sizeof(char *));

    
    python_backslash_char(x2,col_pointers,row_index,data,rhs,n,nz,0); 
    printf("allfine\n");
     printf("%s,%s,%s\n",x2[0],x2[1],x2[2]);
     
     free(x2);*/

     double* x_d=NULL;
     x_d=(double *)malloc(n*sizeof(double));
     char** x_c=NULL;
     x_c=(char **)malloc(n*sizeof(char *));
     bool charOut = false;
     //printf("%d \n", charOut);

     void** x_v=NULL;
     x_v=(void **)malloc(n*sizeof(void *));


     SPEX_python_backslash(x_c,x_d,x_v,charOut,col_pointers,row_index,data,rhs,n,n,nz,2);
     if(charOut)
     {
          printf("STRING SOLUTION\n");
          printf("%s,%s,%s\n",x_c[0],x_c[1],x_c[2]);
          printf("void %s,%s,%s\n",x_v[0],x_v[1],x_v[2]);
     }
     else
     {
          printf("DOUBLE SOLUTION\n");
          printf("%f,%f,%f\n",x_d[0],x_d[1],x_d[2]);
          //printf("void %f,%f,%f\n",*(double *)x_v[0],*(double *)x_v[1],*(double *)x_v[2]);
          printf("void %f,%f,%f\n",x_v[0],x_v[1],x_v[2]);
     }
     free(x_d);
     free(x_c);
     free(x_v);
     return 0;
}

