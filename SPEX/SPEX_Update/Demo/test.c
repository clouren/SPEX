#include "gmp.h"
#include <stdio.h>
		#include <time.h>
void main()
{
    mpz_t a, b, c, r, a_copy, a1;
    mpq_t q, q1;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(r);
    mpz_set_ui(a,2697);
    mpz_set_ui(b,49);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);
    mpz_set_si(a,-2697);
    mpz_set_ui(b,49);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);
    mpz_set_ui(a,2697);
    mpz_set_si(b,-49);
    mpz_fdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd      ",a,b,c,r);
    mpz_fdiv_q(a,a,b);
    gmp_printf("%Zd...?      ",a);
    mpz_set_ui(a,2697);
    mpz_cdiv_qr(c,r,a,b);
    gmp_printf("c=a/b=%Zd/%Zd=%Zd...%Zd\n",a,b,c,r);

    a1->_mp_d = NULL;
    a1->_mp_size = 0;
    a1->_mp_alloc = 0;
    gmp_printf("before swap %Zd\n",a1);
    mpz_swap(a,a1);
    gmp_printf("after swap %Zd %p\n",a1,a1->_mp_d);
    gmp_printf("after swap %Zd %p\n",a,a->_mp_d);
    //mpz_set_si(a,5);
    //gmp_printf("set a = 5 and get a=%Zd\n",a);

    mpq_init(q);
    mpq_set_si(q,2, 3);
    mpq_get_num(a,q);
    mpq_get_den(b,q);
    gmp_printf("%Zd/%Zd %Zd/%Zd\n",mpq_numref(q),mpq_denref(q),a,b);
    mpq_set_num(q,a);
    mpq_set_den(q,c);
    gmp_printf("%Zd/%Zd %Zd/%Zd\n",mpq_numref(q),mpq_denref(q),a,b);
    mpq_canonicalize(q);
    mpq_get_num(a,q);
    mpq_get_den(b,q);
    gmp_printf("%Zd/%Zd %Zd/%Zd\n",mpq_numref(q),mpq_denref(q),a,b);
    mpq_set_si(q,-2, -3);
    mpq_get_num(a,q);
    mpq_get_den(b,q);
    gmp_printf("%Zd/%Zd %Zd/%Zd\n",mpq_numref(q),mpq_denref(q),a,b);
    mpq_set_si(q,-2, 3);
    mpq_get_num(a,q);
    mpq_get_den(b,q);
    gmp_printf("%Zd/%Zd %Zd/%Zd\n",mpq_numref(q),mpq_denref(q),a,b);

 /*   int t[10]= {1,2,3,3,5,7,8,9,3,7};
    int tt[10]={2,4,1,3,1,6,3,6,5,8};
    int max = 10;
    int i;
    for(i = 0;i<max;i++)
    {
        if (tt[i]>t[i])
	{
		printf("yeah i=%d, max = %d\n",i,max);
		max--;
	}
    }
	printf("yeah i=%d, max = %d\n",i,max);

	clock_t start, end;
	double t_fdiv=0, t_cdiv=0, t_tdiv=0, t_divex;

	mpz_set_ui(a,23516871546525);
	mpz_mul(a,a,a);
	mpz_mul(a,a,a);
	mpz_set_si(b,2164466);
	mpz_mul(b,b,b);
	mpz_mul(b,b,b);
	mpz_mul(a,a,b);
	for (int iter=0;iter<100;iter++)
	{
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_cdiv_q(c,a,b);
		}
		end = clock();
		t_cdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_tdiv_q(c,a,b);
		}
		end = clock();
		t_tdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_fdiv_q(c,a,b);
		}
		end = clock();
		t_fdiv += ((double) (end - start)) / CLOCKS_PER_SEC;
		start = clock();
		for(i=0;i<100000;i++)
		{
			mpz_divexact(c,a,b);
		}
		end = clock();
		t_divex += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	printf("time for fdiv: %f, cdiv:%f, tdiv:%f, divex:%f\n",t_fdiv, t_cdiv, t_tdiv, t_divex);*/
}
