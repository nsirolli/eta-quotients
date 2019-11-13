/* This silly code actually proves that a certain Q that
 * fails 4.4 is truly uninteresting.  It just tries to compute
 * p((Q^3*n+1)/24) mod p until we get a non-zero value.  This will
 * eventually come to a crawl so we put a limit on how large n to
 * test.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <flint.h>
#include <arith.h>
#include <fmpz.h>


int main(int argc, char **argv) {
    
    //set up p and Q from the command line
    slong p = atoi(argv[1]);
    slong Q = atoi(argv[2]);
    fmpz_t neg_one, pp, QQ;
    fmpz_init(neg_one);
    fmpz_set_si(neg_one,-1);
    fmpz_init(pp);
    fmpz_init(QQ);
    fmpz_set_si(pp,p);
    fmpz_set_si(QQ,Q);

    slong n;
    for(n=1;n < 10000; n++) {
        //move on if the hypotheses on n are not met
        if(n%24!=1 || n%p==0 || n%Q==0)
            continue;
        fmpz_t negn;
        fmpz_init(negn);
        fmpz_set_si(negn,-1*n);
        //move on if (-n|p) = (-1|p)
        if(fmpz_jacobi(negn,pp)==fmpz_jacobi(neg_one,pp)) {
            fmpz_clear(negn);
            continue;
        }
        //Compute the number of partitions mod p and test
        fmpz_clear(negn);
        slong r = (Q*Q*Q*n+1)/24;
        fmpz_t answer;
        fmpz_init(answer);
        arith_number_of_partitions(answer,(ulong)r);
        fmpz_t mod;
        fmpz_init(mod);
        fmpz_mod(mod,answer,pp);
        fmpz_clear(answer);
        if(fmpz_get_si(mod)%p!=0) {
            flint_printf("Found a failing n of %wd for Q = %wd\n",n,Q);
            break;
        }
    }
    return 0;
}

