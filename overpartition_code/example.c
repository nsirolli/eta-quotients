/* This program tests (4.4) for the overpartition function
 * It uses the fact that the generating function is 
 * eta(q^2) / eta(q)^2.  This requires a lot of memory to
 * compute because we will need to store large truncated
 * power series
*/

#include <stdlib.h>
#include <stdio.h>
#include <flint.h>
#include <arith.h>


slong power(slong x, ulong y, slong p); //modular exponentiation
slong legendre_symbol(slong x, slong y); //for legendre symbol

int main(void) {

    int p = 3; //p = l in the paper
    ulong sturm = 852;//sturm bound
    ulong qlist[] = {47, 191, 239, 383, 431, 479, 719, 863, 911, 1103, 1151, 1439, 1487, 1583, 1823, 1871};// list of Q to test
    ulong n = qlist[sizeof(qlist)/sizeof(qlist[0])-1]*qlist[sizeof(qlist)/sizeof(qlist[0])-1]*sturm;//need power series up to this bound
    flint_printf("Largest value is %wu\n",n);

    //set eta(q)
    nmod_poly_t eta;
    nmod_poly_init2(eta,p,n);
    int k = 1;
    nmod_poly_set_coeff_ui(eta,0,1);
    while(1) {
        int sgn = ((1-2*(k%2))+p)%p;
        if(k*(3*k-1) >= 2*n)
            break;
        nmod_poly_set_coeff_ui(eta,(k*(3*k-1))/2,sgn);
        if(k*(3*k+1) >= 2*n)
            break;
        nmod_poly_set_coeff_ui(eta,(k*(3*k+1))/2,sgn);
        k++;
    }
    //set eta(q^2)
    nmod_poly_t eta2;
    nmod_poly_init2(eta2,p,n);
    k = 1;
    nmod_poly_set_coeff_ui(eta2,0,1);
    while(1) {
        int sgn = ((1-2*(k%2))+p)%p;
        if(k*(3*k-1) >= n)
            break;
        nmod_poly_set_coeff_ui(eta2,k*(3*k-1),sgn);
        if(k*(3*k+1) >= n)
            break;
        nmod_poly_set_coeff_ui(eta2,k*(3*k+1),sgn);
        k++;
    }
    flint_printf("Made it past setting up the two polynomials\n");

    //We've now set up eta(q) and eta2(q) = eta(q^2)
    //Next we need to invert eta and square and then
    //multiply.  Along the way we clear out old results
    //to help free up memory.
    
    nmod_poly_t etainv;
    nmod_poly_init2(etainv,p,n);
    nmod_poly_inv_series(etainv,eta,n);
    nmod_poly_clear(eta);
    flint_printf("Made it past inverting the series\n");
    nmod_poly_t etainvsq;
    nmod_poly_init2(etainvsq,p,n);
    nmod_poly_mullow(etainvsq,etainv,etainv,n);
    nmod_poly_clear(etainv);
    flint_printf("Made it past the squaring\n");
    nmod_poly_t res;
    nmod_poly_init2(res,p,n);
    nmod_poly_mullow(res,etainvsq,eta2,n);
    nmod_poly_clear(eta2);
    nmod_poly_clear(etainvsq);
    flint_printf("Made it past the result phase\n");
    
    
    //so res contains the generating function for the overpartition function.
    //Now we just loop through appropriate n up to the sturm bound and test (4.4)
    
    int i;
    for(i=0;i<sizeof(qlist)/sizeof(qlist[0]);i++) {
        flint_printf("Testing %wu\n",qlist[i]);
        slong Q = qlist[i];
        slong kappa = 71;
        slong Q1 = power(Q,(kappa-3)/2,p);
        slong eps = legendre_symbol(p-1,p);
        slong sgn = 1 - 2*(((kappa-1)/2)%2);
        slong m;
        int is_interesting=1;
        for(m=1;m <= sturm; m++) {
            if(m%p==0 || legendre_symbol(m,p) == eps)
                continue;
            slong fact = legendre_symbol(sgn*m,Q)*Q1;
            //do the test.  We can safely ignore the n/Q^2 part since it
            //is a non-integer
            if((nmod_poly_get_coeff_ui(res,Q*Q*m)+fact*nmod_poly_get_coeff_ui(res,m))%p != 0) {
                is_interesting=0;
                flint_printf("The value Q = %wu fails 4.4 because of %wu\n",Q,m);
                flint_printf("The coefficients are a(n) = %wu and a(Q^2*n) = %wu\n",
                        nmod_poly_get_coeff_ui(res,m),nmod_poly_get_coeff_ui(res,Q*Q*m));
                break;
            }
        }
        if(is_interesting)
            flint_printf("The value Q = %wu is interesting\n",Q); 
    }
    
    return 0;
}

slong power(slong x, ulong y, slong p) 
{ 
    slong res = 1;      // Initialize result 
  
    x = x % p;  // Update x if it is more than or  
                // equal to p 
  
    while (y > 0) 
    { 
        // If y is odd, multiply x with result 
        if (y & 1) 
            res = (res*x) % p; 
  
        // y must be even now 
        y = y>>1; // y = y/2 
        x = (x*x) % p;   
    } 
    return res; 
} 

slong legendre_symbol(slong x, slong y) {

    if(x%y==0)
        return 0;
    slong res = power(x, (y-1)/2, y);
    if((res+y)%y==1)
        return 1;
    return -1;
}
