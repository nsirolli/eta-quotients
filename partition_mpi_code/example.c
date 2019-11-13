/* This code for p(n) is MPI based.  A master thread will run through
 * the n's up to the sturm bound satisfying the hypothesis of Theorem 4.3
 * This master thread then distributes these various values of n to
 * threads which will test (4.4) in parallel. Computing p(n) is done
 * with FLINT.  When run on a hpc this program would be modified to take
 * in the p, Q, and sturm from the command line so that it could be
 * run from a shell script.
*/

#include <stdlib.h>
#include <stdio.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arith.h>
#include <mpi.h>

void master(int p, int QQ, int STURM, int k, int eps_p, int sgn);
void slave(int p, int QQ, int STURM, int k, int eps_p, int sgn);
int powm(int x, unsigned int y, int p); //modular exponentiation

int nextN(int n, int p, int eps_p); //generate next n to test

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int p = 19; //p we're working with.  In paper this is l
    int QQ = 76607; //Q value we're testing
    int STURM = 6548200; //Sturm bound
    int k = p*p-2; //kappa
    int eps_p = powm(-1,(p-1)/2,p); //legendre_symbol (-1 | p)
    int sgn = powm(-1,(k-1)/2,QQ); //sign in front of n in (4.4)
    if (rank == 0)
        master(p, QQ, STURM, k, eps_p,sgn);
    else
        slave(p, QQ, STURM, k, eps_p,sgn);
    MPI_Finalize();
    return 0;
}


void master(int p, int QQ, int STURM, int k, int eps_p, int sgn) {

    int ntasks, rank;
    int work = -1;
    int result;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);	
    /*
     * Seed the slaves.
     */
    for (rank = 1; rank < ntasks; ++rank) {

        work = nextN(work,p,eps_p);
        MPI_Send(&work,		/* message buffer */
                1,		/* one data item */
                MPI_INT,	/* data item is an integer */
                rank,		/* destination process rank */
                1,	/* user chosen message tag */
                MPI_COMM_WORLD);
    }
    /*
     * Receive a result from any slave and dispatch a new work request
     * until work requests have been exhausted.
     */
    work = nextN(work,p,eps_p);

    while (work <= STURM) {
        MPI_Recv(&result,	/* message buffer */
                1,		/* one data item */
                MPI_INT,	/* data item is a double real */
                MPI_ANY_SOURCE,	/* receive from any sender */
                MPI_ANY_TAG,	/* receive any type of message */
                MPI_COMM_WORLD, 
                &status);	/* info about received message */
        MPI_Send(&work, 1, MPI_INT, status.MPI_SOURCE,
                1, MPI_COMM_WORLD);

        work = nextN(work,p,eps_p);
        if(result == -1) {
            printf("The prime %d fails condition 4.4\n",QQ); //exit out of loop.
            work = STURM+1;
        }
    }
    /*
     * Receive results for outstanding work requests.
     */
    for (rank = 1; rank < ntasks; ++rank) {
        MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE,
                MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
    /*
     * Tell all the slaves to exit.
     */
    for (rank = 1; rank < ntasks; ++rank) {
        int g = -1;
        MPI_Send(&g, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
    }
}

void slave(int p, int QQ, int STURM, int k, int eps_p, int sgn) {
    int result;
    int work;
    MPI_Status status;
    int factor = powm(QQ, (k-3)/2, p);//power of Q in (4.4)
    for (;;) {
        MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);
        /*
         * Check the tag of the received message.  Die if another
         * thread found a counterexample or we're done.
         */
        if (work == -1) {
            return;
        }
        slong n = work;
        //Computer the legendre symbol (\pm n | QQ)
        int ls = powm(sgn*n,(QQ-1)/2, QQ);
        if((ls+QQ)%QQ == 1)
            ls=1;
        else
            ls=-1;

        //Compute p((n+1)/24)
        //and p((Q^2*n+1)/24)
        slong Q = QQ;
        fmpz_t P; fmpz_init_set_ui(P,p);
        fmpz_t an; fmpz_init(an);
        fmpz_t aQQn; fmpz_init(aQQn);
        arith_number_of_partitions(an,(n+1)/24);
        arith_number_of_partitions(aQQn,(Q*Q*((slong)n)+1)/24);
        //mod out by p        
        fmpz_t anm; fmpz_init(anm);
        fmpz_mod(anm,an,P);
        fmpz_t aQQnm; fmpz_init(aQQnm);
        fmpz_mod(aQQnm,aQQn,P);
        fmpz_clear(P);
        fmpz_clear(an); fmpz_clear(aQQn);
        int AQQnm = fmpz_get_si(aQQnm);
        int Anm = fmpz_get_si(anm);
        fmpz_clear(anm); fmpz_clear(aQQnm);
        //test (4.4).  We ignore the n/Q^2 part because for us
        //it will always be a non-integer.
        int r = (AQQnm + ls*factor*Anm);
        if(r%p == 0)
            result=0;
        else {
            result=-1;
            flint_printf("4.4 failed because of n = %wd\n",n);
        }
       
        MPI_Send(&result, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

int nextN(int n, int p, int eps_p) {

    int flag=0;
    while(flag == 0) {
        n+=24;
        if((powm(n,(p-1)/2,p) + eps_p)%p==0)
            flag=1;
    }
    return n;

}

int powm(int x, unsigned int y, int p) 
{ 
    int res = 1;      // Initialize result 

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
