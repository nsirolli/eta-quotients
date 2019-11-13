# This code works on sage and utilizes all of the cores
# on a single machine to test 4.4 for the partition function.
# The MPI code is better for an hpc where the number of cores you
# are allowed to use on a single node is limited.

from sage.combinat.partition import *
from multiprocessing import Pool

p = 13
Qlist = [7487,44927,67391]
sturm = 1458909
k = p^2 - 2
eps_p = legendre_symbol(-1,p)
sgn = pow(-1,(k-1)//2)

for Q in Qlist:
    factor = pow(Q,(k-3)//2,p)
    nlist = [n for n in range(sturm+1) if n % 24 == 23 and legendre_symbol(n,p) == -1*eps_p]
    
    def a(m):
        return number_of_partitions((m+1)//24)

    def test(n):
        return (a(Q^2*n) + legendre_symbol(sgn*n,Q)*factor*a(n)) % p == 0

    P = Pool()
    results = P.map(test, nlist)
    flag = reduce(lambda x,y: x and y, results, True)
    print Q, flag

