attach("eta-quotients.sage")

##############################################################################
#Eta quotient of level 1 and weight -12 given by the dictionary {1: -24}
##############################################################################

ex0 = EtaQuotient({1:-24}, desc='1/Delta, https://oeis.org/A006922')

# [(p,ex0.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, -1)), (5, (True, 1)), (7, (True, -1))]

# ex0.interesting_primes(3,100)
# #([2, 5, 11, 17, 23, 29, 41, 47, 53, 59, 71, 83, 89], [])

# a,b = ex0.interesting_primes(3,10^5)
# len(a),len(b)
# # 4807, 0

# a,b = ex0.interesting_primes(5,10^5)
# len(a),len(b)
# # 2390, 0

# a,b = ex0.interesting_primes(7,10^5)
# len(a),len(b)
# # 1596, 0

# a,b = ex0.interesting_primes(11,10^4)
# len(a),len(b)
# # 123, 0

# a,b = ex0.interesting_primes(13,10^4)
# len(a),len(b)
# # 10, 90


###############################################################################
#Eta quotient of level 576 and weight -1/2 given by the dictionary {24: -1}
###############################################################################

ex1 = EtaQuotient({24:-1}, desc='a(n) is the number of partitions of (n+1)/24')

# [(p,ex1.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (False, 0)), (5, (True, 1)), (7, (True, -1))]

# ex1.interesting_primes(5,1000)
# #([], [])

# ex1.interesting_primes(5,3000)
# # Memory error
# # See partition_sage_code for these examples
# # Q = 2879 is the first interesting prime


# #############################################################################
# #Eta quotient of level 1 and weight -24 given by the dictionary {1: -48}
# #############################################################################

ex2 = EtaQuotient({1:-48}, desc='1/Delta^2')

# [(p,ex2.conditionC(p)) for p in prime_range(3,100)]

# #[(3, (False, 0)), (5, (False, 0)), (7, (True, -1)), (11, (False, 0)), (13, (False, 0)), (17, (True, 1)), (19, (False, 0)), (23, (True, -1)), (29, (False, 0)), (31, (True, -1)), (37, (False, 0)), (41, (True, 1)), (43, (False, 0)), (47, (True, -1)), (53, (False, 0)), (59, (False, 0)), (61, (False, 0)), (67, (False, 0)), (71, (True, -1)), (73, (True, 1)), (79, (True, -1)), (83, (False, 0)), (89, (True, 1)), (97, (True, 1))]

# a,b = ex2.interesting_primes(7,10^5)
# len(a),len(b)
# # 1596, 0

# ex2.interesting_primes(17,1000)
# #([67, 101, 271, 373, 509, 577, 883], [])

# ex2.interesting_primes(23,1000)
# #([], [137, 229, 367, 643, 827, 919])


###################################################################################
#Eta quotient of level 16 and weight -1/2 given by the dictionary {1: -2, 2: 1}
##################################################################################

ex3 = EtaQuotient({2:1,1:-2}, desc='a(n) is the number of overpartitions of n')

# [(p,ex3.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, -1)), (5, (True, 1)), (7, (True, -1))]

# ex3.interesting_primes(3,1000)
# #([47, 191, 239, 383, 431, 479, 719, 863, 911], [])
## NCR Memory error now

# ex3.interesting_primes(5,1000)
# #([79, 239, 479, 719], [])
## NCR Memory error now

# ex3.interesting_primes(7,1000)
# #([223], [])
## NCR Memory error now


# ###########################################################################
# #Eta quotient of level 2 and weight 8 given by the dictionary {1: 8, 2: 8}
# ###########################################################################

ex4 = EtaQuotient({1:8,2:8}, desc='A002288')

# [(p,ex4.conditionC(p)) for p in prime_range(3,20)]
# #[(3, (True, 0)), (5, (True, 0)), (7, (True, 0)), (11, (True, 0)), (13, (True, 0)), (17, (True, 0)), (19, (True, 0))]

# a,b = ex4.interesting_primes(3,10^3)
# len(a),len(b)
# # 86, 0

# ex4.interesting_primes(5,1000)
# #([19, 29, 59, 79, 89, 109, 139, 149, 179, 199, 229, 239, 269, 349, 359, 379, 389, 409, 419, 439, 449, 479, 499, 509, 569, 599, 619, 659, 709, 719, 739, 769, 809, 829, 839, 859, 919, 929], [])

# ex4.interesting_primes(7,1000)
# #([139, 349, 503, 587, 643, 769], [13, 41, 83, 97, 167, 181, 223, 251, 293, 307, 419, 433, 461, 601, 727, 797, 811, 839, 853, 881, 937])

# ex4.interesting_primes(11,1000)
# #([131, 439], [43, 109, 197, 241, 263, 307, 373, 461, 571, 593, 659, 769, 857, 967])

# ex4.interesting_primes(13,1000)
# #([337], [103, 181, 233, 311, 389, 467, 571, 701, 727, 857, 883])

# ex4.interesting_primes(17,1000)
# #([67, 101, 271, 373, 509, 577, 883], [])

# ex4.interesting_primes(19,1000)
# #([379], [37, 113, 151, 227, 569, 607, 683, 797, 911])


############################################################################
#Eta quotient of level 3 and weight 3 given by the dictionary {1: -3, 3: 9}
############################################################################

ex5 = EtaQuotient({1:-3,3:9}, desc='A106402,a(n+1) is the number of partition triples of n where each partition is 3-core')

# [(p,ex5.conditionC(p)) for p in prime_range(3,20)]
# #[(3, (True, 0)), (5, (True, 0)), (7, (True, 0)), (11, (True, 0)), (13, (True, 0)), (17, (True, 0)), (19, (True, 0))]

# ex5.interesting_primes(3,1000)
# #([17, 53, 71, 89, 107, 179, 197, 233, 251, 269, 359, 431, 449, 467, 503, 521, 557, 593, 647, 683, 701, 719, 773, 809, 827, 863, 881, 953, 971], [])

# ex5.interesting_primes(7,1000)
# #([41, 83, 167, 251, 293, 419, 461, 503, 587, 797, 839, 881], [])

# ex5.interesting_primes(11,1000)
# #([131, 197, 263, 461, 593, 659, 857], [])

# ex5.interesting_primes(13,1000)
# #([233, 311, 389, 467, 701, 857], [])

# ex5.interesting_primes(17,1000)
# #([101, 509], [])

# ex5.interesting_primes(19,1000)
# #([113, 227, 569, 683, 797, 911], [])
# NCR memory error


####################################################################################
#Eta quotient of level 4 and weight 1 given by the dictionary {1: -4, 2: 10, 4: -4}
####################################################################################

ex6 = EtaQuotient({1:-4,2:10,4:-4}, desc='A004018, number of points in square lattice on the circle of radius sqrt(n). Equivalently, number of Gaussian integers of norm n')

# [(p,ex6.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, 0)), (5, (True, 0)), (7, (True, 0))]

# ex6.interesting_primes(3,1000)
# #([11, 23, 47, 59, 71, 83, 107, 131, 167, 179, 191, 227, 239, 251, 263, 311, 347, 359, 383, 419, 431, 443, 467, 479, 491, 503, 563, 587, 599, 647, 659, 683, 719, 743, 827, 839, 863, 887, 911, 947, 971, 983], [])

# ex6.interesting_primes(5,1000)
# #([19, 59, 79, 139, 179, 199, 239, 359, 379, 419, 439, 479, 499, 599, 619, 659, 719, 739, 839, 859, 919], [])

# ex6.interesting_primes(7,1000)
# #([83, 139, 167, 223, 251, 307, 419, 503, 587, 643, 727, 811, 839], [])


#################################################################################
#Eta quotient of level 144 and weight -2 given by the dictionary {24: -4}
#################################################################################

ex7 = EtaQuotient({24:-4})

# [(p,ex7.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (False, 0)), (5, (False, 0)), (7, (True, -1))]

# ex7.interesting_primes(7,4000)
# #([3023], [])
# NCR now a memory error


#######################################################################################
#Eta quotient of level 8 and weight -3/2 given by the dictionary {1: -2, 2: -3, 4: 2}
#######################################################################################

ex8 = EtaQuotient({1: -2, 2: -3, 4: 2}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex8.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, 1)), (5, (True, -1)), (7, (True, -1))]

# ex8.interesting_primes(3,1000)

# ([], [23, 47, 71, 167, 191, 239, 263, 311, 359, 383, 431, 479, 503, 599, 647, 719, 743, 839, 863, 887, 911, 983])
# NCR memory error


#############################################################################################
#Eta quotient of level 16 and weight -7/2 given by the dictionary {1: 2, 2: -5, 4: -4}
#############################################################################################

ex9 = EtaQuotient({1: 2, 2: -5, 4: -4}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex9.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (False, 0)), (5, (False, 0)), (7, (True, -1))]

# ex9.interesting_primes(7,1000)
# ([], [223])
# NCR memory error


#############################################################################################
#Eta quotient of level 4 and weight -5/2 given by the dictionary {1: -6, 2: -1, 4: 2}
#############################################################################################

ex10 = EtaQuotient({1:-6,2:-1,4:2}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex10.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, -1)), (5, (True, 1)), (7, (True, -1))]

# ex10.interesting_primes(3,1000)

# #([11, 23, 47, 59, 71, 83, 107, 131, 167, 179, 191, 227, 239, 251, 263, 311, 347, 359, 383, 419, 431, 443, 467, 479, 491, 503, 563, 587, 599, 647, 659, 683, 719, 743, 827, 839, 863, 887, 911, 947, 971, 983], [])

# ex10.interesting_primes(5,1000)
# #([19, 59, 79, 139, 179, 199, 239, 359, 379, 419, 439, 479, 499, 599, 619, 659, 719, 739, 839, 859, 919], [])

# #Memory error at p=7


#############################################################################################
#Eta quotient of level 4 and weight -5/2 given by the dictionary {1: 2, 2: -1, 4: -6}
#############################################################################################

ex11 = EtaQuotient({1:2,2:-1,4:-6}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex11.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, -1)), (5, (True, 1)), (7, (True, -1))]

# ex11.interesting_primes(3,100)
# #([11, 23, 47, 59, 71, 83], [])


#############################################################################################
#Eta quotient of level 8 and weight -3/2 given by the dictionary {1: -2, 2: -3, 4: 2}
#############################################################################################

ex12 = EtaQuotient({1:-2,2:-3,4:2}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex12.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, 1)), (5, (True, -1)), (7, (True, -1))]

# ex12.interesting_primes(3,100)
# #([], [23, 47, 71])


#############################################################################################
#Eta quotient of level 4 and weight -1/2 given by the dictionary {1: 10, 2: -5, 4: -6}
#############################################################################################


ex13 = EtaQuotient({2:-5,4:-6,1:10}, desc='an interesting example (half-integral weight, satisfies condC non trivially) of low level')

# [(p,ex13.conditionC(p)) for p in prime_range(3,10)]
# #[(3, (True, -1)), (5, (True, 1)), (7, (True, -1))]

# ex13.interesting_primes(3,100)
# #([11, 23, 47, 59, 71, 83], [])
