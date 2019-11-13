from sage.modular.etaproducts import qexp_eta

class EtaQuotient():
    '''
        INPUT:

        A dictionary {r_delta:delta} satisfying that sum r_delta*delta is zero
        mod 24.

        EXAMPLES:

            sage: ex3 = EtaQuotient({2:1,1:-2})
            sage: ex3
            Eta quotient of level 16 and weight -1/2 given by the dictionary {1: -2, 2: 1}
            sage: [(l,ex3.conditionC(l)) for l in prime_range(3,20)]
            [(3, (True, -1)),
             (5, (True, 1)),
             (7, (True, -1)),
             (11, (True, -1)),
             (13, (True, 1)),
             (17, (True, 1)),
             (19, (True, -1))]
            sage: ex3.interesting_primes(3,200,-1)
            ([47, 191], [])
            sage: ex3.interesting_primes(5,100,1)
            ([79], [])

            sage: ex8 = EtaQuotient({1:-2,2:-3,4:2})
            sage: ex8
            Eta quotient of level 8 and weight -3/2 given by the dictionary {1: -2, 2: -3, 4: 2}
            sage: [(l,ex8.conditionC(l)) for l in prime_range(3,20)]
            [(3, (True, 1)),
             (5, (True, -1)),
             (7, (True, -1)),
             (11, (True, 1)),
             (13, (True, -1)),
             (17, (True, 1)),
             (19, (True, 1))]
            sage: ex8.interesting_primes(3,100,1)
            ([], [23, 47, 71])
            sage: ex8.interesting_primes(5,200,-1)
            ([], [79, 199])

            sage: ex9 = EtaQuotient({1:2,2:-5,4:-4})
            sage: ex9
            Eta quotient of level 16 and weight -7/2 given by the dictionary {1: 2, 2: -5, 4: -4}
            sage: [(l,ex9.conditionC(l)) for l in prime_range(3,20)]
            [(3, (False, 0)),
             (5, (False, 0)),
             (7, (True, -1)),
             (11, (False, 0)),
             (13, (False, 0)),
             (17, (True, 1)),
             (19, (False, 0))]

            sage: ex10 = EtaQuotient({1:-6,2:-1,4:2})
            sage: ex10
            Eta quotient of level 4 and weight -5/2 given by the dictionary {1: -6, 2: -1, 4: 2}
            sage: [(l,ex10.conditionC(l)) for l in prime_range(3,20)]
            [(3, (True, -1)),
             (5, (True, 1)),
             (7, (True, -1)),
             (11, (True, -1)),
             (13, (True, 1)),
             (17, (True, 1)),
             (19, (True, -1))]
            sage: ex10.interesting_primes(3,100,-1)
            ([11, 23, 47, 59, 71, 83], [])
            sage: ex10.interesting_primes(5,100,1)
            ([19, 59, 79], [])
    '''


    def __init__(self,rdict,compute_cusps=True,desc=None,level=None):

        # We first verify that the dictionary satisfies Treneer's (1.4).
        A = sum([rdict[delta]*delta for delta in rdict])
        assert mod(A,24) == 0

        self._N, self._k = rdict_lw(rdict)
        if level: # we can impose a level; otherwise, the level is the minimal one
            assert mod(level,self._N) == 0
            self._N = level
        self._rdict = rdict
        self._group = Gamma0(self._N)
        self._desc = desc

        # When initializing, by default we compute the cusps for Gamma0(N),
        # together with their widths, since they will be used many times.
        if compute_cusps:
            G = self._group
            self._cusps_dict = {cusp:G.cusp_width(cusp) for cusp in G.cusps()}
        else:
            self._cusps_dict = {}

    def __repr__(self):

        k = self._k
        lstr = 'Eta quotient of level %d and weight ' % self._N
        dstr = ' given by the dictionary ' + str(self._rdict)
        return lstr + k.str() + dstr

    def description(self):
        return self._desc

    def level(self):
        return self._N

    def weight(self):
        return self._k

    def dictionary(self):
        return self._rdict

    def cusps_dict(self):
        if self._cusps_dict == {}:
            G = self._group
            return {cusp:G.cusp_width(cusp) for cusp in G.cusps()}
        else:
            return self._cusps_dict

    def character(self):
        rdict = self._rdict
        s = 2 * prod([delta**rdict[delta] for delta in rdict])
        def chi(d):
            return kronecker_symbol(s,d)
        return chi

    def r(self,cusp):
        '''
            Returs the parameter r_s appearing in the Fourier development of
            self at cusp. See Koblitz, p. 182.
            Computed following Treneer's thesis.
        '''
        chi = self.character()
        return cusp_r(cusp,self._N,chi,2*self._k)


    def vanishing_order(self,cusp):
        '''
            Computes the vanishing order at cusp of self. See Ono, Thm 1.65.
        '''
        N = self._N
        c = cusp_level(cusp,N)
        # Beware: Ligozat's formula works for "normalized" cusps.
        # So we cant use simply c = cusp.denominator(), as stated in Treneer's
        # (2.7).
        rdict = self._rdict
        summands = [gcd(c,delta)**2/delta*rdict[delta] for delta in rdict]
        # gcd(c**2,N) = c*gcd(c,N/c) when c divides N. Hence, Treneer's and
        # Ono's formulas agree.
        return N/24/gcd(c**2,N)*sum(summands)

    def is_holomorphic(self):

        cusps_dict = self.cusps_dict()
        for cusp in cusps_dict:
            if self.vanishing_order(cusp) < 0:
                return False
        return True

    def q_expansion(self,cusp,prec=10,coeff_ring=QQ):
        '''
            Letting w denote the width of the cusp, returns a Laurent series
            qexp in q_w, an integer m and a square root rt such that the
            Fourier development of the eta quotient given by rdict at the cusp
            cusp, up to precision prec, is given by qexp * z24**m * rt.
            When cusp == oo, if a coeff_ring is given, we compute the development
            in that ring (to be more efficient in terms of memory, ideally).
        '''
        if cusp == oo:
        # in this case, we use the "native" q-expansion, which is faster
        # than the algorithm which works over arbitrary cusps.
            ring = coeff_ring[['q']]
            q = ring.0
            eta_q = 1
            rdict = self._rdict
            power = sum([delta * rdict[delta] for delta in rdict])
            v = power//24
            if prec <= v: 
                return O(q**prec) * eta_q , 0, 1
            Nprec = prec - v
            for delta in rdict:
                r = rdict[delta]
                ff = qexp_eta_scaled(ring, Nprec, delta) 
                eta_q = eta_q*ff**r 
            return q**v * eta_q, 0, 1 
        else:
            N = self._N
            G = self._group
            cusp = Cusp(cusp)
            w = G.cusp_width(cusp)
            v = self.vanishing_order(cusp)
            Nprec = max([N/w*(prec-1-v) + 1,0])
            KK = CyclotomicField(N) 
            ring = KK[['q_'+str(w)]]
            qw = ring.0
            res = KK(0)
            eta, qN_pow, z24_pow, rt = eta_quotient_at_cusp(self._rdict,cusp,Nprec)
            for t in range(Nprec):
                if eta[t] != 0:
                    res += eta[t]*qw**(w/N*(t+qN_pow))
            return res + O(qw**prec), z24_pow,rt

    def conditionC_data(self):
        '''
            Returns a dictionary whose keys are those cusps where the eta
            quotient is not holomorphic, and whose value at such a cusp
            is a list whose first element is the width of the cusp, and whose
            second element is the support (=list of indexes) of the non 
            holomorphic part of the Fourier development at the cusp.
        '''
        res = {}
        cdict = self.cusps_dict()
        for cusp in cdict:
            v = self.vanishing_order(cusp)
            if v < 0:
                res[cusp] = [cdict[cusp],[]]
                sing = self.q_expansion(cusp,prec=0)[0]
                for n in range(v,0):
                    if sing[n] != 0:
                        res[cusp][1].append(n)
        return res

    def conditionC(self,l):
        '''
            Returns True, eps if condition C holds at l, with eps = eps.
            When it holds trivially (i.e., when self.conditionC_data() is
            empty), returns True, 0.
            Otherwise, returns False, 0.
        '''    
        data = self.conditionC_data()
        eps = 0
        for cusp in data:
            w = data[cusp][0]
            rs = self.r(cusp)
            for n in data[cusp][1]:
                m = 4*n + rs
                if mod(m,l) != 0:
                    kro = kronecker_symbol(m*w,l)
                    if eps == 0:
                        eps = kro
                    else:
                        if eps != kro:
                            return False,0
        return True,eps

    def nontrivialcC(self,lmax):
        '''
            Decides whether there exists a prime l < lmax such that self
            satisfies condicionC at l non trivially, and in that case returns
            the value of l
        '''
        if not self.is_holomorphic():
            for l in prime_range(3,lmax):
                a,b = self.conditionC(l)
                if a:
                    if b !=0:
                        return True,l
        return False,0


    def beta(self,l,j=1,conditionC=False):
        '''
            Computes the minimum value of beta such that g_{l,j} = 1/2 * f_l *
            F_l**(l**beta) is cuspidal. See Treneer, proof of Thm 3.1.
        '''
        if conditionC:
            assert self.conditionC(l)[0], 'Condition C is not satisfied at %d' %l
        beta = j-1
        cdict = self.cusps_dict()
        v0 = min([self.vanishing_order(cusp) for cusp in cdict])
        if v0 >= 0:
            return beta # Beacause in this case f_l is already cuspidal.
        N = self._N
        Gl = Gamma0(N*l**2)
        for cusp in Gl.cusps():
            c = cusp_level(cusp,N*l**2)
            if mod(c,l**2) != 0: # Because otherwise by Treneer's 3.4, f_l is already cuspidal.
                hc = Gl.cusp_width(cusp)
                Fc = orderF(l,N,cusp)
                ll = integer_floor(log(-v0/Fc,l)) + 1
                beta = max([beta,ll])
        return beta

    def kappa(self,l,j=1,conditionC=False):
        k = self._k
        beta = self.beta(l,j,conditionC)
        if l == 3:
            kl = 24
        else: # if l >= 5
            kl = l**2 - 1
        # So F_l has weight kl/2. There is a tiny mistake in Treneer: she doesn't
        # consider separately the case l = 3.
        kappa = 2*k + l**beta*kl # kappa/2 is the weight of g_{l,j}

        # we add 1 to beta until kappa is positive
        while not kappa > 0:
            kappa += l*kl
        return kappa

    def sturm_bound(self,l,j=1,conditionC=False):
        '''
            Returns the Sturm bound for cusp forms of quadratic character, level
            M and weight kappa/2, where M and kappa/2 are the level and weight
            of the form g_{l,j} given in Treneer's Thm 3.1.
            This bound equals 1/Q**e times the number of coefficients at oo of
            self that we need to compute if we want to prove that Q satisfies
            Treneer's Thms 1.1/1.2.
            Here e = 2 in the half-integral weight case, and e = 1 in the
            integral weight case.
        '''
        k = self._k
        N = self._N
        if k in ZZ:
            N = lcm(N,4)
            # this is to be able to see self as a form of level 2k/2 and level
            # divisible by 4, which is needed in Treneer's 3.1.
        Gl = Gamma0(N*l**2)
        m = Gl.index()
        kappa = self.kappa(l,j,conditionC)

        # We return the Sturm bound given by our Prop 4.1
        return integer_floor((kappa*m)/24 - (m-1)/(8*N*l**2)) + 1

    def coeffs_needed(self,l,L,j=1):
        '''
            Returns the number of coefficients at oo of self that we need to
            compute if we want to check whether the primes Q in L are interesting.
        ''' 
        if L == []:
            return 0
        else:
            k = self._k
            e = k.denominator()
            Qmax = max(L)
            sturm = self.sturm_bound(l,j)
            return Qmax**e * sturm

    def satisfies_congruences(self,l,L,j=1,eps=None):
        '''
            Gives the lists of primes Q in L which are interesting and
            potentially uninteresting, according to Theorems 1.2 and 1.3.
            If we know that condition C is satisfied at l with eps = eps, we can
            pass the value eps.
            We assume that the Q's are good candidates (i.e., primes, equal to
            -1 mod N*l**j).
        '''
        if not eps:
            conditionC, eps = self.conditionC(l)
            # we compute this because we need eps!
            assert conditionC, 'Condition C is not satisfied at %d' %l

        # We need to handle this case separately because of max(L) below.
        if L == []:
            return [],[]

        k = self._k
        e = k.denominator()

        Qmax = max(L) 
        # We compute the q-expansion for the bigger Q in L, and use it for
        # every other candidate.
        sturm = self.sturm_bound(l,j)
        nmax = sturm * Qmax**e
        Zl = IntegerModRing(l**j)
        qexp = self.q_expansion(Cusp(oo),prec=nmax+1,coeff_ring=Zl)[0]

        def b(n):
            # this function gives the n-th coefficient of g_{l,j}, mod l**j.
            if not n in ZZ:
                return 0
            if mod(n,l) == 0:
                return 0
            if kronecker_symbol(n,l) == eps:
                return 0
            else:
                return qexp[n]

        kappa = self.kappa(l,j)

        # We now define the function which gives the n-th coefficient of
        # Treneer's (3.17) / (3.21).
        def is_interesting(Q):
            Ql = Zl(Q)
            if k in ZZ:
                Q1 = Ql**(kappa-1)
                def c(n):
                    return b(Q*n) + Q1*b(n/Q)
            else:
                Q2 = Ql**(kappa-2)
                Q3 = Ql**((kappa-3)/2)
                def c(n):
                    return b(Q**2*n) + kronecker_symbol((-1)**((kappa-1)/2)*n,Q)*Q3*b(n) + Q2*b(n/(Q**2))
            for n in range(1,sturm+1):
                if c(n) != 0:
                    return False
            return True
        inter = []
        uninter = []
        for Q in L:
            if is_interesting(Q):
                inter.append(Q)
            else:
                uninter.append(Q)
        return inter,uninter


    def interesting_primes(self,l,Qmax,j=1,eps=None):
        '''
            Returns the lists of interesting and potentially uninteresting primes
            up to Qmax.
        '''
        N = self._N
        L = candidatesQ(N,l,Qmax)
        return self.satisfies_congruences(l,L)

    def check(self,l,L,nmax,j=1,eps=None,Zl=True):
        '''
            Checks if the congruences of Theorems 1.2 and 1.3 are satisfied for
            the candidate primes Q of the list L, up to nmax.
            Returns the list of pairs (Q,n), where Q is a failing prime and n is
            the first failing index.
            By default, we compute the qexpansion of self mod l**j.
        ''' 
        if Zl:
            ring = IntegerModRing(l**j)
        else:
            ring = ZZ

        if not eps:
            conditionC, eps = self.conditionC(l)
            # we compute this because we need eps!
            assert conditionC, 'Condition C is not satisfied at %d' %l

        k = self._k
        e = 2*k.denominator()-1
        qmax = max(L)
        qprec = (nmax-1)*qmax**e + 1
        qexp = self.q_expansion(Cusp(oo),prec=qprec,coeff_ring=ring)[0]

        res = []
        for Q in L:
            for n in range(nmax):
                if mod(n,Q) !=0:
                    if mod(n,l) !=0:
                        if kronecker_symbol(-n,l) != eps:
                            if mod(qexp[n*Q**e],l**j) != 0:
                                res.append((Q,n))
                                break
        return res

    def uninteresting_primes(self,mmax,lmax=None):
        '''
            Returns the list of the uninteresting primes given by the
            q-expansion of self computed with precision mmax.
            A lmax can be given to save memory.
            '''
        res = []
        if not lmax:
            ring = ZZ
        else:
            M = prod(prime_range(3,lmax))
            ring = IntegerModRing(M)
        qexp = self.q_expansion(Cusp(oo),prec=mmax,coeff_ring=ring)[0]
        N = self._N
        k = self._k
        e = 2*k.denominator()-1
        for m in srange(1,mmax):
            for Q in m.prime_divisors():
                if valuation(m,Q) == e:
                    if mod(Q+1,N) == 0:
                        a = ZZ((Q+1)/N)
                        n = m/Q**e
                        for l in a.prime_divisors():
                            if l > 2:
                                if lmax:
                                    if l >= lmax:
                                        continue
                                if mod(n,l) != 0:
                                    ccl = self.conditionC(l)
                                    if ccl[0]:
                                        eps = ccl[1]
                                        if kronecker_symbol(-n,l) != eps:
                                            if mod(qexp[m],l) != 0:
                                                res.append((m,qexp[m],l,Q,eps))
        return sorted(res, key=lambda x: x[2])


### AUXILIARY FUNCTIONS


def cusp_r(cusp,N,chi,k):
    '''
        Computes the r = 0,1,2,3 appearing in the Fourier expansion of weakly
        holomorphic modular forms of level N, character chi and weight k/2.
        We follow Treneer's thesis, Proposition 2.8.
    '''
    a = cusp.numerator()
    c = cusp.denominator()
    h = N / gcd(c^2,N)
    d,b,_,_ = lift_to_sl2z(c,a,1)
    alpha = Matrix(ZZ,2,[a,b,c,d])
    gamma = alpha * Matrix(ZZ,2,[1,h,0,1]) / alpha
    A, B, C, D = gamma.list()
    z = CC(I)
    w = (d*z-b)/(-c*z+a)
    re, im = sqrt(c*w+d) * sqrt(C*z+D) / sqrt(c*(w+h)+d)
    j = round(re)+I*round(im)
    if mod(D,4) == 1:
        epsD = 1
    else: # if mod(d,4) == 3
        epsD = I
    t = kronecker_symbol(C,D) / epsD * j
    e = t**k * chi(1+a*c*h)
    if e == 1:
        return 0
    if e == I:
        return 1
    if e == -1:
        return 2
    if e == -I:
        return 3

def cusp_level(cusp,N):
    '''
        Computes the level of the cusp with respect to Gamma0(N), i.e., the
        unique divisor c of N such that the cusp is equivalent to a/c
        modulo Gamma0(N), with (a,c) = 1.
    '''
    # Such c exists: stated in Ligozat, 3.2.2. It's given by that formula: Cremona,
    # 2.2.3 (there should be a more straightforward proof). Let cusp = p_2/q_2:
    # by the CRT there exists a unit u mod N such that q_2 \equiv u*c mod N. Let
    # a = u*p_2. Then the criterion is satisfied.
    return gcd(cusp.denominator(),N)

def rdict_lw(rdict):
    '''
        Given a dictionary {delta:r_delta}, returns the level (i.e., the minimal
        N divisible by every delta such that Treneer's (4.1) holds) and the
        weight of the corresponding eta quotient.
    '''
    N = lcm(rdict.keys())

    B = sum([rdict[delta]/delta for delta in rdict])
    if B != 0: # otherwise, we're done!
        C = lcm(N*B,24)
        N = ZZ(C/B).abs()
    k = sum(rdict.values())/2
    # in the half integral weight case, the level must be divisible by 4.
    if not k in ZZ:
        N = lcm(4,N)
    return N,k

def orderF(l,N,cusp):
    '''
        Computes the order of vanishing of Treneer's eta-quotient F_l
        at the cusp cusp of Gamma0(Nl^2). 
        We assume that the level of cusp is not divisible by l**2.
    '''
    c = cusp_level(cusp,N*l**2)
    e = valuation(c,l)
    assert e < 2
    a = N/gcd(c**2,N)
    if l == 3:
        if e == 0:
            ordl = 10
        else: # if e == 1
            ordl = 1
    else: # if l >= 5
        if e == 0:
            ordl = (l**4-1)/24
        else: # if e == 1
            ordl = (l**2-1)/24
    return a*ordl



def qexp_eta_scaled(ps_ring, prec, d):
    """
        computes qexp_eta(R,prec)(x^d)
    """
    prec = Integer(prec)
    if not prec > 0:
        raise ValueError("prec must be a positive integer")
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    try:
        n = 1
        while True:
            pm = -pm
            v[d*n*(3*n-1)//2] = pm
            v[d*n*(3*n+1)//2] = pm
            n += 1
    except IndexError:
        pass
    return ps_ring(v, prec=prec)

def find_qexp_data(delta, gamma):
    a,b,c,d = gamma[0][0], gamma[0][1], gamma[1][0], gamma[1][1]
    A = gcd(c,delta)
    D = delta/A
    B = 1
    while (A*d-B*c) % delta !=0:  B = B+1
    gg = matrix(ZZ,2,2,[a,b,c,d])
    dd = matrix(ZZ,2,2,[delta,0,0,1])
    mm = matrix(ZZ,2,2,[A,B,0,D])

    gamma_p = dd*gg*mm.inverse()

    tt = sage.modular.arithgroup.arithgroup_perm.sl2z_word_problem(gamma_p)

    aa, bb = 0,0
    for pair in tt:
        if pair[0] == 0:
            aa = aa + pair[1]
        if pair[0] == 1:
            bb = bb + pair[1]
    return A,B,D,aa,bb



def eta_factor_at_cusp(N,D,B,A,delta, prec, qN,zz):
    '''
        uses Euler's pentagonal formula
    '''
    kk = N/delta
    zzB = zz**(A*kk*B) # zeta_D**B = zeta_N**(A*kk)    

    xx = zzB*qN**(A**2*kk)
    res = 1
    for ii in range(1,prec):
        if A**2*kk*ii*(3*i+1)/2 < prec:
            res = res + (-1)**ii * (xx**(ii*(3*ii+1)/2) + xx**(ii*(3*ii-1)/2))
        else:
            break
    return res + O(qN**prec)

def eta_factor_at_cusp_fast(N,D,B,A,delta,prec,qN,zz):
    ps_ring = qN.parent()
    kk = N/delta
    zzB = zz**(A*kk*B) # zeta_D**B = zeta_N**(A*kk)    
    xx = zzB*qN**(A**2*kk)
    if not prec > 0:
        raise ValueError("prec must be a positive integer")
    v = [Integer(0)] * prec
    pm = zzB
    v[0] = pm
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    try:
        n = 1
        while True:
            pm = -pm
            v[A**2*kk*n*(3*n-1)//2] = pm*zzB**(n*(3*n-1)//2)
            v[A**2*kk*n*(3*n+1)//2] = pm*zzB**(n*(3*n+1)//2)
            n += 1
    except IndexError:
        pass
    return ps_ring(v, prec=prec)

def eta_quotient_at_cusp(rdict, cusp, prec):
    '''
      RETURN:
          a tuple (qexp,m,n,rt) where qexp is a monic power series in ZZ[['qN']]
          (Is ZZ ok? Or should it be QQ?)
          computed with precision prec, m, n are integers and rt is a
          square root, such that the eta quotient is equal, up to the mentioned
          precision, to rt * qN**m * z24**n * qexp.

      EXAMPLES:

          sage: rr = {3:8} #eta(3z)**8
          sage: kappa = Gamma0(9).cusps()[1]
          sage: kappa
          1/3
          sage: eta_quotient_at_cusp(Gamma0(9),kappa,rr,100)
          (1 - 8*q_9^27 + 20*q_9^54 + O(q_9^100), 9, -72, 1)

    '''
    N = rdict_lw(rdict)[0]
    _, gamma = cusp.is_gamma0_equiv(oo,1,'matrix') #SL(2,Z) equivalence
    KK = CyclotomicField(N)
    zz = KK.0
    ring = KK[['q_'+str(N)]]
    res = 1
    qN = ring.0
    qN_pow = 0
    z24_pow = 0
    rt = 1
    for delta in rdict:
        A,B,D,aa,bb = find_qexp_data(delta,gamma)
        k = N/delta
        res = res * (eta_factor_at_cusp(N,D,B,A,delta,prec,qN,zz))**rdict[delta]

        qN_pow = qN_pow +rdict[delta]*A**2*k/24
        z24_pow = z24_pow + rdict[delta]*(9*aa + bb + B/D)
        rt = rt / sqrt(D) 
    return (res + O(qN**prec),qN_pow,z24_pow,rt) 


def candidatesQ(N,l,Qmax,j=1):
    '''
        Gives the list of primes Q such that Q equals -1 mod N*l**j,
        up to Qmax.
    '''
    primes = prime_range(Qmax)
    L = [Q for Q in primes if mod(Q,N*l**j) == -1]
    return L

