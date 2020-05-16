from flint import *
ctx.prec = 128

def Q(n, y, g, b, k):
    print("n = %i, y = %f, g = %f, b = %f, k = %f" % (n, y, g, b, k))

    # Step 1: compute coefficients a_(r,j) used to evaluate residues
    A_cache = {}   # stores a_(r,j), 0 <= r < 7, 0 <= j < n
    for r in range(7):
        # print("n = %i, r = %i, g = %f, b = %f, k = %f" % (n, r, g, b, k))
        ctx.cap = n             # series precision
        t = arb_series([0, 1])  # series generator
        s = (2*r+1) + t
        R = s / ((s+b)*(s+k)**2)
        T1 = (1-t/2).gamma()            # first gamma factor, shifted to 1
        ctx.cap += 1
        T2 = (-r-t/2).rising(r+1)       # rising factorial
        T2 = arb_series(list(T2)[1:])   # divide by t
        ctx.cap -= 1
        T3 = ((1+s+2*g)/2).rgamma()     # second gamma function, nonsingular
        A = R * (T1 / T2 * T3)**n
        for j in range(n):
            A_cache[r,j] = A[j]

    # Step 2: compute the bound M
    # The functions R, T defined in (2), for evaluation
    def RR(s):
        return s / ((s+b)*(s+k)**2)
    def TT(s):
        return (((1-s)/2).gamma() * ((1+s+2*g)/2).rgamma())**n
    pi = arb.pi()
    i = acb(0,1)
    h = arb("1e-3")
    N = 10**4
    rn = 0
    for j in range(1, N+1):
        rn += abs(RR(-g + i*(j-1)*h))
    rn *= h
    rn += (abs(-g + i*N*h) / abs(-g + b + i*N*h)) * (1/(k-g)) * (pi/2 - ((N*h)/(k-g)).atan())
    rn *= 2
    M = ((1/(g-b)) * (rn * (k-b)**2 / (2 * pi * b * TT(-b))).log()).exp() / y
    M = int(M.ceil().unique_fmpz())
    print("M = ", M)

    # Step 3: precompute B, the x-independent factor in RHS of (7)
    q = 14
    B = ((arb(1+q)/2).gamma() * 2**q) * ((1+q+2*g)/2).rgamma()
    for j in range(q):
        B /= abs(q - 1 - 2*j)
    B **= n
    B /= (2*(q+k))

    # Step 4: verify (3) for x = my, 1 <= m <= M
    for m in range(1, M +1):
      #if m <= 1000:  for quick testing
        # Compute s_m = sum of residues
        s = 0
        x = m*y
        lx = x.log()
        for r in range(7):
            resg = 0
            lxpowers = [arb(1), lx]
            for j in range(2, n):
                lxpowers.append(lxpowers[-1] * lx)
            fac = 1 / arb.fac_ui(n-1)  # 1 / (n-1-j)!
            for j in range(n):
                # assert fac.overlaps(1 / arb.fac_ui(n-1-j))
                lxpow = lxpowers[n-1-j]
                # assert lxpow.overlaps(lx**(n-1-j))
                resg += lxpow * A_cache[r,j] * fac
                fac *= (n-1-j)
            resg *= x**(2*r)
            s += resg  #  add this residue

        # Compute bound ER_m(x)
        E = x**(q-1) * B

        # Check bound in column 2 of Table 3 from Lemma 3
        if m == 1:
            print("R_K bound", (s - E) * (1+g).gamma()**n / RR(1) / 2**(n-1))

        if m in [1,10,100,1000,10000,100000,M]:
            print(m, "of", M, "   s = ", s.str(10, radius=False), "   E = ", E.str(3, radius=False))

        # Verify positivity
        assert s > E


Q(n=10, y=arb("0.00021"), g=arb("0.46"), b=arb("0.37"), k=arb("3.8"))
Q(n=11, y=arb("0.00005"), g=arb("0.42"), b=arb("0.33"), k=arb("3.46"))
Q(n=12, y=arb("0.00001"), g=arb("0.375"), b=arb("0.285"), k=arb("3.1"))
Q(n=13, y=arb("0.00001"), g=arb("0.545"), b=arb("0.465"), k=arb("3.5"))
Q(n=14, y=arb("0.00001"), g=arb("0.72"), b=arb("0.64"), k=arb("4.1"))
Q(n=15, y=arb("0.00001"), g=arb("0.885"), b=arb("0.805"), k=arb("4.4"))
Q(n=16, y=arb("0.00001"), g=arb("1.055"), b=arb("0.965"), k=arb("5.3"))
Q(n=17, y=arb("0.00001"), g=arb("1.205"), b=arb("1.115"), k=arb("5.2"))
Q(n=18, y=arb("0.00001"), g=arb("1.355"), b=arb("1.255"), k=arb("5.7"))

