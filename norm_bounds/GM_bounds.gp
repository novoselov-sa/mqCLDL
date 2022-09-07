\\
\\ This set of PARI/gp functions computes bounds for a set of generators of the
\\ class  group of a number field K using the algorithms described in
\\ Grenié-Molteni "Explicit bounds for generators of the class group".
\\
\\ The fundamental instruction is buch_limc_gp(nf):
\\ INPUT: a polinomial or the result of the corresponding nfinit
\\ OUTPUT: a 3x3 matrix M, with the following meaning:
\\ M[,1] = for each row, a bound T,
\\ M[,2] = for each row, the number of ideals below T,
\\ M[,3] = for each row, the time required for the computation.
\\ M[1,] = results for T=T(K)   i.e. Belabas-Diaz y Diaz-Friedman's algorithm
\\ M[2,] = results for T=T_1(K) i.e. Grenié-Molteni's               algorithm
\\ M[3,] = results for T=T_2(K) i.e. Grenié-Molteni's SIMPLIFIED    algorithm
\\
\\ Moreover, one can get the bound for each algorithm separately:
\\ Belabas et al's algorithm using BDyDF(nf),
\\ Grenié-Molteni's algorithm using GRHoptimize(nf), 
\\ Grenié-Molteni's simplified algorithm using newalgo(nf).
\\

if (type(qfeval) != "t_CLOSURE", eval("qfeval(A,v) = v~*A*v;");)

GRHoptimsum(nf, logC)=
{
    my(SA, SB);
    cache_prime_dec(exp(logC + 1e-9), nf);
    for(n = 1, #Sprimes,
	my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
	my(lim = logC\logp);
	if (!lim, break);
	for (j = 1, #fs,
	    my(f = fs[j], M, nb);
	    my(logNP, q, A, B);
	    if (f > lim, next);
	    logNP = f * logp;
	    q = p^(-f/2);
	    A = logNP * q; B = logNP * A; M = lim\f;
	    if (M > 1,
		my(inv1_q = (1 - q)^-1, powqM = q^M);
		A *= (1 - powqM) * inv1_q;
		B *= (1 - powqM * (M+1 - M * q)) * inv1_q^2;
	    );
	    nb = ns[j];
	    SA += nb * A;
	    SB += nb * B;
	);
    );
    [ScD - 2*SA, 2*SB];
}

GRHoptimdifference(nf, N, R1, logC)=
{
    my(cD, cN, cmhalf);
    if (!logC, return(0));
    [cD, cN] = GRHoptimsum(nf, logC);
    cN += ScN;
/*
 * cosh integral:
 * 4*Catalan-4*imag(dilog(I/C^(1/2)))
 *
 * sinh integral:
 * Pi^2/2-4*dilog(C^(-1/2))+dilog(1/C)
 */
    cmhalf = exp(-logC/2);
    if (R1,
	cN -= 4*R1 * imag(dilog(I * cmhalf));
    );
    cN += N * (dilog(cmhalf^2) - dilog(cmhalf)*4);
    cD*logC + cN;
}

cache_prime_dec(LIMC, nf)=
{
    my(limp, SC, i = 0, P = nf.pol,index = nf.index);
    if (Spol != P,
	Sprimes = [];
	Spol = P;
	ScN = 4*nf.r1*Catalan + #nf.zk*Pi^2/2;
	ScD = log(abs(nf.disc)) - #nf.zk*(Euler+log(8*Pi)) - nf.r1*Pi/2;
    );
    limp = if (#Sprimes, Sprimes[#Sprimes][1], 1);
    if (limp >= LIMC, return());
    SC = vector(primepi(LIMC) - primepi(limp));
    forprime(p = limp + 1, LIMC,
	my(Ps, ns, fs, k = 1);
	i++;
	if (index % p,
	    my(fac = factormod(P, p, 1));
	    fs = fac[,1]~;
	,
	    fs = vecsort(apply(P->P.f, idealprimedec(nf, p)));
	);
	ns = vector(#fs, j, 1);
	for(j = 2, #fs,
	    if (fs[j] != fs[k] ,
		k++;
		fs[k] = fs[j];
	    ,
		ns[k]++;
	    );
	);
	fs = fs[1..k];
	ns = ns[1..k];
	SC[i] = [p, log(p), fs, ns];
    );
    Sprimes = concat(Sprimes, SC);
}

listnorms(nf, C)=
{
  my(N = 2, logC = log(C), norms);
  for(n = 1, #Sprimes,
      my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
      my(lim = logC\logp);
      if (!lim, break);
      for (j = 1, #fs, if (fs[j] <= lim, N++));
  );
  norms = vector(N);
  norms[1] = 1;
  norms[N] = ceil(C);
  N = 2;
  for (n = 1, #Sprimes,
      my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
      my(lim = logC\logp);
      if (!lim, break);
      for (j = 1, #fs,
	  my(f = fs[j]);
	  if (f <= lim,
	      norms[N] = p^f;
	      N++;
	  );
      );
  );
  vecsort(norms);
}

GRHsum(logC) =
{
    my(SA=0, SB=0);
    for(n = 1, #Sprimes,
	my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
	my(lim);
	if (logp > logC, break);
	lim = logC\logp;
	for(j = 1, #ns,
	    my(logNP, q, A, B, M, f = fs[j]);
	    if (f > lim, next);
	    logNP = f * logp;
	    q = 1/p^(0.5*f);
	    A = logNP * q; B = logNP * A; M = lim\f;
	    if (M > 1,
		my(inv1_q = 1 / (1-q), qM = q^M);
		A *= (1 - qM) * inv1_q;
		B *= (1 - qM*(M+1 - M*q)) * inv1_q * inv1_q;
	    );
	    nb = ns[j];
	    SA += nb * A;
	    SB += nb * B;
	);
	if (logp == logC, break);
    );
    [ ScD - 2*SA, 2*SB ];
}

GRHdifference(nf, logC)=
{
    my(cD, cN);
    [cD, cN] = GRHsum(logC);
    cD + (ScN + cN)/logC;
}

GRHchk(nf, LIMC)=
{
    cache_prime_dec(LIMC, nf);
    GRHdifference(nf, log(LIMC)) < -1e-8;
}

BDyDF(nf)=
{
    my(high = 2, low = 1, test);
    if (type(nf) == "t_POL", nf = nfinit(nf));
    while(!GRHchk(nf, high),
	low = high;
	high *= 2;
    );
    while (high > low+1,
	test = (high+low)/2;
	if (GRHchk(nf, test),
	    high = test;
	,
	    low = test;
	);
    );
    high;
}

newGRHchk(nf, LIMC, N, R1, LOGD)=
{
    my(dLIMC = LIMC, sLIMC = sqrt(LIMC), isLIMC = 1/sLIMC, logC = log(LIMC));
    my(D);
    D = (sLIMC-1)*(LOGD-(Euler+log(2*Pi*(1-isLIMC)))*N-log(2)*R1)
      - (logC/8-0.5)*logC*(N+R1);
    cache_prime_dec(LIMC, nf);
    for (n = 1, #Sprimes,
	my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
	my(lim);
	if (logp > logC, break);
	lim = logC\logp;
	for(j = 1, #ns,
	    my(logNP, q, A, M, f = fs[j]);
	    if (f > lim, next);
	    logNP = f * logp;
	    M = lim\f;
	    /* sum of x and L-x, partial compensation */
	    A = ((M+1)\2)*(logC-((M+1)\2)*logNP);
	    if (M == 1, /* nothing below L/2 */
	    , M <= 3, /* only one term below L/2 */
		A += -2 + 2*sLIMC/p^f;
	    , /* at least two terms below L/2, this is infrequent */
		my(q = p^(-1.*f));
		A += -2*(M\2) + 2*sLIMC*q*(1-q^(M\2))/(1-q);
	    );
	    D -= ns[j] * A * logNP;
	    if (D < -1e-8, return(1));
	);
	if (logp == logC, break);
    );
    0;
}

/* Give an upper bound for GRHoptimze */
GRHoptimizebound(nf)=
{
    my(TLIM, N, R1, LOGD, low, high);
    /* We have several competing "easy" bounds */
    /* The first one is the "theoretical bound" which is usually worse than
     * the other two, but is here as a safeguard.
     */
    N    = #nf.zk;
    if (N == 1, return(1));
    R1   = nf.r1;
    LOGD = log(abs(nf.disc));
    TLIM = LOGD + log(LOGD) - (Euler+log(2*Pi))*N + 1 + (N+1)*log(7*LOGD)/LOGD;
    TLIM = 4*TLIM^2;
    if (N == 2, TLIM = min(TLIM, 4.01*LOGD^2));
    TLIM = floor(TLIM);
    low = high = 2;
    /* The second and third are T(K) and T_2(K), used together so that
     * the best one will fix the limit */
    while (!GRHchk(nf, high) && !newGRHchk(nf, high, N, R1, LOGD) &&
	   high < TLIM,
	low = high;
	high *= 2;
    );
    if (high > TLIM,
	if (!GRHchk(nf, TLIM) && !newGRHchk(nf, TLIM, N, R1, LOGD),
	    return(TLIM);
	,
	    high = TLIM;
	);
    );
    while (high - low > 1,
	my(test = (high + low)\2);
	if (!GRHchk(nf, test) && !newGRHchk(nf, test, N, R1, LOGD),
	    low = test;
	,
	    high = test;
	);
    );
    high;
}

GRHoptimize(nf)=
{
    my(N, R1, cp, norms, lognorms, LOGD, ohigh, ihigh, LIMD, NDIV);
    my(found);
    if (type(nf) == "t_POL", nf = nfinit(nf));
    N     = #nf.zk;
    if (N == 1, return(1));
    R1    = nf.r1;
    LOGD  = log(abs(nf.disc));
    norms = listnorms(nf, GRHoptimizebound(nf));
    /* Approximate log of next norm *from below* (except last one) */
    lognorms = vector(#norms, i,
			if (i<#norms, (1-1e-9)*log(norms[i+1]), log(norms[i])));
    ohigh = 0;
    ihigh = #norms;
    LIMD  = if (LOGD > 1000, 2*LOGD - 0.000001, 0);
    NDIV  = 1;
    found = 0;
    while (!found || ohigh != ihigh,
	/* This loops on NDIV */
	my(ilow, ieps);
	if (found, ohigh = ihigh);
	NDIV *= 2;
	ilow  = 0;
	ieps  = 1;
	while(ihigh > ilow+1,
	    /* This loops on norms */
	    my(itest, delta, tab, mat, matLDLT);
	    if (!found
	    ,
		/* We are looking for an NDIV for which the upper bound is
		 * satisfied. */
		itest = ihigh;
	    ,
		ilow
	    ,
		itest = (ihigh + ilow)\2;
	    ,
		itest = max(ihigh - ieps, 1);
		ieps *= 2;
	    );
	    delta   = lognorms[itest]/2/NDIV;
	    tab     = vector(2*NDIV);
	    mat     = matid(NDIV);
	    matLDLT = matid(NDIV);
	    for(i = 1, NDIV,
		my(logC = 2*i*delta, elt);
		tab[2*i-1] = GRHoptimdifference(nf, N, R1, logC-delta);
		tab[2*i  ] = GRHoptimdifference(nf, N, R1, logC);
		for(j = 1, i-1,
		    mat[i,j] = mat[j,i] = tab[i+j] - tab[i-j];
		);
		mat[i,i] = tab[2*i];
		for (j = 1, i-1,
		    elt = mat[i,j] - sum(k=1, j-1, matLDLT[i,k]*matLDLT[k,j]);
		    matLDLT[j,i] = elt;
		    matLDLT[i,j] = elt*matLDLT[j,j];
		);
		elt = mat[i,i] - sum(k=1, i-1, matLDLT[i,k]*matLDLT[k,i]);
		if (elt <= 0,
		    /* The couple (NDIV, delta) seems to work, however we want
		     * to check that there were not too much rounding error. */
		    my(col = vectorv(i, j, if (j<i,matLDLT[i,j],-1)));
		    forstep(j = i-1, 1, -1,
			my(piv = col[j]);
			for(k = 1, j-1,
			    col[k] -= piv*matLDLT[j,k];
			);
		    );
		    /* col should now be a vector on which the quadratic form
		     * is negative. */
		    if (qfeval(mat[1..i,1..i], col)<0,
			ihigh = itest;
			found = 1;
		    ,
			/* Below: if !found, we are still looking for the right
			 * NDIV. */
			if (!found, break);
			ilow = itest;
		    );
		    break;
		);
		matLDLT[i,i] = elt^-1;
	    );
	    if (!found, break);
	    if (ihigh != itest,
		/* Failed */
		ilow = itest;
	    );
	);
    );
    norms[ihigh];
}

newalgo(nf)=
{
    my(high = 2, low = 1, test, N, R1, LOGD);
    if (type(nf) == "t_POL", nf = nfinit(nf));
    N = #nf.zk;
    R1 = nf.r1;
    LOGD = log(abs(nf.disc));
    while(!newGRHchk(nf, high, N, R1, LOGD),
	low = high;
	high *= 2;
    );
    while (high > low+1,
	test = (high+low)\2;
	if (newGRHchk(nf, test, N, R1, LOGD),
	    high = test;
	,
	    low = test;
	);
    );
    high;
}

bnfcountideals(nf, LIMC)=
{
    my(logC = log(LIMC), count=0, pC=0, lC);
    lC = isprimepower(LIMC, &pC);
    cache_prime_dec(LIMC, nf);
    for(n = 1, #Sprimes,
	my(pr = Sprimes[n], p = pr[1], logp = pr[2], fs = pr[3], ns = pr[4]);
	my(lim);
	if (logp > logC, break);
	lim = if (p == pC, lC,  logC\logp);
	for (j = 1, #fs, if (fs[j] <= lim, count += ns[j]));
    );
    count;
}

buch_limc_gp(nf)=
{
    my(TK, T1K, T2K, tk, t1k, t2k);
    if (type(nf) == "t_POL", nf = nfinit(nf));
    gettime();
    TK = BDyDF(nf);
    tk = gettime();
    T1K = GRHoptimize(nf);
    t1k = gettime();
    T2K = newalgo(nf);
    t2k = gettime();
    [  TK, bnfcountideals(nf,  TK),  tk;
      T1K, bnfcountideals(nf, T1K), t1k;
      T2K, bnfcountideals(nf, T2K), t2k];
}
