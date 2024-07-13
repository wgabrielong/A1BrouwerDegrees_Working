------------------------
-- Arithmetic operations
------------------------

-- Input: An integer or rational number
-- Output: The smallest magnitude integer in the same integer square class as the input

squarefreePart = method()
squarefreePart (ZZ) := (ZZ) => (n) -> (
    if n == 0 then return 0;
    tableOfPrimeFactors := hashTable(factor(abs(n)));
    return sub(n/abs(n),ZZ)*product(apply(keys(tableOfPrimeFactors),p -> p^(tableOfPrimeFactors#p%2)));
    )

squarefreePart (QQ) := (ZZ) => (n) -> (
    squarefreePart(numerator(n)*denominator(n))
    )

-- Input: An integer or rational number n
-- Output: A list of prime factors of n

primeFactors = method()
primeFactors (ZZ) := List => (n) -> (
    if abs(n) == 1 then return {};
    sort keys(hashTable(factor(abs(n))))
    )

primeFactors (QQ) := List => (n) -> (
    if not liftable(n,ZZ) then error "tried to take prime factors of a rational";
    primeFactors(sub(n,ZZ))   
    )

-- Input: An integer or rational number and a prime number p
-- Output: The p-adic valuation of the integer or rational number

padicValuation = method()
padicValuation (ZZ, ZZ) := (ZZ) => (n, p) -> (
    if n == 0 then error "Trying to find prime factorization of 0";
    H := hashTable(factor(abs(n)));
    if H#?p then (
    	return(H#p);
	)
    else (
	return 0;
	);
    )

padicValuation (QQ, ZZ) := (ZZ) => (q, p) -> (
    padicValuation(numerator(q),p) - padicValuation(denominator(q),p)
    )

-- Input: An element of a finite field
-- Output: Boolean that gives whether the element is a square

legendreBoolean = method()
legendreBoolean (RingElement) := (Boolean) => a -> (
    if not instance(ring(a),GaloisField) then error "legendreBoolean only works for Galois fields";
    q := (ring a).order;
    -- Detects if a is a square in F_q
    a^((q-1)//2) == 1 
    )

-- Input: An integer a and a prime p
-- Output: 1 if a is a unit square, -1 if a = p^(even power) x (non-square unit), 0 otherwise
-- Note: The terminology "Square Symbol" comes from John Voight's Quaternion Algebra book

squareSymbol = method()
squareSymbol(ZZ, ZZ) := (ZZ) => (a, p) -> (
    x := getSymbol "x";
    R := GF(p, Variable => x);
    e1 := padicValuation(a,p);
    if even e1 then (
    	a1 := sub(a/(p^e1), ZZ);
	a2 := sub(a1, R);
	if legendreBoolean a2 then (
	    return 1;
	    ) 
	else (
	    return -1;
	    );
	);
    return 0;
    )

------------------------------
-- P-adic methods
------------------------------

-- Input: Two integers a and b, and a prime number p
-- Output: Boolean that gives whether a and b differ by a square in Q_p

equalUptoPadicSquare = method()
equalUptoPadicSquare (ZZ, ZZ, ZZ) := (Boolean) => (a, b, p) -> (
    
-- One has to separately handle the cases when p is odd and when p = 2

    if odd p then (
        -- p is odd and we need to check that the powers of p have the same parity, and the units
        -- differ by a square in GF(p)
        a1 := squarefreePart a;
        b1 := squarefreePart b;
        if (padicValuation(a1, p) != padicValuation(b1, p)) then (
	    return false;
            )
        else (
    	    -- c1 will be an integer prime to p
	    c1 := squarefreePart(a1*b1);
	    x := getSymbol "x";
	    return (legendreBoolean(sub(c1, GF(p, Variable => x)))); 
	    );
        )
    else (
        -- Case when p=2. Here we have to check that the powers of p have the same parity, and 
        -- that the units agree mod 8.
        a1 = squarefreePart a;
        b1 = squarefreePart b;
        if (padicValuation(a1, p) != padicValuation(b1, p)) then (
	    return false;
            )
        else (
    	    -- c1 will be an integer prime to p
	    c1 = squarefreePart(a1*b1);
	    c1 = c1 % 8;
	    -- If c1 = 1, then the two odd units are congruent mod 8, and are squares in Q_2
	    return (c1 == 1); 
	    );
        );
    )

-- Input: An integer a and a prime number p
-- Output: Boolean that gives whether a is a square in QQ_p

isPadicSquare = method()
isPadicSquare (ZZ, ZZ) := (Boolean) => (a, p) -> (
    equalUptoPadicSquare(a,1,p)
    )

------------------------------
-- Commutative algebra methods
------------------------------

-- Input: A list L of functions f1,...,fn over the same ring R and p is a prime ideal of an isolated zero
-- Output: A list of basis elements of the local k-algebra Q_p(f) = R[x]_p/(f)

localAlgebraBasis = method()
localAlgebraBasis (List, Ideal) := (List) => (L,p) -> (
    
    -- Verify that the ideal p is prime
    if not isPrime(p) then error "ideal is not prime";
    
    -- Ambient ring
    R := ring L#0;
    I := ideal L;
    
    -- Check whether or not the ideal is zero-dimensional
    if dim I > 0 then error "morphism does not have isolated zeroes";
    if not isSubset(I,p) then error "prime is not a zero of function";
    J := I:saturate(I,p);
    A := R/J;
    B := basis A;
    flatten entries(B)
    )

-- Input: A zero-dimensional ideal (f_1,...,f_n) < k[x_1,...,x_n].
-- Output: The rank of the global algebra as a k-vector space.

rankGlobalAlgebra = method()
rankGlobalAlgebra (List) := (ZZ) => (Endo) -> (
    
    -- Get the underlying field    
    kk := coefficientRing(ring(Endo#0));    
    if not isField(kk) then (
    	kk = toField kk;
    	);
    
    -- Let S = k[x_1,...,x_n] be the ambient polynomial ring
    S := ring(Endo#0);
    
    -- First check if the morphism does not have isolated zeroes
    if dim ideal(Endo) > 0  then error "ideal is not zero-dimensional";
    
    -- Get the rank of S/ideal(Endo) as a kk-vector space
    numColumns(basis(S/ideal(Endo)))
    )

