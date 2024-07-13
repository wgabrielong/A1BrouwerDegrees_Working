---------------------------------------
-- Anisotropic dimension
---------------------------------------

-- Input: A Grothendieck-Witt class beta over QQ and a prime number
-- Output: Boolean that indicates whether beta is hyperbolic over the p-adic rationals Q_p
-- Notes: Koprowski/Czogala Algorithm 6

isHyperbolicQp = method()
isHyperbolicQp (GrothendieckWittClass, ZZ) := Boolean => (beta, p) -> (
    B := beta.matrix;
    rankFormBeta := rankForm beta;
    kk := ring B;
    
    if not (kk === QQ) then error "GrothendieckWittClass is not over QQ";
    if not isPrime(p) then error "second argument must be a prime number";
    
    -- Odd rank forms are not hyperbolic
    if odd rankFormBeta then return false; 
    
    -- Hyperbolic forms always have square discriminants
    -- Note that Koprowski and Czogala are using a different, signed, version of the discriminant
    d := (-1)^(rankFormBeta*(rankFormBeta-1)/2) *integralDiscriminant(beta);
    
    -- If this discriminant is not a square in Q_p then return false
    if not isPadicSquare(d,p) then return false;
    
    -- At this stage, the rank and discriminant of our beta agrees with that of a hyperbolic form,
    -- so by e.g. Lam V.3.25 it suffices to check whether their Hasse-Witt invariants agree
    m := sub(rankFormBeta/2,ZZ);
    HasseWittHyperbolicForm := (HilbertSymbol(-1,-1,p))^(m*(m - 1)/2);
    HasseWittBeta := HasseWittInvariant(beta,p);
    HasseWittHyperbolicForm == HasseWittBeta
    )

-- Input: A Grothendieck-Witt class beta over QQ and a prime number
-- Output: An integer, the rank of the anisotropic part of beta over the p-adic rationals Q_p
-- Note that as all quadratic forms over Q_p of dimension >= 5 are isotropic, 
-- this method will always output 0, 1, 2, 3, or 4.
-- This is an implementation of Koprowski/Czogala Algorithm 8

anisotropicDimensionQp = method()
anisotropicDimensionQp (GrothendieckWittClass, ZZ) := ZZ => (beta, p) -> (
    B := beta.matrix;
    rankFormBeta := rankForm beta;
    kk := ring B;
    
    if not (kk === QQ) then error "GrothendieckWittClass is not over QQ";
    if not isPrime(p) then error "second argument must be a prime number";
    
    if even rankFormBeta then (
	-- If the form is hyperbolic it has no anisotropic part
	if isHyperbolicQp(beta,p) then return 0;
       	
	-- Note Koprowski and Czogala use a signed version of the discriminant
	d := (-1)^(rankFormBeta*(rankFormBeta-1)/2) * integralDiscriminant(beta);
	if isPadicSquare(d,p) then return 4;
	return 2;
	);
    
    if odd rankFormBeta then (
	c := (-1)^(rankFormBeta*(rankFormBeta+1)/2) * integralDiscriminant(beta);
	gamma := gwAdd(beta, diagonalForm(QQ,(c)));
	if isHyperbolicQp(gamma,p) then return 1;
	return 3;
	);
    )

-- Input: A Grothendieck-Witt class beta over QQ
-- Output: An integer, the rank of the anisotropic part of beta over QQ
-- Notes: Follows Algorithm 9 of Koprowski/Czogala

anisotropicDimensionQQ = method()
anisotropicDimensionQQ (GrothendieckWittClass) := ZZ => (beta) -> (
    B := beta.matrix;
    rankFormBeta := rankForm beta;
    kk := ring B;
    
    if not (kk === QQ) then error "GrothendieckWittClass is not over QQ";
    
    -- The anisotropic dimension of a form over Q is the maximum of its anisotropic dimensions at any of its completions or over Q_2
    ListOfLocalAnistropicDimensions := {};
    
    -- The anisotropic dimension at RR is the absolute value of the signature of the form
    ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, abs(signature(beta)));
    
    -- We always have to add the anisotropic dimension at the prime 2
    ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, anisotropicDimensionQp(beta,2));
       
    -- For the remaining local fields, we can just look at relevant primes
    for p in relevantPrimes(beta) do (
	ListOfLocalAnistropicDimensions = append(ListOfLocalAnistropicDimensions, anisotropicDimensionQp(beta,p))
	);
    max ListOfLocalAnistropicDimensions
    )

-- Input: A symmetric matrix representing a quadratic form or a GrothendieckWittClass over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: An integer, the rank of the anisotropic part of the quadratic form or GrothendieckWittClass

anisotropicDimension = method()
anisotropicDimension (Matrix) := (ZZ) => (A) -> (
    k := ring A;
    -- Ensure base field is supported
    if not (instance(k,ComplexField) or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    -- Ensure matrix is symmetric
    if not isSquareAndSymmetric(A) then error "Matrix is not symmetric";
    -- Over CC, the anisotropic dimension is 0 or 1 depending on the parity of the rank
    if instance(k,ComplexField) then (
        return (rankForm(A)%2);
        )
    -- Over RR, the anisotropic dimension is the absolute value of the signature
    else if instance(k,RealField) then (
        diagonalA := congruenceDiagonalize A;
        return (abs(numPosDiagEntries(diagonalA) - numNegDiagEntries(diagonalA)));
        )
    -- Over QQ, call anisotropicDimensionQQ
    else if (k === QQ) then (
        return anisotropicDimensionQQ(gwClass(nondegeneratePartDiagonal(A)));
        )
    -- Over a finite field, if the number of nonzero diagonal entries is odd, then the anisotropic dimension is 1
    -- if the number of nonzero diagonal entries is even, then the anisotropic dimension is either 0 or 2
    -- depending on whether the nondegenerate part of the form is totally hyperbolic
    else if (instance(k, GaloisField) and k.char != 2) then (
        diagA := congruenceDiagonalize A;
        if (rankForm(diagA)%2 == 1) then (
            return 1;
            )
        else if (legendreBoolean(det(nondegeneratePartDiagonal(diagA))) == legendreBoolean(sub((-1)^(rankForm(diagA)/2),k))) then (
            return 0;
            )
        else (
            return 2;
            );
        );
    )

anisotropicDimension (GrothendieckWittClass) := (ZZ) => (alpha) -> (
    anisotropicDimension(alpha.matrix)
    )

-- Input: A Grothendieck-Witt class over the complex numbers, the real numbers, the rational numbers, or a finite field of characteristic not 2
-- Output: An integer, the Witt index of the class, i.e. the rank of the maximal totally isotropic subspace

WittIndex = method()
WittIndex (GrothendieckWittClass) := (ZZ) => (alpha) -> (
    sub((rankForm(alpha) - anisotropicDimension(alpha))/2,ZZ)
    )
