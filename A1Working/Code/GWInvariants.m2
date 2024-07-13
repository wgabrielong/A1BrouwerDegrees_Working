---------------------------------------
-- Invariants
---------------------------------------

-- Input: A Grothendieck-Witt class
-- Output: The rank of a quadratic form representing the Grothendieck-Witt class

rankForm = method()
rankForm (GrothendieckWittClass) := (ZZ) => (alpha) -> (
    numRows(alpha.matrix)
    )

rankForm (Matrix) := (ZZ) => (M) -> (
    rank M
    )



-- Input: A symmetric matrix over QQ or RR
-- Output: The number of positive entries on the diagonal of a diagonal matrix to which it is congruent
-- Note: numPosDiagEntries is *not* included as a method in the A1BrowerDegrees package

numPosDiagEntries = method()
numPosDiagEntries (Matrix) := (Matrix) => (A) -> (
    -- Ensure matrix is symmetric
    if not isSquareAndSymmetric(A) then error "Matrix is not symmetric";
    -- Ensure base field is QQ or RR
    k := ring A;
    if not (instance(k,RealField) or k === QQ) then (
        error "Only implemented over QQ and RR";
        );
    if not isDiagonal(A) then (
        A = congruenceDiagonalize(A);
        );
    posDiagEntries := 0;
    for i from 0 to (numRows(A) - 1) do (
        if A_(i,i) > 0 then (
            posDiagEntries = posDiagEntries + 1;
            );
        );
    posDiagEntries
    )

-- Input: A diagonal matrix over QQ or RR
-- Output: The number of positive entries on the diagonal of a diagonal matrix to which it is congruent
-- Note: numPosDiagEntries is *not* included as a method in the A1BrowerDegrees package

numNegDiagEntries = method()
numNegDiagEntries (Matrix) := (Matrix) => (A) -> (
    -- Ensure matrix is symmetric
    if not isSquareAndSymmetric(A) then error "Matrix is not symmetric";
    -- Ensure base field is QQ or RR
    k := ring A;
    if not (instance(k,RealField) or k === QQ) then (
        error "Only implemented over QQ and RR";
        );
    if not isDiagonal(A) then (
        A = congruenceDiagonalize A;
        );
    negDiagEntries := 0;
    for i from 0 to (numRows(A) - 1) do (
        if A_(i,i) < 0 then (
            negDiagEntries = negDiagEntries + 1;
            );
        );
    negDiagEntries
    )

-- Input: A Grothendieck-Witt class defined over QQ or RR
-- Output: The number of positive entries on the diagonal of a diagonal matrix representing the Grothendieck-Witt class
-- Note: numPosEntries is *not* included as a method in the A1BrowerDegrees package

numPosEntries = method()
numPosEntries (GrothendieckWittClass) := ZZ => beta -> (
    numPosDiagEntries(beta.matrix)
    )

-- Input: A Grothendieck-Witt class beta defined over QQ or RR
-- Output: The number of negative entries on the diagonal of a diagonal matrix representing the Grothendieck-Witt class
-- Note: numNegEntries is *not* included as a method in the A1BrowerDegrees package

numNegEntries = method()
numNegEntries (GrothendieckWittClass) := ZZ => beta -> (
    numNegDiagEntries(beta.matrix)
    )

-- Input: A Grothendieck-Witt class beta defined over QQ or RR
-- Output: The signature of beta

signature = method()
signature (GrothendieckWittClass) := ZZ => (beta) -> (
    numPosEntries(beta) - numNegEntries(beta)
    )

---------------------------
-- Comparing forms over QQ
---------------------------

-- Input: A Grothendieck-Witt class defined over QQ
-- Output: A squarefree integral representative of its discriminant

integralDiscriminant = method()
integralDiscriminant (GrothendieckWittClass) := (ZZ) => (beta) -> (
    kk := baseField beta;
    if not kk === QQ then error "GrothendieckWittClass is not over QQ";

    -- Return a squarefree integral representative of the product of diagonal entries of a diagonal representative of the form 
    squarefreePart product(diagonalEntries(beta))
    )

-- Input: A Grothendieck-Witt class defined over QQ
-- Output: The list of primes that divide entries of its diagonal representative

relevantPrimes = method()
relevantPrimes (GrothendieckWittClass) := List => (beta) -> (
    kk := baseField beta;
    if not kk === QQ then error "GrothendieckWittClass is not over QQ";
    
    -- Find the diagonal entries of a diagonal integral representative of the form
    D := diagonalEntries( diagonalClass beta );
    
    -- Make a list of all prime factors of diagonal entries
    L := {};
    for x in D do (
	L = unique(L | primeFactors(sub(x,ZZ)));
	);
    L
    )

-- Input:  A list of the diagonal elements of a quadratic form (assumed to be rational numbers)
-- or a Grothendieck-Witt class over QQ, and a prime number p
-- Output: The Hasse-Witt invariant of the quadratic form or Grothendieck-Witt class over Q_p

HasseWittInvariant = method()
HasseWittInvariant (List, ZZ) := ZZ => (L,p) -> (
    if not isPrime(p) then error "second argument must be a prime number";

    a := 1;
    len := #L;
       
    -- Replace every entry of L by its squarefree part so we can work with integers
    f := {};
    for x in L do (
	f = append(f,squarefreePart(x));
	);
    for i from 0 to len - 2 do (
       	for j from i + 1 to len - 1 do (
	    a = a * HilbertSymbol(f_i, f_j, p);
	    );
	);
    a
    )

HasseWittInvariant(GrothendieckWittClass, ZZ) := ZZ => (beta,p) -> (
    kk := baseField beta;
    if not (kk === QQ) then error "method is only implemented over the rationals";
    HasseWittInvariant(diagonalEntries(beta),p)
    )

