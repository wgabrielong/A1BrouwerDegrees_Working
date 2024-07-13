-- Input: Two Grothendieck-Witt classes or symmetric matrices representing quadratic forms over QQ
-- Output: Boolean that gives whether the Grothendieck-Witt classes/quadratic forms are isomorphic

isIsomorphicFormQ = method()
isIsomorphicFormQ (GrothendieckWittClass, GrothendieckWittClass) := Boolean => (alpha, beta) -> (
    if not (baseField(alpha) === QQ) then error "first input must have base field QQ";
    if not (baseField(beta) === QQ) then error "second input must have base field QQ";
    
    -- If the ranks differ, then the forms are not isomorphic
    if rankForm(alpha) != rankForm(beta) then return false;
    
    -- If the signatures (Hasse-Witt invariants at RR) differ, then the forms are not isomorphic
    if signature(alpha) != signature(beta) then return false;
    
    -- If the discriminants differ, then the forms are not isomorphic
    if integralDiscriminant(alpha) != integralDiscriminant(beta) then return false;
    
    -- Check the Hasse-Witt invariants
    PrimesToCheck := unique(relevantPrimes(alpha) | relevantPrimes(beta));
    for p in PrimesToCheck do (
        -- If any Hasse-Witt invariants differ, then the forms are not isomorphic
	if (HasseWittInvariant(alpha,p) != HasseWittInvariant(beta,p)) then (
	    return false;
	    );
	);
    -- If we get here, then all Hasse-Witt invariants agree and the forms are isomorphic
    true
    )

isIsomorphicFormQ (Matrix, Matrix) := Boolean => (M,N) -> (
    isIsomorphicFormQ(gwClass(M),gwClass(N))
    )

-- Input: Two Grothendieck-Witt classes alpha and beta, defined over CC, RR, QQ, or a finite field of characteristic not 2
-- Output: Boolean that gives whether alpha and beta are the same Grothendieck-Witt class

gwIsomorphic = method()
gwIsomorphic (GrothendieckWittClass,GrothendieckWittClass) := (Boolean) => (alpha,beta) -> (
    isIsometricForm(alpha.matrix,beta.matrix)
    )

-- Input: Two matrices representing symmetric bilinear forms over CC, RR, QQ, or a finite field of characteristic not 2
-- Output: Boolean that gives whether the bilinear forms are isometric

isIsometricForm = method()
isIsometricForm (Matrix,Matrix) := (Boolean) => (A,B) -> (
    k1 := ring A;
    k2 := ring B;
    -- Ensure both base fields are supported
    if not (instance(k1,ComplexField) or instance(k1,RealField) or k1 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    if not (instance(k2,ComplexField) or instance(k2,RealField) or k2 === QQ or (instance(k1, GaloisField) and k1.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    
    -- Ensure both matrices are symmetric
    if not isSquareAndSymmetric(A) then error "Underlying matrix is not symmetric";
    if not isSquareAndSymmetric(B) then error "Underlying matrix is not symmetric";
    
    -----------------------------------
    -- Complex numbers
    -----------------------------------
    
    -- Over CC, forms over spaces of the same dimension are equivalent if and only if they have the same rank
    if (instance(k1,ComplexField) and instance(k2,ComplexField)) then (
        return ((numRows(A) == numRows(B)) and (rankForm(A) == rankForm(B)));
        )
    
    -----------------------------------
    -- Real numbers
    -----------------------------------
    
    -- Over RR, diagonal forms of the same dimension are equivalent if and only if they have the same number of positive and negative entries
    else if (instance(k1,RealField) or instance(k2,RealField)) then (
        diagA := congruenceDiagonalize A;
        diagB := congruenceDiagonalize B;
        return ((numRows(A) == numRows(B)) and (numPosDiagEntries(diagA) == numPosDiagEntries(diagB)) and (numNegDiagEntries(diagA) == numNegDiagEntries(diagB)));
        )
    
    -----------------------------------
    -- Rational numbers
    -----------------------------------
    
    -- Over QQ, if spaces have same dimension, then call isIsomorphicFormQ on their nondegenerate parts
    else if ((k1 === QQ) and (k2 === QQ)) then (
        return ((numRows(A) == numRows(B)) and (isIsomorphicFormQ(nondegeneratePartDiagonal(A),nondegeneratePartDiagonal(B))));
        )
    
    -----------------------------------
    -- Finite fields
    -----------------------------------
    
    -- Over a finite field, diagonal forms over spaces of the same dimension are equivalent if and only if they have the same number of nonzero entries and the product of these nonzero entries is in the same square class
    else if (instance(k1, GaloisField) and instance(k2, GaloisField) and k1.char !=2 and k2.char != 2 and k1.order == k2.order) then (
        return ((numRows(A) == numRows(B)) and (rankForm(A) == rankForm(B)) and (legendreBoolean(det(nondegeneratePartDiagonal(A))) == legendreBoolean(sub(det(nondegeneratePartDiagonal(B)),k1))));
        )
    -- If we get here, the base fields are not the same
    else error "Base fields are not the same"
    )
