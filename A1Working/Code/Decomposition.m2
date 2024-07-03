-- Input: A form q over QQ of anisotropic dimension d >= 4
-- Output: A form < a > so that q + < -a > has anisotropic dimension d - 1
-- Note: This is Koprowski/Rothkegel's Algorithm 5 in the case of QQ

QQanisotropicDimension4 = method()
QQanisotropicDimension4 (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    if not (anisotropicDimensionQQ(beta) >= 4) then error "anisotropic dimension of form is not >= 4";
    
    -- If the signature is non-negative then return < 1 >
    if signature(beta) >= 0 then (
	return diagonalForm(QQ,1);
	);
    
    -- Otherwise return < -1 >
    if signature(beta) < 0 then (
	return diagonalForm(QQ,-1);        
        );	
    )

-- Input: A form q over QQ of anisotropic dimension 3
-- Output: A form < a > so that q + < -a > has anisotropic dimension 2
-- Note: This is Koprowski/Rothkegel's Algorithm 7 in the case of QQ

QQanisotropicDimension3 = method()
QQanisotropicDimension3 (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    if not (anisotropicDimensionQQ(beta) == 3) then error "anisotropic dimension of form is not 3";

    d := integralDiscriminant(beta);
    
    -- Build lists of relevant primes where the p-adic valuation of the discriminant is even or is odd
    L1 := {1};
    L2 := {};
    S1 := {1};
    S2 := {};
    for p in relevantPrimes(beta) do (
	if odd padicValuation(d,p) then (
	    L1 = append(L1,p);
	    S1 = append(S1,d-1);
	    );
	if even padicValuation(d,p) then (
	    L2 = append(L2,p^2);
	    S2 = append(S2,p)
	    );
	);
    
    -- We are looking for an element which is equivalent to d-1 mod p for each p in L1 and equivalent to p mod p^2 for each p in L2
    -- We use the chineseRemainder method from the "Parametrization" package to find such an element
    alpha := chineseRemainder(S1 | S2, L1 | L2);
    a := squarefreePart(alpha);
    diagonalForm(QQ,a)
    )

-- Input: A form q over QQ of anisotropic dimension 2
-- Output: The anisotropic part of q
-- Note: This is Koprowski/Rothkegel's Algorithm 8 in the case of QQ

QQanisotropicDimension2 = method()
QQanisotropicDimension2 (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    if not (anisotropicDimensionQQ(beta) == 2) then error "anisotropic dimension of form is not 2";

    n := numRows beta.matrix;

    -- Shortcut: if the form has anisotropic dimension 2 and the form is dimension 2, return the form itself
    if (n==2) then return beta; 
    
    -- Step 1: We want the Witt index to be 0 mod 4 in their terminology --- note they define the Witt index to be
    -- the integer w so that q = wH + q_a.
    w := WittIndex(beta);
    q := beta;
    if ((w % 4) != 0) then (
	w = w % 4;
	q = gwAdd(q, hyperbolicForm(QQ,2*(4-w)));
	n = n + 2*(4-w);
	);

    -- Step 2: Compute discriminant (note they use a signed version of the discriminant in their algorithm)
    d := ((-1)^(n*(n-1)/2))*integralDiscriminant(q);
    
    -- Step 3: Take relevant primes plus dyadic ones
    S := relevantPrimes(beta);
    if not member(2,S) then (
	S = append(S,2);
	);
    
    -- Start the loop at p=2
    p := 2;
    solnFound := false;
    
    while not solnFound do (
	s := #S;

	-- Step 5a: Make a basis for the group of S-singular elements
	basisES = append(S,-1);
	m := #basisES;

    	-- Step 5c: Make a vector of exponents of Hasse invariants
	W := mutableMatrix(QQ,s,1);
	for i from 0 to (s-1) do (
	    W_(i,0) = (1 - (HasseWittInvariant(q,S_i)))/2;
	    );
       	
	-- Step 5b / 5f: 
	W = matrix(W);
    	if (d < 0) then (
	    if abs(signature(q)) != 2 then error "signature isn't pm 2";
	    if (signature(q) == 2) then (
		W = matrix(QQ,{{0}}) || W;
		);
	    if (signature(q) == -2) then (
		W = matrix(QQ,{{1}}) || W;
	        );
	    );
        
    	-- Step 5e: Make a matrix of Hilbert symbols
    	B := mutableMatrix(QQ,s,m);	
    	for i from 0 to (s-1) do (
	    for j from 0 to (m-1) do (
	    	B_(i,j) = (1 - HilbertSymbol(basisES_j, d, S_i))/2;
	    	);
	    );
	B = matrix(B);
    	
	-- Step 5d: Append a zero column on the front if the discriminant is negative
    	if (d < 0) then (
	    A := mutableMatrix(QQ,1,m);
	    for i from 0 to (m-1) do (
	    	if basisES_i > 0 then (
		    A_(0,i) = 0;
		    )
		else (
		    A_(0,i) = 1;
		    );
	    	);
	    B = matrix(A) || B;
	    );

	-- Step 5f: Try to solve sytem of equations over F_2
        kk := GF(2);
    	W = matrix(kk,entries(W));
    	B = matrix(kk,entries(B));

	if (class(solve(B,W)) === Matrix) then (
	    X := solve(B,W);
	    solnFound = true;
	    break;
	    )
	else (
	    p = nextPrime(p+1);
	    while member(p,S) do (
		p = nextPrime(p+1);
		);
	    S = append(S,p);
	    );
	);
    alpha := sub(1,ZZ);
    for j from 0 to (m-1) do (
	alpha = alpha * ((basisES_j)^(sub(X_(j,0),ZZ)));
	);
    diagonalForm(QQ,(alpha, -squarefreePart(alpha*d)))
    )

-- Input: A Grothendieck-Witt class representing a quadratic form over QQ
-- Output: The Grothendieck-Witt class representing the anisotropic part of this form

QQanisotropicPart = method()
QQanisotropicPart (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (
    beta = diagonalClass(beta);
    n := numRows(beta.matrix);
    
    -- If the form is anisotropic 
    if anisotropicDimension(beta) == n then return beta;
    
    -- Initialize an empty quadratic form
    outputForm := diagonalForm(QQ,());    
    alpha := 1;

    while anisotropicDimension(beta) >= 4 do (
	outputForm = gwAdd(outputForm,QQanisotropicDimension4(beta));
	alpha = ((QQanisotropicDimension4(beta)).matrix)_(0,0);	
	beta = gwAdd(beta, diagonalForm(QQ,((-1)*alpha)));
	);
    
    if anisotropicDimension(beta) == 3 then (
	outputForm = gwAdd(outputForm,QQanisotropicDimension3(beta));
	alpha = ((QQanisotropicDimension3(beta)).matrix)_(0,0);	
	beta = gwAdd(beta, diagonalForm(QQ,((-1)*alpha)));
	);
    
    if anisotropicDimension(beta) == 2 then (
        outputForm = gwAdd(outputForm, QQanisotropicDimension2(beta));
        );
    
    if anisotropicDimension(beta) == 1 then (
	outputForm = gwAdd(outputForm, diagonalForm(QQ,((-1)^((n-1)/2))*integralDiscriminant(beta)));
	);
    
    outputForm
    )

-- Input: A symmetric matrix representing a quadratic form or a GrothendieckWittClass; over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A symmetric matrix or GrothendieckWittClass that is the anisotropic part of the input

anisotropicPart = method()
anisotropicPart (Matrix) := (Matrix) => (A) -> (
    k := ring A;
    -- Ensure base field is supported
    if not (instance(k,ComplexField) or instance(k,RealField) or k === QQ or (instance(k, GaloisField) and k.char != 2)) then (
        error "Base field not supported; only implemented over QQ, RR, CC, and finite fields of characteristic not 2";
        );
    -- Ensure underlying matrix is symmetric
    if not isSquareAndSymmetric(A) then error "Underlying matrix is not symmetric";
    -- Over CC, the anisotropic part is either the rank 0 form or the rank 1 form, depending on the anisotropic dimension
    if instance(k,ComplexField) then (
        if (anisotropicDimension(A) == 0) then (
            return (diagonalMatrix(CC,{}));
            )
        else (
            return (matrix(CC,{{1}}));
            );
        )
    --Over RR, the anisotropic part consists of the positive entries in excess of the number of negative entries, or vice versa
    else if instance(k,RealField) then (
        diagonalA := congruenceDiagonalize(A);
        posEntries := numPosDiagEntries(diagonalA);
        negEntries := numNegDiagEntries(diagonalA);
        if (posEntries > negEntries) then (
            return (id_(RR^(posEntries-negEntries)));
            )
        else if (posEntries < negEntries) then (
            return (-id_(RR^(negEntries-posEntries)));
            )
        else (
            return (diagonalMatrix(RR,{}));
            );
        )
    -- Over QQ, call anisotropicPartQQ
    else if (k === QQ) then (
        return (QQanisotropicPart(gwClass(nondegeneratePartDiagonal(A)))).matrix;
        )
    -- Over a finite field, if the anisotropic dimension is 1, then the form is either < 1 > or < e >, where e is any nonsquare representative,
    -- and if the anisotropic dimension is 2 then the form is <1,-e>
    else if (instance(k, GaloisField) and k.char != 2) then (
        diagA := congruenceDiagonalize(A);
        if (anisotropicDimension(A) == 1) then (
            return (matrix(k,{{sub((-1)^((rank(diagA)-1)/2),k)*det(nondegeneratePartDiagonal(diagA))}}));
            )
        else if (anisotropicDimension(A) == 0) then (
            return (diagonalMatrix(k,{}));
            )
        else (
            return (matrix(k,{{1,0},{0,sub((-1)^((rank(diagA)-2)/2),k)*det(nondegeneratePartDiagonal(diagA))}}));
            );
        );
    )

anisotropicPart (GrothendieckWittClass) := (GrothendieckWittClass) => (alpha) -> (
    gwClass(anisotropicPart(alpha.matrix))
    )

---------------------------------------
-- Simplifying a form
---------------------------------------

-- Input: A Grothendieck-Witt class beta over a field kk
-- Output: A simplified diagonal representative of beta
sumDecompositionVerbose = method()
sumDecompositionVerbose (GrothendieckWittClass) := (GrothendieckWittClass, String) => beta -> (
    -- Get base field of beta
    kk := baseField(beta);

    if numRows(beta.matrix) == 0 then (
	return (gwClass(diagonalMatrix(kk,{})),"empty form");
	);
    
    outputString := "";
    
    -- Get isotropic dimension of beta and construct its isotropic and anistropic parts
    w := WittIndex(beta);
    
    if w > 0 then (
	outputString = outputString | toString(w) | "H";
	);
    
    hyperbolicPart := hyperbolicForm(kk,2*w);
    alpha := anisotropicPart(beta);
    
    if numRows(alpha.matrix) > 0 then (
        D := diagonalEntries(alpha);
        for i from 0 to (length(D)-1) do (
	    outputString = outputString | "+ <" | toString(D_i) | ">";
            );
	);
    
    -- Return a simplified form of beta
    return (gwAdd(alpha,hyperbolicPart),outputString);
    )    

-- Input: A Grothendieck-Witt class beta over a field k
-- Output: A simplified diagonal representative of beta

sumDecomposition = method()
sumDecomposition (GrothendieckWittClass) := (GrothendieckWittClass) => beta -> (
    beta.cache.diagonalClass = (sumDecompositionVerbose(beta))_0;
    (sumDecompositionVerbose(beta))_0
    )

-- Input: A Grothendieck-Witt class beta over a field k
-- Output: The decomposition as a sum of hyperbolic and rank one forms

sumDecompositionString = method()
sumDecompositionString (GrothendieckWittClass) := (String) => beta -> (
    (sumDecompositionVerbose(beta))_1
    )