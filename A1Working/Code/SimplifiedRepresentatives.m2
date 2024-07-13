---------------------
-- Simplifying forms
---------------------

-- Input: A Grothendieck-Witt class beta over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A diagonalized form of beta, with squarefree entries on the diagonal

diagonalClass = method()
diagonalClass (GrothendieckWittClass) := (GrothendieckWittClass) => (beta) -> (

    -- Check if the diagonalClass has already been computed; if so, recall it from the cache
    if beta.cache.?diagonalClass then return beta.cache.diagonalClass;

    diagonalClassOfBetaMatrix := congruenceDiagonalizeSimplify(beta.matrix);

    -- The diagonal form gets cached in the GWClass type
    beta.cache.diagonalClass = gwClass diagonalClassOfBetaMatrix;
    gwClass diagonalClassOfBetaMatrix
    )

-- Input: A Grothendieck-Witt class beta over QQ, RR, CC, or a finite field of characteristic not 2
-- Output: A list of the diagonal entries of a diagonal matrix representing beta

diagonalEntries = method()
diagonalEntries (GrothendieckWittClass) := (List) => (beta) -> (
    
    M := congruenceDiagonalize(beta.matrix);
    n := numRows M;
    L := {};
    
    for i from 0 to (n-1) do (
	L = append(L, M_(i,i));
	);
    L
    )
    
    
