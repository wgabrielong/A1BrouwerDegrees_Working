-- We define GrothendieckWittClass to be a new type, meant to represent the isomorphism class 
-- of a nondegenerate symmetric bilinear form over a field of characteristic not 2

GrothendieckWittClass = new Type of HashTable
GrothendieckWittClass.synonym = "Grothendieck Witt Class"

-- Input: A matrix M representing a nondegenerate symmetric bilinear form over a field of characteristic not 2
-- Output: The GrothendieckWittClass representing the symmetric bilinear form determined by M

gwClass = method()
gwClass (Matrix) := GrothendieckWittClass => M -> (
   if isWellDefined M then (
        new GrothendieckWittClass from {
            symbol matrix => M,
            symbol cache => new CacheTable
            }
        )
    else (
        error "gwClass called on a matrix that does not represent a nondegenerate symmetric bilinear form over a field of characteristic not 2";
        )
    )

-- This allows us to extract the matrix from a Grothendieck-Witt class
matrix GrothendieckWittClass := Matrix => beta -> beta.matrix

-- Input: A matrix
-- Output: Boolean that gives whether the matrix defines a nondegenerate symmetric bilinear form over a field of characteristic not 2

isWellDefined (Matrix) := Boolean => M -> (
    
    -- Return false if the matrix isn't square and symmetric
    if not isSquareAndSymmetric(M) then return false;

    -- Return false if the matrix represents a degenerate form
    if isDegenerate M then return false;

    -- Return false if the matrix isn't defined over a field
    if not isField(ring M) then return false;
    
    -- Returns false if the matrix is defined over a field of characteristic 2
    if char(ring M) == 2 then return false;

    -- Otherwise, return true
    true
    )

-- Input: A Grothendieck-Witt class beta, the isomorphism class of a symmetric bilinear form
-- Output: The base ring of beta

baseField = method()
baseField GrothendieckWittClass := Ring => beta -> (
    ring(beta.matrix)
    )

-- Input: Two Grothendieck-Witt classes beta and gamma
-- Output: The direct sum of beta and gamma

gwAdd = method()
gwAdd(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    Kb := baseField beta;
    Kg := baseField gamma;
    
    -- Galois field case
    if instance(Kb, GaloisField) and instance(Kg, GaloisField) then (
	-- Return an error if the underlying fields of the two classes are different
	if not Kb.order == Kg.order  then error "these classes have different underlying fields";
	return gwClass(beta.matrix ++ substitute(gamma.matrix,Kb));
	);
    
    -- Remaining cases
    if not Kb === Kg then error "these classes have different underlying fields";
    gwClass(beta.matrix ++ gamma.matrix)
    )

-- Input: Two Grothendieck-Witt classes beta and gamma
-- Output: The tensor product of beta and gamma

gwMultiply = method()
gwMultiply(GrothendieckWittClass, GrothendieckWittClass) := GrothendieckWittClass => (beta, gamma) -> (
    Kb := baseField beta;
    Kg := baseField gamma;
    
    -- Galois field case
    if instance(Kb, GaloisField) and instance(Kg, GaloisField) then (
	-- Return an error if the underlying fields of the two classes are different
	if not Kb.order == Kg.order  then error "these classes have different underlying fields";
	return gwClass(beta.matrix ** substitute(gamma.matrix,Kb));
	);
    
    -- Remaining cases
    if not Kb === Kg then error "these classes have different underlying fields";
    gwClass(beta.matrix ** gamma.matrix)
    )
