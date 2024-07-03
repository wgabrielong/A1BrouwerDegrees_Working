-----------------------
-- Producing new forms
-----------------------

-- Input: A field kk, and a list of elements a_1,...,a_n of kk
-- Output: The Grothendieck-Witt class represented by the diagonal form <a_1,...,a_n> 

diagonalForm = method()
diagonalForm(Ring,RingElement) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(Ring,ZZ) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(Ring,QQ) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(Ring, Sequence) := GrothendieckWittClass => (kk,L) -> (
    -- Get the length of the input sequence
    n := #L;
    
    -- Build an n x n mutable identity matrix
    A := mutableIdentity(kk, n);
    
    for i from 0 to (n-1) do (
	A_(i,i) = sub(L_i,kk);
	);
    
    -- A is mutable so we take matrix(A) and form a Grothendieck-Witt class
    gwClass(matrix(A))
    )

diagonalForm(InexactFieldFamily,RingElement) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(InexactFieldFamily,ZZ) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(InexactFieldFamily,QQ) := GrothendieckWittClass => (kk,a) -> (
    gwClass(matrix(kk,{{sub(a,kk)}}))
    )

diagonalForm(InexactFieldFamily, Sequence) := GrothendieckWittClass => (kk,L) -> (
    -- Get the length of the input sequence
    n := #L;
    
    -- Build an n x n mutable identity matrix
    A := mutableIdentity(kk, n);
    
    for i from 0 to (n - 1) do (
	A_(i,i) = sub(L_i,kk);
	);
    
    -- A is mutable so we take matrix(A) and form a Grothendieck-Witt class
    gwClass(matrix(A))
    )

-- Input: A field kk and an optional even rank n (default is n = 2)
-- Output: A Grothendieck-Witt class over kk represented by a totally hyperbolic form or rank n

hyperbolicForm = method()
hyperbolicForm(Ring) := GrothendieckWittClass => (kk) -> (
    gwClass(matrix(kk,{{1,0},{0,-1}}))
    )

hyperbolicForm(Ring,ZZ) := GrothendieckWittClass => (kk,n) -> (
    if odd n then error "entered rank is odd";
    H := matrix(kk,{{1,0},{0,-1}});
    m := sub(n/2,ZZ);
    outputMatrix := diagonalMatrix(kk,{});
    for i from 0 to m - 1 do (
        outputMatrix = outputMatrix ++ H;
        );
    gwClass(outputMatrix)
    )

hyperbolicForm(InexactFieldFamily) := GrothendieckWittClass => (kk) -> (
    gwClass(matrix(kk,{{1,0},{0,-1}}))
    )

hyperbolicForm(InexactFieldFamily,ZZ) := GrothendieckWittClass => (kk,n) -> (
    if odd n then error "entered rank is odd";
    H := matrix(kk,{{1,0},{0,-1}});
    m := sub(n/2,ZZ);
    outputMatrix := diagonalMatrix(kk,{});
    for i from 0 to m - 1 do (
        outputMatrix = outputMatrix ++ H;
        );
    gwClass(outputMatrix)
    )

-- Input: A field kk, and a list of elements a_1,...,a_n of kk
-- Output: The Pfister form <<a_1,...,a_n>>

PfisterForm = method()
PfisterForm(Ring,RingElement) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,(1,sub((-1)*a,kk)))
    )

PfisterForm(Ring,ZZ) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,sub((1,(-1)*a,kk)))
    )

PfisterForm(Ring,QQ) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,(1,sub((-1)*a,kk)))
    )

PfisterForm(Ring,Sequence) := GrothendieckWittClass => (kk,L) -> (
    -- Get the length of the input sequence
    n := #L;
    
    -- Iteratively multiply <1,-L_0>*<1,-L_1>*...
    outputForm := diagonalForm(kk,1);
    for i from 0 to (n-1) do (
	ithPfister := PfisterForm(kk,L_i);
	outputForm = gwMultiply(outputForm,ithPfister);
	);
    outputForm
    )

PfisterForm(InexactFieldFamily,RingElement) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,(1,sub((-1)*a,kk)))
    )

PfisterForm(InexactFieldFamily,ZZ) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,(1,sub((-1)*a,kk)))
    )

PfisterForm(InexactFieldFamily,QQ) := GrothendieckWittClass => (kk,a) -> (
    diagonalForm(kk,(1,sub((-1)*a,kk)))
    )

PfisterForm(InexactFieldFamily,Sequence) := GrothendieckWittClass => (kk,L) -> (
    -- Get the length of the input sequence
    n := #L;
    
    -- Iteratively multiply <1,-L_0>*<1,-L_1>*...
    outputForm := diagonalForm(kk,1);
    for i from 0 to (n-1) do (
	ithPfister := PfisterForm(kk,L_i);
	outputForm = gwMultiply(outputForm,ithPfister);
	);
    outputForm
    )