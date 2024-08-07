document{
    Headline => "outputs the signature of a symmetric bilinear form over the real or rational numbers",
    Key => {getSignature, (getSignature, GrothendieckWittClass)},
    Usage => "getSignature beta",
    Inputs => {
	GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over ", TEX///$\mathbb{Q}$///, " or ", TEX///$\mathbb{R}$///},
	},
    Outputs => {
	ZZ => "n" => {"the ", EM "signature", " of the symmetric bilinear form ", TEX///$\beta$///},
	},
    PARA{"Given a symmetric bilinear form, after diagonalizing it, we can consider the number of positive entries minus the number of negative entries appearing along the diagonal. This is the ", EM "signature", " of a symmetric bilinear form, and is one of the primary invariants we use to classify forms. For more information see ", TO2(isIsomorphicForm,"isIsomorphicForm"), "."},
    EXAMPLE lines ///
    M = matrix(RR, {{0,0,1},{0,1,0},{1,0,0}});
    beta = makeGWClass M;
    getSignature beta
    ///,
    SeeAlso => {"isIsomorphicForm", "getSumDecomposition", "getSumDecompositionString"}
    }


document{
    Key => {getIntegralDiscriminant, (getIntegralDiscriminant, GrothendieckWittClass)},
    Headline => "outputs an integral discriminant for a rational symmetric bilinear form",
    Usage => "getIntegralDiscriminant beta",
    Inputs => {
	GrothendieckWittClass => "beta" => {"denoted by ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///},
	},
    Outputs => {
        ZZ => {"an integral square class representative of ", TEX///$\text{disc}(\beta)$///},
	},
    EXAMPLE lines ///
    beta = makeGWClass matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}});
    getIntegralDiscriminant beta
    getDiagonalClass beta
    ///,
}


document{
    Key => {getHasseWittInvariant, (getHasseWittInvariant, GrothendieckWittClass, ZZ), (getHasseWittInvariant, List, ZZ)},
    Headline => "outputs the Hasse-Witt invariant for a prime p for the quadratic form of the Grothendieck-Witt class",
    Usage => "getHasseWittInvariant(beta, p)",
    Inputs => {
	GrothendieckWittClass => "beta" => {"denoted by ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///},
	ZZ => "p" => {"an integral prime number"},
	},
    Outputs => {
        ZZ => {"the Hasse-Witt invariant for ", TEX///$\beta$///," for the prime ",TEX///$p$///},
	},
    PARA{"The ", EM "Hasse-Witt invariant", " of a diagonal form ", TEX///$\langle a_1,\ldots,a_n\rangle$///, " over a field ", TEX///$K$///, " is defined to be the product ", TEX///$\prod_{i<j} \left((a_i,a_j)_p \right)$///, " where ", TEX///$(-,-)_p$///, " is the ", TO2(getHilbertSymbol,"Hilbert symbol"), "."},
    PARA{"The Hasse-Witt invariant of a form will be equal to 1 for almost all primes. In particular after diagonalizing a form ", TEX///$\beta \cong \left\langle a_1,\ldots,a_n\right\rangle$///, " then the Hasse-Witt invariant at a prime ", TEX///$p$///, " will be 1 automatically if ", TEX///$p\nmid a_i$///, " for all ", TEX///$i$///, ". Thus we only have to compute the invariant at ", TO2(getRelevantPrimes, "primes dividing diagonal entries"),  "."},
    EXAMPLE lines ///
    beta = makeGWClass matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}});
    getHasseWittInvariant(beta, 7)
    ///,
}


document{
    Key => {getRelevantPrimes, (getRelevantPrimes, GrothendieckWittClass)},
    Headline => "outputs a list of primes at which the Hasse-Witt invariants of a symmetric bilinear form may be non-trivial",
    Usage => "getRelevantPrimes beta",
    Inputs => {
	GrothendieckWittClass => "beta" => {"denoted by ", TEX///$\beta\in\text{GW}(\mathbb{Q})$///},
	},
    Outputs => {
        List => {"a finite list of primes ", TEX///$(p_1,\ldots,p_r)$///, " for which the Hasse-Witt invariants ", TEX///$\phi_p(\beta)$///," may be nontrivial"},
	},
    PARA{"It is a classical result that the ", TO2(getHasseWittInvariant,"Hasse-Witt invariants"), " of a quadratic form are equal to 1 for all but finitely many primes (see e.g. [S73, IV Section 3.3]. As the Hasse-Witt invariants are computed as a product of ", TO2(getHilbertSymbol,"Hilbert symbols") , " of the pairwise entries appearing on a diagonalization of the symbol, it suffices to consider primes dividing diagonal entries."},
   EXAMPLE lines ///
   beta = makeDiagonalForm(QQ, (6,7,22));
   getRelevantPrimes beta
   ///,
   PARA{EM "Citations:"},
    UL{
	{"[S73] J.P. Serre, ", EM "A course in arithmetic,", " Springer-Verlag, 1973."},
    },
}

document {
        Key => {getRank, (getRank, GrothendieckWittClass), (getRank, Matrix)},
        Headline => "Calculates the rank of a symmetric bilinear form",
        Usage => "getRank beta",
        Inputs => {
            GrothendieckWittClass => "beta" => {"denoted by ",  TEX///$\beta \in  GW(\mathbb{Q})$///, " or a symmetric matrix ", TEX///$\beta$///}
            },
        Outputs => {
            ZZ => {" the rank of the form ", TEX///$\beta$///}
           },
        PARA {"This computes the rank of the form given by ", TEX///$\beta$/// },
        EXAMPLE lines ///
                 beta = makeDiagonalForm(QQ, (3,5,7,11))
                 getRank beta
		 ///,                
        EXAMPLE lines ///
                 M=matrix(QQ, {{1,4,7},{4,3,-1},{7,-1,5}})
                 getRank M
                 ///,
}
