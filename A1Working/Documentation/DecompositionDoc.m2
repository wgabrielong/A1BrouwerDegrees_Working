document {
    Key => {getSumDecomposition, (getSumDecomposition, GrothendieckWittClass)},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "getSumDecomposition beta",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///},
    },
    Outputs => { 
	GrothendieckWittClass => {"a diagonal representative of the Grothendieck Witt class of the input form"},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we decompose it as a sum of some number of hyperbolic and rank one forms."},
    EXAMPLE lines ///
    M = matrix(RR, {{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = makeGWClass M;
    getSumDecomposition beta
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13), {{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = makeGWClass M;
    getSumDecomposition beta
    ///,
    PARA{EM "Citations:"},
    UL{
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
    SeeAlso => {"getSumDecompositionString"},
}

document {
    Key => {getSumDecompositionString, (getSumDecompositionString, GrothendieckWittClass)},
    Headline => "produces a simplified diagonal representative of a Grothendieck Witt class",
    Usage => "getSumDecompositionString beta",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///},
    },
    Outputs => { 
	String => {"the decomposition as a sum of hyperbolic and rank one forms"},
	},
    PARA {"Given a symmetric bilinear form ", TT"beta", " over a field ", TEX///$k$///, ", we return a simplified diagonal form of ", TT"beta","."},
    EXAMPLE lines ///
    M = matrix(RR, {{2.091,2.728,6.747},{2.728,7.329,6.257},{6.747,6.257,0.294}});
    beta = makeGWClass M;
    getSumDecompositionString beta
    ///,
    PARA {"Over ", TEX///$\mathbb{R}$///, " there are only two square classes and a form is determined uniquely by its rank and signature [L05, II Proposition 3.2]. A form defined by the ", TEX///$3\times 3$///, " Gram matrix ", TT"M", " above is isomorphic to the form ", TEX///$\langle 1,-1,1\rangle $///, "."},
    EXAMPLE lines ///
    M = matrix(GF(13), {{9,1,7,4},{1,10,3,2},{7,3,6,7},{4,2,7,5}});
    beta = makeGWClass M;
    getSumDecompositionString beta
    ///,
    PARA {"Over ", TEX///$\mathbb{F}_{q}$///, " forms can similarly be diagonalized. In this case as ", TEX///$\langle 1,-1,1,-6 \rangle$///, "."},
    PARA{EM "Citations:"},
    UL{
	{"[L05] T.Y. Lam, ", EM "Introduction to quadratic forms over fields,", " American Mathematical Society, 2005."},
	},
    SeeAlso => {"getSumDecomposition"},
}


document {
    Key => {getAnisotropicPart, (getAnisotropicPart, GrothendieckWittClass), (getAnisotropicPart, Matrix)},
    Headline => "returns the anisotropic part of a Grothendieck Witt class",
    Usage => "getAnisotropicPart beta",
    Inputs => {
        GrothendieckWittClass => "beta" => {"a symmetric bilinear form defined over a field ", TEX///$k$///},
    },
    Outputs => { 
	GrothendieckWittClass => {"the anisotropic part of ", TEX///$\beta$///},
    },
    PARA {"Given a form ", TEX///$\beta$///, " we may compute its anisotropic part inductively by reference to its ", TO2(getAnisotropicDimension,"anisotropic dimension"), ". Over the complex numbers and the reals this is trivial, and over finite fields it is a fairly routine computation, however over the rationals some more sophisticated algorithms are needed from the literature. For this methods we implement algorithms developed for number fields by Koprowski and Rothkegel [KR23]."},
    EXAMPLE lines ///
    alpha = makeDiagonalForm(QQ, (3,-3,2,5,1,-9));
    getAnisotropicPart alpha
    ///,
    PARA{EM "Citations:"},
    UL{
	{"[KR23] P. Koprowski and B. Rothkegel, ", EM "The anisotropic part of a quadratic form over a number field,", " Journal of Symbolic Computatoin, 2023."},
	},
    SeeAlso => {"getAnisotropicDimension", "getWittIndex"},
}
