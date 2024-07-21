--A1BrouwerDegrees.m2
newPackage(
    "A1Working",
    Version=>"1.0",
    Date=>"June 5, 2023",
    Authors=>{
        {Name=>"Nikita Borisov",
	 Email=>"nborisov@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/nikita-borisov"},
        {Name=>"Thomas Brazelton",
	 Email=>"brazelton@math.harvard.edu",
	 HomePage=>"https://tbrazel.github.io/"},
        {Name=>"Frenly Espino",
	 Email=>"frenly@sas.upenn.edu",
	 HomePage=>"https://www.math.upenn.edu/people/frenly-espino"},
         {Name=>"Tom Hagedorn",
	 Email=>"hagedorn@tcnj.edu",
	 HomePage=>"https://hagedorn.pages.tcnj.edu/"},
        {Name=>"Zhaobo Han",
	 Email=>"zbtomhan@sas.upenn.edu",
	 HomePage=>"https://www.linkedin.com/in/zhaobo-han-77b1301a2/"},
     	{Name=>"Jordy Lopez Garcia",
	 Email=>"jordy.lopez@tamu.edu",
	 HomePage=>"https://jordylopez27.github.io/"},
        {Name=>"Joel Louwsma",
	 Email=>"jlouwsma@niagara.edu",
	 HomePage=>"https://www.joellouwsma.com/"},
        {Name=>"Andrew Tawfeek",
	 Email=>"atawfeek@uw.edu",
	 HomePage=>"https://www.atawfeek.com/"},
        {Name=>"Wern Juin Gabriel Ong",
	 Email=>"gong@bowdoin.edu",
	 HomePage=>"https://wgabrielong.github.io/"}
	},
    Headline=>"for working with A1Working degree computations",
    PackageImports=>{
	"Parametrization",
	"RealRoots",
	"RationalPoints2"
	},
    PackageExports=>{},
    DebuggingMode=>true
    )


export{
    
    -- ArithmeticMethods.m2
    "getLocalAlgebraBasis",
    "getPadicValuation",
    
    --MatrixMethods.m2
    "diagonalizeViaCongruence",
    
    --GrothendieckWittClasses.m2    
    "GrothendieckWittClass",
    "getMatrix",
    "getBaseField",
    "makeGWClass",
    "addGW",
    "multiplyGW",
    
    --BuildingForms.m2
    "makeDiagonalForm",
    "makeHyperbolicForm",
    "makePfisterForm",
    
    --SimplifiedRepresentatives.m2
    "getDiagonalClass",
    "getDiagonalEntries",
    
    --getHilbertSymbols.m2
    "getHilbertSymbol",
    "getHilbertSymbolReal",
    
    --GWInvariants.m2
    "getSignature",
    "getRank",
    "getIntegralDiscriminant",
    "getRelevantPrimes",
    "getHasseWittInvariant",

    --LocalGlobalDegrees.m2
    "getGlobalA1Degree",
    "getLocalA1Degree",
    
    --IsomorphismOfForms.m2
    "isIsomorphicForm",
    
    --Isotropy.m2
    "isIsotropic",
    "isAnisotropic",

    --AnisotropicDimension.m2
    "getAnisotropicDimensionQQp",
    "getAnisotropicDimension",
    "getWittIndex",
    
    --Decomposition.m2
    "getAnisotropicPart",
    "getSumDecomposition",
    "getSumDecompositionString"
    
    }

-- Basic arithmetic, p-adic, and commutative algebra operations we will use
load "./A1Working/Code/ArithmeticMethods.m2"

-- Basic manipulations of matrices we will use
load "./A1Working/Code/MatrixMethods.m2"

-- Establishing the GrothendieckWittClass type and some basic manipulations
load "./A1Working/Code/GrothendieckWittClasses.m2"

-- For building new symmetric bilinear forms
load "./A1Working/Code/BuildingForms.m2"

-- For providing simplified representatives of symmetric bilinear forms
load "./A1Working/Code/SimplifiedRepresentatives.m2"

-- For Hilbert symbols over p-adic numbers
load "./A1Working/Code/HilbertSymbols.m2"

-- Invariants of symmetric bilinear forms
load "./A1Working/Code/GWInvariants.m2"
    
-- Local and global A1-brouwer degrees
load "./A1Working/Code/LocalGlobalDegrees.m2"

-- Checking if forms are isomorphic
load "./A1Working/Code/IsomorphismOfForms.m2"

-- For verifying (an)isotropy
load "./A1Working/Code/Isotropy.m2"

-- Anisotropic dimension
load "./A1Working/Code/AnisotropicDimension.m2"

-- Finally, decomposing forms
load "./A1Working/Code/Decomposition.m2"

----------------------------
----------------------------
-- DOCUMENTATION
----------------------------
----------------------------

beginDocumentation()

document{
    Key => A1Working,
    Headline => "for working with A1Working degree computations",
    PARA{"This package is intended computing and manipulating ", TO2(getLocalA1Degree,"local"), " and ", TO2(getGlobalA1Degree,"global"), " ", TEX///$\mathbb{A}^1$///, EM "-Brouwer degrees."," Global Brouwer degrees are non-degenerate symmetric bilinear forms valued in the Grothendieck-Witt ring of a field ", TEX///$\text{GW}(k)$///, "."},
    PARA{"In order to simplify the forms produced, this package produces invariants of symmetric bilinear forms, including their ", TO2(getWittIndex,"Witt indices"), ", their ", TO2(getIntegralDiscriminant,"discriminants"), ", and their ", TO2(getHasseWittInvariant, "Hasse Witt invariants"), ". Quadratic forms can furthermore be ", TO2(getSumDecomposition,"decomposed"), " into their isotropic and ", TO2(getAnisotropicPart,"anisotropic parts"), ". Finally, and perhaps most crucially, we can certify whether two symmetric bilinear forms are ", TO2(isIsomorphicForm,"isomorphic") , " in the Grothendieck-Witt ring."},
    }

undocumented {
    }

load "./A1Working/Documentation/ArithmeticMethodsDoc.m2"

load "./A1Working/Documentation/MatrixMethodsDoc.m2"

load "./A1Working/Documentation/GrothendieckWittClassesDoc.m2"

load "./A1Working/Documentation/BuildingFormsDoc.m2"

load "./A1Working/Documentation/SimplifiedRepresentativesDoc.m2"

load "./A1Working/Documentation/HilbertSymbolsDoc.m2"

load "./A1Working/Documentation/GWInvariantsDoc.m2"

load "./A1Working/Documentation/LocalGlobalDegreesDoc.m2"

load "./A1Working/Documentation/IsomorphismOfFormsDoc.m2"

load "./A1Working/Documentation/IsotropyDoc.m2"

load "./A1Working/Documentation/AnisotropicDimensionDoc.m2"

load "./A1Working/Documentation/DecompositionDoc.m2"

----------------------------
----------------------------
-- Testing
----------------------------
----------------------------


-- Diagonal form testing
-- Test 0
TEST ///
M1 = matrix(RR, {{0,1},{1,0}});
G1 = makeGWClass M1;
M2 = getDiagonalClass G1;
assert(getMatrix(M2) === matrix(RR, {{1,0},{0,-1}}));
///

-- Test 1
TEST ///
M3 = matrix(CC, {{1,2,3},{2,4,5},{3,5,7}});
G2 = makeGWClass M3;
M4 = getDiagonalClass G2;
assert(getMatrix(M4) === matrix(CC, {{1,0,0},{0,1,0},{0,0,1}}));
///

--Test 2
TEST ///
M3 = matrix(QQ, {{1,2,3},{2,4,5},{3,5,7}});
G2 = makeGWClass M3;
M4 = getDiagonalClass G2;
assert(getMatrix(M4) === matrix(QQ, {{1,0,0},{0,-2,0},{0,0,2}}));
///

-- Tests for GrothendieckWittClass type
-- Test 3
TEST ///
M = matrix(QQ, {{1,0},{0,1}});
N = matrix(QQ, {{1,2},{2,5}});
beta = makeGWClass M;
gamma = makeGWClass N;
assert(getMatrix(beta) === M);
assert(getBaseField(beta) === QQ);
--Operations within GW-classes
A = addGW(beta, gamma);
B = multiplyGW(beta, gamma);
assert(getMatrix(A) === matrix(QQ, {{1,0,0,0},{0,1,0,0},{0,0,1,2},{0,0,2,5}}));
assert(getMatrix(B) === matrix(QQ, {{1,2,0,0},{2,5,0,0},{0,0,1,2},{0,0,2,5}}));
///

-- Testing for global and local A1 degrees
-- Test 4
TEST ///
T1 = QQ[x];
f = {x^2};
beta = getGlobalA1Degree f;
gamma = makeGWClass matrix(QQ, {{0,1},{1,0}});
assert(isIsomorphicForm(beta, gamma));
///

-- Test 5
TEST ///
T1 = QQ[z_1..z_2];
f1 = {(z_1 - 1)*z_1*z_2, (3/5)*z_1^2 - (17/3)*z_2^2};
f1GD = getGlobalA1Degree f1;
assert(getWittIndex(f1GD) == 3);
q = ideal(z_1, z_2);
r = ideal(z_1 - 1, z_2^2 - 9/85);
f1LDq = getLocalA1Degree(f1, q);
f1LDr = getLocalA1Degree(f1, r);
f1LDsum = addGW(f1LDq, f1LDr);
assert(isIsomorphicForm(f1LDsum, f1GD));
///

-- Test 6
TEST ///
T2 = GF(17)[w];
f2 = {w^4 + w^3 - w^2 - w};
f2GD = getGlobalA1Degree f2;
assert(getWittIndex(f2GD) == 2);
p = ideal(w + 1);
f2LDp = getLocalA1Degree(f2, p);
assert(getWittIndex(f2LDp) == 1);
s = ideal(w - 1);
f2LDs = getLocalA1Degree(f2, s);
t = ideal(w);
f2LDt = getLocalA1Degree(f2, t);
f2LDsum = addGW(addGW(f2LDp, f2LDs), f2LDt);
assert(isIsomorphicForm(f2LDsum, f2GD));
///

-- Testing for building forms
-- Test 7
TEST ///
twoH = makeHyperbolicForm(GF(17), 4);
P = makePfisterForm(GF(17), (2,3));
assert(isIsomorphicForm(P, twoH));
///

-- Test 8
TEST ///
H = makeHyperbolicForm RR;
A = makeDiagonalForm(RR, (1,-1));
B = makeGWClass(matrix(RR, {{0,1},{1,0}}));
assert(isIsomorphicForm(H, A));
assert(isIsomorphicForm(H, B));
///

-- Test for local algebra basis
-- Test 9
TEST ///
QQ[x,y];
f = {x^2 + 1 - y, y};
p = ideal(x^2 + 1, y);
assert(getLocalAlgebraBasis(f, p) == {1,x}); 
///

-- Tests for getDiagonalClass and getDiagonalEntries
-- Test 10
TEST ///
M1 = matrix(CC, {{1,0,0},{0,2,0},{0,0,-3}});
M2 = matrix(CC, {{1,0,0},{0,1,0},{0,0,1}});
G = makeGWClass M1;
assert(getMatrix(getDiagonalClass G) == M2);
assert(getDiagonalEntries(G) == {1,2,-3});
///

-- Test 11
TEST ///
M1 = matrix(RR, {{1,0,0},{0,2,0},{0,0,-3}});
M2 = matrix(RR, {{1,0,0},{0,1,0},{0,0,-1}});
G = makeGWClass M1;
assert(getMatrix(getDiagonalClass G) == M2);
assert(getDiagonalEntries(G) == {1,2,-3});
///

-- Test 12
TEST ///
M = matrix(QQ, {{1,0,0},{0,2,0},{0,0,-3}});
G = makeGWClass M;
assert(getMatrix(getDiagonalClass G) == M);
assert(getDiagonalEntries(G) == {1,2,-3});
///
    
-- Test 13
TEST ///
M = matrix(GF(5), {{1,0,0},{0,2,0},{0,0,-3}});
G = makeGWClass M;
assert(getMatrix(getDiagonalClass G) == M);
assert(getDiagonalEntries(G) == {1,2,-3});
///

-- Test 14
TEST ///
kk = GF(7);
M1 = matrix(kk, {{1,0,0},{0,2,0},{0,0,-3}});
M2 = matrix(kk, {{1,0,0},{0,1,0},{0,0,1}});
G = makeGWClass M1;
assert(getMatrix(getDiagonalClass G) == M2);
assert(getDiagonalEntries(G) == {1,2,-3});
///

-- Test 15
TEST ///
M1 = matrix(QQ, {{18,0,0},{0,125/9,0},{0,0,-8/75}});
M2 = matrix(QQ, {{2,0,0},{0,5,0},{0,0,-6}});
G1 = makeGWClass M1;
assert(getMatrix(getDiagonalClass G1) == M2);
///

-- Test for p-adic valuation
-- Test 16
TEST ///
assert(getPadicValuation(27, 3) == 3);
///

-- Test for diagonalizeViaCongruence
-- Test 17
TEST ///
B = matrix(QQ, {{0/1,1},{1,0}});
eta = makeGWClass B;
assert(getWittIndex(eta) == 1);

P = matrix(QQ, {{0/1, 5,1},{2,2,1},{0,0,1}});
A = matrix(QQ, {{1/1,0,0},{0,-1,0},{0,0,1}});
assert(getWittIndex(makeGWClass(diagonalizeViaCongruence(P*A*transpose(P)))) == 1);
///

-- Test for makeGWClass
-- Test 18
TEST ///
M1 = matrix(QQ, {{1/1,0,0},{0,1,0},{0,0,1}});
M2 = matrix(QQ, {{1/1,24/10,0},{24/10,-5,0},{0,0,69}});
M3 = matrix(GF(7), {{1,0,0},{0,2,0},{0,0,-3}});
assert(class(makeGWClass M1) === GrothendieckWittClass);
assert(class(makeGWClass M2) === GrothendieckWittClass);
assert(class(makeGWClass M3) === GrothendieckWittClass);
///

-- Test for getBaseField
-- Test 19
TEST ///
M = makeGWClass matrix(QQ, {{1/1,0,0},{0,2,3},{0,3,1}});
M1 = makeGWClass matrix(RR, {{1.0,24/10,-2.41},{24/10,-5,0},{-2.41,0,69}});
M2 = makeGWClass matrix(CC, {{1*ii,24/10,-2.41},{24/10,-5,0},{-2.41,0,69+ii}});
M3 = makeGWClass matrix(GF(7), {{1,0,0},{0,2,0},{0,0,-3}});

assert(getBaseField(M) === QQ);
assert(getBaseField(M1) === RR_53);
assert(getBaseField(M2) === CC_53);
assert((getBaseField M3).order == 7);
///

-- Test for gwAdd
-- Test 20
TEST ///
M1 = makeGWClass matrix(QQ, {{1/1,0,-3},{0,23,0},{-3,0,-2/5}});
M2 = makeGWClass matrix(QQ, {{0,1/2,0},{1/2,5/9,0},{0,0,1}});
M3 = makeGWClass matrix(QQ, {{1/1,0,-3,0,0,0},{0,23,0,0,0,0},{-3,0,-2/5,0,0,0},{0,0,0,0,1/2,0},{0,0,0,1/2,5/9,0},{0,0,0,0,0,1}})

G1 = makeGWClass matrix(RR, {{sqrt(2),0,-3},{0,sqrt(5),0},{-3,0,-1/5}});
G2 = makeGWClass matrix(RR, {{1/3}});
G3 = makeGWClass matrix(RR, {{sqrt(2),0,-3,0},{0,sqrt(5),0,0},{-3,0,-1/5,0},{0,0,0,1/3}});

H1 = makeGWClass matrix(CC, {{2*ii,0,0},{0,-2,0},{0,0,-3}});
H2 = makeGWClass matrix(CC, {{1,0,-3+ii,0},{0,-2,0,0},{-3+ii,0,-3,0},{0,0,0,5}});
H3 = makeGWClass matrix(CC, {{2*ii,0,0,0,0,0,0},{0,-2,0,0,0,0,0},{0,0,-3,0,0,0,0},{0,0,0,1,0,-3+ii,0},{0,0,0,0,-2,0,0},{0,0,0,-3+ii,0,-3,0},{0,0,0,0,0,0,5}});

assert(addGW(M1, M2) === M3);
assert(addGW(G1, G2) === G3);
assert(addGW(H1, H2) === H3);
///

-- Test for isIsotropic/isAnisotropic
-- Test 21
TEST ///
A1 = matrix(QQ, {{0,1/1},{1/1,0}});
assert(isIsotropic A1);
assert(not isAnisotropic makeGWClass(A1));

A2 = matrix(RR, {{1,-2,4},{-2,2,0},{4,0,-7}});
assert(not isAnisotropic A2);
assert(isIsotropic makeGWClass A2);

K=GF(13^4);
A3=matrix(K, {{7,81,63},{81,7,55},{63,55,109}});
assert(isIsotropic makeGWClass A3);
--Isotropic by the Chevalley-Warning Theorem.

A4 = matrix(QQ, {{5,0},{0,5}});
assert(isAnisotropic A4);

A5 = matrix(CC, {{3+ii,0},{0,5-ii}});
assert(not isAnisotropic A5);
///

--Tests for isIsomorphicForm
-- Test 22
TEST ///
B1 = matrix(QQ, {{1/1,-2/1,4/1},{-2/1,2/1,0},{4/1,0,-7/1}});
B2 = matrix(QQ, {{-17198/4225,-166126/975,-71771/1560},{-166126/975,-27758641/4050,-251077/135},{-71771/1560,-251077/135,-290407/576}});
assert(isIsomorphicForm(makeGWClass B1, makeGWClass B2));
B3 = matrix(QQ, {{-38/1,-50/1,23/1},{-50/1,-62/1,41/1},{23/1,41/1,29/1}});
assert(isIsomorphicForm(makeGWClass B1, makeGWClass B3));
///

--Test 23

TEST ///
D1 = matrix(QQ, {{1/1,-2/1,4/1},{-2/1,2/1,0},{4/1,0,-7/1}});
D2 = matrix(QQ, {{-38/1,-50/1,23/1},{-50/1,-62/1,41/1},{23/1,41/1,29/1}});
assert(isIsomorphicForm(makeGWClass D1, makeGWClass D2));

C1 = matrix(RR, {{1/1,-2/1,4/1},{-2/1,2/1,0},{4/1,0,-7/1}});
C2 = matrix(RR, {{-38/1,-50/1,23/1},{-50/1,-62/1,41/1},{23/1,41/1,29/1}});
assert(isIsomorphicForm(makeGWClass C1, makeGWClass C2));

M=GF(13^1)
C3 = matrix(M, {{1,11,4},{11,2,0},{4,0,6}});
C4 = matrix(M, {{1,2,10},{2,3,2},{10,2,3}});
assert(isIsomorphicForm(makeGWClass C3, makeGWClass C4));
///

-- Test for GWinvariants
-- Test 24
TEST ///
M1 = makeGWClass matrix(QQ, {{1/1,0,-3},{0,23,0},{-3,0,-2/5}});
M2 = makeGWClass matrix(QQ, {{1/1,0,0},{0, 23,0},{0,0,-2/5}});
M3 = makeGWClass matrix(QQ, {{1/1,0,0},{0,-23,0},{0,0,-2/5}});
M4 = makeGWClass matrix(QQ, {{-1/1,0,0},{0,-23,0},{0,0,-2/5}});

assert(getSignature(M1) == 1);
assert(getSignature(M2) == 1);
assert(getSignature(M3) == -1);
assert(getSignature(M4) == -3);

assert(getIntegralDiscriminant(M1) == -5405);
assert(getRelevantPrimes(M1) == {23, 5, 47});
assert(getHasseWittInvariant(M1, 5) == -1);
assert(getHasseWittInvariant(M1, 23) == 1);
assert(getHasseWittInvariant(M1, 47) == -1);
///

-- Test for HilbertSymbols
-- Test 25
TEST ///
assert(getHilbertSymbol(100, 7, 3) == 1);
assert(getHilbertSymbol(100/121, 7/169, 3) == 1);

assert(getHilbertSymbol(5, 1/9, 7) == 1);
assert(getHilbertSymbol(1/9, 5, 7) == 1);

assert(getHilbertSymbol(3, 11, 3) == -1);
assert(getHilbertSymbol(3, 11, 2) == -1);
assert(getHilbertSymbol(-3, -11, 2) == 1);
assert(getHilbertSymbol(-5, 11, 2) == -1);

assert(getHilbertSymbolReal(-3/1, 5) == 1);
assert(getHilbertSymbolReal(-3, -5/1) == -1);
assert(getHilbertSymbolReal(-3/1, -5) == -1);
assert(getHilbertSymbolReal(3, 5) == 1);
///
