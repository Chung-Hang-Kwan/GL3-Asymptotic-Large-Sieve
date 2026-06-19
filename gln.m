(*
    GL(n)pack: a Mathematica package to go with the book "Automorphic Forms and L
    Functions for the group GL(n,R)" by Dorian Goldfeld with an Appendix by Kevin
    A. Broughan, CUP 2006.





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                    !
    ! IF THIS TEXT APPEARS IN YOUR BROWSER WINDOW, USE   !
    ! THE FILE MENU BUTTON AT THE TOP OF YOUR BROWSER TO !
    ! SAVE THE DATA TO A FILE, WITH NAME gln.m, IN ANY   !
    ! SUBDIRECTORY ON YOUR DIRECTORY TREE IN WHICH YOU   !
    ! ARE ABLE TO SAVE FILES.                            !
    !                                                    !
    ! THEN READ THE FIRST FEW PARAGRAPHS OF THE GL(N)PACK!
    ! DOCUMENTATION, OR THE BOOK APPENDIX, TO SEE HOW TO !
    ! LOAD THE PACKAGE gln.m INTO MATHEMATICA.           !
    !                                                    !
    ! TO FIND OUT ABOUT DIRECTORY TREES CONSULT ONE OF   !
    ! THE TUTORIALS ON THE GL(N)PACK HOME PAGE OR TALK   !
    ! TO YOUR COMPUTER SUPPORT PERSON.                   !
    ! IF THERE ARE PROBLEMS EMAIL kab@waikato.ac.nz .    !
    !                                                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    



    

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should obtain a copy of the GNU General Public License
    along with this program from its website; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA


    Contact details: Kevin Broughan, 
                     Department of Mathematics,
                     University of Waikato, 
                     Hamilton Private Bag 3105,
                     New Zealand.
                     Email: kab@waikato.ac.nz

Package: Name GL(n)pack.
Author: Kevin Broughan
Copyright (c) Kevin Broughan, 2002-2015.
Version: 15th April 2015.
*)
Print[StyleForm["GL(n)pack version 15th April 2015 loading ....\
Copyright(c) 2002-2006 Kevin A. Broughan  \
GL(n)pack comes with ABSOLUTELY NO WARRANTY. This is free software, \
and you are welcome to redistribute it under the terms of the GNU Public Licence. \
For assistance consult kab@waikato.ac.nz .",FontWeight->"Bold"]];


BeginPackage["gln`"]
(********
Global Variables:
        $CartanEpsilon
Functions:
 Public functions:
        GlnVersion
	IwasawaForm (IWF private)
	CartanForm
	BruhatForm
        BruhatCVector
        LanglandsForm
        ParabolicQ
	GetCasimirOperator
        ApplyCasimirOperator

    UnimodularMatrix 
    ZeroMatrix
    MatrixRow, MatrixColumn
    ElementaryMatrix
    InsertMatrixElement
    RemoveMatrixRow, RemoveMatrixColumn
    SwapMatrixRows, SwapMatrixColumns
    BlockMatrix
    MatrixJoinHorizontal, MatrixJoinVertical
    LeadingMatrixBlock, TailingMatrixBlock
    MatrixDiagonal, DiagonalToMatrix
    MatrixUpperTriangle, MatrixLowerTriangle
    UpperTriangleToMatrix, LowerTriagleToMatrix
    WMatrix, WeylGenerator, LongElement
    ConstantMatrix, NRows, NColumns,  
    ModularGenerator, ModularGenerators
    IwasawaQ
    CreQ
   
    EisensteinFourierCoefficient
    EisensteinSeriesTerm
    HeckeEigenvalue
    HermiteFormLower
    HermiteFormUpper
    IFun
    LanglandsIFun
    IwasawaXMatrix
    IwasawaXVariables
    IwasawaYMatrix
    IwasawaYVariables
    RamanujanSum
    SmithForm
    SmithElementaryDivisors
    SmithInvariantFactors
    WeylGroup
    SpecialWeylGroup
  
    MakeXMatrix
    MakeXVariables
    MakeYMatrix
    MakeYVariables
    MakeZMatrix
    MakeZVariables 
    d, Wedge

    VolumeBall
    VolumeSphere
    VolumeHn
    VolumeFormUnimodular
    VolumeFormDiagonal
    VolumeFormHn
    VolumeFormGln

    WhittakerGamma
    Whittaker 
    StadeStar
    KloostermanSum
   
    FunctionalEquation
    KloostermanCompatibility
    KloostermanBruhatCell
    KloostermanSum
    KloostermanSum2
    PluckerCoordinates
    PluckerInverse
    PluckerRelations
    SubscriptedForm
    MPEisensteinLambdas 
    MPEisensteinSeries 
    MPEisensteinGamma 
    MPSymmetricPowerLFun 
    MPExteriorPowerLFun 
    MPSymmetricPowerGamma 
    MPExteriorPowerGamma *)

(*
 Private functions:
        KKSubsets
	SSquareMatrixQ
	NumericSquareMatrixQ
	IntegerTwoOrMoreQ
	TrivialMatrixQ
	UpperTriangular (UT)
	PrintUpperTriangular (PUT)
	UpperTriangularVariables (UTV)
	PrintIwasawaForm (PIF)
	Invert
	DiagonalToMatrix
	PrintCartan (PC)
	NormalizeMatrix
	NormalizeVector
	OneCasimir
	ExtractVariables
	RemoveSmallMatrixElements
	ZeroPad
        ArrayToMatrix (i), MatrixToArray (i) 
        AscendingLengthsQ (i),
        LowerTriangularQ (i) ListOfListsQ (i) MapJoin (i)
        StadeStarR
        SS
        PositiveBallQuad
        ComplexListQ
        IntegerListQ
        RealMatrixQ
        RealQ
        ZeroPadLeft
        cir15
        sph15
        PositiveBallQuad
        ExtractYVariables
        ExtractYVariablesFromIW
        ExtractXMatrixFromIW
        Ifun
        PBesselK
        UpperTriangularMatrixQ

  
      MatrixMinor 
      VectorWedge 
      PluckerRelations 
      PluckerRelationsTop 
      WedgeRelations 
      CramersRelation 
      SignatureSubset 
      PluckerClassic 
      InvertPlucker 
      TestPlucker4 - 7 
      getvector : of mniors for inclusion
       getmatrix : of minors for inclusion
       getminors : for the c_i and zero values
     getciminors : for the c_i only
    fixsignature
    lexibded : to set to 0
    removeinvalid : support
       getzeroelements : the lexicly bded 
       getlamis : indices for the c_i's
       imagemembership : needs revision
    gcdconst : get the sets of minors with gcd 1
       getoneminor : the bottom given col numbers 
    subit : replace integer indices with subscripts for cre's
    GenerallistofintegersQ
    MPEisensteinLambdasInternal
    MPEisensteinGammaInternal
    symmetricandexterior
    getsymmetricindices
    getexteriorindices
  *)


(*=== === === === === === === === === === === === === === === === === === ===*)\

(*------------------------------------------------------*)
Off[General::"spell",General::"spell1",Remove::"rmnsm",UpSet::"write"];

GlnVersion::usage="null -> True"
GlnVersion[]:=
  (Print[StyleForm["GL(n)pack: 15th April 2015",FontWeight->"Bold"]];
   Print[StyleForm[StringJoin["Mathematica: ",$Version],FontWeight->"Bold"]];True)
GlnVersion[argseq__]:= GlnVersion[]

IfPrint::usage="print it if $verbose==True"
$verbose::usage="True|False"

CreQ::usage = "expression -> True|False"

EisensteinFourierCoefficient::usage="expr,expr,natural no -> expr"
EisensteinSeriesTerm::usage="expr,expr,natural no -> expr"
HeckeEigenvalue::usage="natural no,natural no, string -> expr"
HermiteFormLower::usage="non-sing int matrix -> {unimodular matrix, integer matrix}"
HermiteFormUpper::usage="non-sing int matrix -> {unimodular matrix, integer matrix}"
IFun::usage="list of expr,matrix -> expr"
IwasawaXMatrix::usage="matrix -> matrix"
IwasawaXVariables::usage="matrix -> list"
IwasawaYMatrix::usage="matrix -> matrix"
IwasawaYVariables::usage="matrix -> list"


RamanujanSum::usage="natural no, natural no -> integer"
KloostermanSum::usage="list of natural no,list of natural no, \
list of integers, perm matrix->sum"
KloostermanSumQ::usage="list of natural no,list of natural no, \
list of integers, perm matrix->True|False"


SmithForm::usage="integer non-sing matrix -> \
    {unimod mat, diagonal matrix, unimodular mat}"
SmithInvariantFactors::usage="integer matrix -> list"
SmithElementaryDivisors::usage="integer matrix -> list"
WeylGroup::usage= "natural no -> list of matrices"
IwasawaQ::usage = "form -> True|False"
IwasawaForm::usage = "square matrix -> upper triangular matrix"

(*==================================================================================*)

CartanForm::usage = "numeric non-singular matrix -> {orthogonal, exp(symmetric)}"
$CartanEpsilon::usage = "Change to alter the upper bound for the size of matrix entries which are set to zero."

$CartanEpsilon=0.00001

BruhatForm::usage = "invertible -> {upper triangular, diag, perm, upper}"
BruhatCVector::usage = "invertible -> list"
GetCasimirOperator::usage = "natural number,natural number, string, string, string -> expr"

ApplyCasimirOperator::usage = "natural no, natural no, expr, list of symbols 
    -> expr|number"

UnimodularMatrix::usage="natural number, string -> square matrix." 
ZeroMatrix::usage="natural number, natural number -> matrix."
MatrixRow::usage="matrix -> list"
MatrixColumn::usage="matrix -> list"
ElementaryMatrix::usage="natural number, natural number,natural number, expr -> matrix"
InsertMatrixElement::usage="expr, natural number, natural number, matrix -> matrix"
GetMatrixElement::usage="matrix, natural number natural number -> expr"
RemoveMatrixRow::usage="matrix -> matrix"

    RemoveMatrixColumn::usage="matrix -> matrix"
    SwapMatrixRows::usage="matrix -> matrix"
    SwapMatrixColumns::usage="matrix -> matrix"
    BlockMatrix::usage="matrix,list,list -> matrix"
    MatrixJoinHorizontal::usage="matrix, matrix -> matrix"
    MatrixJoinVertical::usage="matrix,matrix -> matrix"
    LeadingMatrixBlock::usage="matrix -> matrix"
    TailingMatrixBlock::usage="matrix -> matrix"
    MatrixDiagonal::usage="matrix -> list"
    DiagonalToMatrix::usage="list -> matrix"
    MatrixUpperTriangle::usage="square matrix -> list of lists"
    MatrixLowerTriangle::usage="square matrix -> list of lists"
    UpperTriangleToMatrix::usage="list of lists -> matrix"
    LowerTriangleToMatrix::usage="list of lists -> matrix"
    WMatrix::usage="natural number -> matrix"
    WeylGenerator::usage="natural no, natural no, natural no -> matrix"
    ModularGenerators::usage="natural number -> {matrix, matrix}"
    LongElement::usage="natural number -> matrix"
    ConstantMatrix::usage="expr, natural number, natural number -> matrix"
    NRows::usage="matrix -> natural number"
    NColumns::usage="matrix -> natural number"  
  d::usage="expr -> expr"
    Wedge::usage="expr, expr, ... -> expr"
reWrite::usage="simplify differential forms expression" 

    VolumeBall::usage="expr, natural number -> expr"
    VolumeSphere::usage="expr, natural number -> expr"
    VolumeHn::usage="natural number -> expr"
    VolumeFormUnimodular::usage="string, natural number -> differential form"
    VolumeFormDiagonal::usage="string, natural number -> differential form"
    VolumeFormHn::usage="string, string, natural number -> differential form"
    VolumeFormGln::usage="string, natural number -> differential form"

    WhittakerGamma::usage="list of CREs -> expr"
    Whittaker::usage="CRE matrix, CRE list, CRE list,symbol -> expression"

    FunctionalEquation::usage="list of exprs, natural no -> list of exprs"

MakeXMatrix::usage = "natural no,string ->matrix"
MakeXVariables::usage = "natural no,string -> list"
MakeYMatrix::usage = "natural no,string ->matrix"
MakeYVariables::usage = "natural no,string -> list"
MakeZMatrix::usage = "natural no,string,string ->matrix"
MakeZVariables::usage = "natural no,string,string -> list"

MakeBlockMatrix::usage = "list of list of matrices of exprs -> matrix of exprs"
MakeMatrix::usage = "string, natural no, natural no -> matrix"

ParabolicQ::usage = "nonsingular matrix,list of positive integers -> True|False"
LanglandsForm::usage = 
 "upper block triangular nonsingular matrix, list of positive integers -> \
  {unipotent,central,block diag}"
LanglandsIFun::usage = "parabolic matrix, list of expr, list of positive integers -> expr"

SchurPolynomial::usage ="list of expr, list of expr -> expression"
HeckeOperator::usage = "natural no, square matrix, string -> expr"
SpecialWeylGroup::usage = "natural no -> list of square integer matrices"
HeckeMultiplicativeSplit::usage = "list of natural nos -> list of list of primes \
and list of powers."

HeckePowerSum::usage = 
    "non-neg integer,list of non-negative integers,string->expression.";


KloostermanCompatibility::usage = 
    "symbol,symbol,symbol,perm matrix,symbol->{{..},{{..},..},{}}";
KloostermanBruhatCell::usage = 
    "symbol,symbol,symbol, perm matrix,symbol -> rules  ";
KloostermanSum::usage = 
    "int list,int list, int list, perm matrix -> exponential sum ";
KloostermanSum2::usage = "int,int,int -> exponential sum ";
PluckerCoordinates::usage = " square CRE matrix -> {{..}, ..} ";
PluckerInverse::usage = " list of lists of integers -> square matrix ";
PluckerRelations::usage = " int >= 2, symbol -> {..,..,..} ";
SubscriptedForm::usage = "expression -> printed expression ";



MPEisensteinLambdas ::usage = "list of CREs -> list of CREs.";
MPEisensteinSeries ::usage = "CRE, list of CREs -> product of zetas.";
MPEisensteinGamma ::usage = "CRE, list of CREs ->  product of gammas.";
MPSymmetricPowerLFun ::usage = 
    "CRE, list of CRE's, nat. no. -> product of zetas.";
MPExteriorPowerLFun ::usage = 
    "CRE, list of CRE's, nat. no. -> product of zetas..";
MPSymmetricPowerGamma ::usage = 
    "CRE, list of CRE's, nat. no. -> product of gammas.";
MPExteriorPowerGamma ::usage = 
    "CRE, list of CRE's, nat. no. -> product of gammas.";




(*** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ********)
glv=GlnVersion;
aco=ApplyCasimirOperator;
gco=GetCasimirOperator;
bru=BruhatForm;
bcv=BruhatCVector;
blm=BlockMatrix;
com=ConstantMatrix;
car=CartanForm;
d2m=DiagonalToMatrix;
cas=CasimirOperator;
efc=EisensteinFourierCoefficient;
est=EisensteinSeriesTerm;
elm=ElementaryMatrix;
feq=FunctionalEquation;
gme=GetMatrixElement;
hms=HeckeMultiplicativeSplit;
hcs=HeckeCoefficientSum;
hps = HeckePowerSum;
hop=HeckeOperator;
hfl=HermiteFormLower;
hfu=HermiteFormUpper;
ifn=IFun;
ime=InsertMatrixElement;
iwf=IwasawaForm;
ixm=IwasawaXMatrix;
ixv=IwasawaXVariables;
iym=IwasawaYMatrix;
iyv=IwasawaYVariables;
lmb=LeadingMatrixBlock;
lel=LongElement;
ltm=LowerTriangleToMatrix;
mcl=MatrixColumn;
mdl=MatrixDiagonal;
mjh=MatrixJoinHorizontal;
mjv=MatrixJoinVertical;
mut=MatrixUpperTriangle;
mlt=MatrixLowerTriangle;
mro=MatrixRow;
mog=ModularGenerators;
nro=NRows;
ncl=NColumns;
rsm=RamanujanSum;
kls=KloostermanSum
klq=KloostermanSumQ
rmr=RemoveMatrixRow;
rmc=RemoveMatrixColumn;
smf=SmithForm;
sed=SmithElementaryDivisors;
sif=SmithInvariantFactors;
smr=SwapMatrixRows;
smc=SwapMatrixColumns;
tmb=TailingMatrixBlock;
utm=UpperTriangleToMatrix;
umm=UnimodularMatrix;
vfd=VolumeFormDiagonal;
vfg=VolumeFormGln;
vfh=VolumeFormHn;
vfu=VolumeFormUnimodular;
vbl=VolumeBall;
vhn=VolumeHn;
vsp=VolumeSphere;
weg=Wedge;
wge=WeylGenerator;
wgr=WeylGroup;
wit=Whittaker;
wig=WhittakerGamma;
wmx=WMatrix;
zmx=ZeroMatrix;

mxm=MakeXMatrix;
mxv=MakeXVariables;
mym=MakeYMatrix;
myv=MakeYVariables;
mzm=MakeZMatrix;
mzv=MakeZVariables;
iwq=IwasawaQ;

paq=ParabolicQ;
llf=LanglandsForm;
lif=LanglandsIFun;
hev=HeckeEigenvalue;
spl=SchurPolynomial;
swg=SpecialWeylGroup;

crq=CreQ;
mbm=MakeBlockMatrix;
mkm=MakeMatrix;

klc = KloostermanCompatibility;
kls = KloostermanSum;
cvm = CVectorToMatrix;
kbc = KloostermanBruhatCell;
plc = PluckerCoordinates;
pli = PluckerInverse;
plr = PluckerRelations;
suf = SubscriptedForm;

eil = MPEisensteinLambdas ;
eis = MPEisensteinSeries ;
eig = MPEisensteinGamma; 
spf = MPSymmetricPowerLFun ;
epl = MPExteriorPowerLFun; 
spg = MPSymmetricPowerGamma;
epg = MPExteriorPowerGamma ;



Begin["`Private`"]


IntegerTwoOrMoreQ[n_] := IntegerQ[n] && n > 1
NumericSquareMatrixQ[M_] := 
  SSquareMatrixQ[M] && NumberQ[Apply[Plus, Flatten[M]]]
SSquareMatrixQ[M_] := MatrixQ[M] && Length[M]==Length[First[M]] 
TrivialMatrixQ[M_]:=Length[M]==1

KKSubsets[lis_,n_]:= Subsets[lis,{n}]

UpperTriangular[n_, v_] := 
  Table[If[i == j, 1, If[i < j, v[i, j], 0]], {i, 1, n}, {j, 1, n}]

PrintUpperTriangular[n_, v_] := TraditionalForm[UpperTriangular[n, v]]

UT = UpperTriangular; PUT = PrintUpperTriangular;

UpperTriangularVariables[n_, v_] := 
  If[n == 2, v[1, 2], Apply[Union, Table[v[i, j], {i, 1, n}, {j, i + 1, n}]]]

UTV = UpperTriangularVariables

(*--------------------------------------------------------------------*)


CreQ::argx = "One argument required.";
CreQ[]:= (Message[CreQ::argx]; Abort[])

CreQ[argseq__] := (If[Length[{argseq}] != 1,(Message[CreQ::argx];Abort[])];
                   creq[First[{argseq}]])

creqold[e_]:=
   MemberQ[{Plus, Times, Power, Integer, Real, Complex, Rational, Symbol}, 
           Head[e]]

ExprQold[x_] := Not[Or[ListQ[x],StringQ[x]]]
ExprQ[x_]:= creq[x]

ExprListQ[x_] := And[ListQ[x], Not[x == {}], Not[MemberQ[ Map[ExprQ, x], False]]]
ExprSquareMatrixQ[x_] := 
   And[MatrixQ[x], Length[x] == Length[First[x]], ExprListQ[Apply[Join, x]]]

ExprMatrixQ[x_] := 
   And[MatrixQ[x], ExprListQ[Apply[Join, x]]]

(*---------------------------------------------------------*)
(*******************************************************************************)
IwasawaForm::argx = "IwasawaForm takes one argument."
IwasawaForm::argtype = "A square matrix of CREs is required."
IwasawaForm::singular = "A non-singular matrix of CREs is required."
IwasawaForm::trivial = "The matrix dimension for IwasawaForm should be 2 or more."

IwasawaForm[]:= (Message[IwasawaForm::argx]; Abort[])


IwasawaForm[argseq__] := Module[{args=List[argseq],M, Eqns, L, u, A, MM, a, n},

    If[Length[args] != 1, (Message[IwasawaForm::argx];Abort[])];
    g=First[args];
    If[Not[ExprSquareMatrixQ[g]], (Message[IwasawaForm::argtype];Abort[])];
    If[TrivialMatrixQ[g], (Message[IwasawaForm::trivial];Abort[])];
    If[Det[g] == 0, (Message[IwasawaForm::singular]; Abort[])];

    n = Length[g];
    M = Simplify[UT[n, \[Lambda]].g.Transpose[g]];
    Eqns = 
      Union[Flatten[ 
          Table[If[i < j, M[[i, j]] == 0, M[[1, 2]] == 0], {i, 1, n}, {j, 1, 
              n}]]];
    L = Flatten[Simplify[Solve[Eqns, UTV[n, \[Lambda]]]]];
    u = Simplify[UT[n, \[Lambda]] /. L];
    M = Simplify[M /. L];
    MM = Simplify[M.Transpose[u]];
    a = Table[If[i == j, Sqrt[MM[[i, j]]], 0], {i, 1, n}, {j, 1, n}];
    A = Simplify[Inverse[u].a];
    Return[Simplify[PowerExpand[Simplify[A/A[[n, n]]]]]] ]

IWF = IwasawaForm

PrintIwasawaForm[g_] := TraditionalForm[IwasawaForm[g]]

PIF = PrintIwasawaForm
(**********************************************************************)
(* nonsingular = upper . diag . permutation . upper where
      Tr(permutation). second upper . permutation is lower *)

BruhatForm::argx = "One argument required."
BruhatForm::argtype = "The argument to BruhatForm should be a \
square matrix of CREs"

BruhatForm::singular = "The argument is a non-singular matrix of CREs."
BruhatForm::trivial = "The matrix should have dimension 2 or more."

BruhatForm[]:= (Message[BruhatForm::argx]; Abort[])

BruhatForm[argseq__] := Module[{args=List[argseq],n, A, B, V, P, U, m, u1, u2, c, w, i, j},
    
    If[Length[args]!= 1, (Message[BruhatForm::argx];Abort[])];
    A=First[args]; n=Length[A];
    If[Not[ExprSquareMatrixQ[A]], (Message[BruhatForm::argtype];Abort[])];
    If[TrivialMatrixQ[A], (Message[BruhatForm::trivial];Abort[])];
    If[Det[A] == 0, (Message[BruhatForm::singular]; Abort[])];

    U = Table[If[i == j, 1, 0], {i, 1, n}, {j, 1, n}];
    P = Table[0, {i, 1, n}, {j, 1, n}];
    V = Table[0, {i, 1, n}, {j, 1, n}];
    B = Table[A[[i, j]], {i, 1, n}, {j, 1, n}];
    Do[{j = n, While[B[[j, i]] == 0, j = j - 1], P[[j, i]] = 1, 
        Do[V[[t, j]] = B[[t, i]], {t, 1, n}], 
        Do[{m = B[[j, k]]/B[[j, i]], U[[i, k]] = m, 
            Do[B[[l, k]] = B[[l, k]] - m  B[[l, i]], {l, 1, j - 1}], 
            B[[j, k]] = 0}, {k, i + 1, n}]}, {i, 1, n}];
    u2 = U; w = P;
    c = Table[0, {i, 1, n}, {j, 1, n}];
    u1 = Table[0, {i, 1, n}, {j, 1, n}];
    Do[c[[i, i]] = V[[i, i]], {i, 1, n}];
    Do[If[i > j, u1[[i, j]] = 0], {i, 1, n}, {j, 1, n}];
    Do[u1[[i, i]] = 1, {i, 1, n}];
    Do[If[i < j, u1[[i, j]] = V[[i, j]]/c[[j, j]]], {i, 1, n}, {j, 1, n}];
    Return[{u1, c, w, u2}]]

(*------------------------------------------------------------------------*)

bcv = BruhatCVector;
BruhatCVector::usage = "matrix,matrix->list"
BruhatCVector::argx = "The function takes two arguments."
BruhatCVector::arg1 = "The first argument is a square non-singular matrix."
BruhatCVector::arg2 = "The second argument is a permutation matrix."

BruhatCVector[]:= (Message[BruhatCVector::argx]; Abort[])

BruhatCVector[argseq__] := 
  Module[{args = List[argseq], n, tw, sigma, i, j, cv, ind,g,w},
    If[Length[args] != 2, (Message[BruhatCVector::argx]; Abort[])];
    g = First[args]; w = Part[args, 2];
    If[Not[And[SSquareMatrixQ[g], Det[g] != 0]], (Message[BruhatCVector::arg1];
         Abort[])];
    If[Not[PermutationMatrixQ[w]], (Message[BruhatCVector::arg2]; 
        Abort[])];

    n = Length[g];
    tw = Transpose[w];
    sigma = Table[0, {i, 1, n}];
    Do[{i = 1, While[tw[[i, j]] == 0, i = i + 1], sigma[[j]] = i}, {j, 1, 
        n}];
    cv = Table[0, {i, 1, n - 1}];
    Do[{z = Take[g, -i], ind = Table[sigma[[j]], {j, n - i + 1, n}], 
        m = Part[Transpose[z], ind], cv[[i]] = Abs[Det[m] ]}, {i, 1, n - 1}];
    Return[cv]]

PermutationMatrixQ[w_]:= And[SSquareMatrixQ[w],
                             Length[w]>1,
                             Abs[Det[w]]==1, 
                             Apply[Plus,Flatten[w]]==Length[w],
                             Or[Union[Flatten[w]]=={0,1},Union[Flatten[w]]=={1,0}]]

(*** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** *)
(* Cartan invertible -> {orthogonal, exp(symmetric)} *)

CartanForm::argx = "CartanForm takes one argument."
CartanForm::argtype = "The argument to CartanForm should be a numeric non-singular matrix."
CartanForm::singular = "The argument to CartanForm should be a non-singular  matrix."
CartanForm::trivial = "The matrix dimension should be greater than 1."

CartanForm[]:= (Message[CartanForm::argx]; Abort[])
CartanForm[argseq__] := Module[{args=List[argseq],n, k1, k2, k, d, l, L, D, S, s,g},

    If[Length[args]!= 1, (Message[CartanForm::argx];Abort[])];
    g=First[args]; n=Length[g];
    If[Not[NumericSquareMatrixQ[g]], (Message[CartanForm::argtype];Abort[])];
    If[TrivialMatrixQ[g], (Message[CartanForm::trivial];Abort[])];
    If[Det[g] == 0, (Message[CartanForm::singular]; Abort[])];

    n = Length[g];
    S = Transpose[g].g;
    l = Re[N[Eigenvalues[S]]];
    k1 = NormalizeMatrix[Re[N[Eigenvectors[S]]]];
    d = DiagonalToMatrix[Map[Invert, Map[Sqrt, l]]];
    L = DiagonalToMatrix[Map[Sqrt, l]];
    k2 = d.k1.Transpose[g]; 
    k = Transpose[k1.k2];
    s = k1.L.k1;
    k=RemoveSmallMatrixElements[k,$CartanEpsilon];
    s=RemoveSmallMatrixElements[s,$CartanEpsilon];
    Return[{k, s}]]


Invert[x_] := 1/x


DiagonalToMatrix::argx = "One argument required."
DiagonalToMatrix::arg1 = "The argument should be a list."



DiagonalToMatrix[argseq__] := Module[{args=List[argseq],d,n, M},
 If[Length[args]!= 1, (Message[DiagonalToMatrix::argx];Abort[])];
    d=First[args];
    If[Not[ExprListQ[d]], 
       (Message[DiagonalToMatrix::arg1];Abort[])];
    n = Length[d];
    M = Table[If[i == j, Part[d, i], 0], {i, 1, n}, {j, 1, n}]; Return[M]]

PrintCartan[g_] := Module[{ans, k, S},
    ans = N[Cartan[g]];
    Print[MatrixForm[First[ans]]];
    Print[MatrixForm[Part[ans, 2]]];
    Return[True]]

PC = PrintCartan

NormalizeMatrix[m_] := Map[NormalizeVector, m]

NormalizeVector[v_] := v/Sqrt[v.v]

RemoveSmallMatrixElements[M_,eps_]:=Module[{Mout=M, n=Length[M], m=Length[First[M]]},
  Do[If[Abs[M[[i,j]]] < eps, Mout[[i,j]]=0], {i,1,n}, {j,1,m}];
  Return[Mout]]

(*** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **)


OneCasimir[n_, is_, x_, y_,f_] := (* is = {i1, ...im} *)
  Module[{X, Y, e = 1, M, m = Length[is], t, Mij, Id, xvars, yvars, i, term, A},
    X = UpperTriangular[n, x];
    (* y[n] = 1;*)
    Y = Table[If[i == j, 1, 0], {i, 1, n}, {j, 1, n}];
    Do[{e = e*y[i], Y[[n - i, n - i]] = e}, {i, 1, n - 1}]; 
    Id = Table[If[i == j, 1, 0], {i, 1, n}, {j, 1, n}];
    Mij = Table[0, {i, 1, n}, {j, 1, n}];
    M = Id;
    Do[{Mij = Id, 
        Mij[[Part[is, i], Part[is, i + 1]]] = 
          t[i] + Mij[[Part[is, i], Part[is, i + 1]]], M = M.Mij}, {i, 1, 
        m - 1}];
    Mij = Id;
    Mij[[Part[is, m], Part[is, 1]]] = Mij[[Part[is, m], Part[is, 1]]] + t[m];
    M = M.Mij;
    A = IWF[X.Y.M];
    xvars=IwasawaXVariables[A];
    yvars=IwasawaYVariables[A];
    fun = Apply[f, Join[xvars,yvars]];
    Do[ fun = D[fun, t[i]], {i, 1, m}];
    term=fun;
    Do[term=term /.t[i]->0, {i,1,m}];
    Return[term]]

ExtractVariables[A_] := Module[{out = {}, n = Length[A]},
    Do[Do[out = Join[{A[[i, j]]}, out], {j, i, n}], {i, 1, n - 1}];
    Return[Reverse[out]]]

GetCasimirOperator::argx = "The function takes 4 arguments."
GetCasimirOperator::arg1 = "The first argument should be an integer 2 or more."
GetCasimirOperator::arg2 = "The second argument should be an integer 2 or more."
GetCasimirOperator::arg3 = "The 3rd/4th/5th arguments should each be a string"

GetCasimirOperator[]:= (Message[GetCasimirOperator::argx]; Abort[])

GetCasimirOperator[argseq__] := 
    Module[{args=List[argseq],op = 0,m,n,x,y,arg3,arg4,arg5,fun},

    If[Length[args] != 5, (Message[GetCasimirOperator::argx];Abort[])];
    m=Part[args,1]; n=Part[args,2]; arg3=Part[args,3]; arg4=Part[args,4]; arg5=Part[args,5];
    If[Not[IntegerTwoOrMoreQ[m]], (Message[GetCasimirOperator::arg1];Abort[])];
    If[Not[IntegerTwoOrMoreQ[n]], (Message[GetCasimirOperator::arg2];Abort[])];
    If[Or[Not[StringQ[arg3]],Not[StringQ[arg4]],Not[StringQ[arg5]]], 
         (Message[GetCasimirOperator::arg3];Abort[])];

    x=Symbol[arg3]; y=Symbol[arg4]; fun=Symbol[arg5];

    Do[{term = OneCasimir[n, ZeroPad[Reverse[IntegerDigits[i, n]], m] + 1, x, y,fun], 
      op = Expand[op + term]}, {i, 0, n^m - 1}];
    Return[Simplify[op, Assumptions->Table[y[i]>0, {i,1,n-1}]]]]

ZeroPad[x_, m_] := 
  If[Length[x] >= m, x, Join[x, Table[0, {i, 1, m - Length[x]}]]]

(*==========================================================================*)


OneApplyCasimir[n_, is_, iwf_, expr_] := (* is = {i1, ...im} *)
  Module[{X, Y, e = 1, M, m = Length[is], t, Mij, Id, vblesx, vblesy,vbles, newvars, 
          newxvars, newyvars, newexpr, i, term, A},
    
   (* X = Table[If[i==j,1,0], {i,1,n},{j,1,n}];
    Do[If[i<j, X[[i,j]]=Part[vbles,n(i-1)-i(i-1)/2+j-i]],{i,1,n-1},{j,2,n}];
   (*  y[n] = 1; *)
    Y = Table[If[i == j, 1, 0], {i, 1, n}, {j, 1, n}];

    Do[{e = e*vbles[[(n^2-n)/2+i]], Y[[n - i, n - i]] = e}, {i, 1, n - 1}]; *)

    Id = Table[If[i == j, 1, 0], {i, 1, n}, {j, 1, n}];
    Mij = Table[0, {i, 1, n}, {j, 1, n}];
    M = Id;
    Do[{Mij = Id, 
        Mij[[Part[is, i], Part[is, i + 1]]] = 
          t[i] + Mij[[Part[is, i], Part[is, i + 1]]], M = M.Mij}, {i, 1, 
        m - 1}];
    Mij = Id;
    Mij[[Part[is, m], Part[is, 1]]] = Mij[[Part[is, m], Part[is, 1]]] + t[m];
    M = M.Mij;

    A = IWF[iwf.M];
    newxvars=IwasawaXVariables[A];
    newyvars=IwasawaYVariables[A];
    newvars=Join[newxvars,newyvars];

    vblesx=IwasawaXVariables[iwf];
    vblesy=IwasawaYVariables[iwf];
    vbles=Join[vblesx,vblesy];

    newexpr=expr;
    Do[newexpr=newexpr/. Part[vbles,i] -> Part[newvars,i],{i,1,Length[vbles]}];

    Do[ newexpr = D[newexpr, t[i]], {i, 1, m}];
    term=newexpr;
    Do[term=term /.t[i]->0, {i,1,m}];
    
    Return[term]]

(*------------------------------------------------------------------------------*)

ApplyCasimirOperator::argx = "The function takes 3 arguments."
ApplyCasimirOperator::arg1 = "The first argument should be an integer 2 or more."
ApplyCasimirOperator::arg2 = "The 2nd argument should be a function or expression in 
a matrix in Iwasawa form."
ApplyCasimirOperator::arg3 = "The 3rd argument should be a matrix of CREs in Iwasawa form"

ApplyCasimirOperator[]:=(Message[ApplyCasimirOperator::argx]; Abort[])

ApplyCasimirOperator[argseq__] := 
    Module[{args=List[argseq],op = 0,m,n,x,y,arg3,expr,vbles,iwf},

    If[Length[args] != 3, (Message[ApplyCasimirOperator::argx];Abort[])];
    m=Part[args,1]; expr=Part[args,2]; iwf=Part[args,3];
    If[Not[IntegerTwoOrMoreQ[m]], (Message[ApplyCasimirOperator::arg1];Abort[])];
   
    If[Not[IwasawaQ[iwf]], 
         (Message[ApplyCasimirOperator::arg3];Abort[])];
    n=Length[iwf]; 
    Do[{term = 
          OneApplyCasimir[n, ZeroPad[Reverse[IntegerDigits[i, n]], m] + 1, iwf, expr], 
      op = Expand[op + term]}, {i, 0, n^m - 1}];
  (*  Return[Simplify[op, Assumptions->Table[y[i]>0, {i,1,n-1}]]]] *)
      Return[op]]

(* ---1-------------------------------------------- *)




NRows::arg1="Argument is a matrix of CREs";
NRows::argx="One argument is required";

NRows[]:= (Message[NRows::argx]; Abort[])

NRows[argseq__]:=Module[{ars=List[argseq],A},
    If[Length[ars]\[NotEqual] 1, (Message[NRows::argx];Abort[])];
    A=First[ars];
    If[Not[ExprMatrixQ[A]],(Message[NRows::arg1];Abort[])];
    Return[Length[A]]]

(* -----2------------------------------------------- *)

NColumns::arg1="Argument is a matrix of CREs";
NColumns::argx="One argument is required";

NColumns[]:= (Message[NColumns::argx]; Abort[])

NColumns[argseq__]:=Module[{ars=List[argseq],A},
    If[Length[ars]\[NotEqual] 1, (Message[NColumns::argx];Abort[])];
    A=First[ars];
    If[Not[ExprMatrixQ[A]],(Message[NColumns::arg1];Abort[])];
    Return[Length[First[A]]]]

(* ----3------------------------------------------------ *)

UnimodularMatrix::argx = "Two arguments required";
UnimodularMatrix::arg1 = "First argument us an integer >= 1";
UnimodularMatrix::arg2="Second argument is a string";

UnimodularMatrix[]:= (Message[UnimodularMatrix::argx]; Abort[])

UnimodularMatrix[argseq__]:=Module[{args=List[argseq],n,ustr,u,i,j,uarray},
    If[Length[args]\[NotEqual]2,(Message[UnimodularMatrix::argx];
        Abort[])];
    n=First[args];ustr=Part[args,2];
    If[!And[IntegerQ[n],n\[GreaterEqual]1],(Message[UnimodularMatrix::arg1];
        Abort[])];
    If[!StringQ[ustr],(Message[UnimodularMatrix::arg2];Abort[])];
    u=Symbol[ustr];
    Do[uarray[i,i]=1, {i,1,n}];
    Do[uarray[i,j]=u[i,j] , {i,1,n-1}, {j,i+1,n}];
    Do[uarray[i,j]=0, {i,2,n},{j,1,i-1}];
    Return[ArrayToMatrix[n,n,uarray]]]

(* ---4------------------------------------------------- *)

LongElement::argx="One argument required";
LongElement::arg1="The argument is an integer >= 1";

LongElement[]:= (Message[LongElement::argx]; Abort[])

LongElement[argseq__]:=Module[{args=List[argseq],n,i,j,marray},
    If[Length[args]\[NotEqual]1,(Message[LongElement::argx];
        Abort[])];
    n=First[args];
    If[!And[IntegerQ[n],n\[GreaterEqual]1],(Message[LongElement::arg1];
        Abort[])];Do[marray[i,j]=0, {i,1,n}, {j,1,n}];
    Do[marray[i,n-i+1]=1, {i,1,n}];
    Return[ArrayToMatrix[n,n,marray]]]

(* ---5---------------------------------------------------- *)

WeylGenerator::argx="Three arguments required";
WeylGenerator::arg1="First argument is an integer >= 2.";
WeylGenerator::arg2="Second argument is an integer <= first.";
WeylGenerator::arg3="Third argument is an integer not equal to the second.";

WeylGenerator[]:= (Message[WeylGenerator::argx]; Abort[])

WeylGenerator[argseq__]:=Module[{args=List[argseq],n,i,j,u,v,marray},
    If[Length[args]\[NotEqual]3,(Message[WeylGenerator::argx];
        Abort[])];
    n=Part[args,1];i=Part[args,2];j=Part[args,3];
    If[!And[IntegerQ[n],n\[GreaterEqual]2],(Message[WeylGenerator::arg1];
        Abort[])];
    If[!And[IntegerQ[i],i\[LessEqual] n,1\[LessEqual] i],(Message[
          WeylGenerator::arg2];Abort[])];
    If[!And[IntegerQ[j],j\[GreaterEqual]1,j\[LessEqual] n,
          j \[NotEqual] i],(Message[WeylGenerator::arg3];Abort[])];
    
    Do[marray[u,v]=0, {u,1,n}, {v,1,n}];
    Do[marray[u,u]=1, {u,1,n}];
    marray[i,j]=-1;marray[j,i]=1;
    marray[i,i]=0; marray[j,j]=0;
    Return[ArrayToMatrix[n,n,marray]]]

(* ---6-------------------------------------------------- *)
ModularGenerators::argx="One argument required";
ModularGenerators::arg1="Argument is an integer >= 2.";

ModularGenerators[]:= (Message[ModularGenerators::argx]; Abort[])

ModularGenerators[argseq__]:=Module[{args=List[argseq],n,i,j,marray,u1,u2},
    If[Length[args]\[NotEqual]1,(Message[ModularGenerators::argx];
        Abort[])];
    n=Part[args,1];
    If[!And[IntegerQ[n],n\[GreaterEqual]2],(Message[ModularGenerators::arg1];
        Abort[])];
    
    Do[marray[i,j]=0, {i,1,n}, {j,1,n}];
    Do[marray[i+1,i]=1, {i,1,n-1}];
    marray[1,n]=(-1)^(n-1);
    u1=ArrayToMatrix[n,n,marray];
    
    Do[marray[i,j]=0, {i,1,n}, {j,1,n}];
    Do[marray[i,i]=1, {i,1,n}];
    marray[1,2]=1;
    u2=ArrayToMatrix[n,n,marray];
    
    Return[{u1,u2}]]

(*---------------------------------------------------------------*)

WMatrix::argx="One argument required";
WMatrix::arg1="Argument is an integer >= 2.";

WMatrix[]:= (Message[WMatrix::argx]; Abort[])

WMatrix[argseq__]:=Module[{args=List[argseq],n,amat, i,j},
    If[Length[args]\[NotEqual]1,(Message[WMatrix::argx];Abort[])];
    n=First[args];
    If[!And[IntegerQ[n],n\[GreaterEqual]2],(Message[WMatrix::arg1];
        Abort[])];Do[amat[i,j]=0, {i,1,n}, {j,1,n}];
    Do[amat[n-i+1,i]=1, {i,1,n-1}];
    amat[1,n]=(-1)^IntegerPart[n/2];
    Return[ArrayToMatrix[n,n,amat]]]
(* ----7---------------------------------------------------  *)
ZeroMatrix::argx="Two arguments required.";
ZeroMatrix::arg1="First Argument is an integer >= 1.";
ZeroMatrix::arg2="Second Argument is an integer >= 1.";

zeroMatrix[]:= (Message[ZeroMatrix::argx]; Abort[])

ZeroMatrix[argseq__]:= Module[{args=List[argseq],m,n,c},
    If[Length[args]\[NotEqual]2,(Message[ZeroMatrix::argx];Abort[])];
    m=Part[args,1];n=Part[args,2];
    If[!And[IntegerQ[m],m\[GreaterEqual]1],(Message[ZeroMatrix::arg1];
        Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual]1],(Message[ZeroMatrix::arg2];
        Abort[])];Return[Table[0,{i,1,m},{j,1,n}]]]

(* ----8-------------------------------------------------- *)
MatrixRow::argx="Two arguments required.";
MatrixRow::arg1="First argument is a matrix of CREs.";
MatrixRow::arg2="Second argument is a valid row index.";

MatrixRow[]:= (Message[MatrixRow::argx]; Abort[])

MatrixRow[argseq__]:= Module[{args=List[argseq],A,i},
    If[Length[args]\[NotEqual]2,(Message[MatrixRow::argx];Abort[])];
    A=Part[args,1];i=Part[args,2];
    If[!ExprMatrixQ[A],(Message[MatrixRow::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NRows[A]],(Message[
          MatrixRow::arg2];Abort[])];
    Return[Part[A,i]]]

(* ---9-------------------------------------------------- *)
MatrixColumn::argx="One argument required";
MatrixColumn::arg1="First argument is a matrix of CREs.";
MatrixColumn::arg2="Second argument is a valid column index.";

MatrixColumn[]:= (Message[MatrixColumn::argx]; Abort[])

MatrixColumn[argseq__]:=Module[{args=List[argseq],A,j,partj},
    If[Length[args]\[NotEqual]2,(Message[MatrixColumn::argx];
        Abort[])];
    A=Part[args,1]; j=Part[args,2];
    If[!ExprMatrixQ[A],(Message[MatrixColumn::arg1];Abort[])];
    If[!And[IntegerQ[j],j\[GreaterEqual]1,j\[LessEqual] NColumns[A]],(Message[
          MatrixColumn::arg2];Abort[])];partj[x_]:=Part[x,j];
    Return[Map[partj,A]]]

(* ---10-------------------------------------------------- *)
ElementaryMatrix::argx="Four arguments required.";
ElementaryMatrix::arg1="First argument is an integer >= 1.";
ElementaryMatrix::arg2="Second argument is a row index.";
ElementaryMatrix::arg3="Third argument is a column index.";
ElementaryMatrix::arg4="Fourth argument is a CRE.";

ElementaryMatrix[]:= (Message[ElementaryMatrix::argx]; Abort[])

ElementaryMatrix[argseq__]:= Module[{args=List[argseq],n,i,j,c},
    If[Length[args]\[NotEqual]4,(Message[ElementaryMatrix::argx];
        Abort[])];
    n=Part[args,1];i=Part[args,2]; j=Part[args,3];
    c=Part[args,4];
    If[!And[IntegerQ[n],n\[GreaterEqual] 1],(Message[ElementaryMatrix::arg1];
        Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] n],(Message[
          ElementaryMatrix::arg2];Abort[])];
    If[!And[IntegerQ[j],1\[LessEqual] j, j\[LessEqual] n],(Message[
          ElementaryMatrix::arg3];Abort[])];
    If[Not[ExprQ[c]],(Message[
          ElementaryMatrix::arg4];Abort[])];
    Return[InsertMatrixElement[c,i,j,IdentityMatrix[n]]]]

(* ----11----------------------------------------------------- *)
GetMatrixElement::argx="Three arguments required.";
GetMatrixElement::arg1="First argument is a matrix of CREs.";
GetMatrixElement::arg2="Second argment is a row index.";
GetMatrixElement::arg3="Third argument is a column index.";

GetMatrixElement[]:= (Message[GetMatrixElement::argx]; Abort[])

GetMatrixElement[argseq__]:= Module[{args=List[argseq],A,i,j},
    If[Length[args]\[NotEqual]3,(Message[GetMatrixElement::argx];
        Abort[])];
    A=Part[args,1];i=Part[args,2]; j=Part[args,3];
    If[!ExprMatrixQ[A],(Message[GetMatrixElement::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NRows[A]],(Message[
          GetMatrixElement::arg2];Abort[])];
    If[!And[IntegerQ[j],1\[LessEqual] j, j\[LessEqual] NColumns[A]],(Message[
          GetMatrixElement::arg3];Abort[])];
    Return[ Part[Part[A,i],j]]]

(* ---12------------------------------------------------------ *)
InsertMatrixElement::argx = "Four arguments required";
InsertMatrixElement::arg1 = "First argument is a CRE.";
InsertMatrixElement::arg2="Second argment is a row index.";
InsertMatrixElement::arg3="Third argument is a column index.";
InsertMatrixElement::arg4="Fourth argument is a matrix of CREs.";

InsertMatrixElement[]:= (Message[InsertMatrixElement::argx]; Abort[])

InsertMatrixElement[argseq__]:=Module[{args=List[argseq],e,i,j,A,marray,m,n},
    If[Length[args]\[NotEqual]4,(Message[InsertMatrixElement::argx];
        Abort[])];
    e=Part[args,1];i=Part[args,2]; j=Part[args,3]; A=Part[args,4];
    If[Not[ExprQ[e]], (Message[InsertMatrixElement::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NRows[A]],(Message[
          InsertMatrixElement::arg2];Abort[])];
    If[!And[IntegerQ[j],1\[LessEqual] j, j\[LessEqual] NColumns[A]],(Message[
          InsertMatrixElement::arg3];Abort[])];
    If[!ExprMatrixQ[A],(Message[InsertMatrixElement::arg4];Abort[])];

    MatrixToArray[A,marray];
    marray[i,j]=e;
    Return[ArrayToMatrix[NRows[A],NColumns[A],marray]]]

(* ----13---------------------------------------------------- *)
RemoveMatrixRow::argx="Two arguments required";
RemoveMatrixRow::arg1="First argument is a matrix of CREs.";
RemoveMatrixRow::arg2="Second argument is a row index.";

RemoveMatrixRow[]:= (Message[RemoveMatrixRow::argx]; Abort[])

RemoveMatrixRow[argseq__]:=
  Module[{args=List[argseq],A, i,above,below, m},
    If[Length[args]\[NotEqual]2,(Message[RemoveMatrixRow::argx];
        Abort[])];
    A=Part[args,1];i=Part[args,2]; 
    If[!ExprMatrixQ[A],(Message[RemoveMatrixRow::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NRows[A]],(Message[
          RemoveMatrixRow::arg2];Abort[])];
    m=Length[A];
    Return[Join[Take[A,i-1],Take[A,i-m]]]]

(* --14----------------------------------------------------- *)

RemoveMatrixColumn::argx="Two arguments required";
RemoveMatrixColumn::arg1="First argument is a matrix of CREs.";
RemoveMatrixColumn::arg2="Second argument is a column index.";

RemoveMatrixColumn[]:= (Message[RemoveMatrixColumn::argx]; Abort[])

RemoveMatrixColumn[argseq__] := Module[{args = List[argseq], A, j, 
  TakeLeft, TakeRight, n, ans = {}}, 
   If[Length[args] 
      \[NotEqual] 2, (Message[RemoveMatrixColumn::argx];
        Abort[])];
    A = Part[args, 1]; j = Part[args, 2]; TakeLeft[x_] := Take[x, j - 1];
    If[Not[ ExprMatrixQ[A]], (Message[RemoveMatrixRow::arg1]; Abort[])];
    If[Not[ And[IntegerQ[j], GreaterEqual[j,1], GreaterEqual[NColumns[A],j]]], 
      (Message[RemoveMatrixColumn::arg2]; Abort[])];
    n=Length[A];
    TakeRight[x_] := Take[x, j - n];
    Do[ans = Join[{Join[TakeLeft[Part[A, i]], TakeRight[Part[A, i]]]}, 
    ans], {i, 1, Length[A]}];
    Return[Reverse[ans]]]
(* ---15-------------------------------------------------- *)

SwapMatrixRows::argx="Three arguments required";
SwapMatrixRows::arg1="First argument is a matrix of CREs.";
SwapMatrixRows::arg2="Second argument is a row index.";
SwapMatrixRows::arg3="Third argument is a distinct row index.";

SwapMatrixRows[]:= (Message[SwapMatrixRows::argx]; Abort[])

SwapMatrixRows[argseq__]:= (* i \[NotEqual] j *) 
  Module[{args=List[argseq],A,i,j,rowi,rowj,above,u,below,between, 
      m},
    If[Length[args]\[NotEqual]3,(Message[SwapMatrixRows::argx];
        Abort[])];
    A=Part[args,1]; i=Part[args,2]; j=Part[args,3];
    If[!ExprMatrixQ[A],(Message[SwapMatrixRows::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NRows[A]],(Message[
          SwapMatrixRows::arg2];Abort[])];
    If[!And[IntegerQ[j],j\[GreaterEqual]1,j\[LessEqual] NRows[A], 
          i\[NotEqual] j],(Message[SwapMatrixRows::arg3];Abort[])];
    If[j<i,{u=j,j=i,i=u}];
    rowi=Part[A,i]; rowj=Part[A,j];
    above=Take[A,i-1]; between=Take[A,{i+1,j-1}];
    m=Length[A];
    below=Take[A,m-j];
    Return[Join[above, {rowj},between, {rowi}, below]]]

(* --16---------------------------------------------------- *)

SwapMatrixColumns::argx="Three arguments required.";
SwapMatrixColumns::arg1="First argument is a matrix of CREs.";
SwapMatrixColumns::arg2="Second argument is a column index.";
SwapMatrixColumns::arg3="Third argument is a distinct column index.";

SwapMatrixColumns[]:= (Message[SwapMatrixColumns::argx]; Abort[])

SwapMatrixColumns[argseq__]:= (* i \[NotEqual] j *) 
  Module[{args=List[argseq],A,i,j, m=Length[A],u},
    If[Length[args]\[NotEqual]3,(Message[SwapMatrixColumns::argx];
        Abort[])];
    A=Part[args,1]; i=Part[args,2]; j=Part[args,3];If[j<i,{u=j,j=i,i=u}];
    If[!ExprMatrixQ[A],(Message[SwapMatrixColumns::arg1];Abort[])];
    If[!And[IntegerQ[i],i\[GreaterEqual]1,i\[LessEqual] NColumns[A]],(Message[
          SwapMatrixColumns::arg2];Abort[])];
    If[!And[IntegerQ[j],j\[GreaterEqual]1,j\[LessEqual] NColumns[A], 
          i\[NotEqual] j],(Message[SwapMatrixColumns::arg3];Abort[])];
    DoSwapEls[x_]:=SwapListElements[x,i,j];
    Return[Map[DoSwapEls,A]]]

(* --17---------------------------------------------------- *)

BlockMatrix::argx="Three arguments required.";
BlockMatrix::arg1="First argument is a matrix of CREs.";
BlockMatrix::arg2=
    "Second argument is a list of two row indices in ascending order.";
BlockMatrix::arg3=
    "Third argument is a list of two column indices in ascending order.";

BlockMatrix[]:= (Message[BlockMatrix::argx]; Abort[])

BlockMatrix[argseq__]:=
  Module[{args=List[argseq],A,rows,cols,subrows,n,m,i1,i2,j1,j2},
    If[Length[args]\[NotEqual]3,(Message[BlockMatrix::argx];
        Abort[])];
    A=Part[args,1]; rows=Part[args,2]; cols=Part[args,3];
    
    If[!ExprMatrixQ[A],(Message[BlockMatrix::arg1];Abort[])];
    
    If[!And[ListQ[rows],Length[rows]\[Equal]2,IntegerQ[First[rows]],
          IntegerQ[Part[rows,2]]],(Message[BlockMatrix::arg2];
        Abort[])];
    
    If[!And[ListQ[cols],Length[cols]\[Equal]2,IntegerQ[First[cols]],
          IntegerQ[Part[cols,2]]],(Message[BlockMatrix::arg3];
        Abort[])];
    
    m=Length[A];n=Length[First[A]];i1=First[rows];i2=Last[rows];
    j1=First[cols];j2=Last[cols];
    
    If[!And[IntegerQ[i1],i1\[GreaterEqual]1,i1\[LessEqual] NRows[A],
          IntegerQ[i2],i2\[GreaterEqual]1,i2\[LessEqual] NRows[A],
          i1\[LessEqual] i2],(Message[BlockMatrix::arg2];Abort[])];
    
    If[!And[IntegerQ[j1],j1\[GreaterEqual]1,j1\[LessEqual] NColumns[A],
          IntegerQ[j2],j2\[GreaterEqual]1,j2\[LessEqual] NColumns[A],
          j1\[LessEqual] j2],(Message[BlockMatrix::arg3];Abort[])];
    
    TakeSubElts[x_]:=Take[x,{j1,j2}];
    Return[Map[TakeSubElts,Take[A,{i1,i2}]]] ]

(* --18---------------------------------------------------- *)

ConstantMatrix::argx="Three arguments required.";
ConstantMatrix::arg2="Second argument is an integer >= 1.";
ConstantMatrix::arg3="Third argument is an integer >= 1.";
ConstantMatrix::arg1="First argument is a CRE.";

ConstantMatrix[]:= (Message[ConstantMatrix::argx]; Abort[])

ConstantMatrix[argseq__]:= Module[{args=List[argseq],m,n,c,i,j},
    If[Length[args]\[NotEqual]3,(Message[ConstantMatrix::argx];
        Abort[])];
    m=Part[args,2];n=Part[args,3]; c=Part[args,1];
    If[Not[ExprQ[c]],(Message[ConstantMatrix::arg1];
        Abort[])];
    If[!And[IntegerQ[m],m\[GreaterEqual]1],(Message[ConstantMatrix::arg2];
        Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual]1],(Message[ConstantMatrix::arg3];
        Abort[])];
    Return[Table[c,{i,1,m},{j,1,n}]]]

(* ---19---------------------------------------------------- *)

MatrixJoinHorizontal::argx="Two arguments requred.";
MatrixJoinHorizontal::arg1="First argument is a matrix of CREs.";
MatrixJoinHorizontal::arg2=
    "Second argment is matrix of CREs with the same number of rows as the first.";

MatrixJoinHorizontal[]:= (Message[MatrixJoinHorizontal::argx]; Abort[])

MatrixJoinHorizontal[argseq__]:= Module[{args=List[argseq],A,B},
    If[Length[args]\[NotEqual]2,(Message[MatrixJoinHorizontal::argx];
        Abort[])];
    A=Part[args,1];B=Part[args,2];
    If[!ExprMatrixQ[A],(Message[MatrixJoinHorizontal::arg1];Abort[])];
    If[!And[ExprMatrixQ[B],NRows[A]\[Equal]NRows[B]],(Message[
          MatrixJoinHorizontal::arg1];Abort[])];
    Return[MapJoin[A,B]]]

(* --20----------------------------------------------------- *)

MatrixJoinVertical::argx="Two arguments requred.";
MatrixJoinVertical::arg1="First argument is a matrix of CREs.";
MatrixJoinVertical::arg2=
    "Second argment is matrix of CREs with the same number of columns as the first.";

MatrixJoinVertical[]:= (Message[MatrixJoinVertical::argx]; Abort[])

MatrixJoinVertical[argseq__]:= Module[{args=List[argseq],A,B},
    If[Length[args]\[NotEqual]2,(Message[MatrixJoinVertical::argx];
        Abort[])];
    A=Part[args,1];B=Part[args,2];
    If[!ExprMatrixQ[A],(Message[MatrixJoinVertical::arg1];Abort[])];
    If[!And[ExprMatrixQ[B],NColumns[A]\[Equal]NColumns[B]],(Message[
          MatrixJoinVertical::arg1];Abort[])];
    Return[
      Join[A,B]]]

(* ---21---------------------------------------------------- *)

LeadingMatrixBlock::argx="Three arguments required.";
LeadingMatrixBlock::arg1="First argument is a matrix of CREs.";
LeadingMatrixBlock::arg2="Second argment is a valid row index.";
LeadingMatrixBlock::arg3="Third argument is a valid column index.";

LeadingMatrixBlock[]:= (Message[LeadingMatrixBlock::argx]; Abort[])

LeadingMatrixBlock[argseq__]:=
  Module[{args=List[argseq],A,a,b,matarray, matblock,i,j},
    If[Length[args]\[NotEqual]3,(Message[LeadingMatrixBlock::argx];
        Abort[])];
    A=Part[args,1]; a=Part[args,2]; b=Part[args,3];
    
    If[!ExprMatrixQ[A],(Message[LeadingMatrixBlock::arg1];Abort[])];
    If[!And[IntegerQ[a],a\[GreaterEqual]1,a\[LessEqual] NRows[A]],(Message[
          LeadingMatrixBlock::arg2];Abort[])];
    If[!And[IntegerQ[b],1\[LessEqual] b, b\[LessEqual] NColumns[A]],(Message[
          LeadingMatrixBlock::arg3];Abort[])];
    
    matarray=MatrixToArray[A,matarray];
    Do[matblock[i,j]=matarray[i,j], {i,1,a}, {j,1,b}];
    Return[ArrayToMatrix[a,b,matblock]]]

(* --22------------------------------------------------- *)

TailingMatrixBlock::argx="Three arguments required.";
TailingMatrixBlock::arg1="First argument is a matrix of CREs.";
TailingMatrixBlock::arg2="Second argment is a valid row index.";
TailingMatrixBlock::arg3="Third argument is a valid column index.";

TailingMatrixBlock[]:= (Message[TailingMatrixBlock::argx]; Abort[])

 TailingMatrixBlock[argseq__]:=
  Module[{args=List[argseq],A,a,b,matarray, matblock,m,n,i,j},
    If[Length[args]\[NotEqual]3,(Message[TailingMatrixBlock::argx];
        Abort[])];
    A=Part[args,1]; a=Part[args,2]; b=Part[args,3];
    m=Length[A];n=Length[First[A]];
    
    If[!ExprMatrixQ[A],(Message[TailingMatrixBlock::arg1];Abort[])];
    If[!And[IntegerQ[a],a\[GreaterEqual]1,a\[LessEqual] NRows[A]],(Message[
          TailingMatrixBlock::arg2];Abort[])];
    If[!And[IntegerQ[b],1\[LessEqual] b, b\[LessEqual] NColumns[A]],(Message[
          TailingMatrixBlock::arg3];Abort[])];
    
    matarray=MatrixToArray[A,matarray];
    Do[matblock[i,j]=matarray[m-a+i,n-b+j], {i,1,a}, {j,1,b}];
    Return[ArrayToMatrix[a,b,matblock]]]

(* ---23-------------------------------------------------- *)
MatrixDiagonal::argx="One argument required.";
MatrixDiagonal::arg1="Argument is a square matrix of CREs.";

MatrixDiagonal[]:= (Message[MatrixDiagonal::argx]; Abort[])

MatrixDiagonal[argseq__]:= Module[{args=List[argseq],A,m,diag={},j},
    If[Length[args]\[NotEqual]1,(Message[MatrixDiagonal::argx];
        Abort[])];
    A=Part[args,1]; 
    If[!ExprSquareMatrixQ[A],(Message[MatrixDiagonal::arg1];Abort[])];
    m=Length[A];
    Do[diag=Join[{Part[Part[A,j],j]},diag], {j,1,m}];
    Return[Reverse[diag]]]

(* --24------------------------------------------------ *)

DiagonalToMatrix::argx="One argument required.";
DiagonalToMatrix::arg1="Argument is a list of CREs.";

DiagonalToMatrix[]:= (Message[DiagonalToMatrix::argx]; Abort[])

DiagonalToMatrix[argseq__]:=Module[{args=List[argseq],diag,n,matarrat,i,j},
    If[Length[args]\[NotEqual]1,(Message[DiagonalToMatrix::argx];
        Abort[])];
    diag=Part[args,1];
    If[!ExprListQ[diag],(Message[DiagonalToMatrix::arg1];Abort[])];
    n=Length[diag];
    Do[matarray[i,j]=0, {i,1,n},{j,1,n}];
    Do[matarray[i,i]=Part[diag,i], {i,1,n}];
    Return[ArrayToMatrix[n,n,matarray]]]

(* ---25----------------------------------------------- *)

MatrixUpperTriangle::argx="One argument required.";
MatrixUpperTriangle::arg1="Argument is a square matrix of CREs.";

MatrixUpperTriangle[]:= (Message[MatrixUpperTriangle::argx]; Abort[])

MatrixUpperTriangle[argseq__]:= 
  Module[{args=List[argseq],A,out={},m,i},
    If[Length[args]\[NotEqual]1,(Message[MatrixUpperTriangle::argx];
        Abort[])];
    A=Part[args,1];
    If[!ExprSquareMatrixQ[A],(Message[MatrixUpperTriangle::arg1];Abort[])];
    m=Length[A];
    Do[out=Join[{Take[Part[A,i],{i,m}]}, out], {i,1,m}];
    Return[Reverse[out]]]

(* ---26----------------------------------------------- *)

MatrixLowerTriangle::argx="One argument required.";
MatrixLowerTriangle::arg1="Argument is a square matrix  of CREs.";

MatrixLowerTriangle[]:= (Message[MatrixLowerTriangle::argx]; Abort[])

MatrixLowerTriangle[argseq__]:=  Module[{args=List[argseq],A,out,m},
    If[Length[args]\[NotEqual]1,(Message[MatrixLowerTriangle::argx];
        Abort[])];
    A=Part[args,1];
    If[!ExprSquareMatrixQ[A],(Message[MatrixLowerTriangle::arg1];Abort[])];
    m=Length[A];out={Part[A,m]};
    Do[out=Join[{Take[Part[A,m-i+1],m-i+1]}, out], {i,2,m}];
    Return[out]]

(* --27---------------------------------------------------- *)

UpperTriangleToMatrix::argx="One argument required.";
UpperTriangleToMatrix::arg1=
    "Argument is a list of lists of CREs in strict descending order of length.";

UpperTriangleToMatrix[]:= (Message[UpperTriangleToMatrix::argx]; Abort[])

UpperTriangleToMatrix[argseq__] := Module[{args = List[argseq], T, m, 
  n, matarray}, 
    If[Length[args] \[NotEqual]1, (Message[UpperTriangleToMatrix::argx];
        Abort[])];
    T = Part[args, 1];
    If[Not[ExprUpperTriangularQ[T]], (Message[UpperTriangleToMatrix::arg1]; 
           Abort[])];

    m = Length[T]; n = Length[First[T]]; Do[matarray[i, j] = 0, {i, 1, n}, {j, 1,
     m}];
    Do[matarray[i, i + j - 1] = Part[Part[T, i], j], {i, 1, m}, {j, 1, n - i + 1}];
    Return[ArrayToMatrix[m, n, matarray]]]

(* --28---------------------------------------------------- *)

LowerTriangleToMatrix::argx="One argument required.";
LowerTriangleToMatrix::arg1="Argument is a list of lists of CREs in strict ascending \
order of length.";

LowerTriangleToMatrix[]:= (Message[LowerTriangleToMatrix::argx]; Abort[])

LowerTriangleToMatrix[argseq__]:=Module[{args=List[argseq],T,m,n,matarray},
    If[Length[args]\[NotEqual]1,(Message[LowerTriangleToMatrix::argx];
        Abort[])];
    T=Part[args,1];
    If[!ExprLowerTriangularQ[T],(Message[LowerTriangleToMatrix::arg1];
        Abort[])];

    m=Length[T];n=Length[Last[T]];Do[matarray[i,j]=0, {i,1,m},{j,1,n}];
    Do[matarray[i,j]=Part[Part[T,i],j],{i,1,m},{j,1,i}];
    Return[ArrayToMatrix[m,n,matarray]]]


(* ============================================================ *)

(* SquareMatrixQ[A_]:=And[MatrixQ[A], NRows[A]\[Equal]NColumns[A]] *)

LowerTriangularQ[A_]:=And[ListOfListsQ[A],AscendingLengthsQ[A]]
ExprLowerTriangularQ[a_]:= And[LowerTriangularQ[a],ExprListQ[Apply[Join,a]]]

UpperTriangularQ[A_]:=LowerTriangularQ[Reverse[A]]
ExprUpperTriangularQ[a_]:= And[UpperTriangularQ[a],ExprListQ[Apply[Join,a]]]

AscendingLengthsQ[A_]:=Module[{n=Length[A], m=Length[Last[A]],ans=True},
    If[n\[NotEqual] m, Return[False]];
    Do[If[
        Length[Part[A,j]]\[NotEqual] j,
          ans=False,Break[]],{j,1,n}];
    Return[ans]]


(*---------------------------------------------------*)
IwasawaQ::argx = "One argument required.";

IwasawaQ[]:= (Message[IwasawaQ::argx]; Abort[])

IwasawaQ[argseq__]:= Module[{args=List[argseq],m},
        If[Length[args] != 1,(Message[IwasawaQ::argx];
        Abort[])];
    m=Part[args,1];
    Return[And[ExprSquareMatrixQ[m],
                   NonZeroDiagonalTermsQ[m],
                   UpperTriangularMatrixQ[m],
                   m[[Length[m],Length[m]]]\[Equal]1]]]
(*---------------------------------------*)

UpperTriangularMatrixQ[m_]:=Module[{n=Length[m], ans=True}, 
    (* must be a sq matrix *)
    Do[If[m[[i,j]]\[NotEqual]0, {ans=False,Break[]}],{i,2,n},{j,1,i-1}];
    Return[ans]]

NonZeroDiagonalTermsQ[m_]:= Module[{n=Length[m], ans=True}, (* must be a sq matrix *)
    Do[If[m[[i,i]]\[Equal]0, {ans=False,Break[]}],{i,1,n}];
    Return[ans]]

ListOfListsQ[x_]:=And[ListQ[x],Apply[And,Map[ListQ,x]]]

ArrayToMatrix[m_,n_,array_]:=Module[{rowi, A={},i,j},
    Do[{rowi={},Do[rowi=Join[{array[i,j]},rowi],{j,1,n}],
        A=Join[{Reverse[rowi]},A]}, {i,1,m}];
    Return[Reverse[A]]]

MatrixToArray[A_,array_]:=Module[{i,j,m=Length[A],n=Length[First[A]]},
    Do[array[i,j]=Part[Part[A,i],j], {i,1,m},{j,1,n}];
    Return[array]]

SwapListElements[x_,i_,j_]:= (* i \[NotEqual] j *)
  Module[{left,right,u,between, m=Length[x],eli, elj},
    If[j<i,{u=j,j=i,i=u}];
    eli=Part[x,i]; elj=Part[x,j];
    left=Take[x,i-1]; between=Take[x,{i+1,j-1}];
    right=Take[x,m-j];
    Return[Join[left, {elj},between, {eli}, right]]]


MapJoin[a_,b_]:= Module[{ans={},n=Length[a]},
    Do[ans=Join[{Join[Part[a,i],Part[b,i]]},ans],{i,1,n}];
    Return[Reverse[ans]]]

(* Wedge, d====================================================== 

Unprotect["Global`*"]; ClearAll["Global`*"]; Remove["Global`*"]; 
Unprotect[In, Out]; Clear[In, Out]; Protect[In, Out]; $Line = 0; 
$RecursionLimit = 256; $IterationLimit = 4096;


Forms[i_] := {}; AllScalForms = {}; AllDifForms = {}; AllSymbols = 
{}; FormVars = {_Wedge, _d}; HeadList = {}; ScalHeadList = {}; 
nodHeads = {Bar, Pattern, Condition, RuleDelayed, 
  SeriesData}; $EDCversion = 360;

zeroQ[0] = True; zeroQ[x_List] := And @@ (zeroQ /@ Union[Flatten[x]]);
 zeroQ[x_SeriesData] := If[x[[3]] === {}, True, False]; 
zeroQ[x_] := False;

SetAttributes[Bar, {Listable}]; Bar[Bar[x_]] = x; 
Bar[Complex[u_, v_]] := Complex[u, -v]; 
Bar[x_Plus | x_Times | x_Wedge | x_Power | x_Rule | x_Equal] := 
 Bar /@ x; 
Bar[x_SeriesData] := 
 x /. {First[x] -> Bar[First[x]], x[[2]] -> Bar[x[[2]]], 
   x[[3]] -> Bar /@ x[[3]]}; 
Bar[DirectedInfinity[x_]] := DirectedInfinity[Bar[x]]; 
Bar[d[x_]] := d[Bar[x]]; 
Bar[Derivative[x__][y_][z__]] := 
 If[Union[Bar[{z}]] === Union[{z}], 
  Derivative[
     Sequence @@ 
      Permutations[{x}][[
       Position[Permutations[{z}], Bar[{z}]][[1, 1]]]]][Bar[y]][z], 
  Derivative[x][Bar[y]][Sequence @@ Bar[{z}]]]; 
Bar[x_] := x /; NumericQ[x];


FormDegree[x_Plus] := FormDegree[First[x]]; 
FormDegree[x_Times] := Plus @@ FormDegree /@ List @@ x;
FormDegree[x_Wedge] := Plus @@ FormDegree /@ List @@ x;
FormDegree[d[x_]] := FormDegree[d[x]] = 1 + FormDegree[x];
FormDegree[x_List] := 
  FormDegree[First[Select[Union[Flatten[x]], ! zeroQ[#] &]]];
FormDegree[Bar[x_]] := FormDegree[Bar[x]] = FormDegree[x]; 
FormDegree[x_SeriesData] := If[x[[3]] === {}, 0, FormDegree[x[[3]]]];
FormDegree[x_] := 0;

DeclareForms[z__] := 
  Block[{h, x = {{z}}, xi, rxi, k, oldHeads, newHeads}, 
   While[Head[x[[1, 1]]] === List, x = First[x]]; 
   Do[xi = x[[i]]; rxi = Rest[xi]; 
    Forms[First[xi]] = Union[Forms[First[xi]], rxi];
    AllScalForms = Union[AllScalForms, rxi];
    Do[h = Head[rxi[[j]]]; 
     If[h === Symbol, FormDegree[rxi[[j]]] = First[xi]; 
      BasicScalFormQ[rxi[[j]]] = True;, FormDegree[_h] = First[xi]; 
      BasicScalFormQ[_h] = True;], {j, Length[rxi]}], {i, Length[x]}];
   AllDifForms = Union[AllScalForms]; oldHeads = ScalHeadList; 
   ScalHeadList = 
    Complement[Union[Head /@ AllScalForms, ScalHeadList], {Symbol}]; 
   newHeads = Complement[ScalHeadList, oldHeads]; 
   AllSymbols = Union[AllDifForms];
   HeadList = Union[Head /@ AllSymbols]; 
   DifFormSymbols = Drop[AllDifForms, -Length[ScalHeadList]];
   k = Thread[Blank[ScalHeadList]]; 
   FormVars = 
    Flatten[{_Wedge, _d, Union[k, Bar[k]], 
      u_ | Bar[u_] /; MemberQ[DifFormSymbols, u], _tr}];];


NoDif[z__] := (nodHeads = Union[nodHeads, Flatten[{z}]];)

BasicScalFormQ[Bar[x_]] := BasicScalFormQ[x]; 
BasicScalFormQ[_d] = True; BasicScalFormQ[x_] := False;

ScalFormQ[x_Times] := Or @@ ScalFormQ /@ List @@ x; 
ScalFormQ[x_Wedge | x_Plus] := And @@ ScalFormQ /@ List @@ x; 
ScalFormQ[x_] := BasicScalFormQ[x];


SetAttributes[d, {Listable}]; 
d[x_Times | x_Wedge] := 
 d[First[x]]\[Wedge]Rest[x] + (-1)^FormDegree[First[x]]*
   First[x]\[Wedge]d[Rest[x]]; d[x_?NumericQ | x_d] = 0; 
d[Power[y_, n_]] := n y^(n - 1) d[y] + y^n Log[y] d[n]; 
d[x_Plus] := d /@ x; 
HoldPattern[d[Bar[x_]]] := 
 With[{evalHx = d[x]}, Bar[evalHx] /; Head[evalHx] =!= d]; 
d[x_Rule | x_Equal] := reWrite[d /@ x]; 
d[x_SeriesData] := (x /. x[[3]] -> d[x[[3]]]) + 
  Wedge[d[First[x]], D[x, First[x]]]; 
d[h_[y__] /; ! MemberQ[nodHeads, h]] := 
 Sum[(Derivative[
        Sequence @@ 
         RotateRight[Append[Table[0, {Length[{y}] - 1}], 1], i]][h][
      y]) d[{y}[[i]]], {i, Length[{y}]}] /; 
  FormDegree[h[y]] === 0 && ! 
    MemberQ[{Integer, Blank, Pattern, Condition}, Head[First[{y}]]];


newSer$[x_SeriesData, k_] := 
 SeriesData[First[x], x[[2]], 
  Flatten[Transpose[
    Prepend[Table[Table[0, {Length[x[[3]]]}], {k - 1}], x[[3]]]]], 
  k x[[4]], k x[[5]], k Last[x]]

Wedge[x_] := x /; Length[{x}] < 2 && Head[{x}[[1]]] =!= Pattern

Default[Wedge] := 1; SetAttributes[Wedge, {Flat, OneIdentity}]; 
Wedge[0, y__] = 0; Wedge[x__, 0] = 0; 
Wedge[x_SeriesData, y_SeriesData] := 
 Block[{x$, y$, r1, r2, res, x3, y3}, 
   If[Last[x] === Last[y], x$ = x; y$ = y, 
    x$ = newSer$[x, LCM[Last[x], Last[y]]/Last[x]]; 
    y$ = newSer$[y, LCM[Last[x], Last[y]]/Last[y]]]; 
   r1 = x$[[-3]] + y$[[-3]]; 
   r2 = Min[x$[[-2]] + y$[[-3]], x$[[-3]] + y$[[-2]]]; 
   If[Length[x$[[3]]] < x$[[-2]] - x$[[-3]], 
    x3 = Join[x$[[3]], 
      Table[0, {x$[[-2]] - x$[[-3]] - Length[x$[[3]]]}]], 
    x3 = x$[[3]]]; 
   If[Length[y$[[3]]] < y$[[-2]] - y$[[-3]], 
    y3 = Join[y$[[3]], 
      Table[0, {y$[[-2]] - y$[[-3]] - Length[y$[[3]]]}]], 
    y3 = y$[[3]]]; 
   Which[r2 === r1, res = {}, r2 - r1 === x$[[-2]] - x$[[-3]], 
    res = Sum[
      Map[Wedge[x3[[i]], #] &, 
       Join[Table[0, {i - 1}], Drop[y3, -i + 1]]], {i, Length[x3]}], 
    True, res = 
     Sum[Map[Wedge[#, y3[[i]]] &, 
       Join[Table[0, {i - 1}], Drop[x3, -i + 1]]], {i, Length[y3]}]]; 
   SeriesData[First[x$], x$[[2]], res, r1, r2, Last[x$]]] /; 
  First[x] === First[y] && x[[2]] === y[[2]]; 
Wedge[y_, x_SeriesData] := x /. x[[3]] -> Map[Wedge[y, #] &, x[[3]]]; 
Wedge[x_SeriesData, y_] := x /. x[[3]] -> Map[Wedge[#, y] &, x[[3]]];
Wedge[x__, y_Plus] := Plus @@ Map[Wedge[x, #] &, List @@ y];
Wedge[x_Plus, y__] := Plus @@ Map[Wedge[#, y] &, List @@ x]; 
Wedge[u__, Times[x_, y_]] := 
 Times[x, Wedge[u, y]] /; NumericQ[x] || FormDegree[x] === 0;
Wedge[Times[x_, y_], z__] := 
 Times[x, Wedge[y, z]] /; NumericQ[x] || FormDegree[x] === 0; 
x_^n_.\[Wedge]y_ := x^n*y /; FormDegree[x] === 0;
y_\[Wedge]x_^n_. := x^n*y /; FormDegree[x] === 0;
Wedge[x_, y___, x_] := 0 /; OddQ[FormDegree[x]] && ScalFormQ[x]; 
Wedge[x__] := 
 Block[{doubL = Transpose[{FormDegree /@ {x}, {x}}]}, 
   Signature[Select[doubL, OddQ[First[#]] &]] Wedge @@ 
     Map[Last[#1] &, Sort[doubL]]] /; 
  Union[BasicScalFormQ /@ {x}] === {True} && 
   Map[Last[#1] &, Sort[Transpose[{FormDegree /@ {x}, {x}}]]] =!= {x};

Wedge[x_] := x /; Length[{x}] < 2 && Head[{x}[[1]]] =!= Pattern

Default[Wedge] := 1; SetAttributes[Wedge, {Flat, OneIdentity}]; 
Wedge[0, y__] = 0; Wedge[x__, 0] = 0; 
Wedge[x_SeriesData, y_SeriesData] := 
 Block[{x$, y$, r1, r2, res, x3, y3}, 
   If[Last[x] === Last[y], x$ = x; y$ = y, 
    x$ = newSer$[x, LCM[Last[x], Last[y]]/Last[x]]; 
    y$ = newSer$[y, LCM[Last[x], Last[y]]/Last[y]]]; 
   r1 = x$[[-3]] + y$[[-3]]; 
   r2 = Min[x$[[-2]] + y$[[-3]], x$[[-3]] + y$[[-2]]]; 
   If[Length[x$[[3]]] < x$[[-2]] - x$[[-3]], 
    x3 = Join[x$[[3]], 
      Table[0, {x$[[-2]] - x$[[-3]] - Length[x$[[3]]]}]], 
    x3 = x$[[3]]]; 
   If[Length[y$[[3]]] < y$[[-2]] - y$[[-3]], 
    y3 = Join[y$[[3]], 
      Table[0, {y$[[-2]] - y$[[-3]] - Length[y$[[3]]]}]], 
    y3 = y$[[3]]]; 
   Which[r2 === r1, res = {}, r2 - r1 === x$[[-2]] - x$[[-3]], 
    res = Sum[
      Map[Wedge[x3[[i]], #] &, 
       Join[Table[0, {i - 1}], Drop[y3, -i + 1]]], {i, Length[x3]}], 
    True, res = 
     Sum[Map[Wedge[#, y3[[i]]] &, 
       Join[Table[0, {i - 1}], Drop[x3, -i + 1]]], {i, Length[y3]}]]; 
   SeriesData[First[x$], x$[[2]], res, r1, r2, Last[x$]]] /; 
  First[x] === First[y] && x[[2]] === y[[2]]; 
Wedge[y_, x_SeriesData] := x /. x[[3]] -> Map[Wedge[y, #] &, x[[3]]]; 
Wedge[x_SeriesData, y_] := x /. x[[3]] -> Map[Wedge[#, y] &, x[[3]]];
Wedge[x__, y_Plus] := Plus @@ Map[Wedge[x, #] &, List @@ y];
Wedge[x_Plus, y__] := Plus @@ Map[Wedge[#, y] &, List @@ x]; 
Wedge[u__, Times[x_, y_]] := 
 Times[x, Wedge[u, y]] /; NumericQ[x] || FormDegree[x] === 0;
Wedge[Times[x_, y_], z__] := 
 Times[x, Wedge[y, z]] /; NumericQ[x] || FormDegree[x] === 0; 
x_^n_.\[Wedge]y_ := x^n*y /; FormDegree[x] === 0;
y_\[Wedge]x_^n_. := x^n*y /; FormDegree[x] === 0;
Wedge[x_, y___, x_] := 0 /; OddQ[FormDegree[x]] && ScalFormQ[x]; 
Wedge[x__] := 
 Block[{doubL = Transpose[{FormDegree /@ {x}, {x}}]}, 
   Signature[Select[doubL, OddQ[First[#]] &]] Wedge @@ 
     Map[Last[#1] &, Sort[doubL]]] /; 
  Union[BasicScalFormQ /@ {x}] === {True} && 
   Map[Last[#1] &, Sort[Transpose[{FormDegree /@ {x}, {x}}]]] =!= {x};

simpRules = {Cos[z_]^2*(x_.) + (x_.)*Sin[z_]^2 -> x}; varList = {};

coll[x_] := Collect[x, {_Wedge, _tr}, Factor];

SetAttributes[reWrite, {Listable}]; reWrite[0] = 0; 
reWrite[x_Equal] := Equal[reWrite[First[x] - Last[x]], 0]; 
reWrite[x_Rule] := Rule[First[x], reWrite[Last[x]]]; 
reWrite[x_SeriesData] := (x /. x[[3]] -> reWrite[x[[3]]]); 
reWrite[x_] := (Collect[coll[x] /. simpRules, FormVars, 
     bestFacRul] /. simpRules) /; FormDegree[x] > 0; 
reWrite[x_] := bestFacRul[x /. simpRules] /. simpRules;

FreeALLQ[x_, y_List] := And @@ Map[FreeQ[x, #] &, y];

bestFacRul[x_] := 
 If[varList === {} || FreeALLQ[x, varList], minFacRul[x], 
  varFacRul[x]]; 
varFacRul[x_] := 
 Collect[Expand[x] /. simpRules, varList, Factor] /. simpRules; 
minFacRul[x_] := Factor[Expand[x] /. simpRules] /. simpRules;

Unprotect[SeriesData]; If[$VersionNumber < 5, 
 Times[x_SeriesData, 0] ^= 0]; 
SeriesData[x1_, x2_, x3_, x4_, x5_, x6_] := 
 Plus @@ If[x2 === Infinity, 
    Map[First[#]/x1^((x4 + Last[#] - 1)/x6) &, 
      Transpose[{x3, Range[Length[x3]]}]] + 
     SeriesData[x1, x2, {}, x4, x5, x6],   
    Map[First[#]*(x1 - x2)^((x4 + Last[#] - 1)/x6) &, 
      Transpose[{x3, Range[Length[x3]]}]] + 
     SeriesData[x1, x2, {}, x4, x5, x6]] /; ! 
   FreeQ[x3, x1]; Protect[SeriesData];

=============================================================================*)


(* Positive definite Unimodular matrices--------------------------*)
VolumeFormUnimodular::argx="Two arguments required.";
VolumeFormUnimodular::arg1="First argument is a string.";
VolumeFormUnimodular::arg2="Second argument is a positive integer >= 2";

VolumeFormUnimodular[]:= (Message[VolumeFormUnimodular::argx]; Abort[])

VolumeFormUnimodular[argseq__]:=
  Module[{ars=List[argseq],astr,n,i,j,x},
    If[Length[ars]\[NotEqual] 2, (Message[VolumeFormUnimodular::argx];
        Abort[])];
    astr=First[ars]; n=Part[ars,2];
    If[!StringQ[astr],(Message[VolumeFormUnimodular::arg1];Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 2],(Message[
          VolumeFormUnimodular::arg2];Abort[])];
    x=Symbol[astr];
    Return[Apply[Wedge,Apply[Join, Table[d[x[i,j]], {i,1,n},{j,i+1,n}]]]]]

(* Upper triangular real matrices with 1's on diagonal: x is a symbol*)
(* Positive definite diagonal matrices--------------------------*)

VolumeFormDiagonal::argx="Two arguments required.";
VolumeFormDiagonal::arg1="First argument is a string.";
VolumeFormDiagonal::arg2="Second argument is a positive integer.";

VolumeFormDiagonal[]:= (Message[VolumeFormDiagonal::argx]; Abort[])

VolumeFormDiagonal[argseq__]:=
  Module[{ars=List[argseq],astr,n,j,a,alis,coefficient},
    If[Length[ars]\[NotEqual] 2, (Message[VolumeFormDiagonal::argx];
        Abort[])];
    astr=First[ars]; n=Part[ars,2];
    If[!StringQ[astr],(Message[VolumeFormDiagonal::arg1];Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 1],(Message[
          VolumeFormDiagonal::arg2];Abort[])];
    a=Symbol[astr];
    (* DeclareForms[{0,a[_]}]; *)
    alis=Table[a[j], {j,1,n}];
    coefficient=1/Apply[Times,alis];
    Return[coefficient*Apply[Wedge,Map[d,alis]]]]

(* ---------------------------------------------------------------*)
VolumeFormGln::argx="Two arguments required.";
VolumeFormGln::arg1="First argument is a string.";
VolumeFormGln::arg2="Second argument is a positive integer n >=  1";

volumeFormGln[]:= (Message[VolumeFormGln::argx]; Abort[])

VolumeFormGln[argseq__]:=
  Module[{ars=List[argseq],astr,n,i,j,g,detg,form,alis,coefficient},
    If[Length[ars]\[NotEqual] 2, (Message[VolumeFormGln::argx];
        Abort[])];
    astr=First[ars]; n=Part[ars,2];
    If[!StringQ[astr],(Message[VolumeFormGln::arg1];Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 1],(Message[VolumeFormGln::arg2];
        Abort[])];
    g=Symbol[astr];
    detg=Det[Table[g[i,j], {i,1,n},{j,1,n}]];
    form=Apply[Wedge,Map[d,Apply[Join,Table[g[i,j], {i,1,n},{j,1,n}]]]];
    Return[(1/detg^n )\[Wedge]form]]

(*------------------------------------------------------------------*)
VolumeFormHn::argx="Three arguments required.";
VolumeFormHn::arg1="First argument is a string.";
VolumeFormHn::arg2="Second argument is a string.";
VolumeFormHn::arg3="Second argument is a positive integer n >=  1";

VolumeFormHn[]:= (Message[VolumeFormHn::argx]; Abort[])

VolumeFormHn[argseq__]:=
  Module[{ars=List[argseq],astr1,astr2,n,i,j,k,x,y,xform,yform,coef},
    If[Length[ars]\[NotEqual] 3, (Message[VolumeFormHn::argx];
        Abort[])];
    astr1=First[ars];astr2=Part[ars,2]; n=Part[ars,3];
    If[!StringQ[astr1],(Message[VolumeFormHn::arg1];Abort[])];
    If[!StringQ[astr2],(Message[VolumeFormHn::arg2];Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 1],(Message[VolumeFormHn::arg2];
        Abort[])];
    x=Symbol[astr1];y=Symbol[astr2];
    xform=Apply[Wedge,Apply[Join, Table[d[x[i,j]], {i,1,n},{j,i+1,n}]]];
    coef=Product[y[k]^(-k(n-k)-1), {k,1,n-1}];
    yform=Apply[Wedge,Table[d[y[i]], {i,1,n-1}]];
    Return[coef \[Wedge]xform \[Wedge] yform]]

(*------------------------------------------------------*)
VolumeHn::argx="One argument is required.";
VolumeHn::arg1="Argument is a positive integer >=2.";

VolumeHn[]:= (Message[VolumeHn::argx]; Abort[])

VolumeHn[argseq__]:=Module[{ars=List[argseq],n,j},
    If[Length[ars]\[NotEqual] 1, (Message[VolumeHn::argx];Abort[])];
    n=First[ars];
    If[!And[IntegerQ[n],n\[GreaterEqual] 2],(Message[VolumeHn::arg1];
        Abort[])];
    Return[Product[Zeta[j]/VolumeSphere[1,j-1], {j,2,n}]]]

(*---------------------------------------------------*)
VolumeSphere::argx="Two arguments are required.";
VolumeSphere::arg1="Argument is a CRE or non-negative number.";
VolumeSphere::arg2="Argument is a positive integer >=1."

VolumeSphere[]:= (Message[VolumeSphere::argx]; Abort[])

VolumeSphere[argseq__]:=Module[{ars=List[argseq],n,r,m},
    If[Length[ars]\[NotEqual] 2, (Message[VolumeSphere::argx];
        Abort[])];
    r=First[ars]; n=Part[ars,2];
    If[And[RealQ[r],r<0],(Message[VolumeSphere::arg1];
        Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 1],(Message[VolumeSphere::arg2];
        Abort[])];m=n+1;
    Return[2 Pi^(m/2)/Gamma[m/2]r^(m-1)]]

(*--------------------------------------------------------*)

VolumeBall::argx="Two arguments are required.";
VolumeBall::arg1="An Expr. If numeric argument is non-negative.";
VolumeBall::arg2="Argument is a positive integer >=2.";

VolumeBall[]:= (Message[VolumeBall::argx]; Abort[])

VolumeBall[argseq__]:=Module[{ars=List[argseq],n,r},
    If[Length[ars]\[NotEqual] 2, (Message[VolumeBall::argx];
        Abort[])];
    r=First[ars]; n=Part[ars,2]; 
    If[And[RealQ[r],r<0],(Message[VolumeBall::arg1];
        Abort[])];
    If[!And[IntegerQ[n],n\[GreaterEqual] 2],(Message[VolumeBall::arg2];
        Abort[])];
    Return[2 Pi^(n/2) /(n Gamma[n/2]) r^n]]
(*==================================================================================*)


RealQ[x_]:= And[NumberQ[x], x \[Element] Reals]

WhittakerGamma::argx="One argument required.";
WhittakerGamma::arg1="Argument is a list of CREs.";
WhittakerGamma::arg2="If numeric the real part of each list element is > 1/n.";

WhittakerGamma[]:= (Message[WhittakerGamma::argx]; Abort[])

WhittakerGamma[argseq__]:=Module[{ars=List[argseq],v,i,j,k,n,gam,nu,ans=True},
    If[Length[ars]\[NotEqual] 1, (Message[WhittakerGamma::argx];
        Abort[])];
    nu=First[ars];
    If[!ExprListQ[nu],(Message[WhittakerGamma::arg1];Abort[])];
    n=Length[nu]+1;
    Do[If[And[NumberQ[Part[nu,j]], Re[Part[nu,j]]<=1/n],ans=False],{j,1,Length[nu]}];

    If[Not[ans],(Message[WhittakerGamma::arg2];Abort[])];
   
    Do[v[j,k]=Sum[( n Part[nu,n-k+i]-1)/2,{i,0,j-1}],{j,1,n-1},{k,j,n-1}];
    gam=Product[Pi^(-1/2-v[j,k])*Gamma[1/2+v[j,k]],{j,1,n-1},{k,j,n-1}];
    Return[Simplify[gam]]]

(* -----------------------------------------------------*)


Whittaker::argx = "Four arguments required.";
Whittaker::arg1 = "First argument is a real non-singular matrix.";
Whittaker::arg2 = "Second argument is a complex list.";
Whittaker::arg3 = "Third argument is a non-zero integer list.";
Whittaker::arg5 = "An argument is incompatible or of incorrect type.";
Whittaker::arg4 = "Fourth argument is a symbol.";

Whittaker[] := (Message[Whittaker::argx]; Abort[])

Whittaker[argseq__] := 
  Module[{ars = List[argseq], z, nu, psi, n, iw, y, x, gam, m, ym, yy = {}, 
      w}, If[Length[ars] != 4, (Message[Whittaker::argx]; Abort[])];
    z = First[ars]; nu = Part[ars, 2]; psi = Part[ars, 3]; w = Part[ars, 4];
    If[Or[Not[MatrixQ[z]], Det[z] == 0], (Message[Whittaker::arg1]; 
        Abort[])];
    If[! ListQ[nu], (Message[Whittaker::arg2]; Abort[])];
    If[! ListQ[psi], (Message[Whittaker::arg3]; Abort[])];
    If[! gln`Private`SymbolQ[w], (Message[Whittaker::arg4]; 
        Abort[])]; (*check*)
    If[! And[NRows[z] == NColumns[z], NRows[z] == Length[nu] + 1, 
          ExprMatrixQ[z], ExprListQ[nu], ExprListQ[psi], 
          Length[psi] == Length[nu]], (Message[Whittaker::arg5];
        Abort[])]; 
    If[MemberQ[psi, 0], (Message[Whittaker::arg3]; Abort[])];
    n = NRows[z];
    iw = IwasawaForm[z];
    y = ExtractYVariablesFromIW[iw];
    x = ExtractXMatrixFromIW[iw];
    gam = WhittakerGamma[nu];
    ym = Abs[psi]*y;
    Return[{WhittakerCoefficient[nu, psi], UnipotentCharacter[psi, x], gam, 
        WhittakerDiagonal[Reverse[ym], nu, w]/gam}]]


WhittakerDiagonal[y_, v_, w_] := 
  StadeQ[Length[y] + 1, v, y]*StadeStarS[Length[y] + 1, v, y, w]


StadeQ[n_, nu_, y_] := 
  Module[{b, mu, p1, p2}, b[i_, j_] := If[i + j <= n, i j, (n - i)(n - j)];
    mu[j_] := Sum[Sum[n nu[[i]]/2, {i, k, k + j - 1}] - j/2, {k, 1, n - j}];
    p1 = Product[y[[i]]^(b[i, j]*nu[[j]]), {i, 1, n - 1}, {j, 1, n - 1}];
    p2 = Product[y[[j]]^(-mu[j]), {j, 1, n - 1}];
    Return[PowerExpand[p1*p2]]]

UnipotentCharacter[m_, x_] := 
  Exp[2 Pi I(Sum[m[[i]]x[[i, i + 1]], {i, 1, Length[m]}])]


WhittakerCoefficient[v_, m_] := Module[{l = Length[m], b, n}, n = l + 1;
    b[i_, j_] := If[i + j <= n, i j, (n - i)(n - j)];
    Return[
      Product[Abs[
              m[[i]]]^(Sum[b[i, j]*v[[j]], {j, 1, l}] - i(n - i)), {i, 1, 
          l}]]]

StadeStarS[n_, v_, y_, w_] :=  (* w is a symbol *)
   Switch[n, 2, SSS[2, v, y, w],
                         3, SSS[3, v, y, w],
                         4, SSS[4, v, y, w],
                         5, SSS[5, v, y, w],
                        _, StadeStarRS[n, v, y, w]]

SSS[2, v_, y_, u_] := 2  K[First[v] - 1/2, 2Pi First[y]]

SSS[3, v_, y_, u_] := 
  Module[{coef, val, eps = 0.000001, integrand3, v1, v2, y1, y2, x}, 
    v1 = Part[v, 1]; v2 = Part[v, 2]; y1 = Part[y, 1]; y2 = Part[y, 2]; 
    coef = 8;
    integrand3[w_] := 
      Module[{s = Sqrt[1 + 1/w^2], a = (3v1 + 3v2 - 2)/2}, 
        Return[K[a, 2Pi y1 s]*
            K[a, 2Pi y2 Sqrt[1 + w^2]]w^(3(v1 - v2)/2 - 1)]];
    val = Integrate[integrand3[u], {u, 0, Infinity}]; Return[coef*val]]


SSS[4, v_, y_, w_] := 
  Module[{coef, value, n = 4, integrand4, h, lam1, rho1, rho2, rho3, wlis, 
      r11, r21, r31, r13, r12, r22, mu1, mu2, mu3, ss, v1, v2, v3, y1, y2, 
      y3}, v1 = Part[v, 1]; v2 = Part[v, 2]; v3 = Part[v, 3]; y1 = Part[y, 1];
     y2 = Part[y, 2];
    y3 = Part[y, 3];
    coef = 2^(2n - 3);(*corrected coefficient not 2^(n - 1)*)lam1 = 2v2;
    rho1 = 2v1; rho2 = 2v2; rho3 = 2v3;
    r11 = rho1 - 1/2;
    r21 = rho1 + rho2 - 1;
    r31 = rho1 + rho2 + rho3 - 3/2;
    r13 = rho3 - 1/2;
    r12 = rho2 - 1/2;
    r22 = rho2 + rho3 - 1;
    mu1 = r11 + r12 + r13;
    mu2 = r21 + r22;
    mu3 = r31;
    integrand4[u1_, u2_] := 
      Module[{bb, uu, s1, s2, s3}, s1 = Sqrt[(1 + 1/u1^2)];
        s2 = Sqrt[(1 + u1^2)(1 + 1/u2^2)];
        s3 = Sqrt[(1 + u2^2)];
        bb = K[mu1, 2Pi y1 s1]*K[mu1, 2Pi y2 s2]*K[mu1, 2Pi y3 s3];
        uu = u1^(r11 - r13 - 1) u2^(r21 - r22 - 1);
        ss = StadeStarS[2, {lam1}, {y2 u1/u2}, {u1, u2}];
        Return[ss*bb*uu]];
    wlis = Table[w[i], {i, 1, 2}];
    Return[
      coef*Integrate[
          Apply[integrand4, wlis], {w[1], 0, Infinity}, {w[2], 0, 
            Infinity}]]]


SSS[5, v_, y_, w_] := 
  Module[{coef, value, eps = 0.000001, n = 5, integrand5, h, lam1, lam2, rho1,
       rho2, rho3, rho4, r11, r21, r31, r41, r13, r12, r22, mu1, mu2, mu3, 
      mu4, ss, v1, v2, v3, v4, y1, y2, y3, y4}, v1 = Part[v, 1]; 
    v2 = Part[v, 2]; v3 = Part[v, 3]; v4 = Part[v, 4];
    y1 = Part[y, 1]; y2 = Part[y, 2]; y3 = Part[y, 3]; y4 = Part[y, 4]; 
    coef = 2^(2n - 3);
    lam1 = n v2/(n - 2); lam2 = n v3/(n - 2);
    rho1 = n v1/2; rho2 = n v2/2; rho3 = n v3/2; rho4 = n v4/2;
    r11 = rho1 - 1/2;
    r21 = rho1 + rho2 - 1;
    r31 = rho1 + rho2 + rho3 - 3/2;
    r41 = rho1 + rho2 + rho3 + rho4 - 2;
    r13 = rho3 - 1/2;
    r12 = rho2 - 1/2;
    r22 = rho2 + rho3 - 1;
    r32 = rho2 + rho3 + rho4 - 3/2;
    r14 = rho4 - 1/2;
    r23 = rho3 + rho4 - 1;
    mu1 = r11 + r12 + r13 + r14;
    mu2 = r21 + r22 + r23;(*again why does Stade define all these?*)
    mu3 = r31 + r32;
    mu4 = r41;
    integrand5[u1_, u2_, u3_] := 
      Module[{bb, uu, s1, s2, s3, s4, ss}, s1 = Sqrt[(1 + 1/u1^2)];
        s2 = Sqrt[(1 + u1^2)(1 + 1/u2^2)];
        s3 = Sqrt[(1 + u2^2)(1 + 1/u3^2)];
        s4 = Sqrt[(1 + u3^2)];
        bb = 
          K[mu1, 2Pi y1 s1]*K[mu1, 2Pi y2 s2]*K[mu1, 2Pi y3 s3]*
            K[mu1, 2Pi y4 s4];
        uu = u1^(r11 - r14 - 1) u2^(r21 - r23 - 1) u3^(r31 - r32 - 1);
        ss = StadeStarS[3, {lam1, lam2}, {y2 u1/u2, y3 u2/u3}, w[1]];
        Return[bb*uu*ss]];
    
    Return[
      coef*Integrate[
          Apply[integrand5, {w[2], w[3], w[4]}], {w[4], 0, Infinity}, {w[3], 
            0, Infinity}, {w[2], 0, Infinity}]]]

StadeStarRS[n_, v_, y_, w_] := 
  Module[{coef, lam, lamlis = {}, rho, j, k, i, r, mu1, integrandn, hh, 
      limits, args, fun, index},
    If[Or[n == 1, n == 2], Return[1.0]];
    coef = 2^(2n - 3);
    Do[lam[j - 1] = n v[[j]]/(n - 1), {j, 2, n - 1}];
    Do[rho[j] = n v[[j]]/2, {j, 1, n - 1}];
    Do[r[j, k] = Sum[rho[i], {i, k, k + j - 1}], {j, 1, n - 1}, {k, 1, 
        n - j}];
    mu1 = Sum[rho[k] - 1/2, {k, 1, n - 1}];
    Do[lamlis = Join[lamlis, {lam[jj]}], {jj, 1, n - 3}];
    
    integrandn[u_] := 
      Module[{s, m = Length[u], bb, uu, ss, jj, ii, vlis}, 
        s[1] = Sqrt[1 + 1/u[[1]]^2];
        Do[s[jj] = Sqrt[(1 + u[[jj - 1]]^2)(1 + 1/u[[jj]]^2)], {jj, 2, m}];
        s[m + 1] = Sqrt[1 + u[[m]]^2];
        bb = Product[K[mu1, 2 Pi y[[jj]] s[jj]], {jj, 1, m + 1}];
        uu = Product[u[[ii]]^(r[ii, 1] - r[ii, m + 2 - ii]), {ii, 1, m}];
        vlis = {};
        Do[
          vlis = Join[vlis, {y[[jj + 1]]u[[jj]]/u[[jj + 1]]}], {jj, 1, 
            m - 1}];
        ss = StadeStarS[m, lamlis, vlis, w];
        
        Return[bb*uu*ss]];
        index = TotalWIndex[n];
        limits = Table[{w[j], 0, Infinity}, {j, index - n + 3, index}];
        args = Table[w[j], {j, index - n + 3, index}];
        fun = integrandn[args];
    Return[coef*Apply[Integrate, Join[{fun}, limits]]]]

TotalWIndex[n_] := 
  Module[{s = 0, j = 2}, While[n - j > 0, {s = s + n - j, j = j + 2}]; 
    Return[s]]

(*-----------------------------------------------------*)

IntegerListQ[x_]:= And[ListQ[x],Apply[And,Map[IntegerQ,x]]]
RealMatrixQ[z_]:=And[MatrixQ[z],Apply[And, Map[RealQ,Apply[Join,z]]]]

RealQ[x_]:=And[NumericQ[x],Im[x]\[Equal]0]
ComplexListQ[x_]:=And[ListQ[x],Apply[And,Map[NumericQ,x]]]

ZeroPadLeft[x_,n_]:=
  Module[{m=Length[x],i},Join[Table[0,{i,1,n-m}],x]]

(*----------------------------------------------------*)
ExtractYVariablesFromIWold[mat_]:= 
  Module[{n=Length[mat], vars={},oldvar=1,newvar=1},
    Do[{newvar=mat[[n-i,n-i]]/oldvar,vars=Join[{newvar},vars],
        oldvar=newvar}, {i,1,n-1}];
    Return[Reverse[vars]]]

ExtractYVariablesFromIW[mat_] := Module[{n = Length[mat], vars = {}, oldd = 1, newd}, 
   Do[(newd = mat[[n - i, n - i]]; 
       vars = Join[{newd/oldd}, vars]; 
       oldd = newd), {i, 1, n - 1}];
    Return[Reverse[vars]]]

ExtractYVariables[mat_] := Module[{n = Length[mat],iw=IwasawaForm[mat],
                                   vars = {}, oldd = 1, newd}, 
   Do[(newd = iw[[n - i, n - i]]; 
       vars = Join[{newd/oldd}, vars]; 
       oldd = newd), {i, 1, n - 1}];
    Return[Reverse[vars]]]

ExtractYVariablesold[mat_]:=
  Module[{n=Length[mat],iw=IwasawaForm[mat], vars={},oldvar=1,newvar=1},
    Do[{newvar=iw[[n-i,n-i]]/oldvar,vars=Join[{newvar},vars],
        oldvar=newvar}, {i,1,n-1}];
    Return[Reverse[vars]]]

ExtractXMatrixFromIW[mat_]:=Module[{n=Length[mat],i,j,out},
    Do[out[i,j]=mat[[i,j]]/mat[[i,i]],{i,1,n},{j,1,n}];
    Return[ArrayToMatrix[n,n,out]]]

(*----------------------------------------------------*)
Ifun[v_,y_]:=Module[{i,j,b, expt,n=Length[y]+1},
    Do[b[i,j]= If[i+j\[LessEqual] n, i j,(n-i)(n-j)],{i,1,n}, {j,1,n}];
    Do[expt[i]=Sum[b[i,j]Part[v,j], {j,1,n-1}],{i,1,n-1}];
    Return[Product[y[[i]]^expt[i], {i,1,n-1}]]]
(*-----------------------------------------------------*)



FunctionalEquation::argx="Two arguments are required.";
FunctionalEquation::arg1="First argument is a list of CREs \
of length n-1.";
FunctionalEquation::arg2="Second argument is a positive integer 1<= i<= \
n-1.";


FunctionalEquation[]:=(Message[FunctionalEquation::argx];Abort[])

FunctionalEquation[argseq__]:=
  Module[{ars=List[argseq],n,nu,i,eqnfirst,eqni, eqnlast,eqnrest={},eqn},
    If[Length[ars]\[NotEqual] 2, (Message[FunctionalEquation::argx];
        Abort[])];
    nu=Part[ars,1];i=Part[ars,2];n=Length[nu]+1;
    If[!And[ExprListQ[nu],n\[GreaterEqual]2],(Message[FunctionalEquation::arg1];
        Abort[])];
    If[!And[IntegerQ[i],i\[LessEqual] n-1,1\[LessEqual] i],(Message[
          FunctionalEquation::arg2];Abort[])];
    If[n\[Equal]2,Return[{1/2-nu[[1]]}]];

    eqn = Table[nu[[l]], {l,1,n-1}];
    Switch[i,1,(eqn[[n-i-1]]=-1/n+nu[[n-i-1]]+nu[[n-i]];
                eqn[[n-i]]=2/n-nu[[n-i]]),
            n-1,(eqn[[n-i]]=2/n-nu[[n-i]];
                 eqn[[n-i+1]]=-1/n+nu[[n-i]]+nu[[n-i+1]]),
            _,(eqn[[n-i-1]]=-1/n+nu[[n-i-1]]+nu[[n-i]];
                  eqn[[n-i]]=2/n-nu[[n-i]];
                  eqn[[n-i+1]]=-1/n+nu[[n-i]]+nu[[n-i+1]])];
    Return[eqn]]
(*=================================================================================*)


(*======================================================*)
EisensteinSeriesTerm::argx="Three arguments required.";
EisensteinSeriesTerm::arg1="First argument is a CRE.";
EisensteinSeriesTerm::arg2="Second argument is a CRE.";
EisensteinSeriesTerm::arg3="Third argument is a list of two integers
with at least one non zero.";


EisensteinSeriesTerm[]:= (Message[EisensteinSeriesTerm::argx]; Abort[])

EisensteinSeriesTerm[argseq__]:=Module[{args=List[argseq],z,s,b,c,d},
    If[Length[args]\[NotEqual]3,(Message[EisensteinSeriesTerm::argx];
        Abort[])];
    
         z=Part[args,1];s=Part[args,2];b=Part[args,3];
    If[Not[ExprQ[z]],(Message[EisensteinSeriesTerm::arg1];Abort[])];
        If[Not[ExprQ[s]],(Message[EisensteinSeriesTerm::arg2];Abort[])];
    If[Not[And[ListQ[b],IntegerQ[c=Part[b,1]],IntegerQ[d=Part[b,2]],c^2+d^2 != 0]], 
       (Message[
          EisensteinSeriesTerm::arg3];Abort[])];
    
    Return[(1/2) Which[GCD[c,d]\[NotEqual]1,0,c^2+d^2\[Equal]0,0,True,
          Im[z]^s/Abs[c z+d]^(2s)]]]

EisensteinFourierCoefficient::argx="Three arguments required."
EisensteinFourierCoefficient::arg1="First argument is a CRE"

EisensteinFourierCoefficient::arg2="Second argument is a CRE."
EisensteinFourierCoefficient::arg3="Third argument is an integer."

EisensteinFourierCoefficient[]:= 
     (Message[EisensteinFourierCoefficient::argx]; Abort[])

EisensteinFourierCoefficient[argseq__]:=
  Module[{args=List[argseq],z,s,b,phi, coef},
    If[Length[args]\[NotEqual]3,(Message[EisensteinFourierCoefficient::argx];
        Abort[])];
         z=Part[args,1];s=Part[args,2];b=Part[args,3];y=Im[z];x=Re[z];
    
    If[Not[ExprQ[z]],(Message[EisensteinFourierCoefficient::arg1];Abort[])];
        If[
      Not[ExprQ[s]],(Message[EisensteinFourierCoefficient::arg2];Abort[])];
    If[Not[IntegerQ[b]], (Message[
          EisensteinFourierCoefficient::arg3];Abort[])];
    
    phi=Sqrt[Pi]Gamma[s-1/2]Zeta[2s-1]/(Gamma[s] Zeta[2s]);
    coef=2 Pi^s Sqrt[y]/(Gamma[s]Zeta[2s]);n=b;
    Return[If[n==0,
      y^s+phi* y^(1-s), 
        coef*DivisorSigma[1-2s,n]Abs[n]^(s-1/2) BesselK[s-1/2,
                  2Pi Abs[n]Im[z]]Exp[2Pi I n Re[z]]]]]
(*----------------------------------------------------------------------*)
HeckeFourierCoefficient::argx="Two arguments required."
HeckeFourierCoefficient::arg1="First argument is a list of functions or \
names."
HeckeFourierCoefficient::arg2="Second argument is a list of whole numbers."
HeckeFourierCoefficient::args="Both argument lists have length n-1."


HeckeFourierCoefficientOld[argseq__]:=Module[{args=List[argseq],fs,ns,n},
    If[Length[args]\[NotEqual]2,(Message[HeckeFourierCoefficient::argx];
        Abort[])];
         fs=Part[args,1];ns=Part[args,2];
    If[Not[
        And[ListQ[fs],ListQ[ns],Length[fs]\[Equal]Length[ns]]], (Message[
          HeckeFourierCoefficient::args];Abort[])];
    If[Not[AllWholeNumbersQ[ns]],(Message[HeckeFourierCoefficient::arg2];
        Abort[])];
    n=Length[fs]+1;
    Which[n\[Equal]2,Return[EvaluateHeckeGL2[First[fs],First[ns]]],n\[Equal]3,
      Return[Apply[EvaluateHeckeGL3, Join[fs,ns]]],Not[MemberQ[{2,3},n]],
      Return[False]]]

EvaluateHeckeGL2[a_,n_]:= Module[{facts, val=1,term,pow,p},
    If[n\[Equal]1, Return[1]];
    facts=FactorInteger[n];
    Do[{term=Part[facts,j],p=First[term],pow=Part[term,2],
        val=val* HeckeRecurseGL2[a,p,pow]},{j,1,Length[facts]}];
    Return[val]]

HeckeRecurseGL2[a_,p_,pow_]:=Which[pow\[Equal]0,1,
    pow==1,a[p],
    pow==2,a[p]^2-1,
    True,a[p]*HeckeRecurseGL2[a,p,pow-1]-HeckeRecurseGL2[a,p,pow-2]]

EvaluateHeckeGL3Left[l_,r_,m_]:= Module[{facts, val=1,term,pow,am1=1,p},
    If[m\[Equal]1, Return[1]];
    facts=FactorInteger[m];
    Do[{term=Part[facts,j],p=First[term],pow=Part[term,2],
        am1=am1* HeckeRecurseGL3Left[l,r,p,pow]},{j,1,Length[facts]}];
    Return[am1]]

EvaluateHeckeGL3Right[l_,r_,n_]:= Module[{facts, val=1,term,pow,a1n=1,p},
    If[n\[Equal]1, Return[1]];
    facts=FactorInteger[n];
    Do[{term=Part[facts,j],p=First[term],pow=Part[term,2],
        a1n=a1n* HeckeRecurseGL3Right[l,r,p,pow]},{j,1,Length[facts]}];
    Return[a1n]]

HeckeRecurseGL3Left[l_,r_,p_,pow_]:=Which[pow\[Equal]0,1,
    pow==1,l[p],
    pow==2,l[p]^2-r[p],
    True,l[p]*HeckeRecurseGL3Left[l,r,p,pow-1]-
      r[p]*HeckeRecurseGL3Left[l,r,p,pow-2]]

HeckeRecurseGL3Right[l_,r_,p_,pow_]:=Which[pow\[Equal]0,1,
    pow==1,r[p],
    pow==2,r[p]^2-l[p],
    True,r[p]*HeckeRecurseGL3Right[l,r,p,pow-1]-
      l[p]*HeckeRecurseGL3Right[l,r,p,pow-2]]

EvaluateHeckeGL3[l_,r_,m_,n_]:=If[1\[Equal]GCD[m,n ],
    EvaluateHeckeGL3Left[l,r,m]*EvaluateHeckeGL3Right[l,r,n],
    DoubleHeckeRecurseGL3[l,r,m,n]]

DoubleHeckeRecurseGL3[l_,r_,m_,n_]:=Module[{g=GCD[m,n],divs},
    divs=Complement[Divisors[g],{1}];
    Return[
      EvaluateHeckeGL3Left[l,r,m]*EvaluateHeckeGL3Right[l,r,n]-
        Sum[EvaluateHeckeGL3[l,r,m/Part[divs,i],n/Part[divs,i]], {i,1,
            Length[divs]}]]]

AllWholeNumbersQ[x_]:=Apply[And,Map[PositiveIntegerQ,x]]
PositiveIntegerQ[x_]:= And[IntegerQ[x],x>0]
(*--------------------------------------------------------------*)
HermiteFormUpper::argx="One argument required."
HermiteFormUpper::arg1="Argument is a square, non-singular integer matrix."

HermiteFormUpper[]:= (Message[HermiteFormUpper::argx]; Abort[])

HermiteFormUpper[argseq__]:=Module[{args=List[argseq],m,r,h,u},
    If[Length[args]\[NotEqual]1,(Message[HermiteFormUpper::argx];
        Abort[])];
    m=Part[args,1];
    If[Not[
        And[IntegerMatrixQ[m],SSquareMatrixQ[m], 
          Det[m]\[NotEqual] 0]],(Message[HermiteFormUpper::arg1];
        Abort[])];
    
    {h,r}=hnf[Change[m],u];
    Return[{Change[Inverse[u]],Change[h]}]]

HermiteFormLower::argx="One argument required."
HermiteFormLower::arg1="Argument is a square, non-singular integer matrix."

HermiteFormLower[]:= (Message[HermiteFormLower::argx]; Abort[])

HermiteFormLower[argseq__]:=Module[{args=List[argseq],m,r,h,u},
    If[Length[args]\[NotEqual]1,(Message[HermiteFormLower::argx];
        Abort[])];
    
    m=Part[args,1];
    If[Not[
        And[IntegerMatrixQ[m], SSquareMatrixQ[m],
          Det[m]\[NotEqual] 0]],(Message[HermiteFormLower::arg1];
        Abort[])];
    {h,r}=hnf[m,u];
    Return[{Inverse[u],h}]]

IntegerMatrixQ[x_]:=MatrixQ[x,IntegerQ]
(* SquareMatrixQ[x_]:=And[MatrixQ[x],Length[x]\[Equal]Length[First[x]]] *)
(*--------------------------------------------------------------------*)
SetAttributes[hnf,HoldAll](* We pass some arguments by name *)

hnf::usage:="{newmatrix,rank}=hnf[mat] gives a lower triangular Hermite \
Normal Form of the integral matrix mat, together with its rank.\n
  {newmatrix,rank}=hnf[matrix,symbol] also assigns to symbol the matrix trans \
so that trans.mat=newmatrix."

hnf[bb_,uu_:False]:=
  Module[
    {u,m=Length[bb],n=Length[bb[[1]]],i,j,l,pp,bb0,rank=0,dimker=0,k=1,maxk=0,
      ans}
    ,
    pp=0 Range[m];(* 
      pp[[i]] will record how far row i sticks out. 
          We initialize pp as a row of zeroes. 
          Any List of length m would do. *)
    If[uu=!=False
      ,(* If the transformation is requested, 
        put an identity matrix in front, 
        to keep track of the transformation. *)
      bb0=Transpose[Join[IdentityMatrix[m],Transpose[bb]]];
      n=n+m
      ,(* else *)
      bb0=bb
      ];
    For[i=1,i\[LessEqual]m,i++
      ,(* Record how far row i sticks out. *)
      j=n;
      While[(j>0)&&(bb0[[i,j]]\[Equal]0),j--];
      pp[[i]]=j;
      ];
    While[k\[LessEqual]m
      ,(* main loop *)
      maxk=k;(* number of new row that is to be processed *)
      While[(k>dimker+1)&&(pp[[k]]\[LessEqual]pp[[k-1]])
        ,
        swap[k,bb0,pp];(* This subroutine uses extended Euclid *)
        k--;
        ];
      If[pp[[k]]\[Equal]0,dimker++,rank++];
      For[Null,k\[LessEqual]maxk,k++
        ,(* Reduce the columns, as required by Hermite. *)
        For[l=k-1,l>dimker,l--
          ,(* 
            We use the order of reduction that is recommended by LLL or Chou \
Collins *)
          j=pp[[l]];
          u=Quotient[bb0[[k,j]],bb0[[l,j]]];
          bb0[[k]]=bb0[[k]]-u bb0[[l]]
          ]
        ];
      k=maxk+1 (* Proceed to next row, 
        now that part above it is in HNF form. *)
      ];(* end of main loop *)
    If[uu=!=False
      ,(* Extract new matrix, transformation and rank. *)
      n=n-m;
      uu=Take[#,m]&/@bb0;
      ans=Take[#,-n]&/@bb0;
      rank=0;
      For[i=1,i\[LessEqual]m,i++,If[pp[[i]]>m,rank++]]
      ,
      ans=bb0
      ];
    {ans,rank}
    ]
SetAttributes[swap,HoldAll];

swap[kk_,bb_,pp_]:=(* Swaps row kk with row kk-1, 
    using extended GCD if needed. Changes bb and pp. *) 
  Module[
    {i,j,l,t1,t2,e,x,y,k=kk}(* So kk is passed by value. *)
    ,
    If[pp[[k]]<pp[[k-1]]
      ,(* Row k-1 sticks out further. *)
      {bb[[k-1]],bb[[k]]}={bb[[k]],
          bb[[k-1]]};
      {pp[[k-1]],pp[[k]]}={pp[[k]],pp[[k-1]]};
      ,(* They stick out equally. *)
      j=pp[[k]];
      {e,{x,y}}=ExtendedGCD[bb[[k-1,j]],bb[[k,j]]];
      t1=Quotient[bb[[k,j]],e];
      t2=-Quotient[bb[[k-1,j]],e];(* This is a doubly extended GCD, 
        so that Det[{{t1,t2},{x,y}}]\[Equal]1. *)
      {bb[[k-1]],
          bb[[k]]}={{t1,t2},{x,y}}.{bb[[k-1]],bb[[k]]};
      If[(j>0)&&(bb[[k,j]]<0),bb[[k]]=-bb[[k]]];(* 
        Keep pseudo-diagonal positive. *)
      
      k--;(* Look again how far it sticks out. *)
      
      While[(j>0)&&(bb[[k,j]]\[Equal]0),j--];
      pp[[k]]=j;
      If[(j>0)&&(bb[[k,j]]<0),bb[[k]]=-bb[[k]]];(* 
        Keep pseudo-diagonal positive. *)
      ]
    ]
Change[x_]:=Reverse[Reverse /@ x]
(*-----------------------------------------------------------------*)
IwasawaYVariables::argx="One argument required."
IwasawaYVariables::arg1="Argument is a non-singular matrix of CREs in Iwasawa form."

IwasawaYVariables[]:= (Message[IwasawaYVariables::argx]; Abort[])

IwasawaYVariables[argseq__]:=Module[{args=List[argseq],a,i,y},
If[Length[args]\[NotEqual]1,(Message[IwasawaYVariables::argx]; Abort[])];
a=Part[args,1];
If[Not[IwasawaQ[a]],(Message[IwasawaYVariables::arg1]; Abort[])];

   n=Length[a]; y=Table[0,{i,1,n-1}];
    Do[y[[n-i]]=a[[i,i]]/a[[i+1,i+1]],{i,1,n-1}];
    Return[y]]


IwasawaYMatrix::argx="One argument required."
IwasawaYMatrix::arg1="Argument is a non-singular matrix of CREs in Iwasawa form."

IwasawaYMatrix[]:= (Message[IwasawaYMatrix::argx]; Abort[])

IwasawaYMatrix[argseq__]:=Module[{args=List[argseq],a,i,n,y},
   If[Length[args]\[NotEqual]1,(Message[IwasawaYMatrix::argx]; Abort[])];
    a=Part[args,1];
If[Not[IwasawaQ[a]],(Message[IwasawaYMatrix::arg1]; Abort[])];
n=Length[a];
y=Table[0,{i,1,n},{j,1,n}];
    Do[y[[i,i]]=a[[i,i]],{i,1,n}];
    Return[y]]

IwasawaXMatrix::argx="One argument required."
IwasawaXMatrix::arg1="Argument is a non-singular matrix of CREs in Iwasawa form."

IwasawaXMatrix[]:= (Message[IwasawaXMatrix::argx]; Abort[])

IwasawaXMatrix[argseq__]:=Module[{args=List[argseq],a,i,j,n,x},
    If[Length[args]\[NotEqual]1,(Message[IwasawaXMatrix::argx]; Abort[])];
a=Part[args,1];
If[Not[IwasawaQ[a]],(Message[IwasawaXMatrix::arg1]; Abort[])];
n=Length[a];
    x=Table[0,{i,1,n},{j,1,n}];
    Do[x[[i,j]]=a[[i,j]]/a[[j,j]],{i,1,n},{j,i,n}];
    Return[x]]

IwasawaXVariables::argx="One argument required."
IwasawaXVariables::arg1="Argument is a non-singular matrix of CREs in Iwasawa form."

IwasawaXVariables[]:= (Message[IwasawaXVariables::argx]; Abort[])

IwasawaXVariables[argseq__]:=
  Module[{args=List[argseq],a,out={},n},
  If[Length[args]\[NotEqual]1,(Message[IwasawaXVariables::argx]; Abort[])];
  a=Part[args,1];
If[Not[IwasawaQ[a]],(Message[IwasawaXVariables::arg1]; Abort[])];
n=Length[a];
 Do[Do[out=Join[{a[[i,j]]/a[[j,j]]},out],{j,i+1,n}],{i,1,n-1}];
    Return[Reverse[out]]]
(*---------------------------------------------------------------------------------------*)

IFun::argx="Two arguments required."
IFun::arg1="First argument is a list of n-1 CREs."
IFun::arg2="Second argument is a non-singular square matrix of CREs of size n"

IFun[]:= (Message[IFun::argx]; Abort[])

  IFun[argseq__]:= Module[{args=List[argseq],nv,z,n1,n,i,j,b,y},
    If[Length[args]\[NotEqual]2,(Message[IFun::argx]; Abort[])];

    nv=Part[args,1];z=Part[args,2];
    If[Not[ExprListQ[nv]],(Message[IFun::arg1]; Abort[])];
    If[Not[And[ExprMatrixQ[z],NonSingularMatrixQ[z]]],
              (Message[IFun::arg2]; Abort[])];
    
    n1=Length[nv];n=n1+1;
    If[n != Length[z],(Message[IFun::arg2]; Abort[])];
    
    Do[b[i,j]=If[i+j\[LessEqual] n, i j, (n-i)(n-j)], {i,1,n1},{j,1,n1}];
    y=IwasawaYVariables[IwasawaForm[z]];
    Return[Product[y[[i]]^(b[i,j]*Part[nv,j]), {j,1,n1},{i,1,n1}]]]

NonSingularMatrixQ[x_]:= And[MatrixQ[x],Det[x] != 0]
(*-----------------------------------------------------------------------*)             

RamanujanSum::argx="Two arguments required."
RamanujanSum::arg1="First argument is a whole number."
RamanujanSum::arg2="Second argument is a whole number."

RamanujanSum[]:= (Message[RamanujanSum::argx]; Abort[])

RamanujanSum[argseq__]:=Module[{args=List[argseq],g,k,n},
If[Length[args] != 2,(Message[RamanujanSum::argx];Abort[])];
k=Part[args,1];n=Part[args,2];
If[Not[And[IntegerQ[k],k>0]],(Message[RamanujanSum::arg1];Abort[])];
If[Not[And[IntegerQ[n],n>0]],(Message[RamanujanSum::arg2];Abort[])];
g=GCD[k,n];
   Return[ Sum[If[Mod[g,d]\[Equal]0,d*MoebiusMu[k/d],0], {d,1,g}]]]

(*---------------------------------------------------*)

SmithForm::argx="One argument required."
SmithElementaryDivisors::argx="One argument required."
SmithInvariantFactors::argx="One argument required."

SmithForm::arg1="Argument should be a non-singular, square, integer matrix."

SmithForm[]:= (Message[SmithForm::argx]; Abort[])


SmithForm[argseq__]:= Module[{args=List[argseq],a},
If[Length[args] != 1,(Message[SmithForm::argx];Abort[])];
a=Part[args,1];

If[Not[And[IntegerMatrixQ[a],Det[a] != 0]],(Message[SmithForm::arg1]; Abort[])];
Return[ExtendedSmithForm[a]]]

SmithElementaryDivisors[]:= (Message[SmithElementaryDivisors::argx]; Abort[])

SmithElementaryDivisors[argseq__]:= Module[{args=List[argseq],a},
If[Length[args] != 1,(Message[SmithElementaryDivisors::argx];Abort[])];
a=Part[args,1];

If[Not[And[IntegerMatrixQ[a],Det[a] != 0]],(Message[SmithForm::arg1]; Abort[])];
Return[ElementaryDivisors[a]]]

SmithInvariantFactors[]:= (Message[SmithInvariantFactors::argx]; Abort[])

SmithInvariantFactors[argseq__]:= Module[{args=List[argseq],a},
If[Length[args] != 1,(Message[SmithInvariantFactors::argx];Abort[])];
a=Part[args,1];

If[Not[And[IntegerMatrixQ[a],Det[a] != 0]],(Message[SmithForm::arg1]; Abort[])];
Return[InvariantFactors[a]]]

(*----------------------------------------------------------------------------*)
NormalizedExtendedGCD[a_,b_]:=
          If[ Mod[b,a]==0, Return[{Abs[a],{Sign[a],0}}],Return[
          ExtendedGCD[a,b]] ]     

R1[i_,u_,A_]:= If[Depth[A]==4, Return[Map[ R1[i,u,#] &, A]],Return[
      Array[Which[ #==i, u A[[i]] , True, A[[#]] ] &, Length[A] ] ]]    ;
C1[i_,u_,A_]:=If[Depth[A]==4, Return[Map[C1[i,u,#] &, A]], 
            Return[Transpose[R1[i,u, Transpose[A]]]]];
R2[c_,i_,j_,A_]:=If[Depth[A]==4, 
               Return[Map[R2[c,i,j,#] &, A]],
              Return[Array[Which[ # == j, A[[j]]+c A[[i]], 
             True, A[[#  ]]] &, Length[A]]] ];
C2[c_, i_,j_,A_]:=If[Depth[A]==4, 
                Return[Map[C2[c,i,j,#] &, A]],
                Return[ Transpose[R2[c,i,j,Transpose[A] ]]]] ;
R3[i_,j_, A_]:= If[Depth[A]==4, Return[Map[R3[i,j, #]&, A]],
              Return[
               Array[Which[ # == i, A[[j]], # ==j, 
                A[[i]], True, A[[#  ]]] &, Length[A] ]] ];              
C3[i_,j_, A_]:= If[Depth[A]==4, Return[Map[C3[i,j, #]&, A]],
                Transpose[R3[i,j,Transpose[A]]]];
R4[i_,j_,d_,a_,b_,x_,y_,A_] := If[Depth[A] ==4, 
                     Return[Map[R4[i,j,d,a,b,x,y,#] &,A]],
                     Return[Array[Which[ # == i, 
                     x A[[i]]+ y A[[j]],
                     # == j, -(b/d) A[[i]] +
                     (a/d) A[[j]], True, A[[#]] ] &,
                     Length[A]] ] ];
                            
C4[i_,j_,d_,a_,b_,x_,y_,A_] := If[Depth[A] ==4, 
                     Return[Map[C4[i,j,d,a,b,x,y, #] &,A]],
                     Transpose[R4[i,j,d,a,b,x,y, Transpose[A]]]];
             
 EliminateCol[A_]:=
            Module[{nA=A,ExtGCD},
            
            Do[If[nA[[m,j]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[m,j]] ];
            nA=C4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[m,j]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], nA]
            ), (Continue[ ])],
             {j,m+1, cols}];
            Return[nA]];   


            
 EliminateRow[A_]:=
            Module[{nA=A, ExtGCD},
            
            Do[If[nA[[j,m]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[j,m]] ];
            nA=R4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[j,m]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], nA]
            ), (Continue[ ])],
             {j,m+1, rows}];
            Return[nA]];            
                       
            
            
ThirdStep[A_]:=Module[{nA=A,InterchangeQ=False},
                (Do[ ( Do[If[ Mod[ A[[i, j]] , A[[m,m]] ] == 0, 
                (Continue[ ]),
                     (nA=R2[1,i, m,nA];
                     
                     nA=EliminateCol[nA];
                     
                     InterchangeQ=True;
                      Break[ ]) ],
                     {j, m+1, Dimensions[A][[2]]}]),
                 
                  {i,m+1,Dimensions[A][[1]]}];
                  
                  If[InterchangeQ , 
                  nA= EliminateRow[nA]
                  ];
                  Return[nA])];             

SmithNormalForm[A_] := Block[ {nA=A,
                        rows, cols,MinEntry=0, Loc,
                        ZeroMatrixQ,NonZeroRowQ=True,NonZeroColQ=
                        True},
          (
                
          (* Some Initialization*)
          {rows, cols}= Dimensions[A];
           
          
                
           (* Now we start the big loop on m*)
           Do[(
           
          (*We first find the min  non-zero value in the remaining
          principal submatrix*)     
          Loc={m,m};ZeroMatrixQ=True;
          Do[ Which[nA[[i,j]] != 0 && ZeroMatrixQ, 
                       (MinEntry=Abs[ nA[[i,j]] ];Loc={i,j};
                        ZeroMatrixQ=False), 
                       nA[[i,j]] != 0 && Abs[ nA[[i,j]] ]< MinEntry,
                       (MinEntry=Abs[ nA[[i,j]] ];
                        Loc={i,j} ), 
                       True, Continue[ ] ], 
                      {i,m,rows},{j, m,cols} 
                      ]; 
           If[ZeroMatrixQ,  Break[ ]];  (*If the remaining
           principal submatrix has all zeros we are done*)        
           
           
           (*First Step: Put min value in (m,m) position *)
           
           If[Loc[[1]] != m, nA=R3[m,Loc[[1]],nA] ];
           If[Loc[[2]] != m, nA=C3[m, Loc[[2]],nA]];
           
           If[nA[[m,m]] < 0, nA=R1[m,-1,nA] ];
           
           
           (*Second Step:  Eliminate row and column *)
           
           If[ m< Min[rows, cols], (
           
           
           
           While[(Do[(
           If[ nA[[i,m]] !=0, (NonZeroRowQ=True;Break[ ]), (Continue[ ])]),
           {i,m+1, rows}];
           If[!NonZeroRowQ,
           Do[(If[ nA[[m,j]] !=0, (NonZeroColQ=True;
           Break[ ]), 
           (Continue[ ])]),
           {j,m+1,cols}]];
           NonZeroRowQ  || NonZeroColQ),
           
           
           If[NonZeroRowQ,nA = EliminateRow[nA]];
           
           
           nA = EliminateCol[ nA];
           
           NonZeroRowQ=False;NonZeroColQ=False];
           
           (*Third Step: Make sure m,m entry divides the rest*)
           nA = ThirdStep[nA];
           
       ),
            (
            Which[rows < cols, 
            
           nA = EliminateCol[nA], 
           rows> cols, 
           nA = EliminateRow[nA]
           , rows==cols, (
           Break [ ] )])]
                
                
                ),{m,1, Min[rows, cols]}];
                
                Return[nA]
                ) 
                ];               

 ExtendedEliminateCol[P_, A_, Q_]:=
            Module[{nA=A,nP=P,tA=A, nQ=Q,ExtGCD},
            
            Do[If[nA[[m,j]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[m,j]] ];
            {nA, nQ}=C4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[m,j]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], {nA, nQ}]
            ), (Continue[ ])],
             {j,m+1, cols}];
            Return[{nP, nA, nQ}]];   


            
 ExtendedEliminateRow[P_, A_, Q_]:=
            Module[{nA=A,nP=P,nQ=Q, ExtGCD},
            
            Do[If[nA[[j,m]] != 0,
            (ExtGCD=NormalizedExtendedGCD[nA[[m,m]], nA[[j,m]] ];
            {nA, nP}=R4[m,j,ExtGCD[[1]],nA[[m,m]], nA[[j,m]],
            ExtGCD[[2,1]], ExtGCD[[2,2]], {nA, nP}]
            ), (Continue[ ])],
             {j,m+1, rows}];
            Return[{nP, nA, nQ}]];            
                       

            
            
ExtendedThirdStep[P_, A_, Q_]:=
                Module[{nA=A,nP=P, nQ=Q,InterchangeQ=False},
                (Do[ ( Do[If[ Mod[ A[[i, j]] , A[[m,m]] ] == 0, 
                (Continue[ ]),
                     ({nA,nP}=R2[1,i, m,{nA,nP}];
                     
                     {nP, nA, nQ}=ExtendedEliminateCol[nP,nA, nQ];
                     
                     InterchangeQ=True;
                      Break[ ]) ],
                     {j, m+1, Dimensions[A][[2]]}]),
                 
                  {i,m+1,Dimensions[A][[1]]}];
                  
                  If[InterchangeQ , 
                  {nP,nA,nQ}= ExtendedEliminateRow[nP,nA,nQ]
                  ];
                  Return[{nP,nA,nQ}])];             

ExtendedSmithForm[A_] := Block[ {nP, nA=A, nQ, 
                        rows, cols,MinEntry=0, Loc,
                        ZeroMatrixQ,NonZeroRowQ=True,NonZeroColQ=
                        True},
          (
                
          (* Some Initialization*)
          {rows, cols}= Dimensions[A];
           nP= IdentityMatrix[rows];
           nQ= IdentityMatrix[cols];
           
          
                
           (* Now we start the big loop on m*)
           Do[(
           
          (*We first find the min  non-zero value in the remaining
          principal submatrix*)     
          Loc={m,m};ZeroMatrixQ=True;
          Do[ Which[nA[[i,j]] != 0 && ZeroMatrixQ, 
                       (MinEntry=Abs[ nA[[i,j]] ];Loc={i,j};
                        ZeroMatrixQ=False), 
                       nA[[i,j]] != 0 && Abs[ nA[[i,j]] ]< MinEntry,
                       (MinEntry=Abs[ nA[[i,j]] ];
                        Loc={i,j} ), 
                       True, Continue[ ] ], 
                      {i,m,rows},{j, m,cols} 
                      ]; 
           If[ZeroMatrixQ,  Break[ ]];  (*If the remaining
           principal submatrix has all zeros we are done*)        
           
           
           (*First Step: Put min value in (m,m) position *)
           
           If[Loc[[1]] != m,{nA,nP}=R3[m,Loc[[1]],{nA, nP}] ];
           If[Loc[[2]] != m, {nA,nQ}=C3[m, Loc[[2]],{nA,nQ}]];
           
           If[nA[[m,m]] < 0, {nA,nP}=R1[m,-1,{nA,nP}] ];
           
           (*Second Step:  Eliminate row and column *)
           
           If[ m< Min[rows, cols], (
           
           
           
           While[(Do[(
           If[ nA[[i,m]] !=0, (NonZeroRowQ=True;Break[ ]), (Continue[ ])]),
           {i,m+1, rows}];
           If[!NonZeroRowQ,
           Do[(If[ nA[[m,j]] !=0, (NonZeroColQ=True;
           Break[ ]), 
           (Continue[ ])]),
           {j,m+1,cols}]];
           NonZeroRowQ  || NonZeroColQ),
           
           
           If[NonZeroRowQ,{nP, nA, nQ} = ExtendedEliminateRow[nP, nA, nQ]];
           
           
           {nP, nA, nQ} = ExtendedEliminateCol[nP, nA, nQ];
           
           NonZeroRowQ=False;NonZeroColQ=False];
           
           (*Third Step: Make sure m,m entry divides the rest*)
           {nP, nA, nQ} = ExtendedThirdStep[nP, nA, nQ]
       ),
            (
            Which[rows < cols, 
            
           {nP, nA, nQ} = ExtendedEliminateCol[nP, nA, nQ], 
           rows> cols, 
           {nP, nA, nQ} = ExtendedEliminateRow[nP, nA, nQ]
           , rows==cols, (
           Break [ ] )])]
                
                
                ),{m,1, Min[rows, cols]}];
                
                Return[{Inverse[nP], nA, Inverse[nQ]}]
                ) 
                ];
               

  InvariantFactors[A_]:=Module[{B,del={1},t,dims=Dimensions[A]},
  If[Max[dims]< 7,

(Do[(t=Apply[GCD, Union[Flatten[Minors[A, k]]]];
AppendTo[del,t]; j=k-1; If[t==0, (Break[ ])] ), {k,1,
Min[dims]}];
Return[Table[ If[k<j+3,del[[k]]/del[[k-1]],0],
{k,2,Min[dims]+1}]]),  
  
(B=SmithNormalForm[A];Return[Table[B[[i,i]],{i,1,Min[dims] }]])] ] ;


ElementaryDivisors[A_]:=Module[{nA=A,T},(T=Sort[Flatten[
 FactorInteger[InvariantFactors[nA]],1]];
  Table [ T[[i,1]]^T[[i,2]], {i, 1, Length[T]}]) ]             


(*-----------------------------------------------------------*)

WeylGroup::argx="One argument required."
WeylGroup::arg1="First argument is a whole number."


WeylGroup[]:= (Message[WeylGroup::argx]; Abort[])

WeylGroup[argseq__]:=
  Module[{args=List[argseq],elements={}, prev, all,newelement,newmatrix, len,
      i,j,k,m,n},
    If[Length[args]\[NotEqual]1,(Message[WeylGroup::argx];Abort[])];
    n=Part[args,1];
    If[Not[And[IntegerQ[n],n>0]],(Message[WeylGroup::arg1];Abort[])];
    If[n\[Equal]1, Return[{{{1}}}]];
    If[n\[Equal]2,Return[{{{1,0},{0,1}},{{0,1},{1,0}}}]];
    all=WeylGroup[n-1];len=Length[all];
    Do[
      {Do[{Do[newelement[i,j]=0, {i,1,n},{j,1,n}], newelement[1,j]=1,
            prev=Part[all,m], 
            {Do[
                Do[newelement[k+1,i]=
                    If[i<j,prev[[k,i]],If[i\[Equal]j,0,prev[[k,i-1]]]],{k,1,
                    n-1}],{i,1,n}]},newmatrix=ArrayToMatrix[n,n,newelement] ,
            elements=Join[{newmatrix},elements]},{m,1,len}]},{j,1,n}];
    Return[Reverse[elements]]]


(*************************************************************)
SpecialWeylGroup::argx = "One argument required."
SpecialWeylGroup::arg1 = "First argument is a whole number."

SpecialWeylGroup[]:= (Message[SpecialWeylGroup::argx]; Abort[])

SpecialWeylGroup[argseq__] := 
  Module[{args = List[argseq], elements = {}, prev, all, newelement, 
      newmatrix, len, i, j, k, m, n, w}, 
    If[Length[args] != 1, (Message[SpecialWeylGroup::argx]; 
        Abort[])];
    n = Part[args, 1];
    If[Not[And[IntegerQ[n], n > 0]], (Message[SpecialWeylGroup::arg1]; 
        Abort[])];
    Return[Select[WholeWeylGroup[n], DetOne]]]
(***************************************************************)

WholeWeylGroup[n_] := 
  Module[{elements = {}, prev, all, newelement, newmatrix, len, i, j, k, m},
    If[n == 1, Return[{{{1}}, {{-1}}}]];
    all = WholeWeylGroup[n - 1]; len = Length[all];
    Do[{Do[{Do[newelement[i, j] = 0, {i, 1, n}, {j, 1, n}], 
            prev = Part[all, m],
            
            newelement[1, j] = 
              1, {Do[Do[
                  newelement[k + 1, i] = 
                    If[i < j, prev[[k, i]], 
                      If[i == j, 0, prev[[k, i - 1]]]], {k, 1, n - 1}], {i, 1,
                   n}]}, newmatrix = ArrayToMatrix[n, n, newelement], 
            elements = Join[{newmatrix}, elements]}, {m, 1, len}]}, {j, 1, 
        n}];
    Do[{Do[{Do[newelement[i, j] = 0, {i, 1, n}, {j, 1, n}], 
            prev = Part[all, m],
            
            newelement[1, j] = -1, {Do[
                Do[newelement[k + 1, i] = 
                    If[i < j, prev[[k, i]], 
                      If[i == j, 0, prev[[k, i - 1]]]], {k, 1, n - 1}], {i, 1,
                   n}]}, newmatrix = ArrayToMatrix[n, n, newelement], 
            elements = Join[{newmatrix}, elements]}, {m, 1, len}]}, {j, 1, 
        n}];
    Return[Reverse[elements]]]

DetOne[m_]:= Det[m]==1

(*=====================================================================*)

a2m[n_,m_,a_]:=ArrayToMatrix[n,m,a]

MakeXMatrix::argx = "Two arguments required";
MakeXMatrix::arg1 = "First argument us an integer >= 1";
MakeXMatrix::arg2 = "Second argument is a string";

MakeXMatrix[]:= (Message[MakeXMatrix::argx]; Abort[])

MakeXMatrix[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, uarray}, 
    If[Length[args] != 2, (Message[MakeXMatrix::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeXMatrix::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeXMatrix::arg2]; Abort[])];
    u = Symbol[ustr];
    Do[uarray[i, i] = 1, {i, 1, n}];
    Do[uarray[i, j] = u[i, j], {i, 1, n - 1}, {j, i + 1, n}];
    Do[uarray[i, j] = 0, {i, 2, n}, {j, 1, i - 1}];
    Return[a2m[n, n, uarray]]]
(*---------------------------------------------------------------*)
MakeMatrix::argx = "Three arguments required";
MakeMatrix::arg2 = "Second argument us an integer >= 1";
MakeMatrix::arg1 = "First argument is a string";
MakeMatrix::arg3 = "Third argument us an integer >= 1";

MakeMatrix[]:=(Message[MakeMatrix::argx];Abort[])

MakeMatrix[argseq__] := 
  Module[{args = List[argseq], n, m, ustr, u, i, j, uarray}, 
    If[Length[args] != 3, (Message[MakeMatrix::argx];
        Abort[])];
    ustr = First[args]; m = Part[args, 2]; n=Part[args,3];
    If[! And[IntegerQ[m], m >= 1], (Message[MakeMatrix::arg2];
        Abort[])];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeMatrix::arg3];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeMatrix::arg1]; Abort[])];

    u = Symbol[ustr];
    Do[uarray[i, j] = u[i, j], {i, 1, m}, {j, 1, n}];
    Return[a2m[m, n, uarray]]]
(*----------------------------------------------------------------*)

MakeBlockMatrix::argx = "One argument required.";
MakeBlockMatrix::arg1 = "Argument is a matrix of matrices.";
MakeBlockMatrix::arg2 = "Each submatrix must have expr elements.";
MakeBlockMatrix::arg3 = 
    "The submatrices must have compatible row and column numbers.";

MakeBlockMatrix[]:=(Message[MakeBlockMatrix::argx];Abort[])

MakeBlockMatrix[argseq__] := 
  Module[{arg = {argseq}, m, n, ms, ns, mb, nb, i, j,ii,ir,jj,jr, A, B, barray},
    If[Length[arg] != 1, (Message[MakeBlockMatrix::argx]; Abort[])];
    A = First[arg];
    If[Or[! MatrixLikeQ[A], MemberQ[Map[MatrixQ, Apply[Join, A]], False]],
      (Message[MakeBlockMatrix::arg1]; Abort[])];

If[MemberQ[Map[ExprMatrixQ, Apply[Join, A]], False],
      (Message[MakeBlockMatrix::arg2]; Abort[])];
    If[! CompatibleBlockMatrixQ[A], (Message[MakeBlockMatrix::arg3]; 
        Abort[])];
    
    m = Length[A]; n = Length[First[A]];
    mb = Apply[Plus, ms = Map[Length, Map[First, A]]];
    nb = Apply[Plus, ns = Map[Length, Map[First, First[A]]]];

    
    Do[{{ii, ir} = Bmatincrement[i, ms], {jj, jr} = Bmatincrement[j, ns], 
        barray[i, j] = (A[[ir, jr]])[[ii, jj]]}, {i, 1, mb}, {j, 1, nb}];
    Return[a2m[mb,nb,barray]]]


MatrixLikeQ[x_]:= And[ListQ[x],!MemberQ[Map[ListQ,x],False],
                    Length[First[x]]==First[Union[Map[Length,x]]]]

Bmatincrement[i_, ms_] := Module[{s = Part[ms, 1], j = 1},
    If[i <= s, Return[{i, 1}]];
    While[s < i, {j = j + 1, s = s + Part[ms, j]}];
    Return[{i - s + Part[ms, j], j}]]


CompatibleBlockMatrixQ[A_] := 
  Module[{m, n, ms, ns, mb, nb, ans = True,ai,bj,i,j}, 
    m = Length[A]; 
    n = Length[First[A]];
    mb = Apply[Plus, ms = Map[Length, Map[First, A]]];
    nb = Apply[Plus, ns = Map[Length, Map[First, First[A]]]];
    Do[If[Or[NRows[A[[i,j]]] != ms[[i]],NColumns[A[[i,j]]] != ns[[j]]],
         {ans=False,Break[]}],{i,1,m},{j,1,n}];
    Return[ans]]
(*--------------------------------------------------------------*)
MakeXVariables::argx = "Two arguments required";
MakeXVariables::arg1 = "First argument us an integer >= 1";
MakeXVariables::arg2 = "Second argument is a string";

MakeXVariables[]:= (Message[MakeXVariables::argx]; Abort[])

MakeXVariables[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, ans={}}, 
    If[Length[args] != 2, (Message[MakeXVariables::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeXVariables::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeXVariables::arg2]; Abort[])];
    u = Symbol[ustr];
    Do[ans = Join[{u[i, j]},ans], {i, 1, n - 1}, {j, i + 1, n}];
    Return[Reverse[ans]]]
    
(*----------------------------------------------------------------*)

MakeYMatrix::argx = "Two arguments required";
MakeYMatrix::arg1 = "First argument us an integer >= 1";
MakeYMatrix::arg2 = "Second argument is a string";

MakeYMatrix[]:= (Message[MakeYMatrix::argx]; Abort[])

MakeYMatrix[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, next=1}, 
    If[Length[args] != 2, (Message[MakeYMatrix::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeYMatrix::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeYMatrix::arg2]; Abort[])];
    u = Symbol[ustr];
    Do[uarray[i,j]=0, {i,1,n},{j,1,n}];
    Do[uarray[i, i] = 1, {i, 1, n}];
    Do[{next=next*u[j],uarray[n-j,n-j] = next}, {j,1,n-1}];
    Return[a2m[n, n, uarray]]]
(*-----------------------------------------------------------------*)

MakeYVariables::argx = "Two arguments required";
MakeYVariables::arg1 = "First argument us an integer >= 1";
MakeYVariables::arg2 = "Second argument is a string";

MakeYVariables[]:= (Message[MakeYVariables::argx]; Abort[])

MakeYVariables[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, ans={}}, 
    If[Length[args] != 2, (Message[MakeYVariables::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeYVariables::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeYVariables::arg2]; Abort[])];

    u = Symbol[ustr];
    Do[ans=Join[{u[i]},ans], {i, 1, n-1}];
    Return[Reverse[ans]]]
(*-----------------------------------------------------------*)

MakeZMatrix::argx = "Three arguments required";
MakeZMatrix::arg1 = "First argument us an integer >= 1";
MakeZMatrix::arg2 = "Second argument is a string";
MakeZMatrix::arg3 = "third argument is a string";

MakeZMatrix[]:= (Message[MakeZMatrix::argx]; Abort[])

MakeZMatrix[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, uarray,y,ystr,next=1}, 
    If[Length[args] != 3, (Message[MakeZMatrix::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2]; ystr=Part[args,3];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeZMatrix::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeZMatrix::arg2]; Abort[])];
    If[! StringQ[ystr], (Message[MakeZMatrix::arg3]; Abort[])];

    u = Symbol[ustr]; y=Symbol[ystr];
    Return[MakeXMatrix[n,ustr].MakeYMatrix[n,ystr]]]

(*-----------------------------------------------------------*)

MakeZVariables::argx = "Two arguments required";
MakeZVariables::arg1 = "First argument us an integer >= 1";
MakeZVariables::arg2 = "Second argument is a string";
MakeZVariables::arg3 = "Third argument is a string";

MakeZVariables[]:= (Message[MakeZVariables::argx]; Abort[])

MakeZVariables[argseq__] := 
  Module[{args = List[argseq], n, ustr, u, i, j, ansx={},ansy={},y}, 
    If[Length[args] != 3, (Message[MakeZVariables::argx];
        Abort[])];
    n = First[args]; ustr = Part[args, 2]; ystr=Part[args,3];
    If[! And[IntegerQ[n], n >= 1], (Message[MakeZVariables::arg1];
        Abort[])];
    If[! StringQ[ustr], (Message[MakeZVariables::arg2]; Abort[])];
    If[! StringQ[ystr], (Message[MakeZVariables::arg3]; Abort[])];

    u = Symbol[ustr];y=Symbol[ystr];
    Do[ansx = Join[{u[i, j]},ansx], {i, 1, n - 1}, {j, i + 1, n}];
     Do[ansy=Join[{y[i]},ansy], {i, 1, n-1}];
    Return[Join[Reverse[ansx], Reverse[ansy]]]]
(*=== === === === === === === === === === === === === === === === === =*)


(* ParabolicQ
    LanglandsForm
    LanglandsIFun
    HeckeEigenvalue
    SchurPolynomial
    HeckeRepresentatives
    HeckeOperator
    SpecialWeylGroup
    HeckeMultiplicativeSplit
       MakeAllLists PositiveIntegerQ, IntegerListQ , GetCall, 
       WholeWeylGroup, DetOne, PrimePower (private)
  *)

hrs = HeckeRepresentatives;
hop = HeckeOperator;
spl = SchurPolynomial;
paq = ParabolicQ;
llf = LanglandsForm;
lif = LanglandsIFun;
hev = HeckeEigenvalue;
mal = MakeAllLists;
swg = SpecialWeylGroup;



(*** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ********************)
  
PositiveIntegerQ[n_] := And[IntegerQ[n], n > 0]
IntegerListQ[x_] := And[ListQ[x], {True} == Union[Map[PositiveIntegerQ, x]]]

ParabolicQ::argx = "Two arguments required."
ParabolicQ::arg1 = "The first argument is a square matrix of CREs of size 2 or more."
ParabolicQ::arg2 = "The second argument is a list of positive integers."

ParabolicQ[]:= (Message[ParabolicQ::argx]; Abort[])

ParabolicQ[argseq__] := 
  Module[{args = List[argseq], a, p, n, r, ans = True, l = 0, h = 0, li, hi, 
      lj, hj, lblocks, hblocks} ,
    If[Length[args] != 2, (Message[ParabolicQ::argx];
        Abort[])];
    a = First[args]; p = Part[args, 2];
    If[Not[ExprSquareMatrixQ[a]], (Message[ParabolicQ::arg1]; Abort[])];
    If[Not[IntegerListQ[p]], (Message[ParabolicQ::arg2]; Abort[])];
    r = Length[p]; n = Length[a];
    Do[If[
        Not[And[IntegerQ[Part[p, j]], Part[p, j] > 0]], {ans = False, 
          Break[]}], {j, 1, r}];
    If[ans == False, Return[ans]];
    If[r == 1, Return[False]];
    If[Not[n == Apply[Plus, p]], Return[False]];
    Do[{l = h + 1, h = h + Part[p, j], 
        If[Det[BlockMatrix[a, {l, h}, {l, h}]] == 0, {ans = False, 
            Break[]}]}, {j, 1, r}];
    If[ans == False, Return[ans]];
    li = 1 + p[[1]]; hi = p[[1]] + p[[2]]; lj = 1; hj = p[[1]];
    lblocks = Table[0, {i, 1, r}];
    hblocks = Table[0, {i, 1, r}];
    lblocks[[1]] = 1;
    Do[lblocks[[i]] = lblocks[[i - 1]] + p[[i - 1]], {i, 2, r}];
    hblocks[[1]] = p[[1]];
    Do[hblocks[[i]] = hblocks[[i - 1]] + p[[i]], {i, 2, r}];
    Do[Do[
        If[Not[SameQ[a[[u, v]], 0]], {ans = False, Break[]}], {u, 
          lblocks[[i]], hblocks[[i]]}, {v, lblocks[[j]], hblocks[[j]]}], {i, 
        2, r}, {j, 1, i - 1}];
    
Return[ans]]

(*** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ***********)

LanglandsForm::argx = "Two arguments required.";
LanglandsForm::arg1 = "The first argument is a nonsingular matrix of CREs.";
LanglandsForm::arg2 = "The second argument is a list of positive integers 
with sum the size of the matrix.";
LanglandsForm::gen = "The first argument is in the parabolic subgroup 
specified by the second argument.";

LanglandsForm[]:= (Message[LanglandsForm::argx]; Abort[])

LanglandsForm[argseq__] := 
  Module[{args = List[argseq], a, p, u, uu, vv, c, cc, cv, s, ss, sv, n, r, 
      l = 0, h = 0, b, d, low, up, ni}, 
    If[Length[args] != 2, (Message[LanglandsForm::argx];
        Abort[])];
    a = First[args]; p = Part[args, 2];
    If[Not[
        And[ExprSquareMatrixQ[a], Det[a] != 0]], (Message[LanglandsForm::arg1];
        Abort[])];
    If[Not[
        And[IntegerListQ[p], Length[a] == Apply[Plus, p]]], (Message[
          LanglandsForm::arg2]; Abort[])];
    If[Not[ParabolicQ[a, p]], (Message[LanglandsForm::gen];
        Abort[])];
    n = Length[a]; r = Length[p];
    low[1] = 1; up[1] = Part[p, 1];
    Do[{ni = Part[p, i], low[i] = up[i - 1] + 1, up[i] = low[i] + ni - 1}, {i,
         2, r}];
    Do[{b[i, i] = BlockMatrix[a, {low[i], up[i]}, {low[i], up[i]}], 
        c[i, i] = Inverse[b[i, i]], d[i, i] = Det[b[i, i]]}, {i, 1, r}];
    Do[If[
        i > j, {b[i, j] = ZeroMatrix[Part[p, i], Part[p, j]], 
          u[i, j] = b[i, j], c[i, j] = b[i, j], s[i, j] = b[i, j]}], {i, 1, 
        r}, {j, 1, r}];
    Do[If[
        i < j, {b[i, j] = ZeroMatrix[Part[p, i], Part[p, j]], 
          c[i, j] = b[i, j], s[i, j] = b[i, j]}], {i, 1, r}, {j, 1, r}];
    Do[u[i, i] = IdentityMatrix[Part[p, i]], {i, 1, r}];
    Do[If[i < j, 
        u[i, j] = 
          BlockMatrix[a, {low[i], up[i]}, {low[j], up[j]}].c[j, j]], {i, 1, 
        r}, {j, 1, r}];
    Do[{uu[i] = u[i, 1], 
        Do[uu[i] = MatrixJoinHorizontal[uu[i], u[i, j]], {j, 2, r}]}, {i, 1, 
        r}];
    vv = uu[1];
    Do[vv = MatrixJoinVertical[vv, uu[i]], {i, 2, r}];
    Do[c[i, i] = (Abs[d[i, i]]^(1/Part[p, i]))*IdentityMatrix[Part[p, i]], {i,
         1, r}];
    Do[If[Not[Equal[i, j]], c[i, j] = ZeroMatrix[Part[p, i], Part[p, j]]], {i,
         1, r}, {j, 1, r}];
    Do[{cc[i] = c[i, 1], 
        Do[cc[i] = MatrixJoinHorizontal[cc[i], c[i, j]], {j, 2, r}]}, {i, 1, 
        r}];
    cv = cc[1];
    Do[cv = MatrixJoinVertical[cv, cc[i]], {i, 2, r}];
    Do[s[i, i] = b[i, i]/Abs[d[i, i]]^(1/Part[p, i]), {i, 1, r}];
    Do[If[Not[Equal[i, j]], s[i, j] = ZeroMatrix[Part[p, i], Part[p, j]]], {i,
         1, r}, {j, 1, r}];
    Do[{ss[i] = s[i, 1], Do[ss[i] = mjh[ss[i], s[i, j]], {j, 2, r}]}, {i, 1, 
        r}];
    sv = ss[1];
    Do[sv = mjv[sv, ss[i]], {i, 2, r}];
    Return[{vv, cv, sv}]]
(*** ** ** ** ** ** ** ******************)
LanglandsIFun::argx = "Three arguments required."
LanglandsIFun::arg1 = "The first argument is a nonsingular matrix of CREs."
LanglandsIFun::arg3 = "The third argument is a list of CREs \
of the same length as the second."
LanglandsIFun::arg2 = "The second argument is a list of positive integers with \
sum the size of the matrix."
LanglandsIFun::gen = "The first argument is in the parabolic subgroup specified \
by the second argument."
LanglandsIFun::gen2 = "The dot product of the second and third arguments should \
be zero."

LanglandsIFun[]:= (Message[LanglandsIFun::argx]; Abort[])

LanglandsIFun[argseq__] := 
  Module[{args = List[argseq], g, s, p, u, c, ss, m, mm, pr = 1, r, low, up, 
      ni},
    If[Length[args] != 3, (Message[LanglandsIFun::argx];
        Abort[])];
    g = First[args]; s = Part[args, 3]; p = Part[args, 2];
    If[Not[And[ExprSquareMatrixQ[g], Det[g] != 0]], (Message[LanglandsIFun::arg1];
         Abort[])];
   
    If[Not[
        And[IntegerListQ[p], Length[g] == Apply[Plus, p]]], (Message[
          LanglandsIFun::arg2]; Abort[])];
     If[Not[And[ListQ[s],Length[s]==Length[p]]], 
           (Message[LanglandsIFun::arg3]; Abort[])];
     If[Not[Apply[Plus,s*p]==0], (Message[LanglandsIFun::gen2]; Abort[])];

    If[Not[ParabolicQ[g,p]], (Message[LanglandsIFun::gen]; Abort[])];
    
    {u, c, ss} = LanglandsForm[g, p];
    r = Length[p];
    mm = c.ss;
    low[1] = 1; up[1] = Part[p, 1];
    Do[{ni = Part[p, i], low[i] = up[i - 1] + 1, up[i] = low[i] + ni - 1}, {i,
         2, r}]; 
    Do[m[j] = BlockMatrix[mm, {low[j], up[j]}, {low[j], up[j]}], {j, 1, r}];
    Do[pr = pr*Abs[Det[m[j]]]^s[[j]], {j, 1, r}];
    Return[Simplify[pr]]]


(*** ** ** ** ** ** ** ** ** ** ** ** ****************)

HeckeEigenvalue::argx = "Three arguments required."
HeckeEigenvalue::arg1 = "The first argument is a natural number."
HeckeEigenvalue::arg2 = "The second argument is an integer > 1."
HeckeEigenvalue::arg3 = "The third argument is a string."

HeckeEigenvalue[]:= (Message[HeckeEigenvalue::argx]; Abort[])


HeckeEigenvalue[argseq__] := 
  Module[{args = List[argseq], m, n, astr, lambdam = 1, k, facs, NoFacts, i, 
      kth, p, r, x, pterm, a, ab}, 
    If[Length[args] != 3, (Message[HeckeEigenvalue::argx];
        Abort[])];
    m = First[args]; n = Part[args, 2]; astr = Part[args, 3];
    If[Not[PositiveIntegerQ[m]], (Message[HeckeEigenvalue::arg1];
        Abort[])];
    If[Not[IntegerTwoOrMoreQ[n]], (Message[HeckeEigenvalue::arg2];
        Abort[])];
    If[Not[StringQ[astr]], (Message[HeckeEigenvalue::arg3];
        Abort[])];
    If[m == 1, Return[lambdam]];
    a = Symbol[astr];
    facs = FactorInteger[m]; NoFacts = Length[facs];
    Do[{kth = Part[facs, i]; p = First[kth]; k = Part[kth, 2];
        pterm = 
          If[k == 1, GetCall[a, n, p, 1], 
            Coefficient[(1 + 
                  Sum[x^j Sum[(-1)^(r + 1)GetCall[a, n, p, r]x^(r - 1), {r, 1,
                             n}]^j, {j, 1, k}]), x, k]]; 
        lambdam = lambdam*pterm}, {i, 1, NoFacts}];
    Return[Simplify[lambdam]]]

GetCall[f_, n_, p_, r_] := 
  Switch[r, n, 1, 1, Apply[f, Join[{p}, Table[1, {i, 1, n - 2}]]], n - 1, 
    Apply[f, Join[Table[1, {i, 1, n - 2}], {p}]], _,
    Apply[f, 
      Join[Table[1, {i, 1, r - 1}], {p}, Table[1, {i, r + 1, n - 1}]]]]
(*************************************************)
HeckeMultiplicativeSplit::argx = "one argument required."
HeckeMultiplicativeSplit::arg1 = "The argument is a list of natural numbers."


HeckeMultiplicativeSplit[]:=(Message[HeckeMultiplicativeSplit::argx];Abort[])

HeckeMultiplicativeSplit[arglist__] := 
  Module[{args = List[arglist], m, ps, p, ans = {}, ansj = {}, i, j},
    If[Length[args] != 1, (Message[HeckeMultiplicativeSplit::argx];
        Abort[])];
    m = First[args]; 
    If[Not[IntegerListQ[m]], (Message[HeckeMultiplicativeSplit::arg1]; 
        Abort[])];

    ps = Union[Map[First, Apply[Join, Map[FactorInteger, Select[m, # > 1 &]]]]];
    If[ps == {}, Return[m]];
    Do[{p = Part[ps, j], ansj = {}, 
        Do[ansj = Join[{PrimePower[p, Part[m, i]]}, ansj], {i, 1, Length[m]}],
         ans = Join[{{p,Reverse[ansj]}}, ans]}, {j, 1, Length[ps]}];
    Return[Reverse[ans]]]


PrimePower[p_, n_] := If[Mod[n, p] != 0, 0, 1 + PrimePower[p, n/p]]
(*----------------------------------------------------------------------*)

ProductSet2[A_, B_] := 
  Module[{i, j}, 
    Return[Apply[Join, 
        Table[{Part[A, i], Part[B, j]}, {i, 1, Length[A]}, {j, 1, 
            Length[B]}]]]]

ProductSet[argseq__] := Module[{args={argseq}},
    Return[If[Length[args] == 2, ProductSet2[Part[args, 1], Part[args, 2]],
        Map[ApplyPrepend, 
          ProductSet2[First[args], Apply[ProductSet, Rest[args]]]]]]]

ApplyPrepend[x_] := Join[{First[x]}, First[Rest[x]]]

CompatibleDivs[m_, ds_] := Module[{ans = {}, j},
    Do[If[Mod[m, Part[ds, j]] == 0, ans = Join[{Part[ds, j]}, ans]], {j, 1, 
        Length[ds]}];
    Return[Reverse[ans]]]

(*----------------------------------------------------------------------------*)

HeckeCoefficientSum::usage = 
    "natural no,list of natural nos, string->expression.";
HeckeCoefficientSum::argx = "Three arguments required.";
HeckeCoefficientSum::arg1 = "First argument is a natural number."; 
HeckeCoefficientSum::arg2 = "Second argument is a list of natural numbers."; 
HeckeCoefficientSum::arg3 = "Third argument is a string.";

HeckeCoefficientSum[]:=(Message[HeckeCoefficientSum::argx];Abort[])

HeckeCoefficientSum[arglist__] := 
  Module[{args = List[arglist], m, ms, astr, n1, i, Di, d, dlist = {}, s = 0, 
      j, a, arg},
    If[Length[args] != 3, (Message[HeckeCoefficientSum::argx]; 
        Abort[])];
    m = First[args]; ms = Part[args, 2]; astr = Part[args, 3];
If[Not[PositiveIntegerQ[m]], (Message[HeckeCoefficientSum::arg1]; 
        Abort[])];
    If[Not[PositiveIntegerListQ[ms]], (Message[HeckeCoefficientSum::arg2]; 
        Abort[])]; 
    If[Not[StringQ[astr]], (Message[HeckeCoefficientSum::arg3]; 
        Abort[])];
    
    n1 = Length[ms];
    a = Symbol[astr];
    Do[{Di = CompatibleDivs[m, Divisors[Part[ms, i]]], 
        dlist = Join[{Di}, dlist]}, {i, 1, n1}];
    Ds = Apply[ProductSet,Reverse[dlist]];
    Do[{cs = Part[Ds, j], d = Apply[Times, cs], 
        If[Mod[m, d] == 0, {cn = m/d, arg = {cn Part[ms, 1]/Part[cs, 1]}, 
            Do[arg = Join[arg, {Part[ms, i]Part[cs, i - 1]/Part[cs, i]}], {i, 
                2, n1}], s = s + Apply[a, arg]}]}, {j, 1, Length[Ds]}];
    Return[s]]

PositiveIntegerQ[n_] := And[IntegerQ[n], n > 0]

(*--------------------------------------------------------------------------------*)

PositiveIntegerListQ[x_] := 
  Apply[And, Join[{ListQ[x]}, Map[PositiveIntegerQ, x]]]


HeckePowerSum::argx = "Three arguments required.";
HeckePowerSum::arg1 = "First argument is a non-negative integer."; 
HeckePowerSum::arg2 = "Second argument is a list of non-negative integers."; 
HeckePowerSum::arg3 = "Third argument is a string.";

HeckePowerSum[]:= (Message[HeckePowerSum::argx];Abort[])

HeckePowerSum[argseq__] := 
  Module[{args = List[argseq], m, ms, astr, n1, i, Di, d, dlist = {}, s = 0, 
      j, a, arg, k},
    If[Length[args] != 3, (Message[HeckePowerSum::argx]; Abort[])];
    m = First[args]; ms = Part[args, 2]; astr = Part[args, 3];
    If[Not[NonNegativeIntegerQ[m]], (Message[HeckePowerSum::arg1]; Abort[])];
    If[Not[NonNegativeIntegerListQ[ms]], (Message[HeckePowerSum::arg2]; 
        Abort[])]; 
    If[Not[StringQ[astr]], (Message[HeckePowerSum::arg3]; Abort[])];
    
    n1 = Length[ms];
    a = Symbol[astr];
    Do[{Di = Table[k, {k, 0, Min[Part[ms, i], m]}], 
        dlist = Join[{Di}, dlist]}, {i, 1, n1}];
    Ds = Apply[ProductSet,Reverse[dlist]]; 
    Do[{cs = Part[Ds, j], d = Apply[Plus, cs], 
        If[d <= m, {cn = m - d, arg = {cn + Part[ms, 1] - Part[cs, 1]}, 
            Do[arg = 
                Join[arg, {Part[ms, i] + Part[cs, i - 1] - Part[cs, i]}], {i, 
                2, n1}], s = s + Apply[a, arg]}]}, {j, 1, Length[Ds]}];
    Return[s]]


NonNegativeIntegerQ[x_]:= And[IntegerQ[x],x>-1]
NonNegativeIntegerListQ[xs_]:= 
        And[ListQ[xs], !MemberQ[Map[NonNegativeIntegerQ,xs],False]]

(********************************************************************)
SchurPolynomial::argx = "Two arguments required."
SchurPolynomial::arg1 = "The first argument is a list of two or more CREs."
SchurPolynomial::arg2 = "The second argument is a list of CREs, \
being one less than that in the first argument."

SchurPolynomial[]:= (Message[SchurPolynomial::argx]; Abort[])

SchurPolynomial[argseq__] := 
  Module[{args = List[argseq], k, x, den, num, mat, n, i, j, index},
    If[Length[args] != 2, (Message[SchurPolynomial::argx];
        Abort[])];
    x = First[args]; k = Part[args, 2];
    If[Not[And[ExprListQ[x], Length[x] > 1]], (Message[SchurPolynomial::arg1]; 
        Abort[])];
    If[Not[
        And[ExprListQ[k], Length[k] + 1 == Length[x]]], (Message[
          SchurPolynomial::arg2]; Abort[])];
    
    n = Length[x];
    den = Product[x[[i]] - x[[j]], {i, 1, n}, {j, i + 1, n}];
    If[den == 0, Return[False]];
    mat = Table[1, {i, 1, n}, {j, 1, n}];
    index = Table[0, {i, 1, n}];
    Do[index[[i + 1]] = index[[i]] + k[[i]], {i, 1, n - 1}];
    Do[mat[[i, j]] = x[[j]]^(index[[n - i + 1]] + n - i), {i, 1, n - 1}, {j, 
        1, n}];
    num = Det[mat];
    Return[Simplify[num/den]]]


(*** ** ** ** ** ** ** ** ** ** **  ** ** ********************)
  
HeckeOperator::argx = "Three arguments required."
HeckeOperator::arg1 = "The first argument should be a natural number."
HeckeOperator::arg2 = "The second argument should be a square matrix of CREs."
HeckeOperator::arg3 = "The third argument should be a string being the name \
of a function of a square matrix or a symbol."

HeckeOperator[]:= (Message[HeckeOperator::argx]; Abort[])

HeckeOperator[argseq__] := 
  Module[{args = List[argseq], N, z, f, fs, coef, n, hr, ans = 0},
    If[Length[args] != 3, (Message[HeckeOperator::argx];
        Abort[])];
    N = First[args]; z = Part[args, 2]; f = Part[args, 3];
    If[Not[PositiveIntegerQ[N]], (Message[HeckeOperator::arg1]; 
        Abort[])];
    If[Not[ExprSquareMatrixQ[z]], (Message[HeckeOperator::arg2]; 
        Abort[])];
    If[Not[StringQ[f]], (Message[HeckeOperator::arg3]; Abort[])];
    
    n = Length[z]; fs = Symbol[f];
    coef = 1/N^((n - 1)/2);
    hr = hrs[N, n];
    Do[ans = ans + fs[Part[hr, j].z], {j, 1, Length[hr]}];
    Return[Simplify[ans]]]



(*** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** ** ** ************)

HeckeRepresentatives[N_,n_] := 
  Module[{SN, i, j, k, old = {}, new = {}, bot, 
      side, divs, d, sided, hrd, zs, hr},
    
    If[n == 1, Return[{{{N}}}]];
    zs = Table[0, {i, 1, n - 1}];
    bot = {Join[zs, {1}]}; 
    side = Transpose[{zs}];
    hr = hrs[N, n - 1];
    Do[old = Join[{mjv[mjh[Part[hr, j], side], bot]}, old], {j, 1, 
        Length[hr]}];
    divs = Complement[Divisors[N], {1}];
    Do[{d = Part[divs, j]; hrd = hrs[N/d, n - 1]; bot = {Join[zs, {d}]};
        sided = Map[Transpose, Map[List, MakeAllLists[d, n - 1]]];
        Do[
          Do[{side = Part[sided, k]; 
              new = Join[{mjv[mjh[Part[hrd, i], side], bot]}, new]}, {k, 1, 
              Length[sided]}],
          {i, 1, Length[hrd]}]}, {j, 1, Length[divs]}];
    Return[Join[new, old]]]


(***********************************************************)



MakeAllLists[d_,n_] := 
  Module[{ans = {}, ind, i, j},
    If[n == 1, Return[Table[{i}, {i, 0, d - 1}]]];
    ind = MakeAllLists[d, n - 1]; 
    Do[Do[ans = Join[{Join[{i}, Part[ind, j]]}, ans], {j, 1, 
          Length[ind]}], {i, 0, d - 1}];
    Return[Reverse[ans]]]

(*= === === === === === === === === === === = === ===*)

(* 
KloostermanSumQ::argx = "Four arguments required."
KloostermanSumQ::arg1 = "The first argument should be a list of integers."
KloostermanSumQ::arg2 = "The second argument should be a list of integers."
KloostermanSumQ::arg3 = "The third argument should be a list of integers."
KloostermanSumQ::arg4 = "The fourth argument should be a permutation matrix."
KloostermanSumQ::argen = "If the matrix is nxn the first three arguments should\
have length n-1."
KloostermanSumQ::abort = "The function KloostermanSumQ is currently under development."

KloostermanSumQ[]:=(Message[KloostermanSumQ::argx];Abort[])

KloostermanSumQ[argseq__] := 
  Module[{args = List[argseq], psi, psip, c, w, n},
  If[True, 
       (Message[KloostermanSumQ::abort];
       Abort[])]; 
   
    If[Length[args] != 3, (Message[KloostermanSumQ::argx];
        Abort[])];
    N = First[args]; z = Part[args, 2]; f = Part[args, 3];
    If[Not[PositiveIntegerQ[N]], (Message[KloostermanSumQ::arg1]; 
        Abort[])];
    If[Not[SSquareMatrixQ[z]], (Message[KloostermanSumQ::arg2]; 
        Abort[])];
    If[Not[StringQ[f]], (Message[KloostermanSumQ::arg3]; Abort[])];
   
  
    n = Length[z]; fs = Symbol[f];
    coef = 1/N^((n - 1)/2);
    hr = hrs[N, n];
    Do[ans = ans + fs[Part[hr, j].z], {j, 1, Length[hr]}];
    Return[Simplify[ans]]]                                      *)
(*====================================================*)
(*------------------------------------------------------*)

SubscriptedForm::argx = "One argument required.";
SubscriptedForm::arg = "Input expression is not suitable for subscripting.";

SubscriptedForm[] := (Message[SubscriptedForm::argx]; Abort[]);

SubscriptedForm[argseq__] := Module[{args = {argseq}},
    If[Length[args] != 1, (Message[SubscriptedForm::argx]; Abort[])];
    If[MemberQ[{Integrate, Sum}, Head[args[[1]]]], 
                 (Message[SubscriptedForm::arg];Abort[])];
    Return[subit[First[args]]]]

subit[x_Integer | x_Rational | x_Real | x_Complex | x_Symbol] := x
subit[x_ + y_] := subit[x] + subit[y]
subit[x_ * y_] := subit[x] * subit[y]
subit[x_ / y_] := subit[x] / subit[y]
subit[-x_] := -subit[x]
subit[x_ - y_] := subit[x] - subit[y]
subit[Rule[x_, y_]] := Rule[subit[x], subit[y]]
subit[x_List] := Map[subit, x]
subit[a_Symbol[n_Integer]] := DisplayForm[Subscript[a, n]]
subit[a_Symbol[n_Integer, seq__Integer]] := 
  DisplayForm[Subscript[a, n, Sequence[seq]]]
subit[a_Symbol[x_]] := a[subit[x]]
subit[a_Symbol[x_, y__]] := a[subit[x], Sequence[subit[y]]]
subit[Derivative[x__][f_][y__]] := 
      Derivative[x][f][DeList[Sequence[subit[{y}]]]]
subit[Integrate[f_,lim__]]:= Integrate[subit[f],Delist[Sequence[subit[{lim}]]]]
subit[x_] := x
DeList[{x__}]:= x
(*------------------------------------------------------------------------------------*)\

creq[_Integer | _Rational | _Real | _Complex | _Symbol] := True
creq[x_ + y_] := creq[x] && creq[y]
creq[x_ * y_] := creq[x] && creq[y]
creq[x_ / y_] := creq[x] && creq[y]
creq[-x_] := creq[x]
creq[x_ - y_] := creq[x] && creq[y]
creq[{}] := False
creq[List[__]] := False
creq[""] := False
creq[String[_]] := False
creq[_Symbol[s_]] := creq[s]
creq[_Symbol[s_, S__]] := And[ creq[s], Apply[And, Map[creq, {S}]]]
creq[_] := False

(*------------------------------------------------------------------------------------*)\


KloostermanSum2::argx = "Four arguments required.";
KloostermanSum2::arg1 = 
    "The first argument should be a list of one integer.";
KloostermanSum2::arg2 = 
    "The second argument should be a list of one integer.";
KloostermanSum2::arg3 = 
    "The third argument should be a list of one integer.";
KloostermanSum2::arg4 = 
    "The fourth argument should be a 2x2 permutation matrix.";
KloostermanSum2::argen = 
    "If the matrix is 2x2 the first three arguments should have length 1.";

KloostermanSum2[] := (Message[KloostermanSum2::argx]; Abort[]);


KloostermanSum2[argseq__] := 
  Module[{args = List[argseq], psi, psip, c, c1, w, u, v, n, s = 0},
    If[Length[args] != 4, (Message[KloostermanSum2::argx];
        Return[$Failed])];
    psi = First[args]; psip = Part[args, 2]; c = Part[args, 3]; 
    w = Part[args, 4];
    
    If[Not[GeneralIntegerListQ[psi]], (Message[KloostermanSum2::arg1];
        Abort[])];
    If[Not[GeneralIntegerListQ[psip]], (Message[KloostermanSum2::arg2];
        Abort[])];
    If[Not[GeneralIntegerListQ[c]], (Message[KloostermanSum2::arg3]; 
        Abort[])];
    If[Not[PermutationMatrixQ[w]], (Message[KloostermanSum2::arg4]; 
        Abort[])];
    If[Not[
        And[TwoOrMoreQ[Length[w]], 
          Length[w] == 
            Union[{Length[psi], Length[psip], Length[c]}]]], (Message[
          KloostermanSum2::argen]; Abort[])];
    
    n = Length[w]; c1 = First[c]; u = First[psi]; v = First[psip];
    If[Not[n == 2], (Message[KloostermanSum2::argen]; Abort[])];
    If[SameQ[w, {{1, 0}, {0, 1}}], Return[1]];
    If[n == 2, 
      Do[If[GCD[a, c1] == 1, 
          s = s + Exp[
                2Pi I((u a + v First[Part[ExtendedGCD[a, c1], 2]])/c1)]], {a, 
          1, c1}]];
    Return[s]]


PermutationMatrixQ[w_] := 
  And[MatrixQ[w], Length[w] == Length[First[w]], 
    Module[{i, t}, t = Table[i, {i, 1, Length[w]}];
      Return[t == Union[t.w]]]]

(*---------------------------------------------------------------------*)

SymbolQ[x_] := SameQ[ Symbol, Head[x]]

KloostermanCompatibility::argx = "Five arguments required.";
KloostermanCompatibility::arg1 = "The first argument should be a symbol.";
KloostermanCompatibility::arg2 = "The second argument should be a symbol.";
KloostermanCompatibility::arg3 = "The third argument should be a symbol.";
KloostermanCompatibility::arg4 = 
    "The fourth argument should be a permutation matrix.";
KloostermanCompatibility::arg5 = "The fifth argument should be a symbol.";


KloostermanCompatibility[] := (Message[KloostermanCompatibility::argx]; 
      Abort[]);

KloostermanCompatibility[argseq__] := 
  Module[{args = {argseq}, n, um, uu, coefs, uuu3u, b1, cm, c, cc, v, t1, t2, 
      chars, divs, minors, psip, psi, w}, 
    If[Length[args] != 5, (Message[KloostermanCompatibility::argx]; 
        Abort[])];
    psi = First[args]; psip = Part[args, 2]; cc = Part[args, 3];
    w = Part[args, 4]; v = Part[args, 5];
    
    If[Not[SymbolQ[psi]], (Message[KloostermanCompatibility::arg1]; 
        Abort[])];
    If[Not[SymbolQ[psip]], (Message[KloostermanCompatibility::arg2]; 
        Abort[])];
    If[Not[SymbolQ[cc]], (Message[KloostermanCompatibility::arg3]; Abort[])];
    If[Not[PermutationMatrixQ[w]], (Message[KloostermanCompatibility::arg4]; 
        Abort[])];
    If[Not[SymbolQ[v]], (Message[KloostermanCompatibility::arg5]; Abort[])];
    
    n = Length[w];
    t1 = Table[psi[i], {i, 1, n - 1}];
    t2 = Table[psip[i], {i, 1, n - 1}];
    c = Table[cc[i], {i, 1, n - 1}];
    
    chars = getcharacterrelns[t1, t2, c, w];
    divs = KloostermanDivisibility2[cc, w, v];
    minors = getminors[w, v, c];
    Return[{chars, divs, minors}]]

CharacterApply[th_, m_] := Sum[th[[i]]*m[[i, i + 1]], {i, 1, Length[th]}]

getcharacterrelns[t1_, t2_, c_, w_] := 
  Module[{args = {argseq}, n, um, uu, coefs, uuu3u, b1, cm}, 
    cm = CVectorToMatrix[c];
    n = Length[t1] + 1;
    um = MakeXMatrix[n, "u"]; 
    uu = UpperTriangleToMatrix[MatrixUpperTriangle[Transpose[w].um.w]];
    coefs = DeleteCases[Apply[Join, uu], _Integer];
     b1 = cm.w.uu.Transpose[w].Inverse[cm]; 
     Return[
      DeleteCases[
        Coefficient[CharacterApply[t1, b1] - CharacterApply[t2, uu], coefs], 
        0]]]

CVectorToMatrix[cv_] := Module[{c, n},
    n = Length[cv] + 1;
    c = ZeroMatrix[n, n];
    c[[1, 1]] = 1/cv[[n - 1]];
    c[[n, n]] = cv[[1]];
    Do[ c[[j, j]] = cv[[n - j + 1]]/cv[[n - j]], {j, 2, n - 1}];
    Return[c]]

(*----------------------------------------------------*)


FindRules[eqns_] := Module[{rules = {}},
    Do[(If[Depth[eqns[[j]]] == 2, 
          rules = Join[{Rule[eqns[[j]], 0]}, rules]]), {j, 1, Length[eqns]}];
    Return[rules]]

productset2[A_, B_] := 
  Module[{i, j}, 
    Return[Apply[Join, 
        Table[{Part[A, i], Part[B, j]}, {i, 1, Length[A]}, {j, 1, 
            Length[B]}]]]]

productset[argseq__] := Module[{args = {argseq}},
    If[Length[args] == 1, Return[Map[List, First[args]]]]; 
    Return[If[Length[args] == 2, productset2[Part[args, 1], Part[args, 2]], 
        Map[applyprepend, 
          productset2[First[args], Apply[productset, Rest[args]]]]]]]


applyprepend[x_] := Join[{First[x]}, First[Rest[x]]]


KloostermanDivisibility2[c_, w_, v_] := 
 Module[{n = Length[w], cirules, relations, cs, divs = {}, vs, 
   lexirules, reduced, zerored = {}, zeroredd, zeroreddd, redd, 
   fourterm = {}, terms, termsij = {}, term, vsij, remainingvs, 
   divisor, dividend, x},
  cs = Table[c[i], {i, 1, n - 1}];
  cirules = getciminors[w, v, cs];
  vs = Map[Function[x, First[First[x]]], cirules];
  (*ci|1*)Do[(If[BiggestQ[vs[[j]], n], 
       divs = Join[{{c[Length[vs[[j]]]], 1}}, divs]] ;;;), {j, 1, 
    Length[vs]}];
  (*lex only v[{}]=0*)
  lexirules = getleximinors[w, v, cs];
  relations = 
   PluckerRelationsInner[v, Table[i, {i, 1, n}]];
  reduced = DeleteCases[relations /. Join[lexirules, cirules], 0];
  (*will have some -1 type duplicates which we ignore*)(*go through \
and find the 0's from the reduced list*)Do[(term = reduced[[j]];
    If[FreeQ[term, Plus], 
      zerored = 
       Join[DeleteCases[Table[D[term, c[i]], {i, 1, n - 1}], 0], 
        zerored]] ;; ;;), {j, 1, Length[reduced]}];
  (*remove minus signs*)
  zeroredd = Union[zerored /. Times[-1, x_] -> x];
  zeroreddd = Map[Function[x, Rule[x, 0]], zeroredd];
  redd = DeleteCases[Union[reduced /. zeroreddd], 0];
  (*now collect all the four term relations of the form ci v_j+ /-cj \
v_k=0 have c<v*)
  Do[term = redd[[j]];
   If[And[FourTermQ[term, cs], Not[SameQ[GetFourTerms[term], False]]],
     fourterm = Join[{GetFourTerms[term]}, fourterm]], {j, 1, 
    Length[redd]}];
  
  (*now go for divisibility for 1 term through max number by checking \
divisibility for ci|cj,all possible distinct pairs i,
  j except when ci=+/-1*)
  Do[If[And[Not[i == j], Not[MemberQ[divs, {c[i], 1}]]], (vsij = {};
     Do[(term = fourterm[[k]];
       If[And[SameQ[term[[1]], c[i]], SameQ[term[[3]], c[j]]], 
        vsij = Union[Join[{term[[4]]}, vsij]]];
       If[And[SameQ[term[[3]], c[i]], SameQ[term[[1]], c[j]]], 
        vsij = Union[Join[{term[[2]]}, vsij]]]), {k, 1, 
       Length[fourterm]}]; 
     remainingvs = 
      Complement[
       getvector[i - 1, n, v] /. 
        Join[lexirules, zeroreddd, cirules], {0, c[i], -c[i]}];
     If[SSubsetQ[remainingvs, vsij], 
      divs = Join[{{c[i], c[j]}}, divs]])], {i, 1, n - 1}, {j, 1, 
    n - 1}];
  Return[Union[divs]]]



BiggestQ[ind_, n_] := 
  Module[{i = Length[ind]}, 
    Return[Apply[Plus, ind] == Sum[j, {j, n - i + 1, n}]]]

FourTermQ[term_, cs_] := 
  And[Length[term] == 2, 
    Apply[Or, Map[Function[x, Not[FreeQ[term, x]]], cs]]]

 GetFourTerms[c1_*v1_ + c2_*v2_] := {c1, v1, c2, v2}  (* order c < v *)
 GetFourTerms[-c1_*v1_ + c2_*v2_] := {c1, v1, c2, v2}
 GetFourTerms[c1_*v1_ - c2_*v2_] := {c1, v1, c2, v2}
 GetFourTerms[-c1_*v1_ - c2_*v2_] := {c1, v1, c2, v2}
 GetFourTerms[x_]:= False

SSubsetQ[a_, b_] := Union[a, b] == b

(*--------------------------------------------------------*)


KloostermanSum::argx = "Four arguments required.";
KloostermanSum::arg1 = "The first argument should be a list of integers.";
KloostermanSum::arg2 = "The second argument should be a list of integers.";
KloostermanSum::arg3 = "The third argument should be a list of integers.";
KloostermanSum::arg4 = 
    "The fourth argument should be a nxn permutation matrix.";
KloostermanSum::argen = 
    "If the matrix is nxn the first three arguments should have length n-1.";
KloostermanSum::char = 
    "The character related and ci arguments are incompatible. Run \
KloostermanCompatibility.";
KloostermanSum::div = 
    "Divisibility relations between the ci's fail. Run \
KloostermanCompatibility.";


KloostermanSum[] := (Message[KloostermanSum::argx]; Abort[]);


KloostermanSum[argseq__] := Module[{args = {argseq}, t1, t2, c, w},
    If[Length[args] != 4, (Message[KloostermanSum::argx]; Abort[])];
    t1 = First[args]; t2 = Part[args, 2]; c = Part[args, 3]; 
    w = Part[args, 4];
    
    If[Not[GeneralIntegerListQ[t1]], (Message[KloostermanSum::arg1]; 
        Abort[])];
    If[Not[GeneralIntegerListQ[t2]], (Message[KloostermanSum::arg2]; 
        Abort[])];
    If[Not[GeneralIntegerListQ[c]], (Message[KloostermanSum::arg3]; 
        Abort[])];
    If[Not[PermutationMatrixQ[w]], (Message[KloostermanSum::arg4]; Abort[])];
    If[Not[
        And[gln`Private`TwoOrMoreQ[Length[w]], {Length[w] - 1} == 
            Union[{Length[t1], Length[t2], Length[c]}]]], (Message[
          KloostermanSum::argen]; Abort[])];
    
    Return[KloostermanSumInner2[t1, t2, c, w]]]




KloostermanSumInner2[t1_, t2_, c_, w_] := 
 Module[{n, sum = 0, b1, b2, x, y, rules, ciconst, crconst = {}, 
   reduced, xmat, xinv, xxmat, ymat, yymat, yyymat, gmat, ggmat, cmat,
    gggmat, ggggmat, rhs, lhs, asubs = {}, newlhs, rmat, brules = {}, 
   eqns = {}, varsx = {}, varsy = {}, singlerules = {}, 
   singleeqns = {}, eqnrules, yeqns = {}, klscomp, lrrules = {}, cc, 
   klsdivis, rij, gij, eqncheck, eqnj, eqj, vbles, solution = {}, 
   eqnss, symmins, symminss, prod, prodd, proddi, proddd, prodddi, 
   reducedbigprod = {}, rbp, x1, y1, c1, w1, start, minors, minorss, 
   a, v, chars, divs, cis, tt1, tt2, ccc, vv, charrules = {}, di, s, 
   relations, bigprod, flatsymmins, bpj, reducedj}, n = Length[c] + 1;
  (*check arguments*)(*check compatibility:
  character and divisibility*)
  {chars, divs, cis} = KloostermanCompatibility[tt1, tt2, ccc, w, vv];
  IfPrint["klc=", {chars, divs, cis}];
  
  Do[charrules = 
    Join[{tt1[i] -> t1[[i]], tt2[i] -> t2[[i]], ccc[i] -> c[[i]]}, 
     charrules], {i, 1, n - 1}];
  IfPrint["charrules = ", charrules];
  If[Not[SameQ[DeleteCases[chars /. charrules, 0], {}]], (Message[
     KloostermanSum::char]; Abort[])];
  Do[(di = divs[[i]];
    If[Not[Mod[di[[2]], di[[1]]] == 0], (Message[KloostermanSum::div];
      Abort[])]), {i, 1, Length[divs]}];
  If[IdentityMatrixQ[w], Return[1]];
  (*initialize*)
  ciconst = getminors[w, v, c];
  IfPrint["ci constraints: ", ciconst];
  symmins = Table[getvector[k, n, v], {k, 0, n - 2}];
  IfPrint["minors = ", symmins];
  symminss = symmins /. ciconst;
  IfPrint["reduced minors = ", symminss];
  Do[(s = symminss[[i]]; 
    prod[i] = Apply[productset, expandvec[s, c[[i]]]];
    prodd[i] = {};
    Do[(If[Apply[GCD, prod[i][[j]]] == 1, 
       prodd[i] = Join[{prod[i][[j]]}, prodd[i]]]), {j, 1, 
      Length[prod[i]]}]), {i, 1, n - 1}];
  IfPrint["# candidates = ", 
   Map[Length, Table[prodd[i], {i, 1, n - 1}]]];
  (*now reduce using plucker*)
  relations = PluckerRelationsInner[v, Table[i, {i, 1, n}]];
  reduced = DeleteCases[relations /. ciconst, 0];
  IfPrint["reduced = ", reduced];
  bigprod = Apply[productset, Table[prodd[i], {i, 1, n - 1}]];
  IfPrint["# bigprod = ", Length[bigprod], " first = ", 
   First[bigprod]];
  If[reduced == {}, 
   reducedbigprod = bigprod, (flatsymmins = Flatten[symmins];
    Do[(bpj = Flatten[bigprod[[j]]];
       rules = MapThread[Rule, {flatsymmins, bpj}];
       reducedj = reduced /. rules;
       If[Union[reducedj] == {0}, 
        reducedbigprod = Join[{bpj}, reducedbigprod]]), {j, 1, 
       Length[bigprod]}] ;;;)];
  IfPrint[" # reduced bigprod = ", Length[reducedbigprod]];
  IfPrint[" reduced bigprod = ", reducedbigprod];
  rbp = reducedbigprod;
  (*produce the sum*)
  Do[(Label[start]; minors = rbp[[j]];
    minorss = RestoreMinors[minors, n];
    IfPrint["minorss=", minorss];
    a = PluckerInverse[minorss];
    IfPrint["a,j = ", {a, j}];
    If[a == False, (j = j + 1; Goto[start])];
    {x1, c1, w1, y1} = BruhatForm[a];
    IfPrint[MatrixForm[c1], MatrixForm[cmat]];
    IfPrint[MatrixForm[w1], MatrixForm[w]];
    b1 = Sum[t1[[i]]*x1[[i, i + 1]], {i, 1, n - 1}];
    b2 = Sum[t2[[i]]*y1[[i, i + 1]], {i, 1, n - 1}];
    sum = sum + Exp[2 Pi I (b1 + b2)]), {j, 1, Length[rbp]}];
  Return[sum]]







IdentityMatrixQ[m_]:= Module[{v=Table[i,{i,1,Length[m]}]}, Return[SameQ[{v},{v}.m]]]



PresentQ[form_, vbles_] := Module[{ans = {}},
    Do[If[Not[FreeQ[form, vbles[[i]]]], ans = Join[{vbles[[i]]}, ans]], {i, 1,
         Length[vbles]}];
    If[ans == {}, Return[False]];
    Return[Reverse[ans]]]

IfPrint[argseq__] := If[$verbose, Print[argseq]]

RestoreMinors[M_, n_] := Module[{minori = {}, ms = {}, start, stop = 0},
    If[n == 2, Return[M]];
    Do[(start = stop + 1; stop = stop + Binomial[n, i]; 
        minori = Part[M, Range[start, stop]];
             ms = Join[{minori}, ms]), {i, 1, n - 1}];
    Return[Reverse[ms]]]

(*=== ==== === === === === === === === === === =*)


KloostermanBruhatCell::argx = "Five arguments required.";
KloostermanBruhatCell::arg1 = "The first argument should be a symbol.";
KloostermanBruhatCell::arg2 = "The second argument should be a symbol.";
KloostermanBruhatCell::arg3 = "The third argument should be a symbol.";
KloostermanBruhatCell::arg4 = 
    "The fourth argument should be a nxn permutation matrix.";
KloostermanBruhatCell::arg5 = "The fifth argument should be a symbol.";


KloostermanBruhatCell[] := (Message[KloostermanBruhatCell::argx]; Abort[]);


KloostermanBruhatCell[argseq__] := Module[{args = {argseq}, a, x, c, w, y},
    If[Not[Length[args] == 5], (Message[KloostermanBruhatCell::argx]; 
        Abort[])];
    a = First[args]; x = Part[args, 2]; c = Part[args, 3]; w = Part[args, 4]; 
    y = Part[args, 5];
    
    If[Not[SymbolQ[a]], (Message[KloostermanBruhatCell::arg1]; Abort[])];
    If[Not[SymbolQ[x]], (Message[KloostermanBruhatCell::arg2]; Abort[])];
    If[Not[SymbolQ[c]], (Message[KloostermanBruhatCell::arg3]; Abort[])];
    If[Not[PermutationMatrixQ[w]], (Message[KloostermanBruhatCell::arg4]; 
        Abort[])];
    If[Not[SymbolQ[y]], (Message[KloostermanBruhatCell::arg5]; Abort[])];
    
    Return[KloostermanBruhatCellInner[a, x, c, w, y]]]

KloostermanBruhatCellInner[a_, x_, c_, w_, y_] := 
  Module[{n, sum = 0, b1, b2, rules, ciconst, crconst = {}, reduced, xmat, 
      xinv, xxmat, ymat, yymat, yyymat, gmat, ggmat, cmat, gggmat, ggggmat, 
      rhs, lhs, asubs = {}, newlhs, rmat, brules = {}, eqns = {}, varsx = {}, 
      varsy = {}, singlerules = {}, singleeqns = {}, eqnrules, yeqns = {}, 
      klscomp, lrrules = {}, cc, klsdivis, rij, gij, eqncheck, eqnj, eqj, 
      vbles, solution = {}, eqnss, symmins, symminss, prod, prodd, proddi, 
      proddd, prodddi, reducedbigprod = {}, rbp, x1, y1, c1, w1, start, 
      minors, minorss, solij, varsij, solutionn,
      singlevars = {}, remainingvars, eqnsss},
    (* remove extra local variables *)
    
    (* initialize *)
    n = Length[w];
    xmat = MakeXMatrix[n, "x0"];
    ymat = MakeXMatrix[n, ToString[y]];
    gmat = MakeMatrix[ToString[a], n, n];
    xinv = MakeXMatrix[n, ToString[x]];
    cmat = CVectorToMatrix[Table[c[i], {i, 1, n - 1}]];
    IfPrint["permutation = ", MatrixForm[w]];
    
    (* reduce y *)
    rmat =  w.ymat.Transpose[w];
    Do[(If[Not[SameQ[rmat[[i, j]] , 0]], 
          brules = Join[{Rule[rmat[[i, j]], 0]}, brules]]), {i, 1, n}, {j, 
        i + 1, n}];
    yymat = ymat /. brules;
    IfPrint["brules = ", brules];
    
    (* form a = x.c.w.y *)
    rhs = xinv.cmat.w.yymat;
    
    IfPrint["xs and ys = ", MatrixForm[subit[xinv]], 
      MatrixForm[subit[yymat]]];
    IfPrint["eqns to solve = ", MatrixForm[subit[gmat]], 
      MatrixForm[subit[rhs]]];
     
    vbles = 
      DeleteCases[
        Apply[Join, 
          Join[ Table[xinv[[i, j]], {i, 1, n - 1}, {j, i + 1, n}], 
            Table[yymat[[i, j]], {i, 1, n - 1}, {j, i + 1, n}]]], 0];
    
    (* solve for xij, yij = f(aij) *)
    Do[(rij = rhs[[i, j]]; gij = gmat[[i, j]];
        If[
          And[Head[rij] == Times, 
            Length[varsij = PresentQ[rij, vbles]] == 1],
              (singlerules = 
              Join[solij = First[Solve[gij == rij, First[varsij]]], 
                singlerules];
            singlevars = Join[{First[solij][[1]]}, singlevars];
                
            singleeqns = 
              Join[{First[varsij] == First[solij][[2]]}, 
                singleeqns])];;;), {i, 1, n}, {j, 1, n}];
    
    IfPrint["singlerules =", singlerules, "singlevars = ", singlevars];
    remainingvars = Complement[vbles, singlevars];
    IfPrint["remainingvars = ", remainingvars];
    IfPrint["vbles = ", vbles]; 
    
    (* While[Not[remainingvars == {}],
          
          Do[(rij = rhs[[i, j]] /. eqns; gij = gmat[[i, j]]; 
              varsij = PresentQ[rij, remainingvars]; 
              If[And[Head[rij] == Plus, 
                  Length[varsij] == 1], (eqns = 
                    Join[First[Solve[gij == rij, First[varsij]]], eqns]; 
                  remainingvars = 
                    Complement[remainingvars, varsij])];;;;;), {i, 1, n}, {j, 
              1, n}]]; *)
    
    Do[(rij = rhs[[i, j]]; gij = gmat[[i, j]]; 
          varsij = PresentQ[rij, remainingvars]; 
          Do[If[And[Head[rij] == Plus, Not[SameQ[varsij, False]], 
                LinearlyPresentQ[rij, varsij[[k]], vbles]], (eqns = 
                  Join[First[Solve[gij == rij, varsij[[k]]]], eqns])], {k, 1, 
              Length[varsij]}];;;;;), {i, 1, n}, {j, 1, n}];
      
      IfPrint["eqns = ", eqns];
    
    eqnss = eqns /. singlerules;
    eqnsss = 
      Map[Function[r, Rule[r[[1]], Simplify[r[[2]] /. eqnss]]], eqnss];
    IfPrint["eqnss = ", eqnss];
    
    solutionn = Sort[Join[brules, singlerules, eqnsss]];
    Return[solutionn]]

LinearlyPresentQ[expr_, var_, vars_] := 
  Not[SameQ[D[expr, var] /. Map[Function[v, Rule[v, 0]], vars], 0]]

(*=== === === === === === === === = === === === === === === =*)


MatrixMinor[a_, is_, js_] := 
  Module[{b, l = Length[is]}, b = Table[0, {i, 1, l}, {j, 1, l}];
    Do[b\[LeftDoubleBracket]i, j\[RightDoubleBracket] = 
        a\[LeftDoubleBracket]is\[LeftDoubleBracket]i\[RightDoubleBracket], 
          js\[LeftDoubleBracket]
            j\[RightDoubleBracket]\[RightDoubleBracket], {i, 1, l}, {j, 1, 
        l}];
    Return[Det[b]]]


VectorWedge[v_, e_, m_, n_] := 
  Module[{sum = 0, subs, ind, minor, term}, 
    subs = KKSubsets[Table[i, {i, 1, n}], m];
    Do[ind = subs\[LeftDoubleBracket]j\[RightDoubleBracket]; minor = v[ind];
      term = minor W[e /@ ind]; sum = sum + term, {j, 1, Length[subs]}]; 
    Return[sum]]

(*-------------------------------------------------------------*)


PluckerRelations::argx = "Two arguments required.";
PluckerRelations::arg1 = "Argument is a positive integer two or more.";
PluckerRelations::arg2 = "Argument is a symbol.";
PluckerRelations[] := (Message[PluckerRelations::argx]; Abort[]);


PluckerRelations[argseq__] := Module[{args = {argseq}, n, v},
    If[Not[Length[args] == 2], (Message[PluckerRelations::argx]; Abort[])];
    n = args[[1]]; v = args[[2]];
    If[Not[And[IntegerQ[n], n > 1]], (Message[PluckerRelations::arg1]; 
        Abort[])]; 
    If[Not[SymbolQ[v]], (Message[PluckerRelations::arg2]; Abort[])];
    Return[Sort[PluckerRelationsInner[v, Table[i, {i, 1, n}]]]]]

PluckerRelationsInner[v_, M_] := 
  Module[{relations = {}, m = Length[M], subsetrelations, subsets, 
      wedgerelations = {}, top,verytop}, 
If[m < 2, Abort[]]; If[m == 2, Return[{}]];
    If[m == 3, 
      Return[{v[{M\[LeftDoubleBracket]1\[RightDoubleBracket]}] v[{M\
\[LeftDoubleBracket]2\[RightDoubleBracket], 
                  M\[LeftDoubleBracket]3\[RightDoubleBracket]}] - 
            v[{M\[LeftDoubleBracket]2\[RightDoubleBracket]}] v[{M\
\[LeftDoubleBracket]1\[RightDoubleBracket], 
                  M\[LeftDoubleBracket]3\[RightDoubleBracket]}] + 
            v[{M\[LeftDoubleBracket]3\[RightDoubleBracket]}] v[{M\
\[LeftDoubleBracket]1\[RightDoubleBracket], 
                  M\[LeftDoubleBracket]2\[RightDoubleBracket]}]}]];
    subsets = KKSubsets[M, m - 1];
    subsetrelations = 
      Join @@ Function[s, PluckerRelationsInner[v, s]] /@ subsets;
    Do[wedgerelations = 
        Join[WedgeRelation[v, M, i, m - i], wedgerelations], {i, 1, m - 1}]; 
    top = PluckerRelationsTop[v, M];
    verytop = prtts[v,M];
    relations = 
      Join[verytop, top, 
           PluckerRelationsTopTwo[v, M],
           {CramersRelation[v, M]}, 
           subsetrelations, 
           wedgerelations];
     Return[Select[Union[relations],SelectPositiveRelation]]]

SelectPositiveRelation[r_] := If[SameQ[r[[1]][[1]], -1], False, True]

PluckerRelationsTop[v_, M_] := 
  Module[{n = Length[M], ans = {}, inds = {}, nlist}, nlist = M; 
    If[n < 4, Return[{}]];
    Do[ans = 
        Join[{PluckerRelationTop[v, 
              M, {M\[LeftDoubleBracket]j\[RightDoubleBracket]}]}, ans], {j, 1,
         n}];
    Return[Reverse[ans]]]

PluckerRelationsTopTwo[v_, M_] := 
  Module[{n = Length[M], ans = {}, inds = {}, nlist}, 
    nlist = M; 
    inds = KKSubsets[M,2];
    If[n < 5, Return[{}]];
    Do[ans = 
        Join[{PluckerRelationTopTwo[v, 
              M, inds[[j]]]}, ans], {j, 1,
         Length[inds]}];
    Return[Reverse[ans]]]

PluckerRelationTop[v_, M_, index_] := 
  Module[{n = Length[M], nlist, i = Length[index], vars, s = 0, k, allvars,L},
    nlist = M;
    vars = Complement[nlist, index];
    allvars = KKSubsets[vars, n - 2];
    Do[(ind = Sort[Join[allvars[[j]], index]]; 
        L= Complement[nlist,ind];
        s = s + SignatureSubset[Join[index, L]]*
                SignatureSubset[Join[ind, L]]* v[Sort[Join[index, L]]]*v[ind]
                ), {j, 1, Length[allvars]}];
    Return[s]]


PluckerRelationTopTwo[v_, M_, index_] := 
  Module[{n = Length[M], nlist, i = Length[index], vars, s = 0, k, allvars,L},
    nlist = M;
    vars = Complement[nlist, index];
    allvars = KKSubsets[vars, n - 3];
    Do[(ind = Sort[Join[allvars[[j]], index]]; 
        L= Complement[nlist,ind];
        s = s + SignatureSubset[Join[index, L]]*
                SignatureSubset[Join[ind, L]]* v[Sort[Join[index, L]]]*v[ind]
                ), {j, 1, Length[allvars]}];
    Return[s]]

prtt[v_,M_,index_]:=Module[{n=Length[M],
  nlist,i=Length[index],vars,s=0,k,allvars,L},nlist=M;
    vars=Complement[nlist,index];
    allvars=KKSubsets[vars,2];
    Do[(ind=Sort[Join[allvars[[j]],index]];
        L=Complement[nlist,ind];
        s=s+SignatureSubset[Join[index,L]]*
      SignatureSubset[Join[ind,L]]*v[Sort[Join[
    index,L]]]*v[ind]),{j,1,Length[allvars]}];
    Return[s]]


prtts[v_,M_]:=Module[{n=Length[M],ans={},inds={},nlist},nlist=M;
    inds=KKSubsets[M,n-3];
    If[n<6,Return[{}]];
    Do[ans=Join[{prtt[v,M,inds[[j]]]},ans],{j,1,Length[inds]}];
    Return[Reverse[ans]]]

WedgeRelation[v_, M_, i_, j_] := 
  Module[{sum = 0, n = i + j, nlist, subi, subsi, subj, subsj}, nlist = M;
    If[! i + j == Length[M], Return[{}]]; subsi = KKSubsets[nlist, i];
    subsj = KKSubsets[nlist, j];
    If[! i == j, Do[subi = subsi\[LeftDoubleBracket]ii\[RightDoubleBracket];
        subj = subsj\[LeftDoubleBracket]jj\[RightDoubleBracket];
        If[subi \[Intersection] subj == {}, 
          sum = sum + SignatureSubset[Join[subi, subj]] v[subi] v[subj]], {ii,
           1, Length[subsi]}, {jj, 1, Length[subsj]}]];
    If[i == j, Return[WedgeRelation[v, M, i]]]; Return[{sum}]]

WedgeRelation[v_, M_, i_] := 
  Module[{sum = 0, n = 2 i, nlist, subi, subsi, subj}, nlist = M;
    If[Not[n == Length[M]], Return[{}]];
    subsi = KKSubsets[nlist, i];
    Do[subi = subsi\[LeftDoubleBracket]ii\[RightDoubleBracket];
      subj = subsi\[LeftDoubleBracket]jj\[RightDoubleBracket];
      If[Intersection[subi, subj] == {}, 
        sum = sum + SignatureSubset[Join[subi, subj]]*v[subi]*v[subj]], {ii, 
        1, Length[subsi]}, {jj, 1, Length[subsi]}];
    If[sum == 0, Return[{}]];
    Return[{Simplify[sum/2]}]]



CramersRelation[v_, M_] := 
  Module[{n = Length[M], s = 0, c}, Do[c = If[EvenQ[j + n], 1, -1];
      s = 
        c v[{M\[LeftDoubleBracket]j\[RightDoubleBracket]}] v[
              Complement[M, {M\[LeftDoubleBracket]j\[RightDoubleBracket]}]] + 
          s, {j, 1, n}]; Return[s]]


SignatureSubset[x_] := Module[{y, i, j, sorted, posns}, sorted = Sort[x];
    posns = 
      Function[j, First[First[j]]] /@ Function[i, Position[sorted, i]] /@ x;
    Return[Signature[posns]]]
(* -----> SignaturePermutation Combinatorica *)

PluckerClassic[v_, m_, n_, as_, bs_] := 
  Module[{t = m - Length[as], q = m - Length[bs], s, ans = 0, nlist, mlist, 
      clist, slist}, s = t + q; nlist = Table[i, {i, 1, n}]; 
    mlist = Table[i, {i, 1, m}];
    slist = Table[i, {i, 1, 2 m - Length[as] - Length[bs]}]; 
    subs = KKSubsets[slist, t];
    Do[perm = 
        Join[Sort[subs\[LeftDoubleBracket]j\[RightDoubleBracket]], 
          clist = Sort[
              Complement[slist, 
                subs\[LeftDoubleBracket]j\[RightDoubleBracket]]]];
      ans = 
        ans + Signature[perm] v[
              Join[as, subs\[LeftDoubleBracket]j\[RightDoubleBracket]]] v[
              Join[clist, bs]], {j, 1, Length[subs]}]; Return[ans]]
(* Combinatorica *)

TestPlucker4[relations_, v_] := 
  Module[{g, rules}, g = MakeMatrix["aa04", 4, 4];
    v[{n1_}] := MatrixMinor[g, {4}, {n1}];
    v[{n1_, n2_}] := MatrixMinor[g, {3, 4}, {n1, n2}];
    v[{n1_, n2_, n3_}] := MatrixMinor[g, {2, 3, 4}, {n1, n2, n3}];
    Return[Expand[relations]]]


TestPlucker5[relations_, v_] := Module[{g}, g = MakeMatrix["aa0", 5, 5];
    v[{n1_}] := MatrixMinor[g, {5}, {n1}];
    v[{n1_, n2_}] := MatrixMinor[g, {4, 5}, {n1, n2}];
    v[{n1_, n2_, n3_}] := MatrixMinor[g, {3, 4, 5}, {n1, n2, n3}];
    v[{n1_, n2_, n3_, n4_}] := 
      MatrixMinor[g, {2, 3, 4, 5}, {n1, n2, n3, n4}];
    Return[Expand[relations]]]

getindices[k_, n_] := Module[{lamk, vec = {}},
    lamk = KKSubsets[Table[i, {i, 1, n}], k];
    Do[vec = Join[vec, {lamk[[i]]}], {i, 1, Binomial[n, k]}];
    Return[vec]]

(*------------------------------------------------------*)

PluckerInverse::argx = "One argument required.";
PluckerInverse::arg1 = 
    "Argument is a list of lists of integers of length n chose i, 1<=i<=n-1";
PluckerInverse[] := (Message[PluckerInverse::argx]; Abort[]);


PluckerInverse[argseq__] := Module[{args = {argseq}, Ms},
    If[Not[Length[args] == 1], (Message[PluckerInverse::argx]; Abort[])];
    Ms = args[[1]];
    If[Not[And[ListQ[Ms], MinorVectorQ[Ms]]], (Message[PluckerInverse::arg1]; 
        Abort[])];  Return[InvertPlucker2[Ms]]]

GeneralIntegerListQ[x_] := And[ListQ[x], {True} == Union[Map[IntegerQ, x]]]


MinorVectorQ[Ms_] := Module[{n = Length[Ms] + 1, ans = True},
    Do[If[
        Not[And[ListQ[Ms[[i]]], GeneralIntegerListQ[Ms[[i]]], 
            Length[Ms[[i]]] == Binomial[n, i]]], (ans = False; Break[])], {i, 
        1, n - 1}];
    Return[ans]]




Ranking[x_, y_] := First[First[Position[y, x]]]

ind2mi[Ms_, ind_, n_] := 
  Ms[[Length[ind]]][[Ranking[ind, getindices[Length[ind], n]]]]

ZeroFreeVbles[sol_] := Module[{ans = {}},
    Do[(If[sol[[j]][[0]] == Element, ans = Join[{sol[[j]][[1]] == 0}, ans]]; 
        If[Not[SameQ[sol[[j]][[0]], Element]], 
          ans = Join[{sol[[j]]}, ans]]), {j, 1, Length[sol]}];
    Return[Apply[And, Reverse[ans]]]]

InvertPlucker2[Ms_] := 
 Module[{a, nlist, lj, u, U, v, V, subsj, coej, sub, rowi, mis, subc, 
   n, C, paramrules, tag = True, coes, x, eqi, solni, solnii, eqnsi, 
   row1, indi},
  n = Length[First[Ms]];
  a = Table[0, {i, 1, n}, {j, 1, n}];
  nlist = Table[i, {i, 1, n}];
  (*bottom row*)
  Do[a[[n, i]] = Ms[[1, i]], {i, 1, n}];
  (*top row*)v = n - 1;
  coes = {};
  Do[(mis = Ms[[v]]; indi = getindices[v, n]; coej = 0;
    lj = Complement[nlist, {j}];
    (V = lj;
     coej = 
      Signature[Join[{j}, V]]*
       mis[[First[First[Position[indi, V]]]]]);
(* combinatorica *)
    coes = Join[coes, {coej}];), {j, 1, n}]; 
  row1 = Part[Apply[ExtendedGCD, coes], 2];
  Do[a[[1, j]] = row1[[j]], {j, 1, n}];
  (*middle rows*)
  paramrules = Table[Rule[C[i], 0], {i, 1, n}];
  Do[(eqnsi = {}; v = n - i + 1;
    indi = getindices[v, n];
    Do[(V = indi[[j]];
      eqi = {ind2mi[Ms, V, n] == 
         Sum[(-1)^(k + 1)*x[V[[k]]]*
           ind2mi[Ms, Complement[V, {V[[k]]}], n], {k, 1, v}]};
      eqnsi = Join[eqi, eqnsi]), {j, 1, Binomial[n, v]}];
    solni = 
     Reduce[eqnsi, Table[x[j], {j, 1, n}], Integers, 
       Backsubstitution -> True, GeneratedParameters -> C] /. 
      paramrules;
    If[solni == False, (tag = False; Break[])];
    solnii = ZeroFreeVbles[solni];
    Do[a[[i, j]] = x[j] /. ToRules[solnii], {j, 1, n}]), {i, 2, 
    n - 1}];
  If[Or[tag == False, a == {}], Return[False]];
  Return[a]]

(*-----------------------------------------------*)
PluckerCoordinates::argx = "One argument required.";
PluckerCoordinates::arg1 = "Argument is a square CRE matrix.";
PluckerCoordinates[] := (Message[PluckerCoordinates::argx]; Abort[])

PluckerCoordinates[argseq__] := 
  Module[{args = {argseq}, a, ans = {}, ansv, inds, ind},
    If[Not[Length[args] == 1], (Message[PluckerCoordinates::argx]; Abort[])];
    a = args[[1]];
    If[Not[
        And[SSquareMatrixQ[a], 
          ExprMatrixQ[a]]], (Message[PluckerCoordinates::arg1]; 
        Abort[])];  
    n = Length[a];
    
    Do[(inds = getindices[v, n]; 
        ansv = Map[
            Function[ind, MatrixMinor[a, Table[i, {i, n - v + 1, n}], ind]], 
            inds]; ans = Join[ans, {ansv}];), {v, 1, n - 1}];
    Return[ans]]


(*=== === === === === === === === === === === === === ==*)

expandvec[s_, ci_] := Module[{ans = {}},
  
    Do[(If[NumberQ[s[[j]]], ans = Join[{{s[[j]]}},
     ans], ans = Join[{Table[i, {i, 0, ci - 1}]}, 
                ans]];;;), {j, 1, Length[s]}];
    Return[Reverse[ans]]]

getvector[k_, n_, v_] := Module[{lamk, vec = {}},
    lamk = KKSubsets[Table[i, {i, 1, n}], k + 1];
    Do[vec = Join[vec, {v[lamk[[i]]]}], {i, 1, Binomial[n, k + 1]}];
    Return[vec]]


getmatrix[k_, n_, v_] := Module[{mat, m = Binomial[n, k + 1], lamk, r},
    lamk = KKSubsets[Table[i, {i, 1, n}], k + 1];
    mat = Table[0, {i, 1, m}, {j, 1, n}];
    Do[  Do[
        mat[[i, j]] = 
          If[MemberQ[lamk[[i]], 
              j], (r = 
                First[First[Position[lamk[[i]], j]]]; (-1)^(r - 1)v[
                  Complement[lamk[[i]], {j}]]), 0], {j, 1, n}], {i, 1, m}];
    Return[mat]]

gcdconst[k_, n_, v_, minors_] := Module[{lamk, vec = {}, redvec},
    lamk = KKSubsets[Table[i, {i, 1, n}], k + 1];
    Do[vec = Join[vec, {v[lamk[[i]]]}], {i, 1, Binomial[n, k + 1]}];
    redvec = DeleteCases[vec /. minors, 0];
    Return[redvec]]

getminors[w_, v_, c_] := Module[{n = Length[w], tw, sigma, ans = {}},
tw = Transpose[w];
sigma = Table[0, {i, 1, n}];
Do[{i = 1, While[tw[[i, j]] == 0, i = i + 1], sigma[[j]] = i}, {j, 1, n}];
Do[{ind = Sort[Table[sigma[[j]], {j, n - i + 1, n}]], 
        ans = Join[ans, {Rule[v[ind], fixsignature[w, i]*c[[i]]]}, 
            Map[Function[x, Rule[x, 0]], Map[v, lexibded[n, ind]]]]}, {i, 1, 
        n - 1}];
Return[ans]]

getciminors[w_, v_, c_] := Module[{n = Length[w], tw, sigma, ans = {}},
tw = Transpose[w];
sigma = Table[0, {i, 1, n}];
Do[{i = 1, While[tw[[i, j]] == 0, i = i + 1], sigma[[j]] = i}, {j, 1, n}];
Do[(ind = Sort[Table[sigma[[j]], {j, n - i + 1, n}]]; 
        ans = Join[ans, {Rule[v[ind], fixsignature[w, i]*c[[i]]]}]), {i, 1, 
        n - 1}];
Return[ans]]

fixsignature[w_, i_] := 
  Module[{n = Length[w], tw, sigma, a = {}, b = {}, ii},
    tw = Transpose[w];
sigma = Table[0, {ii, 1, n}];
Do[{ii = 1, While[tw[[ii, j]] == 0, ii = ii + 1], sigma[[j]] = ii}, {j, 1, 
        n}];
    a = Reverse[Take[sigma, -i]];
    b = Join[Table[j, {j, 1, n - i}], Reverse[Table[j, {j, n - i + 1, n}]]];
    Return[SignatureSubset[a]*SignatureSubset[b]]]
(* combinatorica *)

lexibded[n_, ms_] := 
  Module[{l = Length[ms], m1, nlist, ans = {}, recsets, joinm1, joinj},
    If[l == 0, Return[{}]];
    m1 = First[ms];
    If[l == 1, Return[Map[List, Table[j, {j, 1, m1 - 1}]]]];
    nlist = Table[j, {j, 1, n}];
    Do[(joinj[x_] := Join[{j}, x]; 
        ans = Join[ans, 
            Map[Sort, 
              Map[joinj, KKSubsets[Complement[nlist, {j}], l - 1]]]]), {j, 1, 
        m1 - 1}];
    joinm1[x_] := Join[{m1}, x];
    recsets = Map[Sort, Map[joinm1, lexibded[n, Rest[ms]]]];
    ans = removeinvalid[Join[ans, recsets]];
    Return[Union[ans]]]

getleximinors[w_, v_, c_] := Module[{n = Length[w], tw, sigma, ans = {}},
tw = Transpose[w];
sigma = Table[0, {i, 1, n}];
Do[{i = 1, While[tw[[i, j]] == 0, i = i + 1], sigma[[j]] = i}, {j, 1, n}];
Do[{ind = Sort[Table[sigma[[j]], {j, n - i + 1, n}]], 
        ans = Join[ans, 
            Map[Function[x, Rule[x, 0]], Map[v, lexibded[n, ind]]]]}, {i, 1, 
        n - 1}];
Return[ans]]

removeinvalid[x_] := Module[{ans = {}},
    Do[(ind = x[[j]]; 
        If[Length[ind] == Length[Union[ind]], ans = Join[{ind}, ans]]), {j, 1,
         Length[x]}];
    Return[Reverse[ans]]]

imagemembership[vect_, mat_] := Module[{m, n, x, eqns = {}, vars},
    m = Length[vect];
    n = Length[First[mat]];
    vars = Table[x[i], {j, 1, n}];
    Do[eqns = 
        Join[eqns, {0 == vect[[i]] - Sum[mat[[i, j]]*x[j], {j, 1, n}]}], {i, 
        1, m}];
    Return[Eliminate[eqns, vars]]]

getoneminor[g_, ind_] := Module[{n, i, z},
    n = Length[g]; i = Length[ind];
    z = Take[g, -i]; m = Part[Transpose[z], ind];
    Return[Det[m]]]


getzeroelements[a_, lami_] := 
  Module[{inds, ind, n = Length[a], i = Length[lami], z, ans = {}, m},
    inds = lexibded[n, lami]; Print[inds];
    z = Take[a, -i];
    Do[(ind = inds[[j]]; m = Part[Transpose[z], ind]; 
        ans = Join[{Det[m]}, ans]), {j, 1, Length[inds]}];
    Return[Reverse[ans]]]


getlamis[w_] := Module[{n = Length[w], tw, sigma, ans = {}},
tw = Transpose[w];
sigma = Table[0, {i, 1, n}];
Do[{i = 1, While[tw[[i, j]] == 0, i = i + 1], sigma[[j]] = i}, {j, 1, n}];
Do[{ind = Sort[Table[sigma[[j]], {j, n - i + 1, n}]], 
        ans = Join[ans, {ind}]}, {i, 1, n - 1}];
Return[ans]]

(*=== === === === === === == === === === === === === ==*)

(*=======================================================================*)
MPEisensteinSeriesInternal[s_, v_] := 
  Module[{lams, n = Length[v] + 1},
    lams = MPEisensteinLambdas[v]; 
    Return[Product[Zeta[s - lams[[i]]], {i, 1, n}]]]

MPEisensteinSeries::argx = "Two arguments required.";
MPEisensteinSeries::arg1 = "The first argument is a CRE.";
MPEisensteinSeries::arg2 = "The second argument is a list of n-1 CREs.";
MPEisensteinSeries[] := (Message[MPEisensteinSeries::argx]; Abort[])
MPEisensteinSeries[argseq__] := 
  Module[{args = List[argseq], s, v}, 
    If[Length[args] != 2, (Message[MPEisensteinSeries::argx]; Abort[])];
    s = First[args]; v = Part[args, 2];
    If[Not[creq[s]], (Message[MPEisensteinSeries::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPEisensteinSeries::arg2];
        Abort[])];
     Return[ MPEisensteinSeriesInternal[s, v]]]

MPEisensteinLambdasInternal[v_] := 
  Module[{n = Length[v] + 1, c, m, SLog, cvec, pwers, pwers2, logpwers2, 
      ans},
    SLog[x_^y_] := y SLog[x];
     SLog[x_^y_ * z_] := SLog[x^y] + SLog[z];
    cvec = Join[Table[c[i]/m, {i, 1, n - 1}], {1}];
    pwers = IFun[v - 1/n, DiagonalToMatrix[cvec]];
    pwers2 = Simplify[PowerExpand[pwers]];
    logpwers2 = SLog[pwers2];
    Do[ans[i] = Coefficient[logpwers2, SLog[c[i]]], {i, 1, n - 1}];
    ans[n] = Coefficient[logpwers2, SLog[m]];
    Return[Table[ans[i], {i, 1, n}]]]

MPEisensteinLambdas::argx = "One argument required.";
MPEisensteinLambdas::arg1 = "The argument is a list of n-1 CREs.";
MPEisensteinLambdas[] := (Message[MPEisensteinLambdas::argx]; Abort[])
MPEisensteinLambdas[argseq__] := 
  Module[{args = List[argseq], v}, 
    If[Length[args] != 1, (Message[MPEisensteinLambdas::argx]; Abort[])];
    v = Part[args, 1];
    If[Not[ExprListQ[v]], (Message[MPEisensteinLambdas::arg1];
        Abort[])];
     Return[ MPEisensteinLambdasInternal[v]]]

MPEisensteinGammaInternal[s_, v_] := Module[{n = Length[v] + 1},
    lams = MPEisensteinLambdas[v];
    Return[(Pi^(- n s/2))* Product[Gamma[(s - lams[[i]])/2], {i, 1, n}]]]

MPEisensteinGamma::argx = "Two arguments required.";
MPEisensteinGamma::arg1 = "The first argument is a CRE.";
MPEisensteinGamma::arg2 = "The second argument is a list of n-1 CREs.";
MPEisensteinGamma[] := (Message[MPEisensteinGamma::argx]; Abort[])
MPEisensteinGamma[argseq__] := 
  Module[{args = List[argseq], s, v}, 
    If[Length[args] != 2, (Message[MPEisensteinGamma::argx]; Abort[])];
    s = First[args]; v = Part[args, 2];
    If[Not[creq[s]], (Message[MPEisensteinGamma::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPEisensteinGamma::arg2];
        Abort[])];
     Return[ MPEisensteinGammaInternal[s, v]]]

symmetricandexterior[s_, v_, k_, flag_] := 
  Module[{lmas, n, indices, m, b, ans = 1},
    n = Length[v] + 1;
    lams = MPEisensteinLambdas[v];
    indices = If[flag, getsymmetricindices[n, k], getexteriorindices[n, k]];
    m = Length[indices];
    Do[(indj = indices[[j]]; b[j] = Sum[lams[[indj[[i]]]], {i, 1, k}]; 
        ans = ans*Zeta[s - b[j]]), {j, 1, m}];
    Return[ans]]

MPSymmetricPowerLFun::argx = "Three arguments required.";
MPSymmetricPowerLFun::arg1 = "The first argument is a CRE.";
MPSymmetricPowerLFun::arg2 = "The second argument is a list of n-1 CREs.";
MPSymmetricPowerLFun::arg3 = "The third argument is a positive integer."; 
MPSymmetricPowerLFun[] := (Message[MPSymmetricPowerLFun::argx]; Abort[])
MPSymmetricPowerLFun[argseq__] := 
  Module[{args = List[argseq], s, v, k, n}, 
    If[Length[args] != 3, (Message[MPSymmetricPowerLFun::argx]; Abort[])];
    s = First[args]; v = Part[args, 2]; k = Part[args, 3];
    If[Not[creq[s]], (Message[MPSymmetricPowerLFun::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPSymmetricPowerLFun::arg2];
        Abort[])];
    n = Length[v] + 1;
    If[Not[And[IntegerQ[k], k > 0]], (Message[MPSymmetricPowerLFun::arg3]; 
        Abort[])];
     Return[ symmetricandexterior[s, v, k, True]]]

MPExteriorPowerLFun::argx = "Three arguments required.";
MPExteriorPowerLFun::arg1 = "The first argument is a CRE.";
MPExteriorPowerLFun::arg2 = "The second argument is a list of n-1 CREs.";
MPExteriorPowerLFun::arg3 = "The third argument is a positive integer <= n."; 
MPExteriorPowerLFun[] := (Message[MPExteriorPowerLFun::argx]; Abort[])
MPExteriorPowerLFun[argseq__] := 
  Module[{args = List[argseq], s, v, k, n}, 
    If[Length[args] != 3, (Message[MPExteriorPowerLFun::argx]; Abort[])];
    s = First[args]; v = Part[args, 2]; k = Part[args, 3];
    If[Not[creq[s]], (Message[MPExteriorPowerLFun::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPExteriorPowerLFun::arg2];
        Abort[])];
    n = Length[v] + 1;
    If[Not[
        And[IntegerQ[k], k > 0, k < n + 1]], (Message[
          MPExteriorPowerLFun::arg3]; Abort[])];
     Return[ symmetricandexterior[s, v, k, False]]]

getsymmetricindices[n_, k_] :=
  If[k == 1, Map[List, Table[i, {i, 1, n}]], 
    Union[Map[Sort, 
        Apply[ProductSet, 
          Table[Table[i, {i, 1, n}], {j, 1, k}]]]]]

getexteriorindices[n_, k_] := Module[{ans = {}, sym, symj},
    If[k == 1, Return[getsymmetricindices[n, 1]]];
    sym = getsymmetricindices[n, k];
    Do[(symj = sym[[j]]; 
        If[Length[symj] == Length[Union[symj]], ans = Join[{symj}, ans]]), {j,
         1, Length[sym]}];
    Return[Reverse[ans]]]


MPPowerGamma[s_, v_, k_, flag_] := 
  Module[{lams, lamws, n, indices, w, indj, shifts, eqns, ws, solns}, 
    n = Length[v] + 1;
    lams = MPEisensteinLambdas[v];
    indices = If[flag, getsymmetricindices[n, k], getexteriorindices[n, k]];
    m = Length[indices];
    If[m==1, Return[Pi^(-s/2)Gamma[s/2]]];
    Do[(indj = indices[[j]]; b[j] = Sum[lams[[indj[[i]]]], {i, 1, k}]), {j, 1,
         m - 1}]; 
    ws = Table[w[i], {i, 1, m - 1}];
    shifts = {Take[MPEisensteinLambdasInternal[ws], m - 1], 
        Table[b[j], {j, 1, m - 1}]};
    eqns = MapThread[Equal, shifts];
    solns = First[Solve[eqns, ws]]; 
    lamws = MPEisensteinLambdas[ws];
    Return[
      Simplify[Pi^(-m s/2)Product[Gamma[(s - lamws[[i]])/2], {i, 1, m}] /. 
          solns]]]

MPSymmetricPowerGamma::argx = "Three arguments required.";
MPSymmetricPowerGamma::arg1 = "The first argument is a CRE.";
MPSymmetricPowerGamma::arg2 = "The second argument is a list of n-1 CREs.";
MPSymmetricPowerGamma::arg3 = "The third argument is a positive integer."; 
MPSymmetricPowerGamma[] := (Message[MPSymmetricPowerGamma::argx]; Abort[])
MPSymmetricPowerGamma[argseq__] := 
  Module[{args = List[argseq], s, v, k, n}, 
    If[Length[args] != 3, (Message[MPSymmetricPowerGamma::argx]; Abort[])];
    s = First[args]; v = Part[args, 2]; k = Part[args, 3];
    If[Not[creq[s]], (Message[MPSymmetricPowerGamma::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPSymmetricPowerGamma::arg2];
        Abort[])];
    n = Length[v] + 1;
    If[Not[And[IntegerQ[k], k > 0]], (Message[MPSymmetricPowerGamma::arg3]; 
        Abort[])];
     Return[ MPPowerGamma[s, v, k, True]]]

MPExteriorPowerGamma::argx = "Three arguments required.";
MPExteriorPowerGamma::arg1 = "The first argument is a CRE.";
MPExteriorPowerGamma::arg2 = "The second argument is a list of n-1 CREs.";
MPExteriorPowerGamma::arg3 = "The third argument is a positive integer <= \
n."; 
MPExteriorPowerGamma[] := (Message[MPExteriorPowerGamma::argx]; 
    Abort[])
MPExteriorPowerGamma[argseq__] := 
  Module[{args = List[argseq], s, v, k, n}, 
    If[Length[args] != 3, (Message[MPExteriorPowerGamma::argx]; Abort[])];
    s = First[args]; v = Part[args, 2]; k = Part[args, 3];
    If[Not[creq[s]], (Message[MPExteriorPowerGamma::arg1];
        Abort[])];
    If[Not[ExprListQ[v]], (Message[MPExteriorPowerGamma::arg2];
        Abort[])];
    n = Length[v] + 1;
    If[Not[
        And[IntegerQ[k], k > 0, k < n + 1]], (Message[
          MPExteriorPowerGamma::arg3]; Abort[])];
     Return[ MPPowerGamma[s, v, k, False]]]
(*=========================================================================*)


$RecursionLimit = 256;
$IterationLimit = 4096;

Forms[i_] := {}; AllScalForms = {}; AllDifForms = {}; AllSymbols = \
{}; FormVars = {_Wedge, _d}; ZeroVars = {_Wedge}; HeadList = {}; \
ScalHeadList = {}; ZeroHeads = {}; nodHeads = {Pattern, Blank, 
  Condition, RuleDelayed};


zeroQ[0] = True;
zeroQ[x_List] := And @@ (zeroQ /@ Union[Flatten[x]]); 
zeroQ[x_] := False;
nonZeroL[x_] := Select[Union[Flatten[{x}]], ! zeroQ[#] &]


indepTerms[x_, opt_: {}] := 
 If[NumericQ[x /. _Wedge -> 1], reWrite[indep$cTerms[Factor[x], opt]],
   indep$cTerms[x, opt]]
indep$cTerms[x_, opt_: {}] := redPlus[redTimes[nonZeroL[x], opt]];


redTimes[x_, opt_: {}] := 
 Block[{tmp$ = Select[x, Head[#1] === Times &], rst$}, 
  rst$ = Complement[x, tmp$]; 
  Union[rst$, Map[numb2[#, opt] &, tmp$]]]


numb1[x_, opt_: {}] := 1 /; NumericQ[x]; 
numb1[x_^y_., opt_: {}] := 1 /; MemberQ[opt, x];
numb1[x_, opt_: {}] := x; 
numb2[x_, opt_: {}] := Times @@ Map[numb1[#, opt] &, Level[x, 1]]


redPlus[x_] := 
 Block[{tmp$ = Select[x, Head[#1] === Plus &], cmn$, rst$, i$, negi$},
   rst$ = Complement[x, tmp$]; cmn$ = Intersection[tmp$, -tmp$]; 
  If[cmn$ === {}, Return[x]]; i$ = 1; 
  While[i$ < Length[cmn$], negi$ = Select[cmn$, #1 === -cmn$[[i$]] &];
    If[negi$ =!= {}, cmn$ = Complement[cmn$, negi$]]; i$ = i$ + 1]; 
  Union[rst$, Complement[tmp$, -cmn$]]]


indepCoefs[y_] := 
  indepTerms[
   If[FormDegree[y] > 0, 
    Flatten[(#1\[Wedge]$\[Wedge]$ &) /@ indepTerms[y] //. 
       a_Wedge x_ + (z_: Null) :> {x, z}] /. Null -> Sequence[], y]];


Qvalue[x_List] := 
 Length[Select[
   OrderedQ /@ Transpose[{Drop[x, -1], Rest[x]}], # === False &]]


AllCaseTable[x_List] := 
 With[{rvstr = Reverse[Map[Transpose, x] /. -z_ :> z]}, 
  If[traceOption === 2, 
   Join[Table[{Qvalue[RotateLeft[x, i]], RotateLeft[x, i], i}, {i, 0, 
      Length[x] - 1}], 
    Table[{Qvalue[RotateRight[rvstr, i]], 
      RotateRight[rvstr, i], -i - 1}, {i, 0, Length[x] - 1}]], 
   Table[{Qvalue[RotateLeft[x, i]], RotateLeft[x, i], i}, {i, 0, 
     Length[x] - 1}]]]


bestOrd[x_List] := 
  Block[{tmp$ = Sort[AllCaseTable[x]]}, 
   tmp$ = Select[tmp$, First[#1] === tmp$[[1, 1]]
      &]; 
   If[Length[tmp$] > 1, 
    First[Sort[Select[tmp$, #1[[2]] === tmp$[[1, 2]] &], 
      Abs[Last[#1]] < Abs[Last
          [#2]] &]], First[tmp$]]];


splitFDprod[n_, 
  FDx_List] := (Plus @@ Take[FDx, n]) (Plus @@ Drop[FDx, n]); 
allFDprods[FDx_List] := 
 Sum[FDx[[i]]*(Plus @@ Drop[FDx, i]), {i, Length[FDx] - 1}];


FormDegree[x_Plus] := FormDegree[First[x]]; 
FormDegree[x_Times] := Plus @@ FormDegree /@ List @@ x; 
FormDegree[x_Wedge] := Plus @@ FormDegree /@ List @@ x; 
FormDegree[d[x_]] := FormDegree[d[x]] = 1 + FormDegree[x]; 
FormDegree[x_List] := FormDegree[Last[Union[Flatten[x]]]];
FormDegree[x_] := 0;


FormCoef[x_List, y_] := (FormCoef[#, y] &) /@ x; 
FormCoef[x_Plus, y_] := (FormCoef[#1, y] &) /@ x;
FormCoef[x_.*y_, y_] := x; 
FormCoef[x_, y_Wedge] := Fold[FormCoef, x, Reverse[List @@ y]];
FormCoef[x_, y_] := 0 /; FreeQ[x /. d[y] -> d[$J$], y]; 
FormCoef[x_, 
  y_] := (((x /. d[y] -> d[$J$]) - (x /. d[y] -> d[$J$] /. 
         y -> 0)) /. y -> 1 /. d[$J$] -> d[y]) /; 
  EvenQ[FormDegree[y]]; 
FormCoef[x_, 
  y_] := (((x /. d[y] -> d[$J$]) - (x /. d[y] -> d[$J$] /. 
         y -> 0)) /. 
     d[$J$] -> 
      d[y] /. {HoldPattern[
       Wedge[p___, y, q___]] :> (-1)^Length[Select[{q},
          OddQ[FormDegree[#]] &]] Wedge[p, q]})/OddQ[FormDegree[y]];



DeclareForms[z__] := 
  Module[{h, x = {{z}}, xi, rxi, k, tmp = AllScalForms}, 
   While[Head[x[[1, 1]]] === List, x = First[x]]; 
   Do[xi = x[[i]]; rxi = Rest[xi]; 
    Forms[First[xi]] = Union[Forms[First[xi]], rxi];
    AllScalForms = Union[AllScalForms, rxi];
    Do[h = Head[rxi[[j]]]; 
     If[h === Symbol, FormDegree[rxi[[j]]] = First[xi]; 
      BasicScalFormQ[rxi[[j]]] = True, FormDegree[_h] = First[xi]; 
      BasicScalFormQ[_h] = True], {j, Length[rxi]}], {i, Length[x]}];
   AllDifForms = Union[AllScalForms]; 
   ScalHeadList = 
    Complement[Union[Head /@ AllScalForms, ScalHeadList], {Symbol}]; 
   AllSymbols = Union[AllDifForms];
   HeadList = Union[Head /@ AllSymbols]; 
   DifFormSymbols = Drop[AllDifForms, -Length[Union[ScalHeadList]]]; 
   k = Thread[Blank[Union[ScalHeadList]]]; 
   FormVars = Flatten[{_Wedge, _d, Union[k], Union[DifFormSymbols]}]; 
   If[tmp =!= {}, Print["OK"]]];



If[$VersionNumber < 6., st$pr6$[x__] := StylePrint[x], 
  st$pr6$[x__] := Print[Style[x]]];


stPrint[x_, y_: Null] := 
  If[y === Null, 
   st$pr6$[x, FontSize -> 16, FontSlant -> "Italic", 
    FontColor -> RGBColor[1, 0, 1]], 
   st$pr6$[x, FontSize -> 20, FontWeight -> "Bold", 
    FontColor -> RGBColor[1, 0, 1]]];



NoDif[z__] := (nodHeads = Union[nodHeads, Flatten[{z}]])



BasicScalFormQ[d[x_]] := True;
BasicScalFormQ[x_] := False;


ScalFormQ[x_Times] := Or @@ ScalFormQ /@ List @@ x;
ScalFormQ[x_Wedge | x_Plus] := And @@ ScalFormQ /@ List @@ x; 
ScalFormQ[x_] := BasicScalFormQ[x];

SetAttributes[d, {Listable}]; 
d[x_Times | x_Wedge] := 
 d[First[x]]\[Wedge]Rest[x] + (-1)^FormDegree[First[x]]*
   First[x]\[Wedge]d[Rest[x]]; d[x_?NumericQ | x_d] = 0; 
d[Power[y_, n_]] := n y^(n - 1) d[y] + y^n Log[y] d[n]; 
d[x_Plus] := d /@ x; d[x_Rule | x_Equal] := reWrite[d /@ x]; 
d[h_[y__] /; FreeQ[nodHeads, h]] := 
 Sum[(Derivative[
        Sequence @@ 
         RotateRight[Append[Table[0, {Length[{y}] - 1}], 1], i]][h][
      y]) d[{y}[[i]]], {i, Length[{y}]}] /; 
  FormDegree[h[y]] === 0 && 
   FreeQ[{Integer, Blank, Pattern, Condition}, Head[First[{y}]]];


Wedge[x_] := x /; Length[{x}] < 2 && Head[{x}[[1]]] =!= Pattern
Default[Wedge] := 1; SetAttributes[Wedge, {Flat, OneIdentity}]; 
Wedge[0, y__] = 0; Wedge[x__, 0] = 0; 
Wedge[x__, y_Plus] := Plus @@ Map[Wedge[x, #] &, List @@ y]; 
Wedge[x_Plus, y__] := Plus @@ Map[Wedge[#, y] &, List @@ x]; 
Wedge[z__, Times[x_, y_]] := 
 Times[x, Wedge[z, y]] /; NumericQ[x] || (FormDegree[x] === 0); 
Wedge[Times[x_, y_], z__] := 
 Times[x, Wedge[y, z]] /; NumericQ[x] || (FormDegree[x] === 0); 
Wedge[x_^n_., y_] := x^n*y /; FormDegree[x] === 0; 
Wedge[y_, x_^n_.] := x^n*y /; FormDegree[x] === 0; 
Wedge[x_, y___, x_] := 0 /; OddQ[FormDegree[x]] && ScalFormQ[x]; 
Wedge[x__] := 
 Block[{doubL = Transpose[{FormDegree /@ {x}, {x}}]}, 
   Signature[Select[doubL, OddQ[First[#]] &]] Wedge @@ 
     Map[Last[#1] &, Sort[doubL]]] /; 
  Union[BasicScalFormQ /@ {x}] === {True} && 
   Map[Last[#1] &, Sort[Transpose[{FormDegree /@ {x}, {x}}]]] =!= {x};
 Wedge[x_List, y_List] := 
 If[(FormDegree[x] === 0) || (FormDegree[y] === 0), Dot[x, y], 
  Inner[Wedge, x, y]]; 
Wedge[x_, y_List] := 
 If[FormDegree[x] === 0 || (FormDegree[y] === 0), x *y, 
  Map[Wedge[x, #] &, y]];
Wedge[x_List, y_] := 
  If[FormDegree[y] === 0 || (FormDegree[x] === 0), x*y, 
   Map[Wedge[#, y] &, x]];


TrigRules = {x_.*Cos[h_]^2 + x_.*Sin[h_]^2 :> x, 
  x_.*Cot[h_]^2 + y_.*Csc[h_]^2 :> 
   y /; x + y === 
     0, (Cot[h_] - Csc[h_]) (Cot[h_] + 
      Csc[h_]) -> -1}; HypTrigRules = {x_.*Cosh[h_]^2 + 
    y_.*Sinh[h_]^2 :> x /; x + y === 0, 
  x_.*Coth[h_]^2 + y_.*Csch[h_]^2 :> 
   x /; x + y === 0, (Cosh[x_] - Sinh[x_]) (Cosh[x_] + Sinh[x_]) -> 
   1, (Coth[x_] - Csch[x_]) (Coth[x_] + Csch[x_]) -> 1};


simpRules = {x_.*Cos[h_]^2 + x_.*Sin[h_]^2 :> x}; varList = {};



SetAttributes[reWrite, {Listable}]; reWrite[0] = 0; 
reWrite[x_Equal] := Equal[reWrite[First[x] - Last[x]], 0]; 
reWrite[x_Rule] := Rule[First[x], reWrite[Last[x]]]; 
reWrite[x_RuleDelayed] := 
 Rule[First[x], reWrite[x[[2]]]] /. Rule -> RuleDelayed; 
reWrite[x_Condition] := x /. x[[1]] -> reWrite[x[[1]]]; 
reWrite[x_] := (Collect
     [x /. simpRules, FormVars, bestFacRul] /. simpRules) /; 
  FormDegree[x] > 0; 
reWrite[x_] := bestFacRul[x /. simpRules] /. simpRules;



FreeALLQ[x_, y_List] := And @@ (FreeQ[x, #] & /@ y);
bestFacRul[x_] := 
 If[varList === {} || FreeALLQ[x, varList], minFacRul[x], 
  varFacRul[x]]; 
varFacRul[x_] := 
 Collect[Expand[x] /. simpRules, varList, Factor] /. simpRules; 
minFacRul[x_] := Factor[Expand[x] /. simpRules] /. simpRules;




(*---------------------------------------------------------*)
Print[StyleForm["GL(n)PACK loaded.",FontWeight->"Bold"]];
End[]

On[General::"spell",General::"spell1",UpSet::"write"];

EndPackage[]


