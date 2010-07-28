(* ::Package:: *)

(* ::Title:: *)
(*MathMPS*)
(* -- Tensor Network algorithms in Mathematica*)


(* ::Section:: *)
(*Header*)


Print["MPS manipulation from Mathematica.
MPS Manipulation and creation:
	MPSProductState
	MPSExpandBond
	MPSNormalize
	MPSCanonize
Reading and saving to disk:
	MPSRead
	MPSSave
Expectation values:
	MPSExpectation
	MPSCorrelation
Algorithms:
	MPSApproximate
	MPSMinimizeEnergy
Hamiltonian Creation:
	MPSInitHamiltonian
	MPSHamiltonianAdd

Also available are TEBD algorithms:

Creation and disk interaction:
	TEBDProductState
	TEBDRead
	TEBDSave
Real and imaginary time evolution:
	TEBDEvolve
Expectation values:
	TEBDExpectation
	TEBDCorrelation
	TEBDentropy
	TEBDentropyList
Hamiltonians and evolution operators:
	TEBDInitHamiltonian
	TEBDHamiltonianAdd
	TEBDEvolutionFromHamiltonian
	TEBDInitEvolution
	TEBDEvolutionAddGate
"]


MPSProductState::usage="MPSProductState[numTensors] creates a Matrix Product State of length numTensors. 
Options:
   Spin->s (default = 2)
   Bond->\[Chi] (default = 10)
   Type-> \"Random\" (default), \"Identity\", \"Decaying\". 
As of now the Random option produces an warning but the state is correct." 


MPSExpandBond::usage="MPSExpandBond[mps,new\[Chi]] expands the bond dimension \[Chi] of mps to new\[Chi]. Returns the expanded mps as a new variable"


MPSNormalize::usage="MPSNormalize[mps] normalizes mps so that MPSOverlap[mps,mps]=1. It changes the state."


MPSCanonize::usage="MPSCanonize[mps] transforms mps into its canonical form. 
Options:
      Direction->\"Left\" (default) or \"Right\""


MPSSave::usage="MPSSave[MPS,filename] saves MPS into ASCII files for later retrieval. Right now it produces a lot of files, to be changed in the future. filename needs to be a string"



MPSRead::usage="MPSRead[filename] reads the files produced by MPSSave and returns a new MPS object"


MPSOverlap::usage="MPSOverlap[mps1,mps2] returns the overlap <psi_1|psi_2>, where psi_i is the state represented by mpsi"


MPSExpectation::usage="
MPSExpectation[mps,operator,site]
MPSExpectation[mps,operator,site1,site2]
MPSExpectation[mps,operator1,site1,operator2,site2]

Computes the expectation value of one or two Hermitian operators acting on one or two sites"


MPSCorrelation::usage="MPSCorrelation[mps,operator1,site1,operator2,site2]  Computes all the expectation values of operator1 and operator2 acting on all sites between site1 and site2"


MPSApproximate::usage="MPSApproximate[Bigmps,newBond,Options]  Computes a new MPS that approximates Bigmps but with a smaller bond dimension, specified by newBond. 
Options are:
     UseRandomState\[Rule]bool  (Default is True. It starts the approximation using a random state)
     Ansatz\[Rule]anMPS      (Default Empty List. If UseRandomState->True, it must contain a valid MPS to use as a starting Ansatz)
     Tolerance\[Rule]someNumber  (Default is some small quantity)
     Sweeps\[Rule]someNumber   (Default is 10)
     Verbose\[Rule]bool   (Default is False)

#### WARNING ####
This function is still being debugged"


MPSInitHamiltonian::usage="Ham=MPSInitHamiltonian[length] Initializes the object Ham to be the Hamiltonian matrix for MPS algorithms"


MPSHamiltonianAdd::usage="
MPSHamiltonianAdd[Ham,axis,strength,site]  Adds a field term to the Hamiltonian
MPSHamiltonianAdd[Ham,axis,strength,site1,site2]  Adds an interaction term to the Hamiltonian
Axis is either 1, 2, or 3 (which means X, Y, and Z Pauli matrices)
Ham must be initialized first.
As of now, MPS algorithms can only compute interactions that are of the type XX, YY, or ZZ with an arbitrary strength. This limitation will be lifted in the near future."



MPSMinimizeEnergy::usage="MPSMinimizeEnergy[mps,aHamiltonian,Options]

Computes the ground state of aHamiltonian using a DMRG sweeping technique. It returns the energy of the final state.
Options:

Sweeps\[Rule]MaxNumberOfSweeps
InteractionRange->integer   (The interaction range of aHamiltonian. This will be deprecated in future versions)
Tolerance\[Rule]tol    (Precision tolerance for stopping the algorithm)
MonitorEnergy\[Rule]bool  (Default True, it prints out not only the final energy of the state, but the dynamics of energy as a function of steps)
Verbose\[Rule]bool   (make the routine more talkative)"


TEBDProductState::usage="TEBDProductState[numSites,Options]
Creates a TEBD representation state. 
Options:
	Spin\[Rule]integer (Default is 2)
	BondDimension\[Rule]integer   (default is 20)
	Type\[Rule]{type,list}     (Creates several types of TEBDs: type=\"Random\" , \"Identity\", or \"Specified\".
						If Specified is chosen, then list must be a Table of dimensions {numSites,2}, containing in each element a pair {probAmplitude,phase}"


TEBDSave::usage="TEBDSave[MPS,filename] saves TEBD into ASCII files for later retrieval. Right now it produces a lot of files, to be changed in the future. filename needs to be a string"



TEBDRead::usage="TEBDRead[filename] reads the files produced by TEBDSave and returns a new TEBD object"


TEBDEvolve::usage="
TEBDEvolve[TEBD,OperatorList,repeatTimes]
Evolves the TEBD object with the given OperatorList (see example). It applies this list repeatTimes number of times, but this argument is optional and the default is one. It can take Verbose->True option to print out the steps."


TEBDExpectation::usage="
TEBDExpectation[TEBD,Operator,site]
TEBDExpectation[TEBD,Operator1,site1,Operator2,site2]

Computes the expectation value of one or two operators acting in one or two sites."


TEBDCorrelation::usage="
TEBDCorrelation[TEBD,Operator1,site1,Operator2,site2]

Computes the expectation value of one or two operators acting in one or two sites. The List version returns all the possible correlations between site1 and site2"


TEBDInitHamiltonian::usage="
Ham=TEBDInitHamiltonian[]   Initializes a Hamiltonian object Ham to an empty Hamiltonian"


TEBDInitEvolution::usage="
U=TEBDInitEvolution[]   Initializes an evolution operator object U."


TEBDEvolutionAddGate::usage="
TEBDEvolutionAddGate[U,operator,site]
TEBDEvolutionAddGate[U,operator12,site1,site2]

Adds an operator the the evolution operator U that acts either on a single site, or on two sites. Notice that operator must be a spin x spin matrix, or an arbitrary 2 spin operator (spin^2 x spin^2)."


TEBDHamiltonianAdd::usage="
TEBDHamiltonianAdd[Ham,operator,site]
TEBDHamiltonianAdd[Ham,operator1,site1,,operator2,site2]

Adds an operator the the Hamiltonian Ham that acts either on a single site, or on two sites. Notice that operator must be a spin x spin matrix.
An arbitrary 2 spin operator (spin^2 x spin^2) can also be passed as

TEBDHamiltonianAdd[Ham,operator,site1,site2]"


TEBDEvolutionFromHamiltonian::usage="TEBDEvolutionFromHamiltonian[Ham,time]
Creates an evolution operator List that is ready to be passed to TEBDEvolve. 
NOTICE this function is still being optimized. As of now, it only creates the list with the operators in the same order as given by the Hamiltonian. This means that the evolution list created will very likely create many errors. In the meantime, see examples on how to create a good list by hand."



TEBDEntropyList::usage="TEBDEntropyList[tebd]   Gives a list of all possible bipartite entropies of the TEBD state"


(* ::Section::Closed:: *)
(*MPS Algorithms*)


(* ::Subsection::Closed:: *)
(*Constants*)


DefaultSpinDimension=2;
MPSDefaultBond=10;
DefaultSweeps=10;
DefaultEnergyTolerance=10^(-4);
DefaultApproximationTolerance=10^(-10);
DefaultInteractionRange=1;
MPSMaxBond=100;
ForceUseInternalRoutine=True;


(* ::Subsection::Closed:: *)
(*Auxiliary Functions*)


ClearAll[LProduct,RProduct];
LProduct[A_,B_]:=Plus@@MapThread[ConjugateTranspose[#2].#1&,{A,B}];
LProduct[A_,B_,Lm_]:=Plus@@MapThread[ConjugateTranspose[#2].Lm.#1&,{A,B}];
RProduct[A_,B_]:=Plus@@MapThread[#1.ConjugateTranspose[#2]&,{A,B}];
RProduct[A_,B_,Rm_]:=Plus@@MapThread[#1.Rm.ConjugateTranspose[#2]&,{A,B}];


sigma[0]=SparseArray[{{1.0,0.0},{0.0,1.0}}];
sigma[1]=SparseArray[{{0.0,0.5},{0.5,0.0}}];
sigma[2]=SparseArray[{{0.0,-I 0.5},{I 0.5,0.0}}];
sigma[3]=SparseArray[{{0.5,0.0},{0.0,-0.5}}];


(* ::Subsection::Closed:: *)
(*Creation and basic manipulation*)


Clear[MPSProductState];
Options[MPSProductState]={Spin->DefaultSpinDimension,Bond->MPSDefaultBond,Type->"Random"};
MPSProductState[numTensors_,OptionsPattern[]]:=Module[{\[CapitalGamma]=Array[0&,{numTensors}],type=OptionValue[Type],coeffList,spin=OptionValue[Spin],\[Chi]=OptionValue[Bond]},
Switch[type,
"Identity",
coeffList=Table[Table[If[i==j==1,{1.0}~Join~Table[0.0,{m,1,spin-1}],Table[0.0,{m,spin}]],{k,1,numTensors}],{i,1,\[Chi]},{j,1,\[Chi]}];
"Decaying",
coeffList=Table[Table[If[i==j==1,Normalize[Table[Exp[-0.5 Log[20](m-1)/spin],{m,1,spin}]],0],{k,1,numTensors}],{i,1,\[Chi]},{j,1,\[Chi]}];
"Random",
coeffList=Table[Table[Normalize[RandomComplex[{-1-I,1+I},spin]],{k,1,numTensors}],{i,1,\[Chi]},{j,1,\[Chi]}];
];
{Table[SparseArray[Table[coeffList[[i,j,1,n]],{i,1,1},{j,1,\[Chi]}]],{n,1,spin}]}~Join~Table[
Table[SparseArray[Table[coeffList[[i,j,k,n]],{i,1,\[Chi]},{j,1,\[Chi]}]],{n,1,spin}],
{k,2,numTensors-1}]~Join~{Table[SparseArray[Table[coeffList[[i,j,numTensors,n]],{i,1,\[Chi]},{j,1,1}]],{n,1,spin}]}
];


Clear[MPSExpandBond];
SetAttributes[MPSExpandBond,HoldFirst];
MPSExpandBond[MPS_,new\[Chi]_]:=Module[
{old\[Chi]},
old\[Chi]=Max[Dimensions[#]&/@MPS];
If[old\[Chi]>new\[Chi],Print["WARNING: Trying to expand an MPS to a smaller bond dimension"];MPS,
Print["Grow \[Chi] to "<>ToString[new\[Chi]]];
{SparseArray[PadRight[#,{1,new\[Chi]}]&/@MPS[[1]]]}~Join~Table[SparseArray[PadRight[#,{new\[Chi],new\[Chi]}]&/@MPS[[M]]],{M,2,Length[MPS]-1}]~Join~{SparseArray[PadRight[#,{new\[Chi],1}]&/@MPS[[Length[MPS]]]]}
]
];


ClearAll[MPSOverlap];
SetAttributes[MPSOverlap,HoldAll];
MPSOverlap[mps1_,mps2_]:=Module[{ctemp=0,L},
(* This function also works with Fold, which is appealing because can also give FoldList:
Fold[LProduct[mps1[[#2]],mps2[[#2]],#1]&,LProduct[mps1[[1]],mps2[[1]]],Range[2,Length[mps1]]][[1,1]]
*)
L[1]=LProduct[mps1[[1]],mps2[[1]]];
L[n_]:=LProduct[mps1[[n]],mps2[[n]],L[n-1]];
Chop[L[Length[mps1]][[1,1]]]
];


ClearAll[MPSNormalize];
SetAttributes[MPSNormalize,HoldAll];
MPSNormalize[mps_]:=Module[{norm},
norm=Chop[MPSOverlap[mps,mps]];
mps=Chop[mps/Abs[norm]^(1/(2 Length[mps]))];
norm
];


ClearAll[MPSCanonizeSite];
Options[MPSCanonizeSite]={Direction->"Right",UseMatrix->True};
SetAttributes[MPSCanonizeSite,HoldAll];
MPSCanonizeSite[tensor_,matrix_,OptionsPattern[]]:=Module[{sense=OptionValue[Direction],usematrix=OptionValue[UseMatrix],numTensors,\[Chi]L,\[Chi]R,\[Chi],u,v,t,newTensor},(* Start by multiplying the tensor with the matrix from the previous site *)
If[sense=="Right",
If[usematrix,newTensor=tensor.matrix,newTensor=tensor];
{\[Chi]L,\[Chi]R}=Dimensions[newTensor[[1]]];
\[Chi]=Max[\[Chi]L,\[Chi]R];
(* SVD of the new tensor, putting [chiL, spin*chiR] *)
{u,v,t}=SingularValueDecomposition[Flatten[newTensor,{{2},{1,3}}]];
(* Prepare new right matrix *)
matrix=PadRight[u.v,{Min[\[Chi],\[Chi]L],Min[\[Chi],Length[t],\[Chi]L]}];
(* Form the new tensor with the first row of t^dagger *)
(Partition[ConjugateTranspose[t],{Min[\[Chi],Length[t],\[Chi]L],\[Chi]R}][[1,All]])
,
If[usematrix,newTensor=matrix.#&/@tensor,newTensor=tensor];
{\[Chi]L,\[Chi]R}=Dimensions[newTensor[[1]]];
\[Chi]=Max[\[Chi]L,\[Chi]R];
(* SVD of the new tensor, putting [chiL*spin, chiR] *)
{u,v,t}=SingularValueDecomposition[Flatten[newTensor,{{1,2},{3}}]];
(* Prepare new right matrix *)
matrix=PadRight[v.ConjugateTranspose[t],{Min[\[Chi],Length[u],\[Chi]R],Min[\[Chi],\[Chi]R]}];
(* Form the new tensor with the first column of u *)
(Partition[u,{\[Chi]L,Min[\[Chi],Length[u],\[Chi]R]}][[All,1]])
]
];


ClearAll[MPSCanonize];
Options[MPSCanonize]={Site->1};
SetAttributes[MPSCanonize,HoldAll];
MPSCanonize[mps_,OptionsPattern[]]:=Module[{site=OptionValue[Site],numTensors,xM},
numTensors=Length[mps];
(* First: Right normalization up to site *)
xM={{1.}};
Do[
mps[[s]]=SparseArray[MPSCanonizeSite[mps[[s]],xM]];
,{s,numTensors,site+1,-1}];
(* Now do LEFT normalization *)
xM={{1.}};
Do[
mps[[s]]=SparseArray[MPSCanonizeSite[mps[[s]],xM,Direction->"Left"]];
,{s,1,site-1}];
site
];


ClearAll[MPSCanonizationCheck];
Options[MPSCanonizationCheck]={Site->1};
SetAttributes[MPSCanonizationCheck,HoldAll];
MPSCanonizationCheck[mps_,OptionsPattern[]]:=Module[{checksite=OptionValue[Site],numTensors,spin,norm},
numTensors=Length[mps];
spin=Length[mps[[1]]];
norm=Sum[
Total[
Abs[
Sum[mps[[site,s]].ConjugateTranspose[mps[[site,s]]],{s,1,spin}]-IdentityMatrix[Length[mps[[site,1]]
]
]
],Infinity]
,{site,numTensors,checksite+1,-1}];
norm+=Sum[
Total[
Abs[
Sum[ConjugateTranspose[mps[[site,s]]].mps[[site,s]],{s,1,spin}]-IdentityMatrix[Length[mps[[site,1,1]]
]
]
],Infinity]
,{site,1,checksite-1}];
Chop[norm]==0
];


(* ::Subsection::Closed:: *)
(*Expectation values*)


MPSSiteOperator[tensor_,operator_]:=Module[{},
Flatten[Transpose[tensor,{3,1,2}].operator.Conjugate[tensor],{{3,1},{2,4}}]
];


MPSExpectation[mps_,operator_,site_Integer]:=Module[{L,R,numTensors},
numTensors=Length[mps];
R[numTensors+1]={{1}};
R[n_]:=R[n]=RProduct[mps[[n]],mps[[n]],R[n+1]];
L[0]={{1}};
L[n_]:=L[n]=LProduct[mps[[n]],mps[[n]],L[n-1]];
Chop[
Flatten[L[site-1]].MPSSiteOperator[mps[[site]],operator].Flatten[ R[site+1]]
]
];


MPSExpectation[mps_,operator_,site1_Integer,site2_Integer]:=Module[{L,R,numTensors},
MPSExpectation[mps,operator,site1,operator,site2]
];



MPSExpectation[mps_,operator1_,site1_Integer,operator2_,site2_Integer]:=Module[{L,R,numTensors},
Last[MPSCorrelation[mps,operator1,site1,operator2,site2]]
];


MPSCorrelation[mps_,operator1_,site1_Integer,operator2_,site2_Integer]:=Module[{L,R,Lo,numTensors,siteL,siteR,opL,opR,NormalOrder,corr},
numTensors=Length[mps];
Which[
0<site1<site2<numTensors+1,
{siteL,siteR}={site1,site2};
{opL,opR}={operator1,operator2};
NormalOrder=True;
,0<site2<site1<numTensors+1,
{siteL,siteR}={site2,site1};
{opL,opR}={operator2,operator1};
NormalOrder=False;
,True,
Print["MPSCorrelation called with wrong sites"];Abort[];
];
R[numTensors+1]={{1}};
R[n_]:=R[n]=RProduct[mps[[n]],mps[[n]],R[n+1]];
L[0]={{1}};
L[n_]:=L[n]=LProduct[mps[[n]],mps[[n]],L[n-1]];
Lo[siteL]=LProduct[opL.mps[[siteL]],mps[[siteL]],L[siteL-1]];
Lo[n_]:=Lo[n]=LProduct[mps[[n]],mps[[n]],Lo[n-1]];
corr=Table[
Chop[Flatten[Lo[site-1]].MPSSiteOperator[mps[[site]],opR].Flatten[ R[site+1]]]
,{site,siteL+1,siteR}];
If[NormalOrder,corr,Reverse[corr]]
];


(* ::Subsection::Closed:: *)
(*Disk interaction*)


ClearAll[MPSSave];
SetAttributes[MPSSave,HoldFirst];
MPSSave[MPS_,filename_]:=Module[{numSites,spin},
numSites=Length[MPS];
spin=Length[MPS[[1]]];
Export[filename<>".info",{numSites,spin},"Table"];
Do[
Do[
Export[filename<>"."<>ToString[n]<>"."<>ToString[s]<>".dat",MPS[[n,s]],"Table"]
,{s,1,spin}];
,{n,1,numSites}];
Run["tar -czf "<>filename<>".MPSz "<>filename<>".*.dat "<>filename<>".info"];
Run["rm "<>filename<>"*.dat "<>filename<>".info"]
];


ClearAll[MPSRead];
MPSRead[filename_]:=Module[{MPS,numSites,spin,\[Chi],info},
If[Length[FileNames[filename<>".MPSz"]]=!=1,Return[]];
Run["tar -zxf "<>filename<>".MPSz"];
If[Length[FileNames[filename<>".info"]]=!=1,Return[]];
info=Flatten[Import[filename<>".info","Table"]];
{numSites,spin}=info;
MPS={};
Do[
MPS=Append[MPS,SparseArray[Table[
If[Length[FileNames[filename<>"."<>ToString[n]<>"."<>ToString[s]<>".dat"]]=!=1,Print["Missing File "<>filename<>"."<>ToString[n]<>"."<>ToString[s]<>".dat"];Break[]];
ToExpression[
Import[filename<>"."<>ToString[n]<>"."<>ToString[s]<>".dat","Table"]
]
,{s,1,spin}]]];
,{n,1,numSites}];
Run["rm "<>filename<>"*.dat "<>filename<>".info"];
MPS
];


(* ::Subsection::Closed:: *)
(*Algorithms*)


(* ::Subsubsection::Closed:: *)
(*Approximation*)


SetAttributes[MPSApproximate,HoldAll];
Options[MPSApproximate]={UseRandomState->True,Ansatz->{},Tolerance->DefaultApproximationTolerance,Sweeps->DefaultSweeps,Verbose->False};
MPSApproximate[mps_,new\[Chi]_,OptionsPattern[]]:=Module[{sweeps=OptionValue[Sweeps],sweep=0,tol=OptionValue[Tolerance],L,R,numTensors,defineRight,defineLeft,\[Chi]R,\[Chi]L,new,canon,stillconverging=True,verbose=OptionValue[Verbose],message,info,success,newmps,overlapBIG,overlap,prevoverlap},
If[verbose,message=PrintTemporary["Preparing matrices"]];
(*Preparation assignments*)
(*MPSNormalize[mps];*)
numTensors=Length[mps];
overlapBIG=MPSOverlap[mps,mps];
newmps=MPSProductState[numTensors,Bond->new\[Chi]];
MPSCanonize[newmps];
prevoverlap=1+overlapBIG-2 Re[MPSOverlap[newmps,mps]];
(* These will be used to define left and right matrices *)
defineRight[temp_]:=Module[{fer},
ClearAll[R];
R[numTensors+1]={{1}};
R[n_]:=R[n]=RProduct[mps[[n]],newmps[[n]],R[n+1]];
];
defineLeft[temp_]:=Module[{fer},
ClearAll[L];
L[0]={{1}};
L[n_]:=L[n]=LProduct[mps[[n]],newmps[[n]],L[n-1]];
];
(* Start from the left, so prepare all right matrices *)
defineRight[1];
While[sweep<sweeps&&stillconverging,
(* Sweep to the right clears all left matrices and defines them one by one *)
defineLeft[1];
canon={{1.}};
Do[
If[verbose,NotebookDelete[message];message=PrintTemporary["Right sweep:"<>ToString[sweep]<>", site:"<>ToString[site]<>", overlap:"<>ToString[1+overlapBIG-2 Re[MPSOverlap[newmps,mps]]]];Pause[1]];
success=False;
newmps[[site]]=canon.#&/@newmps[[site]];
new=L[site-1].#.R[site+1]&/@mps[[site]];
newmps[[site]]=MPSCanonizeSite[new,canon,Direction->"Left",UseMatrix->False]; (* This routine changes canon *)
,{site,1,numTensors}];
sweep+=0.5;
(* Sweep to the left clears all right matrices and defines them one by one *)
defineRight[1];
canon={{1.}};
Do[
If[verbose,NotebookDelete[message];message=PrintTemporary["Left sweep:"<>ToString[sweep]<>", site:"<>ToString[site]<>", overlap:"<>ToString[1+overlapBIG-2 Re[MPSOverlap[newmps,mps]]]];Pause[1]];
success=False;
newmps[[site]]=newmps[[site]].canon;
new=L[site-1].#.R[site+1]&/@mps[[site]];
newmps[[site]]=MPSCanonizeSite[new,canon,UseMatrix->False];
,{site,numTensors,1,-1}];
sweep+=0.5;
overlap=1+overlapBIG-2 Re[MPSOverlap[newmps,mps]];
If[Abs[(overlap-prevoverlap)/overlap]<tol,stillconverging=False,prevoverlap=overlap];
];
If[verbose,NotebookDelete[message]];
newmps
];


(* ::Subsubsection::Closed:: *)
(*Energy minimization*)


ClearAll[MPSMinimizeEnergy];
Options[MPSMinimizeEnergy]={Sweeps->DefaultSweeps,Tolerance->DefaultEnergyTolerance,MonitorEnergy->False,Verbose->False,Debug->False};
SetAttributes[MPSMinimizeEnergy,HoldFirst];
MPSMinimizeEnergy[mps_,HObject_,OptionsPattern[]]:=
Module[{energy,prevEnergy=0,energyList={},sweeps=OptionValue[Sweeps],sweep=0,IntRange,monitorenergy=OptionValue[MonitorEnergy],tol=OptionValue[Tolerance],L,R,fieldL,fieldR,operatorsR,operatorsL,interactionsL,interactionsR,Heff,numTensors,defineRight,defineLeft,\[Chi]R,\[Chi]L,new,canon,stillconverging=True,verbose=OptionValue[Verbose],message,info,success,debug=OptionValue[Debug],debuglist={},debugMSG,HMatrix},
If[verbose,message=PrintTemporary["Preparing matrices"]];
(*Preparation assignments*)
(*MPSNormalize[mps];*)
{IntRange,HMatrix}=HObject;
MPSCanonize[mps];
numTensors=Length[mps];
(* These will be used to define left and right matrices *)
defineRight[temp_]:=Module[{fer},
ClearAll[R,fieldR,operatorsR,interactionsR];
R[numTensors+1]={{1}};
R[n_]:=R[n]=RProduct[mps[[n]],mps[[n]],R[n+1]];
fieldR[numTensors+1]=Table[{{0}},{\[Alpha],1,3}];
fieldR[n_]:=fieldR[n]=Table[
RProduct[mps[[n]],mps[[n]],fieldR[n+1][[\[Alpha]]]]+HMatrix[[\[Alpha],n,n]]RProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],R[n+1]],{\[Alpha],1,3}];
operatorsR[numTensors+2]=Table[{},{\[Alpha],1,3}];
operatorsR[numTensors+1]=Table[{},{\[Alpha],1,3}];
operatorsR[n_]:=operatorsR[n]=Table[
Prepend[RProduct[mps[[n]],mps[[n]],#]&/@(Take[operatorsR[n+1][[\[Alpha]]],Min[IntRange-1,Length[operatorsR[n+1][[\[Alpha]]]]]]),RProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],R[n+1]]],{\[Alpha],1,3}];
interactionsR[numTensors+1]=Table[{{0}},{\[Alpha],1,3}];
interactionsR[n_]:=interactionsR[n]=Table[
RProduct[mps[[n]],mps[[n]],interactionsR[n+1][[\[Alpha]]]]+(HMatrix[[\[Alpha],n,n+1;;Min[n+IntRange,numTensors]]].(RProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],#]&/@operatorsR[n+1][[\[Alpha]]])),{\[Alpha],1,3}];
];
defineLeft[temp_]:=Module[{fer},
ClearAll[L,fieldL,operatorsL,interactionsL];
L[0]={{1}};
L[n_]:=L[n]=LProduct[mps[[n]],mps[[n]],L[n-1]];
fieldL[0]=Table[{{0}},{\[Alpha],1,3}];
fieldL[n_]:=fieldL[n]=Table[LProduct[mps[[n]],mps[[n]],fieldL[n-1][[\[Alpha]]]]+HMatrix[[\[Alpha],n,n]]LProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],L[n-1]],{\[Alpha],1,3}];
operatorsL[-1]=Table[{},{\[Alpha],1,3}];
operatorsL[0]=Table[{},{\[Alpha],1,3}];
operatorsL[n_]:=operatorsL[n]=Table[
Prepend[LProduct[mps[[n]],mps[[n]],#]&/@Take[operatorsL[n-1][[\[Alpha]]],Min[IntRange-1,Length[operatorsL[n-1][[\[Alpha]]]]]],LProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],L[n-1]]],{\[Alpha],1,3}];
interactionsL[0]=Table[{{0}},{\[Alpha],1,3}];
interactionsL[n_]:=interactionsL[n]=Table[LProduct[mps[[n]],mps[[n]],interactionsL[n-1][[\[Alpha]]]]+(Reverse[HMatrix[[\[Alpha],Max[n-IntRange,1];;n-1,n]]].(LProduct[sigma[\[Alpha]].mps[[n]],mps[[n]],#]&/@operatorsL[n-1][[\[Alpha]]])),{\[Alpha],1,3}];
];
(* Start from the left, so prepare all right matrices *)

debugMSG[text_,s_,Am_,energ_]:={{text<>", Site:"<>ToString[s] ,
"\nMPS:",Am[[1]],Am[[2]],
"\nInt: ",interactionsL[s-1][[3]],interactionsR[s+1][[3]],
"\nh: ",fieldL[s-1][[1]],fieldR[s+1][[1]],
"\nOps: ",operatorsL[s-1][[3]],operatorsR[s+1][[3]],
"\nelsH: ",HMatrix[[All,Max[1,s-IntRange];;s,Max[1,s-IntRange];;Min[numTensors,s+IntRange]]][[1]],
"\nHeff: ",Sort[ Normal[Diagonal[MPSEffectiveHam[interactionsL[s-1],interactionsR[s+1],fieldL[s-1],fieldR[s+1],operatorsL[s-1],operatorsR[s+1],HMatrix[[All,Max[1,s-IntRange];;s,Max[1,s-IntRange];;Min[numTensors,s+IntRange]]]] ] ]] ,
"\nEnergy: ",energ  }};

defineRight[1];
While[sweep<sweeps&&stillconverging,
(* Sweep to the right clears all left matrices and defines them one by one *)
defineLeft[1];
canon={{1.}};
Do[
If[verbose,NotebookDelete[message];message=PrintTemporary["Right sweep:"<>ToString[sweep]<>", site:"<>ToString[site]<>", Energy:"<>ToString[energy]]];
If[debug,debuglist=debuglist~Join~debugMSG["Before (Right)",site,canon.#&/@ mps[[site]],energy]]; (* *)
success=False;
While[!success,
{energy,new,info}=FindGroundMPSSite[canon.#&/@mps[[site]],interactionsL[site-1],interactionsR[site+1],fieldL[site-1],fieldR[site+1],operatorsL[site-1],operatorsR[site+1],HMatrix[[All,Max[1,site-IntRange];;site,Max[1,site-IntRange];;Min[numTensors,site+IntRange]]]];
If[ForceUseInternalRoutine,
success=True;,
success=(IsLinkActive[]===1);
If[!success,Print["Fallen link..."];ClearLink[link];EstablishLink[link];Print["And we're back."]];
];
];
mps[[site]]=MPSCanonizeSite[new,canon,Direction->"Left",UseMatrix->False]; (* This routine changes canon *)
If[monitorenergy,energyList=energyList~Join~{{energy/(Conjugate[Flatten[new]].Flatten[new]),info}}];
If[debug,debuglist=debuglist~Join~debugMSG["After (Right)",site,mps[[site]].canon,energy/(Conjugate[Flatten[new]].Flatten[new])]];
,{site,1,numTensors}];
sweep+=0.5;
(* Sweep to the left clears all right matrices and defines them one by one *)
defineRight[1];
canon={{1.}};
Do[
If[verbose,NotebookDelete[message];message=PrintTemporary["Left sweep:"<>ToString[sweep]<>", site:"<>ToString[site]<>", Energy:"<>ToString[energy]]];
If[debug,debuglist=debuglist~Join~debugMSG["Before (Left)",site,mps[[site]].canon  ,energy]]; (* *)
success=False;
While[!success,
{energy,new,info}=FindGroundMPSSite[mps[[site]].canon,interactionsL[site-1],interactionsR[site+1],fieldL[site-1],fieldR[site+1],operatorsL[site-1],operatorsR[site+1],HMatrix[[All,Max[1,site-IntRange];;site,Max[1,site-IntRange];;Min[numTensors,site+IntRange]]]];
If[ForceUseInternalRoutine,
success=True;,
success=(IsLinkActive[]===1);
If[!success,Print["Fallen link..."];ClearLink[link];EstablishLink[link];Print["And we're back."]];
];
];
mps[[site]]=MPSCanonizeSite[new,canon,UseMatrix->False];
If[monitorenergy,energyList=energyList~Join~{{energy/(Conjugate[Flatten[new]].Flatten[new]),info}}];
If[debug,debuglist=debuglist~Join~debugMSG["After (Left)",site,canon.#&/@mps[[site]],energy/(Conjugate[Flatten[new]].Flatten[new])]]; (* .canon  *)
,{site,numTensors,1,-1}];
sweep+=0.5;
If[Abs[(energy-prevEnergy)/energy]<tol,stillconverging=False,prevEnergy=energy];
];
If[verbose,NotebookDelete[message]];
If[Total[Abs[#[[2]]]&/@energyList]>0,
Print["Arpack error reported:"];
Print[Cases[#[[2]]&/@energyList,Except[0]]]
];
ClearAll[prevEnergy,sweeps,sweep,IntRange,tol,L,R,fieldL,fieldR,operatorsR,operatorsL,interactionsL,interactionsR,Heff,numTensors,defineRight,defineLeft,\[Chi]R,\[Chi]L,new,canon,stillconverging,verbose,message,info,success];
If[debug,{#[[1]]&/@energyList,debuglist},If[monitorenergy,Chop[energyList],Chop[energy] ] ]
];


(* ::Subsubsection::Closed:: *)
(*Internal use of ARPACK routine*)


(* Define internal versions *)
ClearAll[MPSEffectiveSingleHam];
SetAttributes[MPSEffectiveSingleHam,HoldAll];
MPSEffectiveSingleHam[L_,R_,op_]:=(*KroneckerProduct[op,SparseArray[Flatten[Transpose[{L},{3,1,2}].{R},{{4,1},{3,2}}]] *)
KroneckerProduct[op,L,Transpose[R]];


ClearAll[MPSEffectiveHam];
SetAttributes[MPSEffectiveHam,HoldAll];
MPSEffectiveHam[interactionsL_,interactionsR_,fieldL_,fieldR_,operatorsL_,operatorsR_,Hmatrix_]:=Module[{Htemp=0,\[Chi]L,\[Chi]R,intrangeL,intrangeR,L,R},
(* Preparation of internal matrices and constants *)
intrangeL=Length[operatorsL[[1]]];
intrangeR=Length[operatorsR[[1]]];
\[Chi]L=Length[fieldL[[1]]];
\[Chi]R=Length[fieldR[[1]]];
L=SparseArray[IdentityMatrix[\[Chi]L]];
R=SparseArray[IdentityMatrix[\[Chi]R]];
(* First compute the contribution from the fields h *)
Htemp= Sum[MPSEffectiveSingleHam[L,fieldR[[\[Alpha]]],sigma[0]],{\[Alpha],1,3}];
Htemp+=Sum[MPSEffectiveSingleHam[fieldL[[\[Alpha]]],R,sigma[0]],{\[Alpha],1,3}];
(* Now compute the contribution from the interactions contained in the right and left blocks *)
Htemp+=Sum[MPSEffectiveSingleHam[interactionsL[[\[Alpha]]],R,sigma[0]],{\[Alpha],1,3}];
Htemp+= Sum[MPSEffectiveSingleHam[L,interactionsR[[\[Alpha]]],sigma[0]],{\[Alpha],1,3}];
(* Now compute the interactions between the left and right blocks that do not involve the spin at the site *)
Htemp+=Sum[Sum[Sum[MPSEffectiveSingleHam[operatorsL[[\[Alpha],x]],operatorsR[[\[Alpha],y]],If[(x+y-1)<Max[intrangeL,intrangeR],Hmatrix[[\[Alpha],intrangeL+1-x,intrangeL+1+y]],0]sigma[0]],{\[Alpha],1,3}],{y,1,intrangeR}],{x,1,intrangeL}]; 
(* Now add the term with the field at the site *)
Htemp+= Sum[MPSEffectiveSingleHam[L,R,Hmatrix[[\[Alpha],intrangeL+1,intrangeL+1]]sigma[\[Alpha]] ],{\[Alpha],1,3}];
(* Finally, the interactions between the site and the left and right blocks *) Htemp+=Sum[Sum[MPSEffectiveSingleHam[operatorsL[[\[Alpha],x]]Hmatrix[[\[Alpha],intrangeL+1-x,intrangeL+1]],R,sigma[\[Alpha]]],{\[Alpha],1,3}],{x,1,intrangeL}]; 
Htemp+=Sum[Sum[MPSEffectiveSingleHam[L,operatorsR[[\[Alpha],x]],Hmatrix[[\[Alpha],intrangeL+1,intrangeL+1+x]]sigma[\[Alpha]]],{\[Alpha],1,3}],{x,1,intrangeR}]; 
Return[Htemp]
];


ClearAll[FindGroundMPSSiteManual];
SetAttributes[FindGroundMPSSiteManual,HoldAll];
FindGroundMPSSiteManual[A_,DLeft_,DRight_,hLeft_,hRight_,vLeft_,vRight_,Ham_]:=Module[{H,sol,\[Chi]L,\[Chi]R,spin},
(* Print["Link could not be established, reverting to internal routines..."];*)
H=MPSEffectiveHam[DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham];
sol=Eigensystem[-H,1,Method->{"Arnoldi",MaxIterations->10^5,Criteria->RealPart}];
{spin,\[Chi]L,\[Chi]R}=Dimensions[A];
{Chop[-sol[[1,1]]],Chop[-Partition[Partition[sol[[2,1]],\[Chi]R],\[Chi]L]],0}
];


(*This clears a link and associated variables*)
ClearLink[LINK_]:=Module[{},
Uninstall[LINK];
ClearAll[FindGroundMPSSite,IsLinkActive];
];


(* Now try to establish the link *)
EstablishLink[LINK_]:=(
ClearLink[LINK];
ClearAll[FindGroundMPSSite,IsLinkActive];
If[Length[Links["*arpackformps"]]==0,LINK=Install["arpackformps_"<>$SystemID]];
If[LINK===$Failed||ForceUseInternalRoutine,
SetAttributes[FindGroundMPSSite,HoldAll];
FindGroundMPSSite[A_,DLeft_,DRight_,hLeft_,hRight_,vLeft_,vRight_,Ham_]:=FindGroundMPSSiteManual[A,DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham]
]
);


If[!ForceUseInternalRoutine,
EstablishLink[link];,
ClearAll[FindGroundMPSSite];
SetAttributes[FindGroundMPSSite,HoldAll];
FindGroundMPSSite[A_,DLeft_,DRight_,hLeft_,hRight_,vLeft_,vRight_,Ham_]:=FindGroundMPSSiteManual[A,DLeft,DRight,hLeft,hRight,vLeft,vRight,Ham]
]



(* ::Subsection::Closed:: *)
(*Hamiltonian generation*)


SetAttributes[MPSInitHamiltonian,HoldFirst];
MPSInitHamiltonian[length_]:=({0,Table[Table[0.0,{n,1,length},{m,1,length}],{\[Alpha],1,3}]})


SetAttributes[MPSHamiltonianAdd,HoldFirst];
MPSHamiltonianAdd[Ham_,axis_,strength_,site_]:=Module[{length,dims1,dims2,Numaxes,components,intRange},
components=Length[Ham];
If[components!=2,Print["Hamiltonian seems to be malformed. Please restart it"];Return[];];
intRange=Ham[[1]];
{Numaxes,dims1,dims2}=Dimensions[Ham[[2]]];
length=dims1;
If[Numaxes!=3||dims1!=dims2,Print["Hamiltonian seems to be malformed. Please restart it"];Return[];];
If[1>site||site>length,Print["Site needs to be between 1 and Length."];Return[];];
(Ham[[2,axis,site,site]]=strength)];


MPSHamiltonianAdd[Ham_,axis_,strength_,site1_,site2_]:=Module[{length,dims1,dims2,Numaxes,components,intRange},
components=Length[Ham];
If[components!=2,Print["Hamiltonian seems to be malformed. Please restart it"];Return[];];
intRange=Ham[[1]];
{Numaxes,dims1,dims2}=Dimensions[Ham[[2]]];
If[Numaxes!=3||dims1!=dims2,Print["Hamiltonian seems to be malformed. Please restart it"];Return[];];
length=dims1;
If[1>site1||site1>length||1>site2||site2>length,Print["Site needs to be between 1 and Length."];Return[];];
If[Abs[site1-site2]>intRange,Ham[[1]]=Abs[site1-site2]];
(Ham[[2,axis,site1,site2]]=strength);
(Ham[[2,axis,site2,site1]]=strength);
];


(* ::Section:: *)
(*TEBD Algorithms*)


(* ::Subsection::Closed:: *)
(*Constants*)


TOLERANCE=10^(-10); (* This is the tolerance of PseudoInverse *)


(* ::Subsection::Closed:: *)
(*States*)


Options[TEBDProductState]={Spin->2,BondDimension->20,Type->{"Random",{}}};
TEBDProductState[numSites_,OptionsPattern[]]:=Module[{\[CapitalGamma]=Array[0&,{numSites}],\[Lambda]=Array[0&,{numSites}],tempType=OptionValue[Type],type,phaseList,spin=OptionValue[Spin],\[Chi]=OptionValue[BondDimension]},
{type,phaseList}=tempType;
Switch[type,
"Specified",
If[Dimensions[phaseList]!={numSites,2},Print["TEBDProductState called with wrong list of phases."];Abort[]],
"Identity",
phaseList=Table[{0.0,0.0},{k,1,numSites}];,
"Random",
phaseList=RandomReal[{0,2\[Pi]},{numSites,2}];
];
Table[{
Table[SparseArray[{{i_,j_}/;(i==j==1)->Chop[(2-n)Cos[1.0 phaseList[[k,1]]]+(n-1)Exp[I phaseList[[k,2]]]Sin[1.0 phaseList[[k,1]]]]},{\[Chi],\[Chi]}],{n,1,spin}],
SparseArray[{i_}/;i<=1->1.0 ,{\[Chi]}]},
{k,1,numSites}]
]


(* ::Subsection::Closed:: *)
(*Disk I/O*)


SetAttributes[TEBDSave,HoldFirst];
TEBDSave[TEBD_,filename_]:=Module[{numSites,spin,\[Chi]},
numSites=Length[TEBD];
spin=Length[TEBD[[1,1]]];
\[Chi]=Length[TEBD[[1,1,1]]];
Export[filename<>".info",{numSites,spin,\[Chi]},"Table"];
Do[
Do[
Export[filename<>".G."<>ToString[n]<>"."<>ToString[s]<>".dat",TEBD[[n,1,s]],"Table"]
,{s,1,spin}];
Export[filename<>".L."<>ToString[n]<>".dat",TEBD[[n,2]],"Table"];
,{n,1,numSites}]
]


SetAttributes[TEBDRead,HoldFirst];
TEBDRead[TEBD_,filename_]:=Module[{numSites,spin,\[Chi],info},
info=Flatten[Import[filename<>".info","Table"]];
{numSites,spin,\[Chi]}=info;
TEBD=TEBDProductState[numSites,Spin->spin,BondDimension->\[Chi]];
Do[
TEBD[[n,1]]=SparseArray[Table[
Import[filename<>".G."<>ToString[n]<>"."<>ToString[s]<>".dat","Table"]
,{s,1,spin}]];
TEBD[[n,2]]=SparseArray[Flatten[Import[filename<>".L."<>ToString[n]<>".dat","Table"]]];
,{n,1,numSites}]
]


(* ::Subsection::Closed:: *)
(*Evolution *)


(* ::Subsubsection::Closed:: *)
(*Specific routines*)


Swap={{1,0,0,0},{0,0,1,0},{0,1,0,0},{0,0,0,1}};


SetAttributes[TEBD2QGate,HoldFirst];
TEBD2QGate[TEBD_,leftSite_,H_]:=Module[{\[Chi],\[Chi]R,\[Chi]L,Hint,\[Lambda]LD,\[Lambda]CD,\[Lambda]RD,\[CapitalSigma],X,YT,\[Lambda]temp,\[Lambda]RDi,\[Lambda]LDi,norm,numSites,rightSite,spin},
(*   Print[step=0.1,\[Lambda]L,\[Lambda]C,\[Lambda]R]; Preliminary: prepare matrices for fast products *)
numSites=Length[TEBD];
rightSite=leftSite+1;
If[(1<=leftSite<numSites)!=True,Print["TEBD2QGate called with wrong site parameter."];Abort[]];
{spin,\[Chi]L,\[Chi]R}=Dimensions[TEBD[[leftSite,1]]];
\[Chi]=Max[\[Chi]R,\[Chi]L];
If[Dimensions[H]!={spin^2,spin^2},Print["TEBD2QGate called with wrong operator dimensions."];Abort[]];
Hint=((Partition[#,spin]&/@H));
If[leftSite==1,
\[Lambda]LD=SparseArray[{i_,j_}/;i==j==1->1.0 ,{\[Chi],\[Chi]}];,
\[Lambda]LD=SparseArray[DiagonalMatrix[TEBD[[leftSite-1,2]]]];
];
\[Lambda]CD=SparseArray[DiagonalMatrix[TEBD[[leftSite,2]]]];
\[Lambda]RD=SparseArray[DiagonalMatrix[TEBD[[rightSite,2]]]];
(**O Print["Check 2-1,"];Print[Dimensions[#]&/@{\[Lambda]RD,\[Lambda]CD,TEBD[[leftSite,1]],TEBD[[rightSite,1]]}];DD**)
\[CapitalSigma]=Flatten[Flatten[\[Lambda]LD.Transpose[TEBD[[leftSite,1]],{2,1,3}].\[Lambda]CD.Transpose[TEBD[[rightSite,1]],{2,1,3}].\[Lambda]RD,{{1},{4},{2,3}}].Hint,{{3,1},{4,2}}];
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
(**O Print["Check 2-1B: "<>ToString[leftSite]<>", "<>ToString[Dimensions[#]&/@{X,\[Lambda]temp,YT}]];DD**)
\[Lambda]CD=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); 
norm=Total[\[Lambda]CD^2];
TEBD[[leftSite,2]]=\[Lambda]CD/Sqrt[norm];
(**O Print["Check 2-2"];DD**)
\[Lambda]LDi=SparseArray[PseudoInverse[Normal[\[Lambda]LD],Tolerance->TOLERANCE]];
TEBD[[leftSite,1]]=SparseArray[Chop[(\[Lambda]LDi.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(**O Print["Check 2-3"];DD**)
(**O Print["Check 2-3 B: site "<>ToString[rightSite]<>", "<>ToString[Dimensions[#]&/@{YT,(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]]),Transpose[{PadRight[Conjugate[#],\[Chi]]}]&/@YT}]];DD**)
\[Lambda]RDi=SparseArray[PseudoInverse[Normal[\[Lambda]RD],Tolerance->TOLERANCE]];
TEBD[[rightSite,1]]=SparseArray[Chop[(#.\[Lambda]RDi&)/@(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]])]];
1-norm
]


SetAttributes[TEBDGate,HoldFirst];
TEBDGate[TEBD_,SiteOperator_]:=Module[{error,numSites,numOfSwaps,ActualLeft,ActualRight,leftSite,rightSite,H},
(*   Print[step=0.1,\[Lambda]L,\[Lambda]C,\[Lambda]R]; Preliminary: prepare matrices for fast products *)
numSites=Length[TEBD];
H=SiteOperator[[1]];
leftSite=SiteOperator[[2]];
rightSite=SiteOperator[[3]];
If[!(1<=leftSite<=numSites)||!(1<=rightSite<=numSites),Print["TEBDGate[left,right] called with wrong site parameter:",leftSite,",",rightSite];Abort[]];
error=0;
If[leftSite==rightSite,
(* Equal sites *)
error=TEBD1QGate[TEBD,leftSite,H];
, 
(* Different sites *)
ActualLeft=Min[leftSite,rightSite];
ActualRight=Max[leftSite,rightSite];
numOfSwaps=ActualRight-ActualLeft-If[leftSite<rightSite,1,0];
error=0;
Do[error+=TEBD2QGate[TEBD,ActualRight-n,Swap],{n,1,numOfSwaps,1}];
error+=TEBD2QGate[TEBD,ActualLeft,H];
Do[error+=TEBD2QGate[TEBD,ActualRight-n,Swap],{n,numOfSwaps,1,-1}];
];
error
]


SetAttributes[TEBD3QGate,HoldFirst];
TEBD3QGate[TEBD_,centerSite_,H_]:=Module[{spin,\[Chi],\[Chi]R,\[Chi]L,\[Lambda]LD,\[Lambda]1D,\[Lambda]2D,\[Lambda]RD,\[CapitalSigma],X,YT,\[Lambda]temp,\[Lambda]RDi,\[Lambda]LDi,\[Lambda]1Di,norm1,norm2,numSites,Hint,leftSite,rightSite},
 (*  Print[step=0];  Preliminary: prepare matrices for fast products *)
leftSite=centerSite-1;
rightSite=centerSite+1;
numSites=Length[TEBD];
If[(1<centerSite<numSites)!=True,Print["TEBD3QGate called with wrong site,"<>ToString[centerSite]];Abort[]];
{spin,\[Chi]L,\[Chi]R}=Dimensions[TEBD[[leftSite,1]]];
\[Chi]=Max[\[Chi]R,\[Chi]L];
If[Dimensions[H]!={spin^3,spin^3},Print["TEBD3QGate called with wrong operator dimensions."];Abort[]];
Hint=(Partition[#,spin]&/@(Partition[#,spin]&/@H));
If[leftSite==1,
\[Lambda]LD=SparseArray[{i_,j_}/;i==j==1->1.0 ,{\[Chi],\[Chi]}];,
\[Lambda]LD=SparseArray[DiagonalMatrix[TEBD[[leftSite-1,2]]]];
];
\[Lambda]1D=SparseArray[DiagonalMatrix[TEBD[[leftSite,2]]]];
\[Lambda]2D=SparseArray[DiagonalMatrix[TEBD[[centerSite,2]]]];
\[Lambda]RD=SparseArray[DiagonalMatrix[TEBD[[rightSite,2]]]];
(**O Print["Check 3-1"];DD**)
(* First make 3 spin matrix  *)
\[CapitalSigma]=Flatten[Flatten[\[Lambda]LD.Transpose[TEBD[[leftSite,1]],{2,1,3}].\[Lambda]1D.Transpose[TEBD[[centerSite,1]],{2,1,3}].\[Lambda]2D.Transpose[TEBD[[rightSite,1]],{2,1,3}].\[Lambda]RD,{{1},{5},{2,3,4}}].Hint,{{3,1},{4,5,2}}];
(*  Perform SVD of the three spin matrix as 1-2 spins *)
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
\[Lambda]temp=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); (* Truncated values *)
norm1=Total[\[Lambda]temp^2];
TEBD[[leftSite,2]]=SparseArray[\[Lambda]temp/Sqrt[norm1]];
\[Lambda]1D=DiagonalMatrix[TEBD[[leftSite,2]]];
(* Extract the tensors of the left spin *)
(**O Print["Check 3-2"];DD**)
\[Lambda]LDi=SparseArray[PseudoInverse[Normal[\[Lambda]LD],Tolerance->TOLERANCE]];
TEBD[[leftSite,1]]=SparseArray[Chop[(\[Lambda]LDi.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(* Now rearrange right matrix to represent 2 spins as 1-1, multiplying before by the new \[Lambda]1 *)
(**O Print["Check 3-3"];DD**)
\[CapitalSigma]=Flatten[\[Lambda]1D.Transpose[Partition[Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]],spin],{2,3,1,4}],{{2,1},{3,4}}];
(* And get the SVD as a normal case *)
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
(**O Print[ToString[Dimensions[#]&/@{X,\[Lambda]temp,YT}]];DD**)
\[Lambda]temp=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); (* Truncated values *)
norm2=Total[\[Lambda]temp^2];
TEBD[[centerSite,2]]=SparseArray[\[Lambda]temp/Sqrt[norm2]];
(* The center tensor we obtain taking out \[Lambda]1 *) 
\[Lambda]1Di=SparseArray[PseudoInverse[Normal[\[Lambda]1D],Tolerance->TOLERANCE]];
(**O Print["Check 3-4"];DD**)
TEBD[[centerSite,1]]=SparseArray[Chop[(\[Lambda]1Di.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(* The right tensor needs the inverse right lambda to fit *)
(**O Print["Check 3-5"];DD**)
(**O Print["Check 3-5 B: site "<>ToString[rightSite]<>", "<>ToString[Dimensions[#]&/@{YT,(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]]),Transpose[{PadRight[Conjugate[#],\[Chi]]}]&/@YT}]];DD**)
\[Lambda]RDi=SparseArray[PseudoInverse[Normal[\[Lambda]RD],Tolerance->TOLERANCE]];
TEBD[[rightSite,1]]=SparseArray[Chop[(#.\[Lambda]RDi&)/@(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]])]];
2-norm1-norm2
]


SetAttributes[TEBD4QGate,HoldFirst];
TEBD4QGate[TEBD_,leftSite_,H_]:=Module[{spin,\[Chi],\[Chi]R,\[Chi]L,\[Lambda]LD,\[Lambda]1D,\[Lambda]2D,\[Lambda]3D,\[Lambda]2Di,\[Lambda]RD,\[CapitalSigma],X,YT,\[Lambda]temp,\[Lambda]RDi,\[Lambda]LDi,\[Lambda]1Di,norm1,norm2,norm3,numSites,Hint,centerSite1,centerSite2,rightSite},
 (*  Print[step=0];  Preliminary: prepare matrices for fast products *)
centerSite1=leftSite+1;
centerSite2=leftSite+2;
rightSite=leftSite+3;
numSites=Length[TEBD];
If[(1<=leftSite<rightSite<=numSites)!=True,Print["TEBD4QGate called with wrong site,"<>ToString[leftSite]];Abort[]];
{spin,\[Chi],\[Chi]R}=Dimensions[TEBD[[centerSite1,1]]];
If[Dimensions[H]!={spin^4,spin^4},Print["TEBD4QGate called with wrong operator dimensions."];Abort[]];
Hint=Partition[#,spin]&/@(Partition[#,spin]&/@(Partition[#,spin]&/@H));
If[leftSite==1,
\[Lambda]LD=SparseArray[{i_,j_}/;i==j==1->1.0 ,{\[Chi],\[Chi]}];,
\[Lambda]LD=SparseArray[DiagonalMatrix[TEBD[[leftSite-1,2]]]];
];
\[Lambda]1D=SparseArray[DiagonalMatrix[TEBD[[leftSite,2]]]];
\[Lambda]2D=SparseArray[DiagonalMatrix[TEBD[[centerSite1,2]]]];
\[Lambda]3D=SparseArray[DiagonalMatrix[TEBD[[centerSite2,2]]]];
\[Lambda]RD=SparseArray[DiagonalMatrix[TEBD[[rightSite,2]]]];
(**O Print["Check 3-1"];DD**)
(* First make 3 spin matrix  *)
\[CapitalSigma]=Flatten[Flatten[\[Lambda]LD.Transpose[TEBD[[leftSite,1]],{2,1,3}].\[Lambda]1D.Transpose[TEBD[[centerSite1,1]],{2,1,3}].\[Lambda]2D.Transpose[TEBD[[centerSite2,1]],{2,1,3}].\[Lambda]3D.Transpose[TEBD[[rightSite,1]],{2,1,3}].\[Lambda]RD,{{1},{6},{2,3,4,5}}].Hint,{{3,1},{4,5,2}}];
(*  Perform SVD of the three spin matrix as 1-2 spins *)
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
\[Lambda]temp=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); (* Truncated values *)
norm1=Total[\[Lambda]temp^2];
TEBD[[leftSite,2]]=SparseArray[\[Lambda]temp/Sqrt[norm1]];
\[Lambda]1D=DiagonalMatrix[TEBD[[leftSite,2]]];
(* Extract the tensors of the left spin *)
(**O Print["Check 3-2"];DD**)
\[Lambda]LDi=SparseArray[PseudoInverse[Normal[\[Lambda]LD],Tolerance->TOLERANCE]];
TEBD[[leftSite,1]]=SparseArray[Chop[(\[Lambda]LDi.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(* Now rearrange right matrix to represent 2 spins as 1-1, multiplying before by the new \[Lambda]1 *)
(**O Print["Check 3-3"];DD**)
\[CapitalSigma]=Flatten[\[Lambda]1D.Transpose[Partition[Partition[ConjugateTranspose[YT],{\[Chi],spin^2 \[Chi]}][[1,All]],spin],{2,1,3}],{{2,1},{3}}];
(* And get the SVD as a normal case *)
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
(**O Print[ToString[Dimensions[#]&/@{X,\[Lambda]temp,YT}]];DD**)
\[Lambda]temp=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); (* Truncated values *)
norm2=Total[\[Lambda]temp^2];
TEBD[[centerSite1,2]]=SparseArray[\[Lambda]temp/Sqrt[norm2]];
(* The center tensor we obtain taking out \[Lambda]1 *) 
\[Lambda]1Di=SparseArray[PseudoInverse[Normal[\[Lambda]1D],Tolerance->TOLERANCE]];
(**O Print["Check 3-4"];DD**)
TEBD[[centerSite1,1]]=SparseArray[Chop[(\[Lambda]1Di.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(* The right tensor needs the inverse right lambda to fit *)
(**O Print["Check 3-5"];DD**)
(**O Print["Check 3-5 B: site "<>ToString[rightSite]<>", "<>ToString[Dimensions[#]&/@{YT,(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]]),Transpose[{PadRight[Conjugate[#],\[Chi]]}]&/@YT}]];DD**)
\[CapitalSigma]=Flatten[\[Lambda]2D.Transpose[Partition[Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]],spin],{2,3,1,4}],{{2,1},{3,4}}];
(* And get the SVD as a normal case *)
{X,\[Lambda]temp,YT}=SingularValueDecomposition[\[CapitalSigma]]; 
(**O Print[ToString[Dimensions[#]&/@{X,\[Lambda]temp,YT}]];DD**)
\[Lambda]temp=(PadRight[Diagonal[\[Lambda]temp],\[Chi]]); (* Truncated values *)
norm3=Total[\[Lambda]temp^2];
TEBD[[centerSite2,2]]=SparseArray[\[Lambda]temp/Sqrt[norm2]];
(* The center tensor we obtain taking out \[Lambda]1 *) 
\[Lambda]2Di=SparseArray[PseudoInverse[Normal[\[Lambda]2D],Tolerance->TOLERANCE]];
(**O Print["Check 3-4"];DD**)
TEBD[[centerSite2,1]]=SparseArray[Chop[(\[Lambda]2Di.#&)/@(Partition[X,{\[Chi],\[Chi]}][[All,1]])]];
(* The right tensor needs the inverse right lambda to fit *)
(**O Print["Check 3-5"];DD**)
(**O Print["Check 3-5 B: site "<>ToString[rightSite]<>", "<>ToString[Dimensions[#]&/@{YT,(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]]),Transpose[{PadRight[Conjugate[#],\[Chi]]}]&/@YT}]];DD**)
\[Lambda]RDi=SparseArray[PseudoInverse[Normal[\[Lambda]RD],Tolerance->TOLERANCE]];
TEBD[[rightSite,1]]=SparseArray[Chop[(#.\[Lambda]RDi&)/@(Partition[ConjugateTranspose[YT],{\[Chi],\[Chi]}][[1,All]])]];
3-norm1-norm2-norm3]


SetAttributes[TEBD1QGate,HoldFirst]
TEBD1QGate[TEBD_,site_,H_]:=Module[{error,numSites,spin},
numSites=Length[TEBD];
If[(1<=site<=numSites)!=True,Print["TEBD1QGate called with wrong site."];Abort[]];
spin=Length[TEBD[[site,1]]];
If[Dimensions[H]!={spin,spin},Print["TEBD3QGate called with wrong operator dimensions."];Abort[]];
TEBD[[site,1]]=SparseArray[H.TEBD[[site,1]]];
0
]


SetAttributes[TEBDEvol3BodyPBC,HoldAll]
TEBDEvol3BodyPBC[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin,U,Uedge},
numSites=Length[TEBD];
If[Dimensions[H]!={spin^3,spin^3},Print["TEBDEvol3BodyPBC called with wrong operator dimensions."];Abort[]];
U=MatrixExp[-I H \[Tau]/2];
Uedge=MatrixExp[-I H \[Tau]];
error=0;
(* First do the bulk evolutions *)
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,U];,{n,2,numSites-1,3}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,U];,{n,3,numSites-1,3}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,U];,{n,4,numSites-1,3}];
(* Now the edge interactions are made through swaps. Sites 1 and 2 are taken to the right end *)
Do[
error+=TEBD2QGate[TEBD,n+1,Swap];
error+=TEBD2QGate[TEBD,n,Swap];,
{n,1,numSites-2}];
error+=TEBD3QGate[TEBD,numSites-2,U];
(* Up to here is one evolution operator for half the time. Now do the central *)
error+=TEBD3QGate[TEBD,numSites-1,Uedge];
(* From here do the same transposed *)
error+=TEBD3QGate[TEBD,numSites-2,U];
(* Swap back into position *)
Do[
error+=TEBD2QGate[TEBD,n,Swap];
error+=TEBD2QGate[TEBD,n+1,Swap];,
{n,numSites-2,1,-1}];
(* Do the bulk again *) 
Do[error+=TEBD3QGate[TEBD,n,U];,{n,4,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,U];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,U];,{n,2,numSites-1,3}];
error
]


SetAttributes[TEBDEvol3BodyPBCDisorder,HoldAll]
TEBDEvol3BodyPBCDisorder[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin},
numSites=Length[TEBD];
If[Dimensions[H]!={numSites,spin^3,spin^3},Print["TEBDEvol3BodyPBC called with wrong operator dimensions."];Abort[]];
error=0;
(* First do the bulk evolutions *)
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,2,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,4,numSites-1,3}];
(* Now the edge interactions are made through swaps. Sites 1 and 2 are taken to the right end *)
Do[
error+=TEBD2QGate[TEBD,n+1,Swap];
error+=TEBD2QGate[TEBD,n,Swap];,
{n,1,numSites-2}];
error+=TEBD3QGate[TEBD,numSites-2,MatrixExp[-I H [[numSites-2]]\[Tau]/2]];
(* Up to here is one evolution operator for half the time. Now do the central *)
error+=TEBD3QGate[TEBD,numSites-1,MatrixExp[-I H [[numSites-1]]\[Tau]]];
(* From here do the same transposed *)
error+=TEBD3QGate[TEBD,numSites-2,MatrixExp[-I H [[numSites-2]]\[Tau]/2]];
(* Swap back into position *)
Do[
error+=TEBD2QGate[TEBD,n,Swap];
error+=TEBD2QGate[TEBD,n+1,Swap];,
{n,numSites-2,1,-1}];
(* Do the bulk again *) 
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,4,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,2,numSites-1,3}];
error
]


SetAttributes[TEBDEvol3BodyOBC,HoldAll]
TEBDEvol3BodyOBC[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin,U,Uedge},
numSites=Length[TEBD];
If[Dimensions[H]!={spin^3,spin^3},Print["TEBDEvol3BodyPBC called with wrong operator dimensions."];Abort[]];
U=MatrixExp[-I H \[Tau]/2];
Uedge=MatrixExp[-I H \[Tau]];
error=0;
(* First do the bulk evolutions *)
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,U];,{n,2,numSites-1,3}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,U];,{n,3,numSites-1,3}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD3QGate[TEBD,n,Uedge];,{n,4,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,U];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,U];,{n,2,numSites-1,3}];
error
]


SetAttributes[TEBDEvol3BodyOBCDisorder,HoldAll]
TEBDEvol3BodyOBCDisorder[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin},
numSites=Length[TEBD];
If[Dimensions[H]!={numSites,spin^3,spin^3},Print["TEBDEvol3BodyPBC called with wrong operator dimensions."];Abort[]];
error=0;
(* First do the bulk evolutions *)
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,2,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]]];,{n,4,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,3,numSites-1,3}];
Do[error+=TEBD3QGate[TEBD,n,MatrixExp[-I H [[n]]\[Tau]/2]];,{n,2,numSites-1,3}];
error
]


SetAttributes[TEBDEvol2BodyPBC,HoldAll]
TEBDEvol2BodyPBC[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin,U,Uedge},
numSites=Length[TEBD];
If[Dimensions[H]!={spin^2,spin^2},Print["TEBDEvol2BodyPBC called with wrong operator dimensions."];Abort[]];
U=MatrixExp[-I H \[Tau]/2];
Uedge=MatrixExp[-I H \[Tau]];
error=0;
(* First do the bulk evolutions *)
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD2QGate[TEBD,n,U];,{n,1,numSites-1,2}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD2QGate[TEBD,n,U];,{n,2,numSites-1,2}];
(* Now the edge interactions are made through swaps.  *)
Do[
error+=TEBD2QGate[TEBD,n,Swap];,
{n,1,numSites-1}];
(* Up to here is one evolution operator for half the time. Now do the central *)
error+=TEBD2QGate[TEBD,numSites-1,Uedge];
(* From here do the same transposed *)
(* Swap back into position *)
Do[
error+=TEBD2QGate[TEBD,n,Swap];
,{n,numSites-1,1,-1}];
(* Do the bulk again *) 
Do[error+=TEBD2QGate[TEBD,n,U];,{n,2,numSites-1,2}];
Do[error+=TEBD2QGate[TEBD,n,U];,{n,1,numSites-1,2}];
error
]


SetAttributes[TEBDEvol2BodyOBC,HoldAll]
TEBDEvol2BodyOBC[TEBD_,H_,\[Tau]_]:=Module[{error,numSites,spin,U,Uedge},
numSites=Length[TEBD];
If[Dimensions[H]!={spin^2,spin^2},Print["TEBDEvol3BodyPBC called with wrong operator dimensions."];Abort[]];
U=MatrixExp[-I H \[Tau]/2];
Uedge=MatrixExp[-I H \[Tau]];
error=0;
(* First do the bulk evolutions *)
Do[(**OPrint["n"<>ToString[n]];DD**)
error+=TEBD2QGate[TEBD,n,U],{n,1,numSites-1,2}];
Do[(**OPrint["n"<>ToString[n]];DD**)error+=TEBD2QGate[TEBD,n,Uedge],{n,2,numSites-1,2}];
Do[error+=TEBD2QGate[TEBD,n,U],{n,1,numSites-1,2}];
error
]


(* ::Subsubsection::Closed:: *)
(*Wrappers*)


SetAttributes[TEBDEvolve,HoldFirst];
Options[TEBDEvolve]={Verbose->False};
TEBDEvolve[TEBD_,OperatorList_,steps_:1,OptionsPattern[]]:=Module[{error,verbose=OptionValue[Verbose],message},
error=0;
If[verbose,message=PrintTemporary["Starting up."]];
Do[
If[verbose,NotebookDelete[message];message=PrintTemporary["Evolving step "<>ToString[m]<>" of "<>ToString[steps]]];
Do[
error+=TEBDGate[TEBD,OperatorList[[n]]];
,{n,1,Length[OperatorList],1}];
,{m,1,steps}];
If[verbose,NotebookDelete[message]];
error
]



(* ::Subsection:: *)
(*Expectation values*)


SetAttributes[TEBDExpectation,HoldFirst];
TEBDExpectation[TEBD_,Operator_,site_]:=Module[{Mdown,Mup,value,numSites,\[Lambda]L,\[Lambda]R,spin},
value=0;
numSites=Length[TEBD];
If[(1<=site<=numSites)!=True,Print["TEBDExpectation1O called with wrong site parameter."];Abort[]];
spin=Length[TEBD[[site,1]]];
If[Dimensions[Operator]!={spin,spin},Print["TEBDExpectation1O called with wrong operator dimensions."];Abort[]];
If[site==1,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[numSites,2]]]];,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[site-1,2]]]];
];
\[Lambda]R=SparseArray[DiagonalMatrix[TEBD[[site,2]]]];
Mup=Flatten[\[Lambda]L.Transpose[(Operator.TEBD[[site,1]]),{2,1,3}].\[Lambda]R];
Mdown=Flatten[\[Lambda]L.Transpose[Conjugate[TEBD[[site,1]]],{2,1,3}].\[Lambda]R];
Chop[Mdown.Mup]
]


SetAttributes[TEBDCorrelation,HoldAll];TEBDCorrelation[TEBD_,OperatorL_,site1_,OperatorR_,site2_]:=Module[{numSites,M,temp,spin,\[Lambda]L,siteL,siteR},
numSites=Length[TEBD];
{siteL,siteR}=Sort[{site1,site2}];
If[(1<=siteL<siteR<=numSites)!=True,Print["TEBDExpectation2OList called with wrong site parameter."];Abort[]];
spin=Length[TEBD[[siteL,1]]];
If[Dimensions[OperatorL]!={spin,spin}||Dimensions[OperatorR]!={spin,spin},Print["TEBDExpectation1O called with wrong operator dimensions."];Abort[]];
If[siteL==1,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[numSites,2]]]];,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[siteL-1,2]]]];
];
(* This will do a loop over the different tensors in the TEBD *)

(* First prepare the initial matrix with the action of the Operator1 on the starting tensor *)
M=SparseArray[DiagonalMatrix[TEBD[[siteL,2]]]].Sum[ConjugateTranspose[TEBD[[siteL,1,i]]].(\[Lambda]L^2).((OperatorL.TEBD[[siteL,1]])[[i]]),{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteL,2]]]];
(* The expectation value is added a table as a function of distance *)
Table[{dist,
(* First we calculate the expectation value using Operator2 *)
temp=Tr[SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]].Sum[ConjugateTranspose[TEBD[[siteL+dist,1,i]]].M.((OperatorR.TEBD[[siteL+dist,1]])[[i]]),{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]]];
(* Now almost the same computation but for the following iteration, without the Operator2. The last computation is useless... *)
M=SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]].Sum[ConjugateTranspose[TEBD[[siteL+dist,1,i]]].M.TEBD[[siteL+dist,1,i]],{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]];
Chop[temp]
},{dist,1,siteR-siteL}]
]


SetAttributes[TEBDExpectation,HoldFirst];
TEBDExpectation[TEBD_,OperatorL_,siteL_,OperatorR_,siteR_]:=Last[TEBDCorrelation[TEBD,OperatorL,siteL,OperatorR,siteR]][[2]]


SetAttributes[TEBDExpectation3O,HoldAll];
TEBDExpectation[TEBD_,OperatorL_,site1_,OperatorC_,site2_,OperatorR_,site3_]:=Module[{numSites,M,temp,spin,\[Lambda]L,siteL,siteC,siteR},
numSites=Length[TEBD];
{siteL,siteC,siteR}=Sort[{site1,site2,site3}];
If[(1<=siteL<siteC<siteR<=numSites)!=True,Print["TEBDExpectation3O called with wrong site parameter."];Abort[]];
spin=Length[TEBD[[siteL,1]]];
If[Dimensions[OperatorL]!={spin,spin}||Dimensions[OperatorC]!={spin,spin}||Dimensions[OperatorR]!={spin,spin},Print["TEBDExpectation1O called with wrong operator dimensions."];Abort[]];
If[siteL==1,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[numSites,2]]]];,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[siteL-1,2]]]];
];
(* First prepare the initial matrix with the action of the OperatorL on the starting tensor *)
M=SparseArray[DiagonalMatrix[TEBD[[siteL,2]]]].Sum[ConjugateTranspose[TEBD[[siteL,1,i]]].(\[Lambda]L^2).((OperatorL.TEBD[[siteL,1]])[[i]]),{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteL,2]]]];
(* Now eat tensors away until the center site siteC *)
Do[M=SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]].Sum[ConjugateTranspose[TEBD[[siteL+dist,1,i]]].M.TEBD[[siteL+dist,1,i]],{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteL+dist,2]]]];,{dist,1,siteC-1-siteL}];
(* Now eat the center tensor *)
M=SparseArray[DiagonalMatrix[TEBD[[siteC,2]]]].Sum[ConjugateTranspose[TEBD[[siteC,1,i]]].M.((OperatorR.TEBD[[siteC,1]])[[i]]),{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteC,2]]]];
(* And keep on going until we reach the last tensor *)
Do[M=SparseArray[DiagonalMatrix[TEBD[[siteC+dist,2]]]].Sum[ConjugateTranspose[TEBD[[siteC+dist,1,i]]].M.TEBD[[siteC+dist,1,i]],{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteC+dist,2]]]];,{dist,1,siteR-1-siteC}];
(* Return expectation value *)
Chop[Tr[SparseArray[DiagonalMatrix[TEBD[[siteR,2]]]].Sum[ConjugateTranspose[TEBD[[siteR,1,i]]].M.((OperatorR.TEBD[[siteR,1]])[[i]]),{i,1,spin}].SparseArray[DiagonalMatrix[TEBD[[siteR,2]]]]]]
]

TEBDExpectation3O[TEBD_,Operator_,site_]:=Module[{numSites,M,temp,spin,\[Lambda]L},
numSites=Length[TEBD];
If[(1<site<numSites)!=True,Print["TEBDExpectation3O (B) called with wrong site parameter."];Abort[]];
spin=Length[TEBD[[site,1]]];
If[Dimensions[Operator]!={spin^3,spin^3},Print["TEBDExpectation1O(B) called with wrong operator dimensions."];Abort[]];
If[site==2,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[numSites,2]]]];,
\[Lambda]L=SparseArray[DiagonalMatrix[TEBD[[site-2,2]]]];
];
(* First prepare the initial matrix with the action of the OperatorL on the starting tensor *)
M=SparseArray[DiagonalMatrix[TEBD[[site+1,2]]]].Sum[ConjugateTranspose[Flatten[TEBD[[site-1,1]].SparseArray[DiagonalMatrix[TEBD[[site-1,2]]]].Transpose[TEBD[[site,1]],{2,1,3}].SparseArray[DiagonalMatrix[TEBD[[site,2]]]].Transpose[TEBD[[site+1,1]],{2,1,3}],{{1,3,4},{2},{5}}][[i]]].(\[Lambda]L^2).(Operator.Flatten[TEBD[[site-1,1]].SparseArray[DiagonalMatrix[TEBD[[site-1,2]]]].Transpose[TEBD[[site,1]],{2,1,3}].SparseArray[DiagonalMatrix[TEBD[[site,2]]]].Transpose[TEBD[[site+1,1]],{2,1,3}],{{1,3,4},{2},{5}}]
)[[i]],{i,1,spin^3}].SparseArray[DiagonalMatrix[TEBD[[site+1,2]]]];
Chop[Tr[M]]
]


TEBDEntropy[\[Lambda]_]:=-Sum[If[\[Lambda][[n]]>0,\[Lambda][[n]]^2Log[2,\[Lambda][[n]]^2],0],{n,1,Length[\[Lambda]]}];
TEBDEntropyList[TEBD_]:=Table[{n,TEBDEntropy[TEBD[[n,2]]]},{n,1,Length[TEBD]-1}]


(* ::Subsection::Closed:: *)
(*Hamiltonian and evolution operators*)


SetAttributes[TEBDInitHamiltonian,HoldFirst];
TEBDInitHamiltonian[Ham_]:=({})


SetAttributes[TEBDHamiltonianAdd,HoldFirst];
TEBDHamiltonianAdd[Ham_,operator_,site1_,site2_]:=(Ham=(Ham~Join~{{operator,site1,site2}}))


TEBDHamiltonianAdd[Ham_,operator1_,site1_,operator2_,site2_]:=(Ham=(Ham~Join~{{KroneckerProduct[operator1,operator2],site1,site2}}))


TEBDHamiltonianAdd[Ham_,operator_,site_]:=(Ham=(Ham~Join~{{operator,site,site}}))


TEBDEvolutionFromHamiltonian[Ham_,time_]:=Table[{MatrixExp[-I time Ham[[n,1]]],Ham[[n,2]],Ham[[n,3]]},{n,1,Length[Ham]}]


SetAttributes[TEBDInitEvolution,HoldFirst];
TEBDInitEvolution[]:=({})


SetAttributes[TEBDEvolutionAddGate,HoldFirst];
TEBDEvolutionAddGate[U_,operator_,site1_,site2_]:=(U=(U~Join~{{operator,site1,site2}}))


TEBDEvolutionAddGate[U_,operator_,site_]:=(U=(U~Join~{{operator,site,site}}))



