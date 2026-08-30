(* ::Package:: *)

wkdir = If[
	$InputFileName != "", 
	DirectoryName[AbsoluteFileName[$InputFileName]], 
	NotebookDirectory[]
];


(* ::Section::Closed:: *)
(*Imports*)


(* ::Subsection::Closed:: *)
(*Toolkit*)


(* ::Text:: *)
(*Mino frequencies are related to the Kerr orbital frequencies via \[Omega]r\[Theta]\[Phi] = \[CapitalGamma]r\[Theta]\[Phi]/\[CapitalGamma]t.*)


Get[wkdir<>"KerrGeoOrbit.m"];
Get[wkdir<>"KerrGeodesics-0.9.0/Kernel/OrbitalFrequencies.m"];
Get[wkdir<>"KerrGeodesics-0.9.0/Kernel/SpecialOrbits.m"];


(* ::Subsection::Closed:: *)
(*dELQdtRates*)


(* ::Subsubsection::Closed:: *)
(*Complete expressions*)


{dEdtRate,dLdtRate,dQdtRate} = Get[wkdir<>"Acceleration/dELQdtRate.m"];


Clear[dEdtRate$Fn,dLdtRate$Fn,dQdtRate$Fn]
dEdtRate$Fn[m1_,m2_,b_,\[Epsilon]\[CapitalEpsilon]_][v_,w_][u1_,u2_,u3_,u4_,u_][r_,z_,x_,y_] = dEdtRate;
dLdtRate$Fn[m1_,m2_,b_,\[Epsilon]\[CapitalEpsilon]_][v_,w_][u1_,u2_,u3_,u4_,u_][r_,z_,x_,y_] = dLdtRate;
dQdtRate$Fn[m1_,m2_,b_,\[Epsilon]\[CapitalEpsilon]_][v_,w_][u1_,u2_,u3_,u4_,u_][r_,z_,x_,y_] = dQdtRate;


(* ::Subsubsection::Closed:: *)
(*Compile individual z^m parts*)


(* ::Text:: *)
(*Note: on HPC, need WVM instead of C because is external tool basically and when running things on different node, it can get confused (ask AI).*)
(*The WVM is the default "internal" compiler. When you compile to WVM, Mathematica turns the code into a specialized bytecode that only the Mathematica kernel understands.*)
(*WVM is only slightly slower than C.*)
(*Alternatively, use DistributieDefinitions.*)


Clear[dEdtRate$zmTerm,dLdtRate$zmTerm,dQdtRate$zmTerm]
Do[
dEdtRate$zmTerm[m] = z^m * Simplify@Coefficient[dEdtRate$Fn[m1,m2,b,\[Epsilon]\[CapitalEpsilon]][v,w][u1,u2,u3,u4,u][r,z,x,y], z, m];
dLdtRate$zmTerm[m] = z^m * Simplify@Coefficient[dLdtRate$Fn[m1,m2,b,\[Epsilon]\[CapitalEpsilon]][v,w][u1,u2,u3,u4,u][r,z,x,y], z, m];
dQdtRate$zmTerm[m] = z^m * Simplify@Coefficient[dQdtRate$Fn[m1,m2,b,\[Epsilon]\[CapitalEpsilon]][v,w][u1,u2,u3,u4,u][r,z,x,y], z, m]
,{m,-2,2}]


Do[

dEdtRate$zCoeff$Compiled[m] = With[{rawexpr = dEdtRate$zmTerm[m]},
	Compile[{{m1, _Real}, {m2, _Real}, {b, _Real}, {\[Epsilon]\[CapitalEpsilon], _Real}, {v, _Complex}, {w, _Complex}, {u1, _Real}, {u2, _Real}, {u3, _Real}, {u4, _Real}, {u, _Real}, {r, _Real}, {z, _Complex}, {x, _Complex}, {y, _Complex}},
    rawexpr,
    CompilationTarget -> "WVM",
    RuntimeAttributes -> {Listable},
    Parallelization -> False,
    RuntimeOptions -> "Speed"
   ]
];

dLdtRate$zCoeff$Compiled[m] = With[{rawexpr = dLdtRate$zmTerm[m]},
	Compile[{{m1, _Real}, {m2, _Real}, {b, _Real}, {\[Epsilon]\[CapitalEpsilon], _Real}, {v, _Complex}, {w, _Complex}, {u1, _Real}, {u2, _Real}, {u3, _Real}, {u4, _Real}, {u, _Real}, {r, _Real}, {z, _Complex}, {x, _Complex}, {y, _Complex}},
    rawexpr,
    CompilationTarget -> "WVM",
    RuntimeAttributes -> {Listable},
    Parallelization -> False,
    RuntimeOptions -> "Speed"
   ]
];

dQdtRate$zCoeff$Compiled[m] = With[{rawexpr = dQdtRate$zmTerm[m]},
	Compile[{{m1, _Real}, {m2, _Real}, {b, _Real}, {\[Epsilon]\[CapitalEpsilon], _Real}, {v, _Complex}, {w, _Complex}, {u1, _Real}, {u2, _Real}, {u3, _Real}, {u4, _Real}, {u, _Real}, {r, _Real}, {z, _Complex}, {x, _Complex}, {y, _Complex}},
    rawexpr,
    CompilationTarget -> "WVM",
    RuntimeAttributes -> {Listable},
    Parallelization -> False,
    RuntimeOptions -> "Speed"
   ]
];

,{m,-2,2}]


(* ::Section::Closed:: *)
(*Parallel computing*)


(*ParallelEvaluate[Get[wkdir<>"KerrGeoOrbit.m"]];
ParallelEvaluate[Get[wkdir<>"KerrGeodesics-0.9.0/Kernel/OrbitalFrequencies.m"]];
ParallelEvaluate[Get[wkdir<>"KerrGeodesics-0.9.0/Kernel/SpecialOrbits.m"]];*)


(*ParallelNeeds["KerrGeodesics`KerrGeoOrbit`"];
ParallelNeeds["KerrGeodesics`OrbitalFrequencies`"];
ParallelNeeds["KerrGeodesics`SpecialOrbits`"];
Print["Number of kernels mathematica uses: " <> ToString[$KernelCount]];*)


(* ::Section::Closed:: *)
(*Parameters*)


(* If you want to consider the FULL perturbation (electric + magnetic), choose \[Epsilon]\[CapitalEpsilon]val = 1 *)
(* If you want to consider the MAGNETIC perturbation part ONLY (no electric), choose \[Epsilon]\[CapitalEpsilon]val = 0 *)
\[Epsilon]\[CapitalEpsilon]val = 1;

(* Numerical values *)
au = 1.49597870700 10^11;(*m*)
MSun = 1.98855 10^30; (*kg*)
G = 6.67430 10^-11;(*m^3 kg^-1 s^-2*)
c = 2.99792458 10^8(*m/s*);

(*
m1 = 4 10^6 MSun;(*kg*)(* Central black hole mass *)
m2 = 30 MSun;(*kg*)(* Small EMRI mass *)
m3 = 30 MSun;(*kg*)(* Tidal perturber mass *)
b = 10 au;(*meters*)(* Distance between m1 and m3 *)
*)

m1 = 1;(*kg*)(* Central black hole mass *)
m2 = 1;(*kg*)(* Small EMRI mass *)
m3 = 30 MSun;(*kg*)(* Tidal perturber mass *)
b = 1;(*meters*)(* Distance between m1 and m3 *)

Prec = 50;

fact1 = m1 *G c^-3;(* seconds *)
fact2 = m1*G c^-3/(3.154 10^7(*sec/year*));(* years *)
\[Eta] = m2/m1;(* EMRI mass ratio *)


(* ::Section::Closed:: *)
(*Resonance condition*)


OmegaDynamic[m1_, m3_, b_] := (G*m1)/c^3 Sqrt[(G*(m1+m3))/b^3];


ResConStatic[n_, k_, m_][a_, p_, e_, x_] := Module[{freq, Gammar, Gamma\[Theta], Gamma\[Phi], Gammat},

	freq = KerrGeoFrequencies[a, p, e, x, Time -> "Mino"];
	{Gammar, Gamma\[Theta], Gamma\[Phi], Gammat} = Values@freq;

	n Gammar + k Gamma\[Theta] + m Gamma\[Phi]
  ]


ResConDynamic[n_, k_, m_, s_][a_, p_, e_, x_][m1_, m3_, b_] := Module[{freq, Gammar, Gamma\[Theta], Gamma\[Phi], Gammat},

	freq = KerrGeoFrequencies[a, p, e, x, Time -> "Mino"];
	{Gammar, Gamma\[Theta], Gamma\[Phi], Gammat} = Values@freq;

	n Gammar + k Gamma\[Theta] + m Gamma\[Phi] + s OmegaDynamic[m1, m3, b]
  ]


\[CapitalDelta]sep = 10^(-9);
resConBuffer = 10^(-9);
findp[n_, k_, m_][a_, e_, x_] := 
 Module[{pseparatrix, ResCon, pguess, ResConValue, p},

	pseparatrix = Quiet@KerrGeoSeparatrix[a, e, x];
	ResCon[p_] = ResConStatic[n, k, m][a, p, e, x]; 
   
	(* pguess = p /. FindRoot[Abs@ResCon[p] == 0, {p, pseparatrix+2, pseparatrix+\[CapitalDelta]sep, 100+pseparatrix},WorkingPrecision->100]; *)
	(* pguess = p /. FindRoot[Abs@ResCon[p] == 0, {p, pseparatrix+\[CapitalDelta]sep, 5+pseparatrix},WorkingPrecision->100]; *)
	pguess = p /. Last@Quiet@NMinimize[{Abs@ResCon[p], pseparatrix+\[CapitalDelta]sep < p < 30+pseparatrix}, p,  MaxIterations -> 200, WorkingPrecision -> Prec-3];

     (* check the static resonance condition for the calculated p value *)
	If[Abs@ResCon[pguess] > resConBuffer, "No resonance", pguess]
  ]


(* ::Section::Closed:: *)
(*4 velocity*)


(* ::Text:: *)
(*Fixed the dimensions by introducing factors of m1.*)
(*Note: components do not share the same dimensions!*)


Clear[u1$Fn,u2$Fn,u3$Fn,u4$Fn]
u1$Fn[E_] := -E;
u2$Fn[a_,E_,L_,Q_][r_] := Sqrt@Abs[(E (r^2 + a^2) - a L)^2 - (m1^2 r^2 + Q + (L - a E)^2) (r^2 - 2 m1 r + a^2)]/(r^2 - 2 m1 r + a^2);
u3$Fn[a_,E_,L_,Q_][\[Theta]_] := Sqrt@Abs[Q - L^2 Cot[\[Theta]]^2 - a^2 (m1^2 - E^2) Cos[\[Theta]]^2];
u4$Fn[L_] := L;


(* ::Section::Closed:: *)
(*MinoPeriods and OrbitralFunc*)


GetMinoPeriods[n_, k_, m_][a_, p_, e_, x_] := 
 Module[{freq, Gammat, Gammar, Gamma\[Theta], Gamma\[Phi], \[CapitalLambda]tt, \[CapitalLambda]rr, \[CapitalLambda]\[Theta]\[Theta], \[CapitalLambda]\[Phi]\[Phi]},
	freq = KerrGeoFrequencies[a, p, e, x, Time -> "Mino"];
	{Gammar, Gamma\[Theta], Gamma\[Phi], Gammat} = Values@freq;
	(* equation 2.2.8 in Barts thesis *)
	\[CapitalLambda]tt = 2*Pi / Gammat; (* Temporal Mino Period *)
	\[CapitalLambda]rr = 2*Pi / Gammar; (*Radial Mino Period *)
	\[CapitalLambda]\[Theta]\[Theta] = 2*Pi / Gamma\[Theta]; (* Polar Mino Period *)
	\[CapitalLambda]\[Phi]\[Phi] = 2*Pi / Gamma\[Phi]; (* Azimuthal Mino Period *)
	{\[CapitalLambda]tt, \[CapitalLambda]rr, \[CapitalLambda]\[Theta]\[Theta], \[CapitalLambda]\[Phi]\[Phi]}
  ]


GetOrbitalFunc[n_, k_, m_][a_, p_, e_, x_] := 
 Module[{Case6, ft, fr, f\[Theta], f\[Phi], f\[Phi]1, f\[Phi]2, f\[Phi]3},
	(* Should Print "Reading the changed KerrGeoOrbit package_v5" *)
	Case6 = KerrGeoOrbit[a, p, e, x, Method -> "Analytic"];
	ft[\[Lambda]_] := Re[Case6[\[Lambda]]][[1]]; (* t(\[Lambda]t)  *)
	fr[\[Lambda]_] := Re[Case6[\[Lambda]]][[2]]; (* r(\[Lambda]r) *)
	f\[Theta][\[Lambda]_] := Re[Case6[\[Lambda]]][[3]]; (* \[Theta](\[Lambda]\[Theta]) *)
	f\[Phi][\[Lambda]_] := Re[Case6[\[Lambda]]][[4]]; (* \[Phi](\[Lambda]r,\[Lambda]\[Theta],\[Lambda]\[Phi]); can be split into the thress pieces below; see eq 2.4.8 in Barts thesis *)
	f\[Phi]1[\[Lambda]_] := Re[Case6[\[Lambda]]][[5]]; (* \[Phi](\[Lambda]r) *)
	f\[Phi]2[\[Lambda]_] := Re[Case6[\[Lambda]]][[6]]; (* \[Phi](\[Lambda]\[Theta]) *)
	f\[Phi]3[\[Lambda]_] := Re[Case6[\[Lambda]]][[7]];  (* \[Phi](\[Lambda]\[Phi]) *)
	{ft,fr,f\[Theta],f\[Phi],f\[Phi]1,f\[Phi]2,f\[Phi]3}
  ]


(* ::Text:: *)
(*The below is used when multiplying u2 and u3. We want the sign of u2 (u3) to be positive / negative when r(\[Lambda]) (\[Theta](\[Lambda])) is increasing/decreasing.*)
(*For r(\[Lambda]) for example, because it is periodic with period \[CapitalLambda]r, and symmetric around its half period, that means that \[Lambda] = \[CapitalLambda]r / 2 is where the extrema is reached. So, for 0 < \[Lambda] < \[CapitalLambda]r / 2, we want it positive, and \[CapitalLambda]r / 2 < \[Lambda] < 1 negative.*)


Clear[MinoFreqSign]
MinoFreqSign[\[Lambda]_, \[CapitalLambda]_] := (-1)^Floor[2 * \[Lambda] / \[CapitalLambda]];


(* ::Section::Closed:: *)
(*Jumps*)


(* ::Text:: *)
(*Notes for extreme parameter cases:*)
(*x=1 needs to be exact due to a bug in the KerrGeoSeparatrix package.*)
(*a=0 needs to be exact.*)
(*e=0 does not work if x = 1 as well, due to some bug.*)
(**)
(*a=1 does not work, probably too extreme anyway.*)
(*e=1 does not work because parabolic trajectory.*)
(*x=0 does not work, I think due to coordinate singularity.*)
(**)
(*I recommend 300 <= Nz <= 400, particularly Nz = 300.*)
(*This gets the relative error in all the jumps to be ~10^(-6) or less in most cases. The drop in accuracy happens for a,e ~ 1 and x ~ 0.*)
(*Close to, and beyond I imagine, Nz = 400 the Interpolation starts acting up.*)


(* Reinforcing static case here with u = 0 and y = 1 (\[Omega]=0) *)
ComputeJumpsELQ[n_, k_, m_][a_, e_, x_, incy_, incz_][Nz_] := Module[
{v = Exp[I incz], w = Exp[I incy], u = 0, y = 1,
p, E, L, Q, \[CapitalLambda]t, \[CapitalLambda]r, \[CapitalLambda]\[Theta], \[CapitalLambda]\[Phi], \[CapitalLambda]r\[Theta]Ratio, ft, fr, f\[Theta], f\[Phi], f\[Phi]1, f\[Phi]2, f\[Phi]3, \[Lambda]rGrid, \[Lambda]\[Theta]Grid, \[Lambda]r\[Theta]Grid, fr\[Lambda]rVec, f\[Theta]\[Lambda]\[Theta]Vec, u2Vec, u3Vec, xVec, zMat, z\[Theta], RateFnsList, result, i, dELQdtRateMat, dELQdtRateGrid, dELQdtRateInt},

	p = findp[n, k, m][a, e, x];
	If[NumberQ@p,

		(* Get the conserved quantities E, Lz, Q *)
		E = KerrGeoEnergy[a, p, e, x];
		L = KerrGeoAngularMomentum[a, p, e, x];
		Q = KerrGeoCarterConstant[a, p, e, x];

		{\[CapitalLambda]t, \[CapitalLambda]r, \[CapitalLambda]\[Theta], \[CapitalLambda]\[Phi]} = GetMinoPeriods[n, k, m][a, p, e, x];
		\[CapitalLambda]r\[Theta]Ratio = Round[\[CapitalLambda]r / \[CapitalLambda]\[Theta]];

		{ft, fr, f\[Theta], f\[Phi], f\[Phi]1, f\[Phi]2, f\[Phi]3} =  GetOrbitalFunc[n, k, m][a, p, e, x];

		\[Lambda]rGrid = \[CapitalLambda]r/(\[CapitalLambda]r\[Theta]Ratio * Nz) * Range[0, \[CapitalLambda]r\[Theta]Ratio * Nz];
		\[Lambda]\[Theta]Grid = \[CapitalLambda]\[Theta]/Nz * Range[0, Nz];
		\[Lambda]r\[Theta]Grid = Outer[List, \[Lambda]rGrid, \[Lambda]\[Theta]Grid];

		fr\[Lambda]rVec = fr[\[Lambda]rGrid];
		f\[Theta]\[Lambda]\[Theta]Vec = f\[Theta][\[Lambda]\[Theta]Grid];
		u2Vec = MinoFreqSign[\[Lambda]rGrid, \[CapitalLambda]r] * u2$Fn[a, E, L, Q][fr\[Lambda]rVec];
		u3Vec = MinoFreqSign[\[Lambda]\[Theta]Grid, \[CapitalLambda]\[Theta]] * u3$Fn[a, E, L, Q][f\[Theta]\[Lambda]\[Theta]Vec];

		xVec = Exp[I * f\[Theta]\[Lambda]\[Theta]Vec];
		zMat = KroneckerProduct[Exp[I * f\[Phi]1[\[Lambda]rGrid]], Exp[I * f\[Phi]2[\[Lambda]\[Theta]Grid]]];
 
		result = ConstantArray[Null,3];
		RateFnsList = {dEdtRate$zCoeff$Compiled, dLdtRate$zCoeff$Compiled, dQdtRate$zCoeff$Compiled};

		Do[
		(* Table[dEdtRate$zCoeff$Compiled[m][m1, m2, b, \[Epsilon]\[CapitalEpsilon]val, v, w, u1$Fn[E], u2Vec[[ir+1]], u3Vec[[i\[Theta]+1]], u4$Fn[L], u, fr\[Lambda]rVec[[ir+1]], zVec[[ir+1,i\[Theta]+1]], xVec[[i\[Theta]+1]], y], {ir, 0, (\[CapitalLambda]r\[Theta]Ratio * Nz)}, {i\[Theta], 0, Nz}]; *)
			dELQdtRateMat = MapThread[
				RateFnsList[[i]][m][m1, m2, b, \[Epsilon]\[CapitalEpsilon]val, v, w, u1$Fn[E], #1, u3Vec, u4$Fn[L], u, #3, #2, xVec, y] &,
			{u2Vec, zMat, fr\[Lambda]rVec}];
			dELQdtRateGrid = Transpose[{Flatten[\[Lambda]r\[Theta]Grid, 1], Flatten[dELQdtRateMat]}];
			dELQdtRateInt = Interpolation[dELQdtRateGrid];

			result[[i]] = Quiet[\[CapitalLambda]t/(2*Pi * \[CapitalLambda]r * \[CapitalLambda]\[Theta]) NIntegrate[
				dELQdtRateInt[\[Lambda]r, \[Lambda]\[Theta]] Exp[-I * ((n * 2*Pi)/\[CapitalLambda]r \[Lambda]r + (k * 2*Pi)/\[CapitalLambda]\[Theta] \[Lambda]\[Theta])], 
				{\[Lambda]r, 0, \[CapitalLambda]r}, {\[Lambda]\[Theta], 0, \[CapitalLambda]\[Theta]}]]; (* 2.4.19 of Barts thesis *)

		,{i, 1, Length@RateFnsList}];
	];

    (* If False, p = "No resonance" *)
	If[NumberQ@p, {{n, k, m}, {a, p, e, x}, result}, {{n, k, m}, {a, p, e, x}, p}]
]


(* ::Section::Closed:: *)
(*Computation*)


Nz$a = 20;
Nz$e = 20;
Nz$x = 20;
Nz = 300;
incy = \[Pi]/4;
incz = 0;


klmListProg = {{3,0,-2}, {3,-2,0}, {3,-4,2}, {4,0,-2}, {4,-2,0}, {4,-4,2}, {1,2,-2}, {3,-1,-1}, {3,-3,1}, {4,-1,-1}, {4,-3,1}, {1,1,-1}};
klmListRet = {#1, #2, -#3} & @@@ klmListProg; (* Sends m -> -m *)


(* We remove the first and last elements *)
(* Removing a = 1 case *)
amin = 0;
amax = 1;
aGrid = Most@N[Table[amin + i / Nz$a * (amax - amin),{i, 0, Nz$a}], Prec];
aGrid[[1]] = amin;

(* Removing e = 1 case and replacing e = 0 by small number *)
ebuffer = N[10^(-15), Prec];
emin = 0;
emax = 1;
eGrid = Most@N[Table[emin + i / Nz$e * (emax - emin),{i, 0, Nz$e}], Prec];
eGrid[[1]] = ebuffer;

(* Removing x = 0 case *)
\[Theta]min = 0;
\[Theta]max = Pi;
\[Theta]Grid = N[Table[\[Theta]min + i/ Nz$x * (\[Theta]max - \[Theta]min),{i, 0, Nz$x}], Prec];
xGridProg = Cos[\[Theta]Grid[[;; Floor[(Nz$x+1)/2]]]]; (* 0 \[LessEqual] \[Theta] < \[Pi]/2 *)
xGridProg[[1]] = Cos[\[Theta]min];
xGridRet = Cos[\[Theta]Grid[[Ceiling[(Nz$x+1)/2]+1 ;;]]];  (* \[Pi]/2 < \[Theta] \[LessEqual] \[Pi] *)
xGridRet[[-1]] = Cos[\[Theta]max];


paraGridProg = Tuples[{aGrid, eGrid, xGridProg}];
paraGridRet = Tuples[{aGrid, eGrid, xGridRet}];


(*WriteLine[$Output, "Start computation"];*)


Do[
	{myn, myk, mym} = klmListProg[[i]];
	ELQJumpsGivennkm = ParallelTable[
		{mya, mye, myx} = paraGridProg[[j]];
		ComputeJumpsELQ[myn, myk, mym][mya, mye, myx, incy, incz][Nz]
	,{j,1,Length@paraGridProg}];
	ELQJumpsGivennkm = Prepend[ELQJumpsGivennkm,{{"n","k","m"},{"a","p","e","x"},{"dE","dL","dQ"}}];
	Export[wkdir<>"TidalResonance_Data_Jumps/Prograde/"<>"n"<>ToString[myn]<>"_k"<>ToString[myk]<>"_m"<>ToString[mym]<>".m",N@ELQJumpsGivennkm]
,{i,1,Length@klmListProg}]


(*
Do[
	{myn, myk, mym} = klmListRet[[i]];
	ELQJumpsGivennkm = ParallelTable[
		{mya, mye, myx} = paraGridRet[[j]];
		ComputeJumpsELQ[myn, myk, mym][mya, mye, myx, incy, incz][Nz]
	,{j,1,Length@paraGridRet}];
	Export[wkdir<>"TidalResonance_Data_Jumps/Retrograde/"<>"n"<>ToString[myn]<>"_k"<>ToString[myk]<>"_m"<>ToString[mym]<>".m",N@ELQJumpsGivennkm]
,{i,1,Length@klmListRet}]
*)


(* ::Section::Closed:: *)
(*Ways to improve things*)


(* ::Text:: *)
(*Currently, the largest bottle necks is NIntegrate.*)
(*Currently using Interpolation and then numerically integrate Interpolant. Because we can compute dEdt etc are any point, could simply use that directly in the integrand; respectively use the grid points directly for the numerical integration. I tried this below, but it is slower for some reason?*)
(**)
(*Also, in some cases, integration is actually the bottleneck. For example, use EchoTiming to check for (n,k,m) = (4,-4,2).*)


(* Reinforcing static case here with u = 0 and y = 1 (\[Omega]=0) *)
(*
ComputeJumpsELQDirect[n_, k_, m_][a_, e_, x_, incy_, incz_][Nz_] := Module[
{v = Exp[I incz], w = Exp[I incy], u = 0, y = 1,
p, E, L, Q, \[CapitalLambda]t, \[CapitalLambda]r, \[CapitalLambda]\[Theta], \[CapitalLambda]\[Phi], \[CapitalLambda]r\[Theta]Ratio, ft, fr, f\[Theta], f\[Phi], f\[Phi]1, f\[Phi]2, f\[Phi]3, \[Lambda]rGrid, \[Lambda]\[Theta]Grid, \[Lambda]r\[Theta]Grid, fr\[Lambda]rVec, f\[Theta]\[Lambda]\[Theta]Vec, u2Vec,u3Vec, xVec, zMat, z\[Theta], RateFnsList, result, i, dLdtRateGrid, LrateTab, LrateTabCoeff, LrateInt},

	p = findp[n, k, m][a, e, x];
	If[NumberQ@p,

		(* Get the conserved quantities E, Lz, Q *)
		E = KerrGeoEnergy[a, p, e, x];
		L = KerrGeoAngularMomentum[a, p, e, x];
		Q = KerrGeoCarterConstant[a, p, e, x];

		{\[CapitalLambda]t, \[CapitalLambda]r, \[CapitalLambda]\[Theta], \[CapitalLambda]\[Phi]} = GetMinoPeriods[n, k, m][a, p, e, x];
		\[CapitalLambda]r\[Theta]Ratio = Round[\[CapitalLambda]r / \[CapitalLambda]\[Theta]];

		{ft, fr, f\[Theta], f\[Phi], f\[Phi]1, f\[Phi]2, f\[Phi]3} =  GetOrbitalFunc[n, k, m][a, p, e, x];

		\[Lambda]rGrid = \[CapitalLambda]r/(\[CapitalLambda]r\[Theta]Ratio * Nz) * Range[0, \[CapitalLambda]r\[Theta]Ratio * Nz];
		\[Lambda]\[Theta]Grid = \[CapitalLambda]\[Theta]/Nz * Range[0, Nz];
		\[Lambda]r\[Theta]Grid = Outer[List, \[Lambda]rGrid, \[Lambda]\[Theta]Grid];

		fr\[Lambda]rVec = fr[\[Lambda]rGrid];
		f\[Theta]\[Lambda]\[Theta]Vec = f\[Theta][\[Lambda]\[Theta]Grid];
		u2Vec = MinoFreqSign[\[Lambda]rGrid, \[CapitalLambda]r] * u2$Fn[a, E, L, Q][fr\[Lambda]rVec];
		u3Vec = MinoFreqSign[\[Lambda]\[Theta]Grid, \[CapitalLambda]\[Theta]] * u3$Fn[a, E, L, Q][f\[Theta]\[Lambda]\[Theta]Vec];

		xVec = Exp[I * f\[Theta]\[Lambda]\[Theta]Vec];
		zMat = KroneckerProduct[Exp[I * f\[Phi]1[\[Lambda]rGrid]], Exp[I * f\[Phi]2[\[Lambda]\[Theta]Grid]]];
 
		result = ConstantArray[Null,3];
		RateFnsList = {dEdtRate$zCoeff$Compiled, dLdtRate$zCoeff$Compiled, dQdtRate$zCoeff$Compiled};

		Do[
			result[[i]] = Quiet[\[CapitalLambda]t/(2*Pi * \[CapitalLambda]r * \[CapitalLambda]\[Theta]) NIntegrate[
				RateFnsList[[i]][m][m1, m2, b, \[Epsilon]\[CapitalEpsilon]val, v, w, u1$Fn[E],MinoFreqSign[\[Lambda]r, \[CapitalLambda]r] * u2$Fn[a, E, L, Q][fr[\[Lambda]r]], MinoFreqSign[\[Lambda]\[Theta], \[CapitalLambda]\[Theta]] * u3$Fn[a, E, L, Q][f\[Theta][\[Lambda]\[Theta]]], u4$Fn[L], u, fr[\[Lambda]r],Exp[I * f\[Phi]1[\[Lambda]r]]*Exp[I * f\[Phi]2[\[Lambda]\[Theta]]], Exp[I * f\[Theta][\[Lambda]\[Theta]]], y] * Exp[-I * ((n * 2*Pi)/\[CapitalLambda]r \[Lambda]r + (k * 2*Pi)/\[CapitalLambda]\[Theta] \[Lambda]\[Theta])], 
				{\[Lambda]r, 0, \[CapitalLambda]r}, {\[Lambda]\[Theta], 0, \[CapitalLambda]\[Theta]},Method -> "ClenshawCurtisRule"]]; (* 2.4.19 of Barts thesis *)

		,{i, 1, Length@RateFnsList}];
	];

    (* If False, p = "No resonance" *)
	If[NumberQ@p, {{n, k, m}, {a, p, e, x}, result}, {{n, k, m}, {a, p, e, x}, p}]
]
*)
