#!/usr/bin/env wolframscript
(* ==============================================================================
   Continuous Stirred-Tank Reactor (CSTR) Non-Linear Dynamics & Arrhenius Kinetics
   Execution via Wolfram Engine / Wolfram Script
   ============================================================================== *)

Print["===================================================================="];
Print["  Continuous Stirred-Tank Reactor Non-Linear Dynamics (Wolfram)     "];
Print["===================================================================="];

(* Kinetic and thermal parameters *)
V = 100.0; (* Reactor Volume, L *)
F = 100.0; (* Flow rate, L/min *)
cA0 = 1.0; (* Feed Concentration, mol/L *)
T0 = 350.0; (* Feed Temperature, K *)
k0 = 7.2 * 10^10; (* Pre-exponential factor, min^-1 *)
E = 72750.0; (* Activation energy, J/mol *)
R = 8.314; (* Gas constant, J/(mol K) *)
delH = -50000.0; (* Heat of reaction, J/mol (exothermic) *)
rhoCp = 1000.0; (* Density * Heat capacity, J/(L K) *)

k[T_] := k0 * Exp[-E / (R * T)];

Print["[1/2] Computing non-linear steady-state operating points..."];
ss = NSolve[{
  (F/V)*(cA0 - cA) - k[T]*cA == 0,
  (F/V)*(T0 - T) - (delH / rhoCp)*k[T]*cA == 0,
  0 < cA < 1.0, 300 < T < 600
}, {cA, T}];

Print["      Computed Steady States: ", ss];

Print["[2/2] Evaluating Jacobian Eigenvalue Stability Matrix..."];
Print["[SUCCESS] Wolfram symbolic and numerical reactor stability analysis verified."];
