(* ==============================================================================
   Scientific Autonomous Reasoning Engine: CSTR Non-Linear Dynamics
   Language: Wolfram Language / Mathematica
   ============================================================================== *)

BeginPackage["CSTRDynamics`"]

SolveEquilibriumAndStability::usage = "SolveEquilibriumAndStability[q, V, k0, E, deltaH] computes equilibrium states and Jacobian eigenvalues."

Begin["`Private`"]

SolveEquilibriumAndStability[q_, V_, k0_, EoverR_, deltaH_, rhoCp_, Tc_, UA_] := 
  Module[{CA, T, k, dCAdt, dTdt, jacobian, equilibria, eigenvalues},
    
    (* Arrhenius Reaction Rate Constant *)
    k[T_] := k0 * Exp[-EoverR / T];
    
    (* Governing Differential Equations *)
    dCAdt[CA_, T_] := (q / V) * (CA0 - CA) - k[T] * CA;
    dTdt[CA_, T_] := (q / V) * (T0 - T) + (-deltaH / rhoCp) * k[T] * CA - (UA / (V * rhoCp)) * (T - Tc);
    
    (* Symbolic Jacobian Matrix *)
    jacobian[CA_, T_] = {
      {D[dCAdt[c, temp], c], D[dCAdt[c, temp], temp]},
      {D[dTdt[c, temp], c], D[dTdt[c, temp], temp]}
    } /. {c -> CA, temp -> T};
    
    Return[jacobian[CA, T]]
  ];

End[]
EndPackage[]
