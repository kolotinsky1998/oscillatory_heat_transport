(* Oscillatory heat transport through three nonreciprocally coupled oscillators.
   Run in Mathematica with:  Get["oscillatory_reservoir_demo.wl"]  *)

ClearAll["Global`*"];

(* Physical realization:
   Three optically levitated nanoparticles in one-dimensional traps.
   Optical feedback creates direction-dependent effective spring constants.
   Each particle is coupled to its own Langevin heat reservoir. *)

n = 3;
mass = 1.0;
gamma = 0.1;
xi = 1.0;

M = mass IdentityMatrix[n];
L = gamma IdentityMatrix[n];

(* Force-interaction matrix from the three-site nonreciprocal closed chain.
   The entries F[[i,j]] are the coefficients in
      M q''(t) + F q(t) + L q'(t) == eta(t).
   xi closes the chain; xi = 0 is open, xi = 1 is periodic. *)
F = N[{{1.0, 0.14, 0.06 xi},
       {0.06, 1.0, 0.14},
       {0.14 xi, 0.06, 1.0}}];

Print["Force-interaction matrix F:"];
Print[MatrixForm[F]];

(* First-order deterministic oscillator dynamics:
      d/dt {q, v} = Aosc {q, v}.
   Stability with friction means Re[Eigenvalues[Aosc]] < 0. *)
Aosc = ArrayFlatten[{{ConstantArray[0, {n, n}], IdentityMatrix[n]},
                     {-Inverse[M].F, -Inverse[M].L}}];

oscillatorEigenvalues = Eigenvalues[Aosc];
Print["Eigenvalues of deterministic oscillator matrix Aosc:"];
Print[N[oscillatorEigenvalues, 8]];
Print["Oscillator stable with friction? ",
      Max[Re[oscillatorEigenvalues]] < 0];

ClearAll[CovarianceForTemperatures, EffectiveTemperatureResponse, HeatTransferMatrix];

CovarianceForTemperatures[fmat_, lmat_, tvec_] := Module[
  {nn, aa, noise, cvars, equations, vars, arrays, solution},
  nn = Length[fmat];
  aa = ArrayFlatten[{{ConstantArray[0, {nn, nn}], IdentityMatrix[nn]},
                     {-fmat, -lmat}}];
  noise = ConstantArray[0, {2 nn, 2 nn}];
  Do[
    noise[[nn + i, nn + i]] = 2 lmat[[i, i]] tvec[[i]],
    {i, nn}
  ];
  cvars = Array[c, {2 nn, 2 nn}];
  equations = Flatten[aa.cvars + cvars.Transpose[aa] + noise];
  vars = Flatten[cvars];
  arrays = Normal[CoefficientArrays[equations, vars]];
  solution = LinearSolve[arrays[[2]], -arrays[[1]]];
  Partition[solution, 2 nn]
];

EffectiveTemperatureResponse[fmat_, lmat_] := Module[
  {nn, basis, covariances},
  nn = Length[fmat];
  basis = IdentityMatrix[nn];
  covariances = CovarianceForTemperatures[fmat, lmat, #] & /@ basis;
  Transpose[
    Table[
      Diagonal[covariances[[j]][[nn + 1 ;; 2 nn, nn + 1 ;; 2 nn]]],
      {j, nn}
    ]
  ]
];

HeatTransferMatrix[fmat_, lmat_] := Module[
  {nn, teffResponse},
  nn = Length[fmat];
  teffResponse = EffectiveTemperatureResponse[fmat, lmat];
  IdentityMatrix[nn] - teffResponse
];

S = N[HeatTransferMatrix[F, L]];

Print["Heat-transfer matrix S, defined by J = gamma S T:"];
Print[MatrixForm[S]];

sEigenvalues = Eigenvalues[S];
Print["Eigenvalues of S:"];
Print[N[sEigenvalues, 10]];
Print["Does S contain an oscillatory complex pair? ",
      Max[Abs[Im[sEigenvalues]]] > 10^-8];

(* Reservoir model.
   Without additional heat outflow:
      T'(t) = kappa gamma S T(t).
   Since the nonconservative oscillator chain gives Re[lambda(S)] > 0,
   the reservoir temperatures contain exponentially growing modes.

   Add heat leakage from each reservoir to an ambient thermostat:
      T'(t) = kappa gamma S T(t) - R.(T(t) - Tamb 1).
   Taking R = r I shifts all reservoir eigenvalues by -r and suppresses
   exponential growth while preserving the imaginary parts responsible for
   oscillations. *)

kappa = 20.0;
drainRate = 0.98;
R = drainRate IdentityMatrix[n];
Tamb = 1.0;

Areservoir = kappa gamma S - R;
reservoirEigenvalues = Eigenvalues[Areservoir];

Print["Reservoir evolution matrix kappa gamma S - R:"];
Print[MatrixForm[N[Areservoir, 8]]];
Print["Eigenvalues of stabilized reservoir matrix:"];
Print[N[reservoirEigenvalues, 10]];
Print["Reservoir temperatures bounded after additional outflows? ",
      Max[Re[reservoirEigenvalues]] < 0];
Print["Oscillatory reservoir mode remains? ",
      Max[Abs[Im[reservoirEigenvalues]]] > 10^-8];

tmax = 160;
initialTemperatures = {1.45, 0.95, 0.65};
temperatureVariables = Array[temp, n];
temperatureFunctions = Table[temp[i][t], {i, n}];

temperatureEquations = Join[
  Thread[
    D[temperatureFunctions, t] ==
      Areservoir.temperatureFunctions +
        R.ConstantArray[Tamb, n]
  ],
  Thread[Table[temp[i][0], {i, n}] == initialTemperatures]
];

temperatureSolution = NDSolve[
  temperatureEquations,
  temperatureVariables,
  {t, 0, tmax}
][[1]];

temperaturePlot = Plot[
  Evaluate[temperatureFunctions /. temperatureSolution],
  {t, 0, tmax},
  PlotRange -> All,
  PlotLegends -> Placed[{"T1", "T2", "T3"}, Right],
  Frame -> True,
  FrameLabel -> {"time", "reservoir temperature"},
  PlotLabel -> "Oscillatory heat transport with stabilizing heat outflows",
  ImageSize -> Large
];

deviationPlot = Plot[
  Evaluate[
    (temperatureFunctions /. temperatureSolution) -
      ((temperatureFunctions /. temperatureSolution) /. t -> tmax)
  ],
  {t, 0, tmax},
  PlotRange -> All,
  PlotLegends -> Placed[{"T1 - T1(end)", "T2 - T2(end)", "T3 - T3(end)"}, Right],
  Frame -> True,
  FrameLabel -> {"time", "temperature deviation"},
  PlotLabel -> "Damped oscillatory components",
  ImageSize -> Large
];

Print[temperaturePlot];
Print[deviationPlot];

Export["oscillatory_reservoir_temperatures.png", temperaturePlot];
Export["oscillatory_reservoir_deviations.png", deviationPlot];
