# GraphBasedPredictors

Author: Robin Richter

Institution: Deutsches Zentrum für Neurodegenerative Erkrankungen (DZNE) e.V.

Licence: MIT Licence

Note: Please feel free to ask questions or report bugs.

R functions to compute graph based predictors (GBPs) as baselines for causal structural learning as proposed in the paper RSS23:

- The function oip.R computes observed indegree predictors in all it four forms (OIP, T-OIP, P-OIP and P-T-OIP);
- The function tap.R computes the transitivity assuming predictor (TAP) or its realted q-TAP by Algorithm 1 of [RSS23];
- The function btap.R computes the biased version of the TAP and q-TAP given by Algorithm 2 of [RSS23].

The remaining functions are needed as subroutines.

[RSS23]: Richter, R., Bhamidi, S. & Mukherjee, S. Improved baselines for causal structure learning on interventional data. Stat Comput 33, 93 (2023). https://doi.org/10.1007/s11222-023-10257-9
