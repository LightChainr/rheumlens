# Design-label non-identifiability

## Setup

Let `Y` be a donor-level disease label, `D` a recorded design block and `X` a
donor representation. The scientific target is not merely prediction of `Y`
from `X`, but attribution of the disease-aligned part of `X` to biology rather
than acquisition design.

## Proposition 1: perfect collinearity

Suppose `D = Y` and

`X = beta * Y + delta * D + epsilon`.

The observed model is

`X = (beta + delta) * Y + epsilon`.

For every scalar `c`, the parameter pairs `(beta, delta)` and
`(beta + c, delta - c)` induce the same observed distribution of `(X, D, Y)`.
The biological and technical contributions are therefore not separately
identifiable. This is an observational equivalence result, not a statement
that either contribution is absent.

Residualising `X` on `D` removes the combined design-aligned component.
Restriction within `D` leaves no variation in `Y`. Neither operation can
recover the missing counterfactual cells of the design-by-label table.

## Proposition 2: partial overlap in an additive linear model

Consider the linear model

`X_j = alpha_j + beta_j * Y + D * gamma_j + epsilon_j`

for representation coordinate `j`. Let `M_D` be the residual-maker matrix for
the intercept and encoded design block. By the Frisch-Waugh-Lovell theorem, the
information for `beta_j` is proportional to

`Y' M_D Y`.

Relative to an intercept-only model, define

`I_D = (Y' M_D Y) / (Y' M_1 Y) = 1 - R^2(Y ~ D)`.

Under homoskedastic errors, the conditional variance of the adjusted disease
coefficient is proportional to `1 / (Y' M_D Y)`. Relative variance inflation is
therefore `1 / I_D`. At `I_D = 0`, the disease and design coefficients are not
separately estimable.

This result does not imply that downstream AUC loss is a universal monotone
function of design-only AUC. AUC also depends on feature geometry, nonlinear
signal, technical-effect magnitude, estimator regularisation and finite-sample
variation.

## Empirical interpretation

The paper should report four distinct quantities:

1. design-only disease predictability;
2. representation-to-design predictability;
3. information and support remaining after conditioning on design;
4. source-only performance under an independent target design.

Residualisation asks what remains after deleting design-predictable coordinates.
Restriction asks what is predictable in observed overlap strata. External
transfer asks whether a source-fitted ranking survives a new design. These are
not interchangeable estimands.
