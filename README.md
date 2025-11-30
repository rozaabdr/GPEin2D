# GPEin2D – 2D Schrödinger / Gross–Pitaevskii Split–Step Fourier Solver

MATH 519 – Scientific Computing  
Author: Roza Abdrakhmanova  

This repository contains a compact Python implementation of a 2D split–step Fourier
(Strang splitting) solver for the time–dependent Schrödinger / Gross–Pitaevskii
equation in a harmonic trap. The scripts reproduce the numerical experiments,
figures, LaTeX report and beamer presentation.

The main goal is to verify the numerical method and to study the balance between
accuracy and computational cost in a linear regime and in a strongly nonlinear regime,
using the same test problem.

---

## 1. Model

We solve the 2D Gross–Pitaevskii equation on a square domain with periodic
boundary conditions:

- Domain: Ω = [-L, L] × [-L, L] with L = 10.
- PDE (informally):

  i dψ/dt = -0.5 Δψ + V(x,y) ψ + g |ψ|^2 ψ

  where
  - ψ(t, x, y) is complex–valued,
  - V(x, y) is the harmonic trapping potential,
  - g is the nonlinearity parameter.

- Harmonic trap:

  V(x,y) = 0.5 * ω^2 * (x^2 + y^2) with ω = 0.2.

- Final time: T = 2.

Two regimes are used:

- Linear Schrödinger: g = 0
- Nonlinear GPE (strongly defocusing): g = 50

### Invariants (for verification)

Continuous invariants of the model:

- Mass: integral over the domain of |ψ|^2
- Energy: integral over the domain of

    0.5 |∇ψ|^2 + V |ψ|^2 + 0.5 g |ψ|^4

The code tracks discrete versions of mass and energy at each time step and
uses them as diagnostics to verify the implementation.

### Initial data (shared test for both regimes)

The initial condition is a normalized Gaussian wavepacket with a plane–wave phase:

- Center: x0 = -3, y0 = 0
- Width: σ = 1
- Wave numbers: kx0 = 2, ky0 = 0

The wavepacket is normalized so that the discrete mass at t = 0 is approximately 1.

The same initial data, domain and final time T are used for both g = 0 and g = 50.

---

## 2. Numerical Method

### Spatial discretization

- Grid: Nx = Ny = 64 uniform points on [-L, L] with periodic boundary conditions.
- Spatial derivatives are computed spectrally using 2D FFTs:

  - Compute ψ̂ = FFT(ψ).
  - In Fourier space:
    - -Δψ corresponds to |k|^2 ψ̂,
    - ∂ψ/∂x corresponds to i kx ψ̂,
    - ∂ψ/∂y corresponds to i ky ψ̂.

Discrete diagnostics:

- Discrete mass:

  M^n = sum over grid of |ψ^n_{j,k}|^2 * Δx * Δy

- Discrete energy:

  E^n = sum over grid of
        (0.5 |∇ψ^n|^2 + V |ψ^n|^2 + 0.5 g |ψ^n|^4) * Δx * Δy

where ∇ψ^n is obtained from spectral derivatives.

### Time stepping: split–step Fourier (Strang splitting)

We split the PDE into

- Kinetic part: A ψ = -0.5 Δψ
- Potential/nonlinear part: B(ψ) = (V + g |ψ|^2) ψ

One Strang step from t^n to t^{n+1} = t^n + Δt is:

1. Half kinetic step in Fourier space:

   ψ̂^{n+1/2} = exp(-i (Δt/4) |k|^2) * ψ̂^n

2. Full potential/nonlinear step in physical space:

   ψ* = exp(-i Δt (V + g |ψ^{n+1/2}|^2)) * ψ^{n+1/2}

3. Second half kinetic step:

   ψ̂^{n+1} = exp(-i (Δt/4) |k|^2) * ψ̂*

The method is second order in time (global error proportional to Δt^2 for smooth
solutions).

---

## 3. Experiments and Verification

Two experiment types are implemented.

### 3.1 Experiment A – Invariants and CPU cost

For each regime g in {0, 50} and each time step

Δt in {0.01, 0.04, 0.08}

the solver integrates up to T = 2 and records:

- discrete mass M^n,
- discrete energy E^n,
- CPU time for the run.

Summary of findings:

- Mass conservation:
  - For all runs, max_n |M^n - M^0| is of order 10^{-13} (essentially round–off).
- Energy:
  - g = 0: relative final error in energy is about 10^{-8} for Δt = 0.01.
  - g = 50: relative final error in energy is about 10^{-4} for Δt = 0.01,
    larger for coarser time steps, but improves when Δt is reduced.
- CPU time:
  - Scales approximately like T / Δt (more steps for smaller Δt).
  - Linear and nonlinear runs have almost the same cost, since FFTs dominate.

The plot “cpu_vs_dt.png” should be interpreted as a cost model for
different time resolutions, not as a speedup curve in the sense of
parallel algorithms.

### 3.2 Experiment B – Temporal convergence vs reference solution

For each regime:

1. Compute a refined reference solution with Δt_ref = 0.005 up to T = 2.
2. For Δt in {0.01, 0.02, 0.04, 0.08}, compute the discrete L2 error at T = 2:

   e_L2(Δt) = || ψ_Δt(⋅,T) - ψ_Δt_ref(⋅,T) ||_L2

3. Plot e_L2(Δt) against Δt in log–log scale.

Results:

- g = 0:
  - Errors are roughly:
    e_L2(0.08) ≈ 1e-4,
    e_L2(0.04) ≈ 3e-5,
    e_L2(0.02) ≈ 7e-6,
    e_L2(0.01) ≈ 1e-6.
  - Each halving of Δt reduces the error by about a factor of 4,
    which is consistent with second order in time.

- g = 50:
  - Errors are roughly:
    e_L2(0.08) ≈ 2e-2,
    e_L2(0.04) ≈ 4e-3,
    e_L2(0.02) ≈ 9e-4,
    e_L2(0.01) ≈ 2e-4.
  - Again the ratios are close to 4, so second order is observed,
    but the absolute errors are 1–2 orders of magnitude larger than in
    the linear case.

These tests verify that the code behaves as a second–order method in time
for both regimes, using a refined numerical solution as reference.

---

## 4. Linear vs Nonlinear Behaviour (Illustrations)

Using the same initial data, domain, trap and final time T:

- Linear case (g = 0):
  - The density |ψ|^2 remains localized around the trap centre
    with a bell–shaped profile.
  - The difference between coarse (Δt = 0.08) and fine (Δt = 0.01)
    solutions at T is small and localized near the wavepacket.

- Nonlinear case (g = 50):
  - The density spreads over a larger region with a lower peak
    because of the repulsive interaction term g |ψ|^2 ψ.
  - The difference between coarse and fine solutions is larger
    and more spread out in space.

These plots show that the same time step can be acceptable in the linear regime
and insufficient in the strongly nonlinear regime, even though the method itself
remains second order.
