# GPEin2D – 2D Schrödinger / Gross–Pitaevskii Split–Step Fourier Solver

MATH 519 – Scientific Computing  
Author: Roza Abdrakhmanova  

This repository contains a compact Python implementation of a 2D split–step Fourier
(Strang splitting) solver for the time-dependent Schrödinger / Gross–Pitaevskii
equation in a harmonic trap, together with scripts to reproduce the numerical
experiments, figures, LaTeX report and beamer presentation.

The goal is not only to “simulate something”, but to **verify** the numerical
method and quantify the **accuracy–cost trade-off** in both a linear and a
strongly nonlinear regime on the **same test problem**.

---

## 1. Model

We solve, on a square domain \(\Omega=[-L,L]^2\),
the 2D Gross–Pitaevskii equation with periodic boundary conditions:
\[
i \partial_t \psi(t,x,y)
= -\tfrac12 \Delta \psi(t,x,y)
  + V(x,y)\psi(t,x,y)
  + g|\psi(t,x,y)|^2\psi(t,x,y).
\]

- Domain half-width: \(L = 10\).
- Harmonic trap:
  \[
  V(x,y) = \tfrac12 \omega^2(x^2 + y^2), \quad \omega = 0.2.
  \]
- Final time: \(T = 2\).
- Regimes:
  - Linear Schrödinger: \(g = 0\).
  - Nonlinear GPE: \(g = 50\) (strongly defocusing).

### Invariants

Continuous invariants used for verification:

- Mass:
  \[
  M(t) = \int_{\Omega} |\psi(t,x,y)|^2\,dx\,dy.
  \]
- Energy:
  \[
  E(t) = \int_{\Omega} \Big(
      \tfrac12|\nabla\psi|^2 + V|\psi|^2 + \tfrac12 g|\psi|^4
    \Big)\,dx\,dy.
  \]

The exact dynamics conserves both \(M\) and \(E\) for \(g=0\) and \(g>0\).  
The code monitors discrete versions \(M^n, E^n\) to verify the implementation.

### Initial data (shared test for both regimes)

A normalized Gaussian wavepacket with phase:
\[
\psi_0(x,y)
= C\exp\! \Big(
   -\frac{(x-x_0)^2 + (y-y_0)^2}{2\sigma^2}
  \Big)
  \exp\big(i(k_{x0}x + k_{y0}y)\big),
\]
with
- \(x_0=-3\), \(y_0=0\),
- \(\sigma=1\),
- \(k_{x0}=2\), \(k_{y0}=0\).

The constant \(C\) is chosen so that the **discrete** mass \(M^0\approx 1\).
The same \(\psi_0\), domain and \(T\) are used for \(g=0\) and \(g=50\).

---

## 2. Numerical Method

### Spatial discretization

- Grid: \(N_x=N_y=64\) uniform points on \([-L,L]^2\), periodic.
- Spectral (pseudo-spectral) differentiation via 2D FFT:
  - Compute \(\widehat{\psi}=\mathcal{F}[\psi]\).
  - Apply multipliers in Fourier space:
    - \(-\Delta \psi \leftrightarrow |k|^2\widehat{\psi}\),
    - \(\partial_x\psi \leftrightarrow ik_x\widehat{\psi}\),
    - \(\partial_y\psi \leftrightarrow ik_y\widehat{\psi}\).

Discrete mass:
\[
M^n = \sum_{j,k} |\psi^n_{j,k}|^2\,\Delta x\,\Delta y.
\]

Discrete energy:
\[
E^n = \sum_{j,k}
\Big(
\tfrac12|\nabla_h\psi^n|^2_{j,k}
+ V_{j,k}|\psi^n_{j,k}|^2
+ \tfrac12 g|\psi^n_{j,k}|^4
\Big)\Delta x\Delta y,
\]
where \(\nabla_h\psi^n\) is computed spectrally.

### Time stepping: split–step Fourier (Strang splitting)

We write
\[
i\partial_t\psi = (A + B(\psi))\psi,
\]
with
- \(A\psi = -\tfrac12\Delta\psi\) (kinetic),
- \(B(\psi) = (V+g|\psi|^2)\) (potential + nonlinearity).

One Strang step from \(t^n\) to \(t^{n+1}=t^n+\Delta t\):
\[
\psi^{n+1} \approx
e^{-i\frac{\Delta t}{2}A}\,
e^{-i\Delta t B(\psi^n)}\,
e^{-i\frac{\Delta t}{2}A}\,\psi^n.
\]

Implemented as:
1. Half kinetic step in Fourier space:
   \[
   \widehat{\psi}^{n+1/2}
   = \exp\!\big(-i\tfrac{\Delta t}{4}|k|^2\big)\,\widehat{\psi}^n.
   \]
2. Full potential/nonlinear step in physical space:
   \[
   \psi^{*}(x,y)
   = \exp\!\big(-i\Delta t(V(x,y)+g|\psi^{n+1/2}|^2)\big)\,\psi^{n+1/2}(x,y).
   \]
3. Second half kinetic step:
   \[
   \widehat{\psi}^{n+1}
   = \exp\!\big(-i\tfrac{\Delta t}{4}|k|^2\big)\,\widehat{\psi}^{*}.
   \]

The method is formally **second order in time** for smooth solutions.

---

## 3. Experiments and Verification

Two experiment types are implemented.

### 3.1 Experiment A – Invariants and CPU cost

For each regime \(g\in\{0,50\}\) and
\[
\Delta t\in\{0.01,0.04,0.08\}
\]
we integrate up to \(T=2\) and record:
- discrete mass \(M^n\),
- discrete energy \(E^n\),
- CPU time for one run.

Findings:
- Mass conservation:  
  \(\max_n |M^n-M^0|\lesssim 10^{-13}\) for all runs (round-off level).
- Energy:
  - \(g=0\): relative final energy error \(\sim 10^{-8}\).
  - \(g=50\): relative final energy error \(\sim 10^{-4}\) for \(\Delta t=0.01\),
    larger for coarser \(\Delta t\), but monotone improvement as \(\Delta t\) decreases.
- CPU time:
  - Scales approximately like \(T/\Delta t\).
  - Linear and nonlinear runs have nearly identical cost (FFTs dominate).

**Important:** The CPU vs \(\Delta t\) plot is a **cost model at different time
resolutions**, not an algorithmic speedup curve. No parallelization is used.

### 3.2 Experiment B – Temporal convergence vs reference solution

For each regime:

1. Compute a refined reference solution
   with \(\Delta t_{\text{ref}}=0.005\) up to \(T=2\).
2. For \(\Delta t\in\{0.01,0.02,0.04,0.08\}\), compute
   \[
   e_{L^2}(\Delta t)
   = \big\|\psi_{\Delta t}(\cdot,T)
         -\psi_{\Delta t_{\text{ref}}}(\cdot,T)\big\|_{L^2}.
   \]
3. Plot \(e_{L^2}(\Delta t)\) vs \(\Delta t\) in log–log scale.

Results:

- \(g=0\):  
  \(e_{L^2}(0.08)\approx 10^{-4}\),  
  \(e_{L^2}(0.04)\approx 10^{-5}\),  
  \(e_{L^2}(0.02)\approx 10^{-6}\).  
  Ratios \(\approx 4\) → observed order \(\approx 2\).

- \(g=50\):  
  \(e_{L^2}(0.08)\approx 2\times 10^{-2}\),  
  \(e_{L^2}(0.04)\approx 4\times 10^{-3}\),  
  \(e_{L^2}(0.02)\approx 9\times 10^{-4}\),  
  \(e_{L^2}(0.01)\approx 2\times 10^{-4}\).  
  Again ratios \(\approx 4\) → second order.

Thus, the implementation passes a standard **verification test**:
second-order time convergence against a more refined numerical solution
for both linear and nonlinear regimes.

---

## 4. Linear vs Nonlinear Tests (Illustrations)

Using the **same test problem** (domain, trap, initial data, \(T\)):

- In the linear case (\(g=0\)):
  - The packet stays tightly localized under the trap,
    with a bell-shaped density.
  - Coarse vs fine time-step error is localized near the core
    and very small (\(\sim 10^{-5}\)).

- In the nonlinear case (\(g=50\)):
  - The density spreads over a broader region with lower peak value
    (repulsive interaction).
  - Coarse vs fine error is larger (\(\sim 10^{-3}\)) and more spread out.

This shows that the **same** \(\Delta t\) is much more damaging in the nonlinear case,
even though the scheme remains second order.

