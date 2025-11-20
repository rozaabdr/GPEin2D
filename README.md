This code solves the 2D Schrödinger / Gross–Pitaevskii equation on a square periodic grid with a harmonic trap. The same solver is run for two parameter choices:
- **Linear Schrödinger**: `g = 0` - no interaction
- **Nonlinear Gross–Pitaevskii (GPE)**: `g = 50` - repulsive interaction

In both cases the initial condition is a normalized Gaussian wavepacket.
Time-stepping is performed with a 2D split–step Fourier (Strang splitting) method.


## 1. Background and context

### Schrödinger equation vs Gross–Pitaevskii equation

The time-dependent Schrödinger equation is the basic evolution law of nonrelativistic quantum mechanics. In the simplest form,

```math
i \partial_t \psi = \left(-\tfrac12 \Delta + V(x,y)\right)\psi ,
```

it describes the wavefunction $\psi(t,x,y)$ of a single particle moving in a potential $V(x,y)$.

The Gross–Pitaevskii equation (GPE) is a nonlinear variant used for Bose–Einstein condensates (BECs). At very low temperatures, many bosons occupy the same quantum state and can be modeled by a single macroscopic wavefunction whose dynamics is

```math
i \partial_t \psi = \left(-\tfrac12 \Delta + V(x,y) + g\lvert\psi\rvert^2\right)\psi .
```

The extra term $g\lvert\psi\rvert^2\psi$ encodes mean-field interactions between particles:

- $g = 0$ → purely linear Schrödinger dynamics.
- $g > 0$ → repulsive condensate; high density costs energy and the cloud tends to spread.

This project does not try to be physically exact, but it uses a standard GPE-type model which is widely studied in mathematical physics and computational PDEs.

### Harmonic trap and Gaussian wavepacket

Real BECs are typically held in magnetic or optical harmonic traps, which are well-approximated by

```math
V(x,y) = \tfrac12 \omega^2 (x^2 + y^2).
```

Here we choose a mild trap strength $\omega = 0.2$ so that the wavepacket stays inside the computational domain but still “feels” the confining potential.

The initial condition is a Gaussian wavepacket, i.e. a localized bump with a built-in plane-wave phase:

- localized near $(x_0,y_0) = (-3,0)$,
- width $\sigma = 1$,
- carrier wavenumber $(k_{x0},k_{y0}) = (2,0)$.

This is a classical test case in the Schrödinger/GPE literature because:

- The Gaussian is smooth and rapidly decaying → well-suited to spectral methods.
- The packet is initially displaced from the center of the trap → it moves and oscillates in time, so we see nontrivial dynamics instead of a static picture.
- For $g = 0$, the qualitative behavior (oscillation in the harmonic well) is well understood, providing a sanity check for the numerics.

### Rationale of choosing these parameters:

- **Domain size:** $L = 10$. With the chosen trap and initial packet, essentially all significant mass lies in $[-10,10]^2$; periodic boundary conditions are therefore harmless.
- **Grid:** $64 \times 64$. This is enough to resolve the Gaussian and the potential while keeping run-times below one second on a laptop.
- **Time steps:** $\Delta t \in \{0.01, 0.04, 0.08\}$ up to $T = 2.0$.  
  These values show a clear accuracy–cost trade-off without making the scheme unstable.
- **Nonlinearity strength:** $g = 50$. This is strong enough that the nonlinear case is visibly different from the linear one (the condensate spreads and flattens), but not so strong that the numerics become delicate.

 this setup is a compact prototype of a dispersive nonlinear equation with:
- smooth solution,
- nontrivial potential,
- conserved invariants (mass and energy),
- and parameter regimes where linear and nonlinear dynamics can be compared.



## 2. PDE model

Solve

```math
i \partial_t \psi(t,x,y)
= \left( -\tfrac12 \Delta + V(x,y) + g \lvert\psi\rvert^2 \right) \psi,
\quad (x,y)\in[-L,L]^2,\ L=10,
```

with

```math
V(x,y) = \tfrac12 \omega^2 (x^2 + y^2), \quad \omega = 0.2.
```

Boundary conditions: periodic on the square computational domain.

Initial data: normalized Gaussian wavepacket centered at $(x_0,y_0)=(-3,0)$ with width $\sigma=1$ and carrier wavenumber $(k_{x0},k_{y0})=(2,0)$.

Invariants of interest:

```math
M(t) = \iint \lvert\psi\rvert^2 \, dx\,dy,
```

```math
E(t) = \iint \Big[\tfrac12 \lvert\nabla\psi\rvert^2 + V\lvert\psi\rvert^2
      + \tfrac{g}{2}\lvert\psi\rvert^4\Big] dx\,dy.
```


## 3. Numerical method

- Spatial discretization: uniform grid $64 \times 64$ on $[-L, L]^2$, periodic.
- Derivatives via spectral differentiation using FFT (`numpy.fft`).
- Time discretization: split–step Fourier (also called strang splitting):

  1. Half-step kinetic in Fourier space,
  2. Full-step potential + nonlinearity in physical space,
  3. Half-step kinetic in Fourier space.

- Time steps: $\Delta t = 0.01, 0.04, 0.08$ up to final time $T = 2.0$.

For each $(g, \Delta t)$ we record:

- mass history $M(t)$,
- energy history $E(t)$,
- CPU time,
- maximum mass error $\max_t \lvert M(t) - M(0)\rvert$,
- energy drift $\lvert E(T) - E(0)\rvert$.



## 4.Programming language and libraries:

- Python 3
- `numpy`
- `matplotlib`


## 5. Main files

- `project_gpe2d.py` – main script: grid setup, solver implementation, experiments, plotting, CSV export.
- `results_summary.txt` – human-readable summary table (CPU time, mass/energy errors).
- `results_summary.csv` – same data in CSV format for post-processing.
- `density_*.png`, `mass_*.png`, `energy_*.png` – figures used in the report and Beamer slides.


## 6. Results (short)

- Initial mass: `M(0) ≈ 1.0`.
- For all runs, `max |M(t) - M(0)|` is at machine precision level (`10^{-14}`), so the method is numerically mass-conservative.
- In the **linear** regime (`g = 0`) the energy drift stays below `5.3 × 10^{-6}` even for the largest time step `dt = 0.08`.
- In the **nonlinear** regime (`g = 50`) the energy drift grows with `dt`:
  - ≈ `9.7 × 10^{-4}` for `dt = 0.01`,
  - ≈ `1.6 × 10^{-2}` for `dt = 0.04`,
  - ≈ `7.7 × 10^{-2}` for `dt = 0.08`.

This illustrates the accuracy–cost trade-off: smaller `dt` yields better conservation of invariants, larger `dt` is faster but less accurate.
