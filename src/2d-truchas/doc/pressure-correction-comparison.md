% Pressure correction in free-surface flow
% 2D Truchas working note
% 2026-09-07

# Purpose and scope

This note compares the pressure-correction paths currently used by the
two-dimensional flow code and by mainline Truchas. The comparison is focused
on the cell-centered velocity correction in a fluid/VOID calculation. The
two implementations use closely related pressure systems and face-centered
corrections; the important difference is how the pressure correction is
converted into a cell-centered velocity.

The issue was exposed by a translating fluid plug on a perturbed mesh at a
small Courant number. A mixed fluid/VOID cell retained a very small amount of
real fluid. The pressure solve converged, but the resulting cell-centered
velocity correction became enormous.

# Discrete projection setting

Let \(c\) denote a cell and \(f\) a face. The projection computes a pressure
correction \(\delta p_c\) from a discrete incompressibility equation

\[
  A\,\delta p = b,
\]

where \(A\) is assembled from face coefficients. In both implementations the
coefficient used on a face is the inverse of a face density,

\[
  \alpha_f = \frac{1}{\rho_f}.
\]

Schematically, the matrix is the finite-volume composition

\[
  A \simeq -D\,M_f(\alpha_f)\,G,
\]

with \(G\) a pressure-to-face-gradient operator, \(M_f\) the geometric face
scaling, and \(D\) a face-flux divergence operator. The exact geometric
weights are implementation-specific, but the important point is that the
pressure system is face-based and uses \(\alpha_f\).

The pressure correction is then used in two places:

1. to correct the face-normal velocity used by the divergence constraint; and
2. to update the cell-centered velocity used by the rest of the flow
   algorithm and by output.

The first operation is the one that directly enforces the projection. The
second should represent the same physical correction in a cell-centered form.

# Common face correction

Both codes apply the face correction in the form

\[
  u_f^{n+1} = u_f^* - \Delta t\,\alpha_f\,
      (G_f\,\delta p)_f,
\]

where \(u_f^*\) is the predicted face-normal velocity and \(G_f\) includes
the boundary-condition treatment for pressure correction. In the 2D code
this is the pressure_derivative/inv_density_f path. In mainline it is the
velocity_fc_correct path using gradient_cf and rho_fc.

Consequently, the pressure system and the face velocity correction already
use the same basic design in the two codes. The observed instability is not
evidence by itself that the pressure matrix is nonsymmetric or that the
projection solve failed: in the diagnostic run the projection residual was
about \(10^{-12}\), and the solve reported success.

# The current 2D cell-centered correction

The current 2D projection update stores a cell-centered pressure-gradient-like
quantity at the old and new states. Suppressing boundary and gravity details,
the update has the form

\[
  \boldsymbol u_c^{n+1}
    = \boldsymbol u_c^*
      - \Delta t\,\alpha_c
        \left[\left(\nabla_c p\right)^{n+1}
              -\left(\nabla_c p\right)^n\right],
\]

where

\[
  \alpha_c = \frac{1}{\rho_c}
\]

and \(\nabla_c\) is the 2D cell-centered gradient_cc operator. In the source,
the operation is represented by gradient_pressure followed by a direct
multiplication by inv_density_c.

This is a reasonable approximation for a single-density fluid. If density is
smooth and bounded away from zero, the distinction between

\[
  \frac{1}{\rho_c}\nabla_c p
  \quad\text{and}\quad
  I_{fc}\left(\frac{1}{\rho_f}G_f p\right)
\]

is an ordinary discretization difference. For a fluid/VOID interface,
however, \(\rho_c\) represents the mass density averaged over the active cell
volume. A cell containing a small real-fluid fraction \(\theta\) can therefore
have, schematically,

\[
  \rho_c \simeq \theta\,\rho_{fluid},
  \qquad
  \alpha_c \simeq \frac{1}{\theta\,\rho_{fluid}}.
\]

The cell-centered correction then amplifies an otherwise modest pressure
gradient by \(1/\theta\).

The current free-surface diagnostic gives a concrete example:

\[
  \theta \simeq 3.78\times10^{-8},
  \qquad
  \alpha_c \simeq 2.65\times10^7.
\]

At the first anomalous cell, the projection converged with a relative residual
of approximately \(9.9\times10^{-13}\), but the cell-centered correction was
approximately

\[
  (20.17,\;1.069).
\]

The large velocity is thus a correction-amplification problem after a
successful pressure solve, not a failure of the linear solver to reduce its
residual.

# Mainline's face-density-scaled cell correction

Mainline forms a pressure-gradient quantity already divided by face density.
For each active face it computes, schematically,

\[
  \boldsymbol q_f
    = \frac{1}{\rho_f}
      \left(G_f p + \boldsymbol h_f\right),
\]

where \(\boldsymbol h_f\) represents the gravity-head contribution when
gravity is active. It then interpolates that face quantity to the cell center,

\[
  \boldsymbol q_c = I_{fc}\boldsymbol q_f,
\]

and updates the cell-centered velocity as

\[
  \boldsymbol u_c^{n+1}
    = \boldsymbol u_c^* - \Delta t
      \left(\boldsymbol q_c^{n+1}-\boldsymbol q_c^n\right).
\]

There is no additional multiplication by \(1/\rho_c\) in this final update.
In the source this is the sequence grad_p_rho, interpolate_fc, and
velocity_cc_correct.

Mainline also constructs \(\rho_f\) from neighboring cell properties using
fluid-volume weighting and enforces a positive face-density floor. In
schematic form,

\[
  \rho_f = \max\left(\rho_{f,\mathrm{weighted}},
                       \rho_{f,\min}\right).
\]

The floor is based on the minimum real-fluid density and a configurable face
fraction. Therefore a tiny fluid fragment does not produce an arbitrarily
large inverse density in the pressure correction used by the face and
cell-centered velocity paths.

# Comparison

The essential difference can be summarized as follows.

| Aspect | Current 2D path | Mainline path |
| --- | --- | --- |
| Projection matrix | Face-based; uses \(1/\rho_f\) | Face-based; uses \(1/\rho_f\) |
| Face velocity correction | Face gradient divided by face density | Face gradient divided by face density |
| Cell velocity correction | Cell gradient multiplied by \(1/\rho_c\) | Face gradient divided by \(\rho_f\), then interpolated to cells |
| Small fluid/VOID fraction | Can make \(1/\rho_c\) very large | Controlled by face-density weighting and floor |
| Primary consistency | Cell and face paths use different density/gradient pairings | Cell path is derived from the same face-scaled quantity as the face path |

For a full fluid cell with smooth density, the two cell-centered formulas are
close and may be indistinguishable at the expected discretization error. At
a free surface, they are not equivalent: the current 2D formula applies a
cell-centered inverse density after the gradient has been formed, while the
mainline formula applies the density scaling on faces before interpolation.

# Cell classification is not the root cause

The problematic cell is a mixed fluid/VOID cell adjacent to pure VOID. Both
implementations classify this configuration as an active regular cell rather
than simply deleting it. The reason is that a mixed cell can still exchange
fluid with its neighboring active cells, and the interface reconstruction
needs to retain it.

Thus changing the classification alone would be a poor fix. Promoting the
cell to VOID would discard a geometrically represented fluid fragment; keeping
it active is consistent with the mainline classification contract. The
important question is how pressure and velocity corrections are scaled in
that active cell.

# Implications for the 2D implementation

The most direct alignment with mainline is to replace the current cell
correction by a 2D analogue of the face-density-scaled path:

1. compute the pressure-gradient contribution on faces;
2. divide by the bounded face density, including the free-surface floor;
3. interpolate the resulting vector to cell centers; and
4. subtract the difference between the new and old interpolated quantities.

The face correction and the projection matrix can remain as they are initially.
This change targets the demonstrated amplification mechanism without changing
the material tracker or silently deleting small fluid fragments.

Simply clamping inv_density_c would suppress the symptom, but it would leave
the cell-centered correction mathematically disconnected from the face
correction and would introduce a second, less transparent free-surface
regularization. It may still be useful as a diagnostic comparison, but it is
not the preferred final formulation.

# Limitations and follow-up questions

The face-density-scaled correction resolves the immediate mismatch between the
2D and mainline formulations, but it does not by itself settle every
free-surface question. In particular, future work should verify:

- the appropriate face-density floor for 2D geometries and resolution;
- consistency between the cell-centered velocity used by momentum prediction
  and the face velocity used by the divergence constraint;
- behavior when several disconnected fluid fragments are present;
- conservation and normalization of fluid/VOID volume fractions after
  geometric transport; and
- the treatment of pressure and velocity at free-surface boundaries.

The current evidence supports changing the cell-centered pressure correction,
not changing the active-cell classification or weakening the linear-solver
convergence criterion.
