# blowup-script

`blowup-script` is an interactive Rust CLI for local computations with
projective hypersurfaces over `Q` of dimension 1, 2, and 3.

It works with a homogeneous polynomial in:

- `P^2` for curves, using coordinates `x0, x1, x2`
- `P^3` for surfaces, using coordinates `x0, x1, x2, x3`
- `P^4` for threefolds, using coordinates `x0, x1, x2, x3, x4`

The program builds standard affine charts, computes the Jacobian equations for
the singular locus on each active chart, decomposes the singular locus with an
exact Groebner-basis routine, finds rational intersection points of singular
components, analyzes local point singularity types from multiplicity and tangent
cone data, and reports generic transverse types such as `cA3` along affine
linear singular lines. It lets you blow up a singular component, component
intersection, rational point, or coordinate center, then updates the active
charts. Results can be saved as JSON with schema `blowup-script/v1`, including
active charts, blowup history, current singular-locus components and
intersections, generic line type data, local point type data, and the raw
Jacobian equations with rational point samples.

## Build

```sh
cargo build
```

## Run

```sh
cargo run
```

Example curve:

```text
dimension: 1
polynomial: x0*x2^2 - x1^3
```

Useful commands inside the CLI:

```text
singular
intersections
charts
analyze point 0 0 0
blowup component 1
blowup intersection 1
blowup option 0
blowup point 0 0 0
blowup point 3 x1=0 x2=0 --force
blowup center 0 x1=0 x2=0
resolve
resolve crepant --max-steps 8
set-bound 3
undo
save result.json
quit
```

In an interactive terminal, the command prompt supports Up/Down command history
and Left/Right cursor movement. `undo` restores the state before the previous
state-changing command, currently blowups, automatic resolution, and `set-bound`.

## Automatic Crepant Resolution

Use `resolve` or `resolve crepant` to let the program repeatedly choose and
blow up supported crepant centers for the double cover branched along the input
curve or surface. The resolver exhausts positive-dimensional crepant components
before considering rational component intersections or isolated components. A center is
used only when the double-cover branch discrepancy test is crepant: if the
center has ambient affine codimension `c` and the branch divisor has
multiplicity `m` along it, then `c - 1 - floor(m / 2)` must be zero. Thus, for
branch surfaces in `P^3`, automatic crepant centers are points of multiplicity
4 or 5 and lines of multiplicity 2 or 3.

Manual `blowup ...` commands still compute ordinary strict transforms of the
input hypersurface. The automatic resolver uses the branch-divisor transform:
for odd branch multiplicity it keeps the exceptional divisor in the branch
locus.

The command stops when the transformed branch divisor is smooth, when the step
limit is reached, or when singularities remain but every supported center would
be non-crepant. The default step limit is 32; override it with:

```text
resolve --max-steps 12
```

Manual blowups remain available through the existing `blowup ...` commands.
When the affine-linear curve search is exhausted and
`BLOWUP_SINGULAR_BACKEND=1` is set, the resolver performs an exact
positive-dimensional minimal-prime audit before it considers isolated points.
Nonlinear curve components are therefore reported as unresolved instead of
being miscounted as isolated rational samples. A `save-raw` file written at a
step limit is a checkpoint only; its empty singular-component arrays are not a
smoothness certificate.

## Supported Blowups

The implemented blowup charts are exact symbolic strict transforms for affine
coordinate centers of the form:

```text
var_1 = rational_1, ..., var_k = rational_k
```

This covers rational singular points, arbitrary rational point centers on later
active charts, plus coordinate or affine-linear curves/linear centers in the
current affine chart. `blowup component N` applies the blowup to every affine
chart listed for the displayed singular component, so a projective line or point
that appears in multiple standard charts is handled compatibly. `blowup
intersection N` blows up a displayed rational intersection point of singular
components on its chart. The program validates that each center is contained in
the corresponding chart's singular locus unless `--force` is supplied.

The singular locus display is grouped by prime component:

```text
1: dimension 1, generic multiplicity 2
   Chart 2 poly: f0 = x0, f1 = x1
   Chart 3 poly: f0 = x0, f1 = x1
```

For each component the display includes a dimension, generic multiplicity of the
hypersurface along the component, and one or more chart ideals. If a line has a
finite transverse `A_n` type at a generic rational sample, the component reports
that as `cA_n`. This is a local invariant: after restricting to a transverse
two-variable germ, the program formally solves the nondegenerate critical
direction and reads the order of the resulting one-variable critical germ. It
does not use the dimension of the global affine Jacobian quotient, which can
also count unrelated critical points. Zero-dimensional rational components and component intersections
also show a local type summary: local multiplicity, tangent cone, and quadratic
rank when the tangent cone is quadratic. Labels such as `ordinary double point
(node/A1)` are only emitted when the quadratic tangent cone has full rank;
otherwise the output keeps the more explicit tangent-cone description. When a
component appears in more than three charts, the interactive analysis suppresses
the explicit per-chart polynomials and lists only the chart numbers.

The decomposition uses Buchberger's algorithm with an S-pair product criterion,
plus recursive splitting of monomial products, powers, common monomial factors,
and univariate rational-root factors over `Q`. Affine-linear component
multiplicities and intersections use direct rational row reduction where
possible. Initial projective surfaces in `P^3` also use a bounded projective
linear-subspace search before affine chart decomposition; each candidate line or
point is verified by exact substitution into the partial derivatives. Rational
point options are still detected by exact evaluation over a bounded rational
grid for the legacy `blowup option N` command. Use `set-bound N` to increase the
search height.

For large exact decompositions, an installed
[Singular](https://www.singular.uni-kl.de/) executable can be used as an opt-in
minimal-prime backend:

```sh
BLOWUP_SINGULAR_BACKEND=1 cargo run --release
```

The backend calls `minAss` over `Q` and parses the exact prime generators back
into the Rust model. If it is disabled or unavailable, the in-process
decomposition remains the fallback.
