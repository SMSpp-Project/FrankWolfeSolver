# FrankWolfeSolver

Definition and implementation of the `FrankWolfeSolver` class, which implements
the `CDASolver` interface within the SMS++ framework with a family of
Frank-Wolfe (conditional gradient) type algorithms, ported as closely as
possible from the Julia package
[FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl).

`FrankWolfeSolver` attaches to a "father" `Block` (B) with the following
structure:

- no `Variable` and no `Constraint` of its own: all the `Variable` (and
  `Constraint`) belong to the sub-`Block`;

- (B) has a `FRealObjective` whose `Function` is a `C05Function` (the "linking"
  function), of which the diagonal linearization, i.e., the gradient, is used;

- an arbitrary number of sub-`Block`, each having a `FRealObjective` whose
  `Function` is a `LinearFunction`, a `DQuadFunction` or a `QuadFunction`.

To each sub-`Block` a `:Solver` is registered that acts as a Linear Minimization
Oracle (LMO): each Frank-Wolfe iteration computes the gradient of the father
objective, scatters it into the sub-`Block` objectives, and asks each LMO to
minimize the (modified) sub-`Block` objective over its feasible region. The
overall feasible region is therefore the product of the *convex hulls* of the
sub-`Block` feasible regions, so any integrality in a sub-`Block` is relaxed (the
LMO is free to return integer vertices: `FrankWolfeSolver` minimizes over their
convex hull). The solver is a `CDASolver` because Frank-Wolfe naturally produces,
"for free", a dual (Lagrangian) solution: `get_dual_solution()` simply forwards
to the dual solutions of the sub-`Block` LMO, which at (exact) termination are
valid Lagrangian multipliers for (B).

The solver minimizes the "composite" objective

        F(x) = f_father(x) + sum_j h_j(x_j)

over `X = prod_j conv(X_j)`, where `f_father` is the father objective
(linearized each iteration) and `h_j` is the part of sub-`Block` j's objective
kept exactly in its oracle. Three algorithmic variants are available
(`intAlgorithm`): "vanilla" Frank-Wolfe, Away-step Frank-Wolfe and Blended
Pairwise Conditional Gradient (the last two maintain an *active set* of atoms,
optionally bounded by aggregation via `intMaxAtoms`); the line search
(`intLineSearch`) is the exact one when the father objective is quadratic, the
open-loop `2/(t+2)` rule otherwise.

## Two problems at no extra cost

A qualifying feature of `FrankWolfeSolver` is that, by the choice of a single
algorithmic parameter (`intCvxComb`) and at essentially no extra cost, it solves
either of two distinct problems over the same feasible region `X`:

- **(P1)** `min { f_father(x) + sum_j h_j(x_j) : x in X }` — the *genuine*
  composite objective, with each sub-`Block` cost evaluated at the (generally
  fractional) iterate `x_j`;

- **(P2)** the *Dantzig-Wolfe / simplicial-decomposition* relaxation, where each
  `h_j` is accounted as the convex combination of the costs of the oracle
  vertices, i.e., it is replaced by its convex envelope over those vertices.

The two coincide when the sub-`Block` objectives are linear; they differ for
nonlinear (e.g. quadratic) `h_j`, where convexity gives `(P2) >= (P1)`, so (P2)
is the *stronger* relaxation. When `conv(X_j)` is the convex hull of the integer
solutions of a sub-`Block`, (P2) is exactly the Dantzig-Wolfe bound over that
hull — so `FrankWolfeSolver` can be used as a decomposition alternative to an
explicit Dantzig-Wolfe reformulation. The oracle is identical in the two
modes; only the value / gap / line-search bookkeeping differs, and the
convex-combination value is linear in the step, hence obtained for free. See the
documentation of `intCvxComb` (and the GENERAL NOTES of the class) for the full
discussion, and [frank-wolfe-design.md](frank-wolfe-design.md) for the design.


## Getting started

These instructions will let you build `FrankWolfeSolver`.


### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

It's not a build requirement but you will need a SMS++ `:Solver` to register to
each sub-`Block` as its LMO (e.g. an `MCFSolver` for an `MCFBlock`, or a generic
`:MILPSolver`).


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(FrankWolfeSolver)
target_link_libraries(<my_target> SMS++::FrankWolfeSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Each executable using `FrankWolfeSolver` has to include a "main
makefile" of the module, which typically is either [makefile-c](makefile-c)
including all necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). The makefiles in turn recursively include
all the required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all should
be good to go.

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/frankwolfesolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa


### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
