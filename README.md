# FrankWolfeSolver

Definition and implementation of the `FrankWolfeSolver` class, which implements
the `CDASolver` interface within the SMS++ framework for the Frank-Wolfe
(conditional gradient) family of algorithms, ported as closely as possible from
the Julia package [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl).

`FrankWolfeSolver` attaches to a "father" `Block` having the following
structure:

- a `FRealObjective` whose `Function` is a `C05Function` (the "linking"
  function, of which the diagonal linearization, i.e. the gradient, is used at
  each iteration);

- an arbitrary number of sub-`Block`, which contain *all* the `Variable` and
  `Constraint` of the model (the father `Block` has none of its own);

- each sub-`Block` has a `FRealObjective` whose `Function` is a
  `LinearFunction`, a `DQuadFunction` or a `QuadFunction`.

To each sub-`Block` a `:Solver` acts as a Linear Minimization Oracle (LMO),
i.e., it optimizes a linear function over the sub-`Block` feasible region. As
for `LagrangianDualSolver` and its inner `Solver`, these LMO are not compile-time
dependencies: they are obtained at run time through the `Solver` factory and the
`CDASolver` interface. The overall feasible region that `FrankWolfeSolver`
optimizes over is the product of the convex hulls of the sub-`Block` feasible
regions (any integrality in the sub-`Block` being ignored).

A qualifying feature of this solver is that, by the choice of a single
algorithmic parameter (`intCvxComb`) and at essentially no extra oracle cost, it
can solve either of two distinct problems over the same feasible set:

- the genuine *composite* objective `f_father(x) + sum_j h_j(x_j)` over the
  product of the convex hulls of the sub-`Block`;

- the stronger Dantzig-Wolfe / simplicial-decomposition (perspective-cut)
  relaxation of it, in which each sub-`Block` cost is accounted as the convex
  combination of the vertex costs, i.e. it is replaced by its convex envelope
  over the generated vertices.

The two objectives agree on the vertices of the feasible region and coincide
exactly when all sub-`Block` objectives are linear; otherwise the second one is
the convex envelope of the first, hence a stronger (greater) bound. When the
convex hull of a sub-`Block` is that of its integer solutions and its cost is a
convex-quadratic, the second problem is exactly the perspective reformulation /
perspective-cut bound of that sub-`Block`, so the solver can be used as a
decomposition alternative to an explicit Dantzig-Wolfe / perspective-cut
reformulation.


## Getting started

These instructions will let you build `FrankWolfeSolver`.


### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

It's not a build requirement but you will need one or more SMS++ `Solver`
capable of acting as a Linear Minimization Oracle for the sub-`Block`, obtained
at run time through the `Solver` factory.


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
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `FrankWolfeSolver` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or [makefile-s](makefile-s)
including all necessary libraries but not the "core SMS++" one (for the common
case in which this is used together with other modules that already include
them). The makefiles in turn recursively include all the required other
makefiles, hence one should only need to edit the "main makefile" for
compilation type (C++ compiler and its options) and it all should be good to go.
In case some of the external libraries are not at their default location, it
should only be necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

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
