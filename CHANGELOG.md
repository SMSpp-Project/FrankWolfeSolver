# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `IntegralityBarrierFunction`, a `C05Function` built out of a `Block` of a
  (mixed-)integer linear program: the sum over the integer variables of an
  integrality penalty divided by the product of the normalised slacks of
  the rows where each one appears, raised to exponents that weigh the rows,
  with its gradient in the variables and in the exponents, and two forms of
  the penalty, with or without a square root (`set_psi()`, `intPsi` of the
  `IntegralityBarrierSolver`)

- `IntegralityBarrierSolver`, a heuristic for a feasible solution of the
  integer program of its `Block`, which minimizes that function by a
  `FrankWolfeSolver` whose oracle is a `Solver` of the linear relaxation,
  stops as soon as a rounding is feasible and changes the exponents between
  one run and the next

- `intLineSearch` `LSFixed`, the fixed step `dblFWStep`, and `intInitPoint`,
  which starts the vanilla method from the current values of the `Variable`

- the tester registers to the father another Solver that never computes,
  which has to receive nothing from a `compute()` of many iterations, also
  when the father has a default channel of its own and when the oracle of a
  sub-Block throws halfway, after which the costs have to be the original
  ones and a new `compute()` has to find the optimum

- a tester of the module, `test/`, posed on Blocks of the core alone: a father
  `AbstractBlock` with a separable quadratic `DQuadFunction` over three
  sub-`Block` with box constraints and a `BoxSolver` each, whose optimum is
  known in closed form, solved by every algorithm under both bookkeeping
  modes and both directions, with the bounds checked on both sides of the
  optimum, a coordinate whose box is a single point, a `ColVariable` that the
  `Objective` of its sub-`Block` does not price (added there and taken out
  again when the `Solver` leaves), and a change of the costs of a sub-`Block`
  solved again. The pipeline of the module builds the module alone and runs
  this tester, the suites that pose the decomposition on the `Block` of other
  modules running where those modules are

### Changed

- the checks on what the Block holds ask it for its groups of Variable and of
  Constraint, the vectors of `boost::any` they used to ask for not being
  there any more

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- a Variable of the father Objective that the Objective of its sub-Block does
  not price is no longer refused: what says which sub-Block it belongs to is
  the Block it is of, and the Objective of that sub-Block is where it is
  added, with a zero coefficient, so that the Oracle sees the gradient that is
  scattered onto it. Nothing of the problem changes, and the addition is
  undone when the Solver detaches, so that the Block is left as it was found.
  What is refused is only a Variable of no sub-Block at all, which the
  decomposition cannot move

### Fixed

- a Solver registered to the father next to `FrankWolfeSolver` (e.g., the
  `MILPSolver` of a battery, which solves the Block once and never again)
  received the Modification of every iteration, which rewrites the costs of
  all the sub-Block, and piled them up without bound: 315 KB per iteration
  on `goto10_8`, 25 GB in 4 minutes, which killed the runner of the MCFBlock
  batteries. `compute()` now makes a channel of the father its default one
  while it runs: the sub-Block and their Solver, the LMO first, receive the
  Modification as before, while the channel is emptied at each iteration
  and, at the end, discarded if the costs are back to those at the
  beginning, as `restore_objectives()` leaves them, or else shipped with the
  last Modification, which bring the other Solver to the state of the Block.
  The father gets back the default channel it had, on every way out of the
  method, exceptions comprised

- an exception thrown by an LMO on the sequential path (`intMaxThread` <= 1
  or one sub-Block) was caught and never re-thrown, so that the method went
  on with the result of the previous call, and on the parallel path the
  exceptions after the first were kept and re-thrown by the next call: the
  first one is now re-thrown on both paths, and the others are forgotten

- `scatter()` marks the sub-Block Objectives as changed before changing
  them, so that `restore_objectives()` also undoes a `scatter()` interrupted
  by an exception

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

## [0.2.0] - 2026-09-12

### Changed

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

[Unreleased]: https://gitlab.com/smspp/frankwolfesolver/-/compare/0.2.0...develop
[0.2.0]: https://gitlab.com/smspp/frankwolfesolver/-/compare/0.1.0...0.2.0
