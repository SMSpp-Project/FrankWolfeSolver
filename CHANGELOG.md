# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `intFWDirection` `eDirPrimalDual`, the primal-dual bundle Frank-Wolfe
  method: the master problem of the current gradient and of the aggregate
  minimizes a quadratic upper bound on the gap the direction gives, built
  out of the vertex of the aggregate and of how fast its linear function
  grows away from it on the feasible region of each sub-Block, a Euclidean
  ball or a box, with a linear cost of its own; a step that does not pass
  the test of `dblFWAccept` is not taken, nothing is evaluated anew, and the
  direction refused becomes the aggregate, so that the bound never
  decreases; only with `AlgVanilla`

- the per-iteration log (`intLogVerb` 2) ends with the lower bound of the
  iteration and the seconds since the start of `compute()`

- the makefiles build the C++20 module `smspp.frankwolfesolver` when
  `SMSppMOD` is not empty (e.g., `make SMSppMOD=1`), with clang++ 16 or
  later [see `lib/makefile-mod` of SMS++]

- `intFWCutPolicy` `eCutRelevant`, with `intFWCutPool` and `intFWCutActive`:
  the master problem of `eDirBundleMP` is given, out of a pool of pieces,
  the gradient of the current iterate, the pieces with the largest weights
  in the previous master and those with the smallest linearization error
  at the current iterate, the pool keeping the pieces the master is given
  most often; the default `eCutOldest` keeps the newest pieces, as before

- the per-iteration log (`intLogVerb` 2) also says what the step is taken
  along (the gradient, the direction of the master, the same corrected, or
  the gradient after a refusal), how many calls to the oracle have been made
  so far, how many pieces of the master have a positive weight, the step,
  its norm and its cosine with the previous one in the father `Variable`;
  the final summary also says how many calls to the oracle have been made

- the C++20 module `smspp.frankwolfesolver` (`src/FrankWolfeSolver.cppm`),
  which exports every name that the headers declare at namespace scope and
  re-exports `smspp`, so that a source file can write `import
  smspp.frankwolfesolver;` in place of the `#include` of the headers; it is
  built with the CMake option `SMSpp_BUILD_MODULE`, as the target
  `SMS++::FrankWolfeSolver_module`

### Changed

- the pipeline builds the module also with its makefiles, after CMake, and
  runs on a change of the makefiles too [see `ci/makefile` of the umbrella]

- `SMSpp_BUILD_MODULE` no longer stops the configuration when the C++20
  module cannot be built, i.e., with a CMake older than 3.28, a generator
  other than Ninja or Visual Studio, a compiler other than Clang 16 (with
  `clang-scan-deps`) or MSVC 19.34, or a module it imports that is not
  there: the option is turned off with a warning that says why, and the
  build goes on with the headers [see `cmake/CheckModule.cmake` of SMS++]

- the values in the log have 12 significant digits, the precision of the
  stream being given back at the end of `compute()`: with the default 6 two
  runs that differ by less than one part in a million looked the same

### Fixed

- the makefiles of the tests and of the tools build also with g++: they
  take the option that stops at the first error and the one of the
  debugging information, which clang++ and g++ spell differently
  (`-ferror-limit=1` and `-glldb` against `-fmax-errors=1` and `-ggdb`),
  from `SMSppSWERR` and `SMSppSWDBG` of `lib/makefile-OS` of SMS++

- the pieces of the master problem of `eDirBundleMP` were kept from a
  `compute()` to the next, where a `Modification` may have made them no
  longer linearizations of the linking function, and the direction and the
  bound they gave not valid; they are now dropped at every `compute()`, as
  those of `eDirBundle` and `eDirAggregate` already were

## [0.3.0] - 2026-10-09

### Added

- `intFWDirection` `eDirAggregate`: the master problem of two pieces solved
  in closed form, the older piece being the aggregate of all the previous
  ones (the combination the master of the previous iteration has chosen)
  instead of the gradient of the previous iterate; it needs `dblFWAccept`,
  however small, since without it a direction that does not decrease the
  function makes the method stall at the same iterate

- `intFWOnReject`: a bundle direction refused by `dblFWAccept` that still
  decreases the function can be replaced by the combination of it and of the
  gradient whose step is known to pass the test, at the cost of the one more
  call to the oracle that going back to the gradient costs

- `intFWBestLB`: the test of `dblFWAccept` and the stopping test can use the
  best lower bound on the optimum found so far in the `compute()`, which is
  what `get_lb()` then reports

- `dblFWAccept`, the test a step towards the vertex of a bundle direction
  (`eDirBundle`, `eDirBundleMP`) has to pass: the decrease the gradient
  promises along it has to be at least that fraction of the bound the
  direction gives, or else the oracle is asked again with the gradient in
  the same iteration, which gives the bundle directions the `O(1/k)` rate of
  the classical method; the default 0 makes no test, and the final log
  reports how many directions were taken and how many refused

- `IntegralityBarrierFunction`, a `C05Function` built out of a `Block` of a
  (mixed-)integer linear program: the sum over the integer variables of an
  integrality penalty divided by the product of the normalised slacks of
  the rows where each one appears, raised to exponents that weigh the rows,
  with its gradient in the variables and in the exponents, and two forms of
  the penalty, with or without a square root (`set_psi()`, `intPsi` of the
  `IntegralityBarrierSolver`)

- `IntegralityBarrierSolver`, a heuristic for a feasible solution of the
  integer program of its `Block`, which minimizes that function by a
  `FrankWolfeSolver` whose oracle is a `Solver` of the linear relaxation, stops
  as soon as a rounding is feasible and changes the exponents between one run
  and the next. What Frank-Wolfe does is said by the `ComputeConfig` of
  `strFWCfg`, i.e., where the runs start from (`intInitPoint`, the current
  values of the Variable being `eInitBlock`), which Solver of the copy is the
  oracle (`intLMOSlvr`) and how often the rounding is looked at (`intEverykIt`,
  0 looking at it only at the end of each run); the oracle is registered to the
  copy by the `BlockSolverConfig` of the file `strLMOBSCfg`, `intInitSlvr` says
  the Solver of the Block whose solution is written into the Variable before
  Frank-Wolfe starts (-1, the default, none), and the time of each run is what
  is left of `dblMaxTime` of the class, if smaller than that of `strFWCfg`

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

- the pieces of the bundle directions kept from a previous `compute()` are
  dropped at the start of the next one: a Modification of the father may
  have made them no longer linearizations below the linking function, and
  the direction and the bound built out of them would not have been valid

- `compute()` that throws before it starts, on parameters that do not go
  together or on an oracle that is not there, leaves the `Solver` and the
  father unlocked and still listening to the sub-`Block`: the parameters are
  checked before anything is locked, and what may throw after the locks
  unlocks before the exception goes out

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

[Unreleased]: https://gitlab.com/smspp/frankwolfesolver/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/frankwolfesolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/frankwolfesolver/-/compare/0.1.0...0.2.0
