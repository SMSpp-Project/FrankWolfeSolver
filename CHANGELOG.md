# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
