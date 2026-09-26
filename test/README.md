# test

A tester for `FrankWolfeSolver` that needs nothing but the core SMS++ library.

It builds in memory a father `AbstractBlock` whose `Objective` is a separable
quadratic `DQuadFunction` over the `ColVariable` of its three sub-`Block`, each
of them an `AbstractBlock` with a `BoxConstraint` on every `ColVariable`, a
linear `Objective` and a `BoxSolver` as its oracle. The optimum of the whole
problem is known in closed form (the unconstrained minimum of each coordinate
projected on its box), some coordinates of it lying in the interior of their
box and some on either of its faces.

Every algorithm (vanilla, Away-step, Blended Pairwise) is run under both
bookkeeping modes (`intCvxComb`) and with the direction taken from the gradient
and from the two-piece bundle (`intFWDirection`); each run has to find the
optimal value, and the lower and upper bounds it reports have to hold on the
two sides of it. The last run then changes the costs of a sub-`Block` and
solves again.

The exit code is 0 when every check passes, printing `All tests passed!!`, and
1 otherwise. The `makefile` builds the executable including the
`FrankWolfeSolver` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
