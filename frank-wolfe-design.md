# FrankWolfeSolver — Design document

Porting into SMS++ of the Julia framework [FrankWolfe.jl](https://github.com/ZIB-IOL/FrankWolfe.jl) for
Frank-Wolfe-type algorithms (conditional gradient). The goal is to implement
the FrankWolfe.jl algorithms as *verbatim* as possible, adapted to the
structures of SMS++ (Block tree, `C05Function`, `Solver`, `Solution`,
`Modification`).

Text in English, technical identifiers in English, as per the conventions of the
primer (`SMS++-CONTEXT.md`, §13).

Status: **design approved on the main points; v1 = vanilla Frank-Wolfe.**

---

## 1. What it is and where it attaches

`FrankWolfeSolver` is a `CDASolver` that registers to a "father" Block with
this structure:

- a `FRealObjective` whose `Function` is a `C05Function` (the *linking*
  function, of which we use the *diagonal linearization* = gradient);
- an arbitrary number of **sub-Blocks**, which contain *all* the variables and
  constraints (the father has neither variables nor constraints of its own);
- each sub-Block has a `FRealObjective` containing a `LinearFunction`,
  `DQuadFunction` or `QuadFunction`.

To each sub-Block `FrankWolfeSolver` registers a `:Solver` (the **LMO** of that
Block). The names of the solvers come from the `Configuration` parameters, exactly
as in `LagrangianDualSolver` (LDS): if a sub-Block already has a `:Solver`
registered it is used as is.

The overall feasible region is the **product** of the regions of the children:

```
C = conv(C_1) x conv(C_2) x ... x conv(C_k)
```

where `C_j` is the feasible region of sub-Block `j` (even with integrality
constraints: `FrankWolfeSolver` does not look at them; by combining integer
solutions as a convex combination it effectively minimizes over the convex
hull).

### Analogy with LagrangianDualSolver

| Aspect | LagrangianDualSolver | FrankWolfeSolver |
|---|---|---|
| Linking in the father | linear **constraints** | **objective** `C05Function` |
| What it does with the linking | dualizes (multipliers) | linearizes (gradient) |
| Solver for the children | one per sub-Block | one per sub-Block (= LMO) |
| "Verbatim" direction to copy | solver registration, `Configuration` cache, per-child BSC parameters, `Modification` handling | same |

Most of the *plumbing* (registration/deregistration of the children's solvers,
`Configuration` cache, `str`/`vint` parameters for the per-child
`BlockSolverConfig`, snapshot/restore of the children's objective,
`Modification` handling) is copied as verbatim as possible from LDS.

---

## 2. The composite objective (unification of the three semantics)

The three possible semantics for how the children's objective enters the problem
are special cases of a single scheme of **generalized Frank-Wolfe / composite
conditional gradient**, in which the oracle (= the child's solver) handles
*exactly* a separable part of the objective.

Total objective minimized:

```
F(x) = f_father(x) + Σ_j h_j(x_j),     h_j(y) = α·⟨c_j, y⟩ + β·q_j(y)
```

where `c_j`, `q_j` are the **original** linear and quadratic parts of the
objective of child `j`, and `g = ∇f_father(x)`, restricted to child `j` as `g_j`.
The smooth part `f_father` is linearized; the `h_j` remain exact in the oracle.
The three options are the values of `(α, β)`, selected by the integer parameter
`intLMOObj`:

| `intLMOObj` | name | `(α,β)` | the oracle of child `j` solves | LMO type |
|---|---|---|---|---|
| `0` | `LMOLinear` (default) | `(0,0)` | `min ⟨g_j, v_j⟩` | LP/MILP |
| `1` | `LMOQuad` | `(0,1)` | `min ⟨g_j, v_j⟩ + q_j(v_j)` | QP/MIQP |
| `2` | `LMOFull` | `(1,1)` | `min ⟨c_j + g_j, v_j⟩ + q_j(v_j)` | QP/MIQP |

(The fourth combination `(1,0)` is not required; it is added only if needed.)

**Note (Block-tree semantics)**: the total objective of a Block tree in SMS++ is
the *sum* of all the objectives of the tree; it is exactly what a `:MILPSolver`
on the father minimizes (father objective + Σ children objectives, recursively).
Therefore **`LMOFull` is the mode that coincides with the "true" objective of the
Block tree** and with what the monolithic `:MILPSolver` computes. `LMOLinear`
(children ignored) coincides with the monolithic one only if the children have a
null objective; `LMOQuad` only if they have a purely quadratic objective (null
linear part).

### Mechanics (snapshot + scatter + restore)

In the style of `LagBFunction::cleanup_inner_objective` of LDS:

- in `set_Block`: **snapshot** of the original objective of each child (`c_j`, and
  the quadratic part `q_j` if present);
- each iteration: I write the linear coefficients of the child's objective to
  `α·c_j + g_j` via `LinearFunction::modify_coefficients` /
  `DQuadFunction::modify_linear_coefficients`. The quadratic part is never
  touched: it remains if `β=1`, and is zeroed out once in `set_Block` if `β=0`;
- in `set_Block(nullptr)`/destructor: **restore** the original objective, so the
  Block stays clean.

### The gap is computed in the same way in all modes

The **generalized Frank-Wolfe gap** collapses into a single form:

```
gap(x) = Σ_j [ M_j(x_j) − M_j(v_j) ]   ≥   F(x) − F*
```

where `M_j` is the child's objective *as it has been modified*
(`⟨α·c_j+g_j, ·⟩ + β·q_j`): `M_j(v_j)` is simply **the optimal value that the
child's solver returns**, and `M_j(x_j)` is the same objective evaluated at the
current iterate `x`. In `LMOLinear` mode it reduces exactly to
`⟨g,x⟩ − ⟨g,v⟩`. A single routine for the stopping criterion, independent of the
mode.

### Subtlety for Away-step / Blended Pairwise (post-v1)

In the `LMOQuad`/`LMOFull` modes the oracle is a QP: the *away* direction and the
`active_set_argminmax` must be defined with respect to the smooth part alone
`∇f_father`, with `h_j` entering only in the oracle and in the line search. The
clean classic case is `LMOLinear`. The exact line search remains in closed form if
*both* the father *and* the `q_j` are quadratic (the total Hessian is the sum).

### 2.bis The two problems: value at the iterate vs convex combination (`intCvxComb`)

When the children have a **non-linear** objective (`h_j` with `β=1`, `LMOFull`
mode), the same FW scheme can compute two different quantities — *two different
problems* — depending on how the `h_j` part enters the **value** and the **gap**.
It is a **qualifying feature** of the solver, controlled by `intCvxComb`:

- **(P1) — value at the iterate** (`intCvxComb = eObjAtX = 0`):
  ```
  F_P1(x) = f_father(x) + Σ_j h_j(x_j)
  ```
  one minimizes the *true* composite objective over the convex hull `conv(C)`:
  `h_j` is evaluated **at the current iterate** `x_j` (which is an interior point
  of `conv(C_j)`, not a vertex). It is the bound one would obtain by optimizing
  `F` exactly over the convex closure of the product region.

- **(P2) — convex combination / Dantzig-Wolfe** (`intCvxComb = eObjCvxComb = 1`,
  **default**):
  ```
  F_P2(x) = f_father(x) + Σ_j ( Σ_k λ_{jk} h_j(v_{jk}) )
  ```
  `h_j` is the **convex combination of the costs at the vertices** `v_{jk}` of the
  active set (the `λ_{jk}` are the FW weights), i.e. the **convex hull** of `h_j`
  restricted to the generated vertices. It is exactly the **Dantzig-Wolfe
  decomposition / simplicial decomposition** bound, and — when `conv(C_j)` is the
  integer convex hull of the child — it coincides with the **perspective/Perspective-
  Cut** (P/C) bound that the cut-separated continuous-relaxed DP formulation
  computes.

**Relation.** By Jensen's inequality (`h_j` convex, `x_j = Σ_k λ_{jk} v_{jk}`):
`F_P2(x) ≥ F_P1(x)`, with **equality for linear children** (`q_j=0`): in that case
`h_j` is linear, the convex combination of the vertex costs *is* `h_j(x_j)`, and
the two problems coincide. The difference is purely in the non-linear term of the
children.

**What changes in the algorithm** (everything else is identical — **same LMO**,
same vertices, same direction):
- **bookkeeping of the cost `f_ci`**: each atom stores `c_i = Σ_j h_j(atom_j)`,
  the child-cost *at the vertex*; the child-term of the current value is then
  `cbar = Σ_i λ_i c_i` (P2) instead of `Σ_j h_j(x_j)` (P1). `cost_x = ⟨g..⟩-piece +
  (cvx ? cbar : Σ h_j(x_j))`.
- **line search**: in P2 the child-term is **linear in γ** along `x+γ(v−x)`
  (it interpolates `c_i` between the two atoms), in P1 it is **quadratic** (it is
  `h_j` evaluated at a point that moves). Therefore the exact line search and the
  value along the segment differ, **and consequently the generated iterates
  differ**.
- **gap / stopping criterion**: they use `cost_x` consistent with the mode; in P2
  the gate of the exact line search is relaxed (the model is linear-in-the-weights).

In summary: **same oracle, same scheme, two different value-functions** → P2 ≥
P1, and P2 is the formal bridge between Frank-Wolfe and Dantzig-Wolfe /
Perspective-Cut. This is what allowed validating `FrankWolfeSolver` (P2/`LMOFull`
mode) against `MILPSolver`+DPForm+P/C on `ThermalUnitBlock`: same bound up to
tolerance, even in the case of fractional commitment in which `P1 ≠ P2`.

---

## 3. Gradient extraction and scatter onto the children

The father `C05Function` has as its *active variables* the `ColVariable`s
physically residing in the sub-Blocks. In `set_Block` I build a map

```
gradient position p  →  (sub-Block j, coefficient index i in B_j's objective)
```

by iterating over the active variables of the father and locating their owning
Block (`ColVariable::get_Block()` / climbing up to the direct sub-Block of the
father) and the position in the child's linear objective. This map is static as
long as the structure does not change (cf. §8 Modification).

Each iteration (in the **main thread**, serially, before the fan-out):

1. I write the iterate `x` into the variables of the Block (`Solution::write`);
2. on the father: `compute()` → `compute_new_linearization(true)` →
   `get_linearization_coefficients(buf)` (we assume the father never produces
   *vertical* linearizations);
3. *scatter* of `buf` into the linear coefficients of the children via the map and
   `modify_coefficients`, applying the `(α,β)` mode;
4. (parallel) I run the LMOs, collecting vertices `v_j` and values `M_j(v_j)`.

`min`/`max` handling: the *sense* of the father must coincide with that of the
children (check in `set_Block`, as in LDS). `FrankWolfeSolver` minimizes
internally; for `eMax` the gradient is negated / the sign is handled consistently
with the sense of the children's solvers.

---

## 4. The Frank-Wolfe scheme (v1: vanilla)

The iterate `x` is kept as a `Solution` of the father Block, obtained via
`Block::get_Solution()` (the concrete type depends on the Block: we do not need
to know it). Atoms and active set come with Away-step/BPCG (post-v1); for the
vanilla version the current iterate is enough.

```
x ← initial point (from a first LMO with gradient at a starting point)
for t = 0, 1, 2, ... :
    write x into the Block; compute() father; g ← ∇f_father(x)      # §3
    scatter g onto the children (mode α,β)                          # §2,§3
    parallel for j in 1..k:  v_j, M_j(v_j) ← LMO_j.compute()        # §5
    gap ← Σ_j [ M_j(x_j) − M_j(v_j) ]                               # §2
    if gap ≤ ε  →  STOP (kOK)
    d ← x − v            # v = (v_1,...,v_k) as a Solution of the father Block
    γ ← line_search(...) ∈ [0,1]                                    # §6
    x ← (1−γ)·x + γ·v      # Solution::scale + Solution::sum
```

The combination `(1−γ)x + γv` exploits the **native** support of `Solution` for
(linear/convex) combination — `scale`, `sum` — which automatically recurses into
the sub-Blocks. All `Solution`s have the same interface and here one always and
only combines `Solution`s coming from the **same** Block (the father), hence
mutually compatible whatever the concrete type generated by the Block:
`FrankWolfeSolver` never assumes a concrete type. It is also the "free" basis of
the Away-step (post-v1): the active set will be a set of `Solution`s of the father
with weights `λ_i`.

---

### 4.bis Variants with active set: Away-step and Blended Pairwise (done)

Parameter `intAlgorithm` ∈ {`AlgVanilla`=0, `AlgAwayStep`=1, `AlgBPCG`=2}.
The variants with active set maintain the set of **atoms** (vertices =
`Solution`s of the father) with weights `λ_i`, `x = Σ λ_i a_i`. To choose the
**away vertex** one needs `⟨∇F, a_i⟩` for each atom: the atom is written into the
Block and `Σ_j M_j` (the modified objectives of the children) is evaluated —
`eval_modified_objective()`.

- **Away vertex**: `a = argmax_i ⟨∇F,a_i⟩` (argmin if `eMax`), weight `λ_a`.
- **FW vertex**: `v` from the LMO. `fw_gap = ⟨∇F, x−v⟩`, `away_gap = ⟨∇F, a−x⟩`
  (both ≥0; `fw_gap` certifies optimality → stopping criterion).
- **Away-step**: if `fw_gap ≥ away_gap` → FW step (`d=x−v`, `γ_max=1`); otherwise
  away step (`d=a−x`, `γ_max=λ_a/(1−λ_a)`); update `x ← x − γd`. Drop step
  if `γ=γ_max` (the away atom exits). Weights: FW `λ_i←(1−γ)λ_i`, `λ_v+=γ`; away
  `λ_i←(1+γ)λ_i`, `λ_a−=γ`.
- **BPCG (pairwise)**: `d=a−v`, `γ_max=λ_a`; `x ← x+γ(v−a)`; `λ_a−=γ`, `λ_v+=γ`.
- **Line search**: same exact one, `γ*=clamp(gd/(2Q), 0, γ_max)` with `gd=⟨∇F,d⟩`
  (= the gap of the direction) and `Q=½⟨d,Ad⟩` from `quad_form` (it reuses the
  cached quadratic structure). Agnostic `min(2/(t+2), γ_max)` otherwise.
- **Dedup**: each new vertex `v` is searched for among the atoms (`find_atom`, by
  the values of the father's variables, relative tol); if already present, the
  weight is added instead of creating a duplicate.

**Limited active set (`intMaxAtoms`)**: for large problems, where the active set
can grow and the argmax `O(|aset|)` becomes the bottleneck, a cap is kept.
`intMaxAtoms` (0 = unlimited). When it is exceeded, the **least active** atoms
(counter `f_count` of consecutive iterations with weight > 0 being smaller; in
case of a tie, smaller weight) are **aggregated** into a single atom:
`ā = Σ_{i∈E}(w_i/W) a_i`, weight `W=Σ w_i`. This **preserves `x` exactly**
(analogous to the aggregation of bundle methods). The aggregate is
**indistinguishable** from a vertex (only, it is not extreme): it can itself be an
away-vertex, dropped, or re-aggregated (a recent aggregate may contain previous
ones), and it counts in the budget. **Trade-off**: aggregates are not vertices →
the linear convergence of the away-step is weakened (the "aggregated bundle"
regime), in exchange for limited memory/cost-per-iteration.

**Cached argmax (linear children)**: `⟨∇F(x), a_i⟩ = ⟨g(x), a_i⟩ + c_i`, with
`c_i = ⟨c, a_i⟩` **independent of x** (stored per atom, `Atom::f_ci`, computed at
insertion as `M_j(v)−⟨g,v⟩`, aggregated linearly). Thus the argmax is a **scalar
product** `⟨f_grad, a_i.f_val⟩ + c_i` — no write into the Block nor
`Function::compute` per atom. With quadratic children (`!all_lin`) one goes back to
writing+evaluating. `a_val` (the values of the away atom) is read from
`f_aset[a_idx].f_val` (already stored), no capture.

Status: implemented and validated (away-step + BPCG, dedup, limited active set
with aggregation, cached argmax). 72 base combos + 32 with cap (incl.
`intMaxAtoms=5`) + 32 with caching: `|aset|` at the cap, `x` exact, converges.
TODO perf (postponed): full gram matrix (`active_set_quadratic`) → argmax
`O(|aset|²)` instead of `O(|aset|·G)`, useful for large G if the argmax dominates.

---

## 5. LMO via the children's solvers + parallelism

### 5.1 LMO

**Acquisition of the LMO `:Solver` (v1 choice)**: for each sub-Block the `:Solver`
**already registered** to it is used (registered by the recursive
`BlockSolverConfig` of the tester, `-S`), at index `intLMOSlvr` (default 0) in
`get_registered_solvers()`. The acquisition is **lazy**, at the first `compute()`
(not in `set_Block`), because the registration order of the `:Solver`s to the
various Blocks is not guaranteed when `FrankWolfeSolver` receives `set_Block`. It
is the responsibility of the config that a suitable `:Solver` is registered to
each child (otherwise `compute()` throws an exception). Self-registration from
`FrankWolfeSolver`'s own parameters à la LDS (`str`/`vstr`/`vint` for the
per-child `BlockSolverConfig`) is a v1.1 extension.

Each sub-Block `j` has a `:Solver` registered (the LMO). Running LMO `j`:

1. (already done serially: scatter of the gradient into the child's objective);
2. `solver_j->compute()`;
3. `M_j(v_j) ← solver_j->get_var_value()` (optimal value of the modified objective);
4. the vertex `v_j` stays written in the variables of the sub-Block; it is read
   into a `Solution` (via `Block::get_Solution`) when needed (for `d`, for the
   combination). The `Solution`s of the individual children compose the father's
   `Solution`, which recurses into the sub-Blocks.

### 5.2 Parallelism: persistent bulk-synchronous thread pool

The pattern is bulk-synchronous (barrier after the `k` LMOs), more regular than
the on-demand one of `ParallelBundleSolver`. Proposed scheme (more efficient than
per-iteration `std::async` + active-wait of PBS):

- **Persistent pool** of `min(intMaxThread, k) − 1` worker threads, created lazily
  at the first `compute()`, **reused** for all iterations and all `compute()`s,
  destroyed in `set_Block(nullptr)`/destructor.
- **Atomic counter dispatch**: `std::atomic<Index> next; while( (i =
  next.fetch_add(1)) < k ) run_lmo(i);` → dynamic load balancing (LMOs with
  different times), without static partitioning.
- **The main thread participates** in the work. With `intMaxThread = 1` (default)
  zero extra threads → a purely sequential path (debug, and no cost when not
  needed).
- **Synchronization with `condition_variable` + generation counter** (no
  `sleep`/busy-wait): dispatch barrier (wakes the workers) and join barrier (wakes
  the main on completion).
- **Results in per-index arrays** (`v_value[j]`, `v_status[j]`,
  `v_exception[j]`): no lock; each worker writes only its own slot; each LMO works
  on a distinct Block → no contention on the Block locks.
- **Exceptions**: caught per-slot as `std::exception_ptr` and re-thrown in the
  main after the barrier.

Preconditions (documented): each children's `:Solver` on a distinct Block; the
gradient scatter (which emits `Modification`s) happens in the main *before* the
fan-out; no mutable state shared among workers.

Parameters: `intMaxThread` (default 1).

**Implemented (simple version)**: since for our cases `k` (= number of
sub-Blocks/LMOs) is small (2–3), the cost of creating the threads per-call is
negligible compared to solving the LMOs. So the parallel `run_LMOs` is
load-balanced with an **atomic counter** (`next.fetch_add`) over
`min(intMaxThread, k)` threads (the main participates), no persistent pool. Each
LMO runs on its **distinct** child Block with the identity `f_id` lent to the
`:Solver` (`set_id(f_id)` in `acquire_LMOs`), so `lock(f_id)` on the child already
owned by `f_id` does not contend; no mutable shared state (each LMO writes only in
its own child and in its own `SubBlockData` slot); exceptions caught per-slot
(`SubBlockData::excp`) and re-thrown in the main after the join. Validated: 48
parallel runs (`intMaxThread=4`, all algorithms, with repetitions) identical to
the serial one, no races.

**Persistent pool (postponed)**: the scheme with thread pool +
condition_variable described above is needed only if `k` is large and the LMOs are
tiny (thread creation no longer negligible). To be done if ever needed; likewise
the possible harmonization of `ParallelBundleSolver`.

---

## 6. Line search

From `FrankWolfe.jl/src/linesearch.jl`, for v1:

- **Agnostic** (`2/(t+2)`, i.e. `l/(t+l)`): no extra evaluation, safe default;
  convergence `O(1/t)`.
- **Exact for a quadratic objective**: if `f_father` is `QuadFunction`/
  `DQuadFunction` (plus, in the `LMOQuad`/`LMOFull` modes, the quadratic `q_j`),
  `F` is quadratic along the segment and

  ```
  γ* = clamp( −⟨∇F(x), d⟩ / ⟨d, A d⟩ , 0, γ_max )
  ```

  with `A` the total Hessian (`Quad/DQuadFunction` of the father + Σ `q_j`). For
  `DQuadFunction` `A` is diagonal; for `QuadFunction` the matrix is used (Eigen
  sparse, `get_matrix`).

Post-v1: `Adaptive`/`Backtracking`, `Shortstep`, `Goldenratio`.

`FrankWolfeSolver` checks in `set_Block` whether the father's objective is
`Quad`/`DQuadFunction` to enable the exact path (parameter `intLineSearch` for the
choice; default = auto: exact if quadratic, otherwise Agnostic).

---

## 7. Solver / CDASolver interface

Base class: **`CDASolver`** (as LDS).

- `compute(bool)`: the loop of §4. Lock of the father Block (LDS pattern:
  `is_owned_by` / `lock(f_id)`), `process_outstanding_Modification()`, loop,
  unlock. Return code: `kOK` (gap ≤ ε), `kStopIter`, `kStopTime`,
  `kLowPrecision`, `kInfeasible`/`kUnbounded` (propagated/translated from the
  LMOs).
- `get_var_solution(Configuration*)`: writes the iterate `x` (the maintained
  convex combination) into the variables of the father Block.
- `has_dual_solution()`: `true` **iff** every children's LMO is a `CDASolver` and
  its `has_dual_solution()` is `true`. Otherwise `false`.
- `get_dual_solution(Configuration*)`: **pure forwarding** — for each child,
  `dynamic_cast<CDASolver*>(solver_j)->get_dual_solution(...)`. At termination
  (ideally exact) the optimal duals of the LMOs are valid Lagrange multipliers for
  the overall problem.
- `get_lb()`/`get_ub()`: for minimization, `ub = F(x)` (feasible primal),
  `lb = F(x) − gap` (best seen); at termination `gap → 0` and they coincide. Sign
  inverted for maximization.
- Base events and parameters of `ThinComputeInterface` handled as usual.

---

### 7.bis Children's variables: "leaf" assumption and grandchildren TODO

The gradient→children map and the scatter assume that each *active* variable of
the father's objective belongs to a **direct child** of the father. In general the
variables of a sub-Block "belong" recursively to its ancestors, so a father's
objective could depend on variables defined in a **grandchild** (a child of a
child). v1 **assumes the simple case**: either the children are *leaves* (no
sub-Blocks), or the father's objective depends only on the variables of the direct
children.

**TODO** (not neutral): if the variable is defined in a grandchild, modifying its
coefficient in the objective of that descendant (not of the direct child) might
not be "welcome" to the intermediate child. By analogy with
`LagBFunction::PushCostToOwner`, the possibility is needed of modifying the
coefficient in the Block that *defines* the variable when it is not a direct
child. Postponed.

---

### 7.ter Convergence: guarantee only in the smooth case (father objective)

`FrankWolfeSolver` linearizes the father's objective. **Global convergence to the
minimum is guaranteed only if the father's objective is differentiable** (smooth
with Lipschitz gradient): this is the `[D]QuadFunction` case. If the father's
objective is **nondifferentiable** (e.g. a polyhedral/piecewise-linear function),
Frank-Wolfe with a (sub)gradient **has no guarantee of converging to the
optimum** — it is a known result: F-W is not the subgradient method (it "commits"
to the LMO vertex of *that* subgradient, a potentially ascent direction), and
reducing the step to 0 (DSS) **does not** fix it. Convergent nonsmooth methods
require **smoothing** (Moreau-Yosida / Nesterov), rate $O(1/\sqrt k)$ (tight); see
the papers in `papers/` (White 1993; survey Braun-Pokutta 2022 §3.2.3;
Thekumparampil et al. NeurIPS 2020 — *"FW fails to converge if subgradients are
used instead of gradients"*).

**Decision**: we do NOT implement a dedicated nonsmooth method (no current
application, the fix is non-minor). `FrankWolfeSolver` **makes no direct
reference** to `PolyhedralFunction` or to other nonsmooth types: it treats the
father's objective generically as a `C05Function` and uses its *diagonal
linearization*. It is **explicitly documented** that global convergence is
guaranteed only for the smooth case; with a nonsmooth father objective the
algorithm runs anyway (it may oscillate/stall, Kelley/cutting-plane regime).

Two things hold anyway, for free and correctly:
- the **F-W gap remains a valid bound** even in the nonsmooth case:
  $\text{gap}_t=\langle g_t,x_t-v_t\rangle\ge f(x_t)-f^*$ for every convex
  function and every subgradient $g_t$; therefore `get_lb()`/`get_ub()` are always
  valid — only the bracket **may not close** if the father's objective is
  nonsmooth;
- **design hook** (postponed): an interchangeable "father-gradient provider" + a
  smoothing parameter $\mu$ would allow adding smoothing (for a polyhedral one:
  softmax of the active rows) if ever needed, without the rest of the scheme
  changing.

---

## 8. Modification and warm-start (implemented)

`FrankWolfeSolver` handles **lazily** the `Modification`s coming from the inner
Blocks, keeping the cache structures consistent and **warm-starting** active set /
iterate between one `compute()` and the next. The logic mirrors `LagBFunction`,
but is simpler because `FrankWolfeSolver` is at the top of the chain: it
**consumes** the Modifications, it does not re-translate them into `FunctionMod`
for an external solver.

**Mechanics.**
- No global `inhibit` between two `compute()`s: external Modifications are
  **enqueued** in `v_mod`.
- `process_modifications()`, at the start of `compute()`, drains the queue and
  categorizes it.
- `inhibit_Modification(true)` reused as `play_dumb` **only around the algorithm
  and the restore**: the scatter (and the final restore) change the children's
  objectives and generate "self-inflicted" Modifications to be ignored; the
  inhibit discards only the new ones arriving, **not** the external queue already
  accumulated.
- **Restore-at-end**: at the end of `compute()` the children's objectives return
  to the original `c0` (`restore_objectives`), so an external change between two
  solves acts on a clean state and the re-snapshot of `c0` is correct. At teardown
  the restore uses `eNoMod` (the other `:Solver`s might no longer be there — see
  the note on the teardown crash).

**Categorization** (`guts_of_process_modifications` → bit-mask), which mirrors the
*relax test* of LagBFunction (feasibility is re-checked only when the change can
*restrict* the region):
- **father objective** (bit 1) → re-cache of the quadratic structure
  (`analyze_father`); atoms intact (values and costs do not change).
- **child objective** (bit 2) → re-snapshot `c0` (`snapshot_c0`) + recomputation
  of the atom costs `f_ci` (`recompute_atom_costs`), which are kept.
- **structure** (bit 4) — `C05FunctionModVars` on the objective, `NBModification`,
  variable *fixing* → full re-analysis + drop of the active set.
- **feasible region** (bit 8) — `RowConstraintMod` (LHS/RHS) → check;
  `VariableMod` (sign/integrality), `ConstraintMod` (enforce), `BlockModAD`
  (del-var / add-constraint) → `0` if it *relaxes*, `8`/`4` if it *restricts*;
  generic `BlockMod` → check.
- **Variable fixing**: routed to the **reset** (bit 4), not to 8.
  `ColVariable::is_feasible()` checks sign/integrality but **not** the fixing (a
  *state*, not a constraint), so a feasibility-check would not discard the atoms
  that violate it → the reset is needed.
- *Note* (as in the TODO of LagBFunction): the direction of a `RowConstraintMod`
  (e.g. RHS growing on `≤` → relaxes) is **not** deducible without the previous
  value, so that case is always re-checked.

**Parameter `intHandleMod`** (`eModReset` default / `eModFine`): `eModReset`
always does the worst case (full re-analysis + drop of the active set on any
change, no warm-start); `eModFine` updates only the affected cache and
**preserves the warm-start**.

**Warm-start** (in `compute_active_set` and `compute_vanilla`): if a previous
`compute()` left an active set / iterate still valid (which
`process_modifications` guarantees on objective-only changes), one restarts from
there. Vanilla has no atoms (it cannot reconstruct the exact `cbar`), but it still
reuses the iterate `f_x`, reinitializing `cbar = C_lit(f_x)` (exact for linear
children, "washed away" by the open-loop weights otherwise).

**`feasibility_check`** (eModFine, region change): it writes each atom, calls
`f_Block->is_feasible()`, discards the infeasible ones and rebuilds `f_x` from the
survivors (convex combination of feasible points → feasible); for vanilla it
checks `f_x` directly and keeps it if still feasible, otherwise a cold restart.
**Guard on the aggregates**: the `f_ci` of an aggregated atom is the convex
combination of the costs of the fused (lost) vertices, not `h()` at the point;
with a limited active set (`intMaxAtoms>0`) and quadratic children a
child-objective-change cannot recompute it exactly → in that case the warm-start is
discarded.

**Propagation of infeasibility**: if an LMO returns `kInfeasible` the product
region is empty → `compute()` returns `kInfeasible` (`run_LMOs` sets
`f_lmo_infeas`, checked in the loop and at init).

---

## 9. Module structure and build

Layout that mirrors `LagrangianDualSolver/`:

```
FrankWolfeSolver/
├── include/FrankWolfeSolver.h
├── src/FrankWolfeSolver.cpp
├── obj/
├── makefile        # FWSlvOBJ/INC/H/LIB
├── makefile-s      # + deps (SMS++ excluded)
├── makefile-c      # + SMS++/lib/makefile-c
├── CMakeLists.txt
└── README.md  CHANGELOG.md  LICENSE
```

The testers do **not** live in `FrankWolfeSolver/test`: they live in `tests/`
(see §10), because they depend on the "base" Blocks and on their `:Solver`s,
which are not strict dependencies of `FrankWolfeSolver`.

Factory: `SMSpp_insert_in_factory_h` (header), `SMSpp_insert_in_factory_cpp_0(
FrankWolfeSolver )` (cpp). Makefile macro prefix: `FWSlv` (in the style of
`StcBlk`, `SDDPBk`, `MILP`).

### Parameters (8-point system, §7 of the primer)

The names do not have the `FWSlv` infix: being in the ComputeConfig of a
`FrankWolfeSolver` it is obvious they are its own (unlike LDS, where the names
disambiguate the redirect towards its single inner Solver — a mechanism that does
not apply here).

- `int` (v1, no forwarding): `intLMOObj` (0/1/2), `intLineSearch`
  (auto/agnostic/exact), `intLMOSlvr` (index of the `:Solver` already registered
  to each child to be used as LMO, default 0). Plus the inherited `intMaxThread`
  (default 1), `intMaxIter`.
- `dbl`: inherited `dblRelAcc`/`dblAbsAcc` (tolerance on the gap), `dblMaxTime`.
- v1.1: `intLMOSlvr` → **`vint_LMOSlvr`** (one index per LMO; a vector shorter
  than the number of LMOs → the rest at the default 0; an empty vector [default] →
  all at 0). Plus, à la LDS, `str`/`vstr`/`vint` for the per-child
  `BlockSolverConfig`s and the `Configuration` cache, for the self-registration of
  the LMOs.

---

## 10. Test

**Location**: a single directory `tests/FrankWolfeSolver/` (shared scaffolding in
`fw_test_common.h`, see "Implemented testers"), mirroring
`tests/LagrangianDualSolver_{MMCF,UC,Box}` — from which one copies liberally
(`test.cpp`, `makefile`, config `*.txt`, batch). It lives in `tests/` because it
depends on the base Blocks and on their `:Solver`s, not strict dependencies of
`FrankWolfeSolver`.

Each tester builds a father Block with a simple objective
(`QuadFunction`/`DQuadFunction` or `PolyhedralFunction`) on top of one or more
"base" Blocks of the same type:

- `MCFBlock`, `ThermalUnitBlock` (`UCBlock`), `BinaryKnapsackBlock`.

This way the same Block is solvable by multiple `:Solver`s and they are compared.

**Validity of the cross-check (important)**: F-W minimizes `F` over the **convex
hull** of the children's regions; a monolithic `:MILPSolver` minimizes over the
**true** region. The two optima coincide only when the two sets give the same
value:

| father objective | children | FW vs monolithic MILP | `SolverReading` |
|---|---|---|---|
| linear / `PolyhedralFunction` | even integer | equal (vertices of the hull are integer) | `Exact` (but see caveat below) |
| quadratic (`DQuad`/`Quad`) | **continuous** (e.g. `MCFBlock`) | equal (conv = the region itself) | `Exact` |
| quadratic | generic integer (Knapsack) | FW = relaxation ≤ integer optimum | `LowerBound` (bound, not equality) |
| quadratic | `ThermalUnitBlock` DP form. (P/C) | **equal**: the DP describes the convex hull of the integer solutions | `Exact` |

Order of the tests: **(1) continuous `MCFBlock` + `DQuadFunction` father** → mode
`LMOLinear`, exact line search, cross-check `Exact` (first test, done); (2)
`PolyhedralFunction` father + integer children — **caveat resolved** (see §7.ter):
`F` nondifferentiable ⇒ **no guarantee of global convergence**, F-W runs but may
not reach the optimum; therefore test (2) is not an `Exact` cross-check but an
empirical verification (and `get_lb`/`get_ub` give a valid bracket that may not
close). No smoothing implemented. (3) `ThermalUnitBlock` with DP formulation (P/C)
+ quadratic father → `Exact`, because the DP characterizes the convex hull of the
integer solutions. One step at a time.

**Solver-agnostic, all via configuration**: the C++ code of the tester **makes no
assumption** about which `:Solver`s are attached — neither to the sub-Blocks, nor
to the father Block. The infrastructure of `tests/common_utils.h` is used:

- `process_args` (standard CLI: instance, `-B` BlockConfig, `-S`
  BlockSolverConfig, `-c`/`-p` prefixes); `require_block_config` /
  `require_solver_config`;
- `b_config_Block` applies the `BlockConfig` (-B); `s_config_Block` applies the
  `BlockSolverConfig` (-S) which **registers the `:Solver`s** (to the father and
  to the children) and frees them at the end of the run;
- `SolveAll( block, classify, ref, tol )` runs **all** the `:Solver`s registered
  to the father, prints the uniform line (`print_instance_line`) and does the
  **cross-check** (`cross_check`). For v1: classifier `exact_getter(
  ObjGetter::VarValue )` (both `FrankWolfeSolver` and the comparison solver read as
  the Exact optimum via `get_var_value()`).

The `:Solver`s for the sub-Blocks (the LMOs) and the comparison `:Solver` on the
father are chosen in the config files. The comparison solver on the father
**must** be a `:MILPSolver` (it is the only one that solves the flattened
monolithic problem), but it is the responsibility of whoever writes the config
file that the Solver is correct: otherwise an exception is thrown.
`PolyhedralFunction` as the father's objective is a good non-quadratic case (purely
linear LMO, Agnostic line search).

### Implemented testers

A single directory `tests/FrankWolfeSolver/`, with shared scaffolding in
`fw_test_common.h` (namespace `fwtest`: `build_father`, `make_father_objective`,
`generate_poly`, `collect_vars`, `rnd`/`pos`), and **two** testers that split
along the "independence from the Block type" axis:

- **`test.cpp` (generic, Block-agnostic)** — builds a father with a
  `DQuad`/`Quad`/`Polyhedral` objective on top of `-k` children (of any type, from
  CLI/config) and cross-checks `FrankWolfeSolver` against the comparison solver.
  With **`-M`** it runs random rounds of Modification **of the objective only**
  (father, and possibly of the first child), which are expressed through the
  abstract `FRealObjective` and therefore do **not** require knowing the type of
  the Block. It exercises the *objective-change* branch of `process_modifications`
  (re-snapshot `c0`, recomputation of `f_ci`) and the warm-start, both vanilla and
  active-set, both `eModReset` and `eModFine`.

- **`test-mcf.cpp` (MCF-specific)** — exercises the changes of the **feasible
  region**, which intrinsically require knowing the type: it builds `-k` `MCFBlock`
  under a `DQuad` father, and in random rounds (style of `tests/MCF_MILP`) it
  changes **costs** (`chg_costs` → objective change), **capacities** (`chg_ucaps` →
  region change) or **fixes/unfixes** an arc (closure/reopening = `VariableMod`
  that restricts/relaxes), always with `eModBlck, eModBlck` and clean integer
  capacities ≥ 1 (MCFSimplex is sensitive to many-digit capacities). Cross-check
  against `:MILPSolver` at every round, including the **propagation of
  infeasibility** when the closure of an arc makes the MCF infeasible.

Both run in `regression` (static + `-M` rounds) and in CMake/ctest (`FWS_test`,
`FWS_mcf_test`); the validation matrix (variants × `eModReset`/`eModFine` ×
`DQuad`/`Quad`, and cost/cap/fix × seed × instances) is the one reported in the
tester's README.

### Implementation notes / v1 limits (compute)

The vanilla loop and the active-set variants are implemented in `compute()` and
tested at runtime (see "Implemented testers"). v1 choices/limits:

- **Initialization**: a first LMO at the current point (variables at 0 by
  default) gives the initial vertex `x0 = v0` (always feasible).
- **`f_value = F(x)`** computed as `f_father(x) + (mx_sum − ⟨g,x⟩)` (holds for all
  modes; in `LMOLinear` it reduces to `f_father(x)`). `f_bound = F(x) ∓ gap`.
- **Modification**: full lazy handling with warm-start — see §8. The
  `inhibit_Modification(true)` is no longer global in `set_Block` but reused as
  `play_dumb` only around the algorithm and the restore.
- **`LMOLinear` requires `LinearFunction` children** (no quadratic part to keep in
  the oracle): if a child is DQuad/Quad in `LMOLinear`, `compute()` throws an
  exception (use `LMOQuad`/`LMOFull`). `LMOQuad`/`LMOFull` leave the child's
  quadratic part intact in the oracle.
- **Exact line search** when the father is quadratic (`DQuadFunction` *or*
  `QuadFunction`) and the children are linear (total Hessian = Hessian of the
  father): `γ* = clamp(−gd/(2Q), 0, 1)` with `gd = mv_sum − mx_sum = ⟨∇F,d⟩` and
  `Q = ½⟨d,Ad⟩ = Σ_p a_p d_p² + Σ_(r,c) q_rc d_r d_c` (diagonal + off-diagonal).
  The static quadratic structure of the father (diagonal `a_p` + off-diagonal from
  `mat_nd`) is **cached in `set_Block`**. Non-quadratic father → `Agnostic`.
- **Sequential LMOs** (`run_LMOs`); the persistent thread pool (§5.2) is the next
  step.
- The vertex `v_j` is materialized in the child's variables via
  `lmo->get_var_solution()` after `lmo->compute()`.

---

## 11. Milestone

1. **v1 — vanilla Frank-Wolfe end-to-end** (this doc):
   - module skeleton + factory + base parameters;
   - `set_Block`: structure validation, snapshot of children's objectives,
     gradient map, registration of the children's LMOs (verbatim LDS);
   - gradient + scatter (`LMOLinear` mode first, then `LMOQuad`/`LMOFull`);
   - FW loop + gap + line search (Agnostic + exact quadratic);
   - `get_var_solution` / `get_dual_solution` / `get_lb`-`get_ub`;
   - parallelism (persistent pool; but a first serial version for correctness);
   - tests in `tests/FrankWolfeSolver_<base>` (MCF/Thermal/Knapsack), config-driven,
     cross-check via `SolveAll` against a monolithic `:MILPSolver`.
2. **v2** — active set + Away-step + Blended Pairwise (`afw.jl`,
   `blended_pairwise.jl`, `active_set.jl`, `active_set_quadratic.jl`).
3. **v3** — additional line searches, lazy/cached LMO, fine incremental treatment
   of the `Modification`s.

---

*To be updated as the decisions consolidate.*
