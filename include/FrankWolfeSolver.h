/*--------------------------------------------------------------------------*/
/*----------------------- File FrankWolfeSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the FrankWolfeSolver class, a CDASolver implementing
 * Frank-Wolfe (conditional gradient) type algorithms, ported as closely as
 * possible from the Julia package FrankWolfe.jl.
 *
 * FrankWolfeSolver attaches to a "father" Block having:
 *
 * - a FRealObjective whose Function is a C05Function (the "linking" function,
 *   of which the diagonal linearization, i.e. the gradient, is used);
 *
 * - an arbitrary number of sub-Block, which contain *all* the Variable and
 *   Constraint (the father Block has none of its own);
 *
 * - each sub-Block has a FRealObjective whose Function is a LinearFunction,
 *   a DQuadFunction or a QuadFunction.
 *
 * To each sub-Block a :Solver acts as a Linear Minimization Oracle (LMO).
 * The overall feasible region is the product of the convex hulls of the
 * sub-Block feasible regions (integrality in the sub-Block is ignored).
 *
 * By a single parameter (intCvxComb) the solver can solve, at no extra oracle
 * cost, either the genuine composite objective over that feasible region or the
 * stronger Dantzig-Wolfe / perspective-cut relaxation of it; see the GENERAL
 * NOTES of the class and the intCvxComb parameter.
 *
 * The direction is the one the by-the-book method uses, i.e., the gradient of
 * the linking function at the current iterate is what the LMO is given, which
 * reduces the function to its first-order model and ignores every other piece
 * of information about it. Formulae that instead pass the solution of a small
 * stabilized master problem, built out of the first-order information already
 * at hand and having the by-the-book choice as the special case in which one
 * linearization alone is kept, are those of
 *
 *  A. Frangioni, F. Rinaldi "Bundle-inspired Direction Formulae for
 *  Conditional Gradient Methods", draft
 *
 * Those formulae, in the notation of this solver, are the following. Write
 * \f$ \bar{x}^k \f$ for the current iterate, \f$ \bar{g}^k =
 * \nabla f( \bar{x}^k ) \f$ for the gradient of the linking function there,
 * \f$ \bar{x}^{k-1} \f$ and \f$ \bar{g}^{k-1} \f$ for the same at the
 * previous iterate, and \f$ x^k \f$ and \f$ g^k \f$ for the vertex the
 * oracle has returned and the gradient there. The translated function
 * \f[ h^k( d ) = f( \bar{x}^k + d ) - f( \bar{x}^k ) \f]
 * has \f$ h^k( 0 ) = 0 \f$ and, by convexity, is bounded below by each of
 * the three linear functions
 * \f[ \bar{l}^k( d ) = \bar{g}^k d \quad , \quad
 *     \bar{l}^{k-1}( d ) = \bar{g}^{k-1} d - \bar{\alpha}^{k-1} \quad ,
 *     \quad l^k( d ) = g^k d - \alpha^k \f]
 * whose constants are the linearization errors at the current iterate,
 * \f[ \bar{\alpha}^{k-1} = f( \bar{x}^k ) - f( \bar{x}^{k-1} )
 *      - \bar{g}^{k-1} ( \bar{x}^k - \bar{x}^{k-1} ) \quad , \quad
 *     \alpha^k = f( \bar{x}^k ) - f( x^k ) - g^k ( \bar{x}^k - x^k ) \f]
 * both non-negative, the first one being 0 for the gradient taken at the
 * iterate itself. The by-the-book method keeps \f$ \bar{l}^k \f$ alone; the
 * piecewise-linear model
 * \f$ \check{h}^k( d ) = \max \{ \bar{l}^k( d ) , \bar{l}^{k-1}( d ) ,
 * l^k( d ) \} \leq h^k( d ) \f$ uses all of them, and a stabilizing term
 * makes it into the master problem
 * \f[ \min \{ v + \frac{1}{2t} \| d \|^2 \; : \;
 *     v \geq \bar{g}^k d \; , \;
 *     v \geq \bar{g}^{k-1} d - \bar{\alpha}^{k-1} \; , \;
 *     v \geq g^k d - \alpha^k \} \f]
 * a convex quadratic program in the direction \f$ d \f$ and one further
 * variable. Its dual is the problem of finding the convex multipliers
 * \f$ \theta \f$ that minimize
 * \f[ \bar{\alpha}^{k-1} \bar{\theta}^{k-1} + \alpha^k \theta^k
 *     + \frac{t}{2} \| \bar{g}^{k-1} \bar{\theta}^{k-1} + g^k \theta^k
 *     + \bar{g}^k \bar{\theta}^k \|^2 \f]
 * and the two are tied by \f$ d_*^k = - t z_*^k \f$, with
 * \f[ z_*^k = \bar{g}^{k-1} \bar{\theta}^{k-1} + g^k \theta^k
 *             + \bar{g}^k \bar{\theta}^k \quad , \quad
 *     \alpha_*^k = \bar{\alpha}^{k-1} \bar{\theta}^{k-1}
 *                   + \alpha^k \theta^k \f]
 * so that \f$ z_*^k \in \partial_{\alpha_*^k} f( \bar{x}^k ) \f$: what
 * the oracle is given is an approximate subgradient of the linking function
 * at the current iterate, and a positive multiple of a direction is the same
 * direction for it. The by-the-book choice is the feasible special case
 * \f$ \bar{\theta}^k = 1 \f$, i.e., \f$ z_*^k = \bar{g}^k \f$.
 *
 * The stabilization decides between the two ends: the larger \f$ t \f$ is,
 * the more the norm weighs and the more the direction is the least-norm
 * combination of the gradients at hand, whatever their linearization errors;
 * the smaller it is, the more those errors matter, and as \f$ t \to 0 \f$
 * the gradient is recovered. Since \f$ z_*^k \f$ is an approximate
 * subgradient of the same kind as the pieces it combines, it can take their
 * place in the model of the next iteration without changing the solution of
 * the master, which is what allows keeping a bundle of any size and
 * compressing it to one pair.
 *
 * Of this, intFWDirection == eDirBundle implements the two pieces
 * \f$ ( \bar{g}^k , 0 ) \f$ and
 * \f$ ( \bar{g}^{k-1} , \bar{\alpha}^{k-1} ) \f$, for which the master
 * is solved in closed form: with \f$ d = \bar{g}^{k-1} - \bar{g}^k \f$,
 * \f[ \bar{\theta}^{k-1} = \min \{ 1 , \max \{ 0 ,
 *      - ( t \, \bar{g}^k d + \bar{\alpha}^{k-1} ) /
 *        ( t \| d \|^2 ) \} \} \f]
 * and \f$ z_*^k = \bar{g}^k + \bar{\theta}^{k-1} d \f$. Neither piece
 * costs an evaluation: both are information the method has already paid for.
 *
 * See FrankWolfeSolver/frank-wolfe-design.md for the full design.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __FrankWolfeSolver
 #define __FrankWolfeSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CDASolver.h"

#include "Block.h"

#include "FRealObjective.h"

#include "C05Function.h"

#include "Solution.h"

#include <exception>

#include "MasterProblemBlock.h"

#include <tuple>

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------- FORWARD DECLARATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

 class LinearFunction;   ///< the linear objective case
 class DQuadFunction;    ///< the (diagonal) quadratic objective case (also Quad)
 class QuadFunction;     ///< the full (off-diagonal) quadratic objective case

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS FrankWolfeSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a CDASolver implementing Frank-Wolfe type algorithms
/** The FrankWolfeSolver class derives from CDASolver and implements the
 * Frank-Wolfe (conditional gradient) family of algorithms, ported from the
 * Julia package FrankWolfe.jl. v1 implements the "vanilla" Frank-Wolfe
 * algorithm; Away-step and Blended Pairwise come later.
 *
 * The solver minimizes the "composite" objective
 *
 *     F(x) = f_father(x) + sum_j h_j(x_j) ,
 *            h_j(y) = alpha * <c_j,y> + beta * q_j(y)
 *
 * over X = prod_j conv(X_j), the product of the convex hulls of the sub-Block
 * feasible regions, where f_father is the father Block C05Function (linearized
 * each iteration) and h_j is the part of sub-Block j's objective kept exactly
 * in its oracle (c_j, q_j being its original linear and quadratic terms). The
 * pair (alpha,beta) is selected by intLMOObj; see that parameter.
 *
 * TWO PROBLEMS AT NO EXTRA COST (intCvxComb)
 *
 * A qualifying feature of this solver is that, by the choice of a single
 * algorithmic parameter (intCvxComb) and at essentially no extra cost, it can
 * solve either of two distinct problems over the same feasible set X:
 *
 *   (P1)  min_{x in X}  f_father(x) + sum_j h_j(x_j)
 *
 *   (P2)  min_{lambda}  f_father( sum_k lambda_k v_k )
 *                       + sum_j sum_k lambda_k^j h_j(v_k^j)
 *
 * where the v_k are vertices of X (the oracle outputs) and lambda is a convex
 * combination. (P1) evaluates each sub-Block cost h_j at the (generally
 * fractional) iterate x_j; (P2) accounts it as the convex combination of the
 * vertex costs, i.e. it replaces h_j by its convex envelope over the vertices.
 *
 * Let f_1, f_2 be the two objectives as functions of x in X: they agree on the
 * vertices of X (where the convex combination is trivial) and, by convexity of
 * h_j (Jensen), f_2(x) >= f_1(x) inside X, with f_2 the convex envelope of f_1.
 * Hence val(P2) >= val(P1): (P2) is the *stronger* relaxation. They coincide
 * iff all h_j are linear.
 *
 * (P2) is precisely the Dantzig-Wolfe / simplicial-decomposition bound: FW
 * generates the columns (sub-Block vertices) and the master recombines them
 * with their true costs. When conv(X_j) is the convex hull of the integer
 * solutions of a sub-Block and h_j is a convex-quadratic cost, (P2) is exactly
 * the value of the perspective reformulation / perspective-cut (P/C) bound of
 * that sub-Block -- so the solver can be used as a decomposition alternative to
 * an explicit DW / P-C reformulation.
 *
 * What makes it free: the oracle (LMO) is *identical* in the two modes -- it
 * always returns vertices accounting for the exact h_j. Only the
 * value/gap/line-search bookkeeping differs. In (P2) the sub-Block term is the
 * convex combination of vertex costs, which is *linear* in the step gamma (in
 * the active-set variants it is just sum_i lambda_i c_i, a quantity already
 * cached per atom), so the exact line search reduces to the father-only
 * quadratic case; in (P1) the sub-Block term is h_j(x_j(gamma)), nonlinear in
 * gamma when h_j is nonlinear. Because the *iterates* are driven by the line
 * search, the two modes follow different trajectories and converge to the two
 * different optima -- it is not a mere change of the reported value.
 *
 * intCvxComb selects the mode: eObjCvxComb (the default) solves (P2),
 * eObjAtX solves (P1). See that parameter for the (v1) line-search caveat. */

class FrankWolfeSolver : public CDASolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 using Index = Block::Index;  ///< import Index from Block

 /// how the sub-Block objectives enter the problem (value of int_FWSlv_LMOObj)
 enum lmo_obj_type {
  LMOLinear = 0 ,  ///< (alpha,beta)=(0,0): oracle min <g_j,v_j>; pure LP/MILP
  LMOQuad   = 1 ,  ///< (alpha,beta)=(0,1): oracle min <g_j,v_j> + q_j(v_j)
  LMOFull   = 2    ///< (alpha,beta)=(1,1): oracle min <c_j+g_j,v_j> + q_j(v_j)
  };

 /// which line search to use (value of int_FWSlv_LineSearch)
 enum line_search_type {
  LSAuto     = 0 , ///< exact if the father objective is quadratic, else agnostic
  LSAgnostic = 1 , ///< the open-loop 2/(t+2) rule
  LSExact    = 2   ///< exact line search (requires a quadratic total objective)
  };

 /// which direction the oracle is given (value of intFWDirection)
 enum direction_type {
  eDirGradient = 0 ,  ///< the gradient at the current iterate
  eDirBundle   = 1 ,  ///< the master of two pieces, solved in closed form
  eDirBundleMP = 2    ///< the master of a bundle, a MasterProblemBlock
  };

/*--------------------------------------------------------------------------*/
 /// which Frank-Wolfe variant to run (value of intAlgorithm)
 enum algo_type {
  AlgVanilla  = 0 , ///< vanilla Frank-Wolfe (no active set)
  AlgAwayStep = 1 , ///< Away-step Frank-Wolfe (with active set)
  AlgBPCG     = 2   ///< Blended Pairwise Conditional Gradient (with active set)
  };

 /// how the sub-Block objective is accounted in the value (value of intCvxComb)
 /** Selects which of two distinct (but closely related) problems the solver
  * actually solves; see intCvxComb and the GENERAL NOTES for the full
  * discussion. The two coincide when the sub-Block objectives are linear. */
 enum cvx_comb_type {
  eObjAtX     = 0 , ///< sub-Block cost evaluated at the (fractional) iterate x_j
  eObjCvxComb = 1   ///< sub-Block cost as the convex combination of vertex costs
  };

 /// how Modification from the sub-Block are handled (value of intHandleMod)
 enum handle_mod_type {
  eModReset = 0 , ///< any sub-Block Modification triggers a full re-analysis of
                  ///< the cached structure (simple, always correct)
  eModFine  = 1   ///< categorize each Modification and update only the affected
                  ///< cached information (cheaper, keeps more across re-solves)
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC PARAMETERS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public parameters of FrankWolfeSolver
 *  @{ */

 /// public enum "extending" int_par_type_CDAS to FrankWolfeSolver
 /** Parameter names drop the redundant "FWSlv" infix: being inside the
  * ComputeConfig of a FrankWolfeSolver, it is obvious that they are its own
  * (unlike LagrangianDualSolver, whose names disambiguate the redirect to its
  * single inner Solver, a mechanism that does not apply here). */
 enum int_par_type_FWSlv {
  intLMOObj = intLastParCDAS ,
  ///< how the sub-Block objectives enter the problem
  /**< Selects how the sub-Block FRealObjective enters the problem; one of the
   * lmo_obj_type values: LMOLinear (default), LMOQuad, LMOFull. */

  intLineSearch ,
  ///< which line search to use
  /**< One of the line_search_type values: LSAuto (default), LSAgnostic,
   * LSExact. */

  intLMOSlvr ,
  ///< index of the registered :Solver of each sub-Block to use as its LMO
  /**< For each sub-Block, FrankWolfeSolver uses, as its Linear Minimization
   * Oracle, the :Solver registered to it at this position (in
   * get_registered_solvers()). Default 0, i.e. the first registered :Solver.
   * It is the caller's / configuration's responsibility that an appropriate
   * :Solver is registered to each sub-Block (else compute() throws).
   *
   * v1.1 TODO: turn this into a vint_LMOSlvr parameter (one index per LMO; a
   * vector shorter than the number of LMO leaves the rest at the default 0; an
   * empty vector, the default, leaves all of them at 0). */

  intAlgorithm ,
  ///< which Frank-Wolfe variant to run
  /**< One of the algo_type values: AlgVanilla (default), AlgAwayStep, AlgBPCG. */

  intMaxAtoms ,
  ///< maximum size of the active set (0 = unbounded)
  /**< Only used by the active-set variants (AlgAwayStep, AlgBPCG). When the
   * active set would exceed this size, the least-active atoms (smallest
   * consecutive-active count, ties broken by smallest weight) are merged into a
   * single *aggregate* atom (a convex combination preserving the iterate x).
   * 0 (the default) means no bound. */

  intCvxComb ,
  ///< how the sub-Block objective is accounted: which problem to solve
  /**< Selects how the (kept-exact part of the) sub-Block objective h_j enters
   * the value/gap/line-search, hence which of two distinct problems the solver
   * solves over the same feasible set X = prod_j conv(X_j); one of the
   * cvx_comb_type values:
   *
   * - eObjAtX: each h_j is taken at the current (in general fractional) iterate
   *   x_j. The solver minimizes the *genuine* composite objective
   *       F(x) = f_father(x) + sum_j h_j(x_j)
   *   over the convex hull X.
   *
   * - eObjCvxComb (default): each h_j is replaced by its *convex envelope* over
   *   the sub-Block vertices, i.e. the sub-Block cost is accounted as the convex
   *   combination sum_k lambda_k^j h_j(v_k^j) of the oracle-returned vertices
   *   (the lambda being the very weights with which the iterate is built). This
   *   is the Dantzig-Wolfe / "disaggregated" relaxation: when h_j is a convex
   *   function but conv(X_j) is the convex hull of integer points, this is
   *   exactly the value that an explicit Dantzig-Wolfe reformulation yields, and
   *   for a separable convex-quadratic sub-Block cost it coincides with the
   *   perspective reformulation / perspective-cut (P/C) bound.
   *
   * The two coincide when all sub-Block objectives are linear; they differ only
   * for nonlinear (e.g. quadratic) h_j, where convexity (Jensen) gives
   *   sum_k lambda_k h_j(v_k) >= h_j( sum_k lambda_k v_k ) ,
   * so eObjCvxComb >= eObjAtX, i.e. the DW value is the *stronger* bound.
   *
   * The remarkable point is that the *oracle is identical* in the two modes (it
   * always returns vertices accounting for the exact sub-Block cost): only the
   * value/gap/line-search bookkeeping differs, and the convex-combination value
   * is linear in the step gamma (in the active-set variants it is just
   * sum_i lambda_i c_i, already cached per atom), so it is obtained at no extra
   * oracle cost. This is the simplicial-decomposition / Dantzig-Wolfe reading of
   * Frank-Wolfe: FW generates the columns (sub-Block vertices) and the master
   * recombines them with their true costs. See the GENERAL NOTES.
   *
   * Implementation note (v1): the exact line search of eObjAtX is currently
   * available only when the sub-Block objectives are linear (where eObjAtX and
   * eObjCvxComb coincide); for nonlinear h_j the eObjAtX mode falls back to the
   * agnostic 2/(t+2) rule, so eObjCvxComb (the default) is the mode to use for
   * nonlinear sub-Block objectives. An exact at-iterate line search for the
   * nonlinear case is future work. */

  intHandleMod ,
  ///< how to handle Modification coming from the sub-Block
  /**< Selects how the solver reacts to a Modification coming from a sub-Block
   * (or the father Objective), one of the handle_mod_type values. The solver
   * processes the queued Modification lazily, at the beginning of each
   * compute(); between two compute() the sub-Block objectives are left in their
   * original state (the per-iteration scatter is undone at the end of
   * compute()), so that any external change to them is "clean".
   *
   * - eModReset (default): any such Modification triggers a full re-analysis of
   *   the cached structure (the sub-Block objective snapshots, the
   *   gradient-to-sub-Block scatter map, the father quadratic cache). Simple
   *   and always correct.
   *
   * - eModFine: each Modification is categorized and only the affected cached
   *   information is rebuilt --- a father-Objective change re-caches the father
   *   quadratic structure only; a sub-Block-Objective change re-snapshots its
   *   linear coefficients only; a change of the *variables* (or an
   *   NBModification) still triggers a full re-analysis; a change to the
   *   sub-Block feasible region needs no action here (the active set is rebuilt
   *   from scratch each compute()).
   *
   * Note: the finer handling becomes more valuable once the active set is
   * warm-started across compute() (so that atoms persist and may need a
   * feasibility re-check); v1 cold-starts each compute(), so the two modes
   * differ only in how much of the (cheap) structural cache is rebuilt. */

  intFWDirection ,
  ///< which direction the solver gives the oracle
  /**< Selects what is passed to the LMO in place of the objective, one of the
   * direction_type values:
   *
   * - eDirGradient (default): the gradient of the linking function at the
   *   current iterate, i.e. the by-the-book method, which reduces the
   *   function to its first-order model at that point;
   *
   * - eDirBundle: the solution of a stabilized master problem built out of
   *   the first-order information already at hand, i.e. the current gradient
   *   (whose linearization error at the current iterate is 0) together with
   *   the one of the previous iterate, carried forward with its error. The
   *   direction is then a convex combination of the two, hence an approximate
   *   subgradient of the linking function at the current iterate, and the
   *   weight of the older one grows with dblFWt. eDirGradient is the special
   *   case in which that weight is 0, which is what dblFWt = 0 gives.
   *
   * The information the second mode uses is the one the method has already
   * paid for: no evaluation is added [see the file documentation for where
   * these formulae come from]. */

  intFWBundleSize ,
  ///< how many pieces the master problem of eDirBundleMP keeps
  /**< The size of the bundle of pairs (gradient, linearization error) that
   * the master problem of intFWDirection == eDirBundleMP is built on. With 2
   * it is the master that eDirBundle solves in closed form, which is how the
   * two are checked against each other; the larger it is, the more of the
   * information the method has produced the direction is drawn from. Has no
   * effect under the other directions. Default 10. */

  intLastParFWSlv  ///< first allowed parameter value for derived classes
  };

/*--------------------------------------------------------------------------*/
 /// public enum "extending" dbl_par_type_CDAS to FrankWolfeSolver

 enum dbl_par_type_FWSlv {
  dblFWt = dblLastParCDAS ,
  ///< how much the master problem of eDirBundle is stabilized
  /**< The weight of the squared norm in the master problem that gives the
   * direction when intFWDirection is eDirBundle. The larger it is, the more
   * the direction is the least-norm combination of the gradients at hand,
   * whatever their linearization errors; the smaller it is, the more those
   * errors matter, and with 0 the direction is the gradient at the current
   * iterate, i.e. the by-the-book method. It has no effect under
   * eDirGradient. Default 1. */

  dblLastParFWSlv  ///< first allowed parameter value for derived classes
  };

/*--------------------------------------------------------------------------*/
 /// public enum "extending" str_par_type_CDAS to FrankWolfeSolver

 enum str_par_type_FWSlv {
  strFWMPBCfg = strLastParCDAS ,
  ///< the BlockSolverConfig of the Solver of the master of eDirBundleMP
  /**< The name of the file describing the BlockSolverConfig of the Solver
   * that solves the master problem of intFWDirection == eDirBundleMP, which
   * is a quadratic program. It is required by that direction, the master
   * having no Solver of its own to fall back on, and compute() throws if it
   * is not given; it has no effect under the other directions. */

  strLastParFWSlv  ///< first allowed parameter value for derived classes
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 *  @{ */

 /// constructor: initialises the parameters to their default values
 FrankWolfeSolver( void ) : CDASolver()
 {
  v_events.resize( max_event_number() );  // the three standard event slots
  set_default_parameters();
  }

/*--------------------------------------------------------------------------*/

 /// destructor: detaches from the Block and releases all resources
 ~FrankWolfeSolver() override { set_Block( nullptr ); }

/*--------------------------------------------------------------------------*/

 /// FrankWolfeSolver supports the three standard ThinComputeInterface events
 /** eBeforeTermination (vetoable optimality stop), eEverykIteration (every
  * intEverykIt iterations) and eEveryTTime (every dblEveryTTm seconds). */
 EventID max_event_number( void ) const override { return( 3 ); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) father Block
 /** Validates the structure of the father Block (no own Variable/Constraint,
  * FRealObjective with a C05Function, at least one sub-Block each with a
  * Linear/DQuad/QuadFunction objective of consistent sense), snapshots the
  * sub-Block objectives, builds the gradient-to-sub-Block scatter map and
  * acquires the per-sub-Block LMO :Solver. set_Block( nullptr ) restores the
  * original sub-Block objectives and releases everything. */

 void set_Block( Block * block ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------------ PARAMETER HANDLING ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of FrankWolfeSolver
 *  @{ */

 using ThinComputeInterface::set_par;  // restore the hidden overloaded methods

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastParFWSlv ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE Block ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the problem encoded in the Block
 *  @{ */

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
 *  @{ */

 [[nodiscard]] OFValue get_var_value( void ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] OFValue get_lb( void ) override;

 [[nodiscard]] OFValue get_ub( void ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool has_var_solution( void ) override;

 void get_var_solution( Configuration * solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// tells whether a dual solution is available
 /** True iff every sub-Block LMO :Solver is a CDASolver and each of them has
  * a dual solution available; the dual solution of FrankWolfeSolver is the
  * collection of the LMO dual solutions, which at (exact) termination are
  * valid Lagrangian multipliers for the overall problem. */

 [[nodiscard]] bool has_dual_solution( void ) override;

 /// write the dual solution: pure forwarding to each LMO :Solver
 void get_dual_solution( Configuration * solc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED TYPES -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// per-sub-Block bookkeeping
 struct SubBlockData {
  Block * block = nullptr;          ///< the sub-Block
  FRealObjective * obj = nullptr;   ///< its FRealObjective
  Function * fun = nullptr;         ///< its objective Function (Lin/DQuad/Quad)
  LinearFunction * lin = nullptr;   ///< == fun if a LinearFunction, else nullptr
  DQuadFunction * dq = nullptr;     ///< == fun if a [D]QuadFunction, else nullptr
  Solver * lmo = nullptr;           ///< the :Solver used as LMO

  /// snapshot of the original linear coefficients c_j (objective order)
  std::vector< Function::FunctionValue > c0;

  /// the Variable this Solver added to the Objective of the sub-Block
  /** A Variable of the father Objective that the Objective of its sub-Block
   * does not price is added to it with a zero coefficient, so that the Oracle
   * of that sub-Block sees the gradient this Solver scatters onto it; nothing
   * of the problem changes, and the addition is undone when this Solver
   * detaches [see cleanup()]. */

  std::vector< Variable * > added;

  /// objective-order indices of the variables touched by the father gradient
  std::vector< Block::Index > obj_idx;
  /// father-gradient positions matching obj_idx (parallel vector)
  std::vector< Block::Index > grad_idx;

  /// result slot: optimal value M_j(v_j) of the last LMO call
  Function::FunctionValue value = 0;
  Function::FunctionValue mx = 0;   ///< M_j(x_j) at the current iterate x
  int status = 0;                   ///< return code of the last LMO call
  std::exception_ptr excp;          ///< exception thrown by the LMO (parallel)
  };

/*--------------------------------------------------------------------------*/

 /// an active-set atom: a father Solution with its convex weight
 /** An atom is either a vertex of the product polytope (an LMO solution) or an
  * *aggregate* (a convex combination of evicted atoms, used to keep the active
  * set bounded); the two are indistinguishable except that an aggregate is not
  * an extreme point. f_val caches the father active-variable values (for
  * deduplication and aggregation); f_count is the number of consecutive
  * iterations the atom has been active (weight > 0), used to pick which atoms
  * to evict when the active set exceeds intMaxAtoms. */

 struct Atom {
  Solution * f_sol = nullptr;                       ///< the father Solution
  double f_weight = 0;                              ///< the convex weight
  std::vector< Function::FunctionValue > f_val;     ///< father active-var values
  Function::FunctionValue f_ci = 0;                 ///< the sub-Block cost at the
            ///< atom, sum_j h_j(atom_j) (the x-independent part of
            ///< <grad F, atom>); for an aggregate, the convex combination of the
            ///< merged atoms' costs, so sum_i lambda_i f_ci is the eObjCvxComb
            ///< (Dantzig-Wolfe) sub-Block cost of the iterate (see intCvxComb)
  int f_count = 0;                                  ///< consecutive-active count
  };

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// set every parameter to its default value
 void set_default_parameters( void );

 /// release everything acquired in set_Block and restore the sub-Block objectives
 void cleanup( void );

 /// (re)validate the father Objective and (re)cache its quadratic structure
 void analyze_father( void );

 /// (re)scan the sub-Block: validate, snapshot c0, build the scatter map
 void analyze_subBlocks( void );

 /// re-read the original linear coefficients c0 of every sub-Block Objective
 /// (assumes the variable structure is unchanged)
 void snapshot_c0( void );

 /// recompute the cached sub-Block cost f_ci of every active-set atom after a
 /// sub-Block-Objective change (used to keep the active set warm-started)
 void recompute_atom_costs( void );

 /// drop the active-set atoms that have become infeasible (after a change of a
 /// sub-Block feasible region), rebuilding the iterate from the survivors;
 /// resets the warm start if no atom survives
 void feasibility_check( void );

 /// restore the original sub-Block objective coefficients (c0); quiet == true
 /// issues no Modification (for teardown), else the change is propagated
 void restore_objectives( bool quiet );

 /// process the Modification queued from the sub-Block (lazily, at compute());
 /// categorizes them and rebuilds the affected cached information
 void process_modifications( void );

 /// categorize a single Modification: returns a bit-mask, 0 = harmless,
 /// 1 = father Objective changed, 2 = a sub-Block Objective changed,
 /// 4 = structural change (variables/NBModification) -> full re-analysis
 char guts_of_process_modifications( const Modification * mod ) const;

 /// acquire, for each sub-Block, the LMO :Solver at index intLMOSlvr
 void acquire_LMOs( void );

 /// the direction the oracle is given, out of the information at hand
 /** With intFWDirection == eDirBundle, replaces the direction the oracle is
  * given with the solution of the stabilized master problem built out of the
  * gradient at the current iterate, whose linearization error there is 0, and
  * the one of the previous iterate carried forward with its own error; with
  * two of them the master is solved in closed form. It then records the
  * current pair as the previous one for the next iteration. Takes the value
  * of the linking function at the current iterate, and expects f_xval to hold
  * that iterate. Does nothing under eDirGradient. */

 void bundle_direction( OFValue fx );

/*--------------------------------------------------------------------------*/
 /// the same direction, taken from a MasterProblemBlock [see eDirBundleMP]
 /** Keeps a bundle of the pairs (gradient, linearization error) the method
  * has produced, in a MasterProblemBlock with one component and a proximal
  * stabilization, and takes its aggregated subgradient as the direction.
  * With a bundle of two it is the closed form of bundle_direction(), which
  * is how the two are checked against each other. */

 void master_direction( OFValue fx );

/*--------------------------------------------------------------------------*/
 /// builds the MasterProblemBlock of eDirBundleMP, once

 void build_master( void );

/*--------------------------------------------------------------------------*/
 /// evaluate the father Objective at the current point and fill f_grad
 void evaluate_gradient( void );

 /// scatter the father gradient into the sub-Block linear objectives (mode a,b)
 void scatter( void );

 /// run all the sub-Block LMO and collect v_j and M_j(v_j) into v_sb
 void run_LMOs( bool changedvars );

 /// vanilla Frank-Wolfe loop (no active set)
 int compute_vanilla( bool changedvars );

 /// Away-step / Blended-Pairwise loop (with active set)
 int compute_active_set( bool changedvars );

 /// run all handlers of the given event type and return the decisive action
 /** Invokes, in order, all the registered handlers of event type @p type
  * (an index into v_events) until one returns something other than
  * eContinue; returns that response (eForceContinue / eStopOK / eStopError),
  * or eContinue if there are no handlers or they all return eContinue. */
 int run_event( int type );

 /// sum over the sub-Block of their (modified) objective at the current point
 OFValue eval_modified_objective( void );

 /// read the father active-variable values at the current point into dst
 void capture_father_values( std::vector< Function::FunctionValue > & dst );

 /// 1/2 <d, A d> for the quadratic father, with d_p = b_p - a_p (exact LS)
 OFValue quad_form( const std::vector< Function::FunctionValue > & a ,
                    const std::vector< Function::FunctionValue > & b );

 /// release the active set (deleting the atom Solution)
 void clear_active_set( void );

 /// index of the active-set atom matching the given father values, or
 /// f_aset.size() if none (for deduplicating a new vertex)
 Index find_atom( const std::vector< Function::FunctionValue > & val ) const;

 /// if the active set exceeds intMaxAtoms, merge the least-active atoms into a
 /// single aggregate atom (a convex combination that preserves the iterate x)
 void bound_active_set( void );

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - -

 int f_lmo_obj;        ///< int_FWSlv_LMOObj
 int f_line_search;    ///< int_FWSlv_LineSearch
 int f_lmo_slvr;       ///< int_FWSlv_LMOSlvr
 int f_algorithm;      ///< intAlgorithm
 int f_max_atoms;      ///< intMaxAtoms
 int f_cvx_comb;       ///< intCvxComb (eObjAtX / eObjCvxComb)
 int f_handle_mod;     ///< intHandleMod (eModReset / eModFine)
 int f_direction;      ///< intFWDirection (eDirGradient / eDirBundle / MP)
 double f_t;           ///< dblFWt, the stabilization of that master problem
 int f_bundle_size;    ///< intFWBundleSize, how many pieces the master keeps

 MasterProblemBlock * f_mpb = nullptr;
 ///< the master problem of eDirBundleMP, built once and kept

 std::string f_mpb_cfg;
 ///< strFWMPBCfg, the BlockSolverConfig of the Solver of that master

 int f_next_slot = 0;
 ///< which slot of the bundle the next piece takes when it is full
 int f_max_thread;     ///< intMaxThread
 int f_max_iter;       ///< intMaxIter
 double f_max_time;    ///< dblMaxTime
 double f_rel_acc;     ///< dblRelAcc
 double f_abs_acc;     ///< dblAbsAcc
 int f_log_verb;       ///< intLogVerb (verbosity of the log)
 int f_everyk;         ///< intEverykIt (period of the eEverykIteration events)
 double f_every_t;     ///< dblEveryTTm (period of the eEveryTTime events)

 // statistics of the last compute(), for the final-summary log (intLogVerb 1)
 Index f_niter = 0;        ///< iterations performed by the last compute()
 OFValue f_last_gap = 0;   ///< final Frank-Wolfe gap of the last compute()

 // problem structure - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool f_max = false;           ///< true if the (father) Objective is eMax
 Index f_nsb = 0;              ///< number of sub-Block
 FRealObjective * f_obj = nullptr;  ///< the father FRealObjective
 C05Function * f_fun = nullptr;     ///< the father C05Function
 DQuadFunction * f_dq_father = nullptr;  ///< == f_fun if a [D]QuadFunction
                                         ///< (enables the exact line search)
 QuadFunction * f_quad_father = nullptr; ///< == f_fun if a QuadFunction (off-diag)

 // cached (static) quadratic structure of the father, for the exact line
 // search: diagonal a_p and the off-diagonal terms ( r , c , q ), r > c
 std::vector< Function::FunctionValue > f_father_diag;
 std::vector< std::tuple< Index , Index , Function::FunctionValue > >
                                                            f_father_offdiag;

 std::vector< SubBlockData > v_sb;  ///< per-sub-Block bookkeeping

 std::vector< Function::FunctionValue > f_grad;  ///< father gradient buffer

 std::vector< Function::FunctionValue > f_bdir;
 ///< the direction of eDirBundle, when it is not the gradient itself

 const std::vector< Function::FunctionValue > * p_dir = nullptr;
 ///< what scatter() gives the oracle: f_grad, or f_bdir under eDirBundle

 std::vector< Function::FunctionValue > f_pgrad;  ///< gradient of the previous
 std::vector< Function::FunctionValue > f_pxval;  ///< iterate it was taken at
 OFValue f_pval = 0;                              ///< value of the function there
 std::vector< Function::FunctionValue > f_xval;  ///< father active-var values at x
 std::vector< Function::FunctionValue > f_vval;  ///< father active-var values at v

 // algorithmic state - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Solution * f_x = nullptr;     ///< the current iterate x (a father Solution)

 /// the active set: atoms with convex weights, x = sum_i w_i * atom_i;
 /// used by Away-step and BPCG (bounded by intMaxAtoms via aggregation)
 std::vector< Atom > f_aset;

 OFValue f_value = 0;          ///< F(x), best primal value found
 OFValue f_bound = 0;          ///< F(x) - gap, best bound found
 bool f_has_sol = false;       ///< whether a (primal) solution is available
 bool f_modified = false;      ///< whether the sub-Block objectives were modified
 bool f_lmo_infeas = false;    ///< whether the last run_LMOs found an infeasible
                               ///< sub-Block (=> the father is infeasible)

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class FrankWolfeSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* FrankWolfeSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File FrankWolfeSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
