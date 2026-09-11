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
 * See FrankWolfeSolver/frank-wolfe-design.md for the full design.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
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

  intLastParFWSlv  ///< first allowed parameter value for derived classes
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

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
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
