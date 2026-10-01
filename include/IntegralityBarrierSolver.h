/*--------------------------------------------------------------------------*/
/*-------------------- File IntegralityBarrierSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the IntegralityBarrierSolver class, a heuristic looking
 * for a feasible solution of a (mixed-)integer linear program by minimizing
 * an IntegralityBarrierFunction over the polyhedron of its continuous
 * relaxation with the Frank-Wolfe method, whose oracle solves the linear
 * relaxation, changing the exponents of the barrier between one run of
 * Frank-Wolfe and the next.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IntegralityBarrierSolver
 #define __IntegralityBarrierSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"
#include "FRealObjective.h"
#include "IntegralityBarrierFunction.h"
#include "Solution.h"
#include "Solver.h"

#include <memory>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*--------------------- CLASS IntegralityBarrierSolver ---------------------*/
/*--------------------------------------------------------------------------*/
/// a Frank-Wolfe heuristic for feasible solutions of integer programs
/** The IntegralityBarrierSolver looks for a feasible solution of the
 * (mixed-)integer linear program of the Block it is registered to, which
 * holds ColVariable, FRowConstraint with LinearFunction and bounds, such as
 * an AbstractBlock read out of an MPS file. It minimizes the
 * IntegralityBarrierFunction of the Block [see IntegralityBarrierFunction.h]
 * over the polyhedron of its continuous relaxation by a FrankWolfeSolver,
 * whose oracle is a Solver registered to the Block that solves the linear
 * relaxation: the one of name strLMOSolver, with the ComputeConfig of the
 * file strLMOCfg. Since the integer points are the global minima of the
 * function, the iterates are hoped to get to one of them.
 *
 * To do so, at the first compute() an abstract copy of the Block is made
 * [see AbstractBlock::mirror()], the only sub-Block of an AbstractBlock of
 * its own, whose Objective is the IntegralityBarrierFunction of the copy,
 * and to which a FrankWolfeSolver is registered with the ComputeConfig of
 * the file strFWCfg; the oracle is registered to the copy. Nothing touches
 * the Block or the Solver registered to it, and the copy is made again at
 * the compute() after any Modification of the Block. The point found is
 * mapped back to the Block, and kept only if it is feasible there too, the
 * copy being a relaxation of the Block if mirror() could not reproduce all
 * of it. Then:
 *
 * - the starting point x0 is, with intInitSlvr >= 0, the solution of the
 *   Solver registered to the Block at that position, e.g., the linear
 *   relaxation by a barrier method without crossover, which is a point in
 *   the interior of the polyhedron; with -2, the current values of the
 *   Variable; with -1 (the default) there is none, and each run of
 *   Frank-Wolfe starts from the vertex the oracle gives at the current
 *   values of the Variable;
 *
 * - all the exponents y of the rows start at 0, and then, intMaxIter times
 *   at most: Frank-Wolfe starts from x0 with the current y, stopping as soon
 *   as rounding the integer Variable of the vertex the oracle has given or
 *   of the iterate gives a feasible point; if that never happens, the
 *   exponents are changed by the rule of dblYLevel and x0 is moved towards
 *   the last iterate by dblX0Step.
 *
 * The rule changing the exponents takes, for every row whose normalised
 * slack exceeds dblYLevel, the largest step along the projected gradient of
 * the function with respect to y that keeps y in [ 0 , 1 ], and sets to 0
 * the exponent of every active row (slack 0) that has a fractional
 * Variable, whose slack is then no longer looked at.
 *
 * If the Block has sub-Block, the linear terms of their Objective are
 * moved, in the copy, into the Objective of the copy itself, which is the
 * one FrankWolfeSolver gives the oracle the gradient through, so that the
 * oracle minimizes what it is given whatever the tree.
 *
 * compute() returns kInfeasible if the linear relaxation is empty, and kOK
 * if a feasible point has been found, which
 * get_var_solution() writes into the Variable of the Block and whose value
 * of the Objective of the Block get_ub() returns, and kStopIter or
 * kStopTime otherwise; the lower bound is -INF, and dblRelAcc is +INF,
 * since a heuristic promises no accuracy. */

class IntegralityBarrierSolver : public Solver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*----------------------------- PARAMETERS ---------------------------------*/
/*--------------------------------------------------------------------------*/
 /// the int parameters of IntegralityBarrierSolver, on top of those of
 /// Solver; intMaxIter is the number of runs of Frank-Wolfe

 enum int_par_type_IBSlv {
  intInitSlvr = intLastAlgPar ,  ///< the Solver of the starting point
                                 /**< The position, among the Solver
				  * registered to the Block, of the one whose
				  * solution is the starting point; -2 means
				  * the current values of the Variable, and -1
				  * (the default) no starting point but the
				  * vertex of the oracle [see the class]. */
  intLastAlgParIBSlv             ///< first new int parameter of derived
                                 ///< classes
  };

 /// the double parameters of IntegralityBarrierSolver, on top of those of
 /// Solver; dblMaxTime is the time limit of compute()

 enum dbl_par_type_IBSlv {
  dblYLevel = dblLastAlgPar ,    ///< the slack above which y is moved
                                 /**< The exponent of a row is moved along
				  * the gradient only if its normalised slack
				  * exceeds this value [see the class].
				  * Default 0.9. */
  dblX0Step ,                    ///< how far x0 is moved to the iterate
                                 /**< After each run of Frank-Wolfe the
				  * starting point becomes ( 1 - dblX0Step )
				  * x0 + dblX0Step x. Default 0. */
  dblBarrierEps ,                ///< the epsilon of the function
                                 /**< The epsilon of the
				  * IntegralityBarrierFunction. Default
				  * 1e-8. */
  dblLastAlgParIBSlv             ///< first new double parameter of derived
                                 ///< classes
  };

 /// the string parameters of IntegralityBarrierSolver

 enum str_par_type_IBSlv {
  strFWCfg = strLastAlgPar ,     ///< the ComputeConfig of the Frank-Wolfe
                                 /**< The file of the ComputeConfig of the
				  * FrankWolfeSolver; empty (the default)
				  * means its default parameters, save those
				  * the class sets (intInitPoint, intLMOSlvr
				  * and intEverykIt). */
  strLMOSolver ,                 ///< the name of the Solver of the oracle
                                 /**< The name of the Solver that solves the
				  * linear relaxation [see the class], e.g., a
				  * :MILPSolver; empty by default, which
				  * compute() refuses. */
  strLMOCfg ,                    ///< the ComputeConfig of the oracle
                                 /**< The file of the ComputeConfig of the
				  * Solver of strLMOSolver, which has to make
				  * it solve the linear relaxation (e.g.,
				  * intRelaxIntVars 1 for a :MILPSolver);
				  * empty by default. */
  strLastAlgParIBSlv             ///< first new string parameter of derived
                                 ///< classes
  };

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 IntegralityBarrierSolver( void ) : Solver() {}

 ~IntegralityBarrierSolver() override { set_Block( nullptr ); }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
 /// sets the Block whose integer program is to be solved

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR HANDLING THE PARAMETERS ----------------*/
/*--------------------------------------------------------------------------*/

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 void set_par( idx_type par , std::string && value ) override;

 using Solver::set_par;

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastAlgParIBSlv ) );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastAlgParIBSlv ) );
  }

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastAlgParIBSlv ) );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] OFValue get_lb( void ) override {
  return( - Inf< OFValue >() );
  }

 [[nodiscard]] OFValue get_ub( void ) override { return( f_ub ); }

 [[nodiscard]] bool has_var_solution( void ) override {
  return( f_sol != nullptr );
  }

 /// any Modification of the Block has its copy, the father and the function
 /// made again at the next compute()

 void add_Modification( sp_Mod & mod ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// writes the feasible point found into the Variable of the Block

 void get_var_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the number of runs of Frank-Wolfe of the last compute()

 [[nodiscard]] long get_elapsed_iterations( void ) const override {
  return( f_runs );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the IntegralityBarrierFunction, there after the first compute()

 [[nodiscard]] IntegralityBarrierFunction * get_function( void ) {
  return( f_fun );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 /// builds the father, the function and the FrankWolfeSolver

 void build( void );

 /// frees what build() made

 void unbuild( void );

 /// true if rounding the integer Variable of the Block is feasible, in
 /// which case the point is kept as the solution found

 bool rounding_feasible( void );

 /// changes the exponents y at the current values of the Variable

 void update_y( void );

/*--------------------------- PROTECTED FIELDS -----------------------------*/

 int f_init_slvr = -1;          ///< intInitSlvr
 int f_max_iter = 10;           ///< intMaxIter, the runs of Frank-Wolfe
 int f_log_verb = 0;            ///< intLogVerb
 double f_max_time = Inf< double >();  ///< dblMaxTime
 double f_y_level = 0.9;        ///< dblYLevel
 double f_x0_step = 0;          ///< dblX0Step
 double f_eps = 1e-8;           ///< dblBarrierEps
 std::string f_fw_cfg;          ///< strFWCfg
 std::string f_lmo_name;        ///< strLMOSolver
 std::string f_lmo_cfg;         ///< strLMOCfg
 Solver * f_lmo = nullptr;      ///< the Solver of the oracle, if it is ours

 Block * f_father = nullptr;    ///< the AbstractBlock made by build()
 AbstractBlock * f_copy = nullptr;  ///< the abstract copy of the Block
 IntegralityBarrierFunction * f_fun = nullptr;  ///< the function, which
                                                 ///< f_obj owns
 std::unique_ptr< FRealObjective > f_obj;  ///< the Objective of f_father
 Solver * f_fw = nullptr;       ///< the FrankWolfeSolver

 Solution * f_sol = nullptr;    ///< the feasible point found, if any
 Solution * f_sol_copy = nullptr;  ///< the same, in the copy
 OFValue f_ub = Inf< OFValue >();  ///< its value
 int f_runs = 0;                ///< runs of Frank-Wolfe of the last compute()
 bool f_rebuild = false;        ///< the Block has changed

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( IntegralityBarrierSolver ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* IntegralityBarrierSolver.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File IntegralityBarrierSolver.h -----------------*/
/*--------------------------------------------------------------------------*/
