/*--------------------------------------------------------------------------*/
/*------------------- File IntegralityBarrierSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntegralityBarrierSolver class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Francesca Demelas \n
 *         Laboratoire d'Informatique de Paris Nord \n
 *         Universite' Sorbonne Paris Nord \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli, Francesca Demelas
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>

#include "AbstractBlock.h"
#include "LinearFunction.h"
#include "Modification.h"
#include "FrankWolfeSolver.h"
#include "IntegralityBarrierSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( IntegralityBarrierSolver );

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

// the father made by build(): an AbstractBlock that can let its sub-Block
// and its Objective go without deleting them
class FatherBlock : public AbstractBlock
{
 public:
 FatherBlock( void ) : AbstractBlock( nullptr ) {}
 void release( void ) {
  v_Block.clear();
  reset_objective();
  }
 };

const std::array< std::string , 2 > int_names = { "intInitSlvr" , "intPsi" };
const std::array< std::string , 3 > dbl_names = { "dblYLevel" , "dblX0Step" ,
						   "dblBarrierEps" };
const std::array< std::string , 2 > str_names = { "strFWCfg" ,
						   "strLMOBSCfg" };

// the linear terms of the Objective of the Block inside the copy moved into
// the Objective of the copy itself, which is the one FrankWolfeSolver
// writes the gradient into [see FrankWolfeSolver::intLMOObj], those inside
// being left with zero coefficients: the oracle, which sees all of them,
// then minimizes what FrankWolfeSolver gives it and nothing else
void fold_objectives( AbstractBlock * copy )
{
 auto robj = dynamic_cast< FRealObjective * >( copy->get_objective() );
 auto rlf = robj ? dynamic_cast< LinearFunction * >( robj->get_function() )
		 : nullptr;
 if( ( ! rlf ) && ( ! copy->get_nested_Blocks().empty() ) )
  throw( std::invalid_argument( "IntegralityBarrierSolver::build: the "
				"Objective of a Block with sub-Block is not "
				"linear" ) );

 std::vector< Block * > tree( copy->get_nested_Blocks().begin() ,
			      copy->get_nested_Blocks().end() );
 for( std::size_t t = 0 ; t < tree.size() ; ++t ) {
  for( auto sub : tree[ t ]->get_nested_Blocks() )
   tree.push_back( sub );
  auto obj = dynamic_cast< FRealObjective * >( tree[ t ]->get_objective() );
  if( ! obj )
   continue;
  auto lf = dynamic_cast< LinearFunction * >( obj->get_function() );
  if( ! lf )
   throw( std::invalid_argument( "IntegralityBarrierSolver::build: the "
				 "Objective of a sub-Block is not linear" ) );
  const double sign = ( obj->get_sense() == robj->get_sense() ) ? 1 : -1;
  for( Block::Index i = 0 ; i < lf->get_num_active_var() ; ++i ) {
   const auto c = lf->get_coefficient( i );
   if( c == 0 )
    continue;
   auto x = static_cast< ColVariable * >( lf->get_active_var( i ) );
   const auto k = rlf->is_active( x );
   if( k < rlf->get_num_active_var() )
    rlf->modify_coefficient( k , rlf->get_coefficient( k ) + sign * c ,
			     eNoMod );
   else
    rlf->add_variable( x , sign * c , eNoMod );
   lf->modify_coefficient( i , 0 , eNoMod );
   }
  }
 }

/*--------------------------------------------------------------------------*/

// the ComputeConfig of the file fn, if any, given to the Solver s
void apply_cfg( Solver * s , const std::string & fn )
{
 if( fn.empty() )
  return;
 auto c = Configuration::deserialize( fn );
 auto cc = dynamic_cast< ComputeConfig * >( c );
 if( ! cc ) {
  delete c;
  throw( std::invalid_argument( "IntegralityBarrierSolver::build: " + fn +
				" holds no ComputeConfig" ) );
  }
 s->set_ComputeConfig( cc );
 delete cc;
 }

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 unbuild();
 delete f_sol;
 f_sol = nullptr;
 delete f_sol_copy;
 f_sol_copy = nullptr;
 f_ub = Inf< OFValue >();
 f_rebuild = false;
 Solver::set_Block( block );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::build( void )
{
 if( f_father )
  return;

 if( f_lmo_bscfg.empty() )
  throw( std::invalid_argument( "IntegralityBarrierSolver::build: no "
				"BlockSolverConfig of the oracle "
				"[strLMOBSCfg]" ) );

 // the father, whose only sub-Block is the abstract copy of the Block and
 // whose Objective is the function of the copy: nothing of what follows
 // touches the Block or the Solver registered to it
 auto father = new FatherBlock();
 f_copy = new AbstractBlock();
 f_copy->mirror( f_Block );
 fold_objectives( f_copy );
 father->add_nested_Block( f_copy );

 f_fun = new IntegralityBarrierFunction( f_eps );
 f_fun->set_psi( f_psi , eNoMod );
 f_fun->build( f_copy );
 f_fun->set_y( std::vector< Function::FunctionValue >(
					     f_fun->get_rows().size() , 0 ) ,
	       eNoMod );
 f_obj = std::make_unique< FRealObjective >( father , f_fun );
 f_obj->set_sense( Objective::eMin , eNoMod );
 father->set_objective( f_obj.get() , eNoMod );
 f_father = father;

 // the oracle, on the copy; the BlockSolverConfig, once cleared, is what
 // takes away the Solver it has registered [see unbuild()]
 auto c = Configuration::deserialize( f_lmo_bscfg );
 f_lmo_bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! f_lmo_bsc ) {
  delete c;
  throw( std::invalid_argument( "IntegralityBarrierSolver::build: " +
				f_lmo_bscfg + " holds no "
				"BlockSolverConfig" ) );
  }
 f_lmo_bsc->apply( f_copy );
 f_lmo_bsc->clear();

 f_fw = Solver::new_Solver( "FrankWolfeSolver" );
 apply_cfg( f_fw , f_fw_cfg );
 f_fw_max_time = f_fw->get_dbl_par( dblMaxTime );
 f_father->register_Solver( f_fw );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::unbuild( void )
{
 if( ! f_father )
  return;

 if( f_fw ) {
  f_father->unregister_Solver( f_fw , true );
  f_fw = nullptr;
  }
 if( f_lmo_bsc ) {
  f_lmo_bsc->apply( f_copy );
  delete f_lmo_bsc;
  f_lmo_bsc = nullptr;
  }

 // the Objective, which deletes the function, goes while the Variable of
 // the copy it is active in are still there
 static_cast< FatherBlock * >( f_father )->release();
 f_obj.reset();
 f_fun = nullptr;
 delete f_father;
 f_father = nullptr;
 delete f_copy;
 f_copy = nullptr;
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

bool IntegralityBarrierSolver::rounding_feasible( void )
{
 // the integer Variable rounded, the others as they are; the values are put
 // back if the point is not feasible
 const auto n = f_fun->get_num_active_var();
 std::vector< double > old( n );
 for( Block::Index i = 0 ; i < n ; ++i ) {
  auto x = static_cast< ColVariable * >( f_fun->get_active_var( i ) );
  old[ i ] = x->get_value();
  if( x->is_integer() )
   x->set_value( std::round( old[ i ] ) );
  }

 // the point is kept as a Solution of the copy, which is mapped back to
 // the Block at the end of compute()
 if( f_copy->is_feasible( true ) ) {
  delete f_sol_copy;
  f_sol_copy = f_copy->get_Solution( nullptr , false );
  return( true );
  }

 for( Block::Index i = 0 ; i < n ; ++i )
  static_cast< ColVariable * >( f_fun->get_active_var( i ) )->set_value(
								      old[ i ] );
 return( false );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::update_y( void )
{
 // the gradient with respect to y, normalised, and the largest step along
 // it, among 0, 0.01, ..., 1, that keeps y in [ 0 , 1 ], the components at
 // a bound of y and going out of it being left out
 std::vector< Function::FunctionValue > g;
 f_fun->get_y_gradient( g );
 Function::FunctionValue nrm = 0;
 for( auto & gk : g ) {
  if( ! std::isfinite( gk ) )
   gk = 0;
  nrm += gk * gk;
  }
 nrm = std::sqrt( nrm );
 auto y = f_fun->get_y();
 std::vector< Function::FunctionValue > d( g.size() , 0 );
 for( Block::Index k = 0 ; k < g.size() ; ++k ) {
  const auto gk = nrm > 0 ? g[ k ] / nrm : 0;
  const bool out = ( ( y[ k ] <= 0 ) && ( gk < 0 ) ) ||
		   ( ( y[ k ] >= 1 ) && ( gk > 0 ) );
  if( ! out )
   d[ k ] = gk;
  }
 double step = 0;
 for( int s = 100 ; s > 0 ; --s ) {
  const double t = s / 100.0;
  bool inside = true;
  for( Block::Index k = 0 ; inside && ( k < y.size() ) ; ++k )
   inside = ( y[ k ] + t * d[ k ] >= 0 ) && ( y[ k ] + t * d[ k ] <= 1 );
  if( inside ) {
   step = t;
   break;
   }
  }

 // the rows with a large slack move along the gradient, the active ones
 // with a fractional Variable are no longer looked at
 const auto & rows = f_fun->get_rows();
 for( Block::Index k = 0 ; k < y.size() ; ++k ) {
  const auto s = f_fun->slack( k );
  if( s > f_y_level )
   y[ k ] += step * d[ k ];
  else
   if( std::abs( s ) <= 1e-9 ) {
    bool fractional = false;
    for( auto i : rows[ k ].idx ) {
     auto x = static_cast< ColVariable * >( f_fun->get_active_var( i ) );
     if( x->is_integer() &&
	 ( std::abs( x->get_value() - std::round( x->get_value() ) ) > 1e-9 ) ) {
      fractional = true;
      break;
      }
     }
    if( fractional )
     y[ k ] = 0;
    }
  }
 f_fun->set_y( std::move( y ) );
 }

/*--------------------------------------------------------------------------*/

int IntegralityBarrierSolver::compute( bool changedvars )
{
 const auto start = std::chrono::steady_clock::now();
 auto elapsed = [ & ]() {
  return( std::chrono::duration< double >(
		     std::chrono::steady_clock::now() - start ).count() ); };

 if( ! f_Block )
  throw( std::logic_error( "IntegralityBarrierSolver::compute: no Block" ) );

 if( f_rebuild ) {  // the Block has changed altogether
  unbuild();
  f_rebuild = false;
  }

 delete f_sol;
 f_sol = nullptr;
 delete f_sol_copy;
 f_sol_copy = nullptr;
 f_ub = Inf< OFValue >();
 f_runs = 0;

 int status = kStopIter;
 try {
  // the starting point, in the Variable of the Block
  if( f_init_slvr >= 0 ) {
   const auto & slvrs = f_Block->get_registered_solvers();
   if( f_init_slvr >= int( slvrs.size() ) )
    throw( std::invalid_argument( "IntegralityBarrierSolver::compute: no "
				  "Solver at position intInitSlvr" ) );
   auto s = *std::next( slvrs.begin() , f_init_slvr );
   const int st = s->compute( false );
   if( ( st != kOK ) || ( ! s->has_var_solution() ) )
    throw( std::runtime_error( "IntegralityBarrierSolver::compute: the "
			       "Solver of the starting point has found none" ) );
   s->get_var_solution();
   }

  build();

  // where each run of Frank-Wolfe starts from: x0, written into the copy at
  // the start of each run, or the vertex of the oracle
  const bool init_x0 = ( f_fw->get_int_par( FrankWolfeSolver::intInitPoint )
			 == FrankWolfeSolver::eInitBlock );
  if( init_x0 || ( f_init_slvr >= 0 ) )  // the values, into the copy
   f_Block->map_forward_solution( f_copy );
  Solution * x0 = f_copy->get_Solution( nullptr , false );

  // Frank-Wolfe stops as soon as the rounding of the vertex the oracle has
  // left in the Variable, or of the iterate, is feasible
  bool found = false;
  auto fw = static_cast< FrankWolfeSolver * >( f_fw );
  auto stop = [ this , fw , & found ]( void ) {
   if( rounding_feasible() ) {
    found = true;
    return( int( eStopOK ) );
    }
   fw->get_var_solution();
   if( rounding_feasible() ) {
    found = true;
    return( int( eStopOK ) );
    }
   return( int( eContinue ) );
   };
  auto id = f_fw->set_event_handler( eEverykIteration , stop );

  for( f_runs = 0 ; f_runs < f_max_iter ; ++f_runs ) {
   if( elapsed() >= f_max_time ) {
    status = kStopTime;
    break;
    }
   if( init_x0 )
    x0->write( f_copy );
   // the time of the run, what is left of that of compute() at most
   f_fw->set_par( dblMaxTime , std::min( f_fw_max_time ,
				 std::max( f_max_time - elapsed() , 0.0 ) ) );
   if( f_fw->compute( false ) == kInfeasible ) {  // the linear relaxation is
    status = kInfeasible;                         // empty, so is the program
    ++f_runs;
    break;
    }
   if( found || rounding_feasible() ) {
    status = kOK;
    ++f_runs;
    break;
    }
   if( f_log && ( f_log_verb > 0 ) )
    *f_log << "IntegralityBarrierSolver: run " << f_runs << ", value "
	   << f_fw->get_ub() << std::endl;

   // y changed at the iterate Frank-Wolfe ends at, x0 moved towards it
   f_fw->get_var_solution();
   update_y();
   if( f_x0_step > 0 ) {
    Solution * x = f_copy->get_Solution( nullptr , false );
    Solution * nx = x0->scale( 1 - f_x0_step );
    nx->sum( x , f_x0_step );
    delete x;
    delete x0;
    x0 = nx;
    }
   }

  f_fw->reset_event_handler( eEverykIteration , id );
  delete x0;

  // the point found back into the Block, kept only if it is feasible
  // there too, the copy being a relaxation if mirror() could not reproduce
  // everything; its value is that of the Objective of the copy, which has
  // the terms of the whole tree [see fold_objectives()] and, Frank-Wolfe
  // being over, its own coefficients again
  if( f_sol_copy ) {
   f_sol_copy->write( f_copy );
   f_Block->map_back_solution( f_copy );
   if( f_Block->is_feasible( true ) ) {
    f_sol = f_Block->get_Solution( nullptr , false );
    auto obj = dynamic_cast< RealObjective * >( f_copy->get_objective() );
    if( obj ) {
     obj->compute();
     f_ub = obj->value();
     }
    else
     f_ub = 0;
    }
   else
    status = kStopIter;
   }
  }
 catch( ... ) {
  unbuild();
  throw;
  }

 return( status );
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::add_Modification( sp_Mod & mod )
{
 // the copy no longer is the Block: it is made again at the next compute()
 f_rebuild = true;
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::get_var_solution( Configuration * solc )
{
 if( ! f_sol )
  throw( std::logic_error( "IntegralityBarrierSolver::get_var_solution: no "
			   "feasible point has been found" ) );
 f_sol->write( f_Block );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR HANDLING THE PARAMETERS ----------------*/
/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intInitSlvr ):
   if( value < -1 )
    throw( std::invalid_argument( "IntegralityBarrierSolver::set_par: "
				  "intInitSlvr is -1 or the position of a "
				  "Solver; starting from the current values "
				  "is intInitPoint eInitBlock of the "
				  "FrankWolfeSolver [strFWCfg]" ) );
   f_init_slvr = value;
   return;
  case( intPsi ):
   f_psi = value;
   if( f_fun )
    f_fun->set_psi( value );
   return;
  case( intMaxIter ):  f_max_iter = value;  return;
  case( intLogVerb ):  f_log_verb = value;  return;
  default:             Solver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ):    f_max_time = value; return;
  case( dblYLevel ):     f_y_level = value;  return;
  case( dblX0Step ):     f_x0_step = value;  return;
  case( dblBarrierEps ):
   if( ! ( value > 0 ) )
    throw( std::invalid_argument( "IntegralityBarrierSolver::set_par: "
				  "dblBarrierEps must be positive" ) );
   f_eps = value;
   if( f_fun )
    f_fun->set_eps( value );
   return;
  default:               Solver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void IntegralityBarrierSolver::set_par( idx_type par , std::string && value )
{
 if( ( par >= strFWCfg ) && ( par < strLastAlgParIBSlv ) ) {
  if( f_father )
   throw( std::logic_error( "IntegralityBarrierSolver::set_par: " +
			    str_names[ par - strFWCfg ] + " cannot change "
			    "after the first compute()" ) );
  ( par == strFWCfg ? f_fw_cfg : f_lmo_bscfg ) = std::move( value );
  return;
  }
 Solver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

int IntegralityBarrierSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intInitSlvr ): return( -1 );
  case( intPsi ):      return( 0 );
  case( intMaxIter ):  return( 10 );
  default:             return( Solver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double IntegralityBarrierSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblYLevel ):     return( 0.9 );
  case( dblX0Step ):     return( 0 );
  case( dblBarrierEps ): return( 1e-8 );
  case( dblRelAcc ):     return( Inf< double >() );  // a heuristic promises
                                                     // no accuracy
  default:               return( Solver::get_dflt_dbl_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::string & IntegralityBarrierSolver::get_dflt_str_par( idx_type par )
 const
{
 static const std::string empty;
 if( ( par >= strFWCfg ) && ( par < strLastAlgParIBSlv ) )
  return( empty );
 return( Solver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

int IntegralityBarrierSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intInitSlvr ): return( f_init_slvr );
  case( intPsi ):      return( f_psi );
  case( intMaxIter ):  return( f_max_iter );
  case( intLogVerb ):  return( f_log_verb );
  default:             return( Solver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double IntegralityBarrierSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ):    return( f_max_time );
  case( dblYLevel ):     return( f_y_level );
  case( dblX0Step ):     return( f_x0_step );
  case( dblBarrierEps ): return( f_eps );
  case( dblRelAcc ):     return( Inf< double >() );
  default:               return( Solver::get_dbl_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::string & IntegralityBarrierSolver::get_str_par( idx_type par )
 const
{
 switch( par ) {
  case( strFWCfg ):    return( f_fw_cfg );
  case( strLMOBSCfg ): return( f_lmo_bscfg );
  default:             return( Solver::get_str_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type IntegralityBarrierSolver::int_par_str2idx(
					     const std::string & name ) const
{
 for( idx_type i = 0 ; i < int_names.size() ; ++i )
  if( name == int_names[ i ] )
   return( intInitSlvr + i );
 return( Solver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & IntegralityBarrierSolver::int_par_idx2str( idx_type idx )
 const
{
 if( ( idx >= intInitSlvr ) && ( idx < intLastAlgParIBSlv ) )
  return( int_names[ idx - intInitSlvr ] );
 return( Solver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type IntegralityBarrierSolver::dbl_par_str2idx(
					     const std::string & name ) const
{
 for( idx_type i = 0 ; i < dbl_names.size() ; ++i )
  if( name == dbl_names[ i ] )
   return( dblYLevel + i );
 return( Solver::dbl_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & IntegralityBarrierSolver::dbl_par_idx2str( idx_type idx )
 const
{
 if( ( idx >= dblYLevel ) && ( idx < dblLastAlgParIBSlv ) )
  return( dbl_names[ idx - dblYLevel ] );
 return( Solver::dbl_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type IntegralityBarrierSolver::str_par_str2idx(
					     const std::string & name ) const
{
 for( idx_type i = 0 ; i < str_names.size() ; ++i )
  if( name == str_names[ i ] )
   return( strFWCfg + i );
 return( Solver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & IntegralityBarrierSolver::str_par_idx2str( idx_type idx )
 const
{
 if( ( idx >= strFWCfg ) && ( idx < strLastAlgParIBSlv ) )
  return( str_names[ idx - strFWCfg ] );
 return( Solver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------- End File IntegralityBarrierSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
