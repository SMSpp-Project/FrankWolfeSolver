/*--------------------------------------------------------------------------*/
/*---------------------- File FrankWolfeSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the FrankWolfeSolver class.
 *
 * v1 status: structure (set_Block), parameter system, dual-solution
 * forwarding and result accessors are implemented; the vanilla Frank-Wolfe
 * loop in compute() is not yet implemented.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FrankWolfeSolver.h"

#include "FRealObjective.h"

#include "LinearFunction.h"

#include "DQuadFunction.h"

#include "QuadFunction.h"

#include "ColVariable.h"

#include <algorithm>

#include <atomic>

#include <chrono>

#include <cmath>

#include <thread>

#include <unordered_map>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using FunctionValue = Function::FunctionValue;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( FrankWolfeSolver );

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETER NAME TABLES --------------------------*/
/*--------------------------------------------------------------------------*/

static const std::vector< std::string > FWSlv_int_pars_str = {
 "intLMOObj" , "intLineSearch" , "intLMOSlvr" , "intAlgorithm" , "intMaxAtoms" ,
 "intCvxComb"
 };

/*--------------------------------------------------------------------------*/
/*---------------------------- AUXILIARY ROUTINES --------------------------*/
/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::set_default_parameters( void )
{
 f_lmo_obj     = get_dflt_int_par( intLMOObj );
 f_line_search = get_dflt_int_par( intLineSearch );
 f_lmo_slvr    = get_dflt_int_par( intLMOSlvr );
 f_algorithm   = get_dflt_int_par( intAlgorithm );
 f_max_atoms   = get_dflt_int_par( intMaxAtoms );
 f_cvx_comb    = get_dflt_int_par( intCvxComb );
 f_max_thread  = get_dflt_int_par( intMaxThread );
 f_max_iter    = get_dflt_int_par( intMaxIter );
 f_max_time    = get_dflt_dbl_par( dblMaxTime );
 f_rel_acc     = get_dflt_dbl_par( dblRelAcc );
 f_abs_acc     = get_dflt_dbl_par( dblAbsAcc );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- cleanup --------------------------------*/
/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::cleanup( void )
{
 // if the sub-Block objectives were modified, restore the original linear
 // coefficients so that the Block is left pristine

 // restore the LMO :Solver identities that were lent in acquire_LMOs
 for( auto & d : v_sb )
  if( d.lmo )
   d.lmo->set_id();

 // The restore is issued with eNoMod: the physical (abstract) coefficients are
 // set back to c0, but NO Modification is emitted/propagated. cleanup() runs
 // when FrankWolfeSolver detaches from the Block (typically at teardown), when
 // the other registered :Solver may already be gone; propagating an objective
 // change to them at that point is both useless and unsafe (it would notify a
 // possibly dangling observer). NOTE: for a sub-Block whose abstract objective
 // change is bridged to physical data (e.g. ThermalUnitBlock's per-period
 // costs), the physical data is therefore NOT restored on detach -- re-attach
 // (set_Block) re-snapshots c0, so this only matters if the bare sub-Block is
 // reused by a physical-data solver right after detaching, which the teardown
 // path does not do.
 const bool alpha = ( f_lmo_obj == LMOFull );
 if( f_modified )
  for( auto & d : v_sb ) {
   if( d.c0.empty() )
    continue;
   Index n = Index( d.c0.size() );

   // restore exactly what scatter() touched: when LMOFull updated only a
   // subset (the father-touched variables), restore that subset to c0 and
   // leave the rest (which scatter never changed) alone -- mirroring the
   // Subset path of scatter().
   if( alpha && ( d.obj_idx.size() < std::size_t( n ) ) ) {
    std::vector< Index > srt( d.obj_idx.begin() , d.obj_idx.end() );
    std::sort( srt.begin() , srt.end() );
    Function::Subset nms( srt.begin() , srt.end() );
    Function::Vec_FunctionValue nc( srt.size() );
    for( std::size_t k = 0 ; k < srt.size() ; ++k )
     nc[ k ] = d.c0[ srt[ k ] ];
    if( d.dq )
     d.dq->modify_linear_coefficients( std::move( nc ) , std::move( nms ) ,
                                       true , eNoMod );
    else
     d.lin->modify_coefficients( std::move( nc ) , std::move( nms ) , true ,
                                 eNoMod );
    continue;
    }

   Function::Vec_FunctionValue nc( d.c0 );
   if( d.dq )
    d.dq->modify_linear_coefficients( std::move( nc ) , Function::Range( 0 ,
                                       n ) , eNoMod );
   else
    d.lin->modify_coefficients( std::move( nc ) , Function::Range( 0 , n ) ,
                                eNoMod );
   }

 v_sb.clear();
 f_grad.clear();
 f_xval.clear();
 f_vval.clear();
 f_father_diag.clear();
 f_father_offdiag.clear();
 f_obj = nullptr;
 f_fun = nullptr;
 f_dq_father = nullptr;
 f_quad_father = nullptr;
 f_nsb = 0;
 f_max = false;
 f_modified = false;

 inhibit_Modification( false );

 clear_active_set();

 delete f_x;
 f_x = nullptr;
 f_value = 0;
 f_bound = 0;
 f_has_sol = false;
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::clear_active_set( void )
{
 for( auto & el : f_aset )
  delete el.f_sol;
 f_aset.clear();
 }

/*--------------------------------------------------------------------------*/

FrankWolfeSolver::Index FrankWolfeSolver::find_atom(
                  const std::vector< Function::FunctionValue > & val ) const
{
 // vertices coming from the (deterministic) LMO are bit-identical when equal;
 // a small relative tolerance absorbs any harmless numerical noise
 for( Index i = 0 ; i < f_aset.size() ; ++i ) {
  bool same = true;
  for( Index p = 0 ; p < val.size() ; ++p )
   if( std::abs( val[ p ] - f_aset[ i ].f_val[ p ] ) >
       1e-9 * ( 1 + std::abs( f_aset[ i ].f_val[ p ] ) ) ) {
    same = false;
    break;
    }
  if( same )
   return( i );
  }
 return( Index( f_aset.size() ) );
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::bound_active_set( void )
{
 // keep the active set within intMaxAtoms by merging the least-active atoms
 // into a single aggregate atom (a convex combination that preserves x)

 if( ( f_max_atoms <= 0 ) || ( Index( f_aset.size() ) <= Index( f_max_atoms ) ) )
  return;

 const Index G = Index( f_grad.size() );
 Index n_evict = Index( f_aset.size() ) - Index( f_max_atoms ) + 1;

 // pick the n_evict least-active atoms: smallest f_count, ties by smallest
 // f_weight
 std::vector< Index > idx( f_aset.size() );
 for( Index i = 0 ; i < idx.size() ; ++i )
  idx[ i ] = i;
 std::partial_sort( idx.begin() , idx.begin() + n_evict , idx.end() ,
   [ this ]( Index a , Index b ) {
    if( f_aset[ a ].f_count != f_aset[ b ].f_count )
     return( f_aset[ a ].f_count < f_aset[ b ].f_count );
    return( f_aset[ a ].f_weight < f_aset[ b ].f_weight ); } );

 // build the aggregate: a_bar = sum_{i in E} ( w_i / W ) atom_i , W = sum w_i
 double W = 0;
 for( Index k = 0 ; k < n_evict ; ++k )
  W += f_aset[ idx[ k ] ].f_weight;

 Atom agg;
 agg.f_weight = W;
 agg.f_count = 0;
 agg.f_val.assign( G , 0 );
 for( Index k = 0 ; k < n_evict ; ++k ) {
  Atom & e = f_aset[ idx[ k ] ];
  double w = ( W > 0 ? e.f_weight / W : 1.0 / n_evict );
  if( ! agg.f_sol )
   agg.f_sol = e.f_sol->scale( w );
  else
   agg.f_sol->sum( e.f_sol , w );
  for( Index p = 0 ; p < G ; ++p )
   agg.f_val[ p ] += w * e.f_val[ p ];
  agg.f_ci += w * e.f_ci;       // <c, .> is linear in the atom
  }

 // erase the evicted atoms (largest index first to keep indices valid),
 // deleting their Solution
 std::vector< Index > ev( idx.begin() , idx.begin() + n_evict );
 std::sort( ev.begin() , ev.end() , std::greater< Index >() );
 for( Index e : ev ) {
  delete f_aset[ e ].f_sol;
  f_aset.erase( f_aset.begin() + e );
  }

 f_aset.push_back( std::move( agg ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- set_Block ------------------------------*/
/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::set_Block( Block * block )
{
 if( block == f_Block )  // calling set_Block() twice on the same Block
  return;                // is explicitly allowed: do nothing

 if( f_Block )           // detaching from a previous Block
  cleanup();

 Solver::set_Block( block );  // this sets f_Block

 if( ! block )           // just detaching
  return;

 // the father Block must have no Variable and no Constraint of its own- - - -

 if( ( ! block->get_static_variables().empty() ) ||
     ( ! block->get_dynamic_variables().empty() ) )
  throw( std::invalid_argument(
   "FrankWolfeSolver: the father Block must have no Variable of its own" ) );

 if( ( ! block->get_static_constraints().empty() ) ||
     ( ! block->get_dynamic_constraints().empty() ) )
  throw( std::invalid_argument(
   "FrankWolfeSolver: the father Block must have no Constraint of its own" ) );

 // the father Objective must be a FRealObjective with a C05Function- - - - -

 f_obj = dynamic_cast< FRealObjective * >( block->get_objective() );
 if( ! f_obj )
  throw( std::invalid_argument(
   "FrankWolfeSolver: the father Block must have a FRealObjective" ) );

 f_fun = dynamic_cast< C05Function * >( f_obj->get_function() );
 if( ! f_fun )
  throw( std::invalid_argument(
   "FrankWolfeSolver: the father Objective Function must be a C05Function" ) );

 f_max = ( f_obj->get_sense() == Objective::eMax );

 // there must be at least one sub-Block - - - - - - - - - - - - - - - - - - -

 const auto & sb = block->get_nested_Blocks();
 f_nsb = Index( sb.size() );
 if( ! f_nsb )
  throw( std::invalid_argument(
   "FrankWolfeSolver: the father Block has no sub-Block" ) );

 // scan the sub-Block: validate the objectives, snapshot the original linear
 // coefficients, and record the position of each variable in its objective - -

 v_sb.resize( f_nsb );
 std::unordered_map< Variable * , std::pair< Index , Index > > var2pos;

 for( Index j = 0 ; j < f_nsb ; ++j ) {
  auto & d = v_sb[ j ];
  d.block = sb[ j ];

  d.obj = dynamic_cast< FRealObjective * >( d.block->get_objective() );
  if( ! d.obj )
   throw( std::invalid_argument(
    "FrankWolfeSolver: a sub-Block has no FRealObjective" ) );

  d.fun = d.obj->get_function();
  d.dq  = dynamic_cast< DQuadFunction * >( d.fun );   // covers DQuad and Quad
  d.lin = d.dq ? nullptr : dynamic_cast< LinearFunction * >( d.fun );
  if( ( ! d.dq ) && ( ! d.lin ) )
   throw( std::invalid_argument( "FrankWolfeSolver: a sub-Block Objective "
    "Function is not a Linear/DQuad/QuadFunction" ) );

  if( ( d.obj->get_sense() == Objective::eMax ) != f_max )
   throw( std::invalid_argument( "FrankWolfeSolver: a sub-Block Objective "
    "sense differs from the father one" ) );

  Index n = d.fun->get_num_active_var();
  d.c0.resize( n );
  for( Index i = 0 ; i < n ; ++i ) {
   d.c0[ i ] = d.dq ? d.dq->get_linear_coefficient( i )
                    : d.lin->get_coefficient( i );
   var2pos[ d.fun->get_active_var( i ) ] = { j , i };
   }
  }

 // build the gradient-to-sub-Block scatter map: every active variable of the
 // father Objective must be active in (exactly) one sub-Block Objective - - -

 Index G = f_fun->get_num_active_var();
 f_grad.resize( G );

 for( Index p = 0 ; p < G ; ++p ) {
  auto it = var2pos.find( f_fun->get_active_var( p ) );
  if( it == var2pos.end() )
   throw( std::invalid_argument( "FrankWolfeSolver: a father-Objective "
    "variable is not active in any sub-Block Objective" ) );
  auto & d = v_sb[ it->second.first ];
  d.grad_idx.push_back( p );
  d.obj_idx.push_back( it->second.second );
  }

 f_xval.resize( G );
 f_vval.resize( G );

 // the exact line search is available when the father Objective is quadratic
 // (DQuadFunction or QuadFunction): cache its (static) quadratic structure,
 // i.e. the diagonal a_p and, for a QuadFunction, the off-diagonal terms

 f_dq_father = dynamic_cast< DQuadFunction * >( f_fun );  // DQuad or Quad
 f_quad_father = dynamic_cast< QuadFunction * >( f_fun ); // Quad only

 if( f_dq_father ) {
  f_father_diag.resize( G );
  for( Index p = 0 ; p < G ; ++p )
   f_father_diag[ p ] = f_dq_father->get_quadratic_coefficient( p );

  if( f_quad_father ) {
   auto mat = f_quad_father->get_matrix();   // off-diagonal (lower half)
   for( int k = 0 ; k < mat.outerSize() ; ++k )
    for( QuadFunction::Qmat::InnerIterator it( mat , k ) ; it ; ++it )
     if( it.row() != it.col() )
      f_father_offdiag.emplace_back( Index( it.row() ) , Index( it.col() ) ,
                                     it.value() );
   }
  }

 // v1 assumes a static problem: ignore all Modification

 inhibit_Modification( true );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETER HANDLING -----------------------------*/
/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( intLMOObj ):     f_lmo_obj = value;     return;
  case( intLineSearch ): f_line_search = value; return;
  case( intLMOSlvr ):    f_lmo_slvr = value;    return;
  case( intAlgorithm ):  f_algorithm = value;   return;
  case( intMaxAtoms ):   f_max_atoms = value;   return;
  case( intCvxComb ):    f_cvx_comb = value;    return;
  case( intMaxThread ):         f_max_thread = value;  return;
  case( intMaxIter ):           f_max_iter = value;    return;
  default:                      CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dblMaxTime ): f_max_time = value; return;
  case( dblRelAcc ):  f_rel_acc = value;  return;
  case( dblAbsAcc ):  f_abs_acc = value;  return;
  default:            CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

int FrankWolfeSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( intLMOObj ):     return( LMOLinear );
  case( intLineSearch ): return( LSAuto );
  case( intLMOSlvr ):    return( 0 );
  case( intAlgorithm ):  return( AlgVanilla );
  case( intMaxAtoms ):   return( 0 );
  case( intCvxComb ):    return( eObjCvxComb );
  default:                      return( CDASolver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

int FrankWolfeSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( intLMOObj ):     return( f_lmo_obj );
  case( intLineSearch ): return( f_line_search );
  case( intLMOSlvr ):    return( f_lmo_slvr );
  case( intAlgorithm ):  return( f_algorithm );
  case( intMaxAtoms ):   return( f_max_atoms );
  case( intCvxComb ):    return( f_cvx_comb );
  case( intMaxThread ):         return( f_max_thread );
  case( intMaxIter ):           return( f_max_iter );
  default:                      return( CDASolver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double FrankWolfeSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dblMaxTime ): return( f_max_time );
  case( dblRelAcc ):  return( f_rel_acc );
  case( dblAbsAcc ):  return( f_abs_acc );
  default:            return( CDASolver::get_dbl_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type FrankWolfeSolver::int_par_str2idx( const std::string & name )
 const
{
 for( idx_type i = 0 ; i < FWSlv_int_pars_str.size() ; ++i )
  if( name == FWSlv_int_pars_str[ i ] )
   return( intLastParCDAS + i );

 return( CDASolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FrankWolfeSolver::int_par_idx2str( idx_type idx ) const
{
 if( ( idx >= intLastParCDAS ) && ( idx < intLastParFWSlv ) )
  return( FWSlv_int_pars_str[ idx - intLastParCDAS ] );

 return( CDASolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ SOLVING THE PROBLEM -----------------------------*/
/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::acquire_LMOs( void )
{
 for( Index j = 0 ; j < f_nsb ; ++j ) {
  auto & d = v_sb[ j ];
  const auto & rs = d.block->get_registered_solvers();
  if( int( rs.size() ) <= f_lmo_slvr )
   throw( std::logic_error( "FrankWolfeSolver: sub-Block " +
    std::to_string( j ) + " has no :Solver at index " +
    std::to_string( f_lmo_slvr ) + " to use as LMO" ) );
  auto it = rs.begin();
  std::advance( it , f_lmo_slvr );
  d.lmo = *it;

  // lend FrankWolfeSolver's identity to the LMO :Solver, so that it can
  // lock its sub-Block (which is owned by f_id via the father Block lock)
  d.lmo->set_id( f_id );
  }
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::evaluate_gradient( void )
{
 // evaluate the father Objective at the current point and extract its
 // diagonal linearization (the gradient); the father is assumed never to
 // produce vertical linearizations

 f_fun->compute( true );
 if( ! f_fun->has_linearization( true ) )
  f_fun->compute_new_linearization( true );
 f_fun->get_linearization_coefficients( f_grad.data() );
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::scatter( void )
{
 // write alpha * c_j + g_j into the linear coefficients of each sub-Block
 // Objective; the quadratic part (if any) is left untouched

 const bool alpha = ( f_lmo_obj == LMOFull );

 for( auto & d : v_sb ) {
  Index n = Index( d.c0.size() );

  // LMOFull with the father touching only a proper subset of the sub-Block
  // Objective variables: update *only* those coefficients (a Subset
  // Modification), leaving the father-irrelevant ones at their original value
  // c0. This is both cheaper and necessary when the sub-Block Objective carries
  // extra variables the father does not reach (e.g. formulation auxiliaries):
  // a full-range change would also rewrite those, and some Block reject (or
  // mis-map) a change to them when translating the abstract Objective change
  // into their physical data. The subset is passed sorted (ordered = true), as
  // required by such abstract-to-physical translators.
  if( alpha && ( d.obj_idx.size() < std::size_t( n ) ) ) {
   std::vector< std::pair< Index , FunctionValue > > upd( d.obj_idx.size() );
   for( Index k = 0 ; k < d.obj_idx.size() ; ++k )
    upd[ k ] = { d.obj_idx[ k ] ,
                 d.c0[ d.obj_idx[ k ] ] + f_grad[ d.grad_idx[ k ] ] };
   std::sort( upd.begin() , upd.end() );
   Function::Subset nms( upd.size() );
   Function::Vec_FunctionValue nc( upd.size() );
   for( std::size_t k = 0 ; k < upd.size() ; ++k ) {
    nms[ k ] = upd[ k ].first;
    nc[ k ]  = upd[ k ].second;
    }
   if( d.dq )
    d.dq->modify_linear_coefficients( std::move( nc ) , std::move( nms ) ,
                                      true );
   else
    d.lin->modify_coefficients( std::move( nc ) , std::move( nms ) , true );
   continue;
   }

  Function::Vec_FunctionValue nc( n , 0 );
  if( alpha )
   for( Index i = 0 ; i < n ; ++i )
    nc[ i ] = d.c0[ i ];
  for( Index k = 0 ; k < d.obj_idx.size() ; ++k )
   nc[ d.obj_idx[ k ] ] += f_grad[ d.grad_idx[ k ] ];

  if( d.dq )
   d.dq->modify_linear_coefficients( std::move( nc ) ,
                                     Function::Range( 0 , n ) );
  else
   d.lin->modify_coefficients( std::move( nc ) , Function::Range( 0 , n ) );
  }

 f_modified = true;
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::run_LMOs( bool changedvars )
{
 // run one LMO ( solve + materialize the vertex into the sub-Block variables )
 auto run_one = [ changedvars ]( SubBlockData & d ) {
  try {
   d.status = d.lmo->compute( changedvars );
   d.value  = d.lmo->get_var_value();   // M_j(v_j)
   d.lmo->get_var_solution();           // materialize v_j into the variables
   }
  catch( ... ) { d.excp = std::current_exception(); }
  };

 // sequential path (intMaxThread <= 1 or a single sub-Block)
 if( ( f_max_thread <= 1 ) || ( f_nsb <= 1 ) ) {
  for( auto & d : v_sb )
   run_one( d );
  return;
  }

 // parallel path: each LMO runs on its own (distinct) sub-Block, with the
 // FrankWolfeSolver identity lent to its :Solver (so lock( f_id ) on the
 // already-f_id-owned sub-Block never contends); load-balanced via an atomic
 // counter, the main thread participating. No shared mutable state: each LMO
 // writes only into its own sub-Block and its own SubBlockData slot.

 std::atomic< Index > next( 0 );
 auto chunk = [ & ]() {
  for( ; ; ) {
   Index i = next.fetch_add( 1 );
   if( i >= f_nsb )
    break;
   run_one( v_sb[ i ] );
   }
  };

 int nthreads = std::min< int >( f_max_thread , int( f_nsb ) );
 std::vector< std::thread > pool;
 pool.reserve( nthreads - 1 );
 for( int t = 1 ; t < nthreads ; ++t )
  pool.emplace_back( chunk );
 chunk();                                // the main thread participates
 for( auto & th : pool )
  th.join();

 // re-throw (in the main thread) the first exception, if any
 for( auto & d : v_sb )
  if( d.excp ) {
   auto e = d.excp;
   d.excp = nullptr;
   std::rethrow_exception( e );
   }
 }

/*--------------------------------------------------------------------------*/

int FrankWolfeSolver::compute( bool changedvars )
{
 if( ! f_Block )
  throw( std::logic_error( "FrankWolfeSolver::compute: no Block registered" ) );

 // lock the father Block (LagrangianDualSolver pattern) - - - - - - - - - - -

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  return( kBlockLocked );

 acquire_LMOs();

 // LMOLinear requires purely linear sub-Block objectives (no quadratic term
 // to keep in the oracle); use LMOQuad/LMOFull otherwise

 if( f_lmo_obj == LMOLinear )
  for( auto & d : v_sb )
   if( d.dq )
    throw( std::logic_error( "FrankWolfeSolver: LMOLinear requires "
     "LinearFunction sub-Block objectives; use LMOQuad/LMOFull otherwise" ) );

 int status = ( f_algorithm == AlgVanilla ) ? compute_vanilla( changedvars )
                                            : compute_active_set( changedvars );

 if( ! owned )
  f_Block->unlock( f_id );

 return( status );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue FrankWolfeSolver::eval_modified_objective( void )
{
 // sum_j M_j evaluated at the current Block point ( = <grad F, point> for the
 // linearized total objective, the sub-Block objectives being the modified
 // M_j = <alpha c_j + g_j, .> [+ beta q_j] )

 OFValue s = 0;
 for( auto & d : v_sb ) {
  d.fun->compute( true );
  s += d.fun->get_value();
  }
 return( s );
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::capture_father_values(
                              std::vector< Function::FunctionValue > & dst )
{
 const Index G = Index( f_grad.size() );
 dst.resize( G );
 for( Index p = 0 ; p < G ; ++p )
  dst[ p ] = static_cast< ColVariable * >(
                            f_fun->get_active_var( p ) )->get_value();
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue FrankWolfeSolver::quad_form(
                       const std::vector< Function::FunctionValue > & a ,
                       const std::vector< Function::FunctionValue > & b )
{
 // 1/2 <d, A d> with d = b - a and A the father Hessian = sum_p a_p d_p^2 +
 // sum_{(r,c)} q_rc d_r d_c (off-diagonal, each pair once)

 const Index G = Index( f_grad.size() );
 OFValue Q = 0;
 for( Index p = 0 ; p < G ; ++p ) {
  OFValue dp = b[ p ] - a[ p ];
  Q += f_father_diag[ p ] * dp * dp;
  }
 for( const auto & [ r , c , q ] : f_father_offdiag )
  Q += q * ( b[ r ] - a[ r ] ) * ( b[ c ] - a[ c ] );
 return( Q );
 }

/*--------------------------------------------------------------------------*/

int FrankWolfeSolver::compute_vanilla( bool changedvars )
{
 bool all_lin = true;
 for( auto & d : v_sb )
  if( d.dq ) { all_lin = false; break; }

 // intCvxComb: eObjCvxComb solves the Dantzig-Wolfe / "disaggregated" problem
 // (P2) -- the sub-Block cost is the convex combination of the vertex costs,
 // tracked in cbar below -- while eObjAtX solves the genuine composite problem
 // (P1) -- the sub-Block cost is re-evaluated at the iterate (= mx_sum - gx).
 // The two coincide for linear sub-Block objectives. In (P2) the sub-Block term
 // is linear in the step, so the exact line search is available whenever the
 // father is quadratic (even with quadratic sub-Block objectives); in (P1) it
 // is (for now) available only for linear sub-Block objectives (else agnostic).
 const bool cvx = ( f_cvx_comb == eObjCvxComb );
 bool exact = f_dq_father && ( cvx || all_lin ) &&
              ( ( f_line_search == LSExact ) || ( f_line_search == LSAuto ) );

 auto t_start = std::chrono::steady_clock::now();
 auto elapsed = [ & ]() {
  return( std::chrono::duration< double >(
           std::chrono::steady_clock::now() - t_start ).count() ); };

 const Index G = Index( f_grad.size() );

 // <grad f_father(x), val> over the father active variables
 auto grad_dot = [ this , G ]( const std::vector< FunctionValue > & val ) {
  OFValue s = 0;
  for( Index p = 0 ; p < G ; ++p )
   s += f_grad[ p ] * val[ p ];
  return( s ); };

 // initialization: a first LMO at the current point gives x0 = v0; the convex
 // combination of sub-Block costs starts at the cost of that single vertex,
 // cbar = sum_j h_j(v0_j) = mv_init - <g,v0> (the x-independent part)
 evaluate_gradient();
 scatter();
 run_LMOs( true );
 OFValue mv_init = 0;
 for( auto & d : v_sb )
  mv_init += d.value;
 delete f_x;
 f_x = f_Block->get_Solution( nullptr , false );
 capture_father_values( f_vval );
 OFValue cbar = mv_init - grad_dot( f_vval );    // sum_j h_j(v0_j)
 f_has_sol = true;

 int status = kStopIter;

 for( int t = 0 ; ; ++t ) {
  if( t >= f_max_iter )         { status = kStopIter; break; }
  if( elapsed() >= f_max_time ) { status = kStopTime; break; }

  f_x->write( f_Block );
  evaluate_gradient();
  OFValue father_val = f_fun->get_value();
  scatter();

  OFValue mx_sum = eval_modified_objective();   // sum_j M_j(x_j)
  capture_father_values( f_xval );

  run_LMOs( true );
  OFValue mv_sum = 0;
  for( auto & d : v_sb )
   mv_sum += d.value;                            // sum_j M_j(v_j)
  capture_father_values( f_vval );

  OFValue gx = grad_dot( f_xval );               // <grad f_father(x), x>
  OFValue cv = mv_sum - grad_dot( f_vval );      // sum_j h_j(v_j)

  // cost_x = <grad f_father, x> + ( sub-Block cost at x ): in (P2) the latter
  // is the tracked convex combination cbar; in (P1) it is re-evaluated, i.e.
  // mx_sum - gx, so cost_x = mx_sum. Everything (value, gap, line search) is
  // then identical to the linear case with mx_sum replaced by cost_x.
  OFValue cost_x = cvx ? ( gx + cbar ) : mx_sum;

  f_value = father_val + ( cost_x - gx );
  OFValue gap = f_max ? ( mv_sum - cost_x ) : ( cost_x - mv_sum );
  f_bound = f_max ? ( f_value + gap ) : ( f_value - gap );

  OFValue rel_thr = f_rel_acc * std::max( OFValue( 1 ) , std::abs( f_value ) );
  if( ( gap <= rel_thr ) ||
      ( std::isfinite( f_abs_acc ) && ( gap <= f_abs_acc ) ) ) {
   status = kOK;
   break;
   }

  OFValue gamma;
  if( exact ) {
   OFValue gd = mv_sum - cost_x;                 // d/dgamma at gamma=0 ( <0 good )
   OFValue Q = quad_form( f_xval , f_vval );     // 1/2 <v-x, A(v-x)>
   if( std::abs( Q ) > 1e-12 )
    gamma = std::min( OFValue( 1 ) ,
                      std::max( OFValue( 0 ) , - gd / ( 2 * Q ) ) );
   else
    gamma = ( ( f_max && ( gd > 0 ) ) || ( ( ! f_max ) && ( gd < 0 ) ) )
            ? OFValue( 1 ) : OFValue( 0 );
   }
  else
   gamma = OFValue( 2 ) / OFValue( t + 2 );      // Agnostic open-loop rule

  // x <- ( 1 - gamma ) x + gamma v ; cbar <- ( 1 - gamma ) cbar + gamma cv
  Solution * v_sol = f_Block->get_Solution( nullptr , false );
  Solution * nx = f_x->scale( 1 - gamma );
  nx->sum( v_sol , gamma );
  delete v_sol;
  delete f_x;
  f_x = nx;
  cbar = ( 1 - gamma ) * cbar + gamma * cv;
  }

 f_x->write( f_Block );
 return( status );
 }

/*--------------------------------------------------------------------------*/

int FrankWolfeSolver::compute_active_set( bool changedvars )
{
 // Away-step Frank-Wolfe (AlgAwayStep) and Blended Pairwise / Pairwise CG
 // (AlgBPCG), both maintaining an active set x = sum_i lambda_i atom_i.

 const bool pairwise = ( f_algorithm == AlgBPCG );

 bool all_lin = true;
 for( auto & d : v_sb )
  if( d.dq ) { all_lin = false; break; }

 // see compute_vanilla for intCvxComb: eObjCvxComb (cvx) solves the
 // Dantzig-Wolfe problem (P2), the sub-Block cost being the convex combination
 // cbar = sum_i lambda_i c_i of the active-set atom costs (c_i = Atom::f_ci,
 // already cached); eObjAtX solves the composite problem (P1) (sub-Block cost
 // re-evaluated at x, = mx_sum - gx). cost_x = <grad f_father, x> + sub-cost
 // unifies the two, with cost_x replacing mx_sum throughout.
 const bool cvx = ( f_cvx_comb == eObjCvxComb );
 bool exact = f_dq_father && ( cvx || all_lin ) &&
              ( ( f_line_search == LSExact ) || ( f_line_search == LSAuto ) );

 auto t_start = std::chrono::steady_clock::now();
 auto elapsed = [ & ]() {
  return( std::chrono::duration< double >(
           std::chrono::steady_clock::now() - t_start ).count() ); };

 const Index G = Index( f_grad.size() );
 const double drop_eps = 1e-9;

 // <grad f_father(x), val> over the father active variables
 auto grad_dot = [ this , G ]( const std::vector< FunctionValue > & val ) {
  OFValue s = 0;
  for( Index p = 0 ; p < G ; ++p )
   s += f_grad[ p ] * val[ p ];
  return( s ); };

 // initialization: a first LMO gives x0 = v0; active set = { ( v0 , 1 ) }
 evaluate_gradient();
 scatter();
 run_LMOs( true );
 OFValue mv_init = 0;
 for( auto & d : v_sb )
  mv_init += d.value;
 clear_active_set();
 delete f_x;
 f_x = f_Block->get_Solution( nullptr , false );                 // x = v0
 { Atom a0;
   a0.f_sol = f_Block->get_Solution( nullptr , false );
   a0.f_weight = 1;
   a0.f_count = 1;
   capture_father_values( a0.f_val );
   a0.f_ci = mv_init - grad_dot( a0.f_val );    // <c, v0> ( x-independent )
   f_aset.push_back( std::move( a0 ) ); }
 f_has_sol = true;

 std::vector< FunctionValue > a_val;   // away-atom father values

 int status = kStopIter;

 for( int t = 0 ; ; ++t ) {
  if( t >= f_max_iter )         { status = kStopIter; break; }
  if( elapsed() >= f_max_time ) { status = kStopTime; break; }

  // gradient and linearization at x
  f_x->write( f_Block );
  evaluate_gradient();
  OFValue father_val = f_fun->get_value();
  scatter();
  OFValue mx_sum = eval_modified_objective();    // <grad F, x>
  capture_father_values( f_xval );

  // away vertex: the active-set atom a maximizing (minimizing, if eMax)
  // <grad F, a> — the "worst" atom we want to move weight away from. With
  // linear sub-Block objectives <grad F, a_i> = <grad f_father(x), a_i> + c_i
  // is a cheap cached dot product; otherwise the atom is written and evaluated.
  Index a_idx = 0;
  OFValue La_a = f_max ? Inf< OFValue >() : - Inf< OFValue >();
  for( Index i = 0 ; i < f_aset.size() ; ++i ) {
   OFValue La;
   if( all_lin )
    La = grad_dot( f_aset[ i ].f_val ) + f_aset[ i ].f_ci;
   else {
    f_aset[ i ].f_sol->write( f_Block );
    La = eval_modified_objective();
    }
   if( f_max ? ( La < La_a ) : ( La > La_a ) ) {
    La_a = La;
    a_idx = i;
    }
   }
  double lambda_a = f_aset[ a_idx ].f_weight;
  a_val = f_aset[ a_idx ].f_val;   // away-atom father values (stored)

  // FW vertex: the LMO of grad F
  run_LMOs( true );
  OFValue mv_sum = 0;
  for( auto & d : v_sb )
   mv_sum += d.value;                            // <grad F, v>
  capture_father_values( f_vval );

  // F(x) and the Frank-Wolfe gap (valid optimality bound). cost_x replaces
  // mx_sum: in (P2) the sub-Block cost at x is the convex combination cbar of
  // the atom costs ( cost_x = gx + cbar ); in (P1) it is mx_sum ( = gx + cost
  // re-evaluated at x ). The two coincide for linear sub-Block objectives.
  OFValue gx = grad_dot( f_xval );               // <grad f_father(x), x>
  OFValue cbar = 0;
  if( cvx )
   for( const auto & el : f_aset )
    cbar += el.f_weight * el.f_ci;               // sum_i lambda_i h(atom_i)
  OFValue cost_x = cvx ? ( gx + cbar ) : mx_sum;

  f_value = father_val + ( cost_x - gx );
  OFValue fw_gap   = f_max ? ( mv_sum - cost_x ) : ( cost_x - mv_sum );
  OFValue away_gap = f_max ? ( cost_x - La_a )   : ( La_a - cost_x );
  f_bound = f_max ? ( f_value + fw_gap ) : ( f_value - fw_gap );

  OFValue rel_thr = f_rel_acc * std::max( OFValue( 1 ) , std::abs( f_value ) );
  if( ( fw_gap <= rel_thr ) ||
      ( std::isfinite( f_abs_acc ) && ( fw_gap <= f_abs_acc ) ) ) {
   status = kOK;
   break;
   }

  // choose the step type and its data - - - - - - - - - - - - - - - - - - - -
  // step kinds: 0 = FW ( toward v ), 1 = away ( from a ), 2 = pairwise ( a->v )

  int kind;
  double gamma_max;
  OFValue gd;                  // <grad F, d>, d the (un-normalized) direction
  const std::vector< FunctionValue > * dfrom;  // d = (*dto) - (*dfrom)
  const std::vector< FunctionValue > * dto;

  if( pairwise ) {             // pairwise: move weight from a to v, d = a - v
   kind = 2;
   gamma_max = lambda_a;
   gd = La_a - mv_sum;         // <grad F, a - v>
   dfrom = &f_vval; dto = &a_val;
   }
  else                         // away-step: FW vs away
   if( ( fw_gap >= away_gap ) || ( lambda_a >= 1 - drop_eps ) ) {
    kind = 0;                  // FW step, d = x - v
    gamma_max = 1;
    gd = cost_x - mv_sum;      // d/dgamma along x - v ( cost_x for mx_sum )
    dfrom = &f_vval; dto = &f_xval;
    }
   else {
    kind = 1;                  // away step, d = a - x
    gamma_max = lambda_a / ( 1 - lambda_a );
    gd = La_a - cost_x;        // d/dgamma along a - x ( cost_x for mx_sum )
    dfrom = &f_xval; dto = &a_val;
    }

  // line search: x <- x - gamma d , gamma in [ 0 , gamma_max ] - - - - - - - -
  OFValue gamma;
  if( exact ) {
   OFValue Q = quad_form( *dfrom , *dto );        // 1/2 <d, A d>
   if( std::abs( Q ) > 1e-12 )
    gamma = std::min( OFValue( gamma_max ) ,
                      std::max( OFValue( 0 ) , gd / ( 2 * Q ) ) );
   else
    gamma = ( ( ! f_max && ( gd > 0 ) ) || ( f_max && ( gd < 0 ) ) )
            ? OFValue( gamma_max ) : OFValue( 0 );
   }
  else
   gamma = std::min( OFValue( gamma_max ) , OFValue( 2 ) / OFValue( t + 2 ) );

  // apply the step: update the iterate and the active set - - - - - - - - - -

  // lambda to add ( or merge, by dedup ) a new FW vertex v with weight w
  OFValue v_ci = mv_sum - grad_dot( f_vval );   // <c, v> ( x-independent )
  auto add_vertex = [ & ]( Solution * v_sol , double w ) {
   Index j = find_atom( f_vval );
   if( j < f_aset.size() ) { f_aset[ j ].f_weight += w; delete v_sol; }
   else { Atom a; a.f_sol = v_sol; a.f_weight = w; a.f_count = 0;
          a.f_val = f_vval; a.f_ci = v_ci; f_aset.push_back( std::move( a ) ); }
   };

  if( kind == 0 ) {             // FW step: x <- (1-gamma) x + gamma v
   Solution * v_sol = f_Block->get_Solution( nullptr , false );
   Solution * nx = f_x->scale( 1 - gamma );
   nx->sum( v_sol , gamma );
   delete f_x; f_x = nx;
   for( auto & el : f_aset )
    el.f_weight *= ( 1 - gamma );
   add_vertex( v_sol , gamma );
   }
  else
   if( kind == 1 ) {            // away step: x <- (1+gamma) x - gamma a
    Solution * a_sol = f_aset[ a_idx ].f_sol;
    Solution * nx = f_x->scale( 1 + gamma );
    nx->sum( a_sol , - gamma );
    delete f_x; f_x = nx;
    for( auto & el : f_aset )
     el.f_weight *= ( 1 + gamma );
    f_aset[ a_idx ].f_weight -= gamma;
    }
   else {                       // pairwise step: x <- x + gamma ( v - a )
    Solution * v_sol = f_Block->get_Solution( nullptr , false );
    Solution * a_sol = f_aset[ a_idx ].f_sol;
    Solution * nx = f_x->scale( 1.0 );            // copy of x
    nx->sum( a_sol , - gamma );
    nx->sum( v_sol , gamma );
    delete f_x; f_x = nx;
    f_aset[ a_idx ].f_weight -= gamma;
    add_vertex( v_sol , gamma );
    }

  // drop atoms whose weight has vanished
  for( Index i = 0 ; i < f_aset.size() ; )
   if( f_aset[ i ].f_weight <= drop_eps ) {
    delete f_aset[ i ].f_sol;
    f_aset.erase( f_aset.begin() + i );
    }
   else
    ++i;

  // age the surviving atoms, then bound the active set by aggregation
  for( auto & el : f_aset )
   ++el.f_count;
  bound_active_set();
  }

 f_x->write( f_Block );
 return( status );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- READING THE SOLUTION ------------------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue FrankWolfeSolver::get_var_value( void ) { return( f_value ); }

/*--------------------------------------------------------------------------*/

Solver::OFValue FrankWolfeSolver::get_lb( void )
{
 return( f_max ? f_value : f_bound );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue FrankWolfeSolver::get_ub( void )
{
 return( f_max ? f_bound : f_value );
 }

/*--------------------------------------------------------------------------*/

bool FrankWolfeSolver::has_var_solution( void ) { return( f_has_sol ); }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::get_var_solution( Configuration * solc )
{
 // write the current iterate x (the maintained convex combination) into the
 // father Block variables (which recurse into the sub-Block)

 if( f_x )
  f_x->write( f_Block );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- DUAL SOLUTION ----------------------------------*/
/*--------------------------------------------------------------------------*/

bool FrankWolfeSolver::has_dual_solution( void )
{
 // true iff every sub-Block LMO is a CDASolver with a dual solution available

 if( v_sb.empty() )
  return( false );

 for( auto & d : v_sb ) {
  auto cda = dynamic_cast< CDASolver * >( d.lmo );
  if( ( ! cda ) || ( ! cda->has_dual_solution() ) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/

void FrankWolfeSolver::get_dual_solution( Configuration * solc )
{
 // pure forwarding: at (exact) termination the LMO dual solutions are valid
 // Lagrangian multipliers for the overall problem

 for( auto & d : v_sb )
  if( auto cda = dynamic_cast< CDASolver * >( d.lmo ) )
   cda->get_dual_solution( solc );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File FrankWolfeSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
